//! Benchmark model (bench/pcb model.json with optional `cold_routing` rules) to routing objects and rules.

use crate::cdt::NONE;
use crate::geom::{convex_hull, p2, point_in_poly, poly_edge_dist, P2};
use serde_json::Value;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Class {
    pub name: String,
    pub width: f64,
    pub clear: f64,
    pub via_d: f64,
    pub via_drill: f64,
}

#[derive(Clone, Debug)]
pub struct NetInfo {
    pub name: String,
    pub class: usize,
    /// Index into `Board::layers` of the plane this net is poured on.
    pub plane: Option<usize>,
    /// Signal layers (bit per index into `Board::signal`) the net may use.
    pub allowed: u32,
    pub wpref: Vec<f64>,
    pub wmin: Vec<f64>,
    /// For a ribbon net: its two pair members (first, second).
    pub ribbon: Option<(u32, u32)>,
    /// Via (diameter, drill) overriding the class (ribbon nets use one wide via for the via pair).
    pub via_size: Option<(f64, f64)>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Pad,
    Via,
    Hole,
    Outline,
    Keepout,
}

#[derive(Clone, Debug)]
pub struct ObjSpec {
    pub kind: Kind,
    pub net: u32,
    pub layers: u32,
    /// Convex hull core points (pads), ring points (outline/keepout).
    pub pts: Vec<P2>,
    pub r: f64,
    pub clear: f64,
    pub center: P2,
    pub drill: f64,
    pub name: String,
    /// Pad present on a plane layer (through-hole pads).
    pub through: bool,
}

/// One pour of a plane net on a plane layer.
#[derive(Clone)]
pub struct PlaneRegion {
    pub net: u32,
    pub layer: usize,
    pub priority: i64,
    pub poly: Vec<P2>,
    pub holes: Vec<Vec<P2>>,
}

/// A neck-down stub: a straight track at a narrower width than the net's class allows elsewhere, from a pad
/// out to where the class width fits (a custom rule letting tracks that touch a fine-pitch part narrow).
#[derive(Clone, Debug)]
pub struct Neck {
    pub net: u32,
    /// Signal layer index.
    pub layer: usize,
    pub pad: String,
    /// Export geometry: from inside the pad to the stub's end.
    pub from: P2,
    pub to: P2,
    pub width: f64,
}

#[derive(Clone)]
pub struct Board {
    /// Neck-down stubs; each is also fixed copper of its net (objects named "neck:<pad>").
    pub necks: Vec<Neck>,
    /// Pads that needed a neck-down stub but had no room for one (why).
    pub neck_failed: Vec<String>,
    /// Neck-down width per pad a neck rule covers: a plane net's drop from it may be one straight segment at it.
    pub neck_width: HashMap<String, f64>,
    /// Plane regions; a plane net without any is a full plane (any via of the net joins it).
    pub plane_regions: Vec<PlaneRegion>,
    /// Areas of plane layers kept free of pour (cutouts).
    pub fill_keepouts: Vec<(usize, Vec<P2>)>,
    /// Rule areas that forbid vias but not tracks (e.g. along an RF line): wires pass, vias stay out.
    pub via_keepouts: Vec<Vec<P2>>,
    /// RF frontend nets (bench/pcb/topor_rf.py): each kept on the one layer it starts on, no vias, with why.
    pub rf_nets: Vec<(u32, String)>,
    pub name: String,
    pub size: [f64; 2],
    pub layers: Vec<String>,
    pub signal: Vec<usize>,
    pub classes: Vec<Class>,
    pub nets: Vec<NetInfo>,
    pub net_index: HashMap<String, u32>,
    pub objs: Vec<ObjSpec>,
    pub outline: Vec<P2>,
    pub outline_holes: Vec<Vec<P2>>,
    pub edge_clear: f64,
    pub hole_clear: f64,
    pub hole_to_hole: f64,
    /// Extra clearance between different nets of the same class (custom rules).
    pub class_pair: Vec<f64>,
    /// Least gap between a track and other copper per signal layer (a custom rule such as a fab's inner-layer
    /// spacing), for tracks of every class not exempt (`track_clear_exempt`, per class).
    pub track_clear: Vec<f64>,
    pub track_clear_exempt: Vec<bool>,
    pub routable_nets: usize,
    /// Differential pair partner per net (NONE if unpaired) and the pair's edge-to-edge gap.
    pub partner: Vec<u32>,
    pub pair_gap: Vec<f64>,
    /// Net family: a pair's members and its ribbon share the first member's id (others: own id).
    pub family: Vec<u32>,
    /// Ribbon net of each pair member (NONE if none).
    pub ribbon_of: Vec<u32>,
}

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap_or(0.0)
}

impl Board {
    pub fn nsig(&self) -> usize {
        self.signal.len()
    }
    pub fn sig_name(&self, l: usize) -> &str {
        &self.layers[self.signal[l]]
    }

    pub fn load(m: &Value) -> Result<Board, String> {
        let name = m["name"].as_str().unwrap_or("board").to_string();
        let size = [f(&m["size"][0]), f(&m["size"][1])];
        let layers: Vec<String> = m["layers"].as_array().ok_or("layers")?.iter().map(|x| x.as_str().unwrap().to_string()).collect();
        let cold = &m["cold_routing"];
        let mut plane_of: HashMap<String, usize> = HashMap::new();
        if let Some(pl) = cold["planes"].as_object() {
            for (net, layer) in pl {
                let li = layers.iter().position(|l| l == layer.as_str().unwrap()).ok_or("plane layer")?;
                plane_of.insert(net.clone(), li);
            }
        }
        let mut plane_layers: Vec<usize> = plane_of.values().cloned().collect();
        for l in cold["plane_layers"].as_array().into_iter().flatten() {
            if let Some(li) = layers.iter().position(|x| Some(x.as_str()) == l.as_str()) {
                if !plane_layers.contains(&li) {
                    plane_layers.push(li);
                }
            }
        }
        let signal: Vec<usize> = (0..layers.len()).filter(|i| !plane_layers.contains(i)).collect();
        let gr = &m["global_rules"];
        let min_width = gr["min_track_width"].as_f64().unwrap_or(0.1);
        let edge_clear = gr["min_copper_edge_clearance"].as_f64().unwrap_or(0.5);
        let hole_clear = gr["min_hole_clearance"].as_f64().unwrap_or(0.25);
        let hole_to_hole = gr["min_hole_to_hole"].as_f64().unwrap_or(0.25);

        let rules = &m["rules"];
        let mut classes = Vec::new();
        let mut class_ix: HashMap<String, usize> = HashMap::new();
        if let Some(cl) = m["classes"].as_object() {
            let mut names: Vec<&String> = cl.keys().collect();
            names.sort();
            for n in names {
                let r = &cl[n];
                class_ix.insert(n.clone(), classes.len());
                classes.push(Class { name: n.clone(), width: f(&r[0]), clear: f(&r[1]), via_d: f(&r[2]), via_drill: f(&r[3]) });
            }
        }
        if !class_ix.contains_key("Default") {
            class_ix.insert("Default".into(), classes.len());
            classes.push(Class { name: "Default".into(), width: f(&rules[0]), clear: f(&rules[1]), via_d: f(&rules[2]), via_drill: f(&rules[3]) });
        }
        let mut class_pair = vec![0.0; classes.len()];
        if let Some(cc) = cold["class_clearance"].as_object() {
            for (c, v) in cc {
                if let Some(&i) = class_ix.get(c) {
                    class_pair[i] = f(v);
                }
            }
        }
        let net_class = &m["net_class"];
        let mut nets: Vec<NetInfo> = Vec::new();
        let mut net_index: HashMap<String, u32> = HashMap::new();
        let mut net_of = |name: &str, nets: &mut Vec<NetInfo>| -> u32 {
            if let Some(&i) = net_index.get(name) {
                return i;
            }
            let cname = net_class[name].as_str().unwrap_or("Default");
            let class = *class_ix.get(cname).unwrap_or(&class_ix["Default"]);
            let c = &classes[class];
            let disallow: Vec<String> = cold["disallow"][cname].as_array().map(|a| a.iter().map(|x| x.as_str().unwrap().to_string()).collect()).unwrap_or_default();
            let mut allowed = 0u32;
            let mut wpref = Vec::new();
            let mut wmin = Vec::new();
            for (si, &li) in signal.iter().enumerate() {
                let ln = &layers[li];
                if !disallow.contains(ln) {
                    allowed |= 1 << si;
                }
                let p = cold["layer_widths"][cname][ln].as_f64().unwrap_or(c.width);
                let lm = cold["layer_min_widths"][cname][ln].as_f64();
                let mn = match lm {
                    Some(v) => v,
                    None => min_width.max(p.min(0.15)),
                };
                wpref.push(p);
                wmin.push(mn.min(p));
            }
            let id = nets.len() as u32;
            nets.push(NetInfo { name: name.to_string(), class, plane: plane_of.get(name).cloned(), allowed, wpref, wmin, ribbon: None, via_size: None });
            net_index.insert(name.to_string(), id);
            id
        };

        let fsig = signal.iter().position(|&l| l == 0);
        let bsig = signal.iter().position(|&l| l == layers.len() - 1);
        let all_sig = (1u32 << signal.len()) - 1;
        let mut objs = Vec::new();
        // pads a neck-down rule covers: (object index, neck width (0: the class width), long axis if the pad is
        // elongated, clearance of the stub (0: the class clearance))
        let mut escapes: Vec<(usize, f64, Option<P2>, f64)> = Vec::new();
        for p in m["pads"].as_array().ok_or("pads")? {
            let net = match p["net"].as_str() {
                Some(n) if !n.is_empty() => net_of(n, &mut nets),
                _ => NONE,
            };
            let lay = p["layers"].as_str().unwrap_or("all");
            let npth = p["npth"].as_bool().unwrap_or(false);
            let layers_mask = match lay {
                "F" => fsig.map(|i| 1u32 << i).unwrap_or(0),
                "B" => bsig.map(|i| 1u32 << i).unwrap_or(0),
                _ => all_sig,
            };
            let (x, y, w, h) = (f(&p["x"]), f(&p["y"]), f(&p["w"]), f(&p["h"]));
            let angle = f(&p["angle"]).to_radians();
            let shape = p["shape"].as_str().unwrap_or("rect");
            let rr = f(&p["rr"]);
            let local = f(&p["clearance"]);
            let drill = f(&p["drill_w"]).max(f(&p["drill_h"])).max(f(&p["drill"]));
            let (mut core, r): (Vec<P2>, f64) = if npth {
                let (dw, dh) = (f(&p["drill_w"]).max(f(&p["drill"])), f(&p["drill_h"]).max(f(&p["drill"])));
                if (dw - dh).abs() < 1e-6 {
                    (vec![p2(0., 0.)], dw / 2.)
                } else if dw > dh {
                    (vec![p2(-(dw - dh) / 2., 0.), p2((dw - dh) / 2., 0.)], dh / 2.)
                } else {
                    (vec![p2(0., -(dh - dw) / 2.), p2(0., (dh - dw) / 2.)], dw / 2.)
                }
            } else {
                match shape {
                    "circle" => (vec![p2(0., 0.)], w / 2.),
                    "oval" => {
                        if (w - h).abs() < 1e-6 {
                            (vec![p2(0., 0.)], w / 2.)
                        } else if w > h {
                            (vec![p2(-(w - h) / 2., 0.), p2((w - h) / 2., 0.)], h / 2.)
                        } else {
                            (vec![p2(0., -(h - w) / 2.), p2(0., (h - w) / 2.)], w / 2.)
                        }
                    }
                    "roundrect" => {
                        let r = rr.min(w / 2.).min(h / 2.);
                        let (ax, ay) = (w / 2. - r, h / 2. - r);
                        (vec![p2(-ax, -ay), p2(ax, -ay), p2(ax, ay), p2(-ax, ay)], r)
                    }
                    _ => (vec![p2(-w / 2., -h / 2.), p2(w / 2., -h / 2.), p2(w / 2., h / 2.), p2(-w / 2., h / 2.)], 0.),
                }
            };
            // KiCad orientations turn counter-clockwise as seen, with y pointing down
            let (s, c) = angle.sin_cos();
            for q in core.iter_mut() {
                *q = p2(x + q.x * c + q.y * s, y - q.x * s + q.y * c);
            }
            let core = dedup_hull(&core);
            let kind = if npth { Kind::Hole } else { Kind::Pad };
            // KiCad keeps copper off a non-plated hole by the board-edge clearance as well; a plated hole keeps
            // other copper the hole clearance off its drill, which beyond a thin annular ring is the larger gap
            let clear = if npth {
                hole_clear.max(edge_clear).max(local)
            } else if drill > 0.0 {
                // the thinnest ring, axis by axis (a slot in an oval pad runs along it)
                let (dw, dh) = (f(&p["drill_w"]).max(f(&p["drill"])), f(&p["drill_h"]).max(f(&p["drill"])));
                let ring = ((w - dw).min(h - dh) / 2.0).max(0.0);
                local.max(hole_clear - ring)
            } else {
                local
            };
            let name = format!("{}.{}", p["original_ref"].as_str().unwrap_or("?"), p["original_pin"].as_str().unwrap_or("?"));
            let (e, ec) = (f(&p["escape_width"]), f(&p["escape_clear"]));
            if (e > 0.0 || ec > 0.0) && lay != "all" && !npth && net != NONE {
                let axis = if (w - h).abs() < 0.05 * w.max(h) { None } else if w > h { Some(p2(c, -s)) } else { Some(p2(s, c)) };
                escapes.push((objs.len(), e, axis, ec));
            }
            objs.push(ObjSpec {
                kind,
                net: if npth { NONE } else { net },
                layers: layers_mask,
                pts: core,
                r,
                clear,
                center: p2(x, y),
                drill,
                name,
                through: lay == "all" && !npth,
            });
        }
        // small copper zones that carry connectivity (switch nodes, antenna feeds): each filled island
        // as convex pieces of fixed copper of its net (tied to the pads it reaches, not merged with them);
        // a plain outline (an older model) is one piece
        for z in m["zones"].as_array().into_iter().flatten() {
            let net = match z["net"].as_str() {
                Some(n) if !n.is_empty() => net_of(n, &mut nets),
                _ => continue,
            };
            let li = match signal.iter().position(|&l| layers[l] == z["layer"].as_str().unwrap_or("")) {
                Some(i) => i,
                None => continue,
            };
            let ring = |v: &Value| -> Vec<P2> { v.as_array().map(|a| a.iter().map(|q| p2(f(&q[0]), f(&q[1]))).collect()).unwrap_or_default() };
            let pieces: Vec<Vec<P2>> = match z["pieces"].as_array() {
                Some(a) => a.iter().map(ring).collect(),
                None => vec![ring(&z["points"])],
            };
            let clear = classes[nets[net as usize].class].clear.max(z["clearance"].as_f64().unwrap_or(0.0));
            let name = format!("zone:{}#{}", z["net"].as_str().unwrap_or(""), z["island"].as_u64().unwrap_or(0));
            for pts in pieces {
                let hull = dedup_hull(&pts);
                if hull.len() < 3 {
                    continue;
                }
                let c = hull.iter().fold(p2(0., 0.), |a, &b| a + b) * (1.0 / hull.len() as f64);
                objs.push(ObjSpec { kind: Kind::Pad, net, layers: 1 << li, pts: hull, r: 0.0, clear, center: c, drill: 0.0, name: name.clone(), through: false });
            }
        }
        let sig_of = |ln: &str| -> Option<usize> { signal.iter().position(|&l| layers[l] == ln) };
        let poly = |v: &Value| -> Vec<P2> { v.as_array().map(|a| a.iter().map(|q| p2(f(&q[0]), f(&q[1]))).collect()).unwrap_or_default() };
        let default_clear = classes[class_ix["Default"]].clear;
        let keep = |pts: Vec<P2>, mask: u32, clear: f64, kind: Kind, name: &str, objs: &mut Vec<ObjSpec>| {
            if pts.len() < 3 || mask == 0 {
                return;
            }
            let c = pts.iter().fold(p2(0., 0.), |a, &b| a + b) * (1.0 / pts.len() as f64);
            objs.push(ObjSpec { kind, net: NONE, layers: mask, pts, r: 0., clear, center: c, drill: 0., name: name.into(), through: false });
        };
        for g in m["copper_graphics"].as_array().into_iter().flatten() {
            if let Some(l) = sig_of(g["layer"].as_str().unwrap_or("")) {
                let pg = poly(&g["points"]);
                // a net tie: copper joining pads of different nets. Each pad owns the part of it nearer to
                // that pad than to the others (split 10 µm short of the bisector): copper of the pad's net
                // that its wires may end on, where a keepout would cover the pads themselves
                let touch: Vec<usize> = (0..objs.len())
                    .filter(|&i| {
                        let o = &objs[i];
                        o.kind == Kind::Pad && o.net != NONE && o.layers & (1 << l) != 0 && pg.len() >= 3 && (point_in_poly(&pg, o.center) || o.pts.iter().any(|&q| point_in_poly(&pg, q)) || poly_edge_dist(&pg, o.center) < o.r + 1e-3)
                    })
                    .collect();
                if touch.iter().any(|&i| objs[i].net != objs[touch[0]].net) {
                    for &i in &touch {
                        let mut part = dedup_hull(&pg);
                        for &j in &touch {
                            if objs[j].net == objs[i].net {
                                continue;
                            }
                            let (a, b) = (objs[i].center, objs[j].center);
                            let n = (b - a) * (1.0 / a.dist(b).max(1e-9));
                            let off = (a + b) * 0.5 - n * 0.005;
                            part = clip_half(&part, off, n);
                        }
                        if part.len() >= 3 {
                            let (net, clear) = (objs[i].net, objs[i].clear);
                            let c = part.iter().fold(p2(0., 0.), |s, &q| s + q) * (1.0 / part.len() as f64);
                            objs.push(ObjSpec { kind: Kind::Pad, net, layers: 1 << l, pts: part, r: 0.0, clear, center: c, drill: 0.0, name: "net_tie".into(), through: false });
                        }
                    }
                    continue;
                }
                keep(pg, 1 << l, default_clear, Kind::Keepout, "copper_graphic", &mut objs);
            }
        }
        for g in m["mask_graphics"].as_array().into_iter().flatten() {
            if let Some(l) = sig_of(g["layer"].as_str().unwrap_or("")) {
                keep(poly(&g["points"]), 1 << l, 0.0, Kind::Keepout, "mask_opening", &mut objs);
            }
        }
        for g in m["edge_graphics"].as_array().into_iter().flatten() {
            keep(poly(&g["points"]), all_sig, edge_clear, Kind::Outline, "edge_graphic", &mut objs);
        }
        let mut via_keepouts = Vec::new();
        for a in m["rule_areas"].as_array().into_iter().flatten() {
            // only areas that forbid copper wires keep the router out (a fill-only area does not);
            // one that forbids vias alone keeps out vias only
            let no_tracks = a["tracks"].as_bool().unwrap_or(true);
            if !no_tracks && !a["vias"].as_bool().unwrap_or(true) {
                continue;
            }
            if !no_tracks {
                via_keepouts.extend(a["polygons"].as_array().into_iter().flatten().map(|pg| poly(pg)));
                continue;
            }
            let mut mask = 0;
            for l in a["layers"].as_array().into_iter().flatten() {
                if let Some(i) = sig_of(l.as_str().unwrap_or("")) {
                    mask |= 1 << i;
                }
            }
            for pg in a["polygons"].as_array().into_iter().flatten() {
                keep(poly(pg), mask, 0.0, Kind::Keepout, "rule_area", &mut objs);
            }
        }
        let outline = poly(&m["outline"]);
        let outline_holes: Vec<Vec<P2>> = m["outline_holes"].as_array().into_iter().flatten().map(|h| poly(h)).collect();
        // nets with no pads still get an id via net_class (ignored)
        let routable_nets = m["routable_nets"].as_u64().unwrap_or(0) as usize;
        // plane regions (partial pours): a via joins its net's plane only inside one of them
        let ppoly = |v: &Value| -> Vec<P2> { v.as_array().map(|a| a.iter().map(|q| p2(f(&q[0]), f(&q[1]))).collect()).unwrap_or_default() };
        let mut plane_regions = Vec::new();
        for r in cold["plane_regions"].as_array().into_iter().flatten() {
            let (net, layer) = match (r["net"].as_str().and_then(|n| net_index.get(n)), layers.iter().position(|l| Some(l.as_str()) == r["layer"].as_str())) {
                (Some(&n), Some(l)) => (n, l),
                _ => continue,
            };
            plane_regions.push(PlaneRegion {
                net,
                layer,
                priority: r["priority"].as_i64().unwrap_or(0),
                poly: ppoly(&r["points"]),
                holes: r["holes"].as_array().into_iter().flatten().map(|h| ppoly(h)).collect(),
            });
        }
        let mut track_clear = vec![0.0; signal.len()];
        if let Some(tc) = cold["track_clearance"].as_object() {
            for (ln, v) in tc {
                if let Some(si) = signal.iter().position(|&li| &layers[li] == ln) {
                    track_clear[si] = f(v);
                }
            }
        }
        let exempt: Vec<String> = cold["track_clearance_except"].as_array().map(|a| a.iter().filter_map(|x| x.as_str().map(String::from)).collect()).unwrap_or_default();
        let track_clear_exempt: Vec<bool> = classes.iter().map(|c| exempt.contains(&c.name)).collect();
        let neck_width: HashMap<String, f64> = escapes.iter().filter(|x| x.1 > 0.0).map(|&(i, e, _, _)| (objs[i].name.clone(), e)).collect();
        let (necks, neck_failed) = place_necks(&mut objs, &escapes, &nets, &classes, &track_clear, &track_clear_exempt, &outline, &outline_holes, edge_clear);
        let mut fill_keepouts = Vec::new();
        for k in cold["fill_keepouts"].as_array().into_iter().flatten() {
            if let Some(l) = layers.iter().position(|x| Some(x.as_str()) == k["layer"].as_str()) {
                fill_keepouts.push((l, ppoly(&k["points"])));
            }
        }
        // RF frontend nets keep to their pads' layer: a via in an RF line is an impedance step and a stub
        let mut rf_nets = Vec::new();
        for r in m["rf_nets"].as_array().into_iter().flatten() {
            let (Some(name), Some(ln)) = (r["net"].as_str(), r["layer"].as_str()) else { continue };
            let (Some(&id), Some(si)) = (net_index.get(name), signal.iter().position(|&li| layers[li] == ln)) else { continue };
            let ni = &mut nets[id as usize];
            if ni.allowed & (1 << si) != 0 {
                ni.allowed = 1 << si;
                rf_nets.push((id, format!("{} on {} ({})", name, ln, r["why"].as_str().unwrap_or(""))));
            }
        }
        Ok(Board {
            necks,
            neck_failed,
            neck_width,
            plane_regions,
            fill_keepouts,
            via_keepouts,
            rf_nets,
            name,
            size,
            layers,
            signal,
            classes,
            nets,
            net_index,
            objs,
            outline,
            outline_holes,
            edge_clear,
            hole_clear,
            hole_to_hole,
            class_pair,
            track_clear,
            track_clear_exempt,
            routable_nets,
            partner: Vec::new(),
            pair_gap: Vec::new(),
            family: Vec::new(),
            ribbon_of: Vec::new(),
        })
    }

    /// Pair nets named like ..._DP/_DM, ..._D+/_D-, ..._P/_N within one class. The gap comes from the
    /// project's net class (diff_pair_gap), 0.2 mm if absent.
    pub fn detect_pairs(&mut self, project: Option<&Value>) {
        let n = self.nets.len();
        self.partner = vec![NONE; n];
        self.pair_gap = vec![0.0; n];
        let gap_of = |cls: &str| -> f64 {
            project
                .and_then(|p| p["net_settings"]["classes"].as_array())
                .and_then(|a| a.iter().find(|c| c["name"].as_str() == Some(cls)))
                .and_then(|c| c["diff_pair_gap"].as_f64())
                .unwrap_or(0.2)
        };
        let suffixes = [("_DP", "_DM"), ("_D+", "_D-"), ("D+", "D-"), ("_P", "_N")];
        for i in 0..n {
            let name = self.nets[i].name.clone();
            for (sp, sn) in suffixes {
                if let Some(stem) = name.strip_suffix(sp) {
                    let other = format!("{}{}", stem, sn);
                    if let Some(&j) = self.net_index.get(&other) {
                        if self.nets[j as usize].class == self.nets[i].class && self.partner[i] == NONE {
                            let g = gap_of(&self.classes[self.nets[i].class].name);
                            self.partner[i] = j;
                            self.partner[j as usize] = i as u32;
                            self.pair_gap[i] = g;
                            self.pair_gap[j as usize] = g;
                        }
                    }
                    break;
                }
            }
        }
        // one ribbon net per pair: both members side by side, routed as one wide wire
        self.family = (0..n as u32).collect();
        self.ribbon_of = vec![NONE; n];
        for i in 0..n {
            let j = self.partner[i] as usize;
            if self.partner[i] == NONE || j < i {
                continue;
            }
            let (a, b) = (&self.nets[i], &self.nets[j]);
            let c = &self.classes[a.class];
            // the members' centre spacing beyond their width: the coupling gap, at least the clearance
            let gap = self.pair_gap[i].max(c.clear).max(self.class_pair.get(a.class).cloned().unwrap_or(0.0)) + 0.01;
            let wpref: Vec<f64> = a.wpref.iter().zip(&b.wpref).map(|(x, y)| x.max(*y) * 2.0 + gap).collect();
            let wmin: Vec<f64> = a.wmin.iter().zip(&b.wmin).map(|(x, y)| x.max(*y) * 2.0 + gap).collect();
            let pitch = c.via_d + gap.max(c.clear) + 0.02;
            let r = NetInfo {
                name: format!("ribbon:{}|{}", a.name, b.name),
                class: a.class,
                plane: None,
                allowed: a.allowed & b.allowed,
                wpref,
                wmin,
                ribbon: Some((i as u32, j as u32)),
                via_size: Some((pitch + c.via_d, pitch + c.via_drill)),
            };
            let id = self.nets.len() as u32;
            self.net_index.insert(r.name.clone(), id);
            self.nets.push(r);
            self.partner.push(NONE);
            self.pair_gap.push(gap);
            self.family.push(i as u32);
            self.ribbon_of.push(NONE);
            self.family[j] = i as u32;
            self.ribbon_of[i] = id;
            self.ribbon_of[j] = id;
        }
    }

    /// Whether copper of plane net `net` reaching all layers at pos (radius r) joins its plane: inside one
    /// of its regions, clear of the region's edge and holes, of higher-priority pours of other nets on
    /// that layer and of fill keepouts (the pour keeps its clearance and needs room for a thermal).
    pub fn via_kept_out(&self, pos: P2, r: f64) -> bool {
        self.via_keepouts.iter().any(|k| point_in_poly(k, pos) || poly_edge_dist(k, pos) < r + 0.01)
    }

    pub fn in_plane(&self, net: u32, pos: P2, r: f64) -> bool {
        self.in_plane_on(net, pos, r, None)
    }

    /// in_plane() for copper on one layer only (an SMD pad): a region of its net on that layer. A plane
    /// net without any region joins its plane anywhere, as through copper.
    pub fn in_plane_on(&self, net: u32, pos: P2, r: f64, layer: Option<usize>) -> bool {
        let mut any = false;
        let m = r + 0.25;
        let near = |poly: &[P2]| point_in_poly(poly, pos) || poly_edge_dist(poly, pos) < m;
        for reg in self.plane_regions.iter().filter(|x| x.net == net) {
            any = true;
            if layer.map_or(false, |l| l != reg.layer) {
                continue;
            }
            if !point_in_poly(&reg.poly, pos) || poly_edge_dist(&reg.poly, pos) < m || reg.holes.iter().any(|h| near(h)) {
                continue;
            }
            if self.plane_regions.iter().any(|o| o.layer == reg.layer && o.net != net && o.priority > reg.priority && near(&o.poly)) {
                continue;
            }
            if self.fill_keepouts.iter().any(|(l, k)| *l == reg.layer && near(k)) {
                continue;
            }
            return true;
        }
        !any && layer.is_none()
    }

    /// Via (diameter, drill) for a net.
    pub fn via_size(&self, net: u32) -> (f64, f64) {
        let ni = &self.nets[net as usize];
        ni.via_size.unwrap_or((self.classes[ni.class].via_d, self.classes[ni.class].via_drill))
    }
    /// A via's gap to other copper: its class clearance, or the hole clearance off its drill where the
    /// annular ring is thinner than the difference.
    pub fn via_clear(&self, net: u32) -> f64 {
        let (d, drill) = self.via_size(net);
        self.classes[self.nets[net as usize].class].clear.max(self.hole_clear - (d - drill) / 2.0)
    }
    /// A track's gap to other copper on signal layer `l`.
    pub fn track_clear(&self, net: u32, l: usize) -> f64 {
        let c = self.nets[net as usize].class;
        let t = if self.track_clear_exempt[c] { 0.0 } else { self.track_clear.get(l).cloned().unwrap_or(0.0) };
        self.classes[c].clear.max(t)
    }
}

/// Neck-down stubs for the pads that need one (plane nets' pads drop to vias instead): a pad narrower than its
/// net's class width (or one where that width meets other copper) gets a straight stub at its neck width, outward
/// along the pad (from the part's centre for a square pad), to the first point where a disc of the class width
/// clears the other nets' copper and the class-width starts of the stubs laid before it (so neighbouring stubs
/// stagger), with the stub itself clear.
/// The stub stops 10 µm short of its pad (touching copper, tied to the pad, not hulled with it) and is pushed as
/// fixed copper of the net, named "neck:<pad>"; KiCad gets it as a track from the pad's centre.
#[allow(clippy::too_many_arguments)]
fn place_necks(objs: &mut Vec<ObjSpec>, escapes: &[(usize, f64, Option<P2>, f64)], nets: &[NetInfo], classes: &[Class], track_clear: &[f64], exempt: &[bool], outline: &[P2], holes: &[Vec<P2>], edge_clear: f64) -> (Vec<Neck>, Vec<String>) {
    let mut out = Vec::new();
    let mut failed = Vec::new();
    let mut centre: HashMap<String, (P2, usize)> = HashMap::new();
    for o in objs.iter().filter(|o| o.kind == Kind::Pad) {
        let e = centre.entry(o.name.split('.').next().unwrap_or("").to_string()).or_insert((p2(0., 0.), 0));
        e.0 = e.0 + o.center;
        e.1 += 1;
    }
    // stub ends laid so far: (end, class half-width, clearance, net)
    let mut ends: Vec<(P2, f64, f64, u32)> = Vec::new();
    for &(pi, e, axis, ec) in escapes {
        let pad = objs[pi].clone();
        let net = pad.net;
        let ni = &nets[net as usize];
        // a plane net's pads drop to vias: a straight drop at the neck width needs no stub
        if ni.plane.is_some() {
            continue;
        }
        let Some(l) = (0..32).find(|&l| pad.layers == 1 << l) else { continue };
        let wide = ni.wmin[l].max(ni.wpref[l]);
        let c = ni.class;
        let my_clear = classes[c].clear.max(if exempt[c] { 0.0 } else { track_clear.get(l).cloned().unwrap_or(0.0) });
        // the stub: narrower than the class, or at the class width with a smaller clearance
        let e = if e > 0.0 { e } else { wide };
        let stub_clear = if ec > 0.0 { ec.min(my_clear) } else { my_clear };
        if e >= wide && stub_clear >= my_clear {
            continue;
        }
        // (a pad's own clearance is KiCad's, its class and custom rules resolved; a stub's is its own)
        let need = |o: &ObjSpec, mine: f64| -> f64 {
            match o.kind {
                Kind::Outline => edge_clear,
                Kind::Hole | Kind::Keepout => o.clear,
                _ => mine.max(o.clear),
            }
        };
        // everything within 4 mm (by bounding box: a long keepout's centre may lie far away)
        let near = |o: &ObjSpec| {
            let (lo, hi) = o.pts.iter().fold((p2(1e9, 1e9), p2(-1e9, -1e9)), |(lo, hi), q| (p2(lo.x.min(q.x), lo.y.min(q.y)), p2(hi.x.max(q.x), hi.y.max(q.y))));
            let dx = (lo.x - pad.center.x).max(pad.center.x - hi.x).max(0.0);
            let dy = (lo.y - pad.center.y).max(pad.center.y - hi.y).max(0.0);
            dx.hypot(dy) - o.r < 4.0
        };
        let others: Vec<&ObjSpec> = objs.iter().filter(|o| o.layers & (1 << l) != 0 && (o.net != net || o.kind == Kind::Hole) && near(o)).collect();
        // gap of a capsule a-b of radius r keeping `mine` to everything else (negative: overlap)
        let gap = |a: P2, b: P2, r: f64, mine: f64| -> f64 {
            let mut g = f64::INFINITY;
            for o in &others {
                let d = if o.kind == Kind::Keepout {
                    if point_in_poly(&o.pts, a) || point_in_poly(&o.pts, b) {
                        -1.0
                    } else {
                        (0..o.pts.len()).map(|k| crate::geom::seg_seg_dist(a, b, o.pts[k], o.pts[(k + 1) % o.pts.len()])).fold(f64::INFINITY, f64::min)
                    }
                } else {
                    crate::geom::seg_hull_dist(&o.pts, a, b) - o.r
                };
                g = g.min(d - r - need(o, mine));
            }
            for &(q, rq, cq, n) in &ends {
                if n != net {
                    g = g.min(crate::geom::seg_point_dist(a, b, q) - r - rq - cq.max(mine));
                }
            }
            for ring in std::iter::once(outline).chain(holes.iter().map(|h| h.as_slice())) {
                if ring.len() >= 3 {
                    let d = (0..ring.len()).map(|k| crate::geom::seg_seg_dist(a, b, ring[k], ring[(k + 1) % ring.len()])).fold(f64::INFINITY, f64::min);
                    g = g.min(d - r - edge_clear);
                }
            }
            g
        };
        let part = pad.name.split('.').next().unwrap_or("").to_string();
        let mid = centre.get(&part).map(|&(s, n)| s * (1.0 / n as f64)).unwrap_or(pad.center);
        let out_dir = pad.center - mid;
        let dir = match axis {
            Some(a) => if a.dot(out_dir) >= 0.0 { a } else { a * -1.0 },
            None if out_dir.norm() > 1e-6 => out_dir.unit(),
            None => continue,
        };
        // a pad as wide as the class width, with room for it, needs none
        let across = pad.pts.iter().map(|&q| (q - pad.center).dot(dir.perp()).abs()).fold(0.0, f64::max) + pad.r;
        if 2.0 * across >= wide && gap(pad.center, pad.center, wide / 2.0, my_clear) >= 0.0 {
            continue;
        }
        // where the pad's copper ends along the stub
        let core = pad.pts.iter().map(|&q| (q - pad.center).dot(dir)).fold(0.0, f64::max);
        let reach = core + pad.r;
        // a stub as narrow as the pad runs in from its centre, over its copper, and is modelled from just past
        // the pad (touching copper); one wider than the pad would stand out along its sides, so it starts, in the
        // model and on the board, just past the pad's core, inside its rounded end, and ends round there
        let wider = e > 2.0 * across + 1e-9;
        let start = if wider { pad.center + dir * (core + 0.01) } else { pad.center + dir * (reach + e / 2.0 + 0.01) };
        // (KiCad joins a track to a pad its end lies in: on a square-cornered pad that end stays just inside)
        let from = if wider { pad.center + dir * (core + 0.01).min(reach - 0.005) } else { pad.center };
        let mut found = None;
        let mut why = "no room for the class width within 2.5 mm";
        let mut t = (start - pad.center).dot(dir) + 0.01;
        while t < reach + 2.5 {
            let j = pad.center + dir * t;
            if gap(start, j, e / 2.0, stub_clear) < 0.0 {
                why = "the stub runs into other copper first";
                break;
            }
            if gap(j, j, wide / 2.0, my_clear) >= 0.005 {
                found = Some(j);
                break;
            }
            t += 0.01;
        }
        let Some(j) = found else {
            failed.push(format!("{} ({})", pad.name, why));
            continue;
        };
        let name = format!("neck:{}", pad.name);
        objs.push(ObjSpec { kind: Kind::Pad, net, layers: pad.layers, pts: vec![start, j], r: e / 2.0, clear: stub_clear, center: j, drill: 0.0, name, through: false });
        // (a wire from a stub with its own clearance starts in the courtyard, where that clearance holds too)
        ends.push((j, wide / 2.0, stub_clear, net));
        out.push(Neck { net, layer: l, pad: pad.name.clone(), from, to: j, width: e });
    }
    (out, failed)
}

/// The part of convex polygon `poly` on the side of the line through `o` that normal `n` points away from.
fn clip_half(poly: &[P2], o: P2, n: P2) -> Vec<P2> {
    let side = |q: P2| (q - o).dot(n);
    let mut out = Vec::new();
    for i in 0..poly.len() {
        let (a, b) = (poly[i], poly[(i + 1) % poly.len()]);
        let (sa, sb) = (side(a), side(b));
        if sa <= 0.0 {
            out.push(a);
        }
        if (sa < 0.0) != (sb < 0.0) && sa != sb {
            out.push(a.lerp(b, sa / (sa - sb)));
        }
    }
    dedup_hull(&out)
}

fn dedup_hull(pts: &[P2]) -> Vec<P2> {
    let mut v: Vec<P2> = Vec::new();
    for &p in pts {
        if !v.iter().any(|q| q.dist(p) < 1e-6) {
            v.push(p);
        }
    }
    if v.len() >= 3 {
        convex_hull(&v)
    } else {
        v
    }
}
