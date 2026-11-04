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

#[derive(Clone)]
pub struct Board {
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
            // KiCad keeps copper off a non-plated hole by the board-edge clearance as well
            let clear = if npth { hole_clear.max(edge_clear).max(local) } else { local };
            let name = format!("{}.{}", p["original_ref"].as_str().unwrap_or("?"), p["original_pin"].as_str().unwrap_or("?"));
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
