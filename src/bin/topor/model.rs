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
    /// Signal rise time in ns (`cold_routing.rise_time`, by net or class); 0 for a net that is not high-speed.
    /// A high-speed net's tracks keep off gaps in their reference plane (gaps.rs).
    pub rise: f64,
    /// Cost of a millimetre of wire per signal layer (1 by default; a length group's `layer_cost` raises it, e.g.
    /// to keep its long runs off the crowded pad layer and leave room to tune them).
    pub layer_cost: Vec<f64>,
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
    /// Rule areas that forbid vias: those that forbid vias but not tracks (e.g. along an RF line: wires pass, vias
    /// stay out), and those that forbid both (also keepout objects for the tracks, whose checks reach only near
    /// their edge).
    pub via_keepouts: Vec<Vec<P2>>,
    /// Bounding box of each via keepout (lo, hi).
    pub via_keepout_bb: Vec<(P2, P2)>,
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
    /// A ribbon may touch its own pair's copper (off for pair routing as it was, to compare).
    pub ribbon_touch: bool,
    /// Ribbon net of each pair member (NONE if none).
    pub ribbon_of: Vec<u32>,
    /// Reference planes (indices into `layers`) of each signal layer: the nearest plane layer on each side
    /// with no signal layer between, else the nearest plane layer (`cold_routing.reference_planes` overrides).
    pub ref_planes: Vec<Vec<usize>>,
    /// Pour clearance per layer (`cold_routing.plane_clearance`; 0: the other copper's class clearance) and the
    /// pour's least web width (`plane_min_width`; KiCad's zone minimum width, 0.25 mm by default).
    pub plane_clear: Vec<f64>,
    pub plane_min_w: Vec<f64>,
    /// Gap rule: a slot crossing is allowed while its inductance (Johnson, 5 D ln(D/W) nH) reflects at most
    /// `gap_reflection` of an edge of the net's rise time on a `gap_z0` line.
    pub gap_reflection: f64,
    pub gap_z0: f64,
    /// Length matching (`cold_routing.length_groups`): nets whose routed lengths must agree within a skew.
    pub length_groups: Vec<LengthGroup>,
    /// Meander geometry: edge-to-edge gap between a meander's legs (default the larger of three widths and
    /// the clearance) and the largest bump height (`meander_gap`, `meander_max_amplitude`; mm).
    pub meander_gap: Option<f64>,
    pub meander_max_amp: f64,
    /// Dogbone nets (`cold_routing.dogbone`: true for every net, or a list of nets): their vias go only right
    /// beside their own pads, so a wire either stays on its pads' layer or drops at one end and comes back up
    /// at the other, with no via between. Indexed by net (ribbon nets are not in it).
    pub dogbone: Vec<bool>,
}

/// Nets routed to one length: each within `max_skew` of the longest (or of `target`, when that is longer).
#[derive(Clone, Debug)]
pub struct LengthGroup {
    pub name: String,
    pub nets: Vec<u32>,
    pub max_skew: f64,
    pub target: Option<f64>,
    /// Names given that are not nets of the board.
    pub unknown: Vec<String>,
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
            nets.push(NetInfo { name: name.to_string(), class, plane: plane_of.get(name).cloned(), allowed, wpref, wmin, ribbon: None, via_size: None, rise: 0.0, layer_cost: vec![1.0; signal.len()] });
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
            let no_vias = a["vias"].as_bool().unwrap_or(true);
            for pg in a["polygons"].as_array().into_iter().flatten() {
                keep(poly(pg), mask, 0.0, Kind::Keepout, "rule_area", &mut objs);
                // one that forbids vias as well keeps them out wherever they are in it: the keepout object's checks
                // see only its edge near a via, not a via deep inside a large area
                if no_vias {
                    via_keepouts.push(poly(pg));
                }
            }
        }
        via_keepouts.retain(|k| k.len() >= 3);
        let via_keepout_bb = via_keepouts.iter().map(|k| bbox_of(k)).collect();
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
        // pours on signal layers (e.g. a VDD pour on an inner signal layer under a chip, filled after routing): the
        // net becomes a plane net of that layer, so its pads drop to vias and a via inside the pour joins it, as
        // with plane_regions; the layer stays a signal layer and other nets still route across it (the pour fills
        // round them later). The layer list above is settled first, so the pour does not make it a plane layer.
        for r in cold["signal_pours"].as_array().into_iter().flatten() {
            let (net, layer) = match (r["net"].as_str().and_then(|n| net_index.get(n)), layers.iter().position(|l| Some(l.as_str()) == r["layer"].as_str())) {
                (Some(&n), Some(l)) => (n, l),
                _ => continue,
            };
            if !signal.contains(&layer) {
                return Err(format!("signal_pours: {} is not a signal layer", layers[layer]));
            }
            match nets[net as usize].plane {
                Some(p) if p != layer => return Err(format!("signal_pours: net {} is already a plane net of {}", nets[net as usize].name, layers[p])),
                _ => nets[net as usize].plane = Some(layer),
            }
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
        // high-speed nets by rise time: a net's own entry wins over its class
        if let Some(rt) = cold["rise_time"].as_object() {
            for ni in nets.iter_mut() {
                let v = rt.get(&ni.name).or_else(|| rt.get(&classes[ni.class].name)).and_then(|v| v.as_f64()).unwrap_or(0.0);
                ni.rise = v.max(0.0);
            }
        }
        // per-layer value: one number for every layer or {layer: value}
        let per_layer = |v: &Value, d: f64| -> Vec<f64> {
            layers.iter().map(|l| v.as_f64().or_else(|| v[l.as_str()].as_f64()).unwrap_or(d)).collect()
        };
        let plane_clear = per_layer(&cold["plane_clearance"], 0.0);
        let plane_min_w = per_layer(&cold["plane_min_width"], 0.25);
        let mut ref_planes: Vec<Vec<usize>> = Vec::new();
        for &li in &signal {
            let mut v: Vec<usize> = Vec::new();
            if let Some(list) = cold["reference_planes"][layers[li].as_str()].as_array() {
                v = list.iter().filter_map(|x| layers.iter().position(|l| Some(l.as_str()) == x.as_str())).filter(|i| plane_layers.contains(i)).collect();
            } else {
                // adjacent plane layers above and below (no signal layer between)
                let up = (0..li).rev().take_while(|i| plane_layers.contains(i)).next();
                let down = (li + 1..layers.len()).take_while(|i| plane_layers.contains(i)).next();
                v.extend(up.into_iter().chain(down));
                if v.is_empty() {
                    if let Some(&p) = plane_layers.iter().min_by_key(|&&p| ((p as i64 - li as i64).abs(), p)) {
                        v.push(p);
                    }
                }
            }
            ref_planes.push(v);
        }
        let mut length_groups = Vec::new();
        for (gi, g) in cold["length_groups"].as_array().into_iter().flatten().enumerate() {
            let (mut ids, mut unknown) = (Vec::new(), Vec::new());
            for n in g["nets"].as_array().into_iter().flatten().filter_map(|x| x.as_str()) {
                match net_index.get(n) {
                    Some(&id) if !ids.contains(&id) => ids.push(id),
                    Some(_) => {}
                    None => unknown.push(n.to_string()),
                }
            }
            // {layer: cost of a millimetre} for the group's nets (at least 1)
            if let Some(lc) = g["layer_cost"].as_object() {
                for &id in &ids {
                    for (si, &li) in signal.iter().enumerate() {
                        if let Some(v) = lc.get(layers[li].as_str()).and_then(|v| v.as_f64()) {
                            let c = &mut nets[id as usize].layer_cost[si];
                            *c = c.max(v.max(1.0));
                        }
                    }
                }
            }
            length_groups.push(LengthGroup {
                name: g["name"].as_str().map(String::from).unwrap_or_else(|| format!("group {}", gi + 1)),
                nets: ids,
                max_skew: g["max_skew"].as_f64().unwrap_or(0.1).max(0.0),
                target: g["target"].as_f64(),
                unknown,
            });
        }
        let dogbone: Vec<bool> = match &cold["dogbone"] {
            Value::Bool(b) => vec![*b; nets.len()],
            Value::Array(a) => {
                let names: Vec<&str> = a.iter().filter_map(|x| x.as_str()).collect();
                nets.iter().map(|ni| names.contains(&ni.name.as_str())).collect()
            }
            _ => vec![false; nets.len()],
        };
        Ok(Board {
            dogbone,
            length_groups,
            meander_gap: cold["meander_gap"].as_f64(),
            meander_max_amp: cold["meander_max_amplitude"].as_f64().unwrap_or(1.0),
            ref_planes,
            plane_clear,
            plane_min_w,
            gap_reflection: cold["gap_reflection"].as_f64().unwrap_or(0.05),
            gap_z0: cold["gap_z0"].as_f64().unwrap_or(50.0),
            necks,
            neck_failed,
            neck_width,
            plane_regions,
            fill_keepouts,
            via_keepouts,
            via_keepout_bb,
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
            ribbon_touch: true,
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
                rise: [a.rise, b.rise].into_iter().filter(|&x| x > 0.0).fold(0.0, |m: f64, x| if m == 0.0 { x } else { m.min(x) }),
                layer_cost: a.layer_cost.iter().zip(&b.layer_cost).map(|(x, y)| x.max(*y)).collect(),
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
        let m = r + 0.01;
        self.via_keepouts.iter().zip(&self.via_keepout_bb).any(|(k, &(lo, hi))| pos.x > lo.x - m && pos.x < hi.x + m && pos.y > lo.y - m && pos.y < hi.y + m && (point_in_poly(k, pos) || poly_edge_dist(k, pos) < m))
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
    /// Whether a net keeps off gaps in its reference planes (it has a rise time).
    pub fn high_speed(&self, net: u32) -> bool {
        self.nets.get(net as usize).map_or(false, |n| n.rise > 0.0)
    }
    /// Longest slot (mm) a high-speed net's track on signal layer `l` may cross: Johnson's slot inductance
    /// L = 5 D ln(D/W) nH (D, W in inches; W the track width) reflects L / (2 Z0 tr) of an edge, kept at most
    /// `gap_reflection`. Infinite for a net that is not high-speed.
    pub fn gap_limit(&self, net: u32, l: usize) -> f64 {
        let ni = &self.nets[net as usize];
        if ni.rise <= 0.0 {
            return f64::INFINITY;
        }
        // a ribbon's members carry the signal: their own width
        let w = match ni.ribbon {
            Some((p, _)) => self.nets[p as usize].wpref[l],
            None => ni.wpref[l],
        }
        .max(0.01);
        // D ln(D/W) = L_max (nH) * 25.4 / 5, with L_max = 2 Z0 tr reflection (ohm ns = nH)
        let k = 2.0 * self.gap_z0 * ni.rise * self.gap_reflection * 25.4 / 5.0;
        let (mut lo, mut hi) = (w, 1000.0);
        for _ in 0..60 {
            let m = 0.5 * (lo + hi);
            if m * (m / w).ln() > k {
                hi = m;
            } else {
                lo = m;
            }
        }
        lo
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

fn bbox_of(pts: &[P2]) -> (P2, P2) {
    let lo = pts.iter().fold(p2(f64::INFINITY, f64::INFINITY), |a, q| p2(a.x.min(q.x), a.y.min(q.y)));
    let hi = pts.iter().fold(p2(f64::NEG_INFINITY, f64::NEG_INFINITY), |a, q| p2(a.x.max(q.x), a.y.max(q.y)));
    (lo, hi)
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

#[cfg(test)]
mod signal_pour_tests {
    use super::*;
    use serde_json::json;

    fn model(pours: serde_json::Value) -> serde_json::Value {
        json!({
            "name": "t", "size": [20.0, 20.0], "origin": [0.0, 0.0],
            "outline": [[0.0, 0.0], [20.0, 0.0], [20.0, 20.0], [0.0, 20.0]], "outline_holes": [],
            "layers": ["F.Cu", "In1.Cu", "In2.Cu", "B.Cu"], "rules": [0.15, 0.15, 0.5, 0.3],
            "classes": {"Default": [0.15, 0.15, 0.5, 0.3]}, "net_class": {"VDD": "Default", "GND": "Default"},
            "pads": [
                {"ref": "C1", "num": "1", "x": 2.0, "y": 2.0, "w": 0.5, "h": 0.5, "rot": 0, "shape": "rect", "layers": "F", "net": "VDD"},
                {"ref": "C2", "num": "1", "x": 8.0, "y": 8.0, "w": 0.5, "h": 0.5, "rot": 0, "shape": "rect", "layers": "F", "net": "VDD"}
            ],
            "zones": [], "rule_areas": [], "routable_nets": 1,
            "global_rules": {"min_track_width": 0.1, "min_copper_edge_clearance": 0.3},
            "cold_routing": {"planes": {"GND": "B.Cu"}, "signal_pours": pours}
        })
    }

    /// A pour on a signal layer makes its net a plane net of that layer (vias inside join it) while the layer
    /// stays a signal layer for every other net.
    #[test]
    fn a_signal_layer_pour_joins_vias_inside_it_and_keeps_the_layer_routable() {
        let m = model(json!([{"net": "VDD", "layer": "In2.Cu", "points": [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]}]));
        let b = Board::load(&m).expect("model");
        let in2 = b.layers.iter().position(|l| l == "In2.Cu").unwrap();
        assert!(b.signal.contains(&in2), "In2 stays a signal layer");
        let vdd = b.net_index["VDD"];
        assert_eq!(b.nets[vdd as usize].plane, Some(in2));
        assert!(b.in_plane(vdd, p2(5.0, 5.0), 0.3), "a via inside the pour joins it");
        assert!(!b.in_plane(vdd, p2(15.0, 15.0), 0.3), "outside it does not");
        assert!(!b.in_plane_on(vdd, p2(5.0, 5.0), 0.25, Some(0)), "an F.Cu pad over it still needs a via");
    }

    #[test]
    fn a_signal_pour_on_a_plane_layer_is_refused() {
        let m = model(json!([{"net": "VDD", "layer": "B.Cu", "points": [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]}]));
        assert!(Board::load(&m).is_err());
    }
}
