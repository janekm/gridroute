//! topor: a topological (rubber-band) PCB router engine, inspired by TopoR.
mod cdt;
mod funnel;
mod fx;
mod gaps;
mod geom;
mod model;
mod router;
mod search;
mod spring;
mod svg;
mod tune;
mod world;

use crate::geom::P2;
use crate::model::Kind;
use serde_json::json;
use std::path::PathBuf;
use std::time::Instant;

// the search allocates many small, short-lived vectors: mimalloc serves them far faster than the system allocator
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Move the suggested tracks and vias at these indices into the model's fixed copper: each track segment a zone
/// piece (its capsule outline), each via a fixed via pad, as drawn. Returns the tracks and vias moved, for the
/// layout (they are not in the source board).
fn fix_suggestions(m: &mut serde_json::Value, tracks: &[usize], vias: &[usize]) -> (Vec<serde_json::Value>, Vec<serde_json::Value>) {
    let cls = |m: &serde_json::Value, net: &str| -> Vec<f64> {
        let c = m["net_class"][net].as_str().unwrap_or("Default");
        let r = if m["classes"][c].is_array() { &m["classes"][c] } else { &m["rules"] };
        (0..4).map(|i| r[i].as_f64().unwrap_or(0.0)).collect()
    };
    let sug = m["suggested_copper"].clone();
    // (ids go on from earlier rounds' fixed copper)
    let z0 = m["zones"].as_array().map_or(0, |a| a.iter().filter(|z| z["island"].as_u64().unwrap_or(0) >= 1_000_000).count());
    let v0 = m["pads"].as_array().map_or(0, |a| a.iter().filter(|p| p["ref"].as_str().map_or(false, |r| r.starts_with("via_s"))).count());
    let mut zones = Vec::new();
    for (k, &ti) in tracks.iter().enumerate() {
        let t = &sug["tracks"][ti];
        let pts: Vec<[f64; 2]> = t["pts"].as_array().into_iter().flatten().map(|q| [q[0].as_f64().unwrap_or(0.0), q[1].as_f64().unwrap_or(0.0)]).collect();
        let w = t["width"].as_f64().unwrap_or(0.1);
        for (i, s) in pts.windows(2).enumerate() {
            // the copper as an outline: half discs round both ends, circumscribed (the polygon covers the copper)
            let n = 8;
            let r = w / 2.0 / (std::f64::consts::PI / (2.0 * n as f64)).cos();
            let (a, b) = (s[0], s[1]);
            let ang = if a != b { (b[1] - a[1]).atan2(b[0] - a[0]) } else { 0.0 };
            let mut poly = Vec::new();
            for (c, a0) in [(b, ang - std::f64::consts::FRAC_PI_2), (a, ang + std::f64::consts::FRAC_PI_2)] {
                for j in 0..=n {
                    let u = a0 + std::f64::consts::PI * j as f64 / n as f64;
                    poly.push(json!([c[0] + r * u.cos(), c[1] + r * u.sin()]));
                }
            }
            zones.push(json!({"net": t["net"], "layer": t["layer"], "island": 1_000_000 + (z0 + k) * 1000 + i, "clearance": 0.0, "points": poly}));
        }
    }
    let mut pads = Vec::new();
    let mut out_vias = Vec::new();
    for &vi in vias {
        let v = &sug["vias"][vi];
        let net = v["net"].as_str().unwrap_or("");
        let c = cls(m, net);
        let (d, dr) = (v["d"].as_f64().unwrap_or(c[2]), v["drill"].as_f64().unwrap_or(c[3]));
        let id = v0 + pads.len();
        pads.push(json!({"ref": format!("via_s{}", id), "original_ref": "via", "original_pin": format!("s{}", id), "num": "1", "x": v["x"], "y": v["y"], "w": d, "h": d, "rot": 0, "angle": 0,
                         "shape": "circle", "layers": "all", "npth": false, "drill": dr, "drill_w": dr, "drill_h": dr, "net": net, "clearance": c[1], "rr": 0}));
        out_vias.push(json!({"net": net, "x": v["x"], "y": v["y"], "d": d, "drill": dr}));
    }
    let out_tracks: Vec<serde_json::Value> = tracks.iter().map(|&ti| sug["tracks"][ti].clone()).collect();
    if let Some(a) = m["zones"].as_array_mut() {
        a.extend(zones);
    } else {
        m["zones"] = json!(zones);
    }
    if let Some(a) = m["pads"].as_array_mut() {
        a.extend(pads);
    }
    let keep = |arr: &serde_json::Value, drop: &[usize]| -> serde_json::Value { json!(arr.as_array().into_iter().flatten().enumerate().filter(|(i, _)| !drop.contains(i)).map(|(_, x)| x.clone()).collect::<Vec<_>>()) };
    m["suggested_copper"]["tracks"] = keep(&sug["tracks"], tracks);
    m["suggested_copper"]["vias"] = keep(&sug["vias"], vias);
    (out_tracks, out_vias)
}

/// The suggested tracks and vias at these indices, with the partner's suggested copper beside them (within 1 mm, on
/// the same layer, or any via) for those of a differential pair member.
fn with_partner_copper(sug: &serde_json::Value, tracks: &[usize], vias: &[usize], partner: &dyn Fn(&str) -> Option<String>) -> (Vec<usize>, Vec<usize>) {
    use geom::{p2, seg_point_dist, seg_seg_dist, P2};
    let f = |x: &serde_json::Value| x.as_f64().unwrap_or(0.0);
    let segs = |t: &serde_json::Value| -> Vec<(P2, P2)> {
        let pts: Vec<P2> = t["pts"].as_array().into_iter().flatten().map(|q| p2(f(&q[0]), f(&q[1]))).collect();
        pts.windows(2).map(|w| (w[0], w[1])).collect()
    };
    let at = |v: &serde_json::Value| p2(f(&v["x"]), f(&v["y"]));
    let st = sug["tracks"].as_array().cloned().unwrap_or_default();
    let sv = sug["vias"].as_array().cloned().unwrap_or_default();
    let (mut bt, mut bv) = (tracks.to_vec(), vias.to_vec());
    // fixed copper of a pair member: (partner, layer or None for a via, segments)
    let mut fixed: Vec<(String, Option<String>, Vec<(P2, P2)>)> = Vec::new();
    for &i in tracks {
        if let Some(p) = partner(st[i]["net"].as_str().unwrap_or("")) {
            fixed.push((p, st[i]["layer"].as_str().map(String::from), segs(&st[i])));
        }
    }
    for &i in vias {
        if let Some(p) = partner(sv[i]["net"].as_str().unwrap_or("")) {
            fixed.push((p, None, vec![(at(&sv[i]), at(&sv[i]))]));
        }
    }
    let near = 1.0;
    for (pn, layer, fs) in &fixed {
        for (i, t) in st.iter().enumerate() {
            if bt.contains(&i) || t["net"].as_str() != Some(pn.as_str()) || (layer.is_some() && t["layer"].as_str() != layer.as_deref()) {
                continue;
            }
            if segs(t).iter().any(|&(a, b)| fs.iter().any(|&(c, d)| seg_seg_dist(a, b, c, d) < near)) {
                bt.push(i);
            }
        }
        for (i, v) in sv.iter().enumerate() {
            if !bv.contains(&i) && v["net"].as_str() == Some(pn.as_str()) && fs.iter().any(|&(c, d)| seg_point_dist(c, d, at(v)) < near) {
                bv.push(i);
            }
        }
    }
    bt.sort();
    bv.sort();
    (bt, bv)
}

fn arg(args: &[String], k: &str) -> Option<String> {
    args.iter().position(|a| a == k).and_then(|i| args.get(i + 1).cloned())
}

fn layout_json(r: &router::Router) -> serde_json::Value {
    let w = &r.w;
    let mut tracks = Vec::new();
    for pc in &w.pieces {
        if !pc.alive || pc.poly.len() < 2 {
            continue;
        }
        let mut poly = pc.poly.clone();
        // a wire ends at its pad's centre internally; into a large pad (exposed pad, power pad) it
        // stops just inside the copper instead of crossing it
        for end in [false, true] {
            let t = if end { pc.b } else { pc.a };
            let s = &w.objs[t.obj as usize].s;
            // a merged group's outline spans the gaps between its pads: only a single pad is clipped
            if t.vtx != cdt::NONE || s.kind != Kind::Pad || s.pts.len() < 3 || s.name.contains('+') || world::fixed_copper(&s.name) {
                continue;
            }
            let (lo, hi) = s.pts.iter().fold((P2 { x: 1e9, y: 1e9 }, P2 { x: -1e9, y: -1e9 }), |(lo, hi), q| (P2 { x: lo.x.min(q.x), y: lo.y.min(q.y) }, P2 { x: hi.x.max(q.x), y: hi.y.max(q.y) }));
            if (hi.x - lo.x).min(hi.y - lo.y) + 2.0 * s.r < 1.0 {
                continue;
            }
            if end {
                poly.reverse();
            }
            // walk in from the far end of the wire to where it enters the pad; keep the outside part
            // and end 0.1 mm into the copper
            let inside = |q: P2| geom::point_hull_dist(&s.pts, q) - s.r <= 0.0;
            if let Some(k) = (0..poly.len().saturating_sub(1)).rev().find(|&k| inside(poly[k]) && !inside(poly[k + 1])) {
                let (mut u, mut v) = (poly[k + 1], poly[k]);
                for _ in 0..40 {
                    let m = u.lerp(v, 0.5);
                    if inside(m) {
                        v = m
                    } else {
                        u = m
                    }
                }
                let stop = v + (poly[k] - poly[k + 1]).unit() * 0.1;
                let q = if inside(stop) { stop } else { v };
                let mut kept = vec![q];
                kept.extend_from_slice(&poly[k + 1..]);
                poly = kept;
            }
            if end {
                poly.reverse();
            }
        }
        let pts: Vec<[f64; 2]> = poly.iter().map(|p| [(p.x * 1e6).round() / 1e6, (p.y * 1e6).round() / 1e6]).collect();
        tracks.push(json!({"net": w.b.nets[pc.net as usize].name, "layer": w.b.sig_name(pc.layer as usize), "width": pc.width, "pts": pts}));
    }
    let mut vias = Vec::new();
    for o in &w.objs {
        if o.alive && o.s.kind == Kind::Via {
            vias.push(json!({"net": w.b.nets[o.s.net as usize].name, "x": (o.s.center.x * 1e6).round() / 1e6, "y": (o.s.center.y * 1e6).round() / 1e6, "d": o.s.r * 2.0, "drill": o.s.drill}));
        }
    }
    // neck-down stubs: from the pad's centre (where KiCad joins them) to their end
    for nk in &w.b.necks {
        tracks.push(json!({"net": w.b.nets[nk.net as usize].name, "layer": w.b.sig_name(nk.layer), "width": nk.width, "pts": [[nk.from.x, nk.from.y], [nk.to.x, nk.to.y]]}));
    }
    json!({"outline": w.b.size, "parts": {}, "tracks": tracks, "vias": vias, "zones": [], "keepouts": [], "keepouts_rule": [], "texts": [], "rects": [], "labels": []})
}

/// High-speed tracks over gaps in their reference planes, measured on the copper as exported (`layout`), with
/// a summary on stderr. Null on a board without high-speed nets.
fn plane_gap_report(r: &router::Router, layout: &serde_json::Value) -> serde_json::Value {
    let w = &r.w;
    if !gaps::active(w) {
        return json!(null);
    }
    let g = gaps::Gaps::new(w);
    let sig = |name: &str| (0..w.b.nsig()).find(|&l| w.b.sig_name(l) == name);
    let mut tracks = Vec::new();
    for t in layout["tracks"].as_array().into_iter().flatten() {
        let (Some(&net), Some(l)) = (t["net"].as_str().and_then(|n| w.b.net_index.get(n)), t["layer"].as_str().and_then(sig)) else { continue };
        if !w.b.high_speed(net) {
            continue;
        }
        let pts: Vec<P2> = t["pts"].as_array().into_iter().flatten().map(|q| P2 { x: q[0].as_f64().unwrap_or(0.0), y: q[1].as_f64().unwrap_or(0.0) }).collect();
        tracks.push((net, l, t["width"].as_f64().unwrap_or(0.0), pts));
    }
    let cr = g.report(&tracks);
    let lname = |p: usize| w.b.layers[p].clone();
    let num = |x: f64| if x.is_finite() { json!((x * 1000.0).round() / 1000.0) } else { json!(null) };
    let mut hs: Vec<u32> = (0..w.b.nets.len() as u32).filter(|&n| w.b.high_speed(n) && w.b.nets[n as usize].ribbon.is_none()).collect();
    hs.sort_by(|a, b| w.b.nets[*a as usize].name.cmp(&w.b.nets[*b as usize].name));
    let mut nets_json = Vec::new();
    for &n in &hs {
        let mine: Vec<&gaps::Crossing> = cr.iter().filter(|c| c.net == n && c.spans()).collect();
        let grazed: f64 = cr.iter().filter(|c| c.net == n && !c.spans()).fold(0.0, |s, c| s + c.over);
        let ni = &w.b.nets[n as usize];
        let limits: serde_json::Map<String, serde_json::Value> = (0..w.b.nsig()).filter(|&l| ni.allowed & (1 << l) != 0 && !w.b.ref_planes[l].is_empty()).map(|l| (w.b.sig_name(l).to_string(), num(w.b.gap_limit(n, l)))).collect();
        nets_json.push(json!({
            "net": ni.name, "rise_time_ns": ni.rise, "slot_limit_mm": limits, "crossings": mine.len(),
            "over_gap_mm": num(mine.iter().fold(0.0, |s, c| s + c.over)), "violations": mine.iter().filter(|c| c.violation()).count(),
            "over_void_edge_mm": num(grazed),
        }));
    }
    let viol: Vec<&gaps::Crossing> = cr.iter().filter(|c| c.violation()).collect();
    let spans: Vec<&gaps::Crossing> = cr.iter().filter(|c| c.spans()).collect();
    let grazed: f64 = cr.iter().filter(|c| !c.spans()).fold(0.0, |s, c| s + c.over);
    eprintln!("plane gaps: {} high-speed nets, {} crossings of a gap in their reference plane, {} past the rise-time limit ({:.2} mm over the edge of a via's or pad's clearance){}", hs.len(), spans.len(), viol.len(), grazed, if r.gap_rule { "" } else { "; the board routed without the gap rule completed better, so it was kept" });
    for c in viol.iter().take(40) {
        let slot = if c.slot.is_finite() { format!("{:.2} mm slot", c.slot) } else { "plane cut off".to_string() };
        eprintln!("  gap: {} on {} over {} at ({:.2},{:.2}): {}, {:.2} mm over the void, limit {:.2} mm; {}", w.b.nets[c.net as usize].name, w.b.sig_name(c.layer), lname(c.plane), c.at.x, c.at.y, slot, c.over, c.limit, c.by.join(", "));
    }
    let mut planes = serde_json::Map::new();
    let mut ps: Vec<usize> = w.b.ref_planes.iter().flatten().cloned().collect();
    ps.sort();
    ps.dedup();
    for p in ps {
        let nets: Vec<String> = w.b.nets.iter().filter(|n| n.plane == Some(p)).map(|n| n.name.clone()).chain(w.b.plane_regions.iter().filter(|x| x.layer == p).map(|x| w.b.nets[x.net as usize].name.clone())).collect::<std::collections::BTreeSet<_>>().into_iter().collect();
        planes.insert(lname(p), json!({"nets": nets, "clearance": if w.b.plane_clear[p] > 0.0 { json!(w.b.plane_clear[p]) } else { json!("class clearance") }, "min_width": w.b.plane_min_w[p]}));
    }
    json!({
        "rule": format!("a high-speed track may cross a slot in its reference plane while the slot's inductance, 5 D ln(D/W) nH (Johnson; D the slot across the track, twice the void's reach on its shorter side, W the track width), reflects at most {:.0}% of an edge of the net's rise time on a {:.0} ohm line; a run along a void counts by its length", w.b.gap_reflection * 100.0, w.b.gap_z0),
        "planes": planes,
        "references": (0..w.b.nsig()).map(|l| (w.b.sig_name(l).to_string(), json!(w.b.ref_planes[l].iter().map(|&p| lname(p)).collect::<Vec<_>>()))).collect::<serde_json::Map<String, serde_json::Value>>(),
        "rule_applied": r.gap_rule,
        "violations": viol.len(),
        "crossings_spanning_track": spans.len(),
        "over_void_edge_mm": num(grazed),
        "nets": nets_json,
        "crossings": spans.iter().map(|c| json!({
            "net": w.b.nets[c.net as usize].name, "layer": w.b.sig_name(c.layer), "plane": lname(c.plane), "at": [num(c.at.x), num(c.at.y)],
            "over_mm": num(c.over), "slot_mm": num(c.slot), "cuts_plane_off": !c.slot.is_finite(), "limit_mm": num(c.limit), "violation": c.violation(), "made_by": c.by,
        })).collect::<Vec<_>>(),
    })
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: topor build|route INPUT_DIR [--out DIR] [--svg OUT.svg] [--box x0 y0 x1 y1] [--seconds S] [--via-cost MM] [--single-attempt]");
        std::process::exit(2);
    }
    let cmd = args[1].clone();
    let dir = PathBuf::from(&args[2]);
    let t0 = Instant::now();
    let text = std::fs::read_to_string(dir.join("model.json")).expect("model.json");
    let mut m: serde_json::Value = serde_json::from_str(&text).expect("json");
    let project: Option<serde_json::Value> = std::fs::read_to_string(dir.join("unrouted.kicad_pro")).ok().and_then(|t| serde_json::from_str(&t).ok());
    let no_pairs = args.iter().any(|a| a == "--no-pairs");
    let board_of = |m: &serde_json::Value| {
        let mut b = model::Board::load(m).expect("model");
        b.detect_pairs(if no_pairs { None } else { project.as_ref() });
        if no_pairs {
            b.partner.iter_mut().for_each(|x| *x = cdt::NONE);
        }
        b
    };
    // suggested copper the router cannot take as its own (it breaks a rule, e.g. runs through a keepout drawn round a
    // crystal, or does not thread) stays as drawn: it becomes fixed copper before the board is built
    // (repeated: copper fixed in one round can make other suggested copper fail against it, e.g. of two vias drawn
    // too close the first is taken and the second fixed, and then the first fails against it; both stay as drawn.
    // Later rounds fix only vias and pair members' copper: other tracks failing then mostly overlap the fixed copper
    // at a joint, and fixing them walks the whole drawing of a power net; they are left out for the router as
    // before. A pair member's copper kept as drawn takes the partner's beside it along, so the pair stays together.)
    let mut kept_fixed = (Vec::new(), Vec::new());
    let mut drawn_pairs: Vec<String> = Vec::new();
    if !m["suggested_copper"].is_null() {
        let b0 = board_of(&m);
        let partner = |net: &str| -> Option<String> {
            let &i = b0.net_index.get(net)?;
            let &p = b0.partner.get(i as usize)?;
            (p != cdt::NONE).then(|| b0.nets[p as usize].name.clone())
        };
        for round in 0..8 {
            let w = world::World::new(board_of(&m)).expect("world");
            let cfg = router::Cfg { only: None, via_cost: 4.0, max_exp: 20000, seconds: 60.0, rip_rounds: 0, verbose: false, single_attempt: true, fast: true };
            let mut r = router::Router::new(w, cfg);
            let sg = r.import_suggestions(&m["suggested_copper"]);
            let tracks: Vec<usize> = sg.bad_tracks.iter().cloned().filter(|&i| round == 0 || partner(m["suggested_copper"]["tracks"][i]["net"].as_str().unwrap_or("")).is_some()).collect();
            if tracks.is_empty() && sg.bad_vias.is_empty() {
                break;
            }
            let (bt, bv) = with_partner_copper(&m["suggested_copper"], &tracks, &sg.bad_vias, &partner);
            let (t, v) = fix_suggestions(&mut m, &bt, &bv);
            kept_fixed.0.extend(t);
            kept_fixed.1.extend(v);
        }
        // a pair the suggested copper completes (both members' pads joined as drawn) stays as drawn too: as the
        // router's own copper it is pulled taut, pushed aside and ripped like any wire, and a pair drawn whole is
        // better kept than drawn again; the router takes the pairs drawn in part (TOPOR_ROUTE_DRAWN_PAIRS: all)
        // (TOPOR_SUGGEST_FIXED: all of it as drawn, as fixed copper given as such, to compare)
        if std::env::var("TOPOR_SUGGEST_FIXED").is_ok() {
            let all = |m: &serde_json::Value, k: &str| -> Vec<usize> { (0..m["suggested_copper"][k].as_array().map_or(0, |a| a.len())).collect() };
            let (bt, bv) = (all(&m, "tracks"), all(&m, "vias"));
            let (t, v) = fix_suggestions(&mut m, &bt, &bv);
            kept_fixed.0.extend(t);
            kept_fixed.1.extend(v);
        } else if std::env::var("TOPOR_ROUTE_DRAWN_PAIRS").is_err() {
            let w = world::World::new(board_of(&m)).expect("world");
            let cfg = router::Cfg { only: None, via_cost: 4.0, max_exp: 20000, seconds: 60.0, rip_rounds: 0, verbose: false, single_attempt: true, fast: true };
            let mut r = router::Router::new(w, cfg);
            let _ = r.import_suggestions(&m["suggested_copper"]);
            let names: Vec<String> = r.drawn_pairs.iter().map(|&n| r.w.b.nets[n as usize].name.clone()).collect();
            let sel = |m: &serde_json::Value, k: &str| -> Vec<usize> {
                m["suggested_copper"][k].as_array().into_iter().flatten().enumerate().filter(|(_, x)| names.iter().any(|n| Some(n.as_str()) == x["net"].as_str())).map(|(i, _)| i).collect()
            };
            let (bt, bv) = (sel(&m, "tracks"), sel(&m, "vias"));
            if !bt.is_empty() || !bv.is_empty() {
                let (t, v) = fix_suggestions(&mut m, &bt, &bv);
                kept_fixed.0.extend(t);
                kept_fixed.1.extend(v);
                let mut ps: Vec<String> = r.drawn_pairs.iter().filter(|&&n| n < r.w.b.partner[n as usize]).map(|&n| r.w.b.nets[n as usize].name.clone()).collect();
                ps.sort();
                drawn_pairs = ps;
            }
        }
    }
    let b = board_of(&m);
    let npairs = b.partner.iter().filter(|&&x| x != cdt::NONE).count() / 2;
    eprintln!("{} differential pairs", npairs);
    if std::env::var("TOPOR_PAIRS").is_ok() {
        for (i, &j) in b.partner.iter().enumerate() {
            if j != cdt::NONE && (i as u32) < j {
                eprintln!("  pair {} / {} gap {}", b.nets[i].name, b.nets[j as usize].name, b.pair_gap[i]);
            }
        }
    }
    let w = world::World::new(b).expect("world");
    eprintln!("world built in {:.2}s", t0.elapsed().as_secs_f64());
    if !w.b.rf_nets.is_empty() {
        eprintln!("RF frontend nets kept on one layer (no vias): {}", w.b.rf_nets.iter().map(|x| x.1.as_str()).collect::<Vec<_>>().join("; "));
    }
    if !w.b.necks.is_empty() || !w.b.neck_failed.is_empty() {
        eprintln!("neck-down stubs: {} laid{}", w.b.necks.len(), if w.b.neck_failed.is_empty() { String::new() } else { format!("; no room at {}", w.b.neck_failed.join(", ")) });
    }
    for x in w.warnings.iter().take(10) {
        eprintln!("warning: {}", x);
    }
    let bx = args.iter().position(|a| a == "--box").map(|i| [args[i + 1].parse().unwrap(), args[i + 2].parse().unwrap(), args[i + 3].parse().unwrap(), args[i + 4].parse().unwrap()]);
    let scale = arg(&args, "--scale").map(|s| s.parse().unwrap()).unwrap_or(if bx.is_some() { 60.0 } else { 8.0 });
    if cmd == "build" {
        if let Some(out) = arg(&args, "--svg") {
            let s = svg::render(&w, &svg::Opts { mesh: true, box_: bx, scale, marks: vec![], unrouted: vec![] });
            std::fs::write(&out, s).unwrap();
        }
        return;
    }
    // --fast: a quick look for layout iteration (one routing order, short rip-up, no end-game repairs or
    // polishing); the blocker report then says what keeps each open net out
    let fast = args.iter().any(|a| a == "--fast");
    let cfg = router::Cfg {
        only: arg(&args, "--only"),
        via_cost: arg(&args, "--via-cost").map(|s| s.parse().unwrap()).unwrap_or(4.0),
        // (successful searches rarely need more than 8000 expansions: the quick look stops failing ones there)
        max_exp: arg(&args, "--max-exp").map(|s| s.parse().unwrap()).unwrap_or(if fast { 8000 } else { 20000 }),
        seconds: arg(&args, "--seconds").map(|s| s.parse().unwrap()).unwrap_or(if fast { 60.0 } else { 300.0 }),
        rip_rounds: arg(&args, "--rip-rounds").map(|s| s.parse().unwrap()).unwrap_or(if fast { 12 } else { 40 }),
        verbose: true,
        single_attempt: fast || args.iter().any(|a| a == "--single-attempt"),
        fast,
    };
    let mut r = router::Router::new(w, cfg);
    // hand-drawn copper given as suggestions: routed copper to start from, not fixed copper
    let mut suggested = json!(null);
    if !m["suggested_copper"].is_null() {
        let sg = r.import_suggestions(&m["suggested_copper"]);
        eprintln!(
            "suggested copper: {} vias and {} track chains taken as routed copper; {} tracks and {} vias kept as drawn (fixed: they break a rule as routed copper, or complete a pair); left out {} vias, {} chains",
            sg.vias_kept, sg.tracks_kept, kept_fixed.0.len(), kept_fixed.1.len(), sg.vias_left_out, sg.tracks_left_out
        );
        if !drawn_pairs.is_empty() {
            eprintln!("  pairs the suggested copper completes, kept as drawn: {}", drawn_pairs.join(", "));
        }
        for n in sg.notes.iter().take(20) {
            eprintln!("  suggestion: {}", n);
        }
        suggested = json!({"vias_kept": sg.vias_kept, "chains_kept": sg.tracks_kept, "tracks_fixed": kept_fixed.0.len(), "vias_fixed": kept_fixed.1.len(), "vias_left_out": sg.vias_left_out, "chains_left_out": sg.tracks_left_out, "pairs_kept_as_drawn": drawn_pairs, "notes": sg.notes});
    }
    let t_world = t0.elapsed().as_secs_f64();
    r.route_all();
    // length matching: serpentines on the finished copper (boards with length groups only)
    let mut length_json = json!(null);
    if !r.w.b.length_groups.is_empty() && std::env::var("TOPOR_NO_TUNE").is_err() {
        let t_tune = Instant::now();
        let groups = r.match_lengths();
        r.prof_add("length matching", t_tune);
        let name = |n: u32| r.w.b.nets[n as usize].name.clone();
        let mm = |x: Option<f64>| x.map(|v| json!((v * 1e4).round() / 1e4)).unwrap_or(json!(null));
        for g in &groups {
            eprintln!(
                "length group {}: skew {:.3} -> {:.3} mm (max {:.3}), target {:.3} mm, {}{}",
                g.name,
                g.skew_before,
                g.skew_after,
                g.max_skew,
                g.target,
                if g.met { "met" } else { "NOT met" },
                if g.unknown.is_empty() { String::new() } else { format!("; not nets of the board: {}", g.unknown.join(", ")) }
            );
            for n in g.nets.iter().filter(|n| n.bumps > 0 || n.note.is_some()) {
                eprintln!(
                    "  {}: {} -> {} mm, {} bumps{}",
                    name(n.net),
                    n.before.map_or("-".into(), |v| format!("{:.3}", v)),
                    n.after.map_or("-".into(), |v| format!("{:.3}", v)),
                    n.bumps,
                    n.note.as_ref().map_or(String::new(), |s| format!(" ({})", s))
                );
            }
        }
        length_json = json!(groups
            .iter()
            .map(|g| json!({
                "name": g.name, "max_skew_mm": g.max_skew, "target_mm": (g.target * 1e4).round() / 1e4,
                "skew_before_mm": (g.skew_before * 1e4).round() / 1e4, "skew_after_mm": (g.skew_after * 1e4).round() / 1e4,
                "met": g.met, "not_nets": g.unknown,
                "nets": g.nets.iter().map(|n| json!({"net": name(n.net), "before_mm": mm(n.before), "after_mm": mm(n.after), "bumps": n.bumps, "note": n.note})).collect::<Vec<_>>(),
            }))
            .collect::<Vec<_>>());
    }
    let secs = t0.elapsed().as_secs_f64();
    let t_check = Instant::now();
    let viol = r.w.check_all();
    r.prof_add("final exact check", t_check);
    let t_blk = Instant::now();
    let blockers = r.blockers();
    r.prof_add("blocker report", t_blk);
    let xy = |p: &P2| [(p.x * 100.0).round() / 100.0, (p.y * 100.0).round() / 100.0];
    for b in &blockers {
        let open = b.open_pads.iter().take(3).map(|(n, p)| format!("{} @ ({:.2},{:.2})", n, p.x, p.y)).collect::<Vec<_>>().join(", ");
        let next = format!(
            "next to: {}{}",
            b.near_nets.iter().take(6).map(|(n, k, f)| format!("{} ({}{})", n, if *k > 0 { format!("{} routed", k) } else { String::new() }, if *f > 0 { format!("{}{} fixed", if *k > 0 { ", " } else { "" }, f) } else { String::new() })).collect::<Vec<_>>().join(", "),
            if b.near_parts.is_empty() { String::new() } else { format!("; parts {}", b.near_parts.join(", ")) }
        );
        match &b.way {
            Some(w) => {
                let layers: Vec<&str> = w.legs.iter().map(|x| x.0.as_str()).collect();
                eprintln!(
                    "  blocked: {} ({}): way through {} on {}{}; {}",
                    b.net,
                    open,
                    if w.through.is_empty() { "nothing (room the exact checks reject)".to_string() } else { w.through.iter().map(|(n, k)| format!("{} x{}", n, k)).collect::<Vec<_>>().join(", ") },
                    layers.join(" > "),
                    if w.vias.is_empty() { String::new() } else { format!(" with {} via(s) at {}", w.vias.len(), w.vias.iter().map(|v| format!("({:.2},{:.2})", v.x, v.y)).collect::<Vec<_>>().join(" ")) },
                    next
                );
            }
            None if b.sealed => eprintln!("  blocked: {} ({}): sealed in: no way even through other nets' routed copper (fixed copper, keepouts or the outline close it off; move a part or hand copper); {}", b.net, open, next),
            None => eprintln!("  blocked: {} ({}): no way found within the search budget, even through other nets' copper; {}", b.net, open, next),
        }
    }
    let blockers_json: Vec<serde_json::Value> = blockers
        .iter()
        .map(|b| {
            json!({
                "net": b.net,
                "open_pads": b.open_pads.iter().map(|(n, p)| json!({"pad": n, "at": xy(p)})).collect::<Vec<_>>(),
                "way": b.way.as_ref().map(|w| json!({
                    "through": w.through.iter().map(|(n, k)| json!({"net": n, "wires": k})).collect::<Vec<_>>(),
                    "crossings": w.crossings.iter().map(|(l, p)| json!({"layer": l, "at": xy(p)})).collect::<Vec<_>>(),
                    "legs": w.legs.iter().map(|(l, pts)| json!({"layer": l, "pts": pts.iter().map(|p| xy(p)).collect::<Vec<_>>()})).collect::<Vec<_>>(),
                    "vias": w.vias.iter().map(|p| xy(p)).collect::<Vec<_>>(),
                })),
                "sealed": b.sealed,
                "near_nets": b.near_nets.iter().map(|(n, k, f)| json!({"net": n, "routed": k, "fixed": f})).collect::<Vec<_>>(),
                "near_parts": b.near_parts,
            })
        })
        .collect();
    // where the time went: the phases of the kept attempt (the parallel attempts overlap in wall time), and the
    // searches by kind (summed over every attempt's thread)
    {
        let st = &r.stats;
        let mut ph: Vec<(&str, f64, usize)> = r.prof.iter().map(|(k, v)| (*k, v.0, v.1)).collect();
        ph.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        eprintln!("time: {:.1}s total, {:.1}s to build the board; by phase:", secs, t_world);
        for (k, t, n) in &ph {
            eprintln!("  {:>7.1}s  {:5}x  {}", t, n, k);
        }
        eprintln!(
            "  searches: {} ({} failed), {:.1}s: hard {:.1}s, soft (what to rip) {:.1}s ({} failed, {:.1}s; {} known failures skipped); failed searches {:.1}s ({} of {} expansions)",
            st.searches, st.search_fail, st.search_s, st.search_s - st.soft_s, st.soft_s, st.soft_fail, st.soft_fail_s, st.soft_skipped, st.fail_s, st.fail_exp, st.expansions
        );
    }
    let incomplete = r.incomplete_nets();
    let routable: Vec<u32> = (0..r.w.b.nets.len() as u32).filter(|&n| r.net_pads[n as usize].len() >= 2 || (r.is_plane(n) && !r.net_pads[n as usize].is_empty())).collect();
    let ntracks = r.w.pieces.iter().filter(|p| p.alive).count();
    let nvias = r.w.objs.iter().filter(|o| o.alive && o.s.kind == Kind::Via).count();
    let length: f64 = r.w.pieces.iter().filter(|p| p.alive).map(|p| p.poly.windows(2).map(|s| s[0].dist(s[1])).sum::<f64>()).sum();
    eprintln!(
        "done in {:.1}s: {}/{} nets complete, deficit {}, {} pieces, {} vias, {:.0} mm, {} violations; stats {:?}",
        secs,
        routable.len() - incomplete.iter().filter(|n| routable.contains(n)).count(),
        routable.len(),
        r.total_deficit(),
        ntracks,
        nvias,
        length,
        viol.len(),
        r.stats
    );
    if r.w.thread_drops > 0 {
        eprintln!("  {} wire(s) re-routed after a triangulation change could not re-thread them", r.w.thread_drops);
    }
    let mut by_kind: std::collections::BTreeMap<String, usize> = std::collections::BTreeMap::new();
    for v in &viol {
        *by_kind.entry(v.what.clone()).or_insert(0) += 1;
    }
    if !by_kind.is_empty() {
        eprintln!("  violations by kind: {}", by_kind.iter().map(|(k, n)| format!("{} {}", k, n)).collect::<Vec<_>>().join(", "));
    }
    if std::env::var("TOPOR_VIOLDBG").is_ok() {
        let mut groups: std::collections::BTreeMap<String, (usize, P2)> = std::collections::BTreeMap::new();
        for v in &viol {
            let mut parts: Vec<String> = v.pieces.iter().map(|&q| { let pc = &r.w.pieces[q as usize]; format!("{}{}", r.w.b.nets[pc.net as usize].name, if pc.frozen { "*" } else { "" }) }).collect();
            parts.extend(v.objs.iter().map(|&o| { let ob = &r.w.objs[o as usize]; format!("{:?}:{}", ob.s.kind, r.w.b.nets.get(ob.s.net as usize).map(|n| n.name.as_str()).unwrap_or("-")) }));
            let e = groups.entry(format!("{} L{} {}", v.what, v.layer, parts.join(" | "))).or_insert((0, v.at));
            e.0 += 1;
        }
        for pc in r.w.pieces.iter().filter(|x| x.alive && r.w.b.partner.get(x.net as usize).map_or(false, |&q| q != cdt::NONE)) {
            let len: f64 = pc.poly.windows(2).map(|s| s[0].dist(s[1])).sum();
            let on = |o: u32| { let ob = &r.w.objs[o as usize]; format!("{:?}{}@({:.2},{:.2})", ob.s.kind, if ob.s.name.is_empty() { String::new() } else { format!(" {}", ob.s.name) }, ob.s.center.x, ob.s.center.y) };
            eprintln!("  piece {} L{} {}{}: {} -> {} len {:.2}", r.w.b.nets[pc.net as usize].name, pc.layer, if pc.frozen { "frozen " } else { "" }, pc.width, on(pc.a.obj), on(pc.b.obj), len);
        }
        let mut g: Vec<_> = groups.into_iter().collect();
        g.sort_by(|a, b| b.1 .0.cmp(&a.1 .0));
        for (k, (c, at)) in g.iter().take(60) {
            eprintln!("  vgroup {:4} {} e.g. ({:.3},{:.3})", c, k, at.x, at.y);
        }
    }
    for v in viol.iter().take(10) {
        eprintln!("  violation {} layer {} at ({:.3},{:.3}) d={:.4} need {:.4}", v.what, v.layer, v.at.x, v.at.y, v.dist, v.need);
    }
    for n in incomplete.iter().take(40) {
        eprintln!("  incomplete: {} (deficit {})", r.w.b.nets[*n as usize].name, r.deficit(*n));
    }
    let pairs = r.pair_report_full();
    // per pair: what the router did and why, so the design side can act (move parts, free corridors)
    let mut pair_json = Vec::new();
    for (a, b, lp, ln, cf) in &pairs {
        eprintln!("  pair {} / {}: {:.3} / {:.3} mm, mismatch {:.3}, coupled {:.0}%", a, b, lp, ln, (lp - ln).abs(), cf * 100.0);
        let p = r.w.b.net_index[a].min(r.w.b.net_index[b]);
        let d = r.pair_diag.get(&p).cloned().unwrap_or_default();
        let mut hints: Vec<String> = Vec::new();
        if *cf < 0.6 {
            for (sa, sb, how) in &d.segments {
                if how.starts_with("members routed singly") {
                    hints.push(format!("{} -> {}: no coupled route ({}); more room at those ends (move the parts listed, or the pair's own parts closer/aligned) lets it couple", sa, sb, how));
                }
            }
            if let Some(m) = d.moved.first() {
                if d.coupled_after_pairs >= cf + 0.1 {
                    hints.push(format!("coupling fell from {:.0}% to {:.0}% in rip-up, first at {}; another route for that net (or moving its part off the pair's corridor) keeps the pair coupled", d.coupled_after_pairs * 100.0, cf * 100.0, m));
                }
            }
        }
        for (sa, sb, how) in &d.segments {
            eprintln!("     {} -> {}: {}", sa, sb, how);
        }
        for m in &d.moved {
            eprintln!("     moved: {}", m);
        }
        for x in d.notes.iter() {
            eprintln!("     note: {}", x);
        }
        for h in &hints {
            eprintln!("     hint: {}", h);
        }
        pair_json.push(json!({
            "nets": [a, b], "length_mm": [lp, ln], "mismatch_mm": (lp - ln).abs(), "coupled_fraction": cf,
            "coupled_after_pair_phase": d.coupled_after_pairs,
            "stations": d.stations,
            "segments": d.segments.iter().map(|(x, y, z)| json!({"from": x, "to": y, "routed": z})).collect::<Vec<_>>(),
            "moved_by": d.moved, "notes": d.notes, "hints": hints,
        }));
    }
    for e in &r.events {
        eprintln!("  event: {}", e);
    }
    let mut kinds: std::collections::BTreeMap<String, usize> = std::collections::BTreeMap::new();
    for l in r.log.iter() {
        let k = l.split(": ").nth(2).unwrap_or(l).split(" on layer").next().unwrap_or("").split(" at ").next().unwrap_or("").to_string();
        *kinds.entry(k).or_default() += 1;
    }
    eprintln!("  commit failure kinds: {:?}", kinds);
    if std::env::var("TOPOR_LOG").is_ok() {
        for l in r.log.iter() {
            eprintln!("  log: {}", l);
        }
    }
    let mut layout = layout_json(&r);
    layout["tracks"].as_array_mut().unwrap().extend(kept_fixed.0.iter().cloned());
    layout["vias"].as_array_mut().unwrap().extend(kept_fixed.1.iter().cloned());
    let t_gap = Instant::now();
    let plane_gaps = plane_gap_report(&r, &layout);
    if !plane_gaps.is_null() {
        eprintln!("  plane gap check {:.2}s", t_gap.elapsed().as_secs_f64());
        if std::env::var("TOPOR_GAPDBG").is_ok() {
            let c: Vec<usize> = gaps::CHARGED.iter().map(|x| x.load(std::sync::atomic::Ordering::Relaxed)).collect();
            eprintln!("  gap costs charged in the search (all attempts): {} slots crossed, {} area voids, {} slots closed by a via", c[0], c[1], c[2]);
        }
    }
    if let Some(out) = arg(&args, "--out") {
        let out = PathBuf::from(out);
        std::fs::create_dir_all(&out).unwrap();
        std::fs::write(out.join("layout.json"), serde_json::to_string(&layout).unwrap()).unwrap();
        // incomplete nets with their separate pad groups (where the gap is)
        let missing: Vec<serde_json::Value> = incomplete
            .iter()
            .map(|&n| {
                let groups: Vec<Vec<String>> = r
                    .components(n)
                    .iter()
                    .filter(|c| !r.is_plane(n) || !r.plane_ok(n, c))
                    .map(|c| c.iter().filter(|&&o| r.w.objs[o as usize].s.kind == Kind::Pad).take(6).map(|&o| { let ob = &r.w.objs[o as usize].s; format!("{} @ ({:.2},{:.2})", ob.name, ob.center.x, ob.center.y) }).collect())
                    .collect();
                json!({"net": r.w.b.nets[n as usize].name, "deficit": r.deficit(n), "unconnected_groups": groups})
            })
            .collect();
        let rj = json!({
            "case": r.w.b.name, "algorithm": "topor-engine", "cold": true, "routing_seconds": secs,
            "routable_nets": routable.len(), "complete_nets": routable.len() - incomplete.iter().filter(|n| routable.contains(n)).count(),
            "deficit": r.total_deficit(), "tracks": ntracks, "vias": nvias, "length_mm": length,
            "internal_violations": viol.len(), "internal_violations_by_kind": by_kind, "fast": fast, "blockers": blockers_json, "suggested_copper": suggested,
            "phase_seconds": r.prof.iter().map(|(k, v)| (k.to_string(), json!({"seconds": v.0, "calls": v.1}))).collect::<serde_json::Map<String, serde_json::Value>>(),
            "neck_stubs": r.w.b.necks.iter().map(|n| json!({"pad": n.pad, "net": r.w.b.nets[n.net as usize].name, "width": n.width, "length": n.from.dist(n.to)})).collect::<Vec<_>>(), "neck_failed": r.w.b.neck_failed, "missing": missing,
            "pairs": pair_json,
            "length_groups": length_json,
            "plane_gaps": plane_gaps,
            "events": r.events,
            "rf_nets_single_layer": r.w.b.rf_nets.iter().map(|x| x.1.clone()).collect::<Vec<_>>(),
            "stats": format!("{:?}", r.stats),
            "validation": "internal rules only; KiCad DRC deferred"
        });
        std::fs::write(out.join("routing.json"), serde_json::to_string_pretty(&rj).unwrap()).unwrap();
        let mut unr = Vec::new();
        for &n in &incomplete {
            let comps = r.components(n);
            for c in comps.iter().skip(if r.is_plane(n) { 0 } else { 1 }) {
                if r.is_plane(n) && r.plane_ok(n, c) {
                    continue;
                }
                unr.extend(c.iter().cloned());
            }
        }
        let s = svg::render(&r.w, &svg::Opts { mesh: args.iter().any(|a| a == "--mesh"), box_: bx, scale, marks: viol.iter().map(|v| (v.at, v.what.clone())).collect(), unrouted: unr });
        std::fs::write(out.join("routing.svg"), s).unwrap();
    }
}
