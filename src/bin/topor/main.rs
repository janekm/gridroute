//! topor: a topological (rubber-band) PCB router engine, inspired by TopoR.
mod cdt;
mod funnel;
mod fx;
mod geom;
mod model;
mod router;
mod search;
mod svg;
mod world;

use crate::geom::P2;
use crate::model::Kind;
use serde_json::json;
use std::path::PathBuf;
use std::time::Instant;

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
            if t.vtx != cdt::NONE || s.kind != Kind::Pad || s.pts.len() < 3 || s.name.contains('+') || s.name.starts_with("zone:") {
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
    json!({"outline": w.b.size, "parts": {}, "tracks": tracks, "vias": vias, "zones": [], "keepouts": [], "keepouts_rule": [], "texts": [], "rects": [], "labels": []})
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
    let m: serde_json::Value = serde_json::from_str(&text).expect("json");
    let mut b = model::Board::load(&m).expect("model");
    let project: Option<serde_json::Value> = std::fs::read_to_string(dir.join("unrouted.kicad_pro")).ok().and_then(|t| serde_json::from_str(&t).ok());
    b.detect_pairs(if args.iter().any(|a| a == "--no-pairs") { None } else { project.as_ref() });
    if args.iter().any(|a| a == "--no-pairs") {
        b.partner.iter_mut().for_each(|x| *x = cdt::NONE);
    }
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
    let cfg = router::Cfg {
        only: arg(&args, "--only"),
        via_cost: arg(&args, "--via-cost").map(|s| s.parse().unwrap()).unwrap_or(4.0),
        max_exp: arg(&args, "--max-exp").map(|s| s.parse().unwrap()).unwrap_or(20000),
        seconds: arg(&args, "--seconds").map(|s| s.parse().unwrap()).unwrap_or(300.0),
        rip_rounds: arg(&args, "--rip-rounds").map(|s| s.parse().unwrap()).unwrap_or(40),
        verbose: true,
        single_attempt: args.iter().any(|a| a == "--single-attempt"),
    };
    let mut r = router::Router::new(w, cfg);
    r.route_all();
    let secs = t0.elapsed().as_secs_f64();
    let viol = r.w.check_all();
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
    let pairs = r.pair_report();
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
    if let Some(out) = arg(&args, "--out") {
        let out = PathBuf::from(out);
        std::fs::create_dir_all(&out).unwrap();
        std::fs::write(out.join("layout.json"), serde_json::to_string(&layout_json(&r)).unwrap()).unwrap();
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
            "internal_violations": viol.len(), "missing": missing,
            "pairs": pair_json,
            "events": r.events,
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
