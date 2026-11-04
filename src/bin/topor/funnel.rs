//! Circular funnel: the shortest path through a sleeve of portals whose endpoints are discs.
//!
//! Each portal has a left and a right disc (vertex centre, clearance radius). The path keeps every left
//! disc on its left and every right disc on its right; the result is a list of wrapped discs joined by
//! common tangents, the rubber band of the topological route.

use crate::geom::{ang, p2, tangent, P2};

#[derive(Clone, Copy, Debug)]
pub struct Disc {
    pub c: P2,
    pub r: f64,
    /// Vertex identity (u32::MAX for free points such as terminals).
    pub id: u32,
}

/// A disc wrapped by the path; side +1 = disc on the left (path turns left around it).
#[derive(Clone, Copy, Debug)]
pub struct Wrap {
    pub d: Disc,
    pub side: f64,
}

#[derive(Clone, Copy)]
struct El {
    d: Disc,
    side: f64,
    idx: usize,
}

fn same(a: &El, b: &El) -> bool {
    a.d.id != u32::MAX && a.d.id == b.d.id && a.side == b.side
}

/// Direction of the tangent from apex to x, or None when x is the apex disc itself.
fn dir(apex: &El, x: &El, overlap: &mut Option<P2>) -> Option<P2> {
    if same(apex, x) {
        return None;
    }
    if apex.d.c.dist(x.d.c) < 1e-12 && x.d.id == u32::MAX {
        return None;
    }
    // a disc nested inside the apex disc on the same side is cleared by wrapping the apex
    let dd = apex.d.c.dist(x.d.c);
    if apex.side == x.side && dd + x.d.r <= apex.d.r + 1e-6 {
        return None;
    }
    match tangent(apex.d.c, apex.d.r, apex.side, x.d.c, x.d.r, x.side) {
        Some((u, _, _)) => Some(u),
        None => {
            if std::env::var("TOPOR_DEBUG").is_ok() {
                eprintln!("    tangent fails: apex {:?} r {:.4} s {} id {} -> {:?} r {:.4} s {} id {}", apex.d.c, apex.d.r, apex.side, apex.d.id, x.d.c, x.d.r, x.side, x.d.id);
            }
            if overlap.is_none() {
                *overlap = Some(apex.d.c.lerp(x.d.c, 0.5));
            }
            let v = x.d.c - apex.d.c;
            if v.norm() < 1e-12 {
                None
            } else {
                Some(v.unit())
            }
        }
    }
}

/// Replace each portal disc that lies inside another disc on the same side nearby in the sleeve.
fn dominate(portals: &mut [(Disc, Disc)]) {
    let n = portals.len();
    for side in 0..2 {
        let get = |p: &(Disc, Disc)| if side == 0 { p.0 } else { p.1 };
        for _ in 0..3 {
            let mut changed = false;
            for i in 0..n {
                let a = get(&portals[i]);
                let lo = i.saturating_sub(6);
                let hi = (i + 7).min(n);
                for j in lo..hi {
                    let b = get(&portals[j]);
                    if b.id == a.id {
                        continue;
                    }
                    if a.c.dist(b.c) + a.r <= b.r + 1e-6 && a.r <= b.r {
                        if side == 0 {
                            portals[i].0 = b;
                        } else {
                            portals[i].1 = b;
                        }
                        changed = true;
                        break;
                    }
                }
            }
            if !changed {
                break;
            }
        }
    }
}

/// Length of the tangent from apex to x (distance to where the path would touch x).
fn tlen(apex: &El, x: &El) -> f64 {
    match tangent(apex.d.c, apex.d.r, apex.side, x.d.c, x.d.r, x.side) {
        Some((_, a, b)) => a.dist(b),
        None => apex.d.c.dist(x.d.c),
    }
}

/// Returns the wrapped discs between start and end, and where a tangent was undefined (overlap), if anywhere.
pub fn funnel(start: P2, portals_in: &[(Disc, Disc)], end: P2) -> (Vec<Wrap>, Option<P2>) {
    let mut pv = portals_in.to_vec();
    // a terminal inside a clearance disc: the stack there does not reach the terminal, shrink it
    for (a, b) in pv.iter_mut() {
        for d in [a, b] {
            for t in [start, end] {
                let dd = d.c.dist(t);
                if dd < d.r {
                    d.r = (dd - 1e-6).max(0.0);
                }
            }
        }
    }
    dominate(&mut pv);
    let portals = &pv[..];
    let pt = |p: P2| Disc { c: p, r: 0.0, id: u32::MAX };
    let n = portals.len() + 2;
    // element i: portal i-1 (i = 1..n-2); 0 = start, n-1 = end
    let left_of = |i: usize| -> El {
        if i == 0 {
            El { d: pt(start), side: 0.0, idx: 0 }
        } else if i == n - 1 {
            El { d: pt(end), side: 1.0, idx: i }
        } else {
            El { d: portals[i - 1].0, side: 1.0, idx: i }
        }
    };
    let right_of = |i: usize| -> El {
        if i == 0 {
            El { d: pt(start), side: 0.0, idx: 0 }
        } else if i == n - 1 {
            El { d: pt(end), side: -1.0, idx: i }
        } else {
            El { d: portals[i - 1].1, side: -1.0, idx: i }
        }
    };
    let mut out: Vec<Wrap> = Vec::new();
    let mut overlap: Option<P2> = None;
    let mut apex = left_of(0);
    let mut left = apex;
    let mut right = apex;
    let is_apex = |x: &El, apex: &El| x.idx == apex.idx || same(x, apex);
    let cr = |a: Option<P2>, b: Option<P2>| -> f64 {
        match (a, b) {
            (Some(a), Some(b)) => a.cross(b),
            _ => 0.0,
        }
    };
    let mut i = 1;
    let mut guard = 0usize;
    while i < n {
        guard += 1;
        if guard > 50 * n + 1000 {
            overlap = Some(start);
            break;
        }
        let l = left_of(i);
        let r = right_of(i);
        let da = |x: &El, apex: &El, ov: &mut Option<P2>| if is_apex(x, apex) { None } else { dir(apex, x, ov) };
        // right side
        let d_right = da(&right, &apex, &mut overlap);
        let d_r = da(&r, &apex, &mut overlap);
        if cr(d_right, d_r) >= 0.0 {
            let d_left = da(&left, &apex, &mut overlap);
            if is_apex(&right, &apex) || cr(d_left, d_r) < 0.0 {
                right = r;
            } else {
                // the funnel closed: wrap whichever boundary disc the path reaches first
                let pick = if tlen(&apex, &r) < tlen(&apex, &left) - 1e-9 { r } else { left };
                out.push(Wrap { d: pick.d, side: pick.side });
                apex = pick;
                left = apex;
                right = apex;
                i = apex.idx;
                continue;
            }
        }
        // left side
        let d_left = da(&left, &apex, &mut overlap);
        let d_l = da(&l, &apex, &mut overlap);
        if cr(d_left, d_l) <= 0.0 {
            let d_right = da(&right, &apex, &mut overlap);
            if is_apex(&left, &apex) || cr(d_right, d_l) > 0.0 {
                left = l;
            } else {
                let pick = if tlen(&apex, &l) < tlen(&apex, &right) - 1e-9 { l } else { right };
                out.push(Wrap { d: pick.d, side: pick.side });
                apex = pick;
                left = apex;
                right = apex;
                i = apex.idx;
                continue;
            }
        }
        i += 1;
    }
    // merge repeated wraps of one disc
    let mut merged: Vec<Wrap> = Vec::new();
    for w in out {
        if w.d.id == u32::MAX {
            // start or end point reached as an apex
            continue;
        }
        if let Some(last) = merged.last_mut() {
            if (last.d.id == w.d.id || last.d.c.dist(w.d.c) < 1e-6) && last.side == w.side {
                last.d.r = last.d.r.max(w.d.r);
                continue;
            }
        }
        merged.push(w);
    }
    (merged, overlap)
}

/// Polyline for start -> wraps -> end, arcs discretized outward by `tol`. The flag reports a wrap that
/// the path turns against (the geometry cuts through that disc).
pub fn polyline(start: P2, wraps: &[Wrap], end: P2, tol: f64) -> (Vec<P2>, bool) {
    let mut seq: Vec<Wrap> = Vec::with_capacity(wraps.len() + 2);
    seq.push(Wrap { d: Disc { c: start, r: 0.0, id: u32::MAX }, side: 0.0 });
    seq.extend_from_slice(wraps);
    seq.push(Wrap { d: Disc { c: end, r: 0.0, id: u32::MAX }, side: 0.0 });
    // tangent between consecutive elements: (direction, departure point, arrival point)
    let mut tan: Vec<(P2, P2, P2)> = Vec::with_capacity(seq.len());
    let mut bad = false;
    for i in 0..seq.len() - 1 {
        let (a, b) = (seq[i], seq[i + 1]);
        match tangent(a.d.c, a.d.r, a.side, b.d.c, b.d.r, b.side) {
            Some(t) => tan.push(t),
            None => {
                bad = true;
                let u = (b.d.c - a.d.c).unit();
                tan.push((u, a.d.c + u * a.d.r, b.d.c - u * b.d.r));
            }
        }
    }
    let mut pts = vec![start];
    let push = |pts: &mut Vec<P2>, p: P2| {
        if pts.last().map_or(true, |q| q.dist(p) > 1e-7) {
            pts.push(p);
        }
    };
    for i in 1..seq.len() - 1 {
        let w = seq[i];
        let (u_in, _, arr) = tan[i - 1];
        let (u_out, dep, _) = tan[i];
        let turn = u_in.cross(u_out).atan2(u_in.dot(u_out));
        if turn * w.side < -1e-6 && w.d.r > 1e-3 {
            bad = true;
        }
        if turn * w.side <= 1e-9 || w.d.r <= 1e-9 {
            push(&mut pts, arr.lerp(dep, 0.5));
            continue;
        }
        let c = w.d.c;
        let ro = w.d.r + tol;
        let step = 2.0 * (w.d.r / ro).acos();
        let span = turn.abs();
        let n = ((span / step).ceil() as usize).max(1);
        let a0 = ang(arr - c);
        for k in 0..=n {
            let t = a0 + w.side * span * (k as f64) / (n as f64);
            push(&mut pts, c + p2(t.cos(), t.sin()) * ro);
        }
    }
    push(&mut pts, end);
    if pts.len() == 1 {
        pts.push(end);
    }
    (pts, bad)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn d(x: f64, y: f64, r: f64, id: u32) -> Disc {
        Disc { c: p2(x, y), r, id }
    }

    #[test]
    fn straight_corridor() {
        // portals across a horizontal corridor; straight path, no wraps
        let portals: Vec<(Disc, Disc)> = (1..5).map(|i| (d(i as f64, 1.0, 0.1, i * 2), d(i as f64, -1.0, 0.1, i * 2 + 1))).collect();
        let (w, ov) = funnel(p2(0., 0.), &portals, p2(6., 0.));
        assert!(w.is_empty(), "{:?}", w);
        assert!(ov.is_none());
    }

    #[test]
    fn wraps_obstacle() {
        // path must pass below a disc at (3, 0.2) that is the left endpoint of a portal (disc on left side)
        // portal at x=3: left disc (3,0.2) r .3 ... right disc (3,-3)
        let portals = vec![(d(3.0, 0.2, 0.3, 1), d(3.0, -3.0, 0.1, 2))];
        // travelling +x, left side is +y. Left disc must be on the path's left => path below y=0.2-0.3
        let (w, _) = funnel(p2(0., 0.), &portals, p2(6., 0.));
        assert_eq!(w.len(), 1);
        assert_eq!(w[0].d.id, 1);
        assert!(w[0].side > 0.0);
        let (pl, _) = polyline(p2(0., 0.), &w, p2(6., 0.), 0.001);
        for q in &pl {
            assert!(q.dist(p2(3.0, 0.2)) >= 0.3 - 1e-9, "{:?}", q);
        }
    }

    #[test]
    fn zigzag() {
        // left obstacle low at x=2, right obstacle high at x=4: s-shaped path wrapping both
        let portals = vec![(d(2.0, -0.5, 0.2, 1), d(2.0, -5.0, 0.1, 2)), (d(4.0, 5.0, 0.1, 3), d(4.0, 0.5, 0.2, 4))];
        let (w, _) = funnel(p2(0., 0.), &portals, p2(6., 0.));
        assert_eq!(w.len(), 2, "{:?}", w);
        assert_eq!(w[0].d.id, 1);
        assert_eq!(w[1].d.id, 4);
        let (pl, _) = polyline(p2(0., 0.), &w, p2(6., 0.), 0.001);
        for q in &pl {
            assert!(q.dist(p2(2.0, -0.5)) >= 0.2 - 1e-9);
            assert!(q.dist(p2(4.0, 0.5)) >= 0.2 - 1e-9);
        }
    }
}

#[cfg(test)]
mod regress {
    use super::*;
    use crate::geom::p2;

    fn d(x: f64, y: f64, r: f64, id: u32) -> Disc {
        Disc { c: p2(x, y), r, id }
    }

    fn check(start: P2, portals: &[(Disc, Disc)], end: P2) {
        let (w, ov) = funnel(start, portals, end);
        let (pl, bad) = polyline(start, &w, end, 0.0007);
        for (l, r) in portals {
            for q in pl.windows(2) {
                let dl = crate::geom::seg_point_dist(q[0], q[1], l.c);
                let dr = crate::geom::seg_point_dist(q[0], q[1], r.c);
                assert!(dl >= l.r - 1e-6, "left disc {:?} hit: {} < {} wraps {:?} poly {:?}", l.c, dl, l.r, w, pl);
                assert!(dr >= r.r - 1e-6, "right disc {:?} hit: {} < {} wraps {:?}", r.c, dr, r.r, w);
            }
        }
        assert!(ov.is_none() && !bad, "overlap {:?} bad {} wraps {:?}", ov, bad, w);
    }

    #[test]
    fn s_bend_between_pads() {
        let portals = vec![
            (d(19.05, 16.51, 1.2093, 0x73), d(19.05, 19.05, 1.2093, 0x75)),
            (d(19.05, 16.51, 1.2093, 0x73), d(15.459, 21.1, 0.8195, 0xa8)),
            (d(19.05, 16.51, 1.2093, 0x73), d(12.671159, 16.098515, 0.7593, 0xc3)),
            (d(19.05, 13.97, 1.6695, 0x71), d(12.671159, 16.098515, 0.7593, 0xc3)),
            (d(18.2, 12.28, 0.8195, 0x6f), d(12.671159, 16.098515, 0.7593, 0xc3)),
            (d(11.903, 13.48, 0.0012, 0x8c), d(12.671159, 16.098515, 0.7593, 0xc3)),
            (d(11.903, 13.48, 0.0012, 0x8c), d(11.903, 14.46, 0.0012, 0x8d)),
        ];
        check(p2(21.59, 16.51), &portals, p2(11.438, 13.97));
    }
}
