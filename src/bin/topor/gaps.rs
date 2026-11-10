//! Return-path gaps: high-speed tracks (nets with a rise time) over voids in their reference plane.
//!
//! A plane's voids are the through copper of other nets (vias, through-hole pads) and holes, each grown by the
//! pour clearance; the pull-back from the board edge; copper webs narrower than the pour's least width (they are
//! not poured, so neighbouring clearances merge into one slot); fill keepouts; and, on a plane of partial pours,
//! everywhere outside them. A track that crosses a void sends its return current around it. Taking D as twice
//! the void's reach across the track on its shorter side, the detour adds about 5 D ln(D/W) nH (Johnson, High-Speed
//! Digital Design, slots in ground planes). `Board::gap_limit` gives the longest slot a net's rise time allows.
//!
//! The search charges a high-speed net for each crossing between two through objects whose voids merge into a
//! slot longer than its limit, and for stepping into a fill keepout or across a split between pours; any net pays
//! for a via that would close such a slot across a high-speed wire. The board is also routed without the rule
//! and the rule's board kept only if it is as complete (router::route_all). The report measures the exported
//! copper, sample by sample.
//!
//! Inputs, in model.json `cold_routing`:
//! - `rise_time`: {net or class: ns}. Nets with one are high-speed; nothing changes on a board without.
//! - `plane_clearance`, `plane_min_width`: the pours' clearance and least width (mm; one number or {layer: mm}),
//!   as bench/pcb/topor_planes.py reads them from the KiCad zones. By default each void is the other copper's
//!   own clearance and the least width KiCad's default 0.25 mm.
//! - `reference_planes`: {signal layer: [plane layers]}, where the adjacent planes are not the reference.
//! - `gap_reflection` (0.05) and `gap_z0` (50 ohm): the limit's reflection and line impedance.
//!
//! TOPOR_GAP_COST sets the charge (25 per slot, plus 4 per mm over the limit; 0: report only), TOPOR_GAP_ONLY
//! routes with the rule alone (no board without it), TOPOR_GAPDBG counts what the search charged.

use crate::cdt::NONE;
use crate::fx::FxMap;
use crate::geom::*;
use crate::model::Kind;
use crate::world::{World, IT_OBJ, IT_RING, IT_SEG, MARGIN};
use std::cell::RefCell;

/// Distance between two convex cores (1 point: a disc's centre, 2: a segment, 3 or more: a convex polygon).
fn core_gap(a: &[P2], b: &[P2]) -> f64 {
    let mut d = f64::INFINITY;
    for &q in a {
        d = d.min(point_hull_dist(b, q));
    }
    for &q in b {
        d = d.min(point_hull_dist(a, q));
    }
    if a.len() >= 2 && b.len() >= 2 {
        let na = if a.len() == 2 { 1 } else { a.len() };
        let nb = if b.len() == 2 { 1 } else { b.len() };
        for k in 0..na {
            for m in 0..nb {
                d = d.min(seg_seg_dist(a[k], a[(k + 1) % a.len()], b[m], b[(m + 1) % b.len()]));
            }
        }
    }
    d
}

/// A point of segment a-b inside polygon `poly`: b or a if inside, else just past the first edge crossing.
fn seg_enters(poly: &[P2], a: P2, b: P2) -> Option<P2> {
    if point_in_poly(poly, b) {
        return Some(b);
    }
    if point_in_poly(poly, a) {
        return Some(a);
    }
    let n = poly.len();
    let mut best: Option<f64> = None;
    for i in 0..n {
        if let Some((t, _)) = seg_cross_params(a, b, poly[i], poly[(i + 1) % n]) {
            best = Some(best.map_or(t, |x: f64| x.min(t)));
        }
    }
    best.map(|t| a.lerp(b, (t + 0.01 / a.dist(b).max(0.01)).min(1.0)))
}

/// How often the search charged a slot crossed, an area void, a slot closed by a via (TOPOR_GAPDBG prints them).
pub static CHARGED: [std::sync::atomic::AtomicUsize; 3] = [std::sync::atomic::AtomicUsize::new(0), std::sync::atomic::AtomicUsize::new(0), std::sync::atomic::AtomicUsize::new(0)];
fn charged(k: usize) {
    CHARGED[k].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
}

/// Search cost of crossing a slot `excess` mm longer than the net allows (TOPOR_GAP_COST: the base, 25).
fn pen(excess: f64) -> f64 {
    static BASE: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    let b = *BASE.get_or_init(|| std::env::var("TOPOR_GAP_COST").ok().and_then(|x| x.parse().ok()).unwrap_or(25.0));
    if b <= 0.0 {
        return 0.0;
    }
    b + 4.0 * excess.clamp(0.0, 20.0)
}

/// Nets whose own copper never makes a void for a wire of `net`: the net itself, and for a ribbon its members.
pub fn own_of(w: &World, net: u32) -> Vec<u32> {
    match w.b.nets.get(net as usize).and_then(|n| n.ribbon) {
        Some((p, n)) => vec![net, p, n],
        None => vec![net],
    }
}

/// Whether the board asks for the gap rule: a high-speed net and a signal layer with a reference plane.
pub fn active(w: &World) -> bool {
    w.b.nets.iter().any(|n| n.rise > 0.0) && w.b.ref_planes.iter().any(|r| !r.is_empty())
}

#[derive(Default, Clone)]
struct Clus {
    members: Vec<u32>,
    /// Where the cluster's void meets the pull-back from the board edge.
    edge: Vec<P2>,
}

#[derive(Default)]
struct Clusters {
    of: FxMap<u32, u32>,
    list: Vec<Clus>,
}

/// One crossing of a high-speed track over a void in its reference plane.
#[derive(Clone, Debug)]
pub struct Crossing {
    pub net: u32,
    /// Signal layer index, and the plane's index in `Board::layers`.
    pub layer: usize,
    pub plane: usize,
    pub at: P2,
    /// Track length over the void (mm).
    pub over: f64,
    /// The slot across the track (mm, infinite where the void cuts the plane off); the net's limit.
    pub slot: f64,
    pub limit: f64,
    /// The track's width: a void reaching less far across the track than this only grazes it.
    pub width: f64,
    /// What makes the void (nets of vias, pads, "board edge", "fill keepout", "no pour").
    pub by: Vec<String>,
}

impl Crossing {
    /// The void spans the track (rather than reaching under its centre from one side).
    pub fn spans(&self) -> bool {
        self.slot > self.width
    }
    pub fn violation(&self) -> bool {
        self.slot.max(if self.spans() { self.over } else { 0.0 }) > self.limit
    }
}

pub struct Gaps<'a> {
    w: &'a World,
    /// Nets whose copper joins each plane layer (their through copper makes no void there).
    plane_nets: Vec<Vec<u32>>,
    /// Plane layers some signal layer references.
    planes: Vec<usize>,
    /// Largest growth of a copper core into its void, for grid queries.
    grow_max: f64,
    cache: RefCell<FxMap<(usize, u32), Clusters>>,
    /// Via costs by (net, position).
    vcache: RefCell<FxMap<(u32, i64, i64), f64>>,
}

impl<'a> Gaps<'a> {
    pub fn new(w: &'a World) -> Gaps<'a> {
        let b = &w.b;
        let mut plane_nets = vec![Vec::new(); b.layers.len()];
        for (i, n) in b.nets.iter().enumerate() {
            if let Some(p) = n.plane {
                plane_nets[p].push(i as u32);
            }
        }
        for r in &b.plane_regions {
            if !plane_nets[r.layer].contains(&r.net) {
                plane_nets[r.layer].push(r.net);
            }
        }
        let mut planes: Vec<usize> = b.ref_planes.iter().flatten().cloned().collect();
        planes.sort();
        planes.dedup();
        let mut grow_max = b.hole_clear.max(b.edge_clear);
        for c in &b.classes {
            grow_max = grow_max.max(c.clear);
        }
        for &p in &planes {
            grow_max = grow_max.max(b.plane_clear[p]);
        }
        for o in &w.objs {
            grow_max = grow_max.max(o.s.clear);
        }
        Gaps { w, plane_nets, planes, grow_max: grow_max + 0.05, cache: RefCell::new(FxMap::default()), vcache: RefCell::new(FxMap::default()) }
    }

    fn mw(&self, p: usize) -> f64 {
        self.w.b.plane_min_w[p]
    }

    /// Growth of copper of `net` (clearance `own_clear`) into the void it leaves in plane p.
    fn grow(&self, net: u32, own_clear: f64, p: usize) -> f64 {
        let cc = if net != NONE { self.w.b.classes[self.w.b.nets[net as usize].class].clear } else { 0.0 };
        own_clear.max(cc).max(self.w.b.plane_clear[p])
    }

    /// Radius of object o's void around its core in plane p, if it leaves one there (through copper of a net
    /// other than the plane's and `own`, or a hole).
    pub fn void_r(&self, o: u32, p: usize, own: &[u32]) -> Option<f64> {
        let ob = &self.w.objs[o as usize];
        if !ob.alive {
            return None;
        }
        let s = &ob.s;
        let through = matches!(s.kind, Kind::Via | Kind::Hole) || (s.kind == Kind::Pad && s.through);
        if !through || (s.net != NONE && (own.contains(&s.net) || self.plane_nets[p].contains(&s.net))) {
            return None;
        }
        Some(s.r + self.grow(s.net, s.clear, p))
    }

    /// Objects whose voids in plane p merge with the void (core `pts`, radius r): their gap is under the pour's
    /// least width. Also where it merges with the pull-back from the board edge.
    fn merging(&self, pts: &[P2], r: f64, p: usize, own: &[u32], skip: u32, out: &mut Vec<u32>, edge: &mut Vec<P2>) {
        let w = self.w;
        let mw = self.mw(p);
        let (mut lo, mut hi) = (p2(f64::INFINITY, f64::INFINITY), p2(f64::NEG_INFINITY, f64::NEG_INFINITY));
        for q in pts {
            lo = p2(lo.x.min(q.x), lo.y.min(q.y));
            hi = p2(hi.x.max(q.x), hi.y.max(q.y));
        }
        let m = r + mw + self.grow_max + MARGIN;
        let mut items = Vec::new();
        w.layers[0].grid.query(lo - p2(m, m), hi + p2(m, m), &mut items);
        for it in items {
            if it.id == skip {
                continue;
            }
            match it.kind {
                IT_OBJ => {
                    if out.contains(&it.id) {
                        continue;
                    }
                    if let Some(r2) = self.void_r(it.id, p, own) {
                        if core_gap(pts, &w.objs[it.id as usize].s.pts) - r - r2 < mw {
                            out.push(it.id);
                        }
                    }
                }
                IT_RING => {
                    let ob = &w.objs[it.id as usize];
                    if ob.s.kind != Kind::Outline {
                        continue;
                    }
                    let n = ob.s.pts.len();
                    let (a, c) = (ob.s.pts[it.sub as usize], ob.s.pts[(it.sub as usize + 1) % n]);
                    let d = pts.iter().map(|&q| seg_point_dist(a, c, q)).fold(f64::INFINITY, f64::min);
                    if d - r - w.b.edge_clear < mw {
                        let q = pts.iter().cloned().min_by(|x, y| seg_point_dist(a, c, *x).partial_cmp(&seg_point_dist(a, c, *y)).unwrap()).unwrap();
                        edge.push(a.lerp(c, seg_param(a, c, q)));
                    }
                }
                _ => {}
            }
        }
    }

    /// The void cluster (objects whose voids merge, in turn) of object o in plane p, without `own`'s copper.
    fn cluster(&self, p: usize, own: &[u32], o: u32) -> Clus {
        let key = (p, own[0]);
        if let Some(c) = self.cache.borrow().get(&key).and_then(|cs| cs.of.get(&o).map(|&i| cs.list[i as usize].clone())) {
            return c;
        }
        let mut members = vec![o];
        let mut edge = Vec::new();
        let mut k = 0;
        while k < members.len() && members.len() < 4000 {
            let x = members[k];
            k += 1;
            let r = match self.void_r(x, p, own) {
                Some(r) => r,
                None => continue,
            };
            let mut nb = Vec::new();
            self.merging(&self.w.objs[x as usize].s.pts, r, p, own, x, &mut nb, &mut edge);
            for y in nb {
                if !members.contains(&y) {
                    members.push(y);
                }
            }
        }
        let c = Clus { members, edge };
        let mut cache = self.cache.borrow_mut();
        let cs = cache.entry(key).or_default();
        let i = cs.list.len() as u32;
        for &m in &c.members {
            cs.of.insert(m, i);
        }
        cs.list.push(c.clone());
        c
    }

    /// The slot across a track at m (unit normal nrm) through a void made of `members` (objects), `extra`
    /// elements (core, radius) and edge contacts: twice the void's reach on its shorter side; infinite where the
    /// void meets the board edge on both sides.
    fn slot(&self, p: usize, own: &[u32], members: &[u32], extra: &[(Vec<P2>, f64)], edge: &[P2], m: P2, nrm: P2) -> f64 {
        let (mut lo, mut hi) = (0.0f64, 0.0f64);
        let mut add = |pts: &[P2], r: f64| {
            for &q in pts {
                let s = (q - m).dot(nrm);
                lo = lo.min(s - r);
                hi = hi.max(s + r);
            }
        };
        for &o in members {
            if let Some(r) = self.void_r(o, p, own) {
                add(&self.w.objs[o as usize].s.pts, r);
            }
        }
        for (pts, r) in extra {
            add(pts, *r);
        }
        for &e in edge {
            if (e - m).dot(nrm) < 0.0 {
                lo = f64::NEG_INFINITY;
            } else {
                hi = f64::INFINITY;
            }
        }
        2.0 * (-lo).min(hi)
    }

    /// Search cost for a high-speed wire of `net` on signal layer l crossing the triangulation edge between
    /// objects a and b at `at`, heading `dir`: a slot over its limit in a reference plane.
    pub fn edge_cost(&self, net: u32, own: &[u32], l: usize, a: u32, b: u32, at: P2, dir: P2) -> f64 {
        let w = self.w;
        if a == NONE || b == NONE || a == b || dir.norm2() < 1e-12 {
            return 0.0;
        }
        let nrm = dir.unit().perp();
        let mut cost = 0.0;
        for &p in &w.b.ref_planes[l] {
            let outline = |o: u32| w.objs[o as usize].s.kind == Kind::Outline;
            let (ra, rb) = (self.void_r(a, p, own), self.void_r(b, p, own));
            // a void object and the board edge, or two void objects whose voids merge
            let seed = match (ra, rb) {
                (Some(r1), Some(r2)) => {
                    if core_gap(&w.objs[a as usize].s.pts, &w.objs[b as usize].s.pts) - r1 - r2 >= self.mw(p) {
                        continue;
                    }
                    a
                }
                (Some(_), None) if outline(b) => a,
                (None, Some(_)) if outline(a) => b,
                _ => continue,
            };
            // between a void and the board edge, the void itself must meet the edge's pull-back there
            if outline(a) || outline(b) {
                let (mut nb, mut e) = (Vec::new(), Vec::new());
                let r = self.void_r(seed, p, own).unwrap_or(0.0);
                self.merging(&w.objs[seed as usize].s.pts, r, p, own, seed, &mut nb, &mut e);
                if e.is_empty() {
                    continue;
                }
            }
            let c = self.cluster(p, own, seed);
            // crossing between the void and the board edge: the edge closes the side away from the void
            let mut edge = c.edge.clone();
            if outline(a) || outline(b) {
                edge.push(at + nrm * if (w.objs[seed as usize].s.center - at).dot(nrm) > 0.0 { -1.0 } else { 1.0 });
            }
            let d = self.slot(p, own, &c.members, &[], &edge, at, nrm);
            let lim = w.b.gap_limit(net, l);
            if d > lim {
                cost += pen(d - lim);
                charged(0);
            }
        }
        cost
    }

    /// Search cost for a high-speed wire of `net` on signal layer l stepping from `from` to `to` over an area
    /// void of a reference plane: into a fill keepout whose slot across the step is over the limit, or from one
    /// pour to another (a split).
    pub fn area_cost(&self, net: u32, l: usize, from: P2, to: P2) -> f64 {
        let b = &self.w.b;
        let dir = to - from;
        if dir.norm2() < 1e-12 {
            return 0.0;
        }
        let nrm = dir.unit().perp();
        let mut cost = 0.0;
        for &p in &b.ref_planes[l] {
            let lim = b.gap_limit(net, l);
            for (kl, k) in &b.fill_keepouts {
                if *kl != p {
                    continue;
                }
                if let Some(q) = seg_enters(k, from, to) {
                    let d = self.area_slot(p, q, nrm);
                    if d > lim {
                        cost += pen(d - lim);
                        charged(1);
                    }
                }
            }
            if self.splits(p, from, to) {
                cost += pen(20.0);
                charged(1);
            }
        }
        cost
    }

    /// Whether segment a-b crosses a split in plane p (a plane of partial pours): from one net's pour to another's,
    /// or out of every pour away from the board edge (whose pull-back the outline already keeps). The net under
    /// the segment is read between each pair of its crossings with the pours' outlines.
    fn splits(&self, p: usize, a: P2, b: P2) -> bool {
        let bd = &self.w.b;
        let mut ts: Vec<f64> = vec![0.0, 1.0];
        let mut any = false;
        for r in bd.plane_regions.iter().filter(|r| r.layer == p) {
            any = true;
            for ring in std::iter::once(&r.poly).chain(r.holes.iter()) {
                let n = ring.len();
                for i in 0..n {
                    if let Some((t, _)) = seg_cross_params(a, b, ring[i], ring[(i + 1) % n]) {
                        ts.push(t);
                    }
                }
            }
        }
        if !any || ts.len() == 2 {
            return false;
        }
        ts.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let edge_d = |q: P2| bd.outline_holes.iter().map(|h| poly_edge_dist(h, q)).fold(poly_edge_dist(&bd.outline, q), f64::min);
        let mut last: Option<Option<u32>> = None;
        for w2 in ts.windows(2) {
            if w2[1] - w2[0] < 1e-9 {
                continue;
            }
            let q = a.lerp(b, 0.5 * (w2[0] + w2[1]));
            let n = self.region_at(p, q).map(|i| bd.plane_regions[i].net);
            if n.is_none() && edge_d(q) <= bd.edge_clear + self.mw(p) {
                continue;
            }
            if let Some(l) = last {
                if l != n {
                    return true;
                }
            }
            last = Some(n);
        }
        false
    }

    /// The highest-priority pour of plane p holding q (an index into `Board::plane_regions`), if any.
    fn region_at(&self, p: usize, q: P2) -> Option<usize> {
        let b = &self.w.b;
        (0..b.plane_regions.len())
            .filter(|&i| {
                let r = &b.plane_regions[i];
                r.layer == p && point_in_poly(&r.poly, q) && !r.holes.iter().any(|h| point_in_poly(h, q))
            })
            .max_by_key(|&i| b.plane_regions[i].priority)
    }

    /// The slot across a track at q (unit normal nrm) through the area voids of plane p: twice the void's reach
    /// on its shorter side, stepped (infinite past 50 mm both ways).
    fn area_slot(&self, p: usize, q: P2, nrm: P2) -> f64 {
        let h = self.mw(p) / 2.0;
        let reach = |s: f64| -> f64 {
            let mut t = 0.0;
            while t < 50.0 && self.area_clear(p, q + nrm * (s * t)) < h {
                t += 0.05;
            }
            t
        };
        let (r1, r2) = (reach(1.0), reach(-1.0));
        if r1 >= 50.0 && r2 >= 50.0 {
            f64::INFINITY
        } else {
            2.0 * r1.min(r2)
        }
    }

    /// Search cost of a via of `net` at x: a slot over its limit it would close across a high-speed wire of
    /// another net (its void merging with a void, or the edge pull-back, on the wire's other side).
    pub fn via_cost(&self, net: u32, x: P2) -> f64 {
        let key = (net, (x.x * 2e3).round() as i64, (x.y * 2e3).round() as i64);
        if let Some(&c) = self.vcache.borrow().get(&key) {
            return c;
        }
        let c = self.via_cost_at(net, x);
        self.vcache.borrow_mut().insert(key, c);
        c
    }

    fn via_cost_at(&self, net: u32, x: P2) -> f64 {
        let w = self.w;
        let (d, _) = w.b.via_size(net);
        let mut cost = 0.0;
        for &p in &self.planes {
            if self.plane_nets[p].contains(&net) {
                continue;
            }
            let rv = d / 2.0 + self.grow(net, w.b.via_clear(net), p);
            let mw = self.mw(p);
            let reach = rv + mw + 2.0 * self.grow_max + w.max_hw + 1.0;
            for l in 0..w.nsig() {
                if !w.b.ref_planes[l].contains(&p) {
                    continue;
                }
                let mut segs: Vec<(u32, P2, P2)> = Vec::new();
                w.layers[l].grid.query_each(x - p2(reach, reach), x + p2(reach, reach), |it| {
                    if it.kind != IT_SEG {
                        return;
                    }
                    let pc = &w.pieces[it.id as usize];
                    let k = it.sub as usize;
                    if !pc.alive || k + 1 >= pc.poly.len() || !w.b.high_speed(pc.net) {
                        return;
                    }
                    let (a, b) = (pc.poly[k], pc.poly[k + 1]);
                    if seg_point_dist(a, b, x) < reach && !segs.iter().any(|s| s.0 == pc.net && s.1 == a && s.2 == b) {
                        segs.push((pc.net, a, b));
                    }
                });
                for (hn, a, b) in segs {
                    let own = own_of(w, hn);
                    if own.contains(&net) || (b - a).norm2() < 1e-12 {
                        continue;
                    }
                    let mut nb = Vec::new();
                    let mut edge = Vec::new();
                    self.merging(&[x], rv, p, &own, NONE, &mut nb, &mut edge);
                    // the via closes a slot across the wire if something it merges with lies beyond the wire
                    let across = |q: P2| segs_intersect(x, q, a, b) || seg_point_dist(a, b, q) < 1e-9;
                    let far: Vec<u32> = nb.iter().cloned().filter(|&o| {
                        let s = &w.objs[o as usize].s;
                        let q = s.pts.iter().cloned().min_by(|u, v| u.dist(x).partial_cmp(&v.dist(x)).unwrap()).unwrap_or(s.center);
                        across(q) || across(s.center)
                    }).collect();
                    let edge_far: Vec<P2> = edge.iter().cloned().filter(|&e| across(e)).collect();
                    if far.is_empty() && edge_far.is_empty() {
                        continue;
                    }
                    let mut members: Vec<u32> = Vec::new();
                    let mut ce: Vec<P2> = edge.clone();
                    for &o in &nb {
                        let c = self.cluster(p, &own, o);
                        for m in c.members {
                            if !members.contains(&m) {
                                members.push(m);
                            }
                        }
                        ce.extend(c.edge);
                    }
                    let foot = a.lerp(b, seg_param(a, b, x));
                    let nrm = (b - a).unit().perp();
                    let dd = self.slot(p, &own, &members, &[(vec![x], rv)], &ce, foot, nrm);
                    let lim = w.b.gap_limit(hn, l);
                    if dd > lim {
                        cost += pen(dd - lim);
                        charged(2);
                    }
                }
            }
        }
        cost
    }

    /// Distance from q to the nearest void in plane p (negative inside one), from through copper, the board
    /// edge, fill keepouts and, on a plane of partial pours, the pours' outlines. `els` are the void objects near.
    fn void_clear(&self, p: usize, q: P2, els: &[(u32, f64)], rings: &[(P2, P2)]) -> f64 {
        let w = self.w;
        let mut d = f64::INFINITY;
        for &(o, r) in els {
            d = d.min(point_hull_dist(&w.objs[o as usize].s.pts, q) - r);
        }
        for &(a, c) in rings {
            d = d.min(seg_point_dist(a, c, q) - w.b.edge_clear);
        }
        d.min(self.area_clear(p, q))
    }

    /// Distance from q to the area voids of plane p (fill keepouts, outside its pours); negative inside one.
    fn area_clear(&self, p: usize, q: P2) -> f64 {
        let b = &self.w.b;
        let mut d = f64::INFINITY;
        for (l, k) in &b.fill_keepouts {
            if *l == p {
                let e = poly_edge_dist(k, q);
                d = d.min(if point_in_poly(k, q) { -e } else { e });
            }
        }
        let regs: Vec<&crate::model::PlaneRegion> = b.plane_regions.iter().filter(|r| r.layer == p).collect();
        if !regs.is_empty() {
            // the pour at q: the highest-priority region holding it
            let inside = |r: &crate::model::PlaneRegion| point_in_poly(&r.poly, q) && !r.holes.iter().any(|h| point_in_poly(h, q));
            match regs.iter().filter(|r| inside(r)).max_by_key(|r| r.priority) {
                None => d = d.min(-regs.iter().map(|r| poly_edge_dist(&r.poly, q)).fold(f64::INFINITY, f64::min)),
                Some(r0) => {
                    let mut e = poly_edge_dist(&r0.poly, q);
                    for h in &r0.holes {
                        e = e.min(poly_edge_dist(h, q));
                    }
                    // another net's pour keeps the clearance off this one: a split
                    let g = b.plane_clear[p].max(b.classes.iter().map(|c| c.clear).fold(0.0, f64::max));
                    for r in regs.iter().filter(|r| r.net != r0.net) {
                        e = e.min(poly_edge_dist(&r.poly, q) - g);
                    }
                    d = d.min(e);
                }
            }
        }
        d
    }

    /// Void objects (with void radii) and board-edge segments near q in plane p.
    fn near(&self, p: usize, q: P2, own: &[u32], reach: f64) -> (Vec<(u32, f64)>, Vec<(P2, P2)>) {
        let w = self.w;
        let mut items = Vec::new();
        w.layers[0].grid.query(q - p2(reach, reach), q + p2(reach, reach), &mut items);
        let mut els = Vec::new();
        let mut rings = Vec::new();
        for it in items {
            match it.kind {
                IT_OBJ => {
                    if let Some(r) = self.void_r(it.id, p, own) {
                        if !els.iter().any(|e: &(u32, f64)| e.0 == it.id) {
                            els.push((it.id, r));
                        }
                    }
                }
                IT_RING => {
                    let ob = &w.objs[it.id as usize];
                    if ob.s.kind == Kind::Outline {
                        let n = ob.s.pts.len();
                        rings.push((ob.s.pts[it.sub as usize], ob.s.pts[(it.sub as usize + 1) % n]));
                    }
                }
                _ => {}
            }
        }
        (els, rings)
    }

    /// Whether plane p is void under q: inside a void, or where no disc of the pour's least width fits around q.
    fn is_void(&self, p: usize, q: P2, own: &[u32]) -> bool {
        let mw = self.mw(p);
        let (els, rings) = self.near(p, q, own, mw + self.grow_max + 1.0);
        let c = self.void_clear(p, q, &els, &rings);
        if c < 0.0 {
            return true;
        }
        if c >= mw / 2.0 {
            return false;
        }
        // copper stays only where a disc of diameter mw fits and covers q
        let h = mw / 2.0;
        for f in [0.5, 1.0] {
            for k in 0..12 {
                let an = k as f64 * std::f64::consts::TAU / 12.0;
                let x = q + p2(an.cos(), an.sin()) * (h * f);
                if self.void_clear(p, x, &els, &rings) >= h - 1e-9 {
                    return false;
                }
            }
        }
        true
    }

    /// The crossing at sample m (track direction u) over a void in plane p: its slot across the track and what
    /// makes the void.
    fn crossing_at(&self, p: usize, m: P2, u: P2, own: &[u32]) -> (f64, Vec<String>) {
        let w = self.w;
        let mw = self.mw(p);
        let nrm = u.perp();
        let (els, _) = self.near(p, m, own, mw + self.grow_max + 1.0);
        let mut members: Vec<u32> = Vec::new();
        let mut edge: Vec<P2> = Vec::new();
        for &(o, r) in &els {
            if point_hull_dist(&w.objs[o as usize].s.pts, m) - r < mw {
                let c = self.cluster(p, own, o);
                for x in c.members {
                    if !members.contains(&x) {
                        members.push(x);
                    }
                }
                edge.extend(c.edge);
            }
        }
        let mut by: Vec<String> = Vec::new();
        let mut d = if members.is_empty() && edge.is_empty() { 0.0 } else { self.slot(p, own, &members, &[], &edge, m, nrm) };
        // area voids: their reach along the normal, stepped
        if self.area_clear(p, m) < mw / 2.0 {
            d = d.max(self.area_slot(p, m, nrm));
            let fk = w.b.fill_keepouts.iter().any(|(l, k)| *l == p && point_in_poly(k, m));
            by.push(if fk {
                "fill keepout".into()
            } else if self.region_at(p, m).is_some() {
                // within the clearance between two nets' pours
                let under = self.region_at(p, m).map(|i| w.b.plane_regions[i].net);
                let mut nets: Vec<String> = w.b.plane_regions.iter().filter(|r| r.layer == p && (Some(r.net) == under || poly_edge_dist(&r.poly, m) < 1.0)).map(|r| w.b.nets[r.net as usize].name.clone()).collect();
                nets.sort();
                nets.dedup();
                format!("split between pours ({})", nets.join(" / "))
            } else {
                "no pour".into()
            });
        }
        let mut names: Vec<String> = Vec::new();
        for &o in &members {
            let s = &w.objs[o as usize].s;
            let net = if s.net != NONE { w.b.nets[s.net as usize].name.as_str() } else { "" };
            let n = match s.kind {
                Kind::Via => format!("via {}", net),
                Kind::Hole => format!("hole {}", s.name),
                // a drawn via kept as copper of its own
                _ if s.name.starts_with("via.s") => format!("via {} (drawn)", net),
                _ => format!("pad {} {}", s.name, net),
            };
            if !names.contains(&n) {
                names.push(n);
            }
        }
        names.sort();
        let more = names.len().saturating_sub(8);
        names.truncate(8);
        if more > 0 {
            names.push(format!("{} more", more));
        }
        by.extend(names);
        if !edge.is_empty() {
            by.push("board edge".into());
        }
        (d, by)
    }

    /// Every crossing of the high-speed tracks `tracks` (net, signal layer, width, points) over a void in their
    /// reference planes, sampled every 0.02 mm along the centre line.
    pub fn report(&self, tracks: &[(u32, usize, f64, Vec<P2>)]) -> Vec<Crossing> {
        let w = self.w;
        let mut out = Vec::new();
        for (net, l, width, pts) in tracks {
            if !w.b.high_speed(*net) || pts.len() < 2 {
                continue;
            }
            let own = own_of(w, *net);
            for &p in &w.b.ref_planes[*l] {
                // samples along the whole line, so a crossing over a bend is one crossing
                let mut samp: Vec<(P2, P2, f64)> = Vec::new();
                for s in pts.windows(2) {
                    let len = s[0].dist(s[1]);
                    if len < 1e-9 {
                        continue;
                    }
                    let k = (len / 0.02).ceil().max(1.0) as usize;
                    let u = (s[1] - s[0]).unit();
                    for i in 0..k {
                        samp.push((s[0].lerp(s[1], (i as f64 + 0.5) / k as f64), u, len / k as f64));
                    }
                }
                let void: Vec<bool> = samp.iter().map(|&(q, _, _)| self.is_void(p, q, &own)).collect();
                let mut i = 0;
                while i < samp.len() {
                    if !void[i] {
                        i += 1;
                        continue;
                    }
                    let j0 = i;
                    let mut over = 0.0;
                    while i < samp.len() && void[i] {
                        over += samp[i].2;
                        i += 1;
                    }
                    let (m, u, _) = samp[(j0 + i - 1) / 2];
                    let (slot, by) = self.crossing_at(p, m, u, &own);
                    out.push(Crossing { net: *net, layer: *l, plane: p, at: m, over, slot, limit: w.b.gap_limit(*net, *l), width: *width, by });
                }
            }
        }
        out
    }
}
