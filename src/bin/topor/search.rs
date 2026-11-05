//! Topological A*: states are slots between existing wires on triangulation edges, so a new wire may pass
//! between any wires (pushing them aside) but never crosses another net's wire on its layer. Layer changes
//! place a via at a free point of the current face. A soft mode crosses wires at a price and reports them
//! as victims for rip-up.

use crate::cdt::NONE;
use crate::geom::*;
use crate::model::Kind;
use crate::world::{Attr, Term, World, T_BLOCK};
use std::cmp::Ordering;
use crate::fx::FxMap;
use std::collections::BinaryHeap;
type HashMap<K, V> = FxMap<K, V>;

#[derive(Clone, Copy, PartialEq, Debug)]
enum NK {
    Gap,
    Vert,
    Via,
    GoalHull,
    GoalVert,
    GoalVia,
}

#[derive(Clone, Copy, Debug)]
struct Node {
    layer: u8,
    kind: NK,
    tri: u32,
    k: u8,
    gap: u16,
    vtx: u32,
    pos: P2,
    g: f64,
    parent: u32,
    edge: u32,
    ins: u32,
    left: u32,
    obj: u32,
    vic: u32,
    nvic: u16,
}

#[derive(Clone, Copy)]
pub enum End {
    Obj(Term),
    NewVia(usize),
}

#[derive(Clone)]
pub struct Seg {
    pub layer: usize,
    pub a: End,
    pub b: End,
    pub cx: Vec<(u32, usize, u32)>,
}

#[derive(Clone)]
pub struct PathOut {
    pub segs: Vec<Seg>,
    pub vias: Vec<P2>,
    pub cost: f64,
    pub victims: Vec<u32>,
    pub target: u32,
    /// Triangles where the path crossed (would rip) other wires.
    pub cross_tris: Vec<(u8, u32)>,
}

pub struct Query {
    pub net: u32,
    pub width: Vec<f64>,
    pub allowed: u32,
    pub sources: Vec<u32>,
    pub targets: Vec<u32>,
    pub plane_goal: bool,
    pub via_cost: f64,
    pub soft: bool,
    pub cross_cost: Vec<f64>,
    pub max_exp: usize,
    pub layer_cost: Vec<f64>,
    pub max_len: f64,
    /// Pieces that may be crossed for free (they are about to be ripped).
    pub ignore: Vec<u32>,
    /// Extra cost for entering triangles where an earlier attempt failed to realise.
    pub taboo: HashMap<(u8, u32), f64>,
    /// Negotiated congestion history per triangle (shared by all nets).
    pub history: std::sync::Arc<HashMap<(u8, u32), f64>>,
    /// Via positions to avoid (centre, radius).
    pub via_taboo: Vec<(P2, f64)>,
    /// Only check corners against fixed opposite edges.
    pub wedge_off: bool,
    /// Differential pair partner to run alongside (NONE: none), and the cost factor for coupled steps.
    pub couple: u32,
    pub couple_factor: f64,
    /// Width to reserve (per layer) away from the terminals, so a partner fits alongside; 0 = none.
    pub ribbon: Vec<f64>,
    /// Charge slots in high-speed nets' reference planes (gaps.rs), where the board has such nets.
    pub gap_rule: bool,
    /// Dogbone nets (`cold_routing.dogbone`): the only places a via may go, right beside the net's own pads,
    /// and the layers of those pads (bit per signal layer): each of its vias joins a pad layer to one other
    /// layer, so there is one at each end and no hop between inner layers.
    pub via_spots: Option<std::sync::Arc<Vec<P2>>>,
    pub via_spot_layers: u32,
}

#[derive(Clone, Copy, PartialEq)]
struct HeapE(f64, u32);

/// Weight of the distance heuristic (weighted A*), 1 (plain A*) by default: TOPOR_HARD_W / TOPOR_SOFT_W set
/// one. A greedy soft search (1.5) was faster but ripped other victims, and lost connections on 1BP and EXT.
fn heuristic_weight(soft: bool) -> f64 {
    static W: std::sync::OnceLock<(f64, f64)> = std::sync::OnceLock::new();
    let w = W.get_or_init(|| {
        let env = |k: &str, d: f64| std::env::var(k).ok().and_then(|x| x.parse().ok()).unwrap_or(d);
        (env("TOPOR_HARD_W", 1.0), env("TOPOR_SOFT_W", 1.0))
    });
    if soft {
        w.1
    } else {
        w.0
    }
}
impl Eq for HeapE {}
impl PartialOrd for HeapE {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for HeapE {
    fn cmp(&self, o: &Self) -> Ordering {
        o.0.partial_cmp(&self.0).unwrap_or(Ordering::Equal).then(o.1.cmp(&self.1))
    }
}

#[derive(Clone, Copy, Debug)]
enum CE {
    V(u32),
    C(u32),
}

struct Cyc {
    els: Vec<CE>,
    start: [usize; 4],
    /// (pos a, pos b, net, piece) with a < b
    chords: Vec<(usize, usize, u32, u32)>,
}

pub struct Search<'a> {
    trace: bool,
    w: &'a World,
    q: &'a Query,
    nodes: Vec<Node>,
    heap: BinaryHeap<HeapE>,
    best: HashMap<(u8, u8, u32, u32, u32), f64>,
    cyc: HashMap<(u8, u32), std::rc::Rc<Cyc>>,
    vfit: HashMap<(i64, i64), bool>,
    victims: Vec<u32>,
    tverts: Vec<HashMap<u32, u32>>,
    thull: Vec<HashMap<u32, u32>>,
    tpts: Vec<P2>,
    tlayers: u32,
    spts: Vec<P2>,
    partner_vias: Vec<P2>,
    pub expansions: usize,
    /// Plane gap rule (boards with high-speed nets): costs of slots crossed or closed; whether this net keeps
    /// off them, and the nets whose copper is its own.
    gaps: Option<crate::gaps::Gaps<'a>>,
    gap_hs: bool,
    gap_own: Vec<u32>,
}

fn nx(i: usize) -> usize {
    (i + 1) % 3
}

/// A signal via joins wires on two layers; a wire may end at it on layer l only if it already carries
/// a wire on another layer (otherwise it would be a one-layer, dangling via). Plane vias reach the plane.
fn via_terminal_ok(w: &World, o: u32, l: usize) -> bool {
    let ob = &w.objs[o as usize];
    if ob.s.kind != Kind::Via || ob.plane {
        return true;
    }
    let live: Vec<u32> = w.obj_pieces[o as usize].iter().cloned().filter(|&p| w.pieces[p as usize].alive).collect();
    live.is_empty() || live.iter().any(|&p| w.pieces[p as usize].layer as usize != l)
}

impl<'a> Search<'a> {
    pub fn new(w: &'a World, q: &'a Query) -> Search<'a> {
        let nl = w.nsig();
        let mut tverts: Vec<HashMap<u32, u32>> = vec![HashMap::default(); nl];
        let mut thull: Vec<HashMap<u32, u32>> = vec![HashMap::default(); nl];
        let mut tpts = Vec::new();
        let mut tlayers = 0u32;
        for &o in &q.targets {
            let ob = &w.objs[o as usize];
            tpts.push(ob.s.center);
            for l in 0..nl {
                if !ob.on(l) || q.allowed & (1 << l) == 0 {
                    continue;
                }
                tlayers |= 1 << l;
                if ob.is_poly() {
                    let n = ob.npts();
                    for k in 0..n {
                        if let Some(e) = w.layers[l].cdt.find_edge(ob.v(l, k), ob.v(l, (k + 1) % n)) {
                            thull[l].insert(e, o);
                        }
                    }
                } else {
                    if !via_terminal_ok(w, o, l) {
                        continue;
                    }
                    for k in 0..ob.npts() {
                        let v = ob.v(l, k);
                        if v != NONE {
                            tverts[l].insert(v, o);
                        }
                    }
                }
            }
        }
        let spts: Vec<P2> = q.sources.iter().map(|&o| w.objs[o as usize].s.center).collect();
        let mut partner_vias = Vec::new();
        if q.couple != NONE {
            for o in &w.objs {
                if o.alive && o.s.kind == Kind::Via && o.s.net == q.couple {
                    partner_vias.push(o.s.center);
                }
            }
        }
        Search {
            w,
            q,
            nodes: Vec::new(),
            heap: BinaryHeap::new(),
            best: HashMap::default(),
            cyc: HashMap::default(),
            vfit: HashMap::default(),
            victims: Vec::new(),
            tverts,
            thull,
            tpts,
            tlayers,
            spts,
            partner_vias,
            trace: std::env::var("TOPOR_SEARCHTRACE").is_ok(),
            expansions: 0,
            gaps: if q.gap_rule && crate::gaps::active(w) { Some(crate::gaps::Gaps::new(w)) } else { None },
            gap_hs: w.b.high_speed(q.net),
            gap_own: crate::gaps::own_of(w, q.net),
        }
    }

    fn h(&self, l: usize, p: P2) -> f64 {
        let k = if self.q.couple != NONE { self.q.couple_factor } else { 1.0 };
        k * self.h_raw(l, p) * heuristic_weight(self.q.soft)
    }

    fn h_raw(&self, l: usize, p: P2) -> f64 {
        let mut d = f64::INFINITY;
        for &t in &self.tpts {
            d = d.min(t.dist(p));
        }
        if self.q.plane_goal {
            d = 0.0;
        }
        if self.tlayers & (1 << l) == 0 && !self.q.plane_goal {
            d += self.q.via_cost;
        }
        if d.is_infinite() {
            0.0
        } else {
            d
        }
    }

    fn goal_pt(&self, p: P2) -> P2 {
        let mut best = p;
        let mut d = f64::INFINITY;
        for &t in &self.tpts {
            let x = t.dist(p);
            if x < d {
                d = x;
                best = t;
            }
        }
        best
    }

    fn push(&mut self, n: Node, key: (u8, u8, u32, u32, u32)) {
        if let Some(&g) = self.best.get(&key) {
            if g <= n.g + 1e-9 {
                return;
            }
        }
        self.best.insert(key, n.g);
        let f = n.g + match n.kind {
            NK::GoalHull | NK::GoalVert | NK::GoalVia => 0.0,
            _ => self.h(n.layer as usize, n.pos),
        };
        self.nodes.push(n);
        self.heap.push(HeapE(f, (self.nodes.len() - 1) as u32));
    }

    fn cycle(&mut self, l: usize, t: u32) -> std::rc::Rc<Cyc> {
        if let Some(c) = self.cyc.get(&(l as u8, t)) {
            return c.clone();
        }
        let w = self.w;
        let layer = &w.layers[l];
        let tr = layer.cdt.tris[t as usize];
        let mut els = Vec::new();
        let mut start = [0usize; 4];
        for k in 0..3 {
            start[k] = els.len();
            els.push(CE::V(tr.v[k]));
            let e = tr.e[k];
            let ed = layer.cdt.edges[e as usize];
            let list = &layer.cross[e as usize];
            if tr.v[k] == ed.a {
                for &c in list {
                    els.push(CE::C(c));
                }
            } else {
                for &c in list.iter().rev() {
                    els.push(CE::C(c));
                }
            }
        }
        start[3] = els.len();
        let mut chords = Vec::new();
        for k in 0..3 {
            for i in start[k] + 1..start[k + 1] {
                let c = match els[i] {
                    CE::C(c) => c,
                    _ => continue,
                };
                let cr = &w.crosses[c as usize];
                let pc = &w.pieces[cr.piece as usize];
                let entering = cr.left == tr.v[k];
                let ci = cr.idx as usize;
                let partner: Option<usize> = if entering {
                    if ci + 1 < pc.cx.len() {
                        let pcid = pc.cx[ci + 1];
                        (0..els.len()).find(|&j| matches!(els[j], CE::C(x) if x == pcid))
                    } else if pc.b.vtx != NONE {
                        (0..3).find(|&j| tr.v[j] == pc.b.vtx).map(|j| start[j])
                    } else {
                        None
                    }
                } else if ci > 0 {
                    let pcid = pc.cx[ci - 1];
                    (0..els.len()).find(|&j| matches!(els[j], CE::C(x) if x == pcid))
                } else if pc.a.vtx != NONE {
                    (0..3).find(|&j| tr.v[j] == pc.a.vtx).map(|j| start[j])
                } else {
                    None
                };
                if let Some(j) = partner {
                    let is_v = matches!(els[j], CE::V(_));
                    if is_v || i < j {
                        chords.push((i.min(j), i.max(j), pc.net, cr.piece));
                    }
                }
            }
        }
        let c = std::rc::Rc::new(Cyc { els, start, chords });
        self.cyc.insert((l as u8, t), c.clone());
        c
    }

    fn blocking(&self, net: u32, piece: u32) -> bool {
        net != self.q.net && !self.q.ignore.contains(&piece)
    }

    /// Geometric side of an interior point x relative to chord (a, b) of triangle t: true = inside (a, b).
    fn side(&self, l: usize, t: u32, cy: &Cyc, ch: (usize, usize, u32, u32), x: P2) -> bool {
        let w = self.w;
        let tr = w.layers[l].cdt.tris[t as usize];
        let (a, b, _, p) = ch;
        let endpoint = |i: usize| -> (P2, Option<f64>) {
            match cy.els[i] {
                CE::V(v) => (w.vpos(l, v), None),
                CE::C(c) => {
                    let cr = &w.crosses[c as usize];
                    let par = w.pieces[cr.piece as usize].cparam[cr.idx as usize];
                    (w.point_at(cr.piece, par), Some(par))
                }
            }
        };
        let (pa, ta) = endpoint(a);
        let (pb, tb) = endpoint(b);
        let pc = &w.pieces[p as usize];
        let pa_par = ta.unwrap_or(if pa.dist(w.term_point(l, pc.a)) < 1e-9 { 0.0 } else { (pc.poly.len() - 1) as f64 });
        let pb_par = tb.unwrap_or(if pb.dist(w.term_point(l, pc.a)) < 1e-9 { 0.0 } else { (pc.poly.len() - 1) as f64 });
        let mut poly = vec![pa];
        for k in 0..3 {
            let s = cy.start[k];
            if s > a && s < b {
                poly.push(w.vpos(l, tr.v[k]));
            }
        }
        poly.push(pb);
        // chord from b back to a
        let (lo, hi) = if pa_par < pb_par { (pa_par, pb_par) } else { (pb_par, pa_par) };
        let mut mid = w.sub_poly(p, lo, hi);
        if pa_par < pb_par {
            mid.reverse();
        }
        if mid.len() > 2 {
            poly.extend_from_slice(&mid[1..mid.len() - 1]);
        }
        point_in_poly(&poly, x)
    }

    fn edge_pos(&self, l: usize, e: u32, lo: f64, hi: f64, from: P2) -> P2 {
        let w = self.w;
        let ed = w.layers[l].cdt.edges[e as usize];
        let (a, b) = (w.vpos(l, ed.a), w.vpos(l, ed.b));
        let len = a.dist(b);
        let goal = self.goal_pt(from);
        let d = b - a;
        let r = goal - from;
        let den = r.cross(d);
        let mut t = if den.abs() > 1e-12 { (a - from).cross(r) / den * len } else { ((from - a).dot(d) / len).clamp(0.0, len) };
        if self.q.plane_goal {
            t = (from - a).dot(d) / len;
        }
        let t = t.clamp(lo, hi.max(lo));
        a + d * (t / len)
    }

    fn corner_ok(&self, l: usize, t: u32, cy: &Cyc, k_in: usize, gap_in: usize, k_out: usize, width: f64) -> Result<(), Vec<u32>> {
        let w = self.w;
        let tr = w.layers[l].cdt.tris[t as usize];
        let (shared, opp) = if k_out == nx(k_in) { (nx(k_in), (k_in + 2) % 3) } else { (k_in, nx(k_in)) };
        let v = tr.v[shared];
        let pv = w.vpos(l, v);
        let vobj = w.layers[l].vobj[v as usize];
        let me = w.net_attr(self.q.net, l);
        // wires around the shared vertex: chords between edges k_in and k_out, ordered from the vertex
        let in_lo = cy.start[k_in];
        let in_hi = cy.start[k_in + 1];
        let out_lo = cy.start[k_out];
        let out_hi = cy.start[k_out + 1];
        let mut around: Vec<(usize, Attr, f64)> = Vec::new();
        let mut around_pieces: Vec<u32> = Vec::new();
        for &(a, b, _, p) in &cy.chords {
            let on_in = |x: usize| x > in_lo && x < in_hi;
            let on_out = |x: usize| x > out_lo && x < out_hi;
            let pos_in = if on_in(a) && on_out(b) {
                Some(a)
            } else if on_in(b) && on_out(a) {
                Some(b)
            } else {
                None
            };
            if let Some(pi) = pos_in {
                let pc = &w.pieces[p as usize];
                around.push((pi, w.net_attr(pc.net, pc.layer as usize), pc.width));
                around_pieces.push(p);
            }
        }
        let gap_pos = in_lo + gap_in;
        for x in around.iter_mut() {
            x.0 *= 2;
        }
        around.push((gap_pos * 2 + 1, me, width));
        around.sort_by_key(|x| x.0);
        if shared != k_in {
            around.reverse();
        }
        let ws: Vec<(Attr, f64)> = around.iter().map(|x| (x.1, x.2)).collect();
        let c = w.pack(l, v, &ws);
        let mut reach: f64 = 0.0;
        for (i, (_, wd)) in ws.iter().enumerate() {
            reach = reach.max(c[i] + wd / 2.0);
        }
        let oe = tr.e[opp];
        let fixed = w.layers[l].cdt.edges[oe as usize].fixed;
        let (pa, pb) = (w.vpos(l, tr.v[opp]), w.vpos(l, tr.v[nx(opp)]));
        if !fixed && (self.q.wedge_off || reach + w.max_clear < seg_point_dist(pa, pb, pv) - 1e-9) {
            return Ok(());
        }
        // the stack reaches beyond the triangle: everything else inside the wedge at v must clear it
        let d1 = w.vpos(l, tr.v[nx(shared)]) - pv;
        let d2 = w.vpos(l, tr.v[(shared + 2) % 3]) - pv;
        let in_wedge = |q: P2| d1.cross(q - pv) >= -1e-12 && (q - pv).cross(d2) >= -1e-12;
        let far = reach + w.max_clear + 0.05;
        let ray = |dir: P2| (pv, pv + dir.unit() * (far + 1.0));
        let (r1, r2) = (ray(d1), ray(d2));
        // nearest distance from v to points of segment a-b inside the wedge
        let seg_near = |a: P2, b: P2| -> f64 {
            let mut best = f64::INFINITY;
            let f = a.lerp(b, seg_param(a, b, pv));
            for q in [a, b, f] {
                if in_wedge(q) {
                    best = best.min(q.dist(pv));
                }
            }
            for (ra, rb) in [r1, r2] {
                if let Some((t, _)) = seg_cross_params(ra, rb, a, b) {
                    best = best.min(ra.lerp(rb, t).dist(pv));
                }
            }
            best
        };
        let mut items = Vec::new();
        w.layers[l].grid.query(pv - p2(far, far), pv + p2(far, far), &mut items);
        items.sort_unstable_by_key(|it| (it.kind, it.id, it.sub));
        items.dedup();
        let mut victims: Vec<u32> = Vec::new();
        let mut hard = false;
        for it in &items {
            let (dist, attr, piece_victims): (f64, Attr, Option<Vec<u32>>) = match it.kind {
                crate::world::IT_SEG => {
                    if around_pieces.contains(&it.id) {
                        continue;
                    }
                    let q = &w.pieces[it.id as usize];
                    if !q.alive || q.net == self.q.net {
                        continue;
                    }
                    let k = it.sub as usize;
                    if k + 1 >= q.poly.len() {
                        continue;
                    }
                    (seg_near(q.poly[k], q.poly[k + 1]) - q.width / 2.0, w.net_attr(q.net, q.layer as usize), Some(vec![it.id]))
                }
                _ => {
                    if it.id == vobj {
                        continue;
                    }
                    let ob = &w.objs[it.id as usize];
                    if !ob.alive || (ob.s.net != NONE && ob.s.net == self.q.net) {
                        continue;
                    }
                    let n = ob.s.pts.len();
                    let mut d = f64::INFINITY;
                    if it.kind == crate::world::IT_RING {
                        d = seg_near(ob.s.pts[it.sub as usize], ob.s.pts[(it.sub as usize + 1) % n]);
                    } else if n == 1 {
                        if in_wedge(ob.s.pts[0]) {
                            d = ob.s.pts[0].dist(pv);
                        } else {
                            // outside the wedge: only where the wires cross its sides (at their packed
                            // radius) can come near it; the vertex itself is not wire
                            let q0 = ob.s.pts[0];
                            let oa = w.obj_attr(it.id);
                            let clear_of = ws.iter().enumerate().all(|(i, (a, wd))| match w.clr(a, &oa) {
                                Some(g) => {
                                    let (e1, e2) = (pv + d1.unit() * c[i], pv + d2.unit() * c[i]);
                                    q0.dist(e1).min(q0.dist(e2)) - ob.s.r >= wd / 2.0 + g - 1e-9
                                }
                                None => true,
                            });
                            if clear_of {
                                continue;
                            }
                            d = 0.0;
                        }
                    } else {
                        let m = if n == 2 { 1 } else { n };
                        for k in 0..m {
                            d = d.min(seg_near(ob.s.pts[k], ob.s.pts[(k + 1) % n]));
                        }
                    }
                    let r = if matches!(ob.s.kind, crate::model::Kind::Outline | crate::model::Kind::Keepout) { 0.0 } else { ob.s.r };
                    let pv_ = if ob.s.kind == crate::model::Kind::Via && ob.s.net != self.q.net {
                        let ps: Vec<u32> = w.obj_pieces[it.id as usize].iter().cloned().filter(|&x| w.pieces[x as usize].alive).collect();
                        if ps.is_empty() {
                            None
                        } else {
                            Some(ps)
                        }
                    } else {
                        None
                    };
                    (d - r, w.obj_attr(it.id), pv_)
                }
            };
            if !dist.is_finite() {
                continue;
            }
            let mut ok = true;
            for (i, (a, wd)) in ws.iter().enumerate() {
                if let Some(g) = w.clr(a, &attr) {
                    if c[i] + wd / 2.0 + g > dist + 1e-9 {
                        ok = false;
                        break;
                    }
                }
            }
            if !ok {
                match piece_victims {
                    Some(v) => {
                        for p in v {
                            if !victims.contains(&p) {
                                victims.push(p);
                            }
                        }
                    }
                    None => hard = true,
                }
            }
        }
        if victims.is_empty() && !hard {
            return Ok(());
        }
        if hard {
            return Err(vec![]);
        }
        for p in around_pieces {
            if w.pieces[p as usize].net != self.q.net && !victims.contains(&p) {
                victims.push(p);
            }
        }
        Err(victims)
    }

    pub fn run(&mut self) -> Option<PathOut> {
        let w = self.w;
        let q = self.q;
        // sources
        for &o in &q.sources {
            let ob = &w.objs[o as usize];
            for l in 0..w.nsig() {
                if !ob.on(l) || q.allowed & (1 << l) == 0 {
                    continue;
                }
                let layer = &w.layers[l];
                if ob.is_poly() {
                    let n = ob.npts();
                    for k in 0..n {
                        let e = match layer.cdt.find_edge(ob.v(l, k), ob.v(l, (k + 1) % n)) {
                            Some(e) => e,
                            None => continue,
                        };
                        for side in 0..2 {
                            let t = layer.etri[e as usize][side];
                            if t == NONE || layer.tflag[t as usize] & T_BLOCK != 0 {
                                continue;
                            }
                            let tr = layer.cdt.tris[t as usize];
                            let ki = tr.e.iter().position(|&x| x == e).unwrap();
                            let ed = layer.cdt.edges[e as usize];
                            let nc = layer.cross[e as usize].len();
                            let gap = nc;
                            let ins = if tr.v[ki] == ed.a { gap } else { nc - gap };
                            let (a, b) = (w.vpos(l, tr.v[ki]), w.vpos(l, tr.v[nx(ki)]));
                            let pos = a.lerp(b, 0.5);
                            let n = Node {
                                layer: l as u8,
                                kind: NK::Gap,
                                tri: t,
                                k: ki as u8,
                                gap: gap as u16,
                                vtx: NONE,
                                pos,
                                g: ob.s.center.dist(pos) * q.layer_cost[l],
                                parent: NONE,
                                edge: e,
                                ins: ins as u32,
                                left: tr.v[ki],
                                obj: o,
                                vic: 0,
                                nvic: 0,
                            };
                            self.push(n, (l as u8, 0, e, gap as u32, side as u32));
                        }
                    }
                } else {
                    if !via_terminal_ok(w, o, l) {
                        continue;
                    }
                    for kk in 0..ob.npts() {
                        let v = ob.v(l, kk);
                        if v == NONE {
                            continue;
                        }
                        for t in layer.cdt.star(v) {
                            if layer.tflag[t as usize] & T_BLOCK != 0 {
                                continue;
                            }
                            let n = Node {
                                layer: l as u8,
                                kind: NK::Vert,
                                tri: t,
                                k: 0,
                                gap: 0,
                                vtx: v,
                                pos: w.vpos(l, v),
                                g: 0.0,
                                parent: NONE,
                                edge: NONE,
                                ins: 0,
                                left: NONE,
                                obj: o,
                                vic: 0,
                                nvic: 0,
                            };
                            self.push(n, (l as u8, 1, v, t, 0));
                        }
                    }
                }
            }
        }
        while let Some(HeapE(_, ni)) = self.heap.pop() {
            let n = self.nodes[ni as usize];
            if matches!(n.kind, NK::GoalHull | NK::GoalVert | NK::GoalVia) {
                return Some(self.reconstruct(ni));
            }
            let key = self.key_of(&n);
            if let Some(&g) = self.best.get(&key) {
                if g < n.g - 1e-9 {
                    continue;
                }
            }
            self.expansions += 1;
            if self.expansions > q.max_exp {
                return None;
            }
            if self.trace {
                eprintln!("   expand {:?} layer {} pos ({:.3},{:.3}) g {:.3} tri {}", n.kind, n.layer, n.pos.x, n.pos.y, n.g, n.tri);
            }
            if n.g > q.max_len {
                continue;
            }
            self.expand(ni);
        }
        None
    }

    fn key_of(&self, n: &Node) -> (u8, u8, u32, u32, u32) {
        match n.kind {
            NK::Gap => {
                let ed = self.w.layers[n.layer as usize].cdt.edges[n.edge as usize];
                let side = (n.left == ed.a) as u32;
                (n.layer, 0, n.edge, n.ins, side)
            }
            NK::Vert => (n.layer, 1, n.vtx, n.tri, 0),
            NK::Via => (n.layer, 2, (n.pos.x * 1e4) as i64 as u32, (n.pos.y * 1e4) as i64 as u32, 0),
            _ => (n.layer, 9, n.edge, n.vtx, n.obj),
        }
    }

    fn expand(&mut self, ni: u32) {
        let w = self.w;
        let q = self.q;
        let n = self.nodes[ni as usize];
        let l = n.layer as usize;
        let layer = &w.layers[l];
        let t = n.tri;
        let tr = layer.cdt.tris[t as usize];
        let cy = self.cycle(l, t);
        let mut width = q.width[l];
        if q.ribbon.get(l).cloned().unwrap_or(0.0) > width {
            let near = self.spts.iter().chain(self.tpts.iter()).any(|&c| c.dist(n.pos) < 1.2);
            if !near {
                width = q.ribbon[l];
            }
        }
        // entry position in doubled cycle coordinates, or an interior point
        let entry_x: Option<usize> = match n.kind {
            NK::Gap => Some(2 * (cy.start[n.k as usize] + n.gap as usize) + 1),
            NK::Vert => {
                let j = (0..3).find(|&j| tr.v[j] == n.vtx).unwrap();
                Some(2 * cy.start[j])
            }
            _ => None,
        };
        let blocking: Vec<(usize, usize, u32, u32)> = cy.chords.iter().cloned().filter(|c| self.blocking(c.2, c.3)).collect();
        let entry_side: Vec<bool> = match entry_x {
            Some(x) => blocking.iter().map(|c| 2 * c.0 < x && x < 2 * c.1).collect(),
            None => blocking.iter().map(|&c| self.side(l, t, &cy, c, n.pos)).collect(),
        };
        let sep_cost = |y: usize, this: &Self| -> (f64, Vec<u32>) {
            let mut cost = 0.0;
            let mut v = Vec::new();
            for (i, c) in blocking.iter().enumerate() {
                let ins = 2 * c.0 < y && y < 2 * c.1;
                let on_end = y == 2 * c.0 || y == 2 * c.1;
                if ins != entry_side[i] && !on_end {
                    if !q.soft {
                        return (f64::INFINITY, v);
                    }
                    let pn = this.w.pieces[c.3 as usize].net as usize;
                    cost += q.cross_cost.get(pn).cloned().unwrap_or(10.0);
                    v.push(c.3);
                }
            }
            (cost, v)
        };
        // edge successors
        for k2 in 0..3 {
            if n.kind == NK::Gap && k2 == n.k as usize {
                continue;
            }
            if n.kind == NK::Vert && (tr.v[k2] == n.vtx || tr.v[nx(k2)] == n.vtx) {
                continue;
            }
            let e2 = tr.e[k2];
            let ed = layer.cdt.edges[e2 as usize];
            let nc = layer.cross[e2 as usize].len();
            let tgt_hull = self.thull[l].get(&e2).cloned();
            if ed.fixed && tgt_hull.is_none() {
                continue;
            }
            let al = layer.along[e2 as usize];
            let mut along_victim = None;
            if al != NONE && w.pieces[al as usize].net != q.net {
                if !q.soft {
                    continue;
                }
                along_victim = Some(al);
            }
            let nb = layer.etri[e2 as usize][if tr.v[k2] == ed.a { 1 } else { 0 }];
            if tgt_hull.is_none() && (nb == NONE || layer.tflag[nb as usize] & T_BLOCK != 0) {
                continue;
            }
            // reachable gaps on this edge
            let mut best: Option<(f64, usize, Vec<u32>)> = None;
            let mut run: Vec<usize> = Vec::new();
            let net_of_el = |idx: usize| -> u32 {
                match cy.els[idx] {
                    CE::C(c) => w.pieces[w.crosses[c as usize].piece as usize].net,
                    CE::V(_) => NONE,
                }
            };
            for g in 0..=nc {
                // the slot between the two members of a differential pair belongs to the pair
                if g > 0 && g < nc {
                    let (ia, ib) = (cy.start[k2] + g, cy.start[k2] + g + 1);
                    let (na, nb) = (net_of_el(ia), net_of_el(ib));
                    if na != NONE && na != q.net && nb != q.net && w.b.partner.get(na as usize) == Some(&nb) {
                        // only where the two members actually run coupled across this edge
                        let pt = |idx: usize| -> Option<(P2, f64)> {
                            match cy.els[idx] {
                                CE::C(c) => {
                                    let cr = &w.crosses[c as usize];
                                    let pc = &w.pieces[cr.piece as usize];
                                    pc.cparam.get(cr.idx as usize).map(|&t| (w.point_at(cr.piece, t), pc.width))
                                }
                                _ => None,
                            }
                        };
                        if let (Some((pa, wa)), Some((pb, wb))) = (pt(ia), pt(ib)) {
                            if pa.dist(pb) < (wa + wb) / 2.0 + w.b.pair_gap[na as usize] + 0.08 {
                                continue;
                            }
                        }
                    }
                }
                let y = 2 * (cy.start[k2] + g) + 1;
                let (c, v) = sep_cost(y, self);
                if c == 0.0 {
                    run.push(g);
                } else if c.is_finite() && best.as_ref().map_or(true, |b| c < b.0) {
                    best = Some((c, g, v));
                }
            }
            if run.len() > 1 && run.windows(2).any(|w2| w2[1] != w2[0] + 1) {
                let mut parts: Vec<Vec<usize>> = vec![vec![run[0]]];
                for &g in &run[1..] {
                    if g == parts.last().unwrap().last().unwrap() + 1 {
                        parts.last_mut().unwrap().push(g);
                    } else {
                        parts.push(vec![g]);
                    }
                }
                run = parts.into_iter().max_by_key(|x| x.len()).unwrap();
            }
            let mut cands: Vec<(f64, usize, Vec<u32>)> = Vec::new();
            let mut coupled_gap = NONE as usize;
            if !run.is_empty() {
                let mut pick = run[run.len() / 2];
                if q.couple != NONE {
                    let net_at = |idx: usize| -> u32 {
                        match cy.els[idx] {
                            CE::C(c) => w.pieces[w.crosses[c as usize].piece as usize].net,
                            CE::V(_) => NONE,
                        }
                    };
                    let (g_lo, g_hi) = (run[0], run[run.len() - 1]);
                    if g_lo > 0 && net_at(cy.start[k2] + g_lo) == q.couple {
                        pick = g_lo;
                        coupled_gap = pick;
                    } else if g_hi < nc && net_at(cy.start[k2] + g_hi + 1) == q.couple {
                        pick = g_hi;
                        coupled_gap = pick;
                    }
                }
                cands.push((0.0, pick, vec![]));
            }
            if let Some(av) = along_victim {
                for c in cands.iter_mut() {
                    c.0 += q.cross_cost.get(w.pieces[av as usize].net as usize).cloned().unwrap_or(10.0);
                    c.2.push(av);
                }
            }
            if let Some(b) = best {
                if q.soft {
                    cands.push(b);
                }
            }
            for (ccost, g, vics) in cands {
                if let Some(obj) = tgt_hull {
                    // arrive at a target pad through its hull edge
                    let (a, b) = (w.vpos(l, tr.v[k2]), w.vpos(l, tr.v[nx(k2)]));
                    let c = w.objs[obj as usize].s.center;
                    let p = seg_point_dist(a, b, n.pos);
                    let _ = p;
                    let pos = a.lerp(b, seg_param(a, b, n.pos));
                    let ins = if tr.v[k2] == ed.a { g } else { nc - g };
                    let gg = n.g + (n.pos.dist(pos) + pos.dist(c)) * q.layer_cost[l] + ccost;
                    let vic = self.add_victims(&vics);
                    let node = Node { layer: l as u8, kind: NK::GoalHull, tri: t, k: k2 as u8, gap: g as u16, vtx: NONE, pos: c, g: gg, parent: ni, edge: e2, ins: ins as u32, left: tr.v[nx(k2)], obj, vic, nvic: vics.len() as u16 };
                    let key = self.key_of(&node);
                    self.push(node, key);
                    continue;
                }
                let ins = if tr.v[k2] == ed.a { g } else { nc - g };
                let fit = w.edge_fit(l, e2, ins, q.net, width);
                if self.trace && fit.is_none() {
                    let ed2 = layer.cdt.edges[e2 as usize];
                    eprintln!("      edge {} ({:.3},{:.3})-({:.3},{:.3}) len {:.3} does not fit", e2, w.vpos(l, ed2.a).x, w.vpos(l, ed2.a).y, w.vpos(l, ed2.b).x, w.vpos(l, ed2.b).y, w.edge_len(l, e2));
                }
                let (lo, hi, extra, vics) = match fit {
                    Some((lo, hi)) => (lo, hi, 0.0, vics),
                    None => {
                        if !q.soft {
                            continue;
                        }
                        // overflow: the neighbours at this slot become victims
                        let list = &layer.cross[e2 as usize];
                        let mut v = vics.clone();
                        let mut extra = 0.0;
                        for &j in [ins.wrapping_sub(1), ins].iter() {
                            if j < list.len() {
                                let pc = w.crosses[list[j] as usize].piece;
                                let pn = w.pieces[pc as usize].net;
                                if pn != q.net {
                                    v.push(pc);
                                    extra += q.cross_cost.get(pn as usize).cloned().unwrap_or(10.0);
                                }
                            }
                        }
                        // other nets' vias at the edge ends can be moved: rip the pieces holding them
                        for vx in [ed.a, ed.b] {
                            let o = layer.vobj[vx as usize];
                            if o == NONE {
                                continue;
                            }
                            let ob = &w.objs[o as usize];
                            if ob.alive && ob.s.kind == crate::model::Kind::Via && ob.s.net != q.net {
                                for &pp in &w.obj_pieces[o as usize] {
                                    if w.pieces[pp as usize].alive && !v.contains(&pp) {
                                        v.push(pp);
                                        extra += q.cross_cost.get(ob.s.net as usize).cloned().unwrap_or(10.0);
                                    }
                                }
                            }
                        }
                        if v.len() == vics.len() {
                            continue;
                        }
                        let len = w.edge_len(l, e2);
                        (len * 0.3, len * 0.7, extra, v)
                    }
                };
                let mut vics = vics;
                let mut extra = extra;
                if n.kind == NK::Gap {
                    if let Err(v) = self.corner_ok(l, t, &cy, n.k as usize, n.gap as usize, k2, width) {
                        if self.trace {
                            eprintln!("      corner rejects edge {} ({} victims: {:?})", e2, v.len(), v.iter().map(|&p| format!("{} {:?}", w.b.nets[w.pieces[p as usize].net as usize].name, w.pieces[p as usize].poly.first().zip(w.pieces[p as usize].poly.last()).map(|(a, b)| ((a.x * 100.0).round() / 100.0, (a.y * 100.0).round() / 100.0, (b.x * 100.0).round() / 100.0, (b.y * 100.0).round() / 100.0)))).collect::<Vec<_>>());
                        }
                        if !q.soft || v.is_empty() {
                            continue;
                        }
                        for p in v {
                            extra += q.cross_cost.get(w.pieces[p as usize].net as usize).cloned().unwrap_or(10.0);
                            vics.push(p);
                        }
                    }
                }
                let pos = self.edge_pos(l, e2, lo, hi, n.pos);
                if matches!(n.kind, NK::Via | NK::Vert) {
                    // a hard search only needs to know whether anything is in the way
                    if !q.soft && !w.stretch_clear(l, n.pos, pos, q.net, width, 0.02) {
                        continue;
                    }
                    let blk = if q.soft { w.stretch_block(l, n.pos, pos, q.net, width, 0.02) } else { Vec::new() };
                    if !blk.is_empty() {
                        for p in blk {
                            extra += q.cross_cost.get(w.pieces[p as usize].net as usize).cloned().unwrap_or(10.0);
                            vics.push(p);
                        }
                    }
                }
                let nt = nb;
                let ntr = layer.cdt.tris[nt as usize];
                let nk = ntr.e.iter().position(|&x| x == e2).unwrap();
                let ngap = nc - g;
                let mut pen = q.taboo.get(&(l as u8, nt)).cloned().unwrap_or(0.0) + q.history.get(&(l as u8, nt)).cloned().unwrap_or(0.0);
                if self.gap_hs {
                    if let Some(gs) = &self.gaps {
                        pen += gs.edge_cost(q.net, &self.gap_own, l, layer.vobj[ed.a as usize], layer.vobj[ed.b as usize], pos, pos - n.pos);
                        pen += gs.area_cost(q.net, l, n.pos, pos);
                    }
                }
                let step = if g == coupled_gap { q.couple_factor } else { 1.0 };
                let gg = n.g + n.pos.dist(pos) * q.layer_cost[l] * step + ccost + extra + pen;
                let vic = self.add_victims(&vics);
                let node = Node { layer: l as u8, kind: NK::Gap, tri: nt, k: nk as u8, gap: ngap as u16, vtx: NONE, pos, g: gg, parent: ni, edge: e2, ins: ins as u32, left: ntr.v[nk], obj: NONE, vic, nvic: vics.len() as u16 };
                let key = self.key_of(&node);
                self.push(node, key);
            }
        }
        // target vertices of this triangle
        for j in 0..3 {
            let v = tr.v[j];
            if let Some(&obj) = self.tverts[l].get(&v) {
                if n.kind == NK::Vert && n.vtx == v {
                    continue;
                }
                let y = 2 * cy.start[j];
                let (c, vics) = match entry_x {
                    Some(_) => sep_cost(y, self),
                    None => sep_cost(y, self),
                };
                if c.is_finite() {
                    let p = w.vpos(l, v);
                    let mut c = c;
                    let mut vics = vics;
                    if !q.soft && !w.stretch_clear(l, n.pos, p, q.net, width, 0.02) {
                        continue;
                    }
                    let blk = if q.soft { w.stretch_block(l, n.pos, p, q.net, width, 0.02) } else { Vec::new() };
                    if !blk.is_empty() {
                        for pb in blk {
                            c += q.cross_cost.get(w.pieces[pb as usize].net as usize).cloned().unwrap_or(10.0);
                            vics.push(pb);
                        }
                    }
                    let gg = n.g + n.pos.dist(p) * q.layer_cost[l] + c;
                    let vic = self.add_victims(&vics);
                    let node = Node { layer: l as u8, kind: NK::GoalVert, tri: t, k: j as u8, gap: 0, vtx: v, pos: p, g: gg, parent: ni, edge: NONE, ins: 0, left: NONE, obj, vic, nvic: vics.len() as u16 };
                    let key = self.key_of(&node);
                    self.push(node, key);
                }
            }
        }
        // vias
        if n.kind == NK::Via {
            return;
        }
        let other_layers = q.allowed & !(1 << l);
        if other_layers == 0 && !q.plane_goal {
            return;
        }
        let near_term = n.g < 1.2 || self.h_raw(l, n.pos) < 1.2;
        let mut pts = match &q.via_spots {
            // a dogbone net: only its spots beside its pads that lie in this face
            Some(spots) => {
                let tr = w.layers[l].cdt.tri_ip(t);
                spots.iter().cloned().filter(|&x| { let ix = to_ip(x); orient(tr[0], tr[1], ix) > 0 && orient(tr[1], tr[2], ix) > 0 && orient(tr[2], tr[0], ix) > 0 }).collect()
            }
            None => self.via_candidates(l, t, n.pos, near_term),
        };
        let mut paired_vias: Vec<P2> = Vec::new();
        if q.couple != NONE && q.via_spots.is_none() {
            let c = &w.b.classes[w.b.nets[q.net as usize].class];
            let pitch = c.via_d + w.b.pair_gap[q.net as usize].max(c.clear) + 0.01;
            let tr2 = w.layers[l].cdt.tri_ip(t);
            for &pv in &self.partner_vias {
                if pv.dist(n.pos) > 3.0 {
                    continue;
                }
                for k in 0..12 {
                    let an = k as f64 * std::f64::consts::TAU / 12.0;
                    let x = pv + p2(an.cos(), an.sin()) * pitch;
                    let ix = to_ip(x);
                    if orient(tr2[0], tr2[1], ix) > 0 && orient(tr2[1], tr2[2], ix) > 0 && orient(tr2[2], tr2[0], ix) > 0 {
                        pts.push(x);
                        paired_vias.push(x);
                    }
                }
            }
        }
        for x in pts {
            // same face as the entry
            let mut ok = true;
            let mut vcost = 0.0;
            for (i, &c) in blocking.iter().enumerate() {
                if self.side(l, t, &cy, c, x) != entry_side[i] {
                    if !q.soft {
                        ok = false;
                        break;
                    }
                    vcost += 20.0;
                }
            }
            if !ok || vcost > 0.0 {
                continue;
            }
            if q.via_taboo.iter().any(|&(c, r)| c.dist(x) < r) {
                continue;
            }
            let mut vics: Vec<u32> = Vec::new();
            let mut vc = 0.0;
            if !q.soft && !w.stretch_clear(l, n.pos, x, q.net, width, 0.02) {
                continue;
            }
            if q.soft {
                vics.extend(w.stretch_block(l, n.pos, x, q.net, width, 0.02));
            }
            let key = ((x.x * 2e3).round() as i64, (x.y * 2e3).round() as i64);
            if q.soft {
                match w.via_block(q.net, x) {
                    None => continue,
                    Some(v) => {
                        for p in v {
                            if !vics.contains(&p) {
                                vics.push(p);
                            }
                        }
                    }
                }
                for &p in &vics {
                    vc += q.cross_cost.get(w.pieces[p as usize].net as usize).cloned().unwrap_or(10.0);
                }
            } else {
                let fits = match self.vfit.get(&key) {
                    Some(&f) => f,
                    None => {
                        let f = w.via_fits(q.net, x);
                        self.vfit.insert(key, f);
                        f
                    }
                };
                if !fits {
                    continue;
                }
                if !q.ribbon.is_empty() {
                    let partner = w.b.partner.get(q.net as usize).cloned().unwrap_or(NONE);
                    if partner != NONE {
                        let c = &w.b.classes[w.b.nets[partner as usize].class];
                        let pitch = c.via_d + w.b.pair_gap[partner as usize].max(c.clear) + 0.02;
                        let room = (0..12).any(|k| {
                            let an = k as f64 * std::f64::consts::TAU / 12.0;
                            w.via_fits(partner, x + p2(an.cos(), an.sin()) * pitch)
                        });
                        if !room {
                            continue;
                        }
                    }
                }
            }
            let vcost = if paired_vias.iter().any(|&y| y.dist(x) < 1e-9) { q.via_cost * q.couple_factor } else { q.via_cost };
            let gvc = self.gaps.as_ref().map_or(0.0, |gs| gs.via_cost(q.net, x));
            let base = n.g + n.pos.dist(x) * q.layer_cost[l] + vcost + vc + gvc;
            let vic = self.add_victims(&vics);
            let nvic = vics.len() as u16;
            if q.plane_goal {
                // a drop only reaches the plane inside the net's pour
                if !w.b.in_plane(q.net, x, w.b.via_size(q.net).0 / 2.0) {
                    continue;
                }
                let node = Node { layer: l as u8, kind: NK::GoalVia, tri: t, k: 0, gap: 0, vtx: NONE, pos: x, g: base, parent: ni, edge: NONE, ins: 0, left: NONE, obj: NONE, vic, nvic };
                let key = (l as u8, 8, key.0 as u32, key.1 as u32, 0);
                self.push(node, key);
                continue;
            }
            for l2 in 0..w.nsig() {
                if other_layers & (1 << l2) == 0 {
                    continue;
                }
                // a dogbone via joins its pads' layer to one other layer
                if q.via_spots.is_some() && q.via_spot_layers & ((1 << l) | (1 << l2)) == 0 {
                    continue;
                }
                let lay2 = &w.layers[l2];
                let t2 = match w.locate(l2, x) {
                    Ok(crate::cdt::Loc::Inside(t2)) => t2,
                    _ => continue,
                };
                if lay2.tflag[t2 as usize] & T_BLOCK != 0 {
                    continue;
                }
                let pen = q.taboo.get(&(l2 as u8, t2)).cloned().unwrap_or(0.0);
                let node = Node { layer: l2 as u8, kind: NK::Via, tri: t2, k: 0, gap: 0, vtx: NONE, pos: x, g: base + pen, parent: ni, edge: NONE, ins: 0, left: NONE, obj: NONE, vic, nvic };
                let key = self.key_of(&node);
                self.push(node, key);
            }
        }
    }

    fn add_victims(&mut self, v: &[u32]) -> u32 {
        let s = self.victims.len() as u32;
        self.victims.extend_from_slice(v);
        s
    }

    fn via_candidates(&self, l: usize, t: u32, from: P2, near_term: bool) -> Vec<P2> {
        let w = self.w;
        let tr = w.layers[l].cdt.tris[t as usize];
        let p: Vec<P2> = tr.v.iter().map(|&v| w.vpos(l, v)).collect();
        let (a, b, c) = (p[0], p[1], p[2]);
        let la = b.dist(c);
        let lb = c.dist(a);
        let lc = a.dist(b);
        let per = la + lb + lc;
        let inc = (a * la + b * lb + c * lc) * (1.0 / per);
        let mut out = vec![inc];
        // spread over the triangle, and a ring around the current position (escapes in any direction)
        for &(wa, wb, wc) in &[(0.6, 0.2, 0.2), (0.2, 0.6, 0.2), (0.2, 0.2, 0.6)] {
            out.push(a * wa + b * wb + c * wc);
        }
        if near_term {
            for &rad in &[0.5, 0.8] {
                for k in 0..8 {
                    let an = k as f64 * std::f64::consts::TAU / 8.0 + 0.3;
                    out.push(from + p2(an.cos(), an.sin()) * rad);
                }
            }
        }
        let goal = self.goal_pt(from);
        let dir = (goal - from).unit();
        for k in 1..=3 {
            out.push(from + dir * (0.45 * k as f64));
        }
        if !self.q.plane_goal {
            for &d in &[0.6, 0.9, 1.3, 2.0, 3.0] {
                let back = goal - dir * d;
                out.push(back);
                out.push(back + dir.perp() * 0.5);
                out.push(back - dir.perp() * 0.5);
            }
        } else {
            for k in 0..8 {
                let ang = k as f64 * std::f64::consts::TAU / 8.0;
                out.push(from + p2(ang.cos(), ang.sin()) * 0.55);
            }
        }
        let ia = to_ip(a);
        let ib = to_ip(b);
        let ic = to_ip(c);
        out.retain(|&x| {
            let ix = to_ip(x);
            orient(ia, ib, ix) > 0 && orient(ib, ic, ix) > 0 && orient(ic, ia, ix) > 0
        });
        out
    }

    fn reconstruct(&self, goal: u32) -> PathOut {
        let w = self.w;
        let mut chain = Vec::new();
        let mut i = goal;
        while i != NONE {
            chain.push(i);
            i = self.nodes[i as usize].parent;
        }
        chain.reverse();
        let mut segs: Vec<Seg> = Vec::new();
        let mut vias: Vec<P2> = Vec::new();
        let mut victims: Vec<u32> = Vec::new();
        let mut cur: Option<Seg> = None;
        let mut target = NONE;
        let mut cross_tris = Vec::new();
        for &ni in &chain {
            let n = self.nodes[ni as usize];
            for k in 0..n.nvic as usize {
                victims.push(self.victims[n.vic as usize + k]);
            }
            if n.nvic > 0 && n.parent != NONE {
                let pn = self.nodes[n.parent as usize];
                cross_tris.push((pn.layer, pn.tri));
            }
            let l = n.layer as usize;
            match n.kind {
                NK::Gap => {
                    if n.parent == NONE {
                        cur = Some(Seg { layer: l, a: End::Obj(Term { obj: n.obj, vtx: NONE }), b: End::NewVia(0), cx: vec![(n.edge, n.ins as usize, n.left)] });
                    } else {
                        cur.as_mut().unwrap().cx.push((n.edge, n.ins as usize, n.left));
                    }
                }
                NK::Vert => {
                    cur = Some(Seg { layer: l, a: End::Obj(Term { obj: n.obj, vtx: n.vtx }), b: End::NewVia(0), cx: vec![] });
                }
                NK::Via => {
                    let vi = vias.len();
                    vias.push(n.pos);
                    let mut s = cur.take().unwrap();
                    s.b = End::NewVia(vi);
                    segs.push(s);
                    cur = Some(Seg { layer: l, a: End::NewVia(vi), b: End::NewVia(0), cx: vec![] });
                }
                NK::GoalHull => {
                    let mut s = cur.take().unwrap();
                    s.cx.push((n.edge, n.ins as usize, n.left));
                    s.b = End::Obj(Term { obj: n.obj, vtx: NONE });
                    segs.push(s);
                    target = n.obj;
                }
                NK::GoalVert => {
                    let mut s = cur.take().unwrap();
                    s.b = End::Obj(Term { obj: n.obj, vtx: n.vtx });
                    segs.push(s);
                    target = n.obj;
                }
                NK::GoalVia => {
                    let vi = vias.len();
                    vias.push(n.pos);
                    let mut s = cur.take().unwrap();
                    s.b = End::NewVia(vi);
                    segs.push(s);
                }
            }
        }
        victims.sort();
        victims.dedup();
        let _ = w;
        let _ = Kind::Pad;
        PathOut { segs, vias, cost: self.nodes[goal as usize].g, victims, target, cross_tris }
    }
}
