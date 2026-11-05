//! Spring relaxation of the finished board: all free wires and vias move together (gridroute::relax, on Metal
//! where available). Each wire is a chain of nodes under tension; each via is a free node pulled along each of
//! its wires; every pair of copper elements keeps the least distance its rules require (clearance, hole and
//! drill spacing, board edge, keepouts, via rule areas), so vias pull their wires straight and push other
//! nets' wires aside while those pull tight around them. The relaxed geometry is then threaded back into the
//! topology, realised taut and exact-checked; nets that do not fit are routed again, and the board keeps the
//! result only if it is complete, clean and shorter.
//!
//! Differential pair members, their vias and vias tied into pads stay where they are, as obstacles.

use crate::cdt::NONE;
use crate::geom::*;
use crate::model::{Board, Kind};
use crate::router::{wire_len, Router};
use crate::search::{End, Search};
use crate::world::{clr_in, Attr, Term, MARGIN};
use gridroute::relax::{Params, Relax};
use rayon::prelude::*;

/// Drill "layer": drill-to-drill spacing applies between any nets.
const DRILL: u32 = 1 << 30;
/// Via rule areas: only vias keep out.
const VIARULE: u32 = 1 << 29;
/// Wire nodes per mm at the start (the chains re-space themselves).
const NODE_STEP: f64 = 0.2;

#[derive(Clone)]
struct El {
    layers: u32,
    attr: Attr,
    /// Copper radius around the core (drill radius for a drill element).
    r: f64,
    /// Drill radius of a via's copper (hole-to-copper spacing), else 0.
    rd: f64,
    /// Copper that a via of its own net stays out of (pads, not zone copper or virtual junctions).
    own_pad: bool,
    /// A pad's object (blocked in the topology: even its own net's wires only start or end in it).
    pad: u32,
    /// A wire stretch's terminal objects.
    ends: [u32; 2],
}

struct Chain {
    piece: u32,
    net: u32,
    layer: usize,
    width: f64,
    a: Term,
    b: Term,
    nodes: Vec<u32>,
}

struct Model {
    r: Relax,
    el: Vec<El>,
    el_nodes: Vec<Vec<u32>>,
    mobile: Vec<u32>,
    chains: Vec<Chain>,
    /// (via object, its node)
    vias: Vec<(u32, u32)>,
    /// Same-net (via, pad) pairs that overlap at the start (via in pad): left alone.
    in_pad: std::collections::HashSet<(u32, u32)>,
    max_gap: f64,
    /// Per element: its chain (wire stretches), or NONE.
    el_chain: Vec<u32>,
}

impl Model {
    fn node(&mut self, p: P2, kind: u32) -> u32 {
        self.r.pos.push([p.x as f32, p.y as f32]);
        self.r.kind.push(kind);
        self.r.leash.push([0.0, 0.0, -1.0, 0.0]);
        (self.r.pos.len() - 1) as u32
    }
    fn element(&mut self, nodes: Vec<u32>, e: El, mobile: bool) -> u32 {
        let id = self.el.len() as u32;
        self.el.push(e);
        self.el_nodes.push(nodes);
        if mobile {
            self.mobile.push(id);
        }
        id
    }
    fn pos(&self, n: u32) -> P2 {
        let q = self.r.pos[n as usize];
        p2(q[0] as f64, q[1] as f64)
    }
    /// Core distance of two elements (the first has at most two nodes).
    fn core_dist(&self, e: u32, f: u32) -> f64 {
        let a = &self.el_nodes[e as usize];
        let (pa, qa) = (self.pos(a[0]), self.pos(*a.last().unwrap()));
        let b = &self.el_nodes[f as usize];
        if b.len() <= 2 {
            seg_seg_dist(pa, qa, self.pos(b[0]), self.pos(*b.last().unwrap()))
        } else {
            let pts: Vec<P2> = b.iter().map(|&n| self.pos(n)).collect();
            seg_hull_dist(&pts, pa, qa)
        }
    }
    fn bbox(&self, e: u32) -> (P2, P2) {
        let mut lo = p2(f64::MAX, f64::MAX);
        let mut hi = p2(f64::MIN, f64::MIN);
        for &n in &self.el_nodes[e as usize] {
            let p = self.pos(n);
            lo = p2(lo.x.min(p.x), lo.y.min(p.y));
            hi = p2(hi.x.max(p.x), hi.y.max(p.y));
        }
        let r = self.el[e as usize].r;
        (lo - p2(r, r), hi + p2(r, r))
    }
}

/// The least core distance between movable element e and element f, or None when they may touch.
fn need(bd: &Board, m: &Model, e: u32, f: u32) -> Option<f64> {
    let (a, b) = (&m.el[e as usize], &m.el[f as usize]);
    if a.layers & b.layers == 0 || m.el_nodes[e as usize].iter().any(|n| m.el_nodes[f as usize].contains(n)) {
        return None;
    }
    if a.layers == DRILL || b.layers == DRILL {
        return Some(a.r + b.r + bd.hole_to_hole + MARGIN);
    }
    if b.layers == VIARULE {
        // (the exact check keeps vias 0.01 off a rule area's edge)
        return if a.attr.kind == Kind::Via { Some(a.r + 0.01 + MARGIN) } else { None };
    }
    let gap = match clr_in(bd, &a.attr, &b.attr) {
        Some(g) => g,
        None => {
            // a via stays out of its own net's pads (unless it started in one)
            let via_pad = |x: &El, y: &El| x.attr.kind == Kind::Via && y.own_pad && x.attr.net == y.attr.net;
            if (via_pad(a, b) || via_pad(b, a)) && !m.in_pad.contains(&(e.min(f), e.max(f))) {
                return Some(a.r + b.r);
            }
            // a wire crosses no pad of its own net but those it starts or ends in
            if b.pad != NONE && a.ends[0] != NONE && !a.ends.contains(&b.pad) {
                return Some(0.002);
            }
            return None;
        }
    };
    let mut n = a.r + b.r + gap;
    let copper = |x: &El| x.attr.kind != Kind::Outline && x.attr.kind != Kind::Keepout;
    if a.rd > 0.0 && copper(b) {
        n = n.max(a.rd + b.r + bd.hole_clear + MARGIN);
    }
    if b.rd > 0.0 && copper(a) {
        n = n.max(a.r + b.rd + bd.hole_clear + MARGIN);
    }
    Some(n)
}

impl Router {
    fn spring_model(&self, hold: &std::collections::HashSet<u32>) -> Model {
        let w = &self.w;
        let nsig = w.nsig();
        let sig = (1u32 << nsig) - 1;
        // pair members and held nets stay as they are
        let pair = |net: u32| w.b.partner.get(net as usize).map_or(false, |&q| q != NONE) || hold.contains(&net);
        let mut m = Model { r: Relax::default(), el: Vec::new(), el_nodes: Vec::new(), mobile: Vec::new(), chains: Vec::new(), vias: Vec::new(), in_pad: Default::default(), max_gap: 0.0, el_chain: Vec::new() };
        // movable pieces: every live wire of a net that is not a pair member
        let movable_piece = |p: usize| {
            let pc = &w.pieces[p];
            // (a drop at a neck-down width stays the straight segment it is)
            pc.alive && !pc.frozen && !pair(pc.net) && pc.poly.len() >= 2 && pc.width >= w.b.nets[pc.net as usize].wmin[pc.layer as usize] - 1e-9
        };
        // movable vias: signal and plane vias whose wires all move
        let mut via_node: crate::fx::FxMap<u32, u32> = Default::default();
        for (o, ob) in w.objs.iter().enumerate() {
            // (a dogbone net's via stays beside its pad)
            if !ob.alive || ob.s.kind != Kind::Via || pair(ob.s.net) || self.ties.iter().any(|t| t.0 == o as u32 || t.1 == o as u32) || w.b.dogbone.get(ob.s.net as usize).copied().unwrap_or(false) {
                continue;
            }
            let live: Vec<u32> = w.obj_pieces[o].iter().cloned().filter(|&p| w.pieces[p as usize].alive).collect();
            if live.is_empty() || live.iter().any(|&p| !movable_piece(p as usize)) {
                continue;
            }
            let n = m.node(ob.s.center, 2);
            if ob.plane {
                m.r.leash[n as usize] = [ob.s.center.x as f32, ob.s.center.y as f32, self.plane_leash(ob.s.net, ob.s.center, ob.s.r) as f32, 0.0];
            }
            via_node.insert(o as u32, n);
            m.vias.push((o as u32, n));
        }
        // fixed copper and boundaries
        for (o, ob) in w.objs.iter().enumerate() {
            if !ob.alive || via_node.contains_key(&(o as u32)) {
                continue;
            }
            let s = &ob.s;
            let attr = w.obj_attr(o as u32);
            match s.kind {
                Kind::Outline | Kind::Keepout => {
                    let ring: Vec<u32> = s.pts.iter().map(|&p| m.node(p, 0)).collect();
                    let k = ring.len();
                    for i in 0..k {
                        m.element(vec![ring[i], ring[(i + 1) % k]], El { layers: s.layers & sig, attr, r: 0.0, rd: 0.0, own_pad: false, pad: NONE, ends: [NONE; 2] }, false);
                    }
                }
                _ => {
                    let nodes: Vec<u32> = s.pts.iter().map(|&p| m.node(p, 0)).collect();
                    let rd = if s.kind == Kind::Via { s.drill / 2.0 } else { 0.0 };
                    let own_pad = s.kind == Kind::Pad && !crate::world::fixed_copper(&s.name) && !s.name.starts_with("junction:");
                    m.element(nodes, El { layers: s.layers & sig, attr, r: s.r, rd, own_pad, pad: if s.kind == Kind::Pad && s.pts.len() >= 2 { o as u32 } else { NONE }, ends: [NONE; 2] }, false);
                    if s.drill > 0.0 {
                        let c = m.node(s.center, 0);
                        m.element(vec![c], El { layers: DRILL, attr, r: s.drill / 2.0, rd: 0.0, own_pad: false, pad: NONE, ends: [NONE; 2] }, false);
                    }
                }
            }
        }
        for poly in &w.b.via_keepouts {
            let ring: Vec<u32> = poly.iter().map(|&p| m.node(p, 0)).collect();
            let k = ring.len();
            let attr = Attr { kind: Kind::Keepout, net: NONE, class: NONE, clear: 0.0 };
            for i in 0..k {
                m.element(vec![ring[i], ring[(i + 1) % k]], El { layers: VIARULE, attr, r: 0.0, rd: 0.0, own_pad: false, pad: NONE, ends: [NONE; 2] }, false);
            }
        }
        // the movable vias' copper (on every layer, and in via rule areas) and drills
        for &(o, n) in &m.vias.clone() {
            let s = &w.objs[o as usize].s;
            let c = w.b.nets[s.net as usize].class;
            let attr = Attr { kind: Kind::Via, net: s.net, class: c as u32, clear: s.clear };
            m.element(vec![n], El { layers: (s.layers & sig) | VIARULE, attr, r: s.r, rd: s.drill / 2.0, own_pad: false, pad: NONE, ends: [NONE; 2] }, true);
            m.element(vec![n], El { layers: DRILL, attr, r: s.drill / 2.0, rd: 0.0, own_pad: false, pad: NONE, ends: [NONE; 2] }, true);
        }
        // wires: movable ones as chains, the rest fixed
        for p in 0..w.pieces.len() {
            let pc = &w.pieces[p];
            if !pc.alive || pc.poly.len() < 2 {
                continue;
            }
            let attr = w.net_attr(pc.net, pc.layer as usize);
            let el = El { layers: 1 << pc.layer, attr, r: pc.width / 2.0, rd: 0.0, own_pad: false, pad: NONE, ends: [pc.a.obj, pc.b.obj] };
            if !movable_piece(p) {
                let nodes: Vec<u32> = pc.poly.iter().map(|&q| m.node(q, 0)).collect();
                for k in 0..nodes.len() - 1 {
                    m.element(vec![nodes[k], nodes[k + 1]], el.clone(), false);
                }
                continue;
            }
            let end = |t: Term, at: P2, m: &mut Model| via_node.get(&t.obj).cloned().unwrap_or_else(|| m.node(at, 0));
            let mut nodes = vec![end(pc.a, pc.poly[0], &mut m)];
            for k in 0..pc.poly.len() - 1 {
                let (u, v) = (pc.poly[k], pc.poly[k + 1]);
                let parts = ((u.dist(v) / NODE_STEP).ceil() as usize).max(1);
                for j in 1..parts {
                    let q = u.lerp(v, j as f64 / parts as f64);
                    nodes.push(m.node(q, 1));
                }
                if k + 1 < pc.poly.len() - 1 {
                    nodes.push(m.node(v, 1));
                }
            }
            nodes.push(end(pc.b, *pc.poly.last().unwrap(), &mut m));
            let ci = m.chains.len() as u32;
            for k in 0..nodes.len() - 1 {
                let e = m.element(vec![nodes[k], nodes[k + 1]], el.clone(), true);
                m.el_chain.resize(e as usize + 1, NONE);
                m.el_chain[e as usize] = ci;
            }
            m.chains.push(Chain { piece: p as u32, net: pc.net, layer: pc.layer as usize, width: pc.width, a: pc.a, b: pc.b, nodes });
        }
        m.el_chain.resize(m.el.len(), NONE);
        // tension neighbours: chain neighbours, and each via's first node along each of its wires
        let nn = m.r.pos.len();
        let mut tn: Vec<Vec<u32>> = vec![Vec::new(); nn];
        for ch in &m.chains {
            let k = ch.nodes.len();
            for i in 1..k - 1 {
                tn[ch.nodes[i] as usize] = vec![ch.nodes[i - 1], ch.nodes[i + 1]];
            }
            for (end, next) in [(ch.nodes[0], ch.nodes[1]), (ch.nodes[k - 1], ch.nodes[k - 2])] {
                if m.r.kind[end as usize] == 2 {
                    tn[end as usize].push(next);
                }
            }
        }
        m.r.tn_off = vec![0];
        for t in &tn {
            m.r.tn.extend_from_slice(t);
            m.r.tn_off.push(m.r.tn.len() as u32);
        }
        // elements per node, for the movable nodes
        let mut ne: Vec<Vec<u32>> = vec![Vec::new(); nn];
        for &e in &m.mobile {
            for &n in &m.el_nodes[e as usize] {
                if m.r.kind[n as usize] != 0 {
                    ne[n as usize].push(e);
                }
            }
        }
        m.r.ne_off = vec![0];
        for v in &ne {
            m.r.ne.extend_from_slice(v);
            m.r.ne_off.push(m.r.ne.len() as u32);
        }
        m.r.el_off = vec![0];
        for v in &m.el_nodes {
            m.r.el_node.extend_from_slice(v);
            m.r.el_off.push(m.r.el_node.len() as u32);
        }
        m.max_gap = w.max_clear.max(w.b.hole_clear).max(w.b.hole_to_hole).max(w.b.edge_clear) + MARGIN + 0.01;
        // vias that start inside a pad of their own net stay allowed there
        let pads: Vec<u32> = (0..m.el.len() as u32).filter(|&f| m.el[f as usize].own_pad).collect();
        for &(_, n) in &m.vias {
            let e = m.r.ne[m.r.ne_off[n as usize] as usize];
            let p = m.pos(n);
            for &f in &pads {
                let (lo, hi) = m.bbox(f);
                if p.x < lo.x - 1.0 || p.x > hi.x + 1.0 || p.y < lo.y - 1.0 || p.y > hi.y + 1.0 {
                    continue;
                }
                if m.el[f as usize].attr.net == m.el[e as usize].attr.net && m.core_dist(e, f) < m.el[e as usize].r + m.el[f as usize].r {
                    m.in_pad.insert((e.min(f), e.max(f)));
                }
            }
        }
        m
    }

    /// How far a plane via may move and still reach its pour: the largest radius whose circle (sampled)
    /// stays in the plane.
    fn plane_leash(&self, net: u32, c: P2, r: f64) -> f64 {
        let mut ok = 0.0;
        for k in 1..=10 {
            let rad = k as f64 * 0.1;
            let all = (0..8).all(|j| {
                let an = j as f64 * std::f64::consts::TAU / 8.0;
                self.w.b.in_plane(net, c + p2(an.cos(), an.sin()) * rad, r)
            });
            if !all {
                break;
            }
            ok = rad;
        }
        ok
    }

    /// Constraint lists for the movable elements at the current positions, with `skin` to spare.
    fn spring_constraints(&self, m: &mut Model, skin: f64) {
        const CELL: f64 = 0.5;
        let (lo, hi) = (p2(-1.0, -1.0), p2(self.w.b.size[0] + 1.0, self.w.b.size[1] + 1.0));
        let nx = ((hi.x - lo.x) / CELL).ceil() as usize + 1;
        let ny = ((hi.y - lo.y) / CELL).ceil() as usize + 1;
        let cell = |p: P2| (((p.x - lo.x) / CELL).clamp(0.0, (nx - 1) as f64) as usize, ((p.y - lo.y) / CELL).clamp(0.0, (ny - 1) as f64) as usize);
        let mut grid: Vec<Vec<u32>> = vec![Vec::new(); nx * ny];
        for e in 0..m.el.len() as u32 {
            let (a, b) = m.bbox(e);
            let ((i0, j0), (i1, j1)) = (cell(a), cell(b));
            for j in j0..=j1 {
                for i in i0..=i1 {
                    grid[j * nx + i].push(e);
                }
            }
        }
        let reach = m.max_gap + skin;
        let bd = &self.w.b;
        let mm = &*m;
        let lists: Vec<Vec<(u32, f32)>> = m
            .mobile
            .par_iter()
            .map(|&e| {
                let m = mm;
                let (a, b) = m.bbox(e);
                let ((i0, j0), (i1, j1)) = (cell(a - p2(reach, reach)), cell(b + p2(reach, reach)));
                let mut cand: Vec<u32> = Vec::new();
                for j in j0..=j1 {
                    for i in i0..=i1 {
                        cand.extend_from_slice(&grid[j * nx + i]);
                    }
                }
                cand.sort_unstable();
                cand.dedup();
                let mut out = Vec::new();
                let (qa, qb) = (a - p2(reach, reach), b + p2(reach, reach));
                for f in cand {
                    if f == e {
                        continue;
                    }
                    let (fa, fb) = m.bbox(f);
                    if fa.x > qb.x || fb.x < qa.x || fa.y > qb.y || fb.y < qa.y {
                        continue;
                    }
                    if let Some(need) = need(bd, m, e, f) {
                        if m.core_dist(e, f) < need + skin {
                            // a hair over the rule, so the relaxed geometry passes the exact checks as it is
                            out.push((f, (need + 0.002) as f32));
                        }
                    }
                }
                out
            })
            .collect();
        let ne = m.el.len();
        let mut per: Vec<&[(u32, f32)]> = vec![&[]; ne];
        for (k, &e) in m.mobile.iter().enumerate() {
            per[e as usize] = &lists[k];
        }
        m.r.c_off = Vec::with_capacity(ne + 1);
        m.r.c_el.clear();
        m.r.c_need.clear();
        m.r.c_off.push(0);
        for l in per {
            for &(f, n) in l {
                m.r.c_el.push(f);
                m.r.c_need.push(n);
            }
            m.r.c_off.push(m.r.c_el.len() as u32);
        }
        m.r.c_version += 1;
    }

    /// Chains and nodes the solve left short of a gap (by more than the exact checks' slack).
    fn spring_short(m: &Model) -> (std::collections::HashSet<u32>, std::collections::HashSet<u32>) {
        let mut chains = std::collections::HashSet::new();
        let mut nodes = std::collections::HashSet::new();
        for &e in &m.mobile {
            for c in m.r.c_off[e as usize]..m.r.c_off[e as usize + 1] {
                let f = m.r.c_el[c as usize];
                if m.r.c_need[c as usize] as f64 - m.core_dist(e, f) <= 0.0015 {
                    continue;
                }
                for x in [e, f] {
                    if m.el_chain[x as usize] != NONE {
                        chains.insert(m.el_chain[x as usize]);
                    }
                    for &n in &m.el_nodes[x as usize] {
                        if m.r.kind[n as usize] == 2 {
                            nodes.insert(n);
                        }
                    }
                }
            }
        }
        (chains, nodes)
    }

    fn chain_len(m: &Model) -> f64 {
        m.chains.iter().map(|c| c.nodes.windows(2).map(|s| m.pos(s[0]).dist(m.pos(s[1]))).sum::<f64>()).sum()
    }

    /// The spring relaxation phase (see the module notes). Returns whether the board took the result. Nets
    /// that cannot be joined up again after a relaxation are held as they were in another, up to twice.
    pub fn spring(&mut self) -> bool {
        self.spring_tries(3)
    }

    /// The spring relaxation in at most `tries` solves (each holding the nets the one before could not join up).
    pub fn spring_tries(&mut self, tries: usize) -> bool {
        let mut hold = std::collections::HashSet::new();
        for _ in 0..tries {
            match self.spring_once(&hold) {
                Ok(kept) => return kept,
                Err(nets) => {
                    if self.cfg.verbose {
                        eprintln!("  spring: again with {} more nets held as they are", nets.len());
                    }
                    hold.extend(nets);
                }
            }
        }
        false
    }

    /// One relaxation with `hold` fixed: Ok(kept), or the nets that could not be joined up again.
    fn spring_once(&mut self, hold: &std::collections::HashSet<u32>) -> Result<bool, Vec<u32>> {
        let t0 = self.elapsed();
        let _t = std::time::Instant::now();
        let mut m = self.spring_model(hold);
        self.prof_add("spring: model", _t);
        let t_solve = std::time::Instant::now();
        if m.vias.is_empty() && m.chains.is_empty() {
            return Ok(false);
        }
        let env = |k: &str, d: f64| std::env::var(k).ok().and_then(|s| s.parse().ok()).unwrap_or(d);
        let batch = env("TOPOR_SPRING_BATCH", 32.0) as u32;
        let max_iters = env("TOPOR_SPRING_ITERS", 9600.0) as usize;
        let dmax = env("TOPOR_SPRING_DMAX", 0.006) as f32;
        let skin = env("TOPOR_SPRING_SKIN", 0.4);
        // a wall-clock budget: past it the solve goes straight to the polish, and stops after it
        let budget = env("TOPOR_SPRING_SECONDS", 30.0);
        let mut prm = Params { iters: batch, sub: 2, alpha_wire: dmax * 0.5, alpha_via: dmax, beta: 0.25, dmax, tol: 0.001 };
        let start_len = Self::chain_len(&m);
        let mut last = start_len;
        let mut gpu = false;
        let (mut batches, mut builds, mut redone) = (0, 0, 0);
        let (mut its, mut last_check) = (0usize, 0usize);
        let (mut t_build, mut t_run) = (0.0, 0.0);
        // Verlet lists: rebuilt once some node has moved half the skin since the last build; a batch that
        // went further than that may have missed a pair, so it runs again from where it started on fresh lists
        let moved = |m: &Model, at: &[[f32; 2]]| m.r.pos.iter().zip(at).map(|(a, b)| (a[0] - b[0]).hypot(a[1] - b[1])).fold(0f32, f32::max) as f64;
        let tb = std::time::Instant::now();
        self.spring_constraints(&mut m, skin);
        t_build += tb.elapsed().as_secs_f64();
        builds += 1;
        let mut at_build = m.r.pos.clone();
        let mut fresh = true;
        let mut polish = 0;
        // polish: constraints only, until every gap holds (or a jam will not come apart)
        let residual = |m: &Model| -> f64 {
            m.mobile
                .par_iter()
                .map(|&e| (m.r.c_off[e as usize]..m.r.c_off[e as usize + 1]).map(|c| m.r.c_need[c as usize] as f64 - m.core_dist(e, m.r.c_el[c as usize])).fold(0.0, f64::max))
                .reduce(|| 0.0, f64::max)
        };
        while polish < 40 {
            if self.elapsed() - t0 > budget * 1.5 {
                break;
            }
            if self.elapsed() - t0 > budget && prm.alpha_via != 0.0 {
                prm.alpha_wire = 0.0;
                prm.alpha_via = 0.0;
                prm.beta = 0.0;
                prm.sub = 4;
            }
            let start = m.r.pos.clone();
            let tr = std::time::Instant::now();
            gpu = m.r.run(&prm);
            t_run += tr.elapsed().as_secs_f64();
            if moved(&m, &at_build) > skin / 2.0 && !(fresh && prm.iters == 1) {
                m.r.pos = start;
                if fresh {
                    prm.iters = (prm.iters / 2).max(1);
                }
                let tb = std::time::Instant::now();
                self.spring_constraints(&mut m, skin);
                t_build += tb.elapsed().as_secs_f64();
                builds += 1;
                redone += 1;
                at_build = m.r.pos.clone();
                fresh = true;
                continue;
            }
            fresh = false;
            batches += 1;
            its += prm.iters as usize;
            prm.iters = (prm.iters * 2).min(batch);
            // rebuild before the next batch can outrun the lists, rather than redo it
            if moved(&m, &at_build) > skin / 4.0 {
                let tb = std::time::Instant::now();
                self.spring_constraints(&mut m, skin);
                t_build += tb.elapsed().as_secs_f64();
                builds += 1;
                at_build = m.r.pos.clone();
                fresh = true;
            }
            if prm.alpha_via == 0.0 {
                polish += 1;
                if polish >= 2 && residual(&m) < 0.0015 {
                    break;
                }
            } else if its >= max_iters || its >= last_check + 128 && {
                let l = Self::chain_len(&m);
                let done = last - l < 2e-4 * l;
                last = l;
                last_check = its;
                done
            } {
                prm.alpha_wire = 0.0;
                prm.alpha_via = 0.0;
                prm.beta = 0.0;
                prm.sub = 4;
            }
        }
        let relaxed_len = Self::chain_len(&m);
        self.prof_add("spring: solve (lists + Metal)", t_solve);
        // what the solve left unmet
        let (mut short, mut worst) = (0usize, 0.0f64);
        for &e in &m.mobile {
            for c in m.r.c_off[e as usize]..m.r.c_off[e as usize + 1] {
                let f = m.r.c_el[c as usize];
                let gap = m.r.c_need[c as usize] as f64 - m.core_dist(e, f);
                if gap > 5e-4 {
                    short += 1;
                    worst = worst.max(gap);
                    if gap > 0.01 && std::env::var("TOPOR_SPRINGDBG").is_ok() {
                        let desc = |x: u32| {
                            let el = &m.el[x as usize];
                            let at = m.pos(m.el_nodes[x as usize][0]);
                            format!("{:?} {} n{} layers {:x} r {:.3} at ({:.3},{:.3}) kinds {:?}", el.attr.kind, self.w.b.nets.get(el.attr.net as usize).map_or("-", |n| n.name.as_str()), m.el_nodes[x as usize].len(), el.layers, el.r, at.x, at.y, m.el_nodes[x as usize].iter().map(|&n| m.r.kind[n as usize]).collect::<Vec<_>>())
                        };
                        eprintln!("  short by {:.4}: {} | {}", gap, desc(e), desc(f));
                    }
                }
            }
        }
        if self.cfg.verbose && short > 0 {
            eprintln!("  spring: {} gaps short after the solve (worst by {:.4} mm)", short, worst);
        }
        if self.cfg.verbose {
            eprintln!(
                "[{:.1}s] spring relaxation: {} nodes, {} vias, {} wires, {} iterations in {} batches on {} ({} list builds in {:.2}s, {} batches redone, {:.2}s solve): chains {:.1} -> {:.1} mm",
                self.elapsed(), m.r.pos.len(), m.vias.len(), m.chains.len(), its, batches, if gpu { "Metal" } else { "CPU" }, builds, t_build, redone, t_run, start_len, relaxed_len
            );
        }
        let _t = std::time::Instant::now();
        let ok = self.spring_write_back(&m);
        self.prof_add("spring: write-back", _t);
        if self.cfg.verbose {
            eprintln!("[{:.1}s] spring relaxation {} ({:.1}s)", self.elapsed(), if ok.as_ref().map_or(false, |&k| k) { "kept" } else { "dropped" }, self.elapsed() - t0);
        }
        ok
    }

    /// Thread the relaxed geometry back into the topology: vias to their new places, every movable wire
    /// re-threaded along its relaxed chain and realised taut. Nets that do not fit are routed again; the board
    /// keeps the result only if it is as complete and clean as before and shorter.
    fn spring_write_back(&mut self, m: &Model) -> Result<bool, Vec<u32>> {
        let before_def = self.total_deficit();
        // nets open before the relaxation: the write-back does not try them again (they cost full searches and
        // will not come in now)
        let was_open: std::collections::HashSet<u32> = self.incomplete_nets().into_iter().collect();
        let _t = std::time::Instant::now();
        let before_viol = self.w.check_all().len();
        self.prof_add("spring wb: exact checks", _t);
        let before_len = wire_len(&self.w);
        let _t = std::time::Instant::now();
        let snap = (self.w.clone(), self.ties.clone(), self.history.clone(), self.junctions.clone());
        self.prof_add("spring wb: snapshot", _t);
        let t_thread = std::time::Instant::now();
        let restore = |r: &mut Router, s: (crate::world::World, Vec<(u32, u32)>, Vec<f64>, Vec<u32>)| {
            r.w = s.0;
            r.ties = s.1;
            r.history = s.2;
            r.junctions = s.3;
        };
        // 1. lift the wires that moved
        for ch in &m.chains {
            self.w.remove_piece(ch.piece);
        }
        // what the solve left short: those vias stay where they were, those wires are found again
        let (bad_chain, bad_node) = Self::spring_short(m);
        let dbg = std::env::var("TOPOR_SPRINGDBG").is_ok();
        if dbg {
            for &ci in &bad_chain {
                let ch = &m.chains[ci as usize];
                eprintln!("  left short: {} L{} {} -> {}", self.w.b.nets[ch.net as usize].name, ch.layer, self.w.objs[ch.a.obj as usize].s.name, self.w.objs[ch.b.obj as usize].s.name);
            }
        }
        // 2. vias to their relaxed places
        let mut map: crate::fx::FxMap<u32, u32> = Default::default();
        for &(o, n) in &m.vias {
            let x = m.pos(n);
            let ob = &self.w.objs[o as usize];
            if ob.s.center.dist(x) < 1e-4 || bad_node.contains(&n) {
                map.insert(o, o);
                continue;
            }
            let net = ob.s.net;
            if self.w.remove_via(o).is_err() {
                restore(self, snap);
                return Ok(false);
            }
            match self.w.add_via(net, x) {
                Ok(v) => {
                    map.insert(o, v);
                }
                Err(_) => {
                    restore(self, snap);
                    return Ok(false);
                }
            }
        }
        // 3. the wires along their relaxed chains (not those left short, or ending at a via that stayed)
        let mut open: Vec<u32> = Vec::new();
        // wires to find again between their own ends on their layer: (net, layer, width, ends)
        let mut pending: Vec<(u32, usize, f64, Term, Term)> = Vec::new();
        let mut relaxed: crate::fx::FxMap<u32, Vec<P2>> = Default::default();
        let (mut skipped, mut unthreaded) = (0, 0);
        for (ci, ch) in m.chains.iter().enumerate() {
            let l = ch.layer;
            let term = |t: Term, w: &crate::world::World| match map.get(&t.obj) {
                Some(&v) => Term { obj: v, vtx: w.objs[v as usize].lv[l] },
                None => t,
            };
            let (ta, tb) = (term(ch.a, &self.w), term(ch.b, &self.w));
            let ends_stay = [ch.nodes[0], *ch.nodes.last().unwrap()].iter().any(|n| bad_node.contains(n));
            if bad_chain.contains(&(ci as u32)) || ends_stay {
                skipped += 1;
                pending.push((ch.net, l, ch.width, ta, tb));
                continue;
            }
            let mut pts: Vec<P2> = ch.nodes.iter().map(|&n| m.pos(n)).collect();
            // a polygon pad's terminal point is inside it: drop the chain's nodes that are too (nothing holds
            // them there), so the wire leaves the pad once, straight from its point
            let inside = |t: Term, q: P2| t.vtx == NONE && { let o = &self.w.objs[t.obj as usize].s; o.pts.len() >= 3 && point_in_poly(&o.pts, q) };
            let k0 = (1..pts.len() - 1).take_while(|&k| inside(ta, pts[k])).count();
            let k1 = (1..pts.len() - 1).rev().take_while(|&k| k > k0 && inside(tb, pts[k])).count();
            pts.drain(pts.len() - 1 - k1..pts.len() - 1);
            pts.drain(1..1 + k0);
            // exact ends: the terminals' own points
            pts[0] = self.w.term_point(l, ta);
            *pts.last_mut().unwrap() = self.w.term_point(l, tb);
            let mut poly = simplify(&pts, 0.001);
            let mut cx = self.w.thread_poly(l, &poly, ta, tb);
            if cx.is_none() {
                poly = pts;
                cx = self.w.thread_poly(l, &poly, ta, tb);
            }
            match cx {
                Some(cx) => {
                    let pid = self.w.add_piece(ch.net, l, ch.width, ta, tb, &cx);
                    self.w.seed_poly(pid, poly.clone());
                    relaxed.insert(pid, poly);
                }
                None => {
                    if dbg {
                        let d = |t: Term| { let o = &self.w.objs[t.obj as usize].s; format!("{} {:?} vtx {} c ({:.3},{:.3}) pts {:?}", o.name, o.kind, t.vtx != NONE, o.center.x, o.center.y, o.pts.iter().map(|q| format!("({:.3},{:.3})", q.x, q.y)).collect::<Vec<_>>()) };
                        eprintln!("  unthreaded L{} {} -> {}; poly {:?}", l, d(ta), d(tb), poly.iter().take(6).map(|q| format!("({:.3},{:.3})", q.x, q.y)).collect::<Vec<_>>());
                    }
                    unthreaded += 1;
                    pending.push((ch.net, l, ch.width, ta, tb));
                }
            }
        }
        let added: std::collections::HashSet<u32> = relaxed.keys().cloned().collect();
        self.w.realize_dirty();
        // the taut realisation can cut a clearance disc where a wire crosses an edge obliquely: there the
        // relaxed polyline (clear by construction) stays, and so does that of every wire it then meets
        let mut kept = 0;
        for _ in 0..12 {
            let mut add: Vec<u32> = Vec::new();
            for v in self.w.check_all() {
                for &p in &v.pieces {
                    if relaxed.contains_key(&p) && !add.contains(&p) {
                        add.push(p);
                    }
                }
            }
            if add.is_empty() {
                break;
            }
            for p in add {
                let poly = relaxed.remove(&p).unwrap();
                self.w.set_poly(p, poly);
                kept += 1;
            }
        }
        // 4. what still clashes goes (that wire, or that via with its wires), and its net is joined up again
        let moved_vias: std::collections::HashSet<u32> = map.iter().filter(|(o, v)| o != v).map(|(_, &v)| v).collect();
        let mut ripped = 0;
        for v in self.w.check_all() {
            if dbg {
                let names: Vec<String> = v.pieces.iter().map(|&p| self.w.b.nets[self.w.pieces[p as usize].net as usize].name.clone()).chain(v.objs.iter().map(|&o| format!("{:?}:{}", self.w.objs[o as usize].s.kind, self.w.objs[o as usize].s.name))).collect();
                eprintln!("  spring viol {} L{} at ({:.3},{:.3}) d {:.4} need {:.4}: {:?}", v.what, v.layer, v.at.x, v.at.y, v.dist, v.need, names);
            }
            let mut drop: Vec<u32> = v.pieces.iter().cloned().filter(|p| added.contains(p)).collect();
            if v.what.starts_with("via-") || v.what.starts_with("hole-") {
                if let Some(&o) = moved_vias.iter().find(|&&o| self.w.objs[o as usize].alive && self.w.objs[o as usize].s.center.dist(v.at) < 1e-6) {
                    drop.extend(self.w.obj_pieces[o as usize].iter().cloned().filter(|&p| self.w.pieces[p as usize].alive));
                    let net = self.w.objs[o as usize].s.net;
                    if !open.contains(&net) {
                        open.push(net);
                    }
                }
            }
            for p in drop {
                if self.w.pieces[p as usize].alive {
                    let net = self.w.pieces[p as usize].net;
                    self.rip_piece(p);
                    ripped += 1;
                    if !open.contains(&net) {
                        open.push(net);
                    }
                }
            }
        }
        self.prof_add("spring wb: thread, realise, check", t_thread);
        let t_pending = std::time::Instant::now();
        // the wires left short or not threaded: searched again between their own ends on their own layer
        let mut found = 0;
        for &(net, l, width, ta, tb) in &pending {
            if !self.w.objs[ta.obj as usize].alive || !self.w.objs[tb.obj as usize].alive {
                if !open.contains(&net) {
                    open.push(net);
                }
                continue;
            }
            let mut q = self.query(net, vec![ta.obj], vec![tb.obj], false);
            q.width = vec![width; self.w.nsig()];
            q.allowed = 1 << l;
            let ok = match self.search(&q) {
                Some(path) if path.vias.is_empty() => self.commit(net, &path, &q.width).is_ok(),
                _ => false,
            };
            if ok {
                found += 1;
            } else if !open.contains(&net) {
                open.push(net);
            }
        }
        if dbg {
            eprintln!("  spring: {} wires left short, {} did not thread, {} found again on their layer, {} kept as relaxed, {} ripped for clashes; {} nets to join up", skipped, unthreaded, found, kept, ripped, open.len());
        }
        self.prof_add("spring wb: re-search short wires", t_pending);
        let t_join = std::time::Instant::now();
        // 5. join up what is open: the nets above, and any other the relaxation left apart without a wire of its
        // own to rethread (a wire dropped where the triangles changed under it as its neighbours moved): a next try
        // holds it and the nets beside it
        let mut apart = Vec::new();
        for net in self.incomplete_nets() {
            if !was_open.contains(&net) && !open.contains(&net) {
                if dbg {
                    eprintln!("  spring: {} came apart without a wire to rethread", self.w.b.nets[net as usize].name);
                }
                open.push(net);
                apart.push(net);
            }
        }
        for &net in open.iter().filter(|n| !was_open.contains(n)) {
            for _ in 0..3 {
                if self.route_net(net) == 0 {
                    break;
                }
                for v in self.repair_net(net) {
                    self.route_net(v);
                }
            }
        }
        if self.total_deficit() > before_def {
            self.settle_except(8, &was_open);
        }
        self.cleanup();
        if dbg {
            for n in self.incomplete_nets() {
                eprintln!("  spring: incomplete {} (deficit {})", self.w.b.nets[n as usize].name, self.deficit(n));
            }
        }
        let rerouted = open.len();
        self.prof_add("spring wb: join up, settle, clean-up", t_join);
        let _t = std::time::Instant::now();
        let (def, viol, len) = (self.total_deficit(), self.w.check_all().len(), wire_len(&self.w));
        self.prof_add("spring wb: exact checks", _t);
        let moved = m.vias.iter().filter(|&&(o, n)| snap.0.objs[o as usize].s.center.dist(m.pos(n)) > 1e-3).count();
        if def > before_def || viol > before_viol || len >= before_len {
            if self.cfg.verbose {
                eprintln!("  spring: dropped (deficit {} -> {}, violations {} -> {}, wires {:.1} -> {:.1} mm, {} nets routed again)", before_def, def, before_viol, viol, before_len, len, rerouted);
            }
            // (the nets beside one that came apart, on its layers within 1 mm of its copper as it was)
            let mut hold_next = open.clone();
            for n in nets_beside(&snap.0, &apart, 1.0) {
                if !hold_next.contains(&n) {
                    hold_next.push(n);
                }
            }
            restore(self, snap);
            // left incomplete or unclean by the nets that had to be joined up again: the next try holds them
            return if (def > before_def || viol > before_viol) && !hold_next.is_empty() { Err(hold_next) } else { Ok(false) };
        }
        self.events.push(format!("spring relaxation: {} vias moved, wires {:.1} -> {:.1} mm ({} nets routed again where the relaxed wires did not fit)", moved, before_len, len, rerouted));
        if self.cfg.verbose {
            eprintln!("  spring: {} vias moved, wires {:.1} -> {:.1} mm, {} nets routed again", moved, before_len, len, rerouted);
        }
        Ok(true)
    }
}

/// Douglas-Peucker: drop points within `tol` of the simplified line (ends kept).
fn simplify(poly: &[P2], tol: f64) -> Vec<P2> {
    if poly.len() <= 2 {
        return poly.to_vec();
    }
    let mut keep = vec![false; poly.len()];
    keep[0] = true;
    keep[poly.len() - 1] = true;
    let mut stack = vec![(0usize, poly.len() - 1)];
    while let Some((a, b)) = stack.pop() {
        let mut best = (0.0, 0usize);
        for k in a + 1..b {
            let d = seg_point_dist(poly[a], poly[b], poly[k]);
            if d > best.0 {
                best = (d, k);
            }
        }
        if best.0 > tol {
            keep[best.1] = true;
            stack.push((a, best.1));
            stack.push((best.1, b));
        }
    }
    poly.iter().zip(keep).filter(|(_, k)| *k).map(|(p, _)| *p).collect()
}

impl Router {
    /// Where a movable signal via's wires want it: (length its wires would save there, that point), from the
    /// median of their first bends.
    fn via_pull(&self, o: u32) -> Option<(f64, P2)> {
        let ob = &self.w.objs[o as usize];
        if !ob.alive || ob.s.kind != Kind::Via || ob.plane || self.ties.iter().any(|t| t.0 == o || t.1 == o) {
            return None;
        }
        if self.w.b.partner.get(ob.s.net as usize).map_or(false, |&q| q != NONE) || self.w.b.dogbone.get(ob.s.net as usize).copied().unwrap_or(false) {
            return None;
        }
        let c = ob.s.center;
        let mut qs = Vec::new();
        for &p in &self.w.obj_pieces[o as usize] {
            let pc = &self.w.pieces[p as usize];
            if !pc.alive {
                continue;
            }
            if pc.frozen || pc.poly.len() < 2 {
                return None;
            }
            qs.push(if pc.a.obj == o { pc.poly[1] } else { pc.poly[pc.poly.len() - 2] });
        }
        if qs.len() < 2 {
            return None;
        }
        // Weiszfeld iteration for the geometric median
        let mut m = c;
        for _ in 0..60 {
            let (mut s, mut wsum) = (p2(0.0, 0.0), 0.0);
            for &q in &qs {
                let k = 1.0 / m.dist(q).max(1e-6);
                s = s + q * k;
                wsum += k;
            }
            m = s * (1.0 / wsum);
        }
        let now: f64 = qs.iter().map(|&q| c.dist(q)).sum();
        let best: f64 = qs.iter().map(|&q| m.dist(q)).sum();
        Some((now - best, m))
    }

    /// Stochastic via teleportation (experiment, TOPOR_TELEPORT): the relaxation never takes a via across a
    /// wire, so a via whose wires want it beyond one stays put. Such a via jumps instead, to where its wires
    /// want it or a random spot near there, and its wires are found again by the search on their own layers;
    /// a jump is kept when the board gets shorter. Returns (vias tried, jumps kept, mm saved).
    pub fn teleport_vias(&mut self, tries: usize, seed: u64) -> (usize, usize, f64) {
        let start = wire_len(&self.w);
        let mut rng = seed | 1;
        let mut rnd = move || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            (rng >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut order: Vec<(f64, u32)> = (0..self.w.objs.len() as u32).filter_map(|o| self.via_pull(o).map(|(d, _)| (d, o))).filter(|x| x.0 > 0.2).collect();
        // biggest pulls first, shuffled a little
        for x in order.iter_mut() {
            x.0 *= 0.75 + 0.5 * rnd();
        }
        order.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        let (mut tried, mut kept) = (0, 0);
        for (_, o) in order {
            let Some((d, m)) = self.via_pull(o) else { continue };
            if d <= 0.2 {
                continue;
            }
            tried += 1;
            let c = self.w.objs[o as usize].s.center;
            for k in 0..tries {
                // the median first, then random spots around it, spread growing with the attempts
                let x = if k == 0 {
                    m
                } else {
                    let s = 0.15 + 0.6 * k as f64 / tries as f64;
                    let (u, v) = (rnd().max(1e-12), rnd());
                    let g = (-2.0 * u.ln()).sqrt();
                    let along = c.lerp(m, 0.5 + 0.5 * rnd());
                    along + p2((std::f64::consts::TAU * v).cos(), (std::f64::consts::TAU * v).sin()) * (g * s)
                };
                if self.teleport(o, x) {
                    kept += 1;
                    break;
                }
            }
        }
        (tried, kept, start - wire_len(&self.w))
    }

    /// Via o jumps to x; each of its wires is searched again on its own layer from there. Kept (true) only
    /// when everything checks and the wires get shorter.
    fn teleport(&mut self, o: u32, x: P2) -> bool {
        let net = self.w.objs[o as usize].s.net;
        if !self.w.via_fits(net, x) {
            return false;
        }
        let before = wire_len(&self.w);
        let def0 = self.total_deficit();
        let drops = self.w.thread_drops;
        let old: Vec<(usize, f64, Term)> = self.w.obj_pieces[o as usize]
            .iter()
            .map(|&p| &self.w.pieces[p as usize])
            .filter(|pc| pc.alive)
            .map(|pc| (pc.layer as usize, pc.width, if pc.a.obj == o { pc.b } else { pc.a }))
            .collect();
        self.w.begin();
        let ok = (|| -> Option<()> {
            self.w.remove_via(o).ok()?;
            if !self.w.via_fits(net, x) {
                return None;
            }
            let nv = self.w.add_via(net, x).ok()?;
            for &(l, width, other) in &old {
                let mut q = self.query(net, vec![nv], vec![other.obj], false);
                q.width = vec![width; self.w.nsig()];
                q.allowed = 1 << l;
                let path = Search::new(&self.w, &q).run()?;
                if path.segs.len() != 1 || !path.vias.is_empty() {
                    return None;
                }
                let s = &path.segs[0];
                let (End::Obj(a), End::Obj(b)) = (s.a, s.b) else { return None };
                self.w.add_piece(net, l, width, a, b, &s.cx);
            }
            let mut viol = Vec::new();
            for p in self.w.realize_dirty() {
                if self.w.pieces[p as usize].overlap {
                    return None;
                }
                self.w.check_piece(p, &mut viol);
                if !viol.is_empty() {
                    return None;
                }
            }
            self.w.check_via_at(net, x, nv, &mut viol, true);
            if !viol.is_empty() {
                return None;
            }
            Some(())
        })();
        // a wire the re-triangulation could not re-thread is dropped: that takes the jump back too
        if ok.is_some() && self.w.thread_drops == drops && self.total_deficit() <= def0 && wire_len(&self.w) < before - 0.01 {
            self.w.commit();
            true
        } else {
            self.w.rollback();
            self.w.thread_drops = drops;
            false
        }
    }
}

/// Nets with copper within `d` of the wires of `nets` on the same layer (vias: on any layer), in board w.
fn nets_beside(w: &crate::world::World, nets: &[u32], d: f64) -> Vec<u32> {
    let mut out = Vec::new();
    if nets.is_empty() {
        return out;
    }
    let mine: Vec<(usize, P2, P2, f64)> = w.pieces.iter().filter(|p| p.alive && nets.contains(&p.net)).flat_map(|p| p.poly.windows(2).map(move |s| (p.layer as usize, s[0], s[1], p.width))).collect();
    let near = |l: Option<usize>, a: P2, b: P2, r: f64| mine.iter().any(|&(ml, c, e, mw)| l.map_or(true, |l| l == ml) && seg_seg_dist(a, b, c, e) < d + (r + mw) / 2.0);
    for p in w.pieces.iter().filter(|p| p.alive && !nets.contains(&p.net)) {
        if !out.contains(&p.net) && p.poly.windows(2).any(|s| near(Some(p.layer as usize), s[0], s[1], p.width)) {
            out.push(p.net);
        }
    }
    for o in w.objs.iter().filter(|o| o.alive && o.s.kind == crate::model::Kind::Via && !nets.contains(&o.s.net)) {
        if !out.contains(&o.s.net) && near(None, o.s.center, o.s.center, 2.0 * o.s.r) {
            out.push(o.s.net);
        }
    }
    out
}
