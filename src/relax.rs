//! Elastic relaxation of routed copper: every wire a chain of nodes under tension, every via a free node pulled
//! along each of its wires, every pair of copper elements held apart by the least distance its rules require
//! (position-based, Jacobi), so the whole layout pulls tight at once while wires and vias push each other aside.
//!
//! The caller owns the model: nodes, the convex elements built on them (a via is a 1-node element, a wire
//! stretch a 2-node capsule core, a pad a polygon over fixed nodes) and, per element, the other elements near it
//! with the least core distance between them (radii and clearances folded in). It rebuilds those lists every
//! batch of iterations, with a skin wide enough for the motion a batch allows (`Params::dmax`).
//! Multi-threaded CPU, or Metal on macOS (same arithmetic, f32).

use rayon::prelude::*;

#[derive(Default)]
pub struct Relax {
    /// Node positions (mm), updated in place.
    pub pos: Vec<[f32; 2]>,
    /// Per node: 0 fixed, 1 wire node (tension along its chain), 2 via (pulled along each of its wires).
    pub kind: Vec<u32>,
    /// Tension neighbours per node (CSR): a wire node's two chain neighbours, a via's first node on each wire.
    pub tn_off: Vec<u32>,
    pub tn: Vec<u32>,
    /// Per node (cx, cy, r, 0): stays within r of (cx, cy); r < 0 = free (a plane via keeps to its pour).
    pub leash: Vec<[f32; 4]>,
    /// Elements: convex cores over nodes (CSR; 1 node = disc, 2 = capsule, more = convex polygon, fixed only).
    pub el_off: Vec<u32>,
    pub el_node: Vec<u32>,
    /// Movable elements of each node (CSR): the node takes its share of their constraints.
    pub ne_off: Vec<u32>,
    pub ne: Vec<u32>,
    /// Constraints per element (CSR): the other element and the least core distance to keep from it.
    pub c_off: Vec<u32>,
    pub c_el: Vec<u32>,
    pub c_need: Vec<f32>,
    /// Bumped by the caller whenever the constraint lists change; everything else is fixed after the
    /// first run (the GPU keeps it resident).
    pub c_version: u64,
    /// Per node: step scale (halved while the node's gaps stay short after projection, regrown once they
    /// hold, so tension never drives a node deep into a wedge it cannot leave) and its shortest gap deficit.
    pub scale: Vec<f32>,
    pub viol: Vec<f32>,
    #[cfg(all(target_os = "macos", feature = "metal"))]
    pub(crate) gpu: Option<crate::relax_gpu::State>,
}

#[derive(Clone, Copy)]
#[repr(C)]
pub struct Params {
    /// Iterations to run.
    pub iters: u32,
    /// Constraint projections after each tension step.
    pub sub: u32,
    /// Tension step of wire nodes and of vias (mm per unit pull), the tangential spacing gain, and the most a
    /// node moves in one tension step (mm).
    pub alpha_wire: f32,
    pub alpha_via: f32,
    pub beta: f32,
    pub dmax: f32,
    /// Gap deficit above which a node's step is halved (mm).
    pub tol: f32,
}

type V = [f32; 2];

#[inline]
fn sub(a: V, b: V) -> V {
    [a[0] - b[0], a[1] - b[1]]
}
#[inline]
fn add(a: V, b: V) -> V {
    [a[0] + b[0], a[1] + b[1]]
}
#[inline]
fn mul(a: V, k: f32) -> V {
    [a[0] * k, a[1] * k]
}
#[inline]
fn dot(a: V, b: V) -> f32 {
    a[0] * b[0] + a[1] * b[1]
}
#[inline]
fn unit(a: V) -> V {
    let l = dot(a, a).sqrt();
    if l > 1e-9 { mul(a, 1.0 / l) } else { [0.0, 0.0] }
}

/// Closest points of segments p1-q1 and p2-q2 (either may be a point): parameters (s, t) on each.
#[inline]
fn seg_seg(p1: V, q1: V, p2: V, q2: V) -> (f32, f32) {
    let (d1, d2, r) = (sub(q1, p1), sub(q2, p2), sub(p1, p2));
    let (a, e, f) = (dot(d1, d1), dot(d2, d2), dot(d2, r));
    const EPS: f32 = 1e-12;
    if a <= EPS && e <= EPS {
        return (0.0, 0.0);
    }
    if a <= EPS {
        return (0.0, (f / e).clamp(0.0, 1.0));
    }
    let c = dot(d1, r);
    if e <= EPS {
        return ((-c / a).clamp(0.0, 1.0), 0.0);
    }
    let b = dot(d1, d2);
    let den = a * e - b * b;
    let mut s = if den > EPS { ((b * f - c * e) / den).clamp(0.0, 1.0) } else { 0.0 };
    let mut t = (b * s + f) / e;
    if t < 0.0 {
        t = 0.0;
        s = (-c / a).clamp(0.0, 1.0);
    } else if t > 1.0 {
        t = 1.0;
        s = ((b - c) / a).clamp(0.0, 1.0);
    }
    (s, t)
}

impl Relax {
    #[inline]
    fn minv(&self, n: u32) -> f32 {
        if self.kind[n as usize] != 0 { 1.0 } else { 0.0 }
    }

    /// Node i after one tension step from `p`, and its new step scale.
    fn tension(&self, i: usize, p: &[V], prm: &Params) -> (V, f32) {
        let x = p[i];
        let sc = if self.viol[i] > prm.tol { self.scale[i] * 0.5 } else { (self.scale[i] * 1.25).min(1.0) };
        let (a, b) = (self.tn_off[i] as usize, self.tn_off[i + 1] as usize);
        let d = match self.kind[i] {
            1 if b - a == 2 => {
                let (u, w) = (p[self.tn[a] as usize], p[self.tn[a + 1] as usize]);
                let g = add(unit(sub(u, x)), unit(sub(w, x)));
                let t = unit(sub(w, u));
                let lap = sub(mul(add(u, w), 0.5), x);
                add(mul(g, prm.alpha_wire), mul(t, prm.beta * dot(lap, t)))
            }
            2 => {
                let mut g = [0.0, 0.0];
                for k in a..b {
                    g = add(g, unit(sub(p[self.tn[k] as usize], x)));
                }
                mul(g, prm.alpha_via)
            }
            _ => return (x, sc),
        };
        let d = mul(d, sc);
        let l = dot(d, d).sqrt();
        (add(x, if l > prm.dmax { mul(d, prm.dmax / l) } else { d }), sc)
    }

    fn is_fixed(&self, f: usize) -> bool {
        (self.el_off[f]..self.el_off[f + 1]).all(|k| self.kind[self.el_node[k as usize] as usize] == 0)
    }

    /// The correction for node i (at `xi`, the rest at `p`) from the constraint keeping element e `need`
    /// from element f, and the gap deficit; None when the gap holds.
    fn contact(&self, i: usize, xi: V, e: usize, f: usize, need: f32, p: &[V]) -> Option<(V, f32)> {
        let at = |n: u32| if n as usize == i { xi } else { p[n as usize] };
        let (e0, e1) = (self.el_off[e] as usize, self.el_off[e + 1] as usize);
        let (a0, a1) = (self.el_node[e0], self.el_node[e1 - 1]);
        let (pa, qa) = (at(a0), at(a1));
        let (f0, f1) = (self.el_off[f] as usize, self.el_off[f + 1] as usize);
        // closest pair over the other element's edges (a polygon is fixed: no share for it)
        let mut best = (f32::INFINITY, 0.0f32, 0.0f32, [0.0f32, 0.0], [0.0f32, 0.0], 0u32, 0u32);
        let nf = f1 - f0;
        let edges = if nf <= 2 { 1 } else { nf };
        for j in 0..edges {
            let (b0, b1) = if nf <= 2 { (self.el_node[f0], self.el_node[f1 - 1]) } else { (self.el_node[f0 + j], self.el_node[f0 + (j + 1) % nf]) };
            let (pb, qb) = (at(b0), at(b1));
            let (s, t) = seg_seg(pa, qa, pb, qb);
            let ca = add(pa, mul(sub(qa, pa), s));
            let cb = add(pb, mul(sub(qb, pb), t));
            let dv = sub(ca, cb);
            let d2 = dot(dv, dv);
            if d2 < best.0 {
                best = (d2, s, t, ca, cb, b0, b1);
            }
        }
        let (d2, s, t, ca, cb, b0, b1) = best;
        if d2 >= need * need {
            return None;
        }
        let d = d2.sqrt();
        if d < 1e-7 {
            return None;
        }
        let nrm = mul(sub(ca, cb), 1.0 / d);
        let la = if a0 == a1 { 1.0 } else if i as u32 == a0 { 1.0 - s } else { s };
        let mut w = if a0 == a1 { self.minv(a0) } else { self.minv(a0) * (1.0 - s) * (1.0 - s) + self.minv(a1) * s * s };
        if nf <= 2 {
            w += if b0 == b1 { self.minv(b0) } else { self.minv(b0) * (1.0 - t) * (1.0 - t) + self.minv(b1) * t * t };
        }
        if w < 1e-9 {
            return None;
        }
        // the exact share for a contact near a pinned end is a long lever (penetration / s): capped, so a
        // node never leaps
        Some((mul(nrm, (need - d) * (la / w).min(2.0)), need - d))
    }

    /// Node i after one projection of its elements' constraints, from `p`, and its largest gap deficit. Moving
    /// neighbours share a Jacobi average; fixed obstacles then apply in full from where that leaves the node,
    /// so no number of pushing neighbours can carry a node through a boundary, one after another, so a run of
    /// short boundary segments (a densified curve) pushes once, not once per segment.
    fn project(&self, i: usize, p: &[V]) -> (V, f32) {
        let x = p[i];
        if self.kind[i] == 0 {
            return (x, 0.0);
        }
        let cons = || (self.ne_off[i]..self.ne_off[i + 1]).flat_map(|k| {
            let e = self.ne[k as usize] as usize;
            (self.c_off[e]..self.c_off[e + 1]).map(move |c| (e, self.c_el[c as usize] as usize, self.c_need[c as usize]))
        });
        let (mut acc, mut n, mut worst) = ([0.0f32, 0.0], 0u32, 0.0f32);
        for (e, f, need) in cons() {
            if self.is_fixed(f) {
                continue;
            }
            if let Some((d, pen)) = self.contact(i, x, e, f, need, p) {
                acc = add(acc, d);
                n += 1;
                worst = worst.max(pen);
            }
        }
        let mut y = if n > 0 { add(x, mul(acc, 1.0 / n as f32)) } else { x };
        for (e, f, need) in cons() {
            if !self.is_fixed(f) {
                continue;
            }
            if let Some((d, pen)) = self.contact(i, y, e, f, need, p) {
                y = add(y, d);
                worst = worst.max(pen);
            }
        }
        let lz = self.leash[i];
        if lz[2] >= 0.0 {
            let o = sub(y, [lz[0], lz[1]]);
            let l = dot(o, o).sqrt();
            if l > lz[2] {
                y = add([lz[0], lz[1]], mul(o, lz[2] / l));
            }
        }
        (y, worst)
    }

    /// Run `prm.iters` iterations on the CPU.
    pub fn run_cpu(&mut self, prm: &Params) {
        let n = self.pos.len();
        self.init_scale();
        let mut a = std::mem::take(&mut self.pos);
        let mut b = vec![[0.0f32; 2]; n];
        for _ in 0..prm.iters {
            let ts: Vec<(V, f32)> = (0..n).into_par_iter().map(|i| self.tension(i, &a, prm)).collect();
            for (i, (y, sc)) in ts.into_iter().enumerate() {
                b[i] = y;
                self.scale[i] = sc;
            }
            std::mem::swap(&mut a, &mut b);
            for _ in 0..prm.sub {
                let ps: Vec<(V, f32)> = (0..n).into_par_iter().map(|i| self.project(i, &a)).collect();
                for (i, (y, v)) in ps.into_iter().enumerate() {
                    b[i] = y;
                    self.viol[i] = v;
                }
                std::mem::swap(&mut a, &mut b);
            }
        }
        self.pos = a;
    }

    pub(crate) fn init_scale(&mut self) {
        if self.scale.len() != self.pos.len() {
            self.scale = vec![1.0; self.pos.len()];
            self.viol = vec![0.0; self.pos.len()];
        }
    }

    /// Run on Metal when available (macOS), else on the CPU. Returns whether the GPU ran it.
    pub fn run(&mut self, prm: &Params) -> bool {
        #[cfg(all(target_os = "macos", feature = "metal"))]
        if std::env::var_os("GR_RELAX_CPU").is_none() && crate::relax_gpu::run(self, prm) {
            return true;
        }
        self.run_cpu(prm);
        false
    }
}
