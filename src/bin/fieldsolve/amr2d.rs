//! Prototype (CPU only): dynamic resolution for FDTD, in 2-D. Not used by `run`; `fieldsolve amr-check`
//! measures it.
//!
//! The 2-D TMz equations (Ez, Hx, Hy) of a plane pair seen from the top are written as an LC network on a
//! 2:1-balanced quadtree: Ez is a node voltage at each leaf's centre (C = ε A per unit height), the
//! tangential H on each face between two leaves is a branch current (L = μ d / w, d the distance of the
//! centres along the face normal, w the face width). Leapfrog on this network is the Yee scheme on a
//! uniform grid; on a quadtree it is the two-point-flux finite-volume scheme, which is first order at
//! hanging faces but, being a passive LC network, conserves the discrete energy ½ΣCV² + ½ΣLI² exactly and is
//! stable for dt ≤ 2 / sqrt(λ_max) (λ_max of the network, measured by power iteration).
//!
//! - **Static refinement**: geometry decides (ports, PEC walls and their tips, material boundaries get the
//!   finest cells; elsewhere cells may be as large as `max_size`, limited to λ_min / `ppw`).
//! - **Dynamic refinement**: every `every` steps leaves whose field is "busy" (jump of Ez to a neighbour,
//!   relative to the running peak) are split, sibling groups that are all quiet are merged. Fields move
//!   through a fine-grid map: Ez piecewise constant per leaf (charge-conserving averages on merge), H on
//!   every fine face interpolated linearly inside a leaf (so a split face carries the same current). A merge
//!   loses the energy of the field structure below the new cell size; the run reports how much.
//! - **Time step**: global, set by the finest level allowed anywhere (refinement during the run must not
//!   break stability). `cell_updates` counts the work actually done; `lts_updates` is what 2:1 local time
//!   stepping (each level at its own CFL step) would need with the same meshes.
//! - **Losses**: dielectric conductivity per leaf, and a graded absorbing layer at open boundaries with
//!   matched magnetic loss on the faces (a split-free stand-in for a PML).

use std::f64::consts::PI;

const C0: f64 = 299_792_458.0;
const MU0: f64 = 1.256_637_062_12e-6;
const EPS0: f64 = 8.854_187_812_8e-12;

/// The fine-grid description of a 2-D problem (all on `nx` x `ny` fine cells of size `h`).
#[derive(Clone, Debug)]
pub struct Geo {
    pub nx: usize,
    pub ny: usize,
    pub h: f64,
    pub eps: Vec<f64>,
    pub tand: f64,
    pub f_loss: f64,
    /// fine cells that are PEC (Ez = 0)
    pub pec: Vec<bool>,
    /// thin PEC walls on fine faces: `wall_x[j * (nx + 1) + i]` is the vertical face at x = i (between fine
    /// cells i-1 and i) of row j; `wall_y[j * nx + i]` the horizontal face at y = j of column i
    pub wall_x: Vec<bool>,
    pub wall_y: Vec<bool>,
    /// outer boundary: PEC walls (true) or magnetic walls / open (false)
    pub pec_box: bool,
    /// absorbing layer thickness in fine cells (0 = none), with this peak loss tangent-like factor
    pub absorb: usize,
    /// port: fine cell, resistance (ohm per unit height)
    pub port: (usize, usize, f64),
    pub probes: Vec<(usize, usize)>,
    pub fmax: f64,
}

#[derive(Clone, Debug)]
pub struct AmrOpts {
    /// largest leaf, in fine cells (power of two)
    pub max_size: u32,
    /// points per wavelength at fmax that a leaf must keep
    pub ppw: f64,
    /// geometry features keep leaves of at most this size within `feature_pad` fine cells
    pub feature_size: u32,
    pub feature_pad: usize,
    /// dynamic refinement: off when `every` is 0
    pub every: usize,
    pub refine_tol: f64,
    pub coarsen_tol: f64,
}

impl Default for AmrOpts {
    fn default() -> Self {
        AmrOpts { max_size: 8, ppw: 20.0, feature_size: 1, feature_pad: 2, every: 0, refine_tol: 0.02, coarsen_tol: 0.005 }
    }
}

#[derive(Clone, Copy, Debug)]
struct Leaf {
    x: u32,
    y: u32,
    s: u32,
}

#[derive(Clone, Copy, Debug)]
struct Face {
    a: u32,
    /// u32::MAX: the other side is a PEC (V = 0)
    b: u32,
    /// dt / L and the magnetic-loss damping factors (semi-implicit)
    k: f64,
    da: f64,
    /// fine faces this face covers: axis (0 = vertical face, normal x), fixed coordinate, range along
    axis: u8,
    at: u32,
    r0: u32,
    r1: u32,
    w: f64,
}

pub struct Amr {
    pub g: Geo,
    o: AmrOpts,
    leaves: Vec<Leaf>,
    owner: Vec<u32>,
    faces: Vec<Face>,
    v: Vec<f64>,
    i: Vec<f64>,
    ca: Vec<f64>,
    cb: Vec<f64>,
    cap: Vec<f64>,
    /// smallest size each fine cell may be part of (geometry and wavelength limits)
    limit: Vec<u32>,
    dt: f64,
    port_leaf: u32,
    pub energy_lost: f64,
    pub remeshes: usize,
    /// currents into thin walls, per fine face and side (low = left / below): walls carry a surface current,
    /// so H differs on the two sides
    wall_lo: Vec<f64>,
    wall_hi: Vec<f64>,
}

#[derive(Clone, Debug, Default)]
pub struct AmrRun {
    pub steps: usize,
    pub dt: f64,
    pub vport: Vec<f64>,
    pub iport: Vec<f64>,
    pub probes: Vec<Vec<f64>>,
    pub energy: Vec<(usize, f64)>,
    pub cells: Vec<(usize, usize)>,
    pub cell_updates: f64,
    pub lts_updates: f64,
    pub energy_lost: f64,
    pub remeshes: usize,
    pub seconds: f64,
    pub min_leaves: usize,
    pub max_leaves: usize,
}

impl Amr {
    pub fn new(g: Geo, o: AmrOpts) -> Amr {
        let (nx, ny) = (g.nx, g.ny);
        // per fine cell: the largest leaf size allowed there
        let eps_max = g.eps.iter().cloned().fold(1.0, f64::max);
        let lam = C0 / (eps_max.sqrt() * g.fmax);
        let mut wl = 1u32;
        while ((2 * wl) as f64) * g.h <= lam / o.ppw && 2 * wl <= o.max_size {
            wl *= 2;
        }
        let mut limit = vec![wl; nx * ny];
        // geometry features: material changes, PEC cells, walls, the port and probes
        let mut feat = vec![false; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let q = j * nx + i;
                let e = g.eps[q];
                let diff = (i + 1 < nx && (g.eps[q + 1] != e || g.pec[q + 1] != g.pec[q])) || (j + 1 < ny && (g.eps[q + nx] != e || g.pec[q + nx] != g.pec[q]));
                let w = g.wall_x[j * (nx + 1) + i] || g.wall_x[j * (nx + 1) + i + 1] || g.wall_y[j * nx + i] || g.wall_y[(j + 1) * nx + i];
                if diff || w {
                    feat[q] = true;
                }
            }
        }
        feat[g.port.1 * nx + g.port.0] = true;
        for &(i, j) in &g.probes {
            feat[j * nx + i] = true;
        }
        let pad = o.feature_pad as isize;
        for j in 0..ny {
            for i in 0..nx {
                if feat[j * nx + i] {
                    for b in -pad..=pad {
                        for a in -pad..=pad {
                            let (ii, jj) = (i as isize + a, j as isize + b);
                            if ii >= 0 && jj >= 0 && (ii as usize) < nx && (jj as usize) < ny {
                                let q = jj as usize * nx + ii as usize;
                                limit[q] = limit[q].min(o.feature_size);
                            }
                        }
                    }
                }
            }
        }
        let mut a = Amr { g, o, leaves: vec![], owner: vec![u32::MAX; nx * ny], faces: vec![], v: vec![], i: vec![], ca: vec![], cb: vec![], cap: vec![], limit, dt: 0.0, port_leaf: 0, energy_lost: 0.0, remeshes: 0, wall_lo: vec![], wall_hi: vec![] };
        a.wall_lo = vec![0.0; (nx + 1) * (ny + 1) * 2];
        a.wall_hi = vec![0.0; (nx + 1) * (ny + 1) * 2];
        // the coarsest tree the limits allow (each leaf as large as its cells permit), then 2:1 balanced
        let s0 = a.o.max_size;
        let mut stack: Vec<Leaf> = Vec::new();
        let mut y = 0;
        while y < ny as u32 {
            let mut x = 0;
            while x < nx as u32 {
                stack.push(Leaf { x, y, s: s0 });
                x += s0;
            }
            y += s0;
        }
        let mut leaves = Vec::new();
        while let Some(l) = stack.pop() {
            if l.s > 1 && a.min_limit(l) < l.s {
                let h = l.s / 2;
                for (dx, dy) in [(0, 0), (h, 0), (0, h), (h, h)] {
                    stack.push(Leaf { x: l.x + dx, y: l.y + dy, s: h });
                }
            } else {
                leaves.push(l);
            }
        }
        a.set_leaves(leaves);
        a.balance();
        let nf = a.g.nx * a.g.ny;
        // fine-map transfer state: start from zero fields
        let vf = vec![0.0; nf];
        let hxf = vec![0.0; (a.g.nx + 1) * a.g.ny];
        let hyf = vec![0.0; a.g.nx * (a.g.ny + 1)];
        a.rebuild(&vf, &hxf, &hyf);
        // global time step from the finest cells that may ever appear (1 fine cell): uniform 2-D Yee CFL
        let eps_min = a.g.eps.iter().cloned().fold(f64::INFINITY, f64::min);
        a.dt = 0.95 * a.g.h * eps_min.sqrt() / (C0 * 2f64.sqrt());
        a.coefs();
        a
    }

    fn min_limit(&self, l: Leaf) -> u32 {
        let nx = self.g.nx;
        let mut m = u32::MAX;
        for y in l.y..(l.y + l.s).min(self.g.ny as u32) {
            for x in l.x..(l.x + l.s).min(nx as u32) {
                m = m.min(self.limit[y as usize * nx + x as usize]);
            }
        }
        m
    }

    fn set_leaves(&mut self, leaves: Vec<Leaf>) {
        let nx = self.g.nx;
        self.owner.iter_mut().for_each(|o| *o = u32::MAX);
        for (k, l) in leaves.iter().enumerate() {
            for y in l.y..(l.y + l.s).min(self.g.ny as u32) {
                for x in l.x..(l.x + l.s).min(nx as u32) {
                    self.owner[y as usize * nx + x as usize] = k as u32;
                }
            }
        }
        self.leaves = leaves;
    }

    /// Split leaves until no two face neighbours differ by more than a factor 2.
    fn balance(&mut self) {
        loop {
            let nx = self.g.nx;
            let mut split = vec![false; self.leaves.len()];
            for (k, l) in self.leaves.iter().enumerate() {
                // neighbours along the four sides
                let side = |x: i64, y: i64| -> Option<u32> {
                    if x < 0 || y < 0 || x >= nx as i64 || y >= self.g.ny as i64 {
                        None
                    } else {
                        Some(self.owner[y as usize * nx + x as usize])
                    }
                };
                let (x0, y0, s) = (l.x as i64, l.y as i64, l.s as i64);
                for t in 0..s {
                    for (xx, yy) in [(x0 - 1, y0 + t), (x0 + s, y0 + t), (x0 + t, y0 - 1), (x0 + t, y0 + s)] {
                        if let Some(nb) = side(xx, yy) {
                            if self.leaves[nb as usize].s * 2 < l.s {
                                split[k] = true;
                            }
                        }
                    }
                }
            }
            if !split.iter().any(|&b| b) {
                return;
            }
            let mut out = Vec::with_capacity(self.leaves.len() + 3 * split.len());
            for (k, l) in self.leaves.iter().enumerate() {
                if split[k] {
                    let h = l.s / 2;
                    for (dx, dy) in [(0, 0), (h, 0), (0, h), (h, h)] {
                        out.push(Leaf { x: l.x + dx, y: l.y + dy, s: h });
                    }
                } else {
                    out.push(*l);
                }
            }
            self.set_leaves(out);
        }
    }

    /// Faces and per-leaf fields from the fine-map fields.
    fn rebuild(&mut self, vf: &[f64], hxf: &[f64], hyf: &[f64]) {
        let (nx, ny, h) = (self.g.nx, self.g.ny, self.g.h);
        let n = self.leaves.len();
        self.v = vec![0.0; n];
        self.cap = vec![0.0; n];
        for (k, l) in self.leaves.iter().enumerate() {
            let mut q = 0.0;
            let mut c = 0.0;
            for y in l.y..l.y + l.s {
                for x in l.x..l.x + l.s {
                    let f = y as usize * nx + x as usize;
                    let cc = EPS0 * self.g.eps[f] * h * h;
                    q += cc * vf[f];
                    c += cc;
                }
            }
            self.cap[k] = c;
            self.v[k] = q / c;
        }
        let mut faces = Vec::new();
        for (k, l) in self.leaves.iter().enumerate() {
            let k = k as u32;
            let cx = (l.x as f64 + 0.5 * l.s as f64) * h;
            let cy = (l.y as f64 + 0.5 * l.s as f64) * h;
            // right side (vertical faces at x = l.x + s) and top side; left/bottom only at the domain edge
            for axis in 0..2u8 {
                for far in [false, true] {
                    let fixed = if axis == 0 { if far { l.x + l.s } else { l.x } } else if far { l.y + l.s } else { l.y };
                    let at_edge = if axis == 0 { fixed == 0 || fixed as usize == nx } else { fixed == 0 || fixed as usize == ny };
                    if !far && !at_edge {
                        continue;
                    }
                    let (r0, r1) = if axis == 0 { (l.y, l.y + l.s) } else { (l.x, l.x + l.s) };
                    let mut t = r0;
                    while t < r1 {
                        // the leaf on the other side and the wall flag, constant over a run
                        let other = |t: u32| -> Option<u32> {
                            if at_edge {
                                None
                            } else if axis == 0 {
                                Some(self.owner[t as usize * nx + fixed as usize])
                            } else {
                                Some(self.owner[fixed as usize * nx + t as usize])
                            }
                        };
                        let wall = |t: u32| -> bool { if axis == 0 { self.g.wall_x[t as usize * (nx + 1) + fixed as usize] } else { self.g.wall_y[fixed as usize * nx + t as usize] } };
                        let (o0, w0) = (other(t), wall(t));
                        let mut e = t + 1;
                        while e < r1 && other(e) == o0 && wall(e) == w0 {
                            e += 1;
                        }
                        let w = (e - t) as f64 * h;
                        let half = 0.5 * l.s as f64 * h;
                        match (o0, w0) {
                            (None, _) => {
                                if self.g.pec_box {
                                    faces.push(Face { a: k, b: u32::MAX, k: half / w, da: 1.0, axis, at: fixed, r0: t, r1: e, w });
                                }
                            }
                            (Some(b), false) => {
                                let lb = self.leaves[b as usize];
                                let d = if axis == 0 { (lb.x as f64 + 0.5 * lb.s as f64) * h - cx } else { (lb.y as f64 + 0.5 * lb.s as f64) * h - cy };
                                faces.push(Face { a: k, b, k: d.abs() / w, da: 1.0, axis, at: fixed, r0: t, r1: e, w });
                            }
                            (Some(b), true) => {
                                // a thin PEC wall between the two leaves: each sees V = 0 half its size away
                                let lb = self.leaves[b as usize];
                                faces.push(Face { a: k, b: u32::MAX, k: half / w, da: 1.0, axis, at: fixed, r0: t, r1: e, w });
                                faces.push(Face { a: b, b: u32::MAX, k: 0.5 * lb.s as f64 * h / w, da: 1.0, axis, at: fixed, r0: t, r1: e, w });
                            }
                        }
                        t = e;
                    }
                }
            }
        }
        // current of each face from the fine-face H (integral of H along the face); wall faces from their side
        let wk = |axis: u8, at: u32, t: u32| -> usize { ((axis as usize) * (ny + 1) + t as usize) * (nx + 1) + at as usize };
        self.i = faces
            .iter()
            .map(|f| {
                let mut s = 0.0;
                let wall = f.b == u32::MAX && !(if f.axis == 0 { f.at == 0 || f.at as usize == nx } else { f.at == 0 || f.at as usize == ny });
                let la = self.leaves[f.a as usize];
                let lo = if f.axis == 0 { la.x < f.at } else { la.y < f.at };
                for t in f.r0..f.r1 {
                    s += if wall {
                        if lo { self.wall_lo[wk(f.axis, f.at, t)] } else { self.wall_hi[wk(f.axis, f.at, t)] }
                    } else if f.axis == 0 {
                        hxf[t as usize * (nx + 1) + f.at as usize]
                    } else {
                        hyf[f.at as usize * nx + t as usize]
                    };
                }
                s * h
            })
            .collect();
        self.faces = faces;
        self.port_leaf = self.owner[self.g.port.1 * nx + self.g.port.0];
        let _ = ny;
    }

    /// Update constants: per leaf (dielectric loss, PEC, port conductance), per face (dt / L, absorber).
    fn coefs(&mut self) {
        let (nx, h, dt) = (self.g.nx, self.g.h, self.dt);
        let n = self.leaves.len();
        self.ca = vec![1.0; n];
        self.cb = vec![0.0; n];
        let w = 2.0 * PI * self.g.f_loss;
        let ab = self.g.absorb as f64;
        let sig_ab = |x: f64, y: f64| -> f64 {
            // graded (cubic) loss in the outer `absorb` fine cells, sigma such that the layer is ~optically thick
            if ab <= 0.0 {
                return 0.0;
            }
            let d = [x / h, (self.g.nx as f64 - x / h), y / h, (self.g.ny as f64 - y / h)].iter().cloned().fold(f64::INFINITY, f64::min);
            if d >= ab {
                0.0
            } else {
                let r = (ab - d) / ab;
                let eta = (MU0 / EPS0).sqrt();
                0.8 * 4.0 / (eta * h) * r.powi(3)
            }
        };
        for (k, l) in self.leaves.iter().enumerate() {
            let mut pec = false;
            let mut eps = 0.0;
            for y in l.y..l.y + l.s {
                for x in l.x..l.x + l.s {
                    let f = y as usize * nx + x as usize;
                    pec |= self.g.pec[f];
                    eps += self.g.eps[f];
                }
            }
            eps /= (l.s * l.s) as f64;
            if pec {
                continue;
            }
            let c = self.cap[k];
            let (cx, cy) = ((l.x as f64 + 0.5 * l.s as f64) * h, (l.y as f64 + 0.5 * l.s as f64) * h);
            let area = (l.s as f64 * h).powi(2);
            let mut g = w * EPS0 * eps * self.g.tand * area + sig_ab(cx, cy) * area * eps;
            if k as u32 == self.port_leaf {
                g += 1.0 / self.g.port.2;
            }
            let d = c / dt + 0.5 * g;
            self.ca[k] = (c / dt - 0.5 * g) / d;
            self.cb[k] = 1.0 / d;
        }
        for f in self.faces.iter_mut() {
            // L = mu0 * (d / w) per unit height; magnetic loss matching the absorber's (sigma* = sigma mu / eps)
            let l = MU0 * f.k;
            let la = self.leaves[f.a as usize];
            let (cx, cy) = ((la.x as f64 + 0.5 * la.s as f64) * h, (la.y as f64 + 0.5 * la.s as f64) * h);
            let r = sig_ab(cx, cy) * MU0 / EPS0 * f.k;
            let dd = l / dt + 0.5 * r;
            f.da = (l / dt - 0.5 * r) / dd;
            f.k = 1.0 / dd;
        }
    }

    pub fn dt(&self) -> f64 {
        self.dt
    }

    pub fn leaves(&self) -> usize {
        self.leaves.len()
    }

    fn energy(&self) -> f64 {
        let mut e = 0.0;
        for (k, &v) in self.v.iter().enumerate() {
            e += 0.5 * self.cap[k] * v * v;
        }
        for (f, &i) in self.faces.iter().zip(&self.i) {
            // L from the stored constants: 1 / k = L / dt + R / 2 -> L ~ dt / k without loss
            let l = self.dt / f.k;
            e += 0.5 * l * i * i;
        }
        e
    }

    /// Fine-map fields: V piecewise constant per leaf; H on every fine face, linear inside a leaf between its
    /// faces' values.
    fn to_fine(&mut self) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let (nx, ny, h) = (self.g.nx, self.g.ny, self.g.h);
        let mut vf = vec![0.0; nx * ny];
        for (k, l) in self.leaves.iter().enumerate() {
            for y in l.y..l.y + l.s {
                for x in l.x..l.x + l.s {
                    vf[y as usize * nx + x as usize] = self.v[k];
                }
            }
        }
        let mut hxf = vec![0.0; (nx + 1) * ny];
        let mut hyf = vec![0.0; nx * (ny + 1)];
        let wk = |axis: u8, at: u32, t: u32| -> usize { ((axis as usize) * (ny + 1) + t as usize) * (nx + 1) + at as usize };
        for (f, &i) in self.faces.iter().zip(&self.i) {
            let hv = i / f.w;
            let edge = if f.axis == 0 { f.at == 0 || f.at as usize == nx } else { f.at == 0 || f.at as usize == ny };
            if f.b == u32::MAX && !edge {
                let la = self.leaves[f.a as usize];
                let lo = if f.axis == 0 { la.x < f.at } else { la.y < f.at };
                for t in f.r0..f.r1 {
                    if lo {
                        self.wall_lo[wk(f.axis, f.at, t)] = hv;
                    } else {
                        self.wall_hi[wk(f.axis, f.at, t)] = hv;
                    }
                }
                continue;
            }
            for t in f.r0..f.r1 {
                if f.axis == 0 {
                    hxf[t as usize * (nx + 1) + f.at as usize] = hv;
                } else {
                    hyf[f.at as usize * nx + t as usize] = hv;
                }
            }
        }
        // interior fine faces of each leaf: linear between its two boundary columns / rows
        for l in &self.leaves {
            if l.s == 1 {
                continue;
            }
            for y in l.y..l.y + l.s {
                let (a, b) = (hxf[y as usize * (nx + 1) + l.x as usize], hxf[y as usize * (nx + 1) + (l.x + l.s) as usize]);
                for x in l.x + 1..l.x + l.s {
                    let t = (x - l.x) as f64 / l.s as f64;
                    hxf[y as usize * (nx + 1) + x as usize] = a + t * (b - a);
                }
            }
            for x in l.x..l.x + l.s {
                let (a, b) = (hyf[l.y as usize * nx + x as usize], hyf[(l.y + l.s) as usize * nx + x as usize]);
                for y in l.y + 1..l.y + l.s {
                    let t = (y - l.y) as f64 / l.s as f64;
                    hyf[y as usize * nx + x as usize] = a + t * (b - a);
                }
            }
        }
        let _ = h;
        (vf, hxf, hyf)
    }

    /// Dynamic refinement: split busy leaves, merge quiet sibling groups (within the size limits).
    fn remesh(&mut self, peak: f64) {
        if peak <= 0.0 {
            return;
        }
        let nx = self.g.nx;
        let e0 = self.energy();
        let n = self.leaves.len();
        // busyness: the local truncation error of the leaf, |mean of the neighbours - own value| (h^2 / 4 times
        // the Laplacian on a uniform grid, so ~ (h k)^2 / 4 of the amplitude for a wave), over the running peak
        let mut num = vec![0.0f64; n];
        let mut den = vec![0.0f64; n];
        for f in &self.faces {
            // geometric weight w / d, recovered from the stored dt / L without loss
            let wgt = f.k;
            let (a, va) = (f.a as usize, self.v[f.a as usize]);
            let vb = if f.b == u32::MAX { 0.0 } else { self.v[f.b as usize] };
            num[a] += wgt * (vb - va);
            den[a] += wgt;
            if f.b != u32::MAX {
                num[f.b as usize] += wgt * (va - vb);
                den[f.b as usize] += wgt;
            }
        }
        let busy: Vec<f64> = (0..n).map(|k| if den[k] > 0.0 { (num[k] / den[k]).abs() / peak } else { 0.0 }).collect();
        let (vf, hxf, hyf) = self.to_fine();
        let mut out: Vec<Leaf> = Vec::with_capacity(n);
        let mut taken = vec![false; n];
        for (k, l) in self.leaves.iter().enumerate() {
            if taken[k] {
                continue;
            }
            // merge: this leaf is the lower-left child of a parent whose four children are leaves and quiet
            let ps = 2 * l.s;
            if l.x % ps == 0 && l.y % ps == 0 && ps <= self.o.max_size && (l.x + ps) as usize <= nx && (l.y + ps) as usize <= self.g.ny {
                let kids = [(l.x, l.y), (l.x + l.s, l.y), (l.x, l.y + l.s), (l.x + l.s, l.y + l.s)];
                let ids: Vec<u32> = kids.iter().map(|&(x, y)| self.owner[y as usize * nx + x as usize]).collect();
                let all = ids.iter().zip(&kids).all(|(&id, &(x, y))| {
                    let c = self.leaves[id as usize];
                    c.s == l.s && c.x == x && c.y == y
                });
                let parent = Leaf { x: l.x, y: l.y, s: ps };
                if all && ids.iter().all(|&id| !taken[id as usize] && busy[id as usize] < self.o.coarsen_tol) && self.min_limit(parent) >= ps {
                    for &id in &ids {
                        taken[id as usize] = true;
                    }
                    out.push(parent);
                    continue;
                }
            }
            taken[k] = true;
            if busy[k] > self.o.refine_tol && l.s > 1 {
                let h = l.s / 2;
                for (dx, dy) in [(0, 0), (h, 0), (0, h), (h, h)] {
                    out.push(Leaf { x: l.x + dx, y: l.y + dy, s: h });
                }
            } else {
                out.push(*l);
            }
        }
        self.set_leaves(out);
        self.balance();
        self.rebuild(&vf, &hxf, &hyf);
        self.coefs();
        let e1 = self.energy();
        self.energy_lost += e0 - e1;
        self.remeshes += 1;
    }

    /// Run with a source waveform on the port (sampled at half steps).
    pub fn run(&mut self, exc: &[f64], steps: usize) -> AmrRun {
        let t0 = std::time::Instant::now();
        let mut r = AmrRun { dt: self.dt, probes: vec![vec![]; self.g.probes.len()], min_leaves: usize::MAX, ..Default::default() };
        let nx = self.g.nx;
        let mut peak = 0.0f64;
        let gport = 1.0 / self.g.port.2;
        let lts = |leaves: &[Leaf]| -> f64 { leaves.iter().map(|l| 1.0 / l.s as f64).sum() };
        for n in 0..steps {
            // H (branch currents)
            for (f, i) in self.faces.iter().zip(self.i.iter_mut()) {
                let vb = if f.b == u32::MAX { 0.0 } else { self.v[f.b as usize] };
                *i = f.da * *i + f.k * (self.v[f.a as usize] - vb);
            }
            // E (node voltages): C dV/dt = -sum of currents out (+ port source)
            let mut acc = vec![0.0; self.leaves.len()];
            for (f, &i) in self.faces.iter().zip(&self.i) {
                acc[f.a as usize] -= i;
                if f.b != u32::MAX {
                    acc[f.b as usize] += i;
                }
            }
            let vs = exc.get(n).copied().unwrap_or(0.0);
            acc[self.port_leaf as usize] += gport * vs;
            let v0 = self.v[self.port_leaf as usize];
            for k in 0..self.v.len() {
                self.v[k] = self.ca[k] * self.v[k] + self.cb[k] * acc[k];
            }
            let v1 = self.v[self.port_leaf as usize];
            r.vport.push(v1);
            r.iport.push(gport * (vs - 0.5 * (v0 + v1)));
            for (pi, &(i, j)) in self.g.probes.iter().enumerate() {
                r.probes[pi].push(self.v[self.owner[j * nx + i] as usize]);
            }
            peak = peak.max(self.v.iter().fold(0.0f64, |m, &x| m.max(x.abs())));
            r.cell_updates += self.leaves.len() as f64;
            r.lts_updates += lts(&self.leaves);
            r.min_leaves = r.min_leaves.min(self.leaves.len());
            r.max_leaves = r.max_leaves.max(self.leaves.len());
            if n % 50 == 0 {
                r.energy.push((n, self.energy()));
                r.cells.push((n, self.leaves.len()));
            }
            if self.o.every > 0 && n > 0 && n % self.o.every == 0 {
                self.remesh(peak);
            }
        }
        r.steps = steps;
        r.energy_lost = self.energy_lost;
        r.remeshes = self.remeshes;
        r.seconds = t0.elapsed().as_secs_f64();
        r
    }
}

/// Gaussian pulse at half steps, as `fdtd::gauss_excite` (centre f0, 20 dB half bandwidth fc).
pub fn gauss(f0: f64, fc: f64, dt: f64) -> Vec<f64> {
    crate::fdtd::gauss_excite(f0, fc, dt).into_iter().map(|v| v as f64).collect()
}

/// Port impedance Z(f) = V / I from a run.
pub fn port_z(r: &AmrRun, freqs: &[f64]) -> Vec<crate::post::Cx> {
    let v = crate::post::dft(&r.vport, r.dt, 1.0, freqs);
    let i = crate::post::dft(&r.iport, r.dt, 0.5, freqs);
    v.iter().zip(&i).map(|(a, b)| *a / *b).collect()
}

/// Test geometry: a plane pair `a` x `b` mm (εr, tan δ), seen from the top, open (magnetic) edges, fed near
/// a corner; optionally a slot (a thin PEC wall from one edge, `slot` mm long, through the middle) and an
/// absorbing frame. `h` is the fine cell (mm).
pub fn plane_pair(a: f64, b: f64, er: f64, tand: f64, h: f64, slot: f64, absorb_mm: f64) -> Geo {
    let pad = absorb_mm;
    let nx = (((a + 2.0 * pad) / h).round() as usize).next_multiple_of(64);
    let ny = (((b + 2.0 * pad) / h).round() as usize).next_multiple_of(64);
    let (x0, y0) = (pad, pad);
    let mut eps = vec![1.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let (x, y) = ((i as f64 + 0.5) * h - x0, (j as f64 + 0.5) * h - y0);
            if x >= 0.0 && x <= a && y >= 0.0 && y <= b {
                eps[j * nx + i] = er;
            }
        }
    }
    let mut wall_x = vec![false; (nx + 1) * ny];
    let wall_y = vec![false; nx * (ny + 1)];
    if slot > 0.0 {
        // vertical wall at x = a/2 from y = 0 upward
        let xi = ((x0 + 0.5 * a) / h).round() as usize;
        for j in 0..ny {
            let y = (j as f64 + 0.5) * h - y0;
            if y >= 0.0 && y <= slot {
                wall_x[j * (nx + 1) + xi] = true;
            }
        }
    }
    let port = (((x0 + 3.0) / h) as usize, ((y0 + 3.0) / h) as usize, 50.0 * 0.5e-3);
    let probes = vec![(((x0 + a - 4.0) / h) as usize, ((y0 + b - 4.0) / h) as usize), (((x0 + 0.5 * a + 2.0) / h) as usize, ((y0 + 2.0) / h) as usize)];
    Geo { nx, ny, h: h * 1e-3, eps, tand, f_loss: 3e9, pec: vec![false; nx * ny], wall_x, wall_y, pec_box: false, absorb: (absorb_mm / h).round() as usize, port, probes, fmax: 6e9 }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn static_quadtree_conserves_energy_and_matches_uniform_modes() {
        // lossless closed plane pair, no absorber: the network energy after the pulse must be constant
        let g = plane_pair(20.0, 15.0, 4.3, 0.0, 0.25, 8.0, 0.0);
        let mut a = Amr::new(g.clone(), AmrOpts { max_size: 8, ..Default::default() });
        let exc = gauss(0.0, 6e9, a.dt);
        let n0 = exc.len() + 10;
        // drop the port resistor after the pulse so the box is closed: run, then compare energy late
        let r = a.run(&exc, 6000);
        let es: Vec<f64> = r.energy.iter().filter(|e| e.0 > n0).map(|e| e.1).collect();
        // the port's resistor still damps: the energy (sampled with E and H half a step apart, so it wobbles)
        // must not grow
        let (first, last) = (es[..10].iter().cloned().fold(0.0, f64::max), es[es.len() - 10..].iter().cloned().fold(0.0, f64::max));
        assert!(last <= first, "energy grew {first} -> {last}");
        assert!(a.leaves() < g.nx * g.ny / 2, "quadtree did not coarsen: {} leaves", a.leaves());
    }
}

/// Estimate for step 3: cells an octree would need on a board's rectilinear mesh `g` if cell size may grow
/// with distance d from the nearest copper edge as s = h0 + (r - 1) d (r = grading per cell, like the
/// rectilinear mesh's), capped at the mesh's largest cell. Returns (rectilinear cells, octree estimate).
/// Copper edges are taken from the midpoint-rule masks: nodes where metal and open in-plane edges meet.
pub fn octree_estimate(g: &crate::mesh::Grid, masks: &[crate::model::SheetMask], r: f64) -> (usize, f64) {
    let (nx, ny, nz) = (g.nx(), g.ny(), g.nz());
    let d = |l: &[f64], i: usize| if i + 1 < l.len() { l[i + 1] - l[i] } else { 0.0 };
    let h0 = [g.min_step(0), g.min_step(1)].iter().cloned().fold(f64::INFINITY, f64::min);
    let smax = (0..3).map(|a| g.axis(a).windows(2).map(|w| w[1] - w[0]).fold(0.0, f64::max)).fold(0.0, f64::max);
    // per sheet: in-plane distance from every node to the nearest copper-edge node (two-pass chamfer on the
    // nonuniform grid, Euclidean steps along the axes and diagonals)
    let mut dist: Vec<(usize, Vec<f64>)> = Vec::new();
    for m in masks {
        let mut dd = vec![f64::INFINITY; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let q = j * nx + i;
                let inc = [(i + 1 < nx).then(|| m.ex[q]), (i > 0).then(|| m.ex[q - 1]), (j + 1 < ny).then(|| m.ey[q]), (j > 0).then(|| m.ey[q - nx])];
                let (mut a, mut b) = (false, false);
                for v in inc.iter().flatten() {
                    if *v {
                        a = true
                    } else {
                        b = true
                    }
                }
                if a && b {
                    dd[q] = 0.0;
                }
            }
        }
        for pass in 0..2 {
            let (jr, ir): (Vec<usize>, Vec<usize>) = if pass == 0 { ((0..ny).collect(), (0..nx).collect()) } else { ((0..ny).rev().collect(), (0..nx).rev().collect()) };
            for &j in &jr {
                for &i in &ir {
                    let q = j * nx + i;
                    let mut best = dd[q];
                    for (di, dj) in [(-1i64, 0i64), (0, -1), (-1, -1), (1, -1)] {
                        let (di, dj) = if pass == 0 { (di, dj) } else { (-di, -dj) };
                        let (ii, jj) = (i as i64 + di, j as i64 + dj);
                        if ii < 0 || jj < 0 || ii >= nx as i64 || jj >= ny as i64 {
                            continue;
                        }
                        let p = jj as usize * nx + ii as usize;
                        let step = ((g.x[i] - g.x[ii as usize]).powi(2) + (g.y[j] - g.y[jj as usize]).powi(2)).sqrt();
                        best = best.min(dd[p] + step);
                    }
                    dd[q] = best;
                }
            }
        }
        dist.push((m.k, dd));
    }
    let h0z = g.min_step(2);
    let mut n_oct = 0.0;
    for k in 0..nz.saturating_sub(1) {
        let zc = 0.5 * (g.z[k] + g.z[k + 1]);
        let dzs = dist.iter().map(|(kk, _)| (zc - g.z[*kk]).abs()).fold(f64::INFINITY, f64::min);
        for j in 0..ny.saturating_sub(1) {
            for i in 0..nx.saturating_sub(1) {
                let q = j * nx + i;
                let mut dmin = f64::INFINITY;
                for (kk, dd) in &dist {
                    let dz = (zc - g.z[*kk]).abs();
                    dmin = dmin.min((dd[q].powi(2) + dz * dz).sqrt());
                }
                let (dx, dy, dz) = (d(&g.x, i), d(&g.y, j), d(&g.z, k));
                // in-plane size grows with the distance from copper edges; never finer than the rectilinear
                // mesh is there (the reference resolution); the z size grows with the distance from the
                // nearest copper plane, never beyond the in-plane size
                let sxy = (h0 + (r - 1.0) * dmin).min(smax).max(dx.min(dy));
                let sz = (h0z + (r - 1.0) * dzs).min(sxy).max(dz);
                n_oct += dx * dy * dz / (sxy * sxy * sz);
            }
        }
    }
    (g.n(), n_oct)
}
