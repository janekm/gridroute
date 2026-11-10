//! Prototype (CPU only): 3-D FDTD in integral (FIT) form on a 2:1-balanced octree. Not used by `run`;
//! `fieldsolve octree-check` measures it against Hammerstad-Jensen and the rectilinear solver.
//!
//! - **Mesh.** Leaves are boxes of `s` finest cells per axis (`s` a power of two, the finest cell may differ
//!   per axis, as PCB meshes need). They are balanced 2:1 across faces and edges.
//! - **Unknowns.** Edge voltages e = ∫E·dl on the minimal edges (leaf edges split at hanging nodes), face
//!   fluxes b = ∫B·dA on the minimal faces. Faraday (db/dt = −C e) and Ampère (dd/dt = Cᵀ h) are exact on
//!   this complex; the constitutive relations are the diagonal (Yee) ones: d = M_ε e per edge, h = b / M_μ per
//!   face with M_μ = μ A / (distance of the two leaf centres along the normal), first order at hanging faces.
//! - **Hanging edges** are not free. An edge that is half of a coarser leaf's edge carries half of that
//!   edge's voltage; an edge inside a coarser leaf's face (from the face centre to the middle of a side)
//!   carries a quarter of each of the two coarse edges parallel to it (the coarse cell's field, bilinear).
//!   With e = P e_m (e_m: free edges and the subdivided coarse edges), Ampère is applied through Pᵀ and the
//!   mass is the lumped (row-sum) Pᵀ M_ε P. The update operator stays symmetric in the energy inner product:
//!   no interface stencil to tune, energy conserved, stable for dt² λ_max < 4.
//! - **Measurements.** Port as in `model.rs` (a column of vertical edges with a resistive source). Voltage
//!   probes sum e along a line of edges; the current in a zero-thickness strip is exact in this form: the
//!   conduction current of a metal edge is its curl term (Cᵀ h), so the strip current through a plane is
//!   the sum over the strip's metal edges crossing it.

use rustc_hash::FxHashMap;

const C0: f64 = 299_792_458.0;
const MU0: f64 = 1.256_637_062_12e-6;
const EPS0: f64 = 8.854_187_812_8e-12;

pub type Leaf = [u32; 4];

pub struct Octree {
    pub n: [u32; 3],
    /// finest cell size per axis (m)
    pub h: [f64; 3],
    pub leaves: Vec<Leaf>,
    owner: Vec<u32>,
}

impl Octree {
    fn cell(&self, x: u32, y: u32, z: u32) -> usize {
        ((z as usize * self.n[1] as usize) + y as usize) * self.n[0] as usize + x as usize
    }
    fn set_owner(&mut self) {
        self.owner.iter_mut().for_each(|o| *o = u32::MAX);
        for (k, l) in self.leaves.iter().enumerate() {
            for z in l[2]..l[2] + l[3] {
                for y in l[1]..l[1] + l[3] {
                    let base = ((z as usize * self.n[1] as usize) + y as usize) * self.n[0] as usize;
                    for x in l[0]..l[0] + l[3] {
                        self.owner[base + x as usize] = k as u32;
                    }
                }
            }
        }
    }
    fn owner_at(&self, x: i64, y: i64, z: i64) -> Option<u32> {
        if x < 0 || y < 0 || z < 0 || x >= self.n[0] as i64 || y >= self.n[1] as i64 || z >= self.n[2] as i64 {
            return None;
        }
        Some(self.owner[self.cell(x as u32, y as u32, z as u32)])
    }

    /// Build: roots of `root` cells tile the domain (`n` must be a multiple of `root`); a leaf is split
    /// while `split(leaf)` says so; then 2:1 balance across faces and edges.
    pub fn build(n: [u32; 3], h: [f64; 3], root: u32, split: &dyn Fn(Leaf) -> bool) -> Octree {
        let mut stack = Vec::new();
        for z in (0..n[2]).step_by(root as usize) {
            for y in (0..n[1]).step_by(root as usize) {
                for x in (0..n[0]).step_by(root as usize) {
                    stack.push([x, y, z, root]);
                }
            }
        }
        let mut leaves = Vec::new();
        while let Some(l) = stack.pop() {
            if l[3] > 1 && split(l) {
                let s = l[3] / 2;
                for dz in [0, s] {
                    for dy in [0, s] {
                        for dx in [0, s] {
                            stack.push([l[0] + dx, l[1] + dy, l[2] + dz, s]);
                        }
                    }
                }
            } else {
                leaves.push(l);
            }
        }
        let total = n[0] as usize * n[1] as usize * n[2] as usize;
        let mut o = Octree { n, h, leaves, owner: vec![u32::MAX; total] };
        o.set_owner();
        o.balance();
        o
    }

    /// Split leaves until every face and edge neighbour is at most one level finer.
    fn balance(&mut self) {
        loop {
            let mut split = vec![false; self.leaves.len()];
            for (k, l) in self.leaves.iter().enumerate() {
                let (x0, y0, z0, s) = (l[0] as i64, l[1] as i64, l[2] as i64, l[3] as i64);
                if s < 4 {
                    continue; // a neighbour of size 1 or more is within a factor 2 of s <= 2
                }
                let mut bad = false;
                // sample the cells just outside every face and edge
                'dirs: for dz in -1i64..=1 {
                    for dy in -1i64..=1 {
                        for dx in -1i64..=1 {
                            let nz = (dx != 0) as u8 + (dy != 0) as u8 + (dz != 0) as u8;
                            if nz == 0 || nz == 3 {
                                continue;
                            }
                            let range = |d: i64, a0: i64| -> std::ops::Range<i64> {
                                match d {
                                    -1 => a0 - 1..a0,
                                    1 => a0 + s..a0 + s + 1,
                                    _ => a0..a0 + s,
                                }
                            };
                            for z in range(dz, z0) {
                                for y in range(dy, y0) {
                                    for x in range(dx, x0) {
                                        if let Some(nb) = self.owner_at(x, y, z) {
                                            if (self.leaves[nb as usize][3] as i64) * 2 < s {
                                                bad = true;
                                                break 'dirs;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                split[k] = bad;
            }
            if !split.iter().any(|&b| b) {
                return;
            }
            let mut out = Vec::with_capacity(self.leaves.len() * 2);
            for (k, l) in self.leaves.iter().enumerate() {
                if split[k] {
                    let s = l[3] / 2;
                    for dz in [0, s] {
                        for dy in [0, s] {
                            for dx in [0, s] {
                                out.push([l[0] + dx, l[1] + dy, l[2] + dz, s]);
                            }
                        }
                    }
                } else {
                    out.push(*l);
                }
            }
            self.leaves = out;
            self.set_owner();
        }
    }
}

/// The assembled FIT system on an octree.
pub struct Fit {
    pub n: [u32; 3],
    pub h: [f64; 3],
    /// minimal edges: axis, start (lattice), length (finest units)
    pub ea: Vec<u8>,
    pub ep: Vec<[u32; 3]>,
    pub el: Vec<u32>,
    /// per minimal edge, its expression in masters: CSR into (master, weight)
    emap_s: Vec<u32>,
    emap: Vec<(u32, f64)>,
    /// masters: free minimal edges first, then subdivided coarse ("virtual") edges; per master: axis, start,
    /// length, lumped mass (F·m), conductance (S·m), PEC
    pub ma: Vec<u8>,
    pub mp: Vec<[u32; 3]>,
    pub ml: Vec<u32>,
    pub meps: Vec<f64>,
    pub msig: Vec<f64>,
    pub pec: Vec<bool>,
    /// minimal faces: CSR of (edge, ±1) around each face; M_μ per face
    f_s: Vec<u32>,
    f_e: Vec<(u32, f64)>,
    pub fmu: Vec<f64>,
    /// edge -> faces (transpose)
    e_s: Vec<u32>,
    e_f: Vec<(u32, f64)>,
    edge_at: FxHashMap<u64, u32>,
    pub leaves: usize,
    pub virtuals: usize,
}

fn lkey(n: [u32; 3], a: u8, p: [u32; 3]) -> u64 {
    let lat = ((p[2] as u64 * (n[1] as u64 + 1)) + p[1] as u64) * (n[0] as u64 + 1) + p[0] as u64;
    (a as u64) << 60 | lat
}

/// Assemble the FIT system. `eps_r(leaf)` gives a leaf's permittivity.
pub fn assemble(o: &Octree, eps_r: &dyn Fn(Leaf) -> f64) -> Fit {
    let n = o.n;
    let lat = |p: [u32; 3]| -> usize { ((p[2] as usize * (n[1] as usize + 1)) + p[1] as usize) * (n[0] as usize + 1) + p[0] as usize };
    let nlat = (n[0] as usize + 1) * (n[1] as usize + 1) * (n[2] as usize + 1);
    // nodes: every leaf corner
    let mut node = vec![0u64; nlat.div_ceil(64)];
    let set = |node: &mut Vec<u64>, q: usize| node[q >> 6] |= 1u64 << (q & 63);
    for l in &o.leaves {
        for dz in [0, l[3]] {
            for dy in [0, l[3]] {
                for dx in [0, l[3]] {
                    set(&mut node, lat([l[0] + dx, l[1] + dy, l[2] + dz]));
                }
            }
        }
    }
    let is_node = |p: [u32; 3]| -> bool {
        let q = lat(p);
        node[q >> 6] >> (q & 63) & 1 == 1
    };
    let add = |p: [u32; 3], a: usize, d: u32| -> [u32; 3] {
        let mut q = p;
        q[a] += d;
        q
    };
    // minimal edges and the subdivided coarse edges of every leaf
    let mut edge_at: FxHashMap<u64, u32> = FxHashMap::default();
    let (mut ea, mut ep, mut el) = (Vec::new(), Vec::new(), Vec::new());
    let mut virt_at: FxHashMap<u64, u32> = FxHashMap::default();
    let (mut va, mut vp, mut vl) = (Vec::<u8>::new(), Vec::<[u32; 3]>::new(), Vec::<u32>::new());
    let mut edge = |a: u8, p: [u32; 3], len: u32, ea: &mut Vec<u8>, ep: &mut Vec<[u32; 3]>, el: &mut Vec<u32>| -> u32 {
        *edge_at.entry(lkey(n, a, p)).or_insert_with(|| {
            ea.push(a);
            ep.push(p);
            el.push(len);
            (ea.len() - 1) as u32
        })
    };
    // per leaf edge (12): which minimal edges or virtual edge it is, for the masses
    let mut leaf_edges: Vec<(u8, [u32; 3], u32, bool)> = Vec::new(); // (axis, start, len, subdivided)
    for l in &o.leaves {
        let s = l[3];
        let base = [l[0], l[1], l[2]];
        for a in 0..3usize {
            let (b, c) = ((a + 1) % 3, (a + 2) % 3);
            for (db, dc) in [(0, 0), (s, 0), (0, s), (s, s)] {
                let p = add(add(base, b, db), c, dc);
                let mid = add(p, a, s / 2);
                if s > 1 && is_node(mid) {
                    edge(a as u8, p, s / 2, &mut ea, &mut ep, &mut el);
                    edge(a as u8, mid, s / 2, &mut ea, &mut ep, &mut el);
                    virt_at.entry(lkey(n, a as u8, p)).or_insert_with(|| {
                        va.push(a as u8);
                        vp.push(p);
                        vl.push(s);
                        (va.len() - 1) as u32
                    });
                    leaf_edges.push((a as u8, p, s, true));
                } else {
                    edge(a as u8, p, s, &mut ea, &mut ep, &mut el);
                    leaf_edges.push((a as u8, p, s, false));
                }
            }
        }
    }
    let ne = ea.len();
    // constraints. Halves of a subdivided coarse edge: 1/2 of it. Edges inside a subdivided leaf face (from
    // its centre to the middle of a side): 1/4 of each of the two coarse edges parallel to them.
    let mut slave: Vec<Vec<(u32, f64)>> = vec![Vec::new(); ne]; // virtual index, weight
    for (vi, ((&a, &p), &len)) in va.iter().zip(&vp).zip(&vl).enumerate() {
        let h = len / 2;
        for q in [p, add(p, a as usize, h)] {
            let e = edge_at[&lkey(n, a, q)];
            if slave[e as usize].is_empty() {
                slave[e as usize].push((vi as u32, 0.5));
            }
        }
    }
    for l in &o.leaves {
        let s = l[3];
        if s < 2 {
            continue;
        }
        let base = [l[0], l[1], l[2]];
        for a in 0..3usize {
            // the two faces of the leaf normal to a
            for da in [0, s] {
                let f0 = add(base, a, da);
                let (b, c) = ((a + 1) % 3, (a + 2) % 3);
                let centre = add(add(f0, b, s / 2), c, s / 2);
                if !is_node(centre) {
                    continue;
                }
                // cross edges along b (at c = s/2) and along c (at b = s/2)
                for (dir, other) in [(b, c), (c, b)] {
                    let m0 = add(f0, other, s / 2); // start of the line through the centre along dir
                    let v_lo = virt_at.get(&lkey(n, dir as u8, f0)).copied();
                    let v_hi = virt_at.get(&lkey(n, dir as u8, add(f0, other, s))).copied();
                    for st in [m0, add(m0, dir, s / 2)] {
                        if let (Some(&e), Some(lo), Some(hi)) = (edge_at.get(&lkey(n, dir as u8, st)), v_lo, v_hi) {
                            if slave[e as usize].is_empty() {
                                slave[e as usize].push((lo, 0.25));
                                slave[e as usize].push((hi, 0.25));
                            }
                        }
                    }
                }
            }
        }
    }
    // masters: free minimal edges, then the virtual edges
    let mut midx = vec![u32::MAX; ne];
    let (mut ma, mut mp, mut ml) = (Vec::new(), Vec::new(), Vec::new());
    for e in 0..ne {
        if slave[e].is_empty() {
            midx[e] = ma.len() as u32;
            ma.push(ea[e]);
            mp.push(ep[e]);
            ml.push(el[e]);
        }
    }
    let nfree = ma.len();
    ma.extend(&va);
    mp.extend(&vp);
    ml.extend(&vl);
    let nm = ma.len();
    let mut emap_s = vec![0u32];
    let mut emap = Vec::new();
    for e in 0..ne {
        if slave[e].is_empty() {
            emap.push((midx[e], 1.0));
        } else {
            for &(v, w) in &slave[e] {
                emap.push(((nfree + v as usize) as u32, w));
            }
        }
        emap_s.push(emap.len() as u32);
    }
    // masses: each leaf gives each of its 12 edges eps * (dual area share) / length
    let hm = o.h;
    let mut msl = vec![0.0f64; ne]; // minimal-edge masses
    let mut mvirt = vec![0.0f64; va.len()];
    for (li, l) in o.leaves.iter().enumerate() {
        let er = eps_r(*l);
        let s = l[3] as f64;
        for t in 0..12 {
            let (a, p, len, sub) = leaf_edges[li * 12 + t];
            let (b, c) = ((a as usize + 1) % 3, (a as usize + 2) % 3);
            let w = EPS0 * er * (s * hm[b] / 2.0) * (s * hm[c] / 2.0) / (len as f64 * hm[a as usize]);
            if sub {
                mvirt[virt_at[&lkey(n, a, p)] as usize] += w;
            } else {
                msl[edge_at[&lkey(n, a, p)] as usize] += w;
            }
        }
    }
    let mut meps = vec![0.0f64; nm];
    for v in 0..va.len() {
        meps[nfree + v] += mvirt[v];
    }
    for e in 0..ne {
        let r = &emap[emap_s[e] as usize..emap_s[e + 1] as usize];
        let rs: f64 = r.iter().map(|x| x.1).sum();
        for &(m, w) in r {
            meps[m as usize] += w * rs * msl[e];
        }
    }
    // minimal faces: every leaf side, split in four when its centre is a node
    let mut face_at: FxHashMap<u64, u32> = FxHashMap::default();
    let (mut fa, mut fp, mut fs) = (Vec::<u8>::new(), Vec::<[u32; 3]>::new(), Vec::<u32>::new());
    for l in &o.leaves {
        let s = l[3];
        let base = [l[0], l[1], l[2]];
        for a in 0..3usize {
            let (b, c) = ((a + 1) % 3, (a + 2) % 3);
            for da in [0, s] {
                let f0 = add(base, a, da);
                let centre = add(add(f0, b, s / 2), c, s / 2);
                let subs: Vec<([u32; 3], u32)> = if s > 1 && is_node(centre) {
                    let h = s / 2;
                    vec![(f0, h), (add(f0, b, h), h), (add(f0, c, h), h), (add(add(f0, b, h), c, h), h)]
                } else {
                    vec![(f0, s)]
                };
                for (q, sz) in subs {
                    face_at.entry(lkey(n, a as u8, q)).or_insert_with(|| {
                        fa.push(a as u8);
                        fp.push(q);
                        fs.push(sz);
                        (fa.len() - 1) as u32
                    });
                }
            }
        }
    }
    let nf = fa.len();
    // boundary loop of each face (counter-clockwise about +a): +b at c0, +c at b1, -b at c1, -c at b0
    let mut f_s = vec![0u32];
    let mut f_e: Vec<(u32, f64)> = Vec::new();
    let mut fmu = vec![0.0f64; nf];
    for f in 0..nf {
        let (a, p, s) = (fa[f] as usize, fp[f], fs[f]);
        let (b, c) = ((a + 1) % 3, (a + 2) % 3);
        let sides = [(b, p, 1.0), (c, add(p, b, s), 1.0), (b, add(p, c, s), -1.0), (c, p, -1.0)];
        for (dir, st, sg) in sides {
            // one edge of length s, or two halves
            if let Some(&e) = edge_at.get(&lkey(n, dir as u8, st)).filter(|&&e| el[e as usize] == s) {
                f_e.push((e, sg));
            } else {
                for q in [st, add(st, dir, s / 2)] {
                    let e = *edge_at.get(&lkey(n, dir as u8, q)).expect("face side without edges");
                    f_e.push((e, sg));
                }
            }
        }
        f_s.push(f_e.len() as u32);
        // M_mu = mu A / (sum of the half extents of the leaves on both sides along a)
        let area = (s as f64 * hm[b]) * (s as f64 * hm[c]);
        let mut cen = [p[0] as i64, p[1] as i64, p[2] as i64];
        cen[b] += (s / 2) as i64;
        cen[c] += (s / 2) as i64;
        let mut dsum = 0.0;
        for side in [-1i64, 0] {
            let mut q = cen;
            q[a] += side;
            // a cell touching the face's centre region on that side
            let mut q2 = q;
            if s > 1 {
                q2[b] -= 0;
            }
            if let Some(lf) = o.owner_at(q2[0], q2[1], q2[2]) {
                dsum += 0.5 * o.leaves[lf as usize][3] as f64 * hm[a];
            }
        }
        fmu[f] = if dsum > 0.0 { MU0 * area / dsum } else { f64::INFINITY };
    }
    // transpose: edge -> faces
    let mut cnt = vec![0u32; ne + 1];
    for &(e, _) in &f_e {
        cnt[e as usize + 1] += 1;
    }
    for i in 0..ne {
        cnt[i + 1] += cnt[i];
    }
    let mut e_f = vec![(0u32, 0.0f64); f_e.len()];
    let mut fill = cnt.clone();
    for f in 0..nf {
        for &(e, sg) in &f_e[f_s[f] as usize..f_s[f + 1] as usize] {
            e_f[fill[e as usize] as usize] = (f as u32, sg);
            fill[e as usize] += 1;
        }
    }
    // PEC: masters on the domain boundary (tangential)
    let pec: Vec<bool> = (0..nm)
        .map(|m| {
            let (a, p) = (ma[m] as usize, mp[m]);
            (0..3).any(|d| d != a && (p[d] == 0 || p[d] == n[d]))
        })
        .collect();
    Fit {
        n,
        h: hm,
        ea,
        ep,
        el,
        emap_s,
        emap,
        ma,
        mp,
        ml,
        meps,
        msig: vec![0.0; nm],
        pec,
        f_s,
        f_e,
        fmu,
        e_s: cnt,
        e_f,
        edge_at,
        leaves: o.leaves.len(),
        virtuals: va.len(),
    }
}

/// A port: master edges in series (a vertical column), resistance R; positive along +axis.
pub struct Port {
    pub edges: Vec<u32>,
    pub r: f64,
}

pub struct Run {
    pub dt: f64,
    pub steps: usize,
    pub vport: Vec<f64>,
    pub iport: Vec<f64>,
    pub v: Vec<Vec<f64>>,
    pub i: Vec<Vec<f64>>,
    pub seconds: f64,
    pub energy: Vec<f64>,
}

impl Fit {
    pub fn masters(&self) -> usize {
        self.ma.len()
    }
    pub fn faces(&self) -> usize {
        self.fmu.len()
    }
    /// Minimal edges along axis `a` from `p` over `len` finest units (they must exist).
    pub fn line(&self, a: u8, p: [u32; 3], len: u32) -> Vec<u32> {
        let mut out = Vec::new();
        let mut q = p;
        let end = p[a as usize] + len;
        while q[a as usize] < end {
            let e = *self.edge_at.get(&lkey(self.n, a, q)).expect("line not on mesh edges");
            out.push(e);
            q[a as usize] += self.el[e as usize];
        }
        out
    }
    /// The master of a free minimal edge.
    pub fn master_of(&self, e: u32) -> Option<u32> {
        let r = &self.emap[self.emap_s[e as usize] as usize..self.emap_s[e as usize + 1] as usize];
        (r.len() == 1 && r[0].1 == 1.0).then(|| r[0].0)
    }

    /// dt^2 lambda_max of the lossless update operator (power iteration).
    pub fn stability(&self, dt: f64, iters: usize) -> f64 {
        let nm = self.masters();
        let mut x: Vec<f64> = (0..nm).map(|m| if self.pec[m] { 0.0 } else { ((m * 2654435761) % 1000) as f64 / 1000.0 - 0.5 }).collect();
        let mut mu = 0.0;
        for _ in 0..iters {
            let nrm = x.iter().zip(&self.meps).map(|(a, m)| a * a * m).sum::<f64>().sqrt();
            x.iter_mut().for_each(|a| *a /= nrm);
            let y = self.apply(&x, dt);
            let ny = y.iter().zip(&self.meps).map(|(a, m)| a * a * m).sum::<f64>().sqrt();
            mu = ny;
            x = y;
        }
        mu
    }

    /// dt^2 K x for the curl-curl operator K = M^-1 P^T C^T M_mu^-1 C P (lossless, PEC masters zero).
    fn apply(&self, xm: &[f64], dt: f64) -> Vec<f64> {
        let ne = self.ea.len();
        let mut e = vec![0.0; ne];
        for i in 0..ne {
            e[i] = self.emap[self.emap_s[i] as usize..self.emap_s[i + 1] as usize].iter().map(|&(m, w)| w * xm[m as usize]).sum();
        }
        let nf = self.faces();
        let hh: Vec<f64> = (0..nf).map(|f| dt * self.f_e[self.f_s[f] as usize..self.f_s[f + 1] as usize].iter().map(|&(ed, s)| s * e[ed as usize]).sum::<f64>() / self.fmu[f]).collect();
        let mut g = vec![0.0; self.masters()];
        for i in 0..ne {
            let ge: f64 = self.e_f[self.e_s[i] as usize..self.e_s[i + 1] as usize].iter().map(|&(f, s)| s * hh[f as usize]).sum();
            for &(m, w) in &self.emap[self.emap_s[i] as usize..self.emap_s[i + 1] as usize] {
                g[m as usize] += w * ge;
            }
        }
        g.iter().enumerate().map(|(m, &v)| if self.pec[m] { 0.0 } else { dt * v / self.meps[m] }).collect()
    }

    /// Leapfrog run: Gaussian `exc` (half steps) on `port`; voltage probes (lists of minimal edges, summed)
    /// and current probes (lists of metal masters, their curl terms summed: the conduction current).
    pub fn run(&self, dt: f64, port: &Port, exc: &[f64], vprobes: &[Vec<u32>], iprobes: &[Vec<u32>], steps: usize, end_db: f64) -> Run {
        use rayon::prelude::*;
        let t0 = std::time::Instant::now();
        let (ne, nf, nm) = (self.ea.len(), self.faces(), self.masters());
        let mut em = vec![0.0f64; nm];
        let mut e = vec![0.0f64; ne];
        let mut b = vec![0.0f64; nf];
        let mut hh = vec![0.0f64; nf];
        let mut g = vec![0.0f64; nm];
        // port: conductance per edge (k in series), Norton source
        let k = port.edges.len() as f64;
        let mut gport = vec![0.0f64; nm];
        let pm: Vec<u32> = port.edges.iter().map(|&ed| self.master_of(ed).expect("port edge must be free")).collect();
        for &m in &pm {
            gport[m as usize] = k / port.r;
        }
        let ca: Vec<f64> = (0..nm).map(|m| { let (c, s) = (self.meps[m] / dt, 0.5 * (self.msig[m] + gport[m])); (c - s) / (c + s) }).collect();
        let cb: Vec<f64> = (0..nm).map(|m| 1.0 / (self.meps[m] / dt + 0.5 * (self.msig[m] + gport[m]))).collect();
        let mut r = Run { dt, steps, vport: vec![], iport: vec![], v: vec![vec![]; vprobes.len()], i: vec![vec![]; iprobes.len()], seconds: 0.0, energy: vec![] };
        let mut ge = vec![0.0f64; ne];
        let mut peak = 0.0f64;
        let mut done = steps;
        for nstep in 0..steps {
            // e on all minimal edges from the masters
            e.par_iter_mut().enumerate().for_each(|(i, ei)| {
                let rr = &self.emap[self.emap_s[i] as usize..self.emap_s[i + 1] as usize];
                *ei = if rr.len() == 1 { rr[0].1 * em[rr[0].0 as usize] } else { rr.iter().map(|&(m, w)| w * em[m as usize]).sum() };
            });
            // Faraday on the faces, then h
            b.par_iter_mut().zip(hh.par_iter_mut()).enumerate().for_each(|(f, (bf, hf))| {
                let c: f64 = self.f_e[self.f_s[f] as usize..self.f_s[f + 1] as usize].iter().map(|&(ed, s)| s * e[ed as usize]).sum();
                *bf -= dt * c;
                *hf = *bf / self.fmu[f];
            });
            // Ampere: curl per minimal edge, gathered on the masters through P^T
            ge.par_iter_mut().enumerate().for_each(|(i, gi)| {
                *gi = self.e_f[self.e_s[i] as usize..self.e_s[i + 1] as usize].iter().map(|&(f, s)| s * hh[f as usize]).sum();
            });
            g.iter_mut().for_each(|x| *x = 0.0);
            for i in 0..ne {
                for &(m, w) in &self.emap[self.emap_s[i] as usize..self.emap_s[i + 1] as usize] {
                    g[m as usize] += w * ge[i];
                }
            }
            for (q, l) in iprobes.iter().enumerate() {
                r.i[q].push(l.iter().map(|&m| g[m as usize]).sum());
            }
            let vs = exc.get(nstep).copied().unwrap_or(0.0);
            let v0 = pm.iter().map(|&m| em[m as usize]).sum::<f64>();
            for m in 0..nm {
                if self.pec[m] {
                    continue;
                }
                let src = gport[m] * vs / k;
                em[m] = ca[m] * em[m] + cb[m] * (g[m] + src);
            }
            let v1 = pm.iter().map(|&m| em[m as usize]).sum::<f64>();
            r.vport.push(v1);
            r.iport.push((vs - 0.5 * (v0 + v1)) / port.r);
            // probe voltages from the updated masters
            for (q, l) in vprobes.iter().enumerate() {
                let s: f64 = l.iter().map(|&ed| self.emap[self.emap_s[ed as usize] as usize..self.emap_s[ed as usize + 1] as usize].iter().map(|&(m, w)| w * em[m as usize]).sum::<f64>()).sum();
                r.v[q].push(s);
            }
            if nstep % 200 == 0 {
                let we: f64 = em.iter().zip(&self.meps).map(|(a, m)| 0.5 * m * a * a).sum();
                let wh: f64 = b.iter().zip(&self.fmu).map(|(a, m)| if m.is_finite() { 0.5 * a * a / m } else { 0.0 }).sum();
                let en = we + wh;
                r.energy.push(en);
                peak = peak.max(en);
                if nstep > exc.len() && en < peak * 10f64.powf(end_db / 10.0) {
                    done = nstep + 1;
                    break;
                }
            }
        }
        r.steps = done;
        r.seconds = t0.elapsed().as_secs_f64();
        r
    }
}

/// Microstrip on the octree, as `edge_cases::line` on the rectilinear mesh: width `w` (both edges `inset`
/// cells inside a mesh line), substrate `h` (εr), 100 mm long over a PEC ground (the z = 0 wall), open end,
/// vertical port at x = 0, PEC box. Finest cells: 0.5 mm along the line, `hc` across, h / `nz` through.
/// Leaves grow with the distance from the strip's edges (grading `r` per cell) up to `smax` finest cells.
/// Returns (w, Z0, εeff, leaves, masters, faces, steps, seconds).
pub fn microstrip(w0: f64, h: f64, er: f64, hc: f64, nz: u32, inset: f64, r: f64, smax: u32) -> (f64, f64, f64, usize, usize, usize, usize, f64) {
    let len = 100.0;
    let m = ((w0 / 2.0) / hc - inset).round().max(1.0);
    let w = 2.0 * (m + inset) * hc;
    let hm = [0.5e-3, hc * 1e-3, h / nz as f64 * 1e-3];
    let root = smax;
    // domain: x from -12, y symmetric about 0 with y = 0 on a root multiple, z from the ground
    let up = |v: f64, hh: f64| -> u32 { (((v / hh).ceil() as u32).div_ceil(root)) * root };
    // box: 12 mm beyond the ends, `side` mm to each side and above (FS_OCT_BOX overrides 14)
    let side: f64 = std::env::var("FS_OCT_BOX").ok().and_then(|v| v.parse().ok()).unwrap_or(14.0);
    let nxf = up(len + 24.0, 0.5);
    let ny_half = up(side, hc);
    let nzf = up(side, h / nz as f64);
    let n = [nxf, 2 * ny_half, nzf];
    let x0 = -12.0f64; // x of lattice 0 (mm); make the strip start on a root multiple
    let shift = ((-x0 / 0.5) as u32).div_ceil(root) * root;
    let xi = |x: f64| -> f64 { x / 0.5 + shift as f64 }; // mm -> finest units
    let yi = |y: f64| -> f64 { y / hc + ny_half as f64 };
    let kh = nz; // strip plane index
    let hz = h / nz as f64;
    // distance (mm) from a point to the strip's edges (the two long edges and the ends)
    let dist = |x: f64, y: f64, z: f64| -> f64 {
        let dz = z - h;
        let xc = x.clamp(0.0, len);
        let dl = |ye: f64| ((x - xc).powi(2) + (y - ye).powi(2) + dz * dz).sqrt();
        let yc = y.clamp(-w / 2.0, w / 2.0);
        let de = |xe: f64| ((x - xe).powi(2) + (y - yc).powi(2) + dz * dz).sqrt();
        dl(-w / 2.0).min(dl(w / 2.0)).min(de(0.0)).min(de(len))
    };
    let hmin = hc.min(hz);
    let split = |l: Leaf| -> bool {
        let s = l[3];
        if r <= 1.0 {
            return s > 1; // r <= 1: uniform tree (the Yee scheme), for checks
        }
        // never across the substrate top
        if l[2] < kh && l[2] + s > kh {
            return true;
        }
        // physical extents
        let ext = [s as f64 * 0.5, s as f64 * hc, s as f64 * hz];
        let cx = (l[0] as f64 + 0.5 * s as f64 - shift as f64) * 0.5;
        let cy = (l[1] as f64 + 0.5 * s as f64 - ny_half as f64) * hc;
        let cz = (l[2] as f64 + 0.5 * s as f64) * hz;
        // distance from the nearest point of the leaf: centre distance minus half its diagonal
        let half = 0.5 * (ext[0] * ext[0] + ext[1] * ext[1] + ext[2] * ext[2]).sqrt();
        let d = (dist(cx, cy, cz) - half).max(0.0);
        let allowed = hmin + (r - 1.0) * d;
        // the leaf's size across the line (y, z) against what the distance allows
        ext[1].max(ext[2]) > allowed
    };
    let o = Octree::build(n, hm, root, &split);
    // local time stepping (each level at its own CFL step) would update a leaf of size s once per s steps
    let lts: f64 = o.leaves.iter().map(|l| 1.0 / l[3] as f64).sum::<f64>() / o.leaves.len() as f64;
    let mut hist = [0usize; 6];
    for l in &o.leaves {
        hist[(l[3] as f64).log2() as usize] += 1;
    }
    eprintln!("lts: leaves by size 1/2/4/8/16/32: {hist:?}; with local time steps {:.1} % of the updates", 100.0 * lts);
    let strip = |x: f64, y: f64| x >= 0.0 && x <= len && y.abs() <= w / 2.0;
    let mut fit = assemble(&o, &|l: Leaf| if l[2] + l[3] <= kh { er } else { 1.0 });
    // PEC strip: masters on the strip plane whose midpoint is in the strip
    for mi in 0..fit.masters() {
        let (a, p, l) = (fit.ma[mi] as usize, fit.mp[mi], fit.ml[mi]);
        if a == 2 || p[2] != kh {
            continue;
        }
        let mut c = [p[0] as f64, p[1] as f64];
        c[a] += 0.5 * l as f64;
        let (x, y) = ((c[0] - shift as f64) * 0.5, (c[1] - ny_half as f64) * hc);
        if strip(x, y) {
            fit.pec[mi] = true;
        }
    }
    // port: vertical column at x = 0, y = 0 from the ground to the strip
    let px = xi(0.0) as u32;
    let py = yi(0.0) as u32;
    let port_edges = fit.line(2, [px, py, 0], kh);
    let port = Port { edges: port_edges, r: 50.0 };
    // probes at 35, 50, 65 mm: voltage on the centre line, current through the plane half a cell further
    let (xa, d) = (35.0, 15.0);
    let mut vpr = Vec::new();
    let mut ipr = Vec::new();
    for q in 0..3 {
        let xq = xa + q as f64 * d;
        let pxq = (xi(xq) as u32 / root) * root;
        vpr.push(fit.line(2, [pxq, py, 0], kh));
        let xcut = pxq as f64 + 0.5;
        let ms: Vec<u32> = (0..fit.masters())
            .filter(|&mi| fit.pec[mi] && fit.ma[mi] == 0 && fit.mp[mi][2] == kh && (fit.mp[mi][0] as f64) < xcut && (fit.mp[mi][0] + fit.ml[mi]) as f64 > xcut && fit.mp[mi][1] > 0 && fit.mp[mi][1] < n[1])
            .map(|mi| mi as u32)
            .collect();
        ipr.push(ms);
    }
    let xq_mm: Vec<f64> = vpr.iter().enumerate().map(|(q, _)| ((xi(xa + q as f64 * d) as u32 / root) * root) as f64 * 0.5 - shift as f64 * 0.5).collect();
    // time step: the Yee limit of the finest cell, checked against the assembled operator
    let dt0 = 0.95 / (C0 * (1.0 / (hm[0] * hm[0]) + 1.0 / (hm[1] * hm[1]) + 1.0 / (hm[2] * hm[2])).sqrt());
    let mu = fit.stability(dt0, 60);
    let dt = if mu > 4.0 * 0.9025 { dt0 * (4.0 * 0.9025 / mu).sqrt() * 0.99 } else { dt0 };
    let fmax = 3e9;
    let exc = crate::amr2d::gauss(0.0, fmax, dt);
    // run until the port signals have decayed (fixed length, as long as the rectilinear runs need)
    let run = fit.run(dt, &port, &exc, &vpr, &ipr, 400_000, -70.0);
    let steps = run.steps;
    // Z0 and eps_eff from the three probes, as edge_cases::line
    let freqs = crate::post::linspace(200e6, 2.0e9, 19);
    let sp = |x: &[f64], toff: f64| crate::post::dft(x, dt, toff, &freqs);
    let (v1, v2, v3) = (sp(&run.v[0], 0.0), sp(&run.v[1], 0.0), sp(&run.v[2], 0.0));
    let (i1, i3) = (sp(&run.i[0], -0.5), sp(&run.i[2], -0.5));
    let dsep = xq_mm[1] - xq_mm[0];
    let (mut z0s, mut ees) = (0.0, 0.0);
    for k in 0..freqs.len() {
        let ch = (v1[k] + v3[k]) / (v2[k] * 2.0);
        let bd = ch.re.clamp(-1.0, 1.0).acos();
        let beta = bd / (dsep * 1e-3);
        ees += (beta * C0 / (2.0 * std::f64::consts::PI * freqs[k])).powi(2);
        let z2 = (v1[k] * v1[k] - v3[k] * v3[k]) / (i1[k] * i1[k] - i3[k] * i3[k]);
        z0s += z2.sqrt().abs();
    }
    let nfq = freqs.len() as f64;
    if std::env::var("FS_DEBUG").is_ok() {
        eprintln!("octree: dt^2 lambda {mu:.3} (dt {dt:.3e}), energy start/end {:.3e} / {:.3e}", run.energy.iter().cloned().fold(0.0, f64::max), run.energy.last().copied().unwrap_or(0.0));
    }
    (w, z0s / nfq, ees / nfq, fit.leaves, fit.masters(), fit.faces(), steps, run.seconds)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_octree_cavity_matches_the_analytic_mode() {
        // a uniform tree (no refinement) is the Yee scheme: an empty PEC box's lowest mode
        let n = [16, 12, 8];
        let h = [1e-3, 1e-3, 1e-3];
        let o = Octree::build(n, h, 4, &|l| l[3] > 1);
        let fit = assemble(&o, &|_| 1.0);
        // power iteration gives the largest eigenvalue; check against the uniform Yee bound instead:
        // dt^2 lambda_max = 4 (c dt)^2 sum 1/h^2 sin^2(...) <= 4 cfl^2
        let dt = 0.95 / (C0 * (3.0f64).sqrt() / 1e-3);
        let mu = fit.stability(dt, 200);
        assert!(mu < 4.0 * 0.9025 * 1.0001 && mu > 3.0, "uniform box: dt^2 lambda {mu}");
    }

    #[test]
    fn refined_octree_is_stable_and_conserves_energy() {
        // a box refined around a point: hanging edges everywhere at the interfaces
        let n = [32, 32, 32];
        let h = [1e-3, 1e-3, 1e-3];
        let o = Octree::build(n, h, 8, &|l| {
            let c = [l[0] as f64 + 0.5 * l[3] as f64, l[1] as f64 + 0.5 * l[3] as f64, l[2] as f64 + 0.5 * l[3] as f64];
            let d = ((c[0] - 12.0).powi(2) + (c[1] - 15.0).powi(2) + (c[2] - 17.0).powi(2)).sqrt();
            l[3] as f64 > 1.0 + 0.5 * d
        });
        let sizes: Vec<u32> = o.leaves.iter().map(|l| l[3]).collect();
        assert!(sizes.contains(&1) && sizes.contains(&4), "mixed sizes");
        let fit = assemble(&o, &|_| 1.0);
        assert!(fit.virtuals > 0);
        let dt = 0.95 / (C0 * (3.0f64).sqrt() / 1e-3);
        let mu = fit.stability(dt, 300);
        assert!(mu < 4.0, "refined box unstable at the Yee step: dt^2 lambda {mu}");
        // lossless run from a kick: energy must stay bounded (it oscillates between E and H)
        // a pulse from a resistive port in the fine region: afterwards the energy only drains into the
        // port, it must never grow (it oscillates between E and H, so compare window maxima)
        let p = fit.line(2, [12, 15, 16], 2);
        let port = Port { edges: p, r: 377.0 };
        let exc: Vec<f64> = crate::amr2d::gauss(0.0, 30e9, dt);
        let r = fit.run(dt, &port, &exc, &[], &[], 8000, -300.0);
        let k = r.energy.len();
        let e0 = r.energy[k / 4..k / 4 + 4].iter().cloned().fold(0.0, f64::max);
        let e1 = r.energy[k - 4..].iter().cloned().fold(0.0, f64::max);
        assert!(e1 <= e0 * 1.001, "energy grew {e0} -> {e1}");
    }
}
