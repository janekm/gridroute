//! FDTD time stepping: the multi-threaded CPU reference. The Metal backend (`metal_fdtd`) runs the same
//! updates and must agree with it.
//!
//! Yee layout on a non-uniform grid, index `(k * ny + j) * nx + i` for every component:
//! Ex(i,j,k) at (x_{i+1/2}, y_j, z_k), Hx(i,j,k) at (x_i, y_{j+1/2}, z_{k+1/2}), and cyclically.
//! E updates use `ca`/`cb` per edge (material, lumped loads, metal = 0); H uses dt/mu0.
//! Absorbing walls are CPML (kappa = 1) with recursive-convolution memory only in the boundary slabs.

use crate::model::{Problem, MU0, NONE};
use rayon::prelude::*;

#[derive(Clone, Debug)]
pub struct RunSpec {
    /// Port excited (index into `Problem::ports`), or none (e.g. a pure ring-down).
    pub excite: Option<usize>,
    /// Source voltage per step, sampled at (n + 1/2) dt.
    pub exc: Vec<f32>,
    pub max_steps: usize,
    /// Stop when the field energy falls this far (dB) below its peak, after the excitation.
    pub end_db: f64,
    pub check_every: usize,
    /// Frequencies for field snapshots (running DFT of all six components), sampled every `dft_every` steps.
    pub dft_freqs: Vec<f64>,
    pub dft_every: usize,
    /// Time-domain snapshots of the surface current on mesh planes `snap_k`, every `snap_every` steps
    /// (0 = none).
    pub snap_every: usize,
    pub snap_k: Vec<usize>,
    /// No snapshots after this step (0 = no limit).
    pub snap_until: usize,
}

/// One time-domain snapshot: per requested plane, the sheet current (Jx, Jy) at the nodes, A/m.
#[derive(Clone, Debug)]
pub struct Snap {
    pub step: usize,
    pub j: Vec<[Vec<f32>; 2]>,
}

/// Surface current n x (H+ - H-) on mesh plane k from the H arrays: Jx = Hy(k) - Hy(k-1),
/// Jy = Hx(k-1) - Hx(k), sampled at the nodes (the half-cell offsets of Hx/Hy are ignored).
pub fn sheet_current(nx: usize, ny: usize, nz: usize, hx: &[f32], hy: &[f32], k: usize) -> [Vec<f32>; 2] {
    let mut jx = vec![0.0f32; nx * ny];
    let mut jy = vec![0.0f32; nx * ny];
    if k == 0 || k + 1 >= nz {
        return [jx, jy];
    }
    for j in 0..ny {
        for i in 0..nx {
            let up = (k * ny + j) * nx + i;
            let dn = ((k - 1) * ny + j) * nx + i;
            jx[j * nx + i] = hy[up] - hy[dn];
            jy[j * nx + i] = hx[dn] - hx[up];
        }
    }
    [jx, jy]
}

/// Frequency-domain fields: per frequency, per component (Ex..Ez, Hx..Hz), interleaved (re, im).
#[derive(Clone, Debug, Default)]
pub struct FieldDft {
    pub freqs: Vec<f64>,
    pub data: Vec<[Vec<f32>; 6]>,
}

#[derive(Clone, Debug)]
pub struct RunResult {
    pub steps: usize,
    /// Port voltages at integer steps 0..=steps, per port.
    pub v: Vec<Vec<f32>>,
    /// Field probe samples at steps 0..=steps (E probes at n dt, H probes at (n - 1/2) dt).
    pub probes: Vec<Vec<f32>>,
    pub fields: Option<FieldDft>,
    pub snaps: Vec<Snap>,
    pub seconds: f64,
    pub backend: String,
}

/// Gaussian pulse as openEMS's `SetGaussExcite(f0, fc)`: centre `f0`, 20 dB half bandwidth `fc`.
pub fn gauss_excite(f0: f64, fc: f64, dt: f64) -> Vec<f32> {
    let t0 = 9.0 / (2.0 * std::f64::consts::PI * fc);
    let tau = 3.0 / (2.0 * std::f64::consts::PI * fc);
    let n = (2.0 * t0 / dt).ceil() as usize + 1;
    (0..n)
        .map(|s| {
            let t = (s as f64 + 0.5) * dt - t0;
            ((-(t / tau).powi(2)).exp() * (2.0 * std::f64::consts::PI * f0 * t).cos()) as f32
        })
        .collect()
}

#[derive(Clone, Copy)]
struct Ptr(*mut f32);
unsafe impl Send for Ptr {}
unsafe impl Sync for Ptr {}
impl Ptr {
    #[inline(always)]
    unsafe fn at(self, i: usize) -> *mut f32 {
        self.0.add(i)
    }
}

pub struct CpuState {
    pub e: [Vec<f32>; 3],
    pub h: [Vec<f32>; 3],
    /// CPML memory per axis a: [E comp (a+1)%3, E comp (a+2)%3, H comp (a+1)%3, H comp (a+2)%3].
    psi: [[Vec<f32>; 4]; 3],
}

impl CpuState {
    pub fn new(p: &Problem) -> CpuState {
        let n = p.n();
        let z = || vec![0.0f32; n];
        let nn = [p.nx, p.ny, p.nz];
        let psi = std::array::from_fn(|a| {
            let len = p.pml[a].len() * n / nn[a];
            std::array::from_fn(|_| vec![0.0f32; len])
        });
        CpuState { e: [z(), z(), z()], h: [z(), z(), z()], psi }
    }
}

/// Main H update.
pub fn update_h(p: &Problem, s: &mut CpuState) {
    let (nx, ny, nz) = (p.nx, p.ny, p.nz);
    let k = (p.dt / MU0) as f32;
    let [ex, ey, ez] = &s.e;
    let (hx, hy, hz) = (Ptr(s.h[0].as_mut_ptr()), Ptr(s.h[1].as_mut_ptr()), Ptr(s.h[2].as_mut_ptr()));
    let [idx, idy, idz] = &p.inv_d;
    (0..nz * ny).into_par_iter().for_each(|r| {
        let (j, kk) = (r % ny, r / ny);
        let base = r * nx;
        let jn = j + 1 < ny;
        let kn = kk + 1 < nz;
        let sy = nx;
        let sz = nx * ny;
        for i in 0..nx {
            let id = base + i;
            let inn = i + 1 < nx;
            unsafe {
                if jn && kn {
                    *hx.at(id) -= k * ((ez[id + sy] - ez[id]) * idy[j] - (ey[id + sz] - ey[id]) * idz[kk]);
                }
                if inn && kn {
                    *hy.at(id) -= k * ((ex[id + sz] - ex[id]) * idz[kk] - (ez[id + 1] - ez[id]) * idx[i]);
                }
                if inn && jn {
                    *hz.at(id) -= k * ((ey[id + 1] - ey[id]) * idx[i] - (ex[id + sy] - ex[id]) * idy[j]);
                }
            }
        }
    });
}

/// Main E update (boundary and metal edges have ca = cb = 0 and stay zero).
pub fn update_e(p: &Problem, s: &mut CpuState) {
    let (nx, ny, nz) = (p.nx, p.ny, p.nz);
    let [hx, hy, hz] = &s.h;
    let (ex, ey, ez) = (Ptr(s.e[0].as_mut_ptr()), Ptr(s.e[1].as_mut_ptr()), Ptr(s.e[2].as_mut_ptr()));
    let [idx, idy, idz] = &p.inv_dd;
    let [cax, cay, caz] = &p.ca;
    let [cbx, cby, cbz] = &p.cb;
    (0..nz * ny).into_par_iter().for_each(|r| {
        let (j, kk) = (r % ny, r / ny);
        if j == 0 || kk == 0 || j + 1 >= ny || kk + 1 >= nz {
            // x-walls handled per i below; y/z walls: Ex, Ey/Ez on these rows are tangential or absent,
            // but Ez on a z-wall row (kk == 0) and Ey on a y-wall row exist only as normals: their cb is 0
            // whenever they lie on another wall, so the general path below is still correct for them
        }
        let base = r * nx;
        let sy = nx;
        let sz = nx * ny;
        for i in 0..nx {
            let id = base + i;
            unsafe {
                // the H neighbours below index i-1, j-1, k-1: an edge with any such index at 0 sits on a
                // wall and has cb = 0, so it is skipped
                if cbx[id] != 0.0 && j > 0 && kk > 0 {
                    *ex.at(id) = cax[id] * *ex.at(id) + cbx[id] * ((hz[id] - hz[id - sy]) * idy[j] - (hy[id] - hy[id - sz]) * idz[kk]);
                }
                if cby[id] != 0.0 && i > 0 && kk > 0 {
                    *ey.at(id) = cay[id] * *ey.at(id) + cby[id] * ((hx[id] - hx[id - sz]) * idz[kk] - (hz[id] - hz[id - 1]) * idx[i]);
                }
                if cbz[id] != 0.0 && i > 0 && j > 0 {
                    *ez.at(id) = caz[id] * *ez.at(id) + cbz[id] * ((hy[id] - hy[id - 1]) * idx[i] - (hx[id] - hx[id - sy]) * idy[j]);
                }
            }
        }
    });
}

/// Visit every cell of axis `a`'s CPML slabs: `f(slab_index, s, id, pos_along_axis, psi_index)`.
fn for_slab(p: &Problem, a: usize, f: impl Fn(usize, usize, usize) + Sync + Send) {
    let (nx, ny, nz) = (p.nx, p.ny, p.nz);
    let ax = &p.pml[a];
    let l = ax.len();
    if l == 0 {
        return;
    }
    match a {
        0 => (0..nz * ny).into_par_iter().for_each(|r| {
            for (s, &i) in ax.slab.iter().enumerate() {
                f(s, r * nx + i as usize, r * l + s);
            }
        }),
        1 => (0..nz * l).into_par_iter().for_each(|r| {
            let (s, k) = (r % l, r / l);
            let j = ax.slab[s] as usize;
            for i in 0..nx {
                f(s, (k * ny + j) * nx + i, r * nx + i);
            }
        }),
        _ => (0..l * ny).into_par_iter().for_each(|r| {
            let (j, s) = (r % ny, r / ny);
            let k = ax.slab[s] as usize;
            for i in 0..nx {
                f(s, (k * ny + j) * nx + i, r * nx + i);
            }
        }),
    }
}

/// CPML correction of H after the main H update.
pub fn pml_h(p: &Problem, s: &mut CpuState) {
    let kmu = (p.dt / MU0) as f32;
    let nn = [p.nx, p.ny, p.nz];
    let st = [1usize, p.nx, p.nx * p.ny];
    for a in 0..3 {
        if p.pml[a].len() == 0 {
            continue;
        }
        let (b, c) = ((a + 1) % 3, (a + 2) % 3);
        let ax = &p.pml[a];
        let e = [Ptr(s.e[0].as_mut_ptr()), Ptr(s.e[1].as_mut_ptr()), Ptr(s.e[2].as_mut_ptr())];
        let h = [Ptr(s.h[0].as_mut_ptr()), Ptr(s.h[1].as_mut_ptr()), Ptr(s.h[2].as_mut_ptr())];
        let ps = [Ptr(s.psi[a][2].as_mut_ptr()), Ptr(s.psi[a][3].as_mut_ptr())];
        let invd = &p.inv_d[a];
        let (nb, nc) = (nn[b], nn[c]);
        for_slab(p, a, |sl, id, q| {
            let pos = ax.slab[sl] as usize;
            if pos + 1 >= nn[a] {
                return;
            }
            let coord = [id % p.nx, (id / p.nx) % p.ny, id / (p.nx * p.ny)];
            let (bh, ch) = (ax.bh[sl], ax.ch[sl]);
            unsafe {
                // H_b has +dE_c/da ... sign pattern: H_b -= kmu * (dE_a/dc - dE_c/da); H_c -= kmu * (dE_b/da - dE_a/db)
                // H_b (comp b) uses dE_c/da: valid when coord[c]+1 < nc
                if coord[c] + 1 < nc {
                    let der = (*e[c].at(id + st[a]) - *e[c].at(id)) * invd[pos];
                    let pp = ps[0].at(q);
                    *pp = bh * *pp + ch * der;
                    *h[b].at(id) += kmu * *pp;
                }
                if coord[b] + 1 < nb {
                    let der = (*e[b].at(id + st[a]) - *e[b].at(id)) * invd[pos];
                    let pp = ps[1].at(q);
                    *pp = bh * *pp + ch * der;
                    *h[c].at(id) -= kmu * *pp;
                }
            }
        });
    }
}

/// CPML correction of E after the main E update.
pub fn pml_e(p: &Problem, s: &mut CpuState) {
    let st = [1usize, p.nx, p.nx * p.ny];
    for a in 0..3 {
        if p.pml[a].len() == 0 {
            continue;
        }
        let (b, c) = ((a + 1) % 3, (a + 2) % 3);
        let ax = &p.pml[a];
        let e = [Ptr(s.e[0].as_mut_ptr()), Ptr(s.e[1].as_mut_ptr()), Ptr(s.e[2].as_mut_ptr())];
        let h = [Ptr(s.h[0].as_mut_ptr()), Ptr(s.h[1].as_mut_ptr()), Ptr(s.h[2].as_mut_ptr())];
        let ps = [Ptr(s.psi[a][0].as_mut_ptr()), Ptr(s.psi[a][1].as_mut_ptr())];
        let invdd = &p.inv_dd[a];
        let cb = &p.cb;
        for_slab(p, a, |sl, id, q| {
            let pos = ax.slab[sl] as usize;
            if pos == 0 {
                return;
            }
            let (be, ce) = (ax.be[sl], ax.ce[sl]);
            unsafe {
                // E_b -= cb * dH_c/da ; E_c += cb * dH_b/da
                if cb[b][id] != 0.0 {
                    let der = (*h[c].at(id) - *h[c].at(id - st[a])) * invdd[pos];
                    let pp = ps[0].at(q);
                    *pp = be * *pp + ce * der;
                    *e[b].at(id) -= cb[b][id] * *pp;
                }
                if cb[c][id] != 0.0 {
                    let der = (*h[b].at(id) - *h[b].at(id - st[a])) * invdd[pos];
                    let pp = ps[1].at(q);
                    *pp = be * *pp + ce * der;
                    *e[c].at(id) += cb[c][id] * *pp;
                }
            }
        });
    }
}

pub fn energy(s: &CpuState) -> f64 {
    let eta2 = (crate::model::ETA0 * crate::model::ETA0) as f32;
    let mut t = 0.0f64;
    for a in 0..3 {
        t += s.e[a].par_iter().map(|&v| (v * v) as f64).sum::<f64>();
        t += s.h[a].par_iter().map(|&v| (v * v * eta2) as f64).sum::<f64>();
    }
    t
}

fn probe(p: &Problem, s: &CpuState) -> Vec<f32> {
    p.ports.iter().map(|pt| pt.probe.iter().map(|&(c, id, w)| w * s.e[c as usize][id as usize]).sum()).collect()
}

fn field_probe(p: &Problem, s: &CpuState) -> Vec<f32> {
    p.probes
        .iter()
        .map(|pr| pr.terms.iter().map(|&(c, id, w)| w * if c < 3 { s.e[c as usize][id as usize] } else { s.h[c as usize - 3][id as usize] }).sum())
        .collect()
}

/// Running DFT accumulation of all components at step `n` (E at n dt, H at (n - 1/2) dt).
fn dft_acc(p: &Problem, s: &CpuState, n: usize, every: usize, fd: &mut FieldDft) {
    for (fi, &f) in fd.freqs.iter().enumerate() {
        let w = 2.0 * std::f64::consts::PI * f;
        for c in 0..6 {
            let t = if c < 3 { n as f64 * p.dt } else { (n as f64 - 0.5) * p.dt };
            let (sn, cs) = (-w * t).sin_cos();
            let (wr, wi) = ((cs * p.dt * every as f64) as f32, (sn * p.dt * every as f64) as f32);
            let src = if c < 3 { &s.e[c] } else { &s.h[c - 3] };
            fd.data[fi][c].par_chunks_mut(2).zip(src.par_iter()).for_each(|(o, &v)| {
                o[0] += v * wr;
                o[1] += v * wi;
            });
        }
    }
}

pub fn run_cpu(p: &Problem, spec: &RunSpec) -> RunResult {
    let t0 = std::time::Instant::now();
    let mut s = CpuState::new(p);
    let mut v: Vec<Vec<f32>> = vec![vec![0.0]; p.ports.len()];
    let mut pv: Vec<Vec<f32>> = vec![vec![0.0]; p.probes.len()];
    let mut fd = (!spec.dft_freqs.is_empty()).then(|| FieldDft {
        freqs: spec.dft_freqs.clone(),
        data: spec.dft_freqs.iter().map(|_| std::array::from_fn(|_| vec![0.0f32; 2 * p.n()])).collect(),
    });
    let mut peak = 0.0f64;
    let mut steps = 0;
    let mut snaps = Vec::new();
    for n in 0..spec.max_steps {
        update_h(p, &mut s);
        pml_h(p, &mut s);
        update_e(p, &mut s);
        pml_e(p, &mut s);
        if let Some(pi) = spec.excite {
            if n < spec.exc.len() {
                let vs = spec.exc[n];
                for &(c, id, g) in &p.ports[pi].src {
                    s.e[c as usize][id as usize] += g * vs;
                }
            }
        }
        for (k, x) in probe(p, &s).into_iter().enumerate() {
            v[k].push(x);
        }
        for (k, x) in field_probe(p, &s).into_iter().enumerate() {
            pv[k].push(x);
        }
        steps = n + 1;
        if spec.snap_every > 0 && steps % spec.snap_every == 0 && (spec.snap_until == 0 || steps <= spec.snap_until) {
            snaps.push(Snap { step: steps, j: spec.snap_k.iter().map(|&k| sheet_current(p.nx, p.ny, p.nz, &s.h[0], &s.h[1], k)).collect() });
        }
        if let Some(fd) = fd.as_mut() {
            if steps % spec.dft_every == 0 {
                dft_acc(p, &s, steps, spec.dft_every, fd);
            }
        }
        if steps % spec.check_every == 0 {
            let en = energy(&s);
            peak = peak.max(en);
            if steps > spec.exc.len() && en <= peak * 10f64.powf(spec.end_db / 10.0) {
                break;
            }
        }
    }
    RunResult { steps, v, probes: pv, fields: fd, snaps, seconds: t0.elapsed().as_secs_f64(), backend: "cpu".into() }
}

#[allow(unused)]
pub fn none_used() -> u32 {
    NONE
}
