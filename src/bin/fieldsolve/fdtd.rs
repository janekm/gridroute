//! FDTD time stepping: the multi-threaded CPU reference. The Metal backend (`metal_fdtd`) runs the same
//! updates and must agree with it.
//!
//! Yee layout on a non-uniform grid, index `(k * ny + j) * nx + i` for every component:
//! Ex(i,j,k) at (x_{i+1/2}, y_j, z_k), Hx(i,j,k) at (x_i, y_{j+1/2}, z_{k+1/2}), and cyclically.
//! E updates use `ca`/`cb` per edge (material, lumped loads, metal = 0); H uses dt/mu0.
//! Absorbing walls are CPML (kappa = 1) with recursive-convolution memory only in the boundary slabs.

use crate::model::{Problem, MU0, NONE};
use rayon::prelude::*;

#[derive(Clone, Debug, Default)]
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
    /// Running DFT of selected samples only (e.g. the faces of a near-to-far-field box): `gather[c]` lists
    /// the indices of component c (Ex, Ey, Ez, Hx, Hy, Hz) to accumulate at `gather_freqs`, every
    /// `dft_every` steps.
    pub gather_freqs: Vec<f64>,
    pub gather: [Vec<u32>; 6],
    /// Also stop when every port's signal energy over the last `check_every` steps has fallen `end_db`
    /// below the largest port's peak window energy (the S-parameters only need the port signals to have
    /// decayed; high-Q modes that do not reach the ports can keep the field energy up for much longer).
    pub port_stop: bool,
    /// With the port stop: also require these per-frequency quantities of the gathered DFT (e.g. the
    /// Poynting flux through a far-field box) to have settled, i.e. changed by less than `gather_tol`
    /// (relative) since the check `gather_lag` checks earlier.
    pub gather_metric: Option<GatherMetric>,
    pub gather_tol: f64,
    pub gather_lag: usize,
}

/// A per-frequency scalar of the gathered DFT, evaluated at the energy checks.
#[derive(Clone)]
pub struct GatherMetric(pub std::sync::Arc<dyn Fn(&[[Vec<f32>; 6]]) -> Vec<f64> + Send + Sync>);

impl std::fmt::Debug for GatherMetric {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("GatherMetric")
    }
}

/// Settling test on a history of gathered-DFT metrics (one Vec per check).
pub fn gather_settled(hist: &[Vec<f64>], lag: usize, tol: f64) -> bool {
    if hist.len() <= lag {
        return false;
    }
    let (now, then) = (&hist[hist.len() - 1], &hist[hist.len() - 1 - lag]);
    now.iter().zip(then).all(|(a, b)| (a - b).abs() <= tol * a.abs().max(b.abs()).max(1e-300))
}

/// Port-signal decay test over the last `win` samples of each port: true when every port's window
/// energy is at most `10^(end_db/10)` times the largest window energy seen on any port (`peak`, updated).
pub fn ports_decayed(v: &[Vec<f32>], win: usize, end_db: f64, peak: &mut f64) -> bool {
    let w: Vec<f64> = v.iter().map(|s| s[s.len().saturating_sub(win)..].iter().map(|&x| (x as f64) * (x as f64)).sum()).collect();
    let m = w.iter().cloned().fold(0.0, f64::max);
    *peak = peak.max(m);
    *peak > 0.0 && m <= *peak * 10f64.powf(end_db / 10.0)
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
    /// Gathered DFT: per `gather_freqs` entry, per component, interleaved (re, im) in `gather` order.
    pub gathered: Vec<[Vec<f32>; 6]>,
    /// field energy (arbitrary units) at each energy check: (step, energy)
    pub energy: Vec<(usize, f64)>,
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
    /// Debye currents (compact, `np` per edge) and E before the update over the dispersive planes
    dj: [Vec<f32>; 3],
    e0: [Vec<f32>; 3],
    /// surface-impedance section currents, `nb` per conductor edge
    ci: Vec<f32>,
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
        let (nd, np) = p.disp.as_ref().map(|d| (d.nd(p.nx, p.ny), d.np)).unwrap_or((0, 0));
        CpuState { e: [z(), z(), z()], h: [z(), z(), z()], psi, dj: std::array::from_fn(|_| vec![0.0f32; nd * np]), e0: std::array::from_fn(|_| vec![0.0f32; nd]), ci: vec![0.0f32; p.cond.as_ref().map(|c| c.ids.len() * c.nb).unwrap_or(0)] }
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

/// Copper with a surface impedance (cond.rs): on every listed sheet edge, E = Z_s * J_s with J_s the
/// sheet current from the curl of H across the sheet (after the E update, which left these edges at 0).
pub fn cond_update(p: &Problem, s: &mut CpuState) {
    let Some(c) = &p.cond else { return };
    let (nx, ny) = (p.nx, p.ny);
    let (sy, sz) = (nx, nx * ny);
    let [hx, hy, hz] = &s.h;
    let [idx, idy, idz] = &p.inv_dd;
    let nb = c.nb;
    let ex = Ptr(s.e[0].as_mut_ptr());
    let ey = Ptr(s.e[1].as_mut_ptr());
    s.ci.par_chunks_mut(nb.max(1)).enumerate().for_each(|(t, ik)| {
        let w = c.ids[t];
        let (a, id) = ((w >> 30) as usize, (w & 0x3fff_ffff) as usize);
        let (i, j, k) = (id % nx, (id / nx) % ny, id / sz);
        let curl = if a == 0 {
            (hz[id] - hz[id - sy]) * idy[j] - (hy[id] - hy[id - sz]) * idz[k]
        } else {
            (hx[id] - hx[id - sz]) * idz[k] - (hz[id] - hz[id - 1]) * idx[i]
        };
        let js = curl * c.fac[t];
        let r = &c.tab[c.layer[t] as usize * (nb + 1)..];
        let mut e = r[0] * js;
        for q in 0..nb {
            ik[q] = js + (ik[q] - js) * c.dec[q];
            e += r[q + 1] * (js - ik[q]);
        }
        unsafe {
            *(if a == 0 { ex } else { ey }).at(id) = e;
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

/// Before the E update: keep E over the dispersive planes.
pub fn disp_save(p: &Problem, s: &mut CpuState) {
    if let Some(d) = &p.disp {
        let (a, b) = (d.k0 * p.nx * p.ny, (d.k1 + 1) * p.nx * p.ny);
        for c in 0..3 {
            s.e0[c].copy_from_slice(&s.e[c][a..b]);
        }
    }
}

/// After the E update and its CPML correction: E -= cb sum (1 + k_p)/2 J_p, then
/// J_p = k_p J_p + b_p (E_new - E_old) (disp.rs).
pub fn disp_update(p: &Problem, s: &mut CpuState) {
    let Some(d) = &p.disp else { return };
    let base = d.k0 * p.nx * p.ny;
    let np = d.np;
    for c in 0..3 {
        let e = Ptr(s.e[c].as_mut_ptr());
        let cb = &p.cb[c];
        let b = &d.b[c];
        let e0 = &s.e0[c];
        s.dj[c].par_chunks_mut(np).enumerate().for_each(|(q, jj)| {
            let id = base + q;
            let k = cb[id];
            if k == 0.0 {
                return;
            }
            let bq = &b[q * np..q * np + np];
            if bq.iter().all(|&x| x == 0.0) {
                return;
            }
            let mut sum = 0.0f32;
            for pp in 0..np {
                sum += d.hk[pp] * jj[pp];
            }
            unsafe {
                let e1 = *e.at(id) - k * sum;
                *e.at(id) = e1;
                let de = e1 - e0[q];
                for pp in 0..np {
                    jj[pp] = d.kp[pp] * jj[pp] + bq[pp] * de;
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

/// DFT weights (re, im) of component c at step `n` for angular frequency `w`: E at n dt, H at (n - 1/2) dt.
pub fn dft_weight(p: &Problem, c: usize, n: usize, every: usize, w: f64) -> (f32, f32) {
    let t = if c < 3 { n as f64 * p.dt } else { (n as f64 - 0.5) * p.dt };
    let (sn, cs) = (-w * t).sin_cos();
    ((cs * p.dt * every as f64) as f32, (sn * p.dt * every as f64) as f32)
}

fn gather_acc(p: &Problem, s: &CpuState, n: usize, spec: &RunSpec, acc: &mut [[Vec<f32>; 6]]) {
    for (fi, &f) in spec.gather_freqs.iter().enumerate() {
        let w = 2.0 * std::f64::consts::PI * f;
        for c in 0..6 {
            let (wr, wi) = dft_weight(p, c, n, spec.dft_every, w);
            let src = if c < 3 { &s.e[c] } else { &s.h[c - 3] };
            acc[fi][c].par_chunks_mut(2).zip(spec.gather[c].par_iter()).for_each(|(o, &id)| {
                let v = src[id as usize];
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
    let mut gathered: Vec<[Vec<f32>; 6]> = spec.gather_freqs.iter().map(|_| std::array::from_fn(|c| vec![0.0f32; 2 * spec.gather[c].len()])).collect();
    let mut peak = 0.0f64;
    let mut port_peak = 0.0f64;
    let mut ghist: Vec<Vec<f64>> = Vec::new();
    let mut ehist: Vec<(usize, f64)> = Vec::new();
    let mut steps = 0;
    let mut snaps = Vec::new();
    for n in 0..spec.max_steps {
        update_h(p, &mut s);
        pml_h(p, &mut s);
        disp_save(p, &mut s);
        update_e(p, &mut s);
        pml_e(p, &mut s);
        // the source is part of this step's E change: inject it before the Debye currents see E_new - E_old
        // (with dispersion the order matters; injecting after left the port edges' polarisation without the
        // source's share and made the S-matrix non-passive)
        if let Some(pi) = spec.excite {
            if n < spec.exc.len() {
                let vs = spec.exc[n];
                for &(c, id, g) in &p.ports[pi].src {
                    s.e[c as usize][id as usize] += g * vs;
                }
            }
        }
        disp_update(p, &mut s);
        cond_update(p, &mut s);
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
        if !gathered.is_empty() && steps % spec.dft_every == 0 {
            gather_acc(p, &s, steps, spec, &mut gathered);
        }
        if steps % spec.check_every == 0 {
            let en = energy(&s);
            ehist.push((steps, en));
            peak = peak.max(en);
            // a configured gathered-field metric (e.g. the far field) gates every convergence exit, energy-
            // and port-based alike; only max_steps ends a run without it
            let gs = match &spec.gather_metric {
                Some(m) => {
                    ghist.push((m.0)(&gathered));
                    gather_settled(&ghist, spec.gather_lag, spec.gather_tol)
                }
                None => true,
            };
            if gs && steps > spec.exc.len() && en <= peak * 10f64.powf(spec.end_db / 10.0) {
                break;
            }
            if spec.port_stop {
                let pd = ports_decayed(&v, spec.check_every, spec.end_db, &mut port_peak) && steps > spec.exc.len();
                if pd && gs {
                    break;
                }
            }
        }
    }
    RunResult { steps, v, probes: pv, fields: fd, gathered, energy: ehist, snaps, seconds: t0.elapsed().as_secs_f64(), backend: "cpu".into() }
}

#[allow(unused)]
pub fn none_used() -> u32 {
    NONE
}

/// Late-time decay of the field energy over the last third of the checks: (dB per ns, Q at `f`).
/// A positive rate means the energy grows (an instability, not ringing).
pub fn energy_decay(hist: &[(usize, f64)], dt: f64, f: f64) -> Option<(f64, f64)> {
    let n = hist.len();
    if n < 6 {
        return None;
    }
    let tail = &hist[n - n / 3..];
    let (mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0);
    for &(s, e) in tail {
        let x = s as f64 * dt * 1e9;
        let y = 10.0 * e.max(1e-300).log10();
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    let m = tail.len() as f64;
    let slope = (m * sxy - sx * sy) / (m * sxx - sx * sx);
    // energy ~ exp(-w t / Q): dB/ns = -10 log10(e) w / Q * 1e-9
    let q = -10.0 * std::f64::consts::LOG10_E * 2.0 * std::f64::consts::PI * f * 1e-9 / slope;
    Some((slope, q))
}
