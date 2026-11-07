//! Near-to-far-field transform (Huygens box).
//!
//! A closed box of mesh planes in the air around the structure records the tangential E and H on its six
//! faces as running DFTs (the engines' gathered DFT). On each face cell the equivalent currents are
//! J = n x H and M = -n x E (n outward). With the e^{jwt} convention of the field DFT, and r^ the
//! observation direction (Balanis, Advanced Engineering Electromagnetics, 6.8):
//!   N = sum J e^{+jk r^.r'} dS,   L = sum M e^{+jk r^.r'} dS,
//!   r E_theta e^{jkr} = -jk/(4 pi) (L_phi + eta N_theta),   r E_phi e^{jkr} = +jk/(4 pi) (L_theta - eta N_phi),
//!   U = (|r E_theta|^2 + |r E_phi|^2) / (2 eta)   [W/sr for peak phasors].
//! The box must lie in homogeneous free space (air) inside the absorbing boundary: everything that
//! radiates must be inside it, and nothing may cross its faces.
//!
//! Yee positions (fdtd.rs): E_a at a half index along a, nodes along the other two axes; H_a at a node
//! along a, half indices along the other two. On a face normal to axis a at node plane p, a face cell
//! (m, n) along the tangential axes b = a+1, c = a+2 has its centre at (b_{m+1/2}, c_{n+1/2}, a_p):
//!   E_b: average over c-nodes n, n+1 on plane p;   E_c: average over b-nodes m, m+1 on plane p;
//!   H_b: average over b-nodes m, m+1 and the half planes p-1, p (a_{p-1/2}, a_{p+1/2});
//!   H_c: average over c-nodes n, n+1 and the half planes p-1, p.
//! The DFT already places E at n dt and H at (n - 1/2) dt, so the phases are consistent.

use crate::mesh::Grid;
use crate::model::{C0, ETA0};
use crate::post::Cx;
use rayon::prelude::*;
use std::f64::consts::PI;
use std::path::Path;

/// One rectangular patch of samples of a single component on one mesh plane.
#[derive(Clone, Debug)]
struct Patch {
    comp: usize,
    /// normal axis and node plane index
    a: usize,
    plane: usize,
    /// inclusive index ranges along the tangential axes b = (a+1)%3, c = (a+2)%3
    b: (usize, usize),
    c: (usize, usize),
    /// position of the first sample in the component's gather list
    off: usize,
}

#[derive(Clone, Debug)]
pub struct Plan {
    /// node index bounds of the box (inclusive)
    pub lo: [usize; 3],
    pub hi: [usize; 3],
    pub gather: [Vec<u32>; 6],
    /// per face (axis a, side -1/+1): patches for E_b, E_c, H_b (p-1, p), H_c (p-1, p)
    faces: Vec<(usize, f64, [usize; 6])>,
    patches: Vec<Patch>,
    dims: [usize; 3],
}

impl Plan {
    /// Box with node index bounds `lo..=hi` on a grid of `dims` = [nx, ny, nz] nodes.
    pub fn new(dims: [usize; 3], lo: [usize; 3], hi: [usize; 3]) -> Result<Plan, String> {
        for a in 0..3 {
            if lo[a] < 1 || hi[a] + 1 >= dims[a] || hi[a] <= lo[a] + 1 {
                return Err(format!("far-field box axis {a}: nodes {}..{} do not fit inside 1..{}", lo[a], hi[a], dims[a].saturating_sub(2)));
            }
        }
        let mut plan = Plan { lo, hi, gather: Default::default(), faces: Vec::new(), patches: Vec::new(), dims };
        for a in 0..3 {
            let (b, c) = ((a + 1) % 3, (a + 2) % 3);
            for (side, p) in [(-1.0, lo[a]), (1.0, hi[a])] {
                let full_b = (lo[b], hi[b]);
                let full_c = (lo[c], hi[c]);
                let cell_b = (lo[b], hi[b] - 1);
                let cell_c = (lo[c], hi[c] - 1);
                // E_b at plane p: b cells, c nodes; E_c: b nodes, c cells
                let eb = plan.add(b, a, p, cell_b, full_c);
                let ec = plan.add(c, a, p, full_b, cell_c);
                // H_b: node along b, half along c and a; H_c: node along c, half along a and b
                let hb0 = plan.add(3 + b, a, p - 1, full_b, cell_c);
                let hb1 = plan.add(3 + b, a, p, full_b, cell_c);
                let hc0 = plan.add(3 + c, a, p - 1, cell_b, full_c);
                let hc1 = plan.add(3 + c, a, p, cell_b, full_c);
                plan.faces.push((a, side, [eb, ec, hb0, hb1, hc0, hc1]));
            }
        }
        Ok(plan)
    }

    fn add(&mut self, comp: usize, a: usize, plane: usize, b: (usize, usize), c: (usize, usize)) -> usize {
        let off = self.gather[comp].len();
        let (bb, cc) = ((a + 1) % 3, (a + 2) % 3);
        let [nx, ny, _] = self.dims;
        for jc in c.0..=c.1 {
            for jb in b.0..=b.1 {
                let mut ijk = [0usize; 3];
                ijk[a] = plane;
                ijk[bb] = jb;
                ijk[cc] = jc;
                self.gather[comp].push(((ijk[2] * ny + ijk[1]) * nx + ijk[0]) as u32);
            }
        }
        self.patches.push(Patch { comp, a, plane, b, c, off });
        self.patches.len() - 1
    }

    /// Box from physical bounds (m): the node planes nearest to them.
    pub fn from_bounds(g: &Grid, lo_m: [f64; 3], hi_m: [f64; 3]) -> Result<Plan, String> {
        let snap = |l: &[f64], v: f64| l.iter().enumerate().min_by(|a, b| (a.1 - v).abs().partial_cmp(&(b.1 - v).abs()).unwrap()).unwrap().0;
        let ax = [&g.x, &g.y, &g.z];
        let lo = std::array::from_fn(|a| snap(ax[a], lo_m[a]));
        let hi = std::array::from_fn(|a| snap(ax[a], hi_m[a]));
        Plan::new([g.nx(), g.ny(), g.nz()], lo, hi)
    }

    /// Edges (component, index) of tangential E on the box faces, for a material check: the faces must
    /// lie in free space.
    pub fn face_e_edges(&self) -> Vec<(usize, u32)> {
        let mut v = Vec::new();
        for &(_, _, ids) in &self.faces {
            for &pi in &ids[..2] {
                let pt = &self.patches[pi];
                let n = (pt.b.1 - pt.b.0 + 1) * (pt.c.1 - pt.c.0 + 1);
                v.extend(self.gather[pt.comp][pt.off..pt.off + n].iter().map(|&id| (pt.comp, id)));
            }
        }
        v
    }
}

/// Equivalent sources on the face cells at one frequency.
pub struct Sources {
    pub pos: Vec<[f64; 3]>,
    pub j: Vec<[Cx; 3]>,
    pub m: Vec<[Cx; 3]>,
    pub ds: Vec<f64>,
    /// 1/2 Re of the outward Poynting flux through the box (W, peak phasors)
    pub p_box: f64,
}

/// Face currents from one frequency's gathered DFT (`data[c]` interleaved re/im in gather order).
pub fn sources(plan: &Plan, g: &Grid, data: &[Vec<f32>; 6]) -> Sources {
    let ax = [&g.x, &g.y, &g.z];
    let centre: [f64; 3] = std::array::from_fn(|a| 0.5 * (ax[a][plan.lo[a]] + ax[a][plan.hi[a]]));
    let at = |pi: usize, jb: usize, jc: usize| -> Cx {
        let pt = &plan.patches[pi];
        let nb = pt.b.1 - pt.b.0 + 1;
        let k = pt.off + (jc - pt.c.0) * nb + (jb - pt.b.0);
        let d = &data[pt.comp];
        Cx::new(d[2 * k] as f64, d[2 * k + 1] as f64)
    };
    let mut s = Sources { pos: Vec::new(), j: Vec::new(), m: Vec::new(), ds: Vec::new(), p_box: 0.0 };
    for &(a, side, [eb, ec, hb0, hb1, hc0, hc1]) in &plan.faces {
        let (b, c) = ((a + 1) % 3, (a + 2) % 3);
        let p = plan.patches[eb].plane;
        for n in plan.lo[c]..plan.hi[c] {
            for m in plan.lo[b]..plan.hi[b] {
                let e_b = (at(eb, m, n) + at(eb, m, n + 1)) * 0.5;
                let e_c = (at(ec, m, n) + at(ec, m + 1, n)) * 0.5;
                let h_b = (at(hb0, m, n) + at(hb0, m + 1, n) + at(hb1, m, n) + at(hb1, m + 1, n)) * 0.25;
                let h_c = (at(hc0, m, n) + at(hc0, m, n + 1) + at(hc1, m, n) + at(hc1, m, n + 1)) * 0.25;
                let db = ax[b][m + 1] - ax[b][m];
                let dc = ax[c][n + 1] - ax[c][n];
                let mut r = [0.0; 3];
                r[a] = ax[a][p] - centre[a];
                r[b] = 0.5 * (ax[b][m] + ax[b][m + 1]) - centre[b];
                r[c] = 0.5 * (ax[c][n] + ax[c][n + 1]) - centre[c];
                // n = side e_a:  J = n x H -> J_b = -side H_c, J_c = side H_b;  M = -n x E -> M_b = side E_c, M_c = -side E_b
                let mut j = [Cx::default(); 3];
                let mut mm = [Cx::default(); 3];
                j[b] = h_c * (-side);
                j[c] = h_b * side;
                mm[b] = e_c * side;
                mm[c] = e_b * (-side);
                // (E x H*) . n = side (E_b H_c* - E_c H_b*)
                s.p_box += 0.5 * side * (e_b * h_c.conj() - e_c * h_b.conj()).re * db * dc;
                s.pos.push(r);
                s.j.push(j);
                s.m.push(mm);
                s.ds.push(db * dc);
            }
        }
    }
    s
}

/// Only the Poynting flux through the box (W, peak phasors) from one frequency's gathered DFT: the
/// cheap settling metric for far-field runs.
pub fn box_flux(plan: &Plan, g: &Grid, data: &[Vec<f32>; 6]) -> f64 {
    let ax = [&g.x, &g.y, &g.z];
    let at = |pi: usize, jb: usize, jc: usize| -> Cx {
        let pt = &plan.patches[pi];
        let nb = pt.b.1 - pt.b.0 + 1;
        let k = pt.off + (jc - pt.c.0) * nb + (jb - pt.b.0);
        let d = &data[pt.comp];
        Cx::new(d[2 * k] as f64, d[2 * k + 1] as f64)
    };
    let mut p = 0.0;
    for &(a, side, [eb, ec, hb0, hb1, hc0, hc1]) in &plan.faces {
        let (b, c) = ((a + 1) % 3, (a + 2) % 3);
        for n in plan.lo[c]..plan.hi[c] {
            for m in plan.lo[b]..plan.hi[b] {
                let e_b = (at(eb, m, n) + at(eb, m, n + 1)) * 0.5;
                let e_c = (at(ec, m, n) + at(ec, m + 1, n)) * 0.5;
                let h_b = (at(hb0, m, n) + at(hb0, m + 1, n) + at(hb1, m, n) + at(hb1, m + 1, n)) * 0.25;
                let h_c = (at(hc0, m, n) + at(hc0, m, n + 1) + at(hc1, m, n) + at(hc1, m, n + 1)) * 0.25;
                p += 0.5 * side * (e_b * h_c.conj() - e_c * h_b.conj()).re * (ax[b][m + 1] - ax[b][m]) * (ax[c][n + 1] - ax[c][n]);
            }
        }
    }
    p
}

/// Far field on a theta/phi grid: r E_theta and r E_phi (V, for the e^{-jkr} reference at the box centre).
#[derive(Clone, Debug)]
pub struct Pattern {
    pub theta: Vec<f64>,
    pub phi: Vec<f64>,
    /// index [it * nphi + ip]
    pub eth: Vec<Cx>,
    pub eph: Vec<Cx>,
}

/// theta 0..=180 and phi 0..360 (exclusive) in `step` degrees.
pub fn grid_deg(step: f64) -> (Vec<f64>, Vec<f64>) {
    let nt = (180.0 / step).round() as usize;
    let np = (360.0 / step).round() as usize;
    ((0..=nt).map(|k| (k as f64 * 180.0 / nt as f64).to_radians()).collect(), (0..np).map(|k| (k as f64 * 360.0 / np as f64).to_radians()).collect())
}

pub fn pattern(s: &Sources, f: f64, theta: &[f64], phi: &[f64]) -> Pattern {
    let k = 2.0 * PI * f / C0;
    let dirs: Vec<(f64, f64)> = theta.iter().flat_map(|&t| phi.iter().map(move |&p| (t, p))).collect();
    let out: Vec<(Cx, Cx)> = dirs
        .par_iter()
        .map(|&(th, ph)| {
            let (st, ct) = th.sin_cos();
            let (sp, cp) = ph.sin_cos();
            let rhat = [st * cp, st * sp, ct];
            let mut n = [Cx::default(); 3];
            let mut l = [Cx::default(); 3];
            for q in 0..s.pos.len() {
                let r = s.pos[q];
                let (si, co) = (k * (rhat[0] * r[0] + rhat[1] * r[1] + rhat[2] * r[2])).sin_cos();
                let e = Cx::new(co, si) * s.ds[q];
                for a in 0..3 {
                    n[a] = n[a] + s.j[q][a] * e;
                    l[a] = l[a] + s.m[q][a] * e;
                }
            }
            let nth = n[0] * (ct * cp) + n[1] * (ct * sp) - n[2] * st;
            let nph = n[0] * (-sp) + n[1] * cp;
            let lth = l[0] * (ct * cp) + l[1] * (ct * sp) - l[2] * st;
            let lph = l[0] * (-sp) + l[1] * cp;
            let jk4 = Cx::new(0.0, k / (4.0 * PI));
            let eth = Cx::default() - jk4 * (lph + nth * ETA0);
            let eph = jk4 * (lth - nph * ETA0);
            (eth, eph)
        })
        .collect();
    Pattern { theta: theta.to_vec(), phi: phi.to_vec(), eth: out.iter().map(|x| x.0).collect(), eph: out.iter().map(|x| x.1).collect() }
}

impl Pattern {
    /// Radiation intensity U (W/sr) per direction.
    pub fn u(&self) -> Vec<f64> {
        self.eth.iter().zip(&self.eph).map(|(a, b)| (a.abs().powi(2) + b.abs().powi(2)) / (2.0 * ETA0)).collect()
    }

    /// Total radiated power: integral of U over the sphere (trapezoid in theta, uniform phi).
    pub fn p_rad(&self) -> f64 {
        let u = self.u();
        let (nt, np) = (self.theta.len(), self.phi.len());
        let dph = 2.0 * PI / np as f64;
        let mut p = 0.0;
        for it in 0..nt {
            let w = if it == 0 || it + 1 == nt { 0.5 } else { 1.0 } * (self.theta[nt - 1] - self.theta[0]) / (nt - 1) as f64;
            let row: f64 = (0..np).map(|ip| u[it * np + ip]).sum();
            p += row * self.theta[it].sin() * w * dph;
        }
        p
    }

    /// Maximum of U and its direction (theta, phi) in radians.
    pub fn u_max(&self) -> (f64, f64, f64) {
        let u = self.u();
        let np = self.phi.len();
        let (k, &m) = u.iter().enumerate().max_by(|a, b| a.1.partial_cmp(b.1).unwrap()).unwrap();
        (m, self.theta[k / np], self.phi[k % np])
    }

    /// Linear combination sum_k w_k P_k of patterns on the same grid.
    pub fn combine(parts: &[(&Pattern, Cx)]) -> Pattern {
        let p0 = parts[0].0;
        let n = p0.eth.len();
        let mut eth = vec![Cx::default(); n];
        let mut eph = vec![Cx::default(); n];
        for (p, w) in parts {
            for q in 0..n {
                eth[q] = eth[q] + p.eth[q] * *w;
                eph[q] = eph[q] + p.eph[q] * *w;
            }
        }
        Pattern { theta: p0.theta.clone(), phi: p0.phi.clone(), eth, eph }
    }

    pub fn scale(&self, w: Cx) -> Pattern {
        Pattern::combine(&[(self, w)])
    }
}

/// `"farfield": {"freqs": [2.44e9], "step_deg": 5, "gap": 3}` in a run spec. The box sits `gap` cells
/// inside the absorbing boundary on every side.
#[derive(Clone, Debug)]
pub struct FarOpts {
    pub freqs: Vec<f64>,
    pub step_deg: f64,
    pub gap: usize,
    /// settling tolerance of the box flux for stopping (relative change over 5 checks)
    pub tol: f64,
}

impl FarOpts {
    pub fn from_spec(spec: &serde_json::Value) -> Option<FarOpts> {
        let f = spec.get("farfield")?;
        let freqs: Vec<f64> = f.get("freqs")?.as_array()?.iter().filter_map(|x| x.as_f64()).collect();
        if freqs.is_empty() {
            return None;
        }
        Some(FarOpts {
            freqs,
            step_deg: f.get("step_deg").and_then(|x| x.as_f64()).unwrap_or(5.0),
            gap: f.get("gap").and_then(|x| x.as_u64()).unwrap_or(3) as usize,
            tol: f.get("tol").and_then(|x| x.as_f64()).unwrap_or(0.01),
        })
    }
}

/// The box `gap` cells inside the CPML (`pml` = cells per side: x-, x+, y-, y+, z-, z+).
pub fn auto_plan(p: &crate::model::Problem, pml: [usize; 6], gap: usize) -> Result<Plan, String> {
    let dims = [p.nx, p.ny, p.nz];
    let lo = std::array::from_fn(|a| pml[2 * a] + gap);
    let hi = std::array::from_fn(|a| dims[a] - 1 - pml[2 * a + 1] - gap);
    Plan::new(dims, lo, hi)
}

/// Tangential E edges on the box faces that are not free space (dielectric, loss, metal or a lumped
/// element). The transform assumes free space there.
pub fn non_free_edges(p: &crate::model::Problem, plan: &Plan) -> usize {
    let cb0 = (p.dt / crate::model::EPS0) as f32;
    plan.face_e_edges().iter().filter(|&&(c, id)| {
        let (ca, cb) = (p.ca[c][id as usize], p.cb[c][id as usize]);
        (ca - 1.0).abs() > 1e-5 || (cb - cb0).abs() > 1e-4 * cb0
    }).count()
}

/// Far-field bookkeeping over the excitations of one run set.
pub struct Collector {
    pub opts: FarOpts,
    pub plan: Plan,
    theta: Vec<f64>,
    phi: Vec<f64>,
    /// per excitation port j: per frequency (pattern of the raw run, Poynting flux)
    runs: Vec<Option<Vec<(Pattern, f64)>>>,
    /// port voltage / current spectra at the far-field frequencies: [run j][port k][freq]
    v: Vec<Vec<Vec<Cx>>>,
    i: Vec<Vec<Vec<Cx>>>,
}

impl Collector {
    pub fn new(opts: FarOpts, plan: Plan, nports: usize) -> Collector {
        let (theta, phi) = grid_deg(opts.step_deg);
        Collector { opts, plan, theta, phi, runs: vec![None; nports], v: vec![vec![]; nports], i: vec![vec![]; nports] }
    }

    /// Transform one excitation's gathered DFT; `v`/`i` are all ports' spectra at the far-field frequencies.
    pub fn add_run(&mut self, j: usize, g: &Grid, gathered: &[[Vec<f32>; 6]], v: Vec<Vec<Cx>>, i: Vec<Vec<Cx>>) {
        let per_f = self.opts.freqs.iter().enumerate().map(|(fi, &f)| {
            let s = sources(&self.plan, g, &gathered[fi]);
            (pattern(&s, f, &self.theta, &self.phi), s.p_box)
        }).collect();
        self.runs[j] = Some(per_f);
        self.v[j] = v;
        self.i[j] = i;
    }

    /// Per-port patterns per 1 W incident (other ports matched to z0), power bookkeeping, outputs.
    pub fn finish(&self, names: &[String], z0: f64, out: &Path) -> Result<serde_json::Value, String> {
        use serde_json::json;
        let np = names.len();
        let k = 1.0 / (2.0 * z0.sqrt());
        let all = self.runs.iter().all(|r| r.is_some());
        let s2 = std::f64::consts::SQRT_2; // |a|^2 / 2 = 1 W  <=>  |a| = sqrt 2
        let mut ports_json = vec![json!({}); np];
        let mut summary = Vec::new();
        let mut embedded: Vec<Vec<Option<Pattern>>> = vec![vec![None; self.opts.freqs.len()]; np];
        for (fi, &f) in self.opts.freqs.iter().enumerate() {
            let wave = |j: usize, p: usize, sign: f64| (self.v[j][p][fi] + self.i[j][p][fi] * (sign * z0)) * k;
            // S and the patterns per unit incident wave
            let mut s = vec![vec![Cx::default(); np]; np];
            if all {
                let a: Vec<Vec<Cx>> = (0..np).map(|p| (0..np).map(|j| wave(j, p, 1.0)).collect()).collect();
                let b: Vec<Vec<Cx>> = (0..np).map(|p| (0..np).map(|j| wave(j, p, -1.0)).collect()).collect();
                let ai = crate::post::inv(&a).ok_or("far field: singular wave matrix")?;
                s = crate::post::matmul(&b, &ai);
                for p in 0..np {
                    let parts: Vec<(&Pattern, Cx)> = (0..np).map(|j| (&self.runs[j].as_ref().unwrap()[fi].0, ai[j][p] * s2)).collect();
                    embedded[p][fi] = Some(Pattern::combine(&parts));
                }
            } else {
                for j in (0..np).filter(|&j| self.runs[j].is_some()) {
                    let ajj = wave(j, j, 1.0);
                    for p in 0..np {
                        s[p][j] = wave(j, p, -1.0) / ajj;
                    }
                    embedded[j][fi] = Some(self.runs[j].as_ref().unwrap()[fi].0.scale(Cx::new(s2, 0.0) / ajj));
                }
            }
            for p in 0..np {
                let Some(pat) = &embedded[p][fi] else { continue };
                let p_rad = pat.p_rad();
                let (umax, th, ph) = pat.u_max();
                let p_in = 1.0 - s[p][p].abs().powi(2);
                let p_other: f64 = (0..np).filter(|&q| q != p).map(|q| s[q][p].abs().powi(2)).sum();
                // raw-run check: integrated far field against the Poynting flux through the box
                let (run_rad, run_box) = match &self.runs[p] {
                    Some(r) => {
                        let pinc = wave(p, p, 1.0).abs().powi(2) / 2.0;
                        (r[fi].0.p_rad() / pinc, r[fi].1 / pinc)
                    }
                    None => (f64::NAN, f64::NAN),
                };
                let dbi = |x: f64| 10.0 * x.max(1e-30).log10();
                let row = json!({
                    "port": names[p], "f": f,
                    "p_rad": p_rad, "p_in": p_in, "p_other_ports": p_other, "p_loss": p_in - p_other - p_rad,
                    "directivity_dbi": dbi(4.0 * PI * umax / p_rad),
                    "gain_dbi": dbi(4.0 * PI * umax / p_in),
                    "realized_gain_dbi": dbi(4.0 * PI * umax),
                    "radiation_efficiency": p_rad / (p_in - p_other),
                    "max_dir_deg": [th.to_degrees(), ph.to_degrees()],
                    "check_run_p_far": run_rad, "check_run_p_box": run_box,
                });
                eprintln!(
                    "fieldsolve: far field {} @ {:.3} GHz: P_rad {:.3e} W/W in ({:+.1} dB), D {:.2} dBi, G {:.2} dBi, realized {:.2} dBi, rad. eff. {:.1} %, max at theta {:.0} phi {:.0}; box flux check {:.3e} vs {:.3e}",
                    names[p], f * 1e-9, p_rad, dbi(p_rad), row["directivity_dbi"].as_f64().unwrap(), row["gain_dbi"].as_f64().unwrap(),
                    row["realized_gain_dbi"].as_f64().unwrap(), 100.0 * p_rad / (p_in - p_other), th.to_degrees(), ph.to_degrees(), run_box, run_rad
                );
                summary.push(row);
                let tag = format!("{}_{:.3}GHz", names[p].chars().map(|c| if c.is_ascii_alphanumeric() { c } else { '_' }).collect::<String>(), f * 1e-9);
                polar_svg(&out.join(format!("farfield_{tag}.svg")), pat, &format!("{} at {:.3} GHz: realized gain, 1 W incident", names[p], f * 1e-9)).map_err(|e| e.to_string())?;
            }
        }
        for p in 0..np {
            let per_f: Vec<serde_json::Value> = embedded[p].iter().map(|o| match o {
                Some(pat) => json!({
                    "eth_re": pat.eth.iter().map(|c| c.re as f32).collect::<Vec<_>>(), "eth_im": pat.eth.iter().map(|c| c.im as f32).collect::<Vec<_>>(),
                    "eph_re": pat.eph.iter().map(|c| c.re as f32).collect::<Vec<_>>(), "eph_im": pat.eph.iter().map(|c| c.im as f32).collect::<Vec<_>>(),
                }),
                None => serde_json::Value::Null,
            }).collect();
            ports_json[p] = json!({"name": names[p], "patterns": per_f});
        }
        let ff = json!({
            "freqs": self.opts.freqs,
            "theta_deg": self.theta.iter().map(|t| t.to_degrees()).collect::<Vec<_>>(),
            "phi_deg": self.phi.iter().map(|t| t.to_degrees()).collect::<Vec<_>>(),
            "z0": z0,
            "normalisation": "r*E_theta, r*E_phi (V, peak) per 1 W incident at the port, other ports matched to z0; index [it * nphi + ip]",
            "box_nodes": {"lo": self.plan.lo, "hi": self.plan.hi},
            "ports": ports_json,
            "summary": summary,
        });
        std::fs::write(out.join("farfield.json"), serde_json::to_string(&ff).unwrap()).map_err(|e| e.to_string())?;
        Ok(json!(summary))
    }
}

/// Read the per-port patterns written by `Collector::finish`: (freqs, theta, phi, names, [port][freq]).
pub fn read_farfield(path: &Path) -> Result<(Vec<f64>, Vec<String>, Vec<Vec<Option<Pattern>>>), String> {
    let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?).map_err(|e| e.to_string())?;
    let fl = |x: &serde_json::Value| -> Vec<f64> { x.as_array().map(|a| a.iter().filter_map(|y| y.as_f64()).collect()).unwrap_or_default() };
    let freqs = fl(&v["freqs"]);
    let theta: Vec<f64> = fl(&v["theta_deg"]).iter().map(|d| d.to_radians()).collect();
    let phi: Vec<f64> = fl(&v["phi_deg"]).iter().map(|d| d.to_radians()).collect();
    let mut names = Vec::new();
    let mut pats = Vec::new();
    for p in v["ports"].as_array().ok_or("farfield.json: no ports")? {
        names.push(p["name"].as_str().unwrap_or("").to_string());
        pats.push(p["patterns"].as_array().map(|a| a.iter().map(|q| {
            if q.is_null() {
                return None;
            }
            let c = |re: &str, im: &str| fl(&q[re]).iter().zip(fl(&q[im])).map(|(&r, i)| Cx::new(r, i)).collect::<Vec<_>>();
            Some(Pattern { theta: theta.clone(), phi: phi.clone(), eth: c("eth_re", "eth_im"), eph: c("eph_re", "eph_im") })
        }).collect()).unwrap_or_default());
    }
    Ok((freqs, names, pats))
}

/// Three polar cuts of 4 pi U (realized gain for a 1 W incident pattern), dBi: xz (phi 0/180),
/// yz (phi 90/270) and xy (theta 90).
pub fn polar_svg(path: &Path, pat: &Pattern, title: &str) -> std::io::Result<()> {
    use std::fmt::Write as _;
    let u = pat.u();
    let (nt, np) = (pat.theta.len(), pat.phi.len());
    let g = |it: usize, ip: usize| 10.0 * (4.0 * PI * u[it * np + ip]).max(1e-30).log10();
    let gmax = (0..nt * np).map(|q| g(q / np, q % np)).fold(f64::MIN, f64::max);
    let top = (gmax / 5.0).ceil() * 5.0;
    let bot = top - 40.0;
    let (w, h, rr) = (900.0, 360.0, 120.0);
    let mut o = String::new();
    let _ = write!(o, r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" font-family="sans-serif" font-size="11"><rect width="100%" height="100%" fill="white"/><text x="12" y="18" font-size="13">{title}  (max {gmax:.1} dBi; rings every 10 dB from {top:.0} dBi)</text>"#);
    let pidx = |deg: f64| ((deg / 360.0 * np as f64).round() as usize) % np;
    let tidx = |deg: f64| ((deg / 180.0 * (nt - 1) as f64).round() as usize).min(nt - 1);
    // each cut: list of (angle on the plot, gain); plot angle measured from +vertical (theta cuts) or +x (phi cut)
    let mut cuts: Vec<(&str, Vec<(f64, f64)>, bool)> = Vec::new();
    for (name, p0, p1) in [("xz plane (phi 0/180)", 0.0, 180.0), ("yz plane (phi 90/270)", 90.0, 270.0)] {
        let mut v: Vec<(f64, f64)> = (0..nt).map(|it| (pat.theta[it], g(it, pidx(p0)))).collect();
        v.extend((0..nt).rev().map(|it| (2.0 * PI - pat.theta[it], g(it, pidx(p1)))));
        cuts.push((name, v, true));
    }
    let mut v: Vec<(f64, f64)> = (0..np).map(|ip| (pat.phi[ip], g(tidx(90.0), ip))).collect();
    v.push(v[0]);
    cuts.push(("xy plane (theta 90)", v, false));
    for (ci, (name, pts, vertical)) in cuts.iter().enumerate() {
        let (cx, cy) = (150.0 + 300.0 * ci as f64, 195.0);
        for ring in 0..=4 {
            let r = rr * ring as f64 / 4.0;
            let _ = write!(o, r##"<circle cx="{cx}" cy="{cy}" r="{r}" fill="none" stroke="#ddd"/>"##);
        }
        let _ = write!(o, r##"<line x1="{}" y1="{cy}" x2="{}" y2="{cy}" stroke="#eee"/><line x1="{cx}" y1="{}" x2="{cx}" y2="{}" stroke="#eee"/>"##, cx - rr, cx + rr, cy - rr, cy + rr);
        let path: Vec<String> = pts.iter().map(|&(ang, gd)| {
            let r = rr * ((gd - bot) / (top - bot)).clamp(0.0, 1.0);
            let (x, y) = if *vertical { (cx + r * ang.sin(), cy - r * ang.cos()) } else { (cx + r * ang.cos(), cy - r * ang.sin()) };
            format!("{x:.1},{y:.1}")
        }).collect();
        let _ = write!(o, r##"<polyline fill="none" stroke="#d62728" stroke-width="1.8" points="{}"/><text x="{cx}" y="{}" text-anchor="middle">{name}</text>"##, path.join(" "), cy + rr + 22.0);
        let (lab_a, lab_b) = if *vertical { ("+z", "+x/+y") } else { ("+y", "+x") };
        let _ = write!(o, r##"<text x="{cx}" y="{}" text-anchor="middle" fill="#888">{lab_a}</text><text x="{}" y="{}" fill="#888">{lab_b}</text>"##, cy - rr - 6.0, cx + rr + 4.0, cy + 4.0);
    }
    o.push_str("</svg>");
    std::fs::write(path, o)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A z-directed Hertzian dipole I dl at the origin, sampled on a box: the transform must give
    /// U = eta k^2 |I dl|^2 sin^2(theta) / (32 pi^2), P = eta k^2 |I dl|^2 / (12 pi) and a Poynting flux
    /// equal to P. Exact fields of the dipole fill the gathered samples at their Yee positions.
    #[test]
    fn hertzian_dipole_on_a_box() {
        let f = 3e9;
        let k = 2.0 * PI * f / C0;
        let w = 2.0 * PI * f;
        let idl = Cx::new(1e-3, 0.0);
        // uniform grid, 2 mm cells, box +-30 mm (0.6 lambda), nodes at multiples of 2 mm
        let h = 2e-3;
        let nn = 41;
        let line: Vec<f64> = (0..nn).map(|i| (i as f64 - 20.0) * h).collect();
        let g = Grid { x: line.clone(), y: line.clone(), z: line.clone() };
        let plan = Plan::new([nn, nn, nn], [5, 5, 5], [35, 35, 35]).unwrap();
        // exact fields of a z dipole (e^{jwt}): H_phi, E_r, E_theta
        let field = |r3: [f64; 3]| -> ([Cx; 3], [Cx; 3]) {
            let r = (r3[0] * r3[0] + r3[1] * r3[1] + r3[2] * r3[2]).sqrt();
            let (ct, st) = (r3[2] / r, (r3[0] * r3[0] + r3[1] * r3[1]).sqrt() / r);
            let ph = r3[1].atan2(r3[0]);
            let jkr = Cx::new(0.0, k * r);
            let em = Cx::new((k * r).cos(), -(k * r).sin()); // e^{-jkr}
            let inv = |z: Cx| Cx::new(1.0, 0.0) / z;
            let one = Cx::new(1.0, 0.0);
            // H_phi = j k I dl sin/(4 pi r) (1 + 1/(jkr)) e^{-jkr}
            let hphi = Cx::new(0.0, k) * idl * (st / (4.0 * PI * r)) * (one + inv(jkr)) * em;
            // E_r = eta I dl cos/(2 pi r^2) (1 + 1/(jkr)) e^{-jkr}
            let er = idl * (ETA0 * ct / (2.0 * PI * r * r)) * (one + inv(jkr)) * em;
            // E_theta = j eta k I dl sin/(4 pi r) (1 + 1/(jkr) - 1/(kr)^2) e^{-jkr}
            let eth = Cx::new(0.0, ETA0 * k) * idl * (st / (4.0 * PI * r)) * (one + inv(jkr) - Cx::new(1.0 / (k * r * k * r), 0.0)) * em;
            let (sp, cp) = ph.sin_cos();
            let e = [er * (st * cp) + eth * (ct * cp), er * (st * sp) + eth * (ct * sp), er * ct - eth * st];
            let hh = [hphi * (-sp), hphi * cp, Cx::default()];
            (e, hh)
        };
        let _ = w;
        let pos = |c: usize, id: u32| -> [f64; 3] {
            let id = id as usize;
            let (i, j, kk) = (id % nn, (id / nn) % nn, id / (nn * nn));
            let mut p = [line[i], line[j], line[kk]];
            // E_a: half along a; H_a: half along the other two
            if c < 3 {
                p[c] += h / 2.0;
            } else {
                for a in 0..3 {
                    if a != c - 3 {
                        p[a] += h / 2.0;
                    }
                }
            }
            p
        };
        let data: [Vec<f32>; 6] = std::array::from_fn(|c| {
            plan.gather[c]
                .iter()
                .flat_map(|&id| {
                    let (e, hh) = field(pos(c, id));
                    let v = if c < 3 { e[c] } else { hh[c - 3] };
                    [v.re as f32, v.im as f32]
                })
                .collect()
        });
        let s = sources(&plan, &g, &data);
        let (th, ph) = grid_deg(5.0);
        let pat = pattern(&s, f, &th, &ph);
        let p_exact = ETA0 * k * k * idl.abs().powi(2) / (12.0 * PI);
        let p = pat.p_rad();
        assert!((p / p_exact - 1.0).abs() < 0.03, "P_rad {p} vs {p_exact}");
        assert!((s.p_box / p_exact - 1.0).abs() < 0.03, "Poynting {} vs {p_exact}", s.p_box);
        // directivity 1.5 at the equator, null on the axis
        let (umax, tmax, _) = pat.u_max();
        let d = 4.0 * PI * umax / p;
        assert!((d - 1.5).abs() < 0.03, "D {d}");
        assert!((tmax.to_degrees() - 90.0).abs() < 1e-6);
        let u = pat.u();
        assert!(u[0] < 1e-3 * umax, "on-axis U {} vs max {umax}", u[0]);
    }
}
