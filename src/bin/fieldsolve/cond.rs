//! Conductor loss: a surface impedance on zero-thickness copper.
//!
//! Copper is a sheet on its mesh plane. A lossless sheet forces tangential E to zero; a real conductor of
//! thickness t carries its current in a skin layer and has, for current on one face,
//!   Z_s(w) = k / sigma * coth(k t),  k = (1 + j) / delta,  delta = sqrt(2 / (w mu0 sigma)),
//! which is 1 / (sigma t) at DC and (1 + j) / (sigma delta) once delta << t. Surface roughness raises the
//! loss by the Hammerstad-Jensen factor K = 1 + (2/pi) atan(1.4 (Delta / delta)^2) (Delta the rms
//! roughness), applied to the whole Z_s.
//!
//! On a sheet edge, E = Z_s * J_s with J_s the sheet current, which on the mesh is the curl of H across the
//! sheet times the dual edge length normal to it (Ampere over the dual face; displacement current in the
//! metal is negligible). For a time-domain update Z_s is fitted as a Foster network
//!   Z_s(s) = R_0 + sum_k R_k s tau_k / (1 + s tau_k)          (series chain of parallel R_k || L_k),
//! with R_0 = 1 / (sigma t) and R_k >= 0 by non-negative least squares on Re and Im over the band. Each
//! section's inductor current i_k follows tau_k di_k/dt = J_s - i_k, integrated exactly over a step with
//! J_s held constant (i_k <- J_s + (i_k - J_s) e^{-dt/tau_k}); then
//!   E = R_0 J_s + sum_k R_k (J_s - i_k).
//! The relation is algebraic in E (no eps dE/dt term), so it needs no stiff update; it is dissipative
//! because every section is a positive R in parallel with a positive L. The sheet current is taken to flow
//! on one face, which is right for a trace over a plane and for a plane under a trace or between planes
//! (current on its inner face), and overestimates the loss of an isolated strip carrying current on both
//! faces by up to 2x.

use crate::post::Cx;
use std::f64::consts::PI;

const MU0: f64 = 1.256_637_062_12e-6;

/// Conductor parameters for one copper layer.
#[derive(Clone, Debug)]
pub struct Metal {
    /// conductivity (S/m); annealed copper 5.8e7, electroplated / foil copper about 4.5-5.8e7
    pub sigma: f64,
    /// thickness (m)
    pub t: f64,
    /// rms surface roughness (m); 0 for smooth
    pub rough: f64,
}

/// The fitting band and per-layer metals (indexed like `Layout::sheets`).
#[derive(Clone, Debug)]
pub struct Conductor {
    pub metals: Vec<Metal>,
    pub f_lo: f64,
    pub f_hi: f64,
}

fn cexp(z: Cx) -> Cx {
    let m = z.re.exp();
    Cx::new(m * z.im.cos(), m * z.im.sin())
}

/// Exact surface impedance (ohm per square) at f, current on one face.
pub fn zs(m: &Metal, f: f64) -> Cx {
    let f = f.max(1.0);
    let delta = (1.0 / (PI * f * MU0 * m.sigma)).sqrt();
    let k = Cx::new(1.0 / delta, 1.0 / delta);
    let kt = k * m.t;
    // coth(x) = (1 + e^{-2x}) / (1 - e^{-2x}); for small |x| use 1/x + x/3
    let coth = if kt.abs() < 1e-3 {
        Cx::new(1.0, 0.0) / kt + kt * (1.0 / 3.0)
    } else {
        let e = cexp(kt * -2.0);
        (Cx::new(1.0, 0.0) + e) / (Cx::new(1.0, 0.0) - e)
    };
    let kr = 1.0 + 2.0 / PI * (1.4 * (m.rough / delta).powi(2)).atan();
    k * coth * (kr / m.sigma)
}

/// A fitted Foster network: R_0 plus sections (R_k, corner frequency f_k = 1 / (2 pi tau_k)).
#[derive(Clone, Debug)]
pub struct Fit {
    pub r0: f64,
    pub rk: Vec<f64>,
    pub fk: Vec<f64>,
    /// worst relative error of Re Z (the loss) over the band
    pub worst: f64,
}

impl Fit {
    pub fn z(&self, f: f64) -> Cx {
        let mut z = Cx::new(self.r0, 0.0);
        for (r, fk) in self.rk.iter().zip(&self.fk) {
            let x = f / fk;
            z = z + Cx::new(r * x * x / (1.0 + x * x), r * x / (1.0 + x * x));
        }
        z
    }
}

/// Section corner frequencies shared by all layers: two per decade from below the thickness
/// transition (delta = t at f_t) to 3 f_hi.
pub fn corners(c: &Conductor) -> Vec<f64> {
    let ft = c.metals.iter().map(|m| 1.0 / (PI * MU0 * m.sigma * m.t * m.t)).fold(f64::INFINITY, f64::min);
    let a = (ft / 3.0).min(c.f_lo / 10.0);
    let b = 3.0 * c.f_hi;
    let n = (((b / a).log10() * 2.0).ceil() as usize + 1).clamp(3, 16);
    (0..n).map(|k| a * (b / a).powf(k as f64 / (n - 1) as f64)).collect()
}

/// Fit one metal with the given corners over [min(f_t/3, f_lo/10), f_hi].
pub fn fit(m: &Metal, fk: &[f64], f_lo: f64, f_hi: f64) -> Fit {
    let r0 = 1.0 / (m.sigma * m.t);
    let ft = 1.0 / (PI * MU0 * m.sigma * m.t * m.t);
    let (fa, fb) = ((ft / 3.0).min(f_lo / 10.0), f_hi);
    let nf = 80;
    let fs: Vec<f64> = (0..nf).map(|k| fa * (fb / fa).powf(k as f64 / (nf - 1) as f64)).collect();
    let nb = fk.len();
    // rows: Re and Im of (Z_s - R_0) / |Z_s| at every frequency; columns: sections
    let mut rows: Vec<(Vec<f64>, f64)> = Vec::new();
    for &f in &fs {
        let z = zs(m, f);
        let w = 1.0 / z.abs();
        let re: Vec<f64> = fk.iter().map(|&c| {
            let x = f / c;
            w * x * x / (1.0 + x * x)
        }).collect();
        let im: Vec<f64> = fk.iter().map(|&c| {
            let x = f / c;
            w * x / (1.0 + x * x)
        }).collect();
        // the loss (Re) is what matters; Im is the small internal inductance, and with roughness scaling the
        // whole Z_s the target is not exactly causal, so Im gets a lower weight
        let wi = 0.3;
        rows.push((re, w * (z.re - r0)));
        rows.push((im.iter().map(|v| v * wi).collect(), wi * w * z.im));
    }
    let rk = nnls(&rows, nb);
    let fit = Fit { r0, rk, fk: fk.to_vec(), worst: 0.0 };
    let worst = fs.iter().filter(|&&f| f >= f_lo.min(fb)).map(|&f| (fit.z(f).re / zs(m, f).re - 1.0).abs()).fold(0.0, f64::max);
    Fit { worst, ..fit }
}

/// Non-negative least squares by active set (drop negative coefficients and refit), normal equations.
fn nnls(rows: &[(Vec<f64>, f64)], n: usize) -> Vec<f64> {
    let mut active = vec![true; n];
    let mut x = vec![0.0; n];
    for _ in 0..=n {
        let idx: Vec<usize> = (0..n).filter(|&p| active[p]).collect();
        let m = idx.len();
        if m == 0 {
            return vec![0.0; n];
        }
        let mut a = vec![vec![0.0; m + 1]; m];
        for (r, y) in rows {
            for (i, &p) in idx.iter().enumerate() {
                for (j, &q) in idx.iter().enumerate() {
                    a[i][j] += r[p] * r[q];
                }
                a[i][m] += r[p] * y;
            }
        }
        // a tiny ridge keeps nearly collinear sections solvable
        let tr: f64 = (0..m).map(|i| a[i][i]).sum::<f64>() / m as f64;
        for (i, row) in a.iter_mut().enumerate().take(m) {
            row[i] += 1e-10 * tr;
        }
        for c in 0..m {
            let piv = (c..m).max_by(|&p, &q| a[p][c].abs().partial_cmp(&a[q][c].abs()).unwrap()).unwrap();
            a.swap(c, piv);
            let d = a[c][c];
            if d.abs() < 1e-300 {
                continue;
            }
            for k in c..=m {
                a[c][k] /= d;
            }
            for r in 0..m {
                if r != c {
                    let f = a[r][c];
                    for k in c..=m {
                        a[r][k] -= f * a[c][k];
                    }
                }
            }
        }
        x = vec![0.0; n];
        let mut neg = false;
        for (i, &p) in idx.iter().enumerate() {
            x[p] = a[i][m];
            if x[p] < 0.0 {
                active[p] = false;
                neg = true;
            }
        }
        if !neg {
            break;
        }
    }
    x.iter().map(|v| v.max(0.0)).collect()
}

/// Sheet edges with a surface impedance, as the backends use them.
#[derive(Clone, Debug)]
pub struct Cond {
    /// global edge id | component << 30 (components 0 and 1: in-plane edges)
    pub ids: Vec<u32>,
    /// dual edge length normal to the sheet (m): J_s = curl H * fac
    pub fac: Vec<f32>,
    /// metal (layer) index per edge into `tab`
    pub layer: Vec<u32>,
    /// sections per edge
    pub nb: usize,
    /// per section e^{-dt / tau_k}
    pub dec: Vec<f32>,
    /// per layer [R_0, R_1 .. R_nb]
    pub tab: Vec<f32>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn copper_fit_follows_skin_effect() {
        let m = Metal { sigma: 5.8e7, t: 35e-6, rough: 0.0 };
        // DC and high-frequency limits of the exact form
        assert!((zs(&m, 1.0).re * m.sigma * m.t - 1.0).abs() < 1e-3);
        let f = 2.44e9;
        let rs = (PI * f * MU0 / m.sigma).sqrt();
        let z = zs(&m, f);
        assert!((z.re / rs - 1.0).abs() < 1e-3 && (z.im / rs - 1.0).abs() < 1e-3, "{z:?} vs {rs}");
        let c = Conductor { metals: vec![m.clone()], f_lo: 0.1e9, f_hi: 10e9 };
        let fk = corners(&c);
        let ft = fit(&m, &fk, c.f_lo, c.f_hi);
        assert!(ft.worst < 0.05, "worst {:.3} with {} sections", ft.worst, fk.len());
        // roughness raises the loss towards 2x at high frequency
        let r = Metal { rough: 2e-6, ..m };
        assert!(zs(&r, 10e9).re > 1.7 * zs(&m, 10e9).re);
    }
}
