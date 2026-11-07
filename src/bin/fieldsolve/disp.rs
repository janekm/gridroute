//! Wideband dielectric loss: a sum of Debye poles fitted to a flat loss tangent (the Djordjevic-Sarkar
//! idea of a continuous relaxation spectrum, sampled one pole per decade).
//!
//! A constant conductivity gives tan(delta) = sigma / (w eps), which falls as 1/f: matched at 1 GHz it
//! is 5-7x too small at 5-7 GHz, and plane-cavity modes ring far too long. Here every slab gets
//!   eps(w) = eps_inf + sum_p deps_p / (1 + j w tau_p)
//! with shared poles tau_p and per-slab amplitudes deps_p >= 0 fitted so that eps'' / eps_r follows the
//! target tan(delta)(f) over [f_lo, f_hi]; eps_inf makes eps'(f_ref) = eps_r. Optionally the target adds
//! an equivalent copper loss delta_s(f) / h for a plane pair h apart (Q_c = h / delta_s of a parallel-plate
//! line), a first-order stand-in for a surface-impedance boundary on zero-thickness copper.
//!
//! FDTD (J formulation, Taflove & Hagness 9.4): per pole the polarisation current J_p obeys
//! tau_p dJ_p/dt + J_p = eps0 deps_p dE/dt, discretised as
//!   J_p^{n+1/2} = k_p J_p^{n-1/2} + b_p (E^{n+1} - E^n),
//!   k_p = (1 - dt/(2 tau_p)) / (1 + dt/(2 tau_p)),  b_p = eps0 deps_p / tau_p / (1 + dt/(2 tau_p)).
//! The E update becomes, with e = eps0 eps_inf and B = sum b_p,
//!   E^{n+1} = ca E^n + cb (curl H - sum (1 + k_p)/2 J_p^{n-1/2}),
//!   ca = (e/dt - sigma/2 + B/2) / D,  cb = 1 / D,  D = e/dt + sigma/2 + B/2.

use crate::model::EPS0;
use std::f64::consts::PI;

const MU0: f64 = 1.256_637_062_12e-6;
/// copper resistivity (ohm m) for the equivalent conductor loss
const RHO_CU: f64 = 1.72e-8;

#[derive(Clone, Debug)]
pub struct Wideband {
    /// band over which the loss tangent is fitted flat (Hz)
    pub f_lo: f64,
    pub f_hi: f64,
    /// frequency at which eps' equals the slab's eps_r (datasheet values are usually at 1 GHz)
    pub f_ref: f64,
    /// add the equivalent copper loss of each slab's plane pair
    pub copper: bool,
}

impl Wideband {
    /// One pole per decade from f_lo/3 to 3 f_hi (relaxation frequencies, Hz).
    pub fn poles(&self) -> Vec<f64> {
        let (a, b) = (self.f_lo / 3.0, self.f_hi * 3.0);
        let n = (((b / a).log10().ceil() as usize).max(1) + 1).min(8);
        (0..n).map(|k| a * (b / a).powf(k as f64 / (n - 1) as f64)).collect()
    }
}

/// Equivalent loss tangent of copper planes `h_mm` apart at `f`: delta_s / h.
pub fn copper_tand(f: f64, h_mm: f64) -> f64 {
    if h_mm <= 0.0 {
        return 0.0;
    }
    let ds = (RHO_CU / (PI * f * MU0)).sqrt();
    ds / (h_mm * 1e-3)
}

/// Fitted slab: (eps_inf, deps per pole, worst relative error of tan(delta) over the band).
pub fn fit(eps_r: f64, tand: f64, h_mm: f64, wb: &Wideband) -> (f64, Vec<f64>, f64) {
    let poles = wb.poles();
    let np = poles.len();
    let nf = 60;
    let fs: Vec<f64> = (0..nf).map(|k| wb.f_lo * (wb.f_hi / wb.f_lo).powf(k as f64 / (nf - 1) as f64)).collect();
    let target = |f: f64| tand + if wb.copper { copper_tand(f, h_mm) } else { 0.0 };
    let basis = |f: f64, p: usize| {
        let x = f / poles[p];
        x / (1.0 + x * x)
    };
    // relative least squares on eps'' with non-negative amplitudes (active set: drop negatives, refit)
    let mut active: Vec<bool> = vec![true; np];
    let mut deps = vec![0.0; np];
    for _ in 0..np {
        let idx: Vec<usize> = (0..np).filter(|&p| active[p]).collect();
        let m = idx.len();
        if m == 0 {
            break;
        }
        let mut a = vec![vec![0.0; m + 1]; m];
        for &f in &fs {
            let y = eps_r * target(f);
            let w = 1.0 / (y * y).max(1e-30);
            for (r, &p) in idx.iter().enumerate() {
                for (c, &q) in idx.iter().enumerate() {
                    a[r][c] += w * basis(f, p) * basis(f, q);
                }
                a[r][m] += w * basis(f, p) * y;
            }
        }
        // Gauss-Jordan
        for c in 0..m {
            let piv = (c..m).max_by(|&x, &y| a[x][c].abs().partial_cmp(&a[y][c].abs()).unwrap()).unwrap();
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
        deps = vec![0.0; np];
        let mut neg = false;
        for (r, &p) in idx.iter().enumerate() {
            deps[p] = a[r][m];
            if deps[p] < 0.0 {
                active[p] = false;
                neg = true;
            }
        }
        if !neg {
            break;
        }
    }
    for d in deps.iter_mut() {
        *d = d.max(0.0);
    }
    let eps1 = |f: f64| -> f64 { (0..np).map(|p| deps[p] / (1.0 + (f / poles[p]).powi(2))).sum() };
    let eps_inf = eps_r - eps1(wb.f_ref);
    let mut worst: f64 = 0.0;
    for &f in &fs {
        let e2: f64 = (0..np).map(|p| deps[p] * basis(f, p)).sum();
        let tn = e2 / (eps_inf + eps1(f));
        worst = worst.max((tn / target(f) - 1.0).abs());
    }
    (eps_inf, deps, worst)
}

/// Per-pole FDTD constants for time step `dt`: (k_p, (1 + k_p)/2, b_p per unit deps).
pub fn pole_consts(poles: &[f64], dt: f64) -> Vec<(f64, f64, f64)> {
    poles
        .iter()
        .map(|&fp| {
            let tau = 1.0 / (2.0 * PI * fp);
            let den = 1.0 + dt / (2.0 * tau);
            let k = (1.0 - dt / (2.0 * tau)) / den;
            (k, 0.5 * (1.0 + k), EPS0 / tau / den)
        })
        .collect()
}

/// Complex relative permittivity of a fitted slab at `f`: (eps', eps'').
pub fn eps_at(eps_inf: f64, deps: &[f64], poles: &[f64], f: f64) -> (f64, f64) {
    let mut re = eps_inf;
    let mut im = 0.0;
    for (d, fp) in deps.iter().zip(poles) {
        let x = f / fp;
        re += d / (1.0 + x * x);
        im += d * x / (1.0 + x * x);
    }
    (re, im)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fr4_fit_is_flat_and_hits_eps_r() {
        let wb = Wideband { f_lo: 0.3e9, f_hi: 10e9, f_ref: 1e9, copper: false };
        let (einf, deps, worst) = fit(4.3, 0.015, 0.0, &wb);
        assert!(worst < 0.1, "tan(delta) ripple {worst}");
        let poles = wb.poles();
        let (e1, _) = eps_at(einf, &deps, &poles, 1e9);
        assert!((e1 - 4.3).abs() < 1e-9);
        let mut last = f64::INFINITY;
        for f in [0.5e9, 2.44e9, 4.88e9, 7.32e9] {
            let (re, im) = eps_at(einf, &deps, &poles, f);
            assert!((im / re / 0.015 - 1.0).abs() < 0.1, "tan(delta) at {f}: {}", im / re);
            // eps' falls slowly with frequency (Kramers-Kronig), by a few percent over the band
            assert!(re < last && (re / 4.3 - 1.0).abs() < 0.05, "eps' at {f}: {re}");
            last = re;
        }
        // with copper loss for a 0.21 mm plane pair the target rises at low frequency (~ f^-1/2)
        let wb2 = Wideband { copper: true, ..wb };
        let (_, _, w2) = fit(4.3, 0.015, 0.21, &wb2);
        assert!(w2 < 0.1, "with copper: ripple {w2}");
        assert!((copper_tand(2.44e9, 0.21) - 0.0063).abs() < 0.0005);
    }
}
