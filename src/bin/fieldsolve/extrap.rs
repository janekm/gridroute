//! Late-time extrapolation of a run's running DFTs, so that a run can stop before its slow resonances have
//! rung down.
//!
//! After the source is off, every field quantity is a sum of the structure's modes,
//! f(t) = sum_p a_p e^{s_p t}, with the same complex frequencies s_p everywhere. A running DFT at angular
//! frequency w therefore approaches its limit geometrically:
//!
//!   X(T) = X_inf - sum_p A_p e^{(s_p - i w) T},
//!
//! with one amplitude A_p per accumulator but the same exponents for all of them. The exponents come from
//! the excited port's voltage, sampled every step: linear prediction on its late part, decimated to the
//! running DFT's own sampling. A handful of snapshots of the accumulators, taken at the energy checks, then
//! give X_inf as one fixed linear combination of the snapshots per frequency (least squares, the same
//! weights for every sample of the far-field box). The port voltages are continued with the fitted modes,
//! so the port spectra get the same treatment.
//!
//! The extrapolated objective is only trusted once successive extrapolations, from different snapshot
//! sets and pole fits, agree (`StopState`).

use crate::post::Cx;

/// Extrapolation settings (`"extrapolate"` in a spec's objective).
#[derive(Clone, Debug)]
pub struct Opts {
    /// Accumulator snapshots kept (one per energy check).
    pub snapshots: usize,
    /// Largest number of decimated port samples in the pole fit.
    pub window: usize,
    /// Linear-prediction order (upper bound; also limited by the window).
    pub order: usize,
    /// Time step of the run (s) and the decimation of the port signal (steps; the running DFT's own).
    pub dt: f64,
    pub decimate: usize,
    /// Stop rule: every extrapolated component within `tol_db` over the last `stable_checks` checks.
    pub stop: bool,
    pub tol_db: f64,
    pub stable_checks: usize,
}

impl Default for Opts {
    fn default() -> Opts {
        Opts { snapshots: 10, window: 160, order: 24, dt: 0.0, decimate: 1, stop: false, tol_db: 0.05, stable_checks: 4 }
    }
}

fn cexp(z: Cx) -> Cx {
    let m = z.re.exp();
    Cx::new(m * z.im.cos(), m * z.im.sin())
}

fn cln(z: Cx) -> Cx {
    Cx::new(z.abs().ln(), z.arg())
}

/// Solve A x = b (complex, Gaussian elimination with partial pivoting); None if singular.
pub fn csolve(mut a: Vec<Vec<Cx>>, mut b: Vec<Cx>) -> Option<Vec<Cx>> {
    let n = b.len();
    for c in 0..n {
        let piv = (c..n).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs()))?;
        if a[piv][c].abs() < 1e-300 {
            return None;
        }
        a.swap(c, piv);
        b.swap(c, piv);
        for r in c + 1..n {
            let f = a[r][c] / a[c][c];
            if f.abs() == 0.0 {
                continue;
            }
            for k in c..n {
                let t = a[c][k];
                a[r][k] = a[r][k] - f * t;
            }
            let t = b[c];
            b[r] = b[r] - f * t;
        }
    }
    let mut x = vec![Cx::default(); n];
    for r in (0..n).rev() {
        let mut s = b[r];
        for k in r + 1..n {
            s = s - a[r][k] * x[k];
        }
        x[r] = s / a[r][r];
    }
    Some(x)
}

/// Least squares min |B x - y| through the normal equations with a small relative ridge, on unit-norm
/// columns (mode columns can differ by many orders of magnitude).
fn clstsq(b: &[Vec<Cx>], y: &[Cx]) -> Option<Vec<Cx>> {
    let n = b.first()?.len();
    let norm: Vec<f64> = (0..n).map(|j| b.iter().map(|r| r[j].abs().powi(2)).sum::<f64>().sqrt().max(1e-300)).collect();
    let mut g = vec![vec![Cx::default(); n]; n];
    let mut r = vec![Cx::default(); n];
    for (row, &yv) in b.iter().zip(y) {
        for i in 0..n {
            let ci = row[i].conj() * (1.0 / norm[i]);
            r[i] = r[i] + ci * yv;
            for j in 0..n {
                g[i][j] = g[i][j] + ci * row[j] * (1.0 / norm[j]);
            }
        }
    }
    for (i, gi) in g.iter_mut().enumerate() {
        gi[i].re += 1e-12;
    }
    let x = csolve(g, r)?;
    Some(x.iter().zip(&norm).map(|(v, &nm)| *v * (1.0 / nm)).collect())
}

/// Roots of the monic polynomial z^m + c[0] z^{m-1} + ... + c[m-1] (Aberth-Ehrlich iteration).
pub fn roots(c: &[f64]) -> Vec<Cx> {
    let m = c.len();
    if m == 0 {
        return Vec::new();
    }
    let eval = |z: Cx| -> (Cx, Cx) {
        // p and p' by Horner
        let (mut p, mut d) = (Cx::new(1.0, 0.0), Cx::default());
        for &ck in c {
            d = d * z + p;
            p = p * z + Cx::new(ck, 0.0);
        }
        (p, d)
    };
    // starting points on a circle bounding the roots (Cauchy), rotated off the axes
    let r0 = 1.0 + c.iter().fold(0.0f64, |a, &x| a.max(x.abs()));
    let r0 = r0.min(4.0).max(0.5);
    let mut z: Vec<Cx> = (0..m).map(|k| cexp(Cx::new(r0.ln(), 2.0 * std::f64::consts::PI * (k as f64 + 0.25) / m as f64))).collect();
    for _ in 0..500 {
        let mut moved = 0.0f64;
        for k in 0..m {
            let (p, d) = eval(z[k]);
            if p.abs() == 0.0 {
                continue;
            }
            let ratio = p / d;
            let mut s = Cx::default();
            for j in 0..m {
                if j != k {
                    let diff = z[k] - z[j];
                    if diff.abs() > 0.0 {
                        s = s + Cx::new(1.0, 0.0) / diff;
                    }
                }
            }
            let w = ratio / (Cx::new(1.0, 0.0) - ratio * s);
            z[k] = z[k] - w;
            moved = moved.max(w.abs() / z[k].abs().max(1e-30));
        }
        if moved < 1e-14 {
            break;
        }
    }
    z
}

/// Modes fitted to the late part of a signal: continuous-time exponents `s` (1/s) and complex amplitudes at
/// the last sample, x(t_end + t) = Re sum_p amp_p e^{s_p t}.
#[derive(Clone, Debug, Default)]
pub struct Modes {
    pub s: Vec<Cx>,
    pub amp: Vec<Cx>,
    /// Index of the last sample used.
    pub end: usize,
    /// Residual rms over signal rms in the fitted window.
    pub residual: f64,
}

/// Fit decaying modes to `x[start..]`, decimated by `dec`, with linear prediction of order up to `order`
/// over at most `window` samples.
pub fn fit_modes(x: &[f32], start: usize, dec: usize, window: usize, order: usize, dt: f64) -> Option<Modes> {
    let end = x.len().checked_sub(1)?;
    if end < start {
        return None;
    }
    let dec = dec.max(1);
    let l = ((end - start) / dec + 1).min(window);
    let m = order.min(l / 3);
    if m < 2 {
        return None;
    }
    let y: Vec<f64> = (0..l).map(|k| x[end - (l - 1 - k) * dec] as f64).collect();
    let scale = y.iter().fold(0.0f64, |a, v| a.max(v.abs()));
    if scale == 0.0 {
        return None;
    }
    let y: Vec<f64> = y.iter().map(|v| v / scale).collect();
    // forward linear prediction y[n] = -sum_k c_k y[n-k], ridge-regularised normal equations (the ridge
    // keeps the extra roots of an over-ordered predictor inside the unit circle)
    let mut r = vec![vec![Cx::default(); m]; m];
    let mut rhs = vec![Cx::default(); m];
    for n in m..l {
        for j in 0..m {
            rhs[j].re -= y[n] * y[n - 1 - j];
            for k in 0..m {
                r[j][k].re += y[n - 1 - j] * y[n - 1 - k];
            }
        }
    }
    let tr: f64 = (0..m).map(|i| r[i][i].re).sum::<f64>() / m as f64;
    for (i, ri) in r.iter_mut().enumerate() {
        ri[i].re += 1e-10 * tr.max(1e-300);
    }
    let c: Vec<f64> = csolve(r, rhs)?.iter().map(|v| v.re).collect();
    let z: Vec<Cx> = roots(&c).into_iter().filter(|z| z.abs() < 1.0 && z.abs() > 1e-6).collect();
    if z.is_empty() {
        return None;
    }
    // amplitudes at the last sample: y[l-1-j] = sum_p amp_p z_p^{-j}
    let rows: Vec<Vec<Cx>> = (0..l).map(|k| z.iter().map(|&zp| cexp(cln(zp) * (k as f64 - (l - 1) as f64))).collect()).collect();
    let yc: Vec<Cx> = y.iter().map(|&v| Cx::new(v, 0.0)).collect();
    let amp = clstsq(&rows, &yc)?;
    let mut res = 0.0;
    let mut tot = 0.0;
    for (k, row) in rows.iter().enumerate() {
        let f: Cx = row.iter().zip(&amp).fold(Cx::default(), |a, (&b, &c)| a + b * c);
        res += (f.re - y[k]).powi(2);
        tot += y[k] * y[k];
    }
    let s: Vec<Cx> = z.iter().map(|&zp| cln(zp) * (1.0 / (dec as f64 * dt))).collect();
    Some(Modes { s, amp: amp.iter().map(|a| *a * scale).collect(), end, residual: (res / tot.max(1e-300)).sqrt() })
}

impl Modes {
    /// Continue the signal by `n` samples after `end` (real part of the mode sum).
    pub fn continue_signal(&self, n: usize, dt: f64) -> Vec<f32> {
        let steps: Vec<Cx> = self.s.iter().map(|&s| cexp(s * dt)).collect();
        let mut cur = self.amp.clone();
        (0..n)
            .map(|_| {
                let mut v = 0.0;
                for (c, st) in cur.iter_mut().zip(&steps) {
                    *c = *c * *st;
                    v += c.re;
                }
                v as f32
            })
            .collect()
    }

    /// Samples until every mode has decayed by `nepers` from `end`, capped.
    pub fn ring_length(&self, dt: f64, nepers: f64, cap: usize) -> usize {
        let slow = self.s.iter().map(|s| -s.re).fold(f64::INFINITY, f64::min);
        if !(slow > 0.0) {
            return cap;
        }
        ((nepers / slow / dt).ceil() as usize).min(cap)
    }

    /// The `k` modes with the largest amplitude at `end`, completed to conjugate pairs where the partner
    /// is also present.
    pub fn strongest(&self, k: usize) -> Vec<Cx> {
        let mut idx: Vec<usize> = (0..self.s.len()).collect();
        idx.sort_by(|&a, &b| self.amp[b].abs().total_cmp(&self.amp[a].abs()));
        let mut out: Vec<Cx> = Vec::new();
        for &i in &idx {
            if out.len() >= k {
                break;
            }
            out.push(self.s[i]);
        }
        out
    }
}

/// Weights w_m with X_inf = sum_m w_m X(t_m) for X(t) = X_inf - sum_p A_p e^{lambda_p t} (least squares over
/// the snapshots). `t` relative to the latest snapshot (s).
pub fn limit_weights(t: &[f64], lambda: &[Cx]) -> Option<Vec<Cx>> {
    let n = lambda.len() + 1;
    if t.len() < n {
        return None;
    }
    let b: Vec<Vec<Cx>> = t.iter().map(|&tm| std::iter::once(Cx::new(1.0, 0.0)).chain(lambda.iter().map(|&l| cexp(l * tm))).collect()).collect();
    // column 0 of (B^H B)^-1, g; then w_m = conj(sum_j g_j B_mj)
    let mut gm = vec![vec![Cx::default(); n]; n];
    for row in &b {
        for i in 0..n {
            for j in 0..n {
                gm[i][j] = gm[i][j] + row[i].conj() * row[j];
            }
        }
    }
    let tr: f64 = (0..n).map(|i| gm[i][i].re).sum::<f64>() / n as f64;
    for (i, gi) in gm.iter_mut().enumerate() {
        gi[i].re += 1e-12 * tr;
    }
    let mut e0 = vec![Cx::default(); n];
    e0[0] = Cx::new(1.0, 0.0);
    let g = csolve(gm, e0)?;
    Some(b.iter().map(|row| row.iter().zip(&g).fold(Cx::default(), |a, (&bm, &gj)| a + gj * bm).conj()).collect())
}

/// Apply per-frequency weights to accumulator snapshots (interleaved re, im).
pub fn combine(snaps: &[&[[Vec<f32>; 6]]], w: &[Vec<Cx>]) -> Vec<[Vec<f32>; 6]> {
    use rayon::prelude::*;
    let nf = snaps[0].len();
    (0..nf)
        .map(|fi| {
            std::array::from_fn(|c| {
                let len = snaps[0][fi][c].len();
                let mut out = vec![0f32; len];
                out.par_chunks_mut(2).enumerate().for_each(|(i, o)| {
                    let mut acc = Cx::default();
                    for (m, s) in snaps.iter().enumerate() {
                        acc = acc + w[fi][m] * Cx::new(s[fi][c][2 * i] as f64, s[fi][c][2 * i + 1] as f64);
                    }
                    o[0] = acc.re as f32;
                    o[1] = acc.im as f32;
                });
                out
            })
        })
        .collect()
}

/// Running extrapolation state, kept by the stop rules across checks.
#[derive(Default)]
pub struct State {
    /// (time of the last accumulation in s, accumulators per gather frequency)
    snaps: std::collections::VecDeque<(f64, Vec<[Vec<f32>; 6]>)>,
    /// The latest extrapolation: port voltages continued and the accumulators' limits.
    pub last: Option<(Vec<Vec<f32>>, Vec<[Vec<f32>; 6]>)>,
    pub last_modes: Option<Modes>,
}

impl State {
    /// Record this check's accumulators and, once enough snapshots are kept, extrapolate: continued port
    /// voltages and the accumulators' limits. `src_end`: first step after the source; `port`: the
    /// excited port; `freqs`: the gather frequencies; `dft_every`: accumulation interval in steps.
    #[allow(clippy::too_many_arguments)]
    pub fn update(&mut self, o: &Opts, steps: usize, src_end: usize, port: usize, v: &[Vec<f32>], gathered: &[[Vec<f32>; 6]], freqs: &[f64], dft_every: usize) -> bool {
        self.last = None;
        if steps <= src_end {
            return false;
        }
        let t_acc = (steps / dft_every.max(1) * dft_every.max(1)) as f64 * o.dt;
        if self.snaps.back().is_none_or(|s| s.0 < t_acc) {
            self.snaps.push_back((t_acc, gathered.to_vec()));
            while self.snaps.len() > o.snapshots {
                self.snaps.pop_front();
            }
        }
        let Some(modes) = fit_modes(&v[port], src_end, o.decimate, o.window, o.order, o.dt) else {
            return false;
        };
        // continue the port voltages with the excited port's modes, amplitudes fitted per port
        let n_ext = modes.ring_length(o.dt, 14.0, 1 << 21);
        let vx: Vec<Vec<f32>> = v
            .iter()
            .map(|s| {
                let mut out = s.clone();
                if let Some(mq) = fit_amplitudes(s, src_end, &modes, o.decimate, o.window, o.dt) {
                    out.extend(mq.continue_signal(n_ext, o.dt));
                }
                out
            })
            .collect();
        let mut gx = gathered.to_vec();
        if !freqs.is_empty() && self.snaps.len() >= 3 {
            let k = (self.snaps.len() - 2).min(modes.s.len());
            let lam = modes.strongest(k);
            let latest = self.snaps.back().unwrap().0;
            let t: Vec<f64> = self.snaps.iter().map(|s| s.0 - latest).collect();
            let w: Option<Vec<Vec<Cx>>> = freqs.iter().map(|&f| limit_weights(&t, &lam.iter().map(|&s| s - Cx::new(0.0, 2.0 * std::f64::consts::PI * f)).collect::<Vec<_>>())).collect();
            if let Some(w) = w {
                let refs: Vec<&[[Vec<f32>; 6]]> = self.snaps.iter().map(|s| s.1.as_slice()).collect();
                gx = combine(&refs, &w);
            }
        }
        self.last = Some((vx, gx));
        self.last_modes = Some(modes);
        true
    }

    pub fn snapshots(&self) -> usize {
        self.snaps.len()
    }
}

/// Amplitudes of given modes in another signal (same window and decimation).
fn fit_amplitudes(x: &[f32], start: usize, modes: &Modes, dec: usize, window: usize, dt: f64) -> Option<Modes> {
    let end = x.len().checked_sub(1)?;
    let dec = dec.max(1);
    let l = ((end.checked_sub(start)?) / dec + 1).min(window);
    if l < modes.s.len() + 1 {
        return None;
    }
    let rows: Vec<Vec<Cx>> = (0..l).map(|k| modes.s.iter().map(|&s| cexp(s * (-(((l - 1 - k) * dec) as f64) * dt))).collect()).collect();
    let y: Vec<Cx> = (0..l).map(|k| Cx::new(x[end - (l - 1 - k) * dec] as f64, 0.0)).collect();
    let amp = clstsq(&rows, &y)?;
    Some(Modes { s: modes.s.clone(), amp, end, residual: 0.0 })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two damped sinusoids: linear prediction finds their exponents, and the amplitudes reproduce the
    /// continuation.
    #[test]
    fn modes_of_two_damped_sinusoids() {
        let dt = 1e-12;
        let (s1, s2) = (Cx::new(-2e8, 2.0 * std::f64::consts::PI * 2.6e9), Cx::new(-1e9, 2.0 * std::f64::consts::PI * 5.1e9));
        let f = |n: usize| {
            let t = n as f64 * dt;
            ((s1.re * t).exp() * (s1.im * t + 0.3).cos() + 0.4 * (s2.re * t).exp() * (s2.im * t - 1.0).cos()) as f32
        };
        let x: Vec<f32> = (0..6000).map(f).collect();
        let m = fit_modes(&x, 1000, 37, 120, 12, dt).unwrap();
        assert!(m.residual < 1e-5, "residual {}", m.residual);
        for s in [s1, s2] {
            let best = m.s.iter().map(|q| (*q - s).abs() / s.abs()).fold(f64::INFINITY, f64::min);
            assert!(best < 1e-4, "{s:?} not found in {:?}", m.s);
        }
        let ext = m.continue_signal(3000, dt);
        for (k, v) in ext.iter().enumerate().step_by(97) {
            let want = f(6000 + k);
            assert!((v - want).abs() < 2e-5, "step {k}: {v} vs {want}");
        }
    }

    /// The running DFT of a mode sum: three snapshots of X(T) and the exponents give X_inf.
    #[test]
    fn limit_of_a_running_dft() {
        let dt = 1e-12;
        let w = 2.0 * std::f64::consts::PI * 2.44e9;
        let s = [Cx::new(-3e8, 2.0 * std::f64::consts::PI * 2.67e9), Cx::new(-3e8, -2.0 * std::f64::consts::PI * 2.67e9)];
        let a = [Cx::new(0.5, 0.2), Cx::new(0.5, -0.2)];
        let x_inf = Cx::new(1.234, -0.567);
        // X(T) = X_inf - sum A_p e^{(s_p - iw) T}
        let x_at = |t: f64| (0..2).fold(x_inf, |acc, p| acc - a[p] * cexp((s[p] - Cx::new(0.0, w)) * t));
        let ts: Vec<f64> = (0..6).map(|m| 20000.0 * dt + m as f64 * 1000.0 * dt).collect();
        let latest = *ts.last().unwrap();
        let lam: Vec<Cx> = s.iter().map(|&sp| sp - Cx::new(0.0, w)).collect();
        let wts = limit_weights(&ts.iter().map(|t| t - latest).collect::<Vec<_>>(), &lam).unwrap();
        let est = ts.iter().zip(&wts).fold(Cx::default(), |acc, (&t, &wm)| acc + wm * x_at(t));
        assert!((est - x_inf).abs() < 1e-9, "{est:?} vs {x_inf:?}; raw {:?}", x_at(latest));
    }

    #[test]
    fn roots_of_a_cubic() {
        // (z - 0.5)(z^2 + 0.25) = z^3 - 0.5 z^2 + 0.25 z - 0.125
        let mut r = roots(&[-0.5, 0.25, -0.125]);
        r.sort_by(|a, b| a.im.total_cmp(&b.im));
        assert!((r[0] - Cx::new(0.0, -0.5)).abs() < 1e-12 && (r[1] - Cx::new(0.5, 0.0)).abs() < 1e-12 && (r[2] - Cx::new(0.0, 0.5)).abs() < 1e-12, "{r:?}");
    }
}
