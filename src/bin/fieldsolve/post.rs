//! Port post-processing: voltage and current spectra, impedance, S-parameters, inductance, Touchstone.

use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Cx {
    pub re: f64,
    pub im: f64,
}

impl Cx {
    pub fn new(re: f64, im: f64) -> Cx {
        Cx { re, im }
    }
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }
    pub fn conj(self) -> Cx {
        Cx::new(self.re, -self.im)
    }
    pub fn sqrt(self) -> Cx {
        let r = self.abs().sqrt();
        let t = 0.5 * self.arg();
        Cx::new(r * t.cos(), r * t.sin())
    }
}
impl std::ops::Add for Cx {
    type Output = Cx;
    fn add(self, o: Cx) -> Cx {
        Cx::new(self.re + o.re, self.im + o.im)
    }
}
impl std::ops::Sub for Cx {
    type Output = Cx;
    fn sub(self, o: Cx) -> Cx {
        Cx::new(self.re - o.re, self.im - o.im)
    }
}
impl std::ops::Mul for Cx {
    type Output = Cx;
    fn mul(self, o: Cx) -> Cx {
        Cx::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    }
}
impl std::ops::Mul<f64> for Cx {
    type Output = Cx;
    fn mul(self, o: f64) -> Cx {
        Cx::new(self.re * o, self.im * o)
    }
}
impl std::ops::Div for Cx {
    type Output = Cx;
    fn div(self, o: Cx) -> Cx {
        let d = o.re * o.re + o.im * o.im;
        Cx::new((self.re * o.re + self.im * o.im) / d, (self.im * o.re - self.re * o.im) / d)
    }
}

/// DFT of samples `x[n]` taken at times `(n + t_off) * dt`. The phase advances by a fixed rotation per
/// sample instead of a sin/cos each; it is recomputed exactly every `REANCHOR` samples, so rounding in the
/// recurrence cannot accumulate over long records.
pub fn dft(x: &[f64], dt: f64, t_off: f64, freqs: &[f64]) -> Vec<Cx> {
    const REANCHOR: usize = 1024;
    freqs
        .iter()
        .map(|&f| {
            let w = -2.0 * PI * f * dt;
            let (ss, cs) = w.sin_cos();
            let (mut re, mut im) = (0.0, 0.0);
            for (b, chunk) in x.chunks(REANCHOR).enumerate() {
                let (mut s, mut c) = (w * ((b * REANCHOR) as f64 + t_off)).sin_cos();
                for &v in chunk {
                    re += v * c;
                    im += v * s;
                    (c, s) = (c * cs - s * ss, s * cs + c * ss);
                }
            }
            Cx::new(re * dt, im * dt)
        })
        .collect()
}

#[cfg(test)]
fn dft_direct(x: &[f64], dt: f64, t_off: f64, freqs: &[f64]) -> Vec<Cx> {
    freqs
        .iter()
        .map(|&f| {
            let w = 2.0 * PI * f;
            let (mut re, mut im) = (0.0, 0.0);
            for (n, &v) in x.iter().enumerate() {
                let (s, c) = (-w * (n as f64 + t_off) * dt).sin_cos();
                re += v * c;
                im += v * s;
            }
            Cx::new(re * dt, im * dt)
        })
        .collect()
}

#[cfg(test)]
mod dft_tests {
    use super::*;
    #[test]
    fn recurrence_matches_direct_dft_including_dc_offsets_long_records_and_nulls() {
        let dt = 1e-13;
        // a long record (beyond many re-anchor blocks) of a decaying tone plus a DC step
        let n = 300_000;
        let x: Vec<f64> = (0..n).map(|k| { let t = k as f64 * dt; 0.3 + (2.0 * PI * 2.4e9 * t).sin() * (-t / 2e-9).exp() }).collect();
        let freqs: Vec<f64> = std::iter::once(0.0).chain((1..300).map(|k| k as f64 * 3e7)).collect();
        for t_off in [0.0, 0.5, 1.0] {
            let a = dft(&x, dt, t_off, &freqs);
            let b = dft_direct(&x, dt, t_off, &freqs);
            let peak = b.iter().map(|v| v.abs()).fold(0.0, f64::max);
            for (k, (p, q)) in a.iter().zip(&b).enumerate() {
                let e = (*p - *q).abs();
                // relative to each bin where it is not tiny, and to the peak everywhere (deep nulls)
                assert!(e <= 1e-9 * q.abs().max(1e-6 * peak), "t_off {t_off} bin {k}: {p:?} vs {q:?}");
            }
        }
    }
}

/// Voltage and current spectra of one port. `v` are the probe voltages at integer steps 0..=N, `vs` the
/// source voltage at half steps (empty for a passive port), `r` the port resistance. The current into the
/// structure at step n + 1/2 is `(Vs - (V_n + V_{n+1}) / 2) / R`, exactly as the update applies it.
pub fn port_spectra(v: &[f32], vs: &[f32], r: f64, dt: f64, freqs: &[f64]) -> (Vec<Cx>, Vec<Cx>) {
    let vv: Vec<f64> = v.iter().map(|&x| x as f64).collect();
    let n = v.len().saturating_sub(1);
    let i: Vec<f64> = (0..n)
        .map(|s| {
            let src = vs.get(s).copied().unwrap_or(0.0) as f64;
            (src - 0.5 * (vv[s] + vv[s + 1])) / r
        })
        .collect();
    (dft(&vv, dt, 0.0, freqs), dft(&i, dt, 0.5, freqs))
}

pub fn linspace(a: f64, b: f64, n: usize) -> Vec<f64> {
    if n <= 1 {
        return vec![a];
    }
    (0..n).map(|k| a + (b - a) * k as f64 / (n - 1) as f64).collect()
}

/// Inverse of a small complex matrix (Gauss-Jordan with partial pivoting).
pub fn inv(m: &[Vec<Cx>]) -> Option<Vec<Vec<Cx>>> {
    let n = m.len();
    let mut a: Vec<Vec<Cx>> = m.iter().enumerate().map(|(i, r)| {
        let mut row = r.clone();
        row.extend((0..n).map(|j| if i == j { Cx::new(1.0, 0.0) } else { Cx::default() }));
        row
    }).collect();
    for c in 0..n {
        let piv = (c..n).max_by(|&x, &y| a[x][c].abs().partial_cmp(&a[y][c].abs()).unwrap())?;
        if a[piv][c].abs() < 1e-300 {
            return None;
        }
        a.swap(c, piv);
        let d = a[c][c];
        for k in 0..2 * n {
            a[c][k] = a[c][k] / d;
        }
        for r in 0..n {
            if r != c {
                let f = a[r][c];
                for k in 0..2 * n {
                    let t = a[c][k] * f;
                    a[r][k] = a[r][k] - t;
                }
            }
        }
    }
    Some(a.into_iter().map(|r| r[n..].to_vec()).collect())
}

pub fn matmul(a: &[Vec<Cx>], b: &[Vec<Cx>]) -> Vec<Vec<Cx>> {
    let (n, m, p) = (a.len(), b.len(), b[0].len());
    (0..n).map(|i| (0..p).map(|j| (0..m).fold(Cx::default(), |s, k| s + a[i][k] * b[k][j])).collect()).collect()
}

/// Full S-matrix from all excitations: `v[j][k][f]`, `i[j][k][f]` are port k's voltage and current when
/// port j is excited. With the wave matrices A = (V + Z0 I) / (2 sqrt Z0) and B = (V - Z0 I) / (2 sqrt Z0)
/// (column j = excitation j), S = B A^-1. This holds for any port terminations, unlike b_i / a_j, which
/// assumes the other ports are matched to Z0. Returns `s[i][j][f]`.
pub fn s_matrix(v: &[Vec<Vec<Cx>>], i: &[Vec<Vec<Cx>>], z0: f64) -> Vec<Vec<Vec<Cx>>> {
    let n = v.len();
    let nf = v[0][0].len();
    let k = 1.0 / (2.0 * z0.sqrt());
    let mut s = vec![vec![vec![Cx::default(); nf]; n]; n];
    for f in 0..nf {
        let a: Vec<Vec<Cx>> = (0..n).map(|p| (0..n).map(|j| (v[j][p][f] + i[j][p][f] * z0) * k).collect()).collect();
        let b: Vec<Vec<Cx>> = (0..n).map(|p| (0..n).map(|j| (v[j][p][f] - i[j][p][f] * z0) * k).collect()).collect();
        if let Some(ai) = inv(&a) {
            let sf = matmul(&b, &ai);
            for p in 0..n {
                for j in 0..n {
                    s[p][j][f] = sf[p][j];
                }
            }
        }
    }
    s
}

/// Power-wave reflection/transmission for excitation port `j` with the other ports terminated in their own
/// resistances: `b_i / a_j`. Equal to S only when those terminations are matched to Z0.
pub fn s_column(v: &[Vec<Cx>], i: &[Vec<Cx>], z0: f64, j: usize) -> Vec<Vec<Cx>> {
    let k = 1.0 / (2.0 * z0.sqrt());
    let nf = v[j].len();
    (0..v.len())
        .map(|p| {
            (0..nf)
                .map(|f| {
                    let a = (v[j][f] + i[j][f] * z0) * k;
                    let b = (v[p][f] - i[p][f] * z0) * k;
                    b / a
                })
                .collect()
        })
        .collect()
}

/// Touchstone v1 file, RI format; `s[i][j][f]`. Two ports: S11 S21 S12 S22 on one line. More ports: one
/// matrix row per line (the first on the frequency's line), continued after every four pairs.
pub fn write_touchstone(path: &std::path::Path, freqs: &[f64], s: &[Vec<Vec<Cx>>], z0: f64) -> std::io::Result<()> {
    use std::fmt::Write as _;
    let n = s.len();
    let mut o = String::new();
    let _ = writeln!(o, "! fieldsolve FDTD\n# Hz S RI R {z0}");
    for (fi, f) in freqs.iter().enumerate() {
        let _ = write!(o, "{f:.6e}");
        // Touchstone 2-port order is S11 S21 S12 S22; N > 2 is row-major
        let order: Vec<(usize, usize)> = if n == 2 { vec![(0, 0), (1, 0), (0, 1), (1, 1)] } else { (0..n).flat_map(|i| (0..n).map(move |j| (i, j))).collect() };
        for (q, (i, j)) in order.into_iter().enumerate() {
            if n > 2 && q > 0 && (j == 0 || j % 4 == 0) {
                o.push('\n');
            }
            let x = s[i][j][fi];
            let _ = write!(o, " {:.8e} {:.8e}", x.re, x.im);
        }
        o.push('\n');
    }
    std::fs::write(path, o)
}

/// Least-squares fit of `Z(f) = R + jwL` over the given band (rows with `f` in `[f0, f1]`).
pub fn fit_rl(freqs: &[f64], z: &[Cx], f0: f64, f1: f64) -> (f64, f64) {
    let (mut sr, mut sl, mut n) = (0.0f64, 0.0f64, 0.0f64);
    for (f, z) in freqs.iter().zip(z) {
        if *f >= f0 && *f <= f1 && *f > 0.0 {
            sr += z.re;
            sl += z.im / (2.0 * PI * f);
            n += 1.0;
        }
    }
    (sr / n.max(1.0), sl / n.max(1.0))
}

/// Time-domain reflectometry from S11 on a uniform frequency grid starting at (or near) DC.
/// Returns (t [s], step reflection rho(t), impedance Z0 (1 + rho) / (1 - rho)). The spectrum is shaped
/// by a Gaussian that falls to 0.1 at the top frequency; the step's 10-90 % rise time is then
/// 0.577 / f_g = 0.88 / f_top (returned as the last value).
pub fn tdr(freqs: &[f64], s11: &[Cx], z0: f64, t_max: f64, nt: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, f64) {
    let fmax = *freqs.last().unwrap();
    let fg = fmax / (10f64.ln()).sqrt();
    let df = freqs[1] - freqs[0];
    let w: Vec<f64> = freqs.iter().map(|&f| (-(f / fg).powi(2)).exp()).collect();
    let dt = t_max / nt as f64;
    let mut rho = Vec::with_capacity(nt);
    let mut acc = 0.0;
    let h = |t: f64| -> f64 {
        let mut s = 0.0;
        for (k, (&f, &x)) in freqs.iter().zip(s11).enumerate() {
            let wk = if k == 0 && f == 0.0 { 0.5 } else { 1.0 };
            let (sn, cs) = (2.0 * std::f64::consts::PI * f * t).sin_cos();
            s += wk * w[k] * (x.re * cs - x.im * sn);
        }
        2.0 * df * s
    };
    // the window spreads a reflection at t = 0 symmetrically in time (a Gaussian of width ~ 1/(pi fg)), so
    // the step starts integrating from where that has died away, not from t = 0
    let n_neg = (4.0 / fg / dt).ceil() as usize;
    let mut prev = h(-(n_neg as f64) * dt);
    for k in (0..n_neg).rev() {
        let cur = h(-(k as f64) * dt);
        acc += 0.5 * (prev + cur) * dt;
        prev = cur;
    }
    let ts: Vec<f64> = (0..nt).map(|k| k as f64 * dt).collect();
    for (k, &t) in ts.iter().enumerate() {
        if k > 0 {
            let cur = h(t);
            acc += 0.5 * (prev + cur) * dt;
            prev = cur;
        }
        rho.push(acc);
    }
    let z: Vec<f64> = rho.iter().map(|&r| z0 * (1.0 + r) / (1.0 - r).max(1e-9)).collect();
    (ts, rho, z, 0.577 / fg)
}
