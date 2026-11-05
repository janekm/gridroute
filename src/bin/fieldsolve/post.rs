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

/// DFT of samples `x[n]` taken at times `(n + t_off) * dt`.
pub fn dft(x: &[f64], dt: f64, t_off: f64, freqs: &[f64]) -> Vec<Cx> {
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

/// Power-wave S-parameter column for excitation port `j`: `S_ij = b_i / a_j`.
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

/// Touchstone v1 file, RI format; `s[i][j][f]`.
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
            if n > 2 && q > 0 && q % 4 == 0 {
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
