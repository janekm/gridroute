//! `fieldsolve optimise SPACE.json`: geometry search over a few named parameters (stage 3 of the fast-mode
//! proposal). Cheap screening solves explore; full solves go only to the candidates that matter.
//!
//! The spec is parameterised by expressions: any string `"=expr"` in it (an edit's rectangle, a port position,
//! a cell size) is evaluated with the parameters, e.g. `"rect": [141.8, 74.64, "=141.9 - max(arm, 0)", 74.94]`
//! (numbers, names, `+ - * /`, parentheses, `min`, `max`, `abs`, `clamp`). Values snap to each parameter's
//! `step`, which should not be finer than what the screening mesh resolves.
//!
//! Methods:
//! - `surrogate` (default): an initial design (the start point and a Latin hypercube), then expected
//!   improvement on a Gaussian-process model of the screening objective, inside a trust region around the
//!   incumbent that grows after an improvement and shrinks after repeated failures.
//! - `pattern`: a compass search from the start point (the simple bounded search the proposal asks to
//!   compare against), with the same budget and promotion rules.
//!
//! Promotion: every `promote_every` screening solves the incumbent gets a full solve, and at the end the
//! best screening candidates within `tie_db` of the best (up to `final`) do too. Full-minus-screening
//! differences at the promoted points feed a correction model (a smooth Gaussian process) that is added to the
//! screening values when choosing the incumbent, the next candidate and the final shortlist. The answer is
//! the best full-fidelity candidate.
//!
//! Every solve is cached under `cache/LEVEL/KEY/` (`result.json` plus the resolved spec), so a repeated point is
//! never solved twice and an interrupted search resumes where it stopped.

use crate::{fidelity, gerber, sim, tune::Rng};
use serde_json::{json, Value};
use std::path::{Path, PathBuf};

// ------------------------------------------------------------------------------------------- expressions

/// Evaluate an arithmetic expression over named parameters.
pub fn eval_expr(s: &str, vars: &[(String, f64)]) -> Result<f64, String> {
    let mut p = Expr { b: s.as_bytes(), i: 0, vars };
    let v = p.sum().map_err(|e| format!("expression {s:?}: {e}"))?;
    if p.peek().is_some() {
        return Err(format!("expression {s:?}: unexpected {:?}", &s[p.i..]));
    }
    if !v.is_finite() {
        return Err(format!("expression {s:?} is not finite"));
    }
    Ok(v)
}

struct Expr<'a> {
    b: &'a [u8],
    i: usize,
    vars: &'a [(String, f64)],
}

impl Expr<'_> {
    fn peek(&mut self) -> Option<u8> {
        while self.i < self.b.len() && self.b[self.i].is_ascii_whitespace() {
            self.i += 1;
        }
        self.b.get(self.i).copied()
    }
    fn expect(&mut self, c: u8) -> Result<(), String> {
        if self.peek() == Some(c) {
            self.i += 1;
            Ok(())
        } else {
            Err(format!("expected {:?}", c as char))
        }
    }
    fn sum(&mut self) -> Result<f64, String> {
        let mut v = self.product()?;
        while let Some(c @ (b'+' | b'-')) = self.peek() {
            self.i += 1;
            let r = self.product()?;
            v = if c == b'+' { v + r } else { v - r };
        }
        Ok(v)
    }
    fn product(&mut self) -> Result<f64, String> {
        let mut v = self.unary()?;
        while let Some(c @ (b'*' | b'/')) = self.peek() {
            self.i += 1;
            let r = self.unary()?;
            v = if c == b'*' { v * r } else { v / r };
        }
        Ok(v)
    }
    fn unary(&mut self) -> Result<f64, String> {
        match self.peek() {
            Some(b'-') => {
                self.i += 1;
                Ok(-self.unary()?)
            }
            Some(b'+') => {
                self.i += 1;
                self.unary()
            }
            _ => self.atom(),
        }
    }
    fn atom(&mut self) -> Result<f64, String> {
        let b = self.b;
        match self.peek() {
            Some(b'(') => {
                self.i += 1;
                let v = self.sum()?;
                self.expect(b')')?;
                Ok(v)
            }
            Some(c) if c.is_ascii_digit() || c == b'.' => {
                let s = self.i;
                while self.i < b.len() {
                    let c = b[self.i];
                    let exp_sign = (c == b'-' || c == b'+') && self.i > s && matches!(b[self.i - 1], b'e' | b'E');
                    if c.is_ascii_digit() || c == b'.' || c == b'e' || c == b'E' || exp_sign {
                        self.i += 1;
                    } else {
                        break;
                    }
                }
                let t = std::str::from_utf8(&b[s..self.i]).unwrap();
                t.parse().map_err(|_| format!("bad number {t:?}"))
            }
            Some(c) if c.is_ascii_alphabetic() || c == b'_' => {
                let s = self.i;
                while self.i < b.len() && (b[self.i].is_ascii_alphanumeric() || b[self.i] == b'_') {
                    self.i += 1;
                }
                let name = std::str::from_utf8(&b[s..self.i]).unwrap();
                if self.peek() == Some(b'(') {
                    self.i += 1;
                    let mut a = vec![self.sum()?];
                    while self.peek() == Some(b',') {
                        self.i += 1;
                        a.push(self.sum()?);
                    }
                    self.expect(b')')?;
                    return match (name, a.len()) {
                        ("min", _) => Ok(a.iter().copied().fold(f64::INFINITY, f64::min)),
                        ("max", _) => Ok(a.iter().copied().fold(f64::NEG_INFINITY, f64::max)),
                        ("abs", 1) => Ok(a[0].abs()),
                        ("clamp", 3) => Ok(a[0].max(a[1]).min(a[2])),
                        _ => Err(format!("unknown function {name} with {} arguments", a.len())),
                    };
                }
                self.vars.iter().find(|(n, _)| n == name).map(|v| v.1).ok_or_else(|| format!("unknown parameter {name:?}"))
            }
            c => Err(format!("unexpected {:?}", c.map(|c| c as char))),
        }
    }
}

/// The value with every string `"=expr"` replaced by the expression's value.
pub fn substitute(v: &Value, vars: &[(String, f64)]) -> Result<Value, String> {
    Ok(match v {
        Value::String(s) if s.starts_with('=') => json!(eval_expr(&s[1..], vars)?),
        Value::Array(a) => Value::Array(a.iter().map(|x| substitute(x, vars)).collect::<Result<_, _>>()?),
        Value::Object(m) => Value::Object(m.iter().map(|(k, x)| Ok((k.clone(), substitute(x, vars)?))).collect::<Result<_, String>>()?),
        x => x.clone(),
    })
}

// -------------------------------------------------------------------------------------------- parameters

#[derive(Clone, Debug)]
pub struct Param {
    pub name: String,
    pub min: f64,
    pub max: f64,
    pub step: f64,
    pub start: f64,
}

impl Param {
    /// Number of grid values in [min, max].
    pub fn levels(&self) -> usize {
        ((self.max - self.min) / self.step + 1e-9).floor() as usize + 1
    }
    pub fn snap(&self, x: f64) -> f64 {
        let k = ((x - self.min) / self.step).round().clamp(0.0, (self.levels() - 1) as f64);
        self.at(k as usize)
    }
    pub fn at(&self, k: usize) -> f64 {
        ((self.min + k as f64 * self.step) * 1e9).round() / 1e9 + 0.0
    }
    fn index(&self, x: f64) -> usize {
        ((x - self.min) / self.step).round().max(0.0) as usize
    }
    fn unit(&self, x: f64) -> f64 {
        (x - self.min) / (self.max - self.min)
    }
}

pub fn parse_params(v: &Value) -> Result<Vec<Param>, String> {
    let a = v.as_array().ok_or("\"parameters\" must be a list")?;
    let mut ps = Vec::new();
    for p in a {
        let name = p["name"].as_str().ok_or("parameter needs \"name\"")?.to_string();
        let f = |k: &str| p[k].as_f64().ok_or_else(|| format!("parameter {name}: needs \"{k}\""));
        let (min, max, step) = (f("min")?, f("max")?, f("step")?);
        if !(max > min && step > 0.0) || !name.chars().all(|c| c.is_ascii_alphanumeric() || c == '_') {
            return Err(format!("parameter {name}: needs an identifier name, min < max and step > 0"));
        }
        let mut q = Param { name, min, max, step, start: 0.0 };
        q.start = q.snap(p["start"].as_f64().unwrap_or((min + max) / 2.0));
        ps.push(q);
    }
    if ps.is_empty() {
        return Err("no parameters".into());
    }
    Ok(ps)
}

/// Cache key of a point, e.g. `arm-1.0000_feed+0.3000`.
pub fn key(ps: &[Param], x: &[f64]) -> String {
    ps.iter().zip(x).map(|(p, v)| format!("{}{:+.4}", p.name, v + 0.0)).collect::<Vec<_>>().join("_")
}

fn vars(ps: &[Param], x: &[f64]) -> Vec<(String, f64)> {
    ps.iter().zip(x).map(|(p, v)| (p.name.clone(), *v)).collect()
}

// ------------------------------------------------------------------------------------- Gaussian process

/// Gaussian-process regression with a squared-exponential kernel on unit-scaled inputs.
#[derive(Clone, Debug)]
pub struct Gp {
    xs: Vec<Vec<f64>>,
    l: Vec<f64>,
    alpha: Vec<f64>,
    mean: f64,
    pub ls: Vec<f64>,
    pub sf2: f64,
    pub sn2: f64,
}

fn kern(a: &[f64], b: &[f64], ls: &[f64], sf2: f64) -> f64 {
    let mut s = 0.0;
    for i in 0..a.len() {
        let d = (a[i] - b[i]) / ls[i];
        s += d * d;
    }
    sf2 * (-0.5 * s).exp()
}

/// In-place lower Cholesky factor of an n x n row-major matrix; false if not positive definite.
fn cholesky(a: &mut [f64], n: usize) -> bool {
    for j in 0..n {
        let mut d = a[j * n + j];
        for k in 0..j {
            d -= a[j * n + k] * a[j * n + k];
        }
        if d <= 0.0 {
            return false;
        }
        let d = d.sqrt();
        a[j * n + j] = d;
        for i in j + 1..n {
            let mut s = a[i * n + j];
            for k in 0..j {
                s -= a[i * n + k] * a[j * n + k];
            }
            a[i * n + j] = s / d;
        }
        for k in j + 1..n {
            a[j * n + k] = 0.0;
        }
    }
    true
}

/// Solve L y = b.
fn forward(l: &[f64], n: usize, b: &[f64]) -> Vec<f64> {
    let mut y = b.to_vec();
    for i in 0..n {
        for k in 0..i {
            y[i] -= l[i * n + k] * y[k];
        }
        y[i] /= l[i * n + i];
    }
    y
}

/// Solve L^T x = y.
fn backward(l: &[f64], n: usize, y: &[f64]) -> Vec<f64> {
    let mut x = y.to_vec();
    for i in (0..n).rev() {
        for k in i + 1..n {
            x[i] -= l[k * n + i] * x[k];
        }
        x[i] /= l[i * n + i];
    }
    x
}

impl Gp {
    /// Fit with fixed hyperparameters; returns the model and its log marginal likelihood.
    pub fn fit_with(xs: &[Vec<f64>], ys: &[f64], mean: f64, ls: &[f64], sf2: f64, sn2: f64) -> Option<(Gp, f64)> {
        let n = xs.len();
        let mut k = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..=i {
                let v = kern(&xs[i], &xs[j], ls, sf2) + if i == j { sn2 + 1e-10 * sf2 } else { 0.0 };
                k[i * n + j] = v;
                k[j * n + i] = v;
            }
        }
        if !cholesky(&mut k, n) {
            return None;
        }
        let r: Vec<f64> = ys.iter().map(|y| y - mean).collect();
        let alpha = backward(&k, n, &forward(&k, n, &r));
        let lml = -0.5 * r.iter().zip(&alpha).map(|(a, b)| a * b).sum::<f64>() - (0..n).map(|i| k[i * n + i].ln()).sum::<f64>() - 0.5 * n as f64 * (2.0 * std::f64::consts::PI).ln();
        Some((Gp { xs: xs.to_vec(), l: k, alpha, mean, ls: ls.to_vec(), sf2, sn2 }, lml))
    }

    /// Maximum-likelihood fit over a grid of length scales (per input up to three inputs, shared beyond),
    /// signal variances and the given noise variances.
    pub fn fit(xs: &[Vec<f64>], ys: &[f64], noise: &[f64]) -> Gp {
        let d = xs[0].len();
        let n = ys.len() as f64;
        let mean = ys.iter().sum::<f64>() / n;
        let var = (ys.iter().map(|y| (y - mean).powi(2)).sum::<f64>() / n).max(1e-4);
        const LS: [f64; 5] = [0.1, 0.2, 0.4, 0.8, 1.6];
        let per_dim = d <= 3;
        let combos = if per_dim { LS.len().pow(d as u32) } else { LS.len() };
        let mut best: Option<(Gp, f64)> = None;
        for c in 0..combos {
            let ls: Vec<f64> = if per_dim { (0..d).map(|i| LS[(c / LS.len().pow(i as u32)) % LS.len()]).collect() } else { vec![LS[c]; d] };
            for sf in [0.5, 1.0, 2.0, 4.0] {
                for &sn2 in noise {
                    if let Some((g, lml)) = Gp::fit_with(xs, ys, mean, &ls, var * sf, sn2) {
                        if best.as_ref().is_none_or(|b| lml > b.1) {
                            best = Some((g, lml));
                        }
                    }
                }
            }
        }
        best.map(|b| b.0).unwrap_or_else(|| Gp::fit_with(xs, ys, mean, &vec![0.4; d], var, var).unwrap().0)
    }

    /// Posterior mean and variance at x.
    pub fn predict(&self, x: &[f64]) -> (f64, f64) {
        let n = self.xs.len();
        let ks: Vec<f64> = self.xs.iter().map(|xi| kern(xi, x, &self.ls, self.sf2)).collect();
        let mu = self.mean + ks.iter().zip(&self.alpha).map(|(a, b)| a * b).sum::<f64>();
        let v = forward(&self.l, n, &ks);
        (mu, (self.sf2 - v.iter().map(|a| a * a).sum::<f64>()).max(1e-12))
    }
}

fn norm_cdf(z: f64) -> f64 {
    // Abramowitz and Stegun 7.1.26 (|error| < 1.5e-7)
    let x = z.abs() / std::f64::consts::SQRT_2;
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let erf = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * (-x * x).exp();
    if z >= 0.0 { 0.5 * (1.0 + erf) } else { 0.5 * (1.0 - erf) }
}

/// Expected improvement over `best` (maximisation) of a normal prediction.
pub fn expected_improvement(mu: f64, var: f64, best: f64, xi: f64) -> f64 {
    let sd = var.sqrt();
    let d = mu - best - xi;
    if sd < 1e-9 {
        return d.max(0.0);
    }
    let z = d / sd;
    d * norm_cdf(z) + sd * (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt()
}

// -------------------------------------------------------------------------------------------- the search

pub trait Evaluator {
    /// The objective of point x at a fidelity level, per band frequency (dB, higher is better); the
    /// objective is the worst of them.
    fn eval(&mut self, x: &[f64], level: &str) -> Result<Eval, String>;
}

#[derive(Clone, Debug)]
pub struct Eval {
    /// the objective per band frequency (dB)
    pub per_f: Vec<f64>,
    /// cells x time steps of the solve (summed over excitations)
    pub cost: f64,
    pub seconds: f64,
    pub status: String,
    /// read from the cache instead of solved
    pub cached: bool,
}

fn worst(v: &[f64]) -> f64 {
    v.iter().copied().fold(f64::INFINITY, f64::min)
}

#[derive(Clone, Debug)]
pub struct Opts {
    pub method: String,
    /// screening solves
    pub budget: usize,
    /// initial design size (0: 2 d + 1)
    pub init: usize,
    pub promote_every: usize,
    pub final_n: usize,
    pub tie_db: f64,
    pub screen: String,
    pub full: String,
    pub seed: u64,
    /// initial trust-region half-width, as a fraction of each range
    pub tr0: f64,
    /// screening noise level for the model (dB)
    pub noise_db: f64,
    /// space mapping from screening to full fidelity (off: promote the best screening candidates as they are)
    pub mapping: bool,
}

impl Default for Opts {
    fn default() -> Opts {
        Opts { method: "surrogate".into(), budget: 16, init: 0, promote_every: 6, final_n: 3, tie_db: 1.2, screen: "screening".into(), full: "spec".into(), seed: 1, tr0: 0.5, noise_db: 0.05, mapping: true }
    }
}

impl Opts {
    pub fn from_json(v: &Value) -> Opts {
        let d = Opts::default();
        let n = |k: &str, dv: f64| v.get(k).and_then(|x| x.as_f64()).unwrap_or(dv);
        let s = |k: &str, dv: &str| v.get(k).and_then(|x| x.as_str()).unwrap_or(dv).to_string();
        let method = s("method", &d.method);
        Opts {
            mapping: v.get("mapping").and_then(|x| x.as_bool()).unwrap_or(method == "surrogate"),
            method,
            budget: n("budget", d.budget as f64) as usize,
            init: n("init", 0.0) as usize,
            promote_every: n("promote_every", d.promote_every as f64) as usize,
            final_n: n("final", d.final_n as f64) as usize,
            tie_db: n("tie_db", d.tie_db),
            screen: s("screening", &d.screen),
            full: s("full", &d.full),
            seed: n("seed", d.seed as f64) as u64,
            tr0: n("trust_region", d.tr0),
            noise_db: n("noise_db", d.noise_db),
        }
    }
}

struct Point {
    x: Vec<f64>,
    u: Vec<f64>,
    screen: Option<Vec<f64>>,
    full: Option<Vec<f64>>,
}

/// Screening models per band frequency and the space mapping to full fidelity:
/// full_f(u) ~ screen_f(u + shift) + offset (u in unit coordinates).
struct Model {
    gps: Vec<Gp>,
    shift: Vec<f64>,
    offset: f64,
    rms: f64,
}

impl Model {
    /// Predicted full-fidelity objective (worst band frequency) at u: mean and variance.
    fn predict(&self, u: &[f64]) -> (f64, f64) {
        let v: Vec<f64> = u.iter().zip(&self.shift).map(|(a, s)| a + s).collect();
        let mut best = (f64::INFINITY, 0.0);
        for g in &self.gps {
            let (m, var) = g.predict(&v);
            if m < best.0 {
                best = (m, var);
            }
        }
        (best.0 + self.offset, best.1)
    }
}

/// The shift and offset that best map the screening models onto the full-fidelity responses at the
/// promoted points, with a ridge penalty on the shift (so that one point does not pick an arbitrary one).
fn fit_mapping(gps: &[Gp], pairs: &[(Vec<f64>, Vec<f64>)], d: usize) -> (Vec<f64>, f64, f64) {
    if pairs.is_empty() {
        return (vec![0.0; d], 0.0, 0.0);
    }
    const LAMBDA: f64 = 1.0;
    const SMAX: f64 = 0.35;
    let cost = |s: &[f64]| -> (f64, f64, f64) {
        let mut r = Vec::new();
        for (u, full) in pairs {
            let v: Vec<f64> = u.iter().zip(s).map(|(a, b)| a + b).collect();
            for (g, y) in gps.iter().zip(full) {
                r.push(y - g.predict(&v).0);
            }
        }
        let b = r.iter().sum::<f64>() / r.len() as f64;
        let ss = r.iter().map(|x| (x - b).powi(2)).sum::<f64>();
        (ss + LAMBDA * s.iter().map(|x| x * x).sum::<f64>(), b, (ss / r.len() as f64).sqrt())
    };
    let mut s = vec![0.0; d];
    let mut c = cost(&s);
    let mut step = 0.1;
    while step > 0.004 {
        let mut moved = false;
        for j in 0..d {
            for sg in [1.0, -1.0] {
                let mut t = s.clone();
                t[j] = (t[j] + sg * step).clamp(-SMAX, SMAX);
                let ct = cost(&t);
                if ct.0 < c.0 - 1e-12 {
                    (s, c, moved) = (t, ct, true);
                }
            }
        }
        if !moved {
            step /= 2.0;
        }
    }
    (s, c.1, c.2)
}

struct Search<'a> {
    ps: &'a [Param],
    ev: &'a mut dyn Evaluator,
    o: &'a Opts,
    pts: Vec<Point>,
    hist: Vec<Value>,
    n_screen: usize,
    n_full: usize,
    cached: usize,
    cost: [f64; 2],
    secs: [f64; 2],
    trajectory: Vec<Value>,
    since_promotion: usize,
    model: Option<Model>,
    /// hyperparameters of the last full fit (length scales, signal variance factor, noise), and its age
    hyp: Option<(Vec<f64>, f64, f64)>,
    hyp_age: usize,
    tr: f64,
}

impl Search<'_> {
    fn snap(&self, x: &[f64]) -> Vec<f64> {
        self.ps.iter().zip(x).map(|(p, v)| p.snap(*v)).collect()
    }
    fn from_unit(&self, u: &[f64]) -> Vec<f64> {
        self.ps.iter().zip(u).map(|(p, v)| p.snap(p.min + v * (p.max - p.min))).collect()
    }
    fn unit(&self, x: &[f64]) -> Vec<f64> {
        self.ps.iter().zip(x).map(|(p, v)| p.unit(*v)).collect()
    }
    fn find(&self, x: &[f64]) -> Option<usize> {
        self.pts.iter().position(|p| p.x == x)
    }
    fn screened(&self, x: &[f64]) -> bool {
        self.find(x).is_some_and(|i| self.pts[i].screen.is_some())
    }
    fn fully(&self, x: &[f64]) -> bool {
        self.find(x).is_some_and(|i| self.pts[i].full.is_some())
    }
    fn point(&mut self, x: &[f64]) -> usize {
        if let Some(i) = self.find(x) {
            return i;
        }
        let u = self.unit(x);
        self.pts.push(Point { x: x.to_vec(), u, screen: None, full: None });
        self.pts.len() - 1
    }
    fn needs_model(&self) -> bool {
        self.o.mapping || self.o.method == "surrogate"
    }

    fn solve(&mut self, i: usize, full: bool, why: &str) -> Result<f64, String> {
        let level = if full { self.o.full.clone() } else { self.o.screen.clone() };
        let x = self.pts[i].x.clone();
        // with one level for both, a point solved at one is solved at the other
        let other = if full { &self.pts[i].screen } else { &self.pts[i].full };
        let e = match other {
            Some(v) if self.o.full == self.o.screen => Eval { per_f: v.clone(), cost: 0.0, seconds: 0.0, status: "same level".into(), cached: true },
            _ => self.ev.eval(&x, &level)?,
        };
        let value = worst(&e.per_f);
        let k = full as usize;
        if e.cached {
            self.cached += 1;
        } else {
            self.cost[k] += e.cost;
            self.secs[k] += e.seconds;
        }
        if full {
            self.n_full += 1;
            self.pts[i].full = Some(e.per_f.clone());
        } else {
            self.n_screen += 1;
            self.since_promotion += 1;
            self.pts[i].screen = Some(e.per_f.clone());
        }
        let best_full = self.pts.iter().filter_map(|p| p.full.as_ref().map(|v| worst(v))).fold(f64::NEG_INFINITY, f64::max);
        self.trajectory.push(json!([self.n_screen, self.n_full, if best_full.is_finite() { json!(best_full) } else { Value::Null }]));
        let xs: serde_json::Map<String, Value> = self.ps.iter().zip(&x).map(|(p, v)| (p.name.clone(), json!(v))).collect();
        self.hist.push(json!({"screening_solves": self.n_screen, "full_solves": self.n_full, "level": level, "x": xs, "value_db": value, "per_f_db": e.per_f,
                              "cost_cell_steps": e.cost, "seconds": e.seconds, "status": e.status, "cached": e.cached, "why": why}));
        eprintln!("optimise: [{} screening, {} full] {level:<10} {}: {value:+.3} dB ({why}{})", self.n_screen, self.n_full, key(self.ps, &x), if e.cached { ", cached" } else { "" });
        if self.needs_model() {
            self.refit(!full);
        }
        Ok(value)
    }

    /// Refit the screening models (when `gps`) and the mapping.
    fn refit(&mut self, gps: bool) {
        let idx: Vec<usize> = (0..self.pts.len()).filter(|&i| self.pts[i].screen.is_some()).collect();
        if idx.is_empty() {
            return;
        }
        let d = self.ps.len();
        if gps || self.model.is_none() {
            let xs: Vec<Vec<f64>> = idx.iter().map(|&i| self.pts[i].u.clone()).collect();
            let nf = self.pts[idx[0]].screen.as_ref().unwrap().len();
            let stats = |ys: &[f64]| {
                let m = ys.iter().sum::<f64>() / ys.len() as f64;
                (m, (ys.iter().map(|y| (y - m).powi(2)).sum::<f64>() / ys.len() as f64).max(1e-4))
            };
            // hyperparameters by likelihood on the worst value every few refits, shared by the frequencies
            if self.hyp.is_none() || self.hyp_age >= 3 {
                let ys: Vec<f64> = idx.iter().map(|&i| worst(self.pts[i].screen.as_ref().unwrap())).collect();
                let noise: Vec<f64> = [0.5, 1.0, 2.0].iter().map(|f| (f * self.o.noise_db).powi(2)).collect();
                let g = Gp::fit(&xs, &ys, &noise);
                self.hyp = Some((g.ls.clone(), g.sf2 / stats(&ys).1, g.sn2));
                self.hyp_age = 0;
            }
            self.hyp_age += 1;
            let (ls, sf, sn2) = self.hyp.clone().unwrap();
            let gps: Vec<Gp> = (0..nf)
                .map(|f| {
                    let ys: Vec<f64> = idx.iter().map(|&i| self.pts[i].screen.as_ref().unwrap()[f]).collect();
                    let (m, var) = stats(&ys);
                    Gp::fit_with(&xs, &ys, m, &ls, var * sf, sn2).or_else(|| Gp::fit_with(&xs, &ys, m, &ls, var * sf, sn2 + var)).unwrap().0
                })
                .collect();
            self.model = Some(Model { gps, shift: vec![0.0; d], offset: 0.0, rms: 0.0 });
        }
        if self.o.mapping {
            let pairs: Vec<(Vec<f64>, Vec<f64>)> = self.pts.iter().filter_map(|p| Some((p.u.clone(), p.full.clone()?))).filter(|_| self.o.full != self.o.screen).collect();
            let m = self.model.as_mut().unwrap();
            let (s, b, rms) = fit_mapping(&m.gps, &pairs, d);
            (m.shift, m.offset, m.rms) = (s, b, rms);
        }
    }

    /// The predicted full-fidelity objective of a known point: its full value, else the mapped model (with
    /// mapping), else its screening value.
    fn pred_point(&self, i: usize) -> f64 {
        let p = &self.pts[i];
        if let Some(v) = &p.full {
            return worst(v);
        }
        match (&self.model, self.o.mapping) {
            (Some(m), true) => m.predict(&p.u).0,
            _ => p.screen.as_ref().map(|v| worst(v)).unwrap_or(f64::NEG_INFINITY),
        }
    }
    fn incumbent(&self) -> usize {
        (0..self.pts.len()).max_by(|&a, &b| self.pred_point(a).total_cmp(&self.pred_point(b))).unwrap()
    }

    /// Grid points within the trust region (half-width tr, at least one step) around centre.
    fn region(&self, centre: &[f64], tr: f64, rng: &mut Rng) -> Vec<Vec<f64>> {
        let d = self.ps.len();
        let lo_hi: Vec<(usize, usize)> = self.ps.iter().zip(centre).map(|(p, c)| {
            let r = (tr * (p.max - p.min) / p.step + 1e-9).floor().max(1.0) as usize;
            let k = p.index(*c);
            (k.saturating_sub(r), (k + r).min(p.levels() - 1))
        }).collect();
        let n_box: f64 = lo_hi.iter().map(|(a, b)| (b - a + 1) as f64).product();
        if n_box > 30_000.0 {
            return (0..8000).map(|_| self.ps.iter().zip(&lo_hi).map(|(p, r)| p.at(r.0 + rng.below(r.1 - r.0 + 1))).collect()).collect();
        }
        let mut v = Vec::with_capacity(n_box as usize);
        let mut k: Vec<usize> = lo_hi.iter().map(|r| r.0).collect();
        'outer: loop {
            v.push(self.ps.iter().zip(&k).map(|(p, &i)| p.at(i)).collect());
            for j in 0..d {
                if k[j] < lo_hi[j].1 {
                    k[j] += 1;
                    continue 'outer;
                }
                k[j] = lo_hi[j].0;
            }
            return v;
        }
    }

    /// The point to promote: with mapping, the best predicted full-fidelity grid point near the incumbent
    /// (screened or not) that has no full solve; without, the best screened point that has none.
    fn promotion_candidate(&self, rng: &mut Rng) -> Option<(Vec<f64>, f64)> {
        match (&self.model, self.o.mapping) {
            (Some(m), true) => {
                let centre = self.pts[self.incumbent()].x.clone();
                let mut cands = self.region(&centre, self.tr.max(0.15), rng);
                for p in self.pts.iter().filter(|p| p.screen.is_some()) {
                    let u: Vec<f64> = p.u.iter().zip(&m.shift).map(|(a, s)| a - s).collect();
                    cands.push(self.from_unit(&u));
                    cands.push(p.x.clone());
                }
                cands.into_iter().filter(|x| !self.fully(x)).map(|x| {
                    let v = m.predict(&self.unit(&x)).0;
                    (x, v)
                }).max_by(|a, b| a.1.total_cmp(&b.1))
            }
            _ => self.pts.iter().filter(|p| p.full.is_none() && p.screen.is_some()).map(|p| (p.x.clone(), worst(p.screen.as_ref().unwrap()))).max_by(|a, b| a.1.total_cmp(&b.1)),
        }
    }
    fn promote(&mut self, why: &str, rng: &mut Rng) -> Result<bool, String> {
        let Some((x, pred)) = self.promotion_candidate(rng) else { return Ok(false) };
        let best_full = self.pts.iter().filter_map(|p| p.full.as_ref().map(|v| worst(v))).fold(f64::NEG_INFINITY, f64::max);
        if pred < best_full - self.o.tie_db {
            return Ok(false);
        }
        let i = self.point(&x);
        self.solve(i, true, &format!("{why}, predicted {pred:+.3} dB"))?;
        Ok(true)
    }
    fn screen(&mut self, x: &[f64], why: &str, rng: &mut Rng) -> Result<Option<f64>, String> {
        let i = self.point(x);
        if let Some(v) = &self.pts[i].screen {
            return Ok(Some(worst(v)));
        }
        if self.n_screen >= self.o.budget {
            return Ok(None);
        }
        let v = self.solve(i, false, why)?;
        if self.o.promote_every > 0 && self.since_promotion >= self.o.promote_every && self.n_screen < self.o.budget {
            self.since_promotion = 0;
            self.promote("promote: periodic", rng)?;
        }
        Ok(Some(v))
    }
    /// Final full solves: the best predicted candidates, while they are within tie_db of the best full value.
    fn finish(&mut self, rng: &mut Rng) -> Result<(), String> {
        for _ in 0..self.o.final_n {
            if !self.promote("promote: final", rng)? {
                break;
            }
        }
        Ok(())
    }

    fn surrogate(&mut self, rng: &mut Rng) -> Result<String, String> {
        let d = self.ps.len();
        let start: Vec<f64> = self.ps.iter().map(|p| p.start).collect();
        self.screen(&start, "start", rng)?;
        // initial design: a Latin hypercube on the grid
        let init = if self.o.init > 0 { self.o.init } else { 2 * d + 1 }.min(self.o.budget);
        let m = init.saturating_sub(1);
        let mut perms: Vec<Vec<usize>> = (0..d).map(|_| (0..m).collect()).collect();
        for p in perms.iter_mut() {
            for i in (1..m).rev() {
                p.swap(i, rng.below(i + 1));
            }
        }
        for k in 0..m {
            let u: Vec<f64> = (0..d).map(|j| (perms[j][k] as f64 + rng.unit()) / m as f64).collect();
            let x = self.from_unit(&u);
            self.screen(&x, "initial design", rng)?;
        }
        let total: f64 = self.ps.iter().map(|p| p.levels() as f64).product();
        let min_tr = self.ps.iter().map(|p| p.step / (p.max - p.min)).fold(0.0, f64::max);
        self.tr = self.o.tr0.max(min_tr);
        let mut fails = 0;
        while self.n_screen < self.o.budget {
            let inc = self.incumbent();
            let best = self.pred_point(inc);
            let centre = self.pts[inc].x.clone();
            let cands = self.region(&centre, self.tr, rng);
            let n_box = cands.len() as f64;
            let model = self.model.as_ref().unwrap();
            // expected improvement of the predicted full objective at u; the solve is at u's screening image
            let mut pick: Option<(Vec<f64>, f64)> = None;
            for x in cands {
                let u = self.unit(&x);
                let img = if self.o.mapping { self.from_unit(&u.iter().zip(&model.shift).map(|(a, s)| a + s).collect::<Vec<_>>()) } else { x.clone() };
                if self.screened(&img) {
                    continue;
                }
                let (mu, var) = model.predict(&u);
                let ei = expected_improvement(mu, var, best, 0.01);
                if pick.as_ref().is_none_or(|p| ei > p.1) {
                    pick = Some((img, ei));
                }
            }
            let Some((x, ei)) = pick else {
                if self.tr >= 1.0 || n_box >= total {
                    return Ok("every grid point screened".into());
                }
                self.tr = (self.tr * 2.0).min(1.0);
                continue;
            };
            if ei < 1e-4 && self.tr <= min_tr + 1e-12 {
                return Ok("converged: no expected improvement in the smallest trust region".into());
            }
            let tr = self.tr;
            if self.screen(&x, &format!("expected improvement {ei:.3} dB, trust region {tr:.2}"), rng)?.is_none() {
                break;
            }
            let now = self.pred_point(self.incumbent());
            if now > best + 0.01 {
                self.tr = (self.tr * 2.0).min(1.0);
                fails = 0;
            } else {
                fails += 1;
                if fails >= d.max(2) {
                    self.tr = (self.tr / 2.0).max(min_tr);
                    fails = 0;
                }
            }
        }
        Ok("budget".into())
    }

    fn pattern(&mut self, rng: &mut Rng) -> Result<String, String> {
        let d = self.ps.len();
        let mut x: Vec<f64> = self.ps.iter().map(|p| p.start).collect();
        let Some(mut fx) = self.screen(&x, "start", rng)? else { return Ok("budget".into()) };
        let mut step: Vec<usize> = self.ps.iter().map(|p| (p.levels() / 4).max(1)).collect();
        loop {
            let mut moved = false;
            'poll: for j in 0..d {
                for s in [1i64, -1] {
                    let k = self.ps[j].index(x[j]) as i64 + s * step[j] as i64;
                    if k < 0 || k >= self.ps[j].levels() as i64 {
                        continue;
                    }
                    let mut y = x.clone();
                    y[j] = self.ps[j].at(k as usize);
                    let why = format!("poll {} {:+}", self.ps[j].name, s * step[j] as i64);
                    let Some(fy) = self.screen(&y, &why, rng)? else { return Ok("budget".into()) };
                    if fy > fx {
                        (x, fx, moved) = (y, fy, true);
                        break 'poll;
                    }
                }
            }
            if !moved {
                if step.iter().all(|&s| s == 1) {
                    return Ok("converged: no improving neighbour at the finest step".into());
                }
                for s in step.iter_mut() {
                    *s = (*s / 2).max(1);
                }
            }
        }
    }
}

/// Run a search; returns the record written to optimise.json.
pub fn search(ps: &[Param], ev: &mut dyn Evaluator, o: &Opts) -> Result<Value, String> {
    let mut s = Search {
        ps, ev, o, pts: vec![], hist: vec![], n_screen: 0, n_full: 0, cached: 0, cost: [0.0; 2], secs: [0.0; 2], trajectory: vec![], since_promotion: 0,
        model: None, hyp: None, hyp_age: 0, tr: o.tr0,
    };
    let mut rng = Rng::new(o.seed);
    let stop = match o.method.as_str() {
        "surrogate" => s.surrogate(&mut rng)?,
        "pattern" => s.pattern(&mut rng)?,
        m => return Err(format!("unknown method {m:?}: use surrogate or pattern")),
    };
    s.finish(&mut rng)?;
    let named = |x: &[f64]| -> Value { s.ps.iter().zip(x).map(|(p, v)| (p.name.clone(), json!(v))).collect::<serde_json::Map<_, _>>().into() };
    let val = |v: &Option<Vec<f64>>| v.as_ref().map(|v| worst(v));
    let best_full = (0..s.pts.len()).filter(|&i| s.pts[i].full.is_some()).max_by(|&a, &b| val(&s.pts[a].full).unwrap().total_cmp(&val(&s.pts[b].full).unwrap()));
    let best_screen = (0..s.pts.len()).filter(|&i| s.pts[i].screen.is_some()).max_by(|&a, &b| val(&s.pts[a].screen).unwrap().total_cmp(&val(&s.pts[b].screen).unwrap()));
    let pt = |i: Option<usize>| i.map(|i| json!({"x": named(&s.pts[i].x), "full_db": val(&s.pts[i].full), "screening_db": val(&s.pts[i].screen), "full_per_f_db": s.pts[i].full, "screening_per_f_db": s.pts[i].screen})).unwrap_or(Value::Null);
    let pairs: Vec<f64> = s.pts.iter().filter_map(|p| Some(worst(p.full.as_ref()?) - worst(p.screen.as_ref()?))).collect();
    let mapping = match (&s.model, o.mapping) {
        (Some(m), true) => json!({"shift": s.ps.iter().zip(&m.shift).map(|(p, v)| (p.name.clone(), json!(v * (p.max - p.min)))).collect::<serde_json::Map<_, _>>(),
                                  "offset_db": m.offset, "rms_db": m.rms, "full_points": s.n_full}),
        _ => Value::Null,
    };
    Ok(json!({
        "method": o.method, "stop": stop,
        "parameters": s.ps.iter().map(|p| json!({"name": p.name, "min": p.min, "max": p.max, "step": p.step, "start": p.start})).collect::<Vec<_>>(),
        "opts": {"budget": o.budget, "init": o.init, "promote_every": o.promote_every, "final": o.final_n, "tie_db": o.tie_db,
                 "screening": o.screen, "full": o.full, "seed": o.seed, "trust_region": o.tr0, "noise_db": o.noise_db, "mapping": o.mapping},
        "screening_solves": s.n_screen, "full_solves": s.n_full, "cached": s.cached,
        "cell_steps": {"screening": s.cost[0], "full": s.cost[1], "total": s.cost[0] + s.cost[1]},
        "seconds": {"screening": s.secs[0], "full": s.secs[1], "total": s.secs[0] + s.secs[1]},
        "best": pt(best_full), "best_screening": pt(best_screen), "mapping": mapping,
        "full_minus_screening_db": pairs,
        "trajectory": s.trajectory, "evaluations": s.hist,
    }))
}

// ---------------------------------------------------------------------------------- the EM evaluator

/// Solves of the parameterised spec, cached per level and point.
pub struct SimEval<'a> {
    pub spec: Value,
    pub board: gerber::Board,
    pub kicad: bool,
    pub ps: Vec<Param>,
    pub cache: PathBuf,
    pub backend: sim::Backend<'a>,
}

fn eval_of(r: &Value, cached: bool) -> Result<Eval, String> {
    let per_f: Vec<f64> = r["objective"]["per_f"].as_array().map(|a| a.iter().filter_map(|p| p["total_db"].as_f64()).collect()).unwrap_or_default();
    if per_f.is_empty() {
        return Err("result has no objective (the spec needs \"objective\")".into());
    }
    let steps: f64 = r["runs"].as_array().map(|a| a.iter().filter_map(|x| x["steps"].as_f64()).sum()).unwrap_or(0.0);
    let seconds: f64 = r["runs"].as_array().map(|a| a.iter().filter_map(|x| x["seconds"].as_f64()).sum()).unwrap_or(0.0);
    let status = r["fidelity"]["status"].as_str().unwrap_or("spec").to_string();
    Ok(Eval { per_f, cost: r["cells"].as_f64().unwrap_or(0.0) * steps, seconds, status, cached })
}

impl Evaluator for SimEval<'_> {
    fn eval(&mut self, x: &[f64], level: &str) -> Result<Eval, String> {
        let dir = self.cache.join(level).join(key(&self.ps, x));
        if let Ok(t) = std::fs::read_to_string(dir.join("result.json")) {
            let r: Value = serde_json::from_str(&t).map_err(|e| format!("{}: {e}", dir.display()))?;
            return eval_of(&r, true);
        }
        let spec = fidelity::resolve(&substitute(&self.spec, &vars(&self.ps, x))?, level)?;
        let mut board = self.board.clone();
        sim::apply_edits(&mut board, &spec, self.kicad)?;
        std::fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
        std::fs::write(dir.join("spec.json"), serde_json::to_string_pretty(&spec).unwrap()).map_err(|e| e.to_string())?;
        let mut cx = sim::Ctx { spec, board, kicad: self.kicad, zs: vec![] };
        sim::run(&mut cx, &dir, self.backend)?;
        // score from the file just written, as a cached point is: serde_json's default number parser can be
        // one ulp off the value in memory, and a repeated search must give the same numbers bit for bit
        let t = std::fs::read_to_string(dir.join("result.json")).map_err(|e| format!("{}: {e}", dir.display()))?;
        eval_of(&serde_json::from_str(&t).map_err(|e| format!("{}: {e}", dir.display()))?, false)
    }
}

/// `fieldsolve optimise SPACE.json [--method surrogate|pattern] [--budget N] [--out DIR]`.
pub fn run(space_path: &Path, args: &[String], backend: sim::Backend) -> Result<Value, String> {
    let txt = std::fs::read_to_string(space_path).map_err(|e| format!("{}: {e}", space_path.display()))?;
    let space: Value = serde_json::from_str(&txt).map_err(|e| format!("{}: {e}", space_path.display()))?;
    let base = space_path.parent().unwrap_or(Path::new("."));
    let arg = |k: &str| args.iter().position(|a| a == k).and_then(|i| args.get(i + 1).cloned());
    let mut o = Opts::from_json(&space);
    if let Some(m) = arg("--method") {
        o.method = m;
    }
    if let Some(b) = arg("--budget").and_then(|v| v.parse().ok()) {
        o.budget = b;
    }
    let ps = parse_params(&space["parameters"])?;
    let spec_path = base.join(space["spec"].as_str().ok_or("space needs \"spec\"")?);
    let spec: Value = serde_json::from_str(&std::fs::read_to_string(&spec_path).map_err(|e| format!("{}: {e}", spec_path.display()))?).map_err(|e| format!("{}: {e}", spec_path.display()))?;
    let sdir = spec_path.parent().unwrap_or(Path::new("."));
    let board = gerber::load_dir(&sdir.join(spec["gerbers"].as_str().ok_or("spec needs \"gerbers\"")?))?;
    let kicad = spec.get("coords").and_then(|v| v.as_str()).unwrap_or("kicad") == "kicad";
    let out = arg("--out").map(PathBuf::from).unwrap_or_else(|| base.join(space["out"].as_str().unwrap_or("optimise")));
    let cache = space["cache"].as_str().map(|c| base.join(c)).unwrap_or_else(|| out.join("cache"));
    // every parameter must reach the spec; check the start point resolves before any solve
    let start: Vec<f64> = ps.iter().map(|p| p.start).collect();
    substitute(&spec, &vars(&ps, &start))?;
    let txt = serde_json::to_string(&spec).unwrap();
    for p in &ps {
        if !txt.contains(p.name.as_str()) {
            return Err(format!("parameter {} does not appear in the spec", p.name));
        }
    }
    let t0 = std::time::Instant::now();
    let mut ev = SimEval { spec, board, kicad, ps: ps.clone(), cache, backend };
    let mut r = search(&ps, &mut ev, &o)?;
    r["wall_seconds"] = json!(t0.elapsed().as_secs_f64());
    std::fs::create_dir_all(&out).map_err(|e| e.to_string())?;
    std::fs::write(out.join("optimise.json"), serde_json::to_string_pretty(&r).unwrap()).map_err(|e| e.to_string())?;
    eprintln!(
        "optimise: {} ({}): {} screening and {} full solves ({} cached), {:.3e} cell-steps, best full {} at {}",
        o.method, r["stop"].as_str().unwrap_or(""), r["screening_solves"], r["full_solves"], r["cached"], r["cell_steps"]["total"].as_f64().unwrap_or(0.0),
        r["best"]["full_db"], r["best"]["x"]
    );
    Ok(r)
}

// ------------------------------------------------------------------------- a two-fidelity test problem

/// A synthetic two-fidelity antenna objective on [0, 1]^d: worst in-band efficiency of a resonance whose
/// frequency moves with the parameters (strongly with the first, like a radiator length), with a match
/// factor on the others. The screening level shifts the resonance down by `shift` (the coarse-mesh error
/// measured on the MTH PIFA) and adds repeatable noise.
pub struct Synth {
    c: Vec<f64>,
    u0: Vec<f64>,
    w: Vec<f64>,
    pub shift: f64,
    pub noise: f64,
    pub ps: Vec<Param>,
}

impl Synth {
    pub fn new(d: usize, seed: u64) -> Synth {
        let mut r = Rng::new(seed.wrapping_mul(0x9E37_79B9).wrapping_add(17));
        let sign = |r: &mut Rng| if r.below(2) == 0 { 1.0 } else { -1.0 };
        let c = (0..d).map(|i| sign(&mut r) * if i == 0 { 0.05 + 0.04 * r.unit() } else { 0.01 + 0.02 * r.unit() }).collect();
        let u0 = (0..d).map(|_| 0.15 + 0.7 * r.unit()).collect();
        let w = (0..d).map(|i| if i == 0 { 0.0 } else { 0.4 + 0.8 * r.unit() }).collect();
        let ps = (0..d).map(|i| Param { name: format!("p{i}"), min: 0.0, max: 1.0, step: 0.025, start: 0.5 }).collect();
        Synth { c, u0, w, shift: 0.01, noise: 0.03, ps }
    }
    pub fn per_f(&self, u: &[f64], shift: f64) -> Vec<f64> {
        let du: Vec<f64> = u.iter().zip(&self.u0).map(|(a, b)| a - b).collect();
        let mut s: f64 = self.c.iter().zip(&du).map(|(c, d)| c * d).sum();
        if du.len() > 1 {
            s += 0.5 * self.c[0].abs() * du[0] * du[1];
        }
        let f0 = 2.44 * (1.0 + s - shift);
        let m = (1.0 - self.w.iter().zip(&du).map(|(w, d)| w * d * d).sum::<f64>()).max(0.05);
        [2.40, 2.44, 2.48].iter().map(|f| 10.0 * (m / (1.0 + ((f - f0) / 0.06).powi(2))).log10()).collect()
    }
    pub fn value(&self, u: &[f64], shift: f64) -> f64 {
        worst(&self.per_f(u, shift))
    }
    /// The best full-fidelity value on the parameter grid.
    pub fn optimum(&self) -> f64 {
        let n = self.ps[0].levels();
        let d = self.ps.len();
        let mut best = f64::NEG_INFINITY;
        let mut k = vec![0usize; d];
        loop {
            let u: Vec<f64> = k.iter().map(|&i| self.ps[0].at(i)).collect();
            best = best.max(self.value(&u, 0.0));
            let mut j = 0;
            while j < d {
                k[j] += 1;
                if k[j] < n {
                    break;
                }
                k[j] = 0;
                j += 1;
            }
            if j == d {
                return best;
            }
        }
    }
}

impl Evaluator for Synth {
    fn eval(&mut self, x: &[f64], level: &str) -> Result<Eval, String> {
        if level == "screening" {
            let mut h = Rng::new(x.iter().fold(1469598103934665603u64, |a, v| (a ^ (v * 1e6).round() as i64 as u64).wrapping_mul(1099511628211)));
            let g = (-2.0 * h.unit().max(1e-12).ln()).sqrt() * (2.0 * std::f64::consts::PI * h.unit()).cos();
            Ok(Eval { per_f: self.per_f(x, self.shift).iter().map(|v| v + self.noise * g).collect(), cost: 1.0 / 3.1, seconds: 0.0, status: "screening".into(), cached: false })
        } else {
            Ok(Eval { per_f: self.per_f(x, 0.0), cost: 1.0, seconds: 0.0, status: level.into(), cached: false })
        }
    }
}

/// `fieldsolve optimise-bench [--seeds N]`: both methods and a full-fidelity-only search on the synthetic
/// problem, at several budgets. Cost is in full-solve equivalents (a screening solve = 1 / 3.1, the
/// measured MTH cells x steps saving).
pub fn bench(seeds: u64, base: &Opts) -> Value {
    let mut rows = Vec::new();
    println!("{:<28}{:>7}{:>13}{:>10}{:>11}{:>10}{:>9}{:>12}", "method", "budget", "regret mean", "median", "< 0.1 dB", "screen", "full", "cost (full)");
    // (label, method, screening level, mapping, budget)
    let cases: Vec<(&str, &str, &str, bool, usize)> = vec![
        ("surrogate + mapping", "surrogate", "screening", true, 12), ("surrogate, no mapping", "surrogate", "screening", false, 12), ("pattern", "pattern", "screening", false, 12),
        ("surrogate + mapping", "surrogate", "screening", true, 16), ("surrogate, no mapping", "surrogate", "screening", false, 16), ("pattern", "pattern", "screening", false, 16),
        ("pattern + mapping", "pattern", "screening", true, 16),
        ("surrogate + mapping", "surrogate", "screening", true, 24), ("surrogate, no mapping", "surrogate", "screening", false, 24), ("pattern", "pattern", "screening", false, 24),
        ("surrogate + mapping", "surrogate", "screening", true, 40), ("pattern", "pattern", "screening", false, 40),
        ("surrogate, full only", "surrogate", "spec", false, 8), ("surrogate, full only", "surrogate", "spec", false, 12), ("surrogate, full only", "surrogate", "spec", false, 16),
        ("pattern, full only", "pattern", "spec", false, 16), ("pattern, full only", "pattern", "spec", false, 24),
    ];
    for (name, method, screen, mapping, budget) in cases {
        let mut regret = Vec::new();
        let (mut ns, mut nf, mut cost) = (0.0, 0.0, 0.0);
        for seed in 0..seeds {
            let mut p = Synth::new(3, seed);
            let opt = p.optimum();
            let ps = p.ps.clone();
            let o = Opts { method: method.into(), budget, screen: screen.into(), mapping, seed: seed + 1, ..base.clone() };
            let r = search(&ps, &mut p, &o).unwrap();
            regret.push(opt - r["best"]["full_db"].as_f64().unwrap());
            let full_only = screen == "spec";
            ns += if full_only { 0.0 } else { r["screening_solves"].as_f64().unwrap() };
            nf += if full_only { r["screening_solves"].as_f64().unwrap() } else { r["full_solves"].as_f64().unwrap() };
            cost += r["cell_steps"]["total"].as_f64().unwrap();
        }
        let n = seeds as f64;
        regret.sort_by(f64::total_cmp);
        let mean = regret.iter().sum::<f64>() / n;
        let within = regret.iter().filter(|&&r| r < 0.1).count() as f64 / n;
        println!("{name:<28}{budget:>7}{mean:>13.3}{:>10.3}{:>10.0}%{:>10.1}{:>9.1}{:>12.2}", regret[regret.len() / 2], within * 100.0, ns / n, nf / n, cost / n);
        rows.push(json!({"method": name, "budget": budget, "regret_mean_db": mean, "regret_median_db": regret[regret.len() / 2], "within_0p1_db": within,
                         "screening_solves": ns / n, "full_solves": nf / n, "cost_full_equivalents": cost / n}));
    }
    json!({"seeds": seeds, "promote_every": base.promote_every, "final": base.final_n, "rows": rows})
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(pairs: &[(&str, f64)]) -> Vec<(String, f64)> {
        pairs.iter().map(|(n, x)| (n.to_string(), *x)).collect()
    }

    #[test]
    fn expressions_follow_arithmetic_rules() {
        let e = |s: &str| eval_expr(s, &v(&[("arm", -1.0), ("feed", 0.3)])).unwrap();
        assert_eq!(e("1 + 2 * 3"), 7.0);
        assert_eq!(e("(1 + 2) * 3"), 9.0);
        assert_eq!(e("-(2 - 5)"), 3.0);
        assert_eq!(e("2e-1 * 10"), 2.0);
        assert_eq!(e("8 / 4 / 2"), 1.0);
        assert!((e("141.9 - max(arm, 0)") - 141.9).abs() < 1e-12);
        assert!((e("141.9 + max(-arm, 0)") - 142.9).abs() < 1e-12);
        assert!((e("153.95 + feed") - 154.25).abs() < 1e-12);
        assert_eq!(e("min(3, 1, 2)"), 1.0);
        assert_eq!(e("abs(arm)"), 1.0);
        assert_eq!(e("clamp(5, 0, 2)"), 2.0);
        for bad in ["arm +", "x * 2", "1 2", "sqrt(2)", "1 / 0", "(1"] {
            assert!(eval_expr(bad, &v(&[("arm", 1.0)])).is_err(), "{bad} should fail");
        }
        let spec = json!({"edits": [{"op": "clear", "rect": ["=141.9", 74.5, "=141.9 + max(-arm, 0)", 75.1], "layer": "F.Cu"}], "note": "x=1"});
        let s = substitute(&spec, &v(&[("arm", -0.5)])).unwrap();
        assert_eq!(s["edits"][0]["rect"], json!([141.9, 74.5, 142.4, 75.1]));
        assert_eq!(s["note"], "x=1", "only strings starting with = are expressions");
        assert_eq!(s["edits"][0]["layer"], "F.Cu");
    }

    #[test]
    fn parameters_snap_to_their_grid() {
        let ps = parse_params(&json!([{"name": "arm", "min": -1.9, "max": 1.5, "step": 0.1, "start": 0.02}, {"name": "clr", "min": -0.5, "max": 1.0, "step": 0.1}])).unwrap();
        assert_eq!(ps[0].levels(), 35);
        assert_eq!(ps[0].start, 0.0);
        assert_eq!(ps[0].snap(-3.0), -1.9);
        assert_eq!(ps[0].snap(1.46), 1.5);
        assert_eq!(ps[0].snap(-0.04), 0.0);
        assert_eq!(ps[1].start, 0.3, "default start: the middle, snapped");
        assert_eq!(key(&ps, &[ps[0].snap(-1.0), 0.0]), "arm-1.0000_clr+0.0000");
        assert!(parse_params(&json!([{"name": "a b", "min": 0, "max": 1, "step": 0.1}])).is_err());
        assert!(parse_params(&json!([{"name": "a", "min": 1, "max": 1, "step": 0.1}])).is_err());
    }

    /// The model interpolates a smooth function, is confident at its samples and not far from them.
    #[test]
    fn gaussian_process_interpolates_with_honest_variance() {
        let f = |x: f64| (6.0 * x).sin() - 0.5 * x;
        let xs: Vec<Vec<f64>> = (0..9).map(|i| vec![i as f64 / 8.0]).collect();
        let ys: Vec<f64> = xs.iter().map(|x| f(x[0])).collect();
        let g = Gp::fit(&xs, &ys, &[1e-6, 1e-4]);
        for (x, y) in xs.iter().zip(&ys) {
            assert!((g.predict(x).0 - y).abs() < 0.01);
        }
        for i in 0..8 {
            let x = (i as f64 + 0.5) / 8.0;
            assert!((g.predict(&[x]).0 - f(x)).abs() < 0.05, "midpoint {x}");
        }
        assert!(g.predict(&[2.0]).1 > 100.0 * g.predict(&[0.5]).1, "uncertain away from the data");
        assert!(expected_improvement(1.0, 0.0, 0.5, 0.0) == 0.5 && expected_improvement(0.0, 0.0, 0.5, 0.0) == 0.0);
        assert!(expected_improvement(0.0, 1.0, 0.5, 0.0) > 0.1, "uncertainty alone has expected improvement");
    }

    /// On the synthetic problem the surrogate search finds the full-fidelity optimum closely within its
    /// budget, promotes only as configured, and keeps a consistent record.
    #[test]
    fn surrogate_search_finds_the_full_fidelity_optimum() {
        let mut regrets = Vec::new();
        for seed in 0..12 {
            let mut p = Synth::new(3, seed);
            let opt = p.optimum();
            let ps = p.ps.clone();
            let o = Opts { budget: 24, ..Default::default() };
            let r = search(&ps, &mut p, &o).unwrap();
            assert!(r["screening_solves"].as_u64().unwrap() <= 24);
            let nf = r["full_solves"].as_u64().unwrap();
            assert!((1..=24 / 6 + 2).contains(&nf), "{nf} full solves");
            assert!(r["mapping"]["shift"].is_object());
            let ev = r["evaluations"].as_array().unwrap();
            assert_eq!(ev.len() as u64, r["screening_solves"].as_u64().unwrap() + nf);
            // full solves are promotions
            for e in ev.iter().filter(|e| e["level"] == "spec") {
                assert!(e["why"].as_str().unwrap().starts_with("promote"), "{e}");
            }
            let best = r["best"]["full_db"].as_f64().unwrap();
            assert!(best <= opt + 1e-9);
            // the trajectory's best full value never decreases
            let t: Vec<f64> = r["trajectory"].as_array().unwrap().iter().filter_map(|x| x[2].as_f64()).collect();
            assert!(t.windows(2).all(|w| w[1] >= w[0]));
            regrets.push(opt - best);
        }
        let mean = regrets.iter().sum::<f64>() / regrets.len() as f64;
        assert!(mean < 0.25, "mean regret {mean:.3} dB: {regrets:?}");
    }

    /// The compass search stays on the grid, inside the bounds, and stops when no neighbour improves.
    #[test]
    fn pattern_search_converges_on_a_smooth_problem() {
        let mut p = Synth::new(2, 3);
        p.noise = 0.0;
        p.shift = 0.0;
        let ps = p.ps.clone();
        let o = Opts { method: "pattern".into(), budget: 400, final_n: 1, mapping: false, ..Default::default() };
        let r = search(&ps, &mut p, &o).unwrap();
        assert!(r["stop"].as_str().unwrap().starts_with("converged"), "{}", r["stop"]);
        for e in r["evaluations"].as_array().unwrap() {
            for x in e["x"].as_object().unwrap().values() {
                let x = x.as_f64().unwrap();
                assert!((0.0..=1.0).contains(&x) && ((x / 0.025).round() * 0.025 - x).abs() < 1e-9);
            }
        }
        // a local optimum of a noiseless smooth function: no grid neighbour is better
        let b = &r["best"]["x"];
        let x = [b["p0"].as_f64().unwrap(), b["p1"].as_f64().unwrap()];
        let f = p.value(&x, 0.0);
        for (j, s) in [(0, 0.025), (0, -0.025), (1, 0.025), (1, -0.025)] {
            let mut y = x;
            y[j] = (y[j] + s).clamp(0.0, 1.0);
            assert!(p.value(&y, 0.0) <= f + 1e-12);
        }
    }

    /// A solve that is already in the cache is not repeated, and with one level for screening and full a
    /// promotion costs nothing.
    #[test]
    fn repeated_points_are_not_solved_twice() {
        struct Count(Synth, usize);
        impl Evaluator for Count {
            fn eval(&mut self, x: &[f64], level: &str) -> Result<Eval, String> {
                self.1 += 1;
                self.0.eval(x, level)
            }
        }
        let mut c = Count(Synth::new(2, 5), 0);
        let ps = c.0.ps.clone();
        let o = Opts { screen: "spec".into(), full: "spec".into(), budget: 10, mapping: false, ..Default::default() };
        let r = search(&ps, &mut c, &o).unwrap();
        assert_eq!(c.1, 10, "promotions reuse the screening solve");
        assert_eq!(r["cached"].as_u64().unwrap(), r["full_solves"].as_u64().unwrap());
    }
}
