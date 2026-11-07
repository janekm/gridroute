//! `fieldsolve tune CONFIG.json`: retune the parts on an EM extraction's ports without new EM runs.
//!
//! The EM result (`sparams.sNp`, `result.json`, `farfield.json`) is loaded once. Every allowed part's
//! reflection is cached at the objective frequencies, and the per-port embedded patterns are reduced to one
//! Hermitian radiation matrix per frequency,
//!
//!   R_ij = ∮ (conj Eθ_i Eθ_j + conj Eφ_i Eφ_j) / (2 η0) dΩ,   P_rad(a) = Re(a^H R a),
//!
//! on the far field's own quadrature, so P_rad is algebraically the same as `Pattern::combine(..).p_rad()`
//! (the off-diagonal terms are the ports' interference; adding port powers would be wrong). A candidate then
//! costs one small network solve per frequency and a quadratic form: no pattern is reconstructed.
//!
//! The score is the one `tools/ladder_retune.py` used: per band frequency
//! (1 − |Γ_target|²) × P_rad / P_in, where Γ_target is the power-wave reflection of the impedance seen at the
//! driven port against the `target` network's (else 1 − |S11|² against z0), and P_rad / P_in the share of the
//! power entering that port that is radiated; the objective is the worst of the band, in dB. Optional
//! harmonic constraints bound the radiated power per watt incident at chosen frequencies. Frequency
//! conventions follow `cosim` exactly (S at the EM grid point nearest each frequency; the match from parts
//! at that grid frequency, the radiation from parts at the exact frequency), so the scores equal those
//! computed from `cosim` output.
//!
//! The search is coordinate descent from the as-built parts and from random starts, plus every one- and
//! two-part change from as built. Only the selected candidates get detailed output: they are written as a
//! `cosim` config and run through `cosim` (Touchstone, combined patterns), and the evaluator's scores are
//! checked against that output.
//!
//! ```json
//! { "em": "mth7_board/sparams.s7p", "out": "tune_mth7", "drive": "CHIP",
//!   "target": "nordic_ref_network.s2p", "band": [2.40e9, 2.44e9, 2.48e9],
//!   "harmonics": [{"f": 4.88e9, "max_db": "base"}, {"f": 7.32e9, "max_db": -20}],
//!   "sets": {"caps": {"template": "murata_sparams/lib_grm/GRM1555C1H{code}{tol}A01.s2p", "tols": "BCWDJ",
//!                     "unit": "p", "values": [0.2, 0.3, 1.0, 8.2]}},
//!   "positions": {"C59": ["open", "@caps"], "R70": ["short", {"set": "caps", "min": "0.5p", "suffix": " ser"}],
//!                 "C20": ["open", {"set": "caps", "max": "2p"}, {"0.4 nH": "L:0.4n"}]},
//!   "fixed": {"L4": "L:4.7n"},
//!   "base": {"C59": "1p", "R70": "short", "C20": "open"},
//!   "search": {"restarts": 60, "rounds": 8, "exhaustive": 2, "seed": 7}, "top": 8 }
//! ```

use crate::cosim::{read_touchstone, si, Part};
use crate::nf2ff::Pattern;
use crate::post::Cx;
use rayon::prelude::*;
use serde_json::{json, Value};
use std::path::{Path, PathBuf};
use std::time::Instant;

/// Most EM ports a candidate solve handles (stack buffers).
pub const MAXP: usize = 16;

/// R_ij over the patterns' sphere grid, with `Pattern::p_rad`'s quadrature. Row-major n x n, Hermitian.
pub fn radiation_matrix(pats: &[&Pattern]) -> Vec<Cx> {
    let n = pats.len();
    let p0 = pats[0];
    let (nt, np) = (p0.theta.len(), p0.phi.len());
    let mut r = vec![Cx::default(); n * n];
    let dph = 2.0 * std::f64::consts::PI / np as f64;
    for it in 0..nt {
        let w = if it == 0 || it + 1 == nt { 0.5 } else { 1.0 } * (p0.theta[nt - 1] - p0.theta[0]) / (nt - 1) as f64 * p0.theta[it].sin() * dph / (2.0 * crate::model::ETA0);
        for ip in 0..np {
            let q = it * np + ip;
            for i in 0..n {
                let (ti, pi) = (pats[i].eth[q].conj() * w, pats[i].eph[q].conj() * w);
                for j in i..n {
                    let v = ti * pats[j].eth[q] + pi * pats[j].eph[q];
                    r[i * n + j] = r[i * n + j] + v;
                }
            }
        }
    }
    for i in 0..n {
        r[i * n + i] = Cx::new(r[i * n + i].re, 0.0);
        for j in 0..i {
            r[i * n + j] = r[j * n + i].conj();
        }
    }
    r
}

/// Re(a^H R a).
#[inline]
pub fn quad_power(a: &[Cx], r: &[Cx]) -> f64 {
    let n = a.len();
    let mut p = 0.0;
    for i in 0..n {
        let mut row = Cx::default();
        for j in 0..n {
            row = row + r[i * n + j] * a[j];
        }
        p += (a[i].conj() * row).re;
    }
    p
}

/// Murata-style RKM code: 0.2 -> "R20", 4.7 -> "4N7" (sep 'N'), 10 -> "10N".
pub fn rkm(v: f64, sep: char) -> String {
    if v < 1.0 {
        format!("{sep}{:02}", (v * 100.0).round() as i64)
    } else if v >= 10.0 {
        format!("{}{sep}", v.round() as i64)
    } else {
        let t = (v * 10.0).round() as i64;
        format!("{}{sep}{}", t / 10, t % 10)
    }
}

/// A part spec with vendor files resolved against `base`.
fn resolve_spec(base: &Path, spec: &str) -> String {
    for pre in ["s2p:", "s2p-shunt:", "s1p:"] {
        if let Some(f) = spec.strip_prefix(pre) {
            if Path::new(f).is_relative() {
                return format!("{pre}{}", base.join(f).display());
            }
        }
    }
    spec.to_string()
}

/// One position's allowed parts: (label, spec).
pub type Choices = Vec<(String, String)>;

/// Expand a position's choice list: "open", a part spec, "@set", {"set", "min", "max", "suffix"}, or
/// {"label": "spec", ...}.
pub fn expand_choices(items: &[Value], sets: &serde_json::Map<String, Value>, base: &Path) -> Result<Choices, String> {
    let mut out: Choices = Vec::new();
    let set_members = |name: &str, lo: f64, hi: f64, suffix: &str, out: &mut Choices| -> Result<(), String> {
        let s = sets.get(name).ok_or_else(|| format!("no set {name:?}"))?;
        let unit = s.get("unit").and_then(|v| v.as_str()).unwrap_or("");
        let mul = if unit.is_empty() { 1.0 } else { si(&format!("1{unit}"))? };
        let sep = s.get("code_sep").and_then(|v| v.as_str()).and_then(|c| c.chars().next()).unwrap_or(if unit == "n" { 'N' } else { 'R' });
        let tmpl = s.get("template").and_then(|v| v.as_str()).ok_or_else(|| format!("set {name:?} needs \"template\""))?;
        let tols: Vec<String> = s.get("tols").and_then(|v| v.as_str()).map(|t| t.chars().map(String::from).collect()).unwrap_or_else(|| vec![String::new()]);
        let kind = s.get("kind").and_then(|v| v.as_str()).unwrap_or("s2p");
        for v in s.get("values").and_then(|v| v.as_array()).ok_or_else(|| format!("set {name:?} needs \"values\""))?.iter().filter_map(|x| x.as_f64()) {
            if v * mul < lo * (1.0 - 1e-9) || v * mul > hi * (1.0 + 1e-9) {
                continue;
            }
            let code = rkm(v, sep);
            // the first tolerance grade with a file, as tools/ladder_retune.py
            if let Some(f) = tols.iter().map(|t| base.join(tmpl.replace("{code}", &code).replace("{tol}", t))).find(|f| f.exists()) {
                out.push((format!("{v}{unit}{suffix}"), format!("{kind}:{}", f.display())));
            }
        }
        Ok(())
    };
    for it in items {
        match it {
            Value::String(s) if s.starts_with('@') => set_members(&s[1..], f64::NEG_INFINITY, f64::INFINITY, "", &mut out)?,
            Value::String(s) => out.push((s.clone(), resolve_spec(base, s))),
            Value::Object(o) if o.contains_key("set") => {
                let lim = |k: &str, d: f64| o.get(k).and_then(|v| v.as_str()).map(si).transpose().map(|x| x.unwrap_or(d));
                set_members(o["set"].as_str().unwrap_or(""), lim("min", f64::NEG_INFINITY)?, lim("max", f64::INFINITY)?, o.get("suffix").and_then(|v| v.as_str()).unwrap_or(""), &mut out)?;
            }
            Value::Object(o) => {
                for (k, v) in o {
                    out.push((k.clone(), resolve_spec(base, v.as_str().ok_or("a choice map is label -> part spec")?)));
                }
            }
            _ => return Err(format!("bad choice {it}")),
        }
    }
    for i in 0..out.len() {
        if out[..i].iter().any(|x| x.0 == out[i].0) {
            return Err(format!("duplicate choice label {:?}", out[i].0));
        }
    }
    Ok(out)
}

/// One network solve: the EM S-matrix at a grid point and the terminations' reflections at one frequency.
struct Point {
    /// EM S-matrix, row-major n x n
    s: Vec<Cx>,
    /// reflection per terminated slot: fixed parts' value (searched slots are looked up per choice)
    fixed: Vec<Cx>,
    /// reflections of every choice of every position, `off[p] + c`
    gam: Vec<Cx>,
}

/// The evaluator's view of a candidate: one choice index per position.
pub type Cand = Vec<u16>;

/// Score of one candidate: worst band value (dB) and the total harmonic-constraint violation (dB).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Score {
    pub worst: f64,
    pub violation: f64,
}

impl Score {
    /// One number to rank by: feasible candidates by their worst band value, infeasible ones below all of
    /// them by their violation.
    pub fn rank(&self) -> f64 {
        if self.violation <= 0.0 {
            self.worst
        } else {
            -1e6 - self.violation
        }
    }
}

/// Per-frequency detail of one candidate.
#[derive(Clone, Debug)]
pub struct Detail {
    /// per band frequency: (score dB, |Γ| against the target, radiated / accepted dB, P_in, P_rad)
    pub band: Vec<(f64, f64, f64, f64, f64)>,
    /// per harmonic: radiated power per watt incident, dB
    pub harm: Vec<f64>,
    pub score: Score,
}

pub struct Evaluator {
    pub names: Vec<String>,
    n: usize,
    drive: usize,
    /// terminated EM ports (all but the driven one), and which position (if any) each one is
    term: Vec<usize>,
    slot_pos: Vec<Option<usize>>,
    pub positions: Vec<String>,
    pub choices: Vec<Choices>,
    off: Vec<usize>,
    pub band: Vec<f64>,
    /// per band frequency: (match point, radiation point); per harmonic: radiation point
    band_pts: Vec<(usize, usize)>,
    harm_pts: Vec<usize>,
    pts: Vec<Point>,
    /// radiation matrix per radiation point (n x n); empty for match-only points
    rmat: Vec<Vec<Cx>>,
    /// target impedance at each band frequency's grid point
    target: Vec<Option<Cx>>,
    pub harm_freqs: Vec<f64>,
    pub harm_limit: Vec<f64>,
    z0: f64,
}

/// What `Evaluator::load` needs besides the EM result.
#[derive(Clone)]
pub struct Setup {
    pub em: PathBuf,
    pub drive: String,
    pub target: Option<PathBuf>,
    pub band: Vec<f64>,
    /// (frequency, limit dB or None for "as built")
    pub harmonics: Vec<(f64, Option<f64>)>,
    pub positions: Vec<(String, Choices)>,
    /// other terminated ports: (name, spec)
    pub fixed: Vec<(String, String)>,
    pub deembed: bool,
}

/// Stage timings, seconds.
#[derive(Default, Debug, Clone)]
pub struct Timings {
    pub load: f64,
    pub parts: f64,
    pub radiation: f64,
}

impl Evaluator {
    pub fn load(st: &Setup) -> Result<(Evaluator, Timings), String> {
        let mut tm = Timings::default();
        let t = Instant::now();
        let em = read_touchstone(&st.em)?;
        let res: Value = serde_json::from_str(&std::fs::read_to_string(st.em.with_file_name("result.json")).map_err(|e| format!("result.json: {e}"))?).map_err(|e| e.to_string())?;
        let names: Vec<String> = res["port_names"].as_array().ok_or("result.json has no port_names")?.iter().map(|v| v.as_str().unwrap_or("").to_string()).collect();
        let n = em.s.len();
        if names.len() != n {
            return Err(format!("{} port names for a {n}-port file", names.len()));
        }
        if n > MAXP {
            return Err(format!("{n} ports: the evaluator handles at most {MAXP}"));
        }
        let (ffreq, fnames, fpats) = crate::nf2ff::read_farfield(&st.em.with_file_name("farfield.json"))?;
        let target_ts = st.target.as_ref().map(|p| read_touchstone(p)).transpose()?;
        tm.load = t.elapsed().as_secs_f64();

        let t = Instant::now();
        let idx = |nm: &str| names.iter().position(|x| x == nm).ok_or_else(|| format!("no EM port {nm:?} (have {names:?})"));
        let drive = idx(&st.drive)?;
        let term: Vec<usize> = (0..n).filter(|&k| k != drive).collect();
        let bar_l = |nm: &str| -> f64 {
            res["port_bars"].as_array().and_then(|a| a.iter().find(|b| b["name"].as_str() == Some(nm))).and_then(|b| b["l_gap_nh"].as_f64()).unwrap_or(0.0) * 1e-9
        };
        let part = |nm: &str, spec: &str| -> Result<Part, String> {
            let p = Part::parse(spec)?;
            Ok(if st.deembed && bar_l(nm) != 0.0 { Part::Series { inner: Box::new(p), l: -bar_l(nm) } } else { p })
        };
        let mut slot_pos = vec![None; term.len()];
        let mut fixed_parts: Vec<Option<Part>> = vec![None; term.len()];
        let mut positions = Vec::new();
        let mut choices = Vec::new();
        let mut parts: Vec<Vec<Part>> = Vec::new();
        for (pi, (pn, ch)) in st.positions.iter().enumerate() {
            let k = idx(pn)?;
            let slot = term.iter().position(|&x| x == k).ok_or_else(|| format!("{pn} is the driven port"))?;
            if slot_pos[slot].is_some() {
                return Err(format!("position {pn} given twice"));
            }
            if ch.is_empty() || ch.len() > u16::MAX as usize {
                return Err(format!("position {pn}: {} choices", ch.len()));
            }
            slot_pos[slot] = Some(pi);
            positions.push(pn.clone());
            parts.push(ch.iter().map(|(_, s)| part(pn, s)).collect::<Result<_, _>>()?);
            choices.push(ch.clone());
        }
        for (pn, spec) in &st.fixed {
            let k = idx(pn)?;
            let slot = term.iter().position(|&x| x == k).ok_or_else(|| format!("{pn} is the driven port"))?;
            if slot_pos[slot].is_some() || fixed_parts[slot].is_some() {
                return Err(format!("port {pn} has both a fixed part and a position"));
            }
            fixed_parts[slot] = Some(part(pn, spec)?);
        }
        for (slot, &k) in term.iter().enumerate() {
            if slot_pos[slot].is_none() && fixed_parts[slot].is_none() {
                return Err(format!("EM port {:?} has no part (searched or fixed; use \"open\" to leave it open)", names[k]));
            }
        }
        let mut off = vec![0usize];
        for c in &choices {
            off.push(off.last().unwrap() + c.len());
        }
        let freqs = &em.freqs;
        let at = |f: f64| (0..freqs.len()).min_by(|&a, &b| (freqs[a] - f).abs().partial_cmp(&(freqs[b] - f).abs()).unwrap()).unwrap();
        let mut pts = Vec::new();
        let mut rmat = Vec::new();
        let mut rad_pats: Vec<(usize, usize)> = Vec::new(); // (point, far-field frequency index)
        let mut point = |f_parts: f64, fi: usize, far_q: Option<usize>, pts: &mut Vec<Point>| -> usize {
            let s = (0..n * n).map(|q| em.s[q / n][q % n][fi]).collect();
            let fixed = fixed_parts.iter().map(|p| p.as_ref().map(|p| gamma(p, f_parts, em.z0)).unwrap_or_default()).collect();
            let gam = parts.iter().flat_map(|ps| ps.iter().map(|p| gamma(p, f_parts, em.z0))).collect();
            pts.push(Point { s, fixed, gam });
            if let Some(q) = far_q {
                rad_pats.push((pts.len() - 1, q));
            }
            pts.len() - 1
        };
        let far_q = |f: f64| ffreq.iter().position(|&x| (x - f).abs() <= 1e-6 * f).ok_or_else(|| format!("farfield.json has no pattern at {f:.4e} Hz (have {ffreq:?})"));
        let mut band_pts = Vec::new();
        let mut target = Vec::new();
        for &f in &st.band {
            let fi = at(f);
            let m = point(freqs[fi], fi, None, &mut pts);
            let r = point(f, fi, Some(far_q(f)?), &mut pts);
            band_pts.push((m, r));
            // the target network's input impedance at the grid frequency (port 2 in z0), as cosim
            target.push(target_ts.as_ref().map(|ts| {
                let g = freqs[fi];
                let k = ts.freqs.partition_point(|&x| x < g).clamp(1, ts.freqs.len() - 1);
                let w = ((g - ts.freqs[k - 1]) / (ts.freqs[k] - ts.freqs[k - 1])).clamp(0.0, 1.0);
                crate::cosim::zin(ts.s[0][0][k - 1] * (1.0 - w) + ts.s[0][0][k] * w, ts.z0)
            }));
        }
        let mut harm_pts = Vec::new();
        for &(f, _) in &st.harmonics {
            harm_pts.push(point(f, at(f), Some(far_q(f)?), &mut pts));
        }
        tm.parts = t.elapsed().as_secs_f64();

        let t = Instant::now();
        rmat.resize(pts.len(), Vec::new());
        for &(pi, q) in &rad_pats {
            let pats: Vec<&Pattern> = names
                .iter()
                .map(|nm| {
                    let k = fnames.iter().position(|x| x == nm).ok_or_else(|| format!("farfield.json has no port {nm}"))?;
                    fpats[k].get(q).and_then(|x| x.as_ref()).ok_or_else(|| format!("farfield.json: port {nm} has no pattern at {:.4e} Hz (run with \"excite\": \"all\")", ffreq[q]))
                })
                .collect::<Result<_, String>>()?;
            rmat[pi] = radiation_matrix(&pats);
        }
        tm.radiation = t.elapsed().as_secs_f64();
        let mut ev = Evaluator {
            names,
            n,
            drive,
            term,
            slot_pos,
            positions,
            choices,
            off,
            band: st.band.clone(),
            band_pts,
            harm_pts,
            pts,
            rmat,
            target,
            harm_freqs: st.harmonics.iter().map(|h| h.0).collect(),
            harm_limit: vec![f64::INFINITY; st.harmonics.len()],
            z0: em.z0,
        };
        // "as built" limits are filled in by the caller once the base candidate is known
        for (k, &(_, lim)) in st.harmonics.iter().enumerate() {
            if let Some(l) = lim {
                ev.harm_limit[k] = l;
            }
        }
        Ok((ev, tm))
    }

    pub fn npos(&self) -> usize {
        self.positions.len()
    }

    /// Network solve at point `pi`: incident waves on every port (1 on the driven port) and the wave
    /// reflected at the driven port. `None` if the termination is singular.
    #[inline]
    fn solve(&self, pi: usize, c: &[u16], a: &mut [Cx; MAXP]) -> Option<Cx> {
        let pt = &self.pts[pi];
        let (n, d, nt) = (self.n, self.drive, self.term.len());
        let s = |i: usize, j: usize| pt.s[i * n + j];
        let mut g = [Cx::default(); MAXP];
        for (slot, sp) in self.slot_pos.iter().enumerate() {
            g[slot] = match sp {
                Some(p) => pt.gam[self.off[*p] + c[*p] as usize],
                None => pt.fixed[slot],
            };
        }
        // (I - S_tt Γ) x = S_te
        let mut m = [Cx::default(); MAXP * MAXP];
        let mut x = [Cx::default(); MAXP];
        for r in 0..nt {
            let tr = self.term[r];
            for q in 0..nt {
                let v = Cx::default() - s(tr, self.term[q]) * g[q];
                m[r * nt + q] = if r == q { v + Cx::new(1.0, 0.0) } else { v };
            }
            x[r] = s(tr, d);
        }
        for col in 0..nt {
            let mut piv = col;
            let mut best = m[col * nt + col].abs();
            for r in col + 1..nt {
                let v = m[r * nt + col].abs();
                if v > best {
                    best = v;
                    piv = r;
                }
            }
            if best < 1e-300 {
                return None;
            }
            if piv != col {
                for q in 0..nt {
                    m.swap(col * nt + q, piv * nt + q);
                }
                x.swap(col, piv);
            }
            let inv = Cx::new(1.0, 0.0) / m[col * nt + col];
            for r in col + 1..nt {
                let f = m[r * nt + col] * inv;
                if f.re == 0.0 && f.im == 0.0 {
                    continue;
                }
                for q in col + 1..nt {
                    m[r * nt + q] = m[r * nt + q] - f * m[col * nt + q];
                }
                x[r] = x[r] - f * x[col];
            }
        }
        for r in (0..nt).rev() {
            let mut v = x[r];
            for q in r + 1..nt {
                v = v - m[r * nt + q] * x[q];
            }
            x[r] = v / m[r * nt + r];
        }
        *a = [Cx::default(); MAXP];
        a[d] = Cx::new(1.0, 0.0);
        let mut be = s(d, d);
        for r in 0..nt {
            let ar = g[r] * x[r];
            a[self.term[r]] = ar;
            be = be + s(d, self.term[r]) * ar;
        }
        Some(be)
    }

    /// Worst band score and harmonic violation (allocation-free).
    pub fn score(&self, c: &[u16]) -> Score {
        let mut a = [Cx::default(); MAXP];
        let mut worst = f64::INFINITY;
        for (k, &(mp, rp)) in self.band_pts.iter().enumerate() {
            let (Some(bm), Some(br)) = (self.solve(mp, c, &mut a), self.solve(rp, c, &mut a)) else { return Score { worst: -300.0, violation: 0.0 } };
            let v = self.band_value(k, bm, br, &a).0;
            worst = worst.min(v);
        }
        let mut violation = 0.0;
        for (h, &hp) in self.harm_pts.iter().enumerate() {
            if self.harm_limit[h].is_finite() {
                let Some(_) = self.solve(hp, c, &mut a) else { return Score { worst: -300.0, violation: 0.0 } };
                let db = db10(quad_power(&a[..self.n], &self.rmat[hp]));
                violation += (db - self.harm_limit[h]).max(0.0);
            }
        }
        Score { worst, violation }
    }

    /// (score dB, |Γ_target|, radiated / accepted dB, P_in, P_rad) at band frequency k, from the reflected
    /// wave at the match point, and the waves `a` and reflected wave at the radiation point.
    fn band_value(&self, k: usize, b_match: Cx, b_rad: Cx, a: &[Cx; MAXP]) -> (f64, f64, f64, f64, f64) {
        let rp = self.band_pts[k].1;
        let z = crate::cosim::zin(b_match, self.z0);
        let g = match self.target[k] {
            Some(t) => ((z - t) / (z + t.conj())).abs(),
            None => b_match.abs(),
        };
        let p_rad = quad_power(&a[..self.n], &self.rmat[rp]);
        let p_in = 1.0 - b_rad.abs().powi(2);
        // as tools/ladder_retune.py: radiated / accepted from the dB values (floored at 1e-30)
        let eta_db = db10(p_rad) - db10(p_in);
        let v = 10.0 * ((1.0 - g * g) * 10f64.powf(eta_db / 10.0)).max(1e-12).log10();
        (v, g, eta_db, p_in, p_rad)
    }

    pub fn detail(&self, c: &[u16]) -> Detail {
        let mut a = [Cx::default(); MAXP];
        let mut band = Vec::new();
        for k in 0..self.band_pts.len() {
            let (mp, rp) = self.band_pts[k];
            let bm = self.solve(mp, c, &mut a).unwrap_or_default();
            let br = self.solve(rp, c, &mut a).unwrap_or_default();
            band.push(self.band_value(k, bm, br, &a));
        }
        let harm = self.harm_pts.iter().map(|&hp| {
            self.solve(hp, c, &mut a);
            db10(quad_power(&a[..self.n], &self.rmat[hp]))
        }).collect();
        Detail { band, harm, score: self.score(c) }
    }

    pub fn labels(&self, c: &[u16]) -> Vec<(String, String)> {
        self.positions.iter().zip(c).enumerate().map(|(p, (n, &k))| (n.clone(), self.choices[p][k as usize].0.clone())).collect()
    }

    pub fn cand_of(&self, labels: &serde_json::Map<String, Value>) -> Result<Cand, String> {
        self.positions
            .iter()
            .enumerate()
            .map(|(p, n)| {
                let l = labels.get(n).and_then(|v| v.as_str()).ok_or_else(|| format!("base has no part for {n}"))?;
                self.choices[p].iter().position(|c| c.0 == l).map(|k| k as u16).ok_or_else(|| format!("base {n}={l:?} is not one of its choices"))
            })
            .collect()
    }
}

fn gamma(p: &Part, f: f64, z0: f64) -> Cx {
    match p.z(f) {
        None => Cx::new(1.0, 0.0),
        Some(z) => (z - Cx::new(z0, 0.0)) / (z + Cx::new(z0, 0.0)),
    }
}

fn db10(x: f64) -> f64 {
    10.0 * x.max(1e-30).log10()
}

/// SplitMix64: a small deterministic generator for the random starts.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed)
    }
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    pub fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// The best few distinct candidates by a ranking key.
#[derive(Clone, Default)]
pub struct Top {
    k: usize,
    pub items: Vec<(f64, Score, Cand)>,
}

impl Top {
    pub fn new(k: usize) -> Top {
        Top { k, items: Vec::with_capacity(k + 1) }
    }
    #[inline]
    pub fn offer(&mut self, key: f64, s: Score, c: &[u16]) {
        if self.items.len() >= self.k && key <= self.items[self.items.len() - 1].0 {
            return;
        }
        if self.items.iter().any(|x| x.2 == c) {
            return;
        }
        let at = self.items.partition_point(|x| x.0 >= key);
        self.items.insert(at, (key, s, c.to_vec()));
        self.items.truncate(self.k);
    }
    pub fn merge(mut self, o: Top) -> Top {
        for (k, s, c) in o.items {
            self.offer(k, s, &c);
        }
        self
    }
}

/// Search statistics and the best candidates (by rank, and by worst band value ignoring constraints).
pub struct Found {
    pub top: Top,
    pub unconstrained: Top,
    pub evals: u64,
}

impl Found {
    fn new(k: usize) -> Found {
        Found { top: Top::new(k), unconstrained: Top::new(3), evals: 0 }
    }
    #[inline]
    fn see(&mut self, s: Score, c: &[u16]) {
        self.evals += 1;
        self.top.offer(s.rank(), s, c);
        self.unconstrained.offer(s.worst, s, c);
    }
    fn merge(self, o: Found) -> Found {
        Found { top: self.top.merge(o.top), unconstrained: self.unconstrained.merge(o.unconstrained), evals: self.evals + o.evals }
    }
}

/// Coordinate descent from `start`: positions in random order, each set to its best choice, until a round
/// improves nothing (by more than 1e-4 dB) or `rounds` rounds have run.
pub fn descend(ev: &Evaluator, start: &[u16], rounds: usize, rng: &mut Rng, found: &mut Found) -> (Cand, Score) {
    let mut cur = start.to_vec();
    let mut best = ev.score(&cur);
    found.see(best, &cur);
    let np = ev.npos();
    let mut order: Vec<usize> = (0..np).collect();
    let mut trial = cur.clone();
    for _ in 0..rounds {
        // Fisher-Yates
        for i in (1..np).rev() {
            order.swap(i, rng.below(i + 1));
        }
        let mut improved = false;
        for &p in &order {
            trial.copy_from_slice(&cur);
            let (mut bk, mut bs) = (cur[p], best);
            for k in 0..ev.choices[p].len() as u16 {
                trial[p] = k;
                let s = ev.score(&trial);
                found.see(s, &trial);
                if s.rank() > bs.rank() {
                    bk = k;
                    bs = s;
                }
            }
            if bs.rank() > best.rank() + 1e-4 {
                cur[p] = bk;
                best = bs;
                improved = true;
            }
        }
        if !improved {
            break;
        }
    }
    (cur, best)
}

/// Every candidate differing from `base` in exactly one position, and (if `two`) in exactly two.
pub fn exhaustive(ev: &Evaluator, base: &[u16], two: bool, k: usize) -> Found {
    let np = ev.npos();
    let mut pairs: Vec<(usize, Option<usize>)> = (0..np).map(|p| (p, None)).collect();
    if two {
        for p in 0..np {
            for q in p + 1..np {
                pairs.push((p, Some(q)));
            }
        }
    }
    pairs
        .par_iter()
        .map(|&(p, q)| {
            let mut f = Found::new(k);
            let mut c = base.to_vec();
            for a in (0..ev.choices[p].len() as u16).filter(|&a| a != base[p]) {
                c[p] = a;
                match q {
                    None => f.see(ev.score(&c), &c),
                    Some(q) => {
                        for b in (0..ev.choices[q].len() as u16).filter(|&b| b != base[q]) {
                            c[q] = b;
                            f.see(ev.score(&c), &c);
                        }
                        c[q] = base[q];
                    }
                }
            }
            f
        })
        .reduce(|| Found::new(k), Found::merge)
}

/// Coordinate descent from `base` and `restarts` random starts, in parallel. Start i uses seed + i.
pub fn multistart(ev: &Evaluator, base: &[u16], restarts: usize, rounds: usize, seed: u64, k: usize) -> (Found, Vec<(Cand, Score)>) {
    let res: Vec<(Found, Cand, Score)> = (0..=restarts)
        .into_par_iter()
        .map(|i| {
            let mut rng = Rng::new(seed.wrapping_mul(1_000_003).wrapping_add(i as u64));
            let start: Cand = if i == 0 { base.to_vec() } else { (0..ev.npos()).map(|p| rng.below(ev.choices[p].len()) as u16).collect() };
            let mut f = Found::new(k);
            let (c, s) = descend(ev, &start, rounds, &mut rng, &mut f);
            (f, c, s)
        })
        .collect();
    let mut found = Found::new(k);
    let mut finals = Vec::new();
    for (f, c, s) in res {
        found = found.merge(f);
        finals.push((c, s));
    }
    (found, finals)
}

fn detail_json(ev: &Evaluator, c: &[u16]) -> Value {
    let d = ev.detail(c);
    json!({
        "parts": ev.labels(c).into_iter().map(|(n, l)| (n, Value::String(l))).collect::<serde_json::Map<_, _>>(),
        "worst_db": d.score.worst, "violation_db": d.score.violation, "feasible": d.score.violation <= 0.0,
        "band": ev.band.iter().zip(&d.band).map(|(f, b)| json!({"f": f, "score_db": b.0, "gamma_vs_target": b.1, "rad_over_accepted_db": b.2, "p_in": b.3, "p_rad": b.4})).collect::<Vec<_>>(),
        "harmonics": ev.harm_freqs.iter().zip(&d.harm).zip(&ev.harm_limit).map(|((f, h), l)| json!({"f": f, "p_rad_db": h, "limit_db": if l.is_finite() { json!(l) } else { Value::Null }})).collect::<Vec<_>>(),
    })
}

/// Score a variant from `cosim` output the way tools/ladder_retune.py did: per band frequency
/// (1 - gamma_vs_target^2) * 10^((p_rad_db - p_in_db) / 10), worst in dB; and the harmonics' p_rad_db.
pub fn score_from_cosim(v: &Value, band: &[f64], harm: &[f64]) -> Option<(f64, Vec<f64>)> {
    let near = |rows: &Vec<Value>, f: f64| rows.iter().min_by(|a, b| (a["f"].as_f64().unwrap_or(0.0) - f).abs().partial_cmp(&(b["f"].as_f64().unwrap_or(0.0) - f).abs()).unwrap()).cloned();
    let rep = v["report"].as_array()?;
    let ff = v["farfield"].as_array()?;
    let mut worst = f64::INFINITY;
    for &f in band {
        let r = near(rep, f)?;
        let fr = ff.iter().find(|x| (x["f"].as_f64().unwrap_or(0.0) - f).abs() <= 1e-6 * f)?;
        let g = r["gamma_vs_target"].as_f64().unwrap_or_else(|| 10f64.powf(r["s11_db"].as_f64().unwrap_or(0.0) / 20.0));
        let eta = 10f64.powf((fr["p_rad_db"].as_f64()? - fr["p_in_db"].as_f64()?) / 10.0);
        worst = worst.min(10.0 * ((1.0 - g * g) * eta).max(1e-12).log10());
    }
    let h = harm.iter().map(|&f| ff.iter().find(|x| (x["f"].as_f64().unwrap_or(0.0) - f).abs() <= 1e-6 * f).and_then(|x| x["p_rad_db"].as_f64())).collect::<Option<Vec<_>>>()?;
    Some((worst, h))
}

/// `fieldsolve tune CONFIG.json [--restarts N] [--rounds N] [--exhaustive 0|1|2] [--top N] [--no-outputs]`
pub fn run(cfg_path: &Path, args: &[String]) -> Result<Value, String> {
    let t_all = Instant::now();
    let cfg: Value = serde_json::from_str(&std::fs::read_to_string(cfg_path).map_err(|e| format!("{}: {e}", cfg_path.display()))?).map_err(|e| e.to_string())?;
    // absolute, so the cosim config written under `out` finds the same files
    let parent = cfg_path.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
    let base_dir_abs = std::fs::canonicalize(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    let base_dir = base_dir_abs.as_path();
    let flag = |k: &str| args.iter().position(|a| a == k).and_then(|i| args.get(i + 1)).and_then(|v| v.parse::<f64>().ok());
    let search = cfg.get("search").cloned().unwrap_or(json!({}));
    let sget = |k: &str, d: f64| flag(&format!("--{k}")).or_else(|| search.get(k).and_then(|v| v.as_f64())).unwrap_or(d);
    let (restarts, rounds, exh, seed) = (sget("restarts", 60.0) as usize, sget("rounds", 8.0) as usize, sget("exhaustive", 2.0) as usize, sget("seed", 7.0) as u64);
    let top_k = flag("--top").or_else(|| cfg.get("top").and_then(|v| v.as_f64())).unwrap_or(8.0) as usize;
    let outputs = !args.iter().any(|a| a == "--no-outputs");
    let em = base_dir.join(cfg["em"].as_str().ok_or("tune: needs \"em\"")?);
    let sets = cfg.get("sets").and_then(|v| v.as_object()).cloned().unwrap_or_default();
    let t = Instant::now();
    let mut positions = Vec::new();
    for (pn, items) in cfg["positions"].as_object().ok_or("tune: needs \"positions\"")? {
        positions.push((pn.clone(), expand_choices(items.as_array().ok_or("a position's choices are a list")?, &sets, base_dir)?));
    }
    let t_choices = t.elapsed().as_secs_f64();
    let harmonics: Vec<(f64, Option<f64>)> = cfg.get("harmonics").and_then(|v| v.as_array()).cloned().unwrap_or_default().iter().map(|h| {
        let f = h.get("f").and_then(|v| v.as_f64()).ok_or("harmonic needs \"f\"")?;
        Ok((f, h.get("max_db").and_then(|v| v.as_f64())))
    }).collect::<Result<_, String>>()?;
    let as_built_limit: Vec<bool> = cfg.get("harmonics").and_then(|v| v.as_array()).cloned().unwrap_or_default().iter().map(|h| h.get("max_db").and_then(|v| v.as_str()) == Some("base")).collect();
    let st = Setup {
        em: em.clone(),
        drive: cfg["drive"].as_str().ok_or("tune: needs \"drive\" (the kept, driven port)")?.to_string(),
        target: cfg.get("target").and_then(|v| v.as_str()).map(|t| base_dir.join(t)),
        band: cfg["band"].as_array().ok_or("tune: needs \"band\"")?.iter().filter_map(|v| v.as_f64()).collect(),
        harmonics,
        positions,
        fixed: cfg.get("fixed").and_then(|v| v.as_object()).map(|o| o.iter().map(|(k, v)| (k.clone(), resolve_spec(base_dir, v.as_str().unwrap_or("")))).collect()).unwrap_or_default(),
        deembed: cfg.get("deembed").and_then(|v| v.as_bool()).unwrap_or(false),
    };
    let (mut ev, mut tm) = Evaluator::load(&st)?;
    tm.parts += t_choices;
    let base = ev.cand_of(cfg["base"].as_object().ok_or("tune: needs \"base\" (label per position)")?)?;
    // "max_db": "base": no worse than as built
    let base_detail = ev.detail(&base);
    for (h, &b) in as_built_limit.iter().enumerate() {
        if b {
            ev.harm_limit[h] = base_detail.harm[h];
        }
    }
    let space: f64 = ev.choices.iter().map(|c| c.len() as f64).product();
    eprintln!(
        "tune: {} positions ({}), {:.3e} combinations; load {:.3} s, parts {:.3} s, radiation matrices {:.4} s",
        ev.npos(),
        ev.positions.iter().zip(&ev.choices).map(|(p, c)| format!("{p} {}", c.len())).collect::<Vec<_>>().join(", "),
        space,
        tm.load,
        tm.parts,
        tm.radiation
    );
    let t = Instant::now();
    let (fd, finals) = multistart(&ev, &base, restarts, rounds, seed, top_k);
    let t_descent = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let fe = if exh >= 1 { exhaustive(&ev, &base, exh >= 2, top_k) } else { Found::new(top_k) };
    let t_exh = t.elapsed().as_secs_f64();
    let (n_desc, n_exh) = (fd.evals, fe.evals);
    let found = fd.merge(fe);
    let evals = found.evals;
    let t_search = t_descent + t_exh;
    eprintln!(
        "tune: {evals} evaluations in {t_search:.3} s ({:.3e} per s on {} threads): descent {} starts {n_desc} in {t_descent:.3} s, exhaustive {n_exh} in {t_exh:.3} s",
        evals as f64 / t_search.max(1e-9),
        rayon::current_num_threads(),
        restarts + 1
    );
    let show = |tag: &str, c: &[u16]| {
        let d = ev.detail(c);
        eprintln!(
            "{tag:<12} worst {:+.2} dB{}  {}\n             {}{}",
            d.score.worst,
            if d.score.violation > 0.0 { format!(" (violates by {:.2} dB)", d.score.violation) } else { String::new() },
            ev.labels(c).iter().map(|(n, l)| format!("{n}={l}")).collect::<Vec<_>>().join(" "),
            ev.band.iter().zip(&d.band).map(|(f, b)| format!("{:.2}: {:+.2} dB (|G| {:.2}, rad/acc {:+.2} dB)", f * 1e-9, b.0, b.1, b.2)).collect::<Vec<_>>().join(" | "),
            if ev.harm_freqs.is_empty() { String::new() } else { format!("\n             radiated harmonics per W incident: {}", ev.harm_freqs.iter().zip(&d.harm).map(|(f, h)| format!("{:.2} GHz {h:+.1} dB", f * 1e-9)).collect::<Vec<_>>().join(", ")) }
        );
    };
    show("as built", &base);
    for (i, (_, _, c)) in found.top.items.iter().enumerate().take(5) {
        show(&format!("rank {}", i + 1), c);
    }
    if let Some((_, s, c)) = found.unconstrained.items.first() {
        if s.violation > 0.0 {
            show("unconstr.", c);
        }
    }
    // detailed outputs for the selected candidates only, through cosim
    let out = base_dir.join(cfg.get("out").and_then(|v| v.as_str()).unwrap_or("tune"));
    std::fs::create_dir_all(&out).map_err(|e| e.to_string())?;
    let t = Instant::now();
    let mut selected: Vec<(String, Cand)> = vec![("as built".into(), base.clone())];
    for (i, (_, _, c)) in found.top.items.iter().enumerate() {
        selected.push((format!("rank {}", i + 1), c.clone()));
    }
    if let Some((_, s, c)) = found.unconstrained.items.first() {
        if s.violation > 0.0 && !selected.iter().any(|x| &x.1 == c) {
            selected.push(("best ignoring harmonic limits".into(), c.clone()));
        }
    }
    let mut check = Value::Null;
    if outputs {
        let fixed: serde_json::Map<String, Value> = st.fixed.iter().map(|(k, v)| (k.clone(), json!(v))).collect();
        let variants: serde_json::Map<String, Value> = selected
            .iter()
            .map(|(name, c)| {
                let mut m = fixed.clone();
                for (p, &k) in c.iter().enumerate() {
                    m.insert(ev.positions[p].clone(), json!(ev.choices[p][k as usize].1));
                }
                (name.clone(), Value::Object(m))
            })
            .collect();
        let em_abs = std::fs::canonicalize(&em).unwrap_or(em.clone());
        let mut ccfg = json!({"em": em_abs, "out": "selected", "keep": [st.drive], "variants": variants,
                              "report": ev.band.iter().chain(&ev.harm_freqs).collect::<Vec<_>>(), "plot_ideal": false, "farfield": true, "deembed": st.deembed});
        if let Some(t) = &st.target {
            ccfg["target"] = json!(std::fs::canonicalize(t).unwrap_or(t.clone()));
        }
        let cpath = out.join("selected_cosim.json");
        std::fs::write(&cpath, serde_json::to_string_pretty(&ccfg).unwrap()).map_err(|e| e.to_string())?;
        crate::cosim::run(&cpath)?;
        // the evaluator against cosim on the same candidates
        let cj: Value = serde_json::from_str(&std::fs::read_to_string(out.join("selected/cosim.json")).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
        let mut worst_diff = 0.0f64;
        for (name, c) in &selected {
            let d = ev.detail(c);
            let (w, h) = score_from_cosim(&cj["variants"][name], &ev.band, &ev.harm_freqs).ok_or_else(|| format!("cosim.json: variant {name} incomplete"))?;
            worst_diff = worst_diff.max((w - d.score.worst).abs());
            for (a, b) in h.iter().zip(&d.harm) {
                worst_diff = worst_diff.max((a - b).abs());
            }
        }
        eprintln!("tune: evaluator vs cosim on {} selected candidates: largest difference {worst_diff:.2e} dB", selected.len());
        check = json!({"candidates": selected.len(), "max_abs_diff_db": worst_diff});
    }
    let t_out = t.elapsed().as_secs_f64();
    // "compare_em": the same candidates re-scored on another extraction of the board (e.g. a screening
    // run against the validated one): score differences and whether the ranking holds
    let mut compare = Value::Null;
    if let Some(other) = cfg.get("compare_em").and_then(|v| v.as_str()) {
        let ncmp = cfg.get("compare_top").and_then(|v| v.as_u64()).unwrap_or(30) as usize;
        let mut st2 = Setup { em: base_dir.join(other), ..st.clone() };
        st2.harmonics = st2.harmonics.iter().enumerate().map(|(h, &(f, _))| (f, ev.harm_limit[h].is_finite().then_some(ev.harm_limit[h]))).collect();
        let (ev2, _) = Evaluator::load(&st2)?;
        // a wider pool than `top`: the best distinct candidates of a fresh search
        let (fd, _) = multistart(&ev, &base, restarts, rounds, seed, ncmp);
        let pool = fd.merge(if exh >= 1 { exhaustive(&ev, &base, exh >= 2, ncmp) } else { Found::new(ncmp) }).top;
        let mut cands: Vec<Cand> = vec![base.clone()];
        cands.extend(pool.items.iter().map(|x| x.2.clone()).filter(|c| *c != base));
        let a: Vec<f64> = cands.iter().map(|c| ev.score(c).worst).collect();
        let b: Vec<f64> = cands.iter().map(|c| ev2.score(c).worst).collect();
        let diffs: Vec<f64> = a.iter().zip(&b).map(|(x, y)| y - x).collect();
        let rank = |v: &[f64]| {
            let mut o: Vec<usize> = (0..v.len()).collect();
            o.sort_by(|&i, &j| v[j].partial_cmp(&v[i]).unwrap());
            let mut r = vec![0.0; v.len()];
            for (k, &i) in o.iter().enumerate() {
                r[i] = k as f64;
            }
            r
        };
        let (ra, rb) = (rank(&a), rank(&b));
        let n = a.len() as f64;
        let spearman = 1.0 - 6.0 * ra.iter().zip(&rb).map(|(x, y)| (x - y).powi(2)).sum::<f64>() / (n * (n * n - 1.0)).max(1.0);
        let best_b = (0..cands.len()).max_by(|&i, &j| b[i].partial_cmp(&b[j]).unwrap()).unwrap_or(0);
        eprintln!(
            "tune: {} candidates re-scored on {other}: as built {:+.3} dB vs {:+.3} dB, difference mean {:+.3} max |{:.3}| dB, Spearman {:.3}; its best is rank {} here",
            cands.len(), a[0], b[0], diffs.iter().sum::<f64>() / n, diffs.iter().fold(0.0f64, |m, d| m.max(d.abs())), spearman, ra[best_b] as usize + 1
        );
        compare = json!({"em": other, "candidates": cands.len(), "as_built": [a[0], b[0]], "this": a, "other": b,
                         "diff_mean_db": diffs.iter().sum::<f64>() / n, "diff_max_abs_db": diffs.iter().fold(0.0f64, |m, d| m.max(d.abs())),
                         "spearman": spearman, "other_best_rank_here": ra[best_b] as usize + 1, "this_best_rank_there": rb[0..].iter().zip(&ra).find(|x| *x.1 == 0.0).map(|x| *x.0 as usize + 1),
                         "parts": cands.iter().map(|c| ev.labels(c).into_iter().map(|(n, l)| (n, Value::String(l))).collect::<serde_json::Map<_, _>>()).collect::<Vec<_>>()});
    }
    let res = json!({
        "em": em.to_string_lossy(), "drive": st.drive, "band": ev.band,
        "positions": ev.positions.iter().zip(&ev.choices).map(|(p, c)| json!({"port": p, "choices": c.iter().map(|x| &x.0).collect::<Vec<_>>()})).collect::<Vec<_>>(),
        "combinations": space,
        "search": {"restarts": restarts, "rounds": rounds, "exhaustive": exh, "seed": seed},
        "evaluations": evals, "evaluations_descent": n_desc, "evaluations_exhaustive": n_exh,
        "threads": rayon::current_num_threads(),
        "timings_s": {"load": tm.load, "parts": tm.parts, "radiation_matrices": tm.radiation, "descent": t_descent, "exhaustive": t_exh,
                      "search": t_search, "outputs": t_out, "total": t_all.elapsed().as_secs_f64()},
        "evaluations_per_s": evals as f64 / t_search.max(1e-9),
        "as_built": detail_json(&ev, &base),
        "top": found.top.items.iter().map(|(_, _, c)| detail_json(&ev, c)).collect::<Vec<_>>(),
        "best_ignoring_constraints": found.unconstrained.items.first().map(|(_, _, c)| detail_json(&ev, c)),
        "descent_finals": finals.iter().map(|(c, s)| json!({"worst_db": s.worst, "violation_db": s.violation, "parts": ev.labels(c).into_iter().map(|(n, l)| (n, Value::String(l))).collect::<serde_json::Map<_, _>>()})).collect::<Vec<_>>(),
        "cosim_check": check,
        "compare": compare,
    });
    std::fs::write(out.join("tune.json"), serde_json::to_string_pretty(&res).unwrap()).map_err(|e| e.to_string())?;
    eprintln!("tune: wrote {} (total {:.2} s)", out.display(), t_all.elapsed().as_secs_f64());
    Ok(res)
}

/// Random complex weights with |re|, |im| <= 1.
fn rand_weights(rng: &mut Rng, n: usize) -> Vec<Cx> {
    let mut u = || (rng.next() >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0;
    (0..n).map(|_| Cx::new(u(), u())).collect()
}

/// The weights with a_i = 1 that radiate least: a_rest = -R_rr^-1 R_ri (the deepest cancellation of port
/// i's pattern by the others).
pub fn cancelling_weights(r: &[Cx], n: usize, i: usize) -> Option<Vec<Cx>> {
    let rest: Vec<usize> = (0..n).filter(|&k| k != i).collect();
    let rr: Vec<Vec<Cx>> = rest.iter().map(|&p| rest.iter().map(|&q| r[p * n + q]).collect()).collect();
    let inv = crate::post::inv(&rr)?;
    let mut a = vec![Cx::default(); n];
    a[i] = Cx::new(1.0, 0.0);
    for (x, &p) in rest.iter().enumerate() {
        a[p] = rest.iter().enumerate().fold(Cx::default(), |acc, (y, &q)| acc - inv[x][y] * r[q * n + i]);
    }
    Some(a)
}

/// `fieldsolve radmatrix farfield.json`: the radiation matrix against `Pattern::combine` + `p_rad` on a
/// saved board far field: single ports, `n` random weight vectors (timed both ways) and, per port, the
/// most-cancelling combination of the others. Errors are relative to the incoherent power sum |a_k|^2 R_kk.
pub fn check_radiation_matrix(path: &Path, n: usize) -> Result<Value, String> {
    let (freqs, names, pats) = crate::nf2ff::read_farfield(path)?;
    let np = names.len();
    let mut rows = Vec::new();
    for (q, &f) in freqs.iter().enumerate() {
        let Some(ps) = pats.iter().map(|p| p.get(q).and_then(|x| x.as_ref())).collect::<Option<Vec<&Pattern>>>() else { continue };
        let t = Instant::now();
        let r = radiation_matrix(&ps);
        let t_setup = t.elapsed().as_secs_f64();
        let combine = |a: &[Cx]| Pattern::combine(&ps.iter().zip(a).map(|(p, &w)| (*p, w)).collect::<Vec<_>>()).p_rad();
        let incoherent = |a: &[Cx]| (0..np).map(|k| a[k].abs().powi(2) * r[k * np + k].re).sum::<f64>();
        let (mut max_rel, mut max_abs, mut max_interf) = (0.0f64, 0.0f64, 0.0f64);
        let mut cmp = |a: &[Cx]| -> (f64, f64) {
            let (p0, p1) = (combine(a), quad_power(a, &r));
            let sc = incoherent(a);
            max_abs = max_abs.max((p0 - p1).abs());
            max_rel = max_rel.max((p0 - p1).abs() / sc);
            max_interf = max_interf.max((sc - p0).abs() / p0.abs().max(1e-300));
            (p0, p1)
        };
        for k in 0..np {
            let mut a = vec![Cx::default(); np];
            a[k] = Cx::new(1.0, 0.0);
            cmp(&a);
        }
        let mut rng = Rng::new(0x5eed + q as u64);
        let ws: Vec<Vec<Cx>> = (0..n).map(|_| rand_weights(&mut rng, np)).collect();
        for a in &ws {
            cmp(a);
        }
        // deepest cancellation of each port by the others
        let mut cancel = Vec::new();
        for k in 0..np {
            if let Some(a) = cancelling_weights(&r, np, k) {
                let (p0, p1) = cmp(&a);
                cancel.push(json!({"port": names[k], "alone": r[k * np + k].re, "cancelled_combine": p0, "cancelled_matrix": p1, "depth_db": 10.0 * (p0.abs().max(1e-300) / r[k * np + k].re).log10(),
                                   "incoherent_sum": incoherent(&a)}));
            }
        }
        let t = Instant::now();
        let s0: f64 = ws.iter().map(|a| combine(std::hint::black_box(a))).sum();
        let t_combine = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let s1: f64 = ws.iter().map(|a| quad_power(std::hint::black_box(a), std::hint::black_box(&r))).sum();
        let t_quad = t.elapsed().as_secs_f64();
        std::hint::black_box((s0, s1));
        rows.push(json!({"f": f, "ports": np, "directions": ps[0].eth.len(), "candidates": n,
                         "max_abs_diff_w": max_abs, "max_diff_over_incoherent": max_rel, "max_interference_over_power": max_interf,
                         "cancellation": cancel, "matrix_setup_s": t_setup, "combine_s": t_combine, "matrix_s": t_quad, "speedup": t_combine / t_quad.max(1e-12)}));
    }
    Ok(json!({"farfield": path.to_string_lossy(), "freqs": rows}))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 7-port slice of the MTH R7A extraction (real board patterns, 30° grid): sparams.s7p, result.json and
    /// farfield.json at 2.40/2.44/2.48/4.88 GHz.
    fn fixture() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src/bin/fieldsolve/testdata/mth7_mini")
    }

    /// A temp dir with the fixture under em/, a synthetic target network (series 3 nH, shunt 1 pF) and a
    /// small "vendor library" of series-measured capacitors lib/CAPx{code}{tol}.s2p.
    fn workdir(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("fs_tune_{tag}_{}", std::process::id()));
        std::fs::create_dir_all(d.join("em")).unwrap();
        std::fs::create_dir_all(d.join("lib")).unwrap();
        for f in ["sparams.s7p", "result.json", "farfield.json"] {
            std::fs::copy(fixture().join(f), d.join("em").join(f)).unwrap();
        }
        let freqs = crate::post::linspace(2.0e9, 5.0e9, 61);
        let steps = vec![(false, Part::parse("L:3n").unwrap()), (true, Part::parse("C:1p").unwrap())];
        let net: Vec<[[Cx; 2]; 2]> = freqs.iter().map(|&f| ladder_at(&steps, f)).collect();
        let s: Vec<Vec<Vec<Cx>>> = (0..2).map(|i| (0..2).map(|j| net.iter().map(|m| m[i][j]).collect()).collect()).collect();
        crate::post::write_touchstone(&d.join("target.s2p"), &freqs, &s, 50.0).unwrap();
        for (code, c) in [("R50", 0.5e-12), ("1R0", 1.0e-12), ("1R5", 1.5e-12)] {
            let z0 = Cx::new(50.0, 0.0);
            let s: Vec<Vec<Vec<Cx>>> = (0..2).map(|i| (0..2).map(|j| freqs.iter().map(|&f| {
                let w = 2.0 * std::f64::consts::PI * f;
                let z = Cx::new(0.2, w * 0.35e-9 - 1.0 / (w * c));
                if i == j { z / (z + z0 * 2.0) } else { z0 * 2.0 / (z + z0 * 2.0) }
            }).collect()).collect()).collect();
            crate::post::write_touchstone(&d.join(format!("lib/CAPx{code}B.s2p")), &freqs, &s, 50.0).unwrap();
        }
        d
    }

    fn ladder_at(steps: &[(bool, Part)], f: f64) -> [[Cx; 2]; 2] {
        crate::cosim::ladder(steps, f, 50.0)
    }

    /// On real board patterns the matrix form equals combine + p_rad for single ports, random weights and the
    /// deepest cancellation of each port by the others, where summing port powers would be badly wrong.
    #[test]
    fn radiation_matrix_matches_combine_on_board_patterns() {
        let v = check_radiation_matrix(&fixture().join("farfield.json"), 300).unwrap();
        let rows = v["freqs"].as_array().unwrap();
        assert_eq!(rows.len(), 4);
        let mut deepest = 0.0f64;
        for r in rows {
            let rel = r["max_diff_over_incoherent"].as_f64().unwrap();
            assert!(rel < 1e-12, "f {}: matrix vs combine {rel:e}", r["f"]);
            for c in r["cancellation"].as_array().unwrap() {
                let (a, b, inc) = (c["cancelled_combine"].as_f64().unwrap(), c["cancelled_matrix"].as_f64().unwrap(), c["incoherent_sum"].as_f64().unwrap());
                assert!((a - b).abs() <= 1e-12 * inc, "{c}");
                deepest = deepest.min(c["depth_db"].as_f64().unwrap());
            }
            assert!(r["max_interference_over_power"].as_f64().unwrap() > 10.0, "interference should matter");
        }
        assert!(deepest < -30.0, "the board's ports should cancel deeply somewhere: {deepest} dB");
    }

    fn setup(d: &Path, positions: Vec<(String, Choices)>, harm: Vec<(f64, Option<f64>)>) -> Setup {
        Setup {
            em: d.join("em/sparams.s7p"),
            drive: "CHIP".into(),
            target: Some(d.join("target.s2p")),
            band: vec![2.40e9, 2.44e9, 2.48e9],
            harmonics: harm,
            positions,
            fixed: vec![],
            deembed: false,
        }
    }

    /// The evaluator's scores, |Γ| and harmonics equal those computed from `cosim`'s output for the same parts.
    #[test]
    fn evaluator_matches_cosim() {
        let d = workdir("cosim");
        let lib = |c: &str| format!("s2p:{}", d.join(format!("lib/CAPx{c}B.s2p")).display());
        let variants: Vec<Vec<(&str, String)>> = vec![
            vec![("C59", "C:1p".into()), ("L4", "L:4.7n,q=18@2.4G".into()), ("C60", "C:1.2p".into()), ("L5", "L:2.2n".into()), ("R70", "short".into()), ("C20", "open".into())],
            vec![("C59", "C:1p".into()), ("L4", "L:4.7n,q=18@2.4G".into()), ("C60", "C:1.2p".into()), ("L5", "short".into()), ("R70", "short".into()), ("C20", lib("R50"))],
            vec![("C59", "open".into()), ("L4", "short".into()), ("C60", lib("1R5")), ("L5", "short".into()), ("R70", "L:1.6n,srf=9G".into()), ("C20", "C:0.7p,esr=0.3,esl=0.3n".into())],
            vec![("C59", "open".into()), ("L4", "short".into()), ("C60", "open".into()), ("L5", "L:2.2n".into()), ("R70", "open".into()), ("C20", "short".into())],
        ];
        let mut positions: Vec<(String, Choices)> = Vec::new();
        for v in &variants {
            for (p, s) in v {
                let e = match positions.iter_mut().find(|x| x.0 == *p) {
                    Some(e) => e,
                    None => {
                        positions.push((p.to_string(), vec![]));
                        positions.last_mut().unwrap()
                    }
                };
                if !e.1.iter().any(|c| &c.1 == s) {
                    e.1.push((s.clone(), s.clone()));
                }
            }
        }
        let (ev, _) = Evaluator::load(&setup(&d, positions, vec![(4.88e9, None)])).unwrap();
        let vjson: serde_json::Map<String, Value> = variants.iter().enumerate().map(|(i, v)| (format!("v{i}"), Value::Object(v.iter().map(|(p, s)| (p.to_string(), json!(s))).collect()))).collect();
        let cfg = json!({"em": "em/sparams.s7p", "out": "cosim", "keep": ["CHIP"], "target": "target.s2p", "variants": vjson,
                         "report": [2.40e9, 2.44e9, 2.48e9, 4.88e9], "plot_ideal": false, "farfield": true});
        std::fs::write(d.join("cosim.json"), cfg.to_string()).unwrap();
        crate::cosim::run(&d.join("cosim.json")).unwrap();
        let cj: Value = serde_json::from_str(&std::fs::read_to_string(d.join("cosim/cosim.json")).unwrap()).unwrap();
        let mut spread = Vec::new();
        for (i, v) in variants.iter().enumerate() {
            let labels: serde_json::Map<String, Value> = v.iter().map(|(p, s)| (p.to_string(), json!(s))).collect();
            let c = ev.cand_of(&labels).unwrap();
            let det = ev.detail(&c);
            let cv = &cj["variants"][format!("v{i}")];
            let (w, h) = score_from_cosim(cv, &ev.band, &[4.88e9]).unwrap();
            assert!((w - det.score.worst).abs() < 1e-9, "v{i}: worst {w} (cosim) vs {} (evaluator)", det.score.worst);
            assert!((h[0] - det.harm[0]).abs() < 1e-9, "v{i}: H2 {} vs {}", h[0], det.harm[0]);
            for (k, row) in cv["report"].as_array().unwrap().iter().take(3).enumerate() {
                assert!((row["gamma_vs_target"].as_f64().unwrap() - det.band[k].1).abs() < 1e-12, "v{i} |G| at band {k}");
                let fr = &cv["farfield"][k];
                assert!((fr["p_rad_db"].as_f64().unwrap() - 10.0 * det.band[k].4.log10()).abs() < 1e-9, "v{i} P_rad at band {k}");
                assert!((fr["p_in_db"].as_f64().unwrap() - 10.0 * det.band[k].3.log10()).abs() < 1e-9, "v{i} P_in at band {k}");
            }
            spread.push(det.score.worst);
        }
        assert!(spread.iter().cloned().fold(f64::NEG_INFINITY, f64::max) - spread.iter().cloned().fold(f64::INFINITY, f64::min) > 1.0, "the variants should differ: {spread:?}");
    }

    /// The whole command on the fixture: library sets from a template, a harmonic limit "as built", search,
    /// detailed output for the selected candidates only, and the evaluator-vs-cosim self-check.
    #[test]
    fn tune_command_end_to_end() {
        let d = workdir("run");
        let cfg = json!({
            "em": "em/sparams.s7p", "out": "tuned", "drive": "CHIP", "target": "target.s2p",
            "band": [2.40e9, 2.44e9, 2.48e9], "harmonics": [{"f": 4.88e9, "max_db": "base"}],
            "sets": {"caps": {"template": "lib/CAPx{code}{tol}.s2p", "tols": "CB", "unit": "p", "values": [0.5, 1.0, 1.2, 1.5]}},
            "positions": {"C59": ["open", "@caps"], "L4": ["short", "L:4.7n", "L:3.3n"], "C60": ["open", "@caps", "C:1.2p"],
                          "L5": ["short", "L:2.2n"], "R70": ["short", {"set": "caps", "min": "1p", "suffix": " ser"}], "C20": ["open", {"set": "caps", "max": "1p"}]},
            "base": {"C59": "1p", "L4": "L:4.7n", "C60": "C:1.2p", "L5": "L:2.2n", "R70": "short", "C20": "open"},
            "search": {"restarts": 6, "rounds": 6, "exhaustive": 2, "seed": 3}, "top": 4,
            "compare_em": "em/sparams.s7p", "compare_top": 12
        });
        std::fs::write(d.join("tune.json"), cfg.to_string()).unwrap();
        let r = run(&d.join("tune.json"), &[]).unwrap();
        // 1.2 pF has no library file: the set holds 0.5, 1 and 1.5 only
        let ch = |p: &str| r["positions"].as_array().unwrap().iter().find(|x| x["port"] == p).unwrap()["choices"].clone();
        assert_eq!(ch("C59"), json!(["open", "0.5p", "1p", "1.5p"]));
        assert_eq!(ch("R70"), json!(["short", "1p ser", "1.5p ser"]));
        assert_eq!(ch("C20"), json!(["open", "0.5p", "1p"]));
        assert!(r["cosim_check"]["max_abs_diff_db"].as_f64().unwrap() < 1e-9, "{}", r["cosim_check"]);
        let base = r["as_built"]["worst_db"].as_f64().unwrap();
        let best = &r["top"][0];
        assert!(best["feasible"].as_bool().unwrap());
        assert!(best["worst_db"].as_f64().unwrap() >= base, "search lost to as built");
        assert!(best["harmonics"][0]["p_rad_db"].as_f64().unwrap() <= r["as_built"]["harmonics"][0]["p_rad_db"].as_f64().unwrap() + 1e-12);
        assert!(r["evaluations"].as_u64().unwrap() > 100);
        // re-scoring on the same extraction changes nothing
        let cmp = &r["compare"];
        assert_eq!(cmp["candidates"], 12);
        assert_eq!(cmp["diff_max_abs_db"], 0.0);
        assert_eq!(cmp["spearman"], 1.0);
        // detailed files exist only for the selected candidates
        let n_ts = std::fs::read_dir(d.join("tuned/selected")).unwrap().filter(|e| e.as_ref().unwrap().path().extension().is_some_and(|x| x == "s1p")).count();
        assert!(n_ts >= 2 && n_ts <= 1 + 4 + 1, "{n_ts} Touchstone files");
    }

    /// On a space small enough to enumerate, descent with restarts plus the exhaustive pass finds the optimum,
    /// and an infeasible candidate never outranks a feasible one.
    #[test]
    fn search_finds_brute_force_optimum() {
        let d = workdir("brute");
        let caps: Choices = ["open", "C:0.3p", "C:0.5p", "C:0.8p", "C:1p", "C:1.5p", "C:2.2p"].iter().map(|s| (s.to_string(), s.to_string())).collect();
        let inds: Choices = ["short", "L:1n", "L:2.2n", "L:3.3n", "L:4.7n"].iter().map(|s| (s.to_string(), s.to_string())).collect();
        let positions = vec![("C59".to_string(), caps.clone()), ("L4".to_string(), inds.clone()), ("C60".to_string(), caps.clone()), ("L5".to_string(), inds), ("C20".to_string(), caps)];
        let mut st = setup(&d, positions, vec![(4.88e9, Some(-26.0))]);
        st.fixed = vec![("R70".into(), "short".into())];
        let (ev, _) = Evaluator::load(&st).unwrap();
        let dims: Vec<usize> = ev.choices.iter().map(|c| c.len()).collect();
        let total: usize = dims.iter().product();
        let (mut best, mut bc) = (f64::NEG_INFINITY, vec![]);
        let mut c = vec![0u16; dims.len()];
        let mut n_feasible = 0;
        for k in 0..total {
            let mut r = k;
            for (p, &dp) in dims.iter().enumerate() {
                c[p] = (r % dp) as u16;
                r /= dp;
            }
            let s = ev.score(&c);
            n_feasible += (s.violation <= 0.0) as usize;
            if s.rank() > best {
                best = s.rank();
                bc = c.clone();
            }
        }
        assert!(n_feasible > 0 && n_feasible < total, "the harmonic limit should bind: {n_feasible}/{total}");
        let base = vec![4u16, 4, 4, 2, 0];
        let (fd, _) = multistart(&ev, &base, 12, 8, 7, 5);
        let found = fd.merge(exhaustive(&ev, &base, true, 5));
        let (r0, s0, c0) = &found.top.items[0];
        assert!((r0 - best).abs() < 1e-12, "found {r0} at {c0:?}, optimum {best} at {bc:?}");
        assert!(s0.violation <= 0.0);
        assert!(found.evals < total as u64, "the search should not need the whole space");
    }

    #[test]
    fn rkm_codes() {
        assert_eq!(rkm(0.2, 'R'), "R20");
        assert_eq!(rkm(0.75, 'R'), "R75");
        assert_eq!(rkm(1.0, 'R'), "1R0");
        assert_eq!(rkm(8.2, 'R'), "8R2");
        assert_eq!(rkm(4.7, 'N'), "4N7");
        assert_eq!(rkm(10.0, 'N'), "10N");
        let (a, b) = (Score { worst: -10.0, violation: 0.0 }, Score { worst: -1.0, violation: 0.1 });
        assert!(a.rank() > b.rank());
    }
}
