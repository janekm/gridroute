//! EM–circuit co-simulation (roadmap F1, first part): connect part models to the ports of a board's EM
//! N-port and reduce it to the external ports.
//!
//! The EM run has a port across the pads of every two-terminal part (the part removed, the copper as drawn)
//! plus ports at the chip and connector reference planes. A two-terminal part across its own port is a
//! one-port load, so connecting it terminates that port in the part's impedance Z(f). With reflection
//! Γ = (Z − Z0) / (Z + Z0) on the terminated ports t and the kept ports e,
//!
//!   S'_ee = S_ee + S_et Γ (I − S_tt Γ)^-1 S_te.
//!
//! Parts are ideal or simple parasitic models, or vendor Touchstone files: `s2p:` measured in series through
//! a fixture, Z = 2 Z0 S11 / S21 (Murata); `s2p-shunt:` measured as a shunt to ground on a through line,
//! Z = −Z0 S21 / (2 S11) (Yageo); `s1p:` a one-port, Z = Z0 (1 + S11) / (1 − S11). Both two-port forms are
//! ratios, so a symmetric matched fixture line left in the data cancels.
//! An ideal ladder (no layout) is cascaded with ABCD matrices for comparison, so the layout's own effect
//! can be seen. The kept ports can be renormalised to complex reference impedances (power waves).

use crate::post::{inv, matmul, Cx};
use std::f64::consts::PI;
use std::path::Path;

pub struct Touchstone {
    pub freqs: Vec<f64>,
    pub z0: f64,
    /// `s[i][j][f]`
    pub s: Vec<Vec<Vec<Cx>>>,
}

/// Touchstone v1 reader: `# <unit> S <RI|MA|DB> R <z0>`; N from the file extension (`.sNp`).
pub fn read_touchstone(path: &Path) -> Result<Touchstone, String> {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
    let n: usize = ext.trim_start_matches('s').trim_end_matches('p').parse().map_err(|_| format!("{}: not a .sNp file", path.display()))?;
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let (mut fmul, mut fmt, mut z0) = (1e9, "MA".to_string(), 50.0);
    let mut nums = Vec::new();
    for line in text.lines() {
        let l = line.split('!').next().unwrap_or("").trim();
        if l.is_empty() {
            continue;
        }
        if let Some(opt) = l.strip_prefix('#') {
            let t: Vec<String> = opt.split_whitespace().map(|s| s.to_ascii_uppercase()).collect();
            let mut k = 0;
            while k < t.len() {
                match t[k].as_str() {
                    "HZ" => fmul = 1.0,
                    "KHZ" => fmul = 1e3,
                    "MHZ" => fmul = 1e6,
                    "GHZ" => fmul = 1e9,
                    "RI" | "MA" | "DB" => fmt = t[k].clone(),
                    "R" => {
                        z0 = t.get(k + 1).and_then(|v| v.parse().ok()).unwrap_or(50.0);
                        k += 1;
                    }
                    _ => {}
                }
                k += 1;
            }
            continue;
        }
        for w in l.split_whitespace() {
            nums.push(w.parse::<f64>().map_err(|_| format!("{}: bad number {w:?}", path.display()))?);
        }
    }
    let rec = 1 + 2 * n * n;
    if nums.len() % rec != 0 {
        return Err(format!("{}: {} numbers is not a multiple of {rec}", path.display(), nums.len()));
    }
    let nf = nums.len() / rec;
    let mut s = vec![vec![vec![Cx::default(); nf]; n]; n];
    let mut freqs = Vec::with_capacity(nf);
    for f in 0..nf {
        let r = &nums[f * rec..(f + 1) * rec];
        freqs.push(r[0] * fmul);
        for q in 0..n * n {
            let (i, j) = if n == 2 { [(0, 0), (1, 0), (0, 1), (1, 1)][q] } else { (q / n, q % n) };
            let (a, b) = (r[1 + 2 * q], r[2 + 2 * q]);
            s[i][j][f] = match fmt.as_str() {
                "RI" => Cx::new(a, b),
                "MA" => Cx::new(a * (b * PI / 180.0).cos(), a * (b * PI / 180.0).sin()),
                _ => {
                    let m = 10f64.powf(a / 20.0);
                    Cx::new(m * (b * PI / 180.0).cos(), m * (b * PI / 180.0).sin())
                }
            };
        }
    }
    Ok(Touchstone { freqs, z0, s })
}

/// Number with an SI suffix: `1p`, `4.7n`, `2.4G`, `0.3`, `1e-12`.
pub fn si(s: &str) -> Result<f64, String> {
    let s = s.trim();
    let (body, mul) = match s.chars().last() {
        Some('f') => (&s[..s.len() - 1], 1e-15),
        Some('p') => (&s[..s.len() - 1], 1e-12),
        Some('n') => (&s[..s.len() - 1], 1e-9),
        Some('u') => (&s[..s.len() - 1], 1e-6),
        Some('m') => (&s[..s.len() - 1], 1e-3),
        Some('k') => (&s[..s.len() - 1], 1e3),
        Some('M') => (&s[..s.len() - 1], 1e6),
        Some('G') => (&s[..s.len() - 1], 1e9),
        _ => (s, 1.0),
    };
    body.parse::<f64>().map(|v| v * mul).map_err(|_| format!("bad value {s:?}"))
}

/// A two-terminal part's impedance model.
#[derive(Clone)]
pub enum Part {
    Open,
    Short,
    R(f64),
    /// C with series ESR and ESL
    C { c: f64, esr: f64, esl: f64 },
    /// L with a parallel C (from the SRF) and a series R that grows as sqrt(f) (from Q at f0)
    L { l: f64, cp: f64, q: f64, fq: f64 },
    /// impedance table from a vendor file
    Table { freqs: Vec<f64>, z: Vec<Cx> },
}

impl Part {
    /// `C:1p[,esr=0.3,esl=0.4n]`, `L:4.7n[,srf=6G,q=18@2.4G]`, `R:50`, `open`, `short`, `s2p:FILE` (measured in
    /// series through a fixture), `s2p-shunt:FILE` (measured as a shunt to ground on a through line), `s1p:FILE`.
    pub fn parse(spec: &str) -> Result<Part, String> {
        let (kind, rest) = spec.split_once(':').unwrap_or((spec, ""));
        let mut it = rest.split(',');
        let val = it.next().unwrap_or("");
        let kv: Vec<(String, String)> = it.filter_map(|p| p.split_once('=').map(|(a, b)| (a.trim().to_string(), b.trim().to_string()))).collect();
        let get = |k: &str| kv.iter().find(|x| x.0 == k).map(|x| x.1.clone());
        match kind.to_ascii_lowercase().as_str() {
            "open" => Ok(Part::Open),
            "short" => Ok(Part::Short),
            "r" => Ok(Part::R(si(val)?)),
            "c" => Ok(Part::C { c: si(val)?, esr: get("esr").map(|v| si(&v)).transpose()?.unwrap_or(0.0), esl: get("esl").map(|v| si(&v)).transpose()?.unwrap_or(0.0) }),
            "l" => {
                let l = si(val)?;
                let cp = match get("srf") {
                    Some(f) => 1.0 / ((2.0 * PI * si(&f)?).powi(2) * l),
                    None => 0.0,
                };
                let (q, fq) = match get("q") {
                    Some(q) => {
                        let (a, b) = q.split_once('@').unwrap_or((&q, "1G"));
                        (si(a)?, si(b)?)
                    }
                    None => (f64::INFINITY, 1e9),
                };
                Ok(Part::L { l, cp, q, fq })
            }
            "s2p" | "s2p-shunt" | "s1p" => {
                let ts = read_touchstone(Path::new(val))?;
                let z0 = Cx::new(ts.z0, 0.0);
                let one = Cx::new(1.0, 0.0);
                let k = kind.to_ascii_lowercase();
                let z = (0..ts.freqs.len())
                    .map(|f| {
                        let (s11, s21) = (ts.s[0][0][f], if ts.s.len() > 1 { ts.s[1][0][f] } else { Cx::default() });
                        match k.as_str() {
                            // series through: S11 + S21 = 1 for the bare part; the ratio is unaffected by a
                            // symmetric matched fixture line, whose phase multiplies both
                            "s2p" => z0 * 2.0 * s11 / s21,
                            // shunt to ground on a through line: S21 - S11 = 1 for the bare part
                            "s2p-shunt" => Cx::default() - z0 * s21 / (s11 * 2.0),
                            _ => z0 * (one + s11) / (one - s11),
                        }
                    })
                    .collect();
                Ok(Part::Table { freqs: ts.freqs, z })
            }
            _ => Err(format!("unknown part {spec:?}")),
        }
    }

    /// Impedance at f; `None` for an open.
    pub fn z(&self, f: f64) -> Option<Cx> {
        let w = 2.0 * PI * f.max(1.0);
        match self {
            Part::Open => None,
            Part::Short => Some(Cx::default()),
            Part::R(r) => Some(Cx::new(*r, 0.0)),
            Part::C { c, esr, esl } => Some(Cx::new(*esr, w * esl - 1.0 / (w * c))),
            Part::L { l, cp, q, fq } => {
                let r = if q.is_finite() { 2.0 * PI * fq * l / q * (f / fq).max(0.0).sqrt() } else { 0.0 };
                let zl = Cx::new(r, w * l);
                if *cp > 0.0 {
                    let yc = Cx::new(0.0, w * cp);
                    Some(Cx::new(1.0, 0.0) / (Cx::new(1.0, 0.0) / zl + yc))
                } else {
                    Some(zl)
                }
            }
            Part::Table { freqs, z } => {
                let k = freqs.partition_point(|&x| x < f).clamp(1, freqs.len() - 1);
                let t = ((f - freqs[k - 1]) / (freqs[k] - freqs[k - 1])).clamp(0.0, 1.0);
                Some(z[k - 1] * (1.0 - t) + z[k] * t)
            }
        }
    }

    fn gamma(&self, f: f64, z0: f64) -> Cx {
        match self.z(f) {
            None => Cx::new(1.0, 0.0),
            Some(z) => (z - Cx::new(z0, 0.0)) / (z + Cx::new(z0, 0.0)),
        }
    }
}

/// Terminate ports `t` (with parts) of `s[i][j]` at one frequency; returns S of the kept ports `e`.
pub fn terminate(s: &[Vec<Cx>], e: &[usize], t: &[(usize, Cx)]) -> Vec<Vec<Cx>> {
    let nt = t.len();
    if nt == 0 {
        return e.iter().map(|&i| e.iter().map(|&j| s[i][j]).collect()).collect();
    }
    // M = I - S_tt Γ
    let m: Vec<Vec<Cx>> = (0..nt).map(|a| (0..nt).map(|b| (if a == b { Cx::new(1.0, 0.0) } else { Cx::default() }) - s[t[a].0][t[b].0] * t[b].1).collect()).collect();
    let mi = inv(&m).expect("singular termination");
    let ste: Vec<Vec<Cx>> = t.iter().map(|&(a, _)| e.iter().map(|&j| s[a][j]).collect()).collect();
    let x = matmul(&mi, &ste); // (I - S_tt Γ)^-1 S_te
    e.iter()
        .map(|&i| {
            e.iter()
                .enumerate()
                .map(|(jj, &j)| {
                    let mut v = s[i][j];
                    for (a, &(ta, g)) in t.iter().enumerate() {
                        v = v + s[i][ta] * g * x[a][jj];
                    }
                    v
                })
                .collect()
        })
        .collect()
}

/// Renormalise S (reference z0 on all ports) to power waves with reference impedances `zr` (complex).
pub fn renormalise(s: &[Vec<Cx>], z0: f64, zr: &[Cx]) -> Vec<Vec<Cx>> {
    let n = s.len();
    let one = Cx::new(1.0, 0.0);
    let eye = |i: usize, j: usize| if i == j { one } else { Cx::default() };
    let ims: Vec<Vec<Cx>> = (0..n).map(|i| (0..n).map(|j| eye(i, j) - s[i][j]).collect()).collect();
    let ips: Vec<Vec<Cx>> = (0..n).map(|i| (0..n).map(|j| eye(i, j) + s[i][j]).collect()).collect();
    let z = matmul(&ips, &inv(&ims).expect("S has an eigenvalue 1"));
    let z: Vec<Vec<Cx>> = z.into_iter().map(|r| r.into_iter().map(|v| v * z0).collect()).collect();
    let a: Vec<Vec<Cx>> = (0..n).map(|i| (0..n).map(|j| z[i][j] - if i == j { zr[i].conj() } else { Cx::default() }).collect()).collect();
    let b: Vec<Vec<Cx>> = (0..n).map(|i| (0..n).map(|j| z[i][j] + if i == j { zr[i] } else { Cx::default() }).collect()).collect();
    let sp = matmul(&a, &inv(&b).expect("Z + Zr singular"));
    let fq: Vec<f64> = zr.iter().map(|z| 1.0 / (2.0 * z.re.max(1e-12).sqrt())).collect();
    (0..n).map(|i| (0..n).map(|j| sp[i][j] * (fq[i] / fq[j])).collect()).collect()
}

/// Ideal ladder between two 50 Ω-referenced ends: ("shunt" | "series", part) in order from port 1.
pub fn ladder(steps: &[(bool, Part)], f: f64, z0: f64) -> [[Cx; 2]; 2] {
    let one = Cx::new(1.0, 0.0);
    let mut m = [[one, Cx::default()], [Cx::default(), one]];
    for (shunt, p) in steps {
        let zero = Cx::default();
        let abcd = match (shunt, p.z(f)) {
            (true, None) => [[one, zero], [zero, one]],
            (true, Some(z)) if z.abs() < 1e-15 => [[one, zero], [Cx::new(1e15, 0.0), one]],
            (true, Some(z)) => [[one, zero], [one / z, one]],
            (false, None) => [[one, Cx::new(1e12, 0.0)], [zero, one]],
            (false, Some(z)) => [[one, z], [zero, one]],
        };
        let r = [
            [m[0][0] * abcd[0][0] + m[0][1] * abcd[1][0], m[0][0] * abcd[0][1] + m[0][1] * abcd[1][1]],
            [m[1][0] * abcd[0][0] + m[1][1] * abcd[1][0], m[1][0] * abcd[0][1] + m[1][1] * abcd[1][1]],
        ];
        m = r;
    }
    let (a, b, c, d) = (m[0][0], m[0][1], m[1][0], m[1][1]);
    let z = Cx::new(z0, 0.0);
    let den = a + b / z + c * z + d;
    [[(a + b / z - c * z - d) / den, Cx::new(2.0, 0.0) * (a * d - b * c) / den], [Cx::new(2.0, 0.0) / den, (Cx::default() - a + b / z - c * z + d) / den]]
}

pub fn db(x: Cx) -> f64 {
    20.0 * x.abs().max(1e-15).log10()
}

/// Input impedance at port 1 of a 2-port with port 2 terminated in z0 (both z0-referenced).
pub fn zin(s11: Cx, z0: f64) -> Cx {
    let one = Cx::new(1.0, 0.0);
    Cx::new(z0, 0.0) * (one + s11) / (one - s11)
}

/// |S11| and |S21| in dB against frequency, for several named curves.
pub fn plot_svg(path: &Path, freqs: &[f64], curves: &[(&str, &str, Vec<f64>, Vec<f64>)], marks: &[f64]) -> std::io::Result<()> {
    use std::fmt::Write as _;
    let (w, h, ml, mr, mt, mb) = (760.0, 440.0, 60.0, 20.0, 30.0, 45.0);
    let (f0, f1) = (freqs[0], *freqs.last().unwrap());
    let (y0, y1) = (-50.0, 0.0);
    let px = |f: f64| ml + (f - f0) / (f1 - f0) * (w - ml - mr);
    let py = |d: f64| mt + (y1 - d.clamp(y0, y1)) / (y1 - y0) * (h - mt - mb);
    let mut o = String::new();
    let _ = write!(o, r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" font-family="sans-serif" font-size="12"><rect width="100%" height="100%" fill="white"/>"#);
    for d in (0..=10).map(|k| -5.0 * k as f64) {
        let _ = write!(o, r##"<line x1="{ml}" x2="{}" y1="{y}" y2="{y}" stroke="#ddd"/><text x="{}" y="{}" text-anchor="end">{d}</text>"##, w - mr, ml - 6.0, py(d) + 4.0, y = py(d));
    }
    let mut g = 0.0;
    while g <= f1 / 1e9 + 1e-9 {
        if g * 1e9 >= f0 {
            let x = px(g * 1e9);
            let _ = write!(o, r##"<line x1="{x}" x2="{x}" y1="{mt}" y2="{}" stroke="#eee"/><text x="{x}" y="{}" text-anchor="middle">{g}</text>"##, h - mb, h - mb + 16.0);
        }
        g += 1.0;
    }
    for &m in marks {
        let x = px(m);
        let _ = write!(o, r##"<line x1="{x}" x2="{x}" y1="{mt}" y2="{}" stroke="#bbb" stroke-dasharray="3,3"/>"##, h - mb);
    }
    let _ = write!(o, r#"<text x="{}" y="{}" text-anchor="middle">frequency [GHz]</text><text x="14" y="{}" transform="rotate(-90 14 {})" text-anchor="middle">dB</text>"#, (ml + w - mr) / 2.0, h - 8.0, h / 2.0, h / 2.0);
    for (name, col, s11, s21) in curves {
        for (dash, ys) in [("", s21), ("5,3", s11)] {
            let pts: Vec<String> = freqs.iter().zip(ys.iter()).map(|(&f, &d)| format!("{:.1},{:.1}", px(f), py(d))).collect();
            let _ = write!(o, r#"<polyline fill="none" stroke="{col}" stroke-width="1.6" stroke-dasharray="{dash}" points="{}"/>"#, pts.join(" "));
        }
    }
    // legend on a white box, sized to the longest name (≈ 6.6 px per character at 12 px)
    let longest = curves.iter().map(|c| c.0.chars().count() + 4).max().unwrap_or(0) as f64;
    let lw = longest * 6.6 + 42.0;
    let lx = (w - mr - lw - 6.0).max(ml + 6.0);
    let _ = write!(o, r#"<rect x="{lx}" y="{}" width="{lw}" height="{}" fill="white" fill-opacity="0.85"/>"#, mt + 2.0, 16.0 * 2.0 * curves.len() as f64 + 6.0);
    for (k, (name, col, _, _)) in curves.iter().enumerate() {
        for (j, (dash, lab)) in [("", "S21"), ("5,3", "S11")].into_iter().enumerate() {
            let ly = mt + 16.0 + 16.0 * (2 * k + j) as f64;
            let _ = write!(o, r#"<line x1="{}" x2="{}" y1="{}" y2="{}" stroke="{col}" stroke-width="2" stroke-dasharray="{dash}"/><text x="{}" y="{ly}">{name} {lab}</text>"#, lx + 6.0, lx + 31.0, ly - 4.0, ly - 4.0, lx + 38.0);
        }
    }
    o.push_str("</svg>");
    std::fs::write(path, o)
}

/// `fieldsolve cosim CONFIG.json`:
/// ```json
/// { "em": "run/sparams.s6p", "out": "cosim", "keep": ["CHIP", "J1"],
///   "variants": { "ideal parts": {"C21": "C:1p", "L4": "L:4.7n"},
///                 "typical 0402": {"C21": "C:1p,esr=0.3,esl=0.35n", "L4": "L:4.7n,srf=6G,q=18@2.4G"} },
///   "ladder": [["shunt", "C21"], ["series", "L4"]],
///   "zref": {"CHIP": [50, 0]},
///   "report": [2.44e9, 4.88e9] }
/// ```
/// Port names come from `result.json` next to the EM file. Every EM port that is not kept must have a part
/// in each variant. For each variant this writes `<variant>.s2p` (layout + parts) and the same parts as an
/// ideal ladder for comparison, a plot and `cosim.json`.
pub fn run(cfg_path: &Path) -> Result<(), String> {
    use serde_json::{json, Value};
    let cfg: Value = serde_json::from_str(&std::fs::read_to_string(cfg_path).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
    let base = cfg_path.parent().unwrap_or(Path::new("."));
    let em_path = base.join(cfg["em"].as_str().ok_or("cosim: needs \"em\"")?);
    let em = read_touchstone(&em_path)?;
    let res: Value = serde_json::from_str(&std::fs::read_to_string(em_path.with_file_name("result.json")).map_err(|e| format!("result.json: {e}"))?).map_err(|e| e.to_string())?;
    let names: Vec<String> = res["port_names"].as_array().ok_or("result.json has no port_names")?.iter().map(|v| v.as_str().unwrap_or("").to_string()).collect();
    if names.len() != em.s.len() {
        return Err(format!("{} port names for a {}-port file", names.len(), em.s.len()));
    }
    let out = base.join(cfg["out"].as_str().unwrap_or("cosim"));
    std::fs::create_dir_all(&out).map_err(|e| e.to_string())?;
    let idx = |n: &str| names.iter().position(|x| x == n).ok_or_else(|| format!("no EM port {n:?} (have {names:?})"));
    let keep: Vec<usize> = cfg["keep"].as_array().ok_or("cosim: needs \"keep\"")?.iter().map(|v| idx(v.as_str().unwrap_or(""))).collect::<Result<_, _>>()?;
    let zr: Vec<Cx> = keep
        .iter()
        .map(|&k| cfg["zref"].get(&names[k]).and_then(|v| v.as_array()).map(|a| Cx::new(a[0].as_f64().unwrap_or(em.z0), a.get(1).and_then(|x| x.as_f64()).unwrap_or(0.0))).unwrap_or(Cx::new(em.z0, 0.0)))
        .collect();
    let renorm = zr.iter().any(|z| (z.re - em.z0).abs() > 1e-9 || z.im.abs() > 1e-9);
    let report: Vec<f64> = cfg["report"].as_array().map(|a| a.iter().filter_map(|v| v.as_f64()).collect()).unwrap_or_else(|| vec![2.44e9]);
    let ladder_spec: Vec<(bool, String)> = cfg["ladder"].as_array().map(|a| a.iter().map(|s| (s[0].as_str() == Some("shunt"), s[1].as_str().unwrap_or("").to_string())).collect()).unwrap_or_default();
    let freqs = &em.freqs;
    let plot_ideal = cfg["plot_ideal"].as_bool().unwrap_or(true);
    // target: a measured matching network (port 1 = chip side) loaded with Z0, i.e. the load the chip is
    // designed to see; each case is scored by the power-wave reflection against it
    let target: Option<Vec<Cx>> = match cfg["target"].as_str() {
        Some(t) => {
            let ts = read_touchstone(&base.join(t))?;
            let z = freqs
                .iter()
                .map(|&f| {
                    let k = ts.freqs.partition_point(|&x| x < f).clamp(1, ts.freqs.len() - 1);
                    let (f0, f1) = (ts.freqs[k - 1], ts.freqs[k]);
                    let w = ((f - f0) / (f1 - f0)).clamp(0.0, 1.0);
                    let s11 = ts.s[0][0][k - 1] * (1.0 - w) + ts.s[0][0][k] * w;
                    zin(s11, ts.z0)
                })
                .collect();
            Some(z)
        }
        None => None,
    };
    let gam_t = |z: Cx, fi: usize| target.as_ref().map(|t| ((z - t[fi]) / (z + t[fi].conj())).abs());
    let colours = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"];
    let mut curves = Vec::new();
    let mut out_json = serde_json::Map::new();
    let vars = cfg["variants"].as_object().ok_or("cosim: needs \"variants\"")?;
    let at = |f: f64| freqs.iter().enumerate().min_by(|a, b| (a.1 - f).abs().partial_cmp(&(b.1 - f).abs()).unwrap()).unwrap().0;
    println!("{:<24} {:>8} {:>9} {:>9} {:>18}{}", "case", "f [GHz]", "S11 [dB]", "S21 [dB]", "Zin chip side [Ω]", if target.is_some() { "  |Γ| vs target" } else { "" });
    let tcol = |z: Cx, fi: usize| gam_t(z, fi).map(|g| format!("  {g:13.3}")).unwrap_or_default();
    for (vk, (vname, parts_v)) in vars.iter().enumerate() {
        let mut parts: Vec<(usize, Part)> = Vec::new();
        let mut by_name = std::collections::HashMap::new();
        for (pn, ps) in parts_v.as_object().ok_or("variant must map port -> part")? {
            let mut spec = ps.as_str().unwrap_or("").to_string();
            // vendor files are relative to the config file
            for pre in ["s2p:", "s2p-shunt:", "s1p:"] {
                if let Some(f) = spec.strip_prefix(pre) {
                    if Path::new(f).is_relative() {
                        spec = format!("{pre}{}", base.join(f).display());
                    }
                }
            }
            let p = Part::parse(&spec)?;
            by_name.insert(pn.clone(), p.clone());
            parts.push((idx(pn)?, p));
        }
        for (k, n) in names.iter().enumerate() {
            if !keep.contains(&k) && !parts.iter().any(|p| p.0 == k) {
                return Err(format!("variant {vname:?}: EM port {n:?} has no part (use \"open\" to leave it open)"));
            }
        }
        let mut red: Vec<Vec<Vec<Cx>>> = vec![vec![Vec::with_capacity(freqs.len()); keep.len()]; keep.len()];
        let mut ideal: Vec<[[Cx; 2]; 2]> = Vec::new();
        for (fi, &f) in freqs.iter().enumerate() {
            let sf: Vec<Vec<Cx>> = (0..em.s.len()).map(|i| (0..em.s.len()).map(|j| em.s[i][j][fi]).collect()).collect();
            let t: Vec<(usize, Cx)> = parts.iter().map(|(k, p)| (*k, p.gamma(f, em.z0))).collect();
            let mut r = terminate(&sf, &keep, &t);
            if renorm {
                r = renormalise(&r, em.z0, &zr);
            }
            for i in 0..keep.len() {
                for j in 0..keep.len() {
                    red[i][j].push(r[i][j]);
                }
            }
            if !ladder_spec.is_empty() {
                let steps: Vec<(bool, Part)> = ladder_spec.iter().map(|(sh, n)| (*sh, by_name.get(n).cloned().unwrap_or(Part::Open))).collect();
                let mut l = ladder(&steps, f, em.z0);
                if renorm && keep.len() == 2 {
                    let m = vec![vec![l[0][0], l[0][1]], vec![l[1][0], l[1][1]]];
                    let m = renormalise(&m, em.z0, &zr);
                    l = [[m[0][0], m[0][1]], [m[1][0], m[1][1]]];
                }
                ideal.push(l);
            }
        }
        let safe: String = vname.chars().map(|c| if c.is_ascii_alphanumeric() { c } else { '_' }).collect();
        crate::post::write_touchstone(&out.join(format!("{safe}.s{}p", keep.len())), freqs, &red, em.z0).map_err(|e| e.to_string())?;
        let (col, col2) = if plot_ideal {
            (colours[(2 * vk) % colours.len()], colours[(2 * vk + 1) % colours.len()])
        } else {
            (colours[vk % colours.len()], colours[vk % colours.len()])
        };
        let mut rows = Vec::new();
        for &fr in &report {
            let fi = at(fr);
            let (s11, s21) = (red[0][0][fi], if keep.len() > 1 { red[1][0][fi] } else { Cx::default() });
            let z = zin(s11, zr[0].re);
            println!("{:<24} {:>8.3} {:>9.2} {:>9.2} {:>9.2} {:+8.2}j{}", format!("{vname} + layout"), freqs[fi] / 1e9, db(s11), db(s21), z.re, z.im, tcol(z, fi));
            let mut row = json!({"f": freqs[fi], "s11_db": db(s11), "s21_db": db(s21), "zin": [z.re, z.im]});
            if let Some(g) = gam_t(z, fi) {
                row["gamma_vs_target"] = json!(g);
            }
            if let Some(l) = ideal.get(fi) {
                let zi = zin(l[0][0], zr[0].re);
                if plot_ideal {
                    println!("{:<24} {:>8.3} {:>9.2} {:>9.2} {:>9.2} {:+8.2}j{}", format!("{vname}, ideal ladder"), freqs[fi] / 1e9, db(l[0][0]), db(l[1][0]), zi.re, zi.im, tcol(zi, fi));
                }
                row["ideal"] = json!({"s11_db": db(l[0][0]), "s21_db": db(l[1][0]), "zin": [zi.re, zi.im]});
            }
            rows.push(row);
        }
        out_json.insert(vname.clone(), json!({"touchstone": format!("{safe}.s{}p", keep.len()), "report": rows}));
        if keep.len() == 2 {
            curves.push((format!("{vname} + layout"), col, red[0][0].iter().map(|&x| db(x)).collect::<Vec<_>>(), red[1][0].iter().map(|&x| db(x)).collect::<Vec<_>>()));
            if !ideal.is_empty() && plot_ideal {
                curves.push((format!("{vname}, ideal"), col2, ideal.iter().map(|l| db(l[0][0])).collect(), ideal.iter().map(|l| db(l[1][0])).collect()));
            }
        }
    }
    if !curves.is_empty() {
        let c: Vec<(&str, &str, Vec<f64>, Vec<f64>)> = curves.iter().map(|(n, c, a, b)| (n.as_str(), *c, a.clone(), b.clone())).collect();
        plot_svg(&out.join("cosim.svg"), freqs, &c, &report).map_err(|e| e.to_string())?;
    }
    let kept: Vec<&String> = keep.iter().map(|&k| &names[k]).collect();
    std::fs::write(out.join("cosim.json"), serde_json::to_string_pretty(&json!({"em": em_path.to_string_lossy(), "kept": kept, "zref": zr.iter().map(|z| [z.re, z.im]).collect::<Vec<_>>(), "variants": out_json})).unwrap()).map_err(|e| e.to_string())?;
    println!("cosim: wrote {}", out.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An ideal three-port junction (S = 2/3 - I) with port 3 terminated in Z is a shunt Z between
    /// ports 1 and 2: the termination formula must agree with the ABCD ladder.
    #[test]
    fn termination_matches_ladder_shunt() {
        let s: Vec<Vec<Cx>> = (0..3).map(|i| (0..3).map(|j| Cx::new(2.0 / 3.0 - if i == j { 1.0 } else { 0.0 }, 0.0)).collect()).collect();
        for spec in ["C:1p,esr=0.3,esl=0.4n", "L:4.7n,srf=6G,q=18@2.4G", "R:20", "open", "short"] {
            let p = Part::parse(spec).unwrap();
            for f in [1e8, 2.44e9, 7e9] {
                let r = terminate(&s, &[0, 1], &[(2, p.gamma(f, 50.0))]);
                let l = ladder(&[(true, p.clone())], f, 50.0);
                for (i, j) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                    assert!((r[i][j] - l[i][j]).abs() < 1e-9, "{spec} at {f}: {:?} vs {:?}", r[i][j], l[i][j]);
                }
            }
        }
    }

    #[test]
    fn renormalise_to_own_reference_is_identity_and_matched_load_gives_zero() {
        let s = vec![vec![Cx::new(0.2, -0.1), Cx::new(0.7, 0.3)], vec![Cx::new(0.7, 0.3), Cx::new(-0.1, 0.2)]];
        let r = renormalise(&s, 50.0, &[Cx::new(50.0, 0.0), Cx::new(50.0, 0.0)]);
        for i in 0..2 {
            for j in 0..2 {
                assert!((r[i][j] - s[i][j]).abs() < 1e-12);
            }
        }
        // a one-port load Z referenced to Zr = conj(Z) is conjugate-matched: S11 = 0
        let z = Cx::new(30.0, 12.0);
        let s11 = (z - Cx::new(50.0, 0.0)) / (z + Cx::new(50.0, 0.0));
        let r = renormalise(&[vec![s11]], 50.0, &[z.conj()]);
        assert!(r[0][0].abs() < 1e-12);
    }

    #[test]
    fn touchstone_round_trip() {
        let dir = std::env::temp_dir().join(format!("fs_cosim_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let freqs = vec![1e9, 2e9];
        for n in [2usize, 3] {
            let s: Vec<Vec<Vec<Cx>>> = (0..n).map(|i| (0..n).map(|j| freqs.iter().map(|f| Cx::new(i as f64 + 0.1 * j as f64, f / 1e10)).collect()).collect()).collect();
            let p = dir.join(format!("t.s{n}p"));
            crate::post::write_touchstone(&p, &freqs, &s, 50.0).unwrap();
            let t = read_touchstone(&p).unwrap();
            assert_eq!(t.freqs, freqs);
            for i in 0..n {
                for j in 0..n {
                    for f in 0..2 {
                        assert!((t.s[i][j][f] - s[i][j][f]).abs() < 1e-7);
                    }
                }
            }
        }
    }
    /// A part measured in series or in shunt, with a matched fixture line left in the data, gives back its Z.
    #[test]
    fn vendor_series_and_shunt_files() {
        let dir = std::env::temp_dir().join(format!("fs_cosim_v_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let freqs = vec![1e9, 2e9, 3e9];
        let zpart = |f: f64| Cx::new(0.4, 2.0 * PI * f * 0.4e-9 - 1.0 / (2.0 * PI * f * 1.2e-12));
        let line = |f: f64| { let th = -2.0 * PI * f * 30e-12; Cx::new(th.cos(), th.sin()) };
        for shunt in [false, true] {
            let s: Vec<Vec<Vec<Cx>>> = (0..2).map(|i| (0..2).map(|j| freqs.iter().map(|&f| {
                let z = zpart(f);
                let z0 = Cx::new(50.0, 0.0);
                let (s11, s21) = if shunt { (Cx::default() - z0 / (z * 2.0 + z0), z * 2.0 / (z * 2.0 + z0)) } else { (z / (z + z0 * 2.0), z0 * 2.0 / (z + z0 * 2.0)) };
                (if i == j { s11 } else { s21 }) * line(f)
            }).collect()).collect()).collect();
            let p = dir.join(if shunt { "sh.s2p" } else { "se.s2p" });
            crate::post::write_touchstone(&p, &freqs, &s, 50.0).unwrap();
            let part = Part::parse(&format!("{}:{}", if shunt { "s2p-shunt" } else { "s2p" }, p.display())).unwrap();
            for &f in &freqs {
                assert!((part.z(f).unwrap() - zpart(f)).abs() < 1e-5 * zpart(f).abs(), "shunt={shunt} f={f}");
            }
        }
    }

}
