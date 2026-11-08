//! Analytic validation cases. Each builds a small layout directly (no Gerber), runs it on a backend and
//! compares a derived quantity with a closed-form reference.

use crate::fdtd::{gauss_excite, RunResult, RunSpec};
use crate::geom::{Prim, Shape, P2};
use crate::mesh::{lines_1d, Grid, Refine};
use crate::model::{build, BuildOpts, Layout, Lumped, Problem, Sheet, Slab, C0, ETA0};
use crate::post::{linspace, port_spectra};
use std::f64::consts::PI;

pub type Runner = dyn Fn(&Problem, &RunSpec) -> RunResult;

/// Hammerstad-Jensen microstrip (zero strip thickness): (Z0, eps_eff).
pub fn hammerstad(w: f64, h: f64, er: f64) -> (f64, f64) {
    let u = w / h;
    let a = 1.0 + (1.0 / 49.0) * ((u.powi(4) + (u / 52.0).powi(2)) / (u.powi(4) + 0.432)).ln() + (1.0 / 18.7) * (1.0 + (u / 18.1).powi(3)).ln();
    let b = 0.564 * ((er - 0.9) / (er + 3.0)).powf(0.053);
    let ee = (er + 1.0) / 2.0 + (er - 1.0) / 2.0 * (1.0 + 10.0 / u).powf(-a * b);
    let f = 6.0 + (2.0 * PI - 6.0) * (-(30.666 / u).powf(0.7528)).exp();
    let z01 = ETA0 / (2.0 * PI) * (f / u + (1.0 + 4.0 / (u * u)).sqrt()).ln();
    (z01 / ee.sqrt(), ee)
}

/// Hammerstad open-end extension of a microstrip.
pub fn open_end_dl(w: f64, h: f64, ee: f64) -> f64 {
    let u = w / h;
    0.412 * h * (ee + 0.3) * (u + 0.264) / ((ee - 0.258) * (u + 0.8))
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Prim {
    Prim { shape: Shape::Polygon(vec![P2::new(x0, y0), P2::new(x1, y0), P2::new(x1, y1), P2::new(x0, y1)]), dark: true, net: None, pad: None }
}

fn mm(v: Vec<f64>) -> Vec<f64> {
    v.into_iter().map(|x| x * 1e-3).collect()
}

pub struct CaseResult {
    pub name: String,
    pub lines: Vec<(String, f64, f64)>, // (quantity, simulated, reference)
    pub cells: usize,
    pub steps: usize,
    pub seconds: f64,
    pub backend: String,
}

impl CaseResult {
    pub fn print(&self) {
        println!("{} [{}; {} cells, {} steps, {:.1} s]", self.name, self.backend, self.cells, self.steps, self.seconds);
        for (q, s, r) in &self.lines {
            println!("  {q:<28} sim {s:>12.5}   ref {r:>12.5}   err {:+.2} %", 100.0 * (s / r - 1.0));
        }
    }
    pub fn worst(&self) -> f64 {
        self.lines.iter().map(|(_, s, r)| (s / r - 1.0).abs()).fold(0.0, f64::max)
    }
}

/// Microstrip over a ground wall (the z- boundary, PEC), fed at x = 0, open at x = len. Three voltage probes
/// (ground to strip) and three current loops a distance `d` apart in the middle of the line give, for any
/// standing wave, cosh(gamma d) = (V1 + V3) / (2 V2) and Z0^2 = (V1^2 - V3^2) / (I1^2 - I3^2), free of
/// the feed's parasitics.
pub fn microstrip(run: &Runner, scale: f64) -> CaseResult {
    let (w, h, er, len) = (3.0, 1.6, 4.3, 100.0);
    let hc = 0.2 / scale; // cell across the strip
    let (xa, d) = (35.0, 15.0);
    let x = lines_1d(-12.0, len + 12.0, &[0.0, len, xa, xa + d, xa + 2.0 * d], &[Refine { a: -1.0, b: len + 1.0, h: 0.5 / scale }], 1.3, 2.0, &[]);
    // thirds rule: the strip edge lies a fraction `fe` of a cell outside the last metal line
    let fe = 1.0 / 3.0;
    let ye = vec![-w / 2.0 + fe * hc, -w / 2.0 - (1.0 - fe) * hc, w / 2.0 - fe * hc, w / 2.0 + (1.0 - fe) * hc];
    let y = lines_1d(-14.0, 14.0, &ye, &[Refine { a: -w / 2.0 - 0.6, b: w / 2.0 + 0.6, h: hc }], 1.3, 2.0, &[]);
    let z = lines_1d(0.0, 14.0, &[h], &[Refine { a: 0.0, b: h + 0.4, h: h / (4.0 * scale) }], 1.3, 2.0, &[(0.0, h, (4.0 * scale).round() as usize)]);
    let g = Grid { x: mm(x), y: mm(y), z: mm(z) };
    let strip = Sheet { name: "top".into(), z: h, prims: vec![rect(0.0, -w / 2.0, len, w / 2.0)] };
    let lay = Layout { sheets: vec![strip], slabs: vec![Slab { z0: 0.0, z1: h, eps: er, tand: 0.0 }], barrels: vec![], outline: vec![], copper_window: None, bodies: vec![] };
    let port = Lumped { name: "P1".into(), axis: 2, lo: [0.0, -w / 2.0, 0.0], hi: [0.0, w / 2.0, h], plus_at_hi: true, r: 50.0, c: 0.0, port: true };
    let fmax = 3e9;
    let pml = [8, 8, 8, 8, 0, 8];
    let mut p = build(&g, &lay, &[port], &BuildOpts { pml, fmax, f_loss: 1e9, cfl: 0.95, wideband: None, conductor: None, edge: None });
    for q in 0..3 {
        let xq = xa + q as f64 * d;
        p.add_vprobe(&format!("v{q}"), 2, [xq, 0.0, h], 0.0);
        // the current loop sits half a cell past the voltage line; corrected below by the phase of gamma
        p.add_iprobe(&format!("i{q}"), 0, xq, [0.0, -w / 2.0 - 0.3, h - 0.1], [0.0, w / 2.0 + 0.3, h + 0.1]);
    }
    let exc = gauss_excite(0.0, fmax, p.dt);
    let r = run(&p, &RunSpec { excite: Some(0), exc, max_steps: 200_000, end_db: -70.0, check_every: 500, dft_freqs: vec![], dft_every: 1, snap_every: 0, snap_k: vec![], snap_until: 0, ..Default::default() });
    let freqs = linspace(200e6, 2.0e9, 19);
    let sp = |k: usize, toff: f64| crate::post::dft(&r.probes[k].iter().map(|&v| v as f64).collect::<Vec<f64>>(), p.dt, toff, &freqs);
    let (v1, i1, v2, _i2, v3, i3) = (sp(0, 0.0), sp(1, -0.5), sp(2, 0.0), sp(3, -0.5), sp(4, 0.0), sp(5, -0.5));
    let (z0r, eer) = hammerstad(w, h, er);
    let (mut z0s, mut ees) = (0.0, 0.0);
    for k in 0..freqs.len() {
        let ch = (v1[k] + v3[k]) / (v2[k] * 2.0);
        // lossless: cosh(j beta d) = cos(beta d)
        let bd = ch.re.clamp(-1.0, 1.0).acos();
        let beta = bd / (d * 1e-3);
        ees += (beta * C0 / (2.0 * PI * freqs[k])).powi(2);
        let z2 = (v1[k] * v1[k] - v3[k] * v3[k]) / (i1[k] * i1[k] - i3[k] * i3[k]);
        z0s += z2.sqrt().abs();
        if std::env::var("FS_DEBUG").is_ok() {
            eprintln!("  f={:.2} GHz  eps_eff={:.4}  Z0={:.3}", freqs[k] * 1e-9, (beta * C0 / (2.0 * PI * freqs[k])).powi(2), z2.sqrt().abs());
        }
    }
    let nf = freqs.len() as f64;
    CaseResult {
        name: format!("microstrip w={w} h={h} er={er} (probes, 0.2-2 GHz)"),
        lines: vec![("Z0 [ohm]".into(), z0s / nf, z0r), ("eps_eff".into(), ees / nf, eer)],
        cells: p.n(),
        steps: r.steps,
        seconds: r.seconds,
        backend: r.backend,
    }
}

/// Square patch over a ground wall, fed in its centre: C from Im(1/Z) at low frequency vs the wide-microstrip
/// capacitance with open-end extensions.
pub fn patch_capacitance(run: &Runner, scale: f64) -> CaseResult {
    let (a, h, er) = (10.0, 0.2, 4.3);
    let x = lines_1d(-10.0, a + 10.0, &[0.0, a, a / 2.0], &[Refine { a: -0.5, b: a + 0.5, h: 0.25 / scale }], 1.3, 2.0, &[]);
    let y = x.clone();
    let z = lines_1d(0.0, 10.0, &[h], &[Refine { a: 0.0, b: h + 0.2, h: h / (2.0 * scale) }], 1.3, 2.0, &[(0.0, h, (2.0 * scale).round() as usize)]);
    let g = Grid { x: mm(x), y: mm(y), z: mm(z) };
    let lay = Layout {
        sheets: vec![Sheet { name: "top".into(), z: h, prims: vec![rect(0.0, 0.0, a, a)] }],
        slabs: vec![Slab { z0: 0.0, z1: h, eps: er, tand: 0.0 }],
        barrels: vec![],
        outline: vec![],
        copper_window: None,
        bodies: vec![],
    };
    let port = Lumped { name: "P1".into(), axis: 2, lo: [a / 2.0, a / 2.0, 0.0], hi: [a / 2.0, a / 2.0, h], plus_at_hi: true, r: 50.0, c: 0.0, port: true };
    let p = build(&g, &lay, &[port], &BuildOpts { pml: [8, 8, 8, 8, 0, 8], fmax: 4e9, f_loss: 1e9, cfl: 0.95, wideband: None, conductor: None, edge: None });
    let exc = gauss_excite(0.0, 4e9, p.dt);
    let r = run(&p, &RunSpec { excite: Some(0), exc: exc.clone(), max_steps: 200_000, end_db: -70.0, check_every: 500, dft_freqs: vec![], dft_every: 1, snap_every: 0, snap_k: vec![], snap_until: 0, ..Default::default() });
    let freqs = linspace(50e6, 1000e6, 20);
    let (v, i) = port_spectra(&r.v[0], &exc, 50.0, p.dt, &freqs);
    // well below the first patch resonance (~7 GHz) the feed is a series L-C: Im Z = w L - 1 / (w C);
    // least squares in (L, 1/C)
    let (mut a11, mut a12, mut a22, mut b1, mut b2) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for (&f, (&v, &i)) in freqs.iter().zip(v.iter().zip(&i)) {
        let wv = 2.0 * PI * f;
        let (x1, x2, y) = (wv, -1.0 / wv, (v / i).im);
        a11 += x1 * x1;
        a12 += x1 * x2;
        a22 += x2 * x2;
        b1 += x1 * y;
        b2 += x2 * y;
    }
    let det = a11 * a22 - a12 * a12;
    let inv_c = (a11 * b2 - a12 * b1) / det;
    let c_sim = 1.0 / inv_c;
    let (z0, ee) = hammerstad(a, h, er);
    let dl = open_end_dl(a, h, ee);
    let c_ref = ((a + 2.0 * dl) * 1e-3) * ee.sqrt() / (C0 * z0);
    CaseResult {
        name: format!("patch capacitor {a}x{a} mm, h={h} er={er}"),
        lines: vec![("C [pF]".into(), c_sim * 1e12, c_ref * 1e12)],
        cells: p.n(),
        steps: r.steps,
        seconds: r.seconds,
        backend: r.backend,
    }
}

/// Rectangular plane pair (a x b) over a ground wall, fed near a corner: |Z| peaks at the cavity modes
/// f_mn = c / (2 sqrt(er)) * sqrt((m/a_e)^2 + (n/b_e)^2), with the effective size enlarged by the open edges.
pub fn cavity(run: &Runner, scale: f64) -> CaseResult {
    let (a, b, h, er) = (40.0, 30.0, 0.5, 4.3);
    let x = lines_1d(-8.0, a + 8.0, &[0.0, a, 3.0], &[Refine { a: -0.5, b: a + 0.5, h: 0.5 / scale }], 1.3, 2.0, &[]);
    let y = lines_1d(-8.0, b + 8.0, &[0.0, b, 3.0], &[Refine { a: -0.5, b: b + 0.5, h: 0.5 / scale }], 1.3, 2.0, &[]);
    let z = lines_1d(0.0, 10.0, &[h], &[Refine { a: 0.0, b: h + 0.3, h: h / (2.0 * scale) }], 1.3, 2.0, &[(0.0, h, (2.0 * scale).round() as usize)]);
    let g = Grid { x: mm(x), y: mm(y), z: mm(z) };
    let lay = Layout {
        sheets: vec![Sheet { name: "top".into(), z: h, prims: vec![rect(0.0, 0.0, a, b)] }],
        slabs: vec![Slab { z0: 0.0, z1: h, eps: er, tand: 0.02 }],
        barrels: vec![],
        outline: vec![],
        copper_window: None,
        bodies: vec![],
    };
    let port = Lumped { name: "P1".into(), axis: 2, lo: [3.0, 3.0, 0.0], hi: [3.0, 3.0, h], plus_at_hi: true, r: 50.0, c: 0.0, port: true };
    let p = build(&g, &lay, &[port], &BuildOpts { pml: [8, 8, 8, 8, 0, 8], fmax: 4e9, f_loss: 2e9, cfl: 0.95, wideband: None, conductor: None, edge: None });
    let exc = gauss_excite(0.0, 4e9, p.dt);
    let r = run(&p, &RunSpec { excite: Some(0), exc: exc.clone(), max_steps: 400_000, end_db: -90.0, check_every: 500, dft_freqs: vec![], dft_every: 1, snap_every: 0, snap_k: vec![], snap_until: 0, ..Default::default() });
    let freqs = linspace(1.0e9, 3.5e9, 1001);
    let (v, i) = port_spectra(&r.v[0], &exc, 50.0, p.dt, &freqs);
    let zm: Vec<f64> = v.iter().zip(&i).map(|(&v, &i)| (v / i).abs()).collect();
    // effective dimensions: each open edge adds the open-end extension of a wide microstrip
    let (_, ee_a) = hammerstad(b, h, er);
    let (_, ee_b) = hammerstad(a, h, er);
    let ae = a + 2.0 * open_end_dl(b, h, ee_a);
    let be = b + 2.0 * open_end_dl(a, h, ee_b);
    let fmn = |m: f64, n: f64| C0 / (2.0 * er.sqrt()) * ((m / (ae * 1e-3)).powi(2) + (n / (be * 1e-3)).powi(2)).sqrt();
    let mut lines = Vec::new();
    for (m, n) in [(1.0, 0.0), (0.0, 1.0), (1.0, 1.0)] {
        let fr = fmn(m, n);
        // local maximum of |Z| within +-8 % of the reference
        let (mut best, mut fb) = (0.0, 0.0);
        for (k, &f) in freqs.iter().enumerate() {
            if (f / fr - 1.0).abs() < 0.08 && zm[k] > best {
                best = zm[k];
                fb = f;
            }
        }
        lines.push((format!("f{}{} [GHz]", m as i32, n as i32), fb * 1e-9, fr * 1e-9));
    }
    CaseResult { name: format!("plane-pair cavity {a}x{b} mm, h={h} er={er}"), lines, cells: p.n(), steps: r.steps, seconds: r.seconds, backend: r.backend }
}

/// Empty PEC box (vacuum) on a graded mesh, excited by a small z-directed port off centre: |Z| peaks at the
/// TM_mn0 modes f = c/2 sqrt((m/a)^2 + (n/b)^2). Checks the bare Yee update on a non-uniform grid.
pub fn pec_box(run: &Runner, scale: f64) -> CaseResult {
    let (a, b, c) = (30.0, 20.0, 10.0);
    let x = lines_1d(0.0, a, &[7.0], &[Refine { a: 5.0, b: 9.0, h: 0.25 / scale }], 1.2, 1.0 / scale, &[]);
    let y = lines_1d(0.0, b, &[6.0], &[Refine { a: 4.0, b: 8.0, h: 0.25 / scale }], 1.2, 1.0 / scale, &[]);
    let z = lines_1d(0.0, c, &[], &[], 1.2, 1.0 / scale, &[]);
    let g = Grid { x: mm(x), y: mm(y), z: mm(z) };
    let lay = Layout::default();
    let port = Lumped { name: "P1".into(), axis: 2, lo: [7.0, 6.0, 4.0], hi: [7.0, 6.0, 6.0], plus_at_hi: true, r: 50.0, c: 0.0, port: true };
    let p = build(&g, &lay, &[port], &BuildOpts { pml: [0; 6], fmax: 20e9, f_loss: 1e9, cfl: 0.95, wideband: None, conductor: None, edge: None });
    let exc = gauss_excite(0.0, 20e9, p.dt);
    let r = run(&p, &RunSpec { excite: Some(0), exc: exc.clone(), max_steps: 60_000, end_db: -60.0, check_every: 500, dft_freqs: vec![], dft_every: 1, snap_every: 0, snap_k: vec![], snap_until: 0, ..Default::default() });
    let freqs = linspace(5e9, 16e9, 5501);
    let (v, i) = port_spectra(&r.v[0], &exc, 50.0, p.dt, &freqs);
    let zm: Vec<f64> = v.iter().zip(&i).map(|(&v, &i)| (v / i).abs()).collect();
    let fmn = |m: f64, n: f64| C0 / 2.0 * ((m / (a * 1e-3)).powi(2) + (n / (b * 1e-3)).powi(2)).sqrt();
    let mut lines = Vec::new();
    for (m, n) in [(1.0, 1.0), (2.0, 1.0), (1.0, 2.0)] {
        let fr = fmn(m, n);
        let (mut best, mut fb) = (0.0, 0.0);
        for (k, &f) in freqs.iter().enumerate() {
            if (f / fr - 1.0).abs() < 0.05 && zm[k] > best {
                best = zm[k];
                fb = f;
            }
        }
        lines.push((format!("TM{}{}0 [GHz]", m as i32, n as i32), fb * 1e-9, fr * 1e-9));
    }
    CaseResult { name: format!("PEC box {a}x{b}x{c} mm (graded mesh)"), lines, cells: p.n(), steps: r.steps, seconds: r.seconds, backend: r.backend }
}

/// A single series loop for port checks: ground plane at z = 0, FR4 0.154 mm, a 0.2 mm trace 3 mm long on
/// top between two 0.5 mm square pads. One pad gets a vertical port to the plane, the other a vertical
/// short. With nothing else present, the loop impedance must not depend on which end holds the port, and it
/// must converge with the mesh. `footprint`: port and short span the pad (80 %) instead of one node.
pub fn port_loop_setup(cell: f64, footprint: bool, swap: bool) -> crate::sim::Setup {
    let (h, er, len, pad, w) = (0.154, 4.3, 3.0, 0.5, 0.2);
    let fp = if footprint { 0.8 * pad } else { 0.0 };
    let ends = [0.0, len];
    let mut fx = vec![-pad / 2.0, pad / 2.0, len - pad / 2.0, len + pad / 2.0, 0.0, len];
    let mut fy = vec![-pad / 2.0, pad / 2.0, -w / 2.0, w / 2.0, 0.0];
    if footprint {
        for e in ends {
            fx.extend([e - fp / 2.0, e + fp / 2.0]);
        }
        fy.extend([-fp / 2.0, fp / 2.0]);
    }
    let pml_w = 8.0 * 1.0;
    let x = lines_1d(-6.0 - pml_w, len + 6.0 + pml_w, &fx, &[Refine { a: -0.6, b: len + 0.6, h: cell }], 1.3, 1.0, &[]);
    let y = lines_1d(-6.0 - pml_w, 6.0 + pml_w, &fy, &[Refine { a: -0.6, b: 0.6, h: cell }], 1.3, 1.0, &[]);
    let z = lines_1d(-3.0 - pml_w, h + 3.0 + pml_w, &[0.0, h], &[Refine { a: 0.0, b: h, h: (h / 3.0).min(cell) }], 1.3, 1.0, &[(0.0, h, 3)]);
    let grid = Grid { x: mm(x), y: mm(y), z: mm(z) };
    let top = vec![rect(-pad / 2.0, -pad / 2.0, pad / 2.0, pad / 2.0), rect(len - pad / 2.0, -pad / 2.0, len + pad / 2.0, pad / 2.0), rect(0.0, -w / 2.0, len, w / 2.0)];
    let layout = Layout {
        sheets: vec![Sheet { name: "F.Cu".into(), z: h, prims: top }, Sheet { name: "GND".into(), z: 0.0, prims: vec![rect(-5.0, -5.0, len + 5.0, 5.0)] }],
        slabs: vec![Slab { z0: 0.0, z1: h, eps: er, tand: 0.0 }],
        barrels: vec![],
        outline: vec![],
        copper_window: None,
        bodies: vec![],
    };
    let vert = |name: &str, xc: f64, port: bool| Lumped {
        name: name.into(),
        axis: 2,
        lo: [xc - fp / 2.0, -fp / 2.0, 0.0],
        hi: [xc + fp / 2.0, fp / 2.0, h],
        plus_at_hi: true,
        r: if port { 10.0 } else { 0.0 },
        c: 0.0,
        port,
    };
    let (pa, pb) = if swap { (len, 0.0) } else { (0.0, len) };
    let lumped = vec![vert("P", pa, true), vert("short", pb, false)];
    crate::sim::Setup {
        grid,
        layout,
        lumped,
        opts: BuildOpts { pml: [8; 6], f_loss: 1e9, fmax: 3e9, cfl: 0.95, wideband: None, conductor: None, edge: None },
        freqs: linspace(10e6, 3e9, 300),
        fmax: 3e9,
        end_db: -70.0,
        max_steps: 400_000,
        port_stop: false, bars: vec![], refined: vec![], checks: vec![],
    }
}

/// Width of a 50-ohm microstrip by bisection on Hammerstad-Jensen.
pub fn microstrip_width(z0: f64, h: f64, er: f64) -> f64 {
    let (mut lo, mut hi) = (0.01 * h, 20.0 * h);
    for _ in 0..60 {
        let m = 0.5 * (lo + hi);
        if hammerstad(m, h, er).0 > z0 {
            lo = m
        } else {
            hi = m
        }
    }
    0.5 * (lo + hi)
}

/// A fast edge on a 50-ohm microstrip crossing a slot in its ground plane. Trace 30 mm long on top of
/// 0.2 mm FR4 over a 16 mm wide ground plane; matched 50-ohm ports at both ends (P1 excited). With `slot`,
/// a 1 mm x 10 mm slot is cut out of the plane across the trace's path at x = 15 mm: the return current
/// under the trace has to detour around the slot's ends, which adds series inductance, a reflection and
/// a common-mode / radiating current on the plane.
pub fn gap_setup(slot: bool, cell: f64, fmax: f64) -> crate::sim::Setup {
    let (h, er, len, pw) = (0.2, 4.3, 30.0, 16.0);
    let w = microstrip_width(50.0, h, er);
    let (sx, sg, sl) = (len / 2.0, 1.0, 10.0);
    let fx = vec![-2.0, len + 2.0, 0.0, len, sx - sg / 2.0, sx + sg / 2.0];
    let fy = vec![-pw / 2.0, pw / 2.0, -sl / 2.0, sl / 2.0, -w / 2.0, w / 2.0, 0.0];
    let pml_w = 8.0;
    let x = lines_1d(-6.0 - pml_w, len + 6.0 + pml_w, &fx, &[Refine { a: -2.5, b: len + 2.5, h: cell }], 1.3, 1.0, &[]);
    let y = lines_1d(-pw / 2.0 - 4.0 - pml_w, pw / 2.0 + 4.0 + pml_w, &fy, &[Refine { a: -pw / 2.0 - 0.5, b: pw / 2.0 + 0.5, h: cell }], 1.3, 1.0, &[]);
    let z = lines_1d(-4.0 - pml_w, h + 4.0 + pml_w, &[0.0, h], &[Refine { a: 0.0, b: h, h: h / 4.0 }], 1.3, 1.0, &[(0.0, h, 4)]);
    let grid = Grid { x: mm(x), y: mm(y), z: mm(z) };
    let net = |p: Prim, n: &str| Prim { net: Some(n.into()), ..p };
    let mut gnd = vec![net(rect(-2.0, -pw / 2.0, len + 2.0, pw / 2.0), "GND")];
    if slot {
        gnd.push(Prim { dark: false, ..rect(sx - sg / 2.0, -sl / 2.0, sx + sg / 2.0, sl / 2.0) });
    }
    let layout = Layout {
        sheets: vec![Sheet { name: "TOP".into(), z: h, prims: vec![net(rect(0.0, -w / 2.0, len, w / 2.0), "SIG")] }, Sheet { name: "GND".into(), z: 0.0, prims: gnd }],
        slabs: vec![Slab { z0: 0.0, z1: h, eps: er, tand: 0.0 }],
        barrels: vec![],
        outline: vec![vec![P2::new(-2.0, -pw / 2.0), P2::new(len + 2.0, -pw / 2.0), P2::new(len + 2.0, pw / 2.0), P2::new(-2.0, pw / 2.0)]],
        copper_window: None,
        bodies: vec![],
    };
    let port = |name: &str, x: f64| Lumped { name: name.into(), axis: 2, lo: [x, -w / 2.0, 0.0], hi: [x, w / 2.0, h], plus_at_hi: true, r: 50.0, c: 0.0, port: true };
    crate::sim::Setup {
        grid,
        layout,
        lumped: vec![port("P1", 0.0), port("P2", len)],
        opts: BuildOpts { pml: [8; 6], f_loss: 1e9, fmax, cfl: 0.95, wideband: None, conductor: None, edge: None },
        freqs: linspace(0.0, fmax, 1001),
        fmax,
        end_db: -60.0,
        max_steps: 200_000,
        port_stop: false, bars: vec![], refined: vec![], checks: vec![],
    }
}

/// Centre-fed thin-wire dipole along z in free space (length `len` mm, one-cell feed gap, `cell` mm
/// cells near the wire, graded outwards), with a 50-ohm port. For the far-field checks: a lossless
/// dipole radiates all the power its port accepts; D = 1.5 when short, about 1.64 at resonance.
pub fn dipole_setup(len: f64, cell: f64, fmax: f64, air: f64) -> crate::sim::Setup {
    use crate::model::Barrel;
    let g = cell;
    let pml_w = 8.0 * (C0 / fmax / 15.0 * 1e3).min(4.0);
    let hmax = (C0 / fmax / 15.0 * 1e3).min(4.0);
    let fz = vec![-len / 2.0, -g / 2.0, g / 2.0, len / 2.0, 0.0];
    let x = lines_1d(-air - pml_w, air + pml_w, &[0.0], &[Refine { a: -3.0 * cell, b: 3.0 * cell, h: cell }], 1.25, hmax, &[]);
    let y = x.clone();
    let z = lines_1d(-len / 2.0 - air - pml_w, len / 2.0 + air + pml_w, &fz, &[Refine { a: -len / 2.0 - 3.0 * cell, b: len / 2.0 + 3.0 * cell, h: cell }], 1.25, hmax, &[]);
    let grid = Grid { x: mm(x), y: mm(y), z: mm(z) };
    let wire = |z0: f64, z1: f64| Barrel { c: P2::new(0.0, 0.0), r: 0.01, z0, z1 };
    let layout = Layout { sheets: vec![], slabs: vec![], barrels: vec![wire(-len / 2.0, -g / 2.0), wire(g / 2.0, len / 2.0)], outline: vec![], copper_window: None, bodies: vec![] };
    let port = Lumped { name: "FEED".into(), axis: 2, lo: [0.0, 0.0, -g / 2.0], hi: [0.0, 0.0, g / 2.0], plus_at_hi: true, r: 50.0, c: 0.0, port: true };
    crate::sim::Setup {
        grid,
        layout,
        lumped: vec![port],
        opts: BuildOpts { pml: [8; 6], f_loss: 1e9, fmax, cfl: 0.95, wideband: None, conductor: None, edge: None },
        freqs: linspace(fmax / 50.0, fmax, 300),
        fmax,
        end_db: -60.0,
        max_steps: 100_000,
        port_stop: false, bars: vec![], refined: vec![], checks: vec![],
    }
}

/// Directivity (linear) of a thin dipole of length `len` (m) at `f` with a sinusoidal current:
/// F(theta) = [(cos(kL/2 cos theta) - cos(kL/2)) / sin theta]^2, D = 2 F_max / int F sin theta dtheta.
pub fn dipole_directivity(len: f64, f: f64) -> f64 {
    let k = 2.0 * PI * f / C0;
    let h = k * len / 2.0;
    let ff = |t: f64| {
        let s = t.sin();
        if s.abs() < 1e-9 {
            0.0
        } else {
            ((h * t.cos()).cos() - h.cos()).powi(2) / (s * s)
        }
    };
    let n = 20000;
    let mut int = 0.0;
    let mut fmax: f64 = 0.0;
    for q in 0..n {
        let t = (q as f64 + 0.5) * PI / n as f64;
        let v = ff(t);
        fmax = fmax.max(v);
        int += v * t.sin() * PI / n as f64;
    }
    2.0 * fmax / int
}

/// One far-field check row: (case, f, directivity dBi, theory dBi, radiation efficiency, box/far ratio).
pub type FarCheck = (String, f64, f64, f64, f64, f64);

/// Run a free-space dipole with the far-field box and compare with theory.
pub fn dipole_check(be: &Runner, name: &str, len: f64, cell: f64, fmax: f64, air: f64, ff: &[f64]) -> Vec<FarCheck> {
    let st = dipole_setup(len, cell, fmax, air);
    let ro = crate::sim::RunOpts { farfield: Some(crate::nf2ff::FarOpts { freqs: ff.to_vec(), step_deg: 5.0, gap: 3, tol: 0.01 }), ..Default::default() };
    let dir = std::env::temp_dir().join(format!("fieldsolve-dipole-{name}-{}", std::process::id()));
    let res = crate::sim::run_setup(&st, &dir, be, &ro).unwrap_or_else(|e| panic!("{e}"));
    res["farfield"]
        .as_array()
        .unwrap()
        .iter()
        .map(|r| {
            let f = r["f"].as_f64().unwrap();
            let d = r["directivity_dbi"].as_f64().unwrap();
            let th = 10.0 * dipole_directivity(len * 1e-3, f).log10();
            let eff = r["radiation_efficiency"].as_f64().unwrap();
            let boxr = r["check_run_p_box"].as_f64().unwrap() / r["check_run_p_far"].as_f64().unwrap();
            (name.to_string(), f, d, th, eff, boxr)
        })
        .collect()
}

/// PEC box 30 x 20 x 10 mm filled with eps_r 4, tan(delta) `tand`, weakly probed (10 kohm port):
/// Q of TM110 / TM210 / TM120 from the half-power width of Re Z. With lossless walls Q = 1 / tan(delta)
/// at each mode's frequency, so the wideband model should give a flat Q and the constant conductivity
/// (matched at 1 GHz) a Q rising with frequency. Returns (f0, Q) per mode.
pub fn lossy_cavity(be: &Runner, wideband: bool, tand: f64) -> Vec<(f64, f64)> {
    let (a, b, c) = (30.0, 20.0, 10.0);
    let h = 0.5;
    let x = lines_1d(0.0, a, &[7.0], &[], 1.2, h, &[]);
    let y = lines_1d(0.0, b, &[6.0], &[], 1.2, h, &[]);
    let z = lines_1d(0.0, c, &[], &[], 1.2, h, &[]);
    let g = Grid { x: mm(x), y: mm(y), z: mm(z) };
    let lay = Layout { sheets: vec![], slabs: vec![Slab { z0: 0.0, z1: c, eps: 4.0, tand }], barrels: vec![], outline: vec![], copper_window: None, bodies: vec![] };
    let port = Lumped { name: "P".into(), axis: 2, lo: [7.0, 6.0, 0.0], hi: [7.0, 6.0, c], plus_at_hi: true, r: 1e4, c: 0.0, port: true };
    let wb = wideband.then(|| crate::disp::Wideband { f_lo: 0.3e9, f_hi: 10e9, f_ref: 1e9, copper: false });
    let p = build(&g, &lay, &[port], &BuildOpts { pml: [0; 6], fmax: 10e9, f_loss: 1e9, cfl: 0.95, wideband: wb, conductor: None, edge: None });
    let exc = gauss_excite(6e9, 4e9, p.dt);
    let r = be(&p, &RunSpec { excite: Some(0), exc: exc.clone(), max_steps: 60_000, end_db: -80.0, check_every: 1000, dft_every: 1, ..Default::default() });
    let freqs = linspace(3e9, 9e9, 3001);
    let (v, i) = port_spectra(&r.v[0], &exc, 1e4, p.dt, &freqs);
    let re: Vec<f64> = v.iter().zip(&i).map(|(&v, &i)| (v / i).re).collect();
    let mut out = Vec::new();
    for (lo, hi) in [(3.5e9, 5.4e9), (5.4e9, 7.0e9), (7.0e9, 8.8e9)] {
        let rng: Vec<usize> = (0..freqs.len()).filter(|&k| freqs[k] >= lo && freqs[k] <= hi).collect();
        let &kmax = rng.iter().max_by(|&&p, &&q| re[p].partial_cmp(&re[q]).unwrap()).unwrap();
        let half = re[kmax] / 2.0;
        // linear interpolation of the half-power crossings either side
        let cross = |dir: isize| -> f64 {
            let mut k = kmax as isize;
            while k + dir >= 0 && ((k + dir) as usize) < freqs.len() && re[(k + dir) as usize] > half {
                k += dir;
            }
            let (k0, k1) = (k as usize, (k + dir) as usize);
            freqs[k0] + (freqs[k1] - freqs[k0]) * (re[k0] - half) / (re[k0] - re[k1])
        };
        let (f1, f2) = (cross(-1), cross(1));
        out.push((freqs[kmax], freqs[kmax] / (f2 - f1)));
    }
    out
}

/// Copper loss against theory: a thin cavity between two copper sheets h apart (PEC side walls, lossless
/// dielectric). Its TM_mn0 modes have Q_c = h / delta from the plates alone, independent of the mode.
/// The probe port is weakly coupled (1 Mohm, loading Q ~1e5), so the measured Q is the plates' Q_c.
/// Returns per mode (f, Q with surface-impedance copper, h / delta / K_rough).
pub fn copper_cavity(be: &Runner, rough_um: f64) -> Vec<(f64, f64, f64)> {
    let (a, b, h) = (30.0, 20.0, 0.5);
    let x = lines_1d(0.0, a, &[7.0], &[], 1.2, 0.5, &[]);
    let y = lines_1d(0.0, b, &[6.0], &[], 1.2, 0.5, &[]);
    let z = lines_1d(-1.0, h + 1.0, &[0.0, h], &[Refine { a: 0.0, b: h, h: 0.1 }], 1.3, 0.5, &[]);
    let g = Grid { x: mm(x), y: mm(y), z: mm(z) };
    let plate = |zz: f64| Sheet {
        name: format!("plate{zz}"),
        z: zz,
        prims: vec![Prim { shape: Shape::Polygon(vec![P2::new(-1.0, -1.0), P2::new(a + 1.0, -1.0), P2::new(a + 1.0, b + 1.0), P2::new(-1.0, b + 1.0)]), dark: true, net: None, pad: None }],
    };
    let lay = Layout { sheets: vec![plate(h), plate(0.0)], slabs: vec![Slab { z0: 0.0, z1: h, eps: 4.0, tand: 0.0 }], barrels: vec![], outline: vec![], copper_window: None, bodies: vec![] };
    // a weakly coupled probe port (1 Mohm): its loading Q is ~1e5, so the plates set the Q
    let rp = 1e6;
    let port = Lumped { name: "P".into(), axis: 2, lo: [7.0, 6.0, 0.0], hi: [7.0, 6.0, h], plus_at_hi: true, r: rp, c: 0.0, port: true };
    let sigma = 5.8e7;
    let metal = crate::cond::Metal { sigma, t: 35e-6, rough: rough_um * 1e-6 };
    let q = |conductor: Option<crate::cond::Conductor>, steps: usize| -> Vec<(f64, f64)> {
        let p = build(&g, &lay, &[port.clone()], &BuildOpts { pml: [0; 6], fmax: 10e9, f_loss: 1e9, cfl: 0.95, wideband: None, conductor, edge: None });
        let exc = gauss_excite(6e9, 4e9, p.dt);
        let r = be(&p, &RunSpec { excite: Some(0), exc: exc.clone(), max_steps: steps, end_db: -90.0, check_every: 1000, dft_every: 1, ..Default::default() });
        let freqs = linspace(3e9, 9e9, 12001);
        let (v, i) = port_spectra(&r.v[0], &exc, rp, p.dt, &freqs);
        let re: Vec<f64> = v.iter().zip(&i).map(|(&v, &i)| (v / i).re).collect();
        let mut out = Vec::new();
        for (lo, hi) in [(3.5e9, 5.4e9), (5.4e9, 7.0e9), (7.0e9, 8.8e9)] {
            let rng: Vec<usize> = (0..freqs.len()).filter(|&k| freqs[k] >= lo && freqs[k] <= hi).collect();
            let &kmax = rng.iter().max_by(|&&p, &&q| re[p].partial_cmp(&re[q]).unwrap()).unwrap();
            let half = re[kmax] / 2.0;
            let cross = |dir: isize| -> f64 {
                let mut k = kmax as isize;
                while k + dir >= 0 && ((k + dir) as usize) < freqs.len() && re[(k + dir) as usize] > half {
                    k += dir;
                }
                let (k0, k1) = (k as usize, (k + dir) as usize);
                freqs[k0] + (freqs[k1] - freqs[k0]) * (re[k0] - half) / (re[k0] - re[k1])
            };
            out.push((freqs[kmax], freqs[kmax] / (cross(1) - cross(-1))));
        }
        out
    };
    let q1 = q(Some(crate::cond::Conductor { metals: vec![metal.clone(), metal.clone()], f_lo: 0.3e9, f_hi: 10e9 }), 3_000_000);
    q1.iter()
        .map(|&(f, qc)| {
            let delta = (1.0 / (std::f64::consts::PI * f * 1.256_637_062_12e-6 * sigma)).sqrt();
            // with roughness the plate resistance is K times higher
            let k = crate::cond::zs(&metal, f).re / crate::cond::zs(&crate::cond::Metal { rough: 0.0, ..metal.clone() }, f).re;
            (f, qc, h * 1e-3 / delta / k)
        })
        .collect()
}

