//! Validation cases for the sub-cell copper edge corrections (`edge.rs`): microstrip Z0 and εeff against
//! Hammerstad-Jensen with the strip edges at a chosen position inside a cell, the same line at 45°, and the
//! convergence of a series gap's capacitance.

use crate::cases::{hammerstad, Runner};
use crate::edge::EdgeOpts;
use crate::fdtd::{gauss_excite, RunSpec};
use crate::geom::{Prim, Shape, P2};
use crate::mesh::{lines_1d, Grid, Refine};
use crate::model::{build, BuildOpts, Layout, Lumped, Sheet, Slab, C0};
use crate::post::{dft, linspace, port_spectra, Cx};
use std::f64::consts::PI;

fn mm(v: Vec<f64>) -> Vec<f64> {
    v.into_iter().map(|x| x * 1e-3).collect()
}

fn poly(p: &[(f64, f64)]) -> Prim {
    Prim { shape: Shape::Polygon(p.iter().map(|&(x, y)| P2::new(x, y)).collect()), dark: true, net: None, pad: None }
}

/// Uniform lines `m * h` covering `[lo, hi]` (rounded outwards), as fixed lines for `lines_1d`.
fn uniform(lo: f64, hi: f64, h: f64) -> Vec<f64> {
    let (a, b) = ((lo / h).floor() as i64, (hi / h).ceil() as i64);
    (a..=b).map(|m| m as f64 * h).collect()
}

#[derive(Clone, Debug)]
pub struct LineResult {
    pub w: f64,
    pub z0: f64,
    pub ee: f64,
    pub z0_ref: f64,
    pub ee_ref: f64,
    pub cells: usize,
    pub steps: usize,
    pub seconds: f64,
    pub dt: f64,
}

impl LineResult {
    /// capacitance per length (pF/m) and its Hammerstad reference
    pub fn c(&self) -> (f64, f64) {
        (self.ee.sqrt() / (C0 * self.z0) * 1e12, self.ee_ref.sqrt() / (C0 * self.z0_ref) * 1e12)
    }
}

/// A microstrip of nominal width `w0` on `h` (εr `er`) over a PEC ground wall, 100 mm long, open at the
/// end, fed by a vertical port. The strip is meshed with uniform cells `hc` across it and `nz` cells through
/// the substrate; its width is adjusted (≤ half a cell) so that both edges lie `inset` cells inside the
/// outermost mesh line in the metal (0: edges on mesh lines; 1/3: the thirds rule). `angle` 45 turns the
/// line by 45° in the plane (the mesh stays axis-aligned, so the edges cross the cells at every offset).
/// Three voltage probes and current loops in the middle of the line give Z0 and εeff (as `cases::microstrip`).
pub fn line(run: &Runner, w0: f64, h: f64, er: f64, hc: f64, nz: usize, inset: f64, angle45: bool, edge: Option<EdgeOpts>) -> LineResult {
    let len = if angle45 { 70.0 } else { 100.0 };
    // width with both edges `inset` cells inside a mesh line: w/2 = (m + inset) hc
    let m = ((w0 / 2.0) / hc - inset).round().max(1.0);
    let w = if angle45 { w0 } else { 2.0 * (m + inset) * hc };
    let (xa, d) = (0.35 * len, 0.15 * len);
    let fmax = 3e9;
    let pml = [8, 8, 8, 8, 0, 8];
    let zl = lines_1d(0.0, h + 12.0, &[h], &[Refine { a: 0.0, b: h + 2.0 * h / nz as f64, h: h / nz as f64 }], 1.3, 2.0, &[(0.0, h, nz)]);
    let (g, strip, probes): (Grid, Prim, Vec<([f64; 3], f64, f64, f64)>) = if !angle45 {
        let xl = lines_1d(-12.0, len + 12.0, &[0.0, len, xa, xa + d, xa + 2.0 * d], &[Refine { a: -1.0, b: len + 1.0, h: 0.5 }], 1.3, 2.0, &[]);
        let yr = w / 2.0 + 4.0 * hc + 0.6;
        let yl = lines_1d(-14.0, 14.0, &uniform(-yr, yr, hc), &[Refine { a: -yr, b: yr, h: hc }], 1.3, 2.0, &[]);
        let s = poly(&[(0.0, -w / 2.0), (len, -w / 2.0), (len, w / 2.0), (0.0, w / 2.0)]);
        let pr = (0..3).map(|q| ([xa + q as f64 * d, 0.0, h], xa + q as f64 * d, -w / 2.0 - 0.3, w / 2.0 + 0.3)).collect();
        (Grid { x: mm(xl), y: mm(yl), z: mm(zl) }, s, pr)
    } else {
        // along (1, 1)/sqrt 2 from the origin; both axes uniform over the strip
        let r2 = 0.5f64.sqrt();
        let ext = len * r2;
        let mar = w / 2.0 * r2 + 4.0 * hc + 0.6;
        let fixed = uniform(-mar - 1.0, ext + mar + 1.0, hc);
        let ax = lines_1d(-12.0, ext + 12.0, &fixed, &[Refine { a: -mar - 1.0, b: ext + mar + 1.0, h: hc }], 1.3, 2.0, &[]);
        let (ux, uy) = (r2, r2);
        let (nxv, nyv) = (-r2, r2);
        let p = |s: f64, t: f64| (s * ux + t * nxv, s * uy + t * nyv);
        let s = poly(&[p(0.0, -w / 2.0), p(len, -w / 2.0), p(len, w / 2.0), p(0.0, w / 2.0)]);
        // probes at mesh nodes on the centre line: x = y = multiple of hc
        let pr = (0..3)
            .map(|q| {
                let xq = (((xa + q as f64 * d) * r2) / hc).round() * hc;
                // the loop in the plane x = xq encloses the strip's oblique cut, w sqrt 2 wide in y
                ([xq, xq, h], xq, xq - w * r2 - 0.3, xq + w * r2 + 0.3)
            })
            .collect();
        (Grid { x: mm(ax.clone()), y: mm(ax), z: mm(zl) }, s, pr)
    };
    let lay = Layout { sheets: vec![Sheet { name: "top".into(), z: h, prims: vec![strip] }], slabs: vec![Slab { z0: 0.0, z1: h, eps: er, tand: 0.0 }], barrels: vec![], outline: vec![], copper_window: None };
    // feed: a vertical port across the strip's start
    let port = if !angle45 {
        Lumped { name: "P1".into(), axis: 2, lo: [0.0, -w / 2.0 + hc, 0.0], hi: [0.0, w / 2.0 - hc, h], plus_at_hi: true, r: 50.0, c: 0.0, port: true }
    } else {
        // a single vertical line on the centre line, about 1 mm into the strip
        let x0 = (1.0 / hc).round() * hc * 0.5f64.sqrt().recip() * 0.5f64.sqrt();
        let x0 = (x0 / hc).round() * hc;
        Lumped { name: "P1".into(), axis: 2, lo: [x0, x0, 0.0], hi: [x0, x0, h], plus_at_hi: true, r: 50.0, c: 0.0, port: true }
    };
    let mut p = build(&g, &lay, &[port], &BuildOpts { pml, fmax, f_loss: 1e9, cfl: 0.95, wideband: None, conductor: None, edge });
    for (q, (at, xq, ylo, yhi)) in probes.iter().enumerate() {
        p.add_vprobe(&format!("v{q}"), 2, *at, 0.0);
        p.add_iprobe(&format!("i{q}"), 0, *xq, [0.0, *ylo, h - 0.1], [0.0, *yhi, h + 0.1]);
    }
    let dsep = if angle45 { (probes[1].1 - probes[0].1) * 2f64.sqrt() } else { d };
    let exc = gauss_excite(0.0, fmax, p.dt);
    let r = run(&p, &RunSpec { excite: Some(0), exc, max_steps: 400_000, end_db: -70.0, check_every: 500, dft_freqs: vec![], dft_every: 1, ..Default::default() });
    let freqs = linspace(200e6, 2.0e9, 19);
    let sp = |k: usize, toff: f64| dft(&r.probes[k].iter().map(|&v| v as f64).collect::<Vec<f64>>(), p.dt, toff, &freqs);
    let (v1, i1, v2, v3, i3) = (sp(0, 0.0), sp(1, -0.5), sp(2, 0.0), sp(4, 0.0), sp(5, -0.5));
    let (mut z0s, mut ees) = (0.0, 0.0);
    for k in 0..freqs.len() {
        let ch = (v1[k] + v3[k]) / (v2[k] * 2.0);
        let bd = ch.re.clamp(-1.0, 1.0).acos();
        let beta = bd / (dsep * 1e-3);
        ees += (beta * C0 / (2.0 * PI * freqs[k])).powi(2);
        let z2 = (v1[k] * v1[k] - v3[k] * v3[k]) / (i1[k] * i1[k] - i3[k] * i3[k]);
        z0s += z2.sqrt().abs();
    }
    let nf = freqs.len() as f64;
    let (z0r, eer) = hammerstad(w, h, er);
    LineResult { w, z0: z0s / nf, ee: ees / nf, z0_ref: z0r, ee_ref: eer, cells: p.n(), steps: r.steps, seconds: r.seconds, dt: p.dt }
}

/// A series gap `gap` mm in a `w` mm microstrip on `h` (εr `er`), fed by vertical ports `stub` mm from
/// each side of the gap: the gap capacitance -Im(Y21)/ω, averaged over 0.5-2 GHz (pF), and the cell count.
/// The two-port is symmetric, so one excitation gives Y.
pub fn gap(run: &Runner, w: f64, h: f64, er: f64, gap: f64, hc: f64, nz: usize, edge: Option<EdgeOpts>) -> (f64, usize, usize) {
    let stub = 2.0;
    let (xg0, xg1) = (-gap / 2.0, gap / 2.0);
    let (xp0, xp1) = (xg0 - stub, xg1 + stub);
    let fmax = 6e9;
    let xr = stub + gap / 2.0 + 0.6;
    // the gap's edges stay where they are drawn: the mesh is uniform in x through the gap, not snapped to it
    let xl = lines_1d(-xr - 8.0, xr + 8.0, &[uniform(-xr, xr, hc), vec![xp0, xp1]].concat(), &[Refine { a: -xr, b: xr, h: hc }], 1.3, 1.0, &[]);
    let yr = w / 2.0 + 0.6;
    let yl = lines_1d(-yr - 8.0, yr + 8.0, &uniform(-yr, yr, hc), &[Refine { a: -yr, b: yr, h: hc }], 1.3, 1.0, &[]);
    let zl = lines_1d(0.0, h + 8.0, &[h], &[Refine { a: 0.0, b: h + 2.0 * h / nz as f64, h: h / nz as f64 }], 1.3, 1.0, &[(0.0, h, nz)]);
    let g = Grid { x: mm(xl), y: mm(yl), z: mm(zl) };
    let a = poly(&[(xp0 - 0.3, -w / 2.0), (xg0, -w / 2.0), (xg0, w / 2.0), (xp0 - 0.3, w / 2.0)]);
    let b = poly(&[(xg1, -w / 2.0), (xp1 + 0.3, -w / 2.0), (xp1 + 0.3, w / 2.0), (xg1, w / 2.0)]);
    let lay = Layout { sheets: vec![Sheet { name: "top".into(), z: h, prims: vec![a, b] }], slabs: vec![Slab { z0: 0.0, z1: h, eps: er, tand: 0.0 }], barrels: vec![], outline: vec![], copper_window: None };
    let wp = (w / 2.0 - hc).max(0.0);
    let port = |name: &str, x: f64| Lumped { name: name.into(), axis: 2, lo: [x, -wp, 0.0], hi: [x, wp, h], plus_at_hi: true, r: 50.0, c: 0.0, port: true };
    let p = build(&g, &lay, &[port("P1", xp0), port("P2", xp1)], &BuildOpts { pml: [8, 8, 8, 8, 0, 8], fmax, f_loss: 1e9, cfl: 0.95, wideband: None, conductor: None, edge });
    let exc = gauss_excite(0.0, fmax, p.dt);
    let r = run(&p, &RunSpec { excite: Some(0), exc: exc.clone(), max_steps: 400_000, end_db: -60.0, check_every: 500, port_stop: true, dft_every: 1, ..Default::default() });
    let freqs = linspace(0.5e9, 2.0e9, 7);
    let (v1, i1) = port_spectra(&r.v[0], &exc, 50.0, p.dt, &freqs);
    let (v2, i2) = port_spectra(&r.v[1], &[], 50.0, p.dt, &freqs);
    let mut c = 0.0;
    for k in 0..freqs.len() {
        // symmetric: I1 = a V1 + b V2, I2 = b V1 + a V2
        let den = v1[k] * v1[k] - v2[k] * v2[k];
        let y21: Cx = (i2[k] * v1[k] - i1[k] * v2[k]) / den;
        c += -y21.im / (2.0 * PI * freqs[k]);
    }
    (c / freqs.len() as f64 * 1e12, p.n(), r.steps)
}

/// The edge settings compared in the validation tables.
pub fn variants() -> Vec<(&'static str, Option<EdgeOpts>)> {
    // FS_EDGE='{"reach": 1.5, ...}' overrides the defaults (for parameter scans)
    let base = std::env::var("FS_EDGE").ok().and_then(|s| serde_json::from_str::<serde_json::Value>(&s).ok()).and_then(|v| EdgeOpts::from_json(&v)).unwrap_or_default();
    vec![
        ("none", None),
        ("conformal", Some(EdgeOpts { singular: false, ..base.clone() })),
        ("singular", Some(EdgeOpts { conformal: false, ..base.clone() })),
        ("both", Some(base)),
    ]
}
