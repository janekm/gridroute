//! Tests of the sub-cell copper edge corrections (edge.rs) on both backends.

use crate::edge::EdgeOpts;
use crate::geom::{Prim, Shape, P2};
use crate::{edge_cases, fdtd, mesh, model};

fn poly(p: &[(f64, f64)]) -> Prim {
    Prim { shape: Shape::Polygon(p.iter().map(|&(x, y)| P2::new(x, y)).collect()), dark: true, net: None, pad: None }
}

fn auto() -> Box<dyn Fn(&model::Problem, &fdtd::RunSpec) -> fdtd::RunResult> {
    crate::backend(&[])
}

/// A plane pair whose top copper is a 30 degree rotated rectangle with a notch, over a ground plane, with a
/// via, a port and a lossy (optionally dispersive) slab: every edge kind the corrections touch.
fn problem(edge: Option<EdgeOpts>, wideband: Option<crate::disp::Wideband>, pml: usize) -> model::Problem {
    let l = |a: f64, b: f64, h: f64| mesh::lines_1d(a, b, &[], &[mesh::Refine { a: a + 2.0, b: b - 2.0, h }], 1.3, 1.0, &[]).into_iter().map(|v| v * 1e-3).collect();
    let g = mesh::Grid { x: l(-4.0, 14.0, 0.25), y: l(-4.0, 10.0, 0.25), z: mesh::lines_1d(-2.0, 4.0, &[0.0, 0.4], &[mesh::Refine { a: 0.0, b: 0.4, h: 0.1 }], 1.3, 1.0, &[]).into_iter().map(|v| v * 1e-3).collect() };
    let (c, s) = (30f64.to_radians().cos(), 30f64.to_radians().sin());
    let r = |x: f64, y: f64| (1.0 + x * c - y * s, 0.5 + x * s + y * c);
    let top = poly(&[r(0.0, 0.0), r(7.0, 0.0), r(7.0, 2.0), r(4.13, 2.0), r(4.13, 1.37), r(3.0, 1.37), r(3.0, 2.0), r(0.0, 2.0)]);
    let lay = model::Layout {
        sheets: vec![model::Sheet { name: "top".into(), z: 0.4, prims: vec![top] }, model::Sheet { name: "gnd".into(), z: 0.0, prims: vec![poly(&[(-1.0, -1.0), (11.0, -1.0), (11.0, 7.0), (-1.0, 7.0)])] }],
        slabs: vec![model::Slab { z0: 0.0, z1: 0.4, eps: 4.3, tand: 0.02 }],
        barrels: vec![model::Barrel { c: P2::new(r(6.0, 1.0).0, r(6.0, 1.0).1), r: 0.15, z0: 0.0, z1: 0.4 }],
        outline: vec![vec![P2::new(-2.0, -2.0), P2::new(12.0, -2.0), P2::new(12.0, 8.0), P2::new(-2.0, 8.0)]],
        copper_window: None,
    };
    let (px, py) = r(0.5, 1.0);
    let port = model::Lumped { name: "P1".into(), axis: 2, lo: [px, py, 0.0], hi: [px, py, 0.4], plus_at_hi: true, r: 50.0, c: 0.0, port: true };
    model::build(&g, &lay, &[port], &model::BuildOpts { pml: [pml; 6], fmax: 10e9, f_loss: 2e9, cfl: 0.95, wideband, conductor: None, edge })
}

#[test]
fn edge_corrections_touch_the_rotated_strip() {
    let p0 = problem(None, None, 6);
    let p = problem(Some(EdgeOpts::default()), None, 6);
    let hc = p.hcorr.as_ref().expect("H corrections");
    assert!(hc.ids.len() > 50, "only {} H corrections", hc.ids.len());
    // conformal opens metal edges and changes edge permittivities
    let opened = (0..2).map(|a| (0..p.n()).filter(|&i| p0.cb[a][i] == 0.0 && p.cb[a][i] != 0.0).count()).sum::<usize>();
    assert!(opened > 10, "opened {opened}");
    let changed = (0..3).map(|a| (0..p.n()).filter(|&i| p0.cb[a][i] != 0.0 && (p.cb[a][i] / p0.cb[a][i] - 1.0).abs() > 1e-3).count()).sum::<usize>();
    assert!(changed > 50, "changed {changed}");
    // nothing becomes metal that was not
    let added = (0..3).map(|a| (0..p.n()).filter(|&i| p0.cb[a][i] != 0.0 && p.cb[a][i] == 0.0).count()).sum::<usize>();
    assert_eq!(added, 0);
}

#[test]
fn cpu_and_gpu_agree_with_edge_corrections() {
    for wb in [None, Some(crate::disp::Wideband { f_lo: 0.3e9, f_hi: 10e9, f_ref: 1e9, copper: false })] {
        let p = problem(Some(EdgeOpts::default()), wb, 6);
        let exc = fdtd::gauss_excite(0.0, 10e9, p.dt);
        let spec = fdtd::RunSpec { excite: Some(0), exc, max_steps: 3000, end_db: -300.0, check_every: 1000, dft_freqs: vec![3e9], dft_every: 4, ..Default::default() };
        let a = fdtd::run_cpu(&p, &spec);
        let b = auto()(&p, &spec);
        if b.backend == "cpu" {
            eprintln!("Metal unavailable: nothing to compare");
            return;
        }
        let peak = a.v[0].iter().fold(0.0f32, |m, x| m.max(x.abs()));
        let diff = a.v[0].iter().zip(&b.v[0]).fold(0.0f32, |m, (x, y)| m.max((x - y).abs()));
        assert!(diff <= 1e-4 * peak, "port voltage differs: {diff} vs peak {peak}");
        for c in 0..6 {
            let (fa, fb) = (&a.fields.as_ref().unwrap().data[0][c], &b.fields.as_ref().unwrap().data[0][c]);
            let fpk = fa.iter().fold(0.0f32, |m, x| m.max(x.abs()));
            let fdiff = fa.iter().zip(fb).fold(0.0f32, |m, (x, y)| m.max((x - y).abs()));
            assert!(fdiff <= 1e-3 * fpk, "field DFT c{c} differs: {fdiff} vs peak {fpk}");
        }
    }
}

/// Late / early ratio of the largest field energy in windows of a long lossless run in a PEC box.
fn growth(edge: Option<EdgeOpts>, steps: usize) -> (f64, f64) {
    let mut p = problem(edge, None, 0);
    for a in 0..3 {
        for i in 0..p.n() {
            // lossless: no slab conductivity, no port resistor loss (cb stays: a lossless medium of slightly
            // different permittivity)
            if p.cb[a][i] != 0.0 {
                p.ca[a][i] = 1.0;
            }
        }
    }
    let exc = fdtd::gauss_excite(0.0, 10e9, p.dt);
    let n0 = exc.len() + 200;
    let spec = fdtd::RunSpec { excite: Some(0), exc, max_steps: steps, end_db: -300.0, check_every: 50, ..Default::default() };
    let r = fdtd::run_cpu(&p, &spec);
    let win = |a: usize, b: usize| r.energy.iter().filter(|e| e.0 > a && e.0 <= b).map(|e| e.1).fold(0.0, f64::max);
    let early = win(n0, n0 + 5000);
    let late = win(steps - 5000, steps);
    (late / early, p.dt)
}

#[test]
fn edge_corrections_stay_stable_in_a_closed_box() {
    // lossless, PEC walls, no absorber: the energy oscillates between E and H (the plain sum E^2 + eta^2 H^2
    // is not the conserved discrete energy) but must not grow; an unstable mode grows exponentially
    let steps = 60_000;
    let (g0, _) = growth(None, steps);
    let (g1, _) = growth(Some(EdgeOpts::default()), steps);
    eprintln!("late/early energy: plain {g0:.3}, corrected {g1:.3}");
    assert!(g1 < 1.5 && g1 < 1.5 * g0.max(1.0), "energy grew: {g1:.3} (plain {g0:.3})");
}

#[test]
#[ignore]
fn edge_growth_detects_an_unstable_setting() {
    // with tiny open fractions and no time-step check the scheme does blow up: the test above can see it
    let (g, _) = growth(Some(EdgeOpts { min_open: 0.01, kmin: 0.01, kmax: 100.0, dt_check: false, ..Default::default() }), 60_000);
    eprintln!("late/early energy without the safeguards: {g:.3e}");
    assert!(g > 10.0);
}

#[test]
fn conformal_keeps_a_thin_diagonal_trace_connected() {
    // a 45 degree trace about one cell wide: the staircase connects its ends only through corner nodes,
    // so opening every partly covered edge would cut it
    let l: Vec<f64> = mesh::lines_1d(-3.0, 9.0, &[], &[mesh::Refine { a: -1.0, b: 7.0, h: 0.1 }], 1.3, 1.0, &[]).into_iter().map(|v| v * 1e-3).collect();
    let g = mesh::Grid { x: l.clone(), y: l, z: vec![-1e-3, -0.5e-3, 0.0, 0.2e-3, 0.4e-3, 1e-3, 2e-3] };
    let w = 0.12;
    let r2 = 0.5f64.sqrt();
    let p = |s: f64, t: f64| (s * r2 - t * r2, s * r2 + t * r2);
    let tr = poly(&[p(0.0, -w / 2.0), p(8.0, -w / 2.0), p(8.0, w / 2.0), p(0.0, w / 2.0)]);
    let lay = model::Layout { sheets: vec![model::Sheet { name: "top".into(), z: 0.2, prims: vec![tr] }], slabs: vec![model::Slab { z0: -1.0, z1: 0.2, eps: 4.3, tand: 0.0 }], barrels: vec![], outline: vec![], copper_window: None };
    let pr = model::build(&g, &lay, &[], &model::BuildOpts { pml: [0; 6], edge: Some(EdgeOpts { dt_check: false, ..Default::default() }), ..Default::default() });
    let p0 = model::build(&g, &lay, &[], &model::BuildOpts { pml: [0; 6], ..Default::default() });
    let k = 3;
    let (nx, ny) = (pr.nx, pr.ny);
    let comp = |pp: &model::Problem| -> Vec<u32> {
        // flood the metal in-plane edges of plane k
        let mut lab: Vec<u32> = (0..(nx * ny) as u32).collect();
        fn find(l: &mut [u32], mut a: u32) -> u32 {
            while l[a as usize] != a {
                a = l[a as usize];
            }
            a
        }
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let id = (k * ny + j) * nx + i;
                for (a, o) in [(0usize, 1usize), (1, nx)] {
                    if pp.cb[a][id] == 0.0 {
                        let q = j * nx + i;
                        let (ra, rb) = (find(&mut lab, q as u32), find(&mut lab, (q + o) as u32));
                        lab[ra.max(rb) as usize] = ra.min(rb);
                    }
                }
            }
        }
        (0..(nx * ny) as u32).map(|q| find(&mut lab, q)).collect()
    };
    let node = |x: f64, y: f64| -> usize { mesh::snap(&g.y, y * 1e-3) * nx + mesh::snap(&g.x, x * 1e-3) };
    let (a, b) = (node(p(0.3, 0.0).0, p(0.3, 0.0).1), node(p(7.7, 0.0).0, p(7.7, 0.0).1));
    let (c0, c1) = (comp(&p0), comp(&pr));
    assert_eq!(c0[a], c0[b], "the staircase itself must connect the ends");
    assert_eq!(c1[a], c1[b], "the conformal raster cut the trace");
    let opened = (0..2).map(|aa| (0..pr.n()).filter(|&i| p0.cb[aa][i] == 0.0 && pr.cb[aa][i] != 0.0).count()).sum::<usize>();
    eprintln!("opened {opened} edges on the thin trace");
}

#[test]
fn microstrip_edges_off_the_mesh() {
    // 4 cells through the substrate, edges half a cell inside the outermost metal line: the staircase is
    // about 9 % high in eps_eff; with both corrections within 1 %
    let be = auto();
    let none = edge_cases::line(&*be, 3.0, 1.6, 4.3, 0.4, 4, 0.5, false, None);
    let both = edge_cases::line(&*be, 3.0, 1.6, 4.3, 0.4, 4, 0.5, false, Some(EdgeOpts::default()));
    let e = |r: &edge_cases::LineResult| (r.z0 / r.z0_ref - 1.0, r.ee / r.ee_ref - 1.0);
    let (en, eb) = (e(&none), e(&both));
    eprintln!("none: Z0 {:+.2} % eps {:+.2} %; both: Z0 {:+.2} % eps {:+.2} %", 100.0 * en.0, 100.0 * en.1, 100.0 * eb.0, 100.0 * eb.1);
    assert!(en.1 > 0.05, "the staircase should be visibly off ({:+.3})", en.1);
    assert!(eb.1.abs() < 0.01 && eb.0.abs() < 0.01, "corrected line off: Z0 {:+.3} eps {:+.3}", eb.0, eb.1);
}
