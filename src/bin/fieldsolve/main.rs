//! fieldsolve: a 3-D FDTD field solver for PCB layouts (CPU reference and Metal backend).
#![allow(dead_code)]
mod cases;
mod connect;
mod cosim;
mod fdtd;
mod geom;
mod gerber;
mod mesh;
mod model;
mod post;
mod export;
mod sim;
mod lines;
#[cfg(all(target_os = "macos", feature = "metal"))]
mod metal_fdtd;

fn arg(args: &[String], k: &str) -> Option<String> {
    args.iter().position(|a| a == k).and_then(|i| args.get(i + 1).cloned())
}

/// The backend chosen by `--backend cpu|metal|auto` (auto: Metal when available).
fn backend(args: &[String]) -> Box<dyn Fn(&model::Problem, &fdtd::RunSpec) -> fdtd::RunResult> {
    let want = arg(args, "--backend").unwrap_or_else(|| "auto".into());
    #[cfg(all(target_os = "macos", feature = "metal"))]
    {
        if want != "cpu" && metal_fdtd::available() {
            return Box::new(|p: &model::Problem, s: &fdtd::RunSpec| match metal_fdtd::run_metal(p, s) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("fieldsolve: Metal failed ({e}); running on the CPU");
                    fdtd::run_cpu(p, s)
                }
            });
        }
    }
    if want == "metal" {
        eprintln!("fieldsolve: Metal is not available; running on the CPU");
    }
    Box::new(|p: &model::Problem, s: &fdtd::RunSpec| fdtd::run_cpu(p, s))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cmd = args.get(1).map(|s| s.as_str()).unwrap_or("help");
    match cmd {
        "selftest" => {
            let scale: f64 = arg(&args, "--scale").and_then(|s| s.parse().ok()).unwrap_or(1.0);
            let which = arg(&args, "--case").unwrap_or_else(|| "all".into());
            let be = backend(&args);
            let mut worst: f64 = 0.0;
            for (name, f) in [("microstrip", cases::microstrip as fn(&cases::Runner, f64) -> cases::CaseResult), ("patch", cases::patch_capacitance), ("cavity", cases::cavity), ("box", cases::pec_box)] {
                if which != "all" && which != name {
                    continue;
                }
                let r = f(&*be, scale);
                r.print();
                worst = worst.max(r.worst());
            }
            println!("worst relative error {:.2} %", worst * 100.0);
        }
        "gap" => {
            // fast edge over a slot in the ground plane: TDR, S21 and current frames, with and without slot
            let be = backend(&args);
            let out = std::path::PathBuf::from(arg(&args, "--out").unwrap_or_else(|| "fieldsolve-gap".into()));
            let cell: f64 = arg(&args, "--cell").and_then(|v| v.parse().ok()).unwrap_or(0.1);
            let fmax: f64 = arg(&args, "--fmax").and_then(|v| v.parse().ok()).unwrap_or(20e9);
            let ro = sim::RunOpts { snap_layers: vec!["TOP".into(), "GND".into()], snap_ps: 2.0, snap_until_ps: 500.0, ..Default::default() };
            let mut curves = Vec::new();
            for slot in [false, true] {
                let st = cases::gap_setup(slot, cell, fmax);
                let dir = out.join(if slot { "slot" } else { "solid" });
                let r = sim::run_setup(&st, &dir, &*be, &ro).unwrap_or_else(|e| panic!("{e}"));
                let p0 = &r["ports"][0];
                let freqs: Vec<f64> = p0["freqs"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
                let z: Vec<post::Cx> = p0["z_re"].as_array().unwrap().iter().zip(p0["z_im"].as_array().unwrap()).map(|(a, b)| post::Cx::new(a.as_f64().unwrap(), b.as_f64().unwrap())).collect();
                let s11: Vec<post::Cx> = z.iter().map(|&z| (z - post::Cx::new(50.0, 0.0)) / (z + post::Cx::new(50.0, 0.0))).collect();
                let (t, rho, zt, tr) = post::tdr(&freqs, &s11, 50.0, 600e-12, 1200);
                curves.push((slot, t, rho, zt, tr));
            }
            sim::write_tdr(&out, &curves, cases::hammerstad(cases::microstrip_width(50.0, 0.2, 4.3), 0.2, 4.3).1).unwrap_or_else(|e| panic!("{e}"));
        }
        "ports" => {
            // port model checks on a single series loop: point vs pad-footprint ports, swapped ends, meshes
            let be = backend(&args);
            let out = std::path::PathBuf::from(arg(&args, "--out").unwrap_or_else(|| "fieldsolve-ports".into()));
            let ro = sim::RunOpts { l_band: Some((50e6, 300e6)), ..Default::default() };
            println!("{:<10} {:>6} {:>6} {:>10} {:>10}", "port", "cell", "swap", "L [nH]", "R [mohm]");
            for footprint in [false, true] {
                for cell in [0.1, 0.05, 0.025] {
                    for swap in [false, true] {
                        let st = cases::port_loop_setup(cell, footprint, swap);
                        let dir = out.join(format!("{}_{}_{}", if footprint { "pad" } else { "point" }, cell, if swap { "B" } else { "A" }));
                        let r = sim::run_setup(&st, &dir, &*be, &ro).unwrap_or_else(|e| panic!("{e}"));
                        let lf = &r["ports"][0]["l_fit"];
                        println!("{:<10} {:>6} {:>6} {:>10.4} {:>10.3}", if footprint { "pad" } else { "point" }, cell, swap, lf["l_h"].as_f64().unwrap() * 1e9, lf["r_ohm"].as_f64().unwrap() * 1e3);
                    }
                }
            }
        }
        "run" => {
            let Some(spec) = args.get(2) else {
                eprintln!("usage: fieldsolve run SPEC.json [--out DIR] [--backend cpu|metal|auto]");
                std::process::exit(2);
            };
            let (mut cx, mut out) = match sim::load(std::path::Path::new(spec)) {
                Ok(x) => x,
                Err(e) => {
                    eprintln!("fieldsolve: {e}");
                    std::process::exit(1);
                }
            };
            if let Some(o) = arg(&args, "--out") {
                out = o.into();
            }
            let be = backend(&args);
            let run = |p: &model::Problem, s: &fdtd::RunSpec| be(p, s);
            if let Err(e) = sim::run(&mut cx, &out, &run) {
                eprintln!("fieldsolve: {e}");
                std::process::exit(1);
            }
        }
        "lines" => {
            let Some(dir) = args.get(2) else {
                eprintln!("usage: fieldsolve lines RUN_DIR/fields [--freq TAG] [--field H|E] [--seeds 300] [--phases 1] [--floor-db -40] [--max-len 60] [--crowd 3]");
                std::process::exit(2);
            };
            let dir = std::path::Path::new(dir);
            let tag = arg(&args, "--freq").unwrap_or_else(|| {
                let m: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(dir.join("meta.json")).expect("meta.json")).expect("meta.json");
                m["snapshots"][0]["tag"].as_str().unwrap_or("").to_string()
            });
            let num = |k: &str, d: f64| arg(&args, k).and_then(|v| v.parse().ok()).unwrap_or(d);
            let o = lines::Opts {
                field: arg(&args, "--field").unwrap_or_else(|| "H".into()),
                tag,
                seeds: num("--seeds", 300.0) as usize,
                phases: num("--phases", 1.0) as usize,
                floor_db: num("--floor-db", -40.0),
                max_len_mm: num("--max-len", 60.0),
                crowd: num("--crowd", 3.0) as u16,
            };
            match lines::trace(dir, &o) {
                Ok(p) => eprintln!("fieldsolve lines: wrote {}", p.display()),
                Err(e) => {
                    eprintln!("fieldsolve lines: {e}");
                    std::process::exit(1);
                }
            }
        }
        "cosim" => {
            let Some(cfg) = args.get(2) else {
                eprintln!("usage: fieldsolve cosim CONFIG.json");
                std::process::exit(2);
            };
            if let Err(e) = cosim::run(std::path::Path::new(cfg)) {
                eprintln!("fieldsolve cosim: {e}");
                std::process::exit(1);
            }
        }
        "svg" => {
            let dir = args.get(2).expect("usage: fieldsolve svg GERBER_DIR OUT.svg");
            let b = gerber::load_dir(std::path::Path::new(dir)).unwrap_or_else(|e| panic!("{e}"));
            gerber::render_svg(&b, std::path::Path::new(args.get(3).map(|s| s.as_str()).unwrap_or("board.svg")));
        }
        _ => {
            eprintln!("usage:\n  fieldsolve run SPEC.json [--out DIR] [--backend auto|cpu|metal]\n  fieldsolve selftest [--case microstrip|patch|cavity|box] [--scale 1] [--backend auto|cpu|metal]\n  fieldsolve lines RUN_DIR/fields [--freq TAG] [--field H|E] [--seeds N] [--phases N]\n  fieldsolve cosim CONFIG.json\n  fieldsolve svg GERBER_DIR OUT.svg");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn auto() -> Box<dyn Fn(&model::Problem, &fdtd::RunSpec) -> fdtd::RunResult> {
        backend(&[])
    }

    /// A small lossy plane pair with a port and a short: exercises dielectric, metal, PML and lumped edges.
    fn small_problem() -> model::Problem {
        use geom::{Prim, Shape, P2};
        let rect = |x0: f64, y0: f64, x1: f64, y1: f64| Prim {
            shape: Shape::Polygon(vec![P2::new(x0, y0), P2::new(x1, y0), P2::new(x1, y1), P2::new(x0, y1)]),
            dark: true,
            net: None,
            pad: None,
        };
        let l = |a: f64, b: f64, h: f64| mesh::lines_1d(a, b, &[], &[mesh::Refine { a: a + 2.0, b: b - 2.0, h }], 1.3, 1.0, &[]).into_iter().map(|v| v * 1e-3).collect();
        let g = mesh::Grid { x: l(-4.0, 14.0, 0.3), y: l(-4.0, 10.0, 0.3), z: mesh::lines_1d(-2.0, 4.0, &[0.0, 0.4], &[mesh::Refine { a: 0.0, b: 0.4, h: 0.1 }], 1.3, 1.0, &[]).into_iter().map(|v| v * 1e-3).collect() };
        let lay = model::Layout {
            sheets: vec![model::Sheet { name: "top".into(), z: 0.4, prims: vec![rect(0.0, 0.0, 10.0, 6.0)] }, model::Sheet { name: "gnd".into(), z: 0.0, prims: vec![rect(-1.0, -1.0, 11.0, 7.0)] }],
            slabs: vec![model::Slab { z0: 0.0, z1: 0.4, eps: 4.3, tand: 0.02 }],
            barrels: vec![model::Barrel { c: P2::new(9.0, 5.0), r: 0.15, z0: 0.0, z1: 0.4 }],
            outline: vec![vec![P2::new(-2.0, -2.0), P2::new(12.0, -2.0), P2::new(12.0, 8.0), P2::new(-2.0, 8.0)]],
            copper_window: None,
        };
        let port = model::Lumped { name: "P1".into(), axis: 2, lo: [1.0, 1.0, 0.0], hi: [1.0, 1.0, 0.4], plus_at_hi: true, r: 50.0, c: 0.0, port: true };
        model::build(&g, &lay, &[port], &model::BuildOpts { pml: [6; 6], fmax: 10e9, f_loss: 2e9, cfl: 0.95 })
    }

    #[test]
    fn cpu_and_gpu_agree() {
        let p = small_problem();
        let exc = fdtd::gauss_excite(0.0, 10e9, p.dt);
        let spec = fdtd::RunSpec { excite: Some(0), exc, max_steps: 3000, end_db: -300.0, check_every: 1000, dft_freqs: vec![3e9], dft_every: 4, snap_every: 0, snap_k: vec![], snap_until: 0 };
        let a = fdtd::run_cpu(&p, &spec);
        let b = auto()(&p, &spec);
        if b.backend == "cpu" {
            eprintln!("Metal unavailable: nothing to compare");
            return;
        }
        let peak = a.v[0].iter().fold(0.0f32, |m, x| m.max(x.abs()));
        let diff = a.v[0].iter().zip(&b.v[0]).fold(0.0f32, |m, (x, y)| m.max((x - y).abs()));
        assert!(diff <= 1e-4 * peak, "port voltage differs: {diff} vs peak {peak}");
        let (fa, fb) = (&a.fields.unwrap().data[0][2], &b.fields.unwrap().data[0][2]);
        let fpk = fa.iter().fold(0.0f32, |m, x| m.max(x.abs()));
        let fdiff = fa.iter().zip(fb).fold(0.0f32, |m, (x, y)| m.max((x - y).abs()));
        assert!(fdiff <= 1e-3 * fpk, "field DFT differs: {fdiff} vs peak {fpk}");
    }

    #[test]
    fn pec_box_resonances() {
        let r = cases::pec_box(&*auto(), 1.0);
        r.print();
        assert!(r.worst() < 0.003, "worst {}", r.worst());
    }

    #[test]
    fn microstrip_against_hammerstad() {
        let r = cases::microstrip(&*auto(), 1.0);
        r.print();
        assert!((r.lines[0].1 / r.lines[0].2 - 1.0).abs() < 0.04, "Z0");
        assert!((r.lines[1].1 / r.lines[1].2 - 1.0).abs() < 0.015, "eps_eff");
    }

    #[test]
    fn cavity_modes() {
        let r = cases::cavity(&*auto(), 1.0);
        r.print();
        assert!(r.worst() < 0.015, "worst {}", r.worst());
    }

    #[test]
    fn port_swap_is_symmetric_in_a_series_loop() {
        let be = auto();
        let ro = sim::RunOpts { l_band: Some((50e6, 300e6)), ..Default::default() };
        let dir = std::env::temp_dir().join("fieldsolve-test-ports");
        let l = |swap: bool| {
            let st = cases::port_loop_setup(0.1, true, swap);
            let r = sim::run_setup(&st, &dir.join(format!("{swap}")), &*be, &ro).unwrap();
            r["ports"][0]["l_fit"]["l_h"].as_f64().unwrap()
        };
        let (a, b) = (l(false), l(true));
        assert!((a / b - 1.0).abs() < 1e-3, "{a} vs {b}");
        assert!(a > 0.9e-9 && a < 1.3e-9, "loop L {a}");
    }

    #[test]
    fn net_aware_raster_keeps_nets_apart() {
        use geom::{Prim, Shape, P2};
        // two nets 0.06 mm apart on a 0.1 mm mesh: the plain raster joins them; the check must stay clean
        let rect = |x0: f64, x1: f64, n: &str| Prim {
            shape: Shape::Polygon(vec![P2::new(x0, -1.0), P2::new(x1, -1.0), P2::new(x1, 1.0), P2::new(x0, 1.0)]),
            dark: true,
            net: Some(n.into()),
            pad: None,
        };
        let l = |a: f64, b: f64| mesh::lines_1d(a, b, &[], &[mesh::Refine { a, b, h: 0.1 }], 1.3, 0.1, &[]).into_iter().map(|v| v * 1e-3).collect();
        let g = mesh::Grid { x: l(-2.0, 2.0), y: l(-2.0, 2.0), z: vec![0.0, 0.1e-3, 0.2e-3] };
        let lay = model::Layout {
            sheets: vec![model::Sheet { name: "top".into(), z: 0.1, prims: vec![rect(-1.5, -0.01, "A"), rect(0.05, 1.5, "B")] }],
            barrels: vec![model::Barrel { c: P2::new(1.0, 0.0), r: 0.15, z0: 0.0, z1: 0.1 }],
            ..Default::default()
        };
        assert!(connect::check(&g, &lay).is_empty());
    }

    #[test]
    fn multiport_s_is_reciprocal_with_mismatched_terminations() {
        let be = auto();
        let mut st = cases::gap_setup(true, 0.4, 6e9);
        for l in &mut st.lumped {
            l.r = 10.0; // ports terminated far from the 50-ohm reference
        }
        st.freqs = post::linspace(0.2e9, 6e9, 30);
        let ro = sim::RunOpts { excite_all: true, ..Default::default() };
        let dir = std::env::temp_dir().join("fieldsolve-test-recip");
        sim::run_setup(&st, &dir, &*be, &ro).unwrap();
        let txt = std::fs::read_to_string(dir.join("sparams.s2p")).unwrap();
        for line in txt.lines().filter(|l| !l.starts_with('!') && !l.starts_with('#')) {
            let v: Vec<f64> = line.split_whitespace().map(|x| x.parse().unwrap()).collect();
            let (s21, s12) = (post::Cx::new(v[3], v[4]), post::Cx::new(v[5], v[6]));
            assert!((s21 - s12).abs() <= 0.02 * s21.abs().max(1e-3), "f {} S21 {s21:?} S12 {s12:?}", v[0]);
            // passive: no gain
            assert!(s21.abs() < 1.0 + 1e-3);
        }
    }
}