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
mod nf2ff;
mod disp;
mod cond;
mod post;
mod export;
mod sim;
mod lines;
mod edge;
mod edge_cases;
mod amr2d;
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
            // one device context for every run; a problem's buffers are reused across its excitations
            let session = std::cell::OnceCell::new();
            return Box::new(move |p: &model::Problem, s: &fdtd::RunSpec| {
                match session.get_or_init(metal_fdtd::Session::new).as_ref().map_err(|e| e.clone()).and_then(|ss| ss.run(p, s)) {
                    Ok(r) => r,
                    Err(e) => {
                        eprintln!("fieldsolve: Metal failed ({e}); running on the CPU");
                        fdtd::run_cpu(p, s)
                    }
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
        "edge-diag" => {
            // where the late-time energy sits: edge-diag SPEC [steps] [top]
            let (mut cx, _) = sim::load(std::path::Path::new(&args[2])).unwrap_or_else(|e| { eprintln!("fieldsolve: {e}"); std::process::exit(1) });
            let st = sim::setup(&mut cx).unwrap_or_else(|e| { eprintln!("fieldsolve: {e}"); std::process::exit(1) });
            let p = model::build(&st.grid, &st.layout, &st.lumped, &st.opts);
            let steps = args.get(3).and_then(|v| v.parse().ok()).unwrap_or(20000);
            let top = args.get(4).and_then(|v| v.parse().ok()).unwrap_or(30);
            edge::diag(&p, 0, steps, top);
        }
        "octree-count" => {
            // step-3 estimate: octree cells vs the rectilinear mesh for a spec's board. octree-count SPEC [r]
            let (mut cx, _) = sim::load(std::path::Path::new(&args[2])).unwrap_or_else(|e| { eprintln!("fieldsolve: {e}"); std::process::exit(1) });
            let st = sim::setup(&mut cx).unwrap_or_else(|e| { eprintln!("fieldsolve: {e}"); std::process::exit(1) });
            let (masks, _) = model::sheet_masks(&st.grid, &st.layout);
            for r in [1.3, 1.5, 2.0] {
                let (nr, no) = amr2d::octree_estimate(&st.grid, &masks, r);
                println!("grading {r}: rectilinear {:.2} M cells, octree estimate {:.2} M ({:.1} x fewer)", nr as f64 * 1e-6, no * 1e-6, nr as f64 / no);
            }
        }
        "amr-check" => {
            // dynamic-resolution prototype (amr2d.rs): plane pair with a slot, uniform reference vs static and
            // dynamic quadtrees. amr-check [h_mm] [steps] [absorb_mm]
            let num = |i: usize, d: f64| args.get(i).and_then(|v| v.parse().ok()).unwrap_or(d);
            let (h, steps, ab) = (num(2, 0.125), num(3, 60000.0) as usize, num(4, 0.0));
            let geo = amr2d::plane_pair(40.0, 30.0, 4.3, 0.02, h, 15.0, ab);
            let freqs = post::linspace(0.3e9, 5.0e9, 471);
            let mut cases: Vec<(&str, amr2d::AmrOpts)> = vec![("uniform", amr2d::AmrOpts { max_size: 1, ..Default::default() })];
            for (name, ms, ppw) in [("static ppw20", 8u32, 20.0), ("static ppw10", 16, 10.0), ("static ppw5", 32, 5.0)] {
                cases.push((name, amr2d::AmrOpts { max_size: ms, ppw, ..Default::default() }));
            }
            // refine_tol: (h k)^2 / 4 of a wave at ppw 20 is 0.025, at ppw 10 0.1
            for (name, every, rt, ct) in [("dynamic 0.1/0.025", 10usize, 0.1, 0.025), ("dynamic 0.025/0.006", 10, 0.025, 0.006), ("dynamic 0.01/0.0025", 10, 0.01, 0.0025)] {
                cases.push((name, amr2d::AmrOpts { max_size: 32, ppw: 5.0, every, refine_tol: rt, coarsen_tol: ct, ..Default::default() }));
            }
            let mut zref: Vec<post::Cx> = vec![];
            let mut pref: Vec<f64> = vec![];
            let mut uref = 0.0;
            for (name, o) in cases {
                let mut a = amr2d::Amr::new(geo.clone(), o);
                let n0 = a.leaves();
                // the same dt in every case (set by the finest cell that may appear)
                let exc = amr2d::gauss(0.0, 6e9, a.dt());
                let r = a.run(&exc, steps);
                let z = amr2d::port_z(&r, &freqs);
                let db: Vec<f64> = z.iter().map(|c| 20.0 * c.abs().log10()).collect();
                // the first cavity modes: the |Z| maximum near f10 (1.81), f01 (2.41), f11 (3.01 GHz without the slot)
                let peaks: Vec<f64> = [(1.4e9, 2.1e9), (2.1e9, 2.7e9), (2.7e9, 3.4e9), (3.4e9, 4.2e9)]
                    .iter()
                    .map(|&(lo, hi)| {
                        let k = (0..freqs.len()).filter(|&k| freqs[k] >= lo && freqs[k] <= hi).max_by(|&a, &b| db[a].partial_cmp(&db[b]).unwrap()).unwrap();
                        freqs[k] / 1e9
                    })
                    .collect();
                let probe = &r.probes[0];
                if zref.is_empty() {
                    zref = z.clone();
                    pref = probe.clone();
                    uref = r.cell_updates;
                }
                let rms_db = (db.iter().zip(&zref).map(|(a, b)| (a - 20.0 * b.abs().log10()).powi(2)).sum::<f64>() / db.len() as f64).sqrt();
                let pk = pref.iter().fold(0.0f64, |m, x| m.max(x.abs()));
                let perr = (probe.iter().zip(&pref).map(|(a, b)| (a - b).powi(2)).sum::<f64>() / pref.iter().map(|b| b * b).sum::<f64>().max(1e-300)).sqrt();
                let e_in: f64 = r.vport.iter().zip(&r.iport).zip(&exc).map(|((&v, &i), _)| v * i * r.dt).sum();
                println!(
                    "{name:<16} leaves {n0:>6} (run {}..{})  updates {:.3e} ({:.1} % of uniform; 2:1 LTS {:.1} %)  |Z| rms {:.3} dB  probe rms {:.2e} (peak {:.2e})  peaks {:?}  remeshes {}  energy lost {:.2e} of {:.2e}  {:.1} s",
                    r.min_leaves, r.max_leaves, r.cell_updates, 100.0 * r.cell_updates / uref, 100.0 * r.lts_updates / uref, rms_db, perr, pk,
                    peaks.iter().map(|f| (f * 1000.0).round() / 1000.0).collect::<Vec<_>>(), r.remeshes, r.energy_lost, e_in, r.seconds
                );
            }
        }
        "edge-check" => {
            // sub-cell copper edges (edge.rs): microstrip Z0 / eps_eff against Hammerstad at several cells
            // and edge insets, a 45 degree line, and a series gap; every edge setting
            // edge-check [line|line45|gap] [w h er] [cells,..] [insets,..] [variants,..]
            let be = backend(&args);
            let what = args.get(2).map(|s| s.as_str()).unwrap_or("line");
            let list = |i: usize, d: &str| -> Vec<f64> { args.get(i).map(|s| s.as_str()).unwrap_or(d).split(',').filter_map(|v| v.parse().ok()).collect() };
            let geo = list(3, "3,1.6,4.3");
            let (w, h, er) = (geo[0], geo[1], geo[2]);
            let cells = list(4, "0.4,0.2,0.1");
            let insets = list(5, "0,0.3333,0.5");
            let want: Vec<String> = args.get(6).map(|s| s.split(',').map(|v| v.to_string()).collect()).unwrap_or_default();
            let nzf = |hc: f64| ((h / hc).round() as usize).max(2);
            for (name, eo) in edge_cases::variants() {
                if !want.is_empty() && !want.iter().any(|v| v == name) {
                    continue;
                }
                for &hc in &cells {
                    match what {
                        "gap" => {
                            for &gp in &insets {
                                let (c, n, steps) = edge_cases::gap(&*be, w, h, er, gp, hc, nzf(hc), eo.clone());
                                println!("gap {name:<10} w {w} h {h} gap {gp:.4} cell {hc:.4}: C_gap {c:.5} pF  [{n} cells, {steps} steps]");
                            }
                        }
                        _ => {
                            let a45 = what == "line45";
                            let ins: Vec<f64> = if a45 { vec![0.0] } else { insets.clone() };
                            for &f in &ins {
                                let r = edge_cases::line(&*be, w, h, er, hc, nzf(hc), f, a45, eo.clone());
                                let (cs, cr) = r.c();
                                println!(
                                    "{what} {name:<10} cell {hc:.4} inset {f:.3} w {:.4}: Z0 {:.3} (ref {:.3}, {:+.2} %)  eps_eff {:.4} (ref {:.4}, {:+.2} %)  C {:.2} pF/m ({:+.2} %)  [{} cells, {} steps, dt {:.3e}, {:.0} s]",
                                    r.w, r.z0, r.z0_ref, 100.0 * (r.z0 / r.z0_ref - 1.0), r.ee, r.ee_ref, 100.0 * (r.ee / r.ee_ref - 1.0), cs, 100.0 * (cs / cr - 1.0), r.cells, r.steps, r.dt, r.seconds
                                );
                            }
                        }
                    }
                }
            }
        }
        "debye-port-check" => {
            // two-port microstrip (gap case without the slot) with the wideband dielectric: S must be passive
            let be = backend(&args);
            for (name, wb) in [("constant sigma", None), ("wideband", Some(disp::Wideband { f_lo: 0.3e9, f_hi: 10e9, f_ref: 1e9, copper: false }))] {
                let mut st = cases::gap_setup(false, 0.4, 6e9);
                for s in st.layout.slabs.iter_mut() {
                    s.tand = 0.015;
                }
                st.opts.wideband = wb;
                st.freqs = post::linspace(0.1e9, 6e9, 60);
                if let Some(e) = arg(&args, "--end-db").and_then(|s| s.parse::<f64>().ok()) {
                    st.end_db = e;
                    st.max_steps = 400_000;
                }
                let ro = sim::RunOpts { excite_all: true, ..Default::default() };
                let dir = std::env::temp_dir().join(format!("fieldsolve-debye-{}", name.replace(' ', "_")));
                let r = sim::run_setup(&st, &dir, &*be, &ro).unwrap();
                let txt = std::fs::read_to_string(dir.join("sparams.s2p")).unwrap();
                print!("{name:<16}");
                for line in txt.lines().filter(|l| !l.starts_with('!') && !l.starts_with('#')).step_by(10) {
                    let v: Vec<f64> = line.split_whitespace().map(|x| x.parse().unwrap()).collect();
                    let (s11, s21) = (post::Cx::new(v[1], v[2]), post::Cx::new(v[3], v[4]));
                    print!(" {:.2}GHz:{:+.3}", v[0] * 1e-9, 1.0 - s11.abs().powi(2) - s21.abs().powi(2));
                }
                println!("   steps {:?}", r["runs"].as_array().unwrap().iter().map(|x| x["steps"].as_u64().unwrap()).collect::<Vec<_>>());
            }
        }
        "copper-check" => {
            // thin copper cavity: plate-loss Q against h / delta
            let be = backend(&args);
            for rough in [0.0, 1.0] {
                println!("roughness {rough} um:");
                for (f, qc, th) in cases::copper_cavity(&*be, rough) {
                    println!("  {:.3} GHz  Q {:6.1}  theory h/delta/K {:6.1}  ({:+.1} %)", f * 1e-9, qc, th, 100.0 * (qc / th - 1.0));
                }
            }
        }
        "cond-fit" => {
            // surface impedance of copper: exact against the Foster fit the solver uses
            let num = |k: &str, d: f64| arg(&args, k).and_then(|v| v.parse::<f64>().ok()).unwrap_or(d);
            let m = cond::Metal { sigma: num("--sigma", 5.8e7), t: num("--t-um", 35.0) * 1e-6, rough: num("--rough-um", 1.0) * 1e-6 };
            let c = cond::Conductor { metals: vec![m.clone()], f_lo: num("--f-lo", 0.1e9), f_hi: num("--f-hi", 10e9) };
            let fk = cond::corners(&c);
            let ft = cond::fit(&m, &fk, c.f_lo, c.f_hi);
            println!("{} sections, worst {:.1} % over {:.2e}..{:.2e} Hz", fk.len(), 100.0 * ft.worst, c.f_lo, c.f_hi);
            println!("{:>10} {:>12} {:>12} {:>12} {:>12}", "f [Hz]", "Re exact", "Im exact", "Re fit", "Im fit");
            for e in 6..=10 {
                for m10 in [1.0, 2.5, 5.0] {
                    let f = m10 * 10f64.powi(e);
                    let (z, zf) = (cond::zs(&m, f), ft.z(f));
                    println!("{f:>10.2e} {:>12.4e} {:>12.4e} {:>12.4e} {:>12.4e}", z.re, z.im, zf.re, zf.im);
                }
            }
        }
        "loss-check" => {
            // dielectric-filled PEC cavity: modal Q against 1 / tan(delta), constant conductivity vs wideband
            let be = backend(&args);
            for tand in [0.0, 0.02] {
                for wb in [false, true] {
                    if tand == 0.0 && wb {
                        continue;
                    }
                    let r = cases::lossy_cavity(&*be, wb, tand);
                    let name = if tand == 0.0 { "lossless (port loading only)" } else if wb { "tan 0.02, wideband (Debye)" } else { "tan 0.02, constant sigma @ 1 GHz" };
                    let row: Vec<String> = r.iter().map(|(f, q)| format!("{:.2} GHz Q {:6.1}", f * 1e-9, q)).collect();
                    println!("{name:<34} {}", row.join("   "));
                }
            }
            println!("target for tan 0.02: Q = 1 / tan = 50 at every mode");
        }
        "farfield-check" => {
            // free-space dipoles: directivity against theory, radiated power against port power
            let be = backend(&args);
            println!("{:<12} {:>8} {:>10} {:>10} {:>9} {:>10}", "case", "f [GHz]", "D [dBi]", "theory", "rad. eff", "box/far");
            let mut rows = cases::dipole_check(&*be, "half-wave", 48.0, 1.0, 6e9, 30.0, &[2.6e9, 2.8e9, 3.0e9, 3.2e9]);
            rows.extend(cases::dipole_check(&*be, "short", 10.0, 0.5, 6e9, 30.0, &[1.0e9, 1.5e9]));
            for (n, f, d, th, eff, b) in rows {
                println!("{:<12} {:>8.3} {:>10.2} {:>10.2} {:>8.1}% {:>10.3}", n, f * 1e-9, d, th, eff * 100.0, b);
            }
        }
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
        #[cfg(not(all(target_os = "macos", feature = "metal")))]
        "bench" => {
            eprintln!("fieldsolve bench: times the Metal kernels; this build has no Metal backend (macOS with the \"metal\" feature)");
            std::process::exit(2);
        }
        #[cfg(all(target_os = "macos", feature = "metal"))]
        "bench" => {
            // time the Metal kernels on a spec's real problem: fixed step count, no early stop
            let Some(spec) = args.get(2) else {
                eprintln!("usage: fieldsolve bench SPEC.json [steps]");
                std::process::exit(2);
            };
            let steps: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(3000);
            let (mut cx, _) = sim::load(std::path::Path::new(spec)).unwrap_or_else(|e| { eprintln!("fieldsolve: {e}"); std::process::exit(1) });
            let st = sim::setup(&mut cx).unwrap_or_else(|e| { eprintln!("fieldsolve: {e}"); std::process::exit(1) });
            let p = model::build(&st.grid, &st.layout, &st.lumped, &st.opts);
            let exc = fdtd::gauss_excite(0.0, st.fmax, p.dt);
            let rs = fdtd::RunSpec { excite: Some(0), exc, max_steps: steps, end_db: -1000.0, check_every: 1000, dft_every: 1, ..Default::default() };
            println!("{} x {} x {} = {:.2} M cells, {} steps", p.nx, p.ny, p.nz, p.n() as f64 * 1e-6, steps);
            let session = metal_fdtd::Session::new().unwrap_or_else(|e| { eprintln!("fieldsolve: {e}"); std::process::exit(1) });
            for (fused, tg) in [("0", "32x4x1"), ("1", "32x4x1"), ("1", "64x2x1"), ("1", "32x2x2"), ("1", "16x8x1"), ("1", "128x1x1"), ("1", "32x8x1")] {
                std::env::set_var("FIELDSOLVE_FUSED", fused);
                std::env::set_var("FIELDSOLVE_TG", tg);
                match session.run(&p, &rs) {
                    Ok(r) => println!(
                        "fused {fused} tg {tg:>8}: {:7.1} s  {:6.0} Mcell/s  (set-up {:.2} s, {:.0} MB on the GPU)",
                        r.seconds,
                        p.n() as f64 * r.steps as f64 / r.seconds * 1e-6,
                        r.setup_seconds,
                        session.cached_bytes() as f64 * 1e-6
                    ),
                    Err(e) => println!("fused {fused} tg {tg}: {e}"),
                }
            }
        }
        "prims" => {
            // dump a Gerber file's primitives in drawing order as JSON (mm; dark = additive, else clear)
            let Some(f) = args.get(2) else {
                eprintln!("usage: fieldsolve prims FILE.gbr > prims.json");
                std::process::exit(2);
            };
            let text = std::fs::read_to_string(f).unwrap_or_else(|e| { eprintln!("{f}: {e}"); std::process::exit(1) });
            let prims = gerber::parse_gerber(&text).unwrap_or_else(|e| { eprintln!("{f}: {e}"); std::process::exit(1) });
            let out: Vec<serde_json::Value> = prims
                .iter()
                .map(|p| {
                    fn js(s: &geom::Shape) -> serde_json::Value {
                        match s {
                            geom::Shape::Polygon(v) => serde_json::json!({"polygon": v.iter().map(|q| [q.x, q.y]).collect::<Vec<_>>()}),
                            geom::Shape::Circle { c, r } => serde_json::json!({"circle": [c.x, c.y, *r]}),
                            geom::Shape::Capsule { a, b, r } => serde_json::json!({"capsule": [a.x, a.y, b.x, b.y, *r]}),
                            geom::Shape::Macro(parts) => serde_json::json!({"macro": parts.iter().map(|(s, on)| serde_json::json!({"exposed": on, "shape": js(s)})).collect::<Vec<_>>()}),
                        }
                    }
                    let shape = js(&p.shape);
                    serde_json::json!({"dark": p.dark, "net": p.net, "shape": shape})
                })
                .collect();
            println!("{}", serde_json::to_string(&out).unwrap());
        }
        "connect" => {
            // time the mesh connectivity check on a spec's board and print its report
            let Some(spec) = args.get(2) else {
                eprintln!("usage: fieldsolve connect SPEC.json [repeats]");
                std::process::exit(2);
            };
            let reps: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(5).max(1);
            let (mut cx, _) = sim::load(std::path::Path::new(spec)).unwrap_or_else(|e| { eprintln!("fieldsolve: {e}"); std::process::exit(1) });
            let st = sim::setup(&mut cx).unwrap_or_else(|e| { eprintln!("fieldsolve: {e}"); std::process::exit(1) });
            let (g, lay) = (&st.grid, &st.layout);
            println!("{} x {} x {} = {:.2} M nodes", g.nx(), g.ny(), g.nz(), g.n() as f64 * 1e-6);
            let mut times = Vec::new();
            let mut bridges = Vec::new();
            for _ in 0..reps {
                let t = std::time::Instant::now();
                bridges = connect::check(g, lay);
                times.push(t.elapsed().as_secs_f64());
            }
            times.sort_by(|a, b| a.partial_cmp(b).unwrap());
            println!("connect::check: median {:.1} ms, min {:.1} ms, max {:.1} ms over {reps} runs", times[reps / 2] * 1e3, times[0] * 1e3, times[reps - 1] * 1e3);
            for l in connect::describe(&bridges, cx.kicad) {
                println!("{l}");
            }
        }
        "raster" => {
            // ASCII view of the rasterised model on one copper layer: metal edges and lumped edges
            let (Some(spec), Some(layer)) = (args.get(2), args.get(3)) else {
                eprintln!("usage: fieldsolve raster SPEC.json LAYER x0 y0 x1 y1   (mm, the spec's coordinates)");
                std::process::exit(2);
            };
            let w: Vec<f64> = args[4..8.min(args.len())].iter().filter_map(|s| s.parse().ok()).collect();
            if w.len() != 4 {
                eprintln!("raster: need x0 y0 x1 y1");
                std::process::exit(2);
            }
            let (mut cx, _) = sim::load(std::path::Path::new(spec)).unwrap_or_else(|e| { eprintln!("fieldsolve: {e}"); std::process::exit(1) });
            let st = sim::setup(&mut cx).unwrap_or_else(|e| { eprintln!("fieldsolve: {e}"); std::process::exit(1) });
            let p = model::build(&st.grid, &st.layout, &st.lumped, &st.opts);
            let z = st.layout.sheets.iter().find(|s| &s.name == layer).map(|s| s.z).unwrap_or_else(|| { eprintln!("no layer {layer}"); std::process::exit(1) });
            print!("{}", sim::raster_ascii(&p, &st, z, cx.pt(w[0], w[1]), cx.pt(w[2], w[3]), cx.kicad));
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
mod review_tests;
#[cfg(test)]
mod edge_tests;

#[cfg(test)]
mod tests {
    use super::*;

    fn auto() -> Box<dyn Fn(&model::Problem, &fdtd::RunSpec) -> fdtd::RunResult> {
        backend(&[])
    }

    /// A small lossy plane pair with a port and a short: exercises dielectric, metal, PML and lumped edges.
    fn small_problem() -> model::Problem {
        small_problem_with(None)
    }

    fn small_problem_with(wideband: Option<disp::Wideband>) -> model::Problem {
        small_problem_full(wideband, None)
    }

    fn small_problem_full(wideband: Option<disp::Wideband>, conductor: Option<cond::Conductor>) -> model::Problem {
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
        model::build(&g, &lay, &[port], &model::BuildOpts { pml: [6; 6], fmax: 10e9, f_loss: 2e9, cfl: 0.95, wideband, conductor, edge: None })
    }

    #[test]
    fn cpu_and_gpu_agree() {
        let p = small_problem();
        let exc = fdtd::gauss_excite(0.0, 10e9, p.dt);
        // the far-field box faces as a gather list, at the volume-DFT frequency and one more
        let plan = nf2ff::auto_plan(&p, [6; 6], 2).unwrap();
        let spec = fdtd::RunSpec {
            excite: Some(0), exc, max_steps: 3000, end_db: -300.0, check_every: 1000, dft_freqs: vec![3e9], dft_every: 4, snap_every: 0, snap_k: vec![], snap_until: 0,
            gather_freqs: vec![3e9, 5e9], gather: plan.gather.clone(), ..Default::default()
        };
        let a = fdtd::run_cpu(&p, &spec);
        // the gathered DFT equals the volume DFT at the same samples
        for c in 0..6 {
            let vol = &a.fields.as_ref().unwrap().data[0][c];
            let pk = vol.iter().fold(0.0f32, |m, x| m.max(x.abs())).max(1e-30);
            for (t, &id) in plan.gather[c].iter().enumerate() {
                for q in 0..2 {
                    let (g, v) = (a.gathered[0][c][2 * t + q], vol[2 * id as usize + q]);
                    assert!((g - v).abs() <= 1e-5 * pk, "component {c}: gathered {g} vs volume {v}");
                }
            }
        }
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
        for fi in 0..2 {
            for c in 0..6 {
                let (ga, gb) = (&a.gathered[fi][c], &b.gathered[fi][c]);
                let pk = ga.iter().fold(0.0f32, |m, x| m.max(x.abs())).max(1e-30);
                let d = ga.iter().zip(gb).fold(0.0f32, |m, (x, y)| m.max((x - y).abs()));
                assert!(d <= 1e-3 * pk, "gathered DFT f{fi} c{c} differs: {d} vs peak {pk}");
            }
        }
    }

    #[test]
    fn cpu_and_gpu_agree_with_debye_dielectric() {
        let p = small_problem_with(Some(disp::Wideband { f_lo: 0.3e9, f_hi: 10e9, f_ref: 1e9, copper: true }));
        assert!(p.disp.is_some(), "the lossy slab should be dispersive");
        let exc = fdtd::gauss_excite(0.0, 10e9, p.dt);
        let spec = fdtd::RunSpec { excite: Some(0), exc, max_steps: 3000, end_db: -300.0, check_every: 1000, dft_freqs: vec![3e9], dft_every: 4, ..Default::default() };
        let a = fdtd::run_cpu(&p, &spec);
        let b = auto()(&p, &spec);
        if b.backend == "cpu" {
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
    fn cpu_and_gpu_agree_with_conductor_loss() {
        let m = cond::Metal { sigma: 5.8e7, t: 35e-6, rough: 1e-6 };
        let p = small_problem_full(Some(disp::Wideband { f_lo: 0.3e9, f_hi: 10e9, f_ref: 1e9, copper: false }), Some(cond::Conductor { metals: vec![m.clone(), m], f_lo: 0.3e9, f_hi: 10e9 }));
        assert!(p.cond.as_ref().is_some_and(|c| !c.ids.is_empty()), "the planes should have surface-impedance edges");
        let exc = fdtd::gauss_excite(0.0, 10e9, p.dt);
        let spec = fdtd::RunSpec { excite: Some(0), exc, max_steps: 4000, end_db: -300.0, check_every: 1000, ..Default::default() };
        let a = fdtd::run_cpu(&p, &spec);
        // stable, and lossier than the same problem with perfect copper: the late ringing is weaker
        let lossless = fdtd::run_cpu(&small_problem_full(Some(disp::Wideband { f_lo: 0.3e9, f_hi: 10e9, f_ref: 1e9, copper: false }), None), &spec);
        let peak = a.v[0].iter().fold(0.0f32, |m, x| m.max(x.abs()));
        let tail = |v: &[f32]| v[v.len() - 1000..].iter().map(|x| (x * x) as f64).sum::<f64>();
        let (t1, t0) = (tail(&a.v[0]), tail(&lossless.v[0]));
        assert!(t1.is_finite() && t1 < t0, "late energy with conductor loss {t1:.3e} vs perfect copper {t0:.3e}");
        let b = auto()(&p, &spec);
        if b.backend == "cpu" {
            return;
        }
        let diff = a.v[0].iter().zip(&b.v[0]).fold(0.0f32, |m, (x, y)| m.max((x - y).abs()));
        assert!(diff <= 1e-4 * peak, "port voltage differs: {diff} vs peak {peak}");
    }

    /// A reused Metal problem starts every run from zero: a run without a source after an excited run is
    /// silent, re-running an excitation is bit-identical, another port's excitation matches a fresh context,
    /// and a changed problem is prepared again.
    #[cfg(all(target_os = "macos", feature = "metal"))]
    #[test]
    fn metal_session_reuse_leaks_no_state() {
        if !metal_fdtd::available() {
            return;
        }
        let ss = metal_fdtd::Session::new().unwrap();
        let same = |a: &fdtd::RunResult, b: &fdtd::RunResult, what: &str| {
            assert_eq!(a.steps, b.steps, "{what}: steps");
            assert!(a.v == b.v && a.probes == b.probes, "{what}: probes differ");
            assert!(a.gathered == b.gathered, "{what}: gathered DFT differs");
            assert!(a.fields.as_ref().map(|f| &f.data) == b.fields.as_ref().map(|f| &f.data), "{what}: field DFT differs");
            assert!(a.energy == b.energy, "{what}: energy differs");
        };
        // Debye dielectric and surface-impedance copper: every kind of state buffer is in use
        let m = cond::Metal { sigma: 5.8e7, t: 35e-6, rough: 1e-6 };
        let p = small_problem_full(Some(disp::Wideband { f_lo: 0.3e9, f_hi: 10e9, f_ref: 1e9, copper: false }), Some(cond::Conductor { metals: vec![m.clone(), m], f_lo: 0.3e9, f_hi: 10e9 }));
        let plan = nf2ff::auto_plan(&p, [6; 6], 2).unwrap();
        let exc = fdtd::gauss_excite(0.0, 10e9, p.dt);
        let full = fdtd::RunSpec { excite: Some(0), exc: exc.clone(), max_steps: 1500, end_db: -300.0, check_every: 500, dft_freqs: vec![3e9], dft_every: 4, gather_freqs: vec![3e9], gather: plan.gather.clone(), ..Default::default() };
        let a = ss.run(&p, &full).unwrap();
        let quiet = fdtd::RunSpec { excite: None, ..full.clone() };
        let z = ss.run(&p, &quiet).unwrap();
        assert!(z.v.iter().flatten().all(|&x| x == 0.0), "a run without a source saw the previous run's fields");
        assert!(z.gathered.iter().flat_map(|g| g.iter().flatten()).all(|&x| x == 0.0), "gathered DFT not reset");
        assert!(z.fields.as_ref().unwrap().data.iter().flat_map(|d| d.iter().flatten()).all(|&x| x == 0.0), "field DFT not reset");
        same(&a, &ss.run(&p, &full).unwrap(), "repeated excitation");
        same(&a, &metal_fdtd::run_metal(&p, &full).unwrap(), "fresh context");
        // two ports: port 1 after port 0 (as in a multiport run, the later runs without DFTs)
        let st = cases::gap_setup(true, 0.4, 6e9);
        let q = model::build(&st.grid, &st.layout, &st.lumped, &st.opts);
        assert_eq!(q.ports.len(), 2);
        let qexc = fdtd::gauss_excite(0.0, 6e9, q.dt);
        let r0 = fdtd::RunSpec { excite: Some(0), exc: qexc, max_steps: 1500, end_db: -300.0, check_every: 500, dft_freqs: vec![3e9], dft_every: 4, ..Default::default() };
        let r1 = fdtd::RunSpec { excite: Some(1), dft_freqs: vec![], ..r0.clone() };
        ss.run(&q, &r0).unwrap();
        same(&ss.run(&q, &r1).unwrap(), &metal_fdtd::run_metal(&q, &r1).unwrap(), "second port");
        // a changed problem (same size) is not served from the cache
        let mut q2 = q.clone();
        for x in q2.cb[2].iter_mut().step_by(7) {
            *x *= 0.5;
        }
        same(&ss.run(&q2, &r1).unwrap(), &metal_fdtd::run_metal(&q2, &r1).unwrap(), "changed problem");
        assert!(ss.run(&q2, &r1).unwrap().v != ss.run(&q, &r1).unwrap().v, "the change should matter");
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