//! Tests for the fidelity policies: the objective stop, protected mesh regions, geometry edits and the
//! metadata a run records.
use crate::*;
use geom::{BBox, P2, Prim, Shape};
use serde_json::json;

fn rect(x0: f64, y0: f64, x1: f64, y1: f64, net: &str) -> Prim {
    Prim { shape: Shape::Polygon(vec![P2::new(x0, y0), P2::new(x1, y0), P2::new(x1, y1), P2::new(x0, y1)]), dark: true, net: Some(net.into()), pad: Some((net.into(), "1".into())) }
}

fn ctx(top: Vec<Prim>, spec: serde_json::Value) -> sim::Ctx {
    let plane = rect(-50., -50., 50., 50., "GND");
    let layers = vec![
        gerber::CopperLayer { name: "F.Cu".into(), index: 1, prims: top, bbox: plane.shape.bbox() },
        gerber::CopperLayer { name: "B.Cu".into(), index: 2, prims: vec![plane], bbox: BBox { x0: -50., y0: -50., x1: 50., y1: 50. } },
    ];
    let stack = vec![
        gerber::StackLayer { name: "F.Cu".into(), kind: gerber::StackKind::Copper, thickness_mm: 0.035, eps_r: 1., loss_tan: 0. },
        gerber::StackLayer { name: "FR4".into(), kind: gerber::StackKind::Dielectric, thickness_mm: 0.4, eps_r: 4.3, loss_tan: 0.02 },
        gerber::StackLayer { name: "B.Cu".into(), kind: gerber::StackKind::Copper, thickness_mm: 0.035, eps_r: 1., loss_tan: 0. },
    ];
    sim::Ctx { spec, board: gerber::Board { layers, stack, holes: vec![], outline: vec![] }, kicad: false, zs: vec![] }
}

/// For x(t) = L - A r^t sampled on checks, the tail estimate from two windows is the true remaining
/// change; a component that is not decaying (or oscillates without shrinking) never counts as settled.
#[test]
fn objective_tail_projects_a_geometric_approach() {
    let hist: Vec<(usize, Vec<f64>)> = (0..=40).map(|k| (k * 100, vec![-1.0 - 2.0 * 0.8f64.powi(k as i32), 0.5])).collect();
    let (last, tail) = fdtd::objective_tail(&hist, 500).unwrap();
    let x = |k: i32| -1.0 - 2.0 * 0.8f64.powi(k);
    assert!((last - (x(40) - x(35)).abs()).abs() < 1e-12);
    assert!((tail - (-1.0 - x(40)).abs()).abs() < 1e-9, "tail {tail} vs remaining {}", (-1.0 - x(40)).abs());
    assert!(fdtd::objective_tail(&hist[..8], 500).is_none(), "needs two full windows");
    let osc: Vec<(usize, Vec<f64>)> = (0..=40).map(|k| (k * 100, vec![0.3 * (k as f64 * 0.7).sin()])).collect();
    assert!(fdtd::objective_tail(&osc, 500).unwrap().1 > 0.1);
    // a steady drift is projected at ten times its last change; noise far below the tolerance settles
    let growing: Vec<(usize, Vec<f64>)> = (0..=40).map(|k| (k * 100, vec![k as f64 * 0.01])).collect();
    assert!((fdtd::objective_tail(&growing, 500).unwrap().1 - 0.5).abs() < 1e-9);
    let noise: Vec<(usize, Vec<f64>)> = (0..=40).map(|k| (k * 100, vec![-1.0 + 1e-7 * (k % 3) as f64])).collect();
    assert!(fdtd::objective_tail(&noise, 500).unwrap().1 < 1e-5);
}

/// With the objective stop the run ignores the energy/port rules (an empty problem satisfies the energy
/// rule at once) until the objective has settled after the minimum settling time, and records where the
/// standard rule would have stopped.
#[test]
fn objective_stop_waits_for_the_objective_and_records_the_standard_rule() {
    let g = mesh::Grid { x: vec![0., 1e-3, 2e-3], y: vec![0., 1e-3, 2e-3], z: vec![0., 1e-3, 2e-3] };
    let p = model::build(&g, &model::Layout::default(), &[], &model::BuildOpts { pml: [0; 6], ..Default::default() });
    // converges geometrically: 0.5^(steps / 10) dB away from its limit
    let metric = fdtd::ObjectiveMetric(std::sync::Arc::new(|steps, _, _| vec![-1.0 + 0.5f64.powf(steps as f64 / 10.0)]));
    let base = fdtd::RunSpec { max_steps: 1000, check_every: 5, end_db: -70., port_stop: true, objective: Some(metric), objective_stop: true, objective_tol_db: 1e-3, objective_window: 20, min_steps: 50, ..Default::default() };
    let r = fdtd::run_cpu(&p, &base);
    assert_eq!(r.stop, "objective");
    // settled when the window's spread 3 * 0.5^(t / 10) <= 1e-3: the check at 120, after min_steps
    assert_eq!(r.steps, 120, "stopped at {}", r.steps);
    assert_eq!(r.standard_stop.as_ref().map(|s| s.0), Some(5), "the energy rule fires at the first check");
    assert_eq!(r.objective.len(), r.steps / 5);
    // never settles: max_steps, reported as such
    let never = fdtd::ObjectiveMetric(std::sync::Arc::new(|steps, _, _| vec![steps as f64 * 1e-2]));
    let r = fdtd::run_cpu(&p, &fdtd::RunSpec { objective: Some(never), max_steps: 200, ..base.clone() });
    assert_eq!((r.stop.as_str(), r.steps), ("max_steps", 200));
    // without the objective stop the history is still recorded and the standard rule ends the run
    let r = fdtd::run_cpu(&p, &fdtd::RunSpec { objective_stop: false, ..base });
    assert_eq!((r.stop.as_str(), r.steps), ("energy", 5));
    assert_eq!(r.objective.len(), 1);
    assert!(r.standard_stop.is_none());
}

/// A coarse policy keeps the spec's cell size around every port (and the listed features) and coarsens the
/// rest of the fine region.
#[test]
fn protected_features_keep_their_cell_in_a_coarse_mesh() {
    let pads = vec![rect(-4.2, -0.2, -3.8, 0.2, "A"), rect(-3.2, -0.2, -2.8, 0.2, "B"), rect(3.0, 2.0, 3.4, 2.4, "C")];
    let spec = json!({"roi": [-6, -6, 6, 6], "cell": 0.05, "air": 1, "air_z": 1, "max_cell": 0.5, "pml": 2,
                      "ports": [{"name": "P", "pos": "A.1", "neg": "B.1", "layer": "F.Cu"}],
                      "objective": {"port": "P", "band": [2.4e9]},
                      "fidelity": {"protect": ["C.1"]}});
    let fast = fidelity::resolve(&spec, "screening").unwrap();
    let st = sim::setup(&mut ctx(pads.clone(), fast.clone())).unwrap();
    let full = sim::setup(&mut ctx(pads, spec)).unwrap();
    let step_at = |l: &[f64], x: f64| {
        let i = l.partition_point(|&v| v < x * 1e-3).clamp(1, l.len() - 1);
        (l[i] - l[i - 1]) * 1e3
    };
    for (x, y) in [(-4.0, 0.0), (-3.0, 0.0), (-3.5, 0.0), (3.2, 2.2)] {
        assert!(step_at(&st.grid.x, x) <= 0.05 + 1e-6 && step_at(&st.grid.y, y) <= 0.05 + 1e-6, "protected point ({x}, {y}) coarsened");
    }
    for (x, y) in [(0.5, -4.0), (-5.0, 5.0)] {
        let (sx, sy) = (step_at(&st.grid.x, x), step_at(&st.grid.y, y));
        assert!(sx > 0.08 && sy > 0.08, "({x}, {y}) should be coarse: {sx} {sy}");
    }
    assert!(st.grid.n() * 2 < full.grid.n(), "coarse {} vs full {} cells", st.grid.n(), full.grid.n());
    let what: Vec<&str> = st.refined.iter().map(|r| r.2.as_str()).collect();
    assert_eq!(what, ["port P", "pad C.1"]);
}

/// Geometry edits draw over the Gerber copper in order: a clear cuts it, an add puts copper (with its
/// net) back.
#[test]
fn geometry_edits_clear_and_add_copper() {
    let mut b = ctx(vec![rect(-1.0, -1.0, 1.0, 1.0, "ANT")], json!({})).board;
    let spec = json!({"edits": [{"layer": "F.Cu", "op": "clear", "rect": [0.0, -2.0, 2.0, 2.0]},
                                {"layer": "F.Cu", "op": "add", "rect": [0.5, -0.1, 3.0, 0.1], "net": "ANT"}]});
    sim::apply_edits(&mut b, &spec, true).unwrap();
    let covered = |x: f64, y: f64| {
        let mut on = None;
        for p in &b.layers[0].prims {
            if p.shape.contains(P2::new(x, y)) {
                on = if p.dark { p.net.clone() } else { None };
            }
        }
        on
    };
    assert_eq!(covered(-0.5, 0.5).as_deref(), Some("ANT"));
    assert_eq!(covered(0.3, 0.5), None, "cleared");
    assert_eq!(covered(2.5, 0.0).as_deref(), Some("ANT"), "added (KiCad y is negated: the strip is symmetric)");
    assert!(b.layers[0].bbox.x1 >= 3.0);
    assert!(sim::apply_edits(&mut b, &json!({"edits": [{"layer": "In1.Cu", "op": "add", "rect": [0, 0, 1, 1]}]}), true).is_err());
}

/// An added rectangle with a net, laid over another net's copper, is reported: the net-aware raster would
/// cut it off (an inverted-F's short drawn as net ANT over the GND pour floated). The same net, no
/// overlap, or copper cleared first are not reported.
#[test]
fn an_edit_over_another_nets_copper_is_reported() {
    // board frame (no KiCad flip): the pour spans y 0..2
    let b = ctx(vec![rect(0.0, 0.0, 10.0, 2.0, "GND")], json!({})).board;
    let over = |r: (f64, f64, f64, f64), net: &str| sim::edit_net_overlaps(&b.layers[0].prims, r, net);
    assert_eq!(over((4.0, 1.7, 4.5, 2.5), "ANT"), ["GND"], "the short's foot lies in the pour");
    assert!(over((4.0, 1.7, 4.5, 2.5), "GND").is_empty(), "same net");
    assert!(over((4.0, 2.3, 4.5, 6.0), "ANT").is_empty(), "clear of the pour");
    let mut b2 = ctx(vec![rect(0.0, 0.0, 10.0, 2.0, "GND")], json!({})).board;
    sim::apply_edits(&mut b2, &json!({"edits": [{"layer": "F.Cu", "op": "clear", "rect": [3.0, 1.5, 6.0, 3.0]}]}), false).unwrap();
    assert!(sim::edit_net_overlaps(&b2.layers[0].prims, (4.0, 1.7, 4.5, 2.5), "ANT").is_empty(), "cleared first");
}

/// A screening run records its policy, every override, the actual mesh and time step, the excitation
/// count, each run's stop reason (and where the standard rule would have stopped), the objective with its
/// stability, and earns its level as status only when every excitation converged.
#[test]
fn screening_run_records_policy_metadata() {
    let mut st = cases::gap_setup(false, 0.4, 6e9);
    st.freqs = post::linspace(0.2e9, 6e9, 30);
    st.max_steps = 20_000;
    let spec = json!({"objective": {"port": "P1", "band": [4e9, 5e9]}, "farfield": null});
    let rec = fidelity::resolve(&json!({"cell": 0.4, "objective": spec["objective"]}), "screening").unwrap();
    let mut obj = sim::ObjectiveOpts::from_spec(&rec).unwrap();
    obj.window_periods = 1.0;
    obj.min_settle_periods = 1.0;
    obj.tol_db = 0.2;
    let ro = sim::RunOpts { excite: vec!["P1".into()], objective: Some(obj), no_scene: true, fidelity: rec["fidelity_record"].clone(), ..Default::default() };
    let dir = std::env::temp_dir().join(format!("fieldsolve-test-fidelity-{}", std::process::id()));
    let be = crate::backend(&[]);
    let r = sim::run_setup(&st, &dir, &*be, &ro).unwrap();
    assert!(!dir.join("scene.json").exists(), "no scene export");
    let f = &r["fidelity"];
    assert_eq!(f["level"], "screening");
    assert!(f["overrides"].as_array().unwrap().iter().any(|o| o["key"] == "cell"));
    let a = &f["actual"];
    assert_eq!(a["excitations"], 1);
    assert_eq!(a["mesh"]["dims"].as_array().unwrap().len(), 3);
    assert!(a["mesh"]["min_step_mm"].as_array().unwrap().iter().all(|v| v.as_f64().unwrap() > 0.0));
    assert!(a["dt"].as_f64().unwrap() > 0.0);
    for run in r["runs"].as_array().unwrap() {
        assert!(["objective", "max_steps"].contains(&run["stop"].as_str().unwrap()), "{}", run["stop"]);
        assert!(run["objective_history"].as_array().unwrap().len() > 3);
        assert_eq!(run["objective_components"].as_array().unwrap().len(), 2 * 2);
        assert!(run.get("standard_rule_step").is_some());
    }
    let status = f["status"].as_str().unwrap();
    let all_conv = a["stops"].as_array().unwrap().iter().all(|s| s != "max_steps");
    assert_eq!(status, if all_conv { "screening" } else { "incomplete" });
    let o = &r["objective"];
    assert_eq!(o["port"], "P1");
    assert_eq!(o["per_f"].as_array().unwrap().len(), 2);
    assert!(o["worst_db"].as_f64().unwrap() <= 0.0);
    assert!(o["stability_db"].as_f64().is_some());
    // on a matched line the accepted power is ~ all of it
    assert!(o["per_f"][0]["accepted_db"].as_f64().unwrap() > -0.5, "{o}");
}

/// `optimise` on a parameterised spec: the edit's expression is resolved per point, screening and full
/// solves land in the cache, and a repeated search solves nothing.
#[test]
fn optimise_solves_parameterised_specs_once() {
    let pads = vec![rect(-4.2, -0.2, -3.8, 0.2, "A"), rect(-3.2, -0.2, -2.8, 0.2, "B")];
    let spec = json!({"roi": [-8, -8, 8, 8], "cell": 0.5, "air": 4, "air_z": 4, "max_cell": 2.0, "pml": 4, "fmax": 6e9, "max_steps": 3000,
                      "ports": [{"name": "P", "pos": "A.1", "neg": "B.1", "layer": "F.Cu"}],
                      "farfield": {"freqs": [4e9, 5e9], "step_deg": 30},
                      "objective": {"port": "P", "band": [4e9, 5e9]},
                      "edits": [{"layer": "F.Cu", "op": "add", "rect": [-2.9, -0.2, "=-2.9 + len", 0.2], "net": "B"}]});
    let board = ctx(pads, json!({})).board;
    let ps = optimise::parse_params(&json!([{"name": "len", "min": 2.0, "max": 6.0, "step": 1.0, "start": 4.0}])).unwrap();
    let dir = std::env::temp_dir().join(format!("fieldsolve-test-optimise-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let be = crate::backend(&[]);
    let mut ev = optimise::SimEval { spec, board, kicad: false, ps: ps.clone(), cache: dir.join("cache"), backend: &*be };
    let o = optimise::Opts { budget: 3, init: 3, promote_every: 0, final_n: 1, ..Default::default() };
    let r = optimise::search(&ps, &mut ev, &o).unwrap();
    assert_eq!((r["screening_solves"].as_u64(), r["full_solves"].as_u64(), r["cached"].as_u64()), (Some(3), Some(1), Some(0)), "{}", r["evaluations"]);
    assert!(r["cell_steps"]["screening"].as_f64().unwrap() > 0.0 && r["cell_steps"]["full"].as_f64().unwrap() > 0.0);
    assert_eq!(r["best"]["full_per_f_db"].as_array().unwrap().len(), 2, "one value per band frequency");
    let s: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(dir.join("cache/screening/len+4.0000/spec.json")).unwrap()).unwrap();
    assert!((s["edits"][0]["rect"][2].as_f64().unwrap() - 1.1).abs() < 1e-12, "{}", s["edits"][0]);
    assert_eq!(s["fidelity_record"]["level"], "screening");
    let again = optimise::search(&ps, &mut ev, &o).unwrap();
    assert_eq!(again["cached"].as_u64(), Some(4));
    assert_eq!(again["cell_steps"]["total"].as_f64(), Some(0.0));
    assert_eq!(again["best"], r["best"]);
    let _ = std::fs::remove_dir_all(&dir);
}

/// The scene export's capsule outline is convex: caps on the outer ends, area of a stadium (it once put
/// both caps on the inner side, which bit half-discs out of every track end in drawings and in openEMS's
/// polygon mode).
#[test]
fn capsule_outline_is_a_stadium() {
    let (a, b, r) = (P2::new(1.0, 2.0), P2::new(4.0, 6.0), 0.25);
    let v = sim::capsule_outline(a, b, r);
    let area = 0.5 * (0..v.len()).map(|i| { let (p, q) = (v[i], v[(i + 1) % v.len()]); p.x * q.y - q.x * p.y }).sum::<f64>().abs();
    let want = 5.0 * 2.0 * r + std::f64::consts::PI * r * r;
    assert!((area - want).abs() < 0.01 * want, "area {area} vs {want}");
    // every point lies on the capsule's boundary, and both ends reach out past a and b
    let d = |p: P2| { let t = (((p.x - a.x) * 3.0 + (p.y - a.y) * 4.0) / 25.0).clamp(0.0, 1.0); ((p.x - a.x - 3.0 * t).powi(2) + (p.y - a.y - 4.0 * t).powi(2)).sqrt() };
    assert!(v.iter().all(|&p| (d(p) - r).abs() < 1e-9));
    let along = |p: P2| ((p.x - a.x) * 3.0 + (p.y - a.y) * 4.0) / 5.0;
    assert!(v.iter().any(|&p| along(p) > 5.0 + 0.99 * r) && v.iter().any(|&p| along(p) < -0.99 * r));
}

/// A printed dipole 0.4 mm over a ground plane: a high-Q radiator that rings for many periods.
fn ringing_dipole_spec(extra: serde_json::Value) -> (sim::Ctx, serde_json::Value) {
    let pads = vec![rect(-0.6, -0.3, -0.1, 0.3, "A"), rect(0.1, -0.3, 0.6, 0.3, "B"), rect(-9.0, -0.3, -0.1, 0.3, "A"), rect(0.1, -0.3, 9.0, 0.3, "B")];
    let mut spec = json!({"roi": [-10, -3, 10, 3], "cell": 0.3, "air": 6, "air_z": 6, "max_cell": 2.0, "pml": 6, "fmax": 8e9, "max_steps": 60000,
                          "end_db": -80, "stop": "ports",
                          "ports": [{"name": "P", "pos": "A.1", "neg": "B.1", "layer": "F.Cu"}],
                          "farfield": {"freqs": [4.5e9, 5.0e9], "step_deg": 30},
                          "objective": {"port": "P", "band": [4.5e9, 5.0e9]}});
    for (k, v) in extra.as_object().unwrap() {
        spec[k] = v.clone();
    }
    (ctx(pads, spec.clone()), spec)
}

fn run_ringing(extra: serde_json::Value, edit: impl Fn(&mut serde_json::Value), tag: &str, cpu: bool) -> serde_json::Value {
    let (mut cx, mut spec) = ringing_dipole_spec(extra);
    edit(&mut spec);
    cx.spec = spec.clone();
    let st = sim::setup(&mut cx).unwrap();
    let mut ro = sim::RunOpts::from_spec(&spec);
    ro.no_scene = true;
    let dir = std::env::temp_dir().join(format!("fieldsolve-test-extrap-{tag}-{}", std::process::id()));
    let be = crate::backend(&if cpu { vec!["--backend".to_string(), "cpu".to_string()] } else { vec![] });
    let r = sim::run_setup(&st, &dir, &*be, &ro).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    r
}

/// Late-time extrapolation: on a high-Q radiator the run stops once the extrapolated objective is stable,
/// at a fraction of the steps the port-decay rule needs, and the result (computed from the extrapolated
/// port voltages and far-field accumulators) matches the long run's.
#[test]
fn extrapolated_stop_matches_a_long_run() {
    // FIELDSOLVE_TEST_CPU=1 runs it on the CPU backend (slow, about 3 min)
    let cpu = std::env::var("FIELDSOLVE_TEST_CPU").is_ok();
    let long = run_ringing(json!({}), |_| {}, "long", cpu);
    let fast = run_ringing(json!({"stop": "objective"}), |s| s["objective"]["extrapolate"] = json!(true), "fast", cpu);
    let (rl, rf) = (&long["runs"][0], &fast["runs"][0]);
    assert_eq!(rf["stop"], "extrapolated", "{}", rf["stop"]);
    assert_eq!(rf["extrapolated"], true);
    let (nl, nf) = (rl["steps"].as_u64().unwrap(), rf["steps"].as_u64().unwrap());
    assert!(nf * 2 < nl, "extrapolated run took {nf} steps against {nl}");
    for (a, b) in long["objective"]["per_f"].as_array().unwrap().iter().zip(fast["objective"]["per_f"].as_array().unwrap()) {
        for k in ["total_db", "accepted_db", "s11_db"] {
            let (x, y) = (a[k].as_f64().unwrap(), b[k].as_f64().unwrap());
            eprintln!("{k} at {}: long {x:.4} vs extrapolated {y:.4} (steps {nl} vs {nf})", a["f"]);
            assert!((x - y).abs() < 0.03, "{k} at {}: long {x:.4} vs extrapolated {y:.4} (steps {nl} vs {nf})", a["f"]);
        }
    }
}

#[test]
#[ignore]
fn extrapolation_experiment() {
    let r = run_ringing(json!({}), |s| s["objective"]["extrapolate"] = json!({"stop": false}), "exp", true);
    let run = &r["runs"][0];
    println!("cells {} steps {} stop {} comps {}", r["cells"], run["steps"], run["stop"], run["objective_components"]);
    let raw = run["objective_history"].as_array().unwrap();
    let ex = run["objective_extrap_history"].as_array().unwrap();
    for h in raw.iter() {
        let s = h[0].as_u64().unwrap();
        let x = ex.iter().find(|e| e[0].as_u64() == Some(s)).map(|e| e[1].to_string()).unwrap_or_default();
        println!("{s:6} raw {}  extrap {x}", h[1]);
    }
}
