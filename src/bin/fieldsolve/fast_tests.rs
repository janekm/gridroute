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
