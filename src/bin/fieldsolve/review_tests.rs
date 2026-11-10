use crate::*;
use geom::{BBox, P2, Prim, Shape};
use post::Cx;
use serde_json::json;

fn rect(x0: f64, y0: f64, x1: f64, y1: f64, net: &str) -> Prim {
    Prim { shape: Shape::Polygon(vec![P2::new(x0,y0), P2::new(x1,y0), P2::new(x1,y1), P2::new(x0,y1)]), dark: true, net: Some(net.into()), pad: None }
}

fn ctx(top: Vec<Prim>, spec: serde_json::Value) -> sim::Ctx {
    let plane=rect(-50.,-50.,50.,50.,"GND");
    let layers=vec![gerber::CopperLayer{name:"F.Cu".into(),index:1,prims:top,bbox:plane.shape.bbox()},gerber::CopperLayer{name:"B.Cu".into(),index:2,prims:vec![plane],bbox:BBox{x0:-50.,y0:-50.,x1:50.,y1:50.}}];
    let stack=vec![gerber::StackLayer{name:"F.Cu".into(),kind:gerber::StackKind::Copper,thickness_mm:0.035,eps_r:1.,loss_tan:0.},gerber::StackLayer{name:"FR4".into(),kind:gerber::StackKind::Dielectric,thickness_mm:0.2,eps_r:4.3,loss_tan:0.02},gerber::StackLayer{name:"B.Cu".into(),kind:gerber::StackKind::Copper,thickness_mm:0.035,eps_r:1.,loss_tan:0.}];
    sim::Ctx{spec,board:gerber::Board{layers,stack,holes:vec![],outline:vec![]},kicad:false,zs:vec![]}
}

#[test]
fn review_rotated_port_must_reach_negative_pad() {
    let a=rect(-0.2,-0.2,0.2,0.2,"A");
    let b=rect(0.8,0.8,1.2,1.2,"B");
    let mut cx=ctx(vec![a,b],json!({"roi":[-2,-2,2,2],"cell":0.1,"air":1,"air_z":1,"max_cell":0.5,"pml":2,"ports":[{"name":"P","pos":[0,0],"neg":[1,1],"layer":"F.Cu","width":0.4}]}));
    let st=sim::setup(&mut cx).unwrap();
    let sheet=st.layout.sheet("F.Cu").unwrap();
    let bridges:Vec<_>=sheet.prims.iter().filter(|p|p.net.as_deref()==Some("~P")).collect();
    let touches=bridges.iter().any(|p|p.shape.bbox().overlaps(&BBox{x0:0.8,y0:0.8,x1:1.2,y1:1.2}));
    assert!(touches,"bridge bounds {:?} miss the negative pad at (1,1)",bridges.iter().map(|p|p.shape.bbox()).collect::<Vec<_>>());
}

#[test]
fn review_macro_hole_must_preserve_underlying_copper() {
    let g="%FSLAX26Y26*%\n%MOMM*%\n%AMDONUT*1,1,2,0,0*1,0,1,0,0*%\n%ADD10R,4X4*%\n%ADD11DONUT*%\nD10*X0Y0D03*D11*X0Y0D03*M02*";
    let ps=gerber::parse_gerber(g).unwrap();
    let mut copper=false;
    for p in ps { if p.shape.contains(P2::new(0.,0.)) {copper=p.dark;} }
    assert!(copper,"transparent macro hole erased the pre-existing copper plane");
}

#[test]
fn review_cropped_plane_must_reach_absorber() {
    let mut cx=ctx(vec![rect(-50.,-50.,50.,50.,"TOP")],json!({"roi":[-1,-1,1,1],"cell":0.1,"air":1,"air_z":1,"max_cell":1,"pml":8}));
    let st=sim::setup(&mut cx).unwrap();
    let (masks,_)=model::sheet_masks(&st.grid,&st.layout);
    let g=&st.grid;
    let j=mesh::snap(&g.y,0.);
    let i=st.opts.pml[0];
    assert!(masks[0].ex[j*g.nx()+i],"copper window {:?}, PML inner interface x={} mm, first interior edge is air",st.layout.copper_window,g.x[i]*1e3);
}

fn cosim_fixture(tag: &str) -> std::path::PathBuf {
    let root=std::env::temp_dir().join(format!("fieldsolve-review-tests-{}", std::process::id())).join(tag);
    std::fs::create_dir_all(&root).unwrap();
    std::fs::write(root.join("em.s1p"),"# HZ S RI R 50\n1e9 0 0\n2e9 0 0\n").unwrap();
    std::fs::write(root.join("result.json"),r#"{"port_names":["P"]}"#).unwrap();
    std::fs::write(root.join("cfg.json"),json!({"em":"em.s1p","out":"out","keep":["P"],"variants":{"load":{}},"zref":{"P":[100,30]},"report":[1e9]}).to_string()).unwrap();
    cosim::run(&root.join("cfg.json")).unwrap();
    root
}

#[test]
fn review_complex_reference_must_preserve_reported_impedance() {
    let root=cosim_fixture("complex-zin");
    let r:serde_json::Value=serde_json::from_str(&std::fs::read_to_string(root.join("out/cosim.json")).unwrap()).unwrap();
    let z=&r["variants"]["load"]["report"][0]["zin"];
    assert!((z[0].as_f64().unwrap()-50.).abs()<1e-9 && z[1].as_f64().unwrap().abs()<1e-9,"a 50-ohm resistor is reported as {z}");
}

#[test]
fn review_renormalised_export_must_match_declared_reference() {
    let root=cosim_fixture("complex-export");
    let ts=cosim::read_touchstone(&root.join("out/load.s1p")).unwrap();
    let z=cosim::zin(ts.s[0][0][0],ts.z0);
    assert!((z-Cx::new(50.,0.)).abs()<1e-6,"exported 50-ohm resistor reads as {z:?}, header reference {}",ts.z0);
}

#[test]
fn review_touchstone_3port_rows_must_start_on_new_lines() {
    let root = std::env::temp_dir().join(format!("fieldsolve-review-tests-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let path = root.join("row-layout.s3p");
    let s=vec![vec![vec![Cx::new(0.1,0.)];3];3];
    post::write_touchstone(&path,&[1e9],&s,50.).unwrap();
    let txt=std::fs::read_to_string(path).unwrap();
    let lens:Vec<_>=txt.lines().filter(|l|!l.starts_with('#')&&!l.starts_with('!')).map(|l|l.split_whitespace().count()).collect();
    assert_eq!(lens,vec![7,6,6],"Touchstone v1 3-port matrix rows must be separate");
}

#[test]
fn review_farfield_metric_must_gate_energy_stop() {
    let g=mesh::Grid{x:vec![0.,1e-3,2e-3],y:vec![0.,1e-3,2e-3],z:vec![0.,1e-3,2e-3]};
    let p=model::build(&g,&model::Layout::default(),&[],&model::BuildOpts{pml:[0;6],..Default::default()});
    let calls=std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let cc=calls.clone();
    let spec=fdtd::RunSpec{max_steps:10,check_every:1,end_db:-70.,port_stop:true,gather_lag:5,gather_tol:0.01,gather_metric:Some(fdtd::GatherMetric(std::sync::Arc::new(move |_|vec![cc.fetch_add(1,std::sync::atomic::Ordering::SeqCst) as f64]))),..Default::default()};
    let r=fdtd::run_cpu(&p,&spec);
    assert_eq!(r.steps,10,"stopped before gather metric settled; metric evaluations {}",calls.load(std::sync::atomic::Ordering::SeqCst));
}

#[test]
fn review_tdr_open_must_settle_to_unit_reflection() {
    let fs=post::linspace(0.,10e9,1001);
    let s=vec![Cx::new(1.,0.);fs.len()];
    let (_,rho,_,_)=post::tdr(&fs,&s,50.,1e-9,2000);
    assert!((rho.last().unwrap()-1.).abs()<0.02,"open-circuit reflection plateau {} rather than 1",rho.last().unwrap());
}

// ---- acceptance cases added with the fixes (review of 2097498) ----

fn cosim_fixture_ref(tag: &str, zref: Option<[f64; 2]>) -> std::path::PathBuf {
    let root = std::env::temp_dir().join(format!("fieldsolve-review-tests-{}", std::process::id())).join(tag);
    std::fs::create_dir_all(&root).unwrap();
    std::fs::write(root.join("em.s1p"), "# HZ S RI R 50\n1e9 0 0\n2e9 0 0\n").unwrap();
    std::fs::write(root.join("result.json"), r#"{"port_names":["P"]}"#).unwrap();
    let mut cfg = json!({"em":"em.s1p","out":"out","keep":["P"],"variants":{"load":{}},"report":[1e9]});
    if let Some(z) = zref {
        cfg["zref"] = json!({"P": z});
    }
    std::fs::write(root.join("cfg.json"), cfg.to_string()).unwrap();
    cosim::run(&root.join("cfg.json")).unwrap();
    root
}

#[test]
fn cosim_unchanged_and_changed_real_reference_recover_the_physical_load() {
    for (tag, zr) in [("ref-none", None), ("ref-50", Some([50.0, 0.0])), ("ref-100", Some([100.0, 0.0]))] {
        let root = cosim_fixture_ref(tag, zr);
        let r: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(root.join("out/cosim.json")).unwrap()).unwrap();
        let z = &r["variants"]["load"]["report"][0]["zin"];
        assert!((z[0].as_f64().unwrap() - 50.0).abs() < 1e-9 && z[1].as_f64().unwrap().abs() < 1e-9, "{tag}: reported {z}");
        let ts = cosim::read_touchstone(&root.join("out/load.s1p")).unwrap();
        let ze = cosim::zin(ts.s[0][0][0], ts.z0);
        assert!((ze - Cx::new(50.0, 0.0)).abs() < 1e-6, "{tag}: export reads {ze:?} at its header reference {}", ts.z0);
        // a single real reference is declared and the data renormalised to it
        if let Some(z) = zr {
            assert!((ts.z0 - z[0]).abs() < 1e-9, "{tag}: header reference {}", ts.z0);
        }
    }
}

#[test]
fn tdr_short_matched_and_delayed_open() {
    let fs = post::linspace(0., 10e9, 1001);
    let end = |s: Vec<Cx>| *post::tdr(&fs, &s, 50., 1e-9, 2000).1.last().unwrap();
    assert!((end(vec![Cx::new(-1., 0.); fs.len()]) + 1.).abs() < 0.02, "short");
    assert!(end(vec![Cx::default(); fs.len()]).abs() < 0.02, "matched");
    // an open 0.3 ns away: nothing before the reflection arrives, unity after
    let tau = 0.3e-9;
    let s: Vec<Cx> = fs.iter().map(|&f| { let (sn, cs) = (-std::f64::consts::TAU * f * tau).sin_cos(); Cx::new(cs, sn) }).collect();
    let (t, rho, _, _) = post::tdr(&fs, &s, 50., 1e-9, 2000);
    let at = |x: f64| rho[t.iter().position(|&v| v >= x).unwrap()];
    assert!(at(0.15e-9).abs() < 0.02, "before the delayed open: {}", at(0.15e-9));
    assert!((at(0.6e-9) - 1.).abs() < 0.02, "after the delayed open: {}", at(0.6e-9));
}

#[test]
fn touchstone_4_and_5_port_rows() {
    let root = std::env::temp_dir().join(format!("fieldsolve-review-tests-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    for (n, want) in [(4usize, vec![9, 8, 8, 8]), (5, vec![9, 2, 8, 2, 8, 2, 8, 2, 8, 2])] {
        let path = root.join(format!("rows.s{n}p"));
        let s = vec![vec![vec![Cx::new(0.1, 0.)]; n]; n];
        post::write_touchstone(&path, &[1e9], &s, 50.).unwrap();
        let txt = std::fs::read_to_string(&path).unwrap();
        let lens: Vec<_> = txt.lines().filter(|l| !l.starts_with('#') && !l.starts_with('!')).map(|l| l.split_whitespace().count()).collect();
        assert_eq!(lens, want, "{n}-port rows");
        // and the reader still gets the matrix back
        let ts = cosim::read_touchstone(&path).unwrap();
        assert_eq!(ts.s.len(), n);
    }
}

#[test]
fn clear_macro_flash_keeps_copper_in_its_hole() {
    // a clear donut over a solid rectangle: the ring is removed, the hole keeps the rectangle's copper
    let g = "%FSLAX26Y26*%\n%MOMM*%\n%AMDONUT*1,1,2,0,0*1,0,1,0,0*%\n%ADD10R,4X4*%\n%ADD11DONUT*%\nD10*X0Y0D03*%LPC*%D11*X0Y0D03*M02*";
    let ps = gerber::parse_gerber(g).unwrap();
    let at = |q: P2| ps.iter().fold(false, |c, p| if p.shape.contains(q) { p.dark } else { c });
    assert!(at(P2::new(0., 0.)), "hole centre must stay copper");
    assert!(!at(P2::new(0.75, 0.)), "ring must be cleared");
    assert!(at(P2::new(1.5, 0.)), "outside the donut the rectangle stays");
}

#[test]
fn cropped_plane_reaches_the_absorber_on_all_four_sides() {
    let mut cx = ctx(vec![rect(-50., -50., 50., 50., "TOP")], json!({"roi":[-1,-1,1,1],"cell":0.1,"air":1,"air_z":1,"max_cell":1,"pml":8}));
    let st = sim::setup(&mut cx).unwrap();
    let (masks, _) = model::sheet_masks(&st.grid, &st.layout);
    let g = &st.grid;
    let (nx, ny) = (g.nx(), g.ny());
    let (jm, im) = (mesh::snap(&g.y, 0.), mesh::snap(&g.x, 0.));
    let pml = st.opts.pml;
    // x-edges at the -x / +x interfaces, y-edges at the -y / +y interfaces, and in the absorber itself
    for i in [pml[0], nx - 2 - pml[1], 1, nx - 2] {
        assert!(masks[0].ex[jm * nx + i], "x edge {i} on the centre row is not copper");
    }
    for j in [pml[2], ny - 2 - pml[3], 1, ny - 2] {
        assert!(masks[0].ey[j * nx + im], "y edge {j} on the centre column is not copper");
    }
}

#[test]
fn settled_far_field_metric_still_allows_an_energy_stop() {
    let g = mesh::Grid { x: vec![0., 1e-3, 2e-3], y: vec![0., 1e-3, 2e-3], z: vec![0., 1e-3, 2e-3] };
    let p = model::build(&g, &model::Layout::default(), &[], &model::BuildOpts { pml: [0; 6], ..Default::default() });
    let spec = fdtd::RunSpec { max_steps: 50, check_every: 1, end_db: -70., port_stop: true, gather_lag: 2, gather_tol: 0.01, gather_metric: Some(fdtd::GatherMetric(std::sync::Arc::new(|_| vec![1.0]))), ..Default::default() };
    let r = fdtd::run_cpu(&p, &spec);
    assert!(r.steps < 50, "a constant metric settles, so the (empty) run may stop early; ran {}", r.steps);
}

fn rotated_port_bridges(a: (f64, f64), b: (f64, f64), cell: f64) -> (Vec<geom::BBox>, sim::Setup) {
    let pa = rect(a.0 - 0.2, a.1 - 0.2, a.0 + 0.2, a.1 + 0.2, "A");
    let pb = rect(b.0 - 0.2, b.1 - 0.2, b.0 + 0.2, b.1 + 0.2, "B");
    let mut cx = ctx(vec![pa, pb], json!({"roi":[-2,-2,2,2],"cell":cell,"air":1,"air_z":1,"max_cell":0.5,"pml":2,"ports":[{"name":"P","pos":[a.0,a.1],"neg":[b.0,b.1],"layer":"F.Cu","width":0.4}]}));
    let st = sim::setup(&mut cx).unwrap();
    let bb = st.layout.sheet("F.Cu").unwrap().prims.iter().filter(|p| p.net.as_deref() == Some("~P")).map(|p| p.shape.bbox()).collect();
    (bb, st)
}

#[test]
fn rotated_ports_reach_both_pads_for_every_orientation_and_mesh() {
    // pos -> neg in all four diagonal directions, x- and y-dominant offsets, and a finer mesh
    let cases = [((0., 0.), (1., 1.)), ((1., 1.), (0., 0.)), ((0., 1.), (1., 0.)), ((1., 0.), (0., 1.)), ((0., 0.), (1., 0.4)), ((0., 0.), (0.4, 1.)), ((0., 0.), (-1., -0.6))];
    for cell in [0.1, 0.05] {
        for &(a, b) in &cases {
            let (bb, st) = rotated_port_bridges(a, b, cell);
            let near = |p: (f64, f64)| bb.iter().any(|x| x.overlaps(&geom::BBox { x0: p.0 - 0.2, y0: p.1 - 0.2, x1: p.0 + 0.2, y1: p.1 + 0.2 }));
            assert!(near(a) && near(b), "cell {cell} {a:?} -> {b:?}: bridges {bb:?}");
            // exactly one element edge set, its gap between the pads, polarity following pos
            let l = st.lumped.iter().find(|l| l.name == "P").unwrap();
            let (ax, o) = (l.axis, 1 - l.axis);
            let pc = |p: (f64, f64), k: usize| if k == 0 { p.0 } else { p.1 };
            let (lo, hi) = (pc(a, ax).min(pc(b, ax)), pc(a, ax).max(pc(b, ax)));
            assert!(l.lo[ax] >= lo - 1e-9 && l.hi[ax] <= hi + 1e-9, "cell {cell} {a:?} -> {b:?}: gap {:?}..{:?} outside the pads", l.lo, l.hi);
            assert!(l.hi[o] - l.lo[o] < cell, "gap should be one edge wide across");
            assert_eq!(l.plus_at_hi, pc(a, ax) > pc(b, ax), "cell {cell} {a:?} -> {b:?}: polarity");
        }
    }
}

#[test]
fn nearly_aligned_pads_keep_the_straight_bridge() {
    // an offset below 1.5 cells is treated as aligned: one straight bar each side, as before
    let (bb, st) = rotated_port_bridges((0., 0.), (1., 0.1), 0.1);
    assert_eq!(bb.len(), 2, "straight bridge expected, got {bb:?}");
    assert!(st.lumped.iter().any(|l| l.name == "P" && l.axis == 0));
}

/// The far-field box and the field snapshots are sampled at 3 fmax, so a frequency above fmax aliases.
/// A run asking for one is refused, with the frequency and the fmax to set in the message.
#[test]
fn frequencies_above_fmax_are_rejected() {
    let ok = json!({"fmax": 8e9, "farfield": {"freqs": [2.44e9, 4.88e9, 7.32e9]}, "fields": {"freqs": [2.44e9, 8e9]},
                    "objective": {"band": [2.44e9], "harmonics": [4.88e9, 7.32e9]}});
    assert!(sim::RunOpts::from_spec(&ok).check_freqs(sim::spec_fmax(&ok)).is_ok());
    // the default fmax is 3 GHz: the second harmonic of 2.44 GHz is out of reach
    for (bad, key) in [
        (json!({"farfield": {"freqs": [2.44e9, 4.88e9]}}), "farfield.freqs"),
        (json!({"fields": {"freqs": [4.88e9]}}), "fields.freqs"),
        (json!({"objective": {"band": [2.44e9, 4.88e9]}}), "objective.band"),
        (json!({"objective": {"band": [2.44e9], "harmonics": [4.88e9]}}), "objective.harmonics"),
    ] {
        let e = sim::RunOpts::from_spec(&bad).check_freqs(sim::spec_fmax(&bad)).unwrap_err();
        assert!(e.contains(key) && e.contains("4.880 GHz") && e.contains("3.000 GHz") && e.contains("4.880e9"), "{e}");
        assert!(!e.contains("2.440"), "only the offending frequency is named: {e}");
    }
    // a full run fails before the set-up and before any solver step
    let never: &dyn Fn(&model::Problem, &fdtd::RunSpec) -> fdtd::RunResult = &|_, _| panic!("the solver must not start");
    let spec = json!({"roi": [-5, -5, 5, 5], "cell": 0.5, "ports": [{"name": "P", "pos": "A.1", "neg": "B.1", "layer": "F.Cu"}],
                      "farfield": {"freqs": [2.44e9, 4.88e9]}});
    let dir = std::env::temp_dir().join(format!("fieldsolve-test-fmax-guard-{}", std::process::id()));
    let e = sim::run(&mut ctx(vec![rect(-1., -0.3, -0.1, 0.3, "A"), rect(0.1, -0.3, 1., 0.3, "B")], spec), &dir, never).unwrap_err();
    assert!(e.contains("farfield.freqs") && e.contains("fmax"), "{e}");
    assert!(!dir.exists(), "nothing is written for a refused run");
}

/// Surface-impedance copper on a domain without absorber: the loss must reach the first row of edges.
/// Ex at i = 0 and Ey at j = 0 are interior edges; only Ex on the y = 0 wall, Ey on the x = 0 wall and
/// edges on the z = 0 wall are left out, because the update reads H one cell back across the edge.
#[test]
fn sheet_loss_covers_the_first_row_of_edges() {
    assert!(model::cond_stencil_ok(0, 0, 1, 1) && model::cond_stencil_ok(1, 1, 0, 1));
    assert!(!model::cond_stencil_ok(0, 3, 0, 1) && !model::cond_stencil_ok(1, 0, 3, 1));
    assert!(!model::cond_stencil_ok(0, 3, 3, 0) && !model::cond_stencil_ok(1, 3, 3, 0));

    let lin = |n: usize, h: f64| (0..=n).map(|i| i as f64 * h * 1e-3).collect::<Vec<f64>>();
    let g = mesh::Grid { x: lin(6, 1.0), y: lin(5, 1.0), z: lin(4, 0.25) };
    let plate = |z: f64| model::Sheet {
        name: format!("plate{z}"),
        z,
        prims: vec![Prim { shape: Shape::Polygon(vec![P2::new(-1., -1.), P2::new(7., -1.), P2::new(7., 6.), P2::new(-1., 6.)]), dark: true, net: None, pad: None }],
    };
    let lay = model::Layout { sheets: vec![plate(0.75), plate(0.25)], ..Default::default() };
    let metal = cond::Metal { sigma: 5.8e7, t: 35e-6, rough: 0.0 };
    let conductor = Some(cond::Conductor { metals: vec![metal.clone(), metal], f_lo: 0.3e9, f_hi: 10e9 });
    let p = model::build(&g, &lay, &[], &model::BuildOpts { pml: [0; 6], fmax: 10e9, conductor, ..Default::default() });
    let (nx, ny) = (p.nx, p.ny);
    let edges: Vec<(usize, usize, usize, usize)> = p.cond.as_ref().unwrap().ids.iter().map(|&w| {
        let (a, id) = ((w >> 30) as usize, (w & 0x3fff_ffff) as usize);
        (a, id % nx, (id / nx) % ny, id / (nx * ny))
    }).collect();
    let count = |f: &dyn Fn(&(usize, usize, usize, usize)) -> bool| edges.iter().filter(|e| f(e)).count();
    assert!(count(&|e| e.3 == 1) > 0 && count(&|e| e.3 == 3) > 0 && count(&|e| e.3 != 1 && e.3 != 3) == 0, "both plates, nothing else: {edges:?}");
    // the first row carries as many lossy edges as the second
    assert!(count(&|e| e.0 == 0 && e.1 == 0) > 0, "no Ex edge at i = 0");
    assert_eq!(count(&|e| e.0 == 0 && e.1 == 0), count(&|e| e.0 == 0 && e.1 == 1));
    assert!(count(&|e| e.0 == 1 && e.2 == 0) > 0, "no Ey edge at j = 0");
    assert_eq!(count(&|e| e.0 == 1 && e.2 == 0), count(&|e| e.0 == 1 && e.2 == 1));
    // and the stencil never leaves the grid
    assert_eq!(count(&|e| (e.0 == 0 && e.2 == 0) || (e.0 == 1 && e.1 == 0)), 0);
}

/// Physics review A8: where the excitations are not independent S has no value; the zeros written there
/// must be reported, and the test has to be relative (a weak excitation is not a singular one).
#[test]
fn singular_wave_matrix_is_reported() {
    let z0 = 50.0f64;
    let s_ref = [[Cx::new(0.2, 0.1), Cx::new(0.5, -0.3)], [Cx::new(0.5, -0.3), Cx::new(-0.1, 0.4)]];
    // excitation j with unit incident wave at port j, the other port matched: a = e_j, b = S e_j
    let vi = |j: usize, k: usize, amp: f64| {
        let a = if j == k { Cx::new(1.0, 0.0) } else { Cx::default() };
        let b = s_ref[k][j];
        ((a + b) * (z0.sqrt() * amp), (a - b) * (amp / z0.sqrt()))
    };
    // frequency 0: as is; 1: the second excitation repeats the first; 2: both excitations 1e-9 as strong
    let cell = |j: usize, k: usize, f: usize| match f {
        0 => vi(j, k, 1.0),
        1 => vi(0, k, 1.0),
        _ => vi(j, k, 1e-9),
    };
    let v: Vec<Vec<Vec<Cx>>> = (0..2).map(|j| (0..2).map(|k| (0..3).map(|f| cell(j, k, f).0).collect()).collect()).collect();
    let i: Vec<Vec<Vec<Cx>>> = (0..2).map(|j| (0..2).map(|k| (0..3).map(|f| cell(j, k, f).1).collect()).collect()).collect();
    let (s, singular) = post::s_matrix(&v, &i, z0);
    assert_eq!(singular, vec![1]);
    for p in 0..2 {
        for j in 0..2 {
            assert!((s[p][j][0] - s_ref[p][j]).abs() < 1e-12 && (s[p][j][2] - s_ref[p][j]).abs() < 1e-9, "S{p}{j}: {:?} {:?}", s[p][j][0], s[p][j][2]);
            assert_eq!(s[p][j][1].abs(), 0.0);
        }
    }
}

/// Physics review A8: a conductor one cell across that leaves the far-field box along a face normal sets
/// no tangential edge on the face. The box check has to see it, and the warning has to exist as text (it
/// goes to result.json "checks", not only to the log).
#[test]
fn far_field_box_check_sees_a_wire_along_the_face_normal() {
    let mut st = cases::dipole_setup(10.0, 1.0, 3e9, 10.0);
    let p = model::build(&st.grid, &st.layout, &st.lumped, &st.opts);
    let plan = nf2ff::auto_plan(&p, st.opts.pml, 3).unwrap();
    assert_eq!(nf2ff::non_free_edges(&p, &plan), (0, 0));
    assert!(nf2ff::box_warning(&p, &plan).is_none());
    // the upper arm continued to the top of the domain: it crosses the +z face
    st.layout.barrels[1].z1 = st.grid.z.last().unwrap() * 1e3;
    let p = model::build(&st.grid, &st.layout, &st.lumped, &st.opts);
    let (tangential, normal) = nf2ff::non_free_edges(&p, &plan);
    assert_eq!(tangential, 0, "the wire sets no tangential edge: the old check could not see it");
    assert_eq!(normal, 2, "the edge ending on the face and the one starting on it");
    let w = nf2ff::box_warning(&p, &plan).unwrap();
    assert!(w.contains("2 normal") && w.contains("not valid"), "{w}");
    // a box on the absorber is refused
    assert!(nf2ff::auto_plan(&p, st.opts.pml, 0).unwrap_err().contains("gap"));
}

/// Physics review A8: with a far field the objective is the radiated power, which exists only at the
/// far-field frequencies. A band frequency outside that list used to fall back to the accepted power
/// silently; a harmonic outside it shifted the names of the tracked components.
#[test]
fn objective_band_must_be_in_the_far_field_list() {
    let spec = |band: serde_json::Value, far: serde_json::Value| json!({"fmax": 6e9, "farfield": far, "objective": {"band": band, "harmonics": [4.88e9]}});
    let check = |s: &serde_json::Value| sim::RunOpts::from_spec(s).check_freqs(sim::spec_fmax(s));
    assert!(check(&spec(json!([2.4e9, 2.45e9]), json!({"freqs": [2.4e9, 2.45e9, 2.5e9]}))).is_ok());
    // no far field: the objective is the accepted power, at any frequency of the sweep
    assert!(check(&spec(json!([2.4e9, 2.47e9]), serde_json::Value::Null)).is_ok());
    let e = check(&spec(json!([2.4e9, 2.47e9]), json!({"freqs": [2.4e9, 2.45e9, 2.5e9]}))).unwrap_err();
    assert!(e.contains("objective.band") && e.contains("2.4700") && e.contains("farfield.freqs") && !e.contains("2.4000"), "{e}");
    // the names of the tracked components: a radiated entry only where the box is transformed
    let o = sim::ObjectiveOpts::from_spec(&spec(json!([2.4e9, 2.45e9]), serde_json::Value::Null)).unwrap();
    let names = vec!["P".to_string(), "Q".to_string()];
    let c = o.components(&names, 0, Some(&[2.4e9, 2.45e9, 2.5e9]));
    assert_eq!(c.len(), 3 * 2 + 2, "{c:?}");
    assert_eq!(c.iter().filter(|n| n.starts_with("radiated")).count(), 2);
    assert!(!c.iter().any(|n| n == "radiated 4.880GHz") && c.iter().any(|n| n == "radiated 2.450GHz"), "{c:?}");
    assert_eq!(o.components(&names, 0, None).len(), 3 * 2);
}

/// Physics review A8: settings that contradict each other are refused at set-up, and the wideband fit's
/// error is kept.
#[test]
fn contradictory_loss_models_and_a_zero_gap_are_refused() {
    let base = json!({"roi": [-2, -2, 2, 2], "cell": 0.2, "air": 1, "air_z": 1, "max_cell": 0.5, "pml": 2, "fmax": 6e9});
    let with = |extra: serde_json::Value| {
        let mut s = base.clone();
        for (k, v) in extra.as_object().unwrap() {
            s[k] = v.clone();
        }
        sim::setup(&mut ctx(vec![rect(-0.2, -0.2, 0.2, 0.2, "A")], s))
    };
    let e = with(json!({"loss": {"model": "wideband", "copper": true}, "conductor": true})).err().unwrap();
    assert!(e.contains("counted twice"), "{e}");
    assert!(with(json!({"loss": {"model": "wideband", "copper": true}, "conductor": false})).is_ok());
    let e = with(json!({"farfield": {"freqs": [2.4e9], "gap": 0}})).err().unwrap();
    assert!(e.contains("gap") && e.contains("at least 1"), "{e}");
    assert!(with(json!({"farfield": {"freqs": [2.4e9], "gap": 1}})).is_ok());

    // the fit error: a few percent for a board band, none without the wideband model
    assert!(with(json!({})).unwrap().loss_fit.is_none());
    let st = with(json!({"loss": {"model": "wideband"}, "conductor": true})).unwrap();
    let (poles, worst) = st.loss_fit.unwrap();
    assert!(poles >= 3 && worst > 0.0 && worst < sim::LOSS_FIT_WARN, "{poles} poles, {worst}");
    assert!(!st.checks.iter().any(|c| c.contains("Debye fit")), "{:?}", st.checks);
    // eight poles (the most there are) over nine decades cannot hold the loss tangent flat: say so
    let st = with(json!({"loss": {"model": "wideband", "f_lo": 1e3, "f_hi": 1e12}})).unwrap();
    assert!(st.loss_fit.unwrap().1 > sim::LOSS_FIT_WARN);
    assert!(st.checks.iter().any(|c| c.contains("Debye fit")), "{:?}", st.checks);
}
