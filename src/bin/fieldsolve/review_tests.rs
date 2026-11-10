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
                    "objective": {"band": [2.40e9, 2.48e9], "harmonics": [4.88e9, 7.32e9]}});
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
