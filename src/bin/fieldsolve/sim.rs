//! Board simulations from a JSON spec: Gerber board -> layered layout -> mesh -> ports/elements -> runs ->
//! impedance, S-parameters, loop inductance, field snapshots, and a scene file for the openEMS cross-check.
//!
//! Spec (all lengths mm; `coords: "kicad"` (default) takes x/y as KiCad shows them, Gerber y is negated):
//! ```json
//! { "gerbers": "dir", "out": "dir", "coords": "kicad",
//!   "roi": [x0, y0, x1, y1], "cell": 0.1, "cell_z": 0.077, "min_cells_layer": 2,
//!   "air": 3.0, "air_z": 4.0, "pml": 8, "max_cell": 2.0,
//!   "fmin": 1e6, "fmax": 3e9, "nf": 300, "f_loss": 1e9, "end_db": -70, "max_steps": 400000,
//!   "solder_mask": false, "backend": "auto",
//!   "ports":    [ {"name": "VDD", "pos": "U3.12", "neg": "U3.13", "layer": "F.Cu", "r": 1.0},
//!                 {"name": "ANT", "at": "J1.1", "from": "F.Cu", "to": "In1.Cu", "r": 50} ],
//!   "elements": [ {"kind": "short", "pos": "C5.1", "neg": "C5.2", "layer": "F.Cu"},
//!                 {"kind": "R", "r": 50, "at": [10, 12], "from": "F.Cu", "to": "In1.Cu"} ],
//!   "excite": "all" | ["VDD"], "l_band": [10e6, 100e6],
//!   "fields": {"freqs": [2.4e9]} }
//! ```
//! A pad reference is `REF.PIN` (from the Gerber X2 `.P` attributes); a point is `[x, y]`.

use crate::fdtd::{gauss_excite, RunResult, RunSpec};
use crate::geom::{BBox, Prim, Shape, P2};
use crate::gerber::{Board, StackKind};
use crate::mesh::{lines_1d, Grid, Refine};
use crate::model::{build, BuildOpts, Barrel, Layout, Lumped, Problem, Sheet, Slab, C0};
use crate::post::{fit_rl, linspace, port_spectra, s_column, write_touchstone, Cx};
use serde_json::{json, Value};
use std::path::{Path, PathBuf};

pub struct Ctx {
    pub spec: Value,
    pub board: Board,
    pub kicad: bool,
    /// copper layer name -> z (mm)
    pub zs: Vec<(String, f64)>,
}

fn num(v: &Value, k: &str, d: f64) -> f64 {
    v.get(k).and_then(|x| x.as_f64()).unwrap_or(d)
}

impl Ctx {
    pub fn pt(&self, x: f64, y: f64) -> P2 {
        if self.kicad {
            P2::new(x, -y)
        } else {
            P2::new(x, y)
        }
    }
    pub fn z_of(&self, layer: &str) -> Result<f64, String> {
        self.zs.iter().find(|(n, _)| n == layer).map(|x| x.1).ok_or_else(|| format!("unknown copper layer {layer}"))
    }
    /// A pad reference or point -> (centre, bbox) in Gerber coordinates, searching `layer` (or all layers).
    pub fn locate(&self, v: &Value, layer: Option<&str>) -> Result<(P2, BBox), String> {
        if let Some(a) = v.as_array() {
            let p = self.pt(a[0].as_f64().unwrap_or(0.0), a[1].as_f64().unwrap_or(0.0));
            let mut b = BBox::empty();
            b.add(p);
            return Ok((p, b));
        }
        let s = v.as_str().ok_or("a location is a pad reference \"REF.PIN\" or [x, y]")?;
        let (r, pin) = s.rsplit_once('.').ok_or_else(|| format!("pad reference {s} is not REF.PIN"))?;
        let mut bb = BBox::empty();
        for l in &self.board.layers {
            if layer.is_some_and(|n| n != l.name) {
                continue;
            }
            for p in &l.prims {
                if p.dark && p.pad.as_ref().is_some_and(|(a, b)| a == r && b == pin) {
                    bb.union(&p.shape.bbox());
                }
            }
        }
        if bb.is_empty() {
            return Err(format!("pad {s} not found{}", layer.map(|l| format!(" on {l}")).unwrap_or_default()));
        }
        Ok((P2::new(0.5 * (bb.x0 + bb.x1), 0.5 * (bb.y0 + bb.y1)), bb))
    }
}

/// Copper z positions and dielectric slabs from the stack (bottom copper at z = 0, copper zero-thickness).
pub fn stack_layout(board: &Board, mask: bool) -> (Vec<(String, f64)>, Vec<Slab>) {
    let mut zs = Vec::new();
    let mut slabs = Vec::new();
    let mut z = 0.0;
    let st: Vec<_> = board.stack.iter().rev().collect();
    let first_cu = st.iter().position(|s| s.kind == StackKind::Copper).unwrap_or(0);
    let last_cu = st.iter().rposition(|s| s.kind == StackKind::Copper).unwrap_or(st.len().saturating_sub(1));
    // solder mask below the bottom copper
    if mask {
        for s in st[..first_cu].iter().rev() {
            if s.kind == StackKind::SolderMask {
                slabs.push(Slab { z0: z - s.thickness_mm, z1: z, eps: s.eps_r, tand: s.loss_tan });
            }
        }
    }
    for s in &st[first_cu..=last_cu] {
        match s.kind {
            StackKind::Copper => zs.push((s.name.clone(), z)),
            StackKind::Dielectric => {
                slabs.push(Slab { z0: z, z1: z + s.thickness_mm, eps: s.eps_r, tand: s.loss_tan });
                z += s.thickness_mm;
            }
            _ => {}
        }
    }
    if mask {
        for s in &st[last_cu + 1..] {
            if s.kind == StackKind::SolderMask {
                slabs.push(Slab { z0: z, z1: z + s.thickness_mm, eps: s.eps_r, tand: s.loss_tan });
            }
        }
    }
    (zs, slabs)
}

pub struct Setup {
    pub grid: Grid,
    pub layout: Layout,
    pub lumped: Vec<Lumped>,
    pub opts: BuildOpts,
    pub freqs: Vec<f64>,
    pub fmax: f64,
    pub end_db: f64,
    pub max_steps: usize,
}

/// Build layout, mesh and lumped elements from the spec.
pub fn setup(cx: &mut Ctx) -> Result<Setup, String> {
    let sp = cx.spec.clone();
    let mask = sp.get("solder_mask").and_then(|v| v.as_bool()).unwrap_or(false);
    let (zs, slabs) = stack_layout(&cx.board, mask);
    cx.zs = zs.clone();
    let cell = num(&sp, "cell", 0.1);
    let fmax = num(&sp, "fmax", 3e9);
    let air = num(&sp, "air", 3.0);
    let air_z = num(&sp, "air_z", 4.0);
    let npml = num(&sp, "pml", 8.0) as usize;
    let hmax = num(&sp, "max_cell", (C0 / fmax / 20.0 * 1e3).min(2.0));
    let mut board_bb = BBox::empty();
    cx.board.outline.iter().flatten().for_each(|&p| board_bb.add(p));
    if board_bb.is_empty() {
        cx.board.layers.iter().for_each(|l| board_bb.union(&l.bbox));
    }
    let roi = match sp.get("roi").and_then(|v| v.as_array()) {
        Some(a) => {
            let (p, q) = (cx.pt(a[0].as_f64().unwrap(), a[1].as_f64().unwrap()), cx.pt(a[2].as_f64().unwrap(), a[3].as_f64().unwrap()));
            BBox { x0: p.x.min(q.x), y0: p.y.min(q.y), x1: p.x.max(q.x), y1: p.y.max(q.y) }
        }
        None => board_bb,
    };
    // the PML sits beyond the air margin; it is npml cells of hmax-ish size
    let dom = roi.grow(air);

    // ports and elements -> lumped boxes, bridges on sheets, fixed mesh lines
    let mut sheets: Vec<Sheet> = Vec::new();
    for (name, z) in &zs {
        let l = cx.board.layers.iter().find(|l| &l.name == name).ok_or_else(|| format!("no Gerber for {name}"))?;
        let grown = dom.grow(1.0);
        let prims = l.prims.iter().filter(|p| p.shape.bbox().overlaps(&grown)).cloned().collect();
        sheets.push(Sheet { name: name.clone(), z: *z, prims });
    }
    let mut fx: Vec<f64> = vec![roi.x0, roi.x1];
    let mut fy: Vec<f64> = vec![roi.y0, roi.y1];
    let mut fz: Vec<f64> = zs.iter().map(|x| x.1).collect();
    for s in &slabs {
        fz.push(s.z0);
        fz.push(s.z1);
    }
    struct Pending {
        name: String,
        el: Value,
        port: bool,
    }
    let mut pend = Vec::new();
    for p in sp.get("ports").and_then(|v| v.as_array()).cloned().unwrap_or_default() {
        pend.push(Pending { name: p.get("name").and_then(|v| v.as_str()).unwrap_or("P").to_string(), el: p, port: true });
    }
    for (k, p) in sp.get("elements").and_then(|v| v.as_array()).cloned().unwrap_or_default().into_iter().enumerate() {
        pend.push(Pending { name: format!("el{k}"), el: p, port: false });
    }
    // first pass: geometry of each element, so its lines can be fixed before meshing
    enum Geo {
        Vert { c: P2, half: [f64; 2], z0: f64, z1: f64, plus_hi: bool },
        Plane { a: P2, b: P2, axis: usize, w: f64, z: f64, layer: String },
    }
    let mut geos = Vec::new();
    for pd in &pend {
        let e = &pd.el;
        if let Some(at) = e.get("at") {
            let from = e.get("from").and_then(|v| v.as_str()).ok_or_else(|| format!("{}: vertical element needs from/to layers", pd.name))?;
            let to = e.get("to").and_then(|v| v.as_str()).ok_or_else(|| format!("{}: vertical element needs from/to layers", pd.name))?;
            let (c, _) = cx.locate(at, Some(from))?;
            let (zf, zt) = (cx.z_of(from)?, cx.z_of(to)?);
            let size = e.get("size").and_then(|v| v.as_array()).map(|a| [a[0].as_f64().unwrap_or(0.0), a[1].as_f64().unwrap_or(0.0)]).unwrap_or([0.0, 0.0]);
            fx.push(c.x);
            fy.push(c.y);
            if size[0] > 0.0 {
                fx.extend([c.x - size[0] / 2.0, c.x + size[0] / 2.0]);
            }
            if size[1] > 0.0 {
                fy.extend([c.y - size[1] / 2.0, c.y + size[1] / 2.0]);
            }
            geos.push(Geo::Vert { c, half: [size[0] / 2.0, size[1] / 2.0], z0: zf.min(zt), z1: zf.max(zt), plus_hi: zf > zt });
        } else {
            let layer = e.get("layer").and_then(|v| v.as_str()).unwrap_or("F.Cu").to_string();
            let (a, ba) = cx.locate(e.get("pos").ok_or_else(|| format!("{}: needs pos/neg or at", pd.name))?, Some(&layer))?;
            let (b, bb) = cx.locate(e.get("neg").ok_or_else(|| format!("{}: needs pos/neg or at", pd.name))?, Some(&layer))?;
            let axis = if (b.x - a.x).abs() >= (b.y - a.y).abs() { 0 } else { 1 };
            // bridge width: the smaller pad dimension across the axis
            let wa = if axis == 0 { ba.y1 - ba.y0 } else { ba.x1 - ba.x0 };
            let wb = if axis == 0 { bb.y1 - bb.y0 } else { bb.x1 - bb.x0 };
            let w = num(e, "width", wa.min(wb).max(cell));
            let z = cx.z_of(&layer)?;
            // a straight bridge on the line through pos along the axis
            let b2 = if axis == 0 { P2::new(b.x, a.y) } else { P2::new(a.x, b.y) };
            if axis == 0 {
                fy.extend([a.y - w / 2.0, a.y + w / 2.0]);
                fx.push(0.5 * (a.x + b.x));
            } else {
                fx.extend([a.x - w / 2.0, a.x + w / 2.0]);
                fy.push(0.5 * (a.y + b.y));
            }
            geos.push(Geo::Plane { a, b: b2, axis, w, z, layer });
        }
    }

    // mesh
    let board_t = zs.iter().map(|x| x.1).fold(0.0, f64::max);
    let zlo = slabs.iter().map(|s| s.z0).fold(0.0, f64::min);
    let zhi = slabs.iter().map(|s| s.z1).fold(board_t, f64::max);
    let cell_z = num(&sp, "cell_z", cell.min(slabs.iter().map(|s| s.z1 - s.z0).fold(f64::INFINITY, f64::min) / 2.0));
    let min_layer = num(&sp, "min_cells_layer", 2.0) as usize;
    let pml_w = npml as f64 * hmax.min(1.0).max(cell);
    let mut gx = lines_1d(dom.x0 - pml_w, dom.x1 + pml_w, &fx, &[Refine { a: roi.x0, b: roi.x1, h: cell }], 1.3, hmax, &[]);
    let mut gy = lines_1d(dom.y0 - pml_w, dom.y1 + pml_w, &fy, &[Refine { a: roi.y0, b: roi.y1, h: cell }], 1.3, hmax, &[]);
    let mc: Vec<(f64, f64, usize)> = slabs.iter().map(|s| (s.z0, s.z1, min_layer)).collect();
    let gz = lines_1d(zlo - air_z - pml_w, zhi + air_z + pml_w, &fz, &[Refine { a: zlo, b: zhi, h: cell_z }], 1.3, hmax, &mc);
    // keep only fixed lines that fall inside the domain
    gx.dedup();
    gy.dedup();
    let grid = Grid { x: gx.iter().map(|v| v * 1e-3).collect(), y: gy.iter().map(|v| v * 1e-3).collect(), z: gz.iter().map(|v| v * 1e-3).collect() };

    // second pass: lumped boxes on the mesh; bridges as copper
    let mut lumped = Vec::new();
    for (pd, g) in pend.iter().zip(&geos) {
        let e = &pd.el;
        let kind = e.get("kind").and_then(|v| v.as_str()).unwrap_or(if pd.port { "port" } else { "R" });
        let r = num(e, "r", if pd.port { 50.0 } else { 0.0 });
        let c = num(e, "c", 0.0);
        match g {
            Geo::Vert { c: ctr, half, z0, z1, plus_hi } => {
                lumped.push(Lumped {
                    name: pd.name.clone(),
                    axis: 2,
                    lo: [ctr.x - half[0], ctr.y - half[1], *z0],
                    hi: [ctr.x + half[0], ctr.y + half[1], *z1],
                    plus_at_hi: *plus_hi,
                    r: if kind == "short" { 0.0 } else { r },
                    c: if kind == "short" { 0.0 } else { c },
                    port: pd.port,
                });
            }
            Geo::Plane { a, b, axis, w, z, layer } => {
                let sheet = sheets.iter_mut().find(|s| &s.name == layer).unwrap();
                let (pa, pb) = if *axis == 0 { (a.x, b.x) } else { (a.y, b.y) };
                let oc = if *axis == 0 { a.y } else { a.x };
                // a straight copper bar along the axis between u0 and u1, w wide; bars end exactly on mesh
                // lines so the gap edges between them stay free
                let bar = |u0: f64, u1: f64| -> Prim {
                    let (lo, hi) = (u0.min(u1), u0.max(u1));
                    let pts = if *axis == 0 {
                        vec![P2::new(lo, oc - w / 2.0), P2::new(hi, oc - w / 2.0), P2::new(hi, oc + w / 2.0), P2::new(lo, oc + w / 2.0)]
                    } else {
                        vec![P2::new(oc - w / 2.0, lo), P2::new(oc + w / 2.0, lo), P2::new(oc + w / 2.0, hi), P2::new(oc - w / 2.0, hi)]
                    };
                    Prim { shape: Shape::Polygon(pts), dark: true, net: None, pad: None }
                };
                if kind == "short" && !pd.port {
                    sheet.prims.push(bar(pa, pb));
                    continue;
                }
                let lines = if *axis == 0 { &gx } else { &gy };
                let mid = 0.5 * (pa + pb);
                let mut i = crate::mesh::snap(lines, mid);
                if i + 1 >= lines.len() {
                    i -= 1;
                }
                let (g0, g1) = (lines[i], lines[i + 1]);
                sheet.prims.push(bar(pa.min(pb), g0));
                sheet.prims.push(bar(g1, pa.max(pb)));
                let mut lo = [0.0; 3];
                let mut hi = [0.0; 3];
                lo[*axis] = g0;
                hi[*axis] = g1;
                let o = 1 - *axis;
                lo[o] = oc - w / 2.0;
                hi[o] = oc + w / 2.0;
                lo[2] = *z;
                hi[2] = *z;
                lumped.push(Lumped { name: pd.name.clone(), axis: *axis, lo, hi, plus_at_hi: pa > pb, r: if kind == "short" { 0.0 } else { r }, c, port: pd.port });
            }
        }
    }

    // vias and pins within the domain
    let ncu = zs.len();
    let mut barrels = Vec::new();
    for h in &cx.board.holes {
        if !h.plated || !dom.grow(0.5).contains(h.c) {
            continue;
        }
        let (f, t) = (h.from.clamp(1, ncu) - 1, h.to.clamp(1, ncu) - 1);
        // stack layers are listed top to bottom; zs bottom to top
        let za = zs[ncu - 1 - f].1;
        let zb = zs[ncu - 1 - t].1;
        barrels.push(Barrel { c: h.c, r: h.d / 2.0, z0: za.min(zb), z1: za.max(zb) });
        if let Some(e) = h.end {
            barrels.push(Barrel { c: e, r: h.d / 2.0, z0: za.min(zb), z1: za.max(zb) });
        }
    }
    let layout = Layout { sheets, slabs, barrels, outline: cx.board.outline.clone() };
    let fmin = num(&sp, "fmin", fmax / 300.0);
    let nf = num(&sp, "nf", 300.0) as usize;
    Ok(Setup {
        grid,
        layout,
        lumped,
        opts: BuildOpts { pml: [npml; 6], f_loss: num(&sp, "f_loss", 1e9), fmax, cfl: 0.95 },
        freqs: linspace(fmin, fmax, nf),
        fmax,
        end_db: num(&sp, "end_db", -70.0),
        max_steps: num(&sp, "max_steps", 400_000.0) as usize,
    })
}

/// Scene file: exactly what was simulated, for an independent solver (tools/openems_scene.py).
pub fn scene_json(st: &Setup, exc_port: Option<usize>) -> Value {
    let poly = |p: &Prim| -> Vec<[f64; 2]> {
        let pts: Vec<P2> = match &p.shape {
            Shape::Polygon(v) => v.clone(),
            Shape::Circle { c, r } => (0..32).map(|k| {
                let t = k as f64 / 32.0 * std::f64::consts::TAU;
                P2::new(c.x + r * t.cos(), c.y + r * t.sin())
            }).collect(),
            Shape::Capsule { a, b, r } => {
                let ang = (b.y - a.y).atan2(b.x - a.x);
                let mut v = Vec::new();
                for k in 0..=12 {
                    let t = ang + std::f64::consts::FRAC_PI_2 + k as f64 / 12.0 * std::f64::consts::PI;
                    v.push(P2::new(b.x + r * t.cos(), b.y + r * t.sin()));
                }
                for k in 0..=12 {
                    let t = ang - std::f64::consts::FRAC_PI_2 + k as f64 / 12.0 * std::f64::consts::PI;
                    v.push(P2::new(a.x + r * t.cos(), a.y + r * t.sin()));
                }
                v
            }
        };
        pts.iter().map(|p| [p.x, p.y]).collect()
    };
    json!({
        "unit": 1e-3,
        "grid": {"x": st.grid.x.iter().map(|v| v * 1e3).collect::<Vec<_>>(), "y": st.grid.y.iter().map(|v| v * 1e3).collect::<Vec<_>>(), "z": st.grid.z.iter().map(|v| v * 1e3).collect::<Vec<_>>()},
        "sheets": st.layout.sheets.iter().map(|s| json!({"name": s.name, "z": s.z, "prims": s.prims.iter().map(|p| json!({"dark": p.dark, "pts": poly(p)})).collect::<Vec<_>>()})).collect::<Vec<_>>(),
        "slabs": st.layout.slabs.iter().map(|s| json!({"z0": s.z0, "z1": s.z1, "eps": s.eps, "tand": s.tand})).collect::<Vec<_>>(),
        "outline": st.layout.outline.iter().map(|l| l.iter().map(|p| [p.x, p.y]).collect::<Vec<_>>()).collect::<Vec<_>>(),
        "barrels": st.layout.barrels.iter().map(|b| json!({"c": [b.c.x, b.c.y], "r": b.r, "z0": b.z0, "z1": b.z1})).collect::<Vec<_>>(),
        "lumped": st.lumped.iter().map(|l| json!({"name": l.name, "axis": l.axis, "lo": l.lo, "hi": l.hi, "plus_at_hi": l.plus_at_hi, "r": l.r, "c": l.c, "port": l.port})).collect::<Vec<_>>(),
        "pml": st.opts.pml, "f_loss": st.opts.f_loss, "fmax": st.fmax, "end_db": st.end_db, "max_steps": st.max_steps,
        "freqs": st.freqs, "excite": exc_port,
    })
}

pub type Backend<'a> = &'a dyn Fn(&Problem, &RunSpec) -> RunResult;

/// Run every excitation, write results to `out`.
pub fn run(cx: &mut Ctx, out: &Path, backend: Backend) -> Result<Value, String> {
    let st = setup(cx)?;
    std::fs::create_dir_all(out).map_err(|e| e.to_string())?;
    let t0 = std::time::Instant::now();
    let p = build(&st.grid, &st.layout, &st.lumped, &st.opts);
    let np = p.ports.len();
    eprintln!(
        "fieldsolve: grid {} x {} x {} = {:.2} M cells, dt {:.3e} s, {} ports, build {:.1} s",
        p.nx, p.ny, p.nz, p.n() as f64 * 1e-6, p.dt, np, t0.elapsed().as_secs_f64()
    );
    if np == 0 {
        return Err("no ports".into());
    }
    let names: Vec<String> = p.ports.iter().map(|x| x.name.clone()).collect();
    let excite: Vec<usize> = match cx.spec.get("excite") {
        Some(Value::Array(a)) => a.iter().filter_map(|n| names.iter().position(|m| Some(m.as_str()) == n.as_str())).collect(),
        Some(Value::String(s)) if s == "all" => (0..np).collect(),
        _ => vec![0],
    };
    let fields_f: Vec<f64> = cx.spec.pointer("/fields/freqs").and_then(|v| v.as_array()).map(|a| a.iter().filter_map(|x| x.as_f64()).collect()).unwrap_or_default();
    let exc = gauss_excite(0.0, st.fmax, p.dt);
    let dft_every = ((1.0 / (2.0 * st.fmax * 1.5)) / p.dt).floor().max(1.0) as usize;
    std::fs::write(out.join("scene.json"), serde_json::to_string(&scene_json(&st, excite.first().copied())).unwrap()).map_err(|e| e.to_string())?;
    let freqs = &st.freqs;
    let z0 = 50.0;
    let mut s_mat: Vec<Vec<Vec<Cx>>> = vec![vec![vec![Cx::default(); freqs.len()]; np]; np];
    let mut zself: Vec<Option<Vec<Cx>>> = vec![None; np];
    let mut runs = Vec::new();
    for (q, &j) in excite.iter().enumerate() {
        let spec = RunSpec {
            excite: Some(j),
            exc: exc.clone(),
            max_steps: st.max_steps,
            end_db: st.end_db,
            check_every: 1000,
            dft_freqs: if q == 0 { fields_f.clone() } else { vec![] },
            dft_every,
        };
        let r = backend(&p, &spec);
        eprintln!("fieldsolve: port {} excited: {} steps in {:.1} s on {} ({:.0} Mcell/s)", names[j], r.steps, r.seconds, r.backend, p.n() as f64 * r.steps as f64 / r.seconds * 1e-6);
        let mut vs = Vec::new();
        let mut is = Vec::new();
        for k in 0..np {
            let src: &[f32] = if k == j { &exc } else { &[] };
            let (v, i) = port_spectra(&r.v[k], src, p.ports[k].r, p.dt, freqs);
            vs.push(v);
            is.push(i);
        }
        let col = s_column(&vs, &is, z0, j);
        for k in 0..np {
            s_mat[k][j] = col[k].clone();
        }
        zself[j] = Some(vs[j].iter().zip(&is[j]).map(|(&v, &i)| v / i).collect());
        if let Some(fd) = &r.fields {
            crate::export::write_fields(out, &p, &st, fd).map_err(|e| e.to_string())?;
        }
        runs.push(json!({"port": names[j], "steps": r.steps, "seconds": r.seconds, "backend": r.backend}));
    }
    if excite.len() == np {
        write_touchstone(&out.join(format!("sparams.s{np}p")), freqs, &s_mat, z0).map_err(|e| e.to_string())?;
    }
    let lband = cx.spec.get("l_band").and_then(|v| v.as_array()).map(|a| (a[0].as_f64().unwrap_or(0.0), a[1].as_f64().unwrap_or(0.0)));
    let mut ports_out = Vec::new();
    for j in 0..np {
        let Some(z) = &zself[j] else { continue };
        let mut o = json!({
            "name": names[j],
            "freqs": freqs,
            "z_re": z.iter().map(|c| c.re).collect::<Vec<_>>(),
            "z_im": z.iter().map(|c| c.im).collect::<Vec<_>>(),
            "s11_db": s_mat[j][j].iter().map(|c| 20.0 * c.abs().max(1e-30).log10()).collect::<Vec<_>>(),
        });
        let (f0, f1) = lband.unwrap_or((freqs[0], freqs[freqs.len() / 10].max(freqs[0])));
        let (r_fit, l_fit) = fit_rl(freqs, z, f0, f1);
        o["l_fit"] = json!({"band": [f0, f1], "r_ohm": r_fit, "l_h": l_fit});
        eprintln!("fieldsolve: {}: R = {:.4} ohm, L = {:.4} nH over {:.3e}..{:.3e} Hz", names[j], r_fit, l_fit * 1e9, f0, f1);
        ports_out.push(o);
    }
    let res = json!({
        "grid": [p.nx, p.ny, p.nz], "cells": p.n(), "dt": p.dt, "runs": runs, "ports": ports_out,
        "port_names": names,
    });
    std::fs::write(out.join("result.json"), serde_json::to_string_pretty(&res).unwrap()).map_err(|e| e.to_string())?;
    Ok(res)
}

pub fn load(spec_path: &Path) -> Result<(Ctx, PathBuf), String> {
    let txt = std::fs::read_to_string(spec_path).map_err(|e| format!("{}: {e}", spec_path.display()))?;
    let spec: Value = serde_json::from_str(&txt).map_err(|e| format!("{}: {e}", spec_path.display()))?;
    let base = spec_path.parent().unwrap_or(Path::new("."));
    let gdir = base.join(spec.get("gerbers").and_then(|v| v.as_str()).ok_or("spec needs \"gerbers\"")?);
    let out = base.join(spec.get("out").and_then(|v| v.as_str()).unwrap_or("fieldsolve-out"));
    let board = crate::gerber::load_dir(&gdir)?;
    let kicad = spec.get("coords").and_then(|v| v.as_str()).unwrap_or("kicad") == "kicad";
    Ok((Ctx { spec, board, kicad, zs: vec![] }, out))
}
