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
//!   "fields": {"freqs": [2.4e9]},
//!   "farfield": {"freqs": [2.44e9], "step_deg": 5, "gap": 3},
//!   "extent": "board", "decouple": [{"net": "VCC2V4", "from": "In1.Cu", "to": "In2.Cu"}],
//!   "preset": "antenna" | "rf", "loss": {"model": "wideband"}, "conductor": true, "edges": true,
//!   "mesh": {"no_slivers": true, "z_grade": 3.0}, "refine": [...], "coarse_copper_ok": [[x0, y0, x1, y1]] }
//! ```
//! Which keys to use, and the traps: docs/fieldsolve-settings.md. Every run prints the resolved model and
//! warnings (presets.rs) and records them in result.json "settings" / "checks".
//! A pad reference is `REF.PIN` (from the Gerber X2 `.P` attributes); a point is `[x, y]`.

use crate::fdtd::{gauss_excite, RunResult, RunSpec};
use crate::geom::{BBox, Prim, Shape, P2};
use crate::gerber::{Board, StackKind};
use crate::mesh::{Grid, Refine};
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

    /// Net of the copper at `p` on `layer` (the last Gerber object covering the point decides: dark with its
    /// net, clear with none). `None` for bare board or copper without a net attribute.
    pub fn net_at(&self, p: P2, layer: &str) -> Option<String> {
        let l = self.board.layers.iter().find(|l| l.name == layer)?;
        let mut net = None;
        for q in &l.prims {
            if q.shape.bbox().contains(p) && q.shape.contains(p) {
                net = if q.dark { Some(q.net.clone().unwrap_or_else(|| "copper (no net)".into())) } else { None };
            }
        }
        net
    }
}

/// Port terminal report: `+ NET (layer) / - NET (layer)`, and a warning when a vertical port's ends sit on
/// no ground net (a port referenced to a power plane or a floating island sees that conductor's own
/// resonances; reference it to ground, or tie the plane to ground with its decoupling capacitors).
fn port_nets(name: &str, plus: (Option<String>, &str), minus: (Option<String>, &str), vertical: bool) {
    let show = |n: &Option<String>| n.clone().unwrap_or_else(|| "no copper".into());
    eprintln!("fieldsolve: port {name}: + {} ({}) / - {} ({})", show(&plus.0), plus.1, show(&minus.0), minus.1);
    let gnd = |n: &Option<String>| n.as_deref().is_some_and(|s| s.to_ascii_uppercase().contains("GND"));
    if vertical && !gnd(&plus.0) && !gnd(&minus.0) {
        eprintln!("fieldsolve: WARNING: port {name}: neither end is on a GND net");
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
    /// stop when the port signals have decayed by end_db (spec "stop": "ports", the default for board
    /// runs; "energy" waits for the total field energy)
    pub port_stop: bool,
    /// in-plane ports and elements: (name, length of the bar over no copper (the pad gap) mm, bar width mm,
    /// partial self-inductance of that bare stretch nH). A vendor part model measured on its own fixture
    /// already contains its body's inductance across this gap, so `cosim` can subtract it ("deembed").
    pub bars: Vec<(String, f64, f64, f64)>,
    /// extra refinement boxes from the spec's `"refine"` (fidelity protection): (x0, y0, x1, y1 in the
    /// spec's coordinates, mm; cell mm; what)
    pub refined: Vec<([f64; 4], f64, String)>,
    /// set-up warnings (presets.rs: settings that are known traps, narrow copper in coarse cells)
    pub checks: Vec<String>,
}

/// Partial self-inductance (nH) of a thin flat strip, length l and width w in mm (Grover, t = 0):
/// L = mu0 l / (2 pi) [ln(2 l / w) + 1/2 + w / (3 l)].
pub fn strip_inductance_nh(l: f64, w: f64) -> f64 {
    if l <= 0.0 || w <= 0.0 {
        return 0.0;
    }
    (2e-7 * l * 1e-3 * ((2.0 * l / w).ln() + 0.5 + w / (3.0 * l))).max(0.0) * 1e9
}

/// Build layout, mesh and lumped elements from the spec.
pub fn setup(cx: &mut Ctx) -> Result<Setup, String> {
    let sp = cx.spec.clone();
    let mask = sp.get("solder_mask").and_then(|v| v.as_bool()).unwrap_or(false);
    // "dielectrics": [{"t": mm, "er": .., "tand": ..}, ...] top to bottom between the copper layers,
    // overriding the job file (or its absence)
    if let Some(d) = sp.get("dielectrics").and_then(|v| v.as_array()) {
        use crate::gerber::StackLayer;
        let cu: Vec<String> = cx.board.layers.iter().map(|l| l.name.clone()).collect();
        if d.len() + 1 != cu.len() {
            return Err(format!("\"dielectrics\" needs {} entries for {} copper layers", cu.len().saturating_sub(1), cu.len()));
        }
        let mut st = Vec::new();
        for (k, name) in cu.iter().enumerate() {
            st.push(StackLayer { name: name.clone(), kind: StackKind::Copper, thickness_mm: 0.035, eps_r: 1.0, loss_tan: 0.0 });
            if let Some(e) = d.get(k) {
                st.push(StackLayer {
                    name: format!("{}/{}", name, cu.get(k + 1).cloned().unwrap_or_default()),
                    kind: StackKind::Dielectric,
                    thickness_mm: num(e, "t", 0.2),
                    eps_r: num(e, "er", 4.4),
                    loss_tan: num(e, "tand", 0.02),
                });
            }
        }
        cx.board.stack = st;
    }
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
    // the PML sits beyond the air margin; it is npml cells of hmax-ish size.
    // "extent": "board" takes the whole board (plus `air`) into the domain for far-field runs: the roi
    // stays the fine-meshed region, cells grade up to max_cell outside it, and copper is rasterised
    // everywhere (on coarse cells small gaps merge: fine for planes, not for detail away from the roi)
    let extent_board = sp.get("extent").and_then(|v| v.as_str()) == Some("board");
    let dom = if extent_board {
        let mut b = board_bb;
        b.union(&roi);
        b.grow(air)
    } else {
        roi.grow(air)
    };
    let fine = if extent_board { roi } else { dom };

    // ports and elements -> lumped boxes, bridges on sheets, fixed mesh lines
    let mut sheets: Vec<Sheet> = Vec::new();
    for (name, z) in &zs {
        let l = cx.board.layers.iter().find(|l| &l.name == name).ok_or_else(|| format!("no Gerber for {name}"))?;
        let grown = dom.grow(1.0);
        let prims = l.prims.iter().filter(|p| p.shape.bbox().overlaps(&grown)).cloned().collect();
        sheets.push(Sheet { name: name.clone(), z: *z, prims });
    }
    let mut fx: Vec<f64> = vec![roi.x0, roi.x1, dom.x0, dom.x1];
    let mut fy: Vec<f64> = vec![roi.y0, roi.y1, dom.y0, dom.y1];
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
    // "decouple": [{"net": "VCC2V4", "from": "In1.Cu", "to": "In2.Cu"}]: every two-pad part outside the
    // roi with one pad on `net` and the other on a GND net becomes a vertical short from `from` to `to` at
    // that pad, tying a power plane to ground the way its capacitors do at RF (parts inside the roi are
    // left to explicit elements)
    for d in sp.get("decouple").and_then(|v| v.as_array()).cloned().unwrap_or_default() {
        let net = d.get("net").and_then(|v| v.as_str()).ok_or("decouple entry needs \"net\"")?.to_string();
        let from = d.get("from").and_then(|v| v.as_str()).ok_or("decouple entry needs \"from\"")?.to_string();
        let to = d.get("to").and_then(|v| v.as_str()).ok_or("decouple entry needs \"to\"")?.to_string();
        let mut parts: std::collections::BTreeMap<String, std::collections::BTreeMap<String, (String, BBox)>> = Default::default();
        for l in &cx.board.layers {
            for q in &l.prims {
                if let (true, Some((r, pin)), Some(n)) = (q.dark, &q.pad, &q.net) {
                    let e = parts.entry(r.clone()).or_default().entry(pin.clone()).or_insert((n.clone(), BBox::empty()));
                    e.1.union(&q.shape.bbox());
                }
            }
        }
        let mut added = 0;
        for (r, pins) in &parts {
            if pins.len() != 2 {
                continue;
            }
            let v: Vec<&(String, BBox)> = pins.values().collect();
            let gnd = |s: &str| s.to_ascii_uppercase().contains("GND");
            let hit = if v[0].0 == net && gnd(&v[1].0) { Some(v[0].1) } else if v[1].0 == net && gnd(&v[0].0) { Some(v[1].1) } else { None };
            let Some(bb) = hit else { continue };
            let c = P2::new(0.5 * (bb.x0 + bb.x1), 0.5 * (bb.y0 + bb.y1));
            if roi.contains(c) {
                continue;
            }
            let y = if cx.kicad { -c.y } else { c.y };
            // "esr": a capacitor's series resistance (ohm) instead of a perfect short
            let el = match d.get("esr").and_then(|v| v.as_f64()) {
                Some(esr) if esr > 0.0 => json!({"kind": "R", "r": esr, "at": [c.x, y], "from": from, "to": to, "size": "point"}),
                _ => json!({"kind": "short", "at": [c.x, y], "from": from, "to": to, "size": "point"}),
            };
            pend.push(Pending { name: format!("decouple {r}"), el, port: false });
            added += 1;
        }
        let what = match d.get("esr").and_then(|v| v.as_f64()) {
            Some(e) if e > 0.0 => format!("{e} ohm resistors"),
            _ => "shorts".to_string(),
        };
        eprintln!("fieldsolve: decouple {net}: {added} capacitors outside the roi as {from} -> {to} {what}");
    }
    // first pass: geometry of each element, so its lines can be fixed before meshing
    enum Geo {
        Vert { c: P2, half: [f64; 2], z0: f64, z1: f64, plus_hi: bool },
        Plane { a: P2, b: P2, axis: usize, w: f64, z: f64, layer: String },
        /// pads not on a common mesh line (a part rotated by e.g. 45°): a one-cell staircase bar along the
        /// mesh from each pad centre to a one-cell gap at the part's midpoint
        Diag { a: P2, b: P2, m: P2, axis: usize, z: f64, layer: String },
    }
    let mut geos = Vec::new();
    let mut bars = Vec::new();
    for pd in &pend {
        let e = &pd.el;
        if let Some(at) = e.get("at") {
            let from = e.get("from").and_then(|v| v.as_str()).ok_or_else(|| format!("{}: vertical element needs from/to layers", pd.name))?;
            let to = e.get("to").and_then(|v| v.as_str()).ok_or_else(|| format!("{}: vertical element needs from/to layers", pd.name))?;
            let (c, pbb) = cx.locate(at, Some(from))?;
            let (zf, zt) = (cx.z_of(from)?, cx.z_of(to)?);
            if pd.port {
                port_nets(&pd.name, (cx.net_at(c, from), from), (cx.net_at(c, to), to), true);
            }
            // footprint: explicit [w, h]; "point" for a single mesh line; default for a pad reference is the
            // pad's bounding box shrunk to 80 % (a single node column is a zero-radius wire whose inductance
            // depends on the cell size)
            let size = match e.get("size") {
                Some(Value::Array(a)) => [a[0].as_f64().unwrap_or(0.0), a[1].as_f64().unwrap_or(0.0)],
                Some(Value::String(s)) if s == "point" => [0.0, 0.0],
                _ if at.is_string() => [0.8 * (pbb.x1 - pbb.x0), 0.8 * (pbb.y1 - pbb.y0)],
                _ => [0.0, 0.0],
            };
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
            if pd.port {
                port_nets(&pd.name, (cx.net_at(a, &layer), &layer), (cx.net_at(b, &layer), &layer), false);
            }
            let axis = if (b.x - a.x).abs() >= (b.y - a.y).abs() { 0 } else { 1 };
            // a straight bar along one axis only reaches the far pad when both pads sit on (nearly) one line.
            // For a rotated part the bar would end beside the far pad, usually on the near pad's own copper
            // (shorting the port), so lay a staircase along the line between the pad centres instead.
            let off = if axis == 0 { (b.y - a.y).abs() } else { (b.x - a.x).abs() };
            if off > 1.5 * cell && !e.get("straight").and_then(|v| v.as_bool()).unwrap_or(false) {
                let z = cx.z_of(&layer)?;
                let n = 400;
                let len = ((b.x - a.x).powi(2) + (b.y - a.y).powi(2)).sqrt();
                let free: Vec<bool> = (0..n)
                    .map(|k| {
                        let t = (k as f64 + 0.5) / n as f64;
                        cx.net_at(P2::new(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t), &layer).is_none()
                    })
                    .collect();
                let bare = free.iter().filter(|&&f| f).count();
                let gap = bare as f64 / n as f64 * len;
                bars.push((pd.name.clone(), gap, cell, strip_inductance_nh(gap, cell)));
                // the element's gap goes in the middle of the longest copper-free stretch of the line (the
                // clearance the part body spans), not at the midpoint: that may lie in one pad's copper or
                // pour, which would bridge the gap
                let (mut best, mut run) = ((0usize, 0usize), None::<usize>);
                for k in 0..=n {
                    match (k < n && free[k], run) {
                        (true, None) => run = Some(k),
                        (false, Some(s0)) => {
                            if k - s0 > best.1 - best.0 {
                                best = (s0, k);
                            }
                            run = None;
                        }
                        _ => {}
                    }
                }
                let t = if best.1 > best.0 { 0.5 * (best.0 + best.1) as f64 / n as f64 } else { 0.5 };
                let m = P2::new(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t);
                geos.push(Geo::Diag { a, b, m, axis, z, layer });
                continue;
            }
            // bridge width: the smaller pad dimension across the axis
            let wa = if axis == 0 { ba.y1 - ba.y0 } else { ba.x1 - ba.x0 };
            let wb = if axis == 0 { bb.y1 - bb.y0 } else { bb.x1 - bb.x0 };
            let w = num(e, "width", wa.min(wb).max(cell));
            let z = cx.z_of(&layer)?;
            // a straight bridge on the line through pos along the axis
            let b2 = if axis == 0 { P2::new(b.x, a.y) } else { P2::new(a.x, b.y) };
            // the stretch of the bar over no copper: the gap between the pads that the part body spans
            {
                let n = 400;
                let bare = (0..n).filter(|&k| {
                    let t = (k as f64 + 0.5) / n as f64;
                    cx.net_at(P2::new(a.x + (b2.x - a.x) * t, a.y + (b2.y - a.y) * t), &layer).is_none()
                }).count();
                let gap = bare as f64 / n as f64 * ((b2.x - a.x).abs() + (b2.y - a.y).abs());
                bars.push((pd.name.clone(), gap, w, strip_inductance_nh(gap, w)));
            }
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
    // the absorbing layer's cells are graded up to hmax: with a whole-board domain give it its full width
    // beyond the air margin, so the far-field box keeps the requested distance from the structure
    let pml_w = if extent_board { npml as f64 * hmax.max(cell) } else { npml as f64 * hmax.min(1.0).max(cell) };
    // fixed lines closer than 3/4 of a cell would make a sliver cell that sets the time step for the whole
    // run (e.g. a port footprint edge 0.014 mm from a part's bridge edge); keep the earlier one
    let merge = |v: &[f64]| -> Vec<f64> {
        let tol = 0.75 * cell;
        let mut out: Vec<f64> = Vec::with_capacity(v.len());
        for &x in v {
            if !out.iter().any(|&y| (x - y).abs() < tol) {
                out.push(x);
            }
        }
        out
    };
    let (fx, fy) = (merge(&fx), merge(&fy));
    // "refine": [{"box": [x0, y0, x1, y1]} | {"pads": ["C20.1", ..]} | {"ports": true}, "cell": mm, "margin":
    // mm]: extra fine regions (a coarse fidelity policy protects ports, feeds, gaps and stubs this way). On
    // this tensor-product mesh a box's lines run through the whole domain in x and in y.
    let mut rx = vec![Refine { a: fine.x0, b: fine.x1, h: cell }];
    let mut ry = vec![Refine { a: fine.y0, b: fine.y1, h: cell }];
    let mut refined = Vec::new();
    for r in sp.get("refine").and_then(|v| v.as_array()).cloned().unwrap_or_default() {
        let h = num(&r, "cell", cell);
        let margin = num(&r, "margin", 0.0);
        let mut boxes: Vec<(BBox, String)> = Vec::new();
        if let Some(a) = r.get("box").and_then(|v| v.as_array()) {
            let v: Vec<f64> = a.iter().filter_map(|x| x.as_f64()).collect();
            if v.len() != 4 {
                return Err("refine box needs [x0, y0, x1, y1]".into());
            }
            let (p, q) = (cx.pt(v[0], v[1]), cx.pt(v[2], v[3]));
            boxes.push((BBox { x0: p.x.min(q.x), y0: p.y.min(q.y), x1: p.x.max(q.x), y1: p.y.max(q.y) }, r.get("what").and_then(|v| v.as_str()).unwrap_or("box").to_string()));
        }
        for pad in r.get("pads").and_then(|v| v.as_array()).cloned().unwrap_or_default() {
            let (_, b) = cx.locate(&pad, None)?;
            boxes.push((b, format!("pad {}", pad.as_str().unwrap_or("?"))));
        }
        if r.get("ports").and_then(|v| v.as_bool()).unwrap_or(false) {
            for pd in pend.iter().filter(|pd| pd.port) {
                let layer = pd.el.get("layer").and_then(|v| v.as_str()).or_else(|| pd.el.get("from").and_then(|v| v.as_str()));
                let mut b = BBox::empty();
                for k in ["pos", "neg", "at"] {
                    if let Some(v) = pd.el.get(k) {
                        b.union(&cx.locate(v, layer)?.1);
                    }
                }
                boxes.push((b, format!("port {}", pd.name)));
            }
        }
        for (b, what) in boxes {
            let b = b.grow(margin);
            rx.push(Refine { a: b.x0, b: b.x1, h });
            ry.push(Refine { a: b.y0, b: b.y1, h });
            let (p, q) = (cx.pt(b.x0, b.y0), cx.pt(b.x1, b.y1));
            refined.push(([p.x.min(q.x), p.y.min(q.y), p.x.max(q.x), p.y.max(q.y)], h, what));
        }
    }
    // "mesh": static refinement from the layout (meshpol.rs): no slivers, z grading, copper-edge bands
    let mp = crate::meshpol::MeshPolicy::from_spec(&sp);
    let (bx, by) = mp.copper_bands(&sheets, &fine);
    rx.extend(bx);
    ry.extend(by);
    let mut gx = crate::mesh::lines_1d_slack(dom.x0 - pml_w, dom.x1 + pml_w, &fx, &rx, 1.3, hmax, &[], mp.slack());
    let mut gy = crate::mesh::lines_1d_slack(dom.y0 - pml_w, dom.y1 + pml_w, &fy, &ry, 1.3, hmax, &[], mp.slack());
    let mc: Vec<(f64, f64, usize)> = slabs.iter().map(|s| (s.z0, s.z1, min_layer)).collect();
    let zref = mp.z_refine(zlo, zhi, &fz, cell_z);
    let gz = crate::mesh::lines_1d_slack(zlo - air_z - pml_w, zhi + air_z + pml_w, &fz, &zref, 1.3, hmax, &mc, mp.slack());
    // keep only fixed lines that fall inside the domain
    gx.dedup();
    gy.dedup();
    // known traps in the settings, and narrow copper left in coarse cells next to the fine region
    let mut checks = crate::presets::spec_checks(&sp);
    if extent_board {
        // "coarse_copper_ok": [[x0, y0, x1, y1], ...] marks areas reviewed as not carrying RF current (spec
        // coordinates); clusters inside one are not reported
        let ok: Vec<BBox> = sp
            .get("coarse_copper_ok")
            .and_then(|v| v.as_array())
            .map(|a| {
                a.iter()
                    .filter_map(|r| {
                        let r = r.as_array()?;
                        let (p, q) = (cx.pt(r.first()?.as_f64()?, r.get(1)?.as_f64()?), cx.pt(r.get(2)?.as_f64()?, r.get(3)?.as_f64()?));
                        Some(BBox { x0: p.x.min(q.x), y0: p.y.min(q.y), x1: p.x.max(q.x), y1: p.y.max(q.y) })
                    })
                    .collect()
            })
            .unwrap_or_default();
        let found: Vec<_> = crate::presets::coarse_copper(&sheets, &gx, &gy, &roi, &dom, cell, 2.5 * cell, 1.0, 8.0)
            .into_iter()
            .filter(|c| !ok.iter().any(|o| o.x0 <= c.1.x0 && o.y0 <= c.1.y0 && o.x1 >= c.1.x1 && o.y1 >= c.1.y1))
            .collect();
        if !found.is_empty() {
            let list: Vec<String> = found
                .iter()
                .take(6)
                .map(|(l, b, a, h)| {
                    let (p, q) = (cx.pt(b.x0, b.y0), cx.pt(b.x1, b.y1));
                    format!("{l} x {:.2}-{:.2} y {:.2}-{:.2} ({a:.1} mm², cells up to {h:.2} mm)", p.x.min(q.x), p.x.max(q.x), p.y.min(q.y), p.y.max(q.y))
                })
                .collect();
            checks.push(format!(
                "narrow copper (traces, meanders, antenna arms) in cells coarser than {:.3} mm next to the fine region: {}{}. \
                 If any of it carries RF current (antenna, feed, matching), extend \"roi\" to cover it or add a \"refine\" box; \
                 list areas you have checked in \"coarse_copper_ok\".",
                2.5 * cell,
                list.join("; "),
                if found.len() > 6 { format!("; and {} more", found.len() - 6) } else { String::new() }
            ));
        }
    }
    eprintln!("fieldsolve: {}", crate::presets::summary(&sp).1);
    for c in &checks {
        eprintln!("fieldsolve: WARNING: {c}");
    }
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
                    Prim { shape: Shape::Polygon(pts), dark: true, net: Some(format!("~{}", pd.name)), pad: None }
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
                // clear the gap column over the bar's full width and a cell beyond: copper around the pads
                // (a ground pour reaching into the gap) must not bridge it outside the element's edges, which
                // cover only the middle 80 % (below). Only edges whose midpoint lies in the gap are cleared.
                {
                    let (u0, u1) = (g0 + 1e-6, g1 - 1e-6);
                    let (v0, v1) = (oc - w / 2.0 - cell, oc + w / 2.0 + cell);
                    let pts = if *axis == 0 {
                        vec![P2::new(u0, v0), P2::new(u1, v0), P2::new(u1, v1), P2::new(u0, v1)]
                    } else {
                        vec![P2::new(v0, u0), P2::new(v1, u0), P2::new(v1, u1), P2::new(v0, u1)]
                    };
                    sheet.prims.push(Prim { shape: Shape::Polygon(pts), dark: false, net: None, pad: None });
                }
                let mut lo = [0.0; 3];
                let mut hi = [0.0; 3];
                lo[*axis] = g0;
                hi[*axis] = g1;
                let o = 1 - *axis;
                // the element's edges sit inside the bar: 80 % of its width, as for vertical ports. Mesh
                // lines fall on the bar's (and pad's) edges, and copper edges lying on a polygon boundary
                // may or may not be metal in the raster; an element edge with no metal at one end carries
                // its share of the source and the voltage probe without being part of the circuit, which
                // shows up as a frequency-flat error in the port's resistance (seen as 2-14 % "loss")
                let wp = 0.8 * w;
                lo[o] = oc - wp / 2.0;
                hi[o] = oc + wp / 2.0;
                lo[2] = *z;
                hi[2] = *z;
                lumped.push(Lumped { name: pd.name.clone(), axis: *axis, lo, hi, plus_at_hi: pa > pb, r: if kind == "short" { 0.0 } else { r }, c, port: pd.port });
            }
            Geo::Diag { a, b, m, axis, z, layer } => {
                let sheet = sheets.iter_mut().find(|s| &s.name == layer).unwrap();
                let (ax, o) = (*axis, 1 - *axis);
                let gl = [&gx, &gy];
                let pc = |p: &P2, k: usize| if k == 0 { p.x } else { p.y };
                let node = |p: &P2| [crate::mesh::snap(&gx, p.x), crate::mesh::snap(&gy, p.y)];
                let mut i = crate::mesh::snap(gl[ax], pc(m, ax));
                if i + 1 >= gl[ax].len() {
                    i -= 1;
                }
                let jm = crate::mesh::snap(gl[o], pc(m, o));
                // the gap edge runs along `ax` from line i to i+1 on cross line jm; each pad's staircase ends
                // on its own side of it
                let a_low = pc(a, ax) < pc(b, ax);
                let mut ga = [0usize; 2];
                let mut gb = [0usize; 2];
                ga[ax] = if a_low { i } else { i + 1 };
                gb[ax] = if a_low { i + 1 } else { i };
                ga[o] = jm;
                gb[o] = jm;
                // a mesh edge as a thin strip of copper: only that edge's midpoint falls inside it
                let mut edge_bar = |p: [usize; 2], q: [usize; 2]| {
                    let (x0, x1) = (gx[p[0]].min(gx[q[0]]), gx[p[0]].max(gx[q[0]]));
                    let (y0, y1) = (gy[p[1]].min(gy[q[1]]), gy[p[1]].max(gy[q[1]]));
                    let tx = 0.3 * (gx[(p[0] + 1).min(gx.len() - 1)] - gx[p[0].saturating_sub(1)]).abs() / 2.0;
                    let ty = 0.3 * (gy[(p[1] + 1).min(gy.len() - 1)] - gy[p[1].saturating_sub(1)]).abs() / 2.0;
                    let pts = vec![P2::new(x0 - tx, y0 - ty), P2::new(x1 + tx, y0 - ty), P2::new(x1 + tx, y1 + ty), P2::new(x0 - tx, y1 + ty)];
                    sheet.prims.push(Prim { shape: Shape::Polygon(pts), dark: true, net: Some(format!("~{}", pd.name)), pad: None });
                };
                // 4-connected walk in index space from s to t (alternating steps keep it on the straight line)
                let mut walk = |s: [usize; 2], t: [usize; 2]| {
                    let (mut cur, d) = (s, [t[0] as i64 - s[0] as i64, t[1] as i64 - s[1] as i64]);
                    let n = d[0].abs() + d[1].abs();
                    for k in 1..=n {
                        // ideal position after k unit steps along the straight line
                        let f = k as f64 / n as f64;
                        let want = [s[0] as f64 + d[0] as f64 * f, s[1] as f64 + d[1] as f64 * f];
                        let ex = (want[0] - cur[0] as f64).abs() - (want[1] - cur[1] as f64).abs();
                        let mut nxt = cur;
                        if ex >= 0.0 && cur[0] != t[0] {
                            nxt[0] = (cur[0] as i64 + d[0].signum()) as usize;
                        } else if cur[1] != t[1] {
                            nxt[1] = (cur[1] as i64 + d[1].signum()) as usize;
                        } else {
                            nxt[0] = (cur[0] as i64 + d[0].signum()) as usize;
                        }
                        edge_bar(cur, nxt);
                        cur = nxt;
                    }
                };
                if kind == "short" && !pd.port {
                    walk(node(a), node(b));
                    continue;
                }
                walk(node(a), ga);
                walk(gb, node(b));
                let g0 = gl[ax][i];
                let g1 = gl[ax][i + 1];
                let cl = gl[o][jm];
                let h = 0.5 * (gl[o][(jm + 1).min(gl[o].len() - 1)] - gl[o][jm.saturating_sub(1)]);
                // clear the gap edge itself (and only it) of any copper bridging it
                {
                    let (u0, u1, v0, v1) = (g0 + 1e-6, g1 - 1e-6, cl - 0.5 * h, cl + 0.5 * h);
                    let pts = if ax == 0 {
                        vec![P2::new(u0, v0), P2::new(u1, v0), P2::new(u1, v1), P2::new(u0, v1)]
                    } else {
                        vec![P2::new(v0, u0), P2::new(v1, u0), P2::new(v1, u1), P2::new(v0, u1)]
                    };
                    sheet.prims.push(Prim { shape: Shape::Polygon(pts), dark: false, net: None, pad: None });
                }
                let mut lo = [0.0; 3];
                let mut hi = [0.0; 3];
                lo[ax] = g0;
                hi[ax] = g1;
                lo[o] = cl - 0.4 * h;
                hi[o] = cl + 0.4 * h;
                lo[2] = *z;
                hi[2] = *z;
                lumped.push(Lumped { name: pd.name.clone(), axis: ax, lo, hi, plus_at_hi: !a_low, r: if kind == "short" { 0.0 } else { r }, c, port: pd.port });
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
    let layout = Layout { sheets, slabs, barrels, outline: cx.board.outline.clone(), copper_window: if extent_board { None } else { Some(dom) } };
    let fmin = num(&sp, "fmin", fmax / 300.0);
    let nf = num(&sp, "nf", 300.0) as usize;
    let sheet_names: Vec<String> = layout.sheets.iter().map(|s| s.name.clone()).collect();
    Ok(Setup {
        grid,
        layout,
        lumped,
        opts: BuildOpts {
            pml: [npml; 6],
            f_loss: num(&sp, "f_loss", 1e9),
            fmax,
            cfl: 0.95,
            // "loss": {"model": "wideband", "f_lo": .., "f_hi": .., "f_ref": 1e9, "copper": true}: flat
            // tan(delta) over the band (Debye poles), plus the plane-pair copper loss if asked
            wideband: sp.get("loss").filter(|l| l.get("model").and_then(|m| m.as_str()) == Some("wideband")).map(|l| crate::disp::Wideband {
                f_lo: num(l, "f_lo", fmin.max(0.1e9)),
                f_hi: num(l, "f_hi", fmax),
                f_ref: num(l, "f_ref", 1e9),
                copper: l.get("copper").and_then(|c| c.as_bool()).unwrap_or(false),
            }),
            // "conductor": true or {"sigma": 5.8e7, "roughness_um": 1.0, "thickness_um": 35, "f_lo", "f_hi"}:
            // surface-impedance copper (cond.rs); thickness per layer from the stack-up unless given
            conductor: sp.get("conductor").filter(|c| c.as_bool() != Some(false)).map(|c| {
                let o = |k: &str| c.get(k).and_then(|v| v.as_f64());
                let metals = sheet_names
                    .iter()
                    .map(|name| {
                        let t_mm = o("thickness_um").map(|t| t * 1e-3).unwrap_or_else(|| {
                            cx.board.stack.iter().find(|l| l.kind == crate::gerber::StackKind::Copper && &l.name == name).map(|l| l.thickness_mm).filter(|&t| t > 0.0).unwrap_or(0.035)
                        });
                        crate::cond::Metal { sigma: o("sigma").unwrap_or(5.8e7), t: t_mm * 1e-3, rough: o("roughness_um").unwrap_or(1.0) * 1e-6 }
                    })
                    .collect();
                crate::cond::Conductor { metals, f_lo: o("f_lo").unwrap_or(fmin.max(0.1e9)), f_hi: o("f_hi").unwrap_or(fmax) }
            }),
            // "edges": true or {"conformal": .., "singular": .., ...}: sub-cell copper edges (edge.rs)
            edge: sp.get("edges").and_then(crate::edge::EdgeOpts::from_json),
        },
        freqs: linspace(fmin, fmax, nf),
        fmax,
        end_db: num(&sp, "end_db", -70.0),
        max_steps: num(&sp, "max_steps", 400_000.0) as usize,
        port_stop: sp.get("stop").and_then(|v| v.as_str()).unwrap_or("ports") != "energy",
        bars,
        refined,
        checks,
    })
}

/// The resolved raster of a built problem, so another solver can be given exactly the same structure:
/// - `metal`: runs of perfect-conductor edges `[axis, a0, a1, b, c]`, edges a0..=a1 along `axis` at the
///   transverse node indices (b, c) of the next two axes in cyclic order (x: (j, k), y: (k, i), z: (i, j)).
///   This is every interior edge with cb = 0: rasterised copper after conflict resolution, via barrels and
///   shorted lumped elements. Edges carrying a lumped element are not metal and are not listed.
/// - `cells`: runs of dielectric cells `[slab, i0, i1, j, k0, k1]`, cells i0..=i1 of row j over the z cells
///   k0..=k1 of one slab (cells inside the outline, by centre, as `build` assigns them).
pub fn raster_json(st: &Setup, p: &Problem) -> Value {
    let g = &st.grid;
    let nn = [p.nx, p.ny, p.nz];
    let mut metal: Vec<[usize; 5]> = Vec::new();
    for a in 0..3 {
        let (b, c) = ((a + 1) % 3, (a + 2) % 3);
        for ic in 1..nn[c].saturating_sub(1) {
            for ib in 1..nn[b].saturating_sub(1) {
                let mut run: Option<usize> = None;
                for ia in 0..nn[a] {
                    let mut pos = [0usize; 3];
                    pos[a] = ia;
                    pos[b] = ib;
                    pos[c] = ic;
                    let on = ia + 1 < nn[a] && p.cb[a][g.idx(pos[0], pos[1], pos[2])] == 0.0;
                    match (on, run) {
                        (true, None) => run = Some(ia),
                        (false, Some(s)) => {
                            metal.push([a, s, ia - 1, ib, ic]);
                            run = None;
                        }
                        _ => {}
                    }
                }
            }
        }
    }
    let (inside, slab_of_k) = crate::model::cell_materials(g, &st.layout);
    let (nx, ny) = (p.nx, p.ny);
    // contiguous z ranges of each slab
    let mut zr: Vec<(usize, usize, usize)> = Vec::new();
    for (k, s) in slab_of_k.iter().enumerate() {
        if let Some(s) = *s {
            match zr.last_mut() {
                Some(l) if l.0 == s && l.2 + 1 == k => l.2 = k,
                _ => zr.push((s, k, k)),
            }
        }
    }
    let mut cells: Vec<[usize; 6]> = Vec::new();
    for &(s, k0, k1) in &zr {
        for j in 0..ny.saturating_sub(1) {
            let mut run: Option<usize> = None;
            for i in 0..nx {
                let on = i + 1 < nx && inside[j * nx + i];
                match (on, run) {
                    (true, None) => run = Some(i),
                    (false, Some(a)) => {
                        cells.push([s, a, i - 1, j, k0, k1]);
                        run = None;
                    }
                    _ => {}
                }
            }
        }
    }
    json!({"metal": metal, "cells": cells})
}

/// Outline of a capsule (segment a-b of radius r): the cap around b on the far side from a, then the cap
/// around a on the far side from b (13 points each).
pub fn capsule_outline(a: P2, b: P2, r: f64) -> Vec<P2> {
    let ang = (b.y - a.y).atan2(b.x - a.x);
    let cap = |c: P2, from: f64| (0..=12).map(move |k| {
        let t = from + k as f64 / 12.0 * std::f64::consts::PI;
        P2::new(c.x + r * t.cos(), c.y + r * t.sin())
    });
    cap(b, ang - std::f64::consts::FRAC_PI_2).chain(cap(a, ang + std::f64::consts::FRAC_PI_2)).collect()
}

/// Scene file: exactly what was simulated, for an independent solver (tools/openems_scene.py).
pub fn scene_json(st: &Setup, p: &Problem, exc_port: Option<usize>) -> Value {
    let poly = |p: &Prim| -> Vec<[f64; 2]> {
        let pts: Vec<P2> = match &p.shape.outline_shape() {
            Shape::Macro(_) => Vec::new(), // outline_shape never returns a macro
            Shape::Polygon(v) => v.clone(),
            Shape::Circle { c, r } => (0..32).map(|k| {
                let t = k as f64 / 32.0 * std::f64::consts::TAU;
                P2::new(c.x + r * t.cos(), c.y + r * t.sin())
            }).collect(),
            Shape::Capsule { a, b, r } => capsule_outline(*a, *b, *r),
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
        "raster": raster_json(st, p),
    })
}

pub type Backend<'a> = &'a dyn Fn(&Problem, &RunSpec) -> RunResult;

/// What to run on a set-up: which ports to excite (names, or all), field snapshot frequencies, L-fit band.
#[derive(Clone, Debug, Default)]
pub struct RunOpts {
    /// the resolved model (presets::summary), recorded in result.json "settings"
    pub settings: Value,
    pub excite_all: bool,
    pub excite: Vec<String>,
    pub fields: Vec<f64>,
    pub l_band: Option<(f64, f64)>,
    /// Report positions in KiCad's frame (y negated).
    pub kicad: bool,
    /// Time-domain surface-current snapshots of these copper layers (first excitation only), every
    /// `snap_ps` picoseconds, written to out/frames/.
    pub snap_layers: Vec<String>,
    pub snap_ps: f64,
    pub snap_until_ps: f64,
    /// Near-to-far-field box and frequencies (`"farfield"` in the spec).
    pub farfield: Option<crate::nf2ff::FarOpts>,
    /// The objective tracked during the runs and summarised in result.json (`"objective"` in the spec).
    pub objective: Option<ObjectiveOpts>,
    /// Skip scene.json (the openEMS cross-check export).
    pub no_scene: bool,
    /// The resolved fidelity policy, copied into result.json (`fieldsolve run --fast`).
    pub fidelity: Value,
}

/// `"objective": {"port": "C20", "band": [2.40e9, 2.44e9, 2.48e9], "harmonics": [4.88e9], "tol_db": 0.05,
/// "window_periods": 4, "min_settle_periods": 6}`: the band's radiated power per watt incident at `port` (other ports matched),
/// tracked at every energy check of every excitation together with each port's |S_kj|^2. With `"stop":
/// "objective"` a run ends when all of them have stayed within a spread of `tol_db` over `window_periods`
/// periods of the lowest band frequency, after the source and `min_settle_periods` more periods.
#[derive(Clone, Debug, Default)]
pub struct ObjectiveOpts {
    pub port: Option<String>,
    pub band: Vec<f64>,
    /// further frequencies tracked for the stop (e.g. harmonic constraints), not part of the band score
    pub harmonics: Vec<f64>,
    pub stop: bool,
    pub tol_db: f64,
    pub window_periods: f64,
    pub min_settle_periods: f64,
    /// Late-time extrapolation (`"extrapolate": true` or `{"snapshots": 10, "tol_db": 0.05,
    /// "stable_checks": 4, "stop": true}`): the extrapolated objective is recorded at every check, and with
    /// `stop` (the default) and the objective stop it can end the run early (see `extrap`).
    pub extrapolate: Option<crate::extrap::Opts>,
}

impl ObjectiveOpts {
    pub fn from_spec(spec: &Value) -> Option<ObjectiveOpts> {
        let o = spec.get("objective")?;
        let band: Vec<f64> = o.get("band")?.as_array()?.iter().filter_map(|x| x.as_f64()).collect();
        if band.is_empty() {
            return None;
        }
        Some(ObjectiveOpts {
            port: o.get("port").and_then(|v| v.as_str()).map(String::from),
            band,
            harmonics: o.get("harmonics").and_then(|v| v.as_array()).map(|a| a.iter().filter_map(|x| x.as_f64()).collect()).unwrap_or_default(),
            stop: spec.get("stop").and_then(|v| v.as_str()) == Some("objective"),
            tol_db: num(o, "tol_db", 0.05),
            window_periods: num(o, "window_periods", 4.0),
            min_settle_periods: num(o, "min_settle_periods", 6.0),
            extrapolate: match o.get("extrapolate") {
                Some(Value::Bool(true)) => Some(crate::extrap::Opts { stop: true, ..Default::default() }),
                Some(x @ Value::Object(_)) => {
                    let d = crate::extrap::Opts::default();
                    Some(crate::extrap::Opts {
                        snapshots: num(x, "snapshots", d.snapshots as f64) as usize,
                        window: num(x, "window", d.window as f64) as usize,
                        order: num(x, "order", d.order as f64) as usize,
                        tol_db: num(x, "tol_db", d.tol_db),
                        stable_checks: num(x, "stable_checks", d.stable_checks as f64) as usize,
                        stop: x.get("stop").and_then(|v| v.as_bool()).unwrap_or(true),
                        ..d
                    })
                }
                _ => None,
            },
        })
    }

    /// Floors (dB) below which a component counts as flat: accepted power and radiated power -30 dB of the
    /// incident power, coupling to another port -20 dB.
    const FLOOR_RAD: f64 = 1e-3;
    const FLOOR_COUPLING: f64 = 1e-2;

    /// Component names, in the order `metric` returns them.
    pub fn components(&self, names: &[String], j: usize, radiation: bool) -> Vec<String> {
        let mut c = Vec::new();
        for f in self.band.iter().chain(&self.harmonics) {
            let g = format!("{:.3}GHz", f * 1e-9);
            for (k, n) in names.iter().enumerate() {
                c.push(if k == j { format!("accepted {n} {g}") } else { format!("|S {n},{}|^2 {g}", names[j]) });
            }
            if radiation {
                c.push(format!("radiated {g}"));
            }
        }
        c
    }

    /// The tracked components (dB) of the run exciting port `j`: per band frequency, the power accepted by
    /// port j and coupled to every other port, and the box's Poynting flux, all per watt incident.
    pub fn metric(&self, j: usize, p: &Problem, exc: &[f32], far: Option<(&crate::nf2ff::Plan, &[f64])>) -> crate::fdtd::ObjectiveMetric {
        let band: Vec<f64> = self.band.iter().chain(&self.harmonics).copied().collect();
        let gi: Vec<Option<usize>> = band.iter().map(|f| far.and_then(|(_, ff)| ff.iter().position(|x| (x - f).abs() <= 1e-6 * f))).collect();
        let plan = far.map(|(pl, _)| pl.clone());
        let grid = p.grid.clone();
        let rs: Vec<f64> = p.ports.iter().map(|q| q.r).collect();
        let (exc, dt) = (exc.to_vec(), p.dt);
        let db = |x: f64, floor: f64| 10.0 * x.max(floor).log10();
        crate::fdtd::ObjectiveMetric(std::sync::Arc::new(move |_steps, v, gathered| {
            let z0 = 50.0;
            let k = 1.0 / (2.0 * f64::sqrt(z0));
            let spec: Vec<(Vec<Cx>, Vec<Cx>)> = (0..v.len()).map(|q| port_spectra(&v[q], if q == j { &exc } else { &[] }, rs[q], dt, &band)).collect();
            let mut out = Vec::with_capacity(band.len() * (v.len() + 1));
            for (fi, gfi) in gi.iter().enumerate() {
                let a = (spec[j].0[fi] + spec[j].1[fi] * z0) * k;
                let pinc = (a.abs().powi(2) / 2.0).max(1e-300);
                for q in 0..v.len() {
                    let b = (spec[q].0[fi] - spec[q].1[fi] * z0) * k;
                    let s2 = b.abs().powi(2) / 2.0 / pinc;
                    out.push(if q == j { db(1.0 - s2, Self::FLOOR_RAD) } else { db(s2, Self::FLOOR_COUPLING) });
                }
                if let (Some(g), Some(pl)) = (gfi, &plan) {
                    if let Some(d) = gathered.get(*g) {
                        out.push(db(crate::nf2ff::box_flux(pl, &grid, d) / pinc, Self::FLOOR_RAD));
                    }
                }
            }
            out
        }))
    }
}

impl RunOpts {
    pub fn from_spec(spec: &Value) -> RunOpts {
        RunOpts {
            settings: crate::presets::summary(spec).0,
            excite_all: spec.get("excite").and_then(|v| v.as_str()) == Some("all"),
            excite: spec.get("excite").and_then(|v| v.as_array()).map(|a| a.iter().filter_map(|x| x.as_str().map(String::from)).collect()).unwrap_or_default(),
            fields: spec.pointer("/fields/freqs").and_then(|v| v.as_array()).map(|a| a.iter().filter_map(|x| x.as_f64()).collect()).unwrap_or_default(),
            l_band: spec.get("l_band").and_then(|v| v.as_array()).map(|a| (a[0].as_f64().unwrap_or(0.0), a[1].as_f64().unwrap_or(0.0))),
            kicad: spec.get("coords").and_then(|v| v.as_str()).unwrap_or("kicad") == "kicad",
            snap_layers: spec.pointer("/frames/layers").and_then(|v| v.as_array()).map(|a| a.iter().filter_map(|x| x.as_str().map(String::from)).collect()).unwrap_or_default(),
            snap_ps: spec.pointer("/frames/every_ps").and_then(|v| v.as_f64()).unwrap_or(5.0),
            snap_until_ps: spec.pointer("/frames/until_ps").and_then(|v| v.as_f64()).unwrap_or(0.0),
            farfield: crate::nf2ff::FarOpts::from_spec(spec),
            objective: ObjectiveOpts::from_spec(spec),
            no_scene: spec.get("scene").and_then(|v| v.as_bool()) == Some(false),
            fidelity: spec.get("fidelity_record").cloned().unwrap_or(Value::Null),
        }
    }
}

/// Run every excitation of a board spec, write results to `out`.
pub fn run(cx: &mut Ctx, out: &Path, backend: Backend) -> Result<Value, String> {
    let st = setup(cx)?;
    run_setup(&st, out, backend, &RunOpts::from_spec(&cx.spec))
}

/// Z = Z0 (I + S)(I - S)^-1 for an N-port at one frequency.
pub fn z_from_s(s: &[Vec<Cx>], z0: f64) -> Option<Vec<Vec<Cx>>> {
    let n = s.len();
    let one = Cx::new(1.0, 0.0);
    // a = I - S, b = I + S; solve X a = b  <=>  a^T X^T = b^T (Gauss-Jordan on augmented [a^T | b^T])
    let mut m: Vec<Vec<Cx>> = (0..n)
        .map(|i| {
            let mut row: Vec<Cx> = (0..n).map(|j| if i == j { one - s[j][i] } else { Cx::default() - s[j][i] }).collect();
            row.extend((0..n).map(|j| if i == j { one + s[j][i] } else { s[j][i] }));
            row
        })
        .collect();
    for c in 0..n {
        let piv = (c..n).max_by(|&a, &b| m[a][c].abs().partial_cmp(&m[b][c].abs()).unwrap())?;
        if m[piv][c].abs() < 1e-30 {
            return None;
        }
        m.swap(c, piv);
        let d = m[c][c];
        for k in 0..2 * n {
            m[c][k] = m[c][k] / d;
        }
        for r in 0..n {
            if r != c {
                let f = m[r][c];
                for k in 0..2 * n {
                    let t = m[c][k] * f;
                    m[r][k] = m[r][k] - t;
                }
            }
        }
    }
    // m[:, n..] = X^T
    Some((0..n).map(|i| (0..n).map(|j| m[j][n + i] * z0).collect()).collect())
}

/// Run a prepared set-up, write `scene.json`, `result.json`, Touchstone and fields to `out`.
pub fn run_setup(st: &Setup, out: &Path, backend: Backend, ro: &RunOpts) -> Result<Value, String> {
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
    let bridges = crate::connect::check(&st.grid, &st.layout);
    let bridge_lines = crate::connect::describe(&bridges, ro.kicad);
    for l in bridge_lines.iter().take(12) {
        eprintln!("fieldsolve: WARNING: mesh short: {l}");
    }
    if bridge_lines.len() > 12 {
        eprintln!("fieldsolve: WARNING: ... {} more mesh shorts (all in result.json)", bridge_lines.len() - 12);
    }
    let names: Vec<String> = p.ports.iter().map(|x| x.name.clone()).collect();
    let excite: Vec<usize> = if ro.excite_all {
        (0..np).collect()
    } else if !ro.excite.is_empty() {
        ro.excite.iter().filter_map(|n| names.iter().position(|m| m == n)).collect()
    } else {
        vec![0]
    };
    let fields_f = ro.fields.clone();
    let zmm: Vec<f64> = st.grid.z.iter().map(|v| v * 1e3).collect();
    let snap_sheets: Vec<&crate::model::Sheet> = ro.snap_layers.iter().filter_map(|n| st.layout.sheets.iter().find(|s| &s.name == n)).collect();
    let snap_k: Vec<usize> = snap_sheets.iter().map(|s| crate::mesh::snap(&zmm, s.z)).collect();
    let snap_every = if snap_k.is_empty() { 0 } else { ((ro.snap_ps * 1e-12 / p.dt).round() as usize).max(1) };
    let exc = gauss_excite(0.0, st.fmax, p.dt);
    let dft_every = ((1.0 / (2.0 * st.fmax * 1.5)) / p.dt).floor().max(1.0) as usize;
    if !ro.no_scene {
        std::fs::write(out.join("scene.json"), serde_json::to_string(&scene_json(st, &p, excite.first().copied())).unwrap()).map_err(|e| e.to_string())?;
    }
    let freqs = &st.freqs;
    let z0 = 50.0;
    let mut far = match &ro.farfield {
        Some(fo) => {
            let plan = crate::nf2ff::auto_plan(&p, st.opts.pml, fo.gap)?;
            let bad = crate::nf2ff::non_free_edges(&p, &plan);
            let ax = [&p.grid.x, &p.grid.y, &p.grid.z];
            let mm = |a: usize, i: usize| ax[a][i] * 1e3;
            eprintln!(
                "fieldsolve: far-field box x {:.1}..{:.1}, y {:.1}..{:.1}, z {:.1}..{:.1} mm, {} face samples",
                mm(0, plan.lo[0]), mm(0, plan.hi[0]), mm(1, plan.lo[1]), mm(1, plan.hi[1]), mm(2, plan.lo[2]), mm(2, plan.hi[2]),
                plan.gather.iter().map(|g| g.len()).sum::<usize>()
            );
            if bad > 0 {
                eprintln!("fieldsolve: WARNING: far-field box: {bad} tangential E edges on its faces are not free space (the structure crosses the box); the far field is not valid");
            }
            Some(crate::nf2ff::Collector::new(fo.clone(), plan, np))
        }
        None => None,
    };
    let mut s_mat: Vec<Vec<Vec<Cx>>> = vec![vec![vec![Cx::default(); freqs.len()]; np]; np];
    let mut zself: Vec<Option<Vec<Cx>>> = vec![None; np];
    let mut vall: Vec<Vec<Vec<Cx>>> = vec![vec![]; np];
    let mut iall: Vec<Vec<Vec<Cx>>> = vec![vec![]; np];
    let mut runs = Vec::new();
    // objective stop: window and minimum settling time in periods of the band's lowest frequency
    let obj_stop = ro.objective.as_ref().is_some_and(|o| o.stop);
    let period_steps = |n: f64| ro.objective.as_ref().map(|o| (n / o.band.iter().cloned().fold(f64::INFINITY, f64::min) / p.dt).ceil() as usize).unwrap_or(0);
    let obj_window = period_steps(ro.objective.as_ref().map(|o| o.window_periods).unwrap_or(0.0)).max(3000);
    let obj_min_steps = exc.len() + period_steps(ro.objective.as_ref().map(|o| o.min_settle_periods).unwrap_or(0.0));
    let mut obj_hist: Vec<Option<(Vec<String>, Vec<(usize, Vec<f64>)>)>> = vec![None; np];
    for (q, &j) in excite.iter().enumerate() {
        let spec = RunSpec {
            excite: Some(j),
            exc: exc.clone(),
            max_steps: st.max_steps,
            end_db: st.end_db,
            check_every: 1000,
            dft_freqs: if q == 0 { fields_f.clone() } else { vec![] },
            dft_every,
            snap_every: if q == 0 { snap_every } else { 0 },
            snap_k: snap_k.clone(),
            snap_until: if ro.snap_until_ps > 0.0 { (ro.snap_until_ps * 1e-12 / p.dt) as usize } else { 0 },
            gather_freqs: far.as_ref().map(|c| c.opts.freqs.clone()).unwrap_or_default(),
            gather: far.as_ref().map(|c| c.plan.gather.clone()).unwrap_or_default(),
            port_stop: st.port_stop,
            // far-field runs: the box flux at every far-field frequency must also have settled
            gather_metric: far.as_ref().map(|c| {
                let (plan, grid) = (c.plan.clone(), p.grid.clone());
                crate::fdtd::GatherMetric(std::sync::Arc::new(move |acc: &[[Vec<f32>; 6]]| acc.iter().map(|d| crate::nf2ff::box_flux(&plan, &grid, d)).collect()))
            }),
            gather_tol: far.as_ref().map(|c| c.opts.tol).unwrap_or(0.01),
            gather_lag: 5,
            objective: ro.objective.as_ref().map(|o| o.metric(j, &p, &exc, far.as_ref().map(|c| (&c.plan, c.opts.freqs.as_slice())))),
            objective_stop: obj_stop,
            objective_tol_db: ro.objective.as_ref().map(|o| o.tol_db).unwrap_or(0.0),
            objective_window: obj_window,
            min_steps: obj_min_steps,
            extrap: ro.objective.as_ref().and_then(|o| o.extrapolate.clone()).map(|x| crate::extrap::Opts { dt: p.dt, decimate: dft_every, ..x }),
        };
        let r = backend(&p, &spec);
        eprintln!("fieldsolve: port {} excited: {} steps in {:.1} s (set-up {:.2} s) on {} ({:.0} Mcell/s)", names[j], r.steps, r.seconds, r.setup_seconds, r.backend, p.n() as f64 * r.steps as f64 / r.seconds * 1e-6);
        if let Some((rate, q)) = crate::fdtd::energy_decay(&r.energy, p.dt, st.fmax / 3.0) {
            eprintln!(
                "fieldsolve:   late-time field energy {:+.2} dB/ns{}",
                rate,
                if rate >= 0.0 { " (GROWING: numerical instability, not ringing)".to_string() } else { format!(" (Q {:.0} if at {:.2} GHz)", q, st.fmax / 3e9) }
            );
        }
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
        if let Some(c) = far.as_mut() {
            let ff = c.opts.freqs.clone();
            let (mut fv, mut fi) = (Vec::new(), Vec::new());
            for k in 0..np {
                let src: &[f32] = if k == j { &exc } else { &[] };
                let (v, i) = port_spectra(&r.v[k], src, p.ports[k].r, p.dt, &ff);
                fv.push(v);
                fi.push(i);
            }
            let t1 = std::time::Instant::now();
            c.add_run(j, &p.grid, &r.gathered, fv, fi);
            eprintln!("fieldsolve: far-field transform for {}: {:.1} s", names[j], t1.elapsed().as_secs_f64());
        }
        vall[j] = vs;
        iall[j] = is;
        if let Some(fd) = &r.fields {
            crate::export::write_fields(out, &p, &st, fd).map_err(|e| e.to_string())?;
        }
        if !r.snaps.is_empty() {
            let names: Vec<String> = snap_sheets.iter().map(|s| s.name.clone()).collect();
            crate::export::write_frames(out, &p, &st, &names, &r.snaps, &exc).map_err(|e| e.to_string())?;
        }
        let spread = crate::fdtd::objective_spread(&r.objective, obj_window);
        let tail = crate::fdtd::objective_tail(&r.objective, obj_window);
        if r.stop == "max_steps" {
            eprintln!("fieldsolve: WARNING: port {} reached max_steps ({}) before any stop rule: an incomplete estimate, not converged", names[j], r.steps);
        }
        eprintln!("fieldsolve:   stopped by {}{}", r.stop, tail.map(|(l, t)| format!(", objective change {l:.3} dB over the last {:.2} ns, projected remaining {t:.3} dB", obj_window as f64 * p.dt * 1e9)).unwrap_or_default());
        let mut run_json = json!({"port": names[j], "steps": r.steps, "seconds": r.seconds, "setup_seconds": r.setup_seconds, "backend": r.backend,
                                 "stop": r.stop, "physical_time_ns": r.steps as f64 * p.dt * 1e9,
                                 // where the energy/port rules (the spec's own stop) would have ended this run
                                 "standard_rule_step": r.standard_stop.as_ref().map(|x| x.0), "standard_rule": r.standard_stop.as_ref().map(|x| x.1.clone())});
        if let Some(o) = &ro.objective {
            let comps = o.components(&names, j, far.is_some());
            run_json["objective_spread_db"] = json!(spread);
            run_json["objective_last_window_change_db"] = json!(tail.map(|t| t.0));
            run_json["objective_projected_remaining_db"] = json!(tail.map(|t| if t.1.is_finite() { t.1 } else { -1.0 }));
            run_json["objective_window_steps"] = json!(obj_window);
            run_json["objective_components"] = json!(comps);
            run_json["objective_history"] = json!(r.objective.iter().map(|(s, v)| json!([s, v.iter().map(|x| (x * 1e4).round() / 1e4).collect::<Vec<_>>()])).collect::<Vec<_>>());
            if o.extrapolate.is_some() {
                run_json["objective_extrap_history"] = json!(r.objective_extrap.iter().map(|(s, v)| json!([s, v.iter().map(|x| (x * 1e4).round() / 1e4).collect::<Vec<_>>()])).collect::<Vec<_>>());
                run_json["extrapolated"] = json!(r.stop == "extrapolated");
            }
            obj_hist[j] = Some((comps, r.objective.clone()));
        }
        runs.push(run_json);
    }
    let mut zmat_out = Value::Null;
    let mut balance_out = Value::Null;
    if excite.len() == np {
        s_mat = crate::post::s_matrix(&vall, &iall, z0);
        // power balance per port at the low end of the band, where radiation and dielectric loss are
        // negligible: 1 - sum_k |S_kj|^2 should be ~0; a flat few-percent value means a port's edges are
        // not all part of the circuit (or a port sits across something it should not)
        let nlow = (freqs.len() / 10).max(1);
        let bal: Vec<f64> = (0..np).map(|j| (0..nlow).map(|f| 1.0 - (0..np).map(|k| s_mat[k][j][f].abs().powi(2)).sum::<f64>()).sum::<f64>() / nlow as f64).collect();
        for j in 0..np {
            if bal[j].abs() > 0.02 {
                eprintln!("fieldsolve: WARNING: port {}: power balance {:+.1} % over {:.3e}..{:.3e} Hz (should be ~0): check the port's placement", names[j], 100.0 * bal[j], freqs[0], freqs[nlow - 1]);
            }
        }
        balance_out = json!(names.iter().zip(&bal).map(|(n, b)| json!({"port": n, "low_band_loss": b})).collect::<Vec<_>>());
        write_touchstone(&out.join(format!("sparams.s{np}p")), freqs, &s_mat, z0).map_err(|e| e.to_string())?;
        // impedance matrix per frequency: z[i][j] = {re: [...], im: [...]}
        let zs: Vec<Option<Vec<Vec<Cx>>>> = (0..freqs.len()).map(|f| z_from_s(&(0..np).map(|i| (0..np).map(|j| s_mat[i][j][f]).collect()).collect::<Vec<_>>(), z0)).collect();
        zmat_out = json!((0..np).map(|i| (0..np).map(|j| json!({
            "re": zs.iter().map(|z| z.as_ref().map(|z| z[i][j].re).unwrap_or(f64::NAN)).collect::<Vec<_>>(),
            "im": zs.iter().map(|z| z.as_ref().map(|z| z[i][j].im).unwrap_or(f64::NAN)).collect::<Vec<_>>(),
        })).collect::<Vec<_>>()).collect::<Vec<_>>());
    }
    let lband = ro.l_band;
    let mut ports_out = Vec::new();
    for j in 0..np {
        let Some(z) = &zself[j] else { continue };
        let mut o = json!({
            "name": names[j],
            "freqs": freqs,
            "z_re": z.iter().map(|c| c.re).collect::<Vec<_>>(),
            "z_im": z.iter().map(|c| c.im).collect::<Vec<_>>(),
            "s11_db": s_mat[j][j].iter().map(|c| 20.0 * c.abs().max(1e-30).log10()).collect::<Vec<_>>(),
            // |S_kj| in dB for every port k with this port excited (other ports terminated in their R)
            "s_col_db": (0..np).map(|k| s_mat[k][j].iter().map(|c| 20.0 * c.abs().max(1e-30).log10()).collect::<Vec<_>>()).collect::<Vec<_>>(),
        });
        let (f0, f1) = lband.unwrap_or((freqs[0], freqs[freqs.len() / 10].max(freqs[0])));
        let (r_fit, l_fit) = fit_rl(freqs, z, f0, f1);
        o["l_fit"] = json!({"band": [f0, f1], "r_ohm": r_fit, "l_h": l_fit});
        eprintln!("fieldsolve: {}: R = {:.4} ohm, L = {:.4} nH over {:.3e}..{:.3e} Hz", names[j], r_fit, l_fit * 1e9, f0, f1);
        ports_out.push(o);
    }
    let far_summary = match &far {
        Some(c) => c.finish(&names, z0, out)?,
        None => Value::Null,
    };
    // the objective, from the far-field summary of its port (radiated power per watt incident, other ports
    // matched) or, without a far field, the accepted power
    let objective_out = match &ro.objective {
        Some(o) => {
            let port = o.port.clone().unwrap_or_else(|| names[excite[0]].clone());
            let pi = names.iter().position(|n| *n == port).ok_or_else(|| format!("objective port {port:?} is not a port"))?;
            let db = |x: f64| 10.0 * x.max(1e-30).log10();
            let mut rows = Vec::new();
            for &f in &o.band {
                let row = far_summary.as_array().and_then(|a| a.iter().find(|r| r["port"].as_str() == Some(port.as_str()) && (r["f"].as_f64().unwrap_or(0.0) - f).abs() <= 1e-6 * f));
                let fi = (0..freqs.len()).min_by(|&a, &b| (freqs[a] - f).abs().partial_cmp(&(freqs[b] - f).abs()).unwrap()).unwrap_or(0);
                let (p_in, p_rad) = match row {
                    Some(r) => (r["p_in"].as_f64().unwrap_or(f64::NAN), r["p_rad"].as_f64()),
                    None => (1.0 - s_mat[pi][pi][fi].abs().powi(2), None),
                };
                let total = p_rad.unwrap_or(p_in);
                rows.push(json!({"f": f, "s11_db": db(1.0 - p_in), "accepted_db": db(p_in), "p_rad": p_rad, "total_db": db(total),
                                 "radiated_over_accepted_db": p_rad.map(|x| db(x) - db(p_in))}));
            }
            let worst = rows.iter().map(|r| r["total_db"].as_f64().unwrap_or(f64::NAN)).fold(f64::INFINITY, f64::min);
            let stab = obj_hist[pi].as_ref().and_then(|(_, h)| crate::fdtd::objective_tail(h, obj_window));
            json!({"port": port, "band": o.band, "quantity": if far.is_some() { "radiated power per watt incident (other ports matched)" } else { "accepted power per watt incident" },
                   "per_f": rows, "worst_db": worst,
                   // stability over the objective port's run: change over the last window and the projected
                   // remaining change (-1: not yet decaying geometrically)
                   "stability_db": stab.map(|t| t.0), "projected_remaining_db": stab.map(|t| if t.1.is_finite() { t.1 } else { -1.0 }),
                   "stop_rule": if o.stop { "objective" } else { "energy/ports" },
                   "tol_db": o.tol_db, "window_ns": obj_window as f64 * p.dt * 1e9, "min_steps": obj_min_steps})
        }
        None => Value::Null,
    };
    let ax = [&p.grid.x, &p.grid.y, &p.grid.z];
    let mesh = json!({"dims": [p.nx, p.ny, p.nz], "cells": p.n(), "dt": p.dt,
                      "min_step_mm": (0..3).map(|a| ax[a].windows(2).map(|w| w[1] - w[0]).fold(f64::INFINITY, f64::min) * 1e3).collect::<Vec<_>>(),
                      "refine": st.refined.iter().map(|(b, h, w)| json!({"box": b, "cell": h, "what": w})).collect::<Vec<_>>(),
                      "excitations": excite.len(), "stops": runs.iter().map(|r| r["stop"].clone()).collect::<Vec<_>>()});
    let res = json!({
        "grid": [p.nx, p.ny, p.nz], "cells": p.n(), "dt": p.dt, "runs": runs, "ports": ports_out, "mesh": mesh, "objective": objective_out, "fidelity": ro.fidelity,
        "port_names": names, "port_bars": st.bars.iter().map(|b| json!({"name": b.0, "gap_mm": b.1, "width_mm": b.2, "l_gap_nh": b.3})).collect::<Vec<_>>(), "z_matrix": zmat_out, "mesh_shorts": bridge_lines, "farfield": far_summary, "port_power_balance": balance_out,
    });
    let mut res = res;
    res["settings"] = ro.settings.clone();
    res["checks"] = json!(st.checks);
    res["fidelity"] = crate::fidelity::finish_record(&ro.fidelity, &res);
    if let Some(f) = res["fidelity"].as_object() {
        eprintln!("fieldsolve: fidelity {} ({} cells, dt {:.3e} s, stops {:?}), objective worst {:.2} dB, stability {:?} dB",
            f["status"].as_str().unwrap_or("?"), res["cells"], res["dt"].as_f64().unwrap_or(0.0), f["actual"]["stops"], res.pointer("/objective/worst_db").and_then(|v| v.as_f64()).unwrap_or(f64::NAN), res.pointer("/objective/stability_db"));
    }
    std::fs::write(out.join("result.json"), serde_json::to_string_pretty(&res).unwrap()).map_err(|e| e.to_string())?;
    Ok(res)
}

pub fn load(spec_path: &Path) -> Result<(Ctx, PathBuf), String> {
    let txt = std::fs::read_to_string(spec_path).map_err(|e| format!("{}: {e}", spec_path.display()))?;
    let spec: Value = serde_json::from_str(&txt).map_err(|e| format!("{}: {e}", spec_path.display()))?;
    // "preset": "rf" | "antenna" fills the recommended model settings the spec leaves out (presets.rs)
    let spec = crate::presets::apply(spec)?;
    let base = spec_path.parent().unwrap_or(Path::new("."));
    let gdir = base.join(spec.get("gerbers").and_then(|v| v.as_str()).ok_or("spec needs \"gerbers\"")?);
    let out = base.join(spec.get("out").and_then(|v| v.as_str()).unwrap_or("fieldsolve-out"));
    let mut board = crate::gerber::load_dir(&gdir)?;
    let kicad = spec.get("coords").and_then(|v| v.as_str()).unwrap_or("kicad") == "kicad";
    apply_edits(&mut board, &spec, kicad)?;
    Ok((Ctx { spec, board, kicad, zs: vec![] }, out))
}

/// `"edits": [{"layer": "F.Cu", "op": "add" | "clear", "rect": [x0, y0, x1, y1], "net": "GND"}, ...]`:
/// geometry changes drawn over the Gerber copper in order (a later object wins, as in Gerber), for
/// geometry candidates such as a trimmed radiator, a moved feed or a changed ground clearance.
pub fn apply_edits(board: &mut Board, spec: &Value, kicad: bool) -> Result<(), String> {
    for e in spec.get("edits").and_then(|v| v.as_array()).cloned().unwrap_or_default() {
        let layer = e.get("layer").and_then(|v| v.as_str()).ok_or("edit needs \"layer\"")?;
        let dark = match e.get("op").and_then(|v| v.as_str()) {
            Some("add") => true,
            Some("clear") => false,
            o => return Err(format!("edit op {o:?}: use \"add\" or \"clear\"")),
        };
        let r: Vec<f64> = e.get("rect").and_then(|v| v.as_array()).map(|a| a.iter().filter_map(|x| x.as_f64()).collect()).unwrap_or_default();
        if r.len() != 4 {
            return Err("edit needs \"rect\": [x0, y0, x1, y1]".into());
        }
        let y = |v: f64| if kicad { -v } else { v };
        let (x0, x1, y0, y1) = (r[0].min(r[2]), r[0].max(r[2]), y(r[1]).min(y(r[3])), y(r[1]).max(y(r[3])));
        if x1 - x0 < 1e-9 || y1 - y0 < 1e-9 {
            continue; // empty: lets a parameterised edit vanish (fieldsolve optimise)
        }
        let l = board.layers.iter_mut().find(|l| l.name == layer).ok_or_else(|| format!("edit: no copper layer {layer}"))?;
        let shape = Shape::Polygon(vec![P2::new(x0, y0), P2::new(x1, y0), P2::new(x1, y1), P2::new(x0, y1)]);
        if dark {
            l.bbox.union(&shape.bbox());
        }
        l.prims.push(Prim { shape, dark, net: if dark { e.get("net").and_then(|v| v.as_str()).map(String::from) } else { None }, pad: None });
    }
    Ok(())
}

/// TDR comparison: `tdr.json` and `tdr.svg` (impedance versus distance from the port, two-way time scaled
/// by v = c / sqrt(eps_eff)); prints the impedance excursion of each curve.
pub fn write_tdr(out: &Path, curves: &[(bool, Vec<f64>, Vec<f64>, Vec<f64>, f64)], eps_eff: f64) -> std::io::Result<()> {
    let v = C0 / eps_eff.sqrt();
    let mut js = Vec::new();
    let (w, h, ml, mr, mt, mb) = (900.0, 420.0, 70.0, 20.0, 30.0, 50.0);
    let dmax = 40.0;
    let (zlo, zhi) = (40.0, 120.0);
    let px = |d: f64| ml + (d / dmax) * (w - ml - mr);
    let py = |z: f64| mt + (1.0 - (z - zlo) / (zhi - zlo)) * (h - mt - mb);
    let mut svg = format!("<svg xmlns='http://www.w3.org/2000/svg' width='{w}' height='{h}' font-family='sans-serif' font-size='12'><rect width='100%' height='100%' fill='white'/>");
    for zt in [40, 50, 60, 70, 80, 90, 100, 110, 120] {
        let y = py(zt as f64);
        svg += &format!("<line x1='{ml}' x2='{}' y1='{y}' y2='{y}' stroke='#ddd'/><text x='{}' y='{}' text-anchor='end'>{zt}</text>", w - mr, ml - 6.0, y + 4.0);
    }
    for d in (0..=40).step_by(5) {
        let x = px(d as f64);
        svg += &format!("<line x1='{x}' x2='{x}' y1='{mt}' y2='{}' stroke='#eee'/><text x='{x}' y='{}' text-anchor='middle'>{d}</text>", h - mb, h - mb + 16.0);
    }
    svg += &format!("<text x='{}' y='{}' text-anchor='middle'>distance from port [mm] (v = c / sqrt({eps_eff:.3}))</text>", (ml + w - mr) / 2.0, h - 10.0);
    svg += &format!("<text x='16' y='{}' transform='rotate(-90 16 {})' text-anchor='middle'>TDR impedance [ohm]</text>", h / 2.0, h / 2.0);
    for (slot, t, rho, z, tr) in curves {
        let col = if *slot { "#d62728" } else { "#1f77b4" };
        let pts: Vec<String> = t.iter().zip(z).filter_map(|(&t, &z)| {
            let d = v * t / 2.0 * 1e3;
            (d <= dmax).then(|| format!("{:.1},{:.1}", px(d), py(z.clamp(zlo, zhi))))
        }).collect();
        svg += &format!("<polyline fill='none' stroke='{col}' stroke-width='2' points='{}'/>", pts.join(" "));
        // excursion between the port region and the far termination
        let (mut zmin, mut zmax, mut dmaxz) = (f64::INFINITY, f64::NEG_INFINITY, 0.0);
        for (&t, &zz) in t.iter().zip(z) {
            let d = v * t / 2.0 * 1e3;
            if d > 3.0 && d < 27.0 {
                if zz > zmax {
                    zmax = zz;
                    dmaxz = d;
                }
                zmin = zmin.min(zz);
            }
        }
        let label = if *slot { "slot in ground plane" } else { "solid ground plane" };
        eprintln!("fieldsolve gap: {label}: TDR {zmin:.1} .. {zmax:.1} ohm along the line (peak at {dmaxz:.1} mm), rise time {:.0} ps", tr * 1e12);
        js.push(serde_json::json!({"slot": slot, "t": t, "rho": rho, "z": z, "rise_s": tr, "z_peak": zmax, "d_peak_mm": dmaxz}));
    }
    svg += &format!("<rect x='{}' y='{}' width='12' height='3' fill='#1f77b4'/><text x='{}' y='{}'>solid ground plane</text>", w - 230.0, mt + 10.0, w - 212.0, mt + 14.0);
    svg += &format!("<rect x='{}' y='{}' width='12' height='3' fill='#d62728'/><text x='{}' y='{}'>1 x 10 mm slot at 15 mm</text>", w - 230.0, mt + 28.0, w - 212.0, mt + 32.0);
    svg += "</svg>";
    std::fs::write(out.join("tdr.svg"), svg)?;
    std::fs::write(out.join("tdr.json"), serde_json::to_string(&js).unwrap())
}

/// ASCII raster of the model on the mesh plane nearest `z` (mm) inside the window `a`..`b` (mm, internal
/// coordinates). Per node: '#' x and y metal edges leave the node, '-' only x, '|' only y, 'P' a lumped
/// element edge (port, resistor, short), '.' none. Rows run top to bottom in the board's view.
pub fn raster_ascii(p: &Problem, st: &Setup, z: f64, a: P2, b: P2, kicad: bool) -> String {
    use std::fmt::Write as _;
    let mm = |v: &[f64]| v.iter().map(|x| x * 1e3).collect::<Vec<_>>();
    let (xs, ys, zs) = (mm(&p.grid.x), mm(&p.grid.y), mm(&p.grid.z));
    let k = crate::mesh::snap(&zs, z);
    let (x0, x1, y0, y1) = (a.x.min(b.x), a.x.max(b.x), a.y.min(b.y), a.y.max(b.y));
    let is: Vec<usize> = (0..p.nx).filter(|&i| xs[i] >= x0 && xs[i] <= x1).collect();
    let mut js: Vec<usize> = (0..p.ny).filter(|&j| ys[j] >= y0 && ys[j] <= y1).collect();
    // top row first: the largest internal y (KiCad's internal y is its negated y, so also its top)
    js.reverse();
    let lumped: std::collections::HashSet<(u8, u32)> = p.ports.iter().flat_map(|q| q.src.iter().map(|s| (s.0, s.1))).collect();
    let mut short: std::collections::HashSet<(u8, u32)> = Default::default();
    for l in &st.lumped {
        if l.r == 0.0 && !l.port {
            // shorts are metal edges already; mark them through their box
            for (i, j) in is.iter().flat_map(|&i| js.iter().map(move |&j| (i, j))) {
                let (xm, ym) = (xs[i], ys[j]);
                if xm >= l.lo[0] - 1e-9 && xm <= l.hi[0] + 1e-9 && ym >= l.lo[1] - 1e-9 && ym <= l.hi[1] + 1e-9 {
                    short.insert((l.axis as u8, ((k * p.ny + j) * p.nx + i) as u32));
                }
            }
        }
    }
    let mut o = String::new();
    let _ = writeln!(o, "plane z = {:.3} mm (k = {k}); x {:.3}..{:.3}, {} nodes per row", zs[k], xs[is[0]], xs[*is.last().unwrap()], is.len());
    for &j in &js {
        let yk = if kicad { -ys[j] } else { ys[j] };
        let _ = write!(o, "{yk:9.3} ");
        for &i in &is {
            let id = ((k * p.ny + j) * p.nx + i) as u32;
            let mx = p.cb[0][id as usize] == 0.0 && p.ca[0][id as usize] == 0.0;
            let my = p.cb[1][id as usize] == 0.0 && p.ca[1][id as usize] == 0.0;
            let lp = lumped.contains(&(0, id)) || lumped.contains(&(1, id)) || short.contains(&(0, id)) || short.contains(&(1, id));
            o.push(if lp { 'P' } else if mx && my { '#' } else if mx { '-' } else if my { '|' } else { '.' });
        }
        o.push('\n');
    }
    o
}
