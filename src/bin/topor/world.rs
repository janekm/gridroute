//! Topological routing state: per-layer constrained triangulations of the obstacles, ordered wire
//! crossings on every triangulation edge, rubber-band geometry realised from that topology, a spatial
//! index for exact clearance checks, and an undo journal for transactional commits.

use crate::cdt::{Cdt, Change, NONE};
use crate::funnel::{funnel, polyline, Disc};
use crate::geom::*;
use crate::model::{Board, Kind, ObjSpec};
use std::collections::HashMap;

/// Radius kept around a net's own vertices so wires never pass through them.
pub const EPS_OWN: f64 = 0.0005;
/// Arc discretisation tolerance (vertices sit this far outside the clearance circle).
pub const TOL: f64 = 0.0007;
/// Added to every clearance while routing; the exact check uses the true rules.
pub const MARGIN: f64 = 0.003;

pub const T_BLOCK: u8 = 1;

#[derive(Clone)]
pub struct Obj {
    pub s: ObjSpec,
    pub alive: bool,
    /// Through copper (via, through-hole pad) of a plane net that joins its plane where it stands.
    pub plane: bool,
    /// Vertex id per (layer, point): lv[layer * npts + k].
    pub lv: Vec<u32>,
}

impl Obj {
    pub fn npts(&self) -> usize {
        self.s.pts.len()
    }
    pub fn on(&self, l: usize) -> bool {
        self.s.layers & (1 << l) != 0
    }
    pub fn v(&self, l: usize, k: usize) -> u32 {
        self.lv[l * self.npts() + k]
    }
    pub fn is_poly(&self) -> bool {
        self.s.pts.len() >= 3 && self.s.kind == Kind::Pad
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Term {
    pub obj: u32,
    /// Layer vertex the wire ends at, or NONE for polygon pads (wire ends at the pad centre).
    pub vtx: u32,
}

#[derive(Clone)]
pub struct Crossing {
    pub piece: u32,
    pub edge: u32,
    pub left: u32,
    pub alive: bool,
    /// Index of this crossing in its piece's crossing list.
    pub idx: u32,
}

#[derive(Clone)]
pub struct Piece {
    pub net: u32,
    pub layer: u8,
    pub width: f64,
    pub a: Term,
    pub b: Term,
    pub cx: Vec<u32>,
    pub poly: Vec<P2>,
    pub cparam: Vec<f64>,
    pub alive: bool,
    /// Edge the piece lies along (vertex-to-vertex pieces with no crossings).
    pub along: u32,
    pub overlap: bool,
    /// Where realisation failed (overlapping clearance discs), if it did.
    pub problem: Option<P2>,
    /// Geometry fixed (materialised ribbon members): no topology, never re-realised.
    pub frozen: bool,
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Item {
    pub kind: u8,
    pub id: u32,
    pub sub: u32,
}
pub const IT_OBJ: u8 = 0;
pub const IT_RING: u8 = 1;
pub const IT_SEG: u8 = 2;

#[derive(Clone)]
pub struct Grid {
    pub cell: f64,
    pub nx: usize,
    pub ny: usize,
    pub ox: f64,
    pub oy: f64,
    /// Items per cell, each with the first cell (x, y) of its box: an item spanning several cells is
    /// reported only from the first cell of its box that lies in the query (no sort, no duplicates).
    pub cells: Vec<Vec<(Item, u16, u16)>>,
    /// The wire segments of each cell alone, with their geometry, in the order they have in `cells`: the
    /// clearance scans along a stretch read these without going through the pieces.
    pub segs: Vec<Vec<SegE>>,
}

/// A wire segment as the grid keeps it for the clearance scans.
#[derive(Clone, Copy, Debug)]
pub struct SegE {
    pub id: u32,
    pub sub: u32,
    pub a: P2,
    pub b: P2,
    pub net: u32,
    pub width: f64,
}

impl Grid {
    fn new(lo: P2, hi: P2, cell: f64) -> Grid {
        let nx = (((hi.x - lo.x) / cell).ceil() as usize).max(1);
        let ny = (((hi.y - lo.y) / cell).ceil() as usize).max(1);
        Grid { cell, nx, ny, ox: lo.x, oy: lo.y, cells: vec![Vec::new(); nx * ny], segs: vec![Vec::new(); nx * ny] }
    }
    #[inline]
    pub fn range(&self, lo: P2, hi: P2) -> (usize, usize, usize, usize) {
        let c = |v: f64, o: f64, n: usize| -> usize { (((v - o) / self.cell).floor().max(0.0) as usize).min(n - 1) };
        (c(lo.x, self.ox, self.nx), c(lo.y, self.oy, self.ny), c(hi.x, self.ox, self.nx), c(hi.y, self.oy, self.ny))
    }
    pub fn add(&mut self, lo: P2, hi: P2, it: Item) {
        let (x0, y0, x1, y1) = self.range(lo, hi);
        for y in y0..=y1 {
            for x in x0..=x1 {
                self.cells[y * self.nx + x].push((it, x0 as u16, y0 as u16));
            }
        }
    }
    /// A wire segment (part): in `cells` and, with its geometry, in `segs`.
    pub fn add_seg(&mut self, lo: P2, hi: P2, it: Item, e: SegE) {
        let (x0, y0, x1, y1) = self.range(lo, hi);
        for y in y0..=y1 {
            for x in x0..=x1 {
                self.cells[y * self.nx + x].push((it, x0 as u16, y0 as u16));
                self.segs[y * self.nx + x].push(e);
            }
        }
    }
    pub fn remove(&mut self, lo: P2, hi: P2, it: Item) {
        let (x0, y0, x1, y1) = self.range(lo, hi);
        for y in y0..=y1 {
            for x in x0..=x1 {
                let k = y * self.nx + x;
                let c = &mut self.cells[k];
                if let Some(i) = c.iter().position(|q| q.0 == it) {
                    // `segs` keeps the order the segments have in `cells`: the swap there moves the last item
                    // into the gap, which is the last segment if it is one; otherwise the segments stay in order
                    if it.kind == IT_SEG {
                        let last_seg = c.last().map_or(false, |q| q.0.kind == IT_SEG);
                        let sg = &mut self.segs[k];
                        if let Some(j) = sg.iter().position(|e| e.id == it.id && e.sub == it.sub) {
                            if last_seg {
                                sg.swap_remove(j);
                            } else {
                                sg.remove(j);
                            }
                        }
                    }
                    c.swap_remove(i);
                }
            }
        }
    }
    /// Visit each item whose box meets the query box exactly once.
    #[inline]
    pub fn query_each<F: FnMut(Item)>(&self, lo: P2, hi: P2, mut f: F) {
        let (x0, y0, x1, y1) = self.range(lo, hi);
        for y in y0..=y1 {
            for x in x0..=x1 {
                for &(it, ix, iy) in &self.cells[y * self.nx + x] {
                    if x == (ix as usize).max(x0) && y == (iy as usize).max(y0) {
                        f(it);
                    }
                }
            }
        }
    }
    pub fn query(&self, lo: P2, hi: P2, out: &mut Vec<Item>) {
        out.clear();
        self.query_each(lo, hi, |it| out.push(it));
    }
}

#[derive(Clone)]
pub struct Layer {
    pub cdt: Cdt,
    pub vobj: Vec<u32>,
    pub vk: Vec<u16>,
    pub cross: Vec<Vec<u32>>,
    pub along: Vec<u32>,
    /// Per edge: triangle with the edge directed a->b, and b->a.
    pub etri: Vec<[u32; 2]>,
    pub tflag: Vec<u8>,
    pub grid: Grid,
    /// Recently located triangle per coarse cell, as a starting point for point location.
    pub hint: Vec<std::cell::Cell<u32>>,
}

#[derive(Default, Clone)]
struct LJ {
    nvobj: usize,
    nedge: usize,
    ntri: usize,
    cross: Vec<(u32, Vec<u32>)>,
    along: Vec<(u32, u32)>,
    etri: Vec<(u32, [u32; 2])>,
}

#[derive(Default, Clone)]
struct WJ {
    nobj: usize,
    npiece: usize,
    ncross: usize,
    layers: Vec<LJ>,
    pieces: Vec<(u32, Piece)>,
    crosses: Vec<(u32, Crossing)>,
    objs: Vec<(u32, Obj)>,
    obj_pieces: Vec<(u32, Vec<u32>)>,
}

#[derive(Clone, Copy, Debug)]
pub struct Attr {
    pub kind: Kind,
    pub net: u32,
    pub class: u32,
    pub clear: f64,
}

#[derive(Clone, Debug)]
pub struct Viol {
    pub what: String,
    pub layer: i32,
    pub at: P2,
    pub dist: f64,
    pub need: f64,
    pub pieces: Vec<u32>,
    pub objs: Vec<u32>,
}

#[derive(Clone)]
pub struct World {
    pub b: Board,
    pub objs: Vec<Obj>,
    pub layers: Vec<Layer>,
    pub crosses: Vec<Crossing>,
    pub pieces: Vec<Piece>,
    pub obj_pieces: Vec<Vec<u32>>,
    pub dirty: Vec<u32>,
    jr: Option<WJ>,
    pub max_clear: f64,
    /// Widest wire half-width of any net (segments are indexed by their centre line).
    pub max_hw: f64,
    /// Largest wire-to-wire clearance (classes, class pairs, pair gaps), margin included.
    pub max_wire_clear: f64,
    /// Wires taken out because a triangulation change could not re-thread them consistently.
    pub thread_drops: usize,
    pub warnings: Vec<String>,
    /// Escape lanes in front of fine-pitch chip pads: (net, from, to, half-width); no other net's via
    /// may sit in one (it would seal the pad in). Bucketed by 1 mm cell.
    pub escapes: Vec<(u32, P2, P2, f64, u32)>,
    pub esc_grid: crate::fx::FxMap<(i32, i32), Vec<u32>>,
    pub use_escapes: bool,
    /// Object id for each entry of Board::objs (after merging overlaps).
    pub spec_obj: Vec<u32>,
}

fn tri_edge_index(t: &crate::cdt::Tri, e: u32) -> usize {
    t.e.iter().position(|&x| x == e).unwrap()
}

impl World {
    pub fn nsig(&self) -> usize {
        self.layers.len()
    }

    fn hint_cell(&self, l: usize, p: P2) -> usize {
        let g = &self.layers[l].grid;
        let cx = (((p.x - g.ox) / (2.0 * g.cell)).floor().max(0.0) as usize).min(g.nx / 2);
        let cy = (((p.y - g.oy) / (2.0 * g.cell)).floor().max(0.0) as usize).min(g.ny / 2);
        (cy * (g.nx / 2 + 1) + cx) % self.layers[l].hint.len()
    }

    /// Point location with a per-cell hint (much shorter walks than from an arbitrary triangle).
    pub fn locate(&self, l: usize, p: P2) -> Result<crate::cdt::Loc, String> {
        let lay = &self.layers[l];
        let hc = self.hint_cell(l, p);
        let h = lay.hint[hc].get();
        let start = if h != NONE && (h as usize) < lay.cdt.tris.len() && lay.cdt.tris[h as usize].alive { h } else { lay.cdt.any_alive() };
        let r = lay.cdt.locate(to_ip(p), start)?;
        let t = match r {
            crate::cdt::Loc::Inside(t) | crate::cdt::Loc::OnEdge(t, _) => t,
            crate::cdt::Loc::OnVertex(v) => lay.cdt.vtri[v as usize],
        };
        lay.hint[hc].set(t);
        Ok(r)
    }

    pub fn vpos(&self, l: usize, v: u32) -> P2 {
        from_ip(self.layers[l].cdt.pts[v as usize])
    }

    // ------------------------------------------------------------------ build
    pub fn new(b: Board) -> Result<World, String> {
        let nsig = b.nsig();
        let mut objs: Vec<Obj> = Vec::new();
        // outline rings, densified
        let mut rings: Vec<Vec<P2>> = vec![b.outline.clone()];
        rings.extend(b.outline_holes.iter().cloned());
        for (ri, ring) in rings.into_iter().enumerate() {
            let mut pts = Vec::new();
            let n = ring.len();
            for i in 0..n {
                let (a, c) = (ring[i], ring[(i + 1) % n]);
                let k = ((a.dist(c) / 1.5).ceil() as usize).max(1);
                for j in 0..k {
                    pts.push(a.lerp(c, j as f64 / k as f64));
                }
            }
            let center = pts.iter().fold(p2(0., 0.), |s, &q| s + q) * (1.0 / pts.len() as f64);
            objs.push(Obj {
                s: ObjSpec { kind: Kind::Outline, net: NONE, layers: (1 << nsig) - 1, pts, r: 0., clear: b.edge_clear, center, drill: 0., name: if ri == 0 { "outline".into() } else { "outline_hole".into() }, through: false },
                alive: true,
                plane: false,
                lv: Vec::new(),
            });
        }
        let mut specs = b.objs.clone();
        // densify long keepout rings so triangles near them stay small
        for s in specs.iter_mut() {
            if s.kind == Kind::Keepout || s.kind == Kind::Outline {
                let ring = s.pts.clone();
                let mut pts = Vec::new();
                let n = ring.len();
                for i in 0..n {
                    let (a, c) = (ring[i], ring[(i + 1) % n]);
                    let k = ((a.dist(c) / 1.5).ceil() as usize).max(1);
                    for j in 0..k {
                        pts.push(a.lerp(c, j as f64 / k as f64));
                    }
                }
                s.pts = pts;
            }
        }
        let first = objs.len() as u32;
        let (merged, map, mut warnings) = merge_overlaps(specs, nsig);
        for s in merged {
            // through copper joins its plane in any region of the net; an SMD pad only in one on its layer
            // (an outer-layer pour that is the net's plane, as for a module's ground)
            let smd_layer = (0..b.signal.len()).find(|&l| s.layers == 1 << l).map(|l| b.signal[l]);
            let plane = s.net != NONE && b.nets[s.net as usize].plane.is_some() && s.kind == Kind::Pad && !fixed_copper(&s.name) && if s.through { pad_in_plane(&b, &s) } else { smd_layer.map_or(false, |l| b.in_plane_on(s.net, s.center, 0.0, Some(l))) };
            objs.push(Obj { s, alive: true, plane, lv: Vec::new() });
        }
        let spec_obj: Vec<u32> = map.iter().map(|&k| first + k as u32).collect();
        let lo = p2(-5.0, -5.0);
        let hi = p2(b.size[0] + 5.0, b.size[1] + 5.0);
        let mut layers = Vec::new();
        for o in objs.iter_mut() {
            o.lv = vec![NONE; o.npts() * nsig];
        }
        for l in 0..nsig {
            let mut cdt = Cdt::new(to_ip(lo), to_ip(hi));
            let mut vobj = vec![NONE; 4];
            let mut vk = vec![0u16; 4];
            for (oi, o) in objs.iter_mut().enumerate() {
                if !o.on(l) {
                    continue;
                }
                for k in 0..o.npts() {
                    let p = to_ip(o.s.pts[k]);
                    let v = cdt.insert(p)?;
                    if (v as usize) < vobj.len() {
                        if vobj[v as usize] != oi as u32 {
                            warnings.push(format!("layer {}: point of {} coincides with {}", l, o.s.name, vobj[v as usize]));
                        }
                    } else {
                        vobj.push(oi as u32);
                        vk.push(k as u16);
                    }
                    o.lv[l * o.s.pts.len() + k] = v;
                }
            }
            for o in objs.iter() {
                if !o.on(l) {
                    continue;
                }
                let n = o.npts();
                if n < 2 {
                    continue;
                }
                let m = if n == 2 { 1 } else { n };
                for k in 0..m {
                    let (a, c) = (o.v(l, k), o.v(l, (k + 1) % n));
                    if let Err(e) = cdt.insert_fixed(a, c) {
                        warnings.push(format!("layer {}: fixed edge of {}: {}", l, o.s.name, e));
                    }
                }
            }
            cdt.validate()?;
            let ne = cdt.edges.len();
            let nt = cdt.tris.len();
            let mut layer = Layer {
                cdt,
                vobj,
                vk,
                cross: vec![Vec::new(); ne],
                along: vec![NONE; ne],
                etri: vec![[NONE, NONE]; ne],
                tflag: vec![0; nt],
                grid: Grid::new(lo, hi, 0.6),
                hint: Vec::new(),
            };
            layer.hint = (0..layer.grid.nx * layer.grid.ny / 4 + 8).map(|_| std::cell::Cell::new(NONE)).collect();
            for t in 0..nt as u32 {
                if layer.cdt.tris[t as usize].alive {
                    set_etri(&mut layer, t);
                }
            }
            layers.push(layer);
        }
        let nobj = objs.len();
        let mut w = World {
            b,
            objs,
            layers,
            crosses: Vec::new(),
            pieces: Vec::new(),
            obj_pieces: vec![Vec::new(); nobj],
            dirty: Vec::new(),
            jr: None,
            max_clear: 0.0,
            escapes: Vec::new(),
            esc_grid: Default::default(),
            use_escapes: std::env::var("TOPOR_ESCAPES").is_ok(),
            max_hw: 0.0,
            max_wire_clear: 0.0,
            thread_drops: 0,
            warnings: Vec::new(),
            spec_obj,
        };
        w.warnings.append(&mut warnings);
        let mut mc = w.b.edge_clear.max(w.b.hole_clear);
        for c in &w.b.classes {
            mc = mc.max(c.clear);
        }
        for &c in w.b.class_pair.iter().chain(&w.b.track_clear) {
            mc = mc.max(c);
        }
        for o in &w.objs {
            mc = mc.max(o.s.clear);
        }
        w.max_clear = mc + MARGIN;
        w.max_hw = w.b.nets.iter().flat_map(|n| n.wpref.iter().cloned()).fold(0.1, f64::max) / 2.0;
        w.max_wire_clear = w.b.classes.iter().map(|c| c.clear).chain(w.b.class_pair.iter().cloned()).chain(w.b.pair_gap.iter().cloned()).chain(w.b.track_clear.iter().cloned()).fold(0.0, f64::max) + MARGIN;
        for l in 0..nsig {
            w.mark_blocked(l);
        }
        w.build_escapes();
        for oi in 0..w.objs.len() {
            w.index_obj(oi as u32, true);
        }
        Ok(w)
    }

    fn mark_blocked(&mut self, l: usize) {
        let outline_ids: Vec<u32> = (0..self.objs.len() as u32).filter(|&i| self.objs[i as usize].s.kind == Kind::Outline && self.objs[i as usize].s.name.starts_with("outline")).collect();
        let layer = &mut self.layers[l];
        let nt = layer.cdt.tris.len();
        // outside the board: flood from triangles touching the bounding-box corners, stopping at outline edges
        let is_outline_edge = |layer: &Layer, e: u32| -> bool {
            let ed = layer.cdt.edges[e as usize];
            ed.fixed && outline_ids.contains(&layer.vobj[ed.a as usize]) && outline_ids.contains(&layer.vobj[ed.b as usize])
        };
        let mut seen = vec![false; nt];
        let mut stack: Vec<u32> = Vec::new();
        for t in 0..nt {
            let tr = &layer.cdt.tris[t];
            if tr.alive && tr.v.iter().any(|&v| v < 4) {
                stack.push(t as u32);
            }
        }
        while let Some(t) = stack.pop() {
            if seen[t as usize] {
                continue;
            }
            seen[t as usize] = true;
            layer.tflag[t as usize] |= T_BLOCK;
            let tr = layer.cdt.tris[t as usize];
            for i in 0..3 {
                if is_outline_edge(layer, tr.e[i]) || tr.n[i] == NONE {
                    continue;
                }
                stack.push(tr.n[i]);
            }
        }
        // interiors of objects (all three vertices on one object, centroid inside it)
        for t in 0..nt {
            let tr = layer.cdt.tris[t];
            if !tr.alive {
                continue;
            }
            let o = layer.vobj[tr.v[0] as usize];
            if o == NONE || layer.vobj[tr.v[1] as usize] != o || layer.vobj[tr.v[2] as usize] != o {
                continue;
            }
            let obj = &self.objs[o as usize];
            if obj.npts() < 3 {
                continue;
            }
            let c = (from_ip(layer.cdt.pts[tr.v[0] as usize]) + from_ip(layer.cdt.pts[tr.v[1] as usize]) + from_ip(layer.cdt.pts[tr.v[2] as usize])) * (1.0 / 3.0);
            let inside = point_in_poly(&obj.s.pts, c);
            match obj.s.kind {
                Kind::Outline if obj.s.name == "outline" => {}
                Kind::Pad | Kind::Hole | Kind::Via if obj.npts() >= 3 => {
                    layer.tflag[t] |= T_BLOCK;
                }
                _ => {
                    if inside {
                        layer.tflag[t] |= T_BLOCK;
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------ journal
    pub fn begin(&mut self) {
        let mut j = WJ { nobj: self.objs.len(), npiece: self.pieces.len(), ncross: self.crosses.len(), ..Default::default() };
        for l in self.layers.iter_mut() {
            l.cdt.begin();
            j.layers.push(LJ { nvobj: l.vobj.len(), nedge: l.cross.len(), ntri: l.tflag.len(), ..Default::default() });
        }
        self.jr = Some(j);
    }
    pub fn commit(&mut self) {
        for l in self.layers.iter_mut() {
            l.cdt.commit();
        }
        self.jr = None;
    }
    pub fn rollback(&mut self) {
        let j = match self.jr.take() {
            Some(j) => j,
            None => return,
        };
        // un-index current geometry of touched pieces and new objects
        let mut touched: Vec<u32> = j.pieces.iter().map(|x| x.0).collect();
        touched.extend(j.npiece as u32..self.pieces.len() as u32);
        touched.sort();
        touched.dedup();
        for &p in &touched {
            self.unindex_piece(p);
        }
        for o in j.nobj..self.objs.len() {
            if self.objs[o].alive {
                self.index_obj(o as u32, false);
            }
        }
        let mut touched_objs: Vec<u32> = j.objs.iter().map(|x| x.0).collect();
        touched_objs.sort();
        touched_objs.dedup();
        for &o in &touched_objs {
            if self.objs[o as usize].alive {
                self.index_obj(o, false);
            }
        }
        for (li, (l, lj)) in self.layers.iter_mut().zip(j.layers.into_iter()).enumerate() {
            let _ = li;
            l.cdt.rollback();
            for (e, v) in lj.cross.into_iter().rev() {
                l.cross[e as usize] = v;
            }
            for (e, v) in lj.along.into_iter().rev() {
                l.along[e as usize] = v;
            }
            for (e, v) in lj.etri.into_iter().rev() {
                l.etri[e as usize] = v;
            }
            l.vobj.truncate(lj.nvobj);
            l.vk.truncate(lj.nvobj);
            l.cross.truncate(lj.nedge);
            l.along.truncate(lj.nedge);
            l.etri.truncate(lj.nedge);
            l.tflag.truncate(lj.ntri);
        }
        for (i, p) in j.pieces.into_iter().rev() {
            self.pieces[i as usize] = p;
        }
        for (i, c) in j.crosses.into_iter().rev() {
            self.crosses[i as usize] = c;
        }
        for (i, o) in j.objs.into_iter().rev() {
            self.objs[i as usize] = o;
        }
        for (i, v) in j.obj_pieces.into_iter().rev() {
            self.obj_pieces[i as usize] = v;
        }
        self.pieces.truncate(j.npiece);
        self.crosses.truncate(j.ncross);
        self.objs.truncate(j.nobj);
        self.obj_pieces.truncate(j.nobj);
        for &p in &touched {
            if (p as usize) < self.pieces.len() {
                self.index_piece(p);
            }
        }
        for &o in &touched_objs {
            if (o as usize) < self.objs.len() && self.objs[o as usize].alive {
                self.index_obj(o, true);
            }
        }
        self.dirty.clear();
    }
    fn save_cross_list(&mut self, l: usize, e: u32) {
        if let Some(j) = &mut self.jr {
            if (e as usize) < j.layers[l].nedge {
                j.layers[l].cross.push((e, self.layers[l].cross[e as usize].clone()));
            }
        }
    }
    fn save_along(&mut self, l: usize, e: u32) {
        if let Some(j) = &mut self.jr {
            if (e as usize) < j.layers[l].nedge {
                j.layers[l].along.push((e, self.layers[l].along[e as usize]));
            }
        }
    }
    fn save_etri(&mut self, l: usize, e: u32) {
        if let Some(j) = &mut self.jr {
            if (e as usize) < j.layers[l].nedge {
                j.layers[l].etri.push((e, self.layers[l].etri[e as usize]));
            }
        }
    }
    fn save_piece(&mut self, p: u32) {
        if let Some(j) = &mut self.jr {
            if (p as usize) < j.npiece && !j.pieces.iter().any(|x| x.0 == p) {
                j.pieces.push((p, self.pieces[p as usize].clone()));
            }
        }
    }
    fn save_cross(&mut self, c: u32) {
        if let Some(j) = &mut self.jr {
            if (c as usize) < j.ncross {
                j.crosses.push((c, self.crosses[c as usize].clone()));
            }
        }
    }
    fn save_obj(&mut self, o: u32) {
        if let Some(j) = &mut self.jr {
            if (o as usize) < j.nobj && !j.objs.iter().any(|x| x.0 == o) {
                j.objs.push((o, self.objs[o as usize].clone()));
            }
        }
    }
    fn save_obj_pieces(&mut self, o: u32) {
        if let Some(j) = &mut self.jr {
            if (o as usize) < j.nobj {
                j.obj_pieces.push((o, self.obj_pieces[o as usize].clone()));
            }
        }
    }

    // ------------------------------------------------------------------ spatial index
    fn obj_box(&self, o: u32) -> (P2, P2) {
        let s = &self.objs[o as usize].s;
        let mut lo = p2(f64::INFINITY, f64::INFINITY);
        let mut hi = p2(f64::NEG_INFINITY, f64::NEG_INFINITY);
        for q in &s.pts {
            lo = p2(lo.x.min(q.x), lo.y.min(q.y));
            hi = p2(hi.x.max(q.x), hi.y.max(q.y));
        }
        (lo - p2(s.r, s.r), hi + p2(s.r, s.r))
    }
    fn index_obj(&mut self, o: u32, add: bool) {
        let ring = matches!(self.objs[o as usize].s.kind, Kind::Outline | Kind::Keepout);
        let n = self.objs[o as usize].npts();
        for l in 0..self.nsig() {
            if !self.objs[o as usize].on(l) {
                continue;
            }
            if ring {
                let cell = self.layers[l].grid.cell;
                for k in 0..n {
                    let (a, c) = (self.objs[o as usize].s.pts[k], self.objs[o as usize].s.pts[(k + 1) % n]);
                    let it = Item { kind: IT_RING, id: o, sub: k as u32 };
                    // in parts, like wire segments (an outline edge would otherwise fill its whole box)
                    let m = ((a.dist(c) / cell).ceil() as usize).clamp(1, 1 << 16);
                    for i in 0..m {
                        let (u, w) = (a.lerp(c, i as f64 / m as f64), a.lerp(c, (i + 1) as f64 / m as f64));
                        let lo = p2(u.x.min(w.x), u.y.min(w.y));
                        let hi = p2(u.x.max(w.x), u.y.max(w.y));
                        if add {
                            self.layers[l].grid.add(lo, hi, it);
                        } else {
                            self.layers[l].grid.remove(lo, hi, it);
                        }
                    }
                }
            } else {
                let (lo, hi) = self.obj_box(o);
                let it = Item { kind: IT_OBJ, id: o, sub: 0 };
                if add {
                    self.layers[l].grid.add(lo, hi, it);
                } else {
                    self.layers[l].grid.remove(lo, hi, it);
                }
            }
        }
    }
    /// Grid entries of a piece: each segment in parts no longer than a cell, so a long diagonal does
    /// not fill its whole bounding box (queries may then see a segment once per part).
    fn seg_items(&self, p: u32) -> Vec<(P2, P2, Item)> {
        let pc = &self.pieces[p as usize];
        let cell = self.layers[pc.layer as usize].grid.cell;
        let mut v = Vec::new();
        for k in 0..pc.poly.len().saturating_sub(1) {
            let (a, c) = (pc.poly[k], pc.poly[k + 1]);
            let m = ((a.dist(c) / cell).ceil() as usize).clamp(1, 4096);
            let it = Item { kind: IT_SEG, id: p, sub: k as u32 };
            for i in 0..m {
                let (u, w) = (a.lerp(c, i as f64 / m as f64), a.lerp(c, (i + 1) as f64 / m as f64));
                v.push((p2(u.x.min(w.x), u.y.min(w.y)), p2(u.x.max(w.x), u.y.max(w.y)), it));
            }
        }
        v
    }
    fn index_piece(&mut self, p: u32) {
        if !self.pieces[p as usize].alive {
            return;
        }
        let l = self.pieces[p as usize].layer as usize;
        for (lo, hi, it) in self.seg_items(p) {
            let pc = &self.pieces[p as usize];
            let e = SegE { id: p, sub: it.sub, a: pc.poly[it.sub as usize], b: pc.poly[it.sub as usize + 1], net: pc.net, width: pc.width };
            self.layers[l].grid.add_seg(lo, hi, it, e);
        }
    }
    fn unindex_piece(&mut self, p: u32) {
        if !self.pieces[p as usize].alive {
            return;
        }
        let l = self.pieces[p as usize].layer as usize;
        for (lo, hi, it) in self.seg_items(p) {
            self.layers[l].grid.remove(lo, hi, it);
        }
    }

    // ------------------------------------------------------------------ rules
    pub fn obj_attr(&self, o: u32) -> Attr {
        let s = &self.objs[o as usize].s;
        let class = if s.net != NONE { self.b.nets[s.net as usize].class as u32 } else { NONE };
        Attr { kind: s.kind, net: s.net, class, clear: s.clear }
    }
    /// A wire of `net` on signal layer `l`.
    pub fn net_attr(&self, net: u32, l: usize) -> Attr {
        let c = self.b.nets[net as usize].class;
        Attr { kind: Kind::Pad, net, class: c as u32, clear: self.b.track_clear(net, l) }
    }
    pub fn vert_attr(&self, l: usize, v: u32) -> (Attr, f64) {
        let o = self.layers[l].vobj[v as usize];
        if o == NONE || !self.objs[o as usize].alive {
            return (Attr { kind: Kind::Keepout, net: NONE, class: NONE, clear: -1.0 }, 0.0);
        }
        (self.obj_attr(o), self.objs[o as usize].s.r)
    }
    /// A net tie's copper belongs to no net for KiCad's hole clearance: a via of a net it joins still keeps
    /// its drill that far off (`d`: distance from the via's centre to the copper, `rd`: drill radius).
    fn tie_hole(&self, o: u32, d: f64, rd: f64) -> bool {
        self.objs[o as usize].s.name == "net_tie" && d - rd < self.b.hole_clear + MARGIN
    }
    /// Required gap between two items, or None when they may touch (same net, or free vertex).
    pub fn clr(&self, x: &Attr, y: &Attr) -> Option<f64> {
        clr_in(&self.b, x, y)
    }

    // ------------------------------------------------------------------ topology helpers
    pub fn edge_len(&self, l: usize, e: u32) -> f64 {
        let ed = self.layers[l].cdt.edges[e as usize];
        self.vpos(l, ed.a).dist(self.vpos(l, ed.b))
    }

    /// Pack elements outward from vertex `v` along an edge: returns centre distance of each wire.
    /// `wires` are (attr, width) in order from v.
    pub fn pack(&self, l: usize, v: u32, wires: &[(Attr, f64)]) -> Vec<f64> {
        let mut out = Vec::with_capacity(wires.len());
        self.pack_each(l, v, wires.iter().cloned(), |_, _, _, c| out.push(c));
        out
    }

    /// Centre distance of the last of `wires` packed outward from vertex v.
    pub fn pack_last<I: Iterator<Item = (Attr, f64)>>(&self, l: usize, v: u32, wires: I) -> f64 {
        let mut last = EPS_OWN;
        self.pack_each(l, v, wires, |_, _, _, c| last = c);
        last
    }

    /// Packing without allocation: `emit(k, attr, width, centre)` for each wire in order from v.
    #[inline]
    fn pack_each<I: Iterator<Item = (Attr, f64)>, F: FnMut(usize, Attr, f64, f64)>(&self, l: usize, v: u32, wires: I, mut emit: F) {
        const N: usize = 48;
        let (va, vr) = self.vert_attr(l, v);
        let mut pa: [Attr; N] = [va; N];
        let mut po: [f64; N] = [0.0; N];
        po[0] = vr;
        let mut np = 1;
        let mut spill: Vec<(Attr, f64)> = Vec::new();
        for (k, (a, w)) in wires.enumerate() {
            let mut c = EPS_OWN;
            for j in 0..np {
                if let Some(g) = self.clr(&pa[j], &a) {
                    c = c.max(po[j] + g + w / 2.0);
                }
            }
            for (sa, so) in &spill {
                if let Some(g) = self.clr(sa, &a) {
                    c = c.max(so + g + w / 2.0);
                }
            }
            emit(k, a, w, c);
            if np < N {
                pa[np] = a;
                po[np] = c + w / 2.0;
                np += 1;
            } else {
                spill.push((a, c + w / 2.0));
            }
        }
    }

    /// Wire attributes of the crossings on edge e, in order from edge.a.
    pub fn edge_wires(&self, l: usize, e: u32) -> Vec<(Attr, f64)> {
        self.layers[l].cross[e as usize].iter().map(|&c| self.cross_wire(c)).collect()
    }

    #[inline]
    fn cross_wire(&self, c: u32) -> (Attr, f64) {
        let p = &self.pieces[self.crosses[c as usize].piece as usize];
        (self.net_attr(p.net, p.layer as usize), p.width)
    }

    /// Feasible interval (distance from edge.a) for a new wire inserted at `ins` on edge e, if it fits.
    pub fn edge_fit(&self, l: usize, e: u32, ins: usize, net: u32, width: f64) -> Option<(f64, f64)> {
        let ed = self.layers[l].cdt.edges[e as usize];
        let len = self.edge_len(l, e);
        let list = &self.layers[l].cross[e as usize];
        let n = list.len() + 1;
        let me = (self.net_attr(net, l), width);
        let wire = |k: usize| -> (Attr, f64) {
            if k < ins {
                self.cross_wire(list[k])
            } else if k == ins {
                me
            } else {
                self.cross_wire(list[k - 1])
            }
        };
        // whole edge must fit: every element's packing from a must leave room for b
        let (ba, br) = self.vert_attr(l, ed.b);
        let mut need = 0.0f64;
        let mut lo = 0.0;
        self.pack_each(l, ed.a, (0..n).map(wire), |k, a, w, c| {
            if let Some(g) = self.clr(&a, &ba) {
                need = need.max(c + w / 2.0 + g + br);
            }
            if k == ins {
                lo = c;
            }
        });
        if need > len + 1e-9 {
            return None;
        }
        let fb = self.pack_last(l, ed.b, (ins..n).rev().map(wire));
        let hi = len - fb;
        if lo > hi + 1e-9 {
            return None;
        }
        Some((lo, hi))
    }

    pub fn tri_of_edge(&self, l: usize, e: u32, from_a: bool) -> u32 {
        self.layers[l].etri[e as usize][if from_a { 0 } else { 1 }]
    }
    /// Triangle entered by a crossing (the one where the edge runs left -> right).
    pub fn entered_tri(&self, l: usize, c: u32) -> u32 {
        let cr = &self.crosses[c as usize];
        let ed = self.layers[l].cdt.edges[cr.edge as usize];
        self.tri_of_edge(l, cr.edge, cr.left == ed.a)
    }
    pub fn exited_tri(&self, l: usize, c: u32) -> u32 {
        let cr = &self.crosses[c as usize];
        let ed = self.layers[l].cdt.edges[cr.edge as usize];
        self.tri_of_edge(l, cr.edge, cr.left != ed.a)
    }

    pub fn term_point(&self, l: usize, t: Term) -> P2 {
        if t.vtx != NONE {
            self.vpos(l, t.vtx)
        } else {
            self.objs[t.obj as usize].s.center
        }
    }

    // ------------------------------------------------------------------ pieces
    /// Add a piece with crossings given as (edge, insertion index from edge.a, left vertex).
    pub fn add_piece(&mut self, net: u32, l: usize, width: f64, a: Term, b: Term, cx: &[(u32, usize, u32)]) -> u32 {
        let pid = self.pieces.len() as u32;
        let mut ids = Vec::with_capacity(cx.len());
        let mut done: Vec<(u32, usize)> = Vec::new();
        for &(e, g, left) in cx {
            let shift = done.iter().filter(|(de, dg)| *de == e && *dg <= g).count();
            done.push((e, g));
            let cid = self.crosses.len() as u32;
            self.crosses.push(Crossing { piece: pid, edge: e, left, alive: true, idx: ids.len() as u32 });
            self.save_cross_list(l, e);
            let list = &mut self.layers[l].cross[e as usize];
            let at = (g + shift).min(list.len());
            list.insert(at, cid);
            ids.push(cid);
        }
        let mut along = NONE;
        if cx.is_empty() && a.vtx != NONE && b.vtx != NONE {
            if let Some(e) = self.layers[l].cdt.find_edge(a.vtx, b.vtx) {
                self.save_along(l, e);
                self.layers[l].along[e as usize] = pid;
                along = e;
            }
        }
        self.pieces.push(Piece { net, layer: l as u8, width, a, b, cx: ids, poly: Vec::new(), cparam: Vec::new(), alive: true, along, overlap: false, problem: None, frozen: false });
        for o in [a.obj, b.obj] {
            self.save_obj_pieces(o);
            self.obj_pieces[o as usize].push(pid);
        }
        self.dirty.push(pid);
        for &(e, _, _) in cx {
            for &c in &self.layers[l].cross[e as usize] {
                self.dirty.push(self.crosses[c as usize].piece);
            }
        }
        pid
    }

    pub fn remove_piece(&mut self, p: u32) {
        if !self.pieces[p as usize].alive {
            return;
        }
        self.unindex_piece(p);
        self.save_piece(p);
        let l = self.pieces[p as usize].layer as usize;
        let cx = self.pieces[p as usize].cx.clone();
        for c in cx {
            let e = self.crosses[c as usize].edge;
            self.save_cross_list(l, e);
            self.layers[l].cross[e as usize].retain(|&x| x != c);
            self.save_cross(c);
            self.crosses[c as usize].alive = false;
            for &o in &self.layers[l].cross[e as usize] {
                self.dirty.push(self.crosses[o as usize].piece);
            }
        }
        let al = self.pieces[p as usize].along;
        if al != NONE && self.layers[l].along[al as usize] == p {
            self.save_along(l, al);
            self.layers[l].along[al as usize] = NONE;
        }
        let (a, b) = (self.pieces[p as usize].a.obj, self.pieces[p as usize].b.obj);
        for o in [a, b] {
            self.save_obj_pieces(o);
            self.obj_pieces[o as usize].retain(|&x| x != p);
        }
        self.pieces[p as usize].alive = false;
    }

    /// Topology invariants (debugging): every live crossing is in its edge's list and vice versa.
    pub fn validate_topology(&self) -> Result<(), String> {
        for (pi, pc) in self.pieces.iter().enumerate() {
            if !pc.alive {
                continue;
            }
            let l = pc.layer as usize;
            for (k, &c) in pc.cx.iter().enumerate() {
                let cr = &self.crosses[c as usize];
                if !cr.alive || cr.piece != pi as u32 {
                    return Err(format!("piece {} (layer {}, net {}, terms {:?} {:?}, {} crossings) crossing #{} id {} alive {} owner {} edge {} alive {}", pi, l, self.b.nets[pc.net as usize].name, pc.a, pc.b, pc.cx.len(), k, c, cr.alive, cr.piece, cr.edge, self.layers[l].cdt.edges[cr.edge as usize].alive));
                }
                if cr.idx != k as u32 {
                    return Err(format!("piece {} crossing {} idx {} != {}", pi, c, cr.idx, k));
                }
                if !self.layers[l].cdt.edges[cr.edge as usize].alive {
                    return Err(format!("piece {} crossing {} on dead edge {}", pi, c, cr.edge));
                }
                if !self.layers[l].cross[cr.edge as usize].contains(&c) {
                    return Err(format!("piece {} crossing {} missing from edge {} list", pi, c, cr.edge));
                }
                if k + 1 < pc.cx.len() {
                    let t1 = self.entered_tri(l, c);
                    let t2 = self.exited_tri(l, pc.cx[k + 1]);
                    if t1 != t2 {
                        return Err(format!("piece {} (net {}) crossings #{} and #{} do not share a triangle ({} vs {})", pi, self.b.nets[pc.net as usize].name, k, k + 1, t1, t2));
                    }
                }
            }
        }
        for (l, ly) in self.layers.iter().enumerate() {
            for (e, list) in ly.cross.iter().enumerate() {
                for &c in list {
                    let cr = &self.crosses[c as usize];
                    if !cr.alive || cr.edge != e as u32 || !self.pieces[cr.piece as usize].alive {
                        return Err(format!("layer {} edge {} holds stale crossing {}", l, e, c));
                    }
                }
            }
        }
        Ok(())
    }

    // ------------------------------------------------------------------ realisation
    /// Portal discs for a piece from its topology.
    pub fn portals(&self, p: u32) -> Vec<(Disc, Disc)> {
        let pc = &self.pieces[p as usize];
        let l = pc.layer as usize;
        let me = (self.net_attr(pc.net, pc.layer as usize), pc.width);
        let mut out = Vec::with_capacity(pc.cx.len());
        for &c in &pc.cx {
            let cr = &self.crosses[c as usize];
            let ed = self.layers[l].cdt.edges[cr.edge as usize];
            let list = &self.layers[l].cross[cr.edge as usize];
            let i = list.iter().position(|&x| x == c).unwrap();
            let ra = self.pack_last(l, ed.a, list[..i].iter().map(|&x| self.cross_wire(x)).chain(std::iter::once(me)));
            let rb = self.pack_last(l, ed.b, list[i + 1..].iter().rev().map(|&x| self.cross_wire(x)).chain(std::iter::once(me)));
            let mut da = Disc { c: self.vpos(l, ed.a), r: ra + TOL, id: ed.a };
            let mut db = Disc { c: self.vpos(l, ed.b), r: rb + TOL, id: ed.b };
            if let Some(d) = self.own_vertex_guard(l, ed.a, pc.net, pc.width) {
                if d.r > da.r {
                    da = d;
                }
            }
            if let Some(d) = self.own_vertex_guard(l, ed.b, pc.net, pc.width) {
                if d.r > db.r {
                    db = d;
                }
            }
            if cr.left == ed.a {
                out.push((da, db));
            } else {
                out.push((db, da));
            }
        }
        // a vertex on one side keeps its largest radius over the whole piece (it is wrapped once)
        let mut maxr: HashMap<(u32, bool), f64> = HashMap::new();
        for (a, b) in &out {
            let e = maxr.entry((a.id, true)).or_insert(0.0);
            *e = e.max(a.r);
            let e = maxr.entry((b.id, false)).or_insert(0.0);
            *e = e.max(b.r);
        }
        for (a, b) in out.iter_mut() {
            a.r = maxr[&(a.id, true)];
            b.r = maxr[&(b.id, false)];
        }
        out
    }

    /// A wire may touch its own net's vertices, but other copper reaching into that vertex's surroundings
    /// must still be cleared: return the clearance disc of the copper that reaches furthest past the vertex.
    pub fn own_vertex_guard(&self, l: usize, v: u32, net: u32, width: f64) -> Option<Disc> {
        let o = self.layers[l].vobj[v as usize];
        if o == NONE || !self.objs[o as usize].alive || self.objs[o as usize].s.net != net {
            return None;
        }
        let pv = self.vpos(l, v);
        let me = self.net_attr(net, l);
        let reach = self.max_clear + width;
        let mut items = Vec::new();
        self.layers[l].grid.query(pv - p2(reach, reach), pv + p2(reach, reach), &mut items);
        let mut best: Option<(f64, Disc)> = None;
        for it in &items {
            if it.kind == IT_SEG || it.id == o {
                continue;
            }
            let ob = &self.objs[it.id as usize];
            if !ob.alive {
                continue;
            }
            if let Some(g) = self.clr(&me, &self.obj_attr(it.id)) {
                let (q, r) = if it.kind == IT_RING {
                    let n = ob.s.pts.len();
                    let (a, b) = (ob.s.pts[it.sub as usize], ob.s.pts[(it.sub as usize + 1) % n]);
                    (a.lerp(b, seg_param(a, b, pv)), 0.0)
                } else {
                    (hull_closest(&ob.s.pts, pv), ob.s.r)
                };
                let rad = r + g + width / 2.0 + TOL;
                let over = rad - q.dist(pv);
                if over > 0.0 && best.as_ref().map_or(true, |b| over > b.0) {
                    best = Some((over, Disc { c: q, r: rad, id: 0x8000_0000 | it.id }));
                }
            }
        }
        best.map(|b| b.1)
    }

    /// Geometry of piece p from its topology plus extra clearance discs (insert index, disc, side +1 left).
    fn realize_with(&self, p: u32, extra: &[(usize, Disc, f64)]) -> (Vec<P2>, Vec<f64>, Option<P2>) {
        let l = self.pieces[p as usize].layer as usize;
        let s = self.term_point(l, self.pieces[p as usize].a);
        let t = self.term_point(l, self.pieces[p as usize].b);
        let base = self.portals(p);
        let mut portals = Vec::with_capacity(base.len() + extra.len());
        let endd = Disc { c: t, r: 0.0, id: u32::MAX };
        for i in 0..=base.len() {
            for &(at, d, side) in extra {
                if at == i {
                    let opp = if base.is_empty() { endd } else { base[i.min(base.len() - 1)].clone() .0 };
                    let oppr = if base.is_empty() { endd } else { base[i.min(base.len() - 1)].1 };
                    portals.push(if side > 0.0 { (d, oppr) } else { (opp, d) });
                }
            }
            if i < base.len() {
                portals.push(base[i]);
            }
        }
        let (wraps, overlap) = funnel(s, &portals, t);
        if overlap.is_some() && std::env::var("TOPOR_DEBUG").is_ok() {
            eprintln!("funnel overlap piece {} layer {} start {:?} end {:?} terms {:?} {:?} objs {} {}", p, l, s, t, self.pieces[p as usize].a, self.pieces[p as usize].b, self.objs[self.pieces[p as usize].a.obj as usize].s.name, self.objs[self.pieces[p as usize].b.obj as usize].s.name);
            for (i, (a, b)) in portals.iter().enumerate() {
                eprintln!("  portal {}: L {:?} r {:.4} id {}  R {:?} r {:.4} id {}", i, a.c, a.r, a.id, b.c, b.r, b.id);
            }
        }
        let (poly, bad_turn) = polyline(s, &wraps, t, TOL);
        if bad_turn && std::env::var("TOPOR_DEBUG").is_ok() {
            eprintln!("bad turn piece {} layer {} start {:?} end {:?} {} -> {}", p, l, s, t, self.objs[self.pieces[p as usize].a.obj as usize].s.name, self.objs[self.pieces[p as usize].b.obj as usize].s.name);
            for w in &wraps {
                eprintln!("   wrap {:?} r {:.4} side {} id {:x}", w.d.c, w.d.r, w.side, w.d.id);
            }
            for (i, (a, b)) in portals.iter().enumerate() {
                eprintln!("   portal {}: L {:?} r {:.4} id {:x}  R {:?} r {:.4} id {:x}", i, a.c, a.r, a.id, b.c, b.r, b.id);
            }
        }
        let mut cparam = Vec::with_capacity(base.len());
        let mut seg = 0usize;
        let mut problem = overlap;
        if bad_turn && problem.is_none() {
            problem = Some(s);
        }
        for &c in &self.pieces[p as usize].cx {
            let ed = self.layers[l].cdt.edges[self.crosses[c as usize].edge as usize];
            let (ea, eb) = (self.vpos(l, ed.a), self.vpos(l, ed.b));
            let mut found = None;
            for k in seg..poly.len().saturating_sub(1) {
                if let Some((tt, _)) = seg_cross_params(poly[k], poly[k + 1], ea, eb) {
                    found = Some(k as f64 + tt.min(0.999999));
                    seg = k;
                    break;
                }
            }
            match found {
                Some(v) => cparam.push(v),
                None => {
                    if problem.is_none() {
                        problem = Some(ea.lerp(eb, 0.5));
                    }
                    cparam.push(seg as f64);
                }
            }
        }
        (poly, cparam, problem)
    }

    fn set_geometry(&mut self, p: u32, poly: Vec<P2>, cparam: Vec<f64>, problem: Option<P2>) {
        self.unindex_piece(p);
        self.save_piece(p);
        let pc = &mut self.pieces[p as usize];
        pc.poly = poly;
        pc.cparam = cparam;
        pc.overlap = problem.is_some();
        pc.problem = problem;
        self.index_piece(p);
    }

    pub fn realize(&mut self, p: u32) -> bool {
        if !self.pieces[p as usize].alive || self.pieces[p as usize].frozen {
            return true;
        }
        let (poly, cparam, problem) = self.realize_with(p, &[]);
        let ok = problem.is_none();
        self.set_geometry(p, poly, cparam, problem);
        ok
    }

    /// Triangles a piece passes through: T_0 (before its first crossing) .. T_n (after its last).
    pub fn piece_tris(&self, p: u32) -> Vec<u32> {
        let pc = &self.pieces[p as usize];
        let l = pc.layer as usize;
        let n = pc.cx.len();
        let mut tris = Vec::with_capacity(n + 1);
        if n == 0 {
            return tris;
        }
        tris.push(if pc.a.vtx != NONE { self.exited_tri(l, pc.cx[0]) } else { NONE });
        for k in 0..n {
            tris.push(self.entered_tri(l, pc.cx[k]));
        }
        if pc.b.vtx == NONE {
            tris[n] = NONE;
        }
        tris
    }

    /// Sides (+1 left, -1 right, 0 own terminal) of triangle t's vertices relative to piece p's passage
    /// through t, which enters at crossing `cin` (or starts at vertex) and leaves at `cout` (or ends).
    fn vertex_sides(&self, l: usize, t: u32, cin: Option<u32>, cout: Option<u32>, vstart: u32, vend: u32) -> [f64; 3] {
        let tr = self.layers[l].cdt.tris[t as usize];
        let mut side = [0.0f64; 3];
        let idx = |v: u32| tr.v.iter().position(|&x| x == v);
        if let Some(c) = cin {
            let cr = &self.crosses[c as usize];
            let ed = self.layers[l].cdt.edges[cr.edge as usize];
            let right = if cr.left == ed.a { ed.b } else { ed.a };
            if let Some(i) = idx(cr.left) {
                side[i] = 1.0;
            }
            if let Some(i) = idx(right) {
                side[i] = -1.0;
            }
        }
        if let Some(c) = cout {
            let cr = &self.crosses[c as usize];
            let ed = self.layers[l].cdt.edges[cr.edge as usize];
            let right = if cr.left == ed.a { ed.b } else { ed.a };
            if let Some(i) = idx(cr.left) {
                side[i] = 1.0;
            }
            if let Some(i) = idx(right) {
                side[i] = -1.0;
            }
        }
        for v in [vstart, vend] {
            if v != NONE {
                if let Some(i) = idx(v) {
                    side[i] = 0.0;
                }
            }
        }
        side
    }

    /// Other nets' wires that end inside a triangle p passes through (spokes into a via or pad) bound
    /// p's sleeve there; edge stacking cannot see them. Return clearance discs along them.
    fn spoke_walls(&self, p: u32) -> Vec<(usize, Disc, f64)> {
        let pc = &self.pieces[p as usize];
        let l = pc.layer as usize;
        let n = pc.cx.len();
        let tris = self.piece_tris(p);
        let me = self.net_attr(pc.net, pc.layer as usize);
        let mut out = Vec::new();
        let mut id = 0xC000_0000u32;
        for k in 0..tris.len() {
            let t = tris[k];
            if t == NONE || self.layers[l].tflag[t as usize] & T_BLOCK != 0 {
                continue;
            }
            let tr = self.layers[l].cdt.tris[t as usize];
            let cin = if k > 0 { Some(pc.cx[k - 1]) } else { None };
            let cout = if k < n { Some(pc.cx[k]) } else { None };
            let side = self.vertex_sides(l, t, cin, cout, if k == 0 { pc.a.vtx } else { NONE }, if k == n { pc.b.vtx } else { NONE });
            // approximate direction of travel through t
            let mid = |c: u32| {
                let ed = self.layers[l].cdt.edges[self.crosses[c as usize].edge as usize];
                self.vpos(l, ed.a).lerp(self.vpos(l, ed.b), 0.5)
            };
            let pin = cin.map(mid).unwrap_or_else(|| self.term_point(l, pc.a));
            let pout = cout.map(mid).unwrap_or_else(|| self.term_point(l, pc.b));
            let dir = (pout - pin).unit();
            let mut local: Vec<(f64, Disc, f64)> = Vec::new();
            for ei in 0..3 {
                let e = tr.e[ei];
                for &c in &self.layers[l].cross[e as usize] {
                    let cr = &self.crosses[c as usize];
                    let s_id = cr.piece;
                    if s_id == p {
                        continue;
                    }
                    let sp = &self.pieces[s_id as usize];
                    if !sp.alive || sp.net == pc.net || sp.poly.len() < 2 || sp.cparam.len() != sp.cx.len() {
                        continue;
                    }
                    let i = cr.idx as usize;
                    let entering = cr.left == tr.v[ei];
                    // the spoke's other end inside t: a terminal vertex of t, or a hull crossing that ends it
                    let (wall_side, p0, p1) = if entering {
                        if i + 1 == sp.cx.len() {
                            match tr.v.iter().position(|&v| v == sp.b.vtx) {
                                Some(j) if sp.b.vtx != NONE => (side[j], sp.cparam[i], (sp.poly.len() - 1) as f64),
                                _ => continue,
                            }
                        } else if i + 2 == sp.cx.len() && sp.b.vtx == NONE {
                            let hc = &self.crosses[sp.cx[i + 1] as usize];
                            let hed = self.layers[l].cdt.edges[hc.edge as usize];
                            match tr.v.iter().position(|&v| v == hed.a) {
                                Some(j) if tr.e.contains(&hc.edge) => (side[j], sp.cparam[i], sp.cparam[i + 1]),
                                _ => continue,
                            }
                        } else {
                            continue;
                        }
                    } else if i == 0 {
                        match tr.v.iter().position(|&v| v == sp.a.vtx) {
                            Some(j) if sp.a.vtx != NONE => (side[j], 0.0, sp.cparam[0]),
                            _ => continue,
                        }
                    } else if i == 1 && sp.a.vtx == NONE {
                        let hc = &self.crosses[sp.cx[0] as usize];
                        let hed = self.layers[l].cdt.edges[hc.edge as usize];
                        match tr.v.iter().position(|&v| v == hed.a) {
                            Some(j) if tr.e.contains(&hc.edge) => (side[j], sp.cparam[0], sp.cparam[1]),
                            _ => continue,
                        }
                    } else {
                        continue;
                    };
                    if wall_side == 0.0 {
                        continue;
                    }
                    let need = match self.clr(&me, &self.net_attr(sp.net, sp.layer as usize)) {
                        Some(g) => g + (pc.width + sp.width) / 2.0 + 0.002,
                        None => continue,
                    };
                    let seg = self.sub_poly(s_id, p0.min(p1), p0.max(p1));
                    let mut pts = Vec::new();
                    for w in seg.windows(2) {
                        let m = ((w[0].dist(w[1]) / 0.03).ceil() as usize).clamp(1, 40);
                        for j in 0..m {
                            pts.push(w[0].lerp(w[1], j as f64 / m as f64));
                        }
                    }
                    if let Some(&last) = seg.last() {
                        pts.push(last);
                    }
                    // a wire ending (starting) in t is not walled by spoke copper beyond its end (before
                    // its start): the funnel would have to pass it there
                    let reach = (pout - pin).dot(dir);
                    for q in pts {
                        let along = (q - pin).dot(dir);
                        if (k == n && pc.b.vtx != NONE && along > reach) || (k == 0 && pc.a.vtx != NONE && along < 0.0) {
                            continue;
                        }
                        id += 1;
                        local.push((along, Disc { c: q, r: need, id }, wall_side));
                    }
                }
            }
            local.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            for (_, d, sd) in local {
                out.push((k, d, sd));
            }
        }
        out
    }

    /// Clearance discs for everything piece p violates, on the side where it lies.
    fn feedback_discs(&self, p: u32, next_id: &mut u32) -> Vec<(usize, Disc, f64)> {
        let mut viol = Vec::new();
        self.check_piece(p, &mut viol);
        let pc = &self.pieces[p as usize];
        if viol.is_empty() || pc.poly.len() < 2 {
            return vec![];
        }
        let me = self.net_attr(pc.net, pc.layer as usize);
        let mut out = Vec::new();
        for v in viol.iter().take(12) {
            // closest segment of p to the violation
            let mut best = (f64::INFINITY, 0usize);
            for k in 0..pc.poly.len() - 1 {
                let d = seg_point_dist(pc.poly[k], pc.poly[k + 1], v.at);
                if d < best.0 {
                    best = (d, k);
                }
            }
            let k = best.1;
            let (a, b) = (pc.poly[k], pc.poly[k + 1]);
            let (q, rad) = if let Some(&o) = v.objs.first() {
                let ob = &self.objs[o as usize];
                let g = match self.clr(&me, &self.obj_attr(o)) {
                    Some(g) => g,
                    None => continue,
                };
                // nearest core point of the object to the segment
                let mut bq = ob.s.pts[0];
                let mut bd = f64::INFINITY;
                let n = ob.s.pts.len();
                for i in 0..n.max(1) {
                    let (u, w) = (ob.s.pts[i % n], ob.s.pts[(i + 1) % n]);
                    for t in 0..=8 {
                        let c = u.lerp(w, t as f64 / 8.0);
                        let dd = seg_point_dist(a, b, c);
                        if dd < bd {
                            bd = dd;
                            bq = c;
                        }
                    }
                }
                let r = if matches!(ob.s.kind, Kind::Outline | Kind::Keepout) { 0.0 } else { ob.s.r };
                (bq, r + g + pc.width / 2.0 + TOL + 0.001)
            } else if v.pieces.len() == 2 {
                let other = if v.pieces[0] == p { v.pieces[1] } else { v.pieces[0] };
                let q = &self.pieces[other as usize];
                let g = match self.clr(&me, &self.net_attr(q.net, q.layer as usize)) {
                    Some(g) => g,
                    None => continue,
                };
                let mut bq = q.poly[0];
                let mut bd = f64::INFINITY;
                for w in q.poly.windows(2) {
                    for t in 0..=8 {
                        let c = w[0].lerp(w[1], t as f64 / 8.0);
                        let dd = seg_point_dist(a, b, c);
                        if dd < bd {
                            bd = dd;
                            bq = c;
                        }
                    }
                }
                (bq, g + (pc.width + q.width) / 2.0 + TOL + 0.001)
            } else {
                continue;
            };
            let foot = a.lerp(b, seg_param(a, b, q));
            let side = (b - a).cross(q - foot);
            if side.abs() < 1e-9 {
                continue;
            }
            let param = k as f64 + seg_param(a, b, q);
            let at = pc.cparam.iter().filter(|&&c| c <= param).count();
            *next_id += 1;
            out.push((at, Disc { c: q, r: rad, id: *next_id }, side.signum()));
        }
        out
    }

    /// Pair coupling: where piece p's crossings run beside one partner piece over consecutive edges,
    /// p's geometry there is the partner's path offset by gap + width (the ribbon). Returns whether
    /// anything changed. The partner (lower net id) keeps its own rubber band.
    fn couple_geometry(&mut self, p: u32) -> bool {
        let pc = &self.pieces[p as usize];
        if !pc.alive || (pc.net as usize) >= self.b.partner.len() {
            return false;
        }
        let pn = self.b.partner[pc.net as usize];
        if pn == NONE || pc.net < pn || pc.poly.len() < 2 || pc.cparam.len() != pc.cx.len() {
            return false;
        }
        let l = pc.layer as usize;
        // partner crossing beside each of p's crossings
        let mut nb: Vec<Option<u32>> = Vec::with_capacity(pc.cx.len());
        for &c in &pc.cx {
            let e = self.crosses[c as usize].edge;
            let list = &self.layers[l].cross[e as usize];
            let i = list.iter().position(|&x| x == c).unwrap();
            let pick = |j: usize| -> Option<u32> {
                let q = list[j];
                if self.pieces[self.crosses[q as usize].piece as usize].net == pn {
                    Some(q)
                } else {
                    None
                }
            };
            let mut r = None;
            if i > 0 {
                r = pick(i - 1);
            }
            if r.is_none() && i + 1 < list.len() {
                r = pick(i + 1);
            }
            nb.push(r);
        }
        // runs over consecutive crossings of one partner piece
        let mut runs: Vec<(usize, usize, u32)> = Vec::new();
        let mut k = 0;
        while k < nb.len() {
            if let Some(q0) = nb[k] {
                let qp = self.crosses[q0 as usize].piece;
                let mut j = k;
                let mut dir = 0i64;
                while j + 1 < nb.len() {
                    match nb[j + 1] {
                        Some(q1) if self.crosses[q1 as usize].piece == qp => {
                            let di = self.crosses[q1 as usize].idx as i64 - self.crosses[nb[j].unwrap() as usize].idx as i64;
                            if di.abs() != 1 || (dir != 0 && di != dir) {
                                break;
                            }
                            dir = di;
                            j += 1;
                        }
                        _ => break,
                    }
                }
                if j > k {
                    runs.push((k, j, qp));
                }
                k = j + 1;
            } else {
                k += 1;
            }
        }
        if runs.is_empty() {
            return false;
        }
        let me = self.net_attr(pc.net, pc.layer as usize);
        let width = pc.width;
        let start = self.term_point(l, pc.a);
        let end = self.term_point(l, pc.b);
        let portals = self.portals(p);
        let old = (pc.poly.clone(), pc.cparam.clone());
        let mut runs: Vec<(usize, usize, u32)> = runs;
        for _attempt in 0..16 {
            if runs.is_empty() {
                break;
            }
            // offsets of the partner over each run, extended a little past the run's first and last edge
            let mut offs: Vec<Vec<P2>> = Vec::new();
            let mut bad = false;
            for &(k0, k1, qp) in &runs {
                let q = &self.pieces[qp as usize];
                let qi0 = self.crosses[nb[k0].unwrap() as usize].idx as usize;
                let qi1 = self.crosses[nb[k1].unwrap() as usize].idx as usize;
                let (qa, qb) = (q.cparam[qi0], q.cparam[qi1]);
                let mut qs = self.sub_poly(qp, qa.min(qb), qa.max(qb));
                if qa > qb {
                    qs.reverse();
                }
                let g = self.clr(&me, &self.net_attr(q.net, q.layer as usize));
                if qs.len() < 2 || g.is_none() {
                    bad = true;
                    break;
                }
                let d = g.unwrap() + (width + q.width) / 2.0;
                let pp = self.point_at(p, old.1[k0]);
                let side = (qs[1] - qs[0]).unit().cross(pp - qs[0]).signum();
                let mut off = offset_polyline(&qs, side * d);
                let n = off.len();
                let e0 = off[0] - (off[1] - off[0]).unit() * 0.04;
                let e1 = off[n - 1] + (off[n - 1] - off[n - 2]).unit() * 0.04;
                off.insert(0, e0);
                off.push(e1);
                offs.push(off);
            }
            if bad {
                break;
            }
            // taut stretches between runs (and to the terminals) through the remaining portals
            let mut poly: Vec<P2> = Vec::new();
            let mut from = start;
            let mut next_k = 0usize;
            for (ri, &(k0, k1, _)) in runs.iter().enumerate() {
                let to = offs[ri][0];
                let (w, _) = funnel(from, &portals[next_k..k0], to);
                let (seg, _) = polyline(from, &w, to, TOL);
                poly.extend_from_slice(&seg);
                poly.extend_from_slice(&offs[ri][1..]);
                from = *offs[ri].last().unwrap();
                next_k = k1 + 1;
            }
            let (w, _) = funnel(from, &portals[next_k..], end);
            let (seg, _) = polyline(from, &w, end, TOL);
            poly.extend_from_slice(&seg[1..]);
            poly.dedup_by(|x, y| x.dist(*y) < 1e-7);
            self.unindex_piece(p);
            self.save_piece(p);
            self.pieces[p as usize].poly = poly;
            self.index_piece(p);
            let consistent = self.refresh_cparam(p);
            let mut v = Vec::new();
            if consistent {
                self.check_piece(p, &mut v);
            }
            if consistent && v.is_empty() {
                if std::env::var("TOPOR_PAIRDBG").is_ok() {
                    eprintln!("   couple {} ({}): runs {:?} ok", p, self.b.nets[self.pieces[p as usize].net as usize].name, runs.iter().map(|r| (r.0, r.1)).collect::<Vec<_>>());
                }
                return true;
            }
            // undo; shrink the run nearest the problem (or the longest if the topology broke)
            self.unindex_piece(p);
            self.pieces[p as usize].poly = old.0.clone();
            self.pieces[p as usize].cparam = old.1.clone();
            self.index_piece(p);
            let at = v.first().map(|x| x.at);
            let mut best = 0usize;
            let mut bd = f64::INFINITY;
            for (ri, r) in runs.iter().enumerate() {
                let d = match at {
                    Some(at) => self.point_at(p, old.1[r.0]).dist(at).min(self.point_at(p, old.1[r.1]).dist(at)),
                    None => -((r.1 - r.0) as f64),
                };
                if d < bd {
                    bd = d;
                    best = ri;
                }
            }
            let r = runs[best];
            let near_start = match at {
                Some(at) => self.point_at(p, old.1[r.0]).dist(at) < self.point_at(p, old.1[r.1]).dist(at),
                None => true,
            };
            let nr = if near_start { (r.0 + 1, r.1, r.2) } else { (r.0, r.1 - 1, r.2) };
            if nr.1 > nr.0 {
                runs[best] = nr;
            } else {
                runs.remove(best);
            }
        }
        false
    }

    /// Replace a piece's polyline if the result is clean (final geometry edits such as meanders).
    /// Replace a piece's geometry without checking it (the caller checks the result as a whole).
    pub fn set_poly(&mut self, p: u32, poly: Vec<P2>) {
        self.unindex_piece(p);
        self.save_piece(p);
        self.pieces[p as usize].poly = poly;
        self.index_piece(p);
        self.refresh_cparam(p);
    }

    pub fn try_set_poly(&mut self, p: u32, poly: Vec<P2>) -> bool {
        let old = (self.pieces[p as usize].poly.clone(), self.pieces[p as usize].cparam.clone());
        self.unindex_piece(p);
        self.pieces[p as usize].poly = poly;
        self.index_piece(p);
        let mut v = Vec::new();
        self.check_piece(p, &mut v);
        if v.is_empty() {
            self.refresh_cparam(p);
            return true;
        }
        self.unindex_piece(p);
        self.pieces[p as usize].poly = old.0;
        self.pieces[p as usize].cparam = old.1;
        self.index_piece(p);
        false
    }

    pub fn realize_dirty(&mut self) -> Vec<u32> {
        let mut d = std::mem::take(&mut self.dirty);
        // a pair member follows its partner: when the partner moves, re-offset it too
        let extra: Vec<u32> = d
            .iter()
            .filter_map(|&p| {
                let net = self.pieces[p as usize].net as usize;
                self.b.partner.get(net).cloned().filter(|&x| x != NONE)
            })
            .collect();
        if !extra.is_empty() {
            for (i, pc) in self.pieces.iter().enumerate() {
                if pc.alive && extra.contains(&pc.net) {
                    d.push(i as u32);
                }
            }
        }
        d.sort();
        d.dedup();
        d.retain(|&p| self.pieces[p as usize].alive && !self.pieces[p as usize].frozen);
        for &p in &d {
            self.realize(p);
        }
        // pair followers take their ribbon shape first, so everything else adapts to it
        let mut coupled: Vec<u32> = Vec::new();
        for &p in &d {
            if self.couple_geometry(p) {
                coupled.push(p);
            }
        }
        // spokes (wires ending at a via or pad) bound other wires' sleeves: second pass with walls
        let mut extra: HashMap<u32, Vec<(usize, Disc, f64)>> = HashMap::new();
        for _round in 0..2 {
            for &p in &d {
                if coupled.contains(&p) {
                    continue;
                }
                let walls = self.spoke_walls(p);
                if walls.is_empty() {
                    continue;
                }
                let (poly, cparam, ok) = self.realize_with(p, &walls);
                self.set_geometry(p, poly, cparam, ok);
                extra.insert(p, walls);
            }
        }
        // copper the sleeve cannot see (beyond an edge the wire runs along): exact check, then pull the
        // wire taut around clearance discs of whatever it came too close to
        let mut next_id = 0xE000_0000u32;
        for _round in 0..3 {
            let mut any = false;
            for &p in &d {
                if coupled.contains(&p) {
                    continue;
                }
                let fb = self.feedback_discs(p, &mut next_id);
                if fb.is_empty() {
                    continue;
                }
                any = true;
                let e = extra.entry(p).or_default();
                for f in fb {
                    // one disc per obstacle point and side
                    if let Some(x) = e.iter_mut().find(|x| x.2 == f.2 && x.1.c.dist(f.1.c) < 0.01) {
                        x.1.r = x.1.r.max(f.1.r);
                    } else {
                        e.push(f);
                    }
                }
                let mut ex = e.clone();
                ex.sort_by(|x, y| x.0.cmp(&y.0));
                let (poly, cparam, ok) = self.realize_with(p, &ex);
                self.set_geometry(p, poly, cparam, ok);
            }
            if !any {
                break;
            }
        }
        // followers last, against the leaders' final shape; one whose offset no longer fits is pulled
        // taut on its own with walls and feedback instead
        for &p in &d {
            if self.couple_geometry(p) || !coupled.contains(&p) {
                continue;
            }
            self.realize(p);
            let walls = self.spoke_walls(p);
            let mut ex = walls.clone();
            if !walls.is_empty() {
                let (poly, cparam, ok) = self.realize_with(p, &walls);
                self.set_geometry(p, poly, cparam, ok);
            }
            for _round in 0..3 {
                let fb = self.feedback_discs(p, &mut next_id);
                if fb.is_empty() {
                    break;
                }
                ex.extend(fb);
                ex.sort_by(|x, y| x.0.cmp(&y.0));
                let (poly, cparam, ok) = self.realize_with(p, &ex);
                self.set_geometry(p, poly, cparam, ok);
            }
        }
        d
    }

    pub fn point_at(&self, p: u32, param: f64) -> P2 {
        let pc = &self.pieces[p as usize];
        let k = (param.floor() as usize).min(pc.poly.len().saturating_sub(2));
        let t = param - k as f64;
        pc.poly[k].lerp(pc.poly[k + 1], t.clamp(0.0, 1.0))
    }

    /// Polyline of piece p between two parameters.
    pub fn sub_poly(&self, p: u32, a: f64, b: f64) -> Vec<P2> {
        let pc = &self.pieces[p as usize];
        let mut out = vec![self.point_at(p, a)];
        let k0 = a.floor() as usize + 1;
        let k1 = b.floor() as usize;
        for k in k0..=k1.min(pc.poly.len() - 1) {
            if (k as f64) > a && (k as f64) < b {
                out.push(pc.poly[k]);
            }
        }
        out.push(self.point_at(p, b));
        out
    }

    /// Crossing parameters of a piece's current polyline against its (possibly re-threaded) crossings.
    pub fn refresh_cparam(&mut self, p: u32) -> bool {
        let l = self.pieces[p as usize].layer as usize;
        let mut seg = 0usize;
        let mut out = Vec::with_capacity(self.pieces[p as usize].cx.len());
        let poly = &self.pieces[p as usize].poly;
        let mut all = true;
        for &c in &self.pieces[p as usize].cx {
            let ed = self.layers[l].cdt.edges[self.crosses[c as usize].edge as usize];
            let (ea, eb) = (self.vpos(l, ed.a), self.vpos(l, ed.b));
            let mut found = None;
            for k in seg..poly.len().saturating_sub(1) {
                if let Some((t, _)) = seg_cross_params(poly[k], poly[k + 1], ea, eb) {
                    found = Some(k as f64 + t.min(0.999999));
                    seg = k;
                    break;
                }
            }
            match found {
                Some(v) => out.push(v),
                None => {
                    all = false;
                    out.push(seg as f64);
                }
            }
        }
        self.pieces[p as usize].cparam = out;
        all
    }

    /// A new piece's initial geometry (it is realised afterwards), so later threading can order against it.
    pub fn seed_poly(&mut self, p: u32, poly: Vec<P2>) {
        self.set_geometry(p, poly, Vec::new(), None);
        self.refresh_cparam(p);
    }

    /// Crossings of a polyline running from terminal `ta` to terminal `tb` on layer l, for add_piece:
    /// (edge, insertion index from edge.a among the crossings already there, left vertex). A polygon pad
    /// terminal starts (ends) inside the pad; its hull edge is the first (last) crossing. None when the
    /// polyline passes another fixed edge or blocked triangle, or does not end at tb.
    pub fn thread_poly(&self, l: usize, poly: &[P2], ta: Term, tb: Term) -> Option<Vec<(u32, usize, u32)>> {
        match self.thread_poly_why(l, poly, ta, tb) {
            Ok(v) => Some(v),
            Err(e) => {
                if std::env::var("TOPOR_THREADDBG").is_ok() {
                    eprintln!("   thread {} -> {}: {}", self.objs[ta.obj as usize].s.name, self.objs[tb.obj as usize].s.name, e);
                }
                None
            }
        }
    }

    fn thread_poly_why(&self, l: usize, poly: &[P2], ta: Term, tb: Term) -> Result<Vec<(u32, usize, u32)>, String> {
        if poly.len() < 2 {
            return Err("short".into());
        }
        let lay = &self.layers[l];
        let cdt = &lay.cdt;
        let pos = |v: u32| self.vpos(l, v);
        let blocked = |t: u32| lay.tflag[t as usize] & T_BLOCK != 0;
        let mut t = NONE;
        if ta.vtx != NONE {
            let o = pos(ta.vtx);
            let d0 = poly[1] - poly[0];
            for tt in cdt.star(ta.vtx) {
                let tr = cdt.tris[tt as usize];
                let i = (0..3).find(|&i| tr.v[i] == ta.vtx).ok_or("star")?;
                let (x, y) = (pos(tr.v[(i + 1) % 3]), pos(tr.v[(i + 2) % 3]));
                if (x - o).cross(d0) >= 0.0 && d0.cross(y - o) >= 0.0 {
                    t = tt;
                    break;
                }
            }
        } else {
            // a pad centre often lies on the pad's diagonal: take the side the wire leaves towards
            let d0 = poly[1] - poly[0];
            t = match self.locate(l, poly[0])? {
                crate::cdt::Loc::Inside(t) => t,
                crate::cdt::Loc::OnEdge(t, i) => {
                    let tr = cdt.tris[t as usize];
                    let dir = pos(tr.v[(i + 1) % 3]) - pos(tr.v[i]);
                    if dir.cross(d0) > 0.0 { t } else { tr.n[i] }
                }
                crate::cdt::Loc::OnVertex(_) => return Err("starts on a vertex".into()),
            };
        }
        if t == NONE {
            return Err("no start triangle".into());
        }
        // inside the start pad until its hull is crossed
        let mut inside = ta.vtx == NONE && blocked(t);
        let mut hits: Vec<(u32, u32, P2)> = Vec::new();
        let (mut e_in, mut k, mut s) = (NONE, 0usize, 0.0f64);
        let mut guard = 0;
        let mut arrived = false;
        while k + 1 < poly.len() {
            guard += 1;
            if guard > 200_000 {
                return Err("walk too long".into());
            }
            let tr = cdt.tris[t as usize];
            let (q0, q1) = (poly[k], poly[k + 1]);
            let mut best: Option<(f64, usize)> = None;
            for i in 0..3 {
                if tr.e[i] == e_in {
                    continue;
                }
                if let Some((ts, u)) = seg_cross_params(q0, q1, pos(tr.v[i]), pos(tr.v[(i + 1) % 3])) {
                    if u <= 1e-9 || u >= 1.0 - 1e-9 || ts <= s + 1e-12 {
                        continue;
                    }
                    if best.map_or(true, |b| ts < b.0) {
                        best = Some((ts, i));
                    }
                }
            }
            if std::env::var("TOPOR_THREADTRACE").is_ok() {
                eprintln!("      t {} blocked {} verts {:?} seg {} {:?}->{:?} s {:.4} best {:?} inside {}", t, blocked(t), tr.v.iter().map(|&v| (lay.vobj[v as usize], pos(v))).collect::<Vec<_>>(), k, q0, q1, s, best, inside);
            }
            match best {
                Some((ts, i)) => {
                    let e = tr.e[i];
                    let ed = cdt.edges[e as usize];
                    let nt = tr.n[i];
                    if nt == NONE {
                        return Err("left the triangulation".into());
                    }
                    let left = if (q1 - q0).cross(pos(ed.a) - q0) > 0.0 { ed.a } else { ed.b };
                    let into_block = blocked(nt);
                    if inside {
                        if !into_block {
                            if !ed.fixed {
                                return Err("left the start pad through a free edge".into());
                            }
                            hits.push((e, left, q0.lerp(q1, ts)));
                            inside = false;
                        }
                    } else if into_block || ed.fixed {
                        // the end pad's hull: last crossing
                        let own = |v: u32| lay.vobj[v as usize] == tb.obj;
                        if tb.vtx == NONE && ed.fixed && into_block && own(ed.a) && own(ed.b) {
                            hits.push((e, left, q0.lerp(q1, ts)));
                            arrived = true;
                            break;
                        }
                        return Err(format!("hit a fixed edge or blocked triangle at ({:.3},{:.3})", q0.lerp(q1, ts).x, q0.lerp(q1, ts).y));
                    } else {
                        hits.push((e, left, q0.lerp(q1, ts)));
                    }
                    t = nt;
                    e_in = e;
                    s = ts;
                }
                None => {
                    // the segment must end in this triangle, else the walk lost the wire
                    let v = tr.v.map(pos);
                    let inside_tri = (0..3).all(|i| (v[(i + 1) % 3] - v[i]).cross(q1 - v[i]) >= -1e-9);
                    if !inside_tri {
                        return Err(format!("lost the wire at segment {}", k));
                    }
                    k += 1;
                    s = 0.0;
                    e_in = NONE;
                }
            }
        }
        if inside || (tb.vtx == NONE && !arrived) || (tb.vtx != NONE && !cdt.tris[t as usize].v.contains(&tb.vtx)) {
            return Err(format!("does not end at the terminal (inside {}, arrived {})", inside, arrived));
        }
        // a back-and-forth pair on one edge cancels topologically
        let mut k = 0;
        while k + 1 < hits.len() {
            if hits[k].0 == hits[k + 1].0 {
                hits.drain(k..k + 2);
                k = k.saturating_sub(1);
            } else {
                k += 1;
            }
        }
        // grazing an edge of the terminal vertex itself (nothing between) is no crossing at all
        let grazes = |h: &(u32, u32, P2), v: u32| -> bool {
            let ed = cdt.edges[h.0 as usize];
            if v == NONE || (ed.a != v && ed.b != v) {
                return false;
            }
            let vp = pos(v);
            let dq = h.2.dist(vp);
            !lay.cross[h.0 as usize].iter().any(|&c| {
                let cr = &self.crosses[c as usize];
                let pc = &self.pieces[cr.piece as usize];
                pc.cparam.get(cr.idx as usize).map_or(true, |&t| self.point_at(cr.piece, t).dist(vp) < dq)
            })
        };
        while hits.first().map_or(false, |h| grazes(h, ta.vtx)) {
            hits.remove(0);
        }
        while hits.last().map_or(false, |h| grazes(h, tb.vtx)) {
            hits.pop();
        }
        let mut out = Vec::with_capacity(hits.len());
        for (e, left, q) in hits {
            let ea = pos(cdt.edges[e as usize].a);
            let dq = q.dist(ea);
            let mut g = 0;
            for &c in &lay.cross[e as usize] {
                let cr = &self.crosses[c as usize];
                let pc = &self.pieces[cr.piece as usize];
                let tc = *pc.cparam.get(cr.idx as usize).ok_or("neighbour without geometry")?;
                if self.point_at(cr.piece, tc).dist(ea) < dq {
                    g += 1;
                }
            }
            out.push((e, g, left));
        }
        Ok(out)
    }

    // ------------------------------------------------------------------ cavity re-threading
    /// Apply a triangulation change on layer l over `cavity`, re-threading wires from their geometry.
    /// Drop a piece's first (last) crossings that only graze an edge of its own start (end) vertex, with
    /// no other wire between: they are no crossing at all, and would contradict the taut geometry.
    fn trim_terminal_crossings(&mut self, l: usize, p: u32) {
        for end in [false, true] {
            loop {
                let pc = &self.pieces[p as usize];
                let v = if end { pc.b.vtx } else { pc.a.vtx };
                if v == NONE || pc.cx.is_empty() {
                    break;
                }
                let k = if end { pc.cx.len() - 1 } else { 0 };
                let c = pc.cx[k];
                let e = self.crosses[c as usize].edge;
                let ed = self.layers[l].cdt.edges[e as usize];
                let list = &self.layers[l].cross[e as usize];
                let i = match list.iter().position(|&x| x == c) {
                    Some(i) => i,
                    None => break,
                };
                let next_to = (ed.a == v && i == 0) || (ed.b == v && i + 1 == list.len());
                if !next_to {
                    break;
                }
                self.save_cross_list(l, e);
                self.layers[l].cross[e as usize].remove(i);
                self.save_cross(c);
                self.crosses[c as usize].alive = false;
                self.save_piece(p);
                let pc = &mut self.pieces[p as usize];
                pc.cx.remove(k);
                if k < pc.cparam.len() {
                    pc.cparam.remove(k);
                }
                let cx = pc.cx.clone();
                for (j, c2) in cx.into_iter().enumerate() {
                    if self.crosses[c2 as usize].idx != j as u32 {
                        self.save_cross(c2);
                        self.crosses[c2 as usize].idx = j as u32;
                    }
                }
                self.dirty.push(p);
            }
        }
    }

    fn thread_change<F>(&mut self, l: usize, cavity: &[u32], op: F) -> Result<Change, String>
    where
        F: FnOnce(&mut Cdt) -> Result<Change, String>,
    {
        let blocked = cavity.iter().any(|&t| self.layers[l].tflag[t as usize] & T_BLOCK != 0);
        let flag = if blocked { T_BLOCK } else { 0 };
        let incav = |w: &World, t: u32| cavity.contains(&t) && w.layers[l].tflag[t as usize] & T_BLOCK == 0;
        // affected pieces
        let mut aff: Vec<u32> = Vec::new();
        let mut cav_verts: Vec<u32> = Vec::new();
        for &t in cavity {
            let tr = self.layers[l].cdt.tris[t as usize];
            for i in 0..3 {
                for &c in &self.layers[l].cross[tr.e[i] as usize] {
                    aff.push(self.crosses[c as usize].piece);
                }
                let al = self.layers[l].along[tr.e[i] as usize];
                if al != NONE {
                    aff.push(al);
                }
                cav_verts.push(tr.v[i]);
            }
        }
        cav_verts.sort();
        cav_verts.dedup();
        for &v in &cav_verts {
            let o = self.layers[l].vobj[v as usize];
            if o != NONE && (o as usize) < self.obj_pieces.len() {
                for &p in &self.obj_pieces[o as usize] {
                    let pc = &self.pieces[p as usize];
                    if pc.layer as usize == l && (pc.a.vtx == v || pc.b.vtx == v) {
                        aff.push(p);
                    }
                }
            }
        }
        aff.sort();
        aff.dedup();
        // runs: (piece, first crossing index replaced, end index (exclusive), param from, param to)
        let mut runs: Vec<(u32, usize, usize, f64, f64)> = Vec::new();
        let mut along_pieces: Vec<u32> = Vec::new();
        for &p in &aff {
            let pc = &self.pieces[p as usize];
            if !pc.alive || pc.layer as usize != l {
                continue;
            }
            let n = pc.cx.len();
            if n == 0 {
                if pc.along != NONE {
                    along_pieces.push(p);
                }
                continue;
            }
            let mut tris = Vec::with_capacity(n + 1);
            tris.push(if pc.a.vtx != NONE { self.exited_tri(l, pc.cx[0]) } else { NONE });
            for k in 0..n {
                tris.push(self.entered_tri(l, pc.cx[k]));
            }
            if pc.b.vtx == NONE {
                *tris.last_mut().unwrap() = NONE;
            }
            let mut j = 0;
            while j <= n {
                if tris[j] != NONE && incav(self, tris[j]) {
                    let j0 = j;
                    while j <= n && tris[j] != NONE && incav(self, tris[j]) {
                        j += 1;
                    }
                    let j1 = j - 1;
                    let pa = if j0 > 0 { pc.cparam[j0 - 1] } else { 0.0 };
                    let pb = if j1 < n { pc.cparam[j1] } else { (pc.poly.len() - 1) as f64 };
                    runs.push((p, j0, j1, pa, pb));
                } else {
                    j += 1;
                }
            }
        }
        if std::env::var("TOPOR_VALIDATE").is_ok() {
            // every crossing on an internal cavity edge must be inside a run
            let st: Vec<u32> = cavity.to_vec();
            for &p in &aff {
                let pc = &self.pieces[p as usize];
                if !pc.alive || pc.layer as usize != l {
                    continue;
                }
                for (k, &c) in pc.cx.iter().enumerate() {
                    let e = self.crosses[c as usize].edge;
                    let t0 = self.layers[l].etri[e as usize];
                    if st.contains(&t0[0]) && st.contains(&t0[1]) && !runs.iter().any(|r| r.0 == p && k >= r.1 && k < r.2) {
                        let tris = self.piece_tris(p);
                        eprintln!("thread_change: piece {} crossing #{} on internal edge {} not in a run; etri {:?} tris {:?} flags {:?} blocked {} runs {:?}", p, k, e, t0, tris, tris.iter().map(|&t| if t == NONE { 9 } else { self.layers[l].tflag[t as usize] }).collect::<Vec<_>>(), blocked, runs.iter().filter(|r| r.0 == p).collect::<Vec<_>>());
                    }
                }
            }
        }
        // geometry of runs before the change
        let run_polys: Vec<Vec<P2>> = runs.iter().map(|&(p, _, _, a, b)| self.sub_poly(p, a, b)).collect();
        let along_segs: Vec<(P2, P2)> = along_pieces
            .iter()
            .map(|&p| {
                let pc = &self.pieces[p as usize];
                (self.vpos(l, pc.a.vtx), self.vpos(l, pc.b.vtx))
            })
            .collect();
        let ch = op(&mut self.layers[l].cdt)?;
        // side tables for new edges / triangles
        {
            let ne = self.layers[l].cdt.edges.len();
            let nt = self.layers[l].cdt.tris.len();
            let layer = &mut self.layers[l];
            layer.cross.resize(ne, Vec::new());
            layer.along.resize(ne, NONE);
            layer.etri.resize(ne, [NONE, NONE]);
            layer.tflag.resize(nt, flag);
        }
        for &t in &ch.new_tris {
            let tr = self.layers[l].cdt.tris[t as usize];
            for i in 0..3 {
                let e = tr.e[i];
                self.save_etri(l, e);
            }
            set_etri(&mut self.layers[l], t);
        }
        // kill crossings on removed edges
        for &e in &ch.removed_edges {
            let list = self.layers[l].cross[e as usize].clone();
            for c in list {
                self.save_cross(c);
                self.crosses[c as usize].alive = false;
            }
            self.save_cross_list(l, e);
            self.layers[l].cross[e as usize].clear();
            self.save_along(l, e);
            self.layers[l].along[e as usize] = NONE;
        }
        // re-thread runs, last run of a piece first so indices stay valid
        let new_edges: Vec<(u32, P2, P2, u32, u32)> = ch
            .new_edges
            .iter()
            .map(|&e| {
                let ed = self.layers[l].cdt.edges[e as usize];
                (e, self.vpos(l, ed.a), self.vpos(l, ed.b), ed.a, ed.b)
            })
            .collect();
        let mut edge_new: HashMap<u32, Vec<(f64, u32)>> = HashMap::new();
        let mut order: Vec<usize> = (0..runs.len()).collect();
        order.sort_by(|&x, &y| runs[y].0.cmp(&runs[x].0).then(runs[y].1.cmp(&runs[x].1)));
        for ri in order {
            let (p, j0, j1, _, _) = runs[ri];
            let poly = &run_polys[ri];
            let mut hits: Vec<(f64, u32, u32, P2)> = Vec::new();
            for k in 0..poly.len().saturating_sub(1) {
                let (q0, q1) = (poly[k], poly[k + 1]);
                for &(e, ea, eb, va, vb) in &new_edges {
                    if let Some((t, u)) = seg_cross_params(q0, q1, ea, eb) {
                        if u <= 1e-9 || u >= 1.0 - 1e-9 {
                            continue;
                        }
                        let d = q1 - q0;
                        let left = if d.cross(ea - q0) > 0.0 { va } else { vb };
                        hits.push((k as f64 + t, e, left, q0.lerp(q1, t)));
                    }
                }
            }
            hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            // a back-and-forth pair on one edge cancels topologically
            let mut k = 0;
            while k + 1 < hits.len() {
                if hits[k].1 == hits[k + 1].1 {
                    hits.drain(k..k + 2);
                    k = k.saturating_sub(1);
                } else {
                    k += 1;
                }
            }
            let mut new_ids = Vec::new();
            let mut new_par = Vec::new();
            for (t, e, left, q) in hits {
                let cid = self.crosses.len() as u32;
                self.crosses.push(Crossing { piece: p, edge: e, left, alive: true, idx: 0 });
                let ea = new_edges.iter().find(|x| x.0 == e).unwrap().1;
                edge_new.entry(e).or_default().push((q.dist(ea), cid));
                new_ids.push(cid);
                new_par.push(t);
            }
            self.save_piece(p);
            let pc = &mut self.pieces[p as usize];
            pc.cx.splice(j0..j1, new_ids);
            // parameters are approximate (relative to the run's sub-polyline); realisation recomputes them
            pc.cparam.splice(j0..j1, new_par.iter().map(|_| 0.0));
            self.dirty.push(p);
        }
        // along pieces whose edge died
        for (k, &p) in along_pieces.iter().enumerate() {
            let al = self.pieces[p as usize].along;
            if self.layers[l].cdt.edges[al as usize].alive {
                continue;
            }
            let (q0, q1) = along_segs[k];
            let (va, vb) = (self.pieces[p as usize].a.vtx, self.pieces[p as usize].b.vtx);
            self.save_piece(p);
            if let Some(e) = self.layers[l].cdt.find_edge(va, vb) {
                self.layers[l].along[e as usize] = p;
                self.pieces[p as usize].along = e;
            } else {
                let mut hits: Vec<(f64, u32, u32, P2)> = Vec::new();
                for &(e, ea, eb, a, b) in &new_edges {
                    if let Some((t, u)) = seg_cross_params(q0, q1, ea, eb) {
                        if u <= 1e-9 || u >= 1.0 - 1e-9 {
                            continue;
                        }
                        let left = if (q1 - q0).cross(ea - q0) > 0.0 { a } else { b };
                        hits.push((t, e, left, q0.lerp(q1, t)));
                    }
                }
                hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                let mut ids = Vec::new();
                for (_, e, left, q) in hits {
                    let cid = self.crosses.len() as u32;
                    self.crosses.push(Crossing { piece: p, edge: e, left, alive: true, idx: 0 });
                    let ea = new_edges.iter().find(|x| x.0 == e).unwrap().1;
                    edge_new.entry(e).or_default().push((q.dist(ea), cid));
                    ids.push(cid);
                }
                let n = ids.len();
                self.pieces[p as usize].cx = ids;
                self.pieces[p as usize].cparam = vec![0.0; n];
                self.pieces[p as usize].along = NONE;
            }
            self.dirty.push(p);
        }
        let mut keys: Vec<u32> = edge_new.keys().cloned().collect();
        keys.sort();
        for e in keys {
            let mut v = edge_new.remove(&e).unwrap();
            v.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            self.layers[l].cross[e as usize] = v.into_iter().map(|x| x.1).collect();
        }
        // geometry is unchanged: refresh crossing parameters so later changes can re-thread again
        let mut touched: Vec<u32> = runs.iter().map(|r| r.0).chain(along_pieces.iter().cloned()).collect();
        touched.sort();
        touched.dedup();
        for &p in &touched {
            self.trim_terminal_crossings(l, p);
        }
        for &p in &touched {
            self.refresh_cparam(p);
            let cx = self.pieces[p as usize].cx.clone();
            for (i, c) in cx.into_iter().enumerate() {
                if (c as usize) < self.crosses.len() && self.crosses[c as usize].idx != i as u32 {
                    self.save_cross(c);
                    self.crosses[c as usize].idx = i as u32;
                }
            }
        }
        // a wire whose old geometry no longer threads the new triangles consistently (its crossings
        // out of step with the polyline after an earlier change) is taken out: its net is routed again
        for p in touched {
            if let Some(why) = self.thread_fault(l, p) {
                if std::env::var("TOPOR_VALIDATE").is_ok() || std::env::var("TOPOR_THREADDBG").is_ok() {
                    eprintln!("re-threading left piece {} of {} inconsistent ({}): removed", p, self.b.nets[self.pieces[p as usize].net as usize].name, why);
                }
                self.remove_piece(p);
                self.thread_drops += 1;
            }
        }
        Ok(ch)
    }

    /// What is wrong with a live piece's crossings on layer l, if anything: a dead crossing, one missing
    /// from its edge's list, or two in a row that do not share a triangle.
    fn thread_fault(&self, l: usize, p: u32) -> Option<String> {
        let pc = &self.pieces[p as usize];
        if !pc.alive {
            return None;
        }
        for (k, &c) in pc.cx.iter().enumerate() {
            let cr = &self.crosses[c as usize];
            if !cr.alive || !self.layers[l].cdt.edges[cr.edge as usize].alive || !self.layers[l].cross[cr.edge as usize].contains(&c) {
                return Some(format!("crossing #{} stale", k));
            }
            if k + 1 < pc.cx.len() && self.entered_tri(l, c) != self.exited_tri(l, pc.cx[k + 1]) {
                return Some(format!("crossings #{} and #{} apart", k, k + 1));
            }
        }
        None
    }

    // ------------------------------------------------------------------ vias
    pub fn add_via(&mut self, net: u32, pos: P2) -> Result<u32, String> {
        let (d, drill) = self.b.via_size(net);
        let nsig = self.nsig();
        self.add_point(net, pos, d / 2.0, drill, (1 << nsig) - 1, Kind::Via, format!("via:{}", self.b.nets[net as usize].name))
    }

    /// A virtual single-layer terminal (pair fan-out junction): a wire end, never exported as copper.
    pub fn add_junction(&mut self, net: u32, pos: P2, r: f64, l: usize) -> Result<u32, String> {
        self.add_point(net, pos, r, 0.0, 1 << l, Kind::Pad, format!("junction:{}", self.b.nets[net as usize].name))
    }

    fn add_point(&mut self, net: u32, pos: P2, r: f64, drill: f64, layers: u32, kind: Kind, name: String) -> Result<u32, String> {
        let c = &self.b.classes[self.b.nets[net as usize].class];
        let nsig = self.nsig();
        let oid = self.objs.len() as u32;
        let clear = if kind == Kind::Via { self.b.via_clear(net) } else { c.clear };
        let spec = ObjSpec { kind, net, layers, pts: vec![pos], r, clear, center: pos, drill, name, through: kind == Kind::Via };
        let plane = kind == Kind::Via && self.b.nets[net as usize].plane.is_some() && self.b.in_plane(net, pos, r);
        self.objs.push(Obj { s: spec, alive: true, plane, lv: vec![NONE; nsig] });
        self.obj_pieces.push(Vec::new());
        let ip = to_ip(pos);
        for l in 0..nsig {
            if layers & (1 << l) == 0 {
                continue;
            }
            let hint = match self.locate(l, pos) {
                Ok(crate::cdt::Loc::Inside(t)) | Ok(crate::cdt::Loc::OnEdge(t, _)) => t,
                _ => self.layers[l].cdt.any_alive(),
            };
            let cav = match self.layers[l].cdt.insert_cavity(ip, hint)? {
                Ok(c) => c,
                Err(_) => return Err("point on an existing vertex".into()),
            };
            let mut newv = NONE;
            self.thread_change(l, &cav, |cdt| {
                let (v, ch) = cdt.apply_insert(ip, &cav)?;
                newv = v;
                Ok(ch)
            })?;
            let layer = &mut self.layers[l];
            if layer.vobj.len() != newv as usize {
                return Err("vertex numbering out of step".into());
            }
            layer.vobj.push(oid);
            layer.vk.push(0);
            self.objs[oid as usize].lv[l] = newv;
        }
        self.index_obj(oid, true);
        Ok(oid)
    }

    /// Whether a disc of radius r for `net` at pos on layer l clears other copper (routing margin) and
    /// sits in routable space.
    pub fn point_fits(&self, net: u32, pos: P2, r: f64, l: usize) -> bool {
        let me = self.net_attr(net, l);
        let reach = r + self.max_clear + 0.5;
        let mut items = Vec::new();
        self.layers[l].grid.query(pos - p2(reach, reach), pos + p2(reach, reach), &mut items);
        for it in &items {
            match it.kind {
                IT_OBJ | IT_RING => {
                    let o = &self.objs[it.id as usize];
                    if !o.alive {
                        continue;
                    }
                    if o.s.kind == Kind::Keepout && point_in_poly(&o.s.pts, pos) {
                        return false;
                    }
                    let d = self.obj_dist_seg(it.id, it.sub, it.kind, pos, pos);
                    match self.clr(&me, &self.obj_attr(it.id)) {
                        Some(g) => {
                            if d - r < g {
                                return false;
                            }
                        }
                        None => {
                            if d < EPS_OWN * 4.0 {
                                return false;
                            }
                        }
                    }
                }
                _ => {
                    let q = &self.pieces[it.id as usize];
                    if !q.alive {
                        continue;
                    }
                    let k = it.sub as usize;
                    if k + 1 >= q.poly.len() {
                        continue;
                    }
                    let d = seg_point_dist(q.poly[k], q.poly[k + 1], pos) - q.width / 2.0;
                    let need = self.clr(&me, &self.net_attr(q.net, q.layer as usize)).unwrap_or(0.0);
                    if d - r < need {
                        return false;
                    }
                }
            }
        }
        matches!(self.locate(l, pos), Ok(crate::cdt::Loc::Inside(t)) if self.layers[l].tflag[t as usize] & T_BLOCK == 0)
    }

    /// A piece with fixed geometry and no topology (materialised pair members): never re-realised.
    pub fn add_frozen_piece(&mut self, net: u32, l: usize, width: f64, a: Term, b: Term, poly: Vec<P2>) -> u32 {
        let pid = self.pieces.len() as u32;
        self.pieces.push(Piece { net, layer: l as u8, width, a, b, cx: vec![], poly, cparam: vec![], alive: true, along: NONE, overlap: false, problem: None, frozen: true });
        for o in [a.obj, b.obj] {
            self.save_obj_pieces(o);
            self.obj_pieces[o as usize].push(pid);
        }
        self.index_piece(pid);
        pid
    }

    pub fn remove_via(&mut self, o: u32) -> Result<(), String> {
        if !self.objs[o as usize].alive {
            return Ok(());
        }
        let ps = self.obj_pieces[o as usize].clone();
        for p in ps {
            self.remove_piece(p);
        }
        self.index_obj(o, false);
        self.save_obj(o);
        for l in 0..self.nsig() {
            let v = self.objs[o as usize].lv[l];
            if v == NONE {
                continue;
            }
            let star = self.layers[l].cdt.star(v);
            self.thread_change(l, &star, |cdt| cdt.apply_remove(v, &star))?;
        }
        self.objs[o as usize].alive = false;
        Ok(())
    }

    /// Other nets' pieces a straight stretch a-b of `net` would come too close to.
    pub fn stretch_block(&self, l: usize, a: P2, b: P2, net: u32, width: f64, extra: f64) -> Vec<u32> {
        let mut out = Vec::new();
        self.stretch_scan(l, a, b, net, width, extra, |p| {
            if !out.contains(&p) {
                out.push(p);
            }
            true
        });
        out
    }

    /// Visit other nets' pieces a straight stretch a-b of `net` comes too close to; `f` returns false to
    /// stop. Only grid cells near the stretch are looked at (its box can be large and mostly empty).
    fn stretch_scan<F: FnMut(u32) -> bool>(&self, l: usize, a: P2, b: P2, net: u32, width: f64, extra: f64, mut f: F) {
        let me = self.net_attr(net, l);
        let reach = width / 2.0 + self.max_wire_clear + self.max_hw + extra.max(0.0) + 0.01;
        let g = &self.layers[l].grid;
        let lo = p2(a.x.min(b.x) - reach, a.y.min(b.y) - reach);
        let hi = p2(a.x.max(b.x) + reach, a.y.max(b.y) + reach);
        let (x0, y0, x1, y1) = g.range(lo, hi);
        let mut last_net = NONE;
        let mut last_g: Option<f64> = None;
        let d = b - a;
        for y in y0..=y1 {
            // the part of the stretch within reach of this row, widened by reach: the cells to look at
            let (yl, yh) = (g.oy + y as f64 * g.cell - reach, g.oy + (y + 1) as f64 * g.cell + reach);
            let (t0, t1) = if d.y.abs() < 1e-12 {
                if a.y < yl || a.y > yh {
                    continue;
                }
                (0.0, 1.0)
            } else {
                let (ta, tb) = ((yl - a.y) / d.y, (yh - a.y) / d.y);
                let (ta, tb) = (ta.min(tb).max(0.0), ta.max(tb).min(1.0));
                if ta > tb {
                    continue;
                }
                (ta, tb)
            };
            let (xa, xb) = (a.x + d.x * t0, a.x + d.x * t1);
            let cx = |v: f64| (((v - g.ox) / g.cell).floor().max(0.0) as usize).min(g.nx - 1);
            let (xs, xe) = (cx(xa.min(xb) - reach).max(x0), cx(xa.max(xb) + reach).min(x1));
            for x in xs..=xe {
                // segments are indexed in parts, so an entry may come up more than once: harmless here
                for e in &g.segs[y * g.nx + x] {
                    if e.net == net {
                        continue;
                    }
                    let (s0, s1) = (e.a, e.b);
                    if s0.x.max(s1.x) < lo.x || s0.x.min(s1.x) > hi.x || s0.y.max(s1.y) < lo.y || s0.y.min(s1.y) > hi.y {
                        continue;
                    }
                    if e.net != last_net {
                        last_net = e.net;
                        last_g = self.clr(&me, &self.net_attr(e.net, l));
                    }
                    if let Some(gap) = last_g {
                        if seg_seg_dist(a, b, s0, s1) < gap + (width + e.width) / 2.0 + extra && !f(e.id) {
                            return;
                        }
                    }
                }
            }
        }
    }

    /// What blocks a via of `net` at pos: None when fixed copper does, else the pieces to rip (wires
    /// too close, and the pieces holding other nets' vias in the way).
    pub fn via_block(&self, net: u32, pos: P2) -> Option<Vec<u32>> {
        let (r, rd) = { let (d, dr) = self.b.via_size(net); (d / 2.0, dr / 2.0) };
        if self.in_escape(net, pos, r) {
            return None;
        }
        if self.b.via_kept_out(pos, r) {
            return None;
        }
        let me = Attr { kind: Kind::Via, net, class: self.b.nets[net as usize].class as u32, clear: self.b.via_clear(net) };
        let reach = r + self.max_clear + self.max_hw + 0.01;
        let mut items = Vec::new();
        let mut out: Vec<u32> = Vec::new();
        for l in 0..self.nsig() {
            self.layers[l].grid.query(pos - p2(reach, reach), pos + p2(reach, reach), &mut items);
            for it in &items {
                match it.kind {
                    IT_OBJ | IT_RING => {
                        let o = &self.objs[it.id as usize];
                        if !o.alive {
                            continue;
                        }
                        let d = self.obj_dist_seg(it.id, it.sub, it.kind, pos, pos);
                        let hole = o.s.drill > 0.0 && it.kind == IT_OBJ && o.s.center.dist(pos) - rd - o.s.drill / 2.0 < self.b.hole_to_hole + MARGIN;
                        let inside_keep = o.s.kind == Kind::Keepout && point_in_poly(&o.s.pts, pos);
                        let blocked = match self.clr(&me, &self.obj_attr(it.id)) {
                            Some(g) => d - r < g || (o.s.kind != Kind::Outline && o.s.kind != Kind::Keepout && d - rd < self.b.hole_clear + MARGIN),
                            None => (o.s.kind == Kind::Pad && d < r + EPS_OWN * 4.0) || (o.s.kind == Kind::Via && d < r + o.s.r) || self.tie_hole(it.id, d, rd),
                        };
                        if hole || inside_keep || blocked {
                            if o.s.kind == Kind::Via && o.s.net != net {
                                let ps: Vec<u32> = self.obj_pieces[it.id as usize].iter().cloned().filter(|&p| self.pieces[p as usize].alive).collect();
                                if ps.is_empty() {
                                    return None;
                                }
                                for p in ps {
                                    if !out.contains(&p) {
                                        out.push(p);
                                    }
                                }
                            } else {
                                return None;
                            }
                        }
                    }
                    _ => {
                        let q = &self.pieces[it.id as usize];
                        if !q.alive || out.contains(&it.id) {
                            continue;
                        }
                        let k2 = it.sub as usize;
                        if k2 + 1 >= q.poly.len() {
                            continue;
                        }
                        let d = seg_point_dist(q.poly[k2], q.poly[k2 + 1], pos) - q.width / 2.0;
                        let need = if q.net == net {
                            0.0
                        } else {
                            match self.clr(&me, &self.net_attr(q.net, q.layer as usize)) {
                                Some(g) => g,
                                None => continue,
                            }
                        };
                        if d - r < need || d - rd < self.b.hole_clear + MARGIN {
                            if q.net == net {
                                return None;
                            }
                            out.push(it.id);
                        }
                    }
                }
            }
        }
        Some(out)
    }

    /// Whether a straight wire stretch a-b of `net` clears other nets' wires (with `extra` slack).
    pub fn stretch_clear(&self, l: usize, a: P2, b: P2, net: u32, width: f64, extra: f64) -> bool {
        let mut clear = true;
        self.stretch_scan(l, a, b, net, width, extra, |_| {
            clear = false;
            false
        });
        clear
    }

    // ------------------------------------------------------------------ exact checks
    pub fn obj_dist_seg(&self, o: u32, sub: u32, kind: u8, a: P2, b: P2) -> f64 {
        let s = &self.objs[o as usize].s;
        if kind == IT_RING {
            let n = s.pts.len();
            let (u, w) = (s.pts[sub as usize], s.pts[(sub as usize + 1) % n]);
            return seg_seg_dist(a, b, u, w);
        }
        seg_hull_dist(&s.pts, a, b) - s.r
    }

    /// Clearance violations of one piece against everything else.
    pub fn check_piece(&self, p: u32, out: &mut Vec<Viol>) {
        let pc = &self.pieces[p as usize];
        if !pc.alive {
            return;
        }
        let l = pc.layer as usize;
        let me = self.net_attr(pc.net, pc.layer as usize);
        let hw = pc.width / 2.0;
        let reach = hw + self.max_clear + 0.5;
        let mut items = Vec::new();
        for k in 0..pc.poly.len().saturating_sub(1) {
            let (a, b) = (pc.poly[k], pc.poly[k + 1]);
            let lo = p2(a.x.min(b.x) - reach, a.y.min(b.y) - reach);
            let hi = p2(a.x.max(b.x) + reach, a.y.max(b.y) + reach);
            self.layers[l].grid.query(lo, hi, &mut items);
            // a segment is indexed in parts: look at each once
            items.sort_unstable_by_key(|it| (it.kind, it.id, it.sub));
            items.dedup();
            for it in &items {
                match it.kind {
                    IT_OBJ | IT_RING => {
                        let o = &self.objs[it.id as usize];
                        if !o.alive {
                            continue;
                        }
                        let oa = self.obj_attr(it.id);
                        let need = match self.clr(&me, &oa) {
                            Some(g) => g - MARGIN,
                            None => continue,
                        };
                        let mut d = self.obj_dist_seg(it.id, it.sub, it.kind, a, b) - hw;
                        if o.s.kind == Kind::Keepout && it.kind == IT_RING && (point_in_poly(&o.s.pts, a) || point_in_poly(&o.s.pts, b)) {
                            d = -1.0;
                        }
                        if d < need - 1e-7 {
                            out.push(Viol { what: format!("track-{:?}", o.s.kind), layer: l as i32, at: a.lerp(b, 0.5), dist: d, need, pieces: vec![p], objs: vec![it.id] });
                        }
                    }
                    _ => {
                        if it.id == p {
                            continue;
                        }
                        let q = &self.pieces[it.id as usize];
                        if !q.alive || q.net == pc.net {
                            continue;
                        }
                        let need = match self.clr(&me, &self.net_attr(q.net, q.layer as usize)) {
                            Some(g) => g - MARGIN,
                            None => continue,
                        };
                        let k2 = it.sub as usize;
                        if k2 + 1 >= q.poly.len() {
                            continue;
                        }
                        let d = seg_seg_dist(a, b, q.poly[k2], q.poly[k2 + 1]) - hw - q.width / 2.0;
                        if d < need - 1e-7 {
                            out.push(Viol { what: "track-track".into(), layer: l as i32, at: a.lerp(b, 0.5), dist: d, need, pieces: vec![p, it.id], objs: vec![] });
                        }
                    }
                }
            }
        }
    }

    /// Clearance of a via (copper on every signal layer, drill to holes) at pos for net.
    pub fn check_via_at(&self, net: u32, pos: P2, skip_obj: u32, out: &mut Vec<Viol>, early: bool) {
        let (r, rd) = { let (d, dr) = self.b.via_size(net); (d / 2.0, dr / 2.0) };
        let me = Attr { kind: Kind::Via, net, class: self.b.nets[net as usize].class as u32, clear: self.b.via_clear(net) };
        if self.b.via_kept_out(pos, r) {
            out.push(Viol { what: "via-rule-area".into(), layer: -1, at: pos, dist: -1.0, need: 0.0, pieces: vec![], objs: vec![] });
            if early {
                return;
            }
        }
        let reach = r + self.max_clear + self.max_hw + 0.01;
        let lo = pos - p2(reach, reach);
        let hi = pos + p2(reach, reach);
        let mut items = Vec::new();
        for l in 0..self.nsig() {
            self.layers[l].grid.query(lo, hi, &mut items);
            for it in &items {
                match it.kind {
                    IT_OBJ | IT_RING => {
                        if it.id == skip_obj {
                            continue;
                        }
                        let o = &self.objs[it.id as usize];
                        if !o.alive {
                            continue;
                        }
                        let oa = self.obj_attr(it.id);
                        let d = self.obj_dist_seg(it.id, it.sub, it.kind, pos, pos);
                        // drill to drill, any net
                        if o.s.drill > 0.0 && it.kind == IT_OBJ {
                            let dh = o.s.center.dist(pos) - rd - o.s.drill / 2.0;
                            let need = self.b.hole_to_hole;
                            if dh < need - 1e-7 {
                                out.push(Viol { what: "hole-hole".into(), layer: -1, at: pos, dist: dh, need, pieces: vec![], objs: vec![it.id] });
                                if early {
                                    return;
                                }
                            }
                        }
                        if self.clr(&me, &oa).is_none() && self.tie_hole(it.id, d, rd) {
                            out.push(Viol { what: "hole-net-tie".into(), layer: l as i32, at: pos, dist: d - rd, need: self.b.hole_clear, pieces: vec![], objs: vec![it.id] });
                            if early {
                                return;
                            }
                        }
                        if let Some(g) = self.clr(&me, &oa) {
                            let need = g - MARGIN;
                            let mut dd = d - r;
                            if o.s.kind == Kind::Keepout && point_in_poly(&o.s.pts, pos) {
                                dd = -1.0;
                            }
                            if dd < need - 1e-7 {
                                out.push(Viol { what: format!("via-{:?}", o.s.kind), layer: l as i32, at: pos, dist: dd, need, pieces: vec![], objs: vec![it.id] });
                                if early {
                                    return;
                                }
                            }
                            if o.s.kind != Kind::Outline && o.s.kind != Kind::Keepout {
                                let dh = d - rd;
                                let need = self.b.hole_clear;
                                if dh < need - 1e-7 {
                                    out.push(Viol { what: "hole-copper".into(), layer: l as i32, at: pos, dist: dh, need, pieces: vec![], objs: vec![it.id] });
                                    if early {
                                        return;
                                    }
                                }
                            }
                        }
                    }
                    _ => {
                        let q = &self.pieces[it.id as usize];
                        if !q.alive {
                            continue;
                        }
                        let k2 = it.sub as usize;
                        if k2 + 1 >= q.poly.len() {
                            continue;
                        }
                        let d = seg_point_dist(q.poly[k2], q.poly[k2 + 1], pos) - q.width / 2.0;
                        if q.net == net {
                            continue;
                        }
                        let need = match self.clr(&me, &self.net_attr(q.net, q.layer as usize)) {
                            Some(g) => g - MARGIN,
                            None => continue,
                        };
                        if d - r < need - 1e-7 {
                            out.push(Viol { what: "via-track".into(), layer: l as i32, at: pos, dist: d - r, need, pieces: vec![it.id], objs: vec![] });
                            if early {
                                return;
                            }
                        }
                        if d - rd < self.b.hole_clear - 1e-7 {
                            out.push(Viol { what: "hole-track".into(), layer: l as i32, at: pos, dist: d - rd, need: self.b.hole_clear, pieces: vec![it.id], objs: vec![] });
                            if early {
                                return;
                            }
                        }
                    }
                }
            }
        }
    }

    /// Escape lanes for pads of chips (6+ pads) whose neighbour of another net is within 0.8 mm: from
    /// the pad's outer end 0.8 mm outwards along its long axis, away from the chip centre.
    fn build_escapes(&mut self) {
        let mut per_net: crate::fx::FxMap<u32, usize> = Default::default();
        let mut parts: crate::fx::FxMap<String, Vec<u32>> = Default::default();
        for (i, o) in self.objs.iter().enumerate() {
            if o.s.kind != Kind::Pad || o.s.net == NONE || o.s.name.starts_with("junction:") {
                continue;
            }
            *per_net.entry(o.s.net).or_default() += 1;
            parts.entry(o.s.name.split('.').next().unwrap_or("").to_string()).or_default().push(i as u32);
        }
        for (_, pads) in parts {
            if pads.len() < 6 {
                continue;
            }
            let cc = pads.iter().fold(p2(0.0, 0.0), |a, &o| a + self.objs[o as usize].s.center) * (1.0 / pads.len() as f64);
            for &o in &pads {
                let s = &self.objs[o as usize].s;
                if per_net.get(&s.net).cloned().unwrap_or(0) < 2 && self.b.nets[s.net as usize].plane.is_none() {
                    continue;
                }
                let fine = pads.iter().any(|&q| q != o && self.objs[q as usize].s.net != s.net && self.objs[q as usize].s.center.dist(s.center) < 0.8);
                if !fine || s.pts.len() < 3 {
                    continue;
                }
                // long axis from the hull's longest edge
                let n = s.pts.len();
                let (mut ax, mut best) = (p2(1.0, 0.0), 0.0);
                for k in 0..n {
                    let e = s.pts[(k + 1) % n] - s.pts[k];
                    if e.norm() > best {
                        best = e.norm();
                        ax = e.unit();
                    }
                }
                if ax.dot(s.center - cc) < 0.0 {
                    ax = -ax;
                }
                let (mut far, mut hw) = (0.0f64, 0.0f64);
                for q in &s.pts {
                    far = far.max((*q - s.center).dot(ax));
                    hw = hw.max((*q - s.center).dot(ax.perp()).abs());
                }
                let a = s.center + ax * far;
                let b = a + ax * 0.8;
                let id = self.escapes.len() as u32;
                self.escapes.push((s.net, a, b, hw, o));
                let (lo, hi) = (p2(a.x.min(b.x) - 1.0, a.y.min(b.y) - 1.0), p2(a.x.max(b.x) + 1.0, a.y.max(b.y) + 1.0));
                for gx in lo.x.floor() as i32..=hi.x.floor() as i32 {
                    for gy in lo.y.floor() as i32..=hi.y.floor() as i32 {
                        self.esc_grid.entry((gx, gy)).or_default().push(id);
                    }
                }
            }
        }
    }

    /// Whether a via of `net` (radius r) at pos would sit in another net's pad escape lane.
    pub fn in_escape(&self, net: u32, pos: P2, r: f64) -> bool {
        // experimental (TOPOR_ESCAPES): on DUT it moved failures around rather than removing them
        if !self.use_escapes || self.b.nets.get(net as usize).map_or(true, |n| n.plane.is_some()) {
            return false;
        }
        let key = (pos.x.floor() as i32, pos.y.floor() as i32);
        match self.esc_grid.get(&key) {
            // only pads still waiting for their escape: once a wire leaves the pad the lane is its own
            Some(ids) => ids.iter().any(|&i| {
                let (n, a, b, hw, o) = self.escapes[i as usize];
                n != net && seg_point_dist(a, b, pos) < r + hw + 0.1 && !self.obj_pieces[o as usize].iter().any(|&p| self.pieces[p as usize].alive)
            }),
            None => false,
        }
    }

    /// Whether a via of `net` fits at pos with the routing margin.
    pub fn via_fits(&self, net: u32, pos: P2) -> bool {
        if self.in_escape(net, pos, self.b.via_size(net).0 / 2.0) {
            return false;
        }
        let mut v = Vec::new();
        self.check_via_margin(net, pos, NONE, &mut v);
        v.is_empty()
    }

    /// Whether a via fits entirely inside its own net's pad `pad` (via-in-pad) and clears everything else.
    pub fn via_fits_in_pad(&self, net: u32, pos: P2, pad: u32) -> bool {
        let ob = &self.objs[pad as usize];
        if ob.s.pts.len() < 3 || !point_in_hull(&ob.s.pts, pos) {
            return false;
        }
        // distance from pos to the pad's copper boundary
        let n = ob.s.pts.len();
        let mut edge = f64::INFINITY;
        for i in 0..n {
            edge = edge.min(seg_point_dist(ob.s.pts[i], ob.s.pts[(i + 1) % n], pos));
        }
        if edge + ob.s.r < self.b.via_size(net).0 / 2.0 + 0.02 {
            return false;
        }
        let mut v = Vec::new();
        self.check_via_margin(net, pos, pad, &mut v);
        v.is_empty()
    }

    fn check_via_margin(&self, net: u32, pos: P2, allow: u32, out: &mut Vec<Viol>) {
        // same as check_via_at but demanding the routing margin
        let (r, rd) = { let (d, dr) = self.b.via_size(net); (d / 2.0, dr / 2.0) };
        let me = Attr { kind: Kind::Via, net, class: self.b.nets[net as usize].class as u32, clear: self.b.via_clear(net) };
        if self.b.via_kept_out(pos, r) {
            out.push(Viol { what: "via-rule-area".into(), layer: -1, at: pos, dist: 0., need: 0., pieces: vec![], objs: vec![] });
            return;
        }
        let reach = r + self.max_clear + 0.5;
        let lo = pos - p2(reach, reach);
        let hi = pos + p2(reach, reach);
        let mut items = Vec::new();
        for l in 0..self.nsig() {
            self.layers[l].grid.query(lo, hi, &mut items);
            for it in &items {
                match it.kind {
                    IT_OBJ | IT_RING => {
                        let o = &self.objs[it.id as usize];
                        if !o.alive || it.id == allow {
                            continue;
                        }
                        let d = self.obj_dist_seg(it.id, it.sub, it.kind, pos, pos);
                        if o.s.drill > 0.0 && it.kind == IT_OBJ && o.s.center.dist(pos) - rd - o.s.drill / 2.0 < self.b.hole_to_hole + MARGIN {
                            out.push(Viol { what: "hole-hole".into(), layer: -1, at: pos, dist: 0., need: 0., pieces: vec![], objs: vec![it.id] });
                            return;
                        }
                        if o.s.kind == Kind::Keepout && point_in_poly(&o.s.pts, pos) {
                            out.push(Viol { what: "via-keepout".into(), layer: l as i32, at: pos, dist: 0., need: 0., pieces: vec![], objs: vec![it.id] });
                            return;
                        }
                        let oa = self.obj_attr(it.id);
                        match self.clr(&me, &oa) {
                            Some(g) => {
                                if d - r < g || (o.s.kind != Kind::Outline && o.s.kind != Kind::Keepout && d - rd < self.b.hole_clear + MARGIN) {
                                    out.push(Viol { what: "via-obj".into(), layer: l as i32, at: pos, dist: d, need: g, pieces: vec![], objs: vec![it.id] });
                                    return;
                                }
                            }
                            None => {
                                if self.tie_hole(it.id, d, rd) {
                                    out.push(Viol { what: "hole-net-tie".into(), layer: l as i32, at: pos, dist: d - rd, need: self.b.hole_clear, pieces: vec![], objs: vec![it.id] });
                                    return;
                                }
                                // own-net pad: no via inside small pads; keep the vertex off its hull
                                if d < r + EPS_OWN * 4.0 && o.s.kind == Kind::Pad {
                                    out.push(Viol { what: "via-own-pad".into(), layer: l as i32, at: pos, dist: d, need: 0., pieces: vec![], objs: vec![it.id] });
                                    return;
                                }
                                if o.s.kind == Kind::Via && d < r + o.s.r {
                                    out.push(Viol { what: "via-own-via".into(), layer: l as i32, at: pos, dist: d, need: 0., pieces: vec![], objs: vec![it.id] });
                                    return;
                                }
                            }
                        }
                    }
                    _ => {
                        let q = &self.pieces[it.id as usize];
                        if !q.alive {
                            continue;
                        }
                        let k2 = it.sub as usize;
                        if k2 + 1 >= q.poly.len() {
                            continue;
                        }
                        let d = seg_point_dist(q.poly[k2], q.poly[k2 + 1], pos) - q.width / 2.0;
                        let need = if q.net == net { 0.0 } else {
                            match self.clr(&me, &self.net_attr(q.net, q.layer as usize)) {
                                Some(g) => g,
                                None => continue,
                            }
                        };
                        if d - r < need || d - rd < self.b.hole_clear + MARGIN {
                            out.push(Viol { what: "via-track".into(), layer: l as i32, at: pos, dist: d, need, pieces: vec![it.id], objs: vec![] });
                            return;
                        }
                    }
                }
            }
        }
    }

    pub fn check_all(&self) -> Vec<Viol> {
        let mut out = Vec::new();
        for p in 0..self.pieces.len() as u32 {
            self.check_piece(p, &mut out);
        }
        for (oi, o) in self.objs.iter().enumerate() {
            if o.alive && o.s.kind == Kind::Via {
                self.check_via_at(o.s.net, o.s.center, oi as u32, &mut out, false);
            }
        }
        // a wire under its net's least width is a neck-down: one straight segment from a pad the rule covers
        for (p, pc) in self.pieces.iter().enumerate() {
            if !pc.alive || pc.width >= self.b.nets[pc.net as usize].wmin[pc.layer as usize] - 1e-9 {
                continue;
            }
            let at_neck = |o: u32| self.b.neck_width.get(&self.objs[o as usize].s.name).map_or(false, |&e| pc.width >= e - 1e-9);
            if pc.poly.len() != 2 || !(at_neck(pc.a.obj) || at_neck(pc.b.obj)) {
                out.push(Viol { what: "neck-width".into(), layer: pc.layer as i32, at: pc.poly.first().cloned().unwrap_or(p2(0., 0.)), dist: pc.width, need: self.b.nets[pc.net as usize].wmin[pc.layer as usize], pieces: vec![p as u32], objs: vec![] });
            }
        }
        // each pair reported twice for track-track; dedupe roughly
        out
    }
}

fn set_etri(layer: &mut Layer, t: u32) {
    let tr = layer.cdt.tris[t as usize];
    for i in 0..3 {
        let e = tr.e[i];
        let ed = layer.cdt.edges[e as usize];
        let u = tr.v[i];
        let idx = if u == ed.a { 0 } else { 1 };
        layer.etri[e as usize][idx] = t;
    }
}

pub fn tri_edge_idx(t: &crate::cdt::Tri, e: u32) -> usize {
    tri_edge_index(t, e)
}

/// Merge objects that overlap on a shared layer (conservative convex hulls) so fixed edges never cross.
/// Whether a through pad joins its net's plane: its centre, or for a big pad any part of its copper
/// well inside (a pour that ends across a shunt's terminal still reaches it).
fn pad_in_plane(b: &Board, s: &ObjSpec) -> bool {
    if b.in_plane(s.net, s.center, s.r.min(0.3)) {
        return true;
    }
    // probes 0.35 mm inside the copper: round the core points of a big round or oval pad, in from the
    // corners of a polygon pad
    let mut probes = Vec::new();
    for &q in &s.pts {
        if s.r > 0.35 {
            for k in 0..8 {
                let t = k as f64 * std::f64::consts::FRAC_PI_4;
                probes.push(q + p2(t.cos(), t.sin()) * (s.r - 0.35));
            }
        }
        if s.pts.len() >= 3 && q.dist(s.center) > 0.7 {
            probes.push(q + (s.center - q) * (0.35 / q.dist(s.center)));
        }
    }
    probes.into_iter().any(|q| b.in_plane(s.net, q, 0.3))
}

/// Distance between two objects' cores (point, segment or convex hull); their copper is this less both radii apart.
pub fn core_dist(a: &ObjSpec, b: &ObjSpec) -> f64 {
    let mut d = f64::INFINITY;
    for &q in &a.pts {
        d = d.min(point_hull_dist(&b.pts, q));
    }
    for &q in &b.pts {
        d = d.min(point_hull_dist(&a.pts, q));
    }
    if a.pts.len() >= 2 && b.pts.len() >= 2 {
        let na = if a.pts.len() == 2 { 1 } else { a.pts.len() };
        let nb = if b.pts.len() == 2 { 1 } else { b.pts.len() };
        for k in 0..na {
            for m in 0..nb {
                d = d.min(seg_seg_dist(a.pts[k], a.pts[(k + 1) % a.pts.len()], b.pts[m], b.pts[(m + 1) % b.pts.len()]));
            }
        }
    }
    d
}

fn merge_overlaps(specs: Vec<ObjSpec>, nsig: usize) -> (Vec<ObjSpec>, Vec<usize>, Vec<String>) {
    let mut warnings = Vec::new();
    let n = specs.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut Vec<usize>, x: usize) -> usize {
        let mut r = x;
        while p[r] != r {
            r = p[r];
        }
        let mut y = x;
        while p[y] != r {
            let nx = p[y];
            p[y] = r;
            y = nx;
        }
        r
    }
    // bucket by coarse cells
    let mut buckets: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    let bbox = |s: &ObjSpec| -> (P2, P2) {
        let mut lo = p2(f64::INFINITY, f64::INFINITY);
        let mut hi = p2(f64::NEG_INFINITY, f64::NEG_INFINITY);
        for q in &s.pts {
            lo = p2(lo.x.min(q.x), lo.y.min(q.y));
            hi = p2(hi.x.max(q.x), hi.y.max(q.y));
        }
        (lo - p2(s.r, s.r), hi + p2(s.r, s.r))
    };
    for (i, s) in specs.iter().enumerate() {
        let (lo, hi) = bbox(s);
        for x in (lo.x / 2.0).floor() as i64..=(hi.x / 2.0).floor() as i64 {
            for y in (lo.y / 2.0).floor() as i64..=(hi.y / 2.0).floor() as i64 {
                buckets.entry((x, y)).or_default().push(i);
            }
        }
    }
    let _ = nsig;
    let mut pairs = Vec::new();
    for v in buckets.values() {
        for a in 0..v.len() {
            for b in a + 1..v.len() {
                pairs.push((v[a].min(v[b]), v[a].max(v[b])));
            }
        }
    }
    pairs.sort();
    pairs.dedup();
    for (i, j) in pairs {
        let (a, b) = (&specs[i], &specs[j]);
        if a.layers & b.layers == 0 {
            continue;
        }
        // keepouts touching copper stay separate (a mask opening beside a solder jumper's pads)
        if (a.kind == Kind::Keepout) != (b.kind == Kind::Keepout) {
            continue;
        }
        // core hulls closer than their radii (touching or overlapping copper)?
        let mut d = f64::INFINITY;
        let ring = |s: &ObjSpec| matches!(s.kind, Kind::Keepout | Kind::Outline);
        if ring(a) || ring(b) {
            // polygons: overlap if any point of one is inside the other or edges intersect
            let pa = &a.pts;
            let pb = &b.pts;
            let mut ov = pa.iter().any(|&q| pb.len() >= 3 && point_in_poly(pb, q)) || pb.iter().any(|&q| pa.len() >= 3 && point_in_poly(pa, q));
            if !ov {
                for k in 0..pa.len() {
                    let (u, w) = (pa[k], pa[(k + 1) % pa.len()]);
                    for m in 0..pb.len() {
                        let (x, y) = (pb[m], pb[(m + 1) % pb.len()]);
                        if seg_seg_dist(u, w, x, y) < 1e-4 {
                            ov = true;
                        }
                    }
                }
            }
            if ov {
                d = 0.0;
            }
        } else {
            // core distance; copper touches when it is under ra + rb
            d = core_dist(a, b);
            if d >= a.r + b.r - 1e-6 && d > 1e-4 {
                d = f64::INFINITY;
            } else {
                d = 0.0;
            }
        }
        // zone pieces stay apart: an island's pieces and the pads it reaches are tied, not hulled
        if d == 0.0 && (fixed_copper(&a.name) || fixed_copper(&b.name)) {
            if a.net != b.net {
                warnings.push(format!("zone copper overlaps another net: {} and {}", a.name, b.name));
            }
            continue;
        }
        if d == 0.0 {
            if a.net != b.net {
                warnings.push(format!("overlap of different nets: {} ({:?}) and {} ({:?}); merged as one obstacle", a.name, a.kind, b.name, b.kind));
            }
            let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
            if ri != rj {
                parent[rj] = ri;
            }
        }
    }
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }
    let mut keys: Vec<usize> = groups.keys().cloned().collect();
    keys.sort();
    let mut out = Vec::new();
    let mut map = vec![0usize; n];
    for k in keys {
        let g = &groups[&k];
        for &i in g {
            map[i] = out.len();
        }
        if g.len() == 1 {
            out.push(specs[g[0]].clone());
            continue;
        }
        // merged: convex hull of all copper, as a polygon (radius folded in by sampling discs)
        let mut pts = Vec::new();
        let mut layers = 0;
        let mut clear: f64 = 0.0;
        let mut drill: f64 = 0.0;
        let mut through = false;
        let mut names = Vec::new();
        let net = specs[g[0]].net;
        let mut kind = specs[g[0]].kind;
        for &i in g {
            let s = &specs[i];
            layers |= s.layers;
            clear = clear.max(s.clear);
            drill = drill.max(s.drill);
            through |= s.through;
            names.push(s.name.clone());
            if s.kind == Kind::Pad {
                kind = Kind::Pad;
            }
            for &q in &s.pts {
                if s.r > 0.0 {
                    for m in 0..16 {
                        let t = m as f64 * std::f64::consts::TAU / 16.0;
                        // circumscribed polygon of the disc
                        let rr = s.r / (std::f64::consts::PI / 16.0).cos();
                        pts.push(q + p2(t.cos(), t.sin()) * rr);
                    }
                } else {
                    pts.push(q);
                }
            }
        }
        let hull = convex_hull(&pts);
        let center = specs[g[0]].center;
        let nets_differ = g.iter().any(|&i| specs[i].net != net);
        out.push(ObjSpec {
            kind,
            net: if nets_differ { NONE } else { net },
            layers,
            pts: hull,
            r: 0.0,
            clear,
            center,
            drill,
            name: names.join("+"),
            through,
        });
    }
    (out, map, warnings)
}

/// Copper of a net laid before routing besides its pads: a zone's fill pieces and neck-down stubs. Tied to the pads
/// it touches, never hulled with them.
pub fn fixed_copper(name: &str) -> bool {
    name.starts_with("zone:") || name.starts_with("neck:")
}

/// World::clr on the board's rules alone.
pub fn clr_in(b: &Board, x: &Attr, y: &Attr) -> Option<f64> {
    if x.clear < 0.0 || y.clear < 0.0 {
        return None;
    }
    if x.net != NONE && x.net == y.net {
        return None;
    }
    // a pair's ribbon (its two members as one wire while it is routed) and the pair's own copper may touch: it
    // starts at the members' pads and runs over their fixed copper; the members keep their gap to each other
    if b.ribbon_touch && x.net != NONE && y.net != NONE {
        let (a, c) = (x.net as usize, y.net as usize);
        if a < b.family.len() && c < b.family.len() && b.family[a] == b.family[c] && (b.nets[a].ribbon.is_some() || b.nets[c].ribbon.is_some()) {
            return None;
        }
    }
    let c = if x.kind == Kind::Outline || y.kind == Kind::Outline {
        b.edge_clear
    } else if x.kind == Kind::Hole || y.kind == Kind::Hole {
        if x.kind == Kind::Hole { x.clear } else { y.clear }
    } else if x.kind == Kind::Keepout || y.kind == Kind::Keepout {
        if x.kind == Kind::Keepout { x.clear } else { y.clear }
    } else {
        let mut c = x.clear.max(y.clear);
        if x.class != NONE && x.class == y.class {
            c = c.max(b.class_pair[x.class as usize]);
        }
        // the two members of a differential pair keep exactly their coupling gap
        if x.net != NONE && (x.net as usize) < b.partner.len() && b.partner[x.net as usize] == y.net {
            c = c.max(b.pair_gap[x.net as usize]);
        }
        c
    };
    Some(c + MARGIN)
}
