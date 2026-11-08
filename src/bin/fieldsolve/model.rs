//! Voxel model: turns a layered layout (copper sheets, dielectric slabs, via barrels, board outline) and
//! lumped elements into per-edge FDTD update coefficients on a [`Grid`].
//!
//! Conventions (as in openEMS PCB models): copper is a zero-thickness perfect conductor on the mesh plane
//! of its layer; an in-plane E edge is metal when its midpoint lies in copper. A via barrel is a solid
//! perfect-conductor column (z edges at nodes within its radius). Dielectric fills cells whose centre lies
//! inside the board outline. Lumped elements (ports, resistors, capacitors) occupy a box of edges along one
//! axis: `k` edges in series by `m` in parallel, each edge carrying `R m / k` and `C k / m`.

use crate::geom::{BBox, P2, Prim, Shape};
use crate::mesh::{range_in, snap, Grid};

pub const C0: f64 = 299_792_458.0;
pub const MU0: f64 = 1.256_637_062_12e-6;
pub const EPS0: f64 = 8.854_187_812_8e-12;
pub const ETA0: f64 = 376.730_313_668;

/// Copper on one layer at height `z` (mm).
#[derive(Clone, Debug)]
pub struct Sheet {
    pub name: String,
    pub z: f64,
    pub prims: Vec<Prim>,
}

/// A dielectric layer between `z0` and `z1` (mm).
#[derive(Clone, Copy, Debug)]
pub struct Slab {
    pub z0: f64,
    pub z1: f64,
    pub eps: f64,
    pub tand: f64,
}

/// A plated barrel (via or pin) between `z0` and `z1` (mm).
#[derive(Clone, Copy, Debug)]
pub struct Barrel {
    pub c: P2,
    pub r: f64,
    pub z0: f64,
    pub z1: f64,
}

#[derive(Clone, Debug, Default)]
pub struct Layout {
    pub sheets: Vec<Sheet>,
    pub slabs: Vec<Slab>,
    pub barrels: Vec<Barrel>,
    /// Board outline loops (even-odd); empty means dielectric everywhere in x/y.
    pub outline: Vec<Vec<P2>>,
    /// Copper is rasterised only inside this box (mm): the fine-meshed region. Copper outside it would land
    /// on graded coarse cells and merge across gaps.
    pub copper_window: Option<BBox>,
    /// 3-D dielectric objects (enclosure, FPC substrate): they replace the slab material in the cells they
    /// fill, later ones winning (bodies.rs).
    pub bodies: Vec<crate::bodies::Body>,
}

impl Layout {
    pub fn sheet(&self, name: &str) -> Option<&Sheet> {
        self.sheets.iter().find(|s| s.name == name)
    }
}

/// A lumped two-terminal element on the edges of a box. `lo`/`hi` are box corners in mm; along `axis` the
/// box runs from terminal to terminal, across the other two axes it may be a line, a sheet or a volume.
/// `plus_at_hi`: the + terminal is at the high coordinate end along `axis`.
#[derive(Clone, Debug)]
pub struct Lumped {
    pub name: String,
    pub axis: usize,
    pub lo: [f64; 3],
    pub hi: [f64; 3],
    pub plus_at_hi: bool,
    /// Resistance in ohms (0 = none: then it must have `c` or be a pure probe).
    pub r: f64,
    /// Parallel capacitance in farads.
    pub c: f64,
    /// A port: it is measured, and excited when chosen.
    pub port: bool,
}

/// Runtime form of a port: where its source acts and how its voltage is measured.
#[derive(Clone, Debug)]
pub struct PortRt {
    pub name: String,
    pub r: f64,
    /// (component, edge index, gain): `E += gain * Vs(t)` after each E update.
    pub src: Vec<(u8, u32, f32)>,
    /// (component, edge index, weight): `V = sum weight * E`.
    pub probe: Vec<(u8, u32, f32)>,
}

/// CPML tables for one axis. Positions along the axis map to a slab index `s` (or `NONE`); `be/ce` are the
/// recursive-convolution coefficients at node positions (E derivatives), `bh/ch` at half positions.
#[derive(Clone, Debug, Default)]
pub struct PmlAxis {
    pub map: Vec<u32>,
    pub slab: Vec<u32>,
    pub be: Vec<f32>,
    pub ce: Vec<f32>,
    pub bh: Vec<f32>,
    pub ch: Vec<f32>,
}

impl PmlAxis {
    pub fn len(&self) -> usize {
        self.slab.len()
    }
}

pub const NONE: u32 = u32::MAX;

/// A linear field functional: `sum weight * F[comp][idx]`, comp 0..3 = Ex..Ez, 3..6 = Hx..Hz. A probe is
/// either all-E (a voltage) or all-H (a current loop).
#[derive(Clone, Debug)]
pub struct ProbeRt {
    pub name: String,
    pub h: bool,
    pub terms: Vec<(u8, u32, f32)>,
}

impl Problem {
    fn lines_mm(&self, a: usize) -> Vec<f64> {
        self.grid.axis(a).iter().map(|v| v * 1e3).collect()
    }
    fn dual(&self, a: usize, i: usize) -> f64 {
        let l = self.grid.axis(a);
        let lo = if i > 0 { l[i] - l[i - 1] } else { 0.0 };
        let hi = if i + 1 < l.len() { l[i + 1] - l[i] } else { 0.0 };
        0.5 * (lo + hi)
    }
    /// Voltage probe: integral of E along `axis` from `from` to `to` (mm), on the mesh line nearest the
    /// transverse position of `from`: V = phi(from) - phi(to).
    pub fn add_vprobe(&mut self, name: &str, axis: usize, from: [f64; 3], to: f64) {
        let mut at = [0usize; 3];
        for b in 0..3 {
            at[b] = snap(&self.lines_mm(b), from[b]);
        }
        let l = self.lines_mm(axis);
        let (i0, i1) = (snap(&l, from[axis]), snap(&l, to));
        let sign = if i1 >= i0 { 1.0 } else { -1.0 };
        let mut terms = Vec::new();
        for i in i0.min(i1)..i0.max(i1) {
            let mut p = at;
            p[axis] = i;
            let id = self.grid.idx(p[0], p[1], p[2]);
            terms.push((axis as u8, id as u32, (sign * (l[i + 1] - l[i]) * 1e-3) as f32));
        }
        self.probes.push(ProbeRt { name: name.into(), h: false, terms });
    }
    /// Current probe: the H loop around the E edges along `axis` (current in +axis) at the edge position
    /// nearest `at[axis]`, enclosing the mesh nodes inside the transverse box `lo..hi` (mm).
    pub fn add_iprobe(&mut self, name: &str, axis: usize, at: f64, lo: [f64; 3], hi: [f64; 3]) {
        let (b, c) = ((axis + 1) % 3, (axis + 2) % 3);
        let la = self.lines_mm(axis);
        let mut ia = snap(&la, at);
        if ia + 1 >= la.len() {
            ia -= 1;
        }
        let rb = range_in(&self.lines_mm(b), lo[b], hi[b]);
        let rc = range_in(&self.lines_mm(c), lo[c], hi[c]);
        assert!(!rb.is_empty() && !rc.is_empty() && rb.start > 0 && rc.start > 0, "current probe {name} encloses no nodes");
        let (b0, b1, c0, c1) = (rb.start, rb.end - 1, rc.start, rc.end - 1);
        let idx = |pa: usize, pb: usize, pc: usize| -> u32 {
            let mut p = [0usize; 3];
            p[axis] = pa;
            p[b] = pb;
            p[c] = pc;
            self.grid.idx(p[0], p[1], p[2]) as u32
        };
        let (hb, hc) = ((3 + b) as u8, (3 + c) as u8);
        let mut terms = Vec::new();
        // counter-clockwise seen from +axis: +b along c = c0-1/2, +c along b = b1+1/2, -b along c1+1/2, -c along b0-1/2
        for pb in b0..=b1 {
            let w = self.dual(b, pb) as f32;
            terms.push((hb, idx(ia, pb, c0 - 1), w));
            terms.push((hb, idx(ia, pb, c1), -w));
        }
        for pc in c0..=c1 {
            let w = self.dual(c, pc) as f32;
            terms.push((hc, idx(ia, b1, pc), w));
            terms.push((hc, idx(ia, b0 - 1, pc), -w));
        }
        self.probes.push(ProbeRt { name: name.into(), h: true, terms });
    }
}

/// Everything a backend needs to run.
#[derive(Clone, Debug)]
pub struct Problem {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub dt: f64,
    /// 1 / primary edge length per axis (len n-1, padded to n).
    pub inv_d: [Vec<f32>; 3],
    /// 1 / dual edge length per axis (len n).
    pub inv_dd: [Vec<f32>; 3],
    pub ca: [Vec<f32>; 3],
    pub cb: [Vec<f32>; 3],
    pub pml: [PmlAxis; 3],
    pub ports: Vec<PortRt>,
    /// Field probes recorded every step (E terms at n dt, H terms at (n - 1/2) dt).
    pub probes: Vec<ProbeRt>,
    pub grid: Grid,
    /// Debye polarisation currents of wideband dielectrics (see `disp.rs`), if any.
    pub disp: Option<Disp>,
    /// Surface-impedance copper edges (see `cond.rs`), if conductor loss is modelled.
    pub cond: Option<crate::cond::Cond>,
    /// H components next to copper edges whose update is scaled (see `edge.rs`), if any.
    pub hcorr: Option<HCorr>,
}

/// Sparse H-update corrections: after the H update, `H += fac * dH` on the listed components, where dH is
/// the change the plain update made (`ids`: index | component << 30). `fac = 1/kappa - 1`.
#[derive(Clone, Debug, Default)]
pub struct HCorr {
    pub ids: Vec<u32>,
    pub fac: Vec<f32>,
}

/// Dispersive edges live in the mesh-plane range `k0..=k1`; per component the compact index is
/// q = ((k - k0) * ny + j) * nx + i, and per edge `np` poles.
#[derive(Clone, Debug)]
pub struct Disp {
    pub k0: usize,
    pub k1: usize,
    pub np: usize,
    /// per pole: k_p and (1 + k_p) / 2
    pub kp: Vec<f32>,
    pub hk: Vec<f32>,
    /// per component, compact: b_p of each edge at [q * np + p] (0 where there is no dielectric)
    pub b: [Vec<f32>; 3],
}

impl Disp {
    pub fn nd(&self, nx: usize, ny: usize) -> usize {
        nx * ny * (self.k1 - self.k0 + 1)
    }
}

impl Problem {
    pub fn n(&self) -> usize {
        self.nx * self.ny * self.nz
    }
}

#[derive(Clone, Debug)]
pub struct BuildOpts {
    /// CPML thickness in cells at x-, x+, y-, y+, z-, z+ (0 = PEC wall).
    pub pml: [usize; 6],
    /// Frequency at which `tand` is turned into a conductivity.
    pub f_loss: f64,
    /// Highest frequency of interest (sets the CPML alpha).
    pub fmax: f64,
    pub cfl: f64,
    /// Wideband (Debye) dielectric loss instead of the constant conductivity at `f_loss`.
    pub wideband: Option<crate::disp::Wideband>,
    /// Conductor loss on copper sheets (surface impedance), per sheet in `Layout::sheets` order.
    pub conductor: Option<crate::cond::Conductor>,
    /// Sub-cell copper edges: conformal and edge-singularity corrections (`edge.rs`).
    pub edge: Option<crate::edge::EdgeOpts>,
}

impl Default for BuildOpts {
    fn default() -> Self {
        BuildOpts { pml: [8; 6], f_loss: 1e9, fmax: 6e9, cfl: 0.95, wideband: None, conductor: None, edge: None }
    }
}

/// Per-layer metal flags of in-plane edges: `ex[(j)*(nx) + i]` for Ex(i,j), `ey` likewise for Ey(i,j).
pub struct SheetMask {
    pub k: usize,
    pub ex: Vec<bool>,
    pub ey: Vec<bool>,
    /// Net index (into the list returned with the masks; 0 = no net) of the copper that set each edge.
    pub exn: Vec<u32>,
    pub eyn: Vec<u32>,
}

fn paint(mask: &mut [bool], nets: &mut [u32], net: u32, nx: usize, xs: &[f64], ys: &[f64], p: &Prim) {
    // xs: sample x positions (len nx), ys: sample rows; mask[(j * nx) + i]
    debug_assert_eq!(nx, xs.len());
    paint_with(xs, ys, p, |q| {
        mask[q] = p.dark;
        nets[q] = if p.dark { net } else { 0 };
    });
}

/// Visit the samples `(xs[i], ys[j])` (index `j * xs.len() + i`, both sorted) that lie in `p`'s shape:
/// even-odd scanline spans for polygons, the inside test otherwise.
pub(crate) fn paint_with(xs: &[f64], ys: &[f64], p: &Prim, mut set: impl FnMut(usize)) {
    let nx = xs.len();
    let bb = p.shape.bbox();
    let ir = range_in(xs, bb.x0, bb.x1);
    let jr = range_in(ys, bb.y0, bb.y1);
    if ir.is_empty() || jr.is_empty() {
        return;
    }
    match &p.shape {
        Shape::Polygon(pts) => {
            // scanline: crossings of each row with the polygon edges, even-odd spans
            let mut xc: Vec<f64> = Vec::new();
            for j in jr {
                let y = ys[j];
                xc.clear();
                let n = pts.len();
                for e in 0..n {
                    let (a, b) = (pts[e], pts[(e + 1) % n]);
                    if (a.y > y) != (b.y > y) {
                        xc.push(a.x + (y - a.y) * (b.x - a.x) / (b.y - a.y));
                    }
                }
                xc.sort_by(|a, b| a.partial_cmp(b).unwrap());
                for s in xc.chunks_exact(2) {
                    for i in range_in(xs, s[0], s[1]) {
                        set(j * nx + i);
                    }
                }
            }
        }
        sh => {
            for j in jr {
                for i in ir.clone() {
                    if sh.contains(P2::new(xs[i], ys[j])) {
                        set(j * nx + i);
                    }
                }
            }
        }
    }
}

/// Rasterise copper sheets onto the grid's in-plane edges (coordinates in mm on the layout side). Returns
/// the masks and the net names they index (entry 0 is "" for copper without a net).
pub fn sheet_masks(g: &Grid, lay: &Layout) -> (Vec<SheetMask>, Vec<String>) {
    let mut names: Vec<String> = vec![String::new()];
    let mut index = std::collections::HashMap::new();
    let sheet_nets: Vec<Vec<u32>> = lay
        .sheets
        .iter()
        .map(|s| {
            s.prims
                .iter()
                .map(|p| match &p.net {
                    Some(n) => *index.entry(n.clone()).or_insert_with(|| {
                        names.push(n.clone());
                        (names.len() - 1) as u32
                    }),
                    None => 0,
                })
                .collect()
        })
        .collect();
    use rayon::prelude::*;
    let (nx, ny) = (g.nx(), g.ny());
    let xmm: Vec<f64> = g.x.iter().map(|v| v * 1e3).collect();
    let ymm: Vec<f64> = g.y.iter().map(|v| v * 1e3).collect();
    let zmm: Vec<f64> = g.z.iter().map(|v| v * 1e3).collect();
    // midpoints: Ex samples at (x_{i+1/2}, y_j), Ey at (x_i, y_{j+1/2}); the last entry is a far sentinel
    let mut xh: Vec<f64> = xmm.windows(2).map(|w| 0.5 * (w[0] + w[1])).collect();
    xh.push(f64::INFINITY);
    let mut yh: Vec<f64> = ymm.windows(2).map(|w| 0.5 * (w[0] + w[1])).collect();
    yh.push(f64::INFINITY);
    let masks: Vec<SheetMask> = lay
        .sheets
        .par_iter()
        .zip(sheet_nets.par_iter())
        .map(|(s, sn)| {
            let k = snap(&zmm, s.z);
            let mut ex = vec![false; nx * ny];
            let mut ey = vec![false; nx * ny];
            let mut exn = vec![0u32; nx * ny];
            let mut eyn = vec![0u32; nx * ny];
            for (p, &net) in s.prims.iter().zip(sn) {
                if lay.copper_window.is_some_and(|w| !w.overlaps(&p.shape.bbox())) {
                    continue;
                }
                paint(&mut ex, &mut exn, net, nx, &xh, &ymm, p);
                paint(&mut ey, &mut eyn, net, nx, &xmm, &yh, p);
            }
            if let Some(w) = lay.copper_window {
                // Inside the window the copper is as drawn. Outside it (graded cells and the absorber) copper
                // is not rasterised from the layout, where it would merge across gaps on coarse cells;
                // instead whatever crosses the window's edge continues straight out to the grid boundary,
                // so a cropped plane or track runs into the absorbing layer instead of ending in the domain.
                let span = |v: &[f64], lo: f64, hi: f64| -> Option<(usize, usize)> {
                    let a = v.iter().position(|&x| x >= lo && x <= hi)?;
                    let b = v.iter().rposition(|&x| x >= lo && x <= hi)?;
                    Some((a, b))
                };
                let (exi, exj) = (span(&xh, w.x0, w.x1), span(&ymm, w.y0, w.y1));
                let (eyi, eyj) = (span(&xmm, w.x0, w.x1), span(&yh, w.y0, w.y1));
                let extrude = |m: &mut Vec<bool>, mn: &mut Vec<u32>, si: Option<(usize, usize)>, sj: Option<(usize, usize)>| {
                    let (Some((i0, i1)), Some((j0, j1))) = (si, sj) else {
                        m.iter_mut().for_each(|v| *v = false);
                        return;
                    };
                    let src: Vec<(bool, u32)> = (0..ny * nx)
                        .map(|q| {
                            let (i, j) = (q % nx, q / nx);
                            let s = j.clamp(j0, j1) * nx + i.clamp(i0, i1);
                            (m[s], mn[s])
                        })
                        .collect();
                    for (q, (on, net)) in src.into_iter().enumerate() {
                        m[q] = on;
                        mn[q] = net;
                    }
                };
                extrude(&mut ex, &mut exn, exi, exj);
                extrude(&mut ey, &mut eyn, eyi, eyj);
            }
            SheetMask { k, ex, ey, exn, eyn }
        })
        .collect();
    let mut masks = masks;
    // copper added for lumped elements ("~name") joins nets on purpose: it never conflicts
    let wild: Vec<bool> = names.iter().map(|n| n.starts_with('~')).collect();
    resolve_conflicts(g, lay, &mut masks, &wild);
    (masks, names)
}

/// Keep the mesh's connectivity equal to the layout's where gaps are close to a cell:
/// - at a node where in-plane copper edges of two nets meet, the edges of the net with fewer edges at that
///   node are removed (the other net's copper is the one the node lies in);
/// - at a via barrel's nodes, copper edges of nets other than the barrel's are removed (a too-small antipad
///   in the raster). The barrel's net is the most common net among copper edges at its nodes.
/// Copper without a net never conflicts. Returns the number of edges removed.
pub fn resolve_conflicts(g: &Grid, lay: &Layout, masks: &mut [SheetMask], wild: &[bool]) -> usize {
    let (nx, ny) = (g.nx(), g.ny());
    // element copper and net-less copper are both net 0 here
    let real = |n: u32| if wild[n as usize] { 0 } else { n };
    let mut removed = 0usize;
    // incident in-plane edges of node (i, j): (is_x, q)
    let incident = |i: usize, j: usize| -> [(bool, Option<usize>); 4] {
        [
            (true, (i + 1 < nx).then_some(j * nx + i)),
            (true, (i > 0).then(|| j * nx + i - 1)),
            (false, (j + 1 < ny).then_some(j * nx + i)),
            (false, (j > 0).then(|| (j - 1) * nx + i)),
        ]
    };
    // barrels: their net, from the copper at their nodes on every sheet they cross
    let mut barrel_net: Vec<(Vec<(usize, usize)>, usize, usize, u32)> = Vec::new();
    for b in &lay.barrels {
        let (cells, k0, k1) = barrel_nodes(g, b);
        let mut count: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
        for m in masks.iter().filter(|m| m.k >= k0 && m.k <= k1) {
            for &(i, j) in &cells {
                for (isx, q) in incident(i, j) {
                    if let Some(q) = q {
                        let (on, net) = if isx { (m.ex[q], real(m.exn[q])) } else { (m.ey[q], real(m.eyn[q])) };
                        if on && net != 0 {
                            *count.entry(net).or_default() += 1;
                        }
                    }
                }
            }
        }
        let net = count.into_iter().max_by_key(|&(n, c)| (c, std::cmp::Reverse(n))).map(|x| x.0).unwrap_or(0);
        barrel_net.push((cells, k0, k1, net));
    }
    for m in masks.iter_mut() {
        // barrels first: they are drilled copper of their own net
        for (cells, k0, k1, bnet) in &barrel_net {
            if *bnet == 0 || m.k < *k0 || m.k > *k1 {
                continue;
            }
            for &(i, j) in cells {
                for (isx, q) in incident(i, j) {
                    if let Some(q) = q {
                        let (on, net) = if isx { (&mut m.ex[q], real(m.exn[q])) } else { (&mut m.ey[q], real(m.eyn[q])) };
                        if *on && net != 0 && net != *bnet {
                            *on = false;
                            removed += 1;
                        }
                    }
                }
            }
        }
        // copper of two nets at one node
        for j in 0..ny {
            for i in 0..nx {
                let inc = incident(i, j);
                let mut nets: [(u32, usize); 4] = [(0, 0); 4];
                let mut nn = 0;
                for (isx, q) in inc {
                    if let Some(q) = q {
                        let (on, net) = if isx { (m.ex[q], real(m.exn[q])) } else { (m.ey[q], real(m.eyn[q])) };
                        if on && net != 0 {
                            if let Some(e) = nets[..nn].iter_mut().find(|e| e.0 == net) {
                                e.1 += 1;
                            } else {
                                nets[nn] = (net, 1);
                                nn += 1;
                            }
                        }
                    }
                }
                if nn < 2 {
                    continue;
                }
                let keep = nets[..nn].iter().max_by_key(|e| (e.1, std::cmp::Reverse(e.0))).unwrap().0;
                for (isx, q) in inc {
                    if let Some(q) = q {
                        let (on, net) = if isx { (&mut m.ex[q], real(m.exn[q])) } else { (&mut m.ey[q], real(m.eyn[q])) };
                        if *on && net != 0 && net != keep {
                            *on = false;
                            removed += 1;
                        }
                    }
                }
            }
        }
    }
    removed
}

/// Mesh nodes (i, j) and the z-edge range k0..k1 a barrel occupies: the nodes within its radius, or the
/// nearest node when it is thinner than the mesh.
pub fn barrel_nodes(g: &Grid, b: &Barrel) -> (Vec<(usize, usize)>, usize, usize) {
    let xmm: Vec<f64> = g.x.iter().map(|v| v * 1e3).collect();
    let ymm: Vec<f64> = g.y.iter().map(|v| v * 1e3).collect();
    let zmm: Vec<f64> = g.z.iter().map(|v| v * 1e3).collect();
    let (k0, k1) = (snap(&zmm, b.z0.min(b.z1)), snap(&zmm, b.z0.max(b.z1)));
    let r = b.r;
    let mut cells = Vec::new();
    for j in range_in(&ymm, b.c.y - r, b.c.y + r) {
        for i in range_in(&xmm, b.c.x - r, b.c.x + r) {
            if (xmm[i] - b.c.x).powi(2) + (ymm[j] - b.c.y).powi(2) <= r * r {
                cells.push((i, j));
            }
        }
    }
    if cells.is_empty() {
        cells.push((snap(&xmm, b.c.x), snap(&ymm, b.c.y)));
    }
    (cells, k0, k1)
}

fn cpml_axis(lines: &[f64], nlo: usize, nhi: usize, dt: f64, fmax: f64) -> PmlAxis {
    let n = lines.len();
    let mut ax = PmlAxis { map: vec![NONE; n], ..Default::default() };
    if nlo == 0 && nhi == 0 {
        return ax;
    }
    let m = 3.0f64;
    let alpha_max = 2.0 * std::f64::consts::PI * EPS0 * (fmax / 20.0).max(1e6);
    // depth (0 at the inner PML interface, 1 at the outer wall) of a position on either side
    let mut idx: Vec<usize> = Vec::new();
    if nlo > 0 {
        idx.extend(0..(nlo + 1).min(n));
    }
    if nhi > 0 {
        let start = n.saturating_sub(nhi + 1);
        for i in start..n {
            if !idx.contains(&i) {
                idx.push(i);
            }
        }
    }
    let depth = |x: f64| -> (f64, f64) {
        // returns (rho, local cell size)
        if nlo > 0 && x < lines[nlo] {
            let w = lines[nlo] - lines[0];
            ((lines[nlo] - x) / w, w / nlo as f64)
        } else if nhi > 0 && x > lines[n - 1 - nhi] {
            let w = lines[n - 1] - lines[n - 1 - nhi];
            ((x - lines[n - 1 - nhi]) / w, w / nhi as f64)
        } else {
            (0.0, 1.0)
        }
    };
    let coef = |x: f64| -> (f32, f32) {
        let (rho, h) = depth(x);
        if rho <= 0.0 {
            return (1.0, 0.0);
        }
        let rho = rho.min(1.0);
        let smax = 0.8 * (m + 1.0) / (ETA0 * h);
        let sigma = smax * rho.powf(m);
        let alpha = alpha_max * (1.0 - rho);
        let b = (-(sigma + alpha) * dt / EPS0).exp();
        let c = if sigma > 0.0 { sigma / (sigma + alpha) * (b - 1.0) } else { 0.0 };
        (b as f32, c as f32)
    };
    for (s, &i) in idx.iter().enumerate() {
        ax.map[i] = s as u32;
        ax.slab.push(i as u32);
        let (be, ce) = coef(lines[i]);
        let xh = if i + 1 < n { 0.5 * (lines[i] + lines[i + 1]) } else { lines[i] };
        let (bh, ch) = if i + 1 < n { coef(xh) } else { (1.0, 0.0) };
        ax.be.push(be);
        ax.ce.push(ce);
        ax.bh.push(bh);
        ax.ch.push(ch);
    }
    ax
}

/// Cell materials: whether cell (i, j) lies inside the board outline (by its centre; index j * nx + i), and
/// the slab of each z cell k (by its centre), if any.
pub fn cell_materials(g: &Grid, lay: &Layout) -> (Vec<bool>, Vec<Option<usize>>) {
    let (nx, ny, nz) = (g.nx(), g.ny(), g.nz());
    let mut inside = vec![true; nx * ny];
    if !lay.outline.is_empty() {
        let mut bb = BBox::empty();
        lay.outline.iter().flatten().for_each(|&p| bb.add(p));
        for j in 0..ny.saturating_sub(1) {
            for i in 0..nx.saturating_sub(1) {
                let p = P2::new(0.5e3 * (g.x[i] + g.x[i + 1]), 0.5e3 * (g.y[j] + g.y[j + 1]));
                let mut c = false;
                if bb.contains(p) {
                    for l in &lay.outline {
                        if crate::geom::point_in_polygon(l, p) {
                            c = !c;
                        }
                    }
                }
                inside[j * nx + i] = c;
            }
        }
    }
    let slab_of_k = (0..nz.saturating_sub(1))
        .map(|k| {
            let zc = 0.5e3 * (g.z[k] + g.z[k + 1]);
            lay.slabs.iter().position(|s| zc >= s.z0.min(s.z1) && zc <= s.z0.max(s.z1))
        })
        .collect();
    (inside, slab_of_k)
}

/// Build the update coefficients and runtime ports.
pub fn build(g: &Grid, lay: &Layout, lumped: &[Lumped], o: &BuildOpts) -> Problem {
    let (nx, ny, nz) = (g.nx(), g.ny(), g.nz());
    let n = g.n();
    let dmin = [g.min_step(0), g.min_step(1), g.min_step(2)];
    let dt = o.cfl / (C0 * (1.0 / (dmin[0] * dmin[0]) + 1.0 / (dmin[1] * dmin[1]) + 1.0 / (dmin[2] * dmin[2])).sqrt());

    let d = |a: usize, i: usize| -> f64 {
        let l = g.axis(a);
        if i + 1 < l.len() {
            l[i + 1] - l[i]
        } else {
            0.0
        }
    };
    // dual length at node i: half of each adjacent cell
    let dd = |a: usize, i: usize| -> f64 {
        let l = g.axis(a);
        let lo = if i > 0 { l[i] - l[i - 1] } else { 0.0 };
        let hi = if i + 1 < l.len() { l[i + 1] - l[i] } else { 0.0 };
        0.5 * (lo + hi)
    };
    let mut inv_d: [Vec<f32>; 3] = Default::default();
    let mut inv_dd: [Vec<f32>; 3] = Default::default();
    for a in 0..3 {
        let len = g.axis(a).len();
        inv_d[a] = (0..len).map(|i| if d(a, i) > 0.0 { (1.0 / d(a, i)) as f32 } else { 0.0 }).collect();
        inv_dd[a] = (0..len).map(|i| if dd(a, i) > 0.0 { (1.0 / dd(a, i)) as f32 } else { 0.0 }).collect();
    }

    // cell materials: inside-outline mask per (i,j) cell, slab per k cell
    let (inside, slab_of_k) = cell_materials(g, lay);
    let w_loss = 2.0 * std::f64::consts::PI * o.f_loss;
    // wideband loss: per slab (eps_inf, deps per pole); the copper-loss gap is the slab's thickness
    let poles: Vec<f64> = o.wideband.as_ref().map(|w| w.poles()).unwrap_or_default();
    let npole = poles.len();
    let fits: Vec<Option<(f64, Vec<f64>)>> = lay
        .slabs
        .iter()
        .map(|s| {
            let w = o.wideband.as_ref()?;
            if s.tand <= 0.0 && !w.copper {
                return None;
            }
            let (einf, deps, _) = crate::disp::fit(s.eps, s.tand, (s.z1 - s.z0).abs(), w);
            Some((einf, deps))
        })
        .collect();
    // (eps_r or eps_inf, sigma, slab index for the poles) of cell (i,j,k); out-of-range cells are air
    let body_of = crate::bodies::body_cells(g, &lay.bodies);
    for (bi, b) in lay.bodies.iter().enumerate().filter(|(_, b)| !b.metal) {
        let n = body_of.iter().filter(|&&v| v as usize == bi + 1).count();
        eprintln!("fieldsolve: body {}: eps {} tan d {} in {n} cells", b.name, b.eps, b.tand);
        if n == 0 {
            eprintln!("fieldsolve: WARNING: body {} fills no cells (outside the domain, or thinner than a cell)", b.name);
        }
    }
    let cell = |i: isize, j: isize, k: isize| -> (f64, f64, Option<usize>) {
        if i < 0 || j < 0 || k < 0 || i >= nx as isize - 1 || j >= ny as isize - 1 || k >= nz as isize - 1 {
            return (1.0, 0.0, None);
        }
        if !lay.bodies.is_empty() {
            let b = body_of[(k as usize * (ny - 1) + j as usize) * (nx - 1) + i as usize];
            if b > 0 {
                let b = &lay.bodies[b as usize - 1];
                return (b.eps, w_loss * EPS0 * b.eps * b.tand, None);
            }
        }
        match slab_of_k[k as usize] {
            Some(si) if inside[j as usize * nx + i as usize] => {
                let s = &lay.slabs[si];
                match &fits[si] {
                    Some((einf, _)) => (*einf, 0.0, Some(si)),
                    None => (s.eps, w_loss * EPS0 * s.eps * s.tand, None),
                }
            }
            _ => (1.0, 0.0, None),
        }
    };
    // plane range of edges that touch a dispersive cell
    let (kd0, kd1) = {
        let ks: Vec<usize> = (0..nz.saturating_sub(1)).filter(|&k| slab_of_k[k].is_some_and(|si| fits[si].is_some())).collect();
        match (ks.first(), ks.last()) {
            (Some(&a), Some(&b)) => (a, (b + 1).min(nz - 1)),
            _ => (1, 0),
        }
    };
    let nd = if kd1 >= kd0 { nx * ny * (kd1 - kd0 + 1) } else { 0 };
    let mut dep: [Vec<f32>; 3] = std::array::from_fn(|_| vec![0.0f32; nd * npole]);
    let cl = |a: usize, i: isize| -> f64 { if i < 0 { 0.0 } else { d(a, i as usize) } };

    // per-edge effective eps and sigma: area-weighted mean of the four cells around the edge
    let mut eps: [Vec<f32>; 3] = [vec![1.0; n], vec![1.0; n], vec![1.0; n]];
    let mut sig: [Vec<f32>; 3] = [vec![0.0; n], vec![0.0; n], vec![0.0; n]];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let id = g.idx(i, j, k);
                let (ii, jj, kk) = (i as isize, j as isize, k as isize);
                for a in 0..3 {
                    // the two transverse axes and the cell index along this edge
                    let (b, c) = ((a + 1) % 3, (a + 2) % 3);
                    let pos = [ii, jj, kk];
                    let mut ws = 0.0;
                    let (mut e, mut s) = (0.0, 0.0);
                    let mut dp = [0.0f64; 8];
                    let mut any = false;
                    for db in [-1isize, 0] {
                        for dc in [-1isize, 0] {
                            let mut q = pos;
                            q[b] += db;
                            q[c] += dc;
                            let w = cl(b, q[b]) * cl(c, q[c]);
                            if w <= 0.0 {
                                continue;
                            }
                            let (ce, cs, si) = cell(q[0], q[1], q[2]);
                            ws += w;
                            e += w * ce;
                            s += w * cs;
                            if let Some(si) = si {
                                let deps = &fits[si].as_ref().unwrap().1;
                                for p in 0..npole {
                                    dp[p] += w * deps[p];
                                }
                                any = true;
                            }
                        }
                    }
                    if ws > 0.0 {
                        eps[a][id] = (e / ws) as f32;
                        sig[a][id] = (s / ws) as f32;
                        if any && k >= kd0 && k <= kd1 {
                            let q = ((k - kd0) * ny + j) * nx + i;
                            for p in 0..npole {
                                dep[a][q * npole + p] = (dp[p] / ws) as f32;
                            }
                        }
                    }
                }
            }
        }
    }

    let (masks, mask_nets) = sheet_masks(g, lay);
    // sub-cell copper edges: redistribute the permittivity of edges next to copper (before the lumped
    // elements add theirs)
    let efix = o.edge.as_ref().map(|eo| crate::edge::analyse(g, lay, &masks, lumped, &o.pml, eo));
    if let (Some(fx), Some(eo)) = (&efix, &o.edge) {
        for &(a, id, beta) in &fx.e {
            let (a, id) = (a as usize, id as usize);
            let (i, j, k) = (id % nx, (id / nx) % ny, id / (nx * ny));
            let (b, c) = ((a + 1) % 3, (a + 2) % 3);
            let pos = [i as isize, j as isize, k as isize];
            let mut enew = 0.0;
            let mut t = 0;
            for db in [-1isize, 0] {
                for dc in [-1isize, 0] {
                    let mut q = pos;
                    q[b] += db;
                    q[c] += dc;
                    enew += cell(q[0], q[1], q[2]).0 * beta[t] as f64;
                    t += 1;
                }
            }
            let old = eps[a][id] as f64;
            if old <= 0.0 || enew <= 0.0 {
                continue;
            }
            // quantised (log steps) so that the Metal material table stays small
            let f = ((enew / old).ln() / eo.quant).round() * eo.quant;
            let f = f.exp();
            eps[a][id] = (old * f) as f32;
            sig[a][id] = (sig[a][id] as f64 * f) as f32;
            if npole > 0 && k >= kd0 && k <= kd1 {
                let qd = ((k - kd0) * ny + j) * nx + i;
                for p in 0..npole {
                    dep[a][qd * npole + p] = (dep[a][qd * npole + p] as f64 * f) as f32;
                }
            }
        }
    }

    // lumped elements: extra conductivity / capacitance per edge, and source/probe lists for ports
    let mut ports = Vec::new();
    let mut shorts: Vec<(usize, usize)> = Vec::new();
    let mut loaded: Vec<(usize, usize)> = Vec::new();
    let xs = [&g.x, &g.y, &g.z];
    for l in lumped {
        let a = l.axis;
        let mut rng: [std::ops::Range<usize>; 3] = Default::default();
        for b in 0..3 {
            let (lo, hi) = (l.lo[b] * 1e-3, l.hi[b] * 1e-3);
            if b == a {
                let (i0, i1) = (snap(xs[b], lo), snap(xs[b], hi));
                assert!(i1 > i0, "lumped element {} has no length along its axis", l.name);
                rng[b] = i0..i1; // edges i0..i1-1
            } else {
                let r = range_in(xs[b], lo - 1e-12, hi + 1e-12);
                rng[b] = if r.is_empty() { let s = snap(xs[b], 0.5 * (lo + hi)); s..s + 1 } else { r };
            }
        }
        let k_ser = rng[a].len() as f64;
        let m_par = (rng[(a + 1) % 3].len() * rng[(a + 2) % 3].len()) as f64;
        let sign: f64 = if l.plus_at_hi { -1.0 } else { 1.0 };
        let mut rt = PortRt { name: l.name.clone(), r: l.r, src: vec![], probe: vec![] };
        let mut edges = Vec::new();
        for k in rng[2].clone() {
            for j in rng[1].clone() {
                for i in rng[0].clone() {
                    let id = g.idx(i, j, k);
                    let p = [i, j, k];
                    let len = d(a, p[a]);
                    let area = dd((a + 1) % 3, p[(a + 1) % 3]) * dd((a + 2) % 3, p[(a + 2) % 3]);
                    edges.push((id, len, area));
                }
            }
        }
        if l.r == 0.0 && l.c == 0.0 && !l.port {
            shorts.extend(edges.iter().map(|e| (a, e.0)));
            continue;
        }
        loaded.extend(edges.iter().map(|e| (a, e.0)));
        for &(id, len, area) in &edges {
            if l.r > 0.0 {
                let re = l.r * m_par / k_ser;
                sig[a][id] += (len / (re * area)) as f32;
            }
            if l.c > 0.0 {
                let ce = l.c * k_ser / m_par;
                eps[a][id] += (ce * len / (area * EPS0)) as f32;
            }
        }
        if l.port {
            for &(id, len, area) in &edges {
                let re = l.r * m_par / k_ser;
                // gain multiplies cb later (cb is computed below); keep the raw term here
                rt.src.push((a as u8, id as u32, (sign / (k_ser * re * area)) as f32));
                rt.probe.push((a as u8, id as u32, (sign * len / m_par) as f32));
            }
            ports.push(rt);
        }
    }

    // E update coefficients
    let mut ca: [Vec<f32>; 3] = [vec![0.0; n], vec![0.0; n], vec![0.0; n]];
    let mut cb: [Vec<f32>; 3] = [vec![0.0; n], vec![0.0; n], vec![0.0; n]];
    let pc = crate::disp::pole_consts(&poles, dt);
    // sum of b_p of an edge (Debye poles add B/2 to both sides of the E update)
    let bsum = |a: usize, id: usize| -> f64 {
        if npole == 0 || nd == 0 {
            return 0.0;
        }
        let k = id / (nx * ny);
        if k < kd0 || k > kd1 {
            return 0.0;
        }
        let q = id - kd0 * nx * ny;
        (0..npole).map(|p| pc[p].2 * dep[a][q * npole + p] as f64).sum()
    };
    let coef = |e: f64, s: f64, bb: f64| -> (f32, f32) {
        let d = e / dt + 0.5 * s + 0.5 * bb;
        (((e / dt - 0.5 * s + 0.5 * bb) / d) as f32, (1.0 / d) as f32)
    };
    for a in 0..3 {
        for id in 0..n {
            let (x, y) = coef(eps[a][id] as f64 * EPS0, sig[a][id] as f64, bsum(a, id));
            ca[a][id] = x;
            cb[a][id] = y;
        }
    }
    // edges that do not exist or lie on the outer wall (tangential E = 0) and metal edges
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let id = g.idx(i, j, k);
                let p = [i, j, k];
                let nn = [nx, ny, nz];
                for a in 0..3 {
                    let mut dead = p[a] + 1 >= nn[a];
                    for b in 0..3 {
                        if b != a && (p[b] == 0 || p[b] + 1 >= nn[b]) {
                            dead = true;
                        }
                    }
                    if dead {
                        ca[a][id] = 0.0;
                        cb[a][id] = 0.0;
                    }
                }
            }
        }
    }
    for m in &masks {
        for j in 0..ny {
            for i in 0..nx {
                let id = g.idx(i, j, m.k);
                if m.ex[j * nx + i] {
                    ca[0][id] = 0.0;
                    cb[0][id] = 0.0;
                }
                if m.ey[j * nx + i] {
                    ca[1][id] = 0.0;
                    cb[1][id] = 0.0;
                }
            }
        }
    }
    for b in &lay.barrels {
        let (cells, k0, k1) = barrel_nodes(g, b);
        for (i, j) in cells {
            for k in k0..k1 {
                let id = g.idx(i, j, k);
                ca[2][id] = 0.0;
                cb[2][id] = 0.0;
            }
        }
    }
    // metal boxes (sheets in any axis plane, bodies.rs)
    for (a, id) in crate::bodies::metal_edges(g, &lay.bodies) {
        ca[a][id] = 0.0;
        cb[a][id] = 0.0;
    }
    // conformal edges: metal by the midpoint rule, but mostly open
    if let Some(fx) = &efix {
        for &(a, id) in &fx.open {
            let (a, id) = (a as usize, id as usize);
            let (x, y) = coef(eps[a][id] as f64 * EPS0, sig[a][id] as f64, bsum(a, id));
            ca[a][id] = x;
            cb[a][id] = y;
        }
        for &(a, id) in &fx.close {
            ca[a as usize][id as usize] = 0.0;
            cb[a as usize][id as usize] = 0.0;
        }
    }
    // lumped elements win over copper they overlap (a port drawn across a pad gap that the mesh closed)
    for &(a, id) in &loaded {
        let (x, y) = coef(eps[a][id] as f64 * EPS0, sig[a][id] as f64, bsum(a, id));
        ca[a][id] = x;
        cb[a][id] = y;
    }
    for &(a, id) in &shorts {
        ca[a][id] = 0.0;
        cb[a][id] = 0.0;
    }
    for p in &mut ports {
        for s in &mut p.src {
            s.2 *= cb[s.0 as usize][s.1 as usize];
        }
    }
    let cond = o.conductor.as_ref().map(|c| cond_edges(g, &masks, &mask_nets, &cb, &o.pml, c, dt));

    let pml = [
        cpml_axis(&g.x, o.pml[0], o.pml[1], dt, o.fmax),
        cpml_axis(&g.y, o.pml[2], o.pml[3], dt, o.fmax),
        cpml_axis(&g.z, o.pml[4], o.pml[5], dt, o.fmax),
    ];
    let disp = (npole > 0 && nd > 0).then(|| Disp {
        k0: kd0,
        k1: kd1,
        np: npole,
        kp: pc.iter().map(|c| c.0 as f32).collect(),
        hk: pc.iter().map(|c| c.1 as f32).collect(),
        // b_p per edge; metal edges (cb = 0) keep b but are skipped by the updates
        b: std::array::from_fn(|a| (0..nd * npole).map(|qp| pc[qp % npole].2 as f32 * dep[a][qp]).collect()),
    });
    let hcorr = efix.as_ref().filter(|f| !f.h.is_empty()).map(|f| HCorr {
        ids: f.h.iter().map(|&(a, id, _)| id | (a as u32) << 30).collect(),
        fac: f.h.iter().map(|&(_, _, k)| 1.0 / k - 1.0).collect(),
    });
    let p = Problem { nx, ny, nz, dt, inv_d, inv_dd, ca, cb, pml, ports, probes: vec![], grid: g.clone(), disp, cond, hcorr };
    if let (Some(fx), Some(eo)) = (&efix, &o.edge) {
        let st = &fx.stats;
        eprintln!(
            "fieldsolve: edges ({}{}): {} E conformal, {} E singular, {} H conformal, {} H singular; {} metal edges opened, {} closed, {} staircase nodes kept for connectivity, {} fits rejected, {} shared dropped; only where in-plane cells <= {:.3} mm",
            if eo.conformal { "conformal" } else { "" },
            if eo.singular { if eo.conformal { "+singular" } else { "singular" } } else { "" },
            st.e_conformal, st.e_singular, st.h_conformal, st.h_singular, st.opened, st.closed, st.kept_for_connectivity, st.fit_rejected, st.dropped_shared, st.max_cell
        );
        if eo.dt_check && (!fx.e.is_empty() || !fx.h.is_empty()) {
            let t0 = std::time::Instant::now();
            let mu = crate::edge::stability(&p, eo.dt_iters);
            let lim = 4.0 * o.cfl * o.cfl;
            eprintln!("fieldsolve: edges: dt^2 lambda_max {:.4} (limit {:.4} at cfl {}), {:.1} s", mu, lim, o.cfl, t0.elapsed().as_secs_f64());
            if mu > lim {
                let cfl = o.cfl * (lim / mu).sqrt() * 0.99;
                eprintln!("fieldsolve: edges: time step lowered by {:.3}", cfl / o.cfl);
                let mut o2 = o.clone();
                o2.cfl = cfl;
                o2.edge.as_mut().unwrap().dt_check = false;
                return build(g, lay, lumped, &o2);
            }
        }
    }
    p
}

/// Copper sheet edges that get a surface impedance: metal in-plane edges (cb = 0) of real copper outside
/// the CPML slabs. Element copper ("~name": bars of ports and parts, whose loss is in the part model) and
/// edges inside the absorber (where the curl carries CPML terms) stay perfect conductors.
fn cond_edges(g: &Grid, masks: &[SheetMask], nets: &[String], cb: &[Vec<f32>; 3], pml: &[usize; 6], c: &crate::cond::Conductor, dt: f64) -> crate::cond::Cond {
    use crate::cond::*;
    let (nx, ny, nz) = (g.nx(), g.ny(), g.nz());
    let fk = corners(c);
    let mut tab = Vec::new();
    for (si, m) in c.metals.iter().enumerate() {
        let f = fit(m, &fk, c.f_lo, c.f_hi);
        eprintln!(
            "fieldsolve: conductor on sheet {si}: sigma {:.2e} S/m, t {:.0} um, roughness {:.1} um; R_s {:.1} mohm at 2.44 GHz; fit {} sections, worst {:.1} %",
            m.sigma, m.t * 1e6, m.rough * 1e6, zs(m, 2.44e9).re * 1e3, fk.len(), 100.0 * f.worst
        );
        tab.push(f.r0 as f32);
        tab.extend(f.rk.iter().map(|&r| r as f32));
    }
    let dec: Vec<f32> = fk.iter().map(|&f| (-dt * 2.0 * std::f64::consts::PI * f).exp() as f32).collect();
    let wild: Vec<bool> = nets.iter().map(|n| n.starts_with('~')).collect();
    let inside = |i: usize, j: usize, k: usize| i >= pml[0] && i + pml[1] < nx && j >= pml[2] && j + pml[3] < ny && k >= pml[4] && k + pml[5] < nz;
    let (mut ids, mut fac, mut layer) = (Vec::new(), Vec::new(), Vec::new());
    for (si, m) in masks.iter().enumerate() {
        if si >= c.metals.len() {
            break;
        }
        let k = m.k;
        // dual z length at plane k
        let dz = if k == 0 { g.z[1] - g.z[0] } else if k + 1 >= nz { g.z[k] - g.z[k - 1] } else { 0.5 * (g.z[k + 1] - g.z[k - 1]) };
        for j in 0..ny {
            for i in 0..nx {
                let q = j * nx + i;
                let id = g.idx(i, j, k);
                for (a, on, net) in [(0usize, m.ex[q], m.exn[q]), (1usize, m.ey[q], m.eyn[q])] {
                    if on && !wild[net as usize] && cb[a][id] == 0.0 && inside(i, j, k) && i > 0 && j > 0 && k > 0 {
                        ids.push(id as u32 | (a as u32) << 30);
                        fac.push(dz as f32);
                        layer.push(si as u32);
                    }
                }
            }
        }
    }
    eprintln!("fieldsolve: conductor loss on {} copper edges", ids.len());
    Cond { ids, fac, layer, nb: fk.len(), dec, tab }
}
