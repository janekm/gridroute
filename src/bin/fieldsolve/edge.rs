//! Sub-cell copper edges on zero-thickness sheets: (a) edge-singularity and (b) conformal corrections.
//!
//! Both are written as corrections of the diagonal material factors of the Yee scheme in integral (FIT)
//! form. With e = ∫E·dl on a primary edge, b = ∫B·dA on a primary face, h = ∫H·dl on a dual edge and
//! d = ∫D·dA on a dual face, Faraday (db/dt = -Σ e) and Ampère (dd/dt = Σ h) are exact; the
//! discretisation error sits only in the constitutive factors b = M_μ h and d = M_ε e. The plain Yee scheme
//! takes the field constant over each edge and face: M_μ = μ A / L', M_ε = ε A' / L. Here, next to a
//! copper edge, they come from an assumed local field instead:
//!
//! - **(b) conformal** (Dey–Mittra): the field is uniform on the open part and zero on the metal. An
//!   in-plane E edge that is partly covered carries e = E_open l_open, so in the Yee variable E = e/L its
//!   permittivity becomes ε L / l_open; a face of the sheet plane that is partly covered carries
//!   b = μ H A_open, so its H update is divided by A_open / A. This moves the electrical copper edge from
//!   the mesh line to where the copper is.
//! - **(a) singular**: the field is the static field of a conducting half-plane at the true (conformal) or
//!   staircased edge, φ + iψ = √(u + iw) with u the in-plane distance from the edge (metal at u < 0) and
//!   w the height above the sheet. E = -∇φ and H = -∇ψ (the TEM magnetic field is E turned by 90° in the
//!   cross-section); both vary as r^-1/2. For each E edge and H face near the copper edge the factor is
//!   κ = (field flux through the face) · L / (A · line integral), the ratio of the singular field's face
//!   average to its edge average, which is 1 for a uniform field. Line integrals are potential
//!   differences (exact); fluxes use the potentials along one face direction and Gauss quadrature along the
//!   other. This is the static-field incorporation of Shorthouse and Railton (IEEE MTT 40(5), 1992) and,
//!   for a strip ending on a mesh line, the effect the openEMS "thirds rule" approximates by meshing.
//!
//! The E factors go into the edge's permittivity (and conductivity and Debye weights), so the backends see
//! only a material change. The H factors are a sparse list (`model::HCorr`) applied after the H update.
//! Smaller open fractions make the scheme stiffer: fractions are floored at `min_open`, singular factors
//! clamped to `[kmin, kmax]`, and `dt_check` measures the largest eigenvalue of the corrected operator by
//! power iteration and lowers the time step if needed.
//!
//! The conformal rule never adds metal. It opens a metal edge (by the midpoint rule) only where that cannot
//! change connectivity: every opened edge touches a node outside the copper with no via, port or element on
//! it, and whose other metal neighbours stay connected without it.

use crate::geom::{BBox, P2};
use crate::mesh::{range_in, Grid};
use crate::model::{barrel_nodes, Layout, Lumped, SheetMask};

#[derive(Clone, Debug)]
pub struct EdgeOpts {
    /// (b): partly covered edges and faces get Dey–Mittra fractions; metal edges that only touch copper
    /// are opened
    pub conformal: bool,
    /// (a): edge-singularity factors from the static half-plane field
    pub singular: bool,
    /// supersampling per cell for face coverage and the line fit (edges use 4 times this)
    pub ss: usize,
    /// smallest open fraction used for an E edge (1/f is capped at 1/min_open); a metal edge with less open is
    /// kept metal, and a sheet face with less open area is closed (its H never changes)
    pub min_open: f64,
    /// singular factors are clamped to [kmin, kmax]
    pub kmin: f64,
    pub kmax: f64,
    /// components within this in-plane distance (cells) of a fitted edge line get the full singular
    /// factors; from there they blend linearly back to the uniform-field ones over `taper` cells
    pub reach: f64,
    pub taper: f64,
    /// also give singular factors to sheet-plane components (in-plane E, Hz) sampled over the metal; by
    /// default those keep the uniform-field (conformal) factors
    pub over_metal: bool,
    /// the same for H faces only (experiments)
    pub over_metal_h: bool,
    /// singular factors at full weight only where the nearest copper of another layer is at least this many
    /// in-plane cells away (none at half of it): the half-plane field ignores that conductor
    pub plane_cells: f64,
    /// radius (cells) of the disc used to fit the local edge line, and the misclassified share accepted
    pub fit_r: f64,
    pub fit_tol: f64,
    /// check the corrected operator's stability (power iteration) and lower dt if needed
    pub dt_check: bool,
    pub dt_iters: usize,
    /// log step of the permittivity factors (keeps the Metal material table small)
    pub quant: f64,
}

impl Default for EdgeOpts {
    fn default() -> Self {
        EdgeOpts { conformal: true, singular: true, ss: 8, min_open: 0.15, kmin: 0.25, kmax: 4.0, reach: 1.0, taper: 0.0, over_metal: false, over_metal_h: false, plane_cells: 4.0, fit_r: 1.25, fit_tol: 0.05, dt_check: true, dt_iters: 150, quant: 0.004 }
    }
}

impl EdgeOpts {
    /// From a spec value: `true`, or `{"conformal": .., "singular": .., "ss": .., "min_open": .., ...}`.
    pub fn from_json(v: &serde_json::Value) -> Option<EdgeOpts> {
        if v.as_bool() == Some(false) || v.is_null() {
            return None;
        }
        let mut o = EdgeOpts::default();
        if v.is_object() {
            let b = |k: &str, d: bool| v.get(k).and_then(|x| x.as_bool()).unwrap_or(d);
            let f = |k: &str, d: f64| v.get(k).and_then(|x| x.as_f64()).unwrap_or(d);
            o.conformal = b("conformal", o.conformal);
            o.singular = b("singular", o.singular);
            o.dt_check = b("dt_check", o.dt_check);
            o.dt_iters = f("dt_iters", o.dt_iters as f64) as usize;
            o.ss = f("ss", o.ss as f64).max(2.0) as usize;
            o.min_open = f("min_open", o.min_open);
            o.kmin = f("kmin", o.kmin);
            o.kmax = f("kmax", o.kmax);
            o.reach = f("reach", o.reach);
            o.taper = f("taper", o.taper);
            o.over_metal = b("over_metal", o.over_metal);
            o.over_metal_h = b("over_metal_h", o.over_metal_h);
            o.plane_cells = f("plane_cells", o.plane_cells);
            o.fit_r = f("fit_r", o.fit_r);
            o.fit_tol = f("fit_tol", o.fit_tol);
        }
        (o.conformal || o.singular).then_some(o)
    }
}

/// What `analyse` asks the model builder to change.
#[derive(Clone, Debug, Default)]
pub struct EdgeFix {
    /// E edges whose permittivity is redistributed: (axis, id, β per quadrant). The edge's new permittivity
    /// is Σ ε_q β_q over the four cells around it, in the order of `model::build` (b-, c-), (b-, c),
    /// (b, c-), (b, c) with b = axis + 1, c = axis + 2; a uniform field gives β_q = A_q / A'.
    pub e: Vec<(u8, u32, [f32; 4])>,
    /// metal edges (midpoint rule) that become open
    pub open: Vec<(u8, u32)>,
    /// open edges that become metal: the fourth side of a sliver face whose other three sides are metal
    /// (they already join its two nodes, so no nets merge)
    pub close: Vec<(u8, u32)>,
    /// H components whose update is divided by κ: (axis, id, κ)
    pub h: Vec<(u8, u32, f32)>,
    pub stats: EdgeStats,
}

#[derive(Clone, Debug, Default)]
pub struct EdgeStats {
    pub opened: usize,
    pub closed: usize,
    pub kept_for_connectivity: usize,
    pub e_conformal: usize,
    pub e_singular: usize,
    pub h_conformal: usize,
    pub h_singular: usize,
    pub fit_rejected: usize,
    pub dropped_shared: usize,
}

// ---------------------------------------------------------------------------------------------------
// static field of a conducting half-plane

/// A straight copper edge on the sheet plane z = z0 (mm): through (x0, y0), metal on the side of -n.
#[derive(Clone, Copy, Debug)]
pub struct Line {
    pub x0: f64,
    pub y0: f64,
    pub nx: f64,
    pub ny: f64,
    pub z0: f64,
}

impl Line {
    fn uw(&self, p: [f64; 3]) -> (f64, f64) {
        ((p[0] - self.x0) * self.nx + (p[1] - self.y0) * self.ny, p[2] - self.z0)
    }
    /// in-plane component of the normal along axis a (0 = x, 1 = y; z has none)
    fn n(&self, a: usize) -> f64 {
        match a {
            0 => self.nx,
            1 => self.ny,
            _ => 0.0,
        }
    }
}

/// √(u + iw), principal branch (cut along the metal, u < 0, w = 0); on the cut `side` picks the lip.
fn sqrtz(u: f64, w: f64, side: f64) -> (f64, f64) {
    let r = (u * u + w * w).sqrt();
    let re = (0.5 * (r + u)).max(0.0).sqrt();
    let im = (0.5 * (r - u)).max(0.0).sqrt();
    let s = if w > 0.0 {
        1.0
    } else if w < 0.0 {
        -1.0
    } else if side < 0.0 {
        -1.0
    } else {
        1.0
    };
    (re, s * im)
}

fn phi(u: f64, w: f64) -> f64 {
    sqrtz(u, w, 1.0).0
}

fn psi(u: f64, w: f64, side: f64) -> f64 {
    sqrtz(u, w, side).1
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Field {
    E,
    H,
}

/// Gauss-Legendre nodes and weights on [-1, 1] (16 points).
const GL: [(f64, f64); 8] = [
    (0.0950125098376374, 0.1894506104550685),
    (0.2816035507792589, 0.1826034150449236),
    (0.4580167776572274, 0.1691565193950025),
    (0.6178762444026438, 0.1495959888165767),
    (0.7554044083550030, 0.1246289712555339),
    (0.8656312023878318, 0.0951585116824928),
    (0.9445750230732326, 0.0622535239386479),
    (0.9894009349916499, 0.0271524594117541),
];

/// ∫_a^b f, split at the given points; each piece mapped by a smoothstep (whose derivative vanishes at
/// both ends) so square-root endpoint behaviour integrates accurately with 16 Gauss points.
fn quad(a: f64, b: f64, splits: &[f64], f: &dyn Fn(f64) -> f64) -> f64 {
    if !(b > a) {
        return 0.0;
    }
    let mut pts = vec![a];
    for &s in splits {
        if s > a + 1e-12 * (b - a) && s < b - 1e-12 * (b - a) {
            pts.push(s);
        }
    }
    pts.push(b);
    pts.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let mut sum = 0.0;
    for w in pts.windows(2) {
        let (lo, hi) = (w[0], w[1]);
        let l = hi - lo;
        for &(x, wt) in &GL {
            for t in [0.5 * (1.0 - x), 0.5 * (1.0 + x)] {
                let s = t * t * (3.0 - 2.0 * t);
                let ds = 6.0 * t * (1.0 - t);
                sum += 0.5 * wt * f(lo + l * s) * l * ds;
            }
        }
    }
    sum
}

/// Line integral of the singular field along axis `a` from `p` over length `len` (mm).
pub fn line_integral(fld: Field, ln: &Line, a: usize, p: [f64; 3], len: f64) -> f64 {
    let mut q = p;
    q[a] += len;
    let (ua, wa) = ln.uw(p);
    let (ub, wb) = ln.uw(q);
    match fld {
        // E = -∇φ, φ continuous
        Field::E => phi(ua, wa) - phi(ub, wb),
        // H = -∇ψ; ψ jumps across the metal, so a path through the sheet is split there
        Field::H => {
            if a == 2 && wa * wb < 0.0 && ua < 0.0 {
                let (sa, sb) = (wa.signum(), wb.signum());
                -((psi(ub, wb, sb) - psi(ua, 0.0, sb)) + (psi(ua, 0.0, sa) - psi(ua, wa, sa)))
            } else {
                let side = if wa + wb != 0.0 { (wa + wb).signum() } else { 1.0 };
                -(psi(ub, wb, side) - psi(ua, wa, side))
            }
        }
    }
}

/// Flux of the singular field through the axis-aligned rectangle with normal axis `na` at coordinate `c`,
/// spanning `rb` along axis b = na + 1 and `rc` along c = na + 2 (mm).
pub fn flux(fld: Field, ln: &Line, na: usize, c: f64, rb: (f64, f64), rc: (f64, f64)) -> f64 {
    let (b, cc) = ((na + 1) % 3, (na + 2) % 3);
    let at = |vb: f64, vc: f64| -> [f64; 3] {
        let mut p = [0.0; 3];
        p[na] = c;
        p[b] = vb;
        p[cc] = vc;
        p
    };
    // ranges by axis
    let range = |ax: usize| if ax == b { rb } else { rc };
    if na == 2 {
        // V·z = V_w = ∂G/∂u with G = ψ (E) or -φ (H); analytic along the in-plane axis with the larger
        // normal component, Gauss along the other
        let (pa, qa) = if ln.nx.abs() >= ln.ny.abs() { (0, 1) } else { (1, 0) };
        let np = ln.n(pa);
        let (p0, p1) = range(pa);
        let (q0, q1) = range(qa);
        let g = |pt: [f64; 3]| -> f64 {
            let (u, w) = ln.uw(pt);
            match fld {
                Field::E => psi(u, w, 1.0),
                Field::H => -phi(u, w),
            }
        };
        let pt = |vp: f64, vq: f64| -> [f64; 3] {
            let mut p = [0.0; 3];
            p[2] = c;
            p[pa] = vp;
            p[qa] = vq;
            p
        };
        // the integrand has kinks where the edge line meets the face's p0 / p1 sides
        let mut splits = Vec::new();
        let nq = ln.n(qa);
        if nq.abs() > 1e-12 {
            for &pv in &[p0, p1] {
                // u = 0: np (pv - xp) + nq (vq - xq) = 0
                let x0 = [ln.x0, ln.y0];
                splits.push(x0[qa] - np * (pv - x0[pa]) / nq);
            }
        }
        let f = |vq: f64| (g(pt(p1, vq)) - g(pt(p0, vq))) / np;
        let _ = at;
        quad(q0, q1, &splits, &f)
    } else {
        // in-plane normal: V·n = V_u n_na; the face spans the other in-plane axis p and z
        let nn = ln.n(na);
        if nn.abs() < 1e-12 {
            return 0.0;
        }
        let pa = 1 - na; // the other in-plane axis
        let (p0, p1) = range(pa);
        let (w0, w1) = range(2);
        let np = ln.n(pa);
        let pt = |vp: f64, vz: f64| -> [f64; 3] {
            let mut p = [0.0; 3];
            p[na] = c;
            p[pa] = vp;
            p[2] = vz;
            p
        };
        let s = if np.abs() >= 0.25 {
            // analytic along p: V_u = ∂G/∂u / ... with G = -φ (E) or -ψ (H); Gauss over z split at the sheet
            let g = |q: [f64; 3]| -> f64 {
                let (u, w) = ln.uw(q);
                match fld {
                    Field::E => -phi(u, w),
                    Field::H => -psi(u, w, if w >= 0.0 { 1.0 } else { -1.0 }),
                }
            };
            let f = |vz: f64| (g(pt(p1, vz)) - g(pt(p0, vz))) / np;
            quad(w0, w1, &[ln.z0], &f)
        } else {
            // analytic along z: V_u = ∂G/∂w with G = -ψ (E) or +φ (H); Gauss over p split where u = 0
            let f = |vp: f64| -> f64 {
                let (u, wa) = ln.uw(pt(vp, w0));
                let (_, wb) = ln.uw(pt(vp, w1));
                match fld {
                    Field::E => {
                        if wa < 0.0 && wb > 0.0 && u < 0.0 {
                            -((psi(u, wb, 1.0) - psi(u, 0.0, 1.0)) + (psi(u, 0.0, -1.0) - psi(u, wa, -1.0)))
                        } else {
                            let side = if wa + wb != 0.0 { (wa + wb).signum() } else { 1.0 };
                            -(psi(u, wb, side) - psi(u, wa, side))
                        }
                    }
                    Field::H => phi(u, wb) - phi(u, wa),
                }
            };
            let mut splits = Vec::new();
            if np.abs() > 1e-12 {
                let x0 = [ln.x0, ln.y0];
                let nnv = ln.n(na);
                // u = nn (c - x0[na]) + np (vp - x0[pa]) = 0
                splits.push(x0[pa] - nnv * (c - x0[na]) / np);
            }
            quad(p0, p1, &splits, &f)
        };
        s * nn
    }
}

// ---------------------------------------------------------------------------------------------------
// coverage of one sheet plane

struct Plane {
    k: usize,
    /// fine face samples: centres of an `s` x `s` subdivision of every cell
    s: usize,
    wf: usize,
    xf: Vec<f64>,
    yf: Vec<f64>,
    cov: Vec<bool>,
    /// open fraction of Ex(i, j) / Ey(i, j), index j * nx + i (1 = no copper)
    fx: Vec<f32>,
    fy: Vec<f32>,
    /// open area fraction of the face of cell (i, j)
    ff: Vec<f32>,
    /// node inside copper
    node_in: Vec<bool>,
}

fn sub_centres(l: &[f64], s: usize) -> Vec<f64> {
    let mut v = Vec::with_capacity(l.len().saturating_sub(1) * s);
    for w in l.windows(2) {
        for q in 0..s {
            v.push(w[0] + (q as f64 + 0.5) / s as f64 * (w[1] - w[0]));
        }
    }
    v
}

/// Supersampled true copper of one sheet (inside the copper window, if any).
fn plane_true(g: &Grid, lay: &Layout, si: usize, k: usize, s: usize) -> Plane {
    let (nx, ny) = (g.nx(), g.ny());
    let xmm: Vec<f64> = g.x.iter().map(|v| v * 1e3).collect();
    let ymm: Vec<f64> = g.y.iter().map(|v| v * 1e3).collect();
    let xf = sub_centres(&xmm, s);
    let yf = sub_centres(&ymm, s);
    let se = 4 * s;
    let xe = sub_centres(&xmm, se);
    let ye = sub_centres(&ymm, se);
    let (wf, hf) = (xf.len(), yf.len());
    let mut cov = vec![false; wf * hf];
    let mut ex = vec![false; xe.len() * ny];
    let mut ey = vec![false; nx * ye.len()];
    let mut nd = vec![false; nx * ny];
    let sheet = &lay.sheets[si];
    for p in &sheet.prims {
        if lay.copper_window.is_some_and(|w| !w.overlaps(&p.shape.bbox())) {
            continue;
        }
        let d = p.dark;
        crate::model::paint_with(&xf, &yf, p, |q| cov[q] = d);
        crate::model::paint_with(&xe, &ymm, p, |q| ex[q] = d);
        crate::model::paint_with(&xmm, &ye, p, |q| ey[q] = d);
        crate::model::paint_with(&xmm, &ymm, p, |q| nd[q] = d);
    }
    let mut fx = vec![1.0f32; nx * ny];
    let mut fy = vec![1.0f32; nx * ny];
    let mut ff = vec![1.0f32; nx * ny];
    for j in 0..ny {
        for i in 0..nx.saturating_sub(1) {
            let c = (0..se).filter(|&q| ex[j * xe.len() + i * se + q]).count();
            fx[j * nx + i] = 1.0 - c as f32 / se as f32;
        }
    }
    for j in 0..ny.saturating_sub(1) {
        for i in 0..nx {
            let c = (0..se).filter(|&q| ey[(j * se + q) * nx + i]).count();
            fy[j * nx + i] = 1.0 - c as f32 / se as f32;
        }
    }
    for j in 0..ny.saturating_sub(1) {
        for i in 0..nx.saturating_sub(1) {
            let mut c = 0;
            for b in 0..s {
                for a in 0..s {
                    c += cov[(j * s + b) * wf + i * s + a] as usize;
                }
            }
            ff[j * nx + i] = 1.0 - c as f32 / (s * s) as f32;
        }
    }
    // A straight copper edge that crosses both y edges of a face (or both x edges) covers the mean of
    // their covered fractions: exact, and consistent with the edge fractions (which are sampled 4 times
    // finer than the face). Otherwise keep the sampled area.
    for j in 0..ny.saturating_sub(1) {
        for i in 0..nx.saturating_sub(1) {
            let q = j * nx + i;
            let part = |f: f32| f > 0.0 && f < 1.0;
            let (ya, yb) = (fy[q], fy[q + 1]);
            let (xa, xb) = (fx[q], fx[q + nx]);
            let edge_x = fx[q] == fx[q + nx] && (fx[q] == 0.0 || fx[q] == 1.0);
            let edge_y = fy[q] == fy[q + 1] && (fy[q] == 0.0 || fy[q] == 1.0);
            if part(ya) && part(yb) && !part(xa) && !part(xb) {
                ff[q] = 0.5 * (ya + yb);
            } else if part(xa) && part(xb) && !part(ya) && !part(yb) {
                ff[q] = 0.5 * (xa + xb);
            } else if part(ya) && part(yb) && edge_x {
                ff[q] = 0.5 * (ya + yb);
            } else if part(xa) && part(xb) && edge_y {
                ff[q] = 0.5 * (xa + xb);
            }
        }
    }
    Plane { k, s, wf, xf, yf, cov, fx, fy, ff, node_in: nd }
}

/// The staircased copper the midpoint rule produced: a face is copper when its four edges are metal.
fn plane_stair(g: &Grid, m: &SheetMask, s: usize) -> Plane {
    let (nx, ny) = (g.nx(), g.ny());
    let xmm: Vec<f64> = g.x.iter().map(|v| v * 1e3).collect();
    let ymm: Vec<f64> = g.y.iter().map(|v| v * 1e3).collect();
    let xf = sub_centres(&xmm, s);
    let yf = sub_centres(&ymm, s);
    let wf = xf.len();
    let mut cov = vec![false; wf * yf.len()];
    let mut ff = vec![1.0f32; nx * ny];
    for j in 0..ny.saturating_sub(1) {
        for i in 0..nx.saturating_sub(1) {
            let q = j * nx + i;
            let full = m.ex[q] && m.ex[q + nx] && m.ey[q] && m.ey[q + 1];
            if full {
                ff[q] = 0.0;
                for b in 0..s {
                    for a in 0..s {
                        cov[(j * s + b) * wf + i * s + a] = true;
                    }
                }
            }
        }
    }
    let fx = m.ex.iter().map(|&v| if v { 0.0 } else { 1.0 }).collect();
    let fy = m.ey.iter().map(|&v| if v { 0.0 } else { 1.0 }).collect();
    let node_in = vec![false; nx * ny];
    Plane { k: m.k, s, wf, xf, yf, cov, fx, fy, ff, node_in }
}

/// Fit a straight edge to the copper within radius `r` (mm) of `(px, py)`: the direction from the
/// copper's centroid to the open part's centroid is the normal (exact for a line through a disc), the
/// copper's area share places the line. Returns the line and the misclassified share of the samples.
fn fit(pl: &Plane, xmm: &[f64], ymm: &[f64], px: f64, py: f64, r: f64, z0: f64) -> Option<(Line, f64)> {
    let s = pl.s;
    let ir = range_in(xmm, px - r, px + r);
    let jr = range_in(ymm, py - r, py + r);
    let i0 = ir.start.saturating_sub(1) * s;
    let i1 = (ir.end * s).min(pl.wf);
    let j0 = jr.start.saturating_sub(1) * s;
    let j1 = (jr.end * s).min(pl.yf.len());
    let r2 = r * r;
    let (mut n1, mut n0) = (0usize, 0usize);
    let (mut c1, mut c0) = ([0.0f64; 2], [0.0f64; 2]);
    for jf in j0..j1 {
        let dy = pl.yf[jf] - py;
        for ifn in i0..i1 {
            let dx = pl.xf[ifn] - px;
            if dx * dx + dy * dy > r2 {
                continue;
            }
            if pl.cov[jf * pl.wf + ifn] {
                n1 += 1;
                c1[0] += dx;
                c1[1] += dy;
            } else {
                n0 += 1;
                c0[0] += dx;
                c0[1] += dy;
            }
        }
    }
    if n1 == 0 || n0 == 0 {
        return None;
    }
    let (m1, m0) = ([c1[0] / n1 as f64, c1[1] / n1 as f64], [c0[0] / n0 as f64, c0[1] / n0 as f64]);
    let (dx, dy) = (m0[0] - m1[0], m0[1] - m1[1]);
    let l = (dx * dx + dy * dy).sqrt();
    if l <= 0.0 {
        return None;
    }
    let (nxv, nyv) = (dx / l, dy / l);
    // copper share of the disc -> line offset t along n from the centre: area of {x < t} is
    // r^2 acos(-t/r) + t sqrt(r^2 - t^2)
    let a1 = n1 as f64 / (n1 + n0) as f64;
    let area = |t: f64| -> f64 { ((-t / r).clamp(-1.0, 1.0).acos() + (t / r) * (1.0 - (t / r).powi(2)).max(0.0).sqrt()) / std::f64::consts::PI };
    let (mut lo, mut hi) = (-r, r);
    for _ in 0..60 {
        let m = 0.5 * (lo + hi);
        if area(m) < a1 {
            lo = m;
        } else {
            hi = m;
        }
    }
    let t = 0.5 * (lo + hi);
    let mut mis = 0usize;
    for jf in j0..j1 {
        let dy = pl.yf[jf] - py;
        for ifn in i0..i1 {
            let dx = pl.xf[ifn] - px;
            if dx * dx + dy * dy > r2 {
                continue;
            }
            let u = dx * nxv + dy * nyv;
            if pl.cov[jf * pl.wf + ifn] != (u < t) {
                mis += 1;
            }
        }
    }
    Some((Line { x0: px + t * nxv, y0: py + t * nyv, nx: nxv, ny: nyv, z0 }, mis as f64 / (n1 + n0) as f64))
}

// ---------------------------------------------------------------------------------------------------

/// Nodes (i, j) of plane k that carry something besides sheet copper: via barrels, lumped elements.
fn attached(g: &Grid, lay: &Layout, lumped: &[Lumped], k: usize) -> Vec<bool> {
    let (nx, ny) = (g.nx(), g.ny());
    let mut at = vec![false; nx * ny];
    for b in &lay.barrels {
        let (cells, k0, k1) = barrel_nodes(g, b);
        if k >= k0 && k <= k1 {
            for (i, j) in cells {
                at[j * nx + i] = true;
            }
        }
    }
    let xs = [&g.x, &g.y, &g.z];
    for l in lumped {
        let mut r = [(0usize, 0usize); 3];
        for b in 0..3 {
            let (lo, hi) = (l.lo[b].min(l.hi[b]) * 1e-3, l.lo[b].max(l.hi[b]) * 1e-3);
            r[b] = if b == l.axis {
                (crate::mesh::snap(xs[b], lo), crate::mesh::snap(xs[b], hi))
            } else {
                let rr = range_in(xs[b], lo - 1e-12, hi + 1e-12);
                if rr.is_empty() {
                    let s = crate::mesh::snap(xs[b], 0.5 * (lo + hi));
                    (s, s)
                } else {
                    (rr.start, rr.end - 1)
                }
            };
        }
        if k >= r[2].0 && k <= r[2].1 {
            for j in r[1].0..=r[1].1.min(ny - 1) {
                for i in r[0].0..=r[0].1.min(nx - 1) {
                    at[j * nx + i] = true;
                }
            }
        }
    }
    at
}

struct Uf(Vec<u32>);
impl Uf {
    fn new(n: usize) -> Uf {
        Uf((0..n as u32).collect())
    }
    fn find(&mut self, mut a: u32) -> u32 {
        while self.0[a as usize] != a {
            let p = self.0[self.0[a as usize] as usize];
            self.0[a as usize] = p;
            a = p;
        }
        a
    }
    fn union(&mut self, a: u32, b: u32) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.0[ra.max(rb) as usize] = ra.min(rb);
        }
    }
}

/// Which mask-metal in-plane edges of a plane can be opened: returns per-edge flags (ex, ey) and the
/// nodes that keep a staircase edge (no corrections around them), and counts.
fn openable(g: &Grid, m: &SheetMask, pl: &Plane, att: &[bool], o: &EdgeOpts) -> (Vec<bool>, Vec<bool>, Vec<bool>, usize) {
    let (nx, ny) = (g.nx(), g.ny());
    let cand_x: Vec<bool> = (0..nx * ny).map(|q| m.ex[q] && pl.fx[q] as f64 >= o.min_open).collect();
    let cand_y: Vec<bool> = (0..nx * ny).map(|q| m.ey[q] && pl.fy[q] as f64 >= o.min_open).collect();
    // incident metal edges of node (i, j): (is_x, q)
    let inc = |i: usize, j: usize| -> [(bool, Option<usize>); 4] {
        [
            (true, (i + 1 < nx).then_some(j * nx + i)),
            (true, (i > 0).then(|| j * nx + i - 1)),
            (false, (j + 1 < ny).then_some(j * nx + i)),
            (false, (j > 0).then(|| (j - 1) * nx + i)),
        ]
    };
    let metal = |x: bool, q: usize| if x { m.ex[q] } else { m.ey[q] };
    let cand = |x: bool, q: usize| if x { cand_x[q] } else { cand_y[q] };
    // purely fringe nodes: outside copper, nothing attached, every metal edge a candidate
    let fringe: Vec<bool> = (0..nx * ny)
        .map(|n| {
            let (i, j) = (n % nx, n / nx);
            if pl.node_in[n] || att[n] {
                return false;
            }
            let mut any = false;
            for (x, q) in inc(i, j) {
                if let Some(q) = q {
                    if metal(x, q) {
                        if !cand(x, q) {
                            return false;
                        }
                        any = true;
                    }
                }
            }
            any
        })
        .collect();
    // connectivity without the fringe nodes' edges
    let other = |x: bool, q: usize| -> (usize, usize) {
        let (i, j) = (q % nx, q / nx);
        if x {
            (j * nx + i, j * nx + i + 1)
        } else {
            (j * nx + i, (j + 1) * nx + i)
        }
    };
    let mut uf = Uf::new(nx * ny);
    for q in 0..nx * ny {
        for x in [true, false] {
            if metal(x, q) {
                let (a, b) = other(x, q);
                if b < nx * ny && !fringe[a] && !fringe[b] {
                    uf.union(a as u32, b as u32);
                }
            }
        }
    }
    let mut keep = vec![false; nx * ny];
    let mut kept = 0;
    for n in 0..nx * ny {
        if !fringe[n] {
            continue;
        }
        let (i, j) = (n % nx, n / nx);
        let mut nb: [u32; 4] = [u32::MAX; 4];
        let mut c = 0;
        for (x, q) in inc(i, j) {
            if let Some(q) = q {
                if metal(x, q) {
                    let (a, b) = other(x, q);
                    let o = if a == n { b } else { a };
                    nb[c] = o as u32;
                    c += 1;
                }
            }
        }
        let r0 = uf.find(nb[0]);
        let split = (1..c).any(|t| uf.find(nb[t]) != r0);
        if split {
            keep[n] = true;
            kept += 1;
            for t in 0..c {
                uf.union(n as u32, nb[t]);
            }
        }
    }
    let mut ox = vec![false; nx * ny];
    let mut oy = vec![false; nx * ny];
    for q in 0..nx * ny {
        for x in [true, false] {
            if cand(x, q) {
                let (a, b) = other(x, q);
                let open_a = fringe[a] && !keep[a];
                let open_b = b < nx * ny && fringe[b] && !keep[b];
                if open_a || open_b {
                    if x {
                        ox[q] = true;
                    } else {
                        oy[q] = true;
                    }
                }
            }
        }
    }
    (ox, oy, keep, kept)
}

/// Corrections for every sheet. `masks` are the midpoint-rule masks after net conflict resolution.
pub fn analyse(g: &Grid, lay: &Layout, masks: &[SheetMask], lumped: &[Lumped], pml: &[usize; 6], o: &EdgeOpts) -> EdgeFix {
    let (nx, ny, nz) = (g.nx(), g.ny(), g.nz());
    let xmm: Vec<f64> = g.x.iter().map(|v| v * 1e3).collect();
    let ymm: Vec<f64> = g.y.iter().map(|v| v * 1e3).collect();
    let zmm: Vec<f64> = g.z.iter().map(|v| v * 1e3).collect();
    let d = |l: &[f64], i: usize| if i + 1 < l.len() { l[i + 1] - l[i] } else { 0.0 };
    let mid = |l: &[f64], i: usize| 0.5 * (l[i] + l[i + 1]);
    // dual extent around node i: [l[i] - d(i-1)/2, l[i] + d(i)/2]
    let dual = |l: &[f64], i: usize| -> (f64, f64) {
        let lo = if i > 0 { 0.5 * (l[i] - l[i - 1]) } else { 0.0 };
        let hi = if i + 1 < l.len() { 0.5 * (l[i + 1] - l[i]) } else { 0.0 };
        (l[i] - lo, l[i] + hi)
    };
    let inside = |i: usize, j: usize, k: usize| i > pml[0] && i + pml[1] + 1 < nx && j > pml[2] && j + pml[3] + 1 < ny && k > pml[4] && k + pml[5] + 1 < nz;
    let win: Option<BBox> = lay.copper_window;
    let in_win = |x: f64, y: f64| win.map_or(true, |w| w.contains(P2::new(x, y)));
    // E edges of a lumped element's own axis inside its box keep the element's material
    let in_lumped = |a: usize, p: [f64; 3], tol: [f64; 3]| -> bool {
        lumped.iter().any(|l| {
            l.axis == a
                && (0..3).all(|b| {
                    let (lo, hi) = (l.lo[b].min(l.hi[b]), l.lo[b].max(l.hi[b]));
                    p[b] >= lo - tol[b] && p[b] <= hi + tol[b]
                })
        })
    };
    let mut fix = EdgeFix::default();
    let mut ecomp: Vec<(u8, u32, [f32; 4])> = Vec::new();
    let mut hcomp: Vec<(u8, u32, f32)> = Vec::new();
    for (si, m) in masks.iter().enumerate() {
        let k = m.k;
        if k == 0 || k + 1 >= nz {
            continue;
        }
        let z0 = zmm[k];
        let mut pt = if o.conformal { plane_true(g, lay, si, k, o.ss) } else { plane_stair(g, m, o.ss) };
        let att = attached(g, lay, lumped, k);
        // vertical distance (mm) from each node to the nearest copper of another sheet: the half-plane field
        // ignores that conductor, so the singular factors fade out when it is within a few cells
        // a PEC wall (no absorber on that z side) is a plane everywhere
        let mut wall = f64::INFINITY;
        if pml[4] == 0 {
            wall = wall.min(z0 - zmm[0]);
        }
        if pml[5] == 0 {
            wall = wall.min(zmm[nz - 1] - z0);
        }
        let mut gap_z = vec![wall; nx * ny];
        for m2 in masks.iter().filter(|m2| m2.k != k) {
            let dzk = (zmm[m2.k] - z0).abs();
            for j in 0..ny {
                for i in 0..nx {
                    let q = j * nx + i;
                    let on = m2.ex[q] || m2.ey[q] || (i > 0 && m2.ex[q - 1]) || (j > 0 && m2.ey[q - nx]);
                    if on && dzk < gap_z[q] {
                        gap_z[q] = dzk;
                    }
                }
            }
        }
        let (ox, oy, keep, kept) = if o.conformal { openable(g, m, &pt, &att, o) } else { (vec![false; nx * ny], vec![false; nx * ny], vec![false; nx * ny], 0) };
        fix.stats.kept_for_connectivity += kept;
        let stair_near = |i: usize, j: usize| -> bool {
            // any kept staircase node among (i..=i+1, j..=j+1)
            (0..2).any(|a| (0..2).any(|b| i + a < nx && j + b < ny && keep[(j + b) * nx + i + a]))
        };
        // slivers: a face with less than min_open open area and three metal sides gets the fourth
        let mut cx = vec![false; nx * ny];
        let mut cy = vec![false; nx * ny];
        if o.conformal {
            for j in 0..ny.saturating_sub(1) {
                for i in 0..nx.saturating_sub(1) {
                    let q = j * nx + i;
                    if (pt.ff[q] as f64) >= o.min_open || keep[q] {
                        continue;
                    }
                    let mx = |q: usize| m.ex[q] && !ox[q];
                    let my = |q: usize| m.ey[q] && !oy[q];
                    let sides = [(true, q, mx(q)), (true, q + nx, mx(q + nx)), (false, q, my(q)), (false, q + 1, my(q + 1))];
                    let nm = sides.iter().filter(|s| s.2).count();
                    if nm == 3 {
                        let &(isx, qq, _) = sides.iter().find(|s| !s.2).unwrap();
                        if isx {
                            cx[qq] = true;
                        } else {
                            cy[qq] = true;
                        }
                        // the closed face is copper now: the singular fit must see the edge where the
                        // discrete metal ends
                        let s = pt.s;
                        for b in 0..s {
                            for a in 0..s {
                                pt.cov[(j * s + b) * pt.wf + i * s + a] = true;
                            }
                        }
                        pt.ff[q] = 0.0;
                    }
                }
            }
        }
        // final metal state of the in-plane edges after opening and closing
        let metal_x = |q: usize| (m.ex[q] && !ox[q]) || cx[q];
        let metal_y = |q: usize| (m.ey[q] && !oy[q]) || cy[q];
        for q in 0..nx * ny {
            let (i, j) = (q % nx, q / nx);
            if !inside(i, j, k) || !in_win(xmm[i], ymm[j]) {
                continue;
            }
            if ox[q] {
                fix.open.push((0, g.idx(i, j, k) as u32));
            }
            if oy[q] {
                fix.open.push((1, g.idx(i, j, k) as u32));
            }
            if cx[q] {
                fix.close.push((0, g.idx(i, j, k) as u32));
            }
            if cy[q] {
                fix.close.push((1, g.idx(i, j, k) as u32));
            }
        }
        // candidate cells: faces near a copper boundary (in the fine coverage)
        let mut bnd = vec![false; nx * ny];
        for j in 0..ny.saturating_sub(1) {
            for i in 0..nx.saturating_sub(1) {
                let f = pt.ff[j * nx + i];
                let mixed = f > 0.0 && f < 1.0;
                let diff = (i + 2 < nx && pt.ff[j * nx + i + 1] != f) || (j + 2 < ny && pt.ff[(j + 1) * nx + i] != f);
                if mixed || diff {
                    for b in j.saturating_sub(2)..(j + 3).min(ny) {
                        for a in i.saturating_sub(2)..(i + 3).min(nx) {
                            bnd[b * nx + a] = true;
                        }
                    }
                }
            }
        }
        let hloc = |i: usize, j: usize| 0.5 * (d(&xmm, i.min(nx - 2)) + d(&ymm, j.min(ny - 2)));
        // the singular fit at an in-plane point; None: no edge near, Some(None): fit rejected
        let mut rejected = 0usize;
        let mut line_at = |px: f64, py: f64, h: f64| -> Option<Option<(Line, f64)>> {
            if !o.singular {
                return Some(None);
            }
            let r = (o.fit_r + o.taper).max(o.reach + 0.25) * h;
            match fit(&pt, &xmm, &ymm, px, py, r, z0) {
                None => None,
                Some((ln, err)) => {
                    let u = ((px - ln.x0) * ln.nx + (py - ln.y0) * ln.ny).abs() / h;
                    if u > o.reach + o.taper {
                        None
                    } else if err > o.fit_tol {
                        rejected += 1;
                        Some(None)
                    } else {
                        // weight of the singular factors: 1 within reach, then linear to 0
                        let t = if u <= o.reach { 1.0 } else { 1.0 - (u - o.reach) / o.taper };
                        // fade with the nearest other copper below / above: full weight from `plane_cells`
                        // in-plane cells away, none at half that
                        let (ii, jj) = (crate::mesh::snap(&xmm, px).min(nx - 1), crate::mesh::snap(&ymm, py).min(ny - 1));
                        let hz = gap_z[jj * nx + ii] / h;
                        let pc = o.plane_cells;
                        let tg = ((hz - 0.5 * pc) / (0.5 * pc)).clamp(0.0, 1.0);
                        if t * tg <= 0.0 {
                            Some(None)
                        } else {
                            Some(Some((ln, t * tg)))
                        }
                    }
                }
            }
        };
        let blend = |k: f64, t: f64| 1.0 + t * (k - 1.0);
        let clampk = |v: f64| v.clamp(o.kmin, o.kmax);
        // ---- E edges: Ex, Ey on the sheet; Ez just above and below
        for j in 1..ny.saturating_sub(1) {
            for i in 1..nx.saturating_sub(1) {
                if !bnd[j * nx + i] && !bnd[j * nx + i - 1] && !bnd[(j - 1) * nx + i] && !bnd[(j - 1) * nx + i - 1] {
                    continue;
                }
                let h = hloc(i, j);
                // in-plane edges
                for a in 0..2usize {
                    let q = j * nx + i;
                    let (metal, fo, isz) = if a == 0 { (metal_x(q), pt.fx[q] as f64, i + 1 >= nx) } else { (metal_y(q), pt.fy[q] as f64, j + 1 >= ny) };
                    if metal || isz || !inside(i, j, k) || stair_near(i.saturating_sub(1), j.saturating_sub(1)) || stair_near(i, j) {
                        continue;
                    }
                    let (len, pa) = if a == 0 { (d(&xmm, i), [xmm[i], ymm[j], z0]) } else { (d(&ymm, j), [xmm[i], ymm[j], z0]) };
                    let mut pm = pa;
                    pm[a] += 0.5 * len;
                    if !in_win(pm[0], pm[1]) || in_lumped(a, pm, [0.6 * h, 0.6 * h, 0.6 * d(&zmm, k)]) {
                        continue;
                    }
                    let (b, c) = ((a + 1) % 3, (a + 2) % 3);
                    let ls = [&xmm, &ymm, &zmm];
                    let node = [i, j, k];
                    let (db, dc) = (dual(ls[b], node[b]), dual(ls[c], node[c]));
                    let area = (db.1 - db.0) * (dc.1 - dc.0);
                    let quads = |bq: usize, cq: usize| -> ((f64, f64), (f64, f64)) {
                        let rb = if bq == 0 { (db.0, ls[b][node[b]]) } else { (ls[b][node[b]], db.1) };
                        let rc = if cq == 0 { (dc.0, ls[c][node[c]]) } else { (ls[c][node[c]], dc.1) };
                        (rb, rc)
                    };
                    let mut beta = [0f32; 4];
                    let mut done = false;
                    if let Some(ln) = line_at(pm[0], pm[1], h) {
                        let open_side = |l: &Line| o.over_metal || (pm[0] - l.x0) * l.nx + (pm[1] - l.y0) * l.ny >= 0.0;
                        if let Some((ln, t)) = ln.filter(|l| l.0.n(a).abs() >= 0.25 && open_side(&l.0)) {
                            let e = line_integral(Field::E, &ln, a, pa, len);
                            let mut fl = [0.0f64; 4];
                            for bq in 0..2 {
                                for cq in 0..2 {
                                    let (rb, rc) = quads(bq, cq);
                                    fl[bq * 2 + cq] = flux(Field::E, &ln, a, pm[a], rb, rc);
                                }
                            }
                            let tot: f64 = fl.iter().sum();
                            if e.abs() > 1e-9 && tot / e > 0.0 {
                                let kap = clampk(tot * len / (area * e));
                                for bq in 0..2 {
                                    for cq in 0..2 {
                                        let (rb, rc) = quads(bq, cq);
                                        let q = bq * 2 + cq;
                                        let uni = (rb.1 - rb.0) * (rc.1 - rc.0) / area;
                                        beta[q] = (t * fl[q] / tot * kap + (1.0 - t) * uni) as f32;
                                    }
                                }
                                done = true;
                                fix.stats.e_singular += 1;
                            }
                        }
                        if !done && o.conformal && fo < 1.0 {
                            let kap = 1.0 / fo.max(o.min_open);
                            for bq in 0..2 {
                                for cq in 0..2 {
                                    let (rb, rc) = quads(bq, cq);
                                    beta[bq * 2 + cq] = ((rb.1 - rb.0) * (rc.1 - rc.0) / area * kap) as f32;
                                }
                            }
                            done = true;
                            fix.stats.e_conformal += 1;
                        }
                    } else if o.conformal && fo < 1.0 {
                        let kap = 1.0 / fo.max(o.min_open);
                        for bq in 0..2 {
                            for cq in 0..2 {
                                let (rb, rc) = quads(bq, cq);
                                beta[bq * 2 + cq] = ((rb.1 - rb.0) * (rc.1 - rc.0) / area * kap) as f32;
                            }
                        }
                        done = true;
                        fix.stats.e_conformal += 1;
                    }
                    if done {
                        ecomp.push((a as u8, g.idx(i, j, k) as u32, beta));
                    }
                }
                // Ez above (k) and below (k - 1)
                if !o.singular || stair_near(i.saturating_sub(1), j.saturating_sub(1)) || stair_near(i, j) || !in_win(xmm[i], ymm[j]) {
                    continue;
                }
                let Some(Some((ln, tw))) = line_at(xmm[i], ymm[j], h) else { continue };
                for kk in [k, k - 1] {
                    if !inside(i, j, kk) || kk + 1 >= nz {
                        continue;
                    }
                    let len = d(&zmm, kk);
                    let pa = [xmm[i], ymm[j], zmm[kk]];
                    let pm = [xmm[i], ymm[j], zmm[kk] + 0.5 * len];
                    if in_lumped(2, pm, [0.6 * h, 0.6 * h, 0.0]) {
                        continue;
                    }
                    let (dxr, dyr) = (dual(&xmm, i), dual(&ymm, j));
                    let area = (dxr.1 - dxr.0) * (dyr.1 - dyr.0);
                    let e = line_integral(Field::E, &ln, 2, pa, len);
                    let mut fl = [0.0f64; 4];
                    // quadrant order (b = x, c = y): (x-, y-), (x-, y), (x, y-), (x, y)
                    for bq in 0..2 {
                        for cq in 0..2 {
                            let rb = if bq == 0 { (dxr.0, xmm[i]) } else { (xmm[i], dxr.1) };
                            let rc = if cq == 0 { (dyr.0, ymm[j]) } else { (ymm[j], dyr.1) };
                            fl[bq * 2 + cq] = flux(Field::E, &ln, 2, pm[2], rb, rc);
                        }
                    }
                    let tot: f64 = fl.iter().sum();
                    if e.abs() > 1e-9 && tot / e > 0.0 {
                        let kap = clampk(tot * len / (area * e));
                        let mut beta = [0f32; 4];
                        for t in 0..4 {
                            // the four quadrants of an Ez dual face have equal shares of a uniform field only
                            // on a uniform mesh; take the actual ones
                            let (bq, cq) = (t / 2, t % 2);
                            let wb = if bq == 0 { xmm[i] - dxr.0 } else { dxr.1 - xmm[i] };
                            let wc = if cq == 0 { ymm[j] - dyr.0 } else { dyr.1 - ymm[j] };
                            let uni = wb * wc / area;
                            beta[t] = (tw * fl[t] / tot * kap + (1.0 - tw) * uni) as f32;
                        }
                        ecomp.push((2, g.idx(i, j, kk) as u32, beta));
                        fix.stats.e_singular += 1;
                    }
                }
            }
        }
        // ---- H: Hz on the sheet; Hx, Hy half a cell above and below
        for j in 1..ny.saturating_sub(2) {
            for i in 1..nx.saturating_sub(2) {
                if !bnd[j * nx + i] {
                    continue;
                }
                let h = hloc(i, j);
                let q = j * nx + i;
                // Hz(i, j, k): face of cell (i, j); its edges Ex(j), Ex(j+1), Ey(i), Ey(i+1)
                let nmetal = [metal_x(q), metal_x(q + nx), metal_y(q), metal_y(q + 1)].iter().filter(|&&v| v).count();
                let (xm, ym) = (mid(&xmm, i), mid(&ymm, j));
                if inside(i, j, k) && nmetal > 0 && nmetal < 4 && !stair_near(i, j) && in_win(xm, ym) {
                    let area = d(&xmm, i) * d(&ymm, j);
                    let dz = dual(&zmm, k);
                    let mut kap = None;
                    if let Some(ln) = line_at(xm, ym, h) {
                        if let Some((ln, t)) = ln.filter(|l| o.over_metal || o.over_metal_h || (xm - l.0.x0) * l.0.nx + (ym - l.0.y0) * l.0.ny >= 0.0) {
                            let fl = flux(Field::H, &ln, 2, z0, (xmm[i], xmm[i + 1]), (ymm[j], ymm[j + 1]));
                            let hl = line_integral(Field::H, &ln, 2, [xm, ym, dz.0], dz.1 - dz.0);
                            if hl.abs() > 1e-9 && fl / hl > 0.0 {
                                // the uniform-field factor of this face is its open fraction
                                let uni = if o.conformal && pt.ff[q] > 0.0 { (pt.ff[q] as f64).max(o.min_open) } else { 1.0 };
                                // (faces this small are sampled over the metal and take the branch below)
                                let ks = clampk(fl * (dz.1 - dz.0) / (area * hl));
                                kap = Some((t * ks + (1.0 - t) * uni, true));
                            }
                        }
                        if kap.is_none() && o.conformal && pt.ff[q] < 1.0 && pt.ff[q] > 0.0 {
                            kap = Some((if (pt.ff[q] as f64) < o.min_open { f64::INFINITY } else { pt.ff[q] as f64 }, false));
                        }
                    } else if o.conformal && pt.ff[q] < 1.0 && pt.ff[q] > 0.0 {
                        kap = Some((if (pt.ff[q] as f64) < o.min_open { f64::INFINITY } else { pt.ff[q] as f64 }, false));
                    }
                    if let Some((kv, sing)) = kap {
                        hcomp.push((2, g.idx(i, j, k) as u32, kv as f32));
                        if sing {
                            fix.stats.h_singular += 1;
                        } else {
                            fix.stats.h_conformal += 1;
                        }
                    }
                }
                if !o.singular {
                    continue;
                }
                // Hx(i, j, kk) at (x_i, y_{j+1/2}); Hy(i, j, kk) at (x_{i+1/2}, y_j)
                for (a, px, py, sn) in [(0usize, xmm[i], ym, stair_near(i.saturating_sub(1), j)), (1usize, xm, ymm[j], stair_near(i, j.saturating_sub(1)))] {
                    if sn || !in_win(px, py) {
                        continue;
                    }
                    let Some(Some((ln, tw))) = line_at(px, py, h) else { continue };
                    if ln.n(a).abs() < 0.25 {
                        // the singular H has (almost) no component along a when the edge runs along a
                        continue;
                    }
                    for kk in [k, k - 1] {
                        if !inside(i, j, kk) || kk + 1 >= nz {
                            continue;
                        }
                        let (zl, zh) = (zmm[kk], zmm[kk + 1]);
                        let zm = 0.5 * (zl + zh);
                        let (fl, hl, area, dl) = if a == 0 {
                            let du = dual(&xmm, i);
                            (
                                flux(Field::H, &ln, 0, xmm[i], (ymm[j], ymm[j + 1]), (zl, zh)),
                                line_integral(Field::H, &ln, 0, [du.0, ym, zm], du.1 - du.0),
                                d(&ymm, j) * (zh - zl),
                                du.1 - du.0,
                            )
                        } else {
                            let du = dual(&ymm, j);
                            (
                                flux(Field::H, &ln, 1, ymm[j], (zl, zh), (xmm[i], xmm[i + 1])),
                                line_integral(Field::H, &ln, 1, [xm, du.0, zm], du.1 - du.0),
                                d(&xmm, i) * (zh - zl),
                                du.1 - du.0,
                            )
                        };
                        if hl.abs() > 1e-9 && fl / hl > 0.0 {
                            hcomp.push((a as u8, g.idx(i, j, kk) as u32, blend(clampk(fl * dl / (area * hl)), tw) as f32));
                            fix.stats.h_singular += 1;
                        }
                    }
                }
            }
        }
        fix.stats.fit_rejected += rejected;
    }
    // a component claimed by two sheets (planes a cell apart) gets no correction
    let dedup_e = |v: &mut Vec<(u8, u32, [f32; 4])>| -> usize {
        v.sort_by_key(|t| (t.0, t.1));
        let mut out: Vec<(u8, u32, [f32; 4])> = Vec::with_capacity(v.len());
        let mut dropped = 0;
        let mut i = 0;
        while i < v.len() {
            let mut j = i + 1;
            while j < v.len() && (v[j].0, v[j].1) == (v[i].0, v[i].1) {
                j += 1;
            }
            if j == i + 1 {
                out.push(v[i]);
            } else {
                dropped += j - i;
            }
            i = j;
        }
        *v = out;
        dropped
    };
    let dropped = dedup_e(&mut ecomp);
    hcomp.sort_by_key(|t| (t.0, t.1));
    let mut hout: Vec<(u8, u32, f32)> = Vec::with_capacity(hcomp.len());
    let mut hd = 0;
    let mut i = 0;
    while i < hcomp.len() {
        let mut j = i + 1;
        while j < hcomp.len() && (hcomp[j].0, hcomp[j].1) == (hcomp[i].0, hcomp[i].1) {
            j += 1;
        }
        if j == i + 1 {
            hout.push(hcomp[i]);
        } else {
            hd += j - i;
        }
        i = j;
    }
    fix.open.sort();
    fix.open.dedup();
    fix.close.sort();
    fix.close.dedup();
    fix.stats.opened = fix.open.len();
    fix.stats.closed = fix.close.len();
    fix.stats.dropped_shared = dropped + hd;
    fix.e = ecomp;
    fix.h = hout;
    fix
}

/// Largest eigenvalue of the (lossless, no-CPML) update operator times dt², by power iteration: the
/// leapfrog scheme is stable while this stays below 4. Uses the CPU kernels, H corrections included.
pub fn stability(p: &crate::model::Problem, iters: usize) -> f64 {
    use crate::fdtd::{hcorr_update, update_e, update_h, CpuState};
    let mut s = CpuState::new(p);
    // deterministic pseudo-random start on the live edges
    let mut seed = 0x2545_f491_4f6c_dd1du64;
    for c in 0..3 {
        for id in 0..p.n() {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            s.e[c][id] = if p.cb[c][id] != 0.0 { ((seed >> 40) as f32 / (1u64 << 24) as f32) - 0.5 } else { 0.0 };
        }
    }
    let norm = |e: &[Vec<f32>; 3]| -> f64 { e.iter().map(|v| v.iter().map(|&x| (x as f64) * (x as f64)).sum::<f64>()).sum::<f64>().sqrt() };
    let mut n0 = norm(&s.e);
    let mut mu = 0.0;
    for it in 0..iters {
        for c in 0..3 {
            s.e[c].iter_mut().for_each(|x| *x /= n0 as f32);
            s.h[c].iter_mut().for_each(|x| *x = 0.0);
        }
        update_h(p, &mut s);
        hcorr_update(p, &mut s);
        for c in 0..3 {
            s.e[c].iter_mut().for_each(|x| *x = 0.0);
        }
        update_e(p, &mut s);
        for c in 0..3 {
            s.e[c].iter_mut().for_each(|x| *x = -*x);
        }
        n0 = norm(&s.e);
        mu = n0;
        if n0 == 0.0 {
            break;
        }
        if std::env::var("FS_DEBUG").is_ok() && it % 25 == 0 {
            eprintln!("  power iteration {it}: dt^2 lambda {mu:.4}");
        }
    }
    mu
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn half_plane_line_and_flux_agree_with_brute_force() {
        // a 45 degree edge; compare the flux helpers with a fine midpoint sum of the analytic field
        let s = 0.5f64.sqrt();
        let ln = Line { x0: 0.03, y0: -0.02, nx: s, ny: s, z0: 0.0 };
        let field = |fld: Field, p: [f64; 3]| -> [f64; 3] {
            // E = -grad phi, H = -grad psi by central differences
            let hh = 1e-7;
            let f = |q: [f64; 3]| {
                let (u, w) = ln.uw(q);
                match fld {
                    Field::E => phi(u, w),
                    Field::H => psi(u, w, if w >= 0.0 { 1.0 } else { -1.0 }),
                }
            };
            let mut g = [0.0; 3];
            for a in 0..3 {
                let (mut a1, mut a0) = (p, p);
                a1[a] += hh;
                a0[a] -= hh;
                g[a] = -(f(a1) - f(a0)) / (2.0 * hh);
            }
            g
        };
        for fld in [Field::E, Field::H] {
            for (na, c, rb, rc) in [(2usize, 0.02, (-0.1, 0.1), (-0.08, 0.12)), (0usize, 0.01, (-0.1, 0.1), (0.01, 0.09)), (1usize, -0.03, (0.01, 0.07), (-0.1, 0.1)), (0usize, 0.04, (-0.05, 0.05), (-0.04, -0.01))] {
                let got = flux(fld, &ln, na, c, rb, rc);
                let (b, cc) = ((na + 1) % 3, (na + 2) % 3);
                let n = 400;
                let mut sum = 0.0;
                for ib in 0..n {
                    for ic in 0..n {
                        let mut p = [0.0; 3];
                        p[na] = c;
                        p[b] = rb.0 + (ib as f64 + 0.5) / n as f64 * (rb.1 - rb.0);
                        p[cc] = rc.0 + (ic as f64 + 0.5) / n as f64 * (rc.1 - rc.0);
                        sum += field(fld, p)[na];
                    }
                }
                sum *= (rb.1 - rb.0) * (rc.1 - rc.0) / (n * n) as f64;
                assert!((got - sum).abs() < 2e-3 * sum.abs().max(1e-3), "{fld:?} na {na}: {got} vs {sum}");
            }
            // line integral along z through the metal and along x above it
            for (a, p, len) in [(2usize, [-0.05, -0.05, -0.03], 0.06), (0usize, [-0.1, 0.0, 0.02], 0.2)] {
                let got = line_integral(fld, &ln, a, p, len);
                let n = 20000;
                let mut sum = 0.0;
                for t in 0..n {
                    let mut q = p;
                    q[a] += (t as f64 + 0.5) / n as f64 * len;
                    sum += field(fld, q)[a];
                }
                sum *= len / n as f64;
                assert!((got - sum).abs() < 2e-3 * sum.abs().max(1e-3), "{fld:?} line {a}: {got} vs {sum}");
            }
        }
    }

    #[test]
    fn line_fit_finds_a_straight_edge() {
        // half-plane copper x < 0.0123 + 0.4 y on a 0.05 mm grid
        let l: Vec<f64> = (0..21).map(|i| -0.5 + 0.05 * i as f64).collect();
        let s = 8;
        let xf = sub_centres(&l, s);
        let yf = sub_centres(&l, s);
        let wf = xf.len();
        let mut cov = vec![false; wf * yf.len()];
        for (jf, &y) in yf.iter().enumerate() {
            for (ifn, &x) in xf.iter().enumerate() {
                cov[jf * wf + ifn] = x < 0.0123 + 0.4 * y;
            }
        }
        let n = l.len();
        let pl = Plane { k: 0, s, wf, xf, yf, cov, fx: vec![1.0; n * n], fy: vec![1.0; n * n], ff: vec![1.0; n * n], node_in: vec![false; n * n] };
        let (ln, err) = fit(&pl, &l, &l, 0.02, 0.0, 0.0625, 0.0).unwrap();
        let nn = (1.0f64 + 0.16).sqrt();
        assert!((ln.nx - 1.0 / nn).abs() < 0.03 && (ln.ny + 0.4 / nn).abs() < 0.03, "normal {} {}", ln.nx, ln.ny);
        // distance of the true line from the fitted point
        let dist = (ln.x0 - 0.0123 - 0.4 * ln.y0) / nn;
        assert!(dist.abs() < 0.003, "offset {dist}");
        assert!(err < 0.03, "err {err}");
    }
}

#[cfg(test)]
mod probe_tests {
    use super::*;
    /// Print the factors across a strip edge along x (metal at y < ye) for uniform cells h, dz.
    #[test]
    #[ignore]
    fn print_factors_across_an_edge() {
        let (h, dz) = (1.0, 1.0);
        for inset in [0.05, 0.25, 0.3333, 0.5, 0.6667, 0.75, 0.95] {
            let ye = inset * h; // inner metal line at y = 0
            let ln = Line { x0: 0.0, y0: ye, nx: 0.0, ny: 1.0, z0: 0.0 };
            println!("inset {inset}");
            for m in -2i32..=2 {
                let y0 = m as f64 * h;
                // Ey(m): y0..y0+h on the sheet; dual face y = y0 + h/2, x in [-.5, .5], z in [-dz/2, dz/2]
                let e = line_integral(Field::E, &ln, 1, [0.0, y0, 0.0], h);
                let fl = flux(Field::E, &ln, 1, y0 + 0.5 * h, (-0.5 * dz, 0.5 * dz), (-0.5, 0.5));
                let ke = if e.abs() > 1e-12 { fl * h / (dz * 1.0 * e) } else { f64::NAN };
                // Hz(m): face y0..y0+h, x 0..1 on the sheet; dual edge along z at y0 + h/2
                let fh = flux(Field::H, &ln, 2, 0.0, (0.0, 1.0), (y0, y0 + h));
                let hl = line_integral(Field::H, &ln, 2, [0.5, y0 + 0.5 * h, -0.5 * dz], dz);
                let kh = fh * dz / (h * hl);
                // Ez(m) above at node y0: z 0..dz, dual face z = dz/2, y in [y0 - h/2, y0 + h/2]
                let ez = line_integral(Field::E, &ln, 2, [0.0, y0, 0.0], dz);
                let fz = flux(Field::E, &ln, 2, 0.5 * dz, (-0.5, 0.5), (y0 - 0.5 * h, y0 + 0.5 * h));
                let kz = fz * dz / (h * ez);
                // Hy(m) above: face y = y0, x 0..1, z 0..dz; dual edge along y at z = dz/2 from y0 - h/2
                let fy = flux(Field::H, &ln, 1, y0, (0.0, dz), (0.0, 1.0));
                let hy = line_integral(Field::H, &ln, 1, [0.5, y0 - 0.5 * h, 0.5 * dz], h);
                let ky = fy * h / (dz * hy);
                let fo = ((y0 + h - ye).max(0.0).min(h)) / h;
                println!("  row {m:+}: Ey kappa' {ke:7.3} (DM {:6.3})  Hz kappa {kh:7.3} (DM {fo:5.3})  Ez kappa' {kz:7.3}  Hy kappa {ky:7.3}", if fo > 0.0 { 1.0 / fo } else { f64::NAN });
            }
        }
    }
}

/// Diagnostic: run port `port` for `steps` on the CPU, then list the field components holding the most
/// energy (E^2 and eta^2 H^2), with their position and whether the edge code corrected them.
pub fn diag(p: &crate::model::Problem, port: usize, steps: usize, top: usize) {
    use crate::fdtd::*;
    let mut s = CpuState::new(p);
    let exc = gauss_excite(0.0, 8e9, p.dt);
    for n in 0..steps {
        update_h(p, &mut s);
        hcorr_update(p, &mut s);
        pml_h(p, &mut s);
        disp_save(p, &mut s);
        update_e(p, &mut s);
        pml_e(p, &mut s);
        if n < exc.len() {
            for &(c, id, g) in &p.ports[port].src {
                s.e[c as usize][id as usize] += g * exc[n];
            }
        }
        disp_update(p, &mut s);
        cond_update(p, &mut s);
        if n % 2000 == 0 {
            eprintln!("step {n}: energy {:.3e}", energy(&s));
        }
    }
    let eta2 = (crate::model::ETA0 * crate::model::ETA0) as f32;
    let mut all: Vec<(f32, u8, usize)> = Vec::new();
    for c in 0..3 {
        for id in 0..p.n() {
            let e = s.e[c][id] * s.e[c][id];
            let h = s.h[c][id] * s.h[c][id] * eta2;
            if e > 0.0 {
                all.push((e, c as u8, id));
            }
            if h > 0.0 {
                all.push((h, 3 + c as u8, id));
            }
        }
    }
    all.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    let tot: f64 = all.iter().map(|x| x.0 as f64).sum();
    let hc: Vec<u32> = p.hcorr.as_ref().map(|h| h.ids.clone()).unwrap_or_default();
    let (nx, ny) = (p.nx, p.ny);
    println!("total {tot:.3e}; top {top}:");
    for &(v, c, id) in all.iter().take(top) {
        let (i, j, k) = (id % nx, (id / nx) % ny, id / (nx * ny));
        let g = &p.grid;
        let corr = if c >= 3 { hc.binary_search(&(id as u32 | ((c as u32 - 3) << 30))).is_ok() } else { false };
        let cb = if c < 3 { p.cb[c as usize][id] } else { 0.0 };
        println!(
            "  {:5.2} %  {}{}  ({i}, {j}, {k}) = ({:.3}, {:.3}, {:.3}) mm  {}{}",
            100.0 * v as f64 / tot,
            if c < 3 { "E" } else { "H" },
            ["x", "y", "z"][c as usize % 3],
            g.x[i] * 1e3,
            g.y[j] * 1e3,
            g.z[k] * 1e3,
            if corr { "H corrected " } else { "" },
            if c < 3 { format!("cb {cb:.3e}") } else { String::new() }
        );
    }
}
