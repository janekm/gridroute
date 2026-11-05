//! Layout importer: the copper layers (Gerber RS-274X with X2 attributes), board outline (Edge_Cuts profile),
//! drills (Excellon) and stackup (Gerber job file) of a KiCad fabrication output directory, converted to the
//! shared geometry of `crate::geom` in board millimetres (Gerber orientation: y up).
//!
//! Supported Gerber: %FS (leading or trailing zero omission, absolute), %MO MM/IN (and G70/G71), G01/G02/G03
//! with G74/G75 arcs (linearised to a chord error of `CHORD_TOL`), D01/D02/D03, Dnn, G36/G37 regions (each
//! contour its own polygon), %LPD/%LPC, %LM/%LR/%LS, %SR step and repeat, standard apertures C/R/O/P, aperture
//! macros (primitives 0, 1, 4, 5, 6, 7, 20, 21, variables and arithmetic), X2 object attributes .N and .P.
//! Not supported (rejected with an error): block apertures (%AB), incremental coordinates, negative image
//! polarity (%IPNEG), non-zero %OF offsets.

use crate::geom::{BBox, Hole, Prim, Shape, P2};
use std::collections::HashMap;
use std::f64::consts::{PI, TAU};
use std::path::Path;

/// Maximum chord error when arcs and circles are turned into polylines (mm).
const CHORD_TOL: f64 = 0.002;
/// Endpoint matching tolerance when chaining outline segments into loops (mm).
const OUTLINE_TOL: f64 = 1e-4;
const MM_PER_IN: f64 = 25.4;

// ------------------------------------------------------------------------------------------------------------
// Public API
// ------------------------------------------------------------------------------------------------------------

/// One copper layer of the board.
#[derive(Clone, Debug)]
pub struct CopperLayer {
    /// KiCad layer name ("F.Cu", "In1.Cu", ..., "B.Cu").
    pub name: String,
    /// 1-based Gerber layer number (`%TF.FileFunction,Copper,Lk,...%`).
    pub index: usize,
    /// Objects in file order.
    pub prims: Vec<Prim>,
    /// Bounding box of all objects.
    pub bbox: BBox,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StackKind {
    Copper,
    Dielectric,
    SolderMask,
}

/// One physical layer of the stackup, top to bottom.
#[derive(Clone, Debug)]
pub struct StackLayer {
    pub name: String,
    pub kind: StackKind,
    pub thickness_mm: f64,
    pub eps_r: f64,
    pub loss_tan: f64,
}

/// A whole board.
#[derive(Clone, Debug)]
pub struct Board {
    /// Copper layers top to bottom.
    pub layers: Vec<CopperLayer>,
    /// Closed outline loops, outer loop first (implicitly closed, no repeated end point).
    pub outline: Vec<Vec<P2>>,
    pub holes: Vec<Hole>,
    /// Physical stackup top to bottom (solder masks, copper, dielectrics).
    pub stack: Vec<StackLayer>,
}

/// Parse one RS-274X file into its objects in file order.
pub fn parse_gerber(text: &str) -> Result<Vec<Prim>, String> {
    Ok(interpret(text)?.prims)
}

/// Parse a board profile (Edge_Cuts) file: all draws and arcs chained into closed loops, outer loop first.
pub fn parse_outline(text: &str) -> Result<Vec<Vec<P2>>, String> {
    let it = interpret(text)?;
    let mut paths = it.paths;
    // Region contours in a profile file are closed loops too.
    for p in it.prims {
        if let Shape::Polygon(pts) = p.shape {
            let mut pts = pts;
            pts.push(pts[0]);
            paths.push(pts);
        }
    }
    Ok(chain_loops(paths, OUTLINE_TOL))
}

// ------------------------------------------------------------------------------------------------------------
// Geometry helpers
// ------------------------------------------------------------------------------------------------------------

fn rot(p: P2, deg: f64) -> P2 {
    if deg == 0.0 {
        return p;
    }
    let (s, c) = (deg * PI / 180.0).sin_cos();
    P2::new(p.x * c - p.y * s, p.x * s + p.y * c)
}

fn add(a: P2, b: P2) -> P2 {
    P2::new(a.x + b.x, a.y + b.y)
}

fn poly_area(pts: &[P2]) -> f64 {
    let n = pts.len();
    let mut a = 0.0;
    for i in 0..n {
        let (p, q) = (pts[i], pts[(i + 1) % n]);
        a += p.x * q.y - q.x * p.y;
    }
    0.5 * a
}

/// Number of segments for an arc of radius `r` and sweep `sweep` (radians) with chord error <= CHORD_TOL.
fn arc_segments(r: f64, sweep: f64) -> usize {
    let sweep = sweep.abs();
    let max_step = if r > CHORD_TOL { 2.0 * (1.0 - CHORD_TOL / r).acos() } else { PI / 2.0 };
    let n = (sweep / max_step).ceil().max((sweep / (PI / 2.0)).ceil());
    (n as usize).max(1)
}

/// Points on the circle `c`, `r` from angle `a0` to `a1` (radians, either direction), both ends included.
fn arc_pts(c: P2, r: f64, a0: f64, a1: f64) -> Vec<P2> {
    let n = arc_segments(r, a1 - a0);
    (0..=n)
        .map(|k| {
            let a = a0 + (a1 - a0) * k as f64 / n as f64;
            P2::new(c.x + r * a.cos(), c.y + r * a.sin())
        })
        .collect()
}

/// A full circle as a polygon (no repeated point).
fn circle_poly(c: P2, r: f64) -> Vec<P2> {
    let mut v = arc_pts(c, r, 0.0, TAU);
    v.pop();
    v
}

/// Points of the circular arc from `s` to `e` about `c` (clockwise if `cw`), excluding `s`, ending exactly at
/// `e`. `s == e` is a full circle. The radius is interpolated if start and end radii differ slightly.
fn arc_through(s: P2, e: P2, c: P2, cw: bool) -> Vec<P2> {
    let (r0, r1) = (s.dist(c), e.dist(c));
    if r0 < 1e-9 || r1 < 1e-9 {
        return vec![e];
    }
    let a0 = (s.y - c.y).atan2(s.x - c.x);
    let a1 = (e.y - c.y).atan2(e.x - c.x);
    let mut sweep = if cw { a0 - a1 } else { a1 - a0 };
    while sweep <= 0.0 {
        sweep += TAU;
    }
    while sweep > TAU {
        sweep -= TAU;
    }
    if s.dist(e) < 1e-7 {
        sweep = TAU;
    }
    let n = arc_segments(r0.max(r1), sweep);
    let sgn = if cw { -1.0 } else { 1.0 };
    let mut out = Vec::with_capacity(n);
    for k in 1..n {
        let t = k as f64 / n as f64;
        let a = a0 + sgn * sweep * t;
        let r = r0 + (r1 - r0) * t;
        out.push(P2::new(c.x + r * a.cos(), c.y + r * a.sin()));
    }
    out.push(e);
    out
}

/// Convex hull (Andrew's monotone chain), counter-clockwise, no repeated point.
fn convex_hull(mut pts: Vec<P2>) -> Vec<P2> {
    pts.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap().then(a.y.partial_cmp(&b.y).unwrap()));
    pts.dedup_by(|a, b| (a.x - b.x).abs() < 1e-12 && (a.y - b.y).abs() < 1e-12);
    if pts.len() < 3 {
        return pts;
    }
    let cross = |o: P2, a: P2, b: P2| (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
    let mut h: Vec<P2> = Vec::with_capacity(2 * pts.len());
    for &p in pts.iter() {
        while h.len() >= 2 && cross(h[h.len() - 2], h[h.len() - 1], p) <= 0.0 {
            h.pop();
        }
        h.push(p);
    }
    let lo = h.len() + 1;
    for &p in pts.iter().rev().skip(1) {
        while h.len() >= lo && cross(h[h.len() - 2], h[h.len() - 1], p) <= 0.0 {
            h.pop();
        }
        h.push(p);
    }
    h.pop();
    h
}

/// Apply a point map to a shape; `rs` scales radii.
fn map_shape(s: &Shape, f: &dyn Fn(P2) -> P2, rs: f64) -> Shape {
    match s {
        Shape::Circle { c, r } => Shape::Circle { c: f(*c), r: r * rs },
        Shape::Capsule { a, b, r } => Shape::Capsule { a: f(*a), b: f(*b), r: r * rs },
        Shape::Polygon(p) => Shape::Polygon(p.iter().map(|&q| f(q)).collect()),
    }
}

/// Points describing the outline of a shape (circles sampled), for convex hulls.
fn shape_pts(s: &Shape) -> Vec<P2> {
    match s {
        Shape::Circle { c, r } => circle_poly(*c, *r),
        Shape::Capsule { a, b, r } => {
            let mut v = circle_poly(*a, *r);
            v.extend(circle_poly(*b, *r));
            v
        }
        Shape::Polygon(p) => p.clone(),
    }
}

/// True if the shape has a positive area.
fn non_degenerate(s: &Shape) -> bool {
    match s {
        Shape::Circle { r, .. } | Shape::Capsule { r, .. } => *r > 0.0,
        Shape::Polygon(p) => p.len() >= 3 && poly_area(p).abs() > 1e-14,
    }
}

/// Drop a repeated closing point.
fn open_ring(mut p: Vec<P2>) -> Vec<P2> {
    while p.len() > 1 && p[0].dist(p[p.len() - 1]) < 1e-9 {
        p.pop();
    }
    p
}

// ------------------------------------------------------------------------------------------------------------
// Aperture macro expressions
// ------------------------------------------------------------------------------------------------------------

struct ExprParser<'a> {
    s: &'a [u8],
    i: usize,
    vars: &'a HashMap<u32, f64>,
}

impl ExprParser<'_> {
    fn ws(&mut self) {
        while self.i < self.s.len() && self.s[self.i].is_ascii_whitespace() {
            self.i += 1;
        }
    }
    fn peek(&mut self) -> Option<u8> {
        self.ws();
        self.s.get(self.i).copied()
    }
    fn expr(&mut self) -> Result<f64, String> {
        let mut v = self.term()?;
        while let Some(c) = self.peek() {
            match c {
                b'+' => {
                    self.i += 1;
                    v += self.term()?;
                }
                b'-' => {
                    self.i += 1;
                    v -= self.term()?;
                }
                _ => break,
            }
        }
        Ok(v)
    }
    fn term(&mut self) -> Result<f64, String> {
        let mut v = self.factor()?;
        while let Some(c) = self.peek() {
            match c {
                b'x' | b'X' => {
                    self.i += 1;
                    v *= self.factor()?;
                }
                b'/' => {
                    self.i += 1;
                    v /= self.factor()?;
                }
                _ => break,
            }
        }
        Ok(v)
    }
    fn factor(&mut self) -> Result<f64, String> {
        match self.peek() {
            Some(b'-') => {
                self.i += 1;
                Ok(-self.factor()?)
            }
            Some(b'+') => {
                self.i += 1;
                self.factor()
            }
            Some(b'(') => {
                self.i += 1;
                let v = self.expr()?;
                if self.peek() != Some(b')') {
                    return Err("missing ')' in macro expression".into());
                }
                self.i += 1;
                Ok(v)
            }
            Some(b'$') => {
                self.i += 1;
                let st = self.i;
                while self.i < self.s.len() && self.s[self.i].is_ascii_digit() {
                    self.i += 1;
                }
                let n: u32 = std::str::from_utf8(&self.s[st..self.i])
                    .unwrap()
                    .parse()
                    .map_err(|_| "bad macro variable".to_string())?;
                Ok(self.vars.get(&n).copied().unwrap_or(0.0))
            }
            Some(c) if c.is_ascii_digit() || c == b'.' => {
                let st = self.i;
                while self.i < self.s.len() && (self.s[self.i].is_ascii_digit() || self.s[self.i] == b'.') {
                    self.i += 1;
                }
                let t = std::str::from_utf8(&self.s[st..self.i]).unwrap();
                t.parse::<f64>().map_err(|_| format!("bad number '{t}' in macro expression"))
            }
            Some(c) => Err(format!("unexpected '{}' in macro expression", c as char)),
            None => Err("unexpected end of macro expression".into()),
        }
    }
}

/// Evaluate an aperture macro arithmetic expression (`+ - x / ( )`, unary signs, `$n` variables).
pub fn eval_expr(s: &str, vars: &HashMap<u32, f64>) -> Result<f64, String> {
    let mut p = ExprParser { s: s.as_bytes(), i: 0, vars };
    let v = p.expr()?;
    if p.peek().is_some() {
        return Err(format!("trailing characters in macro expression '{s}'"));
    }
    Ok(v)
}

/// The thermal primitive's quadrant `q` (0..4) in local coordinates (centre at origin, before rotation).
fn thermal_quadrant(ro: f64, ri: f64, gap: f64, q: usize) -> Option<Vec<P2>> {
    let g = gap / 2.0;
    if g * 2f64.sqrt() >= ro {
        return None;
    }
    let ao = (g / ro).asin();
    let mut pts = arc_pts(P2::new(0.0, 0.0), ro, ao, PI / 2.0 - ao);
    if ri > g * 2f64.sqrt() {
        let ai = (g / ri).asin();
        pts.extend(arc_pts(P2::new(0.0, 0.0), ri, PI / 2.0 - ai, ai));
    } else {
        pts.push(P2::new(g, g));
    }
    Some(pts.into_iter().map(|p| rot(p, 90.0 * q as f64)).collect())
}

/// Annulus as a single even-odd polygon (outer ring, seam, inner ring reversed).
fn annulus(c: P2, ro: f64, ri: f64) -> Vec<P2> {
    let mut v = arc_pts(c, ro, 0.0, TAU);
    if ri > 0.0 {
        v.extend(arc_pts(c, ri, TAU, 0.0));
    }
    v
}

/// Evaluate a macro with the given parameters into local shapes (in file units) with their exposure.
fn eval_macro(body: &[String], params: &[f64]) -> Result<Vec<(Shape, bool)>, String> {
    let mut vars: HashMap<u32, f64> = HashMap::new();
    for (k, &v) in params.iter().enumerate() {
        vars.insert(k as u32 + 1, v);
    }
    let mut out = Vec::new();
    for blk in body {
        let b = blk.trim();
        if b.is_empty() || b.starts_with('0') {
            continue; // comment
        }
        if let Some(rest) = b.strip_prefix('$') {
            let (n, e) = rest.split_once('=').ok_or_else(|| format!("bad macro statement '{b}'"))?;
            let n: u32 = n.trim().parse().map_err(|_| format!("bad macro variable in '{b}'"))?;
            let v = eval_expr(e, &vars)?;
            vars.insert(n, v);
            continue;
        }
        let fields: Vec<&str> = b.split(',').collect();
        let code: u32 = fields[0].trim().parse().map_err(|_| format!("bad macro primitive '{b}'"))?;
        let a: Vec<f64> = fields[1..].iter().map(|f| eval_expr(f, &vars)).collect::<Result<_, _>>()?;
        let arg = |i: usize| a.get(i).copied().unwrap_or(0.0);
        let need = |n: usize| -> Result<(), String> {
            if a.len() < n {
                Err(format!("macro primitive {code} needs {n} parameters: '{b}'"))
            } else {
                Ok(())
            }
        };
        let on = |e: f64| e > 0.5;
        match code {
            1 => {
                need(4)?;
                let c = rot(P2::new(arg(2), arg(3)), arg(4));
                out.push((Shape::Circle { c, r: arg(1) / 2.0 }, on(arg(0))));
            }
            20 | 2 => {
                need(7)?;
                let (s, e, w, r) = (P2::new(arg(2), arg(3)), P2::new(arg(4), arg(5)), arg(1), arg(6));
                let l = s.dist(e);
                if l > 0.0 && w > 0.0 {
                    let (nx, ny) = (-(e.y - s.y) / l * w / 2.0, (e.x - s.x) / l * w / 2.0);
                    let pts = vec![
                        P2::new(s.x + nx, s.y + ny),
                        P2::new(e.x + nx, e.y + ny),
                        P2::new(e.x - nx, e.y - ny),
                        P2::new(s.x - nx, s.y - ny),
                    ];
                    out.push((Shape::Polygon(pts.into_iter().map(|p| rot(p, r)).collect()), on(arg(0))));
                }
            }
            21 => {
                need(6)?;
                let (w, h, cx, cy, r) = (arg(1) / 2.0, arg(2) / 2.0, arg(3), arg(4), arg(5));
                let pts = [(-w, -h), (w, -h), (w, h), (-w, h)];
                let pts = pts.iter().map(|&(x, y)| rot(P2::new(cx + x, cy + y), r)).collect();
                out.push((Shape::Polygon(pts), on(arg(0))));
            }
            22 => {
                // deprecated lower-left line: exposure, w, h, x, y, rotation
                need(6)?;
                let (w, h, x, y, r) = (arg(1), arg(2), arg(3), arg(4), arg(5));
                let pts = [(0.0, 0.0), (w, 0.0), (w, h), (0.0, h)];
                let pts = pts.iter().map(|&(dx, dy)| rot(P2::new(x + dx, y + dy), r)).collect();
                out.push((Shape::Polygon(pts), on(arg(0))));
            }
            4 => {
                need(2)?;
                let n = arg(1).round() as usize;
                need(2 + 2 * (n + 1) + 1)?;
                let r = arg(2 + 2 * (n + 1));
                let pts: Vec<P2> = (0..=n).map(|k| rot(P2::new(arg(2 + 2 * k), arg(3 + 2 * k)), r)).collect();
                out.push((Shape::Polygon(open_ring(pts)), on(arg(0))));
            }
            5 => {
                need(5)?;
                let n = (arg(1).round() as usize).max(3);
                let (c, rr, r) = (P2::new(arg(2), arg(3)), arg(4) / 2.0, arg(5));
                let pts = (0..n)
                    .map(|k| {
                        let t = TAU * k as f64 / n as f64;
                        rot(P2::new(c.x + rr * t.cos(), c.y + rr * t.sin()), r)
                    })
                    .collect();
                out.push((Shape::Polygon(pts), on(arg(0))));
            }
            6 => {
                // moire (deprecated): cx, cy, outer d, ring thickness, gap, max rings, cross thickness, cross length, rot
                need(8)?;
                let (c, mut d, th, gap, nmax, ct, cl, r) =
                    (P2::new(arg(0), arg(1)), arg(2), arg(3), arg(4), arg(5) as usize, arg(6), arg(7), arg(8));
                for _ in 0..nmax {
                    if d <= 0.0 {
                        break;
                    }
                    let ri = (d / 2.0 - th).max(0.0);
                    let pts = annulus(c, d / 2.0, ri).into_iter().map(|p| rot(p, r)).collect();
                    out.push((Shape::Polygon(pts), true));
                    d -= 2.0 * (th + gap);
                }
                for (w, h) in [(cl, ct), (ct, cl)] {
                    let pts = [(-w, -h), (w, -h), (w, h), (-w, h)];
                    let pts = pts.iter().map(|&(x, y)| rot(P2::new(c.x + x / 2.0, c.y + y / 2.0), r)).collect();
                    out.push((Shape::Polygon(pts), true));
                }
            }
            7 => {
                need(5)?;
                let (c, ro, ri, gap, r) = (P2::new(arg(0), arg(1)), arg(2) / 2.0, arg(3) / 2.0, arg(4), arg(5));
                for q in 0..4 {
                    if let Some(pts) = thermal_quadrant(ro, ri, gap, q) {
                        out.push((Shape::Polygon(pts.into_iter().map(|p| rot(add(p, c), r)).collect()), true));
                    }
                }
            }
            _ => return Err(format!("unsupported aperture macro primitive {code}")),
        }
    }
    Ok(out)
}

// ------------------------------------------------------------------------------------------------------------
// Gerber interpreter
// ------------------------------------------------------------------------------------------------------------

#[derive(Clone, Debug)]
enum Aperture {
    /// Standard circle, diameter in mm.
    Circle(f64),
    /// Anything else: local shapes in mm (origin = flash point) and their exposure.
    Shapes(Vec<(Shape, bool)>),
}

#[derive(Clone, Copy, PartialEq)]
enum Interp {
    Linear,
    Cw,
    Ccw,
}

/// Aperture transformation (%LM, %LS, %LR).
#[derive(Clone, Copy)]
struct Xform {
    mx: bool,
    my: bool,
    s: f64,
    r: f64,
}

impl Xform {
    fn apply(&self, p: P2) -> P2 {
        let x = if self.mx { -p.x } else { p.x };
        let y = if self.my { -p.y } else { p.y };
        rot(P2::new(x * self.s, y * self.s), self.r)
    }
}

/// Coordinate format from %FS.
#[derive(Clone, Copy)]
struct Fmt {
    int: usize,
    dec: usize,
    trailing_omitted: bool,
}

struct Gerber {
    prims: Vec<Prim>,
    /// Pen-down polylines (each continuous sequence of D01 draws), for outlines.
    paths: Vec<Vec<P2>>,
    cur_path: Vec<P2>,
    fmt: Option<Fmt>,
    /// mm per file unit.
    unit: f64,
    unit_set: bool,
    apertures: HashMap<u32, Aperture>,
    macros: HashMap<String, Vec<String>>,
    cur_ap: Option<u32>,
    pos: P2,
    interp: Interp,
    multi_quadrant: bool,
    /// Some(contour) while in region mode.
    region: Option<Vec<P2>>,
    dark: bool,
    xf: Xform,
    net: Option<String>,
    pad: Option<(String, String)>,
    last_op: Option<u32>,
    /// Open step-and-repeat block: (first prim index, nx, ny, dx, dy).
    sr: Option<(usize, usize, usize, f64, f64)>,
}

/// Decode X2 `\uXXXX` escapes.
fn x2_unescape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut it = s.chars().peekable();
    while let Some(c) = it.next() {
        if c == '\\' && it.peek() == Some(&'u') {
            let hex: String = it.clone().skip(1).take(4).collect();
            if hex.len() == 4 {
                if let Some(ch) = u32::from_str_radix(&hex, 16).ok().and_then(char::from_u32) {
                    out.push(ch);
                    for _ in 0..5 {
                        it.next();
                    }
                    continue;
                }
            }
        }
        out.push(c);
    }
    out
}

/// Split a Gerber data word into (letter, number) tokens: "G01X100Y-20D01" -> [(G,01),(X,100),(Y,-20),(D,01)].
fn tokens(w: &str) -> Result<Vec<(char, &str)>, String> {
    let b = w.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        let l = b[i] as char;
        if !l.is_ascii_alphabetic() {
            return Err(format!("unexpected '{l}' in word '{w}'"));
        }
        i += 1;
        let st = i;
        while i < b.len() && (b[i].is_ascii_digit() || b[i] == b'+' || b[i] == b'-' || b[i] == b'.') {
            i += 1;
        }
        out.push((l.to_ascii_uppercase(), &w[st..i]));
    }
    Ok(out)
}

impl Gerber {
    fn new() -> Gerber {
        Gerber {
            prims: Vec::new(),
            paths: Vec::new(),
            cur_path: Vec::new(),
            fmt: None,
            unit: 1.0,
            unit_set: false,
            apertures: HashMap::new(),
            macros: HashMap::new(),
            cur_ap: None,
            pos: P2::new(0.0, 0.0),
            interp: Interp::Linear,
            multi_quadrant: true,
            region: None,
            dark: true,
            xf: Xform { mx: false, my: false, s: 1.0, r: 0.0 },
            net: None,
            pad: None,
            last_op: None,
            sr: None,
        }
    }

    fn emit(&mut self, shape: Shape, dark: bool) {
        if non_degenerate(&shape) {
            self.prims.push(Prim { shape, dark, net: self.net.clone(), pad: self.pad.clone() });
        }
    }

    fn flush_path(&mut self) {
        if self.cur_path.len() >= 2 {
            self.paths.push(std::mem::take(&mut self.cur_path));
        } else {
            self.cur_path.clear();
        }
    }

    fn flush_contour(&mut self) {
        if let Some(c) = self.region.as_mut() {
            let pts = open_ring(std::mem::take(c));
            self.emit(Shape::Polygon(pts), self.dark);
        }
    }

    fn coord(&self, s: &str) -> Result<f64, String> {
        if s.contains('.') {
            return s.parse::<f64>().map(|v| v * self.unit).map_err(|_| format!("bad coordinate '{s}'"));
        }
        let f = self.fmt.ok_or("coordinate data before %FS")?;
        let (neg, digits) = match s.as_bytes().first() {
            Some(b'-') => (true, &s[1..]),
            Some(b'+') => (false, &s[1..]),
            _ => (false, s),
        };
        if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) {
            return Err(format!("bad coordinate '{s}'"));
        }
        let mut v: f64 = if f.trailing_omitted {
            let mut d = digits.to_string();
            while d.len() < f.int + f.dec {
                d.push('0');
            }
            d.parse::<f64>().unwrap() / 10f64.powi(f.dec as i32)
        } else {
            digits.parse::<f64>().unwrap() / 10f64.powi(f.dec as i32)
        };
        if neg {
            v = -v;
        }
        Ok(v * self.unit)
    }

    /// One `*`-terminated data word. Returns true at end of file (M02).
    fn word(&mut self, raw: &str) -> Result<bool, String> {
        let w: String = raw.chars().filter(|c| !c.is_whitespace()).collect();
        if w.is_empty() {
            return Ok(false);
        }
        let wb = w.as_bytes();
        if w.starts_with("G04") || (w.starts_with("G4") && !wb.get(2).is_some_and(|c| c.is_ascii_digit())) {
            return Ok(false); // comment
        }
        let toks = tokens(&w)?;
        let (mut x, mut y, mut i, mut j, mut d) = (None, None, None, None, None);
        for (l, n) in toks {
            match l {
                'G' => {
                    let g: u32 = n.parse().map_err(|_| format!("bad G code in '{w}'"))?;
                    match g {
                        1 => self.interp = Interp::Linear,
                        2 => self.interp = Interp::Cw,
                        3 => self.interp = Interp::Ccw,
                        36 => {
                            self.flush_path();
                            self.region = Some(Vec::new());
                        }
                        37 => {
                            self.flush_contour();
                            self.region = None;
                        }
                        74 => self.multi_quadrant = false,
                        75 => self.multi_quadrant = true,
                        54 | 55 | 90 => {}
                        70 => self.unit = MM_PER_IN,
                        71 => self.unit = 1.0,
                        91 => return Err("incremental coordinates (G91) are not supported".into()),
                        _ => return Err(format!("unsupported G code G{g:02}")),
                    }
                }
                'X' => x = Some(self.coord(n)?),
                'Y' => y = Some(self.coord(n)?),
                'I' => i = Some(self.coord(n)?),
                'J' => j = Some(self.coord(n)?),
                'D' => d = Some(n.parse::<u32>().map_err(|_| format!("bad D code in '{w}'"))?),
                'M' => {
                    let m: u32 = n.parse().map_err(|_| format!("bad M code in '{w}'"))?;
                    if m == 2 || m == 0 || m == 1 {
                        return Ok(true);
                    }
                    return Err(format!("unsupported M code M{m:02}"));
                }
                'N' => {} // deprecated sequence number
                _ => return Err(format!("unsupported word '{w}'")),
            }
        }
        let has_coord = x.is_some() || y.is_some() || i.is_some() || j.is_some();
        let op = match d {
            Some(n) if n >= 10 => {
                if !self.apertures.contains_key(&n) {
                    return Err(format!("aperture D{n} selected but not defined"));
                }
                self.flush_path();
                self.cur_ap = Some(n);
                None
            }
            Some(n @ 1..=3) => Some(n),
            Some(n) => return Err(format!("unsupported D code D{n:02}")),
            None if has_coord => self.last_op, // deprecated modal operation
            None => None,
        };
        if let Some(op) = op {
            self.last_op = Some(op);
            let to = P2::new(x.unwrap_or(self.pos.x), y.unwrap_or(self.pos.y));
            let off = P2::new(i.unwrap_or(0.0), j.unwrap_or(0.0));
            match op {
                1 => self.draw(to, off)?,
                2 => {
                    if self.region.as_ref().is_some_and(|c| !c.is_empty()) {
                        self.flush_contour();
                    }
                    self.flush_path();
                    self.pos = to;
                }
                _ => {
                    if self.region.is_some() {
                        return Err("D03 flash inside a region".into());
                    }
                    self.flush_path();
                    self.pos = to;
                    self.flash(to)?;
                }
            }
        }
        Ok(false)
    }

    /// The polyline of an interpolation from the current point to `to` (first point = current point).
    fn interpolate(&self, to: P2, off: P2) -> Vec<P2> {
        let s = self.pos;
        let mut pts = vec![s];
        match self.interp {
            Interp::Linear => pts.push(to),
            Interp::Cw | Interp::Ccw => {
                let cw = self.interp == Interp::Cw;
                let c = if self.multi_quadrant {
                    add(s, off)
                } else {
                    // single quadrant: unsigned offsets, pick the centre giving a <= 90 degree arc
                    let mut best: Option<(f64, P2)> = None;
                    for (sx, sy) in [(1.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)] {
                        let c = P2::new(s.x + sx * off.x.abs(), s.y + sy * off.y.abs());
                        let a0 = (s.y - c.y).atan2(s.x - c.x);
                        let a1 = (to.y - c.y).atan2(to.x - c.x);
                        let mut sw = if cw { a0 - a1 } else { a1 - a0 };
                        while sw < 0.0 {
                            sw += TAU;
                        }
                        while sw >= TAU {
                            sw -= TAU;
                        }
                        if sw <= PI / 2.0 + 1e-6 {
                            let err = (s.dist(c) - to.dist(c)).abs();
                            if best.is_none_or(|(e, _)| err < e) {
                                best = Some((err, c));
                            }
                        }
                    }
                    best.map(|b| b.1).unwrap_or(add(s, off))
                };
                if self.multi_quadrant {
                    pts.extend(arc_through(s, to, c, cw));
                } else if s.dist(to) < 1e-7 {
                    pts.push(to); // single quadrant mode: zero-length arc
                } else {
                    pts.extend(arc_through(s, to, c, cw));
                }
            }
        }
        pts
    }

    fn draw(&mut self, to: P2, off: P2) -> Result<(), String> {
        let pts = self.interpolate(to, off);
        self.pos = to;
        if let Some(c) = self.region.as_mut() {
            if c.is_empty() {
                c.push(pts[0]);
            }
            c.extend_from_slice(&pts[1..]);
            return Ok(());
        }
        // pen-down path for outlines
        if self.cur_path.last().is_some_and(|p| p.dist(pts[0]) < 1e-9) {
            self.cur_path.extend_from_slice(&pts[1..]);
        } else {
            self.flush_path();
            self.cur_path = pts.clone();
        }
        let ap = self.cur_ap.ok_or("D01 draw with no aperture selected")?;
        match self.apertures[&ap].clone() {
            Aperture::Circle(d) => {
                let r = d / 2.0 * self.xf.s;
                for w in pts.windows(2) {
                    self.emit(Shape::Capsule { a: w[0], b: w[1], r }, self.dark);
                }
            }
            Aperture::Shapes(shapes) => {
                let xf = self.xf;
                let local: Vec<P2> = shapes
                    .iter()
                    .filter(|s| s.1)
                    .flat_map(|(s, _)| shape_pts(&map_shape(s, &|p| xf.apply(p), xf.s)))
                    .collect();
                for w in pts.windows(2) {
                    let mut hp: Vec<P2> = local.iter().map(|&p| add(p, w[0])).collect();
                    hp.extend(local.iter().map(|&p| add(p, w[1])));
                    self.emit(Shape::Polygon(convex_hull(hp)), self.dark);
                }
            }
        }
        Ok(())
    }

    fn flash(&mut self, at: P2) -> Result<(), String> {
        let ap = self.cur_ap.ok_or("D03 flash with no aperture selected")?;
        let xf = self.xf;
        match self.apertures[&ap].clone() {
            Aperture::Circle(d) => {
                let c = add(at, xf.apply(P2::new(0.0, 0.0)));
                self.emit(Shape::Circle { c, r: d / 2.0 * xf.s }, self.dark);
            }
            Aperture::Shapes(shapes) => {
                for (s, exposed) in shapes {
                    let sh = map_shape(&s, &|p| add(at, xf.apply(p)), xf.s);
                    self.emit(sh, self.dark && exposed);
                }
            }
        }
        Ok(())
    }

    fn close_sr(&mut self) {
        if let Some((start, nx, ny, dx, dy)) = self.sr.take() {
            let block: Vec<Prim> = self.prims[start..].to_vec();
            for iy in 0..ny {
                for ix in 0..nx {
                    if ix == 0 && iy == 0 {
                        continue;
                    }
                    let o = P2::new(ix as f64 * dx, iy as f64 * dy);
                    for p in &block {
                        let mut q = p.clone();
                        q.shape = map_shape(&p.shape, &|v| add(v, o), 1.0);
                        self.prims.push(q);
                    }
                }
            }
        }
    }

    /// The content of one `%...%` extended command.
    fn extended(&mut self, content: &str) -> Result<(), String> {
        let blocks: Vec<String> = content
            .split('*')
            .map(|b| b.replace(['\r', '\n'], "").trim().to_string())
            .filter(|b| !b.is_empty())
            .collect();
        if blocks.is_empty() {
            return Ok(());
        }
        if let Some(name) = blocks[0].strip_prefix("AM") {
            self.macros.insert(name.to_string(), blocks[1..].to_vec());
            return Ok(());
        }
        for b in &blocks {
            self.param(b)?;
        }
        Ok(())
    }

    fn param(&mut self, b: &str) -> Result<(), String> {
        let code = b.get(..2).unwrap_or(b);
        let rest = b.get(2..).unwrap_or("");
        match code {
            "FS" => {
                let r = rest.as_bytes();
                let trailing_omitted = r.first() == Some(&b'T');
                if rest.contains('I') && !rest.starts_with("LA") && !rest.starts_with("TA") {
                    return Err(format!("incremental coordinate format not supported: %{b}%"));
                }
                let xi = rest.find('X').ok_or_else(|| format!("bad %FS: {b}"))?;
                let dig: Vec<usize> = rest[xi + 1..]
                    .bytes()
                    .take(2)
                    .map(|c| (c as char).to_digit(10).map(|d| d as usize))
                    .collect::<Option<_>>()
                    .ok_or_else(|| format!("bad %FS: {b}"))?;
                if dig.len() != 2 {
                    return Err(format!("bad %FS: {b}"));
                }
                self.fmt = Some(Fmt { int: dig[0], dec: dig[1], trailing_omitted });
            }
            "MO" => {
                self.unit = match rest {
                    "MM" => 1.0,
                    "IN" => MM_PER_IN,
                    _ => return Err(format!("bad %MO: {b}")),
                };
                self.unit_set = true;
            }
            "AD" => self.define_aperture(rest)?,
            "LP" => {
                self.dark = match rest {
                    "D" => true,
                    "C" => false,
                    _ => return Err(format!("bad %LP: {b}")),
                }
            }
            "LM" => {
                let (mx, my) = match rest {
                    "N" => (false, false),
                    "X" => (true, false),
                    "Y" => (false, true),
                    "XY" => (true, true),
                    _ => return Err(format!("bad %LM: {b}")),
                };
                self.xf.mx = mx;
                self.xf.my = my;
            }
            "LR" => self.xf.r = rest.parse().map_err(|_| format!("bad %LR: {b}"))?,
            "LS" => self.xf.s = rest.parse().map_err(|_| format!("bad %LS: {b}"))?,
            "SR" => {
                self.close_sr();
                if !rest.is_empty() {
                    let mut v = HashMap::new();
                    for (l, n) in tokens(rest)? {
                        v.insert(l, n.parse::<f64>().map_err(|_| format!("bad %SR: {b}"))?);
                    }
                    let nx = v.get(&'X').copied().unwrap_or(1.0).max(1.0) as usize;
                    let ny = v.get(&'Y').copied().unwrap_or(1.0).max(1.0) as usize;
                    let dx = v.get(&'I').copied().unwrap_or(0.0) * self.unit;
                    let dy = v.get(&'J').copied().unwrap_or(0.0) * self.unit;
                    if nx > 1 || ny > 1 {
                        self.sr = Some((self.prims.len(), nx, ny, dx, dy));
                    }
                }
            }
            "TO" => {
                let body = x2_unescape(rest);
                let mut f = body.split(',');
                let key = f.next().unwrap_or("");
                let vals: Vec<&str> = f.collect();
                match key {
                    ".N" => {
                        let n = vals.first().copied().unwrap_or("");
                        self.net = if n.is_empty() || n == "N/C" { None } else { Some(n.to_string()) };
                    }
                    ".P" => {
                        self.pad =
                            vals.first().map(|r| (r.to_string(), vals.get(1).copied().unwrap_or("").to_string()));
                    }
                    _ => {}
                }
            }
            "TD" => match rest {
                "" => {
                    self.net = None;
                    self.pad = None;
                }
                ".N" => self.net = None,
                ".P" => self.pad = None,
                _ => {}
            },
            "TF" | "TA" | "IN" | "LN" | "AS" => {}
            "IP" => {
                if rest == "NEG" {
                    return Err("negative image polarity (%IPNEG%) is not supported".into());
                }
            }
            "OF" => {
                let zero = tokens(rest)?.iter().all(|(_, n)| n.parse::<f64>().map(|v| v == 0.0).unwrap_or(false));
                if !zero {
                    return Err(format!("image offset %{b}% is not supported"));
                }
            }
            "MI" | "SF" | "IR" => {
                let trivial = matches!(b, "MIA0B0" | "SFA1B1" | "SFA1.0B1.0" | "IR0");
                if !trivial {
                    return Err(format!("deprecated image transformation %{b}% is not supported"));
                }
            }
            "AB" => return Err("block apertures (%AB%) are not supported".into()),
            _ => return Err(format!("unsupported extended command %{b}%")),
        }
        Ok(())
    }

    fn define_aperture(&mut self, rest: &str) -> Result<(), String> {
        // rest = "D10C,0.5X0.2"
        let r = rest.strip_prefix('D').ok_or_else(|| format!("bad %AD: {rest}"))?;
        let nd = r.bytes().take_while(|c| c.is_ascii_digit()).count();
        let code: u32 = r[..nd].parse().map_err(|_| format!("bad aperture number in %AD{rest}"))?;
        let (name, params) = match r[nd..].split_once(',') {
            Some((n, p)) => (n, p),
            None => (&r[nd..], ""),
        };
        let p: Vec<f64> = if params.is_empty() {
            Vec::new()
        } else {
            params
                .split('X')
                .map(|v| v.trim().parse::<f64>())
                .collect::<Result<_, _>>()
                .map_err(|_| format!("bad aperture parameters in %AD{rest}"))?
        };
        let u = self.unit;
        let need = |n: usize| -> Result<(), String> {
            if p.len() < n {
                Err(format!("aperture %AD{rest} needs {n} parameters"))
            } else {
                Ok(())
            }
        };
        let o = P2::new(0.0, 0.0);
        let ap = match name {
            "C" => {
                need(1)?;
                Aperture::Circle(p[0] * u)
            }
            "R" => {
                need(2)?;
                let (w, h) = (p[0] * u / 2.0, p[1] * u / 2.0);
                let pts = vec![P2::new(-w, -h), P2::new(w, -h), P2::new(w, h), P2::new(-w, h)];
                Aperture::Shapes(vec![(Shape::Polygon(pts), true)])
            }
            "O" => {
                need(2)?;
                let (w, h) = (p[0] * u, p[1] * u);
                let pts = if (w - h).abs() < 1e-12 {
                    circle_poly(o, w / 2.0)
                } else if w > h {
                    let (c, r) = ((w - h) / 2.0, h / 2.0);
                    let mut v = arc_pts(P2::new(c, 0.0), r, -PI / 2.0, PI / 2.0);
                    v.extend(arc_pts(P2::new(-c, 0.0), r, PI / 2.0, 1.5 * PI));
                    v
                } else {
                    let (c, r) = ((h - w) / 2.0, w / 2.0);
                    let mut v = arc_pts(P2::new(0.0, c), r, 0.0, PI);
                    v.extend(arc_pts(P2::new(0.0, -c), r, PI, TAU));
                    v
                };
                Aperture::Shapes(vec![(Shape::Polygon(pts), true)])
            }
            "P" => {
                need(2)?;
                let (r, n, a0) = (p[0] * u / 2.0, (p[1].round() as usize).max(3), p.get(2).copied().unwrap_or(0.0));
                let pts = (0..n)
                    .map(|k| {
                        let t = (a0 + 360.0 * k as f64 / n as f64) * PI / 180.0;
                        P2::new(r * t.cos(), r * t.sin())
                    })
                    .collect();
                Aperture::Shapes(vec![(Shape::Polygon(pts), true)])
            }
            m => {
                let body = self.macros.get(m).ok_or_else(|| format!("aperture D{code} uses undefined macro '{m}'"))?;
                let shapes = eval_macro(body, &p)?;
                Aperture::Shapes(
                    shapes.into_iter().map(|(s, e)| (map_shape(&s, &|q| P2::new(q.x * u, q.y * u), u), e)).collect(),
                )
            }
        };
        self.apertures.insert(code, ap);
        Ok(())
    }

    fn finish(&mut self) {
        if self.region.is_some() {
            self.flush_contour();
            self.region = None;
        }
        self.flush_path();
        self.close_sr();
    }
}

fn interpret(text: &str) -> Result<Gerber, String> {
    let mut st = Gerber::new();
    let b = text.as_bytes();
    let mut i = 0;
    let mut line = 1usize;
    let ctx = |e: String, line: usize| format!("Gerber line {line}: {e}");
    while i < b.len() {
        let c = b[i];
        if c == b'%' {
            let end = text[i + 1..].find('%').ok_or_else(|| ctx("unterminated extended command".into(), line))? + i + 1;
            st.extended(&text[i + 1..end]).map_err(|e| ctx(e, line))?;
            line += text[i..end].bytes().filter(|&c| c == b'\n').count();
            i = end + 1;
        } else if c.is_ascii_whitespace() {
            if c == b'\n' {
                line += 1;
            }
            i += 1;
        } else {
            let end = match text[i..].find('*') {
                Some(k) => i + k,
                None => {
                    // tolerate trailing garbage without terminator only if it is whitespace
                    return Err(ctx("unterminated data word".into(), line));
                }
            };
            let done = st.word(&text[i..end]).map_err(|e| ctx(e, line))?;
            line += text[i..end].bytes().filter(|&c| c == b'\n').count();
            if done {
                break;
            }
            i = end + 1;
        }
    }
    st.finish();
    Ok(st)
}

/// The `%TF.FileFunction,...%` value of a Gerber file (also the KiCad `G04 #@! TF.FileFunction,...` form).
fn gerber_file_function(text: &str) -> Option<String> {
    for l in text.lines().take(200) {
        let l = l.trim();
        for pre in ["%TF.FileFunction,", "G04 #@! TF.FileFunction,"] {
            if let Some(r) = l.strip_prefix(pre) {
                return Some(r.trim_end_matches('%').trim_end_matches('*').to_string());
            }
        }
    }
    None
}

// ------------------------------------------------------------------------------------------------------------
// Outline loops
// ------------------------------------------------------------------------------------------------------------

/// Chain polylines into closed loops by endpoint matching within `tol`; outer (largest) loop first.
fn chain_loops(paths: Vec<Vec<P2>>, tol: f64) -> Vec<Vec<P2>> {
    let paths: Vec<Vec<P2>> = paths.into_iter().filter(|p| p.len() >= 2).collect();
    let cell = (tol * 10.0).max(1e-9);
    let key = |p: P2| ((p.x / cell).floor() as i64, (p.y / cell).floor() as i64);
    let mut grid: HashMap<(i64, i64), Vec<(usize, bool)>> = HashMap::new();
    for (k, p) in paths.iter().enumerate() {
        grid.entry(key(p[0])).or_default().push((k, false));
        grid.entry(key(p[p.len() - 1])).or_default().push((k, true));
    }
    let mut used = vec![false; paths.len()];
    let find = |q: P2, used: &[bool]| -> Option<(usize, bool)> {
        let (kx, ky) = key(q);
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(v) = grid.get(&(kx + dx, ky + dy)) {
                    for &(k, is_end) in v {
                        if used[k] {
                            continue;
                        }
                        let p = if is_end { paths[k][paths[k].len() - 1] } else { paths[k][0] };
                        if p.dist(q) <= tol {
                            return Some((k, is_end));
                        }
                    }
                }
            }
        }
        None
    };
    let mut loops = Vec::new();
    for s in 0..paths.len() {
        if used[s] {
            continue;
        }
        used[s] = true;
        let mut chain = paths[s].clone();
        let mut flipped = false;
        let closed = loop {
            if chain.len() >= 3 && chain[0].dist(chain[chain.len() - 1]) <= tol {
                break true;
            }
            match find(chain[chain.len() - 1], &used) {
                Some((k, is_end)) => {
                    used[k] = true;
                    if is_end {
                        chain.extend(paths[k].iter().rev().skip(1));
                    } else {
                        chain.extend(paths[k].iter().skip(1));
                    }
                }
                None if !flipped => {
                    flipped = true;
                    chain.reverse();
                }
                None => break chain.len() >= 3 && chain[0].dist(chain[chain.len() - 1]) <= 0.01,
            }
        };
        if closed {
            let mut c = chain;
            c.pop();
            if c.len() >= 3 {
                loops.push(c);
            }
        } else {
            eprintln!(
                "gerber: outline: dropping open chain of {} points from ({:.4}, {:.4}) to ({:.4}, {:.4})",
                chain.len(),
                chain[0].x,
                chain[0].y,
                chain[chain.len() - 1].x,
                chain[chain.len() - 1].y
            );
        }
    }
    loops.sort_by(|a, b| poly_area(b).abs().partial_cmp(&poly_area(a).abs()).unwrap());
    loops
}

// ------------------------------------------------------------------------------------------------------------
// Excellon
// ------------------------------------------------------------------------------------------------------------

/// Parse an Excellon drill file. Plating and layer span come from the KiCad X2 comments
/// (`; #@! TF.FileFunction,Plated,1,6,PTH`, per-tool `; #@! TA.AperFunction,Plated,...`); without them holes
/// are plated (unless the header mentions NPTH) and span layers 1..=`n_copper`. Slots (G85 or routed
/// G00/M15/G01/M16) become holes with `end` set; a routed path of several segments gives one hole per segment.
pub fn parse_excellon(text: &str, n_copper: usize) -> Result<Vec<Hole>, String> {
    let mut metric = false;
    let mut lz_kept: Option<bool> = None; // Some(true): leading zeros kept (trailing suppressed)
    let mut fmt: Option<(usize, usize)> = None;
    let mut tools: HashMap<u32, (f64, Option<bool>)> = HashMap::new();
    let mut cur_tool: Option<u32> = None;
    let mut pos = P2::new(0.0, 0.0); // file units
    let mut file_plated: Option<bool> = None;
    let mut span = (1usize, n_copper.max(1));
    let mut pending_plated: Option<bool> = None;
    let mut in_header = false;
    let mut rout = false;
    let mut plunged = false;
    let mut absolute = true;
    let mut npth_hint = false;
    let mut holes = Vec::new();

    let num = |s: &str, metric: bool, lz_kept: Option<bool>, fmt: Option<(usize, usize)>| -> Result<f64, String> {
        if s.contains('.') {
            return s.parse::<f64>().map_err(|_| format!("bad Excellon number '{s}'"));
        }
        let (int, dec) = fmt.unwrap_or(if metric { (3, 3) } else { (2, 4) });
        let (neg, d) = match s.as_bytes().first() {
            Some(b'-') => (true, &s[1..]),
            Some(b'+') => (false, &s[1..]),
            _ => (false, s),
        };
        if d.is_empty() || !d.bytes().all(|c| c.is_ascii_digit()) {
            return Err(format!("bad Excellon number '{s}'"));
        }
        let v = if lz_kept == Some(true) {
            let mut t = d.to_string();
            while t.len() < int + dec {
                t.push('0');
            }
            t.parse::<f64>().unwrap() / 10f64.powi(dec as i32)
        } else {
            d.parse::<f64>().unwrap() / 10f64.powi(dec as i32)
        };
        Ok(if neg { -v } else { v })
    };

    for (ln, raw) in text.lines().enumerate() {
        let l = raw.trim();
        let err = |e: String| format!("Excellon line {}: {e}", ln + 1);
        if l.is_empty() {
            continue;
        }
        if let Some(c) = l.strip_prefix(';') {
            let c = c.trim();
            if let Some(r) = c.strip_prefix("#@! TF.FileFunction,") {
                let f: Vec<&str> = r.split(',').map(|s| s.trim()).collect();
                file_plated = match f[0] {
                    "Plated" => Some(true),
                    "NonPlated" => Some(false),
                    _ => None, // MixedPlating: per tool
                };
                if f.len() >= 3 {
                    if let (Ok(a), Ok(b)) = (f[1].parse::<usize>(), f[2].parse::<usize>()) {
                        span = (a.min(b), a.max(b));
                    }
                }
            } else if let Some(r) = c.strip_prefix("#@! TA.AperFunction,") {
                pending_plated = match r.split(',').next().unwrap_or("") {
                    "Plated" => Some(true),
                    "NonPlated" => Some(false),
                    _ => None,
                };
            } else if let Some(r) = c.strip_prefix("FORMAT={") {
                let spec = r.split('/').next().unwrap_or("").trim();
                if let Some((a, b)) = spec.split_once(':') {
                    if let (Ok(a), Ok(b)) = (a.trim().parse(), b.trim().parse()) {
                        fmt = Some((a, b));
                    }
                }
            } else if c.contains("NPTH") && in_header {
                npth_hint = true;
            }
            continue;
        }
        if l == "M48" {
            in_header = true;
            continue;
        }
        if l == "%" || l == "M95" {
            in_header = false;
            continue;
        }
        let units = |l: &str, metric: &mut bool, lz: &mut Option<bool>, fmt: &mut Option<(usize, usize)>| {
            let mut f = l.split(',');
            *metric = f.next() == Some("METRIC");
            for p in f {
                match p {
                    "LZ" => *lz = Some(true),
                    "TZ" => *lz = Some(false),
                    p if p.contains('.') && p.bytes().all(|c| c == b'0' || c == b'.') => {
                        let (a, b) = p.split_once('.').unwrap();
                        *fmt = Some((a.len(), b.len()));
                    }
                    _ => {}
                }
            }
        };
        if l.starts_with("METRIC") || l.starts_with("INCH") {
            units(l, &mut metric, &mut lz_kept, &mut fmt);
            continue;
        }
        let scale = if metric { 1.0 } else { MM_PER_IN };
        if let Some(t) = l.strip_prefix('T') {
            let nd = t.bytes().take_while(|c| c.is_ascii_digit()).count();
            if nd > 0 {
                let n: u32 = t[..nd].parse().unwrap();
                let rest = &t[nd..];
                if let Some(ci) = rest.find('C') {
                    let dstr: String =
                        rest[ci + 1..].chars().take_while(|c| c.is_ascii_digit() || *c == '.' || *c == '-').collect();
                    let d: f64 = dstr.parse().map_err(|_| err(format!("bad tool diameter in '{l}'")))?;
                    tools.insert(n, (d * scale, pending_plated.take()));
                    if !in_header {
                        cur_tool = Some(n);
                    }
                } else if !in_header {
                    cur_tool = if n == 0 { None } else { Some(n) };
                }
                continue;
            }
        }
        if in_header {
            continue; // FMAT, VER, ICI, ... in the header
        }
        match l {
            "G90" => {
                absolute = true;
                continue;
            }
            "G91" => {
                absolute = false;
                continue;
            }
            "G05" | "G81" => {
                rout = false;
                plunged = false;
                continue;
            }
            "M15" => {
                plunged = true;
                continue;
            }
            "M16" | "M17" => {
                plunged = false;
                continue;
            }
            "M30" | "M00" => break,
            "M71" => {
                metric = true;
                continue;
            }
            "M72" => {
                metric = false;
                continue;
            }
            _ => {}
        }
        // coordinate words, possibly prefixed by G00/G01/G02/G03 and with a G85 slot
        let (first, slot_end) = match l.split_once("G85") {
            Some((a, b)) => (a, Some(b)),
            None => (l, None),
        };
        let parse_xy = |s: &str, base: P2| -> Result<(Option<u32>, P2, bool), String> {
            let toks = tokens(s)?;
            let mut g = None;
            let mut p = if absolute { base } else { P2::new(0.0, 0.0) };
            let mut any = false;
            for (c, n) in toks {
                match c {
                    'G' => g = Some(n.parse::<u32>().map_err(|_| format!("bad G code in '{s}'"))?),
                    'X' => {
                        p.x = num(n, metric, lz_kept, fmt)?;
                        any = true;
                    }
                    'Y' => {
                        p.y = num(n, metric, lz_kept, fmt)?;
                        any = true;
                    }
                    'I' | 'J' | 'A' | 'F' | 'S' => {}
                    'M' => {}
                    _ => return Err(format!("unsupported Excellon word '{s}'")),
                }
            }
            if !absolute {
                p = add(p, base);
            }
            Ok((g, p, any))
        };
        let (g, p, any) = parse_xy(first, pos).map_err(err)?;
        let tool = |cur: Option<u32>| -> Result<(f64, bool), String> {
            let t = cur.ok_or("drill hit with no tool selected")?;
            let &(d, pl) = tools.get(&t).ok_or_else(|| format!("tool T{t} not defined"))?;
            Ok((d, pl.or(file_plated).unwrap_or(!npth_hint)))
        };
        let mm = |q: P2| P2::new(q.x * scale, q.y * scale);
        match g {
            Some(0) => {
                rout = true;
                plunged = false;
                pos = p;
            }
            Some(1) | Some(2) | Some(3) if rout => {
                if plunged && any {
                    let (d, plated) = tool(cur_tool).map_err(err)?;
                    holes.push(Hole { c: mm(pos), d, plated, from: span.0, to: span.1, end: Some(mm(p)) });
                }
                pos = p;
            }
            _ => {
                if !any {
                    continue;
                }
                pos = p;
                let (d, plated) = tool(cur_tool).map_err(err)?;
                let end = match slot_end {
                    Some(se) => {
                        let (_, e, _) = parse_xy(se, pos).map_err(err)?;
                        pos = e;
                        Some(mm(e))
                    }
                    None => None,
                };
                holes.push(Hole { c: mm(p), d, plated, from: span.0, to: span.1, end });
            }
        }
    }
    Ok(holes)
}

// ------------------------------------------------------------------------------------------------------------
// Gerber job file stackup
// ------------------------------------------------------------------------------------------------------------

fn json_f64(v: Option<&serde_json::Value>) -> Option<f64> {
    match v? {
        serde_json::Value::Number(n) => n.as_f64(),
        serde_json::Value::String(s) => s.trim().parse().ok(),
        _ => None,
    }
}

/// A generic stackup for `n` copper layers and the given board thickness (used when the job file has none).
fn default_stack(n: usize, board: f64) -> Vec<StackLayer> {
    let n = n.max(1);
    let cu = 0.035;
    let diel = if n > 1 { ((board - n as f64 * cu) / (n - 1) as f64).max(0.05) } else { board };
    let mut v = vec![StackLayer {
        name: "Top Solder Mask".into(),
        kind: StackKind::SolderMask,
        thickness_mm: 0.01,
        eps_r: 3.8,
        loss_tan: 0.025,
    }];
    for k in 1..=n {
        let name = if k == 1 {
            "F.Cu".to_string()
        } else if k == n {
            "B.Cu".to_string()
        } else {
            format!("In{}.Cu", k - 1)
        };
        v.push(StackLayer { name, kind: StackKind::Copper, thickness_mm: cu, eps_r: 1.0, loss_tan: 0.0 });
        if k < n {
            v.push(StackLayer {
                name: format!("dielectric {k}"),
                kind: StackKind::Dielectric,
                thickness_mm: diel,
                eps_r: 4.5,
                loss_tan: 0.02,
            });
        }
    }
    v.push(StackLayer {
        name: "Bottom Solder Mask".into(),
        kind: StackKind::SolderMask,
        thickness_mm: 0.01,
        eps_r: 3.8,
        loss_tan: 0.025,
    });
    v
}

/// The physical stackup of a Gerber job file (`MaterialStackup`, top to bottom; legend and paste layers are
/// skipped). Without a `MaterialStackup` a generic stack for `GeneralSpecs.LayerNumber` copper layers and
/// `BoardThickness` is returned. Copper layers get eps_r 1, loss_tan 0.
pub fn parse_gbrjob(text: &str) -> Result<Vec<StackLayer>, String> {
    let j: serde_json::Value = serde_json::from_str(text).map_err(|e| format!("gbrjob: {e}"))?;
    let specs = j.get("GeneralSpecs");
    let board_t = json_f64(specs.and_then(|s| s.get("BoardThickness"))).unwrap_or(1.6);
    let Some(ms) = j.get("MaterialStackup").and_then(|v| v.as_array()) else {
        let n = json_f64(specs.and_then(|s| s.get("LayerNumber"))).unwrap_or(2.0) as usize;
        return Ok(default_stack(n, board_t));
    };
    let mut out = Vec::new();
    let mut missing = Vec::new();
    for e in ms {
        let ty = e.get("Type").and_then(|v| v.as_str()).unwrap_or("");
        let kind = match ty {
            "Copper" => StackKind::Copper,
            "Dielectric" => StackKind::Dielectric,
            "SolderMask" => StackKind::SolderMask,
            _ => continue,
        };
        let name = e.get("Name").and_then(|v| v.as_str()).unwrap_or(ty).to_string();
        let t = json_f64(e.get("Thickness"));
        let (de, dl) = match kind {
            StackKind::Copper => (1.0, 0.0),
            StackKind::Dielectric => (4.5, 0.02),
            StackKind::SolderMask => (3.8, 0.025),
        };
        let eps_r = if kind == StackKind::Copper { 1.0 } else { json_f64(e.get("DielectricConstant")).unwrap_or(de) };
        let loss_tan = if kind == StackKind::Copper { 0.0 } else { json_f64(e.get("LossTangent")).unwrap_or(dl) };
        let thickness_mm = match (t, kind) {
            (Some(t), _) => t,
            (None, StackKind::Copper) => 0.035,
            (None, StackKind::SolderMask) => 0.01,
            (None, StackKind::Dielectric) => {
                missing.push(out.len());
                0.0
            }
        };
        out.push(StackLayer { name, kind, thickness_mm, eps_r, loss_tan });
    }
    if !missing.is_empty() {
        let used: f64 = out.iter().filter(|l| l.kind != StackKind::SolderMask).map(|l| l.thickness_mm).sum();
        let each = ((board_t - used) / missing.len() as f64).max(0.05);
        for k in missing {
            out[k].thickness_mm = each;
        }
    }
    if !out.iter().any(|l| l.kind == StackKind::Copper) {
        return Err("gbrjob: MaterialStackup has no copper layers".into());
    }
    Ok(out)
}

// ------------------------------------------------------------------------------------------------------------
// Directory loader
// ------------------------------------------------------------------------------------------------------------

/// Classification of one file of an output directory.
enum FileRole {
    /// Copper layer k (1-based) with side token (Top / Inr / Bot).
    Copper(usize, String),
    Profile,
    Drill,
    Job,
    Other,
}

fn classify_function(ff: &str) -> FileRole {
    let f: Vec<&str> = ff.split(',').map(|s| s.trim()).collect();
    match f[0] {
        "Copper" => {
            let k = f.get(1).and_then(|l| l.strip_prefix('L')).and_then(|n| n.parse::<usize>().ok());
            match k {
                Some(k) => FileRole::Copper(k, f.get(2).copied().unwrap_or("").to_string()),
                None => FileRole::Other,
            }
        }
        "Profile" => FileRole::Profile,
        "Plated" | "NonPlated" | "MixedPlating" => FileRole::Drill,
        _ => FileRole::Other,
    }
}

fn is_excellon(text: &str) -> bool {
    text.lines().take(5).any(|l| l.trim() == "M48")
}

/// Load a KiCad Gerber/Excellon output directory. Files are identified by their X2 file function
/// (`%TF.FileFunction%`, `; #@! TF.FileFunction`), then by the job file's `FilesAttributes`, then by the usual
/// KiCad/Protel extensions.
pub fn load_dir(dir: &Path) -> Result<Board, String> {
    let mut entries: Vec<_> = std::fs::read_dir(dir)
        .map_err(|e| format!("{}: {e}", dir.display()))?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.is_file())
        .collect();
    entries.sort();
    let mut files: Vec<(std::path::PathBuf, String, FileRole)> = Vec::new();
    for p in entries {
        let Ok(bytes) = std::fs::read(&p) else {
            continue;
        };
        let text = String::from_utf8_lossy(&bytes).into_owned();
        let ext = p.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
        let role = if ext == "gbrjob" || (text.trim_start().starts_with('{') && text.contains("\"MaterialStackup\"")) {
            FileRole::Job
        } else if is_excellon(&text) {
            FileRole::Drill
        } else if let Some(ff) = gerber_file_function(&text) {
            classify_function(&ff)
        } else {
            FileRole::Other
        };
        files.push((p, text, role));
    }

    // job file: stackup, copper names, file functions of X1 files
    let job_text = files.iter().find(|f| matches!(f.2, FileRole::Job)).map(|f| f.1.clone());
    let mut job_functions: HashMap<String, String> = HashMap::new();
    let mut job_layer_count = None;
    if let Some(jt) = &job_text {
        if let Ok(j) = serde_json::from_str::<serde_json::Value>(jt) {
            if let Some(fa) = j.get("FilesAttributes").and_then(|v| v.as_array()) {
                for a in fa {
                    if let (Some(p), Some(f)) =
                        (a.get("Path").and_then(|v| v.as_str()), a.get("FileFunction").and_then(|v| v.as_str()))
                    {
                        job_functions.insert(p.to_string(), f.to_string());
                    }
                }
            }
            job_layer_count = json_f64(j.get("GeneralSpecs").and_then(|s| s.get("LayerNumber"))).map(|v| v as usize);
        }
    }
    let have_copper = files.iter().any(|f| matches!(f.2, FileRole::Copper(..)));
    let have_profile = files.iter().any(|f| matches!(f.2, FileRole::Profile));
    for f in files.iter_mut() {
        if !matches!(f.2, FileRole::Other) {
            continue;
        }
        let fname = f.0.file_name().and_then(|n| n.to_str()).unwrap_or("").to_string();
        if let Some(ff) = job_functions.get(&fname) {
            f.2 = classify_function(ff);
            continue;
        }
        if have_copper && have_profile {
            continue;
        }
        let ext = f.0.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
        let n = job_layer_count.unwrap_or(2);
        f.2 = match ext.as_str() {
            "gtl" if !have_copper => FileRole::Copper(1, "Top".into()),
            "gbl" if !have_copper => FileRole::Copper(n, "Bot".into()),
            "gm1" | "gko" if !have_profile => FileRole::Profile,
            e if !have_copper && e.len() >= 2 && e.starts_with('g') && e[1..].bytes().all(|c| c.is_ascii_digit()) => {
                FileRole::Copper(e[1..].parse::<usize>().unwrap() + 1, "Inr".into())
            }
            _ => FileRole::Other,
        };
    }

    let stack = match &job_text {
        Some(t) => parse_gbrjob(t)?,
        None => Vec::new(),
    };

    // copper layers (parsed in parallel)
    use rayon::prelude::*;
    let copper: Vec<(usize, String, &std::path::PathBuf, &String)> = files
        .iter()
        .filter_map(|f| match &f.2 {
            FileRole::Copper(k, side) => Some((*k, side.clone(), &f.0, &f.1)),
            _ => None,
        })
        .collect();
    if copper.is_empty() {
        return Err(format!("{}: no copper layer Gerber files found", dir.display()));
    }
    let n_copper = copper.iter().map(|c| c.0).max().unwrap().max(job_layer_count.unwrap_or(0));
    let cu_names: Vec<String> = stack.iter().filter(|l| l.kind == StackKind::Copper).map(|l| l.name.clone()).collect();
    let mut layers: Vec<CopperLayer> = copper
        .par_iter()
        .map(|(k, side, path, text)| -> Result<CopperLayer, String> {
            let prims = parse_gerber(text).map_err(|e| format!("{}: {e}", path.display()))?;
            let mut bbox = BBox::empty();
            for p in &prims {
                bbox.union(&p.shape.bbox());
            }
            let name = if cu_names.len() == n_copper {
                cu_names[k - 1].clone()
            } else if side == "Top" || *k == 1 {
                "F.Cu".to_string()
            } else if side == "Bot" || *k == n_copper {
                "B.Cu".to_string()
            } else {
                format!("In{}.Cu", k - 1)
            };
            Ok(CopperLayer { name, index: *k, prims, bbox })
        })
        .collect::<Result<_, _>>()?;
    layers.sort_by_key(|l| l.index);

    let stack = if stack.is_empty() { default_stack(n_copper, 1.6) } else { stack };

    let mut outline = Vec::new();
    for f in files.iter().filter(|f| matches!(f.2, FileRole::Profile)) {
        outline.extend(parse_outline(&f.1).map_err(|e| format!("{}: {e}", f.0.display()))?);
    }
    outline.sort_by(|a, b| poly_area(b).abs().partial_cmp(&poly_area(a).abs()).unwrap());

    let mut holes = Vec::new();
    for f in files.iter().filter(|f| matches!(f.2, FileRole::Drill)) {
        if !is_excellon(&f.1) {
            continue; // e.g. Gerber X2 drill files: not supported, Excellon expected
        }
        holes.extend(parse_excellon(&f.1, n_copper).map_err(|e| format!("{}: {e}", f.0.display()))?);
    }
    Ok(Board { layers, outline, holes, stack })
}

// ------------------------------------------------------------------------------------------------------------
// SVG debug rendering
// ------------------------------------------------------------------------------------------------------------

fn svg_shape(s: &Shape, out: &mut String) {
    use std::fmt::Write;
    match s {
        Shape::Circle { c, r } => {
            let _ = write!(out, "<circle cx=\"{:.4}\" cy=\"{:.4}\" r=\"{:.4}\"/>", c.x, c.y, r);
        }
        Shape::Capsule { a, b, r } => {
            let _ = write!(
                out,
                "<line x1=\"{:.4}\" y1=\"{:.4}\" x2=\"{:.4}\" y2=\"{:.4}\" stroke-width=\"{:.4}\"/>",
                a.x,
                a.y,
                b.x,
                b.y,
                2.0 * r
            );
        }
        Shape::Polygon(p) => {
            out.push_str("<path d=\"");
            for (k, q) in p.iter().enumerate() {
                let _ = write!(out, "{}{:.4} {:.4}", if k == 0 { "M" } else { "L" }, q.x, q.y);
            }
            out.push_str("Z\"/>");
        }
    }
}

/// Write an SVG of the board for visual debugging: each copper layer a semi-transparent group (toggle by
/// clicking the legend when opened in a browser), holes, and the outline.
pub fn render_svg(board: &Board, path: &Path) {
    use std::fmt::Write;
    let mut bb = BBox::empty();
    for l in &board.outline {
        l.iter().for_each(|&p| bb.add(p));
    }
    if bb.is_empty() {
        for l in &board.layers {
            bb.union(&l.bbox);
        }
    }
    if bb.is_empty() {
        bb = BBox { x0: 0.0, y0: 0.0, x1: 1.0, y1: 1.0 };
    }
    let bb = bb.grow(1.0);
    let legend_w = 12.0;
    let (w, h) = (bb.x1 - bb.x0 + legend_w, bb.y1 - bb.y0);
    let colors = ["#d62728", "#ff7f0e", "#2ca02c", "#17becf", "#9467bd", "#1f4fd6", "#8c564b", "#e377c2"];
    let mut defs = String::new();
    let mut style = String::from(
        "path{fill-rule:evenodd} .clr path,.clr circle{fill:#000;stroke:none} .clr line{stroke:#000;stroke-linecap:round} \
         text{font-family:sans-serif;font-size:1.2px;cursor:pointer}\n",
    );
    let mut body = String::new();
    let mut mask_id = 0;
    for (li, layer) in board.layers.iter().enumerate() {
        let col = colors[li % colors.len()];
        let cls = format!("L{li}");
        let _ = writeln!(
            style,
            ".{cls} path,.{cls} circle{{fill:{col};stroke:none}} .{cls} line{{stroke:{col};stroke-linecap:round;fill:none}}"
        );
        // build content with masks for clear runs
        let mut content = String::new();
        let mut k = 0;
        while k < layer.prims.len() {
            let dark = layer.prims[k].dark;
            let mut run = String::new();
            while k < layer.prims.len() && layer.prims[k].dark == dark {
                svg_shape(&layer.prims[k].shape, &mut run);
                run.push('\n');
                k += 1;
            }
            if dark {
                content.push_str(&run);
            } else {
                mask_id += 1;
                let _ = writeln!(
                    defs,
                    "<mask id=\"m{mask_id}\" maskUnits=\"userSpaceOnUse\" x=\"{:.3}\" y=\"{:.3}\" width=\"{:.3}\" height=\"{:.3}\">\
                     <rect x=\"{:.3}\" y=\"{:.3}\" width=\"{:.3}\" height=\"{:.3}\" fill=\"#fff\"/><g class=\"clr\">{run}</g></mask>",
                    bb.x0, bb.y0, w, h, bb.x0, bb.y0, w, h
                );
                content = format!("<g mask=\"url(#m{mask_id})\">{content}</g>\n");
            }
        }
        let _ = writeln!(body, "<g id=\"layer{li}\" class=\"{cls}\" opacity=\"0.5\">{content}</g>");
    }
    // holes
    body.push_str("<g id=\"holes\">\n");
    for hole in &board.holes {
        let col = if hole.plated { "#222" } else { "#888" };
        match hole.end {
            Some(e) => {
                let _ = writeln!(
                    body,
                    "<line x1=\"{:.4}\" y1=\"{:.4}\" x2=\"{:.4}\" y2=\"{:.4}\" stroke=\"{col}\" stroke-width=\"{:.4}\" stroke-linecap=\"round\" opacity=\"0.8\"/>",
                    hole.c.x, hole.c.y, e.x, e.y, hole.d
                );
            }
            None => {
                let _ = writeln!(
                    body,
                    "<circle cx=\"{:.4}\" cy=\"{:.4}\" r=\"{:.4}\" fill=\"#fff\" stroke=\"{col}\" stroke-width=\"0.03\"/>",
                    hole.c.x,
                    hole.c.y,
                    hole.d / 2.0
                );
            }
        }
    }
    body.push_str("</g>\n<g id=\"outline\" fill=\"none\" stroke=\"#000\" stroke-width=\"0.06\">\n");
    for l in &board.outline {
        let mut d = String::new();
        svg_shape(&Shape::Polygon(l.clone()), &mut d);
        body.push_str(&d.replace("<path ", "<path style=\"fill:none\" "));
        body.push('\n');
    }
    body.push_str("</g>\n");
    // legend (not flipped)
    let mut legend = String::new();
    let lx = bb.x1 + 1.0;
    for (li, layer) in board.layers.iter().enumerate() {
        let col = colors[li % colors.len()];
        let y = -bb.y1 + 2.0 + 2.0 * li as f64;
        let _ = writeln!(
            legend,
            "<g onclick=\"tg('layer{li}')\"><rect x=\"{lx:.3}\" y=\"{:.3}\" width=\"1.2\" height=\"1.2\" fill=\"{col}\" opacity=\"0.6\"/>\
             <text x=\"{:.3}\" y=\"{:.3}\">{} (L{})</text></g>",
            y - 1.0,
            lx + 1.6,
            y,
            layer.name,
            layer.index
        );
    }
    let n = board.layers.len();
    for (k, (id, label)) in [("holes", "holes"), ("outline", "outline")].iter().enumerate() {
        let y = -bb.y1 + 2.0 + 2.0 * (n + k) as f64;
        let _ = writeln!(
            legend,
            "<g onclick=\"tg('{id}')\"><rect x=\"{lx:.3}\" y=\"{:.3}\" width=\"1.2\" height=\"1.2\" fill=\"#444\"/>\
             <text x=\"{:.3}\" y=\"{:.3}\">{label}</text></g>",
            y - 1.0,
            lx + 1.6,
            y
        );
    }
    let svg = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"{:.3} {:.3} {:.3} {:.3}\" width=\"{:.0}mm\" height=\"{:.0}mm\">\n\
         <style>{style}</style>\n<script><![CDATA[function tg(id){{var g=document.getElementById(id);g.style.display=g.style.display=='none'?'':'none';}}]]></script>\n\
         <defs>{defs}</defs>\n<rect x=\"{:.3}\" y=\"{:.3}\" width=\"{:.3}\" height=\"{:.3}\" fill=\"#fff\"/>\n\
         <g transform=\"scale(1,-1)\">\n{body}</g>\n{legend}</svg>\n",
        bb.x0,
        -bb.y1,
        w,
        h,
        w * 4.0,
        h * 4.0,
        bb.x0,
        -bb.y1,
        w,
        h
    );
    if let Err(e) = std::fs::write(path, svg) {
        eprintln!("render_svg: {}: {e}", path.display());
    }
}

// ------------------------------------------------------------------------------------------------------------
// Tests
// ------------------------------------------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const HDR: &str = "%FSLAX46Y46*%\n%MOMM*%\n%LPD*%\nG01*\n";
    const ROUNDRECT: &str = "%AMRoundRect*\n0 Rectangle with rounded corners*\n0 $1 Rounding radius*\n\
        4,1,4,$2,$3,$4,$5,$6,$7,$8,$9,$2,$3,0*\n1,1,$1+$1,$2,$3*\n1,1,$1+$1,$4,$5*\n1,1,$1+$1,$6,$7*\n1,1,$1+$1,$8,$9*\n\
        20,1,$1+$1,$2,$3,$4,$5,0*\n20,1,$1+$1,$4,$5,$6,$7,0*\n20,1,$1+$1,$6,$7,$8,$9,0*\n20,1,$1+$1,$8,$9,$2,$3,0*%\n";

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    fn union_contains(prims: &[Prim], p: P2) -> bool {
        let mut inside = false;
        for q in prims {
            if q.shape.contains(p) {
                inside = q.dark;
            }
        }
        inside
    }

    #[test]
    fn gerber_flash_circle() {
        let g =
            format!("{HDR}%ADD10C,0.500000*%\nD10*\n%TO.P,U1,3*%\n%TO.N,VCC*%\nX1000000Y-2000000D03*\n%TD*%\nM02*\n");
        let p = parse_gerber(&g).unwrap();
        assert_eq!(p.len(), 1);
        match p[0].shape {
            Shape::Circle { c, r } => {
                assert!(close(c.x, 1.0, 1e-9) && close(c.y, -2.0, 1e-9) && close(r, 0.25, 1e-12));
            }
            ref s => panic!("expected circle, got {s:?}"),
        }
        assert!(p[0].dark);
        assert_eq!(p[0].net.as_deref(), Some("VCC"));
        assert_eq!(p[0].pad, Some(("U1".to_string(), "3".to_string())));
    }

    #[test]
    fn gerber_capsule_draw_and_attributes() {
        let g = format!(
            "{HDR}%ADD11C,0.200000*%\n%TO.N,GND*%\nD11*\nX0Y0D02*\nX1000000Y0D01*\nY500000D01*\n%TD*%\n\
             %TO.N,N/C*%\nX2000000Y0D02*\nX3000000D01*\n%TD*%\nX4000000Y0D02*\nX5000000D01*\nM02*\n"
        );
        let p = parse_gerber(&g).unwrap();
        assert_eq!(p.len(), 4);
        assert_eq!(p[0].shape, Shape::Capsule { a: P2::new(0.0, 0.0), b: P2::new(1.0, 0.0), r: 0.1 });
        assert_eq!(p[1].shape, Shape::Capsule { a: P2::new(1.0, 0.0), b: P2::new(1.0, 0.5), r: 0.1 });
        assert_eq!(p[0].net.as_deref(), Some("GND"));
        assert_eq!(p[1].net.as_deref(), Some("GND"));
        assert_eq!(p[2].net, None, "N/C -> None");
        assert_eq!(p[3].net, None, "TD clears");
    }

    #[test]
    fn gerber_rect_draw_is_swept_hull() {
        let g = format!("{HDR}%ADD10R,0.2X0.4*%\nD10*\nX0Y0D02*\nX1000000Y0D01*\nM02*\n");
        let p = parse_gerber(&g).unwrap();
        assert_eq!(p.len(), 1);
        let Shape::Polygon(pts) = &p[0].shape else { panic!() };
        assert_eq!(pts.len(), 4);
        assert!(close(poly_area(pts).abs(), 1.2 * 0.4, 1e-9));
    }

    #[test]
    fn gerber_region_with_arc() {
        // half disc of radius 1: (-1,0) -> (1,0) straight, then CCW arc back to (-1,0) about the origin
        let g = format!(
            "{HDR}G75*\nG36*\nX-1000000Y0D02*\nG01*\nX1000000Y0D01*\nG03*\nX-1000000Y0I-1000000J0D01*\nG01*\nG37*\nM02*\n"
        );
        let p = parse_gerber(&g).unwrap();
        assert_eq!(p.len(), 1);
        let Shape::Polygon(pts) = &p[0].shape else { panic!() };
        assert!(close(poly_area(pts), PI / 2.0, 0.01), "area {}", poly_area(pts));
        for w in pts.windows(2) {
            let (a, b) = (w[0], w[1]);
            if a.y > 1e-9 || b.y > 1e-9 {
                assert!(close(a.dist(P2::new(0.0, 0.0)), 1.0, 1e-9));
                let m = P2::new((a.x + b.x) / 2.0, (a.y + b.y) / 2.0);
                assert!(1.0 - m.dist(P2::new(0.0, 0.0)) <= CHORD_TOL + 1e-12);
            }
        }
        assert!(pts.iter().all(|q| q.y >= -1e-12));
        // two contours in one region -> two polygons
        let g2 = format!(
            "{HDR}G36*\nX0Y0D02*\nX1000000Y0D01*\nX1000000Y1000000D01*\nX0Y0D01*\n\
             X5000000Y0D02*\nX6000000Y0D01*\nX6000000Y1000000D01*\nX5000000Y0D01*\nG37*\nM02*\n"
        );
        let p2 = parse_gerber(&g2).unwrap();
        assert_eq!(p2.len(), 2);
        assert!(p2.iter().all(|p| matches!(&p.shape, Shape::Polygon(v) if v.len() == 3)));
    }

    #[test]
    fn gerber_roundrect_macro_flash() {
        let g = format!(
            "{HDR}{ROUNDRECT}%ADD10RoundRect,0.250000X-0.500000X-0.500000X0.500000X-0.500000X0.500000X0.500000X-0.500000X0.500000X0*%\n\
             D10*\nX10000000Y10000000D03*\nM02*\n"
        );
        let p = parse_gerber(&g).unwrap();
        assert_eq!(p.len(), 9, "1 outline + 4 circles + 4 vector lines");
        let at = |x: f64, y: f64| union_contains(&p, P2::new(10.0 + x, 10.0 + y));
        assert!(at(0.0, 0.0));
        assert!(at(0.74, 0.0) && at(0.0, -0.74) && at(0.6, 0.6));
        assert!(!at(0.74, 0.74), "rounded corner");
        assert!(!at(0.76, 0.0));
        // macro with rotation: a 2 x 0.5 centre line rotated by 90 degrees
        let g = format!("{HDR}%AMROT*\n21,1,2,0.5,0,0,90*%\n%ADD12ROT*%\nD12*\nX0Y0D03*\nM02*\n");
        let p = parse_gerber(&g).unwrap();
        let b = p[0].shape.bbox();
        assert!(close(b.x1 - b.x0, 0.5, 1e-9) && close(b.y1 - b.y0, 2.0, 1e-9));
        // exposure off -> clear prim; variable definitions
        let g =
            format!("{HDR}%AMDONUT*\n$3=$1x0.5*\n1,1,$1,0,0*\n1,0,$3,0,0*%\n%ADD13DONUT,1.0*%\nD13*\nX0Y0D03*\nM02*\n");
        let p = parse_gerber(&g).unwrap();
        assert_eq!(p.len(), 2);
        assert!(p[0].dark && !p[1].dark);
        assert_eq!(p[1].shape, Shape::Circle { c: P2::new(0.0, 0.0), r: 0.25 });
    }

    #[test]
    fn gerber_polarity_and_units() {
        let g = format!(
            "{HDR}%ADD10R,2X2*%\n%ADD11C,1*%\nD10*\nX0Y0D03*\n%LPC*%\nD11*\nX0Y0D03*\n%LPD*%\nX5000000Y0D03*\nM02*\n"
        );
        let p = parse_gerber(&g).unwrap();
        assert_eq!(p.len(), 3);
        assert!(p[0].dark && !p[1].dark && p[2].dark);
        assert!(!union_contains(&p, P2::new(0.0, 0.0)));
        assert!(union_contains(&p, P2::new(0.9, 0.9)));
        // inches, trailing-zero-omitted format
        let g = "%FSLAX24Y24*%\n%MOIN*%\n%ADD10C,0.1*%\nD10*\nX10000Y-5000D03*\nM02*\n";
        let p = parse_gerber(g).unwrap();
        assert_eq!(p[0].shape, Shape::Circle { c: P2::new(25.4, -12.7), r: 1.27 });
        let g = "%FSTAX24Y24*%\n%MOMM*%\n%ADD10C,0.1*%\nD10*\nX01Y-005D03*\nM02*\n";
        let p = parse_gerber(g).unwrap();
        assert_eq!(p[0].shape, Shape::Circle { c: P2::new(1.0, -0.5), r: 0.05 });
    }

    #[test]
    fn gerber_obround_polygon_and_transforms() {
        let g = format!("{HDR}%ADD10O,2X1*%\n%ADD11P,2X6*%\nD10*\nX0Y0D03*\nD11*\n%LR90*%\n%LS2*%\nX0Y0D03*\nM02*\n");
        let p = parse_gerber(&g).unwrap();
        let Shape::Polygon(o) = &p[0].shape else { panic!() };
        assert!(close(poly_area(o), 1.0 + PI * 0.25, 0.005));
        let Shape::Polygon(h) = &p[1].shape else { panic!() };
        assert_eq!(h.len(), 6);
        // hexagon of diameter 2, scaled 2x and rotated 90: vertex on the y axis at distance 2
        assert!(h.iter().any(|q| close(q.x, 0.0, 1e-9) && close(q.y, 2.0, 1e-9)));
    }

    #[test]
    fn gerber_step_repeat_and_errors() {
        let g = format!("{HDR}%ADD10C,0.1*%\nD10*\n%SRX3Y2I1J2*%\nX0Y0D03*\n%SR*%\nM02*\n");
        assert_eq!(parse_gerber(&g).unwrap().len(), 6);
        let g = format!("{HDR}%ABD12*%\n%AB*%\nM02*\n");
        assert!(parse_gerber(&g).unwrap_err().contains("block apertures"));
    }

    #[test]
    fn expressions() {
        let mut v = HashMap::new();
        v.insert(1, 0.15);
        v.insert(2, 2.0);
        assert!(close(eval_expr("$1+$1", &v).unwrap(), 0.3, 1e-12));
        assert!(close(eval_expr("-(2+3)x4/5", &v).unwrap(), -4.0, 1e-12));
        assert!(close(eval_expr("2x-3", &v).unwrap(), -6.0, 1e-12));
        assert!(close(eval_expr("$2-$1X2", &v).unwrap(), 1.7, 1e-12));
        assert!(close(eval_expr("1+2x3", &v).unwrap(), 7.0, 1e-12));
        assert!(close(eval_expr("10/4/5", &v).unwrap(), 0.5, 1e-12));
        assert!(close(eval_expr(" -$9 + .5 ", &v).unwrap(), 0.5, 1e-12));
        assert!(eval_expr("(1+2", &v).is_err());
        assert!(eval_expr("1+", &v).is_err());
    }

    #[test]
    fn excellon_holes_and_slots() {
        let x = "M48\n; #@! TF.FileFunction,Plated,1,4,PTH\nFMAT,2\nMETRIC\n\
                 ; #@! TA.AperFunction,Plated,PTH,ViaDrill\nT1C0.300\nT2C0.600\n%\nG90\nG05\nT1\nX1.0Y-2.0\nY-3.0\n\
                 T2\nX0.875Y-14.732G85X1.925Y-14.732\nG05\nG00X5.0Y5.0\nM15\nG01X5.0Y6.0\nM16\nG05\nM30\n";
        let h = parse_excellon(x, 4).unwrap();
        assert_eq!(h.len(), 4);
        assert!(close(h[0].c.x, 1.0, 1e-12) && close(h[0].d, 0.3, 1e-12) && h[0].end.is_none());
        assert!(close(h[1].c.y, -3.0, 1e-12) && close(h[1].c.x, 1.0, 1e-12), "modal X");
        assert!(h.iter().all(|h| h.plated && h.from == 1 && h.to == 4));
        let e = h[2].end.unwrap();
        assert!(close(h[2].c.x, 0.875, 1e-12) && close(e.x, 1.925, 1e-12) && close(h[2].d, 0.6, 1e-12));
        assert_eq!(h[3].end, Some(P2::new(5.0, 6.0)));
        // inch, leading zeros suppressed (TZ), non-plated
        let x =
            "M48\n; #@! TF.FileFunction,NonPlated,1,2,NPTH\n; FORMAT={2:4/ absolute / inch / suppress leading zeros}\n\
                 INCH,TZ\nT1C0.0394\n%\nT1\nX35071Y-15581\nM30\n";
        let h = parse_excellon(x, 2).unwrap();
        assert_eq!(h.len(), 1);
        assert!(!h[0].plated);
        assert!(close(h[0].c.x, 3.5071 * 25.4, 1e-9) && close(h[0].c.y, -1.5581 * 25.4, 1e-9));
        assert!(close(h[0].d, 0.0394 * 25.4, 1e-12));
        // leading zeros kept (LZ)
        let x = "M48\nINCH,LZ\nT1C0.0394\n%\nT1\nX035071Y-015581\nM30\n";
        let h = parse_excellon(x, 2).unwrap();
        assert!(close(h[0].c.x, 3.5071 * 25.4, 1e-9) && close(h[0].c.y, -1.5581 * 25.4, 1e-9));
    }

    #[test]
    fn outline_chaining() {
        // square drawn with segments out of order and one reversed, plus an inner circle
        let g = format!(
            "{HDR}%ADD10C,0.1*%\nD10*\nX0Y0D02*\nX10000000Y0D01*\nX10000000Y10000000D02*\nX0Y10000000D01*\n\
             X10000000Y10000000D02*\nX10000000Y0D01*\nX0Y0D02*\nX0Y10000000D01*\n\
             G75*\nX6000000Y5000000D02*\nG02*\nX6000000Y5000000I-1000000J0D01*\nG01*\nM02*\n"
        );
        let l = parse_outline(&g).unwrap();
        assert_eq!(l.len(), 2);
        assert_eq!(l[0].len(), 4);
        assert!(close(poly_area(&l[0]).abs(), 100.0, 1e-9));
        assert!(close(poly_area(&l[1]).abs(), PI, 0.02));
    }

    #[test]
    fn gbrjob_stackup() {
        let j = r#"{"GeneralSpecs":{"LayerNumber":2,"BoardThickness":1.6},"MaterialStackup":[
            {"Type":"Legend","Name":"Top Silk Screen"},
            {"Type":"SolderMask","Thickness":0.01,"Name":"Top Solder Mask"},
            {"Type":"Copper","Thickness":0.035,"Name":"F.Cu"},
            {"Type":"Dielectric","Thickness":1.51,"Material":"FR4","Name":"F.Cu/B.Cu","DielectricConstant":4.1,"LossTangent":0.015},
            {"Type":"Copper","Thickness":0.035,"Name":"B.Cu"},
            {"Type":"SolderMask","Thickness":0.01,"Name":"Bottom Solder Mask"}]}"#;
        let s = parse_gbrjob(j).unwrap();
        assert_eq!(s.len(), 5);
        assert_eq!(s[0].kind, StackKind::SolderMask);
        assert!(close(s[0].eps_r, 3.8, 1e-12) && close(s[0].loss_tan, 0.025, 1e-12));
        assert_eq!(s[2].kind, StackKind::Dielectric);
        assert!(close(s[2].eps_r, 4.1, 1e-12) && close(s[2].loss_tan, 0.015, 1e-12));
        assert_eq!(s[3].name, "B.Cu");
    }

    const AMOLED_DIR: &str = "/private/tmp/claude-501/-Users-janekmann-src-gridroute--claude-worktrees-pcb-field-solver-a1d8e5/c00ff31d-4319-416d-ba20-cc0979fdc64c/scratchpad/amoled_gbr";

    #[test]
    fn amoled_board_loads() {
        let dir = Path::new(AMOLED_DIR);
        if !dir.exists() {
            eprintln!("skipping: {AMOLED_DIR} not present");
            return;
        }
        let b = load_dir(dir).unwrap();
        let names: Vec<&str> = b.layers.iter().map(|l| l.name.as_str()).collect();
        assert_eq!(names, ["F.Cu", "In1.Cu", "In2.Cu", "In3.Cu", "In4.Cu", "B.Cu"]);
        assert_eq!(b.layers.iter().map(|l| l.index).collect::<Vec<_>>(), [1, 2, 3, 4, 5, 6]);
        assert!(b.layers.iter().all(|l| !l.prims.is_empty()));
        assert!(b.layers.iter().any(|l| l.prims.iter().any(|p| p.net.as_deref() == Some("GND"))));
        assert!(b.layers[0].prims.iter().any(|p| p.pad.is_some()));
        assert!(!b.outline.is_empty());
        let mut ob = BBox::empty();
        b.outline[0].iter().for_each(|&p| ob.add(p));
        assert!(close(ob.x1 - ob.x0, 33.1, 0.1), "outline width {}", ob.x1 - ob.x0);
        assert!(close(ob.y1 - ob.y0, 40.7, 0.1), "outline height {}", ob.y1 - ob.y0);
        assert!(b.holes.len() > 50, "{} holes", b.holes.len());
        let ogrow = ob.grow(1e-6);
        assert!(b.holes.iter().all(|h| ogrow.contains(h.c) && h.end.is_none_or(|e| ogrow.contains(e))));
        assert!(b.holes.iter().any(|h| h.end.is_some()), "the AMOLED board has slots");
        assert!(b.holes.iter().any(|h| !h.plated));
        let t: f64 = b.stack.iter().filter(|l| l.kind != StackKind::SolderMask).map(|l| l.thickness_mm).sum();
        // 6 x 35 um copper + 5 x 154 um FR4 = 0.98 mm; with the two 10 um masks it is the 1.0 mm board
        assert!(close(t, 1.0, 0.03), "stack thickness {t}");
        let total: f64 = b.stack.iter().map(|l| l.thickness_mm).sum();
        assert!(close(total, 1.0, 0.005), "total stack thickness {total}");
        assert_eq!(b.stack.iter().filter(|l| l.kind == StackKind::Copper).count(), 6);
        for l in &b.layers {
            assert!(
                l.bbox.x0 >= ob.x0 - 0.5
                    && l.bbox.x1 <= ob.x1 + 0.5
                    && l.bbox.y0 >= ob.y0 - 0.5
                    && l.bbox.y1 <= ob.y1 + 0.5
            );
        }
        render_svg(&b, &dir.parent().unwrap().join("amoled.svg"));
    }
}
