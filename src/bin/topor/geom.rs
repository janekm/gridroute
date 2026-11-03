//! Planar geometry in millimetres (f64) and exact predicates on integer nanometre coordinates.

use std::ops::{Add, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct P2 {
    pub x: f64,
    pub y: f64,
}

pub const fn p2(x: f64, y: f64) -> P2 {
    P2 { x, y }
}

impl Add for P2 {
    type Output = P2;
    fn add(self, o: P2) -> P2 {
        p2(self.x + o.x, self.y + o.y)
    }
}
impl Sub for P2 {
    type Output = P2;
    fn sub(self, o: P2) -> P2 {
        p2(self.x - o.x, self.y - o.y)
    }
}
impl Mul<f64> for P2 {
    type Output = P2;
    fn mul(self, k: f64) -> P2 {
        p2(self.x * k, self.y * k)
    }
}
impl Neg for P2 {
    type Output = P2;
    fn neg(self) -> P2 {
        p2(-self.x, -self.y)
    }
}

impl P2 {
    pub fn dot(self, o: P2) -> f64 {
        self.x * o.x + self.y * o.y
    }
    pub fn cross(self, o: P2) -> f64 {
        self.x * o.y - self.y * o.x
    }
    #[inline]
    pub fn norm(self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }
    pub fn norm2(self) -> f64 {
        self.x * self.x + self.y * self.y
    }
    pub fn dist(self, o: P2) -> f64 {
        (self - o).norm()
    }
    /// Left normal (rotated +90 degrees).
    pub fn perp(self) -> P2 {
        p2(-self.y, self.x)
    }
    pub fn unit(self) -> P2 {
        let n = self.norm();
        if n > 0.0 {
            self * (1.0 / n)
        } else {
            p2(1.0, 0.0)
        }
    }
    pub fn lerp(self, o: P2, t: f64) -> P2 {
        self + (o - self) * t
    }
}

/// Nanometre integer coordinates for exact predicates.
pub type IP = [i64; 2];

pub const NM: f64 = 1e6;

pub fn to_ip(p: P2) -> IP {
    [(p.x * NM).round() as i64, (p.y * NM).round() as i64]
}
pub fn from_ip(p: IP) -> P2 {
    p2(p[0] as f64 / NM, p[1] as f64 / NM)
}

/// Sign of the cross product (b-a) x (c-a): >0 when c is left of a->b.
pub fn orient(a: IP, b: IP, c: IP) -> i32 {
    let v = (b[0] - a[0]) as i128 * (c[1] - a[1]) as i128 - (b[1] - a[1]) as i128 * (c[0] - a[0]) as i128;
    v.signum() as i32
}

/// >0 when d lies strictly inside the circumcircle of the counter-clockwise triangle (a, b, c).
pub fn incircle(a: IP, b: IP, c: IP, d: IP) -> i32 {
    let (adx, ady) = ((a[0] - d[0]) as i128, (a[1] - d[1]) as i128);
    let (bdx, bdy) = ((b[0] - d[0]) as i128, (b[1] - d[1]) as i128);
    let (cdx, cdy) = ((c[0] - d[0]) as i128, (c[1] - d[1]) as i128);
    let alift = adx * adx + ady * ady;
    let blift = bdx * bdx + bdy * bdy;
    let clift = cdx * cdx + cdy * cdy;
    let det = alift * (bdx * cdy - bdy * cdx) - blift * (adx * cdy - ady * cdx) + clift * (adx * bdy - ady * bdx);
    det.signum() as i32
}

/// Distance from p to the boundary of polygon `poly` (closed ring).
pub fn poly_edge_dist(poly: &[P2], p: P2) -> f64 {
    let n = poly.len();
    (0..n).map(|i| seg_point_dist(poly[i], poly[(i + 1) % n], p)).fold(f64::INFINITY, f64::min)
}

pub fn seg_point_dist(a: P2, b: P2, p: P2) -> f64 {
    let d = b - a;
    let l2 = d.norm2();
    let t = if l2 > 0.0 { ((p - a).dot(d) / l2).clamp(0.0, 1.0) } else { 0.0 };
    (a + d * t).dist(p)
}

/// Closest parameter on segment a-b to p (clamped).
pub fn seg_param(a: P2, b: P2, p: P2) -> f64 {
    let d = b - a;
    let l2 = d.norm2();
    if l2 > 0.0 {
        ((p - a).dot(d) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    }
}

fn orient_f(a: P2, b: P2, c: P2) -> f64 {
    (b - a).cross(c - a)
}

/// Proper or touching intersection of segments a-b and c-d.
pub fn segs_intersect(a: P2, b: P2, c: P2, d: P2) -> bool {
    let d1 = orient_f(a, b, c);
    let d2 = orient_f(a, b, d);
    let d3 = orient_f(c, d, a);
    let d4 = orient_f(c, d, b);
    ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0)) && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0))
}

pub fn seg_seg_dist(a: P2, b: P2, c: P2, d: P2) -> f64 {
    if segs_intersect(a, b, c, d) {
        return 0.0;
    }
    seg_point_dist(a, b, c)
        .min(seg_point_dist(a, b, d))
        .min(seg_point_dist(c, d, a))
        .min(seg_point_dist(c, d, b))
}

/// Parameter t on a-b and u on c-d where the segments cross, if they do (non-parallel).
pub fn seg_cross_params(a: P2, b: P2, c: P2, d: P2) -> Option<(f64, f64)> {
    let r = b - a;
    let s = d - c;
    let den = r.cross(s);
    if den.abs() < 1e-18 {
        return None;
    }
    let t = (c - a).cross(s) / den;
    let u = (c - a).cross(r) / den;
    if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
        Some((t, u))
    } else {
        None
    }
}

/// Distance from p to a convex polygon (counter-clockwise or clockwise points); 0 inside.
/// One point is a disc centre, two points a segment.
pub fn point_hull_dist(pts: &[P2], p: P2) -> f64 {
    match pts.len() {
        0 => f64::INFINITY,
        1 => pts[0].dist(p),
        2 => seg_point_dist(pts[0], pts[1], p),
        n => {
            let mut inside = true;
            let mut sign = 0.0;
            let mut best = f64::INFINITY;
            for i in 0..n {
                let a = pts[i];
                let b = pts[(i + 1) % n];
                let o = orient_f(a, b, p);
                if sign == 0.0 && o != 0.0 {
                    sign = o.signum();
                }
                if o * sign < 0.0 {
                    inside = false;
                }
                best = best.min(seg_point_dist(a, b, p));
            }
            if inside {
                0.0
            } else {
                best
            }
        }
    }
}

pub fn point_in_hull(pts: &[P2], p: P2) -> bool {
    pts.len() >= 3 && point_hull_dist(pts, p) == 0.0
}

/// Closest point of a convex hull (or segment / point) to p.
pub fn hull_closest(pts: &[P2], p: P2) -> P2 {
    match pts.len() {
        0 => p,
        1 => pts[0],
        2 => pts[0].lerp(pts[1], seg_param(pts[0], pts[1], p)),
        n => {
            if point_in_hull(pts, p) {
                return p;
            }
            let mut best = pts[0];
            let mut bd = f64::INFINITY;
            for i in 0..n {
                let (a, b) = (pts[i], pts[(i + 1) % n]);
                let q = a.lerp(b, seg_param(a, b, p));
                let d = q.dist(p);
                if d < bd {
                    bd = d;
                    best = q;
                }
            }
            best
        }
    }
}

/// Distance from segment a-b to a convex hull (as for point_hull_dist).
pub fn seg_hull_dist(pts: &[P2], a: P2, b: P2) -> f64 {
    match pts.len() {
        0 => f64::INFINITY,
        1 => seg_point_dist(a, b, pts[0]),
        2 => seg_seg_dist(a, b, pts[0], pts[1]),
        n => {
            if point_in_hull(pts, a) || point_in_hull(pts, b) {
                return 0.0;
            }
            let mut best = f64::INFINITY;
            for i in 0..n {
                best = best.min(seg_seg_dist(a, b, pts[i], pts[(i + 1) % n]));
            }
            best
        }
    }
}

/// Point-in-polygon (even-odd) for simple polygons.
pub fn point_in_poly(poly: &[P2], p: P2) -> bool {
    let mut inside = false;
    let n = poly.len();
    let mut j = n.wrapping_sub(1);
    for i in 0..n {
        let (a, b) = (poly[i], poly[j]);
        if (a.y > p.y) != (b.y > p.y) && p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x {
            inside = !inside;
        }
        j = i;
    }
    inside
}

pub fn poly_area(poly: &[P2]) -> f64 {
    let n = poly.len();
    let mut s = 0.0;
    for i in 0..n {
        s += poly[i].cross(poly[(i + 1) % n]);
    }
    s * 0.5
}

/// Convex hull (counter-clockwise), removing duplicate and collinear points.
pub fn convex_hull(points: &[P2]) -> Vec<P2> {
    let mut pts: Vec<P2> = points.to_vec();
    pts.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap().then(a.y.partial_cmp(&b.y).unwrap()));
    pts.dedup_by(|a, b| (a.x - b.x).abs() < 1e-9 && (a.y - b.y).abs() < 1e-9);
    if pts.len() < 3 {
        return pts;
    }
    let mut lower: Vec<P2> = Vec::new();
    for &p in &pts {
        while lower.len() >= 2 && orient_f(lower[lower.len() - 2], lower[lower.len() - 1], p) <= 1e-12 {
            lower.pop();
        }
        lower.push(p);
    }
    let mut upper: Vec<P2> = Vec::new();
    for &p in pts.iter().rev() {
        while upper.len() >= 2 && orient_f(upper[upper.len() - 2], upper[upper.len() - 1], p) <= 1e-12 {
            upper.pop();
        }
        upper.push(p);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

/// Unit direction of the common tangent leaving circle (c0, r0) and arriving at circle (c1, r1).
/// `s0`/`s1` are +1 when the circle lies on the left of the travel direction, -1 on the right.
/// Returns (direction, departure point, arrival point), or None if the circles leave no such tangent.
pub fn tangent(c0: P2, r0: f64, s0: f64, c1: P2, r1: f64, s1: f64) -> Option<(P2, P2, P2)> {
    let d = c1 - c0;
    let dd = d.norm2();
    let k = s1 * r1 - s0 * r0;
    let t2 = dd - k * k;
    if t2 <= -1e-7 * dd.max(1e-6) || dd < 1e-18 {
        return None;
    }
    let t = t2.max(0.0).sqrt();
    let u = (d * t - d.perp() * k) * (1.0 / dd);
    let n = u.perp();
    let p0 = c0 - n * (s0 * r0);
    let p1 = c1 - n * (s1 * r1);
    Some((u, p0, p1))
}

/// Polyline offset by d to the left (negative: right), mitred at vertices.
pub fn offset_polyline(pts: &[P2], d: f64) -> Vec<P2> {
    let n = pts.len();
    if n < 2 {
        return pts.to_vec();
    }
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let a = if i == 0 { pts[0] } else { pts[i - 1] };
        let b = if i + 1 == n { pts[n - 1] } else { pts[i + 1] };
        let (d0, d1) = (if i == 0 { (pts[1] - pts[0]).unit() } else { (pts[i] - a).unit() }, if i + 1 == n { (pts[i] - pts[i - 1]).unit() } else { (b - pts[i]).unit() });
        let t = (d0 + d1).unit();
        let nrm = t.perp();
        let c = nrm.dot(d0.perp()).max(0.5);
        out.push(pts[i] + nrm * (d / c));
    }
    out
}

/// Angle of v in (-pi, pi].
pub fn ang(v: P2) -> f64 {
    v.y.atan2(v.x)
}

/// Discretize an arc around c of radius r from point a to point b, turning left (ccw) when s>0.
/// Vertices are placed slightly outside the circle so chords never come nearer than r.
pub fn arc_points(c: P2, r: f64, a: P2, b: P2, s: f64, tol: f64, out: &mut Vec<P2>) {
    let a0 = ang(a - c);
    let a1 = ang(b - c);
    let mut span = if s > 0.0 { a1 - a0 } else { a0 - a1 };
    while span < 0.0 {
        span += std::f64::consts::TAU;
    }
    while span >= std::f64::consts::TAU {
        span -= std::f64::consts::TAU;
    }
    if r <= 1e-9 || span < 1e-9 {
        return;
    }
    let ro = r + tol;
    let step = 2.0 * (r / ro).acos();
    let n = ((span / step).ceil() as usize).max(1);
    let dir = if s > 0.0 { 1.0 } else { -1.0 };
    for k in 1..n {
        let t = a0 + dir * span * (k as f64) / (n as f64);
        out.push(c + p2(t.cos(), t.sin()) * ro);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tangent_outer_and_inner() {
        // equal circles, same side: parallel line offset by r
        let (u, p0, p1) = tangent(p2(0., 0.), 1.0, 1.0, p2(10., 0.), 1.0, 1.0).unwrap();
        assert!((u.x - 1.0).abs() < 1e-12 && u.y.abs() < 1e-12);
        assert!((p0.y + 1.0).abs() < 1e-12 && (p1.y + 1.0).abs() < 1e-12);
        // opposite sides: crossing tangent through midpoint
        let (_, p0, p1) = tangent(p2(0., 0.), 1.0, 1.0, p2(10., 0.), 1.0, -1.0).unwrap();
        let m = p0.lerp(p1, 0.5);
        assert!(m.dist(p2(5., 0.)) < 1e-9);
        assert!((p0.dist(p2(0., 0.)) - 1.0).abs() < 1e-12);
        assert!((p1.dist(p2(10., 0.)) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn predicates() {
        assert_eq!(orient([0, 0], [10, 0], [5, 3]), 1);
        assert_eq!(orient([0, 0], [10, 0], [5, -3]), -1);
        assert_eq!(orient([0, 0], [10, 0], [5, 0]), 0);
        assert_eq!(incircle([0, 0], [10, 0], [0, 10], [5, 5]), 1);
        assert_eq!(incircle([0, 0], [10, 0], [0, 10], [20, 20]), -1);
        assert_eq!(incircle([0, 0], [10, 0], [0, 10], [10, 10]), 0);
    }
}
