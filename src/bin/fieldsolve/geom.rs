//! 2-D layout geometry shared by the importers (Gerber, Excellon) and the mesher: copper shapes in board
//! millimetres, with the inside test the rasteriser uses.

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct P2 {
    pub x: f64,
    pub y: f64,
}

impl P2 {
    pub fn new(x: f64, y: f64) -> P2 {
        P2 { x, y }
    }
    pub fn dist(self, o: P2) -> f64 {
        ((self.x - o.x).powi(2) + (self.y - o.y).powi(2)).sqrt()
    }
}

/// Axis-aligned box `[x0, x1] x [y0, y1]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BBox {
    pub x0: f64,
    pub y0: f64,
    pub x1: f64,
    pub y1: f64,
}

impl BBox {
    pub fn empty() -> BBox {
        BBox { x0: f64::INFINITY, y0: f64::INFINITY, x1: f64::NEG_INFINITY, y1: f64::NEG_INFINITY }
    }
    pub fn add(&mut self, p: P2) {
        self.x0 = self.x0.min(p.x);
        self.y0 = self.y0.min(p.y);
        self.x1 = self.x1.max(p.x);
        self.y1 = self.y1.max(p.y);
    }
    pub fn union(&mut self, o: &BBox) {
        self.x0 = self.x0.min(o.x0);
        self.y0 = self.y0.min(o.y0);
        self.x1 = self.x1.max(o.x1);
        self.y1 = self.y1.max(o.y1);
    }
    pub fn grow(&self, d: f64) -> BBox {
        BBox { x0: self.x0 - d, y0: self.y0 - d, x1: self.x1 + d, y1: self.y1 + d }
    }
    pub fn is_empty(&self) -> bool {
        self.x0 > self.x1 || self.y0 > self.y1
    }
    pub fn contains(&self, p: P2) -> bool {
        p.x >= self.x0 && p.x <= self.x1 && p.y >= self.y0 && p.y <= self.y1
    }
    pub fn overlaps(&self, o: &BBox) -> bool {
        self.x0 <= o.x1 && o.x0 <= self.x1 && self.y0 <= o.y1 && o.y0 <= self.y1
    }
}

/// A filled copper shape.
#[derive(Clone, Debug, PartialEq)]
pub enum Shape {
    /// Disc.
    Circle { c: P2, r: f64 },
    /// Stadium: all points within `r` of the segment `a`-`b` (a track drawn with a round aperture).
    Capsule { a: P2, b: P2, r: f64 },
    /// Simple polygon, implicitly closed, even-odd fill (arcs already linearised).
    Polygon(Vec<P2>),
}

fn seg_dist2(p: P2, a: P2, b: P2) -> f64 {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let l2 = dx * dx + dy * dy;
    let t = if l2 > 0.0 { (((p.x - a.x) * dx + (p.y - a.y) * dy) / l2).clamp(0.0, 1.0) } else { 0.0 };
    let (qx, qy) = (a.x + t * dx - p.x, a.y + t * dy - p.y);
    qx * qx + qy * qy
}

pub fn point_in_polygon(pts: &[P2], p: P2) -> bool {
    let mut inside = false;
    let n = pts.len();
    if n < 3 {
        return false;
    }
    let mut j = n - 1;
    for i in 0..n {
        let (a, b) = (pts[i], pts[j]);
        if (a.y > p.y) != (b.y > p.y) && p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x {
            inside = !inside;
        }
        j = i;
    }
    inside
}

impl Shape {
    pub fn contains(&self, p: P2) -> bool {
        match self {
            Shape::Circle { c, r } => (p.x - c.x).powi(2) + (p.y - c.y).powi(2) <= r * r,
            Shape::Capsule { a, b, r } => seg_dist2(p, *a, *b) <= r * r,
            Shape::Polygon(pts) => point_in_polygon(pts, p),
        }
    }
    pub fn bbox(&self) -> BBox {
        let mut b = BBox::empty();
        match self {
            Shape::Circle { c, r } => {
                b.add(*c);
                b = b.grow(*r);
            }
            Shape::Capsule { a, b: e, r } => {
                b.add(*a);
                b.add(*e);
                b = b.grow(*r);
            }
            Shape::Polygon(pts) => pts.iter().for_each(|&p| b.add(p)),
        }
        b
    }
}

/// One object of a copper layer, in file order. A clear (`dark == false`) object removes copper painted
/// before it (Gerber `%LPC%`).
#[derive(Clone, Debug)]
pub struct Prim {
    pub shape: Shape,
    pub dark: bool,
    /// X2 object attribute `.N` (net name), if present.
    pub net: Option<String>,
    /// X2 object attribute `.P` (reference designator, pin number), if present.
    pub pad: Option<(String, String)>,
}

/// A drilled hole. `plated` holes are via or pin barrels between copper layers `from`..=`to` (1-based Gerber
/// layer numbers as in the drill file's FileFunction; through holes span all layers).
#[derive(Clone, Debug)]
pub struct Hole {
    pub c: P2,
    pub d: f64,
    pub plated: bool,
    pub from: usize,
    pub to: usize,
    /// For slots: the far end of the slot (the hole is the capsule `c`-`end` of diameter `d`).
    pub end: Option<P2>,
}
