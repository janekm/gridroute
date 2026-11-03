//! SVG rendering of the routing state, one panel per signal layer.

use crate::cdt::NONE;
use crate::geom::P2;
use crate::model::Kind;
use crate::world::{World, T_BLOCK};
use std::fmt::Write;

pub struct Opts {
    pub mesh: bool,
    pub box_: Option<[f64; 4]>,
    pub scale: f64,
    pub marks: Vec<(P2, String)>,
    pub unrouted: Vec<u32>,
}

const COLORS: [&str; 6] = ["#e8695e", "#6f86e8", "#5fbf6c", "#e2b347", "#b07ee0", "#4fc1d6"];

pub fn render(w: &World, o: &Opts) -> String {
    let (x0, y0, x1, y1) = match o.box_ {
        Some(b) => (b[0], b[1], b[2], b[3]),
        None => (-1.0, -1.0, w.b.size[0] + 1.0, w.b.size[1] + 1.0),
    };
    let nl = w.nsig();
    let cols = if nl > 2 { 2 } else { nl };
    let rows = (nl + cols - 1) / cols;
    let pw = (x1 - x0) * o.scale;
    let ph = (y1 - y0) * o.scale + 24.0;
    let mut s = String::new();
    let _ = write!(s, r#"<svg xmlns="http://www.w3.org/2000/svg" width="{:.0}" height="{:.0}" viewBox="0 0 {:.0} {:.0}">"#, pw * cols as f64, ph * rows as f64, pw * cols as f64, ph * rows as f64);
    let _ = write!(s, r##"<rect width="100%" height="100%" fill="#0d1621"/>"##);
    for l in 0..nl {
        let ox = (l % cols) as f64 * pw;
        let oy = (l / cols) as f64 * ph;
        let tx = |p: P2| ((p.x - x0) * o.scale + ox, (p.y - y0) * o.scale + oy + 24.0);
        let _ = write!(s, r##"<g><text x="{}" y="{}" fill="#e8f0f6" font-size="16" font-family="monospace">{}</text>"##, ox + 6.0, oy + 17.0, w.b.sig_name(l));
        let _ = write!(s, r##"<clipPath id="c{}"><rect x="{}" y="{}" width="{}" height="{}"/></clipPath><g clip-path="url(#c{})">"##, l, ox, oy + 24.0, pw, ph - 24.0, l);
        // board
        let mut d = String::new();
        for (i, q) in w.b.outline.iter().enumerate() {
            let (x, y) = tx(*q);
            let _ = write!(d, "{}{:.2},{:.2} ", if i == 0 { "M" } else { "L" }, x, y);
        }
        let _ = write!(s, r##"<path d="{}Z" fill="#182c34"/>"##, d);
        let layer = &w.layers[l];
        if o.mesh {
            let mut d = String::new();
            for t in 0..layer.cdt.tris.len() {
                let tr = &layer.cdt.tris[t];
                if !tr.alive || layer.tflag[t] & T_BLOCK != 0 {
                    continue;
                }
                let p: Vec<(f64, f64)> = tr.v.iter().map(|&v| tx(w.vpos(l, v))).collect();
                let _ = write!(d, "M{:.2},{:.2}L{:.2},{:.2}L{:.2},{:.2}Z", p[0].0, p[0].1, p[1].0, p[1].1, p[2].0, p[2].1);
            }
            let _ = write!(s, r##"<path d="{}" fill="none" stroke="#2f4a55" stroke-width="0.6"/>"##, d);
        }
        // objects
        for (oi, ob) in w.objs.iter().enumerate() {
            if !ob.alive || !ob.on(l) {
                continue;
            }
            let fill = match ob.s.kind {
                Kind::Pad => {
                    if ob.s.net == NONE {
                        "#6b6650"
                    } else if w.b.nets[ob.s.net as usize].plane.is_some() {
                        "#9fb48a"
                    } else {
                        "#e5d8a4"
                    }
                }
                Kind::Via => "#d0d0d0",
                Kind::Hole => "#444444",
                Kind::Keepout => "#5a2a2a",
                Kind::Outline => continue,
            };
            let mut d = String::new();
            for (i, q) in ob.s.pts.iter().enumerate() {
                let (x, y) = tx(*q);
                let _ = write!(d, "{}{:.2},{:.2} ", if i == 0 { "M" } else { "L" }, x, y);
            }
            let sw = (2.0 * ob.s.r * o.scale).max(0.0);
            let unr = o.unrouted.contains(&(oi as u32));
            let stroke = if unr { "#ff3030" } else { fill };
            let _ = write!(s, r##"<path d="{}Z" fill="{}" stroke="{}" stroke-width="{:.2}" stroke-linejoin="round" stroke-linecap="round" fill-opacity="0.95"/>"##, d, fill, stroke, sw.max(if unr { 2.0 } else { 0.0 }));
            if ob.s.kind == Kind::Via {
                let (x, y) = tx(ob.s.center);
                let _ = write!(s, r##"<circle cx="{:.2}" cy="{:.2}" r="{:.2}" fill="#0d1621"/>"##, x, y, ob.s.drill / 2.0 * o.scale);
            }
            if unr {
                let (x, y) = tx(ob.s.center);
                let _ = write!(s, r##"<circle cx="{:.2}" cy="{:.2}" r="{:.2}" fill="none" stroke="#ff3030" stroke-width="2"/>"##, x, y, 0.6 * o.scale);
            }
        }
        // pieces
        let col = COLORS[if l == 0 { 0 } else if l == nl - 1 { 1 } else { 2 + (l - 1) % 4 }];
        for pc in &w.pieces {
            if !pc.alive || pc.layer as usize != l || pc.poly.len() < 2 {
                continue;
            }
            let mut d = String::new();
            for (i, q) in pc.poly.iter().enumerate() {
                let (x, y) = tx(*q);
                let _ = write!(d, "{}{:.2},{:.2} ", if i == 0 { "M" } else { "L" }, x, y);
            }
            let c = if pc.overlap { "#ff00ff" } else { col };
            let _ = write!(s, r##"<path d="{}" fill="none" stroke="{}" stroke-width="{:.2}" stroke-linejoin="round" stroke-linecap="round" stroke-opacity="0.9"/>"##, d, c, (pc.width * o.scale).max(0.8));
        }
        for (p, label) in &o.marks {
            let (x, y) = tx(*p);
            let _ = write!(s, r##"<circle cx="{:.2}" cy="{:.2}" r="{:.2}" fill="none" stroke="#ff40ff" stroke-width="2"/><text x="{:.2}" y="{:.2}" fill="#ff80ff" font-size="10">{}</text>"##, x, y, 0.4 * o.scale, x + 6.0, y - 6.0, label);
        }
        s.push_str("</g></g>");
    }
    s.push_str("</svg>");
    s
}
