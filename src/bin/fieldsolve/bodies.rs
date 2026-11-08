//! 3-D dielectric objects (roadmap A2, first step): enclosure shells from STL, boxes such as an FPC antenna's
//! substrate. Each object is voxelised by cell centre (a cell belongs to the last object containing its
//! centre) and replaces the board's slab material there, with a constant conductivity from `tand` at the
//! build's `f_loss`. Closed triangle meshes are filled by vertical ray parity per cell column.

use crate::mesh::Grid;

#[derive(Clone, Debug)]
pub enum BodyShape {
    /// Axis-aligned box, mm: `lo` and `hi` corners (x, y, z).
    Box { lo: [f64; 3], hi: [f64; 3] },
    /// Closed (watertight) triangle mesh, mm.
    Mesh { tris: Vec<[[f64; 3]; 3]> },
}

#[derive(Clone, Debug)]
pub struct Body {
    pub name: String,
    pub eps: f64,
    pub tand: f64,
    pub shape: BodyShape,
}

impl Body {
    /// Bounding box corners, mm.
    pub fn bounds(&self) -> ([f64; 3], [f64; 3]) {
        match &self.shape {
            BodyShape::Box { lo, hi } => (*lo, *hi),
            BodyShape::Mesh { tris } => {
                let mut lo = [f64::INFINITY; 3];
                let mut hi = [f64::NEG_INFINITY; 3];
                for t in tris {
                    for v in t {
                        for a in 0..3 {
                            lo[a] = lo[a].min(v[a]);
                            hi[a] = hi[a].max(v[a]);
                        }
                    }
                }
                (lo, hi)
            }
        }
    }
}

/// Triangles of an STL file (binary or ASCII), mm, each vertex mapped by `f`.
pub fn read_stl(path: &std::path::Path, f: impl Fn([f64; 3]) -> [f64; 3]) -> Result<Vec<[[f64; 3]; 3]>, String> {
    let b = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut tris = Vec::new();
    let binary_len = if b.len() >= 84 { 84 + 50 * u32::from_le_bytes([b[80], b[81], b[82], b[83]]) as usize } else { 0 };
    if b.len() >= 84 && binary_len == b.len() {
        let n = (b.len() - 84) / 50;
        tris.reserve(n);
        for t in 0..n {
            let o = 84 + 50 * t + 12;
            let v = |k: usize| f32::from_le_bytes([b[o + 4 * k], b[o + 4 * k + 1], b[o + 4 * k + 2], b[o + 4 * k + 3]]) as f64;
            tris.push([f([v(0), v(1), v(2)]), f([v(3), v(4), v(5)]), f([v(6), v(7), v(8)])]);
        }
    } else {
        let txt = String::from_utf8_lossy(&b);
        let mut cur: Vec<[f64; 3]> = Vec::with_capacity(3);
        for line in txt.lines() {
            let mut it = line.split_whitespace();
            if it.next() == Some("vertex") {
                let p: Vec<f64> = it.filter_map(|s| s.parse().ok()).collect();
                if p.len() != 3 {
                    return Err(format!("{}: bad vertex line {line:?}", path.display()));
                }
                cur.push(f([p[0], p[1], p[2]]));
                if cur.len() == 3 {
                    tris.push([cur[0], cur[1], cur[2]]);
                    cur.clear();
                }
            }
        }
    }
    if tris.is_empty() {
        return Err(format!("{}: no triangles", path.display()));
    }
    Ok(tris)
}

/// Body of each cell (index into `bodies` plus one, 0 = none), cell (i, j, k) at `(k * (ny-1) + j) * (nx-1) + i`.
/// Later bodies win where they overlap.
pub fn body_cells(g: &Grid, bodies: &[Body]) -> Vec<u16> {
    let (nx, ny, nz) = (g.nx(), g.ny(), g.nz());
    let (cx, cy, cz) = (nx.saturating_sub(1), ny.saturating_sub(1), nz.saturating_sub(1));
    let centres = |v: &[f64]| -> Vec<f64> { v.windows(2).map(|w| 0.5e3 * (w[0] + w[1])).collect() };
    let (xc, yc, zc) = (centres(&g.x), centres(&g.y), centres(&g.z));
    let mut out = vec![0u16; cx * cy * cz];
    if bodies.is_empty() {
        return out;
    }
    // cells whose centre lies in [a, b]
    let span = |c: &[f64], a: f64, b: f64| -> (usize, usize) { (c.partition_point(|&x| x < a), c.partition_point(|&x| x <= b)) };
    for (bi, body) in bodies.iter().enumerate() {
        let tag = (bi + 1).min(u16::MAX as usize) as u16;
        match &body.shape {
            BodyShape::Box { lo, hi } => {
                let (i0, i1) = span(&xc, lo[0], hi[0]);
                let (j0, j1) = span(&yc, lo[1], hi[1]);
                let (k0, k1) = span(&zc, lo[2], hi[2]);
                for k in k0..k1 {
                    for j in j0..j1 {
                        let row = (k * cy + j) * cx;
                        out[row + i0..row + i1].fill(tag);
                    }
                }
            }
            BodyShape::Mesh { tris } => {
                // z of every crossing of each column's vertical centre line with the surface; the line is
                // nudged off the cell centre so that it never runs exactly through a shared edge or vertex
                let (ex, ey) = (1.234_567e-7, 7.654_321e-8);
                let mut hits: Vec<Vec<f64>> = vec![Vec::new(); cx * cy];
                for t in tris {
                    let (x0, x1) = (t[0][0].min(t[1][0]).min(t[2][0]), t[0][0].max(t[1][0]).max(t[2][0]));
                    let (y0, y1) = (t[0][1].min(t[1][1]).min(t[2][1]), t[0][1].max(t[1][1]).max(t[2][1]));
                    let (i0, i1) = span(&xc, x0 - ex, x1 - ex);
                    let (j0, j1) = span(&yc, y0 - ey, y1 - ey);
                    if i0 >= i1 || j0 >= j1 {
                        continue;
                    }
                    let (a, b, c) = (t[0], t[1], t[2]);
                    let det = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
                    if det.abs() < 1e-18 {
                        continue; // vertical facet: parallel to the ray
                    }
                    for j in j0..j1 {
                        let py = yc[j] + ey;
                        for i in i0..i1 {
                            let px = xc[i] + ex;
                            let u = ((px - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (py - a[1])) / det;
                            let v = ((b[0] - a[0]) * (py - a[1]) - (px - a[0]) * (b[1] - a[1])) / det;
                            if u >= 0.0 && v >= 0.0 && u + v <= 1.0 {
                                hits[j * cx + i].push(a[2] + u * (b[2] - a[2]) + v * (c[2] - a[2]));
                            }
                        }
                    }
                }
                for j in 0..cy {
                    for i in 0..cx {
                        let h = &mut hits[j * cx + i];
                        if h.len() < 2 {
                            continue;
                        }
                        h.sort_by(|p, q| p.partial_cmp(q).unwrap());
                        for pair in h.chunks_exact(2) {
                            let (k0, k1) = span(&zc, pair[0], pair[1]);
                            for k in k0..k1 {
                                out[(k * cy + j) * cx + i] = tag;
                            }
                        }
                    }
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid(n: usize, h: f64) -> Grid {
        let l: Vec<f64> = (0..=n).map(|i| (i as f64 * h - n as f64 * h / 2.0) * 1e-3).collect();
        Grid { x: l.clone(), y: l.clone(), z: l }
    }

    /// The 12 triangles of an axis-aligned box.
    fn box_mesh(lo: [f64; 3], hi: [f64; 3]) -> Vec<[[f64; 3]; 3]> {
        let p = |i: usize| [if i & 1 > 0 { hi[0] } else { lo[0] }, if i & 2 > 0 { hi[1] } else { lo[1] }, if i & 4 > 0 { hi[2] } else { lo[2] }];
        let quads = [[0, 1, 3, 2], [4, 6, 7, 5], [0, 4, 5, 1], [2, 3, 7, 6], [0, 2, 6, 4], [1, 5, 7, 3]];
        quads.iter().flat_map(|q| [[p(q[0]), p(q[1]), p(q[2])], [p(q[0]), p(q[2]), p(q[3])]]).collect()
    }

    #[test]
    fn mesh_box_fills_the_same_cells_as_the_box() {
        let g = grid(20, 0.5);
        let (lo, hi) = ([-2.1, -1.3, -0.7], [1.6, 2.4, 3.1]);
        let a = body_cells(&g, &[Body { name: "b".into(), eps: 3.0, tand: 0.0, shape: BodyShape::Box { lo, hi } }]);
        let b = body_cells(&g, &[Body { name: "m".into(), eps: 3.0, tand: 0.0, shape: BodyShape::Mesh { tris: box_mesh(lo, hi) } }]);
        let n = a.iter().filter(|&&v| v > 0).count();
        assert!(n > 0);
        assert_eq!(a, b, "ray parity must reproduce the box");
    }

    #[test]
    fn hollow_shell_leaves_its_inside_empty_and_later_bodies_win() {
        let g = grid(24, 0.5);
        // a closed shell: outer box minus inner box (inner faces reversed do not matter for parity)
        let mut tris = box_mesh([-5.0, -5.0, -5.0], [5.0, 5.0, 5.0]);
        tris.extend(box_mesh([-3.0, -3.0, -3.0], [3.0, 3.0, 3.0]));
        let fill = Body { name: "core".into(), eps: 2.0, tand: 0.0, shape: BodyShape::Box { lo: [-1.0; 3], hi: [1.0; 3] } };
        let v = body_cells(&g, &[Body { name: "shell".into(), eps: 3.0, tand: 0.0, shape: BodyShape::Mesh { tris } }, fill]);
        let (cx, cy) = (24, 24);
        let at = |x: f64, y: f64, z: f64| {
            let i = ((x + 6.0) / 0.5) as usize;
            let j = ((y + 6.0) / 0.5) as usize;
            let k = ((z + 6.0) / 0.5) as usize;
            v[(k * cy + j) * cx + i]
        };
        assert_eq!(at(4.1, 0.1, 0.1), 1, "in the wall");
        assert_eq!(at(2.1, 0.1, 0.1), 0, "in the cavity");
        assert_eq!(at(0.1, 0.1, 0.1), 2, "the later body wins");
        assert_eq!(at(5.6, 0.1, 0.1), 0, "outside");
    }
}
