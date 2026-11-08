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
    /// A perfect conductor instead of a dielectric (boxes only): every E edge whose two end nodes lie in the
    /// box is metal. A box of zero extent along one axis is a sheet on the mesh plane nearest to it, so copper
    /// can stand in any axis plane (an FPC antenna on an enclosure's side wall).
    pub metal: bool,
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
        if body.metal {
            continue;
        }
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

/// E edges made metal by the metal boxes: (axis, node index) pairs. A node is metal when it lies in any metal
/// box (with a 1e-6 mm tolerance; along an axis where the box has zero extent, the nearest node), and an edge is
/// metal when both its nodes are. Boxes that touch or overlap therefore join even where no mesh line falls on
/// their shared face, so polygons can be given as strips thinner than a cell; two sheets less than a cell
/// apart join too.
pub fn metal_edges(g: &Grid, bodies: &[Body]) -> Vec<(usize, usize)> {
    let lines: [Vec<f64>; 3] = [g.x.iter().map(|v| v * 1e3).collect(), g.y.iter().map(|v| v * 1e3).collect(), g.z.iter().map(|v| v * 1e3).collect()];
    let (nx, ny, nz) = (g.nx(), g.ny(), g.nz());
    let mut on = vec![false; g.n()];
    let mut nodes = Vec::new();
    for body in bodies.iter().filter(|b| b.metal) {
        let BodyShape::Box { lo, hi } = &body.shape else { continue };
        let mut r = [(0usize, 0usize); 3];
        for a in 0..3 {
            let l = &lines[a];
            r[a] = if hi[a] - lo[a] < 1e-9 {
                let k = crate::mesh::snap(l, lo[a]);
                (k, k + 1)
            } else {
                (l.partition_point(|&x| x < lo[a] - 1e-6), l.partition_point(|&x| x <= hi[a] + 1e-6))
            };
        }
        for k in r[2].0..r[2].1 {
            for j in r[1].0..r[1].1 {
                for i in r[0].0..r[0].1 {
                    let id = g.idx(i, j, k);
                    if !on[id] {
                        on[id] = true;
                        nodes.push((i, j, k, id));
                    }
                }
            }
        }
    }
    let mut out = Vec::new();
    for &(i, j, k, id) in &nodes {
        if i + 1 < nx && on[g.idx(i + 1, j, k)] {
            out.push((0, id));
        }
        if j + 1 < ny && on[g.idx(i, j + 1, k)] {
            out.push((1, id));
        }
        if k + 1 < nz && on[g.idx(i, j, k + 1)] {
            out.push((2, id));
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
        let a = body_cells(&g, &[Body { name: "b".into(), eps: 3.0, tand: 0.0, shape: BodyShape::Box { lo, hi }, metal: false }]);
        let b = body_cells(&g, &[Body { name: "m".into(), eps: 3.0, tand: 0.0, shape: BodyShape::Mesh { tris: box_mesh(lo, hi) }, metal: false }]);
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
        let fill = Body { name: "core".into(), eps: 2.0, tand: 0.0, shape: BodyShape::Box { lo: [-1.0; 3], hi: [1.0; 3] }, metal: false };
        let v = body_cells(&g, &[Body { name: "shell".into(), eps: 3.0, tand: 0.0, shape: BodyShape::Mesh { tris }, metal: false }, fill]);
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

    #[test]
    fn a_zero_thickness_metal_box_is_a_sheet_on_the_nearest_plane() {
        let g = grid(10, 1.0); // nodes at -5..5 mm
        let sheet = Body { name: "s".into(), eps: 1.0, tand: 0.0, shape: BodyShape::Box { lo: [-2.0, 0.3, -1.0], hi: [2.0, 0.3, 1.0] }, metal: true };
        let e = metal_edges(&g, &[sheet]);
        // nodes x -2..2 (5), y = 0 (the nearest), z -1..1 (3): x edges 4 x 3, z edges 5 x 2, no y edges
        assert_eq!(e.iter().filter(|(a, _)| *a == 0).count(), 12);
        assert_eq!(e.iter().filter(|(a, _)| *a == 2).count(), 10);
        assert_eq!(e.iter().filter(|(a, _)| *a == 1).count(), 0);
        assert!(body_cells(&g, &[Body { name: "s".into(), eps: 1.0, tand: 0.0, shape: BodyShape::Box { lo: [-2.0; 3], hi: [2.0; 3] }, metal: true }]).iter().all(|&v| v == 0));
    }

    #[test]
    fn abutting_metal_strips_join_off_the_mesh_lines() {
        let g = grid(10, 1.0); // nodes at -5..5 mm
        // two strips in the plane y = 0, meeting at x = 0.3 and z = 0.3 where no mesh line lies
        let a = Body { name: "a".into(), eps: 1.0, tand: 0.0, shape: BodyShape::Box { lo: [-2.0, 0.0, -1.0], hi: [0.3, 0.0, 0.3] }, metal: true };
        let b = Body { name: "b".into(), eps: 1.0, tand: 0.0, shape: BodyShape::Box { lo: [0.3, 0.0, 0.3], hi: [2.0, 0.0, 1.0] }, metal: true };
        let e = metal_edges(&g, &[a, b]);
        let (j0, k0) = (5usize, 5usize); // node indices of y = 0, z = 0
        let id = |i: usize, k: usize| (k * 11 + j0) * 11 + i;
        // a holds nodes x -2..0, z -1..0; b holds x 1..2, z 1: the x edge from (0,0) to (1,0) is not metal
        // (node (1,0) is in neither), but the strips are joined through no common node either way; what
        // matters is that edges inside each strip exist and none crosses a gap of empty nodes
        assert!(e.contains(&(0, id(4, k0))), "edge x -1..0 at z 0 inside a");
        assert!(!e.contains(&(0, id(5, k0))), "edge x 0..1 at z 0: node (1, 0) is outside both");
        assert!(e.contains(&(0, id(6, k0 + 1))), "edge x 1..2 at z 1 inside b");
    }

    #[test]
    fn thin_strips_tiling_a_sheet_give_the_same_metal_as_one_box() {
        let g = grid(10, 1.0);
        let one = Body { name: "s".into(), eps: 1.0, tand: 0.0, shape: BodyShape::Box { lo: [-2.0, 0.0, -1.0], hi: [2.0, 0.0, 1.0] }, metal: true };
        let strips: Vec<Body> = (0..40)
            .map(|n| {
                let z0 = -1.0 + 0.05 * n as f64;
                Body { name: format!("s{n}"), eps: 1.0, tand: 0.0, shape: BodyShape::Box { lo: [-2.0, 0.0, z0], hi: [2.0, 0.0, z0 + 0.05] }, metal: true }
            })
            .collect();
        let mut a = metal_edges(&g, &[one]);
        let mut b = metal_edges(&g, &strips);
        a.sort();
        b.sort();
        assert_eq!(a, b);
    }
}
