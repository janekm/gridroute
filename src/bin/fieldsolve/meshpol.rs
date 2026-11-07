//! Static (a-priori) mesh refinement on the rectilinear mesh, from the layout: `"mesh"` in the spec.
//!
//! ```json
//! "mesh": {"no_slivers": true, "z_grade": 3.0,
//!          "copper_bands": {"layers": ["F.Cu"], "cell": 0.05, "margin": 0.1}}
//! ```
//!
//! - `no_slivers`: fill the intervals between fixed lines with a rounding slack (`mesh::lines_1d_slack`), so
//!   that two fixed lines a little over a cell apart do not leave two half cells. The smallest cell sets the
//!   time step on the whole mesh.
//! - `z_grade` (a factor f): through the board, `cell_z` only within two cells of each copper plane and of
//!   each dielectric interface; in between cells grow (1.3 per cell) up to f · `cell_z`. Thick cores then
//!   carry a few graded cells instead of a uniform stack. The minimum cells per layer still holds.
//! - `copper_bands`: the copper window gets the spec's `cell` as a base, and the listed layers' copper
//!   edges inside it get bands of the finer `cell` (± `margin` mm). The edges are found on a supersampled
//!   raster of the copper: x bands where the copper changes along x (edges running in y), y bands where it
//!   changes along y. On a tensor-product mesh a band runs through the whole domain, so this helps where
//!   copper edges are sparse in x or in y.

use crate::geom::BBox;
use crate::mesh::Refine;
use crate::model::Sheet;
use serde_json::Value;

#[derive(Clone, Debug, Default)]
pub struct MeshPolicy {
    pub no_slivers: bool,
    pub z_grade: Option<f64>,
    pub bands: Option<Bands>,
}

#[derive(Clone, Debug)]
pub struct Bands {
    pub layers: Vec<String>,
    pub cell: f64,
    pub margin: f64,
}

impl MeshPolicy {
    pub fn from_spec(sp: &Value) -> MeshPolicy {
        let Some(m) = sp.get("mesh") else { return MeshPolicy::default() };
        let bands = m.get("copper_bands").map(|b| Bands {
            layers: b.get("layers").and_then(|v| v.as_array()).map(|a| a.iter().filter_map(|x| x.as_str().map(String::from)).collect()).unwrap_or_default(),
            cell: b.get("cell").and_then(|v| v.as_f64()).unwrap_or(0.05),
            margin: b.get("margin").and_then(|v| v.as_f64()).unwrap_or(0.1),
        });
        MeshPolicy {
            no_slivers: m.get("no_slivers").and_then(|v| v.as_bool()).unwrap_or(false),
            z_grade: m.get("z_grade").and_then(|v| v.as_f64().or_else(|| v.as_bool().filter(|&b| b).map(|_| 3.0))),
            bands,
        }
    }

    pub fn slack(&self) -> f64 {
        if self.no_slivers {
            0.35
        } else {
            0.0
        }
    }

    /// Refinement regions for copper edges in `win` (mm): (x regions, y regions).
    pub fn copper_bands(&self, sheets: &[Sheet], win: &BBox) -> (Vec<Refine>, Vec<Refine>) {
        let Some(b) = &self.bands else { return (vec![], vec![]) };
        let s = b.cell / 4.0;
        let nx = ((win.x1 - win.x0) / s).ceil() as usize + 1;
        let ny = ((win.y1 - win.y0) / s).ceil() as usize + 1;
        let xs: Vec<f64> = (0..nx).map(|i| win.x0 + (i as f64 + 0.5) * s).collect();
        let ys: Vec<f64> = (0..ny).map(|j| win.y0 + (j as f64 + 0.5) * s).collect();
        // 1-D edge positions, per layer: x where its copper changes along a row, y along a column
        let mut ex = vec![false; nx];
        let mut ey = vec![false; ny];
        for sh in sheets.iter().filter(|sh| b.layers.is_empty() || b.layers.iter().any(|l| l == &sh.name)) {
            let mut c = vec![false; nx * ny];
            for p in &sh.prims {
                if !p.shape.bbox().overlaps(win) {
                    continue;
                }
                let d = p.dark;
                crate::model::paint_with(&xs, &ys, p, |q| c[q] = d);
            }
            for j in 0..ny {
                for i in 0..nx {
                    let q = j * nx + i;
                    if i + 1 < nx && c[q] != c[q + 1] {
                        ex[i] = true;
                    }
                    if j + 1 < ny && c[q] != c[q + nx] {
                        ey[j] = true;
                    }
                }
            }
        }
        let regions = |e: &[bool], v: &[f64]| -> Vec<Refine> {
            let mut out: Vec<Refine> = Vec::new();
            for (k, &on) in e.iter().enumerate() {
                if !on {
                    continue;
                }
                let x = v[k] + 0.5 * s;
                let (a, bb) = (x - b.margin, x + b.margin);
                match out.last_mut() {
                    Some(r) if a <= r.b => r.b = r.b.max(bb),
                    _ => out.push(Refine { a, b: bb, h: b.cell }),
                }
            }
            out
        };
        (regions(&ex, &xs), regions(&ey, &ys))
    }

    /// z refinement: either the whole board at `cell_z`, or (z_grade) `cell_z` near every copper plane and
    /// dielectric interface with graded cells capped at z_grade · cell_z in between.
    pub fn z_refine(&self, zlo: f64, zhi: f64, planes: &[f64], cell_z: f64) -> Vec<Refine> {
        match self.z_grade {
            None => vec![Refine { a: zlo, b: zhi, h: cell_z }],
            Some(f) => {
                let mut v = vec![Refine { a: zlo, b: zhi, h: f * cell_z }];
                for &z in planes {
                    v.push(Refine { a: z - 2.0 * cell_z, b: z + 2.0 * cell_z, h: cell_z });
                }
                v
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slack_removes_slivers() {
        // fixed lines 1.02 cells apart: ceil gives two 0.51 cells, the slack one
        let a = crate::mesh::lines_1d(0.0, 1.0, &[0.3, 0.351], &[Refine { a: 0.0, b: 1.0, h: 0.05 }], 1.3, 1.0, &[]);
        let b = crate::mesh::lines_1d_slack(0.0, 1.0, &[0.3, 0.351], &[Refine { a: 0.0, b: 1.0, h: 0.05 }], 1.3, 1.0, &[], 0.35);
        let min = |l: &[f64]| l.windows(2).map(|w| w[1] - w[0]).fold(f64::INFINITY, f64::min);
        assert!(min(&a) < 0.03, "plain fill has the sliver ({})", min(&a));
        assert!(min(&b) > 0.04, "slack fill has none ({})", min(&b));
        let max = b.windows(2).map(|w| w[1] - w[0]).fold(0.0, f64::max);
        assert!(max <= 0.05 * 1.36, "cells stay within the slack ({max})");
    }
}
