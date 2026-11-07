//! Non-uniform rectilinear (Yee) mesh: node lines per axis, in metres.
//!
//! A 1-D line set is generated from fixed lines (which always become mesh lines: interfaces, region edges)
//! and a target cell size that is `h` inside refinement regions and grows geometrically (`ratio` per cell)
//! with distance from them, capped at `hmax`.

#[derive(Clone, Debug)]
pub struct Grid {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
}

impl Grid {
    pub fn nx(&self) -> usize {
        self.x.len()
    }
    pub fn ny(&self) -> usize {
        self.y.len()
    }
    pub fn nz(&self) -> usize {
        self.z.len()
    }
    pub fn n(&self) -> usize {
        self.nx() * self.ny() * self.nz()
    }
    #[inline]
    pub fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        (k * self.ny() + j) * self.nx() + i
    }
    pub fn axis(&self, a: usize) -> &[f64] {
        match a {
            0 => &self.x,
            1 => &self.y,
            _ => &self.z,
        }
    }
    pub fn min_step(&self, a: usize) -> f64 {
        self.axis(a).windows(2).map(|w| w[1] - w[0]).fold(f64::INFINITY, f64::min)
    }
}

/// A refinement region `[a, b]` with cell size `h`.
#[derive(Clone, Copy, Debug)]
pub struct Refine {
    pub a: f64,
    pub b: f64,
    pub h: f64,
}

/// Lines covering `[lo, hi]`, through every `fixed` line inside it. `min_cells` optionally forces a
/// minimum cell count in given intervals (e.g. at least two cells across each dielectric layer).
pub fn lines_1d(lo: f64, hi: f64, fixed: &[f64], refine: &[Refine], ratio: f64, hmax: f64, min_cells: &[(f64, f64, usize)]) -> Vec<f64> {
    let size = |x: f64| -> f64 {
        let mut s = hmax;
        for r in refine {
            let d = if x < r.a { r.a - x } else if x > r.b { x - r.b } else { 0.0 };
            // geometric growth: after n cells beyond the region the size is h * ratio^n, i.e. at distance d
            // the size is about h + (ratio - 1) * d
            s = s.min(r.h + (ratio - 1.0) * d);
        }
        s.max(1e-12)
    };
    let mut fx: Vec<f64> = fixed.iter().copied().filter(|&f| f > lo && f < hi).collect();
    fx.push(lo);
    fx.push(hi);
    fx.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let tol = 1e-9 * (hi - lo).max(1e-6);
    fx.dedup_by(|a, b| (*a - *b).abs() < tol);
    let mut out = vec![fx[0]];
    for w in fx.windows(2) {
        let (a, b) = (w[0], w[1]);
        // cumulative integral of 1/size over [a, b]
        const M: usize = 4000;
        let dx = (b - a) / M as f64;
        let mut cum = vec![0.0f64; M + 1];
        for s in 0..M {
            let xm = a + (s as f64 + 0.5) * dx;
            cum[s + 1] = cum[s] + dx / size(xm);
        }
        let mut n = cum[M].ceil().max(1.0) as usize;
        for &(c0, c1, m) in min_cells {
            if a >= c0 - tol && b <= c1 + tol {
                // the interval lies inside a min-cells range: give it its share
                let share = ((b - a) / (c1 - c0) * m as f64).ceil() as usize;
                n = n.max(share.max(1));
            }
        }
        let tot = cum[M];
        let mut s = 0usize;
        for q in 1..n {
            let target = tot * q as f64 / n as f64;
            while s < M && cum[s + 1] < target {
                s += 1;
            }
            let f = if cum[s + 1] > cum[s] { (target - cum[s]) / (cum[s + 1] - cum[s]) } else { 0.0 };
            out.push(a + (s as f64 + f) * dx);
        }
        out.push(b);
    }
    out
}

/// Index of the line nearest `v`.
pub fn snap(lines: &[f64], v: f64) -> usize {
    match lines.binary_search_by(|p| p.partial_cmp(&v).unwrap()) {
        Ok(i) => i,
        Err(i) => {
            if i == 0 {
                0
            } else if i >= lines.len() {
                lines.len() - 1
            } else if v - lines[i - 1] < lines[i] - v {
                i - 1
            } else {
                i
            }
        }
    }
}

/// Indices `i` with `lo <= lines[i] <= hi` (a half-open range `[first, last+1)`).
pub fn range_in(lines: &[f64], lo: f64, hi: f64) -> std::ops::Range<usize> {
    let a = lines.partition_point(|&p| p < lo);
    let b = lines.partition_point(|&p| p <= hi);
    a..b.max(a)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn graded_lines_hit_fixed_and_grow() {
        let l = lines_1d(-10.0, 10.0, &[0.123], &[Refine { a: -1.0, b: 1.0, h: 0.05 }], 1.3, 1.0, &[]);
        assert!(l.iter().any(|&x| (x - 0.123).abs() < 1e-12));
        assert_eq!(l[0], -10.0);
        assert_eq!(*l.last().unwrap(), 10.0);
        let d: Vec<f64> = l.windows(2).map(|w| w[1] - w[0]).collect();
        assert!(d.iter().all(|&s| s > 0.0 && s <= 1.0 + 1e-9));
        // fine inside the region
        let mid: Vec<f64> = l.windows(2).filter(|w| w[0] > -0.9 && w[1] < 0.9).map(|w| w[1] - w[0]).collect();
        assert!(mid.iter().all(|&s| s <= 0.0501), "{mid:?}");
        // neighbouring cells never jump by much more than the ratio
        for w in d.windows(2) {
            let r = w[1] / w[0];
            assert!(r < 1.6 && r > 1.0 / 1.6, "jump {r}");
        }
    }

    #[test]
    fn min_cells_per_layer() {
        let l = lines_1d(0.0, 1.0, &[0.1, 0.2], &[], 1.3, 1.0, &[(0.1, 0.2, 3)]);
        let n = l.iter().filter(|&&x| x > 0.1 + 1e-12 && x < 0.2 - 1e-12).count();
        assert_eq!(n, 2);
    }
}
