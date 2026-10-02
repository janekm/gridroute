//! Obstacle dilation: which cells of a window lie within a structuring element of an obstacle.
//!
//! The structuring element is given as row spans: `span[di + rc]` is the half-width in columns of the element on
//! row offset `di` (`-1` leaves the row out). A disc of radius r gives `span[k] = floor(sqrt(r² - di²))`; the
//! caller computes the spans, so the exact rounding rules of its own raster stay with it and every backend gives
//! bit-identical masks.
//!
//! A source cell is an obstacle when its value is non-zero and not in the exclusion list (the router's own net),
//! or when the optional extra mask is set there. Cells outside the source grid are never obstacles.
//!
//! The CPU algorithm is separable: for every source row it first finds, per output column, the horizontal
//! distance to the nearest obstacle (saturating at 255); a cell is then blocked when any row offset `di` has that
//! distance within `span[di]`. For an OR over several layers the per-row distances are combined with `min` before
//! the vertical pass, because every layer uses the same element.

use rayon::prelude::*;

/// Several equally sized 2-D planes stored back to back (a `[nL][ny][nx]` array, or a single plane).
#[derive(Clone, Copy)]
pub struct Grid3<T> {
    pub ptr: *const T,
    pub nl: usize,
    pub ny: usize,
    pub nx: usize,
}
unsafe impl<T> Send for Grid3<T> {}
unsafe impl<T> Sync for Grid3<T> {}

impl<T: Copy> Grid3<T> {
    #[inline(always)]
    pub fn row(&self, l: usize, i: usize) -> &[T] {
        unsafe { std::slice::from_raw_parts(self.ptr.add((l * self.ny + i) * self.nx), self.nx) }
    }
}

#[derive(Clone, Copy)]
pub struct Window {
    pub i0: usize,
    pub j0: usize,
    pub h: usize,
    pub w: usize,
}

pub struct Element<'a> {
    pub span: &'a [i32],
}

impl Element<'_> {
    pub fn rc(&self) -> usize {
        self.span.len() / 2
    }
    pub fn wmax(&self) -> usize {
        self.span.iter().copied().max().unwrap_or(-1).max(0) as usize
    }
}

#[inline(always)]
fn is_obstacle(v: i16, excl: &[i16]) -> bool {
    v != 0 && !excl.contains(&v)
}

/// Per output column x of the window, the distance (columns) from source column j0 + x to the nearest obstacle in
/// source row `i` of layer `l`, saturated at 255; written to `hd` (len w). `acc_min`: combine with `min` into `hd`.
fn row_distance(
    src: &Grid3<i16>,
    extra: Option<&Grid3<u8>>,
    excl: &[i16],
    l: usize,
    i: usize,
    win: &Window,
    wmax: usize,
    hd: &mut [u8],
    obs: &mut Vec<u8>,
    acc_min: bool,
) {
    // column range of the source that can matter: [j0 - wmax, j0 + w - 1 + wmax] clipped
    let c0 = win.j0.saturating_sub(wmax);
    let c1 = (win.j0 + win.w - 1 + wmax).min(src.nx - 1);
    let row = src.row(l, i);
    let n = c1 - c0 + 1;
    obs.clear();
    obs.resize(n, 0);
    match extra {
        Some(e) => {
            let er = e.row(0, i);
            for k in 0..n {
                obs[k] = (is_obstacle(row[c0 + k], excl) || er[c0 + k] != 0) as u8;
            }
        }
        None => {
            for k in 0..n {
                obs[k] = is_obstacle(row[c0 + k], excl) as u8;
            }
        }
    }
    // forward / backward passes over the column range: distance to the last obstacle seen
    let mut left = vec![255u8; n];
    let mut d: u32 = 255;
    for k in 0..n {
        d = if obs[k] != 0 { 0 } else { (d + 1).min(255) };
        left[k] = d as u8;
    }
    d = 255;
    for k in (0..n).rev() {
        d = if obs[k] != 0 { 0 } else { (d + 1).min(255) };
        let v = (d as u8).min(left[k]);
        left[k] = v;
    }
    let off = win.j0 - c0;
    let src_row = &left[off..off + win.w];
    if acc_min {
        for (h, &v) in hd.iter_mut().zip(src_row) {
            *h = (*h).min(v);
        }
    } else {
        hd.copy_from_slice(src_row);
    }
}

/// CPU dilation. `layers`: source layers to test. `reduce_or`: one output plane, the OR over the layers;
/// otherwise one output plane per layer (in order). `extra` (single plane, ORed into the obstacles) applies to
/// the layers whose flag in `extra_on` is set. `out` holds `planes * h * w` bytes (1 = blocked).
pub fn dilate_cpu(
    src: &Grid3<i16>,
    layers: &[usize],
    excl: &[i16],
    extra: Option<&Grid3<u8>>,
    extra_on: &[bool],
    win: &Window,
    el: &Element,
    reduce_or: bool,
    out: &mut [u8],
) {
    let rc = el.rc();
    let wmax = el.wmax();
    assert!(wmax < 255, "structuring element too wide");
    // source rows touched: [i0 - rc, i0 + h - 1 + rc] clipped to the grid
    let r0 = win.i0 as isize - rc as isize;
    let nrows = win.h + 2 * rc;
    let hw = win.h * win.w;
    let groups: Vec<Vec<usize>> = if reduce_or {
        vec![layers.to_vec()]
    } else {
        layers.iter().map(|&l| vec![l]).collect()
    };
    for (g, ls) in groups.iter().enumerate() {
        // horizontal distances for every touched row (255 for rows outside the grid)
        let mut hd = vec![255u8; nrows * win.w];
        hd.par_chunks_mut(win.w).enumerate().for_each_init(Vec::new, |obs, (k, hrow)| {
            let i = r0 + k as isize;
            if i < 0 || i >= src.ny as isize {
                return;
            }
            for (n, &l) in ls.iter().enumerate() {
                let ex = if extra_on.get(l).copied().unwrap_or(false) { extra } else { None };
                row_distance(src, ex, excl, l, i as usize, win, wmax, hrow, obs, n > 0);
            }
        });
        let plane = &mut out[g * hw..(g + 1) * hw];
        plane.par_chunks_mut(win.w).enumerate().for_each(|(y, orow)| {
            orow.fill(0);
            for (k, &s) in el.span.iter().enumerate() {
                if s < 0 {
                    continue;
                }
                let s = s as u8;
                let hrow = &hd[(y + k) * win.w..(y + k + 1) * win.w];
                for (o, &h) in orow.iter_mut().zip(hrow) {
                    *o |= (h <= s) as u8;
                }
            }
        });
    }
}

/// Reference implementation (direct scan of the element), used by the tests.
pub fn dilate_naive(
    src: &Grid3<i16>,
    layers: &[usize],
    excl: &[i16],
    extra: Option<&Grid3<u8>>,
    extra_on: &[bool],
    win: &Window,
    el: &Element,
    reduce_or: bool,
    out: &mut [u8],
) {
    let rc = el.rc() as isize;
    let hw = win.h * win.w;
    out.fill(0);
    for (n, &l) in layers.iter().enumerate() {
        let g = if reduce_or { 0 } else { n };
        for y in 0..win.h {
            for x in 0..win.w {
                let mut hit = false;
                'el: for di in -rc..=rc {
                    let s = el.span[(di + rc) as usize];
                    if s < 0 {
                        continue;
                    }
                    let i = (win.i0 + y) as isize + di;
                    if i < 0 || i >= src.ny as isize {
                        continue;
                    }
                    for dj in -(s as isize)..=(s as isize) {
                        let j = (win.j0 + x) as isize + dj;
                        if j < 0 || j >= src.nx as isize {
                            continue;
                        }
                        let v = src.row(l, i as usize)[j as usize];
                        let e = extra_on.get(l).copied().unwrap_or(false)
                            && extra.map_or(false, |e| e.row(0, i as usize)[j as usize] != 0);
                        if is_obstacle(v, excl) || e {
                            hit = true;
                            break 'el;
                        }
                    }
                }
                if hit {
                    out[g * hw + y * win.w + x] = 1;
                }
            }
        }
    }
}
