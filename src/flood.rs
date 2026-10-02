//! Copper islands: the cells of one net connected to seed cells, within a window.
//!
//! Connectivity is 4-neighbour within a layer, plus a jump between all layers at cells where the through-hole
//! plane (vias, plated holes) carries the net.

use crate::dilate::{Grid3, Window};

/// `core`: [nl][ny][nx] copper net ids, `thru`: [ny][nx] hole net ids. `seeds`: (layer, y, x) relative to the
/// window; seeds outside the window or not on the net's copper are ignored. `out`: nl * h * w bytes.
pub fn island(core: &Grid3<i16>, thru: &Grid3<i16>, nid: i16, win: &Window, seeds: &[[i32; 3]], out: &mut [u8]) {
    island_crop(core, thru, nid, win, seeds, win, out);
}

/// `island` traced over `win` but written only for the sub-window `crop` (out: nl * crop.h * crop.w bytes). The
/// trace keeps its own visited set, so the large trace window costs only the net's copper, not a full mask.
pub fn island_crop(core: &Grid3<i16>, thru: &Grid3<i16>, nid: i16, win: &Window, seeds: &[[i32; 3]], crop: &Window,
                   out: &mut [u8]) {
    let (h, w, nl) = (win.h, win.w, core.nl);
    let hw = h * w;
    out.fill(0);
    let mut seen = std::collections::HashSet::<u32, std::hash::BuildHasherDefault<Fx>>::default();
    let mut q: Vec<u32> = Vec::new();
    let own = |l: usize, y: usize, x: usize| core.row(l, win.i0 + y)[win.j0 + x] == nid;
    let mut visit = |k: usize, q: &mut Vec<u32>| {
        if seen.insert(k as u32) {
            q.push(k as u32);
        }
    };
    for s in seeds {
        let (l, y, x) = (s[0], s[1], s[2]);
        if l < 0 || l as usize >= nl || y < 0 || x < 0 || y as usize >= h || x as usize >= w {
            continue;
        }
        let (l, y, x) = (l as usize, y as usize, x as usize);
        if own(l, y, x) {
            visit(l * hw + y * w + x, &mut q);
        }
    }
    let mut head = 0;
    while head < q.len() {
        let s = q[head] as usize;
        head += 1;
        let l = s / hw;
        let c = s - l * hw;
        let y = c / w;
        let x = c - y * w;
        let (gy, gx) = (win.i0 + y, win.j0 + x);
        if gy >= crop.i0 && gx >= crop.j0 && gy < crop.i0 + crop.h && gx < crop.j0 + crop.w {
            out[(l * crop.h + gy - crop.i0) * crop.w + gx - crop.j0] = 1;
        }
        if x > 0 && own(l, y, x - 1) {
            visit(s - 1, &mut q);
        }
        if x + 1 < w && own(l, y, x + 1) {
            visit(s + 1, &mut q);
        }
        if y > 0 && own(l, y - 1, x) {
            visit(s - w, &mut q);
        }
        if y + 1 < h && own(l, y + 1, x) {
            visit(s + w, &mut q);
        }
        if thru.row(0, gy)[gx] == nid {
            for l2 in 0..nl {
                if l2 != l && own(l2, y, x) {
                    visit(l2 * hw + c, &mut q);
                }
            }
        }
    }
}

/// Multiplicative hasher for the u32 cell keys of the visited set.
#[derive(Default)]
pub struct Fx(u64);

impl std::hash::Hasher for Fx {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, b: &[u8]) {
        for &x in b {
            self.0 = (self.0.rotate_left(5) ^ x as u64).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
        }
    }
    fn write_u32(&mut self, x: u32) {
        self.0 = (self.0.rotate_left(5) ^ x as u64).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
    }
}

/// Copper components of net `nid` at query cells: label[k] = the component of cell k (all cells reached from the
/// same seed share a label), -1 where the cell is not the net's copper. One pass over the net's copper reachable
/// from the queries (4-neighbour, layer jumps at the net's holes), over the whole grid.
pub fn label_cells(core: &Grid3<i16>, thru: &Grid3<i16>, nid: i16, cells: &[[i32; 3]], label: &mut [i32]) {
    let (ny, nx, nl) = (core.ny, core.nx, core.nl);
    let hw = ny * nx;
    let own = |l: usize, y: usize, x: usize| core.row(l, y)[x] == nid;
    let mut seen = std::collections::HashMap::<u32, i32, std::hash::BuildHasherDefault<Fx>>::default();
    let mut q: Vec<u32> = Vec::new();
    for (k, c) in cells.iter().enumerate() {
        let (l, y, x) = (c[0], c[1], c[2]);
        label[k] = -1;
        if l < 0 || y < 0 || x < 0 || l as usize >= nl || y as usize >= ny || x as usize >= nx {
            continue;
        }
        let (l, y, x) = (l as usize, y as usize, x as usize);
        if !own(l, y, x) {
            continue;
        }
        let s = (l * hw + y * nx + x) as u32;
        if let Some(&lab) = seen.get(&s) {
            label[k] = lab;
            continue;
        }
        let lab = k as i32;
        label[k] = lab;
        seen.insert(s, lab);
        q.clear();
        q.push(s);
        let mut head = 0;
        while head < q.len() {
            let s = q[head] as usize;
            head += 1;
            let l = s / hw;
            let c = s - l * hw;
            let (y, x) = (c / nx, c % nx);
            let mut nb = [usize::MAX; 4];
            if x > 0 { nb[0] = s - 1; }
            if x + 1 < nx { nb[1] = s + 1; }
            if y > 0 { nb[2] = s - nx; }
            if y + 1 < ny { nb[3] = s + nx; }
            for t in nb {
                if t != usize::MAX {
                    let (tl, tc) = (t / hw, t % hw);
                    if own(tl, tc / nx, tc % nx) && !seen.contains_key(&(t as u32)) {
                        seen.insert(t as u32, lab);
                        q.push(t as u32);
                    }
                }
            }
            if thru.row(0, y)[x] == nid {
                for l2 in 0..nl {
                    let t = l2 * hw + c;
                    if l2 != l && own(l2, y, x) && !seen.contains_key(&(t as u32)) {
                        seen.insert(t as u32, lab);
                        q.push(t as u32);
                    }
                }
            }
        }
    }
}
