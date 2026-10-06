//! Cost-to-target field: for every state, the cheapest cost to any target state over the search graph without
//! turn penalties (a relaxation of the A* graph, so the field is an admissible and consistent heuristic).
//!
//! `field_cpu` is a backward Dijkstra (exact). The GPU version (gpu_field.rs) computes the same field by parallel
//! relaxation sweeps.

use crate::astar::{Problem, DI, DJ};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[derive(PartialEq)]
struct E(f32, u32);
impl Eq for E {}
impl PartialOrd for E {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for E {
    fn cmp(&self, o: &Self) -> Ordering {
        o.0.partial_cmp(&self.0).unwrap_or(Ordering::Equal).then(o.1.cmp(&self.1))
    }
}

/// Exact field by backward Dijkstra; unreachable states get +inf. `out`: nl * h * w.
pub fn field_cpu(p: &Problem, out: &mut [f32]) {
    let (h, w) = (p.h as i64, p.w as i64);
    let hw = h * w;
    let n = p.nl * p.h * p.w;
    out[..n].fill(f32::INFINITY);
    let mut q = BinaryHeap::new();
    for s in 0..n {
        if p.tgt[s] != 0 && p.blk[s] == 0 {
            out[s] = 0.0;
            q.push(E(0.0, s as u32));
        }
    }
    while let Some(E(g, s2)) = q.pop() {
        let s2 = s2 as i64;
        if g > out[s2 as usize] {
            continue;
        }
        let l = s2 / hw;
        let c = s2 - l * hw;
        let i2 = c / w;
        let j2 = c - i2 * w;
        let lb = l * hw;
        // predecessors s = s2 - (DI, DJ): the forward move s -> s2 in direction d
        for d in 0..8 {
            let (i, j) = (i2 - DI[d] as i64, j2 - DJ[d] as i64);
            if i < 0 || j < 0 || i >= h || j >= w {
                continue;
            }
            let s = (lb + i * w + j) as usize;
            if p.blk[s] != 0 {
                continue;
            }
            if d >= 4 && (p.blk[(lb + i * w + j2) as usize] != 0 || p.blk[(lb + i2 * w + j) as usize] != 0) {
                continue;
            }
            let g1 = g + p.mcost[l as usize * 8 + d];
            if g1 < out[s] {
                out[s] = g1;
                q.push(E(g1, s as u32));
            }
        }
        // a via from (l1, c) to (l, c) needs vok[c] and lay_ok[l]
        if p.vok[c as usize] != 0 && p.lay_ok[l as usize] != 0 {
            for l1 in 0..p.nl as i64 {
                if l1 == l {
                    continue;
                }
                let s = (l1 * hw + c) as usize;
                if p.blk[s] != 0 {
                    continue;
                }
                let g1 = g + p.vcost;
                if g1 < out[s] {
                    out[s] = g1;
                    q.push(E(g1, s as u32));
                }
            }
        }
    }
}

/// A coarse copy of the search graph, `f` x `f` cells per coarse cell, for an approximate field: a coarse cell is
/// blocked only if all its cells are (optimistic), allows a via if any cell does, is a target if any free cell is,
/// and a coarse move costs `f` times the fine move. The field on it, looked up per fine cell, is a cheap
/// heuristic that still knows the walls; it can underestimate through gaps narrower than a coarse cell.
pub struct Coarse {
    pub f: usize,
    pub h: usize,
    pub w: usize,
    pub blk: Vec<u8>,
    pub vok: Vec<u8>,
    pub tgt: Vec<u8>,
    pub mcost: Vec<f32>,
}

pub fn coarsen(p: &Problem, f: usize) -> Coarse {
    use rayon::prelude::*;
    let (h, w) = (p.h.div_ceil(f), p.w.div_ceil(f));
    let (nl, hf, wf) = (p.nl, p.h, p.w);
    let mut blk = vec![1u8; nl * h * w];
    let mut tgt = vec![0u8; nl * h * w];
    let mut vok = vec![0u8; h * w];
    // one coarse row (of one layer) per task
    blk.par_chunks_mut(w).zip(tgt.par_chunks_mut(w)).enumerate().for_each(|(r, (brow, trow))| {
        let (l, ci) = (r / h, r % h);
        for i in ci * f..((ci + 1) * f).min(hf) {
            let row = (l * hf + i) * wf;
            let (bf, tf) = (&p.blk[row..row + wf], &p.tgt[row..row + wf]);
            for j in 0..wf {
                if bf[j] == 0 {
                    brow[j / f] = 0;
                    if tf[j] != 0 {
                        trow[j / f] = 1;
                    }
                }
            }
        }
    });
    vok.par_chunks_mut(w).enumerate().for_each(|(ci, vrow)| {
        for i in ci * f..((ci + 1) * f).min(hf) {
            let vf = &p.vok[i * wf..(i + 1) * wf];
            for j in 0..wf {
                if vf[j] != 0 {
                    vrow[j / f] = 1;
                }
            }
        }
    });
    let mcost = p.mcost.iter().map(|&c| c * f as f32).collect();
    Coarse { f, h, w, blk, vok, tgt, mcost }
}

impl Coarse {
    pub fn problem<'a>(&'a self, p: &Problem<'a>) -> Problem<'a> {
        Problem {
            nl: p.nl,
            h: self.h,
            w: self.w,
            blk: &self.blk,
            vok: &self.vok,
            tgt: &self.tgt,
            src: &[],
            mcost: &self.mcost,
            lay_ok: p.lay_ok,
            vcost: p.vcost,
            turn: 0.0,
            hmul: 1.0,
            tbox: (0, 0, 0, 0),
            max_exp: 0,
            touched: std::ptr::null_mut(),
            cost: None,
        }
    }
}
