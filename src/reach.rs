//! Reachability with exactly the A* move set: is any target state connected to a source state?
//!
//! A search that cannot reach its target fails only after exhausting its open set or its expansion budget, which
//! for a board-sized window can take seconds. This flood fill answers the same question in one linear pass, so
//! `astar` runs it alongside the search and stops as soon as it proves the target unreachable (the result, "no
//! path", is the one the search would have returned).

use crate::astar::{Problem, DI, DJ};
use std::sync::atomic::{AtomicBool, Ordering};

/// Some(true) if a target is reachable, Some(false) if not, None if `cancel` was raised first.
pub fn reachable(p: &Problem, cancel: &AtomicBool) -> Option<bool> {
    let (h, w) = (p.h as i64, p.w as i64);
    let hw = h * w;
    let n = p.nl * p.h * p.w;
    let mut seen = vec![0u64; n.div_ceil(64)];
    let test = |seen: &Vec<u64>, k: usize| seen[k >> 6] >> (k & 63) & 1 != 0;
    let mut q: Vec<u32> = Vec::new();
    for &s in p.src {
        if s < 0 || s as usize >= n || p.blk[s as usize] != 0 {
            continue;
        }
        let k = s as usize;
        if !test(&seen, k) {
            seen[k >> 6] |= 1 << (k & 63);
            q.push(s as u32);
        }
    }
    let mut head = 0;
    while head < q.len() {
        if head & 0xffff == 0 && cancel.load(Ordering::Relaxed) {
            return None;
        }
        let s = q[head] as i64;
        head += 1;
        if p.tgt[s as usize] != 0 {
            return Some(true);
        }
        let l = s / hw;
        let c = s - l * hw;
        let i = c / w;
        let j = c - i * w;
        let lb = l * hw;
        for d in 0..8 {
            let (i2, j2) = (i + DI[d] as i64, j + DJ[d] as i64);
            if i2 < 0 || j2 < 0 || i2 >= h || j2 >= w {
                continue;
            }
            let s2 = (lb + i2 * w + j2) as usize;
            if p.blk[s2] != 0 || test(&seen, s2) {
                continue;
            }
            if d >= 4 && (p.blk[(lb + i * w + j2) as usize] != 0 || p.blk[(lb + i2 * w + j) as usize] != 0) {
                continue;
            }
            seen[s2 >> 6] |= 1 << (s2 & 63);
            q.push(s2 as u32);
        }
        if p.vok[c as usize] != 0 {
            for l2 in 0..p.nl as i64 {
                if l2 == l || p.lay_ok[l2 as usize] == 0 {
                    continue;
                }
                let s2 = (l2 * hw + c) as usize;
                if p.blk[s2] != 0 || test(&seen, s2) {
                    continue;
                }
                seen[s2 >> 6] |= 1 << (s2 & 63);
                q.push(s2 as u32);
            }
        }
    }
    Some(false)
}
