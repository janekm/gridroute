//! Multi-layer grid A* with per-layer direction costs, a turn penalty and vias.
//!
//! State = (layer, row, col) flattened as `L*H*W + i*W + j`. Moves: 8 neighbours on a layer (cost per layer and
//! direction from `mcost[L*8+d]`); diagonal moves only when both orthogonal neighbours are free; a layer change
//! (via) at cells where `vok[i*W+j]` is set, to layers with `lay_ok` set. A turn penalty is added when the
//! direction changes. The heuristic is the octile distance to the target bounding box times `hmul`.
//!
//! This is a faithful port of the original C search (same heap, same tie-breaking, same fused multiply-adds as
//! clang emits for it), so it returns the same paths; the gain is in buffer reuse and the call overhead.

use std::cell::RefCell;
use std::sync::atomic::{AtomicBool, Ordering};

/// From this many states, `astar` runs the reachability flood fill (reach.rs) alongside the search.
pub const REACH_MIN_STATES: usize = 1 << 21;

#[derive(Clone, Copy)]
struct Node {
    f: f32,
    g: f32,
    s: i32,
}

fn hpush(h: &mut Vec<Node>, x: Node) {
    let mut i = h.len();
    h.push(x);
    while i > 0 {
        let p = (i - 1) / 2;
        if h[p].f <= x.f {
            break;
        }
        h[i] = h[p];
        i = p;
    }
    h[i] = x;
}

fn hpop(h: &mut Vec<Node>) -> Node {
    let top = h[0];
    let x = h.pop().unwrap();
    let n = h.len();
    if n == 0 {
        return top;
    }
    let mut i = 0;
    loop {
        let mut c = 2 * i + 1;
        if c >= n {
            break;
        }
        if c + 1 < n && h[c + 1].f < h[c].f {
            c += 1;
        }
        if h[c].f >= x.f {
            break;
        }
        h[i] = h[c];
        i = c;
    }
    h[i] = x;
    top
}


/// Binary heap with the same comparisons and moves as hpush / hpop (so the same pop order, ties included). Keys are
/// stored apart from the payload (a sift-down level touches one cache line), followed by two +inf sentinels so the
/// larger-child choice is branch-free (`f[n] = inf` never wins, exactly like the `c + 1 < n` test it replaces).
struct Heap {
    f: Vec<f32>,
    v: Vec<(f32, i32)>,
    n: usize,
}

impl Heap {
    const fn new() -> Self {
        Heap { f: Vec::new(), v: Vec::new(), n: 0 }
    }
    fn clear(&mut self) {
        self.n = 0;
        self.f.clear();
        self.f.extend([f32::INFINITY; 2]);
        self.v.clear();
    }
    fn is_empty(&self) -> bool {
        self.n == 0
    }
    #[inline(always)]
    fn push(&mut self, x: Node) {
        let mut i = self.n;
        self.n += 1;
        self.f.push(f32::INFINITY); // keep two sentinels past the end
        self.v.push((x.g, x.s));
        let (f, v) = (self.f.as_mut_ptr(), self.v.as_mut_ptr());
        unsafe {
            while i > 0 {
                let p = (i - 1) / 2;
                if *f.add(p) <= x.f {
                    break;
                }
                *f.add(i) = *f.add(p);
                *v.add(i) = *v.add(p);
                i = p;
            }
            *f.add(i) = x.f;
            *v.add(i) = (x.g, x.s);
        }
    }
    #[inline(always)]
    fn pop(&mut self) -> Node {
        let top = Node { f: self.f[0], g: self.v[0].0, s: self.v[0].1 };
        self.n -= 1;
        let n = self.n;
        let xf = self.f[n];
        let xv = self.v.pop().unwrap();
        self.f[n] = f32::INFINITY;
        self.f.pop();
        if n == 0 {
            return top;
        }
        let (f, v) = (self.f.as_mut_ptr(), self.v.as_mut_ptr());
        let mut i = 0;
        unsafe {
            loop {
                let mut c = 2 * i + 1;
                if c >= n {
                    break;
                }
                c += (*f.add(c + 1) < *f.add(c)) as usize;
                if *f.add(c) >= xf {
                    break;
                }
                *f.add(i) = *f.add(c);
                *v.add(i) = *v.add(c);
                i = c;
            }
            *f.add(i) = xf;
            *v.add(i) = xv;
        }
        top
    }
}

pub struct Problem<'a> {
    pub nl: usize,
    pub h: usize,
    pub w: usize,
    pub blk: &'a [u8],
    pub vok: &'a [u8],
    pub tgt: &'a [u8],
    pub src: &'a [i32],
    pub mcost: &'a [f32],
    pub lay_ok: &'a [u8],
    pub vcost: f32,
    pub turn: f32,
    pub hmul: f32,
    /// target bounding box (x0, y0, x1, y1) in cells
    pub tbox: (i32, i32, i32, i32),
    pub max_exp: i64,
    /// Optional read set: [ceil(h/32)][ceil(w/32)] bytes, set to 1 for every 32x32 tile holding a cell the search
    /// read (expanded cells and their neighbours, all layers). Null when not wanted.
    pub touched: *mut u8,
    /// Optional extra cost of entering each state ([nl*h*w], >= 0), e.g. negotiated congestion. `astar_ref`
    /// ignores it; the reachability check and cost-to-target fields stay admissible because costs are >= 0.
    pub cost: Option<&'a [f32]>,
}

pub const TOUCH_TILE: usize = 32;

// the read-set pointer is only written by the searching thread
unsafe impl Sync for Problem<'_> {}
unsafe impl Send for Problem<'_> {}

#[inline(always)]
fn touch(p: &Problem, i: i64, j: i64) {
    if p.touched.is_null() {
        return;
    }
    let tw = p.w.div_ceil(TOUCH_TILE) as i64;
    let (h, w) = (p.h as i64, p.w as i64);
    let (r0, r1) = (((i - 1).max(0)) >> 5, ((i + 1).min(h - 1)) >> 5);
    let (c0, c1) = (((j - 1).max(0)) >> 5, ((j + 1).min(w - 1)) >> 5);
    for r in r0..=r1 {
        for c in c0..=c1 {
            unsafe { *p.touched.add((r * tw + c) as usize) = 1 };
        }
    }
}

pub enum Outcome {
    Found(Vec<i32>),
    NotFound,
}

struct Scratch {
    dist: Vec<f32>,
    prev: Vec<i32>,
    dirn: Vec<u8>,
    heap: Vec<Node>,
}

thread_local! {
    static SCRATCH: RefCell<Scratch> = RefCell::new(Scratch { dist: Vec::new(), prev: Vec::new(), dirn: Vec::new(), heap: Vec::new() });
}

pub(crate) const DI: [i32; 8] = [0, 0, 1, -1, 1, 1, -1, -1];
pub(crate) const DJ: [i32; 8] = [1, -1, 0, 0, 1, -1, 1, -1];

/// Reference implementation (the original C search transcribed); `astar` must return the same results.
pub fn astar_ref(p: &Problem) -> (Outcome, i64) {
    SCRATCH.with(|sc| {
        let mut sc = sc.borrow_mut();
        let sc = &mut *sc;
        let (h, w) = (p.h as i64, p.w as i64);
        let hw = h * w;
        let n = p.nl * p.h * p.w;
        sc.dist.clear();
        sc.dist.resize(n, f32::INFINITY);
        sc.prev.clear();
        sc.prev.resize(n, -1);
        sc.dirn.clear();
        sc.dirn.resize(n, 8);
        let (dist, prev, dirn, heap) = (&mut sc.dist, &mut sc.prev, &mut sc.dirn, &mut sc.heap);
        heap.clear();
        let (blk, vok, tgt) = (p.blk, p.vok, p.tgt);
        let (tx0, ty0, tx1, ty1) = p.tbox;
        for &s in p.src {
            if s < 0 || s as usize >= n || blk[s as usize] != 0 {
                continue;
            }
            dist[s as usize] = 0.0;
            hpush(heap, Node { f: 0.0, g: 0.0, s });
        }
        let mut found: i32 = -1;
        let mut exp: i64 = 0;
        while !heap.is_empty() {
            let cur = hpop(heap);
            let s = cur.s as i64;
            let su = s as usize;
            if cur.g > dist[su] {
                continue;
            }
            if tgt[su] != 0 {
                found = cur.s;
                break;
            }
            exp += 1;
            if exp > p.max_exp {
                break;
            }
            let l = s / hw;
            let c = s - l * hw;
            let i = c / w;
            let j = c - i * w;
            let din = dirn[su];
            let lbase = l * hw;
            for d in 0..8usize {
                let i2 = i + DI[d] as i64;
                let j2 = j + DJ[d] as i64;
                if i2 < 0 || j2 < 0 || i2 >= h || j2 >= w {
                    continue;
                }
                let s2 = (lbase + i2 * w + j2) as usize;
                if blk[s2] != 0 {
                    continue;
                }
                if DI[d] != 0 && DJ[d] != 0 {
                    if blk[(lbase + i * w + j2) as usize] != 0 || blk[(lbase + i2 * w + j) as usize] != 0 {
                        continue;
                    }
                }
                let t = if din != 8 && din as usize != d { p.turn } else { 0.0 };
                let g2 = cur.g + p.mcost[l as usize * 8 + d] + t;
                if g2 < dist[s2] {
                    dist[s2] = g2;
                    prev[s2] = s as i32;
                    dirn[s2] = d as u8;
                    let (j2, i2) = (j2 as i32, i2 as i32);
                    let dx = if j2 < tx0 { tx0 - j2 } else if j2 > tx1 { j2 - tx1 } else { 0 };
                    let dy = if i2 < ty0 { ty0 - i2 } else if i2 > ty1 { i2 - ty1 } else { 0 };
                    let hh = (dx.min(dy) as f32).mul_add(-0.5858f32, (dx + dy) as f32);
                    hpush(heap, Node { f: hh.mul_add(p.hmul, g2), g: g2, s: s2 as i32 });
                }
            }
            if vok[c as usize] != 0 {
                for l2 in 0..p.nl as i64 {
                    if l2 == l || p.lay_ok[l2 as usize] == 0 {
                        continue;
                    }
                    let s2 = (l2 * hw + c) as usize;
                    if blk[s2] != 0 {
                        continue;
                    }
                    let g2 = cur.g + p.vcost;
                    if g2 < dist[s2] {
                        dist[s2] = g2;
                        prev[s2] = s as i32;
                        dirn[s2] = 8;
                        hpush(heap, Node { f: cur.f - cur.g + g2, g: g2, s: s2 as i32 });
                    }
                }
            }
        }
        if found < 0 {
            return (Outcome::NotFound, exp);
        }
        let mut path = Vec::new();
        let mut s = found;
        while s >= 0 {
            path.push(s);
            s = prev[s as usize];
        }
        path.reverse();
        (Outcome::Found(path), exp)
    })
}

/// Search state: g (−∞ marks a blocked state, so the relaxation test `g2 < d` rejects it without a separate load)
/// and the predecessor (low 28 bits, +1) with the arrival direction (high 4 bits, 8 = none / via).
#[derive(Clone, Copy)]
#[repr(C)]
struct St {
    d: f32,
    pd: u32,
}

const PMASK: u32 = (1 << 28) - 1;

thread_local! {
    static STATES: RefCell<(Vec<St>, Heap)> = RefCell::new((Vec::new(), Heap::new()));
}

/// Returns the path (source first) and the number of expansions. Same results as `astar_ref`. Large searches
/// run a reachability check on a second thread and stop early when it proves the target unreachable.
pub fn astar(p: &Problem) -> (Outcome, i64) {
    let n = p.nl * p.h * p.w;
    if n >= PMASK as usize {
        return astar_ref(p);
    }
    if n < REACH_MIN_STATES {
        return astar_fast(p, None);
    }
    let cancel = AtomicBool::new(false);
    let unreachable = AtomicBool::new(false);
    std::thread::scope(|sc| {
        sc.spawn(|| {
            if crate::reach::reachable(p, &cancel) == Some(false) {
                unreachable.store(true, Ordering::Relaxed);
            }
        });
        let r = astar_fast(p, Some(&unreachable));
        cancel.store(true, Ordering::Relaxed);
        r
    })
}

/// `astar` without the concurrent reachability check.
pub fn astar_plain(p: &Problem) -> (Outcome, i64) {
    astar_fast(p, None)
}

fn astar_fast(p: &Problem, stop: Option<&AtomicBool>) -> (Outcome, i64) {
    let (tx0, ty0, tx1, ty1) = p.tbox;
    astar_core(p, stop, p.hmul, false, |_s2, i2, j2| {
        let dx = if j2 < tx0 { tx0 - j2 } else if j2 > tx1 { j2 - tx1 } else { 0 };
        let dy = if i2 < ty0 { ty0 - i2 } else if i2 > ty1 { i2 - ty1 } else { 0 };
        (dx.min(dy) as f32).mul_add(-0.5858f32, (dx + dy) as f32)
    })
}

/// A* guided by a cost-to-target field (field.rs) instead of the octile distance: f = g + weight * field.
/// States with an infinite field value cannot reach a target and are never queued. Not result-identical to
/// `astar` (different heuristic), but with weight 1 it returns a cheapest path and expands far fewer states.
pub fn astar_guided(p: &Problem, field: &[f32], weight: f32) -> (Outcome, i64) {
    astar_core(p, None, weight, true, |s2, _i2, _j2| field[s2])
}

/// `astar_guided` with the field given as a lookup (state -> cost to target), e.g. a coarse field.
pub fn astar_guided_by<F: Fn(usize) -> f32>(p: &Problem, field: F, weight: f32) -> (Outcome, i64) {
    astar_core(p, None, weight, true, |s2, _i2, _j2| field(s2))
}

#[inline(always)]
fn astar_core<HF: Fn(usize, i32, i32) -> f32>(
    p: &Problem,
    stop: Option<&AtomicBool>,
    hmul: f32,
    via_h: bool,
    heur: HF,
) -> (Outcome, i64) {
    let n = p.nl * p.h * p.w;
    STATES.with(|cell| {
        let mut guard = cell.borrow_mut();
        let (st, heap) = &mut *guard;
        // (re)initialise in parallel: blocked states get −∞
        if st.len() < n {
            st.resize(n, St { d: 0.0, pd: 0 });
        }
        {
            use rayon::prelude::*;
            const CH: usize = 1 << 16;
            st[..n].par_chunks_mut(CH).zip(p.blk[..n].par_chunks(CH)).for_each(|(a, b)| {
                for (x, &k) in a.iter_mut().zip(b) {
                    *x = St { d: if k != 0 { f32::NEG_INFINITY } else { f32::INFINITY }, pd: 8 << 28 };
                }
            });
        }
        let st = &mut st[..n];
        heap.clear();
        let (h, w) = (p.h as i64, p.w as i64);
        let hw = h * w;
        let (vok, tgt) = (p.vok, p.tgt);
        let ecost = |s2: usize| p.cost.map_or(0.0, |c| c[s2]);
        let mut off = [0i64; 8];
        for d in 0..8 {
            off[d] = DI[d] as i64 * w + DJ[d] as i64;
        }
        for &s in p.src {
            if s < 0 || s as usize >= n || p.blk[s as usize] != 0 {
                continue;
            }
            st[s as usize].d = 0.0;
            heap.push(Node { f: 0.0, g: 0.0, s });
        }
        let mut found: i32 = -1;
        let mut exp: i64 = 0;
        while !heap.is_empty() {
            let cur = heap.pop();
            let s = cur.s as i64;
            let su = s as usize;
            if cur.g > st[su].d {
                continue;
            }
            if tgt[su] != 0 {
                found = cur.s;
                break;
            }
            exp += 1;
            if exp > p.max_exp {
                break;
            }
            if exp & 0xfff == 0 && stop.is_some_and(|f| f.load(Ordering::Relaxed)) {
                break;
            }
            let l = s / hw;
            let c = s - l * hw;
            let i = c / w;
            let j = c - i * w;
            touch(p, i, j);
            let din = st[su].pd >> 28;
            let mrow = &p.mcost[l as usize * 8..l as usize * 8 + 8];
            for d in 0..8usize {
                let i2 = i + DI[d] as i64;
                let j2 = j + DJ[d] as i64;
                if i2 < 0 || j2 < 0 || i2 >= h || j2 >= w {
                    continue;
                }
                let s2 = (s + off[d]) as usize;
                let d2 = st[s2].d;
                if d2 == f32::NEG_INFINITY {
                    continue;
                }
                if d >= 4 && (st[(s + DJ[d] as i64) as usize].d == f32::NEG_INFINITY
                    || st[(s + DI[d] as i64 * w) as usize].d == f32::NEG_INFINITY)
                {
                    continue;
                }
                let t = if din != 8 && din as usize != d { p.turn } else { 0.0 };
                let g2 = cur.g + mrow[d] + t + ecost(s2);
                if g2 < d2 {
                    st[s2] = St { d: g2, pd: ((d as u32) << 28) | (s as u32 + 1) };
                    let hh = heur(s2, i2 as i32, j2 as i32);
                    if hh == f32::INFINITY {
                        continue;
                    }
                    heap.push(Node { f: hh.mul_add(hmul, g2), g: g2, s: s2 as i32 });
                }
            }
            if vok[c as usize] != 0 {
                for l2 in 0..p.nl as i64 {
                    if l2 == l || p.lay_ok[l2 as usize] == 0 {
                        continue;
                    }
                    let s2 = (l2 * hw + c) as usize;
                    let g2 = cur.g + p.vcost + ecost(s2);
                    if g2 < st[s2].d {
                        st[s2] = St { d: g2, pd: (8 << 28) | (s as u32 + 1) };
                        let f = if via_h {
                            let hh = heur(s2, i as i32, j as i32);
                            if hh == f32::INFINITY {
                                continue;
                            }
                            hh.mul_add(hmul, g2)
                        } else {
                            cur.f - cur.g + g2
                        };
                        heap.push(Node { f, g: g2, s: s2 as i32 });
                    }
                }
            }
        }
        if found < 0 {
            return (Outcome::NotFound, exp);
        }
        let mut path = Vec::new();
        let mut s = found as u32;
        loop {
            path.push(s as i32);
            let pr = st[s as usize].pd & PMASK;
            if pr == 0 {
                break;
            }
            s = pr - 1;
        }
        path.reverse();
        (Outcome::Found(path), exp)
    })
}
