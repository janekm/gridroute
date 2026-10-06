//! gridroute: kernels for grid-based PCB routers, with a C ABI (see include/gridroute.h).
//!
//! - `gr_dilate`: obstacle dilation with net exclusions (clearance masks), multi-threaded CPU or Metal GPU.
//! - `gr_island`: copper island flood fill through vias / plated holes.
//! - `gr_astar`: multi-layer A* with direction costs, turn penalty and vias (same paths as the original C search).
//! - `gr_field`: cost-to-target field (Metal tiled relaxation or CPU Dijkstra); `gr_astar_hybrid`: A* that switches
//!   to field guidance when a search gets expensive.
//! - `gr_hash_*`: content hashes for memoising searches across runs.
//!
//! All grids are dense row-major arrays owned by the caller; nothing is retained between calls.

pub mod astar;
pub mod dilate;
pub mod field;
pub mod flood;
pub mod hash;
pub mod reach;
#[cfg(all(target_os = "macos", feature = "metal"))]
pub mod gpu;

use dilate::{Element, Grid3, Window};
use std::sync::atomic::{AtomicI64, Ordering};

pub const BACKEND_AUTO: i32 = 0;
pub const BACKEND_CPU: i32 = 1;
pub const BACKEND_METAL: i32 = 2;
pub const BACKEND_NAIVE: i32 = 3;

/// Work size (output cells x element rows x layers) from which `BACKEND_AUTO` uses the GPU.
static METAL_THRESHOLD: AtomicI64 = AtomicI64::new(i64::MAX);

/// Cells per side of a field cell in gr_astar_hybrid (1 = exact field at search resolution).
static FIELD_FACTOR: AtomicI64 = AtomicI64::new(1);

/// Resolution of the hybrid search's field: `f` x `f` search cells per field cell (1 = exact; 2 or 3 trade a
/// slightly weaker heuristic for f^2 less field work).
#[no_mangle]
pub extern "C" fn gr_set_field_factor(f: i32) {
    FIELD_FACTOR.store(f.max(1) as i64, Ordering::Relaxed);
}

/// Field for `p` at the configured factor; returns (values, coarse w, coarse h*w, factor).
fn compute_field(p: &astar::Problem) -> (Vec<f32>, usize, usize, usize) {
    let f = FIELD_FACTOR.load(Ordering::Relaxed) as usize;
    let co;
    let q;
    let tc = std::time::Instant::now();
    let (pp, w, hw) = if f > 1 {
        co = field::coarsen(p, f);
        if std::env::var_os("GR_TIME").is_some() {
            eprintln!("coarsen {:.1} ms ({} -> {} cells)", tc.elapsed().as_secs_f64() * 1e3, p.nl * p.h * p.w, p.nl * co.h * co.w);
        }
        q = co.problem(p);
        (&q, co.w, co.h * co.w)
    } else {
        (p, p.w, p.h * p.w)
    };
    let mut fv = vec![0f32; pp.nl * hw];
    #[cfg(all(target_os = "macos", feature = "metal"))]
    let gpu_ok = std::env::var_os("GR_FIELD_CPU").is_none() && gpu::field_metal(pp, &mut fv).is_some();
    #[cfg(not(all(target_os = "macos", feature = "metal")))]
    let gpu_ok = false;
    if !gpu_ok {
        field::field_cpu(pp, &mut fv);
    }
    (fv, w, hw, f)
}

#[no_mangle]
pub extern "C" fn gr_version() -> u32 {
    1
}

#[no_mangle]
pub extern "C" fn gr_metal_available() -> i32 {
    #[cfg(all(target_os = "macos", feature = "metal"))]
    {
        gpu::available() as i32
    }
    #[cfg(not(all(target_os = "macos", feature = "metal")))]
    {
        0
    }
}

/// Copies the GPU name (NUL-terminated) into buf; returns its length (0 without a GPU).
#[no_mangle]
pub unsafe extern "C" fn gr_device_name(buf: *mut u8, cap: usize) -> usize {
    #[cfg(all(target_os = "macos", feature = "metal"))]
    let name = gpu::device_name();
    #[cfg(not(all(target_os = "macos", feature = "metal")))]
    let name = String::new();
    let b = name.as_bytes();
    let n = b.len().min(cap.saturating_sub(1));
    std::ptr::copy_nonoverlapping(b.as_ptr(), buf, n);
    if cap > 0 {
        *buf.add(n) = 0;
    }
    n
}

/// Work size from which BACKEND_AUTO dilations run on the GPU (negative: never).
#[no_mangle]
pub extern "C" fn gr_set_metal_threshold(work: i64) {
    METAL_THRESHOLD.store(if work < 0 { i64::MAX } else { work }, Ordering::Relaxed);
}

/// Number of CPU worker threads (call before the first kernel; 0 = one per core).
#[no_mangle]
pub extern "C" fn gr_set_threads(n: i32) -> i32 {
    rayon::ThreadPoolBuilder::new().num_threads(n.max(0) as usize).build_global().is_ok() as i32
}

/// Obstacle dilation (see dilate.rs).
///
/// src: [nl][ny][nx] int16 net ids (0 = free). layers[nlayers]: source layers. excl[nexcl]: values that are not
/// obstacles. extra: optional [ny][nx] uint8 mask ORed into the obstacles of the layers whose extra_on[l] is set
/// (extra_on has nl entries; may be NULL when extra is NULL). Window (i0, j0, h, w) must lie inside the grid.
/// span[nspan] (odd length 2 rc + 1): half-width per element row, -1 = row unused. reduce_or: 1 = a single output
/// plane (OR over the layers), else one plane per layer. out: planes * h * w bytes.
/// Returns the backend used (1 CPU, 2 Metal, 3 naive) or a negative error.
#[no_mangle]
pub unsafe extern "C" fn gr_dilate(
    src: *const i16,
    nl: i32,
    ny: i32,
    nx: i32,
    layers: *const i32,
    nlayers: i32,
    excl: *const i16,
    nexcl: i32,
    extra: *const u8,
    extra_on: *const u8,
    i0: i32,
    j0: i32,
    h: i32,
    w: i32,
    span: *const i32,
    nspan: i32,
    reduce_or: i32,
    backend: i32,
    out: *mut u8,
) -> i32 {
    if src.is_null() || out.is_null() || nlayers <= 0 || nspan <= 0 || nspan % 2 == 0 || h <= 0 || w <= 0 {
        return -1;
    }
    if i0 < 0 || j0 < 0 || i0 + h > ny || j0 + w > nx {
        return -2;
    }
    let g = Grid3 { ptr: src, nl: nl as usize, ny: ny as usize, nx: nx as usize };
    let ls: Vec<usize> = std::slice::from_raw_parts(layers, nlayers as usize).iter().map(|&l| l as usize).collect();
    if ls.iter().any(|&l| l >= g.nl) {
        return -3;
    }
    let ex: &[i16] = if nexcl > 0 { std::slice::from_raw_parts(excl, nexcl as usize) } else { &[] };
    let eg = Grid3 { ptr: extra, nl: 1, ny: g.ny, nx: g.nx };
    let extra_g = if extra.is_null() { None } else { Some(&eg) };
    let on: Vec<bool> = if extra_on.is_null() {
        vec![!extra.is_null(); g.nl]
    } else {
        std::slice::from_raw_parts(extra_on, g.nl).iter().map(|&v| v != 0).collect()
    };
    let win = Window { i0: i0 as usize, j0: j0 as usize, h: h as usize, w: w as usize };
    let el = Element { span: std::slice::from_raw_parts(span, nspan as usize) };
    let planes = if reduce_or != 0 { 1 } else { ls.len() };
    let o = std::slice::from_raw_parts_mut(out, planes * win.h * win.w);
    let work = (win.h * win.w * el.span.len() * ls.len()) as i64;
    let want_gpu = backend == BACKEND_METAL || (backend == BACKEND_AUTO && work >= METAL_THRESHOLD.load(Ordering::Relaxed));
    if want_gpu {
        #[cfg(all(target_os = "macos", feature = "metal"))]
        if gpu::dilate_metal(&g, &ls, ex, extra_g, &on, &win, &el, reduce_or != 0, o) {
            return BACKEND_METAL;
        }
        if backend == BACKEND_METAL {
            return -4;
        }
    }
    if backend == BACKEND_NAIVE {
        dilate::dilate_naive(&g, &ls, ex, extra_g, &on, &win, &el, reduce_or != 0, o);
        return BACKEND_NAIVE;
    }
    dilate::dilate_cpu(&g, &ls, ex, extra_g, &on, &win, &el, reduce_or != 0, o);
    BACKEND_CPU
}

/// Copper island of net `nid` connected to the seeds (see flood.rs). core: [nl][ny][nx], thru: [ny][nx].
/// seeds: nseeds x (layer, y, x) relative to the window. out: nl * h * w bytes. Returns the cell count.
#[no_mangle]
pub unsafe extern "C" fn gr_island(
    core: *const i16,
    nl: i32,
    ny: i32,
    nx: i32,
    thru: *const i16,
    nid: i16,
    i0: i32,
    j0: i32,
    h: i32,
    w: i32,
    seeds: *const i32,
    nseeds: i32,
    out: *mut u8,
) -> i64 {
    if i0 < 0 || j0 < 0 || h <= 0 || w <= 0 || i0 + h > ny || j0 + w > nx {
        return -2;
    }
    let c = Grid3 { ptr: core, nl: nl as usize, ny: ny as usize, nx: nx as usize };
    let t = Grid3 { ptr: thru, nl: 1, ny: ny as usize, nx: nx as usize };
    let win = Window { i0: i0 as usize, j0: j0 as usize, h: h as usize, w: w as usize };
    let sd: &[[i32; 3]] =
        if nseeds > 0 { std::slice::from_raw_parts(seeds as *const [i32; 3], nseeds as usize) } else { &[] };
    let o = std::slice::from_raw_parts_mut(out, nl as usize * win.h * win.w);
    flood::island(&c, &t, nid, &win, sd, o);
    o.iter().map(|&v| v as i64).sum()
}

/// A* search (see astar.rs); the same arguments as the original astar.c plus `nexp` (expansions, may be NULL).
/// Returns the path length written to out (source first), -1 if no path was found, -2 if out is too small.
#[no_mangle]
pub unsafe extern "C" fn gr_astar(
    nl: i32,
    h: i32,
    w: i32,
    blk: *const u8,
    vok: *const u8,
    tgt: *const u8,
    src: *const i32,
    nsrc: i32,
    mcost: *const f32,
    lay_ok: *const u8,
    vcost: f32,
    turn: f32,
    hmul: f32,
    tx0: i32,
    ty0: i32,
    tx1: i32,
    ty1: i32,
    max_exp: i64,
    out: *mut i32,
    out_cap: i32,
    nexp: *mut i64,
) -> i32 {
    gr_astar_impl(nl, h, w, blk, vok, tgt, src, nsrc, mcost, lay_ok, vcost, turn, hmul, tx0, ty0, tx1, ty1, max_exp, out, out_cap, nexp, false)
}

/// gr_astar with the reference (unoptimised) search; for equivalence tests and benchmarks.
#[no_mangle]
pub unsafe extern "C" fn gr_astar_ref(
    nl: i32, h: i32, w: i32, blk: *const u8, vok: *const u8, tgt: *const u8, src: *const i32, nsrc: i32,
    mcost: *const f32, lay_ok: *const u8, vcost: f32, turn: f32, hmul: f32, tx0: i32, ty0: i32, tx1: i32, ty1: i32,
    max_exp: i64, out: *mut i32, out_cap: i32, nexp: *mut i64,
) -> i32 {
    gr_astar_impl(nl, h, w, blk, vok, tgt, src, nsrc, mcost, lay_ok, vcost, turn, hmul, tx0, ty0, tx1, ty1, max_exp, out, out_cap, nexp, true)
}

#[allow(clippy::too_many_arguments)]
unsafe fn gr_astar_impl(
    nl: i32, h: i32, w: i32, blk: *const u8, vok: *const u8, tgt: *const u8, src: *const i32, nsrc: i32,
    mcost: *const f32, lay_ok: *const u8, vcost: f32, turn: f32, hmul: f32, tx0: i32, ty0: i32, tx1: i32, ty1: i32,
    max_exp: i64, out: *mut i32, out_cap: i32, nexp: *mut i64, reference: bool,
) -> i32 {
    let n = nl as usize * h as usize * w as usize;
    let p = astar::Problem {
        nl: nl as usize,
        h: h as usize,
        w: w as usize,
        blk: std::slice::from_raw_parts(blk, n),
        vok: std::slice::from_raw_parts(vok, h as usize * w as usize),
        tgt: std::slice::from_raw_parts(tgt, n),
        src: std::slice::from_raw_parts(src, nsrc.max(0) as usize),
        mcost: std::slice::from_raw_parts(mcost, nl as usize * 8),
        lay_ok: std::slice::from_raw_parts(lay_ok, nl as usize),
        vcost,
        turn,
        hmul,
        tbox: (tx0, ty0, tx1, ty1),
        max_exp,
        touched: std::ptr::null_mut(),
        cost: None,
    };
    let (res, e) = if reference { astar::astar_ref(&p) } else { astar::astar(&p) };
    if !nexp.is_null() {
        *nexp = e;
    }
    match res {
        astar::Outcome::NotFound => -1,
        astar::Outcome::Found(path) => {
            if path.len() > out_cap as usize {
                return -2;
            }
            std::ptr::copy_nonoverlapping(path.as_ptr(), out, path.len());
            path.len() as i32
        }
    }
}

/// Field-guided A* (experimental, not result-identical to gr_astar): computes the cost-to-target field (backend
/// GR_CPU: exact backward Dijkstra) and searches with f = g + weight * field. Same arguments as gr_astar, plus
/// `weight` and `field_secs` (time spent on the field, may be NULL).
#[no_mangle]
pub unsafe extern "C" fn gr_astar_guided(
    nl: i32, h: i32, w: i32, blk: *const u8, vok: *const u8, tgt: *const u8, src: *const i32, nsrc: i32,
    mcost: *const f32, lay_ok: *const u8, vcost: f32, turn: f32, weight: f32, max_exp: i64, out: *mut i32,
    out_cap: i32, nexp: *mut i64, field_secs: *mut f64,
) -> i32 {
    let n = nl as usize * h as usize * w as usize;
    let p = astar::Problem {
        nl: nl as usize,
        h: h as usize,
        w: w as usize,
        blk: std::slice::from_raw_parts(blk, n),
        vok: std::slice::from_raw_parts(vok, h as usize * w as usize),
        tgt: std::slice::from_raw_parts(tgt, n),
        src: std::slice::from_raw_parts(src, nsrc.max(0) as usize),
        mcost: std::slice::from_raw_parts(mcost, nl as usize * 8),
        lay_ok: std::slice::from_raw_parts(lay_ok, nl as usize),
        vcost,
        turn,
        hmul: weight,
        tbox: (0, 0, 0, 0),
        max_exp,
        touched: std::ptr::null_mut(),
        cost: None,
    };
    let t = std::time::Instant::now();
    let mut field = vec![0f32; n];
    #[cfg(all(target_os = "macos", feature = "metal"))]
    let gpu_ok = std::env::var("GR_FIELD_CPU").is_err() && gpu::field_metal(&p, &mut field).is_some();
    #[cfg(not(all(target_os = "macos", feature = "metal")))]
    let gpu_ok = false;
    if !gpu_ok {
        field::field_cpu(&p, &mut field);
    }
    if !field_secs.is_null() {
        *field_secs = t.elapsed().as_secs_f64();
    }
    let (res, e) = astar::astar_guided(&p, &field, weight);
    if !nexp.is_null() {
        *nexp = e;
    }
    match res {
        astar::Outcome::NotFound => -1,
        astar::Outcome::Found(path) => {
            if path.len() > out_cap as usize {
                return -2;
            }
            std::ptr::copy_nonoverlapping(path.as_ptr(), out, path.len());
            path.len() as i32
        }
    }
}

/// Cost-to-target field (field.rs): out[nl*h*w] = cheapest cost from each state to a target without turn
/// penalties (+inf if none). backend GR_CPU: exact Dijkstra; GR_METAL / GR_AUTO: tiled GPU relaxation (CPU if no
/// GPU). Returns the backend used (or < 0); `rounds` (may be NULL) gets the GPU round count.
#[no_mangle]
pub unsafe extern "C" fn gr_field(
    nl: i32, h: i32, w: i32, blk: *const u8, vok: *const u8, tgt: *const u8, mcost: *const f32, lay_ok: *const u8,
    vcost: f32, backend: i32, out: *mut f32, rounds: *mut i64,
) -> i32 {
    let n = nl as usize * h as usize * w as usize;
    let p = astar::Problem {
        nl: nl as usize,
        h: h as usize,
        w: w as usize,
        blk: std::slice::from_raw_parts(blk, n),
        vok: std::slice::from_raw_parts(vok, h as usize * w as usize),
        tgt: std::slice::from_raw_parts(tgt, n),
        src: &[],
        mcost: std::slice::from_raw_parts(mcost, nl as usize * 8),
        lay_ok: std::slice::from_raw_parts(lay_ok, nl as usize),
        vcost,
        turn: 0.0,
        hmul: 1.0,
        tbox: (0, 0, 0, 0),
        max_exp: 0,
        touched: std::ptr::null_mut(),
        cost: None,
    };
    let o = std::slice::from_raw_parts_mut(out, n);
    if backend != BACKEND_CPU {
        #[cfg(all(target_os = "macos", feature = "metal"))]
        if let Some(r) = gpu::field_metal(&p, o) {
            if !rounds.is_null() {
                *rounds = r as i64;
            }
            return BACKEND_METAL;
        }
        if backend == BACKEND_METAL {
            return -4;
        }
    }
    field::field_cpu(&p, o);
    BACKEND_CPU
}

/// Hybrid search: gr_astar with an expansion `budget`; a search that is still open then switches to the
/// cost-to-target field (GPU when available) and continues as field-guided A* (f = g + weight * field) within
/// max_exp. Easy searches return exactly what gr_astar returns; hard ones return a near-cheapest path much faster
/// (and may succeed where gr_astar's budget runs out). info (may be NULL): [phase (1 = octile, 2 = guided,
/// 3 = no path by field), field seconds, total expansions].
#[no_mangle]
pub unsafe extern "C" fn gr_astar_hybrid(
    nl: i32, h: i32, w: i32, blk: *const u8, vok: *const u8, tgt: *const u8, src: *const i32, nsrc: i32,
    mcost: *const f32, lay_ok: *const u8, vcost: f32, turn: f32, hmul: f32, tx0: i32, ty0: i32, tx1: i32, ty1: i32,
    max_exp: i64, out: *mut i32, out_cap: i32, budget: i64, weight: f32, info: *mut f64, touched: *mut u8,
) -> i32 {
    gr_astar_hybrid_cost(nl, h, w, blk, vok, tgt, src, nsrc, mcost, lay_ok, vcost, turn, hmul, tx0, ty0, tx1, ty1,
                         max_exp, out, out_cap, budget, weight, info, touched, std::ptr::null())
}

/// gr_astar_hybrid with `cost` (may be NULL): float [nl*h*w], an extra cost (>= 0) for entering each state.
#[no_mangle]
pub unsafe extern "C" fn gr_astar_hybrid_cost(
    nl: i32, h: i32, w: i32, blk: *const u8, vok: *const u8, tgt: *const u8, src: *const i32, nsrc: i32,
    mcost: *const f32, lay_ok: *const u8, vcost: f32, turn: f32, hmul: f32, tx0: i32, ty0: i32, tx1: i32, ty1: i32,
    max_exp: i64, out: *mut i32, out_cap: i32, budget: i64, weight: f32, info: *mut f64, touched: *mut u8,
    cost: *const f32,
) -> i32 {
    let n = nl as usize * h as usize * w as usize;
    let mut p = astar::Problem {
        nl: nl as usize,
        h: h as usize,
        w: w as usize,
        blk: std::slice::from_raw_parts(blk, n),
        vok: std::slice::from_raw_parts(vok, h as usize * w as usize),
        tgt: std::slice::from_raw_parts(tgt, n),
        src: std::slice::from_raw_parts(src, nsrc.max(0) as usize),
        mcost: std::slice::from_raw_parts(mcost, nl as usize * 8),
        lay_ok: std::slice::from_raw_parts(lay_ok, nl as usize),
        vcost,
        turn,
        hmul,
        tbox: (tx0, ty0, tx1, ty1),
        max_exp: max_exp.min(budget),
        touched,
        cost: if cost.is_null() { None } else { Some(std::slice::from_raw_parts(cost, n)) },
    };
    let put = |path: Vec<i32>| -> i32 {
        if path.len() > out_cap as usize {
            return -2;
        }
        std::ptr::copy_nonoverlapping(path.as_ptr(), out, path.len());
        path.len() as i32
    };
    let mut inf = [1.0f64, 0.0, 0.0];
    let (res, e1) = astar::astar_plain(&p);
    inf[2] = e1 as f64;
    let r = match res {
        astar::Outcome::Found(path) => put(path),
        astar::Outcome::NotFound if e1 <= p.max_exp || budget >= max_exp => -1,
        astar::Outcome::NotFound => {
            let t = std::time::Instant::now();
            let (fv, cw, chw, f) = compute_field(&p);
            inf[1] = t.elapsed().as_secs_f64();
            let (fw, fhw) = (p.w, p.h * p.w);
            let look = |s: usize| {
                if f == 1 {
                    return fv[s];
                }
                let l = s / fhw;
                let c = s - l * fhw;
                fv[l * chw + (c / fw) / f * cw + (c % fw) / f]
            };
            if !p.src.iter().any(|&s| s >= 0 && (s as usize) < n && look(s as usize).is_finite()) {
                inf[0] = 3.0;
                -1
            } else {
                inf[0] = 2.0;
                p.max_exp = max_exp;
                let (res, e2) = astar::astar_guided_by(&p, look, weight);
                inf[2] += e2 as f64;
                match res {
                    astar::Outcome::Found(path) => put(path),
                    astar::Outcome::NotFound => -1,
                }
            }
        }
    };
    if !info.is_null() {
        std::ptr::copy_nonoverlapping(inf.as_ptr(), info, 3);
    }
    r
}

/// 128-bit content hash (out[2]) of the window rows i0..=i1, cols j0..=j1 of each of several [nl][ny][nx] arrays
/// (arrays[k], dims[k] = (nl, ny, nx, element size)), chained from `seed`. Windows are clipped to each array.
#[no_mangle]
pub unsafe extern "C" fn gr_hash_windows(
    narr: i32, arrays: *const *const u8, dims: *const i64, i0: i64, j0: i64, i1: i64, j1: i64, seed: u64,
    out: *mut u64,
) {
    let mut h = hash::Hasher::new(seed);
    for k in 0..narr as usize {
        let d = std::slice::from_raw_parts(dims.add(4 * k), 4);
        let (nl, ny, nx, es) = (d[0] as usize, d[1] as usize, d[2] as usize, d[3] as usize);
        let win = (i0.max(0) as usize, j0.max(0) as usize, (i1 as usize).min(ny - 1), (j1 as usize).min(nx - 1));
        hash::hash_window(*arrays.add(k), nl, ny, nx, es, win, &mut h);
    }
    let (a, b) = h.finish();
    *out = a;
    *out.add(1) = b;
}

/// 128-bit hash (out[2]) of several contiguous byte buffers (ptrs[k], lens[k]), in order; large buffers are hashed
/// in parallel 1 MiB chunks (the result does not depend on the thread count).
#[no_mangle]
pub unsafe extern "C" fn gr_hash_buffers(n: i32, ptrs: *const *const u8, lens: *const i64, seed: u64, out: *mut u64) {
    use rayon::prelude::*;
    const CH: usize = 1 << 20;
    let mut h = hash::Hasher::new(seed);
    for k in 0..n as usize {
        let b = std::slice::from_raw_parts(*ptrs.add(k), *lens.add(k) as usize);
        h.bytes(&(b.len() as u64).to_le_bytes());
        let parts: Vec<(u64, u64)> = b
            .par_chunks(CH)
            .map(|c| {
                let mut hc = hash::Hasher::new(seed ^ 0x5151);
                hc.bytes(c);
                hc.finish()
            })
            .collect();
        for (a, bb) in parts {
            h.bytes(&a.to_le_bytes());
            h.bytes(&bb.to_le_bytes());
        }
    }
    let (a, b) = h.finish();
    *out = a;
    *out.add(1) = b;
}

/// 128-bit hash (out[2]) of the listed 32x32 tiles (tiles: ntiles x (ty, tx)) of several [nl][h][w] byte arrays
/// (arrays[k], nls[k] planes; all h x w): the read-set fingerprint of a search.
#[no_mangle]
pub unsafe extern "C" fn gr_hash_tiles(
    narr: i32, arrays: *const *const u8, nls: *const i32, h: i32, w: i32, ntiles: i32, tiles: *const i32, out: *mut u64,
) {
    let (h, w) = (h as usize, w as usize);
    let t = astar::TOUCH_TILE;
    let mut hs = hash::Hasher::new(7);
    let tl = std::slice::from_raw_parts(tiles, 2 * ntiles.max(0) as usize);
    for k in 0..narr as usize {
        let a = *arrays.add(k);
        let nl = *nls.add(k) as usize;
        for q in tl.chunks_exact(2) {
            let (ty, tx) = (q[0] as usize, q[1] as usize);
            let (r0, c0) = (ty * t, tx * t);
            let (r1, c1) = ((r0 + t).min(h), (c0 + t).min(w));
            for l in 0..nl {
                for r in r0..r1 {
                    hs.bytes(std::slice::from_raw_parts(a.add((l * h + r) * w + c0), c1 - c0));
                }
            }
        }
    }
    let (x, y) = hs.finish();
    *out = x;
    *out.add(1) = y;
}

/// gr_island traced over window (i0, j0, h, w) but written only for the crop window (ci0, cj0, ch, cw), which
/// must lie inside it. out: nl * ch * cw bytes. Returns the number of crop cells set.
#[no_mangle]
pub unsafe extern "C" fn gr_island_crop(
    core: *const i16, nl: i32, ny: i32, nx: i32, thru: *const i16, nid: i16, i0: i32, j0: i32, h: i32, w: i32,
    seeds: *const i32, nseeds: i32, ci0: i32, cj0: i32, ch: i32, cw: i32, out: *mut u8,
) -> i64 {
    if i0 < 0 || j0 < 0 || h <= 0 || w <= 0 || i0 + h > ny || j0 + w > nx || ci0 < i0 || cj0 < j0 || ch <= 0 || cw <= 0
        || ci0 + ch > i0 + h || cj0 + cw > j0 + w
    {
        return -2;
    }
    let c = Grid3 { ptr: core, nl: nl as usize, ny: ny as usize, nx: nx as usize };
    let t = Grid3 { ptr: thru, nl: 1, ny: ny as usize, nx: nx as usize };
    let win = Window { i0: i0 as usize, j0: j0 as usize, h: h as usize, w: w as usize };
    let crop = Window { i0: ci0 as usize, j0: cj0 as usize, h: ch as usize, w: cw as usize };
    let sd: &[[i32; 3]] =
        if nseeds > 0 { std::slice::from_raw_parts(seeds as *const [i32; 3], nseeds as usize) } else { &[] };
    let o = std::slice::from_raw_parts_mut(out, nl as usize * crop.h * crop.w);
    flood::island_crop(&c, &t, nid, &win, sd, &crop, o);
    o.iter().map(|&v| v as i64).sum()
}

/// Indices of the set bytes of `a AND NOT b` (b may be NULL) over n bytes, in order, written to out (up to cap);
/// returns the count (which may exceed cap). With bbox (may be NULL) also its (x0, y0, x1, y1) over [.., h, w]
/// planes ((-1, ..) when empty). Parallel: the scan of a board-sized mask takes a few ms.
#[no_mangle]
pub unsafe extern "C" fn gr_mask_scan(
    a: *const u8, b: *const u8, n: i64, h: i32, w: i32, out: *mut i32, cap: i64, bbox: *mut i32,
) -> i64 {
    use rayon::prelude::*;
    let n = n as usize;
    let a = std::slice::from_raw_parts(a, n);
    let b = if b.is_null() { None } else { Some(std::slice::from_raw_parts(b, n)) };
    const CH: usize = 1 << 18;
    let parts: Vec<Vec<i32>> = a
        .par_chunks(CH)
        .enumerate()
        .map(|(k, c)| {
            let base = k * CH;
            let mut v = Vec::new();
            match b {
                Some(b) => {
                    let bc = &b[base..base + c.len()];
                    for (i, (&x, &y)) in c.iter().zip(bc).enumerate() {
                        if x != 0 && y == 0 {
                            v.push((base + i) as i32);
                        }
                    }
                }
                None => {
                    for (i, &x) in c.iter().enumerate() {
                        if x != 0 {
                            v.push((base + i) as i32);
                        }
                    }
                }
            }
            v
        })
        .collect();
    let total: usize = parts.iter().map(|p| p.len()).sum();
    let (h, w) = (h.max(1) as usize, w.max(1) as usize);
    let (mut x0, mut y0, mut x1, mut y1) = (i32::MAX, i32::MAX, -1, -1);
    let mut k = 0usize;
    for p in &parts {
        for &s in p {
            if !bbox.is_null() {
                let c = s as usize % (h * w);
                let (i, j) = ((c / w) as i32, (c % w) as i32);
                x0 = x0.min(j);
                x1 = x1.max(j);
                y0 = y0.min(i);
                y1 = y1.max(i);
            }
            if (k as i64) < cap {
                *out.add(k) = s;
            }
            k += 1;
        }
    }
    if !bbox.is_null() {
        let v = if total == 0 { [-1, -1, -1, -1] } else { [x0, y0, x1, y1] };
        std::ptr::copy_nonoverlapping(v.as_ptr(), bbox, 4);
    }
    total as i64
}

/// Per-tile 128-bit hashes: for each listed tile (tiles: ntiles x (ty, tx), `tile` cells per side) the hash of that
/// tile's rows in every array (arrays[k] with dims[k] = (nl, ny, nx, element size), clipped), written to
/// out[2 * k..2 * k + 2]. Tiles are hashed in parallel; used to keep incremental content hashes of a board.
#[no_mangle]
pub unsafe extern "C" fn gr_tile_hashes(
    narr: i32, arrays: *const *const u8, dims: *const i64, tile: i32, ntiles: i32, tiles: *const i32, out: *mut u64,
) {
    use rayon::prelude::*;
    let t = tile as usize;
    let tl = std::slice::from_raw_parts(tiles, 2 * ntiles.max(0) as usize);
    let arrs: Vec<(usize, [usize; 4])> = (0..narr as usize)
        .map(|k| {
            let d = std::slice::from_raw_parts(dims.add(4 * k), 4);
            (*arrays.add(k) as usize, [d[0] as usize, d[1] as usize, d[2] as usize, d[3] as usize])
        })
        .collect();
    let res: Vec<(u64, u64)> = tl
        .par_chunks_exact(2)
        .map(|q| {
            let (ty, tx) = (q[0] as usize, q[1] as usize);
            let mut h = hash::Hasher::new(11);
            for &(p, [nl, ny, nx, es]) in &arrs {
                let (i0, j0) = (ty * t, tx * t);
                if i0 >= ny || j0 >= nx {
                    continue;
                }
                let win = (i0, j0, (i0 + t - 1).min(ny - 1), (j0 + t - 1).min(nx - 1));
                hash::hash_window(p as *const u8, nl, ny, nx, es, win, &mut h);
            }
            h.finish()
        })
        .collect();
    for (k, (a, b)) in res.into_iter().enumerate() {
        *out.add(2 * k) = a;
        *out.add(2 * k + 1) = b;
    }
}

/// Paint a capsule (segment (x0, y0)-(x1, y1) of radius r, coordinates in the grid's units with pitch g) into plane
/// `layer` of arr [nl][ny][nx] int16: cells within r (centre to cell point, +1e-9) that are 0 take `nid`; with
/// `keep` set, cells already holding another value become `keep_val` (a clearance clash). The window is
/// (floor(min/g), ceil(max/g)) clipped, the same arithmetic as numpy, so the result is bit-identical.
/// win_out (may be NULL) receives the window (i0, j0, i1, j1). Returns the number of cells in the capsule.
#[no_mangle]
pub unsafe extern "C" fn gr_paint_capsule(
    arr: *mut i16, nl: i32, ny: i32, nx: i32, layer: i32, x0: f64, y0: f64, x1: f64, y1: f64, r: f64, g: f64,
    nid: i16, keep: i32, keep_val: i16, win_out: *mut i32,
) -> i64 {
    let (ny, nx) = (ny as i64, nx as i64);
    let i0 = ((y0.min(y1) - r) / g).floor().max(0.0) as i64;
    let j0 = ((x0.min(x1) - r) / g).floor().max(0.0) as i64;
    let i1 = (((y0.max(y1) + r) / g).ceil() as i64).min(ny - 1);
    let j1 = (((x0.max(x1) + r) / g).ceil() as i64).min(nx - 1);
    if !win_out.is_null() {
        std::ptr::copy_nonoverlapping([i0 as i32, j0 as i32, i1 as i32, j1 as i32].as_ptr(), win_out, 4);
    }
    if layer < 0 || layer >= nl {
        return 0;
    }
    let (vx, vy) = (x1 - x0, y1 - y0);
    let l2 = vx * vx + vy * vy;
    let rr = r * r + 1e-9;
    let base = arr.add((layer as i64 * ny * nx) as usize);
    let mut n = 0i64;
    for i in i0..=i1 {
        let ys = i as f64 * g;
        let row = base.add((i * nx) as usize);
        for j in j0..=j1 {
            let xs = j as f64 * g;
            let t = if l2 > 0.0 { (((xs - x0) * vx + (ys - y0) * vy) / l2).clamp(0.0, 1.0) } else { 0.0 };
            let dx = xs - (x0 + t * vx);
            let dy = ys - (y0 + t * vy);
            if dx * dx + dy * dy <= rr {
                n += 1;
                let c = row.add(j as usize);
                let v = *c;
                if v == 0 {
                    *c = nid;
                } else if keep != 0 && v != nid {
                    *c = keep_val;
                }
            }
        }
    }
    n
}

/// Copper component labels of net nid at ncells query cells ((layer, y, x) in grid coordinates): see
/// flood::label_cells. Cells with equal labels are connected; -1 = not the net's copper.
#[no_mangle]
pub unsafe extern "C" fn gr_label_cells(
    core: *const i16, nl: i32, ny: i32, nx: i32, thru: *const i16, nid: i16, cells: *const i32, ncells: i32,
    label: *mut i32,
) {
    let c = Grid3 { ptr: core, nl: nl as usize, ny: ny as usize, nx: nx as usize };
    let t = Grid3 { ptr: thru, nl: 1, ny: ny as usize, nx: nx as usize };
    let q = std::slice::from_raw_parts(cells as *const [i32; 3], ncells.max(0) as usize);
    flood::label_cells(&c, &t, nid, q, std::slice::from_raw_parts_mut(label, ncells.max(0) as usize));
}
