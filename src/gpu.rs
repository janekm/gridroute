//! Metal backend for the obstacle dilation (macOS). Shaders are compiled from source at first use.
//!
//! The source planes are wrapped without a copy when the touched rows lie in page-aligned memory (large numpy
//! arrays do); otherwise the touched rows are copied into a shared buffer.

use crate::dilate::{Element, Grid3, Window};
use metal::*;
use metal::foreign_types::ForeignType;
use objc::{sel, sel_impl};
use std::ffi::c_void;
use std::sync::{Mutex, OnceLock};

const SRC: &str = include_str!("../shaders/dilate.metal");
const FIELD_SRC: &str = include_str!("../shaders/field.metal");

#[repr(C)]
#[derive(Clone, Copy, Default)]
struct Params {
    nx: u32,
    ny: u32,
    i0: u32,
    j0: u32,
    h: u32,
    w: u32,
    r0: i32,
    nrows: u32,
    wmax: u32,
    nspan: u32,
    nlayers: u32,
    nexcl: u32,
    has_extra: u32,
    base_layer: u32,
    base_row: u32,
    plane_rows: u32,
    src_off: u32,
    extra_off: u32,
    hd_off: u32,
    out_off: u32,
    layer: [u32; 16],
    extra_on: [u32; 16],
    excl: [i16; 8],
}

struct Ctx {
    dev: Device,
    queue: CommandQueue,
    p_row: ComputePipelineState,
    p_comb: ComputePipelineState,
    p_finit: ComputePipelineState,
    p_fround: ComputePipelineState,
    p_freset: ComputePipelineState,
    p_fcollect: ComputePipelineState,
    dummy: Buffer,
}

unsafe impl Send for Ctx {}

static CTX: OnceLock<Option<Mutex<Ctx>>> = OnceLock::new();

fn ctx() -> Option<&'static Mutex<Ctx>> {
    CTX.get_or_init(|| {
        let dev = Device::system_default()?;
        let lib = dev.new_library_with_source(SRC, &CompileOptions::new()).ok()?;
        let f_row = lib.get_function("row_distance", None).ok()?;
        let f_comb = lib.get_function("combine", None).ok()?;
        let p_row = dev.new_compute_pipeline_state_with_function(&f_row).ok()?;
        let p_comb = dev.new_compute_pipeline_state_with_function(&f_comb).ok()?;
        let fsrc = FIELD_SRC.replace("constant constexpr int TS = 32;", &format!("constant constexpr int TS = {};", field_ts()));
        let flib = dev.new_library_with_source(&fsrc, &CompileOptions::new()).map_err(|e| eprintln!("gridroute: field.metal: {e}")).ok()?;
        let f_init = flib.get_function("field_init", None).ok()?;
        let f_round = flib.get_function("field_round", None).ok()?;
        let p_finit = dev.new_compute_pipeline_state_with_function(&f_init).ok()?;
        let p_fround = dev.new_compute_pipeline_state_with_function(&f_round).ok()?;
        let f_reset = flib.get_function("field_reset", None).ok()?;
        let f_collect = flib.get_function("field_collect", None).ok()?;
        let p_freset = dev.new_compute_pipeline_state_with_function(&f_reset).ok()?;
        let p_fcollect = dev.new_compute_pipeline_state_with_function(&f_collect).ok()?;
        let queue = dev.new_command_queue();
        let dummy = dev.new_buffer(16, MTLResourceOptions::StorageModeShared);
        Some(Mutex::new(Ctx { dev, queue, p_row, p_comb, p_finit, p_fround, p_freset, p_fcollect, dummy }))
    })
    .as_ref()
}

pub fn available() -> bool {
    ctx().is_some()
}

pub fn device_name() -> String {
    ctx().map(|c| c.lock().unwrap().dev.name().to_string()).unwrap_or_default()
}

const PAGE: usize = 16384;

/// A buffer covering `len` bytes at `ptr`: wrapped in place when possible, else copied. Returns the buffer and
/// the byte offset of `ptr` inside it.
fn wrap(dev: &Device, ptr: *const u8, len: usize) -> (Buffer, usize) {
    let a = ptr as usize;
    let lo = a & !(PAGE - 1);
    let hi = (a + len + PAGE - 1) & !(PAGE - 1);
    unsafe {
        let raw: *mut MTLBuffer = objc::msg_send![dev.as_ref(), newBufferWithBytesNoCopy: lo as *const c_void
            length: (hi - lo) as NSUInteger
            options: MTLResourceOptions::StorageModeShared
            deallocator: std::ptr::null::<c_void>()];
        if !raw.is_null() {
            return (Buffer::from_ptr(raw), a - lo);
        }
    }
    (dev.new_buffer_with_data(ptr as *const c_void, len as u64, MTLResourceOptions::StorageModeShared), 0)
}

fn size2(w: usize, h: usize) -> MTLSize {
    MTLSize { width: w as u64, height: h as u64, depth: 1 }
}

/// Same contract as `dilate::dilate_cpu`. Returns false when Metal is unavailable or the request does not fit
/// the shader limits (the caller then uses the CPU).
pub fn dilate_metal(
    src: &Grid3<i16>,
    layers: &[usize],
    excl: &[i16],
    extra: Option<&Grid3<u8>>,
    extra_on: &[bool],
    win: &Window,
    el: &Element,
    reduce_or: bool,
    out: &mut [u8],
) -> bool {
    let Some(c) = ctx() else { return false };
    if layers.len() > 16 || excl.len() > 8 || el.wmax() >= 255 || layers.is_empty() {
        return false;
    }
    let c = c.lock().unwrap();
    let rc = el.rc();
    let r0 = win.i0 as isize - rc as isize;
    let rs0 = r0.max(0) as usize;
    let rs1 = ((win.i0 + win.h - 1 + rc).min(src.ny - 1)) as usize;
    let lmin = *layers.iter().min().unwrap();
    let lmax = *layers.iter().max().unwrap();
    // source bytes from (lmin, rs0, 0) to the end of (lmax, rs1, nx - 1)
    let first = (lmin * src.ny + rs0) * src.nx;
    let last = (lmax * src.ny + rs1 + 1) * src.nx;
    let (sbuf, soff) = wrap(&c.dev, unsafe { src.ptr.add(first) } as *const u8, (last - first) * 2);
    let (ebuf, eoff) = match extra {
        Some(e) => wrap(&c.dev, unsafe { e.ptr.add(rs0 * e.nx) }, (rs1 + 1 - rs0) * e.nx),
        None => (c.dummy.clone(), 0),
    };
    let nrows = win.h + 2 * rc;
    let hw = win.h * win.w;
    let groups: Vec<Vec<usize>> = if reduce_or { vec![layers.to_vec()] } else { layers.iter().map(|&l| vec![l]).collect() };
    let hd = c.dev.new_buffer((groups.len() * nrows * win.w).max(1) as u64, MTLResourceOptions::StorageModePrivate);
    let obuf = c.dev.new_buffer((groups.len() * hw).max(1) as u64, MTLResourceOptions::StorageModeShared);
    let span: Vec<i32> = el.span.to_vec();
    let cmd = c.queue.new_command_buffer();
    for (g, ls) in groups.iter().enumerate() {
        let mut p = Params {
            nx: src.nx as u32,
            ny: src.ny as u32,
            i0: win.i0 as u32,
            j0: win.j0 as u32,
            h: win.h as u32,
            w: win.w as u32,
            r0: r0 as i32,
            nrows: nrows as u32,
            wmax: el.wmax() as u32,
            nspan: span.len() as u32,
            nlayers: ls.len() as u32,
            nexcl: excl.len() as u32,
            has_extra: extra.is_some() as u32,
            base_layer: lmin as u32,
            base_row: rs0 as u32,
            plane_rows: src.ny as u32,
            src_off: (soff / 2) as u32,
            extra_off: eoff as u32,
            hd_off: (g * nrows * win.w) as u32,
            out_off: (g * hw) as u32,
            ..Default::default()
        };
        for (n, &l) in ls.iter().enumerate() {
            p.layer[n] = l as u32;
            p.extra_on[n] = extra_on.get(l).copied().unwrap_or(false) as u32;
        }
        p.excl[..excl.len()].copy_from_slice(excl);
        let pp = &p as *const Params as *const c_void;
        let psz = std::mem::size_of::<Params>() as u64;
        let enc = cmd.new_compute_command_encoder();
        enc.set_compute_pipeline_state(&c.p_row);
        enc.set_buffer(0, Some(&sbuf), 0);
        enc.set_buffer(1, Some(&ebuf), 0);
        enc.set_bytes(2, psz, pp);
        enc.set_buffer(3, Some(&hd), 0);
        enc.dispatch_threads(size2(win.w, nrows), size2(32, 8));
        enc.set_compute_pipeline_state(&c.p_comb);
        enc.set_buffer(0, Some(&hd), 0);
        enc.set_bytes(1, psz, pp);
        enc.set_bytes(2, (span.len() * 4) as u64, span.as_ptr() as *const c_void);
        enc.set_buffer(3, Some(&obuf), 0);
        enc.dispatch_threads(size2(win.w, win.h), size2(32, 8));
        enc.end_encoding();
    }
    cmd.commit();
    cmd.wait_until_completed();
    unsafe {
        std::ptr::copy_nonoverlapping(obuf.contents() as *const u8, out.as_mut_ptr(), groups.len() * hw);
    }
    true
}

#[repr(C)]
#[derive(Clone, Copy)]
struct FieldParams {
    nl: u32,
    h: u32,
    w: u32,
    tiles_x: u32,
    tiles_y: u32,
    lay_ok: u32,
    vcost: f32,
    maxit: u32,
    round: u32,
}

const MAX_ROUNDS: usize = 1 << 16;

fn env_usize(k: &str, d: usize) -> usize {
    std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d)
}

/// Field tile size (cells per side), local relaxation iterations per round, rounds per command buffer.
fn field_ts() -> usize {
    env_usize("GR_FIELD_TS", 8)
}

/// Cost-to-target field (field.rs definition) by tiled relaxation on the GPU. Returns the number of rounds, or
/// None when Metal is unavailable or the problem has more than 6 layers.
pub fn field_metal(p: &crate::astar::Problem, out: &mut [f32]) -> Option<usize> {
    let c = ctx()?;
    if p.nl > 6 {
        return None;
    }
    let c = c.lock().unwrap();
    let n = p.nl * p.h * p.w;
    let opt = MTLResourceOptions::StorageModeShared;
    let (blk, boff) = wrap(&c.dev, p.blk.as_ptr(), n);
    let (tgt, toff) = wrap(&c.dev, p.tgt.as_ptr(), n);
    let (vok, voff) = wrap(&c.dev, p.vok.as_ptr(), p.h * p.w);
    if boff != 0 || toff != 0 || voff != 0 {
        // the kernels index from 0: fall back to tight copies
        return field_metal_copied(&c, p, out);
    }
    run_field(&c, p, &blk, &tgt, &vok, out, opt)
}

fn field_metal_copied(c: &Ctx, p: &crate::astar::Problem, out: &mut [f32]) -> Option<usize> {
    let opt = MTLResourceOptions::StorageModeShared;
    let mk = |s: &[u8]| c.dev.new_buffer_with_data(s.as_ptr() as *const c_void, s.len().max(1) as u64, opt);
    let (blk, tgt, vok) = (mk(p.blk), mk(p.tgt), mk(p.vok));
    run_field(c, p, &blk, &tgt, &vok, out, opt)
}

fn run_field(c: &Ctx, p: &crate::astar::Problem, blk: &Buffer, tgt: &Buffer, vok: &Buffer, out: &mut [f32],
             opt: MTLResourceOptions) -> Option<usize> {
    let n = p.nl * p.h * p.w;
    let fbuf = c.dev.new_buffer((n * 4) as u64, opt);
    let ts = field_ts();
    let (tx, ty) = (p.w.div_ceil(ts), p.h.div_ceil(ts));
    let last = c.dev.new_buffer((tx * ty * 4) as u64, opt);
    let changed = c.dev.new_buffer((MAX_ROUNDS * 4) as u64, opt);
    let list = c.dev.new_buffer((tx * ty * 4) as u64, MTLResourceOptions::StorageModePrivate);
    let args = c.dev.new_buffer(16, opt);
    unsafe {
        let a = args.contents() as *mut u32;
        *a = 0;
        *a.add(1) = 1;
        *a.add(2) = 1;
    }
    unsafe {
        std::ptr::write_bytes(last.contents() as *mut u8, 0, tx * ty * 4);
        std::ptr::write_bytes(changed.contents() as *mut u8, 0, MAX_ROUNDS * 4);
    }
    let mut lay = 0u32;
    for (l, &k) in p.lay_ok.iter().enumerate() {
        if k != 0 {
            lay |= 1 << l;
        }
    }
    let mut fp = FieldParams {
        nl: p.nl as u32,
        h: p.h as u32,
        w: p.w as u32,
        tiles_x: tx as u32,
        tiles_y: ty as u32,
        lay_ok: lay,
        vcost: p.vcost,
        maxit: env_usize("GR_FIELD_MAXIT", 32) as u32,
        round: 0,
    };
    let psz = std::mem::size_of::<FieldParams>() as u64;
    let nn = n as u32;
    let mut round = 0usize;
    let batch = env_usize("GR_FIELD_BATCH", 24);
    {
        let cmd = c.queue.new_command_buffer();
        let enc = cmd.new_compute_command_encoder();
        enc.set_compute_pipeline_state(&c.p_finit);
        enc.set_buffer(0, Some(&fbuf), 0);
        enc.set_buffer(1, Some(blk), 0);
        enc.set_buffer(2, Some(tgt), 0);
        enc.set_bytes(3, 4, &nn as *const u32 as *const c_void);
        enc.dispatch_threads(MTLSize { width: n as u64, height: 1, depth: 1 }, MTLSize { width: 256, height: 1, depth: 1 });
        enc.end_encoding();
        cmd.commit();
        cmd.wait_until_completed();
    }
    let mut src_reached_at: Option<usize> = None;
    loop {
        let cmd = c.queue.new_command_buffer();
        let enc = cmd.new_compute_command_encoder();
        let end = (round + batch).min(MAX_ROUNDS);
        let ntiles = (tx * ty) as u64;
        for r in round..end {
            fp.round = r as u32;
            let pp = &fp as *const FieldParams as *const c_void;
            enc.set_compute_pipeline_state(&c.p_freset);
            enc.set_buffer(0, Some(&args), 0);
            enc.dispatch_threads(MTLSize { width: 1, height: 1, depth: 1 }, MTLSize { width: 1, height: 1, depth: 1 });
            enc.memory_barrier_with_resources(&[&args]);
            enc.set_compute_pipeline_state(&c.p_fcollect);
            enc.set_bytes(0, psz, pp);
            enc.set_buffer(1, Some(&last), 0);
            enc.set_buffer(2, Some(&args), 0);
            enc.set_buffer(3, Some(&list), 0);
            enc.dispatch_threads(MTLSize { width: ntiles, height: 1, depth: 1 }, MTLSize { width: 64, height: 1, depth: 1 });
            enc.memory_barrier_with_resources(&[&args, &list]);
            enc.set_compute_pipeline_state(&c.p_fround);
            enc.set_buffer(0, Some(&fbuf), 0);
            enc.set_buffer(1, Some(blk), 0);
            enc.set_buffer(2, Some(vok), 0);
            enc.set_bytes(3, (p.mcost.len() * 4) as u64, p.mcost.as_ptr() as *const c_void);
            enc.set_bytes(4, psz, pp);
            enc.set_buffer(5, Some(&last), 0);
            enc.set_buffer(6, Some(&changed), 0);
            enc.set_buffer(7, Some(&list), 0);
            enc.dispatch_thread_groups_indirect(&args, 0, MTLSize { width: ts as u64, height: ts as u64, depth: 1 });
            enc.memory_barrier_with_resources(&[&fbuf, &last, &changed]);
        }
        enc.end_encoding();
        cmd.commit();
        cmd.wait_until_completed();
        let ch = unsafe { std::slice::from_raw_parts(changed.contents() as *const u32, MAX_ROUNDS) };
        let mut done = ch[end - 1] == 0;
        // early stop: once a source state has a finite cost, `extra` more batches (the wave has reached the
        // source; values behind it are still settling, so the field is then an approximate heuristic)
        if !done && p.src.len() > 0 {
            if let Some(extra) = early_extra() {
                if src_reached_at.is_none() {
                    let f = unsafe { std::slice::from_raw_parts(fbuf.contents() as *const f32, n) };
                    if p.src.iter().any(|&s| s >= 0 && (s as usize) < n && f[s as usize].is_finite()) {
                        src_reached_at = Some(round);
                    }
                }
                if let Some(r0) = src_reached_at {
                    if end >= r0 + batch * (extra + 1) {
                        done = true;
                    }
                }
            }
        }
        round = end;
        if done || round >= MAX_ROUNDS {
            // trim: the first round without change
            round = (0..end).find(|&r| ch[r] == 0).map_or(end, |r| r + 1);
            let _ = &mut round;
            break;
        }
    }
    unsafe {
        std::ptr::copy_nonoverlapping(fbuf.contents() as *const f32, out.as_mut_ptr(), n);
    }
    Some(round)
}

/// Batches (of 24 rounds) to keep relaxing after the wave reaches a source; None = run to convergence.
/// Set with gr_set_field_early_stop (or $GR_FIELD_EXTRA for experiments).
static EARLY_EXTRA: std::sync::atomic::AtomicI64 = std::sync::atomic::AtomicI64::new(-2);

pub fn set_early_extra(v: i64) {
    EARLY_EXTRA.store(v, std::sync::atomic::Ordering::Relaxed);
}

fn early_extra() -> Option<usize> {
    let mut v = EARLY_EXTRA.load(std::sync::atomic::Ordering::Relaxed);
    if v == -2 {
        v = std::env::var("GR_FIELD_EXTRA").ok().and_then(|s| s.parse().ok()).unwrap_or(-1);
        EARLY_EXTRA.store(v, std::sync::atomic::Ordering::Relaxed);
    }
    if v < 0 { None } else { Some(v as usize) }
}
