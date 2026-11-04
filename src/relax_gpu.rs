//! Metal backend of the elastic relaxation (src/relax.rs). Shaders are compiled from source at first use.

use crate::relax::{Params, Relax};
use metal::*;
use std::ffi::c_void;
use std::sync::{Mutex, OnceLock};

const SRC: &str = include_str!("../shaders/relax.metal");

#[repr(C)]
#[derive(Clone, Copy)]
struct GpuParams {
    iters: u32,
    sub: u32,
    alpha_wire: f32,
    alpha_via: f32,
    beta: f32,
    dmax: f32,
    tol: f32,
    n: u32,
}

struct Ctx {
    dev: Device,
    queue: CommandQueue,
    p_tension: ComputePipelineState,
    p_project: ComputePipelineState,
}

unsafe impl Send for Ctx {}

static CTX: OnceLock<Option<Mutex<Ctx>>> = OnceLock::new();

fn ctx() -> Option<&'static Mutex<Ctx>> {
    CTX.get_or_init(|| {
        let dev = Device::system_default()?;
        let lib = dev.new_library_with_source(SRC, &CompileOptions::new()).map_err(|e| eprintln!("gridroute: relax.metal: {e}")).ok()?;
        let f_tension = lib.get_function("relax_tension", None).ok()?;
        let f_project = lib.get_function("relax_project", None).ok()?;
        let p_tension = dev.new_compute_pipeline_state_with_function(&f_tension).ok()?;
        let p_project = dev.new_compute_pipeline_state_with_function(&f_project).ok()?;
        let queue = dev.new_command_queue();
        Some(Mutex::new(Ctx { dev, queue, p_tension, p_project }))
    })
    .as_ref()
}

fn buf<T>(dev: &Device, v: &[T]) -> Buffer {
    let len = (std::mem::size_of_val(v)).max(16) as u64;
    if v.is_empty() {
        return dev.new_buffer(len, MTLResourceOptions::StorageModeShared);
    }
    dev.new_buffer_with_data(v.as_ptr() as *const c_void, len, MTLResourceOptions::StorageModeShared)
}

/// Buffers kept on the GPU between runs: the fixed model once, the constraint lists per version.
pub struct State {
    pos: [Buffer; 2],
    kind: Buffer,
    tn_off: Buffer,
    tn: Buffer,
    leash: Buffer,
    el_off: Buffer,
    el_node: Buffer,
    ne_off: Buffer,
    ne: Buffer,
    scale: Buffer,
    viol: Buffer,
    c: Option<(u64, Buffer, Buffer, Buffer)>,
}

unsafe impl Send for State {}

/// Run `prm.iters` iterations; false when Metal is unavailable (the caller then uses the CPU).
pub fn run(r: &mut Relax, prm: &Params) -> bool {
    let Some(c) = ctx() else { return false };
    let c = c.lock().unwrap();
    let n = r.pos.len();
    if n == 0 {
        return true;
    }
    let gp = GpuParams { iters: prm.iters, sub: prm.sub, alpha_wire: prm.alpha_wire, alpha_via: prm.alpha_via, beta: prm.beta, dmax: prm.dmax, tol: prm.tol, n: n as u32 };
    if r.gpu.is_none() {
        r.init_scale();
        r.gpu = Some(State {
            scale: buf(&c.dev, &r.scale),
            viol: buf(&c.dev, &r.viol),
            pos: [buf(&c.dev, &r.pos), buf(&c.dev, &r.pos)],
            kind: buf(&c.dev, &r.kind),
            tn_off: buf(&c.dev, &r.tn_off),
            tn: buf(&c.dev, &r.tn),
            leash: buf(&c.dev, &r.leash),
            el_off: buf(&c.dev, &r.el_off),
            el_node: buf(&c.dev, &r.el_node),
            ne_off: buf(&c.dev, &r.ne_off),
            ne: buf(&c.dev, &r.ne),
            c: None,
        });
    }
    let st = r.gpu.as_mut().unwrap();
    if st.c.as_ref().map_or(true, |x| x.0 != r.c_version) {
        st.c = Some((r.c_version, buf(&c.dev, &r.c_off), buf(&c.dev, &r.c_el), buf(&c.dev, &r.c_need)));
    }
    // positions: the caller may have moved them (a batch taken back)
    unsafe {
        std::ptr::copy_nonoverlapping(r.pos.as_ptr(), st.pos[0].contents() as *mut [f32; 2], n);
    }
    let (_, c_off, c_el, c_need) = st.c.as_ref().unwrap();
    let psz = std::mem::size_of::<GpuParams>() as u64;
    let pp = &gp as *const GpuParams as *const c_void;
    let grid = MTLSize { width: n as u64, height: 1, depth: 1 };
    let tg = MTLSize { width: 128, height: 1, depth: 1 };
    let cmd = c.queue.new_command_buffer();
    let enc = cmd.new_compute_command_encoder();
    // ping-pong: `cur` holds the positions after each dispatch
    let mut cur = 0usize;
    for _ in 0..prm.iters {
        enc.set_compute_pipeline_state(&c.p_tension);
        enc.set_buffer(0, Some(&st.pos[cur]), 0);
        enc.set_buffer(1, Some(&st.pos[1 - cur]), 0);
        enc.set_buffer(2, Some(&st.kind), 0);
        enc.set_buffer(3, Some(&st.tn_off), 0);
        enc.set_buffer(4, Some(&st.tn), 0);
        enc.set_bytes(5, psz, pp);
        enc.set_buffer(6, Some(&st.scale), 0);
        enc.set_buffer(7, Some(&st.viol), 0);
        enc.dispatch_threads(grid, tg);
        enc.memory_barrier_with_resources(&[&st.pos[1 - cur]]);
        cur = 1 - cur;
        for _ in 0..prm.sub {
            enc.set_compute_pipeline_state(&c.p_project);
            enc.set_buffer(0, Some(&st.pos[cur]), 0);
            enc.set_buffer(1, Some(&st.pos[1 - cur]), 0);
            enc.set_buffer(2, Some(&st.kind), 0);
            enc.set_buffer(3, Some(&st.leash), 0);
            enc.set_buffer(4, Some(&st.el_off), 0);
            enc.set_buffer(5, Some(&st.el_node), 0);
            enc.set_buffer(6, Some(&st.ne_off), 0);
            enc.set_buffer(7, Some(&st.ne), 0);
            enc.set_buffer(8, Some(c_off), 0);
            enc.set_buffer(9, Some(c_el), 0);
            enc.set_buffer(10, Some(c_need), 0);
            enc.set_bytes(11, psz, pp);
            enc.set_buffer(12, Some(&st.viol), 0);
            enc.dispatch_threads(grid, tg);
            enc.memory_barrier_with_resources(&[&st.pos[1 - cur]]);
            cur = 1 - cur;
        }
    }
    enc.end_encoding();
    cmd.commit();
    cmd.wait_until_completed();
    unsafe {
        std::ptr::copy_nonoverlapping(st.pos[cur].contents() as *const [f32; 2], r.pos.as_mut_ptr(), n);
    }
    true
}
