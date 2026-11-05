//! Metal backend: the CPU reference updates (`fdtd.rs`) as compute kernels (`shaders/fdtd.metal`). Steps are
//! batched into command buffers; probes are written to a per-batch buffer and read back after each batch.

use crate::fdtd::{FieldDft, RunResult, RunSpec};
use crate::model::{Problem, MU0};
use metal::*;
use std::ffi::c_void;

const SRC: &str = include_str!("../../../shaders/fdtd.metal");

#[repr(C)]
#[derive(Clone, Copy, Default)]
struct P {
    nx: u32,
    ny: u32,
    nz: u32,
    kmu: f32,
    axis: u32,
    slab_len: u32,
    step: u32,
    vs: f32,
    nsrc: u32,
    nprobe: u32,
    w: [f32; 12],
}

struct Pipes {
    h: ComputePipelineState,
    e: ComputePipelineState,
    ph: ComputePipelineState,
    pe: ComputePipelineState,
    src: ComputePipelineState,
    probe: ComputePipelineState,
    dft: ComputePipelineState,
    energy: ComputePipelineState,
}

fn buf_f32(dev: &Device, v: &[f32]) -> Buffer {
    let len = (v.len().max(1) * 4) as u64;
    if v.is_empty() {
        return dev.new_buffer(len, MTLResourceOptions::StorageModeShared);
    }
    dev.new_buffer_with_data(v.as_ptr() as *const c_void, len, MTLResourceOptions::StorageModeShared)
}
fn buf_u32(dev: &Device, v: &[u32]) -> Buffer {
    let len = (v.len().max(1) * 4) as u64;
    if v.is_empty() {
        return dev.new_buffer(len, MTLResourceOptions::StorageModeShared);
    }
    dev.new_buffer_with_data(v.as_ptr() as *const c_void, len, MTLResourceOptions::StorageModeShared)
}
fn zeros(dev: &Device, n: usize) -> Buffer {
    let b = dev.new_buffer((n.max(1) * 4) as u64, MTLResourceOptions::StorageModeShared);
    unsafe { std::ptr::write_bytes(b.contents() as *mut u8, 0, n.max(1) * 4) };
    b
}
fn read(b: &Buffer, n: usize) -> Vec<f32> {
    unsafe { std::slice::from_raw_parts(b.contents() as *const f32, n).to_vec() }
}

pub fn available() -> bool {
    Device::system_default().is_some()
}

pub fn run_metal(p: &Problem, spec: &RunSpec) -> Result<RunResult, String> {
    let t0 = std::time::Instant::now();
    let dev = Device::system_default().ok_or("no Metal device")?;
    let lib = dev.new_library_with_source(SRC, &CompileOptions::new()).map_err(|e| format!("fdtd.metal: {e}"))?;
    let pipe = |name: &str| -> Result<ComputePipelineState, String> {
        let f = lib.get_function(name, None).map_err(|e| e.to_string())?;
        dev.new_compute_pipeline_state_with_function(&f).map_err(|e| e.to_string())
    };
    let pp = Pipes { h: pipe("update_h")?, e: pipe("update_e")?, ph: pipe("pml_h")?, pe: pipe("pml_e")?, src: pipe("add_src")?, probe: pipe("probe")?, dft: pipe("dft")?, energy: pipe("energy")? };
    let queue = dev.new_command_queue();
    let (nx, ny, nz) = (p.nx, p.ny, p.nz);
    let n = p.n();
    let e: Vec<Buffer> = (0..3).map(|_| zeros(&dev, n)).collect();
    let h: Vec<Buffer> = (0..3).map(|_| zeros(&dev, n)).collect();
    let ca: Vec<Buffer> = p.ca.iter().map(|v| buf_f32(&dev, v)).collect();
    let cb: Vec<Buffer> = p.cb.iter().map(|v| buf_f32(&dev, v)).collect();
    let cat = |a: &[Vec<f32>; 3]| -> Vec<f32> { a.iter().flatten().copied().collect() };
    let inv_d = buf_f32(&dev, &cat(&p.inv_d));
    let inv_dd = buf_f32(&dev, &cat(&p.inv_dd));
    let nn = [nx, ny, nz];
    struct Pml {
        slab: Buffer,
        be: Buffer,
        ce: Buffer,
        bh: Buffer,
        ch: Buffer,
        psi: Vec<Buffer>,
        len: usize,
    }
    let pml: Vec<Pml> = (0..3)
        .map(|a| {
            let ax = &p.pml[a];
            let len = ax.len();
            let psize = len * n / nn[a];
            Pml {
                slab: buf_u32(&dev, &ax.slab),
                be: buf_f32(&dev, &ax.be),
                ce: buf_f32(&dev, &ax.ce),
                bh: buf_f32(&dev, &ax.bh),
                ch: buf_f32(&dev, &ax.ch),
                psi: (0..4).map(|_| zeros(&dev, psize)).collect(),
                len,
            }
        })
        .collect();
    // sources
    let (sc, si, sg): (Vec<u32>, Vec<u32>, Vec<f32>) = match spec.excite {
        Some(k) => (p.ports[k].src.iter().map(|s| s.0 as u32).collect(), p.ports[k].src.iter().map(|s| s.1).collect(), p.ports[k].src.iter().map(|s| s.2).collect()),
        None => (vec![], vec![], vec![]),
    };
    let nsrc = sc.len();
    let (sc, si, sg) = (buf_u32(&dev, &sc), buf_u32(&dev, &si), buf_f32(&dev, &sg));
    // probes: ports first, then field probes
    let mut start = vec![0u32];
    let (mut tc, mut ti, mut tw) = (Vec::new(), Vec::new(), Vec::new());
    for pt in &p.ports {
        for &(c, id, w) in &pt.probe {
            tc.push(c as u32);
            ti.push(id);
            tw.push(w);
        }
        start.push(tc.len() as u32);
    }
    for pr in &p.probes {
        for &(c, id, w) in &pr.terms {
            tc.push(c as u32);
            ti.push(id);
            tw.push(w);
        }
        start.push(tc.len() as u32);
    }
    let nprobe = p.ports.len() + p.probes.len();
    let (startb, tcb, tib, twb) = (buf_u32(&dev, &start), buf_u32(&dev, &tc), buf_u32(&dev, &ti), buf_f32(&dev, &tw));
    const BATCH: usize = 256;
    let pout = zeros(&dev, BATCH * nprobe.max(1));
    let eout = zeros(&dev, ny * nz);
    let dftb: Vec<Vec<Buffer>> = spec.dft_freqs.iter().map(|_| (0..6).map(|_| zeros(&dev, 2 * n)).collect()).collect();

    let prm = P { nx: nx as u32, ny: ny as u32, nz: nz as u32, kmu: (p.dt / MU0) as f32, nsrc: nsrc as u32, nprobe: nprobe as u32, ..Default::default() };
    let psz = std::mem::size_of::<P>() as u64;
    let tg3 = MTLSize { width: 32, height: 4, depth: 1 };
    let grid3 = MTLSize { width: nx as u64, height: ny as u64, depth: nz as u64 };
    let mut v: Vec<Vec<f32>> = vec![vec![0.0]; p.ports.len()];
    let mut pv: Vec<Vec<f32>> = vec![vec![0.0]; p.probes.len()];
    let mut peak = 0.0f64;
    let mut steps = 0usize;
    let mut done = false;
    let mut snaps = Vec::new();
    while !done && steps < spec.max_steps {
        let mut nb = BATCH.min(spec.max_steps - steps);
        if spec.snap_every > 0 && (spec.snap_until == 0 || steps < spec.snap_until) {
            nb = nb.min(spec.snap_every - steps % spec.snap_every);
        }
        let cmd = queue.new_command_buffer();
        let enc = cmd.new_compute_command_encoder();
        for b in 0..nb {
            let n_step = steps + b;
            let set_p = |prm: &P| enc.set_bytes(7, psz, prm as *const P as *const c_void);
            // H
            enc.set_compute_pipeline_state(&pp.h);
            for c in 0..3 {
                enc.set_buffer(c as u64, Some(&h[c]), 0);
                enc.set_buffer(3 + c as u64, Some(&e[c]), 0);
            }
            enc.set_buffer(6, Some(&inv_d), 0);
            set_p(&prm);
            enc.dispatch_threads(grid3, tg3);
            for a in 0..3 {
                let pm = &pml[a];
                if pm.len == 0 {
                    continue;
                }
                enc.set_compute_pipeline_state(&pp.ph);
                // f0..f2 = H comps (a, a+1, a+2 indexed absolutely), e0..e2 = E comps
                for c in 0..3 {
                    enc.set_buffer(c as u64, Some(&h[c]), 0);
                    enc.set_buffer(3 + c as u64, Some(&e[c]), 0);
                }
                enc.set_buffer(6, Some(&inv_d), 0);
                enc.set_buffer(8, Some(&pm.slab), 0);
                enc.set_buffer(9, Some(&pm.bh), 0);
                enc.set_buffer(10, Some(&pm.ch), 0);
                enc.set_buffer(11, Some(&pm.psi[2]), 0);
                enc.set_buffer(12, Some(&pm.psi[3]), 0);
                let mut q = prm;
                q.axis = a as u32;
                q.slab_len = pm.len as u32;
                set_p(&q);
                let mut g = [nx as u64, ny as u64, nz as u64];
                g[a] = pm.len as u64;
                enc.dispatch_threads(MTLSize { width: g[0], height: g[1], depth: g[2] }, tg3);
            }
            // E
            enc.set_compute_pipeline_state(&pp.e);
            for c in 0..3 {
                enc.set_buffer(c as u64, Some(&e[c]), 0);
                enc.set_buffer(3 + c as u64, Some(&h[c]), 0);
                enc.set_buffer(8 + c as u64, Some(&ca[c]), 0);
                enc.set_buffer(11 + c as u64, Some(&cb[c]), 0);
            }
            enc.set_buffer(6, Some(&inv_dd), 0);
            set_p(&prm);
            enc.dispatch_threads(grid3, tg3);
            for a in 0..3 {
                let pm = &pml[a];
                if pm.len == 0 {
                    continue;
                }
                enc.set_compute_pipeline_state(&pp.pe);
                for c in 0..3 {
                    enc.set_buffer(c as u64, Some(&e[c]), 0);
                    enc.set_buffer(3 + c as u64, Some(&h[c]), 0);
                    enc.set_buffer(13 + c as u64, Some(&cb[c]), 0);
                }
                enc.set_buffer(6, Some(&inv_dd), 0);
                enc.set_buffer(8, Some(&pm.slab), 0);
                enc.set_buffer(9, Some(&pm.be), 0);
                enc.set_buffer(10, Some(&pm.ce), 0);
                enc.set_buffer(11, Some(&pm.psi[0]), 0);
                enc.set_buffer(12, Some(&pm.psi[1]), 0);
                let mut q = prm;
                q.axis = a as u32;
                q.slab_len = pm.len as u32;
                set_p(&q);
                let mut g = [nx as u64, ny as u64, nz as u64];
                g[a] = pm.len as u64;
                enc.dispatch_threads(MTLSize { width: g[0], height: g[1], depth: g[2] }, tg3);
            }
            // source
            if nsrc > 0 && n_step < spec.exc.len() {
                enc.set_compute_pipeline_state(&pp.src);
                for c in 0..3 {
                    enc.set_buffer(c as u64, Some(&e[c]), 0);
                }
                enc.set_buffer(3, Some(&sc), 0);
                enc.set_buffer(4, Some(&si), 0);
                enc.set_buffer(5, Some(&sg), 0);
                let mut q = prm;
                q.vs = spec.exc[n_step];
                set_p(&q);
                enc.dispatch_threads(MTLSize { width: nsrc as u64, height: 1, depth: 1 }, MTLSize { width: 64.min(nsrc as u64).max(1), height: 1, depth: 1 });
            }
            // probes
            if nprobe > 0 {
                enc.set_compute_pipeline_state(&pp.probe);
                for c in 0..3 {
                    enc.set_buffer(c as u64, Some(&e[c]), 0);
                    enc.set_buffer(3 + c as u64, Some(&h[c]), 0);
                }
                enc.set_buffer(6, Some(&startb), 0);
                enc.set_buffer(8, Some(&tcb), 0);
                enc.set_buffer(9, Some(&tib), 0);
                enc.set_buffer(10, Some(&twb), 0);
                enc.set_buffer(11, Some(&pout), 0);
                let mut q = prm;
                q.step = b as u32;
                set_p(&q);
                enc.dispatch_threads(MTLSize { width: nprobe as u64, height: 1, depth: 1 }, MTLSize { width: 32.min(nprobe as u64), height: 1, depth: 1 });
            }
            // field DFT
            let s1 = n_step + 1;
            if !spec.dft_freqs.is_empty() && s1 % spec.dft_every == 0 {
                enc.set_compute_pipeline_state(&pp.dft);
                for (fi, &f) in spec.dft_freqs.iter().enumerate() {
                    let w = 2.0 * std::f64::consts::PI * f;
                    let mut q = prm;
                    for c in 0..6 {
                        let t = if c < 3 { s1 as f64 * p.dt } else { (s1 as f64 - 0.5) * p.dt };
                        let (sn, cs) = (-w * t).sin_cos();
                        q.w[2 * c] = (cs * p.dt * spec.dft_every as f64) as f32;
                        q.w[2 * c + 1] = (sn * p.dt * spec.dft_every as f64) as f32;
                    }
                    for c in 0..6 {
                        q.axis = c as u32;
                        enc.set_buffer(0, Some(if c < 3 { &e[c] } else { &h[c - 3] }), 0);
                        enc.set_buffer(1, Some(&dftb[fi][c]), 0);
                        set_p(&q);
                        enc.dispatch_threads(MTLSize { width: n as u64, height: 1, depth: 1 }, MTLSize { width: 256, height: 1, depth: 1 });
                    }
                }
            }
        }
        let end = steps + nb;
        let check = end / spec.check_every > steps / spec.check_every;
        if check {
            enc.set_compute_pipeline_state(&pp.energy);
            for c in 0..3 {
                enc.set_buffer(c as u64, Some(&e[c]), 0);
                enc.set_buffer(3 + c as u64, Some(&h[c]), 0);
            }
            enc.set_buffer(6, Some(&eout), 0);
            enc.set_bytes(7, psz, &prm as *const P as *const c_void);
            enc.dispatch_threads(MTLSize { width: (ny * nz) as u64, height: 1, depth: 1 }, MTLSize { width: 64, height: 1, depth: 1 });
        }
        enc.end_encoding();
        cmd.commit();
        cmd.wait_until_completed();
        let pr = read(&pout, nb * nprobe);
        for b in 0..nb {
            for k in 0..p.ports.len() {
                v[k].push(pr[b * nprobe + k]);
            }
            for k in 0..p.probes.len() {
                pv[k].push(pr[b * nprobe + p.ports.len() + k]);
            }
        }
        steps = end;
        if spec.snap_every > 0 && steps % spec.snap_every == 0 && (spec.snap_until == 0 || steps <= spec.snap_until) {
            // unified memory: read H in place
            let hx = unsafe { std::slice::from_raw_parts(h[0].contents() as *const f32, n) };
            let hy = unsafe { std::slice::from_raw_parts(h[1].contents() as *const f32, n) };
            snaps.push(crate::fdtd::Snap { step: steps, j: spec.snap_k.iter().map(|&k| crate::fdtd::sheet_current(nx, ny, nz, hx, hy, k)).collect() });
        }
        if check {
            let en: f64 = read(&eout, ny * nz).iter().map(|&x| x as f64).sum();
            peak = peak.max(en);
            if steps > spec.exc.len() && en <= peak * 10f64.powf(spec.end_db / 10.0) {
                done = true;
            }
            if !en.is_finite() {
                return Err(format!("fields diverged at step {steps}"));
            }
        }
    }
    let fields = (!spec.dft_freqs.is_empty()).then(|| FieldDft {
        freqs: spec.dft_freqs.clone(),
        data: dftb.iter().map(|bs| std::array::from_fn(|c| read(&bs[c], 2 * n))).collect(),
    });
    Ok(RunResult { steps, v, probes: pv, fields, snaps, seconds: t0.elapsed().as_secs_f64(), backend: format!("metal ({})", dev.name()) })
}
