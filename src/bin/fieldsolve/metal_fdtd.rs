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
    ng: u32,
    pl: [u32; 3],
    kd: [u32; 3],
    w: [f32; 12],
}

struct Pipes {
    h: ComputePipelineState,
    e: ComputePipelineState,
    ei: ComputePipelineState,
    hf: ComputePipelineState,
    ef: ComputePipelineState,
    ph: ComputePipelineState,
    pe: ComputePipelineState,
    src: ComputePipelineState,
    probe: ComputePipelineState,
    dft: ComputePipelineState,
    gather: ComputePipelineState,
    energy: ComputePipelineState,
    cond: ComputePipelineState,
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
    let pp = Pipes { h: pipe("update_h")?, e: pipe("update_e")?, ei: pipe("update_e_idx")?, hf: pipe("update_h_f")?, ef: pipe("update_e_f")?, ph: pipe("pml_h")?, pe: pipe("pml_e")?, src: pipe("add_src")?, probe: pipe("probe")?, dft: pipe("dft")?, gather: pipe("dft_gather")?, energy: pipe("energy")?, cond: pipe("cond_e")? };
    let queue = dev.new_command_queue();
    let (nx, ny, nz) = (p.nx, p.ny, p.nz);
    let n = p.n();
    let e: Vec<Buffer> = (0..3).map(|_| zeros(&dev, n)).collect();
    let h: Vec<Buffer> = (0..3).map(|_| zeros(&dev, n)).collect();
    let ca: Vec<Buffer> = p.ca.iter().map(|v| buf_f32(&dev, v)).collect();
    let cb: Vec<Buffer> = p.cb.iter().map(|v| buf_f32(&dev, v)).collect();
    // material indices: the distinct (ca, cb) pairs over all edges, if they fit in u16
    let mat = {
        let mut tab: Vec<[f32; 2]> = Vec::new();
        let mut seen: std::collections::HashMap<(u32, u32), u16> = std::collections::HashMap::new();
        let mut ids: Vec<Vec<u16>> = vec![Vec::with_capacity(n); 3];
        let mut ok = true;
        'outer: for c in 0..3 {
            for q in 0..n {
                let key = (p.ca[c][q].to_bits(), p.cb[c][q].to_bits());
                let m = match seen.get(&key) {
                    Some(&m) => m,
                    None => {
                        if tab.len() >= u16::MAX as usize {
                            ok = false;
                            break 'outer;
                        }
                        tab.push([p.ca[c][q], p.cb[c][q]]);
                        let m = (tab.len() - 1) as u16;
                        seen.insert(key, m);
                        m
                    }
                };
                ids[c].push(m);
            }
        }
        ok.then(|| {
            let flat: Vec<f32> = tab.iter().flatten().copied().collect();
            let mk = |v: &[u16]| dev.new_buffer_with_data(v.as_ptr() as *const c_void, (v.len().max(1) * 2) as u64, MTLResourceOptions::StorageModeShared);
            (ids.iter().map(|v| mk(v)).collect::<Vec<Buffer>>(), buf_f32(&dev, &flat))
        })
    };
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
        map: Buffer,
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
                map: buf_u32(&dev, &ax.map),
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
    let gidx: Vec<Buffer> = (0..6).map(|c| buf_u32(&dev, &spec.gather[c])).collect();
    let gacc: Vec<Vec<Buffer>> = spec.gather_freqs.iter().map(|_| (0..6).map(|c| zeros(&dev, 2 * spec.gather[c].len())).collect()).collect();

    // Debye currents (compact over the dispersive planes): [component][q][pole]
    let (djb, dbb, khb, kd) = match &p.disp {
        Some(d) => {
            let nd = d.nd(nx, ny);
            let b: Vec<f32> = (0..3).flat_map(|c| d.b[c].iter().copied()).collect();
            let kh: Vec<f32> = (0..d.np).flat_map(|q| [d.kp[q], d.hk[q]]).collect();
            (zeros(&dev, 3 * nd * d.np), buf_f32(&dev, &b), buf_f32(&dev, &kh), [d.k0 as u32, d.k1 as u32, d.np as u32])
        }
        None => (zeros(&dev, 1), zeros(&dev, 1), zeros(&dev, 2), [0, 0, 0]),
    };
    // surface-impedance copper edges
    let ncond = p.cond.as_ref().map(|c| c.ids.len()).unwrap_or(0);
    let (cidb, cfacb, clayb, ccib, ctabb, cdecb, cnb) = match &p.cond {
        Some(c) if ncond > 0 => (buf_u32(&dev, &c.ids), buf_f32(&dev, &c.fac), buf_u32(&dev, &c.layer), zeros(&dev, ncond * c.nb), buf_f32(&dev, &c.tab), buf_f32(&dev, &c.dec), c.nb),
        _ => (zeros(&dev, 1), zeros(&dev, 1), zeros(&dev, 1), zeros(&dev, 1), zeros(&dev, 1), zeros(&dev, 1), 0),
    };
    let prm = P { nx: nx as u32, ny: ny as u32, nz: nz as u32, kmu: (p.dt / MU0) as f32, nsrc: nsrc as u32, nprobe: nprobe as u32, pl: [pml[0].len as u32, pml[1].len as u32, pml[2].len as u32], kd, ..Default::default() };
    // fused H/E kernels (CPML inline, material indices) unless disabled for comparison
    let fused = mat.is_some() && std::env::var("FIELDSOLVE_FUSED").map(|v| v != "0").unwrap_or(true);
    if p.disp.is_some() && !fused {
        return Err("wideband (Debye) dielectrics need the fused Metal kernel (material table too large, or FIELDSOLVE_FUSED=0)".into());
    }
    let psz = std::mem::size_of::<P>() as u64;
    // threadgroup shape for the 3-D kernels: FIELDSOLVE_TG="32x4x1" overrides (for benchmarking)
    let tg3 = std::env::var("FIELDSOLVE_TG").ok().and_then(|s| {
        let v: Vec<u64> = s.split('x').filter_map(|x| x.parse().ok()).collect();
        (v.len() == 3).then(|| MTLSize { width: v[0], height: v[1], depth: v[2] })
    }).unwrap_or(MTLSize { width: 32, height: 4, depth: 1 });
    let grid3 = MTLSize { width: nx as u64, height: ny as u64, depth: nz as u64 };
    let mut v: Vec<Vec<f32>> = vec![vec![0.0]; p.ports.len()];
    let mut pv: Vec<Vec<f32>> = vec![vec![0.0]; p.probes.len()];
    let mut peak = 0.0f64;
    let mut port_peak = 0.0f64;
    let mut ghist: Vec<Vec<f64>> = Vec::new();
    let mut ehist: Vec<(usize, f64)> = Vec::new();
    let mut steps = 0usize;
    let mut done = false;
    let mut snaps = Vec::new();
    let setup_seconds = t0.elapsed().as_secs_f64();
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
            if fused {
                let (mi, tab) = mat.as_ref().unwrap();
                enc.set_compute_pipeline_state(&pp.hf);
                for c in 0..3 {
                    enc.set_buffer(c as u64, Some(&h[c]), 0);
                    enc.set_buffer(3 + c as u64, Some(&e[c]), 0);
                }
                enc.set_buffer(6, Some(&inv_d), 0);
                set_p(&prm);
                for a in 0..3 {
                    let pm = &pml[a];
                    let base = 8 + 5 * a as u64;
                    enc.set_buffer(base, Some(&pm.map), 0);
                    enc.set_buffer(base + 1, Some(&pm.bh), 0);
                    enc.set_buffer(base + 2, Some(&pm.ch), 0);
                    enc.set_buffer(base + 3, Some(&pm.psi[2]), 0);
                    enc.set_buffer(base + 4, Some(&pm.psi[3]), 0);
                }
                enc.dispatch_threads(grid3, tg3);
                enc.set_compute_pipeline_state(&pp.ef);
                for c in 0..3 {
                    enc.set_buffer(c as u64, Some(&e[c]), 0);
                    enc.set_buffer(3 + c as u64, Some(&h[c]), 0);
                    enc.set_buffer(8 + c as u64, Some(&mi[c]), 0);
                }
                enc.set_buffer(6, Some(&inv_dd), 0);
                enc.set_buffer(11, Some(tab), 0);
                set_p(&prm);
                for a in 0..3 {
                    let pm = &pml[a];
                    let base = 12 + 5 * a as u64;
                    enc.set_buffer(base, Some(&pm.map), 0);
                    enc.set_buffer(base + 1, Some(&pm.be), 0);
                    enc.set_buffer(base + 2, Some(&pm.ce), 0);
                    enc.set_buffer(base + 3, Some(&pm.psi[0]), 0);
                    enc.set_buffer(base + 4, Some(&pm.psi[1]), 0);
                }
                enc.set_buffer(27, Some(&djb), 0);
                enc.set_buffer(28, Some(&dbb), 0);
                enc.set_buffer(29, Some(&khb), 0);
                enc.dispatch_threads(grid3, tg3);
            } else {
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
            if let Some((mi, tab)) = &mat {
                enc.set_compute_pipeline_state(&pp.ei);
                for c in 0..3 {
                    enc.set_buffer(c as u64, Some(&e[c]), 0);
                    enc.set_buffer(3 + c as u64, Some(&h[c]), 0);
                    enc.set_buffer(8 + c as u64, Some(&mi[c]), 0);
                }
                enc.set_buffer(11, Some(tab), 0);
            } else {
                enc.set_compute_pipeline_state(&pp.e);
                for c in 0..3 {
                    enc.set_buffer(c as u64, Some(&e[c]), 0);
                    enc.set_buffer(3 + c as u64, Some(&h[c]), 0);
                    enc.set_buffer(8 + c as u64, Some(&ca[c]), 0);
                    enc.set_buffer(11 + c as u64, Some(&cb[c]), 0);
                }
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
                enc.set_buffer(8, Some(&djb), 0);
                enc.set_buffer(9, Some(&dbb), 0);
                let mut q = prm;
                q.vs = spec.exc[n_step];
                set_p(&q);
                enc.dispatch_threads(MTLSize { width: nsrc as u64, height: 1, depth: 1 }, MTLSize { width: 64.min(nsrc as u64).max(1), height: 1, depth: 1 });
            }
            // conductor loss: E on copper edges from the sheet current (after the source, before the probes)
            if ncond > 0 {
                enc.set_compute_pipeline_state(&pp.cond);
                enc.set_buffer(0, Some(&e[0]), 0);
                enc.set_buffer(1, Some(&e[1]), 0);
                for c in 0..3 {
                    enc.set_buffer(3 + c as u64, Some(&h[c]), 0);
                }
                enc.set_buffer(6, Some(&inv_dd), 0);
                enc.set_buffer(8, Some(&cidb), 0);
                enc.set_buffer(9, Some(&cfacb), 0);
                enc.set_buffer(10, Some(&clayb), 0);
                enc.set_buffer(11, Some(&ccib), 0);
                enc.set_buffer(12, Some(&ctabb), 0);
                enc.set_buffer(13, Some(&cdecb), 0);
                let mut q = prm;
                q.axis = cnb as u32;
                q.slab_len = ncond as u32;
                set_p(&q);
                enc.dispatch_threads(MTLSize { width: ncond as u64, height: 1, depth: 1 }, MTLSize { width: 256, height: 1, depth: 1 });
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
            // gathered DFT (near-to-far-field box faces)
            if !spec.gather_freqs.is_empty() && s1 % spec.dft_every == 0 {
                enc.set_compute_pipeline_state(&pp.gather);
                for (fi, &f) in spec.gather_freqs.iter().enumerate() {
                    let w = 2.0 * std::f64::consts::PI * f;
                    let mut q = prm;
                    for c in 0..6 {
                        let (wr, wi) = crate::fdtd::dft_weight(p, c, s1, spec.dft_every, w);
                        q.w[2 * c] = wr;
                        q.w[2 * c + 1] = wi;
                    }
                    for c in 0..6 {
                        let cnt = spec.gather[c].len();
                        if cnt == 0 {
                            continue;
                        }
                        q.axis = c as u32;
                        q.ng = cnt as u32;
                        enc.set_buffer(0, Some(if c < 3 { &e[c] } else { &h[c - 3] }), 0);
                        enc.set_buffer(1, Some(&gacc[fi][c]), 0);
                        enc.set_buffer(2, Some(&gidx[c]), 0);
                        set_p(&q);
                        enc.dispatch_threads(MTLSize { width: cnt as u64, height: 1, depth: 1 }, MTLSize { width: 256, height: 1, depth: 1 });
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
            ehist.push((steps, en));
            peak = peak.max(en);
            if !en.is_finite() {
                return Err(format!("fields diverged at step {steps}"));
            }
            // as on the CPU: a configured gathered-field metric gates every convergence exit
            let gs = match &spec.gather_metric {
                Some(m) => {
                    // unified memory: read the accumulators in place
                    let cur: Vec<[Vec<f32>; 6]> = gacc.iter().map(|bs| std::array::from_fn(|c| read(&bs[c], 2 * spec.gather[c].len()))).collect();
                    ghist.push((m.0)(&cur));
                    crate::fdtd::gather_settled(&ghist, spec.gather_lag, spec.gather_tol)
                }
                None => true,
            };
            if gs && steps > spec.exc.len() && en <= peak * 10f64.powf(spec.end_db / 10.0) {
                done = true;
            }
            if spec.port_stop {
                let pd = crate::fdtd::ports_decayed(&v, spec.check_every, spec.end_db, &mut port_peak) && steps > spec.exc.len();
                if pd && gs {
                    done = true;
                }
            }
        }
    }
    let fields = (!spec.dft_freqs.is_empty()).then(|| FieldDft {
        freqs: spec.dft_freqs.clone(),
        data: dftb.iter().map(|bs| std::array::from_fn(|c| read(&bs[c], 2 * n))).collect(),
    });
    let gathered = gacc.iter().map(|bs| std::array::from_fn(|c| read(&bs[c], 2 * spec.gather[c].len()))).collect();
    Ok(RunResult { steps, v, probes: pv, fields, gathered, energy: ehist, snaps, seconds: t0.elapsed().as_secs_f64(), setup_seconds, backend: format!("metal ({})", dev.name()) })
}
