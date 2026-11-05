//! Field lines: streamlines of the instantaneous field Re[F e^{j phase}] from a snapshot's complex node
//! vectors (`fields/E_vec_<f>.f32`, `H_vec_<f>.f32`).
//!
//! Seeds are drawn with probability proportional to |F| (deterministic), outside the PML. Each line is
//! traced both ways with RK4 on the unit direction field (trilinear interpolation on the non-uniform node
//! grid), with a step of half the local cell. A line stops where the field falls below `floor_db` of the
//! maximum, when it leaves the interior, closes on itself, or runs into a voxel already crowded by other
//! lines (an occupancy grid keeps the density even, after Jobard and Lefer).
//!
//! Output `fields/lines_<F>_<f>.json`: `{"field", "freq", "phases": [{"phase", "lines": [{"p": [x,y,z,...],
//! "m": [|F|...]}]}], "max"}` in mm, Gerber orientation.

use serde_json::{json, Value};
use std::path::Path;

struct Field {
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
    v: Vec<[f32; 3]>,
    nx: usize,
    ny: usize,
    nz: usize,
}

impl Field {
    fn cell(lines: &[f64], u: f64) -> Option<(usize, f64)> {
        if u < lines[0] || u > lines[lines.len() - 1] {
            return None;
        }
        let i = (lines.partition_point(|&p| p <= u)).clamp(1, lines.len() - 1) - 1;
        Some((i, (u - lines[i]) / (lines[i + 1] - lines[i])))
    }
    fn at(&self, p: [f64; 3]) -> Option<[f64; 3]> {
        let (i, fx) = Self::cell(&self.x, p[0])?;
        let (j, fy) = Self::cell(&self.y, p[1])?;
        let (k, fz) = Self::cell(&self.z, p[2])?;
        let mut o = [0.0f64; 3];
        for (dk, wz) in [(0, 1.0 - fz), (1, fz)] {
            for (dj, wy) in [(0, 1.0 - fy), (1, fy)] {
                for (di, wx) in [(0, 1.0 - fx), (1, fx)] {
                    let w = wx * wy * wz;
                    if w == 0.0 {
                        continue;
                    }
                    let id = ((k + dk) * self.ny + j + dj) * self.nx + i + di;
                    let v = self.v[id];
                    for c in 0..3 {
                        o[c] += w * v[c] as f64;
                    }
                }
            }
        }
        Some(o)
    }
    fn local_h(&self, p: [f64; 3]) -> f64 {
        let mut h = f64::INFINITY;
        for (l, u) in [(&self.x, p[0]), (&self.y, p[1]), (&self.z, p[2])] {
            if let Some((i, _)) = Self::cell(l, u) {
                h = h.min(l[i + 1] - l[i]);
            }
        }
        h
    }
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

pub struct Opts {
    pub field: String,
    pub tag: String,
    pub seeds: usize,
    pub phases: usize,
    pub floor_db: f64,
    pub max_len_mm: f64,
    pub crowd: u16,
}

pub fn trace(dir: &Path, o: &Opts) -> Result<std::path::PathBuf, String> {
    let meta: Value = serde_json::from_str(&std::fs::read_to_string(dir.join("meta.json")).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
    let sh: Vec<usize> = meta["shape"].as_array().ok_or("meta.json: shape")?.iter().map(|v| v.as_u64().unwrap() as usize).collect();
    let (nx, ny, nz) = (sh[0], sh[1], sh[2]);
    let ax = |k: &str| -> Vec<f64> { meta[k].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect() };
    let (x, y, z) = (ax("x"), ax("y"), ax("z"));
    let pml: Vec<usize> = meta.get("pml").and_then(|v| v.as_array()).map(|a| a.iter().map(|v| v.as_u64().unwrap_or(0) as usize).collect()).unwrap_or(vec![0; 6]);
    let raw = std::fs::read(dir.join(format!("{}_vec_{}.f32", o.field, o.tag))).map_err(|e| format!("{}_vec_{}.f32: {e} (re-run fieldsolve to export vector fields)", o.field, o.tag))?;
    let fl: Vec<f32> = raw.chunks_exact(4).map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]])).collect();
    let n = nx * ny * nz;
    if fl.len() != 6 * n {
        return Err("vector file does not match meta.json shape".into());
    }
    // interior box (outside the PML)
    let lo = [x[pml[0].min(nx - 1)], y[pml[2].min(ny - 1)], z[pml[4].min(nz - 1)]];
    let hi = [x[nx - 1 - pml[1].min(nx - 1)], y[ny - 1 - pml[3].min(ny - 1)], z[nz - 1 - pml[5].min(nz - 1)]];
    let inside = |p: [f64; 3]| (0..3).all(|c| p[c] >= lo[c] && p[c] <= hi[c]);
    let hmin = [&x, &y, &z].iter().map(|l| l.windows(2).map(|w| w[1] - w[0]).fold(f64::INFINITY, f64::min)).fold(f64::INFINITY, f64::min);
    // occupancy voxels: 4x the finest cell
    let vox = 4.0 * hmin;
    let (ox, oy, oz) = (((hi[0] - lo[0]) / vox).ceil() as usize + 1, ((hi[1] - lo[1]) / vox).ceil() as usize + 1, ((hi[2] - lo[2]) / vox).ceil() as usize + 1);
    let vid = |p: [f64; 3]| -> usize {
        let i = (((p[0] - lo[0]) / vox) as usize).min(ox - 1);
        let j = (((p[1] - lo[1]) / vox) as usize).min(oy - 1);
        let k = (((p[2] - lo[2]) / vox) as usize).min(oz - 1);
        (k * oy + j) * ox + i
    };
    let mut phases_out = Vec::new();
    let mut gmax = 0.0f64;
    for ph in 0..o.phases.max(1) {
        let phase = ph as f64 / o.phases.max(1) as f64 * std::f64::consts::TAU;
        let (c, s) = (phase.cos() as f32, phase.sin() as f32);
        let v: Vec<[f32; 3]> = (0..n).map(|id| {
            let q = &fl[6 * id..6 * id + 6];
            [q[0] * c - q[3] * s, q[1] * c - q[4] * s, q[2] * c - q[5] * s]
        }).collect();
        let f = Field { x: x.clone(), y: y.clone(), z: z.clone(), v, nx, ny, nz };
        // magnitude at interior nodes; seed CDF
        let mut mags = Vec::with_capacity(n);
        let mut fmax = 0.0f64;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let p = [x[i], y[j], z[k]];
                    let m = if inside(p) { norm(f.v[(k * ny + j) * nx + i].map(|a| a as f64)) } else { 0.0 };
                    fmax = fmax.max(m);
                    mags.push(m);
                }
            }
        }
        gmax = gmax.max(fmax);
        let floor = fmax * 10f64.powf(o.floor_db / 20.0);
        let mut cdf = Vec::with_capacity(n);
        let mut acc = 0.0;
        for &m in &mags {
            acc += if m > floor { m } else { 0.0 };
            cdf.push(acc);
        }
        let mut occ = vec![0u16; ox * oy * oz];
        let mut rng = 0x9E3779B97F4A7C15u64 ^ (ph as u64 + 1);
        let mut next = || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            (rng >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut lines = Vec::new();
        let mut tries = 0;
        while lines.len() < o.seeds && tries < o.seeds * 20 && acc > 0.0 {
            tries += 1;
            let r = next() * acc;
            let id = cdf.partition_point(|&c| c < r).min(n - 1);
            let (i, j, k) = (id % nx, (id / nx) % ny, id / (nx * ny));
            // jitter inside the node's neighbourhood
            let jit = |l: &Vec<f64>, a: usize, t: f64| -> f64 {
                let h = if a + 1 < l.len() { l[a + 1] - l[a] } else { l[a] - l[a - 1] };
                (l[a] + (t - 0.5) * h).clamp(l[0], l[l.len() - 1])
            };
            let seed = [jit(&x, i, next()), jit(&y, j, next()), jit(&z, k, next())];
            if !inside(seed) || occ[vid(seed)] >= o.crowd {
                continue;
            }
            let mut halves: Vec<Vec<([f64; 3], f64)>> = Vec::new();
            let mut closed = false;
            for sgn in [1.0f64, -1.0] {
                if closed {
                    halves.push(Vec::new());
                    continue;
                }
                let mut pts = Vec::new();
                let mut p = seed;
                let mut len = 0.0;
                let mut own: Vec<usize> = Vec::new();
                for _ in 0..20000 {
                    let Some(fv) = f.at(p) else { break };
                    let m = norm(fv);
                    if m < floor || !inside(p) {
                        break;
                    }
                    pts.push((p, m));
                    let vv = vid(p);
                    if !own.contains(&vv) {
                        if occ[vv] >= o.crowd && len > 2.0 * vox {
                            break;
                        }
                        own.push(vv);
                    }
                    let h = 0.5 * f.local_h(p);
                    let dirf = |q: [f64; 3]| -> Option<[f64; 3]> {
                        let v = f.at(q)?;
                        let m = norm(v);
                        if m == 0.0 {
                            return None;
                        }
                        Some([sgn * v[0] / m, sgn * v[1] / m, sgn * v[2] / m])
                    };
                    let add = |a: [f64; 3], b: [f64; 3], t: f64| [a[0] + t * b[0], a[1] + t * b[1], a[2] + t * b[2]];
                    let Some(k1) = dirf(p) else { break };
                    let Some(k2) = dirf(add(p, k1, h / 2.0)) else { break };
                    let Some(k3) = dirf(add(p, k2, h / 2.0)) else { break };
                    let Some(k4) = dirf(add(p, k3, h)) else { break };
                    let d = [(k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0, (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0, (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]) / 6.0];
                    p = add(p, d, h);
                    len += h;
                    if len > o.max_len_mm {
                        break;
                    }
                    // closed loop: back near the seed after going some way
                    if len > 6.0 * h && norm([p[0] - seed[0], p[1] - seed[1], p[2] - seed[2]]) < 0.6 * h {
                        pts.push((seed, m));
                        closed = true;
                        break;
                    }
                }
                halves.push(pts);
            }
            let mut pts: Vec<([f64; 3], f64)> = halves[1].iter().rev().cloned().collect();
            if !pts.is_empty() {
                pts.pop(); // the seed appears in both halves
            }
            pts.extend(halves[0].iter().cloned());
            if pts.len() < 4 {
                continue;
            }
            for (q, _) in &pts {
                let vv = vid(*q);
                occ[vv] = occ[vv].saturating_add(1);
            }
            let flat: Vec<f64> = pts.iter().flat_map(|(q, _)| q.iter().map(|v| (v * 1e4).round() / 1e4)).collect();
            let mm: Vec<f64> = pts.iter().map(|(_, m)| (*m as f32) as f64).collect();
            lines.push(json!({"p": flat, "m": mm}));
        }
        eprintln!("fieldsolve lines: phase {:.0} deg: {} lines", phase.to_degrees(), lines.len());
        phases_out.push(json!({"phase": phase, "lines": lines}));
    }
    let out = dir.join(format!("lines_{}_{}.json", o.field, o.tag));
    std::fs::write(&out, serde_json::to_string(&json!({"field": o.field, "freq": o.tag, "max": gmax, "floor_db": o.floor_db, "phases": phases_out})).unwrap()).map_err(|e| e.to_string())?;
    Ok(out)
}
