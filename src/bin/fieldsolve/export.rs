//! Field and geometry export.
//!
//! For each snapshot frequency `f` (in `out/fields/`):
//! - `fields_<f>.vtr`: VTK rectilinear grid (ParaView), node-interpolated complex E and H as magnitude
//!   arrays `E_abs`, `H_abs` plus the real parts of each component at phase 0 (`E_re`, `H_re` vectors).
//! - `E_abs_<f>.f32`, `H_abs_<f>.f32`: raw little-endian f32 volumes, x fastest (for Blender).
//! - `J_<layer>_<f>.f32`: surface current density magnitude |n x (H+ - H-)| on each copper layer (A/m),
//!   a 2-D (nx, ny) raw array, x fastest.
//! - `E_vec_<f>.f32`, `H_vec_<f>.f32`: complex node vectors, 6 f32 per node (re x, y, z, im x, y, z), for
//!   `fieldsolve lines` (field-line tracing).
//! - `meta.json`: grid lines (mm), layer heights, array names.
//! - `copper_<layer>.obj`: the layer's copper as n-gon faces at its height (mm), for context in Blender.

use crate::fdtd::FieldDft;
use crate::geom::{Shape, P2};
use crate::model::Problem;
use crate::sim::Setup;
use serde_json::json;
use std::io::Write;
use std::path::Path;

fn write_f32(path: &Path, v: &[f32]) -> std::io::Result<()> {
    let mut b = Vec::with_capacity(v.len() * 4);
    for x in v {
        b.extend_from_slice(&x.to_le_bytes());
    }
    std::fs::write(path, b)
}

/// Complex component `c` (0..6) at node (i,j,k), averaged from its staggered samples: E_a from its two
/// samples along a, H_a from its four samples in the plane across a (offsets fixed, no allocation).
fn node_val(p: &Problem, fd: &[Vec<f32>; 6], c: usize, i: usize, j: usize, k: usize) -> (f32, f32) {
    let a = c % 3;
    let mut offs = [[0isize; 3]; 4];
    let len = if c >= 3 {
        let (b, d) = ((a + 1) % 3, (a + 2) % 3);
        offs[1][b] = -1;
        offs[2][d] = -1;
        offs[3][b] = -1;
        offs[3][d] = -1;
        4
    } else {
        offs[1][a] = -1;
        2
    };
    let (mut re, mut im, mut n) = (0.0f32, 0.0f32, 0.0f32);
    for o in &offs[..len] {
        let (ii, jj, kk) = (i as isize + o[0], j as isize + o[1], k as isize + o[2]);
        if ii < 0 || jj < 0 || kk < 0 || ii >= p.nx as isize || jj >= p.ny as isize || kk >= p.nz as isize {
            continue;
        }
        let id = ((kk as usize * p.ny + jj as usize) * p.nx + ii as usize) * 2;
        re += fd[c][id];
        im += fd[c][id + 1];
        n += 1.0;
    }
    if n > 0.0 {
        (re / n, im / n)
    } else {
        (0.0, 0.0)
    }
}

/// The previous allocating version, kept to check the fixed-offset one against.
#[cfg(test)]
fn node_val_reference(p: &Problem, fd: &[Vec<f32>; 6], c: usize, i: usize, j: usize, k: usize) -> (f32, f32) {
    let (nx, ny, nz) = (p.nx, p.ny, p.nz);
    let a = c % 3;
    let mut offs: Vec<[isize; 3]> = vec![[0, 0, 0]];
    let stag: Vec<usize> = if c >= 3 { vec![(a + 1) % 3, (a + 2) % 3] } else { vec![a] };
    for s in stag {
        let more: Vec<[isize; 3]> = offs.iter().map(|o| { let mut q = *o; q[s] = -1; q }).collect();
        offs.extend(more);
    }
    let (mut re, mut im, mut n) = (0.0f32, 0.0f32, 0.0f32);
    for o in offs {
        let (ii, jj, kk) = (i as isize + o[0], j as isize + o[1], k as isize + o[2]);
        if ii < 0 || jj < 0 || kk < 0 || ii >= nx as isize || jj >= ny as isize || kk >= nz as isize {
            continue;
        }
        let id = ((kk as usize * ny + jj as usize) * nx + ii as usize) * 2;
        re += fd[c][id];
        im += fd[c][id + 1];
        n += 1.0;
    }
    if n > 0.0 { (re / n, im / n) } else { (0.0, 0.0) }
}

#[cfg(test)]
mod node_val_tests {
    use super::*;
    #[test]
    fn fixed_offsets_match_the_allocating_version_exactly() {
        let g = crate::mesh::Grid { x: (0..7).map(|v| v as f64 * 1e-3).collect(), y: (0..6).map(|v| v as f64 * 1e-3).collect(), z: (0..5).map(|v| v as f64 * 1e-3).collect() };
        let p = crate::model::build(&g, &crate::model::Layout::default(), &[], &crate::model::BuildOpts { pml: [0; 6], ..Default::default() });
        let n = p.nx * p.ny * p.nz * 2;
        let mut seed = 12345u64;
        let mut rnd = || { seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); ((seed >> 33) as f32) / (1u64 << 31) as f32 - 0.5 };
        let fd: [Vec<f32>; 6] = std::array::from_fn(|_| (0..n).map(|_| rnd()).collect());
        for c in 0..6 {
            for k in 0..p.nz {
                for j in 0..p.ny {
                    for i in 0..p.nx {
                        assert_eq!(node_val(&p, &fd, c, i, j, k), node_val_reference(&p, &fd, c, i, j, k), "c {c} at {i},{j},{k}");
                    }
                }
            }
        }
    }
}

fn vtr(path: &Path, p: &Problem, arrays: &[(&str, usize, &[f32])]) -> std::io::Result<()> {
    let (nx, ny, nz) = (p.nx, p.ny, p.nz);
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    let ext = format!("0 {} 0 {} 0 {}", nx - 1, ny - 1, nz - 1);
    writeln!(f, "<?xml version=\"1.0\"?>\n<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">")?;
    writeln!(f, "<RectilinearGrid WholeExtent=\"{ext}\"><Piece Extent=\"{ext}\">")?;
    let mut off = 0u64;
    let mut blobs: Vec<Vec<u8>> = Vec::new();
    let mut push = |v: Vec<u8>, blobs: &mut Vec<Vec<u8>>| -> u64 {
        let o = off;
        off += 8 + v.len() as u64;
        blobs.push(v);
        o
    };
    writeln!(f, "<PointData>")?;
    for (name, ncomp, data) in arrays {
        let mut b = Vec::with_capacity(data.len() * 4);
        data.iter().for_each(|x| b.extend_from_slice(&x.to_le_bytes()));
        let o = push(b, &mut blobs);
        writeln!(f, "<DataArray type=\"Float32\" Name=\"{name}\" NumberOfComponents=\"{ncomp}\" format=\"appended\" offset=\"{o}\"/>")?;
    }
    writeln!(f, "</PointData><Coordinates>")?;
    for a in 0..3 {
        let mut b = Vec::new();
        p.grid.axis(a).iter().for_each(|x| b.extend_from_slice(&((x * 1e3) as f32).to_le_bytes()));
        let o = push(b, &mut blobs);
        writeln!(f, "<DataArray type=\"Float32\" Name=\"{}\" format=\"appended\" offset=\"{o}\"/>", ["x", "y", "z"][a])?;
    }
    writeln!(f, "</Coordinates></Piece></RectilinearGrid>\n<AppendedData encoding=\"raw\">")?;
    f.write_all(b"_")?;
    for b in blobs {
        f.write_all(&(b.len() as u64).to_le_bytes())?;
        f.write_all(&b)?;
    }
    writeln!(f, "\n</AppendedData></VTKFile>")?;
    Ok(())
}

pub fn write_copper_obj(path: &Path, prims: &[crate::geom::Prim], z: f64) -> std::io::Result<()> {
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    let mut nv = 0usize;
    for pr in prims.iter().filter(|p| p.dark) {
        let pts: Vec<P2> = match &pr.shape.outline_shape() {
            Shape::Macro(_) => Vec::new(), // outline_shape never returns a macro
            Shape::Polygon(v) => v.clone(),
            Shape::Circle { c, r } => (0..24).map(|k| {
                let t = k as f64 / 24.0 * std::f64::consts::TAU;
                P2::new(c.x + r * t.cos(), c.y + r * t.sin())
            }).collect(),
            Shape::Capsule { a, b, r } => {
                let ang = (b.y - a.y).atan2(b.x - a.x);
                let mut v = Vec::new();
                for (c, base) in [(b, ang - std::f64::consts::FRAC_PI_2), (a, ang + std::f64::consts::FRAC_PI_2)] {
                    for k in 0..=8 {
                        let t = base + k as f64 / 8.0 * std::f64::consts::PI;
                        v.push(P2::new(c.x + r * t.cos(), c.y + r * t.sin()));
                    }
                }
                v
            }
        };
        if pts.len() < 3 {
            continue;
        }
        for q in &pts {
            writeln!(f, "v {:.5} {:.5} {:.5}", q.x, q.y, z)?;
        }
        write!(f, "f")?;
        for k in 0..pts.len() {
            write!(f, " {}", nv + k + 1)?;
        }
        writeln!(f)?;
        nv += pts.len();
    }
    Ok(())
}

pub fn write_fields(out: &Path, p: &Problem, st: &Setup, fd: &FieldDft) -> std::io::Result<()> {
    let dir = out.join("fields");
    std::fs::create_dir_all(&dir)?;
    let (nx, ny, nz) = (p.nx, p.ny, p.nz);
    let n = p.n();
    let zmm: Vec<f64> = p.grid.z.iter().map(|v| v * 1e3).collect();
    let mut layers = Vec::new();
    for s in &st.layout.sheets {
        let k = crate::mesh::snap(&zmm, s.z);
        let tag = s.name.replace('.', "_");
        write_copper_obj(&dir.join(format!("copper_{tag}.obj")), &s.prims, s.z)?;
        layers.push((s.name.clone(), tag, k, s.z));
    }
    let mut snaps = Vec::new();
    for (fi, &f) in fd.freqs.iter().enumerate() {
        let d = &fd.data[fi];
        let ftag = format!("{:.4}GHz", f * 1e-9);
        let mut eabs = vec![0.0f32; n];
        let mut habs = vec![0.0f32; n];
        let mut ere = vec![0.0f32; 3 * n];
        let mut hre = vec![0.0f32; 3 * n];
        // complex node vectors (re x y z, im x y z) for field-line tracing
        let mut evec = vec![0.0f32; 6 * n];
        let mut hvec = vec![0.0f32; 6 * n];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let id = (k * ny + j) * nx + i;
                    let (mut se, mut sh) = (0.0f32, 0.0f32);
                    for c in 0..3 {
                        let (r, im) = node_val(p, d, c, i, j, k);
                        se += r * r + im * im;
                        ere[3 * id + c] = r;
                        evec[6 * id + c] = r;
                        evec[6 * id + 3 + c] = im;
                        let (r, im) = node_val(p, d, c + 3, i, j, k);
                        sh += r * r + im * im;
                        hre[3 * id + c] = r;
                        hvec[6 * id + c] = r;
                        hvec[6 * id + 3 + c] = im;
                    }
                    eabs[id] = se.sqrt();
                    habs[id] = sh.sqrt();
                }
            }
        }
        vtr(&dir.join(format!("fields_{ftag}.vtr")), p, &[("E_abs", 1, &eabs), ("H_abs", 1, &habs), ("E_re", 3, &ere), ("H_re", 3, &hre)])?;
        write_f32(&dir.join(format!("E_abs_{ftag}.f32")), &eabs)?;
        write_f32(&dir.join(format!("H_abs_{ftag}.f32")), &habs)?;
        write_f32(&dir.join(format!("E_vec_{ftag}.f32")), &evec)?;
        write_f32(&dir.join(format!("H_vec_{ftag}.f32")), &hvec)?;
        // surface current on each sheet: Jx = Hy(k) - Hy(k-1), Jy = -(Hx(k) - Hx(k-1)) at the node
        for (_, tag, k, _) in &layers {
            let k = *k;
            if k == 0 || k + 1 >= nz {
                continue;
            }
            let mut js = vec![0.0f32; nx * ny];
            for j in 0..ny {
                for i in 0..nx {
                    let g = |c: usize, kk: usize, ii: usize, jj: usize| -> (f32, f32) {
                        let id = ((kk * ny + jj) * nx + ii) * 2;
                        (d[c][id], d[c][id + 1])
                    };
                    // Hy at (x_{i+1/2}, y_j, z_{k+-1/2}); Hx at (x_i, y_{j+1/2}, z_{k+-1/2})
                    let i1 = i.min(nx - 2);
                    let j1 = j.min(ny - 2);
                    let (hy_u, hy_ui) = g(4, k, i1, j);
                    let (hy_d, hy_di) = g(4, k - 1, i1, j);
                    let (hx_u, hx_ui) = g(3, k, i, j1);
                    let (hx_d, hx_di) = g(3, k - 1, i, j1);
                    let (jxr, jxi) = (hy_u - hy_d, hy_ui - hy_di);
                    let (jyr, jyi) = (hx_d - hx_u, hx_di - hx_ui);
                    js[j * nx + i] = (jxr * jxr + jxi * jxi + jyr * jyr + jyi * jyi).sqrt();
                }
            }
            write_f32(&dir.join(format!("J_{tag}_{ftag}.f32")), &js)?;
        }
        snaps.push(json!({"f": f, "tag": ftag}));
    }
    let meta = json!({
        "shape": [nx, ny, nz],
        "x": p.grid.x.iter().map(|v| v * 1e3).collect::<Vec<_>>(),
        "y": p.grid.y.iter().map(|v| v * 1e3).collect::<Vec<_>>(),
        "z": zmm,
        "layers": layers.iter().map(|(n, t, k, z)| json!({"name": n, "tag": t, "k": k, "z": z})).collect::<Vec<_>>(),
        "snapshots": snaps,
        "pml": st.opts.pml,
        "units": {"length": "mm", "E": "V/m per V of source spectrum", "H": "A/m", "J": "A/m"},
        "kicad_y_flipped": true,
    });
    std::fs::write(dir.join("meta.json"), serde_json::to_string_pretty(&meta).unwrap())
}

/// Time-domain surface-current frames: `out/frames/J_<layer>_<NNNN>.f32` (Jx plane then Jy plane, nx * ny
/// f32 each, x fastest) and `frames.json` (grid lines, layer heights, frame steps and times, the source
/// voltage at each frame, and the largest |J| per layer for a fixed colour scale).
pub fn write_frames(out: &Path, p: &Problem, st: &Setup, layers: &[String], snaps: &[crate::fdtd::Snap], exc: &[f32]) -> std::io::Result<()> {
    let dir = out.join("frames");
    std::fs::create_dir_all(&dir)?;
    let (nx, ny) = (p.nx, p.ny);
    let mut jmax = vec![0.0f32; layers.len()];
    for (fi, sn) in snaps.iter().enumerate() {
        for (li, name) in layers.iter().enumerate() {
            let [jx, jy] = &sn.j[li];
            let mut buf = Vec::with_capacity(2 * nx * ny);
            buf.extend_from_slice(jx);
            buf.extend_from_slice(jy);
            for q in 0..nx * ny {
                jmax[li] = jmax[li].max((jx[q] * jx[q] + jy[q] * jy[q]).sqrt());
            }
            write_f32(&dir.join(format!("J_{}_{:04}.f32", name.replace('.', "_"), fi)), &buf)?;
        }
    }
    let zmm: Vec<f64> = p.grid.z.iter().map(|v| v * 1e3).collect();
    let meta = json!({
        "shape": [nx, ny],
        "x": p.grid.x.iter().map(|v| v * 1e3).collect::<Vec<_>>(),
        "y": p.grid.y.iter().map(|v| v * 1e3).collect::<Vec<_>>(),
        "layers": layers.iter().zip(&jmax).map(|(n, m)| {
            let s = st.layout.sheets.iter().find(|s| &s.name == n).unwrap();
            json!({"name": n, "tag": n.replace('.', "_"), "z": s.z, "k": crate::mesh::snap(&zmm, s.z), "j_max": m})
        }).collect::<Vec<_>>(),
        "frames": snaps.iter().map(|s| json!({"step": s.step, "t_ps": s.step as f64 * p.dt * 1e12, "vs": exc.get(s.step).copied().unwrap_or(0.0)})).collect::<Vec<_>>(),
        "dt": p.dt,
        "pml": st.opts.pml,
    });
    std::fs::write(dir.join("frames.json"), serde_json::to_string_pretty(&meta).unwrap())
}
