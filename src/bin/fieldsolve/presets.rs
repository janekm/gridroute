//! Settings presets and set-up checks (docs/fieldsolve-settings.md).
//!
//! `"preset": "rf"` or `"antenna"` fills in the recommended model settings that a spec leaves out (keys the
//! spec sets itself always win). Every run then prints one line naming the model it actually uses and
//! warns about settings that are known traps:
//! - `"loss": {"copper": true}` is a plane-pair stand-in (an extra tan delta in the dielectric), not
//!   conductor loss;
//! - a far-field run without `"conductor"` treats copper as lossless, so its efficiency is optimistic;
//! - narrow copper (traces, meanders, antenna arms) in coarse cells next to the fine region (`roi`).
//!   This happened on the MTH PIFA: the open end sat in 0.35-0.5 mm cells and the resonance came out
//!   about 3 % low.
use crate::geom::{BBox, P2};
use serde_json::{json, Value};

pub const PRESETS: &[&str] = &["rf", "antenna"];

/// Fill the keys `spec["preset"]` recommends and record what was filled in `spec["preset_record"]`.
pub fn apply(mut spec: Value) -> Result<Value, String> {
    let Some(name) = spec.get("preset").and_then(|v| v.as_str()).map(str::to_string) else {
        return Ok(spec);
    };
    let mut fill: Vec<(&str, Value)> = vec![
        ("loss", json!({"model": "wideband"})),
        ("conductor", json!(true)),
        ("edges", json!(true)),
        ("mesh", json!({"no_slivers": true, "z_grade": 3.0})),
    ];
    match name.as_str() {
        "rf" => {}
        "antenna" => fill.push(("extent", json!("board"))),
        o => return Err(format!("unknown \"preset\" {o:?} (have {PRESETS:?})")),
    }
    let mut filled = Vec::new();
    for (k, v) in fill {
        if spec.get(k).is_none() {
            spec[k] = v.clone();
            filled.push(json!({"key": k, "value": v}));
        } else if k == "loss" && spec["loss"].get("model").is_none() {
            spec["loss"]["model"] = json!("wideband");
            filled.push(json!({"key": "loss.model", "value": "wideband"}));
        }
    }
    spec["preset_record"] = json!({"preset": name, "filled": filled});
    Ok(spec)
}

fn on(v: Option<&Value>) -> bool {
    v.is_some_and(|v| v.as_bool() != Some(false) && !v.is_null())
}

/// The model a spec resolves to, as JSON (result.json "settings") and as one printable line.
pub fn summary(spec: &Value) -> (Value, String) {
    let num = |v: &Value, k: &str| v.get(k).and_then(|x| x.as_f64());
    let loss = spec.get("loss");
    let wide = loss.and_then(|l| l.get("model")).and_then(|m| m.as_str()) == Some("wideband");
    let copper_standin = loss.and_then(|l| l.get("copper")).and_then(|c| c.as_bool()).unwrap_or(false);
    let dielectric = if wide {
        let l = loss.unwrap();
        format!("wideband tan d ({:.2}-{:.1} GHz)", num(l, "f_lo").unwrap_or(f64::NAN) / 1e9, num(l, "f_hi").unwrap_or(f64::NAN) / 1e9)
            .replace("NaN", "default")
    } else {
        format!("tan d exact only at f_loss {:.2} GHz", num(spec, "f_loss").unwrap_or(1e9) / 1e9)
    };
    let cond = spec.get("conductor");
    let copper = if on(cond) {
        let c = cond.unwrap();
        format!("surface impedance (sigma {:.2e} S/m, roughness {} um)", num(c, "sigma").unwrap_or(5.8e7), num(c, "roughness_um").unwrap_or(1.0))
    } else if copper_standin {
        "perfect conductor + plane-pair stand-in in the dielectric (\"loss.copper\")".to_string()
    } else {
        "perfect conductor (lossless)".to_string()
    };
    let edges = match spec.get("edges") {
        Some(Value::Bool(true)) => "conformal + singular".to_string(),
        Some(v @ Value::Object(_)) => format!("custom {}", v),
        _ => "off (staircase)".to_string(),
    };
    let mesh = match spec.get("mesh") {
        Some(v) if !v.is_null() => v.to_string(),
        _ => "default".to_string(),
    };
    let extent = spec.get("extent").and_then(|v| v.as_str()).unwrap_or("roi");
    let roi = spec.get("roi").cloned().unwrap_or(Value::Null);
    let cell = num(spec, "cell").unwrap_or(0.1);
    let preset = spec.pointer("/preset_record/preset").and_then(|v| v.as_str()).unwrap_or("none");
    let line = format!(
        "model: preset {preset} | dielectric: {dielectric} | copper: {copper} | edges: {edges} | mesh: {mesh} | fine region {roi} at {cell} mm, extent {extent}"
    );
    let v = json!({"preset": preset, "preset_filled": spec.pointer("/preset_record/filled").cloned().unwrap_or(Value::Null),
                   "dielectric": dielectric, "copper": copper, "edges": edges, "mesh": mesh, "extent": extent, "roi": roi, "cell": cell});
    (v, line)
}

/// Warnings about settings that are known traps (independent of the geometry).
pub fn spec_checks(spec: &Value) -> Vec<String> {
    let mut w = Vec::new();
    let far = spec.get("farfield").is_some_and(|v| !v.is_null());
    let loss = spec.get("loss");
    let wide = loss.and_then(|l| l.get("model")).and_then(|m| m.as_str()) == Some("wideband");
    let standin = loss.and_then(|l| l.get("copper")).and_then(|c| c.as_bool()).unwrap_or(false);
    let cond = on(spec.get("conductor"));
    if standin {
        w.push(
            "\"loss\": {\"copper\": true} is a plane-pair stand-in: it adds about delta_s/h to the dielectric's tan delta \
             and leaves the copper lossless. It is not conductor loss. Use \"conductor\": true (surface impedance) \
             for traces, antennas and anything that radiates."
                .to_string(),
        );
    }
    if far && !cond {
        w.push(
            "far-field run with lossless copper: radiation efficiency and gain will be optimistic. Add \
             \"conductor\": true, or use \"preset\": \"antenna\"."
                .to_string(),
        );
    }
    if far && !wide {
        w.push(
            "far-field run without \"loss\": {\"model\": \"wideband\"}: tan delta is right only near f_loss, so the \
             substrate loss is wrong across the band."
                .to_string(),
        );
    }
    if far && spec.get("roi").is_some() && spec.get("extent").and_then(|v| v.as_str()) != Some("board") {
        w.push(
            "far-field run with the domain cropped to \"roi\": the rest of the board's ground (which radiates on a \
             small board) is missing. Use \"extent\": \"board\"."
                .to_string(),
        );
    }
    if far && spec.get("edges").is_none() {
        w.push("far-field run without \"edges\": copper edges are staircased (see docs/fieldsolve-settings.md).".to_string());
    }
    if spec.pointer("/preset_record/preset").and_then(|v| v.as_str()) == Some("antenna") && !far {
        w.push("\"preset\": \"antenna\" without \"farfield\": no radiation pattern or efficiency will be computed.".to_string());
    }
    w
}

/// Narrow copper (features under `narrow` mm across in x or y) that lies outside `roi`, within `reach` mm
/// of it, in mesh cells coarser than `coarse` mm. Returns (layer, bbox, area mm², largest local cell mm),
/// nearest to the roi first. `gx`/`gy` are the mesh lines in mm.
pub fn coarse_copper(
    sheets: &[crate::model::Sheet],
    gx: &[f64],
    gy: &[f64],
    roi: &BBox,
    dom: &BBox,
    cell: f64,
    coarse: f64,
    narrow: f64,
    reach: f64,
) -> Vec<(String, BBox, f64, f64)> {
    let pitch = cell.max(0.025);
    let g = roi.grow(reach);
    let reg = BBox { x0: g.x0.max(dom.x0), y0: g.y0.max(dom.y0), x1: g.x1.min(dom.x1), y1: g.y1.min(dom.y1) };
    if reg.is_empty() || gx.len() < 2 || gy.len() < 2 {
        return Vec::new();
    }
    let nx = ((reg.x1 - reg.x0) / pitch).ceil() as usize;
    let ny = ((reg.y1 - reg.y0) / pitch).ceil() as usize;
    if nx == 0 || ny == 0 || nx * ny > 40_000_000 {
        return Vec::new();
    }
    let step = |l: &[f64], v: f64| {
        let i = l.partition_point(|&a| a <= v).clamp(1, l.len() - 1);
        l[i] - l[i - 1]
    };
    let xs: Vec<f64> = (0..nx).map(|i| reg.x0 + (i as f64 + 0.5) * pitch).collect();
    let ys: Vec<f64> = (0..ny).map(|j| reg.y0 + (j as f64 + 0.5) * pitch).collect();
    let hx: Vec<f64> = xs.iter().map(|&x| step(gx, x)).collect();
    let hy: Vec<f64> = ys.iter().map(|&y| step(gy, y)).collect();
    let mut out = Vec::new();
    for sh in sheets {
        // primitives binned into 1 mm buckets, kept in file order so clears apply in sequence
        let bs = 1.0;
        let (bx, by) = (((reg.x1 - reg.x0) / bs).ceil() as usize + 1, ((reg.y1 - reg.y0) / bs).ceil() as usize + 1);
        let mut bins: Vec<Vec<u32>> = vec![Vec::new(); bx * by];
        for (k, p) in sh.prims.iter().enumerate() {
            let b = p.shape.bbox();
            if !b.overlaps(&reg) {
                continue;
            }
            let i0 = (((b.x0 - reg.x0) / bs).floor().max(0.0) as usize).min(bx - 1);
            let i1 = (((b.x1 - reg.x0) / bs).floor().max(0.0) as usize).min(bx - 1);
            let j0 = (((b.y0 - reg.y0) / bs).floor().max(0.0) as usize).min(by - 1);
            let j1 = (((b.y1 - reg.y0) / bs).floor().max(0.0) as usize).min(by - 1);
            for j in j0..=j1 {
                for i in i0..=i1 {
                    bins[j * bx + i].push(k as u32);
                }
            }
        }
        let mut cu = vec![false; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let p = P2::new(xs[i], ys[j]);
                let b = &bins[((ys[j] - reg.y0) / bs) as usize * bx + ((xs[i] - reg.x0) / bs) as usize];
                let mut v = false;
                for &k in b {
                    let pr = &sh.prims[k as usize];
                    if pr.dark != v && pr.shape.contains(p) {
                        v = pr.dark;
                    }
                }
                cu[j * nx + i] = v;
            }
        }
        // run lengths of copper along x and y
        let mut rx = vec![0u32; nx * ny];
        let mut ry = vec![0u32; nx * ny];
        for j in 0..ny {
            let mut i = 0;
            while i < nx {
                if !cu[j * nx + i] {
                    i += 1;
                    continue;
                }
                let s = i;
                while i < nx && cu[j * nx + i] {
                    i += 1;
                }
                (s..i).for_each(|q| rx[j * nx + q] = (i - s) as u32);
            }
        }
        for i in 0..nx {
            let mut j = 0;
            while j < ny {
                if !cu[j * nx + i] {
                    j += 1;
                    continue;
                }
                let s = j;
                while j < ny && cu[j * nx + i] {
                    j += 1;
                }
                (s..j).for_each(|q| ry[q * nx + i] = (j - s) as u32);
            }
        }
        let lim = (narrow / pitch).round().max(1.0) as u32;
        let flag: Vec<bool> = (0..nx * ny)
            .map(|q| {
                let (i, j) = (q % nx, q / nx);
                cu[q] && rx[q].min(ry[q]) <= lim && hx[i].max(hy[j]) > coarse && !roi.contains(P2::new(xs[i], ys[j]))
            })
            .collect();
        // 8-connected clusters
        let mut seen = vec![false; nx * ny];
        let mut stack = Vec::new();
        for q0 in 0..nx * ny {
            if !flag[q0] || seen[q0] {
                continue;
            }
            seen[q0] = true;
            stack.push(q0);
            let (mut n, mut bb, mut hmax) = (0usize, BBox::empty(), 0f64);
            while let Some(q) = stack.pop() {
                let (i, j) = (q % nx, q / nx);
                n += 1;
                bb.add(P2::new(xs[i] - pitch / 2.0, ys[j] - pitch / 2.0));
                bb.add(P2::new(xs[i] + pitch / 2.0, ys[j] + pitch / 2.0));
                hmax = hmax.max(hx[i].max(hy[j]));
                for dj in -1i64..=1 {
                    for di in -1i64..=1 {
                        let (a, b) = (i as i64 + di, j as i64 + dj);
                        if a < 0 || b < 0 || a >= nx as i64 || b >= ny as i64 {
                            continue;
                        }
                        let r = b as usize * nx + a as usize;
                        if flag[r] && !seen[r] {
                            seen[r] = true;
                            stack.push(r);
                        }
                    }
                }
            }
            let area = n as f64 * pitch * pitch;
            if area >= 0.25 {
                out.push((sh.name.clone(), bb, area, hmax));
            }
        }
    }
    let gap = |b: &BBox| {
        let dx = (roi.x0 - b.x1).max(b.x0 - roi.x1).max(0.0);
        let dy = (roi.y0 - b.y1).max(b.y0 - roi.y1).max(0.0);
        dx.hypot(dy)
    };
    out.sort_by(|a, b| gap(&a.1).partial_cmp(&gap(&b.1)).unwrap().then(b.2.partial_cmp(&a.2).unwrap()));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Prim, Shape};

    #[test]
    fn preset_fills_only_missing_keys() {
        let s = apply(json!({"preset": "antenna", "edges": false, "loss": {"f_lo": 0.5e9}})).unwrap();
        assert_eq!(s["conductor"], json!(true));
        assert_eq!(s["edges"], json!(false), "a key the spec sets wins");
        assert_eq!(s["loss"]["model"], json!("wideband"));
        assert_eq!(s["loss"]["f_lo"], json!(0.5e9));
        assert_eq!(s["extent"], json!("board"));
        assert!(apply(json!({"preset": "nope"})).is_err());
        let r = apply(json!({"preset": "rf"})).unwrap();
        assert!(r.get("extent").is_none());
    }

    #[test]
    fn copper_standin_and_lossless_far_field_warn() {
        let w = spec_checks(&json!({"farfield": {"freqs": [2.4e9]}, "loss": {"model": "wideband", "copper": true}, "roi": [0, 0, 1, 1]}));
        assert!(w.iter().any(|m| m.contains("plane-pair stand-in")));
        assert!(w.iter().any(|m| m.contains("lossless copper")));
        assert!(w.iter().any(|m| m.contains("cropped")));
        let ok = spec_checks(&apply(json!({"preset": "antenna", "farfield": {"freqs": [2.4e9]}, "roi": [0, 0, 1, 1]})).unwrap());
        assert!(ok.is_empty(), "{ok:?}");
        let (_, line) = summary(&apply(json!({"preset": "antenna"})).unwrap());
        assert!(line.contains("surface impedance") && line.contains("conformal"), "{line}");
    }

    #[test]
    fn narrow_copper_in_coarse_cells_is_found_and_planes_are_not() {
        // roi 0..10 x 0..10 at 0.05 mm; outside it the mesh grades to 0.5 mm. A 0.4 mm wide arm at x -4..-1
        // and a big plane at y -12..-6 lie outside the roi.
        let rect = |x0: f64, y0: f64, x1: f64, y1: f64| Prim {
            shape: Shape::Polygon(vec![P2::new(x0, y0), P2::new(x1, y0), P2::new(x1, y1), P2::new(x0, y1)]),
            dark: true,
            net: None,
            pad: None,
        };
        let arm = rect(-4.0, 4.8, -1.0, 5.2);
        let plane = rect(0.0, -12.0, 10.0, -6.0);
        let sheets = vec![crate::model::Sheet { name: "F.Cu".into(), z: 0.0, prims: vec![arm, plane] }];
        let mut gx: Vec<f64> = (0..=200).map(|i| i as f64 * 0.05).collect();
        let mut gy = gx.clone();
        for k in 1..=40 {
            gx.insert(0, -0.5 * k as f64);
            gy.insert(0, -0.5 * k as f64);
        }
        gx.extend((1..=20).map(|k| 10.0 + 0.5 * k as f64));
        gy.extend((1..=20).map(|k| 10.0 + 0.5 * k as f64));
        let roi = BBox { x0: 0.0, y0: 0.0, x1: 10.0, y1: 10.0 };
        let dom = BBox { x0: -20.0, y0: -20.0, x1: 20.0, y1: 20.0 };
        let c = coarse_copper(&sheets, &gx, &gy, &roi, &dom, 0.05, 0.125, 1.0, 8.0);
        assert_eq!(c.len(), 1, "{c:?}");
        let (_, b, area, h) = &c[0];
        assert!(b.x0 < -3.9 && b.x1 > -1.1 && (area - 1.2).abs() < 0.2 && *h > 0.4, "{c:?}");
        // the same arm inside a roi that covers it is not reported
        let roi2 = BBox { x0: -5.0, y0: 0.0, x1: 10.0, y1: 10.0 };
        assert!(coarse_copper(&sheets, &gx, &gy, &roi2, &dom, 0.05, 0.125, 1.0, 8.0).is_empty());
    }
}
