//! Mesh connectivity check: does the voxel model connect copper of different nets?
//!
//! Copper is rasterised onto mesh edges, and two copper features that are closer than about a cell can end
//! up sharing a mesh node. A via barrel also connects to plane copper whose metal edges reach its nodes
//! through a too-small antipad. Either short silently changes every impedance the solver reports. This
//! check labels each metal edge with the Gerber net of the copper that painted it and floods the labels
//! through the conductor graph (copper edges and via barrels, not the intended lumped shorts). Wherever two
//! labels meet, the model has a short between those nets; it is reported with its position.

use crate::geom::P2;
use crate::mesh::Grid;
use crate::model::{barrel_nodes, sheet_masks, Layout};
use std::collections::VecDeque;

#[derive(Clone, Debug, PartialEq)]
pub struct Bridge {
    pub a: String,
    pub b: String,
    /// Number of meeting nodes found for this pair.
    pub count: usize,
    /// One meeting node, mm, Gerber frame.
    pub at: [f64; 3],
    /// Through a barrel (true) or copper touching in the raster (false).
    pub via: bool,
}

/// No metal edge here (in the per-plane net tables).
const ABSENT: u32 = u32::MAX;

pub fn check(g: &Grid, lay: &Layout) -> Vec<Bridge> {
    labels(g, lay).1
}

/// The flood: the net label of every mesh node (0 = none) and the bridges, most meeting nodes first (ties
/// by net pair). Deterministic: edges seed in increasing node index, x edges before y edges.
///
/// Everything is indexed by dense node ids: the in-plane metal edges as net tables over the few mesh planes
/// that carry copper, barrel z-edges as a bitset, the node labels (which double as the visited set) as a
/// vector; only the bridges (a handful of net pairs) are hashed.
pub(crate) fn labels(g: &Grid, lay: &Layout) -> (Vec<u32>, Vec<Bridge>) {
    let (nx, ny, nz) = (g.nx(), g.ny(), g.nz());
    let n = g.n();
    let plane = nx * ny;
    let (masks, names) = sheet_masks(g, lay);
    // in-plane edges through metal, per copper plane: the edge's net label (0 = none) or ABSENT
    let mut slot = vec![ABSENT; nz];
    let mut nplanes = 0usize;
    for m in &masks {
        if slot[m.k] == ABSENT {
            slot[m.k] = nplanes as u32;
            nplanes += 1;
        }
    }
    let mut exn = vec![ABSENT; nplanes * plane];
    let mut eyn = vec![ABSENT; nplanes * plane];
    let wild: Vec<bool> = names.iter().map(|n| n.starts_with('~')).collect();
    for m in &masks {
        let base = slot[m.k] as usize * plane;
        for j in 0..ny {
            for i in 0..nx {
                let q = j * nx + i;
                // element copper (an intended short) is not part of the layout's connectivity; a later sheet
                // on the same plane overrides an earlier one's label
                if m.ex[q] && !wild[m.exn[q] as usize] && i + 1 < nx {
                    exn[base + q] = m.exn[q];
                }
                if m.ey[q] && !wild[m.eyn[q] as usize] && j + 1 < ny {
                    eyn[base + q] = m.eyn[q];
                }
            }
        }
    }
    // the +x / +y edge at node id
    let edge = |t: &[u32], id: usize| -> u32 {
        match slot[id / plane] {
            ABSENT => ABSENT,
            s => t[s as usize * plane + id % plane],
        }
    };
    // barrel z-edges (node id to id + plane), unlabelled
    let mut ez = vec![0u64; n.div_ceil(64)];
    for b in &lay.barrels {
        let (cells, k0, k1) = barrel_nodes(g, b);
        for (i, j) in cells {
            for k in k0..k1 {
                let id = g.idx(i, j, k);
                ez[id / 64] |= 1 << (id % 64);
            }
        }
    }
    let has_ez = |id: usize| ez[id / 64] >> (id % 64) & 1 != 0;
    let mut label = vec![0u32; n];
    let mut pairs: rustc_hash::FxHashMap<(u32, u32), Bridge> = Default::default();
    let pos = |id: usize| -> [f64; 3] {
        let (i, j, k) = (id % nx, (id / nx) % ny, id / plane);
        [g.x[i] * 1e3, g.y[j] * 1e3, g.z[k] * 1e3]
    };
    let mut meet = |a: u32, b: u32, id: usize, via: bool| {
        let key = (a.min(b), a.max(b));
        let e = pairs.entry(key).or_insert_with(|| Bridge { a: names[key.0 as usize].clone(), b: names[key.1 as usize].clone(), count: 0, at: pos(id), via });
        e.count += 1;
    };
    let mut queue = VecDeque::new();
    // seed: both end nodes of every labelled edge; copper of two nets meeting at a node is a raster short
    for (t, step) in [(&exn, 1usize), (&eyn, nx)] {
        for k in (0..nz).filter(|&k| slot[k] != ABSENT) {
            let base = slot[k] as usize * plane;
            for q in 0..plane {
                let net = t[base + q];
                if net == ABSENT || net == 0 {
                    continue;
                }
                let id = k * plane + q;
                for node in [id, id + step] {
                    if label[node] == 0 {
                        label[node] = net;
                        queue.push_back(node);
                    } else if label[node] != net {
                        meet(label[node], net, node, false);
                    }
                }
            }
        }
    }
    // flood through unlabelled conductors (net-less copper and barrels)
    let mut nb = [(0usize, false); 6];
    while let Some(id) = queue.pop_front() {
        let l = label[id];
        let (i, j, k) = (id % nx, (id / nx) % ny, id / plane);
        let mut cnt = 0;
        let mut push = |m: usize, via: bool| {
            nb[cnt] = (m, via);
            cnt += 1;
        };
        if i + 1 < nx && edge(&exn, id) != ABSENT {
            push(id + 1, false);
        }
        if i > 0 && edge(&exn, id - 1) != ABSENT {
            push(id - 1, false);
        }
        if j + 1 < ny && edge(&eyn, id) != ABSENT {
            push(id + nx, false);
        }
        if j > 0 && edge(&eyn, id - nx) != ABSENT {
            push(id - nx, false);
        }
        if k + 1 < nz && has_ez(id) {
            push(id + plane, true);
        }
        if k > 0 && has_ez(id - plane) {
            push(id - plane, true);
        }
        for &(m, via) in &nb[..cnt] {
            if label[m] == 0 {
                label[m] = l;
                queue.push_back(m);
            } else if label[m] != l {
                meet(l, label[m], m, via);
            }
        }
    }
    let mut v: Vec<((u32, u32), Bridge)> = pairs.into_iter().collect();
    v.sort_by(|a, b| b.1.count.cmp(&a.1.count).then(a.0.cmp(&b.0)));
    (label, v.into_iter().map(|x| x.1).collect())
}

/// Report lines for the bridges, coordinates in the user's frame.
pub fn describe(b: &[Bridge], kicad: bool) -> Vec<String> {
    b.iter()
        .map(|x| {
            let p = P2::new(x.at[0], if kicad { -x.at[1] } else { x.at[1] });
            format!(
                "nets {:?} and {:?} meet at {} mesh node(s), e.g. ({:.3}, {:.3}) mm, z = {:.3} mm{}",
                x.a,
                x.b,
                x.count,
                p.x,
                p.y,
                x.at[2],
                if x.via { " (through a via barrel)" } else { " (copper touching in the raster)" }
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Prim, Shape};
    use crate::model::{Barrel, Sheet};
    use std::collections::{BTreeMap, BTreeSet, HashMap};

    /// The original traversal (hash maps of node ids, a `Vec` of neighbours per node), with ordered maps so
    /// that it seeds in increasing node index like `labels`: the reference the dense version must match.
    fn reference(g: &Grid, lay: &Layout) -> (Vec<u32>, Vec<Bridge>) {
        let (nx, ny, nz) = (g.nx(), g.ny(), g.nz());
        let n = g.n();
        let (masks, names) = sheet_masks(g, lay);
        let mut ex: BTreeMap<usize, u32> = BTreeMap::new();
        let mut ey: BTreeMap<usize, u32> = BTreeMap::new();
        let wild: Vec<bool> = names.iter().map(|n| n.starts_with('~')).collect();
        for m in &masks {
            for j in 0..ny {
                for i in 0..nx {
                    let q = j * nx + i;
                    if m.ex[q] && wild[m.exn[q] as usize] || m.ey[q] && wild[m.eyn[q] as usize] {
                        if m.ex[q] && !wild[m.exn[q] as usize] && i + 1 < nx {
                            ex.insert(g.idx(i, j, m.k), m.exn[q]);
                        }
                        if m.ey[q] && !wild[m.eyn[q] as usize] && j + 1 < ny {
                            ey.insert(g.idx(i, j, m.k), m.eyn[q]);
                        }
                        continue;
                    }
                    if m.ex[q] && i + 1 < nx {
                        ex.insert(g.idx(i, j, m.k), m.exn[q]);
                    }
                    if m.ey[q] && j + 1 < ny {
                        ey.insert(g.idx(i, j, m.k), m.eyn[q]);
                    }
                }
            }
        }
        let mut ez: BTreeSet<usize> = BTreeSet::new();
        for b in &lay.barrels {
            let (cells, k0, k1) = barrel_nodes(g, b);
            for (i, j) in cells {
                for k in k0..k1 {
                    ez.insert(g.idx(i, j, k));
                }
            }
        }
        let mut label = vec![0u32; n];
        let mut pairs: HashMap<(u32, u32), Bridge> = HashMap::new();
        let pos = |id: usize| -> [f64; 3] {
            let (i, j, k) = (id % nx, (id / nx) % ny, id / (nx * ny));
            [g.x[i] * 1e3, g.y[j] * 1e3, g.z[k] * 1e3]
        };
        let meet = |a: u32, b: u32, id: usize, via: bool, pairs: &mut HashMap<(u32, u32), Bridge>| {
            let key = (a.min(b), a.max(b));
            let e = pairs.entry(key).or_insert_with(|| Bridge { a: names[key.0 as usize].clone(), b: names[key.1 as usize].clone(), count: 0, at: pos(id), via });
            e.count += 1;
        };
        let mut queue = VecDeque::new();
        for (map, step) in [(&ex, 1usize), (&ey, nx)] {
            for (&id, &net) in map.iter() {
                if net == 0 {
                    continue;
                }
                for node in [id, id + step] {
                    if label[node] == 0 {
                        label[node] = net;
                        queue.push_back(node);
                    } else if label[node] != net {
                        meet(label[node], net, node, false, &mut pairs);
                    }
                }
            }
        }
        let nbrs = |id: usize| -> Vec<(usize, bool)> {
            let (i, j, k) = (id % nx, (id / nx) % ny, id / (nx * ny));
            let mut v = Vec::with_capacity(6);
            if i + 1 < nx && ex.contains_key(&id) {
                v.push((id + 1, false));
            }
            if i > 0 && ex.contains_key(&(id - 1)) {
                v.push((id - 1, false));
            }
            if j + 1 < ny && ey.contains_key(&id) {
                v.push((id + nx, false));
            }
            if j > 0 && ey.contains_key(&(id - nx)) {
                v.push((id - nx, false));
            }
            if k + 1 < nz && ez.contains(&id) {
                v.push((id + nx * ny, true));
            }
            if k > 0 && ez.contains(&(id - nx * ny)) {
                v.push((id - nx * ny, true));
            }
            v
        };
        while let Some(id) = queue.pop_front() {
            let l = label[id];
            for (m, via) in nbrs(id) {
                if label[m] == 0 {
                    label[m] = l;
                    queue.push_back(m);
                } else if label[m] != l {
                    meet(l, label[m], m, via, &mut pairs);
                }
            }
        }
        let mut v: Vec<((u32, u32), Bridge)> = pairs.into_iter().collect();
        v.sort_by(|a, b| b.1.count.cmp(&a.1.count).then(a.0.cmp(&b.0)));
        (label, v.into_iter().map(|x| x.1).collect())
    }

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64, net: Option<&str>) -> Prim {
        Prim { shape: Shape::Polygon(vec![P2::new(x0, y0), P2::new(x1, y0), P2::new(x1, y1), P2::new(x0, y1)]), dark: true, net: net.map(String::from), pad: None }
    }

    /// Three planes: a top layer with pads of nets A, B, C, net-less copper and element copper, a GND plane
    /// and a bottom layer of net D, joined by barrels. The raster keeps labelled nets apart (antipads and
    /// touching copper are resolved in `sheet_masks`); net-less copper joins them: A and B through a strip
    /// across their gap, and C through a patch and its barrel to the planes.
    fn shorted_board() -> (Grid, Layout) {
        let l = |a: f64, b: f64, h: f64| crate::mesh::lines_1d(a, b, &[], &[crate::mesh::Refine { a, b, h }], 1.3, h, &[]).into_iter().map(|v| v * 1e-3).collect();
        let g = Grid { x: l(-3.0, 3.0, 0.1), y: l(-2.0, 2.0, 0.1), z: vec![-0.4e-3, -0.2e-3, 0.0, 0.2e-3, 0.4e-3, 0.6e-3] };
        let top = vec![
            rect(-2.5, -1.5, -1.0, 1.5, Some("A")),
            rect(-0.95, -0.5, 0.5, 0.5, Some("B")),
            rect(-1.3, -0.2, -0.6, 0.2, None), // net-less strip across the A-B gap
            rect(0.6, 0.5, 2.5, 1.5, Some("C")),
            rect(0.6, -1.5, 2.5, 0.6, None),        // net-less patch overlapping C
            rect(2.0, -0.6, 2.2, 0.6, Some("~E1")), // element copper (an intended short): no connectivity
        ];
        let gnd = vec![rect(-3.0, -2.0, 3.0, 2.0, Some("GND"))];
        let bot = vec![rect(-2.5, -1.5, 2.5, 1.5, Some("D"))];
        let lay = Layout {
            sheets: vec![Sheet { name: "top".into(), z: 0.4, prims: top }, Sheet { name: "gnd".into(), z: 0.0, prims: gnd }, Sheet { name: "bot".into(), z: -0.4, prims: bot }],
            barrels: vec![
                Barrel { c: P2::new(-2.0, 0.0), r: 0.15, z0: -0.4, z1: 0.4 }, // A through GND to D (resolved: no short)
                Barrel { c: P2::new(1.2, -1.0), r: 0.1, z0: -0.4, z1: 0.4 },  // net-less patch to GND and D
                Barrel { c: P2::new(-1.6, -1.0), r: 0.1, z0: 0.0, z1: 0.4 },  // A to GND
                Barrel { c: P2::new(1.5, 0.4), r: 0.04, z0: 0.0, z1: 0.4 },   // patch next to C, to GND: meets in the barrel
            ],
            ..Default::default()
        };
        (g, lay)
    }

    #[test]
    fn dense_traversal_matches_the_reference() {
        let (g, lay) = shorted_board();
        let (la, ba) = labels(&g, &lay);
        let (lb, bb) = reference(&g, &lay);
        // A-B through the strip; C-GND through the patch, its barrels and (counted) meetings inside a barrel
        let pairs: Vec<(&str, &str)> = ba.iter().map(|b| (b.a.as_str(), b.b.as_str())).collect();
        assert!(pairs.contains(&("A", "B")) && pairs.contains(&("C", "GND")), "the test board's shorts: {ba:?}");
        assert_eq!(la, lb, "node labels differ");
        assert_eq!(ba, bb, "bridges differ");
        // and without the barrels (raster shorts only), and with no copper at all
        let mut lay2 = lay.clone();
        lay2.barrels.clear();
        assert_eq!(labels(&g, &lay2), reference(&g, &lay2));
        lay2.sheets.clear();
        assert_eq!(labels(&g, &lay2), reference(&g, &lay2));
    }

    /// The same comparison on a real board: `FIELDSOLVE_CONNECT_SPEC=board.json cargo test --release
    /// --bin fieldsolve real_board_matches -- --ignored`.
    #[test]
    #[ignore]
    fn real_board_matches_the_reference() {
        let spec = std::env::var("FIELDSOLVE_CONNECT_SPEC").expect("set FIELDSOLVE_CONNECT_SPEC");
        let (mut cx, _) = crate::sim::load(std::path::Path::new(&spec)).unwrap();
        let st = crate::sim::setup(&mut cx).unwrap();
        let t = std::time::Instant::now();
        let (la, ba) = labels(&st.grid, &st.layout);
        let t1 = t.elapsed();
        let (lb, bb) = reference(&st.grid, &st.layout);
        eprintln!("{} nodes, {} labelled, {} bridges; dense {:.1} ms", la.len(), la.iter().filter(|&&l| l != 0).count(), ba.len(), t1.as_secs_f64() * 1e3);
        assert_eq!(la, lb, "node labels differ");
        assert_eq!(ba, bb, "bridges differ");
    }
}
