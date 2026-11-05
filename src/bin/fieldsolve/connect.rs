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
use std::collections::{HashMap, VecDeque};

#[derive(Clone, Debug)]
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

pub fn check(g: &Grid, lay: &Layout) -> Vec<Bridge> {
    let (nx, ny, nz) = (g.nx(), g.ny(), g.nz());
    let n = g.n();
    let (masks, names) = sheet_masks(g, lay);
    // adjacency through metal: in-plane edges carry a net label (0 = none); barrel z-edges are unlabelled
    let mut ex: HashMap<usize, u32> = HashMap::new();
    let mut ey: HashMap<usize, u32> = HashMap::new();
    let wild: Vec<bool> = names.iter().map(|n| n.starts_with('~')).collect();
    for m in &masks {
        for j in 0..ny {
            for i in 0..nx {
                let q = j * nx + i;
                // element copper (an intended short) is not part of the layout's connectivity
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
    let mut ez: std::collections::HashSet<usize> = std::collections::HashSet::new();
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
    // seed: both end nodes of every labelled edge; copper of two nets meeting at a node is a raster short
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
    // flood through unlabelled conductors (net-less copper and barrels)
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
    let mut v: Vec<Bridge> = pairs.into_values().collect();
    v.sort_by(|a, b| b.count.cmp(&a.count));
    v
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
