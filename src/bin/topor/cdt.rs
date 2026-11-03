//! Constrained Delaunay triangulation on exact integer coordinates, with point insertion and removal
//! (Bowyer-Watson cavities), fixed segments and an undo journal.
//!
//! Triangle edge `i` runs `v[i] -> v[(i+1)%3]` (counter-clockwise); `n[i]` is the neighbour across it
//! and `e[i]` the undirected edge id. Triangle and edge ids are never reused, so callers can keep
//! per-id side tables and treat dead ids as stale.

use crate::geom::{incircle, orient, IP};
use std::collections::HashMap;

pub const NONE: u32 = u32::MAX;

#[derive(Clone, Copy, Debug)]
pub struct Tri {
    pub v: [u32; 3],
    pub n: [u32; 3],
    pub e: [u32; 3],
    pub alive: bool,
}

#[derive(Clone, Copy, Debug)]
pub struct Edge {
    pub a: u32,
    pub b: u32,
    pub fixed: bool,
    pub alive: bool,
}

#[derive(Debug)]
pub enum Loc {
    Inside(u32),
    OnEdge(u32, usize),
    OnVertex(u32),
}

#[derive(Debug, Default, Clone)]
pub struct Change {
    pub removed_tris: Vec<u32>,
    pub new_tris: Vec<u32>,
    pub removed_edges: Vec<u32>,
    pub new_edges: Vec<u32>,
}

#[derive(Default, Clone)]
struct Journal {
    npts: usize,
    ntris: usize,
    nedges: usize,
    tris: Vec<(u32, Tri)>,
    edges: Vec<(u32, Edge)>,
    vtri: Vec<(u32, u32)>,
    last: u32,
}

#[derive(Clone)]
pub struct Cdt {
    pub pts: Vec<IP>,
    pub tris: Vec<Tri>,
    pub edges: Vec<Edge>,
    pub vtri: Vec<u32>,
    last: u32,
    jr: Option<Journal>,
    stamp: Vec<u32>,
    stamp_id: u32,
}

#[inline]
fn nx(i: usize) -> usize {
    (i + 1) % 3
}
#[inline]
fn pv(i: usize) -> usize {
    (i + 2) % 3
}

impl Cdt {
    /// Bounding box corners become vertices 0..4; everything is inserted inside.
    pub fn new(lo: IP, hi: IP) -> Cdt {
        let pts = vec![[lo[0], lo[1]], [hi[0], lo[1]], [hi[0], hi[1]], [lo[0], hi[1]]];
        let mut c = Cdt {
            pts,
            tris: Vec::new(),
            edges: Vec::new(),
            vtri: vec![0, 0, 0, 0],
            last: 0,
            jr: None,
            stamp: Vec::new(),
            stamp_id: 0,
        };
        // tri 0: 0,1,2  tri 1: 0,2,3 ; shared edge 0-2
        let e01 = c.new_edge(0, 1);
        let e12 = c.new_edge(1, 2);
        let e20 = c.new_edge(0, 2);
        let e23 = c.new_edge(2, 3);
        let e30 = c.new_edge(0, 3);
        c.tris.push(Tri { v: [0, 1, 2], n: [NONE, NONE, 1], e: [e01, e12, e20], alive: true });
        c.tris.push(Tri { v: [0, 2, 3], n: [0, NONE, NONE], e: [e20, e23, e30], alive: true });
        c.vtri = vec![0, 0, 0, 1];
        c.stamp = vec![0, 0];
        c
    }

    fn new_edge(&mut self, a: u32, b: u32) -> u32 {
        self.edges.push(Edge { a: a.min(b), b: a.max(b), fixed: false, alive: true });
        (self.edges.len() - 1) as u32
    }

    // ---------------- journal ----------------
    pub fn begin(&mut self) {
        self.jr = Some(Journal {
            npts: self.pts.len(),
            ntris: self.tris.len(),
            nedges: self.edges.len(),
            last: self.last,
            ..Default::default()
        });
    }
    pub fn commit(&mut self) {
        self.jr = None;
    }
    pub fn rollback(&mut self) {
        if let Some(j) = self.jr.take() {
            for (id, t) in j.tris.into_iter().rev() {
                self.tris[id as usize] = t;
            }
            for (id, e) in j.edges.into_iter().rev() {
                self.edges[id as usize] = e;
            }
            for (id, t) in j.vtri.into_iter().rev() {
                self.vtri[id as usize] = t;
            }
            self.tris.truncate(j.ntris);
            self.edges.truncate(j.nedges);
            self.pts.truncate(j.npts);
            self.vtri.truncate(j.npts);
            self.stamp.truncate(j.ntris);
            self.last = j.last;
        }
    }
    fn save_tri(&mut self, t: u32) {
        if let Some(j) = &mut self.jr {
            if (t as usize) < j.ntris {
                j.tris.push((t, self.tris[t as usize]));
            }
        }
    }
    fn save_edge(&mut self, e: u32) {
        if let Some(j) = &mut self.jr {
            if (e as usize) < j.nedges {
                j.edges.push((e, self.edges[e as usize]));
            }
        }
    }
    fn save_vtri(&mut self, v: u32) {
        if let Some(j) = &mut self.jr {
            if (v as usize) < j.npts {
                j.vtri.push((v, self.vtri[v as usize]));
            }
        }
    }

    fn next_stamp(&mut self) -> u32 {
        self.stamp_id = self.stamp_id.wrapping_add(1);
        if self.stamp_id == 0 {
            self.stamp.iter_mut().for_each(|s| *s = 0);
            self.stamp_id = 1;
        }
        if self.stamp.len() < self.tris.len() {
            self.stamp.resize(self.tris.len(), 0);
        }
        self.stamp_id
    }

    // ---------------- queries ----------------
    pub fn tri_ip(&self, t: u32) -> [IP; 3] {
        let v = self.tris[t as usize].v;
        [self.pts[v[0] as usize], self.pts[v[1] as usize], self.pts[v[2] as usize]]
    }

    pub fn any_alive(&self) -> u32 {
        if (self.last as usize) < self.tris.len() && self.tris[self.last as usize].alive {
            return self.last;
        }
        self.tris.iter().rposition(|t| t.alive).unwrap() as u32
    }

    pub fn locate(&self, p: IP, hint: u32) -> Result<Loc, String> {
        let mut t = if (hint as usize) < self.tris.len() && self.tris[hint as usize].alive { hint } else { self.any_alive() };
        let mut off = 0usize;
        for _ in 0..(4 * self.tris.len() + 64) {
            let tr = &self.tris[t as usize];
            let mut moved = false;
            for k in 0..3 {
                let i = (k + off) % 3;
                let a = self.pts[tr.v[i] as usize];
                let b = self.pts[tr.v[nx(i)] as usize];
                if orient(a, b, p) < 0 {
                    if tr.n[i] == NONE {
                        return Err("point outside triangulation".into());
                    }
                    t = tr.n[i];
                    moved = true;
                    break;
                }
            }
            off = (off + 1) % 3;
            if !moved {
                let mut zeros = [false; 3];
                for i in 0..3 {
                    let a = self.pts[tr.v[i] as usize];
                    let b = self.pts[tr.v[nx(i)] as usize];
                    zeros[i] = orient(a, b, p) == 0;
                }
                let nz = zeros.iter().filter(|z| **z).count();
                if nz >= 2 {
                    for i in 0..3 {
                        if self.pts[tr.v[i] as usize] == p {
                            return Ok(Loc::OnVertex(tr.v[i]));
                        }
                    }
                }
                for i in 0..3 {
                    if self.pts[tr.v[i] as usize] == p {
                        return Ok(Loc::OnVertex(tr.v[i]));
                    }
                }
                if nz >= 1 {
                    let i = zeros.iter().position(|z| *z).unwrap();
                    return Ok(Loc::OnEdge(t, i));
                }
                return Ok(Loc::Inside(t));
            }
        }
        Err("locate did not converge".into())
    }

    /// Triangles incident to vertex v.
    pub fn star(&self, v: u32) -> Vec<u32> {
        let mut out = Vec::new();
        let t0 = self.vtri[v as usize];
        if t0 == NONE {
            return out;
        }
        let mut stack = vec![t0];
        while let Some(t) = stack.pop() {
            if out.contains(&t) {
                continue;
            }
            let tr = &self.tris[t as usize];
            debug_assert!(tr.alive && tr.v.contains(&v));
            out.push(t);
            let i = tr.v.iter().position(|&x| x == v).unwrap();
            for &k in &[i, pv(i)] {
                let nb = tr.n[k];
                if nb != NONE && !out.contains(&nb) {
                    stack.push(nb);
                }
            }
        }
        out
    }

    /// Edge id between a and b if it exists.
    pub fn find_edge(&self, a: u32, b: u32) -> Option<u32> {
        for t in self.star(a) {
            let tr = &self.tris[t as usize];
            for i in 0..3 {
                let (x, y) = (tr.v[i], tr.v[nx(i)]);
                if (x == a && y == b) || (x == b && y == a) {
                    return Some(tr.e[i]);
                }
            }
        }
        None
    }

    /// The triangle that has directed edge u->v, if any.
    pub fn tri_with_dir(&self, u: u32, v: u32) -> Option<(u32, usize)> {
        for t in self.star(u) {
            let tr = &self.tris[t as usize];
            for i in 0..3 {
                if tr.v[i] == u && tr.v[nx(i)] == v {
                    return Some((t, i));
                }
            }
        }
        None
    }

    // ---------------- region replacement ----------------
    /// Replace the `removed` triangles by `new` (counter-clockwise vertex triples covering the same region).
    pub fn replace_region(&mut self, removed: &[u32], new: &[[u32; 3]], fixed: &[(u32, u32)]) -> Result<Change, String> {
        let st = self.next_stamp();
        for &t in removed {
            self.stamp[t as usize] = st;
        }
        // directed boundary edges -> (outside tri, index in outside tri, edge id)
        let mut boundary: HashMap<(u32, u32), (u32, usize, u32)> = HashMap::new();
        let mut internal: Vec<u32> = Vec::new();
        for &t in removed {
            let tr = self.tris[t as usize];
            for i in 0..3 {
                let nb = tr.n[i];
                if nb == NONE || self.stamp[nb as usize] != st {
                    let k = if nb == NONE {
                        0
                    } else {
                        self.tris[nb as usize].n.iter().position(|&x| x == t).ok_or("broken adjacency")?
                    };
                    boundary.insert((tr.v[i], tr.v[nx(i)]), (nb, k, tr.e[i]));
                } else if !internal.contains(&tr.e[i]) {
                    internal.push(tr.e[i]);
                }
            }
        }
        let mut ch = Change { removed_tris: removed.to_vec(), removed_edges: internal.clone(), ..Default::default() };
        for &t in removed {
            self.save_tri(t);
            self.tris[t as usize].alive = false;
        }
        for &e in &internal {
            self.save_edge(e);
            self.edges[e as usize].alive = false;
        }
        let base = self.tris.len() as u32;
        let mut pending: HashMap<(u32, u32), (u32, usize, u32)> = HashMap::new();
        let mut used = 0usize;
        for (j, tv) in new.iter().enumerate() {
            let id = base + j as u32;
            let mut tri = Tri { v: *tv, n: [NONE; 3], e: [NONE; 3], alive: true };
            for i in 0..3 {
                let (u, w) = (tv[i], tv[nx(i)]);
                if let Some(&(nb, k, eid)) = boundary.get(&(u, w)) {
                    tri.n[i] = nb;
                    tri.e[i] = eid;
                    used += 1;
                    if nb != NONE {
                        self.save_tri(nb);
                        self.tris[nb as usize].n[k] = id;
                    }
                } else if let Some((ot, oi, eid)) = pending.remove(&(w, u)) {
                    tri.n[i] = ot;
                    tri.e[i] = eid;
                    if ot < base + j as u32 {
                        self.tris[ot as usize].n[oi] = id;
                    }
                } else {
                    let eid = self.new_edge(u, w);
                    ch.new_edges.push(eid);
                    tri.e[i] = eid;
                    pending.insert((u, w), (id, i, eid));
                }
            }
            self.tris.push(tri);
            ch.new_tris.push(id);
        }
        // pending entries whose partner came later were linked from the later side; fix the earlier side
        if !pending.is_empty() {
            return Err(format!("region retriangulation left {} unmatched edges", pending.len()));
        }
        if used != boundary.len() {
            return Err(format!("region boundary mismatch: used {} of {}", used, boundary.len()));
        }
        // second pass: links where the later triangle found the earlier one already set n on itself only
        for j in 0..new.len() {
            let id = base + j as u32;
            for i in 0..3 {
                let ot = self.tris[id as usize].n[i];
                if ot != NONE && ot >= base {
                    let (u, w) = (self.tris[id as usize].v[i], self.tris[id as usize].v[nx(i)]);
                    let o = &mut self.tris[ot as usize];
                    for k in 0..3 {
                        if o.v[k] == w && o.v[nx(k)] == u {
                            o.n[k] = id;
                        }
                    }
                }
            }
        }
        for &(a, b) in fixed {
            for &id in &ch.new_tris {
                let tr = self.tris[id as usize];
                for i in 0..3 {
                    let (u, w) = (tr.v[i], tr.v[nx(i)]);
                    if (u == a && w == b) || (u == b && w == a) {
                        self.edges[tr.e[i] as usize].fixed = true;
                    }
                }
            }
        }
        for &id in &ch.new_tris {
            for k in 0..3 {
                let v = self.tris[id as usize].v[k];
                self.save_vtri(v);
                self.vtri[v as usize] = id;
            }
        }
        self.stamp.resize(self.tris.len(), 0);
        if let Some(&t) = ch.new_tris.last() {
            self.last = t;
        }
        Ok(ch)
    }

    // ---------------- insertion ----------------
    /// Triangles to remove for inserting p, or the existing vertex at p.
    pub fn insert_cavity(&mut self, p: IP, hint: u32) -> Result<Result<Vec<u32>, u32>, String> {
        let loc = self.locate(p, hint)?;
        let (t0, extra) = match loc {
            Loc::OnVertex(v) => return Ok(Err(v)),
            Loc::Inside(t) => (t, NONE),
            Loc::OnEdge(t, i) => {
                let tr = self.tris[t as usize];
                if self.edges[tr.e[i] as usize].fixed {
                    return Err("point on a fixed edge".into());
                }
                (t, tr.n[i])
            }
        };
        let st = self.next_stamp();
        let mut cav = vec![t0];
        self.stamp[t0 as usize] = st;
        if extra != NONE {
            cav.push(extra);
            self.stamp[extra as usize] = st;
        }
        let mut k = 0;
        while k < cav.len() {
            let t = cav[k];
            k += 1;
            let tr = self.tris[t as usize];
            for i in 0..3 {
                if self.edges[tr.e[i] as usize].fixed {
                    continue;
                }
                let nb = tr.n[i];
                if nb == NONE || self.stamp[nb as usize] == st {
                    continue;
                }
                let [a, b, c] = self.tri_ip(nb);
                if incircle(a, b, c, p) > 0 {
                    self.stamp[nb as usize] = st;
                    cav.push(nb);
                }
            }
        }
        // keep the cavity star-shaped from p
        loop {
            let st2 = self.next_stamp();
            for &t in &cav {
                self.stamp[t as usize] = st2;
            }
            let mut bad = None;
            for (ci, &t) in cav.iter().enumerate() {
                if t == t0 || t == extra {
                    continue;
                }
                let tr = self.tris[t as usize];
                for i in 0..3 {
                    let nb = tr.n[i];
                    if nb == NONE || self.stamp[nb as usize] != st2 {
                        let (a, b) = (self.pts[tr.v[i] as usize], self.pts[tr.v[nx(i)] as usize]);
                        if orient(a, b, p) <= 0 {
                            bad = Some(ci);
                        }
                    }
                }
                if bad.is_some() {
                    break;
                }
            }
            match bad {
                Some(ci) => {
                    cav.remove(ci);
                }
                None => break,
            }
        }
        Ok(Ok(cav))
    }

    /// Insert p replacing `cavity` (from insert_cavity). Returns the new vertex.
    pub fn apply_insert(&mut self, p: IP, cavity: &[u32]) -> Result<(u32, Change), String> {
        let v = self.pts.len() as u32;
        self.pts.push(p);
        self.vtri.push(NONE);
        let st = self.next_stamp();
        for &t in cavity {
            self.stamp[t as usize] = st;
        }
        let mut new = Vec::new();
        for &t in cavity {
            let tr = self.tris[t as usize];
            for i in 0..3 {
                let nb = tr.n[i];
                if nb == NONE || self.stamp[nb as usize] != st {
                    new.push([tr.v[i], tr.v[nx(i)], v]);
                }
            }
        }
        for tv in &new {
            if orient(self.pts[tv[0] as usize], self.pts[tv[1] as usize], p) <= 0 {
                return Err("degenerate insertion".into());
            }
        }
        let ch = self.replace_region(cavity, &new, &[])?;
        Ok((v, ch))
    }

    pub fn insert(&mut self, p: IP) -> Result<u32, String> {
        let hint = self.last;
        match self.insert_cavity(p, hint)? {
            Err(v) => Ok(v),
            Ok(cav) => Ok(self.apply_insert(p, &cav)?.0),
        }
    }

    // ---------------- removal ----------------
    /// Counter-clockwise link polygon of v, from its star.
    pub fn link(&self, v: u32, star: &[u32]) -> Result<Vec<u32>, String> {
        let mut next: HashMap<u32, u32> = HashMap::new();
        for &t in star {
            let tr = &self.tris[t as usize];
            let i = tr.v.iter().position(|&x| x == v).ok_or("star")?;
            next.insert(tr.v[nx(i)], tr.v[pv(i)]);
        }
        let start = *next.keys().min().ok_or("empty star")?;
        let mut poly = vec![start];
        let mut cur = start;
        loop {
            let n = *next.get(&cur).ok_or("open link (hull vertex?)")?;
            if n == start {
                break;
            }
            poly.push(n);
            cur = n;
            if poly.len() > next.len() {
                return Err("bad link".into());
            }
        }
        if poly.len() != next.len() {
            return Err("link is not a single cycle".into());
        }
        Ok(poly)
    }

    /// Remove vertex v (kept as a dead slot) and retriangulate its link.
    pub fn apply_remove(&mut self, v: u32, star: &[u32]) -> Result<Change, String> {
        let mut poly = self.link(v, star)?;
        let mut new = Vec::new();
        while poly.len() > 3 {
            let n = poly.len();
            let mut fallback = None;
            let mut chosen = None;
            for i in 0..n {
                let (a, b, c) = (poly[(i + n - 1) % n], poly[i], poly[(i + 1) % n]);
                let (pa, pb, pc) = (self.pts[a as usize], self.pts[b as usize], self.pts[c as usize]);
                if orient(pa, pb, pc) <= 0 {
                    continue;
                }
                let mut ok = true;
                let mut delaunay = true;
                for &q in &poly {
                    if q == a || q == b || q == c {
                        continue;
                    }
                    let pq = self.pts[q as usize];
                    if orient(pa, pb, pq) >= 0 && orient(pb, pc, pq) >= 0 && orient(pc, pa, pq) >= 0 {
                        ok = false;
                        break;
                    }
                    if incircle(pa, pb, pc, pq) > 0 {
                        delaunay = false;
                    }
                }
                if !ok {
                    continue;
                }
                if delaunay {
                    chosen = Some(i);
                    break;
                }
                if fallback.is_none() {
                    fallback = Some(i);
                }
            }
            let i = chosen.or(fallback).ok_or("no ear in link polygon")?;
            let (a, b, c) = (poly[(i + n - 1) % n], poly[i], poly[(i + 1) % n]);
            new.push([a, b, c]);
            poly.remove(i);
        }
        if orient(self.pts[poly[0] as usize], self.pts[poly[1] as usize], self.pts[poly[2] as usize]) <= 0 {
            return Err("degenerate final ear".into());
        }
        new.push([poly[0], poly[1], poly[2]]);
        let ch = self.replace_region(star, &new, &[])?;
        self.save_vtri(v);
        self.vtri[v as usize] = NONE;
        Ok(ch)
    }

    // ---------------- constraints ----------------
    /// Make segment a-b an edge and mark it fixed. Vertices exactly on the segment split it.
    pub fn insert_fixed(&mut self, a: u32, b: u32) -> Result<(), String> {
        if a == b {
            return Ok(());
        }
        if let Some(e) = self.find_edge(a, b) {
            self.save_edge(e);
            self.edges[e as usize].fixed = true;
            return Ok(());
        }
        let pa = self.pts[a as usize];
        let pb = self.pts[b as usize];
        let between = |c: IP| -> bool {
            let d = [(pb[0] - pa[0]) as i128, (pb[1] - pa[1]) as i128];
            let q = [(c[0] - pa[0]) as i128, (c[1] - pa[1]) as i128];
            let dot = d[0] * q[0] + d[1] * q[1];
            dot > 0 && dot < d[0] * d[0] + d[1] * d[1]
        };
        // first triangle around a that the segment enters
        let mut start = None;
        for t in self.star(a) {
            let tr = self.tris[t as usize];
            let i = tr.v.iter().position(|&x| x == a).unwrap();
            let (u, w) = (tr.v[nx(i)], tr.v[pv(i)]);
            for &c in &[u, w] {
                let pc = self.pts[c as usize];
                if orient(pa, pb, pc) == 0 && between(pc) {
                    self.insert_fixed(a, c)?;
                    return self.insert_fixed(c, b);
                }
            }
            let ou = orient(pa, pb, self.pts[u as usize]);
            let ow = orient(pa, pb, self.pts[w as usize]);
            if ou < 0 && ow > 0 {
                start = Some((t, nx(i)));
            }
        }
        let (mut t, mut i) = start.ok_or("fixed segment: no start triangle")?;
        let mut crossed = vec![t];
        let mut left: Vec<u32> = Vec::new();
        let mut right: Vec<u32> = Vec::new();
        {
            let tr = self.tris[t as usize];
            right.push(tr.v[i]);
            left.push(tr.v[nx(i)]);
        }
        loop {
            let tr = self.tris[t as usize];
            if self.edges[tr.e[i] as usize].fixed {
                return Err("fixed segments cross".into());
            }
            let (u, w) = (tr.v[i], tr.v[nx(i)]);
            let nb = tr.n[i];
            if nb == NONE {
                return Err("fixed segment leaves triangulation".into());
            }
            let ntr = self.tris[nb as usize];
            let k = ntr.n.iter().position(|&x| x == t).unwrap();
            let q = ntr.v[pv(k)];
            crossed.push(nb);
            if q == b {
                break;
            }
            let pq = self.pts[q as usize];
            let oq = orient(pa, pb, pq);
            if oq == 0 {
                // collinear vertex: split the constraint there (no changes made yet)
                self.insert_fixed(a, q)?;
                return self.insert_fixed(q, b);
            }
            // in nb, edge k is w->u. Its other edges: k+1: u->q, k+2: q->w
            if oq > 0 {
                left.push(q);
                // segment continues between u (right) and q (left): edge u->q is index nx(k) in nb
                t = nb;
                i = nx(k);
                let _ = w;
            } else {
                right.push(q);
                t = nb;
                i = pv(k);
                let _ = u;
            }
        }
        let mut new = Vec::new();
        tri_pseudo(&self.pts, a, b, &left, &mut new);
        let rrev: Vec<u32> = right.iter().rev().cloned().collect();
        tri_pseudo(&self.pts, b, a, &rrev, &mut new);
        self.replace_region(&crossed, &new, &[(a, b)])?;
        Ok(())
    }

    /// Check structural invariants (tests and debugging).
    pub fn validate(&self) -> Result<(), String> {
        for (ti, t) in self.tris.iter().enumerate() {
            if !t.alive {
                continue;
            }
            let [a, b, c] = self.tri_ip(ti as u32);
            if orient(a, b, c) <= 0 {
                return Err(format!("tri {} not ccw", ti));
            }
            for i in 0..3 {
                let e = self.edges[t.e[i] as usize];
                if !e.alive {
                    return Err(format!("tri {} uses dead edge", ti));
                }
                let (u, w) = (t.v[i], t.v[nx(i)]);
                if !((e.a == u.min(w)) && (e.b == u.max(w))) {
                    return Err(format!("tri {} edge {} endpoints mismatch", ti, i));
                }
                let nb = t.n[i];
                if nb != NONE {
                    let o = &self.tris[nb as usize];
                    if !o.alive {
                        return Err(format!("tri {} neighbour dead", ti));
                    }
                    let k = o.n.iter().position(|&x| x == ti as u32).ok_or(format!("tri {} not back-linked", ti))?;
                    if o.v[k] != w || o.v[nx(k)] != u || o.e[k] != t.e[i] {
                        return Err(format!("tri {} neighbour edge mismatch", ti));
                    }
                }
            }
        }
        Ok(())
    }

    /// Delaunay check ignoring fixed edges (tests).
    pub fn is_delaunay(&self) -> bool {
        for (ti, t) in self.tris.iter().enumerate() {
            if !t.alive {
                continue;
            }
            let [a, b, c] = self.tri_ip(ti as u32);
            for i in 0..3 {
                if self.edges[t.e[i] as usize].fixed || t.n[i] == NONE {
                    continue;
                }
                let o = &self.tris[t.n[i] as usize];
                let k = o.n.iter().position(|&x| x == ti as u32).unwrap();
                let q = self.pts[o.v[pv(k)] as usize];
                if incircle(a, b, c, q) > 0 {
                    return false;
                }
            }
        }
        true
    }
}

/// Triangulate the pseudo-polygon with base a->b and chain (ordered from a's side to b's side) on its left.
fn tri_pseudo(pts: &[IP], a: u32, b: u32, chain: &[u32], out: &mut Vec<[u32; 3]>) {
    if chain.is_empty() {
        return;
    }
    let (pa, pb) = (pts[a as usize], pts[b as usize]);
    let mut ci = 0;
    for i in 1..chain.len() {
        if incircle(pa, pb, pts[chain[ci] as usize], pts[chain[i] as usize]) > 0 {
            ci = i;
        }
    }
    let c = chain[ci];
    tri_pseudo(pts, a, c, &chain[..ci], out);
    tri_pseudo(pts, c, b, &chain[ci + 1..], out);
    out.push([a, b, c]);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rng(seed: &mut u64) -> u64 {
        *seed ^= *seed << 13;
        *seed ^= *seed >> 7;
        *seed ^= *seed << 17;
        *seed
    }

    #[test]
    fn random_insert_remove_constraints() {
        let mut c = Cdt::new([-1000, -1000], [101_000, 101_000]);
        let mut s = 12345u64;
        let mut vs = Vec::new();
        for _ in 0..400 {
            let x = (rng(&mut s) % 100_000) as i64;
            let y = (rng(&mut s) % 100_000) as i64;
            vs.push(c.insert([x, y]).unwrap());
        }
        c.validate().unwrap();
        assert!(c.is_delaunay());
        // a grid of collinear points and constraints through them
        let mut g = Vec::new();
        for i in 0..10 {
            g.push(c.insert([5000 + i * 3000, 50_500]).unwrap());
        }
        c.insert_fixed(g[0], g[9]).unwrap();
        c.validate().unwrap();
        for i in 0..9 {
            assert!(c.edges[c.find_edge(g[i], g[i + 1]).unwrap() as usize].fixed);
        }
        // journal: insert then roll back
        let before = c.tris.iter().filter(|t| t.alive).count();
        c.begin();
        for _ in 0..30 {
            let x = (rng(&mut s) % 100_000) as i64;
            let y = (rng(&mut s) % 100_000) as i64;
            let _ = c.insert([x, y]);
        }
        c.rollback();
        c.validate().unwrap();
        assert_eq!(before, c.tris.iter().filter(|t| t.alive).count());
        // remove random vertices
        for k in 0..100 {
            let v = vs[k * 3];
            let star = c.star(v);
            if star.is_empty() {
                continue;
            }
            c.apply_remove(v, &star).unwrap();
        }
        c.validate().unwrap();
    }
}
