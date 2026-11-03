//! Routing control: plane connections, net-by-net tree growth with transactional commits, and negotiated
//! rip-up and reroute for connections the plain search cannot place.

use crate::cdt::NONE;
use crate::geom::*;
use crate::model::Kind;
use crate::search::{End, PathOut, Query, Search};
use crate::world::{Term, Viol, World};
use crate::fx::FxMap as HashMap;
use std::time::Instant;

#[derive(Clone)]
pub struct Cfg {
    /// One routing order only (no parallel fan-out-first attempt): faster, sometimes less complete.
    pub single_attempt: bool,
    pub only: Option<String>,
    pub via_cost: f64,
    pub max_exp: usize,
    pub seconds: f64,
    pub rip_rounds: usize,
    pub verbose: bool,
}

/// What happened to one differential pair, for the design side to act on (move a part, free a corridor).
#[derive(Clone, Default, Debug)]
pub struct PairDiag {
    /// Stations in path order: "U1202.8/U1202.9 @ (89.25,46.94)".
    pub stations: Vec<String>,
    /// Per segment: (from station, to station, how it was routed and why).
    pub segments: Vec<(String, String, String)>,
    /// Coupled fraction right after the pair phase (before other nets were routed).
    pub coupled_after_pairs: f64,
    /// Rip-up disturbances: which net moved the pair, where, and how the pair was rebuilt.
    pub moved: Vec<String>,
    pub notes: Vec<String>,
}

#[derive(Clone)]
pub struct Router {
    pub w: World,
    pub cfg: Cfg,
    pub net_pads: Vec<Vec<u32>>,
    pub history: Vec<f64>,
    pub log: Vec<String>,
    pub t0: Instant,
    pub stats: Stats,
    pub dumps: usize,
    /// Location of the violation that rejected the last commit.
    pub last_viol: Option<Viol>,
    pub last_exp: usize,
    /// Objects joined without a piece (a via inside its pad).
    pub ties: Vec<(u32, u32)>,
    /// Wire-end objects of unzipped pair members (virtual, never exported).
    pub junctions: Vec<u32>,
    /// Unzip mode: members leave the ribbon to reach stations with spread-out pads by stubs.
    pub split_wide: bool,
    /// Other nets may not rip differential pair members (relaxed for the last rip-up rounds).
    pub protect_pairs: bool,
    /// Pad-to-pad span of each pair member (0 for other nets), for the cost of ripping it.
    pub pair_len: Vec<f64>,
    /// Last rip-up phase: any net may move pairs at the plain cost.
    pub relax_pairs: bool,
    /// Per pair (first member net): diagnostics reported in routing.json.
    pub pair_diag: HashMap<u32, PairDiag>,
    /// Flow events worth reporting (rip-up phases, pair rescues and rebuilds).
    pub events: Vec<String>,
    /// Current rip-up round (-1: first pass / pair phase).
    pub round: i32,
    /// Search budget multiplier: searches fail fast while many nets compete, and dig deeper once
    /// rip-up stalls on the few that are left.
    pub budget_scale: f64,
    /// Set when a parallel attempt is called off.
    pub stop: Option<std::sync::Arc<std::sync::atomic::AtomicBool>>,
    /// Negotiated congestion history per (layer, triangle).
    pub tri_hist: std::sync::Arc<HashMap<(u8, u32), f64>>,
}

#[derive(Default, Debug, Clone)]
pub struct Stats {
    pub searches: usize,
    pub search_fail: usize,
    pub commits: usize,
    pub commit_fail: usize,
    pub drc_reject: usize,
    pub rips: usize,
    pub expansions: usize,
    pub search_s: f64,
    pub commit_s: f64,
    /// Expansions and time spent in searches that found nothing; time in soft (rip-up) searches.
    pub fail_exp: usize,
    pub fail_s: f64,
    pub soft_s: f64,
}

/// A segment failure in words for the design side.
fn plain_reason(e: &str) -> String {
    if let Some(rest) = e.strip_prefix("unzip failed: ") {
        if let Some(i) = rest.find(" at (") {
            return format!("the two members collide with copper at {} where the ribbon is split", &rest[i + 4..]);
        }
        return format!("splitting the ribbon into its members failed ({})", rest);
    }
    if e.starts_with("ribbon route failed") {
        return "no corridor wide enough for the coupled pair".into();
    }
    if e.starts_with("no room for a ribbon junction") {
        return format!("no room next to the pads to bring the members together{}", e.trim_start_matches("no room for a ribbon junction"));
    }
    if e.starts_with("stub ") {
        return format!("a member could not reach its pad from the pair ({})", e);
    }
    e.to_string()
}

fn tr_centroid(w: &World, l: usize, t: u32) -> P2 {
    let tr = w.layers[l].cdt.tris[t as usize];
    (w.vpos(l, tr.v[0]) + w.vpos(l, tr.v[1]) + w.vpos(l, tr.v[2])) * (1.0 / 3.0)
}

fn find(p: &mut Vec<usize>, x: usize) -> usize {
    let mut r = x;
    while p[r] != r {
        r = p[r];
    }
    let mut y = x;
    while p[y] != r {
        let n = p[y];
        p[y] = r;
        y = n;
    }
    r
}

impl Router {
    pub fn new(w: World, cfg: Cfg) -> Router {
        let mut net_pads = vec![Vec::new(); w.b.nets.len()];
        for (oi, o) in w.objs.iter().enumerate() {
            if o.s.kind == Kind::Pad && o.s.net != NONE {
                net_pads[o.s.net as usize].push(oi as u32);
            }
        }
        // a zone's pieces join the pads and pieces of its net they touch (the fill's own connections)
        let mut ties = Vec::new();
        for (zi, z) in w.objs.iter().enumerate() {
            if !z.s.name.starts_with("zone:") || z.s.net == NONE {
                continue;
            }
            for &oi in &net_pads[z.s.net as usize] {
                let o = &w.objs[oi as usize];
                let once = !o.s.name.starts_with("zone:") || oi as usize > zi;
                if once && o.s.layers & z.s.layers != 0 && crate::world::core_dist(&z.s, &o.s) - o.s.r < 0.03 {
                    ties.push((zi as u32, oi));
                }
            }
        }
        let n = w.b.nets.len();
        Router { w, cfg, net_pads, history: vec![0.0; n], log: Vec::new(), t0: Instant::now(), stats: Stats::default(), dumps: 0, last_viol: None, last_exp: 0, ties, junctions: Vec::new(), split_wide: true, protect_pairs: false, pair_len: Vec::new(), relax_pairs: false, pair_diag: HashMap::default(), events: Vec::new(), round: -1, budget_scale: 1.0, stop: None, tri_hist: std::sync::Arc::new(HashMap::default()) }
    }

    pub fn elapsed(&self) -> f64 {
        self.t0.elapsed().as_secs_f64()
    }

    pub fn is_plane(&self, net: u32) -> bool {
        self.w.b.nets[net as usize].plane.is_some()
    }

    /// Objects of a net (pads and live vias).
    pub fn net_objs(&self, net: u32) -> Vec<u32> {
        let mut v = self.net_pads[net as usize].clone();
        for (oi, o) in self.w.objs.iter().enumerate() {
            if o.alive && o.s.kind == Kind::Via && o.s.net == net {
                v.push(oi as u32);
            }
        }
        for &j in &self.junctions {
            let o = &self.w.objs[j as usize];
            if o.alive && o.s.net == net {
                v.push(j);
            }
        }
        v
    }

    /// Connected components of a net's objects through its pieces; plane objects join through the plane.
    pub fn components(&self, net: u32) -> Vec<Vec<u32>> {
        let objs = self.net_objs(net);
        let idx = |o: u32| objs.iter().position(|&x| x == o);
        let mut par: Vec<usize> = (0..objs.len()).collect();
        let plane = self.is_plane(net);
        let mut plane_root: Option<usize> = None;
        for (i, &o) in objs.iter().enumerate() {
            for &p in &self.w.obj_pieces[o as usize] {
                let pc = &self.w.pieces[p as usize];
                if !pc.alive {
                    continue;
                }
                for other in [pc.a.obj, pc.b.obj] {
                    if let Some(j) = idx(other) {
                        let (a, b) = (find(&mut par, i), find(&mut par, j));
                        if a != b {
                            par[b] = a;
                        }
                    }
                }
            }
            for &(x, y) in &self.ties {
                if x == o || y == o {
                    let other = if x == o { y } else { x };
                    if !self.w.objs[other as usize].alive {
                        continue;
                    }
                    if let Some(j) = idx(other) {
                        let (a, b) = (find(&mut par, i), find(&mut par, j));
                        if a != b {
                            par[b] = a;
                        }
                    }
                }
            }
            if plane {
                let ob = &self.w.objs[o as usize];
                if ob.plane {
                    match plane_root {
                        None => plane_root = Some(i),
                        Some(r) => {
                            let (a, b) = (find(&mut par, r), find(&mut par, i));
                            if a != b {
                                par[b] = a;
                            }
                        }
                    }
                }
            }
        }
        let mut groups: Vec<(usize, Vec<u32>)> = Vec::new();
        for (i, &o) in objs.iter().enumerate() {
            let r = find(&mut par, i);
            match groups.iter_mut().find(|g| g.0 == r) {
                Some(g) => g.1.push(o),
                None => groups.push((r, vec![o])),
            }
        }
        // drop groups made only of vias (dangling vias are not terminals)
        let mut out: Vec<Vec<u32>> = groups.into_iter().map(|g| g.1).filter(|g| g.iter().any(|&o| self.w.objs[o as usize].s.kind == Kind::Pad) || plane && g.iter().any(|&o| self.w.objs[o as usize].s.kind == Kind::Via)).collect();
        out.sort_by_key(|g| std::cmp::Reverse(g.len()));
        out
    }

    pub fn plane_ok(&self, net: u32, comp: &[u32]) -> bool {
        comp.iter().any(|&o| self.w.objs[o as usize].plane) && self.is_plane(net)
    }

    /// Unconnected (pad component) count for a net: components - 1, or plane-less pad groups for plane nets.
    pub fn deficit(&self, net: u32) -> usize {
        let comps = self.components(net);
        if self.is_plane(net) {
            comps.iter().filter(|c| !self.plane_ok(net, c)).count()
        } else {
            comps.len().saturating_sub(1)
        }
    }

    fn query(&self, net: u32, sources: Vec<u32>, targets: Vec<u32>, narrow: bool) -> Query {
        let ni = &self.w.b.nets[net as usize];
        let nl = self.w.nsig();
        let mut cross_cost = vec![0.0; self.w.b.nets.len()];
        let plane_net = ni.plane.is_some();
        for (i, h) in self.history.iter().enumerate() {
            let paired = self.w.b.partner.get(i).map_or(false, |&x| x != NONE);
            // differential pairs are placed first and rebuilt only as pairs: other nets go round them
            // (a plane drop can always go elsewhere; a long pair costs more to move than a short one)
            let pair_cost = if self.protect_pairs {
                f64::INFINITY
            } else if self.relax_pairs {
                40.0
            } else if plane_net {
                f64::INFINITY
            } else {
                40.0 + 0.5 * self.pair_len.get(i).cloned().unwrap_or(0.0)
            };
            cross_cost[i] = if self.w.b.nets[i].plane.is_some() { 2.5 } else if paired { pair_cost } else { 6.0 } + h;
        }
        Query {
            net,
            width: if narrow { ni.wmin.clone() } else { ni.wpref.clone() },
            allowed: ni.allowed & ((1 << nl) - 1),
            sources,
            targets,
            plane_goal: false,
            via_cost: self.cfg.via_cost,
            soft: false,
            cross_cost,
            max_exp: (self.cfg.max_exp as f64 * self.budget_scale) as usize,
            layer_cost: vec![1.0; nl],
            max_len: f64::INFINITY,
            ignore: vec![],
            taboo: HashMap::default(),
            via_taboo: vec![],
            history: self.tri_hist.clone(),
            wedge_off: std::env::var("TOPOR_WEDGE_OFF").is_ok(),
            couple: NONE,
            couple_factor: 0.3,
            ribbon: vec![],
        }
    }

    pub fn search(&mut self, q: &Query) -> Option<PathOut> {
        let t = Instant::now();
        let mut s = Search::new(&self.w, q);
        let r = s.run();
        self.last_exp = s.expansions;
        self.stats.searches += 1;
        self.stats.expansions += s.expansions;
        self.stats.search_s += t.elapsed().as_secs_f64();
        if std::env::var("TOPOR_EXPHIST").is_ok() {
            eprintln!("exphist {} {} {} {} {}", if r.is_some() { "ok" } else { "fail" }, if q.soft { "soft" } else { "hard" }, s.expansions, q.max_exp, self.round);
        }
        if r.is_none() && std::env::var("TOPOR_SEARCHDBG").is_ok() {
            let names = |v: &[u32]| v.iter().take(6).map(|&o| self.w.objs[o as usize].s.name.clone()).collect::<Vec<_>>().join(",");
            eprintln!("search fail {}: {} expansions (cap {}), {} -> {}", self.w.b.nets[q.net as usize].name, s.expansions, q.max_exp, names(&q.sources), names(&q.targets));
        }
        if r.is_none() {
            self.stats.search_fail += 1;
            self.stats.fail_exp += s.expansions;
            self.stats.fail_s += t.elapsed().as_secs_f64();
        }
        if q.soft {
            self.stats.soft_s += t.elapsed().as_secs_f64();
        }
        r
    }

    /// Commit a found path transactionally. Returns false (state unchanged) if it does not realise cleanly.
    pub fn commit(&mut self, net: u32, path: &PathOut, width: &[f64]) -> Result<Vec<u32>, String> {
        self.last_viol = None;
        let t = Instant::now();
        let r = self.commit_inner(net, path, width);
        self.stats.commit_s += t.elapsed().as_secs_f64();
        if std::env::var("TOPOR_VALIDATE").is_ok() {
            if let Err(e) = self.w.validate_topology() {
                panic!("topology broken after commit of {} ({:?}): {}", self.w.b.nets[net as usize].name, r.as_ref().map(|_| ()), e);
            }
        }
        match &r {
            Ok(_) => self.stats.commits += 1,
            Err(_) => self.stats.commit_fail += 1,
        }
        r
    }

    fn commit_inner(&mut self, net: u32, path: &PathOut, width: &[f64]) -> Result<Vec<u32>, String> {
        self.w.begin();
        let mut vias = Vec::new();
        for &p in &path.vias {
            match self.w.add_via(net, p) {
                Ok(o) => vias.push(o),
                Err(e) => {
                    self.w.rollback();
                    return Err(format!("via: {}", e));
                }
            }
        }
        let mut pieces = Vec::new();
        for s in &path.segs {
            let l = s.layer;
            let term = |e: End| -> Term {
                match e {
                    End::Obj(t) => t,
                    End::NewVia(i) => Term { obj: vias[i], vtx: self.w.objs[vias[i] as usize].lv[l] },
                }
            };
            let (a, b) = (term(s.a), term(s.b));
            // crossings on edges removed by the via insertions: allowed only next to a via end
            let alive: Vec<bool> = s.cx.iter().map(|c| self.w.layers[l].cdt.edges[c.0 as usize].alive).collect();
            let mut lo = 0;
            let mut hi = s.cx.len();
            if matches!(s.a, End::NewVia(_)) {
                while lo < hi && !alive[lo] {
                    lo += 1;
                }
            }
            if matches!(s.b, End::NewVia(_)) {
                while hi > lo && !alive[hi - 1] {
                    hi -= 1;
                }
            }
            let cx: Vec<(u32, usize, u32)> = s.cx[lo..hi].to_vec();
            // crossings must still chain through shared triangles after the via insertions
            let lay = &self.w.layers[l];
            let tri_of = |e: u32, left: u32, entering: bool| -> u32 {
                let ed = lay.cdt.edges[e as usize];
                let from_a = (left == ed.a) == entering;
                lay.etri[e as usize][if from_a { 0 } else { 1 }]
            };
            let has_v = |t: u32, v: u32| t != NONE && lay.cdt.tris[t as usize].v.contains(&v);
            let mut chained = alive[lo..hi].iter().all(|&x| x);
            if chained {
                for k in 0..cx.len().saturating_sub(1) {
                    if tri_of(cx[k].0, cx[k].2, true) != tri_of(cx[k + 1].0, cx[k + 1].2, false) {
                        chained = false;
                    }
                }
                if let (Some(f), true) = (cx.first(), a.vtx != NONE) {
                    if !has_v(tri_of(f.0, f.2, false), a.vtx) {
                        chained = false;
                    }
                }
                if let (Some(f), true) = (cx.last(), b.vtx != NONE) {
                    if !has_v(tri_of(f.0, f.2, true), b.vtx) {
                        chained = false;
                    }
                }
                if cx.is_empty() && a.vtx != NONE && b.vtx != NONE && lay.cdt.find_edge(a.vtx, b.vtx).is_none() {
                    chained = false;
                }
            }
            if !chained {
                // re-search this segment between its fixed ends on its layer
                let mut q = self.query(net, vec![a.obj], vec![b.obj], false);
                q.width = width.to_vec();
                q.allowed = 1 << l;
                let mut se = Search::new(&self.w, &q);
                match se.run() {
                    Some(p2) if p2.segs.len() == 1 => {
                        let s2 = &p2.segs[0];
                        let a2 = match s2.a {
                            End::Obj(t) => t,
                            _ => a,
                        };
                        let b2 = match s2.b {
                            End::Obj(t) => t,
                            _ => b,
                        };
                        let pid = self.w.add_piece(net, l, width[l], a2, b2, &s2.cx);
                        pieces.push(pid);
                        continue;
                    }
                    _ => {
                        self.w.rollback();
                        return Err("segment crosses a via cavity".into());
                    }
                }
            }
            let pid = self.w.add_piece(net, l, width[l], a, b, &cx);
            pieces.push(pid);
        }
        let realized = self.w.realize_dirty();
        let mut viol: Vec<Viol> = Vec::new();
        for &p in &realized {
            if self.w.pieces[p as usize].overlap {
                let at = self.w.pieces[p as usize].problem.unwrap_or(self.w.pieces[p as usize].poly[0]);
                viol.push(Viol { what: "overflow".into(), layer: self.w.pieces[p as usize].layer as i32, at, dist: 0., need: 0., pieces: vec![p], objs: vec![] });
            }
            self.w.check_piece(p, &mut viol);
            if !viol.is_empty() {
                break;
            }
        }
        if viol.is_empty() {
            for &o in &vias {
                let c = self.w.objs[o as usize].s.center;
                self.w.check_via_at(net, c, o, &mut viol, true);
            }
        }
        if !viol.is_empty() {
            let v = &viol[0];
            let msg = format!("{} on layer {} at ({:.3},{:.3}) d={:.4} need {:.4}", v.what, v.layer, v.at.x, v.at.y, v.dist, v.need);
            let want = std::env::var("TOPOR_DUMP_KIND").map(|k| v.what.contains(&k)).unwrap_or(true);
            if let (Ok(dir), true) = (std::env::var("TOPOR_DUMP"), want) {
                self.dumps += 1;
                if self.dumps <= 6 {
                    let b = [v.at.x - 2.5, v.at.y - 2.0, v.at.x + 2.5, v.at.y + 2.0];
                    let svg = crate::svg::render(&self.w, &crate::svg::Opts { mesh: true, box_: Some(b), scale: 70.0, marks: vec![(v.at, format!("{} d={:.3}", v.what, v.dist))], unrouted: vec![] });
                    let _ = std::fs::write(format!("{}/dump{}.svg", dir, self.dumps), svg);
                    eprintln!("dump{}: {} net {} pieces {:?} objs {:?}", self.dumps, msg, self.w.b.nets[net as usize].name, v.pieces, v.objs.iter().map(|&o| self.w.objs[o as usize].s.name.clone()).collect::<Vec<_>>());
                    if let Some(&pp) = v.pieces.first() {
                        let pc = &self.w.pieces[pp as usize];
                        let l = pc.layer as usize;
                        eprintln!("   piece {} {} -> {} width {}", pp, self.w.objs[pc.a.obj as usize].s.name, self.w.objs[pc.b.obj as usize].s.name, pc.width);
                        for (i, (a, b)) in self.w.portals(pp).iter().enumerate() {
                            let on = |id: u32| if id < 0x8000_0000 { let o = self.w.layers[l].vobj[id as usize]; if o == NONE { "?".to_string() } else { self.w.objs[o as usize].s.name.clone() } } else { format!("x{:x}", id) };
                            eprintln!("     portal {} L {} ({:.3},{:.3}) r {:.4} | R {} ({:.3},{:.3}) r {:.4}", i, on(a.id), a.c.x, a.c.y, a.r, on(b.id), b.c.x, b.c.y, b.r);
                        }
                        eprintln!("     poly {:?}", pc.poly.iter().map(|q| format!("({:.3},{:.3})", q.x, q.y)).collect::<Vec<_>>());
                        for &c in &pc.cx {
                            let cr = &self.w.crosses[c as usize];
                            let ed = self.w.layers[l].cdt.edges[cr.edge as usize];
                            let (ea, eb) = (self.w.vpos(l, ed.a), self.w.vpos(l, ed.b));
                            let on = |v: u32| { let o = self.w.layers[l].vobj[v as usize]; if o == NONE { "?".to_string() } else { self.w.objs[o as usize].s.name.clone() } };
                            eprintln!("     cross edge {} {}({:.3},{:.3}) - {}({:.3},{:.3}) left {} list {:?}", cr.edge, on(ed.a), ea.x, ea.y, on(ed.b), eb.x, eb.y, on(cr.left), self.w.layers[l].cross[cr.edge as usize].iter().map(|&x| self.w.b.nets[self.w.pieces[self.w.crosses[x as usize].piece as usize].net as usize].name.clone()).collect::<Vec<_>>());
                        }
                        // the new pieces of this commit
                        for &np in &pieces {
                            let q = &self.w.pieces[np as usize];
                            eprintln!("     new piece {} {} -> {} cx {} poly {:?}", np, self.w.objs[q.a.obj as usize].s.name, self.w.objs[q.b.obj as usize].s.name, q.cx.len(), q.poly.iter().map(|z| format!("({:.3},{:.3})", z.x, z.y)).collect::<Vec<_>>());
                        }
                    }
                }
            }
            self.last_viol = Some(v.clone());
            self.w.rollback();
            self.stats.drc_reject += 1;
            return Err(msg);
        }
        self.w.commit();
        Ok(pieces)
    }

    /// Penalise the neighbourhood of the last rejected commit's violation for the next attempt.
    fn learn(&self, q: &mut Query) {
        let v = match &self.last_viol {
            Some(v) => v.clone(),
            None => return,
        };
        if v.layer < 0 || v.what.starts_with("hole") || v.what.starts_with("via") {
            q.via_taboo.push((v.at, 0.6));
        }
        let layers: Vec<usize> = if v.layer >= 0 { vec![v.layer as usize] } else { (0..self.w.nsig()).collect() };
        for l in layers {
            let lay = &self.w.layers[l];
            if let Ok(loc) = self.w.locate(l, v.at) {
                let t0 = match loc {
                    crate::cdt::Loc::Inside(t) | crate::cdt::Loc::OnEdge(t, _) => t,
                    crate::cdt::Loc::OnVertex(vx) => lay.cdt.vtri[vx as usize],
                };
                let mut stack = vec![t0];
                let mut seen = vec![t0];
                while let Some(t) = stack.pop() {
                    *q.taboo.entry((l as u8, t)).or_insert(0.0) += 4.0;
                    let tr = lay.cdt.tris[t as usize];
                    for &nb in &tr.n {
                        if nb == NONE || seen.contains(&nb) || !lay.cdt.tris[nb as usize].alive {
                            continue;
                        }
                        let c = tr_centroid(&self.w, l, nb);
                        if c.dist(v.at) < 0.8 && seen.len() < 40 {
                            seen.push(nb);
                            stack.push(nb);
                        }
                    }
                }
            }
        }
    }

    /// Remove a piece and any via left without pieces.
    pub fn rip_piece(&mut self, p: u32) {
        if !self.w.pieces[p as usize].alive {
            return;
        }
        let (a, b) = (self.w.pieces[p as usize].a.obj, self.w.pieces[p as usize].b.obj);
        if std::env::var("TOPOR_PAIRDBG").is_ok() && self.w.b.partner.get(self.w.pieces[p as usize].net as usize).map_or(false, |&x| x != NONE) {
            eprintln!("rip pair piece {} of {}: {}", p, self.w.b.nets[self.w.pieces[p as usize].net as usize].name, std::backtrace::Backtrace::force_capture().to_string().lines().filter(|l| l.contains("topor::router")).take(3).collect::<Vec<_>>().join(" <- "));
        }
        self.w.remove_piece(p);
        if std::env::var("TOPOR_VALIDATE").is_ok() {
            if let Err(e) = self.w.validate_topology() {
                panic!("topology broken after removing piece {}: {}", p, e);
            }
        }
        for o in [a, b] {
            let ob = &self.w.objs[o as usize];
            if ob.alive && ob.s.kind == Kind::Via && self.w.obj_pieces[o as usize].iter().all(|&q| !self.w.pieces[q as usize].alive) {
                if std::env::var("TOPOR_VALIDATE").is_ok() {
                    if let Err(e) = self.w.validate_topology() {
                        panic!("topology broken before removing via {}: {}", o, e);
                    }
                    eprintln!("removing via {} at {:?} pieces {:?}", o, self.w.objs[o as usize].s.center, self.w.obj_pieces[o as usize]);
                }
                if let Err(e) = self.w.remove_via(o) {
                    self.log.push(format!("via removal failed: {}", e));
                }
                if std::env::var("TOPOR_VALIDATE").is_ok() {
                    if let Err(e) = self.w.validate_topology() {
                        panic!("topology broken after removing via {}: {}", o, e);
                    }
                }
            }
        }
        self.stats.rips += 1;
    }

    /// Vias inside a large pad of the group (exposed pads): no track needed. A pad big enough for several
    /// gets an array (about 1 mm pitch) for low inductance and heat; a small one a single via.
    fn via_in_pad(&mut self, net: u32, comp: &[u32]) -> bool {
        let (vd, _) = self.w.b.via_size(net);
        let mut placed = 0;
        for &o in comp {
            let ob = &self.w.objs[o as usize];
            if ob.s.kind != Kind::Pad || ob.s.pts.len() < 3 {
                continue;
            }
            let (lo, hi) = ob.s.pts.iter().fold((p2(1e9, 1e9), p2(-1e9, -1e9)), |(lo, hi), q| (p2(lo.x.min(q.x), lo.y.min(q.y)), p2(hi.x.max(q.x), hi.y.max(q.y))));
            let c = ob.s.center;
            // grid over the pad, inset by the via and a little copper margin
            let pitch = (vd + 0.5).max(1.0);
            let span = |a: f64, b: f64| -> Vec<f64> {
                let free = (b - a) - vd - 0.2;
                let n = if free > 0.0 { (free / pitch).floor() as usize + 1 } else { 1 };
                let n = n.min(4);
                let mid = (a + b) / 2.0;
                (0..n).map(|i| mid + (i as f64 - (n - 1) as f64 / 2.0) * pitch).collect()
            };
            let mut cands: Vec<P2> = Vec::new();
            for &y in &span(lo.y, hi.y) {
                for &x in &span(lo.x, hi.x) {
                    cands.push(p2(x, y));
                }
            }
            // centre-most first, so a crowded pad still gets its best via
            cands.sort_by(|a, b| a.dist(c).partial_cmp(&b.dist(c)).unwrap());
            if cands.len() == 1 {
                for &(dx, dy) in &[(0.6, 0.0), (-0.6, 0.0), (0.0, 0.6), (0.0, -0.6), (0.6, 0.6), (-0.6, -0.6), (0.6, -0.6), (-0.6, 0.6)] {
                    cands.push(c + p2(dx, dy));
                }
            }
            let array = cands.len() > 1 && span(lo.x, hi.x).len() * span(lo.y, hi.y).len() > 1;
            for x in cands {
                if !self.w.via_fits_in_pad(net, x, o) || !self.w.b.in_plane(net, x, vd / 2.0) {
                    continue;
                }
                self.w.begin();
                match self.w.add_via(net, x) {
                    Ok(v) => {
                        let mut viol = Vec::new();
                        self.w.check_via_at(net, x, v, &mut viol, true);
                        if viol.is_empty() {
                            self.w.commit();
                            self.ties.push((o, v));
                            placed += 1;
                            if !array {
                                return true;
                            }
                        } else {
                            self.w.rollback();
                        }
                    }
                    Err(_) => self.w.rollback(),
                }
            }
        }
        placed > 0
    }

    /// Smaller side of a pad's copper box (0 for vias and point pads).
    fn pad_extent(&self, o: u32) -> f64 {
        let s = &self.w.objs[o as usize].s;
        if s.kind != Kind::Pad || s.pts.len() < 3 {
            return 0.0;
        }
        let (lo, hi) = s.pts.iter().fold((p2(1e9, 1e9), p2(-1e9, -1e9)), |(lo, hi), q| (p2(lo.x.min(q.x), lo.y.min(q.y)), p2(hi.x.max(q.x), hi.y.max(q.y))));
        (hi.x - lo.x).min(hi.y - lo.y) + 2.0 * s.r
    }

    /// Connect one plane-net pad group to its plane (an existing plane object nearby or a new via).
    fn route_plane_comp(&mut self, net: u32, comp: &[u32]) -> bool {
        if comp.iter().any(|&o| {
            let s = &self.w.objs[o as usize].s;
            s.pts.len() >= 3 && {
                let (lo, hi) = s.pts.iter().fold((p2(1e9, 1e9), p2(-1e9, -1e9)), |(lo, hi), q| (p2(lo.x.min(q.x), lo.y.min(q.y)), p2(hi.x.max(q.x), hi.y.max(q.y))));
                hi.x - lo.x + 2.0 * s.r >= 1.6 && hi.y - lo.y + 2.0 * s.r >= 1.6
            }
        }) && self.via_in_pad(net, comp)
        {
            return true;
        }
        let center = self.w.objs[comp[0] as usize].s.center;
        let mut targets = Vec::new();
        for o in self.net_objs(net) {
            if comp.contains(&o) {
                continue;
            }
            let ob = &self.w.objs[o as usize];
            if ob.plane && ob.s.center.dist(center) < 2.0 {
                targets.push(o);
            }
        }
        // or tie to a neighbouring pad of the net that already has its drop
        let comps = self.components(net);
        for c in &comps {
            if c.as_slice() == comp || !self.plane_ok(net, c) {
                continue;
            }
            for &o in c {
                let ob = &self.w.objs[o as usize];
                // by distance to the pad's copper: a chip's exposed pad is big, its centre far away
                if ob.s.kind == Kind::Pad && point_hull_dist(&ob.s.pts, center) - ob.s.r < 2.0 && !targets.contains(&o) {
                    targets.push(o);
                }
            }
        }
        // a short tie to a neighbouring pad that is already on the plane (a chip's exposed pad with its
        // vias) beats a via of its own
        let ties: Vec<u32> = targets.iter().cloned().filter(|&o| self.w.objs[o as usize].s.kind == Kind::Pad).collect();
        if !ties.is_empty() {
            let mut q = self.query(net, comp.to_vec(), ties, false);
            q.max_len = 1.5;
            let mut layers = 0u32;
            for &o in comp {
                layers |= self.w.objs[o as usize].s.layers;
            }
            q.allowed &= layers;
            q.via_cost = 1e6;
            let w = q.width.clone();
            for _attempt in 0..2 {
                let path = match self.search(&q) {
                    Some(p) => p,
                    None => break,
                };
                if path.vias.is_empty() && self.commit(net, &path, &w).is_ok() {
                    return true;
                }
                self.learn(&mut q);
            }
        }
        for narrow in [false, true] {
            let mut q = self.query(net, comp.to_vec(), targets.clone(), narrow);
            q.plane_goal = true;
            q.max_len = 4.0;
            // plane drops stay on the pad's own layer
            let mut layers = 0u32;
            for &o in comp {
                layers |= self.w.objs[o as usize].s.layers;
            }
            q.allowed &= layers;
            let w = q.width.clone();
            for _attempt in 0..4 {
                let path = match self.search(&q) {
                    Some(p) => p,
                    None => break,
                };
                match self.commit(net, &path, &w) {
                    Ok(_) => return true,
                    Err(e) => {
                        if self.cfg.verbose {
                            self.log.push(format!("plane {} {}: {}", self.w.b.nets[net as usize].name, self.w.objs[comp[0] as usize].s.name, e));
                        }
                        self.learn(&mut q);
                    }
                }
            }
        }
        // no pour of the net within reach (a partial plane): wire it to a part of the net that is on it
        if !self.w.b.plane_regions.is_empty() {
            let reach: Vec<u32> = self.components(net).iter().filter(|c| c.as_slice() != comp && self.plane_ok(net, c)).flatten().cloned().collect();
            if !reach.is_empty() && self.connect(net, comp[0], reach[0], NONE) {
                return true;
            }
            if !reach.is_empty() {
                let mut q = self.query(net, comp.to_vec(), reach, false);
                q.max_exp = (q.max_exp / 2).max(5000);
                let w = q.width.clone();
                for _attempt in 0..3 {
                    let path = match self.search(&q) {
                        Some(p) => p,
                        None => break,
                    };
                    if self.commit(net, &path, &w).is_ok() {
                        return true;
                    }
                    self.learn(&mut q);
                }
            }
        }
        false
    }

    /// Grow a net's tree until connected or stuck. Returns remaining deficit.
    pub fn route_net(&mut self, net: u32) -> usize {
        if self.is_plane(net) {
            // groups with a big pad first (vias in the pad): small pads beside it then tie to it
            let mut comps = self.components(net);
            let big = |me: &Router, c: &Vec<u32>| -> f64 { c.iter().map(|&o| me.pad_extent(o)).fold(0.0, f64::max) };
            comps.sort_by(|a, b| big(self, b).partial_cmp(&big(self, a)).unwrap());
            for c in comps {
                if !self.plane_ok(net, &c) && self.components(net).contains(&c) {
                    self.route_plane_comp(net, &c);
                }
            }
            return self.deficit(net);
        }
        let mut stuck = 0;
        loop {
            let comps = self.components(net);
            if comps.len() <= 1 {
                return 0;
            }
            if stuck >= comps.len() {
                return comps.len() - 1;
            }
            let src = comps[stuck % comps.len()].clone();
            let targets: Vec<u32> = comps.iter().filter(|c| **c != src).flatten().cloned().collect();
            let mut done = false;
            'widths: for narrow in [false, true] {
                let mut q = self.query(net, src.clone(), targets.clone(), narrow);
                // successful searches rarely need more than a few thousand expansions; only the first
                // attempt gets the full budget, retries (narrow, other sources) a sixth of it
                if narrow || stuck > 0 {
                    q.max_exp = (q.max_exp / 6).max(5000);
                }
                let w = q.width.clone();
                for _attempt in 0..4 {
                    let path = match self.search(&q) {
                        Some(p) => p,
                        None => break,
                    };
                    match self.commit(net, &path, &w) {
                        Ok(_) => {
                            done = true;
                            break 'widths;
                        }
                        Err(e) => {
                            if self.cfg.verbose {
                                self.log.push(format!("{}: commit failed: {}", self.w.b.nets[net as usize].name, e));
                            }
                            self.learn(&mut q);
                        }
                    }
                }
            }
            if done {
                stuck = 0;
            } else {
                stuck += 1;
            }
        }
    }

    /// Negotiated repair of one net: find the cheapest path through other nets' wires, rip them, route.
    pub fn repair_net(&mut self, net: u32) -> Vec<u32> {
        let comps = self.components(net);
        let (src, targets): (Vec<u32>, Vec<u32>) = if self.is_plane(net) {
            match comps.iter().find(|c| !self.plane_ok(net, c)) {
                Some(c) => (c.clone(), vec![]),
                None => return vec![],
            }
        } else {
            if comps.len() <= 1 {
                return vec![];
            }
            (comps[0].clone(), comps[1..].iter().flatten().cloned().collect())
        };
        let mut q = self.query(net, src.clone(), targets.clone(), true);
        q.soft = true;
        q.max_exp = (self.cfg.max_exp as f64 * self.budget_scale) as usize * 7;
        if self.is_plane(net) {
            q.plane_goal = true;
            q.max_len = 5.0;
            let mut layers = 0u32;
            for &o in &src {
                layers |= self.w.objs[o as usize].s.layers;
            }
            q.allowed &= layers;
        }
        let dbg = std::env::var("TOPOR_REPAIR").is_ok();
        let mut ripped_nets: Vec<u32> = Vec::new();
        let mut path = None;
        for _attempt in 0..5 {
            let p = match self.search(&q) {
                Some(p) => p,
                None => {
                    if dbg {
                        eprintln!("repair {}: soft search failed (src {:?} targets {:?}) after {} expansions", self.w.b.nets[net as usize].name, src.iter().map(|&o| self.w.objs[o as usize].s.name.clone()).collect::<Vec<_>>(), targets.iter().map(|&o| self.w.objs[o as usize].s.name.clone()).collect::<Vec<_>>(), self.last_exp);
                    }
                    return ripped_nets;
                }
            };
            if dbg {
                eprintln!("repair {}: soft path cost {:.2}, {} victims {:?}", self.w.b.nets[net as usize].name, p.cost, p.victims.len(), p.victims.iter().map(|&x| self.w.b.nets[self.w.pieces[x as usize].net as usize].name.clone()).collect::<Vec<_>>());
            }

            if !p.cost.is_finite() {
                return ripped_nets;
            }
            if !p.victims.is_empty() {
                path = Some(p);
                break;
            }
            // the soft search found room the hard checks rejected: let the exact check decide, and rip
            // whatever other net's copper it names
            let w = q.width.clone();
            match self.commit(net, &p, &w) {
                Ok(_) => return ripped_nets,
                Err(e) => {
                    if dbg {
                        eprintln!("   victim-free soft commit: {}", e);
                    }
                    let v = self.last_viol.clone();
                    let mut culprits: Vec<u32> = Vec::new();
                    if let Some(v) = &v {
                        for &pp in &v.pieces {
                            if (pp as usize) < self.w.pieces.len() && self.w.pieces[pp as usize].alive && self.w.pieces[pp as usize].net != net {
                                culprits.push(pp);
                            }
                        }
                        for &o in &v.objs {
                            if (o as usize) < self.w.objs.len() && self.w.objs[o as usize].alive && self.w.objs[o as usize].s.kind == Kind::Via && self.w.objs[o as usize].s.net != net {
                                culprits.extend(self.w.obj_pieces[o as usize].iter().cloned().filter(|&x| self.w.pieces[x as usize].alive));
                            }
                        }
                    }
                    // differential pairs are rebuilt only as pairs: do not rip them for a single net
                    if self.protect_pairs {
                        culprits.retain(|&pp| self.w.b.partner.get(self.w.pieces[pp as usize].net as usize).map_or(true, |&x| x == NONE));
                    }

                    if culprits.is_empty() {
                        self.learn(&mut q);
                        continue;
                    }
                    for pp in culprits {
                        let pn = self.w.pieces[pp as usize].net;
                        self.note_pair_moved(pn, pp, net);
                        if !ripped_nets.contains(&pn) {
                            ripped_nets.push(pn);
                        }
                        self.history[pn as usize] += 1.0;
                        self.rip_piece(pp);
                    }
                    if let Some(v) = &v {
                        self.learn(&mut q);
                        let _ = v;
                    }
                }
            }
        }
        let path = match path {
            Some(p) => p,
            None => return ripped_nets,
        };
        let mut victim_nets = ripped_nets;
        {
            let h = std::sync::Arc::make_mut(&mut self.tri_hist);
            for &t in &path.cross_tris {
                *h.entry(t).or_insert(0.0) += 3.0;
            }
        }
        // the soft path itself, with its crossing slots shifted past the wires about to go: valid as
        // long as ripping removes no via (that would change the triangulation under it)
        let removes_via = path.victims.iter().any(|&p| {
            let pc = &self.w.pieces[p as usize];
            [pc.a.obj, pc.b.obj].iter().any(|&o| {
                let ob = &self.w.objs[o as usize];
                ob.s.kind == Kind::Via && self.w.obj_pieces[o as usize].iter().all(|&q| !self.w.pieces[q as usize].alive || path.victims.contains(&q))
            })
        });
        let mut direct = path.clone();
        if !removes_via {
            for sg in direct.segs.iter_mut() {
                let l = sg.layer;
                for c in sg.cx.iter_mut() {
                    let list = &self.w.layers[l].cross[c.0 as usize];
                    let before = list.iter().take(c.1).filter(|&&x| path.victims.contains(&self.w.crosses[x as usize].piece)).count();
                    c.1 -= before;
                }
            }
        }
        for &p in &path.victims {
            let pn = self.w.pieces[p as usize].net;
            self.note_pair_moved(pn, p, net);
            if !victim_nets.contains(&pn) {
                victim_nets.push(pn);
            }
            self.history[pn as usize] += 1.0;
            self.rip_piece(p);
        }
        self.history[net as usize] += 0.5;
        if !removes_via && std::env::var("TOPOR_DIRECT").is_ok() {
            let w = self.query(net, Vec::new(), Vec::new(), true).width;
            let _ = self.commit(net, &direct, &w);
        }
        victim_nets
    }

    /// Record that `by` displaced pair member `net` (piece pc) during rip-up.
    fn note_pair_moved(&mut self, net: u32, pc: u32, by: u32) {
        let pn = match self.w.b.partner.get(net as usize) {
            Some(&x) if x != NONE => x,
            _ => return,
        };
        let key = net.min(pn);
        let poly = &self.w.pieces[pc as usize].poly;
        let target: Vec<P2> = self.net_pads[by as usize].iter().map(|&o| self.w.objs[o as usize].s.center).collect();
        let at = poly.iter().cloned().min_by(|a, b| {
            let da = target.iter().map(|t| t.dist(*a)).fold(f64::INFINITY, f64::min);
            let db = target.iter().map(|t| t.dist(*b)).fold(f64::INFINITY, f64::min);
            da.partial_cmp(&db).unwrap()
        });
        let phase = if self.relax_pairs { "pairs released" } else { "pairs costly to move" };
        let msg = match at {
            Some(q) => format!("rip round {} ({}): {} needed the room at ({:.2},{:.2}) on {}", self.round, phase, self.w.b.nets[by as usize].name, q.x, q.y, self.w.b.sig_name(self.w.pieces[pc as usize].layer as usize)),
            None => format!("rip round {} ({}): moved by {}", self.round, phase, self.w.b.nets[by as usize].name),
        };
        let d = self.pair_diag.entry(key).or_default();
        if d.moved.last() != Some(&msg) {
            d.moved.push(msg);
        }
    }

    /// Up to `max` fitting junction positions for a station, best first, spread at least 0.25 apart.
    fn pair_junctions(&self, pp: u32, nn: u32, r: u32, toward: P2, max: usize) -> Vec<(P2, usize, f64)> {
        let mut out: Vec<(P2, usize, f64)> = Vec::new();
        let (op, on) = (&self.w.objs[pp as usize], &self.w.objs[nn as usize]);
        let ri = &self.w.b.nets[r as usize];
        let layers = op.s.layers & on.s.layers & ri.allowed;
        let l = match (0..self.w.nsig()).find(|&l| layers & (1 << l) != 0) {
            Some(l) => l,
            None => return out,
        };
        let (cp, cn) = (op.s.center, on.s.center);
        let a = (cn - cp).unit();
        let m = cp.lerp(cn, 0.5);
        let refname = op.s.name.split('.').next().unwrap_or("").to_string();
        let mut cc = p2(0., 0.);
        let mut k = 0.0;
        for o in &self.w.objs {
            if o.s.kind == Kind::Pad && o.s.name.split('.').next() == Some(refname.as_str()) {
                cc = cc + o.s.center;
                k += 1.0;
            }
        }
        let cc = if k > 0.0 { cc * (1.0 / k) } else { m };
        let tw = (toward - m).unit();
        let score = |e: P2| e.dot(tw) + if e.dot(m - cc) < -0.05 { -10.0 } else { 0.0 };
        let mut sides = vec![a.perp(), -a.perp()];
        if score(sides[1]) > score(sides[0]) {
            sides.swap(0, 1);
        }
        // pads apart (a protection IC, series resistors): the members split there anyway; the ribbon
        // ends short of them, towards the other station
        let split = self.split_wide && self.wide_station(pp, nn);
        if split {
            sides.insert(0, tw);
        }
        let rr = ri.wpref[l] / 2.0;
        let room = if split { rr + 0.08 } else { rr };
        for e in sides {
            let mut ext: f64 = 0.0;
            for o in [op, on] {
                for q in &o.s.pts {
                    ext = ext.max((*q - m).dot(e) + o.s.r);
                }
            }
            for &extra in &[0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 2.8] {
                for &tilt in &[0.0f64, 20.0, -20.0, 40.0, -40.0, 60.0, -60.0] {
                    let (sn, cs) = tilt.to_radians().sin_cos();
                    let e2 = p2(e.x * cs - e.y * sn, e.x * sn + e.y * cs);
                    let j = m + e2 * (ext + extra);
                    // the fan-out from here to the pads must not cut across other copper
                    let fan = |o: &crate::world::Obj| self.w.stretch_clear(l, j, o.s.center, o.s.net, self.w.b.nets[o.s.net as usize].wpref[l], 0.0);
                    if out.iter().all(|x| x.0.dist(j) > 0.25) && self.w.point_fits(r, j, room, l) && fan(op) && fan(on) {
                        out.push((j, l, rr));
                        if out.len() >= max {
                            return out;
                        }
                    }
                }
            }
        }
        out
    }

    /// "U1202.8/U1202.9 @ (89.25,46.94)"
    pub fn station_label(&self, st: (u32, u32)) -> String {
        let (a, b) = (&self.w.objs[st.0 as usize].s, &self.w.objs[st.1 as usize].s);
        let m = a.center.lerp(b.center, 0.5);
        format!("{}/{} @ ({:.2},{:.2})", a.name, b.name, m.x, m.y)
    }

    /// Other nets with copper (pieces or pads) within r of `at`, nearest first: who crowds a spot.
    pub fn crowd(&self, at: P2, r: f64, except: &[u32]) -> Vec<String> {
        let mut found: Vec<(f64, u32)> = Vec::new();
        let mut note = |d: f64, net: u32| {
            if net == NONE || except.contains(&net) || d > r {
                return;
            }
            match found.iter_mut().find(|x| x.1 == net) {
                Some(x) => x.0 = x.0.min(d),
                None => found.push((d, net)),
            }
        };
        for pc in self.w.pieces.iter().filter(|x| x.alive) {
            let d = pc.poly.windows(2).map(|sg| seg_point_dist(sg[0], sg[1], at)).fold(f64::INFINITY, f64::min);
            note(d, pc.net);
        }
        for o in self.w.objs.iter().filter(|o| o.alive && o.s.kind == Kind::Pad) {
            note(o.s.center.dist(at), o.s.net);
        }
        found.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        found.into_iter().take(5).map(|(_, n)| self.w.b.nets[n as usize].name.clone()).collect()
    }

    /// Nets with pins on fine-pitch parts (16+ pads, neighbours of other nets under 0.7 mm apart), in
    /// fan-out order: pins deep inside a part (inner rows) first, as their lanes out run between the
    /// outer pins; then outer pins round each part. Plane nets are left to the plane drops.
    pub fn fanout_order(&self) -> Vec<u32> {
        let mut parts: HashMap<String, Vec<u32>> = HashMap::default();
        for (i, o) in self.w.objs.iter().enumerate() {
            if o.alive && o.s.kind == Kind::Pad && o.s.net != NONE && !o.s.name.starts_with("junction:") && !o.s.name.starts_with("zone:") {
                parts.entry(o.s.name.split('.').next().unwrap_or("").to_string()).or_default().push(i as u32);
            }
        }
        let mut ranked: Vec<(i64, f64, u32)> = Vec::new();
        for (_, pads) in parts {
            if pads.len() < 16 {
                continue;
            }
            let c = |o: u32| self.w.objs[o as usize].s.center;
            let fine = pads.iter().any(|&a| pads.iter().any(|&b| a != b && self.w.objs[a as usize].s.net != self.w.objs[b as usize].s.net && c(a).dist(c(b)) < 0.7));
            if !fine {
                continue;
            }
            // pin field: the small pads (an exposed pad is not a pin row)
            let small: Vec<u32> = pads.iter().cloned().filter(|&o| {
                let s = &self.w.objs[o as usize].s;
                s.pts.iter().map(|q| q.dist(s.center)).fold(0.0, f64::max) + s.r < 0.8
            }).collect();
            if small.is_empty() {
                continue;
            }
            let (mut lo, mut hi) = (p2(f64::INFINITY, f64::INFINITY), p2(f64::NEG_INFINITY, f64::NEG_INFINITY));
            for &o in &small {
                let q = c(o);
                lo = p2(lo.x.min(q.x), lo.y.min(q.y));
                hi = p2(hi.x.max(q.x), hi.y.max(q.y));
            }
            for &o in &small {
                let net = self.w.objs[o as usize].s.net;
                if self.is_plane(net) || self.w.b.partner.get(net as usize).map_or(false, |&x| x != NONE) || self.net_pads[net as usize].len() < 2 {
                    continue;
                }
                let q = c(o);
                let depth = (q.x - lo.x).min(hi.x - q.x).min(q.y - lo.y).min(hi.y - q.y);
                ranked.push(((depth / 0.2).round() as i64, self.net_span(net), net));
            }
        }
        // deepest first; then short nets first
        ranked.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.partial_cmp(&b.1).unwrap()));
        let mut out = Vec::new();
        for (_, _, n) in ranked {
            if !out.contains(&n) {
                out.push(n);
            }
        }
        out
    }

    fn net_span(&self, net: u32) -> f64 {
        let (mut lo, mut hi) = (p2(f64::INFINITY, f64::INFINITY), p2(f64::NEG_INFINITY, f64::NEG_INFINITY));
        for &o in &self.net_pads[net as usize] {
            let q = self.w.objs[o as usize].s.center;
            lo = p2(lo.x.min(q.x), lo.y.min(q.y));
            hi = p2(hi.x.max(q.x), hi.y.max(q.y));
        }
        (hi.x - lo.x) + (hi.y - lo.y)
    }

    /// A station whose two pads are not neighbours: the members leave the ribbon to reach them singly.
    pub fn wide_station(&self, pp: u32, nn: u32) -> bool {
        self.w.objs[pp as usize].s.center.dist(self.w.objs[nn as usize].s.center) > 1.2
    }

    /// A pair's stations, in path order: (p pad, n pad) couples where both members meet one part
    /// (connector, protection IC, chip, series resistors). The members run through every station.
    pub fn pair_stations(&self, p: u32, n: u32) -> Vec<(u32, u32)> {
        let c = |o: u32| self.w.objs[o as usize].s.center;
        let mut cand: Vec<(f64, u32, u32, P2)> = Vec::new();
        let mut cen = p2(0.0, 0.0);
        let mut k = 0.0;
        for &a in &self.net_pads[p as usize] {
            cen = cen + c(a);
            k += 1.0;
            for &b in &self.net_pads[n as usize] {
                if c(a).dist(c(b)) < 3.0 {
                    cand.push((c(a).dist(c(b)), a, b, c(a).lerp(c(b), 0.5)));
                }
            }
        }
        for &b in &self.net_pads[n as usize] {
            cen = cen + c(b);
            k += 1.0;
        }
        if k == 0.0 {
            return Vec::new();
        }
        let cen = cen * (1.0 / k);
        // closest couples first (an interleaved connector row offers several), ties towards the middle
        cand.sort_by(|x, y| ((x.0 / 0.05).round(), x.3.dist(cen)).partial_cmp(&((y.0 / 0.05).round(), y.3.dist(cen))).unwrap());
        // one station per part (a connector's second row is joined to the first afterwards)
        let part = |o: u32| self.w.objs[o as usize].s.name.split('.').next().unwrap_or("").to_string();
        let mut st: Vec<(u32, u32, P2)> = Vec::new();
        for (_, a, b, m) in cand {
            if st.iter().any(|s| s.0 == a || s.1 == b || s.2.dist(m) < 1.0 || (part(s.0) == part(a) && part(s.1) == part(b))) {
                continue;
            }
            st.push((a, b, m));
        }
        if st.len() < 2 {
            return Vec::new();
        }
        // the two farthest apart are the ends; the others in between along the shortest path
        let (mut i0, mut i1, mut bd) = (0, 1, -1.0);
        for i in 0..st.len() {
            for j in i + 1..st.len() {
                if st[i].2.dist(st[j].2) > bd {
                    bd = st[i].2.dist(st[j].2);
                    i0 = i;
                    i1 = j;
                }
            }
        }
        let mut order = vec![i0];
        let mut rest: Vec<usize> = (0..st.len()).filter(|&i| i != i0 && i != i1).collect();
        while !rest.is_empty() {
            let last = st[*order.last().unwrap()].2;
            let (ri, _) = rest.iter().enumerate().min_by(|x, y| (st[*x.1].2.dist(last) - st[*x.1].2.dist(st[i1].2) * 0.0).partial_cmp(&(st[*y.1].2.dist(last))).unwrap()).unwrap();
            order.push(rest.remove(ri));
        }
        order.push(i1);
        order.into_iter().map(|i| (st[i].0, st[i].1)).collect()
    }

    /// Route one net point to point (a pad or junction to another), optionally coupled to a partner.
    fn connect(&mut self, net: u32, a: u32, b: u32, couple: u32) -> bool {
        self.connect_within(net, a, b, couple, f64::INFINITY)
    }

    /// connect() with a cap on the path cost (a stub must not wander off round the board).
    fn connect_within(&mut self, net: u32, a: u32, b: u32, couple: u32, max_len: f64) -> bool {
        let mut q = self.query(net, vec![a], vec![b], false);
        q.max_len = max_len;
        if couple != NONE {
            q.couple = couple;
            q.couple_factor = 0.6;
        }
        let w = q.width.clone();
        for _ in 0..4 {
            let path = match self.search(&q) {
                Some(x) => x,
                None => return false,
            };
            if self.commit(net, &path, &w).is_ok() {
                return true;
            }
            self.learn(&mut q);
        }
        false
    }

    /// Route a pair station to station: per segment one ribbon (a wide wire with room for both members)
    /// between junctions just off the two stations, unzipped into the members, which end in the station
    /// pads (so a protection IC sits in series on the path). Members that cannot reach a pad straight
    /// from their lane get a point-to-point stub; a segment without a ribbon routes its members singly.
    pub fn route_ribbon(&mut self, p: u32, n: u32) -> bool {
        let r = match self.w.b.ribbon_of.get(p as usize) {
            Some(&r) if r != NONE => r,
            _ => return false,
        };
        let dbg = std::env::var("TOPOR_PAIRDBG").is_ok();
        let st = self.pair_stations(p, n);
        if st.len() < 2 {
            return false;
        }
        let labels: Vec<String> = st.iter().map(|&x| self.station_label(x)).collect();
        let mut segs: Vec<(String, String, String)> = Vec::new();
        for k in 0..st.len() - 1 {
            match self.route_segment(r, p, n, st[k], st[k + 1]) {
                Ok(()) => {
                    if dbg {
                        eprintln!("ribbon {} segment {}: ok", self.w.b.nets[p as usize].name, k);
                    }
                    let split = self.wide_station(st[k].0, st[k].1) || self.wide_station(st[k + 1].0, st[k + 1].1);
                    segs.push((labels[k].clone(), labels[k + 1].clone(), if split { "ribbon (members split to reach spread-out pads)".into() } else { "ribbon".into() }));
                }
                Err(e) => {
                    if dbg {
                        eprintln!("ribbon {} segment {}: {}; members singly", self.w.b.nets[p as usize].name, k, e);
                    }
                    let a = self.connect(p, st[k].0, st[k + 1].0, NONE);
                    let b = self.connect(n, st[k].1, st[k + 1].1, p);
                    // what crowds the two ends, for whoever can move parts
                    let me = [p, n, r];
                    let at = |me2: &Router, x: (u32, u32)| me2.w.objs[x.0 as usize].s.center.lerp(me2.w.objs[x.1 as usize].s.center, 0.5);
                    let ca = self.crowd(at(self, st[k]), 1.5, &me);
                    let cb = self.crowd(at(self, st[k + 1]), 1.5, &me);
                    let first = plain_reason(e.split("; ").next().unwrap_or(""));
                    segs.push((
                        labels[k].clone(),
                        labels[k + 1].clone(),
                        format!("members routed singly{} (ribbon: {}); copper near the ends: [{}] / [{}]", if a && b { "" } else { ", incomplete" }, first, ca.join(", "), cb.join(", ")),
                    ));
                }
            }
        }
        let d = self.pair_diag.entry(p).or_default();
        d.stations = labels;
        d.segments = segs;
        self.route_net(p);
        self.route_net(n);
        let ok = self.deficit(p) == 0 && self.deficit(n) == 0;
        if dbg && !ok {
            eprintln!("ribbon {}: members incomplete", self.w.b.nets[p as usize].name);
            for x in [p, n] {
                let comps: Vec<String> = self.components(x).iter().map(|c| c.iter().map(|&o| { let ob = &self.w.objs[o as usize]; format!("{}@({:.2},{:.2})", ob.s.name, ob.s.center.x, ob.s.center.y) }).collect::<Vec<_>>().join(",")).collect();
                eprintln!("   {}: {}", self.w.b.nets[x as usize].name, comps.join(" | "));
            }
        }
        ok
    }

    /// One ribbon segment between stations sa and sb, stubs included: tried over a few junction
    /// placements until both members join their two station pads.
    fn route_segment(&mut self, r: u32, p: u32, n: u32, sa: (u32, u32), sb: (u32, u32)) -> Result<(), String> {
        let c = |me: &Router, o: u32| me.w.objs[o as usize].s.center;
        let (ma, mb) = (c(self, sa.0).lerp(c(self, sa.1), 0.5), c(self, sb.0).lerp(c(self, sb.1), 0.5));
        // straight into the pads first; at stations whose pads are apart, split with stubs next
        let wide = self.wide_station(sa.0, sa.1) || self.wide_station(sb.0, sb.1);
        let mut tries: Vec<(P2, usize, f64, P2, usize, f64, bool)> = Vec::new();
        for split in [false, true] {
            if split && !wide {
                break;
            }
            self.split_wide = split;
            let ca = self.pair_junctions(sa.0, sa.1, r, mb, 3);
            let cb = self.pair_junctions(sb.0, sb.1, r, ma, 3);
            let mut combos: Vec<(usize, usize)> = Vec::new();
            for i in 0..ca.len() {
                for j in 0..cb.len() {
                    combos.push((i, j));
                }
            }
            combos.sort_by_key(|&(i, j)| (i + j, i.max(j)));
            combos.truncate(4);
            for (i, j) in combos {
                tries.push((ca[i].0, ca[i].1, ca[i].2, cb[j].0, cb[j].1, cb[j].2, split));
            }
        }
        if tries.is_empty() {
            self.split_wide = true;
            let na = self.pair_junctions(sa.0, sa.1, r, mb, 3).len();
            let nb = self.pair_junctions(sb.0, sb.1, r, ma, 3).len();
            return Err(format!("no room for a ribbon junction ({} at {}, {} at {})", na, self.w.objs[sa.0 as usize].s.name, nb, self.w.objs[sb.0 as usize].s.name));
        }
        let mut why = String::new();
        for (ja, la, ra, jb, lb, rb, split) in tries {
            self.split_wide = split;
            let snap = (self.w.clone(), self.ties.clone(), self.history.clone(), self.junctions.clone());
            match self.try_segment(r, p, n, sa, sb, [(ja, la, ra), (jb, lb, rb)]) {
                Ok(()) => {
                    self.split_wide = true;
                    return Ok(());
                }
                Err(e) => {
                    why = if why.is_empty() { e } else { format!("{}; {}", why, e) };
                    self.w = snap.0;
                    self.ties = snap.1;
                    self.history = snap.2;
                    self.junctions = snap.3;
                    self.net_pads[r as usize].clear();
                }
            }
        }
        self.split_wide = true;
        Err(why)
    }

    fn try_segment(&mut self, r: u32, p: u32, n: u32, sa: (u32, u32), sb: (u32, u32), ends: [(P2, usize, f64); 2]) -> Result<(), String> {
        self.w.begin();
        let mut js = [NONE; 2];
        for (k, &(j, l, rr)) in ends.iter().enumerate() {
            match self.w.add_junction(r, j, rr, l) {
                Ok(a) => js[k] = a,
                Err(e) => {
                    self.w.rollback();
                    return Err(format!("junction: {}", e));
                }
            }
        }
        self.w.realize_dirty();
        self.w.commit();
        self.net_pads[r as usize] = js.to_vec();
        if self.route_net(r) > 0 {
            return Err("ribbon route failed".into());
        }
        let stubs = self.unzip(r, p, n, js, [sa.0, sb.0], [sa.1, sb.1]).map_err(|e| format!("unzip failed: {}", e))?;
        self.net_pads[r as usize].clear();
        // stubs from lane ends to their pads, either way round, failures retried after the others
        let mut left = stubs;
        for _round in 0..2 {
            let mut failed = Vec::new();
            for &(net, j, pad) in &left {
                let cap = 3.0 * self.w.objs[j as usize].s.center.dist(self.w.objs[pad as usize].s.center) + 2.0 + 2.0 * self.cfg.via_cost;
                if !(self.connect_within(net, j, pad, NONE, cap) || self.connect_within(net, pad, j, NONE, cap)) {
                    failed.push((net, j, pad));
                }
            }
            if failed.is_empty() || failed.len() == left.len() {
                left = failed;
                break;
            }
            left = failed;
        }
        if let Some(&(_, j, pad)) = left.first() {
            return Err(format!("stub {} -> {} failed", self.w.objs[j as usize].s.name, self.w.objs[pad as usize].s.name));
        }
        let joined = |me: &Router, net: u32, a: u32, b: u32| me.components(net).iter().any(|g| g.contains(&a) && g.contains(&b));
        if !joined(self, p, sa.0, sb.0) || !joined(self, n, sa.1, sb.1) {
            return Err("members not joined station to station".into());
        }
        Ok(())
    }

    /// Replace a routed ribbon by its two members: each member threaded through the triangulation along
    /// the ribbon's offset (so the two run side by side in every crossing list), ending at member
    /// junctions at the ribbon's ends and at side-by-side vias where the ribbon changes layer. The members
    /// are ordinary pieces afterwards (realised, coupled, rubber-banded like any other wire).
    fn unzip(&mut self, r: u32, p: u32, n: u32, js: [u32; 2], ppads: [u32; 2], npads: [u32; 2]) -> Result<Vec<(u32, u32, u32)>, String> {
        let chain = self.piece_chain(r, js[0], js[1]).filter(|c| !c.is_empty()).ok_or("no ribbon chain")?;
        let oriented = |w: &World, pc: u32, from: u32| -> Vec<P2> {
            let piece = &w.pieces[pc as usize];
            let mut v = piece.poly.clone();
            if piece.a.obj != from {
                v.reverse();
            }
            v
        };
        let polys: Vec<(Vec<P2>, usize)> = chain.iter().map(|&(pc, from)| (oriented(&self.w, pc, from), self.w.pieces[pc as usize].layer as usize)).collect();
        if polys.iter().any(|x| x.0.len() < 2) {
            return Err("degenerate ribbon".into());
        }
        let gap = self.w.b.pair_gap[r as usize];
        let (j0, j1) = (self.w.objs[js[0] as usize].s.center, self.w.objs[js[1] as usize].s.center);
        // lanes: p takes the ribbon's left lane when its pad comes first turning counter-clockwise from
        // the ribbon's start direction (so the fan-out does not cross); at the far end, seen from the
        // ribbon arriving, the order must agree or the members swap sides there
        let ang = |t: P2, v: P2| -> f64 {
            let a = t.cross(v).atan2(t.dot(v));
            if a < 0.0 { a + std::f64::consts::TAU } else { a }
        };
        let pc = |o: u32| self.w.objs[o as usize].s.center;
        let first = &polys[0].0;
        let t0 = (first[1] - first[0]).unit();
        let side = if ang(t0, pc(ppads[0]) - j0) < ang(t0, pc(npads[0]) - j0) { 1.0 } else { -1.0 };
        let last = &polys[polys.len() - 1].0;
        let tb = (last[last.len() - 2] - last[last.len() - 1]).unit();
        // arriving, the forward right lane is the reverse left lane
        let (right_pad, left_pad) = if side > 0.0 { (npads[1], ppads[1]) } else { (ppads[1], npads[1]) };
        let swapped = ang(tb, pc(right_pad) - j1) > ang(tb, pc(left_pad) - j1);
        // vias between ribbon pieces: (ribbon via, p position, n position)
        let (vd, _) = self.w.b.via_size(p);
        let half = (vd + gap + 0.02) / 2.0;
        let rv = self.w.b.via_size(r).0 / 2.0;
        let mut vias: Vec<(u32, P2, P2)> = Vec::new();
        for k in 0..chain.len() - 1 {
            let x = {
                let piece = &self.w.pieces[chain[k].0 as usize];
                if piece.a.obj == chain[k].1 { piece.b.obj } else { piece.a.obj }
            };
            let xc = self.w.objs[x as usize].s.center;
            let (pin, pout) = (&polys[k].0, &polys[k + 1].0);
            let din = (pin[pin.len() - 1] - pin[pin.len() - 2]).unit();
            let dout = (pout[1] - pout[0]).unit();
            if din.dot(dout) < -0.3 {
                return Err("ribbon turns back at a via".into());
            }
            let u = (din + dout).unit();
            vias.push((x, xc + u.perp() * (side * half), xc - u.perp() * (side * half)));
        }
        // member polylines: offsets of each ribbon piece, entering via discs straight to the member via
        let mut members: Vec<[Vec<P2>; 2]> = Vec::new();
        for (k, (poly, l)) in polys.iter().enumerate() {
            let wd = self.w.b.nets[p as usize].wpref[*l].max(self.w.b.nets[n as usize].wpref[*l]);
            let d = (wd + gap) / 2.0;
            let mut out: [Vec<P2>; 2] = [Vec::new(), Vec::new()];
            for (mi, sd) in [(0usize, side), (1usize, -side)] {
                let mut m = offset_polyline(poly, sd * d);
                if k > 0 {
                    let c = poly[0];
                    let i = m.iter().position(|q| q.dist(c) > rv + 0.02).unwrap_or(m.len() - 1);
                    m.drain(..i);
                    m.insert(0, if mi == 0 { vias[k - 1].1 } else { vias[k - 1].2 });
                }
                if k + 1 < polys.len() {
                    let c = *poly.last().unwrap();
                    let i = m.iter().rposition(|q| q.dist(c) > rv + 0.02).unwrap_or(0);
                    m.truncate(i + 1);
                    m.push(if mi == 0 { vias[k].1 } else { vias[k].2 });
                }
                // at a split station the lanes end a little apart, so each member's stub has room to
                // leave past the other's end
                let spread = 0.06;
                if k == 0 && self.split_wide && self.wide_station(ppads[0], npads[0]) && m.len() >= 2 {
                    let t = (m[1] - m[0]).unit();
                    m[0] = m[0] + t.perp() * (sd * spread);
                }
                if k + 1 == polys.len() && self.split_wide && self.wide_station(ppads[1], npads[1]) && m.len() >= 2 {
                    let q = m.len();
                    let t = (m[q - 1] - m[q - 2]).unit();
                    m[q - 1] = m[q - 1] + t.perp() * (sd * spread);
                }
                m.dedup_by(|x, y| x.dist(*y) < 1e-6);
                if m.len() < 2 {
                    return Err("member offset collapsed".into());
                }
                out[mi] = m;
            }
            members.push(out);
        }
        self.w.begin();
        for &(pc, _) in &chain {
            self.w.remove_piece(pc);
        }
        for &(x, _, _) in &vias {
            self.w.remove_via(x).map_err(|e| format!("ribbon via: {}", e))?;
        }
        for &j in &js {
            self.w.remove_via(j).map_err(|e| format!("ribbon junction: {}", e))?;
        }
        let mut made: Vec<u32> = Vec::new();
        let mut stubs: Vec<(u32, u32, u32)> = Vec::new();
        let res = (|| -> Result<(), String> {
            // via terminals per piece and member: [start, end] objects (NONE: a pad or junction end)
            let np = polys.len();
            let mut ends: Vec<[[u32; 2]; 2]> = vec![[[NONE; 2]; 2]; np];
            for (k, &(_, vp, vn)) in vias.iter().enumerate() {
                let a = self.w.add_via(p, vp)?;
                let b = self.w.add_via(n, vn)?;
                ends[k][0][1] = a;
                ends[k + 1][0][0] = a;
                ends[k][1][1] = b;
                ends[k + 1][1][0] = b;
            }
            let pads = [[ppads[0], ppads[1]], [npads[0], npads[1]]];
            for k in 0..np {
                let l = polys[k].1;
                for (mi, net) in [(0usize, p), (1usize, n)] {
                    let lane = &members[k][mi];
                    let wd = self.w.b.nets[net as usize].wpref[l];
                    // pad ends first; a junction at the lane end (connected afterwards) where the direct
                    // fan-out does not thread, and at a swapped far end
                    let mut tries: Vec<[bool; 2]> = Vec::new();
                    let split = |me: &Router, a: u32, b: u32| me.split_wide && me.wide_station(a, b);
                    let can = [k == 0 && !split(self, ppads[0], npads[0]), k + 1 == np && !swapped && !split(self, ppads[1], npads[1])];
                    for t in [[true, true], [false, true], [true, false], [false, false]] {
                        let t = [t[0] && can[0], t[1] && can[1]];
                        if !tries.contains(&t) {
                            tries.push(t);
                        }
                    }
                    let mut done = false;
                    for t in tries {
                        let mut poly = lane.clone();
                        let mut terms = [Term { obj: NONE, vtx: NONE }; 2];
                        let mut fresh: Vec<u32> = Vec::new();
                        let mut bad = false;
                        for side_i in 0..2 {
                            let at = if side_i == 0 { lane[0] } else { *lane.last().unwrap() };
                            let via = ends[k][mi][side_i];
                            if via != NONE {
                                terms[side_i] = Term { obj: via, vtx: self.w.objs[via as usize].lv[l] };
                                continue;
                            }
                            if t[side_i] {
                                match self.pad_term(pads[mi][side_i], l, at) {
                                    Some((tm, q)) => {
                                        terms[side_i] = tm;
                                        // a polygon pad is entered across the hull edge nearest the lane
                                        let mut add = vec![q];
                                        if tm.vtx == NONE {
                                            // through the hull edge facing the ribbon's end, off its corners
                                            // (not between neighbouring pads; a line through a corner runs
                                            // along the pad's own diagonal)
                                            let hull = &self.w.objs[pads[mi][side_i] as usize].s.pts;
                                            let jc = if side_i == 0 { j0 } else { j1 };
                                            let toward = (jc - q).unit();
                                            let mut best: Option<(f64, usize)> = None;
                                            for i in 0..hull.len() {
                                                let (a, b) = (hull[i], hull[(i + 1) % hull.len()]);
                                                let mut nrm = (b - a).perp().unit();
                                                if nrm.dot(a.lerp(b, 0.5) - q) < 0.0 {
                                                    nrm = -nrm;
                                                }
                                                let f = nrm.dot(toward);
                                                if best.map_or(true, |x| f > x.0) {
                                                    best = Some((f, i));
                                                }
                                            }
                                            if let Some((_, i)) = best {
                                                let (a, b) = (hull[i], hull[(i + 1) % hull.len()]);
                                                let d = b - a;
                                                let u = ((at - a).dot(d) / d.dot(d)).clamp(0.25, 0.75);
                                                let hq = a.lerp(b, u);
                                                let mut nrm = d.perp().unit();
                                                if nrm.dot(hq - q) < 0.0 {
                                                    nrm = -nrm;
                                                }
                                                add.push(hq + nrm * 0.02);
                                            }
                                        }
                                        if side_i == 0 {
                                            poly.splice(0..0, add);
                                        } else {
                                            poly.extend(add.into_iter().rev());
                                        }
                                    }
                                    None => bad = true,
                                }
                            } else {
                                match self.w.add_junction(net, at, wd / 2.0, l) {
                                    Ok(j) => {
                                        fresh.push(j);
                                        terms[side_i] = Term { obj: j, vtx: self.w.objs[j as usize].lv[l] };
                                    }
                                    Err(_) => bad = true,
                                }
                            }
                        }
                        if !bad {
                            for (side_i, tm) in terms.iter().enumerate() {
                                if tm.vtx != NONE {
                                    let q = self.w.vpos(l, tm.vtx);
                                    if side_i == 0 {
                                        poly[0] = q;
                                    } else {
                                        *poly.last_mut().unwrap() = q;
                                    }
                                }
                            }
                            if let Some(cx) = self.w.thread_poly(l, &poly, terms[0], terms[1]) {
                                let pid = self.w.add_piece(net, l, wd, terms[0], terms[1], &cx);
                                self.w.seed_poly(pid, poly);
                                made.extend(fresh.iter().cloned());
                                for side_i in 0..2 {
                                    if ends[k][mi][side_i] == NONE && !t[side_i] {
                                        stubs.push((net, terms[side_i].obj, pads[mi][side_i]));
                                    }
                                }
                                done = true;
                                break;
                            }
                        }
                        for j in fresh {
                            let _ = self.w.remove_via(j);
                        }
                    }
                    if !done {
                        return Err(format!("threading member {} piece {}", mi, k));
                    }
                }
            }
            if std::env::var("TOPOR_VALIDATE").is_ok() {
                self.w.validate_topology()?;
            }
            let re = self.w.realize_dirty();
            let mut v = Vec::new();
            for &q in &re {
                self.w.check_piece(q, &mut v);
            }
            if let Some(x) = v.first() {
                if let Ok(dir) = std::env::var("TOPOR_UNZIPDUMP") {
                    let b = [x.at.x - 2.5, x.at.y - 2.0, x.at.x + 2.5, x.at.y + 2.0];
                    let svg = crate::svg::render(&self.w, &crate::svg::Opts { mesh: true, box_: Some(b), scale: 70.0, marks: v.iter().map(|y| (y.at, y.what.clone())).collect(), unrouted: vec![] });
                    let _ = std::fs::write(format!("{}/unzip-{}.svg", dir, self.w.b.nets[p as usize].name.replace('/', "_")), svg);
                    for y in v.iter().take(4) {
                        eprintln!("   unzip viol {} at ({:.3},{:.3}) d={:.3} need {:.3} pieces {:?} objs {:?}", y.what, y.at.x, y.at.y, y.dist, y.need, y.pieces.iter().map(|&q| self.w.b.nets[self.w.pieces[q as usize].net as usize].name.clone()).collect::<Vec<_>>(), y.objs.iter().map(|&o| self.w.objs[o as usize].s.name.clone()).collect::<Vec<_>>());
                    }
                }
                return Err(format!("{} violations, first {} at ({:.3},{:.3})", v.len(), x.what, x.at.x, x.at.y));
            }
            Ok(())
        })();
        match res {
            Ok(()) => {
                self.w.commit();
                self.junctions.extend(made);
                Ok(stubs)
            }
            Err(e) => {
                self.w.rollback();
                Err(e)
            }
        }
    }

    /// A wire terminal on pad o on layer l: a polygon pad's centre, else its vertex nearest `toward`.
    fn pad_term(&self, o: u32, l: usize, toward: P2) -> Option<(Term, P2)> {
        let ob = &self.w.objs[o as usize];
        if !ob.alive || !ob.on(l) {
            return None;
        }
        if ob.is_poly() {
            return Some((Term { obj: o, vtx: NONE }, ob.s.center));
        }
        let k = (0..ob.npts()).min_by(|&a, &b| ob.s.pts[a].dist(toward).partial_cmp(&ob.s.pts[b].dist(toward)).unwrap())?;
        let v = ob.v(l, k);
        if v == NONE {
            return None;
        }
        Some((Term { obj: o, vtx: v }, self.w.vpos(l, v)))
    }

    /// Remove a pair entirely: both members' pieces and vias and its unzip junctions.
    pub fn clear_pair(&mut self, p: u32, n: u32) {
        for x in [p, n] {
            let ps: Vec<u32> = (0..self.w.pieces.len() as u32).filter(|&q| self.w.pieces[q as usize].alive && self.w.pieces[q as usize].net == x).collect();
            for q in ps {
                self.rip_piece(q);
            }
        }
        let js: Vec<u32> = self.junctions.iter().cloned().filter(|&j| {
            let o = &self.w.objs[j as usize];
            o.alive && (o.s.net == p || o.s.net == n)
        }).collect();
        for j in js {
            let _ = self.w.remove_via(j);
        }
        self.junctions.retain(|&j| self.w.objs[j as usize].alive);
    }

    /// Reconnect a net after it lost pieces; a differential pair member is rebuilt with its partner.
    pub fn reroute(&mut self, net: u32) {
        let pn = self.w.b.partner.get(net as usize).cloned().unwrap_or(NONE);
        if pn == NONE || self.deficit(net) == 0 {
            self.route_net(net);
            return;
        }
        let (p, n) = (net.min(pn), net.max(pn));
        self.clear_pair(p, n);
        self.route_pair(p, n);
        if self.deficit(p) > 0 {
            self.route_net(p);
        }
        if self.deficit(n) > 0 {
            self.route_net(n);
        }
        let cf = self.pair_report().into_iter().find(|x| self.w.b.net_index.get(&x.0) == Some(&p)).map(|x| x.4).unwrap_or(0.0);
        let done = self.deficit(p) + self.deficit(n) == 0;
        if let Some(m) = self.pair_diag.entry(p).or_default().moved.last_mut() {
            m.push_str(&format!("; rebuilt {}, {:.0}% coupled", if done { "complete" } else { "incomplete" }, cf * 100.0));
        }
    }

    /// The end stations of a differential pair: ((p_a, n_a), (p_b, n_b)).
    pub fn pair_ends(&self, p: u32, n: u32) -> Option<((u32, u32), (u32, u32))> {
        let st = self.pair_stations(p, n);
        if st.len() < 2 {
            return None;
        }
        Some((st[0], st[st.len() - 1]))
    }

    /// Route differential pairs first: one member reserving room for the ribbon, the other alongside.
    pub fn route_pairs(&mut self) {
        let nn = self.w.b.nets.len() as u32;
        let mut pairs: Vec<(f64, u32, u32)> = Vec::new();
        for p in 0..nn {
            let n = self.w.b.partner.get(p as usize).cloned().unwrap_or(NONE);
            if n == NONE || n < p {
                continue;
            }
            if let Some(((pa, _), (pb, _))) = self.pair_ends(p, n) {
                let d = self.w.objs[pa as usize].s.center.dist(self.w.objs[pb as usize].s.center);
                pairs.push((d, p, n));
            }
        }
        // which pair goes first decides who has to make way where pairs crowd a chip's pad row or
        // cross: try short-first and long-first, keep the better coupled result
        let start = (self.w.clone(), self.ties.clone(), self.history.clone(), self.junctions.clone());
        let mut best: Option<(f64, (World, Vec<(u32, u32)>, Vec<f64>, Vec<u32>), HashMap<u32, PairDiag>)> = None;
        let orders: Vec<&str> = match std::env::var("TOPOR_PAIRORDER").as_deref() {
            Ok("long") => vec!["long"],
            Ok("short") => vec!["short"],
            _ => vec!["short", "long"],
        };
        for (oi, order) in orders.iter().enumerate() {
            if oi > 0 {
                self.w = start.0.clone();
                self.ties = start.1.clone();
                self.history = start.2.clone();
                self.junctions = start.3.clone();
            }
            let mut ps = pairs.clone();
            if *order == "long" {
                ps.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
            } else {
                ps.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            }
            for &(_, p, n) in &ps {
                self.route_pair(p, n);
                self.route_net(p);
                self.route_net(n);
            }
            let score: f64 = ps.iter().map(|&(_, p, n)| if self.deficit(p) + self.deficit(n) > 0 { -1.0 } else { 0.0 }).sum::<f64>() + self.pair_report().iter().map(|x| x.4).sum::<f64>();
            if std::env::var("TOPOR_PAIRDBG").is_ok() {
                eprintln!("pair order {}: score {:.2}", order, score);
            }
            if best.as_ref().map_or(true, |b| score > b.0) {
                best = Some((score, (self.w.clone(), self.ties.clone(), self.history.clone(), self.junctions.clone()), self.pair_diag.clone()));
            }
        }
        let mut score = 0.0;
        if let Some((sc, st, dg)) = best {
            score = sc;
            self.w = st.0;
            self.ties = st.1;
            self.history = st.2;
            self.junctions = st.3;
            self.pair_diag = dg;
        }
        for (a, _, _, _, cf) in self.pair_report() {
            if let Some(&p) = self.w.b.net_index.get(&a) {
                self.pair_diag.entry(p).or_default().coupled_after_pairs = cf;
            }
        }
        if std::env::var("TOPOR_NEGOTIATE").is_ok() {
            self.negotiate_pairs(&pairs, score);
        }
    }

    fn pairs_score(&self, pairs: &[(f64, u32, u32)]) -> f64 {
        pairs.iter().map(|&(_, p, n)| if self.deficit(p) + self.deficit(n) > 0 { -1.0 } else { 0.0 }).sum::<f64>() + self.pair_report().iter().map(|x| x.4).sum::<f64>()
    }

    /// A poorly coupled pair and a pair crowding its stations swap who goes first, kept if the pairs
    /// come out better coupled overall.
    fn negotiate_pairs(&mut self, pairs: &[(f64, u32, u32)], mut score: f64) {
        let dbg = std::env::var("TOPOR_PAIRDBG").is_ok();
        let report = self.pair_report();
        let coupled = |me: &Router, p: u32| report.iter().find(|x| me.w.b.net_index.get(&x.0) == Some(&p) || me.w.b.net_index.get(&x.1) == Some(&p)).map(|x| x.4).unwrap_or(1.0);
        let weak: Vec<(u32, u32)> = pairs.iter().filter(|x| coupled(self, x.1) < 0.6).map(|x| (x.1, x.2)).collect();
        for (wp, wn) in weak {
            let st: Vec<P2> = self.pair_stations(wp, wn).iter().flat_map(|&(a, b)| [self.w.objs[a as usize].s.center, self.w.objs[b as usize].s.center]).collect();
            let near = |me: &Router, net: u32| me.w.pieces.iter().any(|pc| pc.alive && pc.net == net && pc.poly.windows(2).any(|sg| st.iter().any(|&c| seg_point_dist(sg[0], sg[1], c) < 2.0)));
            let blockers: Vec<(u32, u32)> = pairs.iter().filter(|x| x.1 != wp && (near(self, x.1) || near(self, x.2))).map(|x| (x.1, x.2)).collect();
            for (op, on) in blockers {
                let snap = (self.w.clone(), self.ties.clone(), self.history.clone(), self.junctions.clone());
                self.clear_pair(wp, wn);
                self.clear_pair(op, on);
                for (p, n) in [(wp, wn), (op, on)] {
                    self.route_pair(p, n);
                    self.route_net(p);
                    self.route_net(n);
                }
                let sc = self.pairs_score(pairs);
                if dbg {
                    eprintln!("negotiate {} before {}: score {:.2} -> {:.2}", self.w.b.nets[wp as usize].name, self.w.b.nets[op as usize].name, score, sc);
                }
                if sc > score + 0.05 {
                    score = sc;
                    break;
                }
                self.w = snap.0;
                self.ties = snap.1;
                self.history = snap.2;
                self.junctions = snap.3;
            }
        }
    }

    /// One member reserving room for the ribbon, the other following it.
    pub fn route_pair(&mut self, p: u32, n: u32) {
        if self.pair_stations(p, n).len() < 2 {
            return;
        }
        if self.w.b.ribbon_of.get(p as usize).map_or(false, |&r| r != NONE) {
            self.route_ribbon(p, n);
            return;
        }
        {
            let ((pa, na), (pb, nb)) = self.pair_ends(p, n).unwrap();
            let gap = self.w.b.pair_gap[p as usize];
            let mut q = self.query(p, vec![pa], vec![pb], false);
            let w = q.width.clone();
            q.ribbon = w.iter().map(|&x| 2.0 * x + gap).collect();
            let mut ok = false;
            for ribbon in [true, false] {
                if !ribbon {
                    q.ribbon.clear();
                }
                for _ in 0..4 {
                    let path = match self.search(&q) {
                        Some(x) => x,
                        None => break,
                    };
                    if self.commit(p, &path, &w).is_ok() {
                        ok = true;
                        break;
                    }
                    self.learn(&mut q);
                }
                if ok {
                    break;
                }
            }
            if !ok {
                self.log.push(format!("pair {}: first member failed", self.w.b.nets[p as usize].name));
                return;
            }
            let followed = self.follow_partner(p, n, pa, na, pb, nb);
            if std::env::var("TOPOR_PAIRDBG").is_ok() {
                // crossings of n adjacent to a crossing of p on the same edge
                let (mut adj, mut tot) = (0, 0);
                for pc in self.w.pieces.iter().filter(|x| x.alive && x.net == n) {
                    let l = pc.layer as usize;
                    for &c in &pc.cx {
                        tot += 1;
                        let e = self.w.crosses[c as usize].edge;
                        let list = &self.w.layers[l].cross[e as usize];
                        let i = list.iter().position(|&x| x == c).unwrap();
                        let nb = |j: usize| j < list.len() && self.w.pieces[self.w.crosses[list[j] as usize].piece as usize].net == p;
                        if (i > 0 && nb(i - 1)) || nb(i + 1) {
                            adj += 1;
                        }
                    }
                }
                let pcx: usize = self.w.pieces.iter().filter(|x| x.alive && x.net == p).map(|x| x.cx.len()).sum();
                eprintln!("pair {} / {}: follow {}; n crossings adjacent to p: {}/{} (p has {})", self.w.b.nets[p as usize].name, self.w.b.nets[n as usize].name, followed, adj, tot, pcx);
            }
            if !followed {
                // independent route alongside where it can
                let mut q2 = self.query(n, vec![na], vec![nb], false);
                q2.couple = p;
                q2.couple_factor = 0.6;
                let w2 = q2.width.clone();
                for _ in 0..4 {
                    let path = match self.search(&q2) {
                        Some(x) => x,
                        None => break,
                    };
                    if self.commit(n, &path, &w2).is_ok() {
                        break;
                    }
                    self.learn(&mut q2);
                }
            }
        }
    }

    /// Pieces of a net forming the path from object a to b, each with the object it starts from.
    pub fn piece_chain(&self, net: u32, a: u32, b: u32) -> Option<Vec<(u32, u32)>> {
        let objs = self.net_objs(net);
        let mut prev: Vec<Option<(u32, u32)>> = vec![None; objs.len()];
        let idx = |o: u32| objs.iter().position(|&x| x == o);
        let ia = idx(a)?;
        let ib = idx(b)?;
        let mut seen = vec![false; objs.len()];
        seen[ia] = true;
        let mut queue = std::collections::VecDeque::from(vec![ia]);
        while let Some(u) = queue.pop_front() {
            if u == ib {
                break;
            }
            for &pc in &self.w.obj_pieces[objs[u] as usize] {
                let piece = &self.w.pieces[pc as usize];
                if !piece.alive {
                    continue;
                }
                let other = if piece.a.obj == objs[u] { piece.b.obj } else { piece.a.obj };
                if let Some(j) = idx(other) {
                    if !seen[j] {
                        seen[j] = true;
                        prev[j] = Some((pc, objs[u]));
                        queue.push_back(j);
                    }
                }
            }
        }
        if !seen[ib] {
            return None;
        }
        let mut chain = Vec::new();
        let mut cur = ib;
        while cur != ia {
            let (pc, from) = prev[cur]?;
            chain.push((pc, from));
            cur = idx(from)?;
        }
        chain.reverse();
        Some(chain)
    }

    /// Route pair member n alongside p: vias beside p's vias, then each layer segment between those
    /// terminals with a coupled single-layer search, committed as one transaction.
    fn follow_partner(&mut self, p: u32, n: u32, pa: u32, na: u32, pb: u32, nb: u32) -> bool {
        let dbg = std::env::var("TOPOR_PAIRDBG").is_ok();
        let why = |m: &str| {
            if dbg {
                eprintln!("   follow: {}", m);
            }
        };
        let chain = match self.piece_chain(p, pa, pb) {
            Some(c) if !c.is_empty() => c,
            _ => {
                why("no chain for the first member");
                return false;
            }
        };
        let oriented = |w: &World, pc: u32, from: u32| -> Vec<P2> {
            let piece = &w.pieces[pc as usize];
            let mut v = piece.poly.clone();
            if piece.a.obj != from {
                v.reverse();
            }
            v
        };
        let p0 = oriented(&self.w, chain[0].0, chain[0].1);
        if p0.len() < 2 {
            return false;
        }
        let d0 = (p0[p0.len().min(3) - 1] - p0[0]).unit();
        let side = d0.cross(self.w.objs[na as usize].s.center - self.w.objs[pa as usize].s.center).signum();
        {
            let last = &chain[chain.len() - 1];
            let pl = oriented(&self.w, last.0, last.1);
            let d1 = (pl[pl.len() - 1] - pl[pl.len().saturating_sub(3).min(pl.len() - 2)]).unit();
            let side_end = d1.cross(self.w.objs[nb as usize].s.center - self.w.objs[pb as usize].s.center).signum();
            if side_end != side {
                why(&format!("pin order swapped between the ends ({} vs {})", side, side_end));
            }
        }
        let c = &self.w.b.classes[self.w.b.nets[n as usize].class];
        let pitch = c.via_d + self.w.b.pair_gap[n as usize].max(c.clear) + 0.02;
        let mut via_pos = Vec::new();
        for k in 0..chain.len() - 1 {
            let x = {
                let piece = &self.w.pieces[chain[k].0 as usize];
                if piece.a.obj == chain[k].1 { piece.b.obj } else { piece.a.obj }
            };
            let xc = self.w.objs[x as usize].s.center;
            let pin = oriented(&self.w, chain[k].0, chain[k].1);
            let pout = oriented(&self.w, chain[k + 1].0, chain[k + 1].1);
            let din = (pin[pin.len() - 1] - pin[pin.len().saturating_sub(3).min(pin.len() - 2)]).unit();
            let dout = (pout[pout.len().min(3) - 1] - pout[0]).unit();
            let u = (din + dout).unit();
            let nrm = u.perp() * side;
            let mut found = None;
            'cand: for &scale in &[1.0, 1.15, 1.3, 1.6, 2.0] {
                for &ang in &[0.0f64, 15.0, -15.0, 30.0, -30.0, 45.0, -45.0, 60.0, -60.0, 75.0, -75.0, 90.0, -90.0, 120.0, -120.0, 180.0] {
                    let (sn, cs) = ang.to_radians().sin_cos();
                    let dir = p2(nrm.x * cs - nrm.y * sn, nrm.x * sn + nrm.y * cs);
                    let cand = xc + dir * (pitch * scale);
                    if self.w.via_fits(n, cand) {
                        found = Some(cand);
                        break 'cand;
                    }
                }
            }
            match found {
                Some(v) => via_pos.push(v),
                None => {
                    why(&format!("no room for a partner via beside {:?}", xc));
                    return false;
                }
            }
        }
        let width = self.w.b.nets[n as usize].wpref.clone();
        self.w.begin();
        let mut terms = vec![na];
        for &v in &via_pos {
            match self.w.add_via(n, v) {
                Ok(o) => terms.push(o),
                Err(_) => {
                    self.w.rollback();
                    return false;
                }
            }
        }
        terms.push(nb);
        let mut new_vias: Vec<u32> = terms[1..terms.len() - 1].to_vec();
        for k in 0..chain.len() {
            let l = self.w.pieces[chain[k].0 as usize].layer as usize;
            let mut q = self.query(n, vec![terms[k]], vec![terms[k + 1]], false);
            q.allowed &= 1 << l;
            q.couple = p;
            // follow, do not wander: a detour much longer than the partner's stretch means the sides swap
            let plen: f64 = self.w.pieces[chain[k].0 as usize].poly.windows(2).map(|s| s[0].dist(s[1])).sum();
            q.max_len = (1.3 * plen + 2.0) * q.couple_factor.max(0.3) + 2.0;
            let path = {
                let mut se = Search::new(&self.w, &q);
                se.run()
            };
            let path = match path {
                Some(x) if x.segs.len() == 1 => x,
                other => {
                    why(&format!("segment {} on layer {} search {:?}", k, l, other.map(|x| x.segs.len())));
                    self.w.rollback();
                    return false;
                }
            };
            let sg = &path.segs[0];
            let (a, b) = match (sg.a, sg.b) {
                (End::Obj(a), End::Obj(b)) => (a, b),
                _ => {
                    self.w.rollback();
                    return false;
                }
            };
            self.w.add_piece(n, l, width[l], a, b, &sg.cx);
        }
        let realized = self.w.realize_dirty();
        let mut viol: Vec<Viol> = Vec::new();
        for &pp in &realized {
            if self.w.pieces[pp as usize].overlap {
                viol.push(Viol { what: "overflow".into(), layer: 0, at: P2::default(), dist: 0., need: 0., pieces: vec![pp], objs: vec![] });
                break;
            }
            self.w.check_piece(pp, &mut viol);
            if !viol.is_empty() {
                break;
            }
        }
        for &o in &new_vias {
            if viol.is_empty() {
                let c = self.w.objs[o as usize].s.center;
                self.w.check_via_at(n, c, o, &mut viol, true);
            }
        }
        new_vias.clear();
        if !viol.is_empty() {
            why(&format!("{} at {:?}", viol[0].what, viol[0].at));
            self.w.rollback();
            return false;
        }
        self.w.commit();
        true
    }

    /// Path length between two objects of a net through its pieces (Dijkstra over objects).
    pub fn net_path_len(&self, net: u32, a: u32, b: u32) -> Option<f64> {
        let objs = self.net_objs(net);
        let mut dist: Vec<f64> = vec![f64::INFINITY; objs.len()];
        let idx = |o: u32| objs.iter().position(|&x| x == o);
        let ia = idx(a)?;
        let ib = idx(b)?;
        dist[ia] = 0.0;
        let mut done = vec![false; objs.len()];
        loop {
            let mut u = None;
            for i in 0..objs.len() {
                if !done[i] && dist[i].is_finite() && u.map_or(true, |j: usize| dist[i] < dist[j]) {
                    u = Some(i);
                }
            }
            let u = u?;
            if u == ib {
                return Some(dist[u]);
            }
            done[u] = true;
            for &pc in &self.w.obj_pieces[objs[u] as usize] {
                let piece = &self.w.pieces[pc as usize];
                if !piece.alive {
                    continue;
                }
                let other = if piece.a.obj == objs[u] { piece.b.obj } else { piece.a.obj };
                if let Some(j) = idx(other) {
                    let len: f64 = piece.poly.windows(2).map(|s| s[0].dist(s[1])).sum();
                    if dist[u] + len < dist[j] {
                        dist[j] = dist[u] + len;
                    }
                }
            }
        }
    }

    /// Per pair: member path lengths between the pair's ends, mismatch, and the coupled fraction.
    pub fn pair_report(&self) -> Vec<(String, String, f64, f64, f64)> {
        let mut out = Vec::new();
        for p in 0..self.w.b.nets.len() as u32 {
            let n = self.w.b.partner.get(p as usize).cloned().unwrap_or(NONE);
            if n == NONE || n < p {
                continue;
            }
            let ((pa, na), (pb, nb)) = match self.pair_ends(p, n) {
                Some(x) => x,
                None => continue,
            };
            let lp = self.net_path_len(p, pa, pb).unwrap_or(f64::NAN);
            let ln = self.net_path_len(n, na, nb).unwrap_or(f64::NAN);
            // coupled fraction: length of n's pad-to-pad path within the coupling distance of p's copper
            // on the same layer (sampled every 50 um)
            let gap = self.w.b.pair_gap[p as usize];
            let mut tot = 0.0;
            let mut cpl = 0.0;
            let chain: Vec<u32> = self.piece_chain(n, na, nb).map(|c| c.iter().map(|x| x.0).collect()).unwrap_or_default();
            for &pi in &chain {
                let pc = &self.w.pieces[pi as usize];
                for sgm in pc.poly.windows(2) {
                    let len = sgm[0].dist(sgm[1]);
                    let k = (len / 0.05).ceil().max(1.0) as usize;
                    for i in 0..k {
                        let m = sgm[0].lerp(sgm[1], (i as f64 + 0.5) / k as f64);
                        tot += len / k as f64;
                        let near = self.w.pieces.iter().filter(|x| x.alive && x.net == p && x.layer == pc.layer).any(|q| {
                            q.poly.windows(2).any(|t| seg_point_dist(t[0], t[1], m) < (q.width + pc.width) / 2.0 + gap + 0.03)
                        });
                        if near {
                            cpl += len / k as f64;
                        }
                    }
                }
            }
            out.push((self.w.b.nets[p as usize].name.clone(), self.w.b.nets[n as usize].name.clone(), lp, ln, if tot > 0.0 { cpl / tot } else { 0.0 }));
        }
        out
    }

    fn trace_pairs(&self, stage: &str) {
        if std::env::var("TOPOR_PAIRTRACE").is_err() {
            return;
        }
        let r = self.pair_report();
        let s: Vec<String> = r.iter().map(|x| format!("{}:{:.0}%", x.0.split('/').last().unwrap_or("").trim_end_matches("_DM").trim_end_matches("_DP"), x.4 * 100.0)).collect();
        eprintln!("pairs after {}: {}", stage, s.join(" "));
    }

    pub fn route_all(&mut self) {
        let nn = self.w.b.nets.len() as u32;
        if self.cfg.only.as_deref() == Some("pairs") {
            self.route_pairs();
            return;
        }
        if let Some(name) = self.cfg.only.clone() {
            if let Some(&net) = self.w.b.net_index.get(&name) {
                let d = self.route_net(net);
                eprintln!("only {}: deficit {}", name, d);
            }
            return;
        }
        // differential pairs claim their corridors first, then plane drops (TOPOR_FANOUT: the nets of
        // fine-pitch parts before the planes; on the test boards it did not pay off)
        self.pair_len = vec![0.0; nn as usize];
        for p in 0..nn {
            if let Some(((pa, _), (pb, _))) = self.w.b.partner.get(p as usize).filter(|&&x| x != NONE).and_then(|&x| self.pair_ends(p.min(x), p.max(x))) {
                self.pair_len[p as usize] = self.w.objs[pa as usize].s.center.dist(self.w.objs[pb as usize].s.center);
            }
        }
        self.route_pairs();
        if self.cfg.verbose {
            eprintln!("[{:.1}s] pairs done", self.elapsed());
        }
        self.trace_pairs("pairs");
        // the rest of the board; routing is order-sensitive where fine-pitch parts crowd, so other
        // orders run at once (the fine-pitch parts fanned out first), the best one kept. An attempt that
        // completes the board calls off the ones after it; the earlier ones run on, so the first complete
        // attempt in this order always wins.
        let has_fanout = !self.fanout_order().is_empty();
        let forced: Option<usize> = std::env::var("TOPOR_ATTEMPT").ok().and_then(|x| x.parse().ok());
        // (long nets first, with or without the fan-out, did worse on every test board: TOPOR_ATTEMPT=2, 3)
        let variants: Vec<(bool, bool)> = if forced.is_some() { vec![(false, false), (true, false), (false, true), (true, true)] } else if has_fanout { vec![(false, false), (true, false)] } else { vec![(false, false)] };
        if let Some(k) = forced {
            let (fan, long) = variants[k.min(variants.len() - 1)];
            self.route_rest(fan, long);
        } else if std::env::var("TOPOR_FANOUT").is_ok() {
            self.route_rest(true, false);
        } else if self.cfg.single_attempt || variants.len() == 1 {
            self.route_rest(false, false);
        } else {
            let stops: Vec<std::sync::Arc<std::sync::atomic::AtomicBool>> = variants.iter().map(|_| std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false))).collect();
            let mut others: Vec<Router> = (1..variants.len())
                .map(|k| {
                    let mut o = self.clone();
                    o.cfg.verbose = false;
                    o.stop = Some(stops[k].clone());
                    o
                })
                .collect();
            let call_off_after = |k: usize, stops: &[std::sync::Arc<std::sync::atomic::AtomicBool>]| {
                for s in &stops[k + 1..] {
                    s.store(true, std::sync::atomic::Ordering::Relaxed);
                }
            };
            let done: Vec<(usize, Router)> = std::thread::scope(|sc| {
                let hs: Vec<_> = others
                    .drain(..)
                    .enumerate()
                    .map(|(i, mut o)| {
                        let (fan, long) = variants[i + 1];
                        let stops = &stops;
                        sc.spawn(move || {
                            o.route_rest(fan, long);
                            if !o.stopped() && o.total_deficit() == 0 {
                                call_off_after(i + 1, stops);
                            }
                            o
                        })
                    })
                    .collect();
                self.route_rest(false, false);
                if self.total_deficit() == 0 {
                    call_off_after(0, &stops);
                }
                hs.into_iter().enumerate().map(|(i, h)| (i + 1, h.join().expect("routing attempt"))).collect()
            });
            let name = |k: usize| match variants[k] {
                (false, false) => "plain order",
                (true, false) => "fine-pitch parts fanned out first",
                (false, true) => "long nets first",
                (true, true) => "fan-out first, then long nets first",
            };
            let first = self.total_deficit();
            let mut best: Option<(usize, usize)> = None;
            let mut summary = vec![format!("{}: deficit {}", name(0), first)];
            for (k, o) in &done {
                if o.stopped() {
                    summary.push(format!("{}: called off", name(*k)));
                    continue;
                }
                let d = o.total_deficit();
                summary.push(format!("{}: deficit {}", name(*k), d));
                if d < best.map_or(first, |b| b.1) {
                    best = Some((*k, d));
                }
            }
            if self.cfg.verbose {
                eprintln!("[{:.1}s] {}", self.elapsed(), summary.join("; "));
            }
            if let Some((k, d)) = best {
                let mut o = done.into_iter().find(|x| x.0 == k).unwrap().1;
                o.cfg.verbose = self.cfg.verbose;
                o.stop = None;
                o.t0 = self.t0;
                o.events.push(format!("the board completed better with the {} (deficit {} instead of {} in plain order)", name(k), d, first));
                *self = o;
            }
        }
        self.improve_pairs();
        self.cleanup();
        self.trace_pairs("cleanup");
        self.match_pairs(0.1);
    }

    fn stopped(&self) -> bool {
        self.stop.as_ref().map_or(false, |s| s.load(std::sync::atomic::Ordering::Relaxed))
    }

    /// Planes, signals, rip-up and the end-game repairs, after the pairs.
    fn route_rest(&mut self, fanout: bool, long_first: bool) {
        let nn = self.w.b.nets.len() as u32;
        if fanout {
            let fan = self.fanout_order();
            for &net in &fan {
                self.route_net(net);
            }
            if self.cfg.verbose {
                eprintln!("[{:.1}s] fan-out of fine-pitch parts: {} nets, deficit {}", self.elapsed(), fan.len(), fan.iter().map(|&n| self.deficit(n)).sum::<usize>());
            }
        }
        for net in 0..nn {
            if self.is_plane(net) {
                self.route_net(net);
            }
        }
        if self.cfg.verbose {
            eprintln!("[{:.1}s] planes done", self.elapsed());
        }
        self.trace_pairs("planes");
        let mut order: Vec<(f64, u32)> = Vec::new();
        for net in 0..nn {
            if self.is_plane(net) || self.net_pads[net as usize].len() < 2 {
                continue;
            }
            let mut lo = p2(f64::INFINITY, f64::INFINITY);
            let mut hi = p2(f64::NEG_INFINITY, f64::NEG_INFINITY);
            for &o in &self.net_pads[net as usize] {
                let c = self.w.objs[o as usize].s.center;
                lo = p2(lo.x.min(c.x), lo.y.min(c.y));
                hi = p2(hi.x.max(c.x), hi.y.max(c.y));
            }
            order.push(((hi.x - lo.x) + (hi.y - lo.y), net));
        }
        order.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        if long_first {
            order.reverse();
        }
        for &(_, net) in &order {
            self.route_net(net);
            if self.elapsed() > self.cfg.seconds || self.stopped() {
                break;
            }
        }
        if self.stopped() {
            return;
        }
        if self.cfg.verbose {
            eprintln!("[{:.1}s] first pass: deficit {}", self.elapsed(), self.total_deficit());
        }
        self.trace_pairs("first pass");
        // rip-up and reroute
        let mut best_def = usize::MAX;
        let mut stale = 0;
        // the best board seen: a round can make things worse, the flow then goes on from the best
        type Snap = (World, Vec<(u32, u32)>, Vec<f64>, Vec<u32>, HashMap<u32, PairDiag>);
        let mut best_state: Option<(usize, Snap)> = None;
        for round in 0..self.cfg.rip_rounds {
            if self.stopped() {
                return;
            }
            self.round = round as i32;
            // rounds that no longer make progress go to the pair rescue instead
            let d = self.total_deficit();
            if best_state.as_ref().map_or(true, |b| d < b.0) {
                best_state = Some((d, (self.w.clone(), self.ties.clone(), self.history.clone(), self.junctions.clone(), self.pair_diag.clone())));
            }
            if d < best_def {
                best_def = d;
                stale = 0;
            } else {
                stale += 1;
                // with pairs released, rounds are cheap and the last few nets need patience
                if stale >= if self.relax_pairs { 10 } else { 4 } {
                    // stuck with pairs held: let the remaining nets move them, once
                    if self.relax_pairs {
                        break;
                    }
                    self.relax_pairs = true;
                    self.budget_scale = 3.0;
                    stale = 0;
                    let stuck: Vec<String> = (0..nn).filter(|&x| self.deficit(x) > 0).map(|x| self.w.b.nets[x as usize].name.clone()).collect();
                    let ev = format!("rip-up stalled at round {} with deficit {} ({}) while moving differential pairs was costly (and barred to plane drops); from here on any net may move pairs", round, d, stuck.join(", "));
                    if self.cfg.verbose {
                        eprintln!("[{:.1}s] {}", self.elapsed(), ev);
                    }
                    self.events.push(ev);
                }
            }
            // pairs stand firm first; whatever is still blocked later may move them
            self.protect_pairs = round < std::env::var("TOPOR_PROTECT").ok().and_then(|x| x.parse().ok()).unwrap_or(0);
            let failing: Vec<u32> = (0..nn).filter(|&n| self.deficit(n) > 0 && self.net_pads[n as usize].len() >= 1).collect();
            if failing.is_empty() || self.elapsed() > self.cfg.seconds {
                break;
            }
            for net in failing {
                if self.elapsed() > self.cfg.seconds {
                    break;
                }
                if self.deficit(net) == 0 {
                    continue;
                }
                let victims = self.repair_net(net);
                // a ripped pair is rebuilt first: it needs a corridor, the net only a track
                let paired = |me: &Router, v: u32| me.w.b.partner.get(v as usize).map_or(false, |&x| x != NONE);
                if std::env::var("TOPOR_PAIRFIRST").is_ok() {
                    for &v in &victims {
                        if paired(self, v) {
                            self.reroute(v);
                        }
                    }
                }
                self.route_net(net);
                for &v in &victims {
                    if !paired(self, v) || std::env::var("TOPOR_PAIRFIRST").is_err() {
                        self.reroute(v);
                    }
                }
            }
            if self.cfg.verbose {
                eprintln!("[{:.1}s] rip round {}: deficit {}", self.elapsed(), round, self.total_deficit());
            }
            self.trace_pairs(&format!("rip round {}", round));
        }
        if let Some((d, st)) = best_state {
            if d < self.total_deficit() {
                if self.cfg.verbose {
                    eprintln!("[{:.1}s] back to the best rip-up state (deficit {} instead of {})", self.elapsed(), d, self.total_deficit());
                }
                self.w = st.0;
                self.ties = st.1;
                self.history = st.2;
                self.junctions = st.3;
                self.pair_diag = st.4;
            }
        }
        if self.stopped() {
            return;
        }
        self.clear_area_repairs();
        self.rescue_pairs();
    }

    /// Nets still open after rip-up: a pad sealed in by other nets' vias and wires the soft search cannot
    /// pass. Clear other copper round the open pads, route the net first, then the displaced nets;
    /// kept only if the board ends up less incomplete.
    fn clear_area_repairs(&mut self) {
        let nn = self.w.b.nets.len() as u32;
        for radius in [1.2, 2.5] {
            let failing: Vec<u32> = (0..nn).filter(|&x| self.deficit(x) > 0 && self.net_pads[x as usize].len() >= 1).collect();
            for net in failing {
                if self.deficit(net) == 0 || self.elapsed() > self.cfg.seconds * 1.5 {
                    continue;
                }
                // pads of the groups that are not connected (all but the largest; every group for planes)
                let comps = self.components(net);
                let open: Vec<P2> = comps
                    .iter()
                    .skip(if self.is_plane(net) { 0 } else { 1 })
                    .filter(|c| !self.is_plane(net) || !self.plane_ok(net, c))
                    .flatten()
                    .filter(|&&o| self.w.objs[o as usize].s.kind == Kind::Pad)
                    .map(|&o| self.w.objs[o as usize].s.center)
                    .collect();
                if open.is_empty() {
                    continue;
                }
                let near = |q: P2| open.iter().any(|c| c.dist(q) < radius);
                let mut victims: Vec<u32> = Vec::new();
                for (i, pc) in self.w.pieces.iter().enumerate() {
                    if pc.alive && pc.net != net && pc.poly.windows(2).any(|sg| open.iter().any(|&c| seg_point_dist(sg[0], sg[1], c) < radius)) {
                        victims.push(i as u32);
                    }
                }
                for (oi, o) in self.w.objs.iter().enumerate() {
                    if o.alive && o.s.kind == Kind::Via && o.s.net != net && near(o.s.center) {
                        victims.extend(self.w.obj_pieces[oi].iter().cloned().filter(|&p| self.w.pieces[p as usize].alive));
                    }
                }
                victims.sort();
                victims.dedup();
                if victims.is_empty() {
                    continue;
                }
                let before = self.total_deficit();
                let snap = (self.w.clone(), self.ties.clone(), self.history.clone(), self.junctions.clone());
                let dsnap = self.pair_diag.clone();
                let mut vnets: Vec<u32> = victims.iter().map(|&p| self.w.pieces[p as usize].net).collect();
                vnets.sort();
                vnets.dedup();
                for p in victims {
                    self.rip_piece(p);
                }
                self.route_net(net);
                for v in &vnets {
                    self.reroute(*v);
                }
                for v in &vnets {
                    if self.deficit(*v) > 0 {
                        let more = self.repair_net(*v);
                        self.route_net(*v);
                        for m in more {
                            self.reroute(m);
                        }
                    }
                }
                let after = self.total_deficit();
                let kept = after < before;
                let ev = format!("{} was sealed in near {}: cleared {:.1} mm round its open pads ({} nets moved), deficit {} -> {}{}", self.w.b.nets[net as usize].name, open.iter().take(2).map(|c| format!("({:.2},{:.2})", c.x, c.y)).collect::<Vec<_>>().join(" "), radius, vnets.len(), before, after, if kept { "" } else { ", undone" });
                if self.cfg.verbose {
                    eprintln!("[{:.1}s] {}", self.elapsed(), ev);
                }
                if kept {
                    self.events.push(ev);
                } else {
                    self.w = snap.0;
                    self.ties = snap.1;
                    self.history = snap.2;
                    self.junctions = snap.3;
                    self.pair_diag = dsnap;
                }
            }
        }
    }

    /// On the finished board, rebuild poorly coupled pairs as ribbons; a rebuild is kept only if the
    /// board stays as complete and the pair comes out better coupled.
    fn improve_pairs(&mut self) {
        let report = self.pair_report();
        let total = self.total_deficit();
        for (a, _, lp0, ln0, cf) in report {
            if cf >= 0.6 || self.elapsed() > self.cfg.seconds * 1.5 {
                continue;
            }
            let len0 = if lp0.is_finite() && ln0.is_finite() { lp0.max(ln0) } else { f64::INFINITY };
            let p = self.w.b.net_index[&a];
            let n = self.w.b.partner[p as usize];
            let (p, n) = (p.min(n), p.max(n));
            let snap = (self.w.clone(), self.ties.clone(), self.history.clone(), self.junctions.clone());
            let dsnap = self.pair_diag.clone();
            self.clear_pair(p, n);
            self.route_pair(p, n);
            self.route_net(p);
            self.route_net(n);
            let rep = self.pair_report().into_iter().find(|x| x.0 == a);
            let after = rep.as_ref().map(|x| x.4).unwrap_or(0.0);
            let len1 = rep.as_ref().map(|x| x.2.max(x.3)).unwrap_or(f64::INFINITY);
            // better coupled, not bought with a detour
            let keep = self.total_deficit() <= total && after > cf + 0.05 && len1 <= len0 * 1.15 + 1.0;
            if self.cfg.verbose {
                eprintln!("[{:.1}s] improve pair {}: {:.0}% -> {:.0}% {}", self.elapsed(), a, cf * 100.0, after * 100.0, if keep { "kept" } else { "reverted" });
            }
            if keep {
                let note = format!("rebuilt as ribbon on the finished board: {:.0}% -> {:.0}% coupled", cf * 100.0, after * 100.0);
                self.pair_diag.entry(p).or_default().notes.push(note);
            }
            if !keep {
                self.w = snap.0;
                self.ties = snap.1;
                self.history = snap.2;
                self.junctions = snap.3;
                self.pair_diag = dsnap;
            }
        }
    }

    /// A pair still incomplete after rip-up: route its members as plain nets (pads in any order), kept
    /// only if the board ends up less incomplete; then a few more rip-up rounds.
    fn rescue_pairs(&mut self) {
        let nn = self.w.b.nets.len() as u32;
        let mut any = false;
        for p in 0..nn {
            let n = self.w.b.partner.get(p as usize).cloned().unwrap_or(NONE);
            if n == NONE || n < p || self.deficit(p) + self.deficit(n) == 0 {
                continue;
            }
            let before = self.total_deficit();
            let snap = (self.w.clone(), self.ties.clone(), self.history.clone(), self.junctions.clone());
            let dsnap = self.pair_diag.clone();
            self.clear_pair(p, n);
            self.route_net(p);
            self.route_net(n);
            for x in [p, n] {
                if self.deficit(x) > 0 {
                    let victims = self.repair_net(x);
                    self.route_net(x);
                    for v in victims {
                        if v != p && v != n {
                            self.reroute(v);
                        }
                    }
                }
            }
            let after = self.total_deficit();
            if self.cfg.verbose {
                eprintln!("[{:.1}s] rescue pair {}: deficit {} -> {}", self.elapsed(), self.w.b.nets[p as usize].name, before, after);
            }
            if after < before {
                any = true;
                let ev = format!("pair {} / {} still incomplete after rip-up: routed as two plain nets (board deficit {} -> {}), not coupled", self.w.b.nets[p as usize].name, self.w.b.nets[n as usize].name, before, after);
                self.pair_diag.entry(p).or_default().notes.push(ev.clone());
                self.events.push(ev);
            } else {
                self.w = snap.0;
                self.ties = snap.1;
                self.history = snap.2;
                self.junctions = snap.3;
                self.pair_diag = dsnap;
            }
        }
        if any {
            for _ in 0..3 {
                let failing: Vec<u32> = (0..nn).filter(|&x| self.deficit(x) > 0 && self.net_pads[x as usize].len() >= 1).collect();
                if failing.is_empty() {
                    break;
                }
                for net in failing {
                    if self.deficit(net) == 0 {
                        continue;
                    }
                    let victims = self.repair_net(net);
                    self.route_net(net);
                    for v in victims {
                        self.reroute(v);
                    }
                }
            }
        }
    }

    /// Remove vias that connect nothing (one piece or none) and their lone pieces.
    pub fn cleanup(&mut self) {
        loop {
            let mut changed = false;
            for o in 0..self.w.objs.len() as u32 {
                let ob = &self.w.objs[o as usize];
                if !ob.alive || ob.s.kind != Kind::Via || self.ties.iter().any(|t| t.1 == o) {
                    continue;
                }
                let live: Vec<u32> = self.w.obj_pieces[o as usize].iter().cloned().filter(|&p| self.w.pieces[p as usize].alive).collect();
                // a plane drop has one wire; a via of a plane net off its pour is just a dangling via
                let plane = ob.plane;
                if live.is_empty() || (live.len() == 1 && !plane) {
                    self.w.begin();
                    for p in live {
                        self.w.remove_piece(p);
                    }
                    match self.w.remove_via(o) {
                        Ok(_) => {
                            let re = self.w.realize_dirty();
                            let mut v = Vec::new();
                            for &q in &re {
                                self.w.check_piece(q, &mut v);
                            }
                            if v.is_empty() {
                                self.w.commit();
                                changed = true;
                            } else {
                                self.w.rollback();
                            }
                        }
                        Err(_) => self.w.rollback(),
                    }
                }
            }
            if !changed {
                break;
            }
        }
        // signal vias whose wires all sit on one layer: rip them and reconnect the net without them
        for o in 0..self.w.objs.len() as u32 {
            let ob = &self.w.objs[o as usize];
            if !ob.alive || ob.s.kind != Kind::Via || ob.plane {
                continue;
            }
            let live: Vec<u32> = self.w.obj_pieces[o as usize].iter().cloned().filter(|&p| self.w.pieces[p as usize].alive).collect();
            let mut layers: Vec<u8> = live.iter().map(|&p| self.w.pieces[p as usize].layer).collect();
            layers.sort();
            layers.dedup();
            if layers.len() >= 2 {
                continue;
            }
            let net = ob.s.net;
            for p in live {
                self.rip_piece(p);
            }
            for _ in 0..3 {
                if self.route_net(net) == 0 {
                    break;
                }
                for v in self.repair_net(net) {
                    self.route_net(v);
                }
            }
        }
    }

    /// Length matching: serpentine bumps on the shorter member of each pair, on a straight stretch, on the
    /// side away from its partner, exact-checked. Final geometry only (nothing is routed afterwards).
    pub fn match_pairs(&mut self, tolerance: f64) {
        for p in 0..self.w.b.nets.len() as u32 {
            let n = self.w.b.partner.get(p as usize).cloned().unwrap_or(NONE);
            if n == NONE || n < p {
                continue;
            }
            let ((pa, na), (pb, nb)) = match self.pair_ends(p, n) {
                Some(x) => x,
                None => continue,
            };
            for _round in 0..6 {
                let lp = match self.net_path_len(p, pa, pb) {
                    Some(x) => x,
                    None => break,
                };
                let ln = match self.net_path_len(n, na, nb) {
                    Some(x) => x,
                    None => break,
                };
                let delta = (lp - ln).abs();
                if delta <= tolerance {
                    break;
                }
                let (short, sa, sb, other) = if lp < ln { (p, pa, pb, n) } else { (n, na, nb, p) };
                let chain = match self.piece_chain(short, sa, sb) {
                    Some(c) => c,
                    None => break,
                };
                if !self.add_meander(&chain, other, delta) {
                    break;
                }
            }
        }
    }

    /// Insert bumps adding `delta` of length into one straight stretch of the chain's pieces.
    fn add_meander(&mut self, chain: &[(u32, u32)], other: u32, delta: f64) -> bool {
        let mut cands: Vec<(f64, u32, usize)> = Vec::new();
        for &(pc, _) in chain {
            let pl = &self.w.pieces[pc as usize].poly;
            for k in 0..pl.len().saturating_sub(1) {
                let len = pl[k].dist(pl[k + 1]);
                if len > 0.8 {
                    cands.push((len, pc, k));
                }
            }
        }
        cands.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        for (len, pc, k) in cands.into_iter().take(12) {
            let piece = &self.w.pieces[pc as usize];
            let l = piece.layer as usize;
            let w = piece.width;
            let (a, b) = (piece.poly[k], piece.poly[k + 1]);
            let dir = (b - a).unit();
            // bump away from the partner
            let mid = a.lerp(b, 0.5);
            let mut toward = 0.0;
            for q in self.w.pieces.iter().filter(|x| x.alive && x.net == other && x.layer as usize == l) {
                for sgm in q.poly.windows(2) {
                    let f = sgm[0].lerp(sgm[1], seg_param(sgm[0], sgm[1], mid));
                    if f.dist(mid) < 1.0 {
                        toward += dir.cross(f - mid).signum();
                    }
                }
            }
            let sides: Vec<f64> = if toward > 0.0 { vec![-1.0] } else if toward < 0.0 { vec![1.0] } else { vec![1.0, -1.0] };
            let pitch = 4.0 * w + 0.3;
            for side in sides {
                for nb in 1..=((len - 0.6) / pitch).floor().max(0.0) as usize {
                    let amp = delta / (2.0 * nb as f64);
                    if amp < w || amp > 1.2 {
                        continue;
                    }
                    let nrm = dir.perp() * side;
                    let span = nb as f64 * pitch;
                    let t0 = (len - span) / 2.0;
                    let mut pts = vec![a];
                    for i in 0..nb {
                        let s0 = t0 + i as f64 * pitch;
                        let base = a + dir * s0;
                        // trapezoid bump: up, across, down (adds 2 * amp)
                        pts.push(base + dir * (0.25 * pitch));
                        pts.push(base + dir * (0.25 * pitch) + nrm * amp);
                        pts.push(base + dir * (0.75 * pitch) + nrm * amp);
                        pts.push(base + dir * (0.75 * pitch));
                    }
                    pts.push(b);
                    let mut np = self.w.pieces[pc as usize].poly.clone();
                    np.splice(k..=k + 1, pts);
                    if self.w.try_set_poly(pc, np) {
                        return true;
                    }
                }
            }
        }
        false
    }

    pub fn total_deficit(&self) -> usize {
        (0..self.w.b.nets.len() as u32).map(|n| self.deficit(n)).sum()
    }

    pub fn incomplete_nets(&self) -> Vec<u32> {
        (0..self.w.b.nets.len() as u32).filter(|&n| self.deficit(n) > 0).collect()
    }
}
