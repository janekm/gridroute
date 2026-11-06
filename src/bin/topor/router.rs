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
    /// Quick look for layout iteration: one routing order, short rip-up, no end-game repairs, no polishing
    /// (pairs, spring, via bridges); the blocker report says what keeps each open net out.
    pub fast: bool,
}

/// What keeps an open group of a net out, for the design side: the cheapest way to it through other nets'
/// copper (the nets it would cross, and where it runs: a route to draw by hand, or nets to move), or none if
/// fixed copper and the outline seal it in (then a part must move); and what is next to its open pads.
#[derive(Clone, Debug)]
pub struct Blocker {
    pub net: String,
    pub open_pads: Vec<(String, P2)>,
    /// None: no way even through other nets' copper.
    pub way: Option<Way>,
    /// With no way: true if the search ran out of places to go (sealed in by fixed copper, keepouts or the
    /// outline), false if it ran out of budget first.
    pub sealed: bool,
    /// Other nets' copper within 1 mm of the open pads: (net, routed wires and vias, fixed copper: pads, hand
    /// drawn tracks and vias).
    pub near_nets: Vec<(String, usize, usize)>,
    /// Other parts with pads within 1 mm of the open pads.
    pub near_parts: Vec<String>,
}

/// What became of the suggested copper (see import_suggestions).
#[derive(Clone, Debug, Default)]
pub struct Suggested {
    pub vias_kept: usize,
    pub vias_left_out: usize,
    pub tracks_kept: usize,
    pub tracks_left_out: usize,
    pub notes: Vec<String>,
    /// Indices into the model's suggested vias and tracks that broke a rule (not those reaching no pad).
    pub bad_vias: Vec<usize>,
    pub bad_tracks: Vec<usize>,
}

#[derive(Clone, Debug)]
pub struct Way {
    /// The nets it crosses, with how many of their wires.
    pub through: Vec<(String, usize)>,
    /// Where it crosses them (triangle centres, layer name).
    pub crossings: Vec<(String, P2)>,
    /// Its course per layer: from its start over the edges it crosses to its end.
    pub legs: Vec<(String, Vec<P2>)>,
    pub vias: Vec<P2>,
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
    /// Per dogbone net: where its vias may go (right beside its pads, clear of the fixed copper).
    pub dogbone_spots: Vec<Option<std::sync::Arc<Vec<P2>>>>,
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
    /// Set when the other routing attempts are done and this one has had its grace: wrap up (keep the best
    /// state so far, start no more repairs); the result still counts.
    pub hurry: Option<std::sync::Arc<std::sync::atomic::AtomicBool>>,
    /// Soft searches that found nothing, by net and what they depend on (the net's source group and targets,
    /// whether pairs may be crossed, the budget): other nets' copper is passable to them, so the same search
    /// fails again until one of those changes, and is not run again.
    pub soft_dead: std::collections::HashSet<(u32, u64)>,
    /// Pieces imported from the model's suggested copper (hand-drawn tracks routed as the router's own), and
    /// their nets: rip-up moves them at a small premium, pair stations treat them like fixed copper.
    pub suggested: std::collections::HashSet<u32>,
    pub suggested_nets: std::collections::HashSet<u32>,
    /// Pair members whose pads the suggested copper joins already (a pair drawn by hand, as found by the import):
    /// the model keeps their suggested copper as drawn (see main).
    pub drawn_pairs: std::collections::HashSet<u32>,
    /// Per net: the coupled length (mm) of the differential pair it belongs to, as last measured; crossing it
    /// costs in proportion, so rip-up does not break good coupling for a short cut.
    pub pair_coupled: Vec<f64>,
    /// Pair stations that fixed copper already joins are merged into one (see pair_stations); off for the
    /// attempts that route the pairs as before, to compare.
    pub merge_stations: bool,
    /// Pairs routed as before: no merged stations, no ribbon start at the pad tips, no ribbon over the pair's own
    /// copper, no lane retreat for a member's stub, no cost for lost coupling (see route_all).
    pub legacy_pairs: bool,
    /// Wall time and calls per phase of the run (the attempt that is kept), for the time report.
    pub prof: std::collections::BTreeMap<&'static str, (f64, usize)>,
    /// High-speed nets keep off gaps in their reference planes (gaps.rs); off for the attempts that route
    /// the board without the rule, to compare (see route_all).
    pub gap_rule: bool,
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
    /// Soft searches that found nothing, and their time; repeats of one that is known to fail, skipped.
    pub soft_fail: usize,
    pub soft_fail_s: f64,
    pub soft_skipped: usize,
}

/// Extra cost per millimetre of coupling a rip of a pair would lose (TOPOR_COUPLE_COST).
fn couple_cost() -> f64 {
    static C: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *C.get_or_init(|| std::env::var("TOPOR_COUPLE_COST").ok().and_then(|x| x.parse().ok()).unwrap_or(2.0))
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
        // and a neck-down stub its pad
        for (zi, z) in w.objs.iter().enumerate() {
            if !crate::world::fixed_copper(&z.s.name) || z.s.net == NONE {
                continue;
            }
            for &oi in &net_pads[z.s.net as usize] {
                let o = &w.objs[oi as usize];
                let once = !crate::world::fixed_copper(&o.s.name) || oi as usize > zi;
                if once && o.s.layers & z.s.layers != 0 && crate::world::core_dist(&z.s, &o.s) - o.s.r - z.s.r < 0.03 {
                    ties.push((zi as u32, oi));
                }
            }
        }
        let n = w.b.nets.len();
        let dogbone_spots = (0..n).map(|net| if w.b.dogbone.get(net).copied().unwrap_or(false) { Some(std::sync::Arc::new(dogbone_spots(&w, net as u32, &net_pads[net]))) } else { None }).collect();
        Router { w, cfg, net_pads, dogbone_spots, history: vec![0.0; n], log: Vec::new(), t0: Instant::now(), stats: Stats::default(), dumps: 0, last_viol: None, last_exp: 0, ties, junctions: Vec::new(), split_wide: true, protect_pairs: false, pair_len: Vec::new(), relax_pairs: false, pair_diag: HashMap::default(), events: Vec::new(), round: -1, budget_scale: 1.0, stop: None, tri_hist: std::sync::Arc::new(HashMap::default()), hurry: None, soft_dead: Default::default(), suggested: Default::default(), suggested_nets: Default::default(), drawn_pairs: Default::default(), pair_coupled: Vec::new(), merge_stations: true, legacy_pairs: false, prof: Default::default(), gap_rule: true }
    }

    pub fn elapsed(&self) -> f64 {
        self.t0.elapsed().as_secs_f64()
    }

    /// Add the time since `since` to phase `name` of the time report.
    pub fn prof_add(&mut self, name: &'static str, since: Instant) {
        let e = self.prof.entry(name).or_insert((0.0, 0));
        e.0 += since.elapsed().as_secs_f64();
        e.1 += 1;
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
        // (an index and one pass over the ties: the groups come out as before, in the same order, without the
        // quadratic lookups that big nets made expensive)
        let mut index: HashMap<u32, usize> = HashMap::default();
        for (i, &o) in objs.iter().enumerate() {
            index.entry(o).or_insert(i);
        }
        let idx = |o: u32| index.get(&o).cloned();
        let mut par: Vec<usize> = (0..objs.len()).collect();
        let plane = self.is_plane(net);
        let mut plane_root: Option<usize> = None;
        for &(x, y) in &self.ties {
            if let (Some(i), Some(j)) = (idx(x), idx(y)) {
                // as seen from either end: the other end must be alive
                if self.w.objs[y as usize].alive || self.w.objs[x as usize].alive {
                    let (a, b) = (find(&mut par, i), find(&mut par, j));
                    if a != b {
                        par[b] = a;
                    }
                }
            }
        }
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

    /// Whether pad `o` has a neck-down stub (tied to it).
    pub fn has_neck(&self, o: u32) -> bool {
        self.ties.iter().any(|&(x, y)| (x == o && self.w.objs[y as usize].s.name.starts_with("neck:")) || (y == o && self.w.objs[x as usize].s.name.starts_with("neck:")))
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

    pub(crate) fn query(&self, net: u32, sources: Vec<u32>, targets: Vec<u32>, narrow: bool) -> Query {
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
            // an RF line can only go back on its one layer: as costly to move as a pair
            let rf = self.w.b.rf_nets.iter().any(|x| x.0 as usize == i);
            // (a coupled pair costs more for every millimetre of coupling a rip would lose)
            let coupled = if paired && !self.legacy_pairs { couple_cost() * self.pair_coupled.get(i).cloned().unwrap_or(0.0) } else { 0.0 };
            // (suggested copper stays unless the net in the way needs it: a small premium)
            let sug = if self.suggested_nets.contains(&(i as u32)) { 10.0 } else { 0.0 };
            cross_cost[i] = if self.w.b.nets[i].plane.is_some() { 2.5 } else if paired || rf { pair_cost + coupled } else { 6.0 } + h + sug;
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
            layer_cost: if ni.layer_cost.len() == nl { ni.layer_cost.clone() } else { vec![1.0; nl] },
            max_len: f64::INFINITY,
            ignore: vec![],
            taboo: HashMap::default(),
            via_taboo: vec![],
            history: self.tri_hist.clone(),
            wedge_off: std::env::var("TOPOR_WEDGE_OFF").is_ok(),
            couple: NONE,
            couple_factor: 0.3,
            ribbon: vec![],
            gap_rule: self.gap_rule,
            via_spots: self.dogbone_spots.get(net as usize).cloned().flatten(),
            via_spot_layers: self.net_pads.get(net as usize).map_or(0, |ps| ps.iter().fold(0, |m, &o| m | self.w.objs[o as usize].s.layers)),
        }
    }

    /// High-speed wires (as routed) crossing a slot in a reference plane past their net's limit.
    pub fn gap_violations(&self) -> usize {
        if !crate::gaps::active(&self.w) {
            return 0;
        }
        let tracks: Vec<(u32, usize, f64, Vec<P2>)> = self.w.pieces.iter().filter(|p| p.alive && p.poly.len() >= 2 && self.w.b.high_speed(p.net)).map(|p| (p.net, p.layer as usize, p.width, p.poly.clone())).collect();
        crate::gaps::Gaps::new(&self.w).report(&tracks).iter().filter(|c| c.violation()).count()
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
            if r.is_none() {
                self.stats.soft_fail += 1;
                self.stats.soft_fail_s += t.elapsed().as_secs_f64();
            }
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
            let msg = format!("{} on layer {} at ({:.3},{:.3}) d={:.4} need {:.4}{}", v.what, v.layer, v.at.x, v.at.y, v.dist, v.need, v.objs.first().map(|&o| format!(" ({})", self.w.objs[o as usize].s.name)).unwrap_or_default());
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
        // a pad a neck-down rule covers drops at the neck width, if the drop is one straight segment from it
        // (a track that touches the part's courtyard)
        let neck = comp.iter().filter_map(|&o| self.w.b.neck_width.get(&self.w.objs[o as usize].s.name).cloned()).fold(f64::INFINITY, f64::min);
        if neck.is_finite() {
            let mut q = self.query(net, comp.to_vec(), targets.clone(), true);
            q.plane_goal = true;
            q.max_len = 2.0;
            let mut layers = 0u32;
            for &o in comp {
                layers |= self.w.objs[o as usize].s.layers;
            }
            q.allowed &= layers;
            for x in q.width.iter_mut() {
                *x = x.min(neck);
            }
            let w = q.width.clone();
            for _attempt in 0..4 {
                let path = match self.search(&q) {
                    Some(p) => p,
                    None => break,
                };
                match self.commit(net, &path, &w) {
                    Ok(ps) => {
                        let straight = ps.iter().all(|&p| {
                            let pc = &self.w.pieces[p as usize];
                            pc.poly.len() == 2 && (comp.contains(&pc.a.obj) || comp.contains(&pc.b.obj))
                        });
                        if straight {
                            return true;
                        }
                        for p in ps {
                            self.rip_piece(p);
                        }
                        self.learn(&mut q);
                    }
                    Err(_) => self.learn(&mut q),
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
        let dead_key = {
            use std::hash::{Hash, Hasher};
            let mut h = std::collections::hash_map::DefaultHasher::new();
            let mut sorted = src.clone();
            sorted.sort();
            (sorted, targets.len(), self.relax_pairs, self.protect_pairs, self.budget_scale.to_bits(), self.w.b.nets[net as usize].allowed).hash(&mut h);
            (net, h.finish())
        };
        if self.soft_dead.contains(&dead_key) {
            self.stats.soft_skipped += 1;
            return vec![];
        }
        let mut q = self.query(net, src.clone(), targets.clone(), true);
        q.soft = true;
        // (the quick look takes a smaller soft budget: what still needs more is reported as blocked)
        q.max_exp = (self.cfg.max_exp as f64 * self.budget_scale) as usize * if self.cfg.fast { 3 } else { 7 };
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
        // a plane pad a neck-down rule covers (a power pin of a fine-pitch part) may leave its pad only at the neck
        // width (route_plane_comp's last try): when the class width finds no way out, the soft search looks for
        // what is in the way of that drop instead
        let neck = if self.is_plane(net) && std::env::var("TOPOR_NO_PLANE_NECK_REPAIR").is_err() {
            src.iter().filter_map(|&o| self.w.b.neck_width.get(&self.w.objs[o as usize].s.name).cloned()).fold(f64::INFINITY, f64::min)
        } else {
            f64::INFINITY
        };
        let mut necked = false;
        let mut ripped_nets: Vec<u32> = Vec::new();
        let mut path = None;
        for _attempt in 0..5 {
            let p = match self.search(&q) {
                Some(p) => p,
                None if neck.is_finite() && !necked => {
                    necked = true;
                    for x in q.width.iter_mut() {
                        *x = x.min(neck);
                    }
                    q.max_len = 2.0;
                    continue;
                }
                None => {
                    if dbg {
                        eprintln!("repair {}: soft search failed (src {:?} targets {:?}) after {} expansions", self.w.b.nets[net as usize].name, src.iter().map(|&o| self.w.objs[o as usize].s.name.clone()).collect::<Vec<_>>(), targets.iter().map(|&o| self.w.objs[o as usize].s.name.clone()).collect::<Vec<_>>(), self.last_exp);
                    }
                    if ripped_nets.is_empty() {
                        self.soft_dead.insert(dead_key);
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
                // (a drop at the neck width is one straight segment from its pad, else it is taken up again and
                // the drop left to route_net)
                Ok(ps) if necked && !ps.iter().all(|&x| self.w.pieces[x as usize].poly.len() == 2 && (src.contains(&self.w.pieces[x as usize].a.obj) || src.contains(&self.w.pieces[x as usize].b.obj))) => {
                    for x in ps {
                        self.rip_piece(x);
                    }
                    return ripped_nets;
                }
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
                    // nothing this net may not cross at any price (a pair while pairs are protected, or for a plane
                    // net, whose drop can go elsewhere): the search was kept off it, and the exact check is not
                    // a way round that
                    culprits.retain(|&pp| q.cross_cost.get(self.w.pieces[pp as usize].net as usize).map_or(true, |c| c.is_finite()));

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
        // pads already about as far apart as the coupled members (fine-pitch pins side by side): the members are
        // together at the pad tips, so the ribbon may start right there, without room to fan in
        let coupled = self.w.b.nets[op.s.net as usize].wpref[l].max(self.w.b.nets[on.s.net as usize].wpref[l]) + self.w.b.pair_gap[op.s.net as usize];
        let close = !self.legacy_pairs && !split && cp.dist(cn) <= coupled + 0.1;
        let extras: &[f64] = if close { &[0.05, 0.15, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 2.8] } else { &[0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 2.8] };
        for e in sides {
            let mut ext: f64 = 0.0;
            for o in [op, on] {
                for q in &o.s.pts {
                    ext = ext.max((*q - m).dot(e) + o.s.r);
                }
            }
            for &extra in extras {
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

    /// Pin rows of fine-pitch parts (16+ pads; elongated pads, rectangular or oval, with a neighbour of
    /// another net closer than 0.7 mm), each row ordered along the part's side: (pad, pad's outer end,
    /// outward direction).
    fn escape_rows(&self) -> Vec<Vec<(u32, P2, P2)>> {
        let mut parts: HashMap<String, Vec<u32>> = HashMap::default();
        for (i, o) in self.w.objs.iter().enumerate() {
            if o.alive && o.s.kind == Kind::Pad && o.s.net != NONE && !o.s.name.starts_with("junction:") && !crate::world::fixed_copper(&o.s.name) {
                parts.entry(o.s.name.split('.').next().unwrap_or("").to_string()).or_default().push(i as u32);
            }
        }
        let mut names: Vec<String> = parts.keys().cloned().collect();
        names.sort();
        let mut out = Vec::new();
        for name in names {
            let pads = &parts[&name];
            if pads.len() < 16 {
                continue;
            }
            let c = |o: u32| self.w.objs[o as usize].s.center;
            let cc = pads.iter().fold(p2(0.0, 0.0), |a, &o| a + c(o)) * (1.0 / pads.len() as f64);
            let mut rows: HashMap<i32, Vec<(f64, u32, P2, P2)>> = HashMap::default();
            for &o in pads {
                let s = &self.w.objs[o as usize].s;
                // polygon pads, or oval ones (a segment core with a radius)
                if s.pts.len() < 2 || (s.pts.len() == 2 && s.r <= 0.0) || !pads.iter().any(|&q| q != o && self.w.objs[q as usize].s.net != s.net && c(q).dist(s.center) < 0.7) {
                    continue;
                }
                let n = s.pts.len();
                let (mut ax, mut best) = (p2(1.0, 0.0), 0.0);
                for k in 0..n {
                    let e = s.pts[(k + 1) % n] - s.pts[k];
                    if e.norm() > best {
                        best = e.norm();
                        ax = e.unit();
                    }
                }
                if ax.dot(s.center - cc) < 0.0 {
                    ax = -ax;
                }
                let (mut far, mut hw) = (0.0f64, 0.0f64);
                for q in &s.pts {
                    far = far.max((*q - s.center).dot(ax) + s.r);
                    hw = hw.max((*q - s.center).dot(ax.perp()).abs() + s.r);
                }
                // not a pin (an exposed pad), or no clear long axis
                if far > 0.75 || far < 1.5 * hw {
                    continue;
                }
                let side = (ax.y.atan2(ax.x) / std::f64::consts::FRAC_PI_2).round() as i32 & 3;
                rows.entry(side).or_default().push(((s.center - cc).dot(ax.perp()), o, s.center + ax * far, ax));
            }
            let mut sides: Vec<i32> = rows.keys().cloned().collect();
            sides.sort();
            for side in sides {
                let mut row = rows[&side].clone();
                row.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                out.push(row.into_iter().map(|(_, o, a, ax)| (o, a, ax)).collect());
            }
        }
        out
    }

    /// Fan-out of fine-pitch parts, planned per row of pins the way a 0.4 mm QFN is fanned out by hand.
    /// A pin whose net has another pad within 3.5 mm stays on its layer (it goes to the passives round
    /// the part) unless its way there must cross a neighbour's; every other pin drops to a via in front
    /// of the row. In row order, each via goes where it costs least: a dogbone right at its pin, which
    /// the neighbours' wires then pass over, or further out, clear of the straight way out of every pin
    /// that stays on its layer; always clear of the vias and stubs already planned, with the stubs in pin
    /// order so that none crosses another. Vias and stubs are laid as the nets' copper; once every row has
    /// them, the pins that stay on their layer route, row by row in nested lanes; the other nets go on
    /// from their vias later. (Short straight stubs for the pins that stay on their layer took the room the power pins'
    /// plane drops need: worse on EXT.)
    fn plan_fanout(&mut self) {
        let mut made = 0;
        let mut wanted = 0;
        for row in self.escape_rows() {
            let (m, w, _, _) = self.fan_row(&row, 0, row.len(), false);
            made += m;
            wanted += w;
        }
        if self.cfg.verbose {
            eprintln!("[{:.1}s] planned fan-out: {} of {} far-going fine-pitch pins dropped to vias", self.elapsed(), made, wanted);
        }
    }

    /// The fan-out plan of the stretch lo..hi of a pin row: lays the vias (and their stubs) of the
    /// stretch's far-going pins that have no copper yet. With `dogbones` a via may also sit right at its
    /// pin between pins that stay on their layer, whose wires then pass over it (planned on the whole
    /// board at once, before the rest routes, that cost other nets their room on EXT; it pays when a sealed
    /// stretch is rebuilt). Returns the vias laid, the vias wanted and the nets of the stretch's pins that
    /// stay on their layer in nesting order.
    fn fan_row(&mut self, row: &[(u32, P2, P2)], lo: usize, hi: usize, dogbones: bool) -> (usize, usize, Vec<u32>, Vec<(u32, u32, f64)>) {
        let mut made = 0;
        let mut wanted = 0;
        // dogbones that wires pass over get their stubs once those wires are in: (pad, via, stub reach)
        let mut pending: Vec<(u32, u32, f64)> = Vec::new();
        let (a0, ax) = (row[0].1, row[0].2);
        let perp = ax.perp();
        let u: Vec<f64> = row.iter().map(|r| (r.1 - a0).dot(perp)).collect();
        // 0: routes otherwise (plane, pair, single pad), 1: stays on its layer, 2: drops to a via
        // where each pin's nearest pad of its net elsewhere lies along the row
        let mut target: Vec<f64> = vec![0.0; row.len()];
        let mut kind: Vec<u8> = row
            .iter()
            .enumerate()
            .map(|(i, &(o, _, _))| {
                let s = &self.w.objs[o as usize].s;
                let net = s.net;
                let paired = self.w.b.partner.get(net as usize).map_or(false, |&x| x != NONE);
                // a pin with a neck-down stub routes on from the stub's end
                if self.is_plane(net) || paired || self.net_pads[net as usize].len() < 2 || self.has_neck(o) {
                    return 0;
                }

                let part = s.name.split('.').next().unwrap_or("");
                let nearest = self.net_pads[net as usize]
                    .iter()
                    .filter(|&&q| q != o && self.w.objs[q as usize].s.name.split('.').next().unwrap_or("") != part)
                    .map(|&q| {
                        let t = &self.w.objs[q as usize].s;
                        (point_hull_dist(&t.pts, s.center) - t.r, t.center)
                    })
                    .fold((f64::INFINITY, s.center), |a, b| if b.0 < a.0 { b } else { a });
                target[i] = (nearest.1 - a0).dot(perp);
                if nearest.0 < 3.5 {
                    1
                } else {
                    2
                }
            })
            .collect();
        // a net kept to its pad's layer (an RF line) never drops to a via
        let locked: Vec<bool> = row.iter().map(|&(o, _, _)| { let s = &self.w.objs[o as usize].s; self.w.b.nets[s.net as usize].allowed & !s.layers == 0 }).collect();
        for i in 0..row.len() {
            if locked[i] && kind[i] == 2 {
                kind[i] = 1;
            }
        }
        // two pins that stay on their layer with their targets the other way round along the row must
        // cross there: the one in most such crossings drops to a via instead, until none is left
        loop {
            let inv = |i: usize, me: &Vec<u8>| (0..row.len()).filter(|&j| j != i && me[j] == 1 && ((u[j] - u[i]) * (target[j] - target[i]) < 0.0) && (target[j] - target[i]).abs() > 0.2).count();
            let worst = (0..row.len()).filter(|&i| kind[i] == 1 && !locked[i]).map(|i| (inv(i, &kind), i)).max();
            match worst {
                Some((n, i)) if n > 0 => kind[i] = 2,
                _ => break,
            }
        }
        let w_of = |me: &Router, o: u32| -> (f64, f64, f64) {
            let net = me.w.objs[o as usize].s.net;
            let (vd, _) = me.w.b.via_size(net);
            let ni = &me.w.b.nets[net as usize];
            let l = (0..me.w.nsig()).find(|&l| me.w.objs[o as usize].s.layers & (1 << l) != 0).unwrap_or(0);
            (vd, ni.wpref[l], me.w.b.classes[ni.class].clear)
        };
        // planned so far: (row index, lateral position, distance out)
        let mut placed: Vec<(usize, f64, f64)> = Vec::new();
        for i in lo..hi {
            // a pin already connected (outside the stretch being rebuilt) keeps its copper
            if kind[i] != 2 || self.w.obj_pieces[row[i].0 as usize].iter().any(|&p| self.w.pieces[p as usize].alive) {
                continue;
            }
            wanted += 1;
            let (o, a, _) = row[i];
            let net = self.w.objs[o as usize].s.net;
            let (vd, w, clr) = w_of(self, o);
            let gap_track = vd / 2.0 + w / 2.0 + clr + 0.03;
            let gap_via = vd + clr + 0.03;
            let mut cands: Vec<(f64, f64, f64)> = Vec::new();
            for &sd in &[0.0, 0.3, 0.6, 0.9, 1.2] {
                let s_out = vd / 2.0 + clr + 0.05 + sd;
                for &d in &[0.0, 0.15, -0.15, 0.3, -0.3, 0.45, -0.45, 0.6, -0.6, 0.8, -0.8] {
                    cands.push((s_out + 1.5 * f64::abs(d), s_out, d));
                }
            }
            // finer dogbone spots: a little sideways, a little further out
            if dogbones {
                for &sd in &[0.0, 0.075] {
                    let s_out = vd / 2.0 + clr + 0.05 + sd;
                    for &d in &[0.05, -0.05, 0.1, -0.1] {
                        cands.push((s_out + 1.5 * f64::abs(d), s_out, d));
                    }
                }
            }
            // neighbours that stay on their layer and must pass over this pin's dogbone, from which side
            let over_left = (0..row.len()).any(|j| kind[j] == 1 && u[j] < u[i] && target[j] > u[i] + 0.1);
            let over_right = (0..row.len()).any(|j| kind[j] == 1 && u[j] > u[i] && target[j] < u[i] - 0.1);
            // a dogbone that wires pass over first tries the spot that leaves them room round it (as by
            // hand: a quarter pitch sideways, away from them, and a little further out than the tightest)
            if dogbones && over_left != over_right {
                cands.insert(0, (0.0, vd / 2.0 + clr + 0.05 + 0.075, if over_left { 0.1 } else { -0.1 }));
            }
            cands.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
            let mut done = false;
            for &(_, s_out, d) in &cands {
                let v = u[i] + d;
                let (lo, hi) = (u[i].min(v), u[i].max(v));
                // the ways out of pins that stay on their layer: a via close in (a dogbone at its pin)
                // only has to leave them room to leave their pads, as their wires then run on over it;
                // further out it keeps off their straight way (vias there blocked their lanes on EXT)
                let dogbone = dogbones && s_out <= vd / 2.0 + clr + 0.05 + 0.15 && d.abs() <= 0.15;
                let blocks_way = (0..row.len()).any(|j| {
                    if kind[j] != 1 {
                        return false;
                    }
                    if dogbone {
                        // a dogbone leans away from the side the wires passing over it come up on
                        if (over_left && !over_right && d <= 0.0) || (over_right && !over_left && d >= 0.0) {
                            return true;
                        }
                        let (f0, f1) = (p2(u[j], 0.0), p2(u[j], 0.1));
                        seg_point_dist(f0, f1, p2(v, s_out)) < gap_track - 0.01 || seg_seg_dist(f0, f1, p2(u[i], 0.0), p2(v, s_out)) < w + clr + 0.02
                    } else {
                        (v - u[j]).abs() < gap_track || (u[j] > lo - (w + clr + 0.03) && u[j] < hi + (w + clr + 0.03))
                    }
                });
                if blocks_way {
                    continue;
                }
                // vias and stubs already planned in this row: apart, in pin order, uncrossed
                let stub = |k: usize, vk: f64, sk: f64| (p2(u[k], 0.0), p2(vk, sk));
                let me_site = p2(v, s_out);
                let (m0, m1) = stub(i, v, s_out);
                let clash = placed.iter().any(|&(k, vk, sk)| {
                    let (k0, k1) = stub(k, vk, sk);
                    let site = p2(vk, sk);
                    site.dist(me_site) < gap_via || (k < i) != (vk < v) || seg_point_dist(k0, k1, me_site) < gap_track || seg_point_dist(m0, m1, site) < gap_track || seg_seg_dist(k0, k1, m0, m1) < w + clr + 0.03
                });
                if clash {
                    continue;
                }
                let pos = a + ax * s_out + perp * d;
                if !self.w.via_fits(net, pos) {
                    continue;
                }
                self.w.begin();
                let vid = match self.w.add_via(net, pos) {
                    Ok(v) => v,
                    Err(_) => {
                        self.w.rollback();
                        continue;
                    }
                };
                let mut viol = Vec::new();
                self.w.check_via_at(net, pos, vid, &mut viol, true);
                if !viol.is_empty() {
                    self.w.rollback();
                    continue;
                }
                self.w.commit();
                // a stub to a dogbone that wires pass over would stand in their way round the via: later
                if dogbone && over_left != over_right {
                    pending.push((o, vid, s_out + d.abs() + 0.6));
                    placed.push((i, v, s_out));
                    made += 1;
                    done = true;
                    break;
                }
                if self.connect_within(net, o, vid, NONE, s_out + d.abs() + 0.6) {
                    placed.push((i, v, s_out));
                    made += 1;
                    done = true;
                    break;
                }
                let _ = self.w.remove_via(vid);
            }
            if !done && std::env::var("TOPOR_FANDBG").is_ok() {
                eprintln!("fan-out: no via site for {} ({})", self.w.objs[o as usize].s.name, self.w.b.nets[net as usize].name);
            }
        }
        // the stretch's pins that stay on their layer in nesting order: those whose pad lies straight
        // ahead first, then of those going one way along the row the one nearest that side first, so
        // that routed in this order each next wire runs round the ones already there (over the
        // dogbones) to its pad
        let mut fl: Vec<usize> = (lo..hi).filter(|&i| kind[i] == 1).collect();
        let key = |i: usize| -> (u8, f64) {
            let off = target[i] - u[i];
            if off.abs() <= 0.2 {
                (0, 0.0)
            } else if off > 0.0 {
                (1, -u[i])
            } else {
                (2, u[i])
            }
        };
        fl.sort_by(|&a, &b| key(a).partial_cmp(&key(b)).unwrap());
        if std::env::var("TOPOR_FANDBG").is_ok() {
            let tag: String = kind.iter().map(|k| match k { 0 => '.', 1 => 'F', _ => 'V' }).collect();
            eprintln!("fan-out row {} .. {}: {}", self.w.objs[row[0].0 as usize].s.name, self.w.objs[row[row.len() - 1].0 as usize].s.name, tag);
        }
        (made, wanted, fl.iter().map(|&i| row[i].0).collect(), pending)
    }

    /// End game for pins of fine-pitch rows still open, freed the way it is done by hand: rip all copper in
    /// front of the pin's stretch of the row (three pins either side, out to the passives round the part),
    /// re-plan the stretch with dogbones allowed (vias right at their pins), route its pins that stay on
    /// their layer in nested lanes, then everything else that was ripped; a few rip-up rounds settle it.
    /// Kept only if the board ends up less incomplete.
    fn re_fan_repairs(&mut self) {
        let rows = self.escape_rows();
        let nn = self.w.b.nets.len() as u32;
        // each stretch rebuild is a trial on its own copy of the board, all at once; the best that pays is kept
        // and the rest try again from there
        loop {
            if self.elapsed() > self.cfg.seconds * 1.5 || self.hurried() {
                return;
            }
            let failing: Vec<u32> = (0..nn).filter(|&x| self.deficit(x) > 0 && !self.is_plane(x)).collect();
            let mut jobs: Vec<(u32, usize, usize)> = Vec::new();
            for &net in &failing {
                let pads = &self.net_pads[net as usize];
                for (ri, row) in rows.iter().enumerate() {
                    if let Some(k) = row.iter().position(|x| pads.contains(&x.0)) {
                        jobs.push((net, ri, k));
                    }
                }
            }
            if !self.adopt_best_trial(&jobs, |r, &(net, ri, k)| r.re_fan_trial(net, &rows[ri], k)) {
                break;
            }
        }
    }

    /// One stretch rebuild for `net` round pin k of `row`, on this copy of the board: its event line.
    fn re_fan_trial(&mut self, net: u32, row: &[(u32, P2, P2)], k: usize) -> Option<String> {
        if self.deficit(net) == 0 {
            return None;
        }
        let (lo, hi) = (k.saturating_sub(3), (k + 4).min(row.len()));
        let before = self.total_deficit();
        // the window in row coordinates: along the stretch and 0.6 mm beyond, out to 1.8 mm, and as
        // far as the nearby pads its pins go to (a lane's target left outside kept its old wires)
        let (a0, ax) = (row[0].1, row[0].2);
        let perp = ax.perp();
        let (mut ul, mut uh, mut sh) = ((row[lo].1 - a0).dot(perp) - 0.6, (row[hi - 1].1 - a0).dot(perp) + 0.6, 1.8f64);
        for i in lo..hi {
            let (o, a, _) = row[i];
            let s0 = &self.w.objs[o as usize].s;
            let part = s0.name.split('.').next().unwrap_or("").to_string();
            for &q in &self.net_pads[s0.net as usize] {
                let t = &self.w.objs[q as usize].s;
                if q == o || t.name.split('.').next().unwrap_or("") == part || t.center.dist(a) > 3.5 {
                    continue;
                }
                let (tu, ts) = ((t.center - a0).dot(perp), (t.center - a0).dot(ax));
                ul = ul.min(tu - 0.6);
                uh = uh.max(tu + 0.6);
                sh = sh.max(ts + 0.4);
            }
        }
        let inside = |q: P2| {
            let d = q - a0;
            let (uu, ss) = (d.dot(perp), d.dot(ax));
            uu > ul && uu < uh && ss > -0.05 && ss < sh
        };
        let mut rip: Vec<u32> = Vec::new();
        for (pi, pc) in self.w.pieces.iter().enumerate() {
            if pc.alive && (pc.poly.iter().any(|&q| inside(q)) || pc.poly.windows(2).any(|sg| inside(sg[0].lerp(sg[1], 0.5)))) {
                rip.push(pi as u32);
            }
        }
        for (oi, o) in self.w.objs.iter().enumerate() {
            if o.alive && o.s.kind == Kind::Via && inside(o.s.center) {
                rip.extend(self.w.obj_pieces[oi].iter().cloned().filter(|&p| self.w.pieces[p as usize].alive));
            }
        }
        rip.sort();
        rip.dedup();
        let mut moved: Vec<u32> = rip.iter().map(|&p| self.w.pieces[p as usize].net).collect();
        moved.push(net);
        moved.sort();
        moved.dedup();
        for p in rip {
            self.rip_piece(p);
        }
        let dbg = std::env::var("TOPOR_FANDBG").is_ok();
        let open = |me: &Router| -> Vec<String> { moved.iter().filter(|&&x| me.deficit(x) > 0).map(|&x| me.w.b.nets[x as usize].name.clone()).collect() };
        let (made, wanted, nested, pending) = self.fan_row(row, lo, hi, true);
        if dbg {
            eprintln!("re-fan {}: moved {:?}; vias {} of {}; nested {:?}", self.w.b.nets[net as usize].name, moved.iter().map(|&x| self.w.b.nets[x as usize].name.clone()).collect::<Vec<_>>(), made, wanted, nested.iter().map(|&o| self.w.objs[o as usize].s.name.clone()).collect::<Vec<_>>());
        }
        // the lanes stay on the pins' layer (left free, a lane took a via wherever that was cheaper)
        for &o in &nested {
            let x = self.w.objs[o as usize].s.net;
            if self.deficit(x) > 0 {
                let ok = self.route_on_layer(x, o);
                if dbg {
                    eprintln!("   lane {} ({}): {}", self.w.objs[o as usize].s.name, self.w.b.nets[x as usize].name, if ok { "on its layer".to_string() } else { format!("failed after {} expansions", self.last_exp) });
                }
                if !ok {
                    self.route_net(x);
                }
            }
        }
        for (o, vid, reach) in pending {
            let x = self.w.objs[o as usize].s.net;
            if !self.connect_within(x, o, vid, NONE, reach) {
                let _ = self.w.remove_via(vid);
            }
        }
        if dbg {
            eprintln!("   after nested lanes, open: {:?}", open(self));
            if let Ok(dir) = std::env::var("TOPOR_FANSVG") {
                let c = row[k].1;
                let svg = crate::svg::render(&self.w, &crate::svg::Opts { mesh: false, box_: Some([c.x - 2.8, c.y - 2.8, c.x + 2.8, c.y + 1.2]), scale: 120.0, marks: vec![], unrouted: vec![] });
                let _ = std::fs::write(format!("{}/refan-{}.svg", dir, self.w.objs[row[k].0 as usize].s.name), svg);
            }
        }
        for &x in &moved {
            if self.deficit(x) > 0 {
                self.reroute(x);
            }
        }
        if dbg {
            eprintln!("   after rerouting the rest, open: {:?}", open(self));
        }
        self.settle(8);
        Some(format!("{} was sealed in at {} by its neighbours' escapes: rebuilt that stretch of the row with dogbone vias and nested lanes ({} nets moved), deficit {} -> {}", self.w.b.nets[net as usize].name, self.w.objs[row[k].0 as usize].s.name, moved.len(), before, self.total_deficit()))
    }

    /// Connect pad o to the rest of its net without leaving the pad's layer (no vias).
    fn route_on_layer(&mut self, net: u32, o: u32) -> bool {
        let comps = self.components(net);
        let (src, targets): (Vec<u32>, Vec<u32>) = match comps.iter().position(|c| c.contains(&o)) {
            Some(i) => (comps[i].clone(), comps.iter().enumerate().filter(|&(j, _)| j != i).flat_map(|(_, c)| c.iter().cloned()).collect()),
            None => return false,
        };
        if targets.is_empty() {
            return true;
        }
        let mut q = self.query(net, src, targets, false);
        q.allowed &= self.w.objs[o as usize].s.layers;
        q.via_cost = 1e6;
        let w = q.width.clone();
        for _ in 0..3 {
            let path = match self.search(&q) {
                Some(p) => p,
                None => return false,
            };
            if path.vias.is_empty() && self.commit(net, &path, &w).is_ok() {
                return true;
            }
            self.learn(&mut q);
        }
        false
    }

    /// A few rip-up rounds over whatever is still open (after a local rebuild), ending on the best board
    /// they reach.
    pub(crate) fn settle(&mut self, rounds: usize) {
        self.settle_except(rounds, &Default::default());
    }

    /// `settle`, leaving out the nets in `skip` (ones that were open already and will not come in now).
    pub(crate) fn settle_except(&mut self, rounds: usize, skip: &std::collections::HashSet<u32>) {
        let nn = self.w.b.nets.len() as u32;
        let mut best: Option<(usize, (World, Vec<(u32, u32)>, Vec<f64>, Vec<u32>))> = None;
        for _ in 0..rounds {
            let d = self.total_deficit();
            if best.as_ref().map_or(true, |b| d < b.0) {
                best = Some((d, (self.w.clone(), self.ties.clone(), self.history.clone(), self.junctions.clone())));
            }
            let failing: Vec<u32> = (0..nn).filter(|&n| !skip.contains(&n) && self.deficit(n) > 0 && self.net_pads[n as usize].len() >= 1).collect();
            if failing.is_empty() || self.elapsed() > self.cfg.seconds * 1.5 || self.stopped() {
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
        if let Some((d, st)) = best {
            if d < self.total_deficit() {
                self.w = st.0;
                self.ties = st.1;
                self.history = st.2;
                self.junctions = st.3;
            }
        }
    }

    /// Nets with pins on fine-pitch parts (16+ pads, neighbours of other nets under 0.7 mm apart), in
    /// fan-out order: pins deep inside a part (inner rows) first, as their lanes out run between the
    /// outer pins; then outer pins round each part. Plane nets are left to the plane drops.
    pub fn fanout_order(&self) -> Vec<u32> {
        let mut parts: HashMap<String, Vec<u32>> = HashMap::default();
        for (i, o) in self.w.objs.iter().enumerate() {
            if o.alive && o.s.kind == Kind::Pad && o.s.net != NONE && !o.s.name.starts_with("junction:") && !crate::world::fixed_copper(&o.s.name) {
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
        self.pair_stations_with(p, n, self.merge_stations)
    }

    /// pair_stations, with stations that fixed copper joins merged or not.
    fn pair_stations_with(&self, p: u32, n: u32, merge: bool) -> Vec<(u32, u32)> {
        let c = |o: u32| self.w.objs[o as usize].s.center;
        // (merging: the vias of the members' suggested copper count too, as a fixed via pair would: the pair then
        // starts where the hand-drawn escapes leave the layer)
        let ends = |x: u32| -> Vec<u32> {
            let mut v = self.net_pads[x as usize].clone();
            if merge {
                let mut vias: Vec<u32> = self
                    .suggested
                    .iter()
                    .map(|&q| &self.w.pieces[q as usize])
                    .filter(|pc| pc.alive && pc.net == x)
                    .flat_map(|pc| [pc.a.obj, pc.b.obj])
                    .filter(|&o| self.w.objs[o as usize].alive && self.w.objs[o as usize].s.kind == Kind::Via)
                    .collect();
                vias.sort();
                vias.dedup();
                v.extend(vias);
            }
            v
        };
        let (ends_p, ends_n) = (ends(p), ends(n));
        let mut cand: Vec<(f64, u32, u32, P2)> = Vec::new();
        let mut cen = p2(0.0, 0.0);
        let mut k = 0.0;
        for &a in &self.net_pads[p as usize] {
            cen = cen + c(a);
            k += 1.0;
        }
        for &a in &ends_p {
            for &b in &ends_n {
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
        // stations the members' fixed copper (hand-drawn tracks and vias, zones) already joins, both members, are one:
        // the pair needs no ribbon between them, only from the copper's tip that faces the other stations
        let (gp, gn) = (self.fixed_groups(p), self.fixed_groups(n));
        let key = |x: &(u32, u32, P2)| (gp.get(&x.0).cloned(), gn.get(&x.1).cloned());
        let mut keys: Vec<(Option<usize>, Option<usize>)> = Vec::new();
        for x in &st {
            if !keys.contains(&key(x)) {
                keys.push(key(x));
            }
        }
        if merge && keys.len() < st.len() {
            let mut merged: Vec<(u32, u32, P2)> = Vec::new();
            for k in &keys {
                let mine: Vec<&(u32, u32, P2)> = st.iter().filter(|x| key(x) == *k).collect();
                let others: Vec<P2> = st.iter().filter(|x| key(x) != *k).map(|x| x.2).collect();
                // the copper's own end before a part's pads it starts from: its via pair (where it leaves the
                // layer), else its pieces
                let rank = |x: &(u32, u32, P2)| {
                    let kinds: Vec<u8> = [x.0, x.1]
                        .iter()
                        .map(|&o| {
                            let s = &self.w.objs[o as usize].s;
                            if s.name.starts_with("via.") || s.kind == Kind::Via {
                                0
                            } else if crate::world::fixed_copper(&s.name) {
                                1
                            } else {
                                2
                            }
                        })
                        .collect();
                    kinds[0].max(kinds[1])
                };
                let best_rank = mine.iter().map(|x| rank(x)).min().unwrap_or(2);
                let pool: Vec<&(u32, u32, P2)> = mine.iter().cloned().filter(|x| rank(x) == best_rank).collect();
                let pick = if others.is_empty() {
                    pool[0]
                } else {
                    let c = others.iter().fold(p2(0., 0.), |a, &b| a + b) * (1.0 / others.len() as f64);
                    pool.iter().min_by(|a, b| a.2.dist(c).partial_cmp(&b.2.dist(c)).unwrap()).unwrap()
                };
                merged.push(*pick);
            }
            st = merged;
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

    /// The groups a net's pads and fixed copper form through their ties alone (not routed wires): which of them
    /// hand-drawn copper or zones already join. Object -> group.
    fn fixed_groups(&self, net: u32) -> HashMap<u32, usize> {
        let objs = &self.net_pads[net as usize];
        let mut index: HashMap<u32, usize> = HashMap::default();
        for (i, &o) in objs.iter().enumerate() {
            index.entry(o).or_insert(i);
        }
        let mut par: Vec<usize> = (0..objs.len()).collect();
        let join = |x: u32, y: u32, par: &mut Vec<usize>| {
            if let (Some(&i), Some(&j)) = (index.get(&x), index.get(&y)) {
                let (a, b) = (find(par, i), find(par, j));
                if a != b {
                    par[b] = a;
                }
            }
        };
        for &(x, y) in &self.ties {
            join(x, y, &mut par);
        }
        // suggested copper counts as drawn: chains of its pieces through vias and junctions join what they reach
        let mut sug: Vec<(u32, u32)> = self.suggested.iter().filter(|&&q| self.w.pieces[q as usize].alive && self.w.pieces[q as usize].net == net).map(|&q| (self.w.pieces[q as usize].a.obj, self.w.pieces[q as usize].b.obj)).collect();
        sug.sort();
        // a via or junction between two pieces is not a station object itself: join through it
        let mut through: HashMap<u32, Vec<u32>> = HashMap::default();
        for &(a, b) in &sug {
            through.entry(a).or_default().push(b);
            through.entry(b).or_default().push(a);
        }
        for (_, ends) in through.iter() {
            for w in ends.windows(2) {
                join(w[0], w[1], &mut par);
            }
        }
        for &(a, b) in &sug {
            join(a, b, &mut par);
        }
        objs.iter().enumerate().map(|(i, &o)| (o, find(&mut par, i))).collect()
    }

    /// Whether objects a and b of a net are joined by its copper (wires, and ties through fixed copper).
    fn joined(&self, net: u32, a: u32, b: u32) -> bool {
        self.components(net).iter().any(|c| c.contains(&a) && c.contains(&b))
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
            // both members joined station to station already by hand copper given as suggestions: a ribbon or
            // member wires there would only double it (where fixed copper joins them, the segment is routed as
            // before)
            let sug = self.suggested_nets.contains(&p) || self.suggested_nets.contains(&n);
            if sug && self.joined(p, st[k].0, st[k + 1].0) && self.joined(n, st[k].1, st[k + 1].1) {
                segs.push((labels[k].clone(), labels[k + 1].clone(), "joined already (suggested copper)".into()));
                continue;
            }
            // a dogbone pair: no ribbon (its via pair would sit wherever the ribbon changes layer); each member
            // drops at its own pads' dogbones, the second alongside the first
            let dogbone = self.w.b.dogbone.get(p as usize).copied().unwrap_or(false) && self.w.b.dogbone.get(n as usize).copied().unwrap_or(false);
            let res = if dogbone { Err("dogbone vias".to_string()) } else { self.route_segment(r, p, n, st[k], st[k + 1]) };
            match res {
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
                if !(self.connect_within(net, j, pad, NONE, cap) || self.connect_within(net, pad, j, NONE, cap) || self.stub_after_retreat(net, j, pad, cap)) {
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

    /// A member's stub from its lane end `j` to `pad` that did not fit: the lane ends earlier instead (0.3, 0.6 or
    /// 1 mm back along it, this member only), which leaves its stub room to turn. Kept if the stub then joins.
    fn stub_after_retreat(&mut self, net: u32, j: u32, pad: u32, cap: f64) -> bool {
        if self.legacy_pairs {
            return false;
        }
        for d in [0.3, 0.6, 1.0] {
            let snap = (self.w.clone(), self.junctions.clone());
            if let Some(nj) = self.retreat_lane_end(net, j, d) {
                if self.connect_within(net, nj, pad, NONE, cap + d) || self.connect_within(net, pad, nj, NONE, cap + d) {
                    self.junctions.retain(|&x| x != j);
                    self.junctions.push(nj);
                    return true;
                }
            }
            self.w = snap.0;
            self.junctions = snap.1;
        }
        false
    }

    /// End the member lane of `net` that ends at junction `j` a length `d` earlier, at a new junction (returned).
    /// The board may be left half changed on failure: the caller restores it.
    fn retreat_lane_end(&mut self, net: u32, j: u32, d: f64) -> Option<u32> {
        let pid = self.w.obj_pieces[j as usize].iter().cloned().find(|&q| self.w.pieces[q as usize].alive && self.w.pieces[q as usize].net == net)?;
        let pc = self.w.pieces[pid as usize].clone();
        let l = pc.layer as usize;
        // the lane from its other end to j
        let (other, mut poly) = if pc.b.obj == j { (pc.a, pc.poly.clone()) } else { (pc.b, pc.poly.iter().rev().cloned().collect::<Vec<_>>()) };
        let total: f64 = poly.windows(2).map(|x| x[0].dist(x[1])).sum();
        if total < d + 0.2 {
            return None;
        }
        let mut left = d;
        while poly.len() >= 2 {
            let k = poly.len();
            let seg = poly[k - 1].dist(poly[k - 2]);
            if seg > left + 1e-9 {
                poly[k - 1] = poly[k - 1].lerp(poly[k - 2], left / seg);
                break;
            }
            left -= seg;
            poly.pop();
        }
        let cut = *poly.last()?;
        let r = self.w.objs[j as usize].s.r;
        self.w.remove_piece(pid);
        self.w.remove_via(j).ok()?;
        let nj = self.w.add_junction(net, cut, r, l).ok()?;
        let tb = Term { obj: nj, vtx: self.w.objs[nj as usize].lv[l] };
        if other.vtx != NONE {
            poly[0] = self.w.vpos(l, other.vtx);
        }
        *poly.last_mut().unwrap() = self.w.vpos(l, tb.vtx);
        let cx = self.w.thread_poly(l, &poly, other, tb)?;
        let np = self.w.add_piece(net, l, pc.width, other, tb, &cx);
        self.w.seed_poly(np, poly);
        self.w.realize_dirty();
        Some(nj)
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
        self.note_rebuild(p, n, "rebuilt");
    }

    /// Add how a rip-up rebuild of pair (p, n) ended to its last displacement note.
    fn note_rebuild(&mut self, p: u32, n: u32, how: &str) {
        let cf = self.pair_measure_routed(p, n).map_or(0.0, |x| x.2);
        let done = self.deficit(p) + self.deficit(n) == 0;
        if let Some(m) = self.pair_diag.entry(p).or_default().moved.last_mut() {
            m.push_str(&format!("; {} {}, {:.0}% coupled", how, if done { "complete" } else { "incomplete" }, cf * 100.0));
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
        let mut best: Option<(f64, (World, Vec<(u32, u32)>, Vec<f64>, Vec<u32>), HashMap<u32, PairDiag>)> = None;
        let orders: Vec<&str> = match std::env::var("TOPOR_PAIRORDER").as_deref() {
            Ok("long") => vec!["long"],
            Ok("short") => vec!["short"],
            _ => vec!["short", "long"],
        };
        // (the two orders are independent: each runs on its own copy, at once)
        let runs: Vec<(f64, Router)> = std::thread::scope(|sc| {
            let hs: Vec<_> = orders
                .iter()
                .map(|order| {
                    let mut r = self.clone();
                    let mut ps = pairs.clone();
                    sc.spawn(move || {
                        if *order == "long" {
                            ps.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
                        } else {
                            ps.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                        }
                        for &(_, p, n) in &ps {
                            r.route_pair(p, n);
                            r.route_net(p);
                            r.route_net(n);
                        }
                        let score: f64 = ps.iter().map(|&(_, p, n)| if r.deficit(p) + r.deficit(n) > 0 { -1.0 } else { 0.0 }).sum::<f64>() + r.pair_report().iter().map(|x| x.4).sum::<f64>();
                        if std::env::var("TOPOR_PAIRDBG").is_ok() {
                            eprintln!("pair order {}: score {:.2}", order, score);
                        }
                        (score, r)
                    })
                })
                .collect();
            hs.into_iter().map(|h| h.join().expect("pair order")).collect()
        });
        for (score, r) in runs {
            if best.as_ref().map_or(true, |b| score > b.0) {
                best = Some((score, (r.w, r.ties, r.history, r.junctions), r.pair_diag));
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
    /// Measure each pair's coupled length again (for the rip-up costs).
    pub fn update_pair_coupling(&mut self) {
        let mut v = vec![0.0; self.w.b.nets.len()];
        for (a, b, lp, ln, cf) in self.pair_report() {
            let (Some(&p), Some(&n)) = (self.w.b.net_index.get(&a), self.w.b.net_index.get(&b)) else { continue };
            let fixed = || self.w.objs.iter().any(|o| o.alive && (o.s.net == p || o.s.net == n) && crate::world::fixed_copper(&o.s.name));
            let c = if lp.is_finite() && ln.is_finite() {
                cf * lp.min(ln)
            } else if !fixed() {
                0.0
            } else {
                // no chain of the router's own wires between the pair's ends (hand copper kept as drawn in
                // between): the coupled length of n's own wires beside p's copper, fixed copper too, which is
                // what a rip can take away
                let steps: Vec<(Vec<P2>, usize, f64)> = self.w.pieces.iter().filter(|x| x.alive && x.net == n).map(|q| (q.poly.clone(), q.layer as usize, q.width)).collect();
                let len: f64 = steps.iter().map(|s| s.0.windows(2).map(|t| t[0].dist(t[1])).sum::<f64>()).sum();
                self.coupled_fraction(p, &steps, true) * len
            };
            v[p as usize] = c;
            v[n as usize] = c;
            if std::env::var("TOPOR_COUPLEDBG").is_ok() {
                eprintln!("pair coupling {}: lp {:.2} ln {:.2} cf {:.2} -> {:.2} mm", a, lp, ln, cf, c);
            }
        }
        self.pair_coupled = v;
    }

    /// Per pair (p, n): member path lengths between the pair's ends and the coupled fraction, on the router's own
    /// wires: what routing decisions go by (pair order, attempts, rip-up costs, rebuilds).
    pub fn pair_report(&self) -> Vec<(String, String, f64, f64, f64)> {
        self.pairs_with(|me, p, n| me.pair_measure_routed(p, n))
    }

    /// The same as reported: fixed copper counts too (a pair kept as drawn is measured), between the pads.
    pub fn pair_report_full(&self) -> Vec<(String, String, f64, f64, f64)> {
        self.pairs_with(|me, p, n| me.pair_measure(p, n))
    }

    fn pairs_with(&self, f: impl Fn(&Router, u32, u32) -> Option<(f64, f64, f64)>) -> Vec<(String, String, f64, f64, f64)> {
        let mut out = Vec::new();
        for p in 0..self.w.b.nets.len() as u32 {
            let n = self.w.b.partner.get(p as usize).cloned().unwrap_or(NONE);
            if n == NONE || n < p {
                continue;
            }
            if let Some((lp, ln, cf)) = f(self, p, n) {
                out.push((self.w.b.nets[p as usize].name.clone(), self.w.b.nets[n as usize].name.clone(), lp, ln, cf));
            }
        }
        out
    }

    /// One pair on the router's wires: (p's length, n's length, coupled fraction: the length of n's pad-to-pad path
    /// within the coupling distance of p's wires on the same layer, sampled every 50 um).
    pub fn pair_measure_routed(&self, p: u32, n: u32) -> Option<(f64, f64, f64)> {
        let ((pa, na), (pb, nb)) = self.pair_ends(p, n)?;
        let lp = self.net_path_len(p, pa, pb).unwrap_or(f64::NAN);
        let ln = self.net_path_len(n, na, nb).unwrap_or(f64::NAN);
        let chain: Vec<u32> = self.piece_chain(n, na, nb).map(|c| c.iter().map(|x| x.0).collect()).unwrap_or_default();
        let steps: Vec<(Vec<P2>, usize, f64)> = chain.iter().map(|&pi| { let q = &self.w.pieces[pi as usize]; (q.poly.clone(), q.layer as usize, q.width) }).collect();
        Some((lp, ln, self.coupled_fraction(p, &steps, false)))
    }

    /// One pair as reported: lengths between its outermost pads (stations not merged through fixed copper) and the
    /// coupled fraction, both through fixed copper too (the lines between the centres of tied pieces).
    pub fn pair_measure(&self, p: u32, n: u32) -> Option<(f64, f64, f64)> {
        let st = self.pair_stations_with(p, n, false);
        if st.len() < 2 {
            return None;
        }
        let ((pa, na), (pb, nb)) = (st[0], st[st.len() - 1]);
        let path_p = self.report_path(p, pa, pb);
        let path_n = self.report_path(n, na, nb);
        let len = |x: &Option<Vec<(Vec<P2>, usize, f64)>>| x.as_ref().map_or(f64::NAN, |v| v.iter().map(|s| s.0.windows(2).map(|t| t[0].dist(t[1])).sum::<f64>()).sum());
        if std::env::var("TOPOR_MEASUREDBG").is_ok() {
            let steps: Vec<String> = path_n.iter().flatten().map(|s| format!("L{} w{:.2} {:.1}mm {}pts", s.1, s.2, s.0.windows(2).map(|t| t[0].dist(t[1])).sum::<f64>(), s.0.len())).collect();
            eprintln!("measure {}: routed {:?}; full ends {}/{} {}/{}; n path: {}", self.w.b.nets[p as usize].name, self.pair_measure_routed(p, n), self.w.objs[pa as usize].s.name, self.w.objs[na as usize].s.name, self.w.objs[pb as usize].s.name, self.w.objs[nb as usize].s.name, steps.join(", "));
        }
        Some((len(&path_p), len(&path_n), path_n.as_ref().map_or(0.0, |st| self.coupled_fraction(p, st, true))))
    }

    /// The part of the path `steps` (polyline, layer, width) of p's partner that runs beside p's copper on the same
    /// layer (its wires, and with `fixed` its fixed copper ties) at the coupled pitch (the two half widths and the
    /// pair gap, within 0.05 mm) and parallel (within 15 degrees), sampled every 50 um.
    fn coupled_fraction(&self, p: u32, steps: &[(Vec<P2>, usize, f64)], fixed: bool) -> f64 {
        // p's copper by layer, indexed once
        let mut pcu: Vec<Vec<(P2, P2, f64)>> = vec![Vec::new(); self.w.nsig()];
        for q in self.w.pieces.iter().filter(|x| x.alive && x.net == p) {
            for t in q.poly.windows(2) {
                pcu[q.layer as usize].push((t[0], t[1], q.width));
            }
        }
        if fixed {
            for (x, y, l, w) in self.tie_links(p) {
                pcu[l].push((self.w.objs[x as usize].s.center, self.w.objs[y as usize].s.center, w));
            }
        }
        let gap = self.w.b.pair_gap[p as usize];
        let (mut tot, mut cpl) = (0.0, 0.0);
        let sin_max = 15f64.to_radians().sin();
        for (poly, l, w) in steps {
            for sgm in poly.windows(2) {
                let len = sgm[0].dist(sgm[1]);
                if len < 1e-9 {
                    continue;
                }
                let u = (sgm[1] - sgm[0]) * (1.0 / len);
                let k = (len / 0.05).ceil().max(1.0) as usize;
                for i in 0..k {
                    let m = sgm[0].lerp(sgm[1], (i as f64 + 0.5) / k as f64);
                    tot += len / k as f64;
                    let beside = pcu[*l].iter().any(|&(a, b, qw)| {
                        let lv = a.dist(b);
                        let parallel = lv < 1e-9 || ((b - a) * (1.0 / lv)).cross(u).abs() <= sin_max;
                        parallel && (seg_point_dist(a, b, m) - ((qw + w) / 2.0 + gap)).abs() <= 0.05
                    });
                    if beside {
                        cpl += len / k as f64;
                    }
                }
            }
        }
        if tot > 0.0 { cpl / tot } else { 0.0 }
    }

    /// A net's fixed copper ties (a line between the tied objects' centres): (from, to, layer, width).
    fn tie_links(&self, net: u32) -> Vec<(u32, u32, usize, f64)> {
        let width = |s: &crate::model::ObjSpec| -> f64 {
            if s.pts.len() < 3 {
                return 2.0 * s.r;
            }
            let m = s.pts.len();
            let inner = (0..m).map(|i| seg_point_dist(s.pts[i], s.pts[(i + 1) % m], s.center)).fold(f64::INFINITY, f64::min);
            2.0 * (inner + s.r)
        };
        let mut out = Vec::new();
        for &(x, y) in &self.ties {
            let (sx, sy) = (&self.w.objs[x as usize].s, &self.w.objs[y as usize].s);
            if sx.net != net || sy.net != net || sx.layers & sy.layers == 0 {
                continue;
            }
            let w = width(sx).min(width(sy));
            let mut ls = sx.layers & sy.layers;
            while ls != 0 {
                let l = ls.trailing_zeros() as usize;
                ls &= ls - 1;
                if l < self.w.nsig() {
                    out.push((x, y, l, w));
                }
            }
        }
        out
    }

    /// The shortest path of a net between two of its objects through its wires and its fixed copper ties, as steps
    /// (polyline, layer, width).
    fn report_path(&self, net: u32, a: u32, b: u32) -> Option<Vec<(Vec<P2>, usize, f64)>> {
        let objs = self.net_objs(net);
        let ix: HashMap<u32, usize> = objs.iter().enumerate().map(|(i, &o)| (o, i)).collect();
        let mut adj: Vec<Vec<(usize, Vec<P2>, usize, f64)>> = vec![Vec::new(); objs.len()];
        for (i, &o) in objs.iter().enumerate() {
            for &pc in &self.w.obj_pieces[o as usize] {
                let q = &self.w.pieces[pc as usize];
                if !q.alive || q.net != net || q.poly.len() < 2 {
                    continue;
                }
                let other = if q.a.obj == o { q.b.obj } else { q.a.obj };
                if let Some(&j) = ix.get(&other) {
                    let mut poly = q.poly.clone();
                    if q.a.obj != o {
                        poly.reverse();
                    }
                    adj[i].push((j, poly, q.layer as usize, q.width));
                }
            }
        }
        for (x, y, l, w) in self.tie_links(net) {
            if let (Some(&i), Some(&j)) = (ix.get(&x), ix.get(&y)) {
                let (cx, cy) = (self.w.objs[x as usize].s.center, self.w.objs[y as usize].s.center);
                adj[i].push((j, vec![cx, cy], l, w));
                adj[j].push((i, vec![cy, cx], l, w));
            }
        }
        let (ia, ib) = (*ix.get(&a)?, *ix.get(&b)?);
        let mut dist = vec![f64::INFINITY; objs.len()];
        let mut prev: Vec<Option<(usize, usize)>> = vec![None; objs.len()];
        let mut done = vec![false; objs.len()];
        dist[ia] = 0.0;
        loop {
            let u = (0..objs.len()).filter(|&i| !done[i] && dist[i].is_finite()).min_by(|&x, &y| dist[x].partial_cmp(&dist[y]).unwrap())?;
            if u == ib {
                break;
            }
            done[u] = true;
            for (k, (j, poly, _, _)) in adj[u].iter().enumerate() {
                let d = dist[u] + poly.windows(2).map(|t| t[0].dist(t[1])).sum::<f64>();
                if d < dist[*j] {
                    dist[*j] = d;
                    prev[*j] = Some((u, k));
                }
            }
        }
        let mut steps = Vec::new();
        let mut cur = ib;
        while cur != ia {
            let (u, k) = prev[cur]?;
            let (_, poly, l, w) = &adj[u][k];
            steps.push((poly.clone(), *l, *w));
            cur = u;
        }
        steps.reverse();
        Some(steps)
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
            // (TOPOR_LEGACYPAIRS: the pairs as routed before, to compare)
            if std::env::var("TOPOR_LEGACYPAIRS").is_ok() {
                self.merge_stations = false;
                self.legacy_pairs = true;
                self.w.b.ribbon_touch = false;
            }
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
        // the pairs are routed both as they are now (stations fixed copper joins merged, ribbons from the pad tips,
        // member stubs that may end the lane earlier, rip-up charged for lost coupling) and as before, at once;
        // the routing attempts then run on both, the new ones first, and the board keeps whichever comes out more
        // complete: the coupling never costs a connection (TOPOR_NOSPLITPAIRS: the new ones alone)
        // likewise, on a board with high-speed nets, the board is also routed without the plane gap rule: the rule
        // never costs a connection either (TOPOR_GAP_ONLY: with the rule alone)
        let mut split_pairs: Option<Router> = None;
        let mut no_gap: Option<Router> = None;
        let mut split_no_gap: Option<Router> = None;
        {
            let _t = Instant::now();
            let pairs = !self.cfg.fast && !self.cfg.single_attempt && std::env::var("TOPOR_NOSPLITPAIRS").is_err() && (0..nn).any(|p| {
                let n = self.w.b.partner.get(p as usize).cloned().unwrap_or(NONE);
                n != NONE && p < n && self.pair_stations_with(p, n, false).len() >= 2
            });
            let gap_alt = !self.cfg.fast && !self.cfg.single_attempt && self.gap_rule && crate::gaps::active(&self.w) && std::env::var("TOPOR_GAP_ONLY").is_err();
            let mut other = if pairs {
                let mut o = self.clone();
                o.merge_stations = false;
                o.legacy_pairs = true;
                o.w.b.ribbon_touch = false;
                o.cfg.verbose = false;
                Some(o)
            } else {
                None
            };
            let off = |r: &Router| {
                let mut o = r.clone();
                o.gap_rule = false;
                o.cfg.verbose = false;
                o
            };
            let mut ng = if gap_alt { Some(off(self)) } else { None };
            let mut sng = if gap_alt { other.as_ref().map(off) } else { None };
            std::thread::scope(|sc| {
                let spawn = |r: Option<Router>| r.map(|mut o| sc.spawn(move || {
                    o.route_pairs();
                    o
                }));
                let (h1, h2, h3) = (spawn(other.take()), spawn(ng.take()), spawn(sng.take()));
                self.route_pairs();
                split_pairs = h1.map(|h| h.join().expect("pairs routed as before"));
                no_gap = h2.map(|h| h.join().expect("pairs routed without the gap rule"));
                split_no_gap = h3.map(|h| h.join().expect("pairs routed as before, without the gap rule"));
            });
            self.update_pair_coupling();
            for o in split_pairs.iter_mut().chain(no_gap.iter_mut()).chain(split_no_gap.iter_mut()) {
                o.update_pair_coupling();
            }
            self.prof_add("pairs", _t);
        }
        if self.cfg.verbose {
            eprintln!("[{:.1}s] pairs done", self.elapsed());
        }
        self.trace_pairs("pairs");
        // the rest of the board; routing is order-sensitive where fine-pitch parts crowd, so other
        // orders run at once (the fine-pitch parts fanned out first, or given a planned fan-out first),
        // the best one kept. An attempt that
        // completes the board calls off the ones after it; the earlier ones run on, so the first complete
        // attempt in this order always wins.
        let t_attempts = Instant::now();
        let has_fanout = !self.fanout_order().is_empty();
        let forced: Option<usize> = std::env::var("TOPOR_ATTEMPT").ok().and_then(|x| x.parse().ok());
        // (long nets first, with or without the fan-out, did worse on every test board: TOPOR_ATTEMPT=3, 4)
        let all = vec![(false, false, false), (true, false, false), (true, false, true), (false, true, false), (true, true, false)];
        let variants: Vec<(bool, bool, bool)> = if forced.is_some() { all } else if has_fanout { all[..3].to_vec() } else { vec![all[0]] };
        if let Some(k) = forced {
            let (fan, long, esc) = variants[k.min(variants.len() - 1)];
            self.route_rest(fan, long, esc);
        } else if std::env::var("TOPOR_FANOUT").is_ok() {
            self.route_rest(true, false, false);
        } else if (self.cfg.single_attempt || variants.len() == 1) && split_pairs.is_none() && no_gap.is_none() {
            self.route_rest(false, false, false);
        } else {
            // the orders on these pairs, then (if they were routed so too) on the pairs routed as before, and on both
            // routed without the gap rule (0: these, 1: as before, 2 and 3: the same without the gap rule)
            let mut atts: Vec<(u8, (bool, bool, bool))> = variants.iter().map(|&v| (0, v)).collect();
            for (k, r) in [(1u8, &split_pairs), (2, &no_gap), (3, &split_no_gap)] {
                if r.is_some() {
                    atts.extend(variants.iter().map(|&v| (k, v)));
                }
            }
            let stops: Vec<std::sync::Arc<std::sync::atomic::AtomicBool>> = atts.iter().map(|_| std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false))).collect();
            let mut others: Vec<Router> = (1..atts.len())
                .map(|k| {
                    let mut o = match atts[k].0 {
                        1 => split_pairs.as_ref().unwrap().clone(),
                        2 => no_gap.as_ref().unwrap().clone(),
                        3 => split_no_gap.as_ref().unwrap().clone(),
                        _ => self.clone(),
                    };
                    o.cfg.verbose = false;
                    o.stop = Some(stops[k].clone());
                    o
                })
                .collect();
            // a complete board calls off the later attempts on the same pairs (the other pairs still finish, to
            // compare their coupling)
            let nv = variants.len();
            let group_end = |k: usize| (k / nv + 1) * nv;
            let call_off_after = |k: usize, stops: &[std::sync::Arc<std::sync::atomic::AtomicBool>]| {
                for s in &stops[k + 1..group_end(k)] {
                    s.store(true, std::sync::atomic::Ordering::Relaxed);
                }
            };
            // a race: once the first attempt is done, the others get as long again as it took (at least 5 s;
            // TOPOR_GRACE: another fraction of it), then wrap up with what they have, which still counts. Once one
            // is done complete, the others can only win by coupling more, complete too: they get half as long
            // again, then those still running are called off (the first, on this thread, wraps up)
            let hurry = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
            type Fin = (Vec<bool>, Option<Instant>, Option<Instant>);
            let finished = std::sync::Arc::new(std::sync::Mutex::new((vec![false; atts.len()], None::<Instant>, None::<Instant>)));
            let grace_frac: f64 = std::env::var("TOPOR_GRACE").ok().and_then(|x| x.parse().ok()).unwrap_or(1.0);
            let finish = |f: &std::sync::Mutex<Fin>, k: usize, complete: bool| {
                let mut g = f.lock().unwrap();
                g.0[k] = true;
                g.1.get_or_insert_with(Instant::now);
                if complete {
                    g.2.get_or_insert_with(Instant::now);
                }
            };
            for o in others.iter_mut() {
                o.hurry = Some(hurry.clone());
            }
            self.hurry = Some(hurry.clone());
            let mut at0 = 0.0;
            // (whether an attempt was called off is read when it finishes: a later call-off of the attempts after a
            // complete one does not undo one that has finished already)
            let done: Vec<(usize, Router, f64, bool)> = std::thread::scope(|sc| {
                let hs: Vec<_> = others
                    .drain(..)
                    .enumerate()
                    .map(|(i, mut o)| {
                        let (fan, long, esc) = atts[i + 1].1;
                        let stops = &stops;
                        let finished = &finished;
                        sc.spawn(move || {
                            o.route_rest(fan, long, esc);
                            let off = o.stopped();
                            let complete = !off && o.total_deficit() == 0;
                            if complete {
                                call_off_after(i + 1, stops);
                            }
                            finish(finished, i + 1, complete);
                            let at = o.elapsed();
                            (o, at, off)
                        })
                    })
                    .collect();
                {
                    let (hurry, finished, stops) = (hurry.clone(), finished.clone(), &stops);
                    sc.spawn(move || loop {
                        std::thread::sleep(std::time::Duration::from_millis(100));
                        let g = finished.lock().unwrap();
                        if g.0.iter().all(|&x| x) {
                            break;
                        }
                        let late = |t: Option<Instant>, frac: f64| t.map_or(false, |t| t.elapsed().as_secs_f64() > (frac * grace_frac * t.duration_since(t_attempts).as_secs_f64()).max(5.0));
                        if late(g.2, 0.5) {
                            for k in 1..g.0.len() {
                                if !g.0[k] {
                                    stops[k].store(true, std::sync::atomic::Ordering::Relaxed);
                                }
                            }
                        }
                        if late(g.1, 1.0) || late(g.2, 0.5) {
                            hurry.store(true, std::sync::atomic::Ordering::Relaxed);
                            break;
                        }
                    });
                }
                self.route_rest(false, false, false);
                let complete = !self.stopped() && self.total_deficit() == 0;
                if complete {
                    call_off_after(0, &stops);
                }
                finish(&finished, 0, complete);
                at0 = self.elapsed();
                hs.into_iter().enumerate().map(|(i, h)| { let (o, at, off) = h.join().expect("routing attempt"); (i + 1, o, at, off) }).collect()
            });
            let hurried = hurry.load(std::sync::atomic::Ordering::Relaxed);
            self.hurry = None;
            let name = |k: usize| {
                let order = match atts[k].1 {
                    (true, _, true) => "fine-pitch parts given a planned fan-out first",
                    (false, false, _) => "plain order",
                    (true, false, _) => "fine-pitch parts fanned out first",
                    (false, true, _) => "long nets first",
                    (true, true, _) => "fan-out first, then long nets first",
                };
                match atts[k].0 {
                    1 => format!("{} (pairs as before)", order),
                    2 => format!("{} (gap rule off)", order),
                    3 => format!("{} (pairs as before, gap rule off)", order),
                    _ => order.to_string(),
                }
            };
            // the most complete board wins; between equally complete ones, the one with fewer high-speed wires over
            // slots past their limit (boards with such nets), then the one with the most coupled length less the
            // pairs' skews (a member's detour costs what it is longer than its partner)
            let coupled_len = |r: &Router| -> f64 { r.pair_report().iter().map(|x| if x.2.is_finite() && x.3.is_finite() { x.4 * x.2.min(x.3) - (x.2 - x.3).abs() } else { 0.0 }).sum() };
            let gaps_on = crate::gaps::active(&self.w);
            let gap_note = |g: usize| if gaps_on { format!(", {} past the gap limit", g) } else { String::new() };
            let first = self.total_deficit();
            let first_c = coupled_len(self);
            let first_g = self.gap_violations();
            let mut best: Option<(usize, usize, usize, f64)> = None;
            let mut summary = vec![format!("{}: deficit {}{}, {:.0} mm coupled less skew (done at {:.1}s)", name(0), first, gap_note(first_g), first_c, at0)];
            for (k, o, at, off) in &done {
                if *off {
                    summary.push(format!("{}: called off (at {:.1}s)", name(*k), at));
                    continue;
                }
                let d = o.total_deficit();
                let c = coupled_len(o);
                let g = o.gap_violations();
                summary.push(format!("{}: deficit {}{}, {:.0} mm coupled less skew (done at {:.1}s)", name(*k), d, gap_note(g), c, at));
                let (bd, bg, bc) = best.map_or((first, first_g, first_c), |b| (b.1, b.2, b.3));
                if d < bd || (d == bd && (g < bg || (g == bg && c > bc + 1e-6))) {
                    best = Some((*k, d, g, c));
                }
            }
            if self.cfg.verbose {
                eprintln!("[{:.1}s] {}{}", self.elapsed(), summary.join("; "), if hurried { " (the slower ones wrapped up after their grace)" } else { "" });
            }
            if let Some((k, d, _, _)) = best {
                let mut o = done.into_iter().find(|x| x.0 == k).unwrap().1;
                o.cfg.verbose = self.cfg.verbose;
                o.stop = None;
                o.hurry = None;
                o.t0 = self.t0;
                o.events.push(format!("the board completed better with the {} (deficit {} instead of {} in plain order)", name(k), d, first));
                *self = o;
            }
        }
        self.prof_add("routing attempts (in parallel; their phases below)", t_attempts);
        if self.cfg.fast {
            return;
        }
        { let _t = Instant::now(); self.improve_pairs(); self.prof_add("improve pairs", _t); }
        let d0 = self.total_deficit();
        { let _t = Instant::now(); self.cleanup(); self.prof_add("clean-up", _t); }
        if self.cfg.verbose {
            eprintln!("[{:.1}s] after the pair improvements: deficit {}, after clean-up {}", self.elapsed(), d0, self.total_deficit());
        }
        // the spring relaxation: every free wire and via pulls tight together (TOPOR_NOSPRING: off)
        if std::env::var("TOPOR_NOSPRING").is_err() {
            { let _t = Instant::now(); self.spring(); self.prof_add("spring", _t); }
            // (experiment: TOPOR_SPRING2 relaxes twice, the control for teleportation's second pass)
            if std::env::var("TOPOR_SPRING2").is_ok() {
                self.spring();
            }
            if let Some(tries) = std::env::var("TOPOR_TELEPORT").ok().and_then(|s| s.parse::<usize>().ok()) {
                let (tried, kept, saved) = self.teleport_vias(tries, 0x9e37_79b9_7f4a_7c15);
                if self.cfg.verbose {
                    eprintln!("[{:.1}s] via teleportation: {} vias tried, {} jumps kept, wires {:.1} mm shorter", self.elapsed(), tried, kept, saved);
                }
                if kept > 0 {
                    self.events.push(format!("via teleportation: {} jumps kept of {} vias tried, wires {:.1} mm shorter", kept, tried, saved));
                    self.spring();
                }
            }
        }
        // via bridges the relaxed board no longer needs; then the board relaxes again around the change
        let _t = Instant::now();
        let (bridges, removed) = self.remove_bridges();
        self.prof_add("via bridges", _t);
        if removed > 0 {
            self.events.push(format!("{} of {} via bridges removed (a net hopping to another layer and straight back)", removed, bridges));
            if self.cfg.verbose {
                eprintln!("[{:.1}s] via bridges: {} of {} removed", self.elapsed(), removed, bridges);
            }
            // (one solve: the board has relaxed once already, and a second try's gain is a few mm; the retries
            // holding the nets it could not join up were dropped every time on the test boards)
            if std::env::var("TOPOR_NOSPRING").is_err() {
                { let _t = Instant::now(); self.spring_tries(1); self.prof_add("spring", _t); }
            }
        } else if self.cfg.verbose && bridges > 0 {
            eprintln!("[{:.1}s] via bridges: none of {} removable", self.elapsed(), bridges);
        }
        self.trace_pairs("cleanup");
        self.match_pairs(0.1);
    }

    fn stopped(&self) -> bool {
        self.stop.as_ref().map_or(false, |s| s.load(std::sync::atomic::Ordering::Relaxed))
    }

    /// Called off, or asked to wrap up (see `hurry`).
    fn hurried(&self) -> bool {
        self.stopped() || self.hurry.as_ref().map_or(false, |s| s.load(std::sync::atomic::Ordering::Relaxed))
    }

    /// Planes, signals, rip-up and the end-game repairs, after the pairs.
    fn route_rest(&mut self, fanout: bool, long_first: bool, escapes: bool) {
        let nn = self.w.b.nets.len() as u32;
        // RF frontends first: kept to one layer, short, beside their radios; they take the straight lines
        let rf: Vec<u32> = self.w.b.rf_nets.iter().map(|x| x.0).collect();
        let _t = Instant::now();
        for &net in &rf {
            self.route_net(net);
        }
        self.prof_add("RF frontends", _t);
        if self.cfg.verbose && !rf.is_empty() {
            eprintln!("[{:.1}s] RF frontends: {} nets, deficit {}", self.elapsed(), rf.len(), rf.iter().map(|&n| self.deficit(n)).sum::<usize>());
        }
        if escapes {
            { let _t = Instant::now(); self.plan_fanout(); self.prof_add("planned fan-out", _t); }
        }
        if fanout {
            let fan = self.fanout_order();
            let _t = Instant::now();
            for &net in &fan {
                self.route_net(net);
            }
            self.prof_add("fan-out", _t);
            if self.cfg.verbose {
                eprintln!("[{:.1}s] fan-out of fine-pitch parts: {} nets, deficit {}", self.elapsed(), fan.len(), fan.iter().map(|&n| self.deficit(n)).sum::<usize>());
            }
        }
        let _t = Instant::now();
        for net in 0..nn {
            if self.is_plane(net) {
                self.route_net(net);
            }
        }
        self.prof_add("plane drops", _t);
        if self.cfg.verbose {
            eprintln!("[{:.1}s] planes done", self.elapsed());
        }
        self.trace_pairs("planes");
        let t_first = Instant::now();
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
            if self.elapsed() > self.cfg.seconds || self.hurried() {
                break;
            }
        }
        if self.stopped() {
            return;
        }
        self.prof_add("first pass", t_first);
        if self.cfg.verbose {
            eprintln!("[{:.1}s] first pass: deficit {}", self.elapsed(), self.total_deficit());
        }
        self.trace_pairs("first pass");
        // rip-up and reroute
        let t_rip = Instant::now();
        let mut best_def = usize::MAX;
        let mut stale = 0;
        // rip-up stops when it no longer pays: no new best for half the time it has run (and 4 rounds), the
        // end-game repairs do better from there than more rounds of the same
        let rip_start = self.elapsed();
        let mut best_at = rip_start;
        // the best board seen: a round can make things worse, the flow then goes on from the best
        type Snap = (World, Vec<(u32, u32)>, Vec<f64>, Vec<u32>, HashMap<u32, PairDiag>);
        let mut best_state: Option<(usize, Snap)> = None;
        // a repair that leaves the board more incomplete is undone (the history it added stays, so the next try goes
        // another way): rounds then only keep what pays, where before a round could undo two rounds' progress and
        // the deficit swung by a fifth from round to round on a crowded board (TOPOR_RIP_UNDO_OFF: kept as before)
        let rip_undo = std::env::var("TOPOR_RIP_UNDO_OFF").is_err();
        let mut undone = 0usize;
        for round in 0..self.cfg.rip_rounds {
            if self.stopped() {
                return;
            }
            if self.hurried() {
                break;
            }
            self.round = round as i32;
            self.update_pair_coupling();
            // rounds that no longer make progress go to the pair rescue instead
            let d = self.total_deficit();
            if best_state.as_ref().map_or(true, |b| d < b.0) {
                let _t = Instant::now();
                best_state = Some((d, (self.w.clone(), self.ties.clone(), self.history.clone(), self.junctions.clone(), self.pair_diag.clone())));
                self.prof_add("rip-up: snapshots", _t);
            }
            if d < best_def {
                best_def = d;
                stale = 0;
                best_at = self.elapsed();
            } else {
                stale += 1;
                if stale >= 4 && self.elapsed() - best_at > 0.5 * (self.elapsed() - rip_start) && std::env::var("TOPOR_NOSTALL").is_err() {
                    if self.cfg.verbose {
                        eprintln!("[{:.1}s] rip-up stops: no gain since {:.1}s (deficit {})", self.elapsed(), best_at, best_def);
                    }
                    break;
                }
                // with pairs released, rounds are cheap and the last few nets need patience
                if stale >= if self.cfg.fast { 2 } else if self.relax_pairs { 10 } else { 4 } {
                    // stuck with pairs held: let the remaining nets move them, once (the quick look stops here:
                    // its blocker report searches with pairs crossable, so a pair in the way shows there)
                    if self.relax_pairs || self.cfg.fast {
                        break;
                    }
                    self.relax_pairs = true;
                    // (TOPOR_RELAX_BUDGET: another factor; 1 was faster, no success needed more on dut-pi or EXT,
                    // but the changed searches cost EXT its last net)
                    self.budget_scale = std::env::var("TOPOR_RELAX_BUDGET").ok().and_then(|x| x.parse().ok()).unwrap_or(if self.cfg.fast { 1.0 } else { 3.0 });
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
                let undo = if rip_undo {
                    let _t = Instant::now();
                    let s = (self.total_deficit(), self.w.clone(), self.ties.clone(), self.junctions.clone(), self.pair_diag.clone());
                    self.prof_add("rip-up: undo snapshots", _t);
                    Some(s)
                } else {
                    None
                };
                let _t = Instant::now();
                let victims = self.repair_net(net);
                self.prof_add("rip-up: find what to rip (soft search)", _t);
                // a ripped pair is rebuilt first: it needs a corridor, the net only a track
                let paired = |me: &Router, v: u32| me.w.b.partner.get(v as usize).map_or(false, |&x| x != NONE);
                if std::env::var("TOPOR_PAIRFIRST").is_ok() {
                    for &v in &victims {
                        if paired(self, v) {
                            self.reroute(v);
                        }
                    }
                }
                let _t = Instant::now();
                self.route_net(net);
                self.prof_add("rip-up: route the net", _t);
                let _t = Instant::now();
                for &v in &victims {
                    if !paired(self, v) || std::env::var("TOPOR_PAIRFIRST").is_err() {
                        self.reroute(v);
                    }
                }
                self.prof_add("rip-up: reroute the ripped", _t);
                if let Some((d0, w, ties, junctions, pair_diag)) = undo {
                    if self.total_deficit() > d0 {
                        self.w = w;
                        self.ties = ties;
                        self.junctions = junctions;
                        self.pair_diag = pair_diag;
                        undone += 1;
                    }
                }
            }
            if self.cfg.verbose {
                eprintln!("[{:.1}s] rip round {}: deficit {}{}", self.elapsed(), round, self.total_deficit(), if rip_undo { format!(" ({} repairs undone so far)", undone) } else { String::new() });
            }
            self.trace_pairs(&format!("rip round {}", round));
        }
        self.prof_add("rip-up rounds", t_rip);
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
        if self.hurried() || self.cfg.fast {
            return;
        }
        { let _t = Instant::now(); self.clear_area_repairs(); self.prof_add("end game: clear areas", _t); }
        { let _t = Instant::now(); self.re_fan_repairs(); self.prof_add("end game: re-fan rows", _t); }
        { let _t = Instant::now(); self.rescue_pairs(); self.prof_add("end game: rescue pairs", _t); }
    }

    /// Nets still open after rip-up: a pad sealed in by other nets' vias and wires the soft search cannot
    /// pass. Clear other copper round the open pads, route the net first, then the displaced nets;
    /// kept only if the board ends up less incomplete. The trials are independent (nearly all are undone), so
    /// they run at once, each on its own copy of the board: the best one that pays is kept and the others try
    /// again from there, until a wave brings nothing.
    fn clear_area_repairs(&mut self) {
        let nn = self.w.b.nets.len() as u32;
        for radius in [1.2, 2.5] {
            loop {
                if self.elapsed() > self.cfg.seconds * 1.5 || self.hurried() {
                    return;
                }
                let failing: Vec<u32> = (0..nn).filter(|&x| self.deficit(x) > 0 && self.net_pads[x as usize].len() >= 1).collect();
                let before = self.total_deficit();
                if !self.adopt_best_trial(&failing, |r, &net| r.clear_area_trial(net, radius, before)) {
                    break;
                }
            }
        }
    }

    /// Run `trial` for every job, each on its own copy of the board, in parallel (a chunk of copies per round of
    /// threads, so memory stays bounded), and keep the copy that leaves the board least incomplete if that is
    /// better than now (ties: the earlier job). The trials' event lines are printed in job order. Returns
    /// whether a trial was kept.
    fn adopt_best_trial<T: Sync, F: Fn(&mut Router, &T) -> Option<String> + Sync>(&mut self, jobs: &[T], trial: F) -> bool {
        use rayon::prelude::*;
        let before = self.total_deficit();
        let chunk = rayon::current_num_threads().max(1);
        let mut best: Option<(usize, String, Router)> = None;
        for part in jobs.chunks(chunk) {
            if self.elapsed() > self.cfg.seconds * 1.5 || self.hurried() {
                break;
            }
            let _t = Instant::now();
            let copies: Vec<Router> = part
                .iter()
                .map(|_| {
                    let mut r = self.clone();
                    r.cfg.verbose = false;
                    r
                })
                .collect();
            self.prof_add("end game: board copies", _t);
            let done: Vec<(Option<String>, usize, Router)> = copies
                .into_par_iter()
                .zip(part.par_iter())
                .map(|(mut r, job)| {
                    let ev = trial(&mut r, job);
                    let d = r.total_deficit();
                    (ev, d, r)
                })
                .collect();
            for (ev, d, r) in done {
                let Some(ev) = ev else { continue };
                let better = d < before && best.as_ref().map_or(true, |b| d < b.0);
                if self.cfg.verbose {
                    eprintln!("[{:.1}s] {}{}", self.elapsed(), ev, if better { "" } else if d < before { ", not the best of its wave" } else { ", undone" });
                }
                if better {
                    best = Some((d, ev, r));
                }
            }
        }
        let Some((_, ev, r)) = best else { return false };
        let (verbose, stop, t0) = (self.cfg.verbose, self.stop.clone(), self.t0);
        *self = r;
        self.cfg.verbose = verbose;
        self.stop = stop;
        self.t0 = t0;
        self.events.push(ev);
        true
    }

    /// One clear-area trial for `net` on this copy of the board: the event line if there was copper to move.
    fn clear_area_trial(&mut self, net: u32, radius: f64, before: usize) -> Option<String> {
        if self.deficit(net) == 0 {
            return None;
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
            return None;
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
            return None;
        }
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
        Some(format!("{} was sealed in near {}: cleared {:.1} mm round its open pads ({} nets moved), deficit {} -> {}", self.w.b.nets[net as usize].name, open.iter().take(2).map(|c| format!("({:.2},{:.2})", c.x, c.y)).collect::<Vec<_>>().join(" "), radius, vnets.len(), before, self.total_deficit()))
    }

    /// On the finished board, rebuild poorly coupled pairs as ribbons; a rebuild is kept only if the
    /// board stays as complete and the pair comes out better coupled. The rebuilds are tried at once, each on its
    /// own copy of the board; taken in order, the first one kept is adopted and those after it are tried again on
    /// that board: the outcome of trying them one after another.
    fn improve_pairs(&mut self) {
        use rayon::prelude::*;
        let total = self.total_deficit();
        let mut todo: Vec<(String, f64, f64, f64)> = self.pair_report().into_iter().filter(|x| x.4 < 0.6).map(|(a, _, lp, ln, cf)| (a, lp, ln, cf)).collect();
        while !todo.is_empty() && self.elapsed() <= self.cfg.seconds * 1.5 {
            let copies: Vec<Router> = todo
                .iter()
                .map(|_| {
                    let mut r = self.clone();
                    r.cfg.verbose = false;
                    r
                })
                .collect();
            let tried: Vec<(Option<String>, bool, Router)> = copies
                .into_par_iter()
                .zip(todo.par_iter())
                .map(|(mut r, c)| {
                    let (msg, keep) = r.improve_pair(c, total);
                    (msg, keep, r)
                })
                .collect();
            let mut rest = Vec::new();
            for (k, (msg, keep, r)) in tried.into_iter().enumerate() {
                if let (Some(m), true) = (&msg, self.cfg.verbose) {
                    eprintln!("[{:.1}s] {}", self.elapsed(), m);
                }
                if keep {
                    let (verbose, t0) = (self.cfg.verbose, self.t0);
                    *self = r;
                    self.cfg.verbose = verbose;
                    self.t0 = t0;
                    rest = todo[k + 1..].to_vec();
                    break;
                }
            }
            todo = rest;
        }
    }

    /// Rebuild pair (named by its member a; its report before: lengths, coupled fraction) as a ribbon on this
    /// board: the line for the log (none when out of time) and whether the board is better for it.
    fn improve_pair(&mut self, c: &(String, f64, f64, f64), total: usize) -> (Option<String>, bool) {
        let (a, lp0, ln0, cf) = (c.0.as_str(), c.1, c.2, c.3);
        if self.elapsed() > self.cfg.seconds * 1.5 {
            return (None, false);
        }
        let len0 = if lp0.is_finite() && ln0.is_finite() { lp0.max(ln0) } else { f64::INFINITY };
        let p = self.w.b.net_index[a];
        let n = self.w.b.partner[p as usize];
        let (p, n) = (p.min(n), p.max(n));
        self.clear_pair(p, n);
        self.route_pair(p, n);
        self.route_net(p);
        self.route_net(n);
        let rep = self.pair_measure_routed(p, n);
        let after = rep.map(|x| x.2).unwrap_or(0.0);
        let len1 = rep.map(|x| x.0.max(x.1)).unwrap_or(f64::INFINITY);
        // better coupled, not bought with a detour
        let keep = self.total_deficit() <= total && after > cf + 0.05 && len1 <= len0 * 1.15 + 1.0;
        if keep {
            let note = format!("rebuilt as ribbon on the finished board: {:.0}% -> {:.0}% coupled", cf * 100.0, after * 100.0);
            self.pair_diag.entry(p).or_default().notes.push(note);
        }
        (Some(format!("improve pair {}: {:.0}% -> {:.0}% {}", a, cf * 100.0, after * 100.0, if keep { "kept" } else { "reverted" })), keep)
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

    /// Remove vias that connect nothing (one piece or none) and their lone pieces (unless the net comes apart: a
    /// via on fixed copper of its net, such as suggested copper kept as drawn, connects through it).
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
                    let net = ob.s.net;
                    let (before, drops) = (self.deficit(net), self.w.thread_drops);
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
                            // (a wire of another net that no longer threads the new triangles is dropped: not worth it)
                            if v.is_empty() && self.deficit(net) <= before && self.w.thread_drops == drops {
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
        // signal vias whose wires all sit on one layer: rip them and reconnect the net without them (not a pair
        // member's: it would come back as a lone wire, its coupling lost)
        for o in 0..self.w.objs.len() as u32 {
            let ob = &self.w.objs[o as usize];
            if !ob.alive || ob.s.kind != Kind::Via || ob.plane || self.w.b.partner.get(ob.s.net as usize).map_or(false, |&q| q != NONE) {
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
            // kept only if the net reconnects without it and nothing else comes apart
            let before = self.total_deficit();
            let snap = (self.w.clone(), self.ties.clone(), self.history.clone(), self.junctions.clone());
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
            if self.total_deficit() > before {
                self.w = snap.0;
                self.ties = snap.1;
                self.history = snap.2;
                self.junctions = snap.3;
            }
        }
    }

    /// Via bridges: a net that drops through a via to another layer for one wire and comes straight back
    /// through a second via (P -L- A =M= B -L- Q). Once the board has relaxed there is often room to stay on
    /// L: the three wires and both vias go, and P and Q are joined on L alone. Kept when everything checks and
    /// the new way is no longer than the bridge plus 1 mm for each via saved. Returns (bridges, removed).
    pub fn remove_bridges(&mut self) -> (usize, usize) {
        let (mut tried, mut removed) = (0, 0);
        let mut o = 0;
        while o < self.w.objs.len() as u32 {
            if let Some((pa, pb, pq)) = self.bridge_at(o) {
                tried += 1;
                if self.unbridge(pa, pb, pq) {
                    removed += 1;
                }
            }
            o += 1;
        }
        (tried, removed)
    }

    /// The bridge whose first via is o: (outer wire at o, the wire across, outer wire at the second via).
    fn bridge_at(&self, a: u32) -> Option<(u32, u32, u32)> {
        let w = &self.w;
        let via = |o: u32| {
            let ob = &w.objs[o as usize];
            ob.alive && ob.s.kind == Kind::Via && !ob.plane && !self.ties.iter().any(|t| t.0 == o || t.1 == o) && !w.b.partner.get(ob.s.net as usize).map_or(false, |&q| q != NONE)
        };
        let live = |o: u32| -> Vec<u32> { w.obj_pieces[o as usize].iter().cloned().filter(|&p| w.pieces[p as usize].alive && !w.pieces[p as usize].frozen).collect() };
        if !via(a) {
            return None;
        }
        let pa = live(a);
        if pa.len() != 2 || w.pieces[pa[0] as usize].layer == w.pieces[pa[1] as usize].layer {
            return None;
        }
        for (outer, inner) in [(pa[0], pa[1]), (pa[1], pa[0])] {
            let ip = &w.pieces[inner as usize];
            let b = if ip.a.obj == a { ip.b.obj } else { ip.a.obj };
            // each bridge once, from its lower via
            if b <= a || !via(b) {
                continue;
            }
            let pb = live(b);
            if pb.len() != 2 {
                continue;
            }
            let other = if pb[0] == inner { pb[1] } else { pb[0] };
            if w.pieces[other as usize].layer == w.pieces[outer as usize].layer {
                return Some((outer, inner, other));
            }
        }
        None
    }

    fn unbridge(&mut self, pa: u32, pb: u32, pq: u32) -> bool {
        let len = |w: &World, p: u32| w.pieces[p as usize].poly.windows(2).map(|s| s[0].dist(s[1])).sum::<f64>();
        let (net, l, width) = { let pc = &self.w.pieces[pa as usize]; (pc.net, pc.layer as usize, pc.width) };
        let vias: Vec<u32> = [pb].iter().flat_map(|&p| [self.w.pieces[p as usize].a.obj, self.w.pieces[p as usize].b.obj]).collect();
        let end = |w: &World, p: u32| { let pc = &w.pieces[p as usize]; if vias.contains(&pc.a.obj) { pc.b.obj } else { pc.a.obj } };
        let (from, to) = (end(&self.w, pa), end(&self.w, pq));
        let bridge = len(&self.w, pa) + len(&self.w, pb) + len(&self.w, pq);
        let (before, def0, drops) = (wire_len(&self.w), self.total_deficit(), self.w.thread_drops);
        let snap = (self.w.clone(), self.ties.clone(), self.history.clone(), self.junctions.clone());
        for p in [pa, pb, pq] {
            self.rip_piece(p);
        }
        let mut q = self.query(net, vec![from], vec![to], false);
        q.allowed = 1 << l;
        q.width = vec![width; self.w.nsig()];
        let ok = match self.search(&q) {
            Some(path) if path.vias.is_empty() => self.commit(net, &path, &q.width).is_ok(),
            _ => false,
        };
        let added = wire_len(&self.w) - (before - bridge);
        if std::env::var("TOPOR_BRIDGEDBG").is_ok() {
            let at: Vec<String> = vias.iter().map(|&v| { let c = snap.0.objs[v as usize].s.center; format!("({:.2},{:.2})", c.x, c.y) }).collect();
            eprintln!("  bridge {} L{} {} -> {} vias {:?} hop {:.2} mm: {} (bridge {:.2} mm, new {:.2} mm)", self.w.b.nets[net as usize].name, l, self.w.objs[from as usize].s.name, self.w.objs[to as usize].s.name, at, len(&snap.0, pb), if !ok { "no way on the layer" } else if added > bridge + 2.0 { "longer" } else { "ok" }, bridge, added);
        }
        if ok && self.w.thread_drops == drops && self.total_deficit() <= def0 && added <= bridge + 2.0 {
            return true;
        }
        self.w = snap.0;
        self.ties = snap.1;
        self.history = snap.2;
        self.junctions = snap.3;
        false
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

    /// Hand-drawn copper given as suggestions (model `suggested_copper`: `tracks` [{net, layer, width, pts}],
    /// `vias` [{net, x, y}]), imported before routing as the router's own copper of its nets: it bootstraps a board
    /// that is hard to route, while rip-up may still move it (at a small premium) and the spring relaxation improves
    /// it like any wire. Vias become router vias; track segments are joined into polylines between pads, vias and
    /// branch points and threaded into the triangulation along their course. What breaks a rule is left out (each
    /// noted) for the router to do itself. Returns (vias kept, vias left out, polylines kept, left out, notes, and
    /// the indices of the vias and tracks left out for a rule, which the caller may keep as fixed copper instead).
    pub fn import_suggestions(&mut self, v: &serde_json::Value) -> Suggested {
        let f = |x: &serde_json::Value| x.as_f64().unwrap_or(0.0);
        let mut notes = Vec::new();
        let (mut vok, mut vbad, mut tok, mut tbad) = (0, 0, 0, 0);
        let (mut bad_vias, mut bad_tracks): (Vec<usize>, Vec<usize>) = (Vec::new(), Vec::new());
        let mut vias: Vec<(u32, P2, u32)> = Vec::new();
        for (vi, x) in v["vias"].as_array().into_iter().flatten().enumerate() {
            let Some(&net) = x["net"].as_str().and_then(|n| self.w.b.net_index.get(n)) else {
                vbad += 1;
                continue;
            };
            let pos = p2(f(&x["x"]), f(&x["y"]));
            self.w.begin();
            let res = self.w.add_via(net, pos).map_err(|e| e.to_string()).and_then(|o| {
                let mut viol = Vec::new();
                self.w.check_via_at(net, pos, o, &mut viol, true);
                match viol.first() {
                    Some(x) => {
                        let with = x.objs.first().map(|&q| self.w.objs[q as usize].s.name.clone())
                            .or(x.pieces.first().map(|&q| format!("a wire of {}", self.w.b.nets[self.w.pieces[q as usize].net as usize].name)));
                        Err(match with {
                            Some(wn) if x.dist >= 0.0 => format!("{} ({:.3} mm from {}, needs {:.3})", x.what, x.dist, wn, x.need),
                            Some(wn) => format!("{} ({})", x.what, wn),
                            None => x.what.clone(),
                        })
                    }
                    None => Ok(o),
                }
            });
            match res {
                Ok(o) => {
                    self.w.commit();
                    vias.push((net, pos, o));
                    vok += 1;
                }
                Err(e) => {
                    self.w.rollback();
                    vbad += 1;
                    bad_vias.push(vi);
                    notes.push(format!("via of {} at ({:.2},{:.2}) left out: {}", self.w.b.nets[net as usize].name, pos.x, pos.y, e));
                }
            }
        }
        // a via drawn on copper kept as drawn (fixed: zone pieces of its net) joins it, as the fill joins a pad
        // (else the net comes apart there: the wire drawn up to the via is a dead end, and the router goes round)
        for &(net, _, o) in &vias {
            for z in 0..self.w.objs.len() as u32 {
                let (zs, vs) = (&self.w.objs[z as usize].s, &self.w.objs[o as usize].s);
                if zs.net == net && self.w.objs[z as usize].alive && crate::world::fixed_copper(&zs.name) && zs.layers & vs.layers != 0 && crate::world::core_dist(zs, vs) - zs.r - vs.r < 0.03 && !self.ties.contains(&(z, o)) {
                    self.ties.push((z, o));
                }
            }
        }
        // segments by (net, layer, width), their end points joined at 1 um
        let key = |p: P2| ((p.x * 1000.0).round() as i64, (p.y * 1000.0).round() as i64);
        let mut segs: Vec<(u32, usize, f64, P2, P2)> = Vec::new();
        let mut seg_track: Vec<usize> = Vec::new();
        for (ti, t) in v["tracks"].as_array().into_iter().flatten().enumerate() {
            let net = t["net"].as_str().and_then(|n| self.w.b.net_index.get(n)).cloned();
            let l = (0..self.w.nsig()).find(|&l| Some(self.w.b.sig_name(l)) == t["layer"].as_str());
            let (Some(net), Some(l)) = (net, l) else {
                tbad += 1;
                continue;
            };
            let w = f(&t["width"]).max(self.w.b.nets[net as usize].wmin[l]);
            let pts: Vec<P2> = t["pts"].as_array().into_iter().flatten().map(|q| p2(f(&q[0]), f(&q[1]))).collect();
            for k in 0..pts.len().saturating_sub(1) {
                if pts[k].dist(pts[k + 1]) > 1e-6 {
                    segs.push((net, l, w, pts[k], pts[k + 1]));
                    seg_track.push(ti);
                }
            }
        }
        // a track that ends on another track of its net and layer, inside its copper but off its points (a T drawn
        // by overlap, as an ESD tap off a line running past it), joins it there: the other track is split at the
        // foot of that end, and the end moved onto it (else the tap stays apart, and the router draws the line again
        // through the part)
        {
            let mut ends: HashMap<(u32, usize, (i64, i64)), usize> = HashMap::default();
            for &(net, l, _, a, b) in &segs {
                *ends.entry((net, l, key(a))).or_insert(0) += 1;
                *ends.entry((net, l, key(b))).or_insert(0) += 1;
            }
            let at_via = |net: u32, p: P2| vias.iter().any(|&(n, q, _)| n == net && q.dist(p) < 0.02);
            let at_pad = |me: &Router, net: u32, l: usize, p: P2| {
                me.net_pads[net as usize].iter().any(|&o| {
                    let s = &me.w.objs[o as usize].s;
                    s.layers & (1 << l) != 0 && point_hull_dist(&s.pts, p) - s.r <= 1e-3
                })
            };
            for i in 0..segs.len() {
                for side in 0..2 {
                    let (net, l, w, a, b) = segs[i];
                    let e = if side == 0 { a } else { b };
                    if ends.get(&(net, l, key(e))).cloned().unwrap_or(0) != 1 || at_via(net, e) || at_pad(self, net, l, e) {
                        continue;
                    }
                    let foot = (0..segs.len()).filter(|&j| j != i && segs[j].0 == net && segs[j].1 == l).find_map(|j| {
                        let (_, _, wj, c, d) = segs[j];
                        let t = seg_param(c, d, e);
                        let q = c.lerp(d, t);
                        (q.dist(e) < (w + wj) / 2.0 && q.dist(c) > 1e-3 && q.dist(d) > 1e-3).then_some((j, q))
                    });
                    let Some((j, q)) = foot else { continue };
                    let (_, _, wj, c, d) = segs[j];
                    segs[j] = (net, l, wj, c, q);
                    segs.push((net, l, wj, q, d));
                    seg_track.push(seg_track[j]);
                    if side == 0 {
                        segs[i].3 = q;
                    } else {
                        segs[i].4 = q;
                    }
                    ends.remove(&(net, l, key(e)));
                    *ends.entry((net, l, key(q))).or_insert(0) += 3;
                }
            }
        }
        // stops: points where a chain must end (a pad or via of the net there, a branch, an end, a width change)
        let mut deg: HashMap<(u32, usize, (i64, i64)), usize> = HashMap::default();
        for &(net, l, _, a, b) in &segs {
            *deg.entry((net, l, key(a))).or_insert(0) += 1;
            *deg.entry((net, l, key(b))).or_insert(0) += 1;
        }
        let on_terminal = |me: &Router, net: u32, l: usize, p: P2| -> bool {
            vias.iter().any(|&(n, q, _)| n == net && q.dist(p) < 0.02)
                || me.net_pads[net as usize].iter().any(|&o| {
                    let s = &me.w.objs[o as usize].s;
                    s.layers & (1 << l) != 0 && point_hull_dist(&s.pts, p) - s.r <= 1e-3
                })
        };
        let mut used = vec![false; segs.len()];
        let mut chains: Vec<(u32, usize, f64, Vec<P2>, Vec<usize>)> = Vec::new();
        for i in 0..segs.len() {
            if used[i] {
                continue;
            }
            used[i] = true;
            let (net, l, w, a, b) = segs[i];
            let mut poly = vec![a, b];
            let mut from = vec![seg_track[i]];
            // grow both ways through points of degree two that are no terminal, at the same width
            for dir in 0..2 {
                loop {
                    let end = if dir == 0 { *poly.last().unwrap() } else { poly[0] };
                    if deg.get(&(net, l, key(end))).cloned().unwrap_or(0) != 2 || on_terminal(self, net, l, end) {
                        break;
                    }
                    let next = (0..segs.len()).find(|&j| !used[j] && segs[j].0 == net && segs[j].1 == l && (segs[j].2 - w).abs() < 1e-6 && (key(segs[j].3) == key(end) || key(segs[j].4) == key(end)));
                    let Some(j) = next else { break };
                    used[j] = true;
                    from.push(seg_track[j]);
                    let other = if key(segs[j].3) == key(end) { segs[j].4 } else { segs[j].3 };
                    if dir == 0 {
                        poly.push(other);
                    } else {
                        poly.insert(0, other);
                    }
                }
            }
            chains.push((net, l, w, poly, from));
        }
        // thread each chain between its end terminals (pad, via, or a junction shared by chains meeting there)
        let mut joints: Vec<(u32, usize, P2, u32)> = Vec::new();
        for (net, l, w, mut poly, from) in chains {
            let name = self.w.b.nets[net as usize].name.clone();
            self.w.begin();
            let mut fresh: Vec<u32> = Vec::new();
            let mut terms = [Term { obj: NONE, vtx: NONE }; 2];
            let mut bad: Option<String> = None;
            for side in 0..2 {
                let at = if side == 0 { poly[0] } else { *poly.last().unwrap() };
                let via = vias.iter().find(|&&(n, q, _)| n == net && q.dist(at) < 0.02).map(|x| x.2);
                let pad = self.net_pads[net as usize].iter().cloned().find(|&o| {
                    let s = &self.w.objs[o as usize].s;
                    s.layers & (1 << l) != 0 && point_hull_dist(&s.pts, at) - s.r <= 1e-3
                });
                let joint = joints.iter().find(|x| x.0 == net && x.1 == l && x.2.dist(at) < 1e-3).map(|x| x.3);
                let (tm, q) = if let Some(o) = via {
                    (Term { obj: o, vtx: self.w.objs[o as usize].lv[l] }, None)
                } else if let Some(o) = pad {
                    match self.pad_term(o, l, at) {
                        Some((tm, q)) => (tm, Some(q)),
                        None => {
                            bad = Some("its pad end".into());
                            break;
                        }
                    }
                } else if let Some(j) = joint {
                    (Term { obj: j, vtx: self.w.objs[j as usize].lv[l] }, None)
                } else {
                    match self.w.add_junction(net, at, w / 2.0, l) {
                        Ok(j) => {
                            fresh.push(j);
                            (Term { obj: j, vtx: self.w.objs[j as usize].lv[l] }, None)
                        }
                        Err(e) => {
                            bad = Some(format!("an end: {}", e));
                            break;
                        }
                    }
                };
                terms[side] = tm;
                // the wire ends at its terminal's vertex, or at a polygon pad's centre
                let p = if tm.vtx != NONE { self.w.vpos(l, tm.vtx) } else { q.unwrap_or(at) };
                if side == 0 {
                    if tm.vtx == NONE { poly.insert(0, p) } else { poly[0] = p }
                } else if tm.vtx == NONE {
                    poly.push(p)
                } else {
                    *poly.last_mut().unwrap() = p
                }
            }
            poly.dedup_by(|a, b| a.dist(*b) < 1e-6);
            let res = match bad {
                Some(e) => Err(e),
                None if poly.len() < 2 || terms[0].obj == terms[1].obj => Err("too short".into()),
                None => match self.w.thread_poly(l, &poly, terms[0], terms[1]) {
                    None => Err("does not thread through the triangulation".into()),
                    Some(cx) => {
                        let pid = self.w.add_piece(net, l, w, terms[0], terms[1], &cx);
                        self.w.seed_poly(pid, poly.clone());
                        self.w.realize_dirty();
                        let mut viol = Vec::new();
                        self.w.check_piece(pid, &mut viol);
                        match viol.first() {
                            Some(x) => Err(format!("{} at ({:.2},{:.2})", x.what, x.at.x, x.at.y)),
                            None => Ok(pid),
                        }
                    }
                },
            };
            match res {
                Ok(pid) => {
                    self.w.commit();
                    for &j in &fresh {
                        joints.push((net, l, self.w.objs[j as usize].s.center, j));
                        self.junctions.push(j);
                    }
                    self.suggested.insert(pid);
                    self.suggested_nets.insert(net);
                    tok += 1;
                }
                Err(e) => {
                    self.w.rollback();
                    tbad += 1;
                    bad_tracks.extend(from.iter().cloned());
                    notes.push(format!("track of {} on {} from ({:.2},{:.2}) left out: {}", name, self.w.b.sig_name(l), poly[0].x, poly[0].y, e));
                }
            }
        }
        // copper that reaches none of its net's pads or vias is of no use as a start: dropped
        for net in self.suggested_nets.clone() {
            for g in self.components(net) {
                if g.iter().all(|&o| self.w.objs[o as usize].s.name.starts_with("junction:")) {
                    for &o in &g {
                        for p in self.w.obj_pieces[o as usize].clone() {
                            if self.w.pieces[p as usize].alive {
                                self.rip_piece(p);
                                tok -= 1;
                                tbad += 1;
                            }
                        }
                        let _ = self.w.remove_via(o);
                    }
                    notes.push(format!("suggested copper of {} reaching none of its pads left out", self.w.b.nets[net as usize].name));
                }
            }
        }
        self.junctions.retain(|&j| self.w.objs[j as usize].alive);
        // pairs drawn whole: both members complete as suggested
        for &net in &self.suggested_nets {
            let pn = self.w.b.partner.get(net as usize).cloned().unwrap_or(NONE);
            if pn != NONE && self.deficit(net) == 0 && self.deficit(pn) == 0 {
                self.drawn_pairs.insert(net);
                self.drawn_pairs.insert(pn);
            }
        }
        bad_tracks.sort();
        bad_tracks.dedup();
        Suggested { vias_kept: vok, vias_left_out: vbad, tracks_kept: tok, tracks_left_out: tbad, notes, bad_vias, bad_tracks }
    }

    /// The blocker report: for each open group of each incomplete net (two per net, forty in all), the cheapest
    /// way to the rest of the net through other nets' copper (pairs included, the normal budget), and what is next
    /// to its open pads. Nothing on the board changes.
    pub fn blockers(&mut self) -> Vec<Blocker> {
        let mut out = Vec::new();
        let saved = (self.relax_pairs, self.protect_pairs, self.budget_scale);
        self.relax_pairs = true;
        self.protect_pairs = false;
        self.budget_scale = 1.0;
        for net in self.incomplete_nets() {
            let comps = self.components(net);
            let plane = self.is_plane(net);
            let groups: Vec<Vec<u32>> = if plane { comps.iter().filter(|c| !self.plane_ok(net, c)).cloned().collect() } else { comps.iter().skip(1).cloned().collect() };
            for g in groups.into_iter().take(2) {
                if out.len() >= 40 {
                    break;
                }
                let pad = |o: u32| {
                    let s = &self.w.objs[o as usize].s;
                    s.kind == Kind::Pad && !crate::world::fixed_copper(&s.name) && !s.name.starts_with("junction:")
                };
                let open_pads: Vec<(String, P2)> = g.iter().filter(|&&o| pad(o)).map(|&o| (self.w.objs[o as usize].s.name.clone(), self.w.objs[o as usize].s.center)).collect();
                if open_pads.is_empty() {
                    continue;
                }
                let targets: Vec<u32> = if plane { vec![] } else { comps.iter().filter(|c| **c != g).flatten().cloned().collect() };
                let mut q = self.query(net, g.clone(), targets, true);
                q.soft = true;
                q.max_exp = self.cfg.max_exp * 3;
                if plane {
                    q.plane_goal = true;
                    q.max_len = 5.0;
                    q.allowed &= g.iter().fold(0u32, |m, &o| m | self.w.objs[o as usize].s.layers);
                }
                let mut found = self.search(&q);
                // (a plane pad a neck-down rule covers drops at the neck width: see repair_net)
                let neck = if plane { g.iter().filter_map(|&o| self.w.b.neck_width.get(&self.w.objs[o as usize].s.name).cloned()).fold(f64::INFINITY, f64::min) } else { f64::INFINITY };
                if found.is_none() && neck.is_finite() {
                    for x in q.width.iter_mut() {
                        *x = x.min(neck);
                    }
                    q.max_len = 2.0;
                    found = self.search(&q);
                }
                let sealed = found.is_none() && self.last_exp < q.max_exp;
                let way = found.filter(|p| p.cost.is_finite()).map(|p| self.way_of(&p));
                // what is next to the open pads: other nets' wires and vias, their fixed copper, other parts
                let near = |q: P2| open_pads.iter().any(|(_, c)| c.dist(q) < 1.0);
                let mut nets: HashMap<u32, (usize, usize)> = HashMap::default();
                for pc in self.w.pieces.iter().filter(|p| p.alive && p.net != net) {
                    if pc.poly.windows(2).any(|sg| open_pads.iter().any(|(_, c)| seg_point_dist(sg[0], sg[1], *c) < 1.0)) {
                        nets.entry(pc.net).or_insert((0, 0)).0 += 1;
                    }
                }
                let mut parts: Vec<String> = Vec::new();
                let own: Vec<String> = open_pads.iter().map(|(n, _)| n.split('.').next().unwrap_or("").to_string()).collect();
                for o in self.w.objs.iter().filter(|o| o.alive && o.s.net != net) {
                    if o.s.kind == Kind::Via && near(o.s.center) {
                        nets.entry(o.s.net).or_insert((0, 0)).0 += 1;
                        continue;
                    }
                    if o.s.kind != Kind::Pad || o.s.name.starts_with("junction:") || !open_pads.iter().any(|(_, c)| point_hull_dist(&o.s.pts, *c) - o.s.r < 1.0) {
                        continue;
                    }
                    if o.s.net != NONE {
                        nets.entry(o.s.net).or_insert((0, 0)).1 += 1;
                    }
                    // a real part (hand-drawn copper comes as pads of a "via"/"track" part, or as zone pieces)
                    let part = o.s.name.split('.').next().unwrap_or("").to_string();
                    if !crate::world::fixed_copper(&o.s.name) && !matches!(part.as_str(), "via" | "track" | "copper") && !own.contains(&part) && !parts.contains(&part) {
                        parts.push(part);
                    }
                }
                let mut near_nets: Vec<(String, usize, usize)> = nets.into_iter().map(|(n, (k, f))| (self.w.b.nets[n as usize].name.clone(), k, f)).collect();
                near_nets.sort_by(|a, b| (b.1 + b.2).cmp(&(a.1 + a.2)).then(a.0.cmp(&b.0)));
                parts.sort();
                out.push(Blocker { net: self.w.b.nets[net as usize].name.clone(), open_pads, way, sealed, near_nets, near_parts: parts });
            }
        }
        (self.relax_pairs, self.protect_pairs, self.budget_scale) = saved;
        out
    }

    /// A found (soft) path as a way for the blocker report.
    fn way_of(&self, p: &PathOut) -> Way {
        let mut by: HashMap<u32, usize> = HashMap::default();
        for &v in &p.victims {
            *by.entry(self.w.pieces[v as usize].net).or_insert(0) += 1;
        }
        let mut through: Vec<(String, usize)> = by.into_iter().map(|(n, k)| (self.w.b.nets[n as usize].name.clone(), k)).collect();
        through.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        let tri_at = |l: usize, t: u32| {
            let tr = self.w.layers[l].cdt.tris[t as usize];
            (self.w.vpos(l, tr.v[0]) + self.w.vpos(l, tr.v[1]) + self.w.vpos(l, tr.v[2])) * (1.0 / 3.0)
        };
        let crossings = p.cross_tris.iter().map(|&(l, t)| (self.w.b.sig_name(l as usize).to_string(), tri_at(l as usize, t))).collect();
        let end_pos = |l: usize, e: &End| match e {
            End::Obj(t) => if t.vtx != NONE { self.w.vpos(l, t.vtx) } else { self.w.objs[t.obj as usize].s.center },
            End::NewVia(i) => p.vias[*i],
        };
        let legs = p
            .segs
            .iter()
            .map(|sg| {
                let l = sg.layer;
                let mut pts = vec![end_pos(l, &sg.a)];
                for &(e, _, _) in &sg.cx {
                    let ed = self.w.layers[l].cdt.edges[e as usize];
                    pts.push((self.w.vpos(l, ed.a) + self.w.vpos(l, ed.b)) * 0.5);
                }
                pts.push(end_pos(l, &sg.b));
                (self.w.b.sig_name(l).to_string(), pts)
            })
            .collect();
        Way { through, crossings, legs, vias: p.vias.clone() }
    }

    pub fn incomplete_nets(&self) -> Vec<u32> {
        (0..self.w.b.nets.len() as u32).filter(|&n| self.deficit(n) > 0).collect()
    }
}

/// Total length of the live wires.
pub fn wire_len(w: &World) -> f64 {
    w.pieces.iter().filter(|p| p.alive).map(|p| p.poly.windows(2).map(|s| s[0].dist(s[1])).sum::<f64>()).sum()
}

/// Where a dogbone net's vias may go: on rings just outside each of its pads (the via's copper clear of the pad
/// and out to about a pitch beyond it), wherever a via fits among the fixed copper. In a BGA these are the gaps
/// between the balls round each of the net's own balls; at a part's edge also the room outside it.
fn dogbone_spots(w: &World, net: u32, pads: &[u32]) -> Vec<P2> {
    let rv = w.b.via_size(net).0 / 2.0;
    let mut out: Vec<P2> = Vec::new();
    for &o in pads {
        let s = &w.objs[o as usize].s;
        if s.layers == (1 << w.nsig()) - 1 {
            continue; // a through pad is its own via
        }
        let ext = s.pts.iter().map(|q| q.dist(s.center)).fold(0.0, f64::max) + s.r;
        let r0 = ext + rv + 0.02;
        for k in 0..10 {
            let r = r0 + 0.04 * k as f64;
            for a in 0..24 {
                let an = a as f64 * std::f64::consts::TAU / 24.0;
                let x = s.center + p2(an.cos(), an.sin()) * r;
                if w.via_fits(net, x) && !out.iter().any(|y| y.dist(x) < 0.02) {
                    out.push(x);
                }
            }
        }
    }
    out
}
