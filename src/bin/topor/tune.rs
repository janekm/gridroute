//! Length matching: after routing, serpentines lengthen the shorter nets of each length group
//! (`cold_routing.length_groups`: `[{name, nets, max_skew, target?}]`) until every net is within its group's skew
//! of the longest (or of `target`, when that is longer). A differential pair is a group of its two nets.
//!
//! A net's length is its pads' distance along its own copper (two-pad nets; vias count as nothing). Each step
//! puts a row of chamfered bumps into one straight stretch of the net's wires and keeps it only if the piece
//! still passes the exact clearance check; the stretches are tried longest first, on both sides (a pair member
//! bumps away from its partner), with low bumps first (they need the least room beside the wire). A net is aimed
//! at the middle of the upper half of its group's window, so the group ends within half its skew where the room
//! allows. Final geometry only: nothing is routed or relaxed afterwards.
//!
//! Meander geometry (`cold_routing`): `meander_gap` (edge-to-edge gap between legs; default the larger of three
//! track widths and the clearance) and `meander_max_amplitude` (bump height, default 1 mm).
//!
//! Room to tune: taut wires packed between pads leave none, so a group may give `layer_cost` ({layer: cost of a
//! millimetre, at least 1}) for its nets' routing, e.g. {"F.Cu": 2.5} to take the long runs off a BGA pad layer to
//! inner layers, where serpentines fit; the pad escapes stay short.

use crate::cdt::NONE;
use crate::geom::*;
use crate::router::Router;

/// One net of a length group: its length before and after, the bumps added, and why it is not tuned (if not).
#[derive(Clone, Debug)]
pub struct NetTune {
    pub net: u32,
    pub before: Option<f64>,
    pub after: Option<f64>,
    pub bumps: usize,
    pub note: Option<String>,
}

#[derive(Clone, Debug)]
pub struct GroupTune {
    pub name: String,
    pub max_skew: f64,
    /// The length the group was matched to (its longest net, or the given target when longer).
    pub target: f64,
    pub nets: Vec<NetTune>,
    pub skew_before: f64,
    pub skew_after: f64,
    pub met: bool,
    pub unknown: Vec<String>,
}

/// Length a run of `nb` bumps of height `amp` adds over the straight line (chamfer `c` at every corner).
fn bump_gain(amp: f64, c: f64) -> f64 {
    2.0 * amp - (8.0 - 4.0 * std::f64::consts::SQRT_2) * c
}

/// Chamfer of a bump: a quarter of the leg pitch, or half its height when it is low.
fn chamfer(g: f64, amp: f64) -> f64 {
    (g / 4.0).min(amp / 2.0)
}

/// Bump height adding `d` per bump at leg pitch g.
fn amp_for(d: f64, g: f64) -> f64 {
    let c = g / 4.0;
    let amp = (d + (8.0 - 4.0 * std::f64::consts::SQRT_2) * c) / 2.0;
    if amp >= 2.0 * c {
        amp
    } else {
        // low bumps: chamfers only (c = amp / 2), gain 2 amp - 2.34 amp / 2
        d / (2.0 - (4.0 - 2.0 * std::f64::consts::SQRT_2))
    }
}

/// The stretch a..b with `nb` bumps of height amp on `side`, the first leg at `start` along it.
fn serpentine(a: P2, b: P2, start: f64, nb: usize, g: f64, amp: f64, side: f64) -> Vec<P2> {
    let u = (b - a).unit();
    let n = u.perp() * side;
    let c = chamfer(g, amp);
    let at = |x: f64, y: f64| a + u * x + n * y;
    let mut pts = vec![a];
    for i in 0..nb {
        let x = start + i as f64 * 2.0 * g;
        pts.extend([at(x - c, 0.0), at(x, c), at(x, amp - c), at(x + c, amp), at(x + g - c, amp), at(x + g, amp - c), at(x + g, c), at(x + g + c, 0.0)]);
    }
    pts.push(b);
    // (low bumps: drop the zero-length legs)
    pts.dedup_by(|p, q| p.dist(*q) < 1e-9);
    pts
}

impl Router {
    /// A two-pad net's length between its pads along its own wires; None for another net or an open one.
    pub fn net_length(&self, net: u32) -> Option<f64> {
        let pads = &self.net_pads[net as usize];
        if pads.len() != 2 {
            return None;
        }
        self.net_path_len(net, pads[0], pads[1])
    }

    /// Match the lengths of every length group (loosest first, so a tighter group inside one finishes last).
    pub fn match_lengths(&mut self) -> Vec<GroupTune> {
        let mut groups = self.w.b.length_groups.clone();
        groups.sort_by(|a, b| b.max_skew.partial_cmp(&a.max_skew).unwrap());
        let mut out = Vec::new();
        for g in &groups {
            let before: Vec<Option<f64>> = g.nets.iter().map(|&n| self.net_length(n)).collect();
            let known: Vec<f64> = before.iter().flatten().cloned().collect();
            let longest = known.iter().cloned().fold(0.0, f64::max);
            let target = g.target.map_or(longest, |t| t.max(longest));
            let lo = target - g.max_skew / 2.0;
            let goal = target - g.max_skew / 4.0;
            let base_notes: Vec<Option<String>> = before
                .iter()
                .zip(&g.nets)
                .map(|(l, &n)| match l {
                    Some(_) => None,
                    None if self.net_pads[n as usize].len() != 2 => Some(format!("{} pads: only two-pad nets are measured", self.net_pads[n as usize].len())),
                    None => Some("not routed".to_string()),
                })
                .collect();
            // the shortest first (it needs the most room); when some come up short, again from the start with
            // those first, keeping the attempt that leaves the least length missing
            let mut order: Vec<usize> = (0..g.nets.len()).filter(|&i| before[i].is_some()).collect();
            order.sort_by(|&i, &j| before[i].partial_cmp(&before[j]).unwrap());
            let start = self.w.clone();
            let mut best: Option<(f64, crate::world::World, Vec<usize>, Vec<Option<String>>)> = None;
            for _attempt in 0..4 {
                let mut bumps = vec![0usize; g.nets.len()];
                let mut notes = base_notes.clone();
                let mut short: Vec<usize> = Vec::new();
                let mut missing = 0.0;
                for &i in &order {
                    let net = g.nets[i];
                    let mut len = self.net_length(net).unwrap_or(before[i].unwrap());
                    // a pair member whose partner is in the group too: while both are short, both together
                    let partner = self.w.b.partner.get(net as usize).cloned().unwrap_or(NONE);
                    let j = g.nets.iter().position(|&x| x == partner && partner != NONE);
                    for _ in 0..24 {
                        if len >= lo - 1e-9 {
                            break;
                        }
                        if let Some(j) = j {
                            let lj = self.net_length(partner).unwrap_or(f64::INFINITY);
                            if lj < lo - 1e-9 {
                                if let Some(k) = self.pair_meander_step(net, partner, (goal - len).min(goal - lj)) {
                                    bumps[i] += k;
                                    bumps[j] += k;
                                    len = self.net_length(net).unwrap_or(len);
                                    continue;
                                }
                            }
                        }
                        match self.meander_step(net, goal - len) {
                            Some(k) => bumps[i] += k,
                            None => break,
                        }
                        len = self.net_length(net).unwrap_or(len);
                    }
                    if len < lo - 1e-9 {
                        notes[i] = Some(format!("{:.3} mm short of the window: no room for more bumps", lo - len));
                        short.push(i);
                        missing += lo - len;
                    }
                }
                let better = best.as_ref().map_or(true, |b| missing < b.0 - 1e-6);
                if better {
                    best = Some((missing, self.w.clone(), bumps, notes));
                }
                if short.is_empty() || short.iter().all(|i| order.iter().position(|x| x == i) < Some(short.len())) {
                    break;
                }
                let rest: Vec<usize> = order.iter().cloned().filter(|i| !short.contains(i)).collect();
                order = short.into_iter().chain(rest).collect();
                self.w = start.clone();
            }
            let (_, w, bumps, notes) = best.unwrap();
            self.w = w;
            let after: Vec<Option<f64>> = g.nets.iter().map(|&n| self.net_length(n)).collect();
            let skew = |v: &[Option<f64>]| {
                let k: Vec<f64> = v.iter().flatten().cloned().collect();
                if k.is_empty() {
                    0.0
                } else {
                    k.iter().cloned().fold(f64::MIN, f64::max) - k.iter().cloned().fold(f64::MAX, f64::min)
                }
            };
            let (sb, sa) = (skew(&before), skew(&after));
            let met = after.iter().all(|x| x.is_some()) && sa <= g.max_skew + 1e-9 && g.unknown.is_empty();
            out.push(GroupTune {
                name: g.name.clone(),
                max_skew: g.max_skew,
                target,
                nets: (0..g.nets.len()).map(|i| NetTune { net: g.nets[i], before: before[i], after: after[i], bumps: bumps[i], note: notes[i].clone() }).collect(),
                skew_before: sb,
                skew_after: sa,
                met,
                unknown: g.unknown.clone(),
            });
        }
        out
    }

    /// One serpentine on both members of a pair where they run side by side, adding up to `want` to each: the
    /// outer member bumps away from the other, which follows it at the pair's spacing (the offset of the outer
    /// member's bumps), so both gain the same and stay coupled. The bumps added, or None when nothing fits.
    fn pair_meander_step(&mut self, p: u32, n: u32, want: f64) -> Option<usize> {
        if want < 0.05 {
            return None;
        }
        let ends = |r: &Router, net: u32| -> Option<Vec<(u32, u32)>> {
            let pads = &r.net_pads[net as usize];
            if pads.len() != 2 {
                return None;
            }
            r.piece_chain(net, pads[0], pads[1])
        };
        let (cp, cn) = (ends(self, p)?, ends(self, n)?);
        let w0 = cp.iter().chain(&cn).map(|&(x, _)| self.w.pieces[x as usize].width).fold(0.0, f64::max);
        let clear = self.w.b.classes[self.w.b.nets[p as usize].class].clear;
        let gap_m = self.w.b.meander_gap.unwrap_or((3.0 * w0).max(clear));
        let max_amp = self.w.b.meander_max_amp;
        // side-by-side straight stretches: (overlap, outer piece, its segment, inner piece, its segment, overlap
        // start and end along the outer segment, spacing); either member may be the outer one
        let mut cands: Vec<(f64, u32, usize, u32, usize, f64, f64, f64)> = Vec::new();
        for (ca, cb) in [(&cp, &cn), (&cn, &cp)] {
            for &(po, _) in ca.iter() {
                for &(pi, _) in cb.iter() {
                    let (o, i) = (&self.w.pieces[po as usize], &self.w.pieces[pi as usize]);
                    if o.layer != i.layer {
                        continue;
                    }
                    for ko in 0..o.poly.len().saturating_sub(1) {
                        let (a, b) = (o.poly[ko], o.poly[ko + 1]);
                        let len = a.dist(b);
                        if len < 2.0 {
                            continue;
                        }
                        let u = (b - a).unit();
                        for ki in 0..i.poly.len().saturating_sub(1) {
                            let (c, e) = (i.poly[ki], i.poly[ki + 1]);
                            if c.dist(e) < 2.0 || u.cross((e - c).unit()).abs() > 0.01 {
                                continue;
                            }
                            let (dc, de) = (u.cross(c - a), u.cross(e - a));
                            if (dc - de).abs() > 0.01 || dc.abs() > 3.0 * (w0 + clear) + 0.1 {
                                continue;
                            }
                            let (tc, te) = ((c - a).dot(u), (e - a).dot(u));
                            let (t0, t1) = (tc.min(te).max(0.0), tc.max(te).min(len));
                            if t1 - t0 > 2.0 {
                                cands.push((t1 - t0, po, ko, pi, ki, t0, t1, dc));
                            }
                        }
                    }
                }
            }
        }
        cands.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        cands.truncate(24);
        for &(_, po, ko, pi, ki, t0, t1, dc) in &cands {
            // the outer member bumps away from the inner one
            let side = -dc.signum();
            let d = dc.abs();
            let g = ((3.0 * w0).max(clear) + w0).max(2.0 * d + w0 + gap_m);
            let margin = g;
            let room = t1 - t0 - 2.0 * margin;
            let mut amps = Vec::new();
            let mut a = (g / 2.0).max(2.0 * w0);
            while a < max_amp {
                amps.push(a);
                a *= 1.4;
            }
            amps.push(max_amp.max(g / 2.0));
            let mut tries: Vec<(usize, f64)> = Vec::new();
            for &amp in &amps {
                let c = chamfer(g, amp);
                let nb_fit = (((room - g - 2.0 * c) / (2.0 * g)).floor() + 1.0).max(0.0) as usize;
                let per = bump_gain(amp, c);
                if nb_fit == 0 || per <= 1e-3 {
                    continue;
                }
                let nb = (want / per).ceil().max(1.0) as usize;
                if nb <= nb_fit {
                    tries.push((nb, amp_for(want / nb as f64, g).min(max_amp)));
                } else {
                    for x in [nb_fit, nb_fit / 2, 1] {
                        if x >= 1 {
                            tries.push((x, amp));
                        }
                    }
                }
            }
            for &(nb, amp) in &tries {
                let c = chamfer(g, amp);
                let span = (nb - 1) as f64 * 2.0 * g + g + 2.0 * c;
                let free = room - span;
                if free < -1e-9 {
                    continue;
                }
                for off in [0.5, 0.0, 1.0] {
                    let (a, b) = { let pl = &self.w.pieces[po as usize].poly; (pl[ko], pl[ko + 1]) };
                    let start = t0 + margin + c + free * off;
                    let outer = serpentine(a, b, start, nb, g, amp, side);
                    // the inner member: the bumps' offset onto its line, in its own direction
                    let mut follow = offset_polyline(&outer[1..outer.len() - 1], -side * d);
                    let (ci, ei) = { let pl = &self.w.pieces[pi as usize].poly; (pl[ki], pl[ki + 1]) };
                    if (ei - ci).dot(b - a) < 0.0 {
                        follow.reverse();
                    }
                    let li = ci.dist(ei);
                    let ui = (ei - ci).unit();
                    let (f0, f1) = ((follow[0] - ci).dot(ui), (follow[follow.len() - 1] - ci).dot(ui));
                    if f0 < 0.05 || f1 > li - 0.05 || f0 >= f1 {
                        continue;
                    }
                    let old = self.w.pieces[po as usize].poly.clone();
                    let mut np = old.clone();
                    np.splice(ko..=ko + 1, outer);
                    if !self.w.try_set_poly(po, np) {
                        continue;
                    }
                    let mut ni = self.w.pieces[pi as usize].poly.clone();
                    let mut seg = vec![ci];
                    seg.extend(follow);
                    seg.push(ei);
                    ni.splice(ki..=ki + 1, seg);
                    if self.w.try_set_poly(pi, ni) {
                        return Some(nb);
                    }
                    self.w.set_poly(po, old);
                }
            }
        }
        None
    }

    /// One serpentine on the net, adding up to `want`: the bumps added, or None when nothing fits.
    fn meander_step(&mut self, net: u32, want: f64) -> Option<usize> {
        let pads = self.net_pads[net as usize].clone();
        let chain = self.piece_chain(net, pads[0], pads[1])?;
        let partner = self.w.b.partner.get(net as usize).cloned().unwrap_or(NONE);
        let w0 = chain.iter().map(|&(p, _)| self.w.pieces[p as usize].width).fold(0.0, f64::max);
        let clear = self.w.b.classes[self.w.b.nets[net as usize].class].clear;
        let g = self.w.b.meander_gap.unwrap_or((3.0 * w0).max(clear)) + w0; // leg pitch, centre to centre
        let max_amp = self.w.b.meander_max_amp.max(g);
        // the straight stretches, longest first; a bump stays a leg pitch off a stretch's ends (its corners)
        let margin = g;
        let mut segs: Vec<(f64, u32, usize)> = Vec::new();
        for &(pc, _) in &chain {
            let pl = &self.w.pieces[pc as usize].poly;
            for k in 0..pl.len().saturating_sub(1) {
                let len = pl[k].dist(pl[k + 1]);
                if len >= 2.0 * margin + 1.5 * g {
                    segs.push((len, pc, k));
                }
            }
        }
        segs.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        segs.truncate(40);
        // bump heights from low (little room needed) to the limit
        let mut amps = Vec::new();
        let mut a = (g / 2.0).max(2.0 * w0);
        while a < max_amp {
            amps.push(a);
            a *= 1.4;
        }
        amps.push(max_amp);
        // a pair member bumps away from its partner where that runs beside it
        let sides_of = |r: &Router, pc: u32, k: usize| -> Vec<f64> {
            if partner == NONE {
                return vec![1.0, -1.0];
            }
            let pcs = &r.w.pieces[pc as usize];
            let (a, b) = (pcs.poly[k], pcs.poly[k + 1]);
            let (dir, mid) = ((b - a).unit(), a.lerp(b, 0.5));
            let mut toward = 0.0;
            for q in r.w.pieces.iter().filter(|x| x.alive && x.net == partner && x.layer == pcs.layer) {
                for s in q.poly.windows(2) {
                    let f = s[0].lerp(s[1], seg_param(s[0], s[1], mid));
                    if f.dist(mid) < 1.0 {
                        toward += dir.cross(f - mid).signum();
                    }
                }
            }
            if toward > 0.0 {
                vec![-1.0]
            } else if toward < 0.0 {
                vec![1.0]
            } else {
                vec![1.0, -1.0]
            }
        };
        // (whole amount in one stretch first; else as much as one stretch takes, the most first)
        for whole in [true, false] {
            for &(len, pc, k) in &segs {
                let room = len - 2.0 * margin;
                let mut tries: Vec<(usize, f64)> = Vec::new();
                for &amp in &amps {
                    let c = chamfer(g, amp);
                    let nb_fit = (((room - g - 2.0 * c) / (2.0 * g)).floor() + 1.0).max(0.0) as usize;
                    if nb_fit == 0 {
                        continue;
                    }
                    let per = bump_gain(amp, c);
                    if per <= 1e-3 {
                        continue;
                    }
                    if whole {
                        let nb = (want / per).ceil().max(1.0) as usize;
                        if nb <= nb_fit {
                            let exact = amp_for(want / nb as f64, g).min(max_amp);
                            tries.push((nb, exact));
                        }
                    } else {
                        // shorter runs too: they fit between what crowds a part of the stretch
                        let mut nbs = vec![nb_fit, nb_fit * 2 / 3, nb_fit / 2, nb_fit / 3, 2, 1];
                        nbs.retain(|&x| x >= 1 && x <= nb_fit && x as f64 * per <= want + 1e-6);
                        tries.extend(nbs.into_iter().map(|x| (x, amp)));
                    }
                }
                if !whole {
                    // the most length first
                    tries.sort_by(|x, y| (y.0 as f64 * bump_gain(y.1, chamfer(g, y.1))).partial_cmp(&(x.0 as f64 * bump_gain(x.1, chamfer(g, x.1)))).unwrap());
                }
                tries.dedup_by(|x, y| x.0 == y.0 && (x.1 - y.1).abs() < 1e-6);
                for side in sides_of(self, pc, k) {
                    for &(nb, amp) in &tries {
                        let c = chamfer(g, amp);
                        let span = (nb - 1) as f64 * 2.0 * g + g + 2.0 * c;
                        let free = room - span;
                        if free < -1e-9 {
                            continue;
                        }
                        let (a, b) = {
                            let pl = &self.w.pieces[pc as usize].poly;
                            (pl[k], pl[k + 1])
                        };
                        for off in [0.5, 0.0, 1.0, 0.25, 0.75] {
                            let start = margin + c + free * off;
                            let pts = serpentine(a, b, start, nb, g, amp, side);
                            let mut np = self.w.pieces[pc as usize].poly.clone();
                            np.splice(k..=k + 1, pts);
                            if self.w.try_set_poly(pc, np) {
                                return Some(nb);
                            }
                        }
                    }
                }
            }
        }
        None
    }
}
