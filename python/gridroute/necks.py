"""Narrowing of a power net: how wide the copper is on the way from a source to every pad of the net.

A pour that looks wide can be pinched to a sliver by the antipads of the vias and holes that pass through it, and a
trunk can meet its pour through less copper than either has. A clearance check does not see this: the copper is all
connected. This module rasterises the net's copper on every layer (zones less the clearance round other nets'
copper and holes, plus the net's own tracks, vias and pads), and finds for each pad the widest path from the source:
the largest width w such that a disc of diameter w can slide from the source to the pad without leaving the copper.

    from gridroute import necks
    r = necks.analyse(bd, '+5V_DUT', source=['C1059', 'C1060'], layers=CU, planes={'+3V3': 'In4.Cu'}, outline=poly)
    print(necks.report(r, below=1.0))
    necks.render(r, 'dut5.png')

Widths are found in steps (STEPS), so a result reads "at least this wide". Layers join at the net's vias and plated
holes; a via counts as a conductor as wide as its barrel's circumference. Zone connections are taken as solid (no
thermal spokes). The neck position is where the copper that only this width reaches hangs on the wider copper nearer the
source (the contact nearest the pad, if there are several). Raster estimate: about one pixel (pitch) of error on each edge.

Needs numpy, scipy and Pillow.
"""
import math

import numpy as np

STEPS = (3.0, 2.5, 2.0, 1.5, 1.2, 1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1)


def _pads(bd):
    for part in bd.parts.values():
        if getattr(part, 'placed', True):
            yield from part.pads()


def _pad_layers(p, layers):
    """Copper layers of a pad: gridroute pads carry 'F', 'B' or 'all' (layers: top first, bottom last)."""
    if p.get('npth'):
        return []
    pl = p['layers']
    if p.get('drill', 0) > 0 or pl == 'all':
        return list(layers)
    if pl in ('F', 'B'):
        return [layers[0] if pl == 'F' else layers[-1]]
    return [l for l in pl if l in layers]


class _Raster:
    def __init__(self, x0, y0, x1, y1, pitch):
        from PIL import Image, ImageDraw
        self.x0, self.y0, self.s = x0, y0, 1.0 / pitch
        self.nx, self.ny = int(math.ceil((x1 - x0) * self.s)) + 1, int(math.ceil((y1 - y0) * self.s)) + 1
        self.im = Image.new('L', (self.nx, self.ny), 0)
        self.d = ImageDraw.Draw(self.im)

    def P(self, x, y):
        return ((x - self.x0) * self.s, (y - self.y0) * self.s)

    def disc(self, x, y, r, fill):
        if r > 0:
            cx, cy = self.P(x, y)
            r *= self.s
            self.d.ellipse([cx - r, cy - r, cx + r, cy + r], fill=fill)

    def line(self, pts, w, fill):
        """Polyline with round ends and joints, w wide."""
        q = [self.P(*p) for p in pts]
        if len(q) > 1:
            self.d.line(q, fill=fill, width=max(int(round(w * self.s)), 1))
        for p in pts:
            self.disc(p[0], p[1], w / 2, fill)

    def poly(self, pts, fill, grow=0.0):
        if len(pts) >= 3:
            self.d.polygon([self.P(*p) for p in pts], fill=fill)
        if grow > 0:
            self.line(list(pts) + [pts[0]], 2 * grow, fill)

    def pad(self, p, fill, grow=0.0):
        if p.get('polygons'):
            for g in p['polygons']:
                self.poly(g['points'], fill, grow)
            return
        w, h = p['w'] + 2 * grow, p['h'] + 2 * grow
        m = min(p['w'], p['h'])
        r = {'circle': m / 2, 'oval': m / 2, 'roundrect': p.get('rr', 0.25) * m}.get(p['shape'], 0.0) + grow
        a = math.radians(p.get('angle', 0.0))
        if abs(a) > 1e-9 or r <= 0:
            c, s = math.cos(a), math.sin(a)
            self.poly([(p['x'] + dx * c - dy * s, p['y'] + dx * s + dy * c)
                       for dx, dy in ((-w / 2, -h / 2), (w / 2, -h / 2), (w / 2, h / 2), (-w / 2, h / 2))], fill)
            return
        x0, y0 = self.P(p['x'] - w / 2, p['y'] - h / 2)
        x1, y1 = self.P(p['x'] + w / 2, p['y'] + h / 2)
        self.d.rounded_rectangle([x0, y0, x1, y1], radius=min(r * self.s, (x1 - x0) / 2, (y1 - y0) / 2), fill=fill)

    def array(self):
        return np.asarray(self.im) > 0


def analyse(bd, net, source, layers, planes=None, outline=None, pitch=0.05, clearance=0.2, hole_clearance=0.25,
            edge_clearance=0.3, steps=STEPS, window=None, reach=1.5, via_width=None):
    """Widest path from the source to every pad of net.

    bd: a gridroute Board (zones, tracks, vias, parts). source: refs, or (ref, pin) pairs, of the pads the supply
    comes from. layers: every copper layer name. planes: {net: layer} of the nets that fill a whole layer inside
    outline (a polygon). clearance / hole_clearance: zone fill to other nets' copper / to holes. window: (x0, y0,
    x1, y1) to limit the area; default the net's copper. reach: a pad (source or load) counts as joined once copper of the width
    comes this near it on any layer, so that the short drop from a pour to each pad (a via and a stub) does not hide
    the pours' and trunks' own widths; 0 follows the copper right into every pad. via_width: the width a via or
    plated hole counts as (default: its barrel's circumference); a large value shows the pours alone, as if every
    drop had vias enough.

    Returns a dict: pads (one per pad: ref, pin, x, y, width, neck (x, y, layer) or None, via (the neck is a via),
    source), layers, and the raster (for render)."""
    from scipy import ndimage
    planes = planes or {}
    layers = list(layers)
    pads = [p for p in _pads(bd) if not p.get('npth')]
    mine = [p for p in pads if p['net'] == net]
    tracks = [t for t in bd.tracks if t['layer'] in layers]
    zones = [z for z in getattr(bd, 'zones', ()) if z['layer'] in layers]
    plane_layer = planes.get(net)
    if plane_layer and not outline:
        raise ValueError('a plane net needs the board outline')
    # window: everything of the net
    xs, ys = [], []
    for p in mine:
        xs += [p['x'] - p['w'], p['x'] + p['w']]
        ys += [p['y'] - p['h'], p['y'] + p['h']]
    for t in tracks:
        if t['net'] == net:
            xs += [q[0] for q in t['pts']]
            ys += [q[1] for q in t['pts']]
    for z in zones:
        if z['net'] == net:
            xs += [q[0] for q in z['pts']]
            ys += [q[1] for q in z['pts']]
    if plane_layer:
        xs += [q[0] for q in outline]
        ys += [q[1] for q in outline]
    if not xs:
        raise ValueError('no copper on net %r' % net)
    if window is None:
        window = (min(xs) - 1.0, min(ys) - 1.0, max(xs) + 1.0, max(ys) + 1.0)
    x0, y0, x1, y1 = window
    inside = lambda x, y, m=3.0: x0 - m <= x <= x1 + m and y0 - m <= y <= y1 + m

    masks = {}
    for L in layers:
        R = _Raster(x0, y0, x1, y1, pitch)
        own = [z for z in zones if z['layer'] == L and z['net'] == net]
        others = [z for z in zones if z['layer'] == L and z['net'] != net]
        is_plane = plane_layer == L
        if is_plane:
            R.poly(outline, 255)
            R.line(list(outline) + [outline[0]], 2 * edge_clearance, 0)
        for z in own:
            R.poly(z['pts'], 255)
        if is_plane or own:
            top = max([z.get('priority', 1) for z in own], default=-1)
            for z in others:                       # a plane yields to every pour; a pour to pours of higher priority
                if is_plane or z.get('priority', 1) > top:
                    R.poly(z['pts'], 0, clearance)
            for v in bd.vias:
                if v['net'] != net and inside(v['x'], v['y']):
                    R.disc(v['x'], v['y'], max(v['d'] / 2 + clearance, v['drill'] / 2 + hole_clearance), 0)
            for p in _pads(bd):
                if p['net'] == net and not p.get('npth') or not inside(p['x'], p['y'], 6.0):
                    continue
                if p.get('npth') or p.get('drill', 0) > 0:
                    dw, dh = p.get('drill_w', p['drill']), p.get('drill_h', p['drill'])
                    hole = dict(p, w=dw, h=dh, shape='oval', polygons=None)
                    R.pad(hole, 0, hole_clearance)
                if not p.get('npth') and L in _pad_layers(p, layers):
                    R.pad(p, 0, clearance)
            for t in tracks:
                if t['net'] != net and t['layer'] == L:
                    R.line(t['pts'], t['width'] + 2 * clearance, 0)
        for t in tracks:
            if t['net'] == net and t['layer'] == L:
                R.line(t['pts'], t['width'], 255)
        for v in bd.vias:
            if v['net'] == net:
                R.disc(v['x'], v['y'], v['d'] / 2, 255)
        for p in mine:
            if L in _pad_layers(p, layers):
                R.pad(p, 255)
        m = R.array()
        if m.any():
            masks[L] = m
    used = [L for L in layers if L in masks]
    edt = {L: ndimage.distance_transform_edt(masks[L]) * pitch for L in used}
    ny, nx = next(iter(masks.values())).shape
    s = 1.0 / pitch

    def box(x, y, hw, hh):
        cx, cy = (x - x0) * s, (y - y0) * s
        return (slice(max(int(cy - hh * s) - 1, 0), min(int(cy + hh * s) + 3, ny)),
                slice(max(int(cx - hw * s) - 1, 0), min(int(cx + hw * s) + 3, nx)))

    # layer links: the net's vias and plated holes
    links = [(v['x'], v['y'], v['d'] / 2, via_width or math.pi * v['drill'], used) for v in bd.vias if v['net'] == net]
    for p in mine:
        pl = [L for L in _pad_layers(p, layers) if L in masks]
        if p.get('drill', 0) > 0 and len(pl) > 1:
            links.append((p['x'], p['y'], min(p['w'], p['h']) / 2, via_width or math.pi * p['drill'], pl))
    links = [(box(x, y, r, r), cap, ls, (x, y)) for x, y, r, cap, ls in links if inside(x, y, 0.0)]

    def is_source(p):
        for q in source:
            if isinstance(q, str):
                if p['ref'] == q:
                    return True
            elif (p['ref'], str(p['num'])) == (q[0], str(q[1])):
                return True
        return False

    terms = []
    for p in mine:
        pl = [L for L in _pad_layers(p, layers) if L in masks]
        if pl and inside(p['x'], p['y'], 0.0):
            terms.append(dict(ref=p['ref'], pin=str(p['num']), x=p['x'], y=p['y'], width=0.0, neck=None, via=False,
                              source=is_source(p), _layers=pl if reach <= 0 else used,
                              _box=box(p['x'], p['y'], p['w'] / 2 + reach, p['h'] / 2 + reach)))
    if not any(t['source'] for t in terms):
        raise ValueError('no source pad of %r among %r' % (net, source))

    prev_mask, prev_src, prev_w, pending = None, None, None, [t for t in terms if not t['source']]
    for t in terms:
        if t['source']:
            t['width'] = float('inf')
    for w in sorted(steps, reverse=True):
        if not pending:
            break
        lab, base, n = {}, {}, 1
        cur = {L: edt[L] >= w / 2 - pitch * 0.5 for L in used}
        for L in used:
            lab[L], k = ndimage.label(cur[L])
            base[L] = n
            n += k
        parent = np.arange(n)

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        def ids(L, bx):
            u = np.unique(lab[L][bx])
            return [int(i) + base[L] - 1 for i in u if i]

        via_of = {}
        for bx, cap, ls, xy in links:
            if w > cap + 1e-9:
                continue
            got = [i for L in ls for i in ids(L, bx)]
            for i in got[1:]:
                a, b = find(got[0]), find(i)
                if a != b:
                    parent[a] = b
            for i in got:
                via_of.setdefault(i, xy)
        roots = set()
        for t in terms:
            if t['source']:
                roots |= {find(i) for L in t['_layers'] for i in ids(L, t['_box'])}
        if not roots:
            prev_mask, prev_src, prev_w = cur, {L: np.zeros_like(cur[L]) for L in used}, w
            continue
        root_of = np.array([find(i) for i in range(n)])
        in_src = np.isin(root_of, list(roots))
        src_px = {L: in_src[np.where(lab[L] > 0, lab[L] + base[L] - 1, 0)] & (lab[L] > 0) for L in used}
        still, joined = [], []
        for t in pending:
            if any(src_px[L][t['_box']].any() for L in t['_layers']):
                t['width'] = w
                joined.append(t)
            else:
                still.append(t)
        if joined and prev_src is not None:
            # what this step added to the source's copper, in pieces; each piece hangs on the wider copper of the
            # step above through its contacts: those are the necks of this width
            rl, rb, m = {}, {}, 1
            # the step also grows the wider copper by a band of half the width difference all round: not a neck
            band = (prev_w - w) / 2 / pitch + 1.5
            near = {L: (ndimage.distance_transform_edt(~prev_src[L]) <= band) if prev_src[L].any()
                    else np.zeros_like(prev_src[L]) for L in used}
            for L in used:
                rl[L], k = ndimage.label(src_px[L] & ~near[L])
                rb[L] = m
                m += k
            par = list(range(m))

            def fnd(i):
                while par[i] != i:
                    par[i] = par[par[i]]
                    i = par[i]
                return i

            contacts = {}
            for bx, cap, ls, xy in links:
                if w > cap + 1e-9:
                    continue
                got = [int(i) + rb[L] - 1 for L in ls for i in np.unique(rl[L][bx]) if i]
                for i in got[1:]:
                    par[fnd(got[0])] = fnd(i)
                if got and any(prev_src[L][bx].any() for L in ls):
                    contacts.setdefault(fnd(got[0]), []).append((xy[0], xy[1], 'via'))
            for L in used:
                yy, xx = np.nonzero((rl[L] > 0) & ndimage.binary_dilation(near[L], iterations=2))
                for i, x_, y_ in zip(rl[L][yy, xx], xx, yy):
                    contacts.setdefault(int(i) + rb[L] - 1, []).append((x_ / s + x0, y_ / s + y0, L))
            byroot = {}
            for i, c in contacts.items():
                byroot.setdefault(fnd(i), []).extend(c)
            for t in joined:
                best = None
                for L in t['_layers']:
                    for i in np.unique(rl[L][t['_box']]):
                        for c in byroot.get(fnd(int(i) + rb[L] - 1), ()) if i else ():
                            d2 = (c[0] - t['x']) ** 2 + (c[1] - t['y']) ** 2
                            if best is None or d2 < best[0]:
                                best = (d2, c)
                if best:
                    c = best[1]
                    t['neck'] = (round(float(c[0]), 2), round(float(c[1]), 2), c[2])
                    t['via'] = c[2] == 'via'
        pending = still
        prev_mask, prev_src, prev_w = cur, src_px, w
    for t in terms:
        t.pop('_box')
        t.pop('_layers')
    return dict(net=net, pads=terms, layers=used, window=window, pitch=pitch, masks=masks, edt=edt, steps=tuple(steps))


def report(result, below=None, loads=None):
    """Text table: the pads grouped by neck, narrowest first. below: only the necks narrower than this. loads: only
    the pads of these refs (the parts that draw the current)."""
    groups = {}
    for t in result['pads']:
        if t['source'] or (below is not None and t['width'] >= below) or (loads and t['ref'] not in loads):
            continue
        key = (t['width'], t['neck'])
        groups.setdefault(key, []).append('%s.%s' % (t['ref'], t['pin']))
    out = ['%s: %d pads, %d source' % (result['net'], len(result['pads']), sum(t['source'] for t in result['pads']))]
    for (w, neck), names in sorted(groups.items(), key=lambda kv: (kv[0][0], str(kv[0][1]))):
        where = 'at (%.2f, %.2f) %s' % neck if neck else ''
        head = 'open' if w == 0 else '>= %.2f mm' % w
        out.append('  %-10s %-28s %d pads: %s%s' % (head, where, len(names), ' '.join(sorted(names)[:8]),
                                                    ' ...' if len(names) > 8 else ''))
    return '\n'.join(out)


def local_width(result, layer):
    """Per pixel, the width of the widest disc of STEPS that fits in the copper there (0 off the copper)."""
    from scipy import ndimage
    e, pitch = result['edt'][layer], result['pitch']
    out = np.zeros(e.shape, np.float32)
    out[result['masks'][layer]] = 0.05
    for w in sorted(result['steps']):
        core = e >= w / 2 - pitch * 0.5
        if not core.any():
            break
        out[ndimage.distance_transform_edt(~core) * pitch <= w / 2] = w
    return out


def render(result, path, layers=None, scale=None, good=2.0, bad=0.3, window=None):
    """PNG of the net's copper coloured by its local width (red at bad and below, green at good and above), the
    layers side by side, with each neck ringed and labelled with its width and each source pad marked."""
    from PIL import Image, ImageDraw
    layers = [L for L in (layers or result['layers']) if L in result['masks']]
    x0, y0, x1, y1 = result['window']
    pitch = result['pitch']
    wx0, wy0, wx1, wy1 = window or result['window']
    c0, c1 = int((wx0 - x0) / pitch), int((wx1 - x0) / pitch)
    r0, r1 = int((wy0 - y0) / pitch), int((wy1 - y0) / pitch)
    tiles = []
    for L in layers:
        w = local_width(result, L)[max(r0, 0):r1, max(c0, 0):c1]
        t = np.clip((np.log(np.maximum(w, 1e-3)) - math.log(bad)) / (math.log(good) - math.log(bad)), 0, 1)
        rgb = np.zeros(w.shape + (3,), np.uint8)
        rgb[..., 0] = (235 * np.minimum(1, 2 - 2 * t)).astype(np.uint8)
        rgb[..., 1] = (200 * np.minimum(1, 2 * t)).astype(np.uint8)
        rgb[..., 2] = 40
        rgb[w == 0] = (20, 20, 26)
        im = Image.fromarray(rgb)
        if scale:
            im = im.resize((max(int(im.size[0] * scale), 1), max(int(im.size[1] * scale), 1)), Image.LANCZOS)
        k = (scale or 1.0) / pitch
        d = ImageDraw.Draw(im)
        seen = set()
        for p in result['pads']:
            if p['source']:
                cx, cy = (p['x'] - wx0) * k, (p['y'] - wy0) * k
                d.rectangle([cx - 4, cy - 4, cx + 4, cy + 4], outline=(90, 170, 255), width=2)
            n = p['neck']
            if n and (n[2] == L or n[2] == 'via') and (n, p['width']) not in seen and p['width'] < good:
                seen.add((n, p['width']))
                cx, cy = (n[0] - wx0) * k, (n[1] - wy0) * k
                r = max(0.8 * k, 7)
                d.ellipse([cx - r, cy - r, cx + r, cy + r], outline=(255, 255, 255), width=2)
                d.text((cx + r + 2, cy - 6), '%.2f' % p['width'], fill=(255, 255, 255))
        d.text((4, 2), '%s  %s' % (result['net'], L), fill=(255, 255, 255))
        tiles.append(im)
    W = sum(t.size[0] for t in tiles) + 6 * (len(tiles) - 1)
    out = Image.new('RGB', (W, max(t.size[1] for t in tiles)), (60, 60, 70))
    x = 0
    for t in tiles:
        out.paste(t, (x, 0))
        x += t.size[0] + 6
    out.save(path)
    return path
