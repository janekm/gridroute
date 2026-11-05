#!/usr/bin/env python3
"""topor: read a board's custom design rules (unrouted.kicad_dru) into a prepared model. Run on a case input dir
(plain Python; KiCad Python also works).

The engine takes the rules as data (model.json cold_routing and pad fields); the forms read here are the ones
seen on real boards, each matched exactly, and every other rule is listed as not modelled:

- class track width: `track_width (min W)` on `A.NetClass == 'C'` (optionally `(layer outer|inner)`): W is a hard
  minimum for C on those layers (layer_min_widths, and the preferred width at least W);
- layer track width for all classes but some: `(layer inner|outer)`, `A.NetClass != 'X' && ...`;
- track spacing on a layer: `clearance (min G)` with `(layer ...)` on `(A.Type == 'Track' || B.Type == 'Track')`,
  optionally `&& A.NetClass != 'X' && B.NetClass != 'X'` (track_clearance, track_clearance_except);
- neck-down: `track_width (min W)` on `(A.NetClass == 'C' || ...) && (A.intersectsCourtyard('REF') || ...)`
  (optionally `(layer outer)`): every outer-layer pad of those parts on a net of those classes gets
  escape_width W; the engine lays a straight stub at W from it to where the class width fits (a stub from a pad
  of the part always touches its courtyard). Later rules win, as in KiCad.
- escape clearance: `clearance (min G)` on the same kind of condition: those pads get escape_clear G, and the
  engine's stub from them keeps G (at the class width) out to where the class clearance holds.
  (`insideCourtyard` is KiCad's old name of `intersectsCourtyard`: any part of the item in the courtyard.)

A rule that relaxes a constraint where it is not modelled is left out: the class value then holds there, which is
conservative.
"""
import json, pathlib, re, sys

case = pathlib.Path(sys.argv[1])
m = json.loads((case / 'model.json').read_text())
dru = case / 'unrouted.kicad_dru'
cr = m.setdefault('cold_routing', {})
layers = m['layers']
outer, inner = [layers[0], layers[-1]], layers[1:-1]


def forms(text):
    text = re.sub(r'(?m)^\s*#.*$', '', text)
    tokens = iter(re.findall(r'"(?:\\.|[^"\\])*"|[()]|[^\s()]+', text))

    def parse(t):
        if t == '(':
            out = []
            for t in tokens:
                if t == ')':
                    return out
                out.append(parse(t))
            raise ValueError('unclosed rule')
        return json.loads(t) if t.startswith('"') else t
    return [parse(t) for t in tokens]


def mm(v):
    v = str(v)
    return float(v[:-2]) if v.endswith('mm') else float(v)


def field(rule, key):
    return [x for x in rule if isinstance(x, list) and x and x[0] == key]


classes = m['classes']
lw = cr.setdefault('layer_widths', {})
lm = cr.setdefault('layer_min_widths', {})


def raise_width(c, ls, w):
    pref = classes[c][0] if c in classes else 0.0
    for l in ls:
        lm.setdefault(c, {})[l] = max(lm.get(c, {}).get(l, 0.0), w)
        lw.setdefault(c, {})[l] = max(lw.get(c, {}).get(l, pref), w)


done, skipped = [], []
neck = {}  # (ref, class) -> width, later rules win
esc_clear = {}  # (ref, class) -> clearance
rules = forms(dru.read_text()) if dru.exists() else []
for r in rules:
    if not (isinstance(r, list) and r and r[0] == 'rule'):
        continue
    name = r[1]
    cond = (field(r, 'condition') or [[None, '']])[0][1].strip()
    lay = (field(r, 'layer') or [[None, None]])[0][1]
    ls = {'outer': outer, 'inner': inner, None: layers}.get(lay, [lay] if lay in layers else None)
    cons = field(r, 'constraint')
    if ls is None or len(cons) != 1:
        skipped.append(name)
        continue
    kind = cons[0][1]
    mins = [mm(x[1]) for x in cons[0][2:] if isinstance(x, list) and x[0] == 'min']
    if not mins:
        skipped.append(name)
        continue
    v = mins[0]
    one = re.fullmatch(r"A\.NetClass == '(\w+)'", cond)
    excl = re.fullmatch(r"A\.NetClass != '\w+'( && A\.NetClass != '\w+')*", cond)
    nk = re.fullmatch(r"\(?((?:A\.NetClass == '\w+'(?: \|\| )?)+)\)? && \(?((?:A\.(?:intersects|inside)Courtyard\('\w+'\)(?: \|\| )?)+)\)?", cond)
    sp = re.fullmatch(r"\(A\.Type == 'Track' \|\| B\.Type == 'Track'\)((?: && [AB]\.NetClass != '\w+')*)", cond)
    if kind == 'track_width' and one:
        if one.group(1) in classes:
            raise_width(one.group(1), ls, v)
        done.append(name)
    elif kind == 'track_width' and excl:
        out = set(re.findall(r"'(\w+)'", cond))
        for c in classes:
            if c not in out:
                raise_width(c, ls, v)
        done.append(name)
    elif kind == 'track_width' and nk and set(ls) <= set(outer):
        for c in re.findall(r"'(\w+)'", nk.group(1)):
            for ref in re.findall(r"'(\w+)'", nk.group(2)):
                neck[(ref, c)] = v
        done.append(name)
    elif kind == 'clearance' and nk:
        for c in re.findall(r"'(\w+)'", nk.group(1)):
            for ref in re.findall(r"'(\w+)'", nk.group(2)):
                esc_clear[(ref, c)] = v
        done.append(name)
    elif kind == 'clearance' and sp and lay is not None:
        tc = cr.setdefault('track_clearance', {})
        for l in ls:
            tc[l] = max(tc.get(l, 0.0), v)
        cr['track_clearance_except'] = sorted(set(cr.get('track_clearance_except', [])) | set(re.findall(r"'(\w+)'", sp.group(1))))
        done.append(name)
    else:
        skipped.append(name)

net_class = m.get('net_class', {})
pads = 0
for p in m['pads']:
    key = (p.get('original_ref'), net_class.get(p.get('net') or '', 'Default'))
    if p.get('layers') not in ('F', 'B'):
        continue
    if key in neck:
        p['escape_width'] = neck[key]
        pads += 1
    if key in esc_clear:
        p['escape_clear'] = esc_clear[key]
        pads += 1
m['custom_rules'] = dict(modelled=done, not_modelled=skipped)
(case / 'model.json').write_text(json.dumps(m, indent=2))
print(json.dumps(dict(modelled=done, not_modelled=skipped, neck_pads=pads, track_clearance=cr.get('track_clearance'),
                      layer_min_widths={c: sorted(set(v.values())) for c, v in lm.items()})))
