"""MTH R7A PIFA geometry candidates for the fast-mode calibration (docs/fieldsolve-fast-mode.md): antenna
alone at the feed (port C20, R70 and the ladder left open), whole board, roi widened to cover the radiator's
open end. Writes one spec per candidate into the current directory.

    python3 make_specs.py MTH_7PORT_SPEC.json GERBER_DIR

MTH_7PORT_SPEC is the all-port MTH run spec of docs/mth-antenna-findings.md; GERBER_DIR is relative to the
output directory."""
import json, copy, sys
base = json.load(open(sys.argv[1]))
base['gerbers'] = sys.argv[2]
base['roi'] = [141.5, 70.0, 159.6, 87.5]
base['ports'] = [p for p in base['ports'] if p['name'] == 'C20']
base.pop('fields', None)
base['farfield'] = {"freqs": [2.40e9, 2.44e9, 2.48e9], "step_deg": 5, "gap": 3, "tol": 0.01}
base['objective'] = {"port": "C20", "band": [2.40e9, 2.44e9, 2.48e9]}
# feed, feed-to-short gap and the short's junction with the ground edge keep the full resolution
base['fidelity'] = {"protect": [[153.8, 75.4, 156.3, 77.0]]}
F, B = 'F.Cu', 'B.Cu'
def ed(layer, op, r): return {"layer": layer, "op": op, "rect": r, "net": "GND"}
def feed(dx):
    x0 = 153.95 + dx
    lo, hi = min(x0, 153.95), max(x0 + 0.5, 154.45)
    return [ed(F, 'clear', [153.9, 71.05, 154.5, 74.45]), ed(F, 'add', [x0, 70.9, x0 + 0.5, 74.95]), ed(F, 'add', [lo, 74.45, hi, 74.95])]
def gnd_clear(y):
    return [ed(F, 'clear', [139.9, 75.9, 152.3, y]), ed(F, 'clear', [154.7, 75.9, 155.2, y]), ed(F, 'clear', [156.2, 75.9, 160.1, y]), ed(B, 'clear', [139.9, 75.9, 160.1, y])]
cands = {
    'base': [],
    'trim1': [ed(F, 'clear', [141.8, 74.5, 142.9, 75.1])],
    'trim2': [ed(F, 'clear', [141.8, 74.5, 143.84, 75.1])],
    'ext1': [ed(F, 'add', [140.9, 74.64, 141.95, 74.94])],
    'ext15': [ed(F, 'add', [140.4, 74.64, 141.95, 74.94])],
    'feedL05': feed(-0.5),
    'feedL10': feed(-1.0),
    'feedR03': feed(0.3),
    'gnd+05': gnd_clear(76.5),
    'gnd+10': gnd_clear(77.0),
    'gnd-05': [ed(F, 'add', [140.2, 75.5, 153.3, 76.05]), ed(F, 'add', [156.15, 75.5, 159.8, 76.05]), ed(B, 'add', [140.2, 75.5, 159.8, 76.05])],
    'trim1_feedL05': [ed(F, 'clear', [141.8, 74.5, 142.9, 75.1])] + feed(-0.5),
}
for name, edits in cands.items():
    s = copy.deepcopy(base)
    if edits:
        s['edits'] = edits
    s['out'] = f'runs/{name}'
    json.dump(s, open(f'{name}.json', 'w'), indent=1)
print(' '.join(cands))
