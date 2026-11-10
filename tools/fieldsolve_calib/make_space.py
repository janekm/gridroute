"""MTH R7A PIFA geometry search space for `fieldsolve optimise` (docs/fieldsolve-fast-mode.md, Stage 3).

    python3 make_space.py base.json        (base.json from make_specs.py, in the current directory)

Writes mth_param.json (base.json with parameterised edits) and space.json. Parameters (mm): arm lengthens
(> 0) or trims (< 0) the radiator's open end at x 141.9; feed moves the feed strip (x 153.95-154.45) towards
the open end (< 0) or the short (> 0); clr widens (> 0) or narrows (< 0) the ground clearance under the
radiator (edge at y 76.0). Each expression gives an empty rectangle at 0, so the start point is as built."""
import json, sys
b = json.load(open(sys.argv[1]))
b.pop('out', None)
b.pop('edits', None)
F, B = 'F.Cu', 'B.Cu'
def ed(layer, op, r): return {"layer": layer, "op": op, "rect": r, "net": "GND"}
yc0, yc1 = "=76.0 - min(max(clr, 0), 0.1)", "=76.0 + max(clr, 0)"
ya0, ya1 = "=76.0 + min(clr, 0)", "=76.0 + min(max(-clr, 0), 0.05)"
b['edits'] = [
    ed(F, 'clear', ["=141.9 - min(max(-arm, 0), 0.1)", 74.5, "=141.9 + max(-arm, 0)", 75.1]),
    ed(F, 'add', ["=141.9 - max(arm, 0)", 74.64, "=141.9 + min(max(arm, 0), 0.05)", 74.94]),
    ed(F, 'clear', [153.9, 71.05, 154.5, 74.45]),
    ed(F, 'add', ["=153.95 + feed", 70.9, "=154.45 + feed", 74.95]),
    ed(F, 'add', ["=min(153.95 + feed, 153.95)", 74.45, "=max(154.45 + feed, 154.45)", 74.95]),
    ed(F, 'clear', [139.9, yc0, 152.3, yc1]), ed(F, 'clear', [154.7, yc0, 155.2, yc1]), ed(F, 'clear', [156.2, yc0, 160.1, yc1]),
    ed(B, 'clear', [139.9, yc0, 160.1, yc1]),
    ed(F, 'add', [140.2, ya0, 153.3, ya1]), ed(F, 'add', [156.15, ya0, 159.8, ya1]), ed(B, 'add', [140.2, ya0, 159.8, ya1]),
]
json.dump(b, open('mth_param.json', 'w'), indent=1)
json.dump({"spec": "mth_param.json", "out": "surrogate", "cache": "cache", "method": "surrogate", "budget": 16, "seed": 1,
           "parameters": [{"name": "arm", "min": -1.9, "max": 1.5, "step": 0.1, "start": 0},
                          {"name": "feed", "min": -1.0, "max": 0.5, "step": 0.1, "start": 0},
                          {"name": "clr", "min": -0.5, "max": 1.0, "step": 0.1, "start": 0}]}, open('space.json', 'w'), indent=1)
