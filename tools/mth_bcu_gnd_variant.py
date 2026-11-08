"""MTH R7A layout variant for fieldsolve: all VCC3V3 copper off the bottom layer. The U-shaped VCC3V3 zone along
the board edges and the VCC3V3 B.Cu jogs become GND pour (spec "edits"), and the VCC3V3 vias plus the H5.3
header hole are dropped from a copy of the drill file, so nothing shorts 3.3 V to the new ground. The 3.3 V
distribution that replaces them (top-layer tracks, crossings away from the RF corner) is not modelled.
See docs/mth-antenna-findings.md, "Efficiency".

    PCB=.../MTH_R7A.kicad_pcb GBR=mth7_gbr BASE=mth7p_bom.json python3 mth_bcu_gnd_variant.py

Writes <GBR>_gnd3v3/ and mth7r_gnd.json (BASE with the new Gerbers, out and edits), plus a 20-step dry-run
spec mth7r_gnd_dry.json."""
import json, math, os, re, shutil
PCB = os.path.expanduser(os.environ.get('PCB', '~/src/p-hw/easyeda2kicad/7A/MTH/MTH_R7A.kicad_pcb'))
pcb = open(PCB).read()
blocks = re.split(r'\n\t\(', pcb)
# 1. GND fill: the VCC3V3 B.Cu zone (ring) plus its clearance, keeping off the H5 header pads
rects = [(140.0, 81.4, 141.95, 130.0), (140.0, 128.15, 157.1, 130.0), (157.95, 86.35, 159.95, 105.0),
         (158.45, 105.0, 159.95, 119.8), (157.95, 119.8, 159.95, 127.15), (141.6, 92.6, 142.7, 94.1)]
segs = []
for b in blocks:
    if b.startswith('segment') and re.search(r'\(layer "?B.Cu"?\)', b):
        x0, y0, x1, y1 = map(float, re.search(r'\(start ([-\d.]+) ([-\d.]+)\)\s*\(end ([-\d.]+) ([-\d.]+)\)', b).groups())
        w = float(re.search(r'\(width ([\d.]+)\)', b).group(1)); n = re.search(r'\(net "?([^")]*)"?\)', b).group(1)
        segs.append((x0, y0, x1, y1, w, n))
# 2. every VCC3V3 B.Cu track: squares along it, track width + clearance
sq = []
for x0, y0, x1, y1, w, n in segs:
    if n != 'VCC3V3':
        continue
    L = math.hypot(x1 - x0, y1 - y0); k = max(1, int(L / 0.15))
    h = w / 2 + 0.3
    for i in range(k + 1):
        x, y = x0 + (x1 - x0) * i / k, y0 + (y1 - y0) * i / k
        sq.append((x - h, y - h, x + h, y + h))
allr = rects + sq
# conflict check: other-net B.Cu tracks, vias and TH pads touching the fill
def d2r(x, y, rs):
    return min(math.hypot(max(r[0] - x, 0, x - r[2]), max(r[1] - y, 0, y - r[3])) for r in rs)
bad = []
for x0, y0, x1, y1, w, n in segs:
    if n not in ('VCC3V3', 'GND') and min(d2r(x0 + (x1 - x0) * t / 20, y0 + (y1 - y0) * t / 20, allr) for t in range(21)) < w / 2 + 0.05:
        bad.append(('seg', n, x0, y0))
drop = []
for b in blocks:
    if b.startswith('via'):
        x, y = map(float, re.search(r'\(at ([-\d.]+) ([-\d.]+)', b).groups()); n = re.search(r'\(net "?([^")]*)"?\)', b).group(1)
        if n == 'VCC3V3':
            drop.append((x, y))
        elif n != 'GND' and d2r(x, y, allr) < 0.25:
            bad.append(('via', n, x, y))
drop.append((157.57, 116.21))  # H5 pad 3 (VCC3V3), its B.Cu pad joined the ring
print('conflicts:', bad)
# drill file without the VCC3V3 holes
src = os.environ.get('GBR', 'mth7_gbr'); dst = src.rstrip('/') + '_gnd3v3'
os.makedirs(dst, exist_ok=True)
for f in os.listdir(src):
    shutil.copy(os.path.join(src, f), dst)
lines = open(os.path.join(src, 'MTH_R7A-PTH.DRL')).read().split('\n')
out, nd = [], 0
for ln in lines:
    m = re.match(r'X([-\d.]+)Y([-\d.]+)$', ln)
    if m and any(abs(float(m.group(1)) - x) < 0.02 and abs(-float(m.group(2)) - y) < 0.02 for x, y in drop):
        nd += 1
        continue
    out.append(ln)
open(os.path.join(dst, 'MTH_R7A-PTH.DRL'), 'w').write('\n'.join(out))
print(f'dropped {nd} of {len(drop)} VCC3V3 holes; {len(allr)} GND fill rects')
spec = json.load(open(os.environ.get('BASE', 'mth7p_bom.json')))
spec['gerbers'] = dst; spec['out'] = 'mth7r_gnd'
spec['edits'] = [{'layer': 'B.Cu', 'op': 'add', 'net': 'GND', 'rect': [round(v, 4) for v in r]} for r in allr]
json.dump(spec, open('mth7r_gnd.json', 'w'), indent=1)
d = dict(spec, out='mth7r_gnd_dry', max_steps=20); d['scene'] = True; d.pop('frames', None)
json.dump(d, open('mth7r_gnd_dry.json', 'w'), indent=1)
