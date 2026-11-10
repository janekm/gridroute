#!/usr/bin/env python3
"""Draw a layout.json over its model's pads as SVG (tracks coloured by layer, vias as rings; with a check.json, open
nets' pads outlined and overlapping copper ringed in yellow).

    render_layout.py MODEL.json LAYOUT.json OUT.svg [--check CHECK.json] [--title TEXT] [--px 900]
"""
import argparse
import json
import math
import pathlib

COLOURS = ['#d62728', '#1f77b4', '#2ca02c', '#9467bd', '#ff7f0e', '#17becf', '#8c564b', '#e377c2']


def main():
    p = argparse.ArgumentParser()
    p.add_argument('model', type=pathlib.Path)
    p.add_argument('layout', type=pathlib.Path)
    p.add_argument('out', type=pathlib.Path)
    p.add_argument('--check', type=pathlib.Path)
    p.add_argument('--title', default='')
    p.add_argument('--px', type=float, default=900, help='pixels along the longer side')
    a = p.parse_args()
    m = json.loads(a.model.read_text())
    lay = json.loads(a.layout.read_text())
    open_nets, overlaps = set(), []
    if a.check and a.check.exists():
        chk = json.loads(a.check.read_text())
        open_nets = {x['net'] for x in chk['incomplete']}
        overlaps = chk.get('overlap_list', [])
    layers = m['layers']
    col = {l: COLOURS[i % len(COLOURS)] for i, l in enumerate(layers)}
    xs, ys = [q[0] for q in m['outline']], [q[1] for q in m['outline']]
    fs = max(max(xs) - min(xs) + 2, max(ys) - min(ys) + 5) * 0.025  # title size, from the longer side
    x0, y0, x1, y1 = min(xs) - 1, min(ys) - max(4, 2.8 * fs), max(xs) + 1, max(ys) + 1
    s = a.px / max(x1 - x0, y1 - y0)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x0} {y0} {x1 - x0} {y1 - y0}" '
           f'width="{(x1 - x0) * s:.0f}" height="{(y1 - y0) * s:.0f}" style="background:#111">',
           f'<polygon points="{" ".join(f"{x},{y}" for x, y in m["outline"])}" fill="#1b1b1b" stroke="#888" stroke-width="{1.5 / s}"/>']
    for pd in m['pads']:
        w, h, ang = pd['w'], pd['h'], pd.get('angle', 0)
        fill = '#555' if pd['layers'] != 'B' else '#3a3a3a'
        stroke = ' stroke="#ff0" stroke-width="%.3f"' % (2.5 / s) if pd.get('net') in open_nets else ''
        if pd['shape'] == 'circle':
            out.append(f'<circle cx="{pd["x"]}" cy="{pd["y"]}" r="{w / 2}" fill="{fill}"{stroke}/>')
        else:
            rx = min(w, h) / 2 if pd['shape'] == 'oval' else pd.get('rr', 0) if pd['shape'] == 'roundrect' else 0
            out.append(f'<rect x="{-w / 2}" y="{-h / 2}" width="{w}" height="{h}" rx="{rx}" fill="{fill}"{stroke} '
                       f'transform="translate({pd["x"]} {pd["y"]}) rotate({-ang})"/>')
    for t in sorted(lay['tracks'], key=lambda t: -layers.index(t['layer'])):
        pts = ' '.join(f'{x:.4f},{y:.4f}' for x, y in t['pts'])
        out.append(f'<polyline points="{pts}" fill="none" stroke="{col[t["layer"]]}" stroke-width="{t["width"]}" '
                   f'stroke-linecap="round" stroke-linejoin="round" opacity="0.85"/>')
    for v in lay['vias']:
        out.append(f'<circle cx="{v["x"]}" cy="{v["y"]}" r="{v["d"] / 2}" fill="#ddd"/>'
                   f'<circle cx="{v["x"]}" cy="{v["y"]}" r="{v["drill"] / 2}" fill="#111"/>')
    for o in overlaps:
        out.append(f'<circle cx="{o["at"][0]}" cy="{o["at"][1]}" r="{18 / s}" fill="none" stroke="#ff0" '
                   f'stroke-width="{3 / s}"/>')
    legend = '  '.join(f'<tspan fill="{col[l]}">{l}</tspan>' for l in layers)
    out.append(f'<text x="{x0 + 0.5}" y="{y0 + fs * 1.2}" font-size="{fs}" fill="#eee" font-family="sans-serif">{a.title}</text>')
    out.append(f'<text x="{x0 + 0.5}" y="{y0 + fs * 2.5}" font-size="{fs * 0.8}" font-family="sans-serif">{legend}</text>')
    out.append('</svg>')
    a.out.write_text('\n'.join(out))


if __name__ == '__main__':
    main()
