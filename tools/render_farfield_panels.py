"""Render far-field patterns of several boards side by side on shared scales (one row per frequency).

    python3 tools/render_farfield_panels.py PANELS.json OUT_PREFIX

Needs Blender (`blender` on PATH) and ImageMagick (`magick`). Each panel's patterns come from a
`fieldsolve cosim` output directory written with "farfield": true (farfield_<variant>_<f>GHz.json).
Outputs land next to PANELS.json: <OUT_PREFIX>_<f>MHz.png per frequency, <OUT_PREFIX>_all.png, and a .blend
per panel and frequency.

PANELS.json: [{"dir": cosim output dir, "label": "...", "board": board_top.png,
               "outline": [x0, y0, x1, y1], "thickness": mm, "centre": [x, y, z]}, ...]
"""
import glob
import json
import pathlib
import subprocess
import sys

TOOL = str(pathlib.Path(__file__).with_name('blender_farfield.py'))
# a font file ImageMagick can read (macOS system font; override with FIELDSOLVE_FONT)
FONT = __import__('os').environ.get('FIELDSOLVE_FONT', '/System/Library/Fonts/SFNS.ttf')
RANGE = 20.0
STOPS = ['#000004', '#57106e', '#bc3754', '#f98e09', '#fcffa4']

panels = json.loads(pathlib.Path(sys.argv[1]).read_text())
prefix = sys.argv[2]
here = pathlib.Path(sys.argv[1]).resolve().parent


def patterns(d):
    out = {}
    for f in sorted(glob.glob(f'{d}/farfield_*GHz.json')):
        p = json.loads(pathlib.Path(f).read_text())
        out[round(p['f'] / 1e6)] = (f, p)
    return out


pats = [patterns(p['dir']) for p in panels]
freqs = sorted(set.intersection(*[set(p) for p in pats]))
W = 1400 * len(panels)
bar = here / f'{prefix}_bar.png'
subprocess.run(['magick', *[f'xc:{c}' for c in STOPS], '+append', '-filter', 'triangle', '-resize', f'{W // 2}x30!', str(bar)], check=True)
rows = []
for fm in freqs:
    top = max(max(pp[fm][1]['eirp_per_w_db']) for pp in pats)
    tiles = []
    for k, (pn, pp) in enumerate(zip(panels, pats)):
        path, p = pp[fm]
        raw = here / f'{prefix}_{fm}MHz_{k}_raw.png'
        blend = here / f'{prefix}_{fm}MHz_{k}.blend'
        subprocess.run(['blender', '--background', '--python', TOOL, '--', path, '--board', pn['board'],
                        '--outline', *[str(v) for v in pn['outline']], '--thickness', str(pn['thickness']),
                        '--centre', *[str(v) for v in pn['centre']], '--top', f'{top}', '--range', f'{RANGE}',
                        '--render', str(raw), '--save', str(blend)], check=True, capture_output=True)
        gmax = max(p['eirp_per_w_db'])
        cap = (f"{pn['label']}\n{p['f'] / 1e9:.2f} GHz (H{round(p['f'] / 2.44e9)}): radiated {p['p_rad_db']:+.1f} dB, "
               f"to connector {p['p_conducted_db']:+.1f} dB, peak EIRP {gmax:+.1f} dB (per W into the chip port)")
        tile = here / f'{prefix}_{fm}MHz_{k}.png'
        subprocess.run(['magick', str(raw), '-gravity', 'north', '-background', '#202024', '-splice', '0x90', '-fill', 'white',
                        '-font', FONT, '-pointsize', '27', '-annotate', '+0+12', cap, str(tile)], check=True)
        tiles.append(str(tile))
    legend = here / f'{prefix}_{fm}MHz_legend.png'
    subprocess.run(['magick', '-size', f'{W}x120', 'xc:#202024', str(bar), '-gravity', 'north', '-geometry', '+0+42', '-composite',
                    '-fill', 'white', '-font', FONT, '-pointsize', '24', '-gravity', 'northwest', '-annotate', f'+{W // 4}+10', f'{top - RANGE:+.1f} dB',
                    '-gravity', 'northeast', '-annotate', f'+{W // 4}+10', f'{top:+.1f} dB', '-gravity', 'north', '-annotate', '+0+80',
                    'EIRP per W into the chip port; colour and radius share this scale in every panel of the row (radius 0 below the bottom)', str(legend)], check=True)
    row = here / f'{prefix}_{fm}MHz.png'
    subprocess.run(['magick', *tiles, '+append', str(legend), '-append', str(row)], check=True)
    rows.append(str(row))
subprocess.run(['magick', *rows, '-append', '-resize', '50%', str(here / f'{prefix}_all.png')], check=True)
print('\n'.join(rows + [str(here / f'{prefix}_all.png')]))
