#!/usr/bin/env python3
"""Route a prepared case with the topological engine (target/release/topor), optionally validate in KiCad.

The case directory (bench/pcb/cases/CASE, or --input) is copied to OUT/input, so the KiCad export and
acceptance steps read exactly the model the engine routed. Internal rules come first; --kicad adds the
native export, DRC and baseline comparison. RF frontend nets are marked on the copy first (topor_rf.py, KiCad
Python) unless the model already lists them, so the engine keeps each on one layer.
"""
import argparse
import json
import pathlib
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser()
    p.add_argument('case')
    p.add_argument('--out', type=pathlib.Path, required=True)
    p.add_argument('--input', type=pathlib.Path, help='prepared input directory (default bench/pcb/cases/CASE)')
    p.add_argument('--seconds', type=float, default=300)
    p.add_argument('--engine-args', default='', help='extra arguments for topor route')
    p.add_argument('--kicad', action='store_true', help='export, run KiCad DRC and compare with the baseline')
    p.add_argument('--no-rf', action='store_true', help='do not mark RF frontend nets (they may then change layers)')
    p.add_argument('--kicad-python', default='/Applications/KiCad.app/Contents/Frameworks/Python.framework/Versions/Current/bin/python3')
    p.add_argument('--kicad-cli', default='/Applications/KiCad.app/Contents/MacOS/kicad-cli')
    a = p.parse_args()
    if a.out.exists():
        p.error('use a fresh output directory')
    src = a.input or ROOT / 'cases' / a.case
    a.out.mkdir(parents=True)
    shutil.copytree(src, a.out / 'input')
    model = json.loads((a.out / 'input' / 'model.json').read_text())
    if not a.no_rf and 'rf_nets' not in model and pathlib.Path(a.kicad_python).exists():
        r = subprocess.run([a.kicad_python, str(ROOT / 'topor_rf.py'), str(a.out / 'input')], capture_output=True, text=True)
        if r.returncode:
            print('RF frontend detection failed:', r.stderr.strip().splitlines()[-1:], flush=True)
    exe = ROOT.parents[1] / 'target/release/topor'
    cmd = [str(exe), 'route', str(a.out / 'input'), '--out', str(a.out), '--seconds', str(a.seconds)] + a.engine_args.split()
    with (a.out / 'engine.log').open('w') as log:
        code = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
    print((a.out / 'engine.log').read_text().splitlines()[-1] if code == 0 else f'engine failed ({code})', flush=True)
    if code or not a.kicad:
        return code
    steps = [
        [a.kicad_python, str(ROOT / 'kicad_io.py'), 'export', a.case, '--out', str(a.out)],
        [a.kicad_cli, 'pcb', 'drc', '--format', 'json', '--severity-all', '--all-track-errors',
         '-o', str(a.out / 'routed-drc.json'), str(a.out / 'routed.kicad_pcb')],
        [sys.executable, str(ROOT / 'accept.py'), a.case, '--out', str(a.out)],
    ]
    for s in steps:
        r = subprocess.run(s)
        if r.returncode and s is not steps[-1]:
            return r.returncode
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
