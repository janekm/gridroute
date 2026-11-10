#!/usr/bin/env python3
"""Route a prepared case with the topological engine (target/release/topor), optionally validate in KiCad.

The case directory (bench/pcb/cases/CASE, or --input) is copied to OUT/input, so the KiCad export and
acceptance steps read exactly the model the engine routed. Internal rules come first; --kicad adds the
native export, DRC and baseline comparison. RF frontend nets are marked on the copy first (topor_rf.py, KiCad
Python) unless the model already lists them, so the engine keeps each on one layer. Overlapping rule areas
are merged into their union first (topor_keepouts.py), so the engine does not hull them, and custom design
rules (unrouted.kicad_dru) are read into the model (topor_rules.py). $TOPOR overrides the engine binary.

topor is the default engine of run.py and suite.py, which call route() below; this script is the same route
with its own options.
"""
import argparse
import hashlib
import json
import os
import pathlib
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent


KICAD_PYTHON = '/Applications/KiCad.app/Contents/Frameworks/Python.framework/Versions/Current/bin/python3'
KICAD_CLI = '/Applications/KiCad.app/Contents/MacOS/kicad-cli'


def route(case, out, seconds, src=None, engine_args='', no_rf=False, kicad_python=KICAD_PYTHON):
    """Copy the case to OUT/input, prepare it and route it with topor. Returns the engine's exit code.

    OUT may exist (suite.py creates it and keeps its log there) but must not hold an earlier run.
    """
    out = pathlib.Path(out)
    if any((out / n).exists() for n in ['input', 'routing.json', 'layout.json', 'routed.kicad_pcb', 'acceptance.json']):
        raise FileExistsError('use a fresh output directory; benchmark evidence must not be overwritten')
    exe = pathlib.Path(os.environ.get('TOPOR') or ROOT.parents[1] / 'target/release/topor')
    if not exe.exists():
        raise FileNotFoundError(f'{exe} not found: run cargo build --release (or set $TOPOR)')
    src = pathlib.Path(src) if src else ROOT / 'cases' / case
    if not (src / 'model.json').exists():
        raise FileNotFoundError(f'{src / "model.json"} not found')
    out.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, out / 'input')
    hashes = {q.name: hashlib.sha256(q.read_bytes()).hexdigest() for q in (out / 'input').iterdir() if q.is_file() and q.name != 'sha256.json'}
    (out / 'input' / 'sha256.json').write_text(json.dumps(hashes, indent=2))
    model = json.loads((out / 'input' / 'model.json').read_text())
    if 'custom_rules' not in model and (out / 'input' / 'unrouted.kicad_dru').exists():
        r = subprocess.run([sys.executable, str(ROOT / 'topor_rules.py'), str(out / 'input')], capture_output=True, text=True)
        print('custom rules:', r.stdout.strip().splitlines()[-1:] if r.returncode == 0 else r.stderr.strip().splitlines()[-1:], flush=True)
    if not model.get('rule_areas_merged') and model.get('rule_areas') and pathlib.Path(kicad_python).exists():
        r = subprocess.run([kicad_python, str(ROOT / 'topor_keepouts.py'), str(out / 'input')], capture_output=True, text=True)
        print('rule areas:', r.stdout.strip().splitlines()[-1:] if r.returncode == 0 else r.stderr.strip().splitlines()[-1:], flush=True)
    if not no_rf and 'rf_nets' not in model and pathlib.Path(kicad_python).exists():
        r = subprocess.run([kicad_python, str(ROOT / 'topor_rf.py'), str(out / 'input')], capture_output=True, text=True)
        if r.returncode:
            print('RF frontend detection failed:', r.stderr.strip().splitlines()[-1:], flush=True)
    cmd = [str(exe), 'route', str(out / 'input'), '--out', str(out), '--seconds', str(seconds)] + engine_args.split()
    with (out / 'engine.log').open('w') as log:
        code = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
    lines = (out / 'engine.log').read_text().splitlines()
    print(next((l for l in reversed(lines) if l.startswith('done in')), lines[-1] if lines else '') if code == 0 else f'engine failed ({code})', flush=True)
    return code


def main():
    p = argparse.ArgumentParser()
    p.add_argument('case')
    p.add_argument('--out', type=pathlib.Path, required=True)
    p.add_argument('--input', type=pathlib.Path, help='prepared input directory (default bench/pcb/cases/CASE)')
    p.add_argument('--seconds', type=float, default=300)
    p.add_argument('--engine-args', default='', help='extra arguments for topor route')
    p.add_argument('--kicad', action='store_true', help='export, run KiCad DRC and compare with the baseline')
    p.add_argument('--no-rf', action='store_true', help='do not mark RF frontend nets (they may then change layers)')
    p.add_argument('--kicad-python', default=KICAD_PYTHON)
    p.add_argument('--kicad-cli', default=KICAD_CLI)
    a = p.parse_args()
    if a.out.exists():
        p.error('use a fresh output directory')
    code = route(a.case, a.out, a.seconds, a.input, a.engine_args, a.no_rf, a.kicad_python)
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
