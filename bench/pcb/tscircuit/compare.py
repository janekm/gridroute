#!/usr/bin/env python3
"""Route one prepared case with topor and with tscircuit's autorouter, and judge both the same way.

    compare.py CASE --input DIR --out OUT [--engines topor,tscircuit] [--seconds S] [--tscircuit-seconds S]

OUT/CASE/topor: bench/pcb/topor_engine.py (engine binary from $TOPOR; KiCad export, DRC and acceptance when the
input has an unrouted.kicad_pcb). OUT/CASE/tscircuit: to_srj.py, route.mjs (package from $TSCIRCUIT_AUTOROUTER),
from_srj.py, then the same KiCad export, DRC and acceptance. A case made from one of tscircuit's own boards
(srj_to_model.py: srj.json and srj-meta.json in the input) gives tscircuit its original SRJ instead of to_srj.py. Both get check_layout.py (check.json), which also
covers a case without a KiCad board. Topor scripts are taken from $TOPOR_BENCH (a bench/pcb directory),
default this one.
"""
import argparse
import json
import os
import pathlib
import shutil
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
KICAD_PY = '/Applications/KiCad.app/Contents/Frameworks/Python.framework/Versions/Current/bin/python3'
KICAD_CLI = '/Applications/KiCad.app/Contents/MacOS/kicad-cli'


def run(cmd, log, **kw):
    with open(log, 'a') as f:
        f.write('$ ' + ' '.join(map(str, cmd)) + '\n')
        f.flush()
        return subprocess.run(list(map(str, cmd)), stdout=f, stderr=subprocess.STDOUT, **kw).returncode


def kicad_gate(case, out, bench, log):
    if not (out / 'input' / 'unrouted.kicad_pcb').exists():
        return
    run([KICAD_PY, bench / 'kicad_io.py', 'export', case, '--out', out], log)
    run([KICAD_CLI, 'pcb', 'drc', '--format', 'json', '--severity-all', '--all-track-errors',
         '-o', out / 'routed-drc.json', out / 'routed.kicad_pcb'], log)
    run([sys.executable, bench / 'accept.py', case, '--out', out], log)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('case')
    p.add_argument('--input', type=pathlib.Path, required=True)
    p.add_argument('--out', type=pathlib.Path, required=True)
    p.add_argument('--engines', default='topor,tscircuit')
    p.add_argument('--seconds', type=float, default=300, help='topor time limit')
    p.add_argument('--tscircuit-seconds', type=float, default=1800)
    p.add_argument('--pipeline', default='AutoroutingPipelineSolver')
    p.add_argument('--lift-pair-cap', action='store_true', help="raise tscircuit's 1000-connection pairing cap")
    a = p.parse_args()
    bench = pathlib.Path(os.environ.get('TOPOR_BENCH', HERE.parent))
    engines = a.engines.split(',')
    if 'topor' in engines:
        out = a.out / a.case / 'topor'
        out.parent.mkdir(parents=True, exist_ok=True)
        log = a.out / a.case / 'topor.log'
        run([sys.executable, bench / 'topor_engine.py', a.case, '--input', a.input, '--out', out,
             '--seconds', a.seconds] + (['--kicad'] if (a.input / 'unrouted.kicad_pcb').exists() else []), log)
        run([sys.executable, HERE / 'check_layout.py', out / 'input' / 'model.json', out / 'layout.json',
             '--out', out / 'check.json'], log)
    if 'tscircuit' in engines:
        out = a.out / a.case / 'tscircuit'
        out.mkdir(parents=True)
        log = a.out / a.case / 'tscircuit.log'
        shutil.copytree(a.input, out / 'input')
        if (a.input / 'srj.json').exists() and (a.input / 'srj-meta.json').exists():
            # one of tscircuit's own boards (srj_to_model.py): tscircuit routes its original SRJ
            shutil.copyfile(a.input / 'srj.json', out / 'srj.json')
            shutil.copyfile(a.input / 'srj-meta.json', out / 'srj-meta.json')
        else:
            run([sys.executable, HERE / 'to_srj.py', out / 'input', '--out', out], log)
        run(['node', '--max-old-space-size=16000', HERE / 'route.mjs', out / 'srj.json', out / 'tscircuit-out.json',
             '--seconds', a.tscircuit_seconds, '--pipeline', a.pipeline] + (['--lift-pair-cap'] if a.lift_pair_cap else []), log)
        run([sys.executable, HERE / 'from_srj.py', out, out / 'tscircuit-out.json', '--out', out], log)
        kicad_gate(a.case, out, bench, log)
        run([sys.executable, HERE / 'check_layout.py', out / 'input' / 'model.json', out / 'layout.json',
             '--out', out / 'check.json'], log)
    for e in engines:
        f = a.out / a.case / e / 'check.json'
        if f.exists():
            r = json.loads(f.read_text())
            acc = a.out / a.case / e / 'acceptance.json'
            kicad = json.loads(acc.read_text()) if acc.exists() else {}
            print(e, json.dumps(dict(complete=f"{r['complete_nets']}/{r['nets']}", open=r['open_pad_groups'],
                                     clearance=r['clearance_total'], vias=r['vias'], length_mm=r['track_length_mm'],
                                     kicad_passed=kicad.get('passed'), kicad_unconnected=kicad.get('native_unconnected_items'))))


if __name__ == '__main__':
    main()
