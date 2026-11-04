#!/usr/bin/env python3
"""Turn tscircuit's routed traces (route.mjs output) into a gridroute layout.json for kicad_io.py export.

Wires become tracks on their KiCad layer and net, vias through vias of the size tscircuit gave them (else the
SimpleRouteJson minimum). Coordinates go back to y down. Also writes routing.json with what tscircuit reported
and simple counts: connections with copper, length, vias, signal copper on plane layers.
"""
import argparse
import json
import math
import pathlib


def convert(srj, meta, out):
    back = {v: k for k, v in meta['layer_names'].items()}
    conns = meta['connections']
    plane_layers = set(meta['plane_layers'])
    tracks, vias, skipped = [], [], {}
    routed = set()
    length = 0.0
    on_plane = 0.0
    dropped_on_plane = 0.0
    for t in out['traces']:
        name = t.get('connection_name') or t['pcb_trace_id']
        info = conns.get(name) or conns.get(name.split('__')[0])
        if info is None:
            root = next((c for c in conns if name.startswith(c + '_')), None)
            info = conns.get(root)
        if info is None:
            skipped['trace of unknown connection'] = skipped.get('trace of unknown connection', 0) + 1
            continue
        net = info['net']
        routed.add(name)
        run = None
        layer_prev_xy = None

        def close():
            nonlocal run, length, on_plane
            if run and len(run['pts']) >= 2:
                tracks.append(run)
                seg = sum(math.dist(a, b) for a, b in zip(run['pts'], run['pts'][1:]))
                length += seg
                if run['layer'] in plane_layers and info['kind'] == 'signal':
                    on_plane += seg
            run = None

        for p in t['route']:
            kind = p.get('route_type')
            if kind == 'wire':
                layer, w, xy = back[p['layer']], p['width'], [p['x'], -p['y']]
                if meta['planes'].get(net) == layer:
                    # a plane net's wire on its own plane: the pour carries it (a drop ends at its via)
                    close()
                    if layer_prev_xy is not None:
                        dropped_on_plane += math.dist(layer_prev_xy, xy)
                    layer_prev_xy = xy
                    continue
                layer_prev_xy = None
                if run and (run['layer'] != layer or abs(run['width'] - w) > 1e-9):
                    last = run['pts'][-1]
                    close()
                    run = dict(layer=layer, net=net, width=w, pts=[last])
                if run is None:
                    run = dict(layer=layer, net=net, width=w, pts=[])
                if not run['pts'] or run['pts'][-1] != xy:
                    run['pts'].append(xy)
            elif kind == 'via':
                close()
                layer_prev_xy = None
                vias.append(dict(x=p['x'], y=-p['y'], net=net,
                                 d=p.get('via_diameter') or srj.get('minViaPadDiameter') or srj['minViaDiameter'],
                                 drill=p.get('via_hole_diameter') or srj['minViaHoleDiameter']))
            else:
                skipped[kind] = skipped.get(kind, 0) + 1
        close()
    layout = dict(keepouts=[], keepouts_rule=[], labels=[], outline=[], parts={}, rects=[], texts=[],
                  tracks=tracks, vias=vias, zones=[])
    by_kind = {}
    for name, info in conns.items():
        k = by_kind.setdefault(info['kind'], [0, 0])
        k[0] += 1
        k[1] += name in routed
    report = dict(algorithm='tscircuit capacity-autorouter', version=out.get('version'), pipeline=out.get('pipeline'),
                  solver=out.get('solver'), routing_seconds=out.get('route_seconds'), solved=out.get('solved'),
                  failed=out.get('failed'), timed_out=out.get('timed_out'), solver_error=out.get('solver_error'),
                  exception=(out.get('exception') or '').split('\n')[0] or None, output_error=out.get('output_error'),
                  connections_with_copper={k: dict(total=v[0], with_copper=v[1]) for k, v in by_kind.items()},
                  tracks=len(tracks), vias=len(vias), length_mm=length, signal_copper_on_plane_layers_mm=on_plane,
                  plane_net_wire_left_to_pour_mm=dropped_on_plane,
                  skipped_route_kinds=skipped, unsupported_inputs=meta.get('unsupported'))
    return layout, report


def main():
    p = argparse.ArgumentParser()
    p.add_argument('srj_dir', type=pathlib.Path, help='directory with srj.json and srj-meta.json')
    p.add_argument('routed', type=pathlib.Path, help='route.mjs output')
    p.add_argument('--out', type=pathlib.Path, required=True)
    a = p.parse_args()
    srj = json.loads((a.srj_dir / 'srj.json').read_text())
    meta = json.loads((a.srj_dir / 'srj-meta.json').read_text())
    layout, report = convert(srj, meta, json.loads(a.routed.read_text()))
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / 'layout.json').write_text(json.dumps(layout))
    (a.out / 'routing.json').write_text(json.dumps(report, indent=1))
    print(json.dumps({k: v for k, v in report.items() if k != 'unsupported_inputs'}))


if __name__ == '__main__':
    main()
