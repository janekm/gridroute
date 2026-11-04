#!/usr/bin/env python3
"""topor: mark a board's RF frontend nets in a prepared model, so the engine keeps each on one layer (no
vias). Run with KiCad Python on a case input dir.

An RF frontend starts at an RF pin: a pin whose function names it (ANT, RF, RF_IN, LNA_IN, PA_OUT, ...), or a
signal pin of an antenna or coaxial connector footprint (U.FL/IPEX/MHF, SMA/BNC/MMCX connectors), or a net
labelled that way. From there it runs on through the matching network: across each two-pad R, L, C or ferrite
to the net on its other pad, and stops at ground, planes and power, and at anything that is not a two-pad part
(the radio, the connector). Each RF net keeps to the copper layer its SMD pads are on.

These are usually routed by hand; this is so the engine does not hop them through vias when they are not.
"""
import json, pathlib, re, sys
import pcbnew

case = pathlib.Path(sys.argv[1])
m = json.loads((case / 'model.json').read_text())
b = pcbnew.LoadBoard(str(case / 'unrouted.kicad_pcb'))
planes = set(m.get('cold_routing', {}).get('planes', {}))

RF_TOKENS = {'ANT', 'ANTENNA', 'RF', 'RFIO', 'RFIN', 'RFOUT', 'RFP', 'RFN', 'LNA', 'LNAIN', 'PAOUT', 'TXRX', 'RFX'}
# coaxial and antenna footprints (not "D_SMA": that SMA is a diode package)
RF_FP = re.compile(r'antenna|u\.?fl|ipex|mhf\d?|coax|rf[_ ]?conn|(^|[^a-z])(sma|bnc|mmcx)[_-]?(jack|plug|edge|conn|female|male|vertical|horizontal)', re.I)
POWER = re.compile(r'(^|[^A-Z])(A|D|P|S)?GND|^(\+|-)?\d+V\d*|(^|[^A-Z])(VCC|VDD|VSS|VBAT|VBUS|VSYS|VIN|VREG|VUSB)', re.I)
PASSIVE = re.compile(r'^(R|L|C|FB|FL)\d')


def tokens(s):
    return [t for t in re.split(r'[^A-Z0-9]+', s.upper()) if t]


def rf_name(s, pin=''):
    """A pin function or net label that names an RF signal (a trailing pin number ignored)."""
    t = tokens(s)
    if t and pin and t[-1] == pin.upper():
        t = t[:-1]
    return any(x in RF_TOKENS or (x.startswith('ANT') and len(x) <= 6) or re.fullmatch(r'RF\d*[A-Z]?', x) for x in t)


npads = {}
for p in b.GetPads():
    if p.GetNetCode():
        npads.setdefault(p.GetNetname(), []).append(p)


def signal_net(net):
    """A net an RF frontend may run on: a routable signal, not ground, a plane or power."""
    if not net or net.startswith('unconnected-') or net in planes or len(npads.get(net, [])) < 2:
        return False
    return not POWER.search(net.split('/')[-1]) and len(npads[net]) <= 6


rf = {}  # net -> why
for f in b.GetFootprints():
    fp = '%s:%s' % (f.GetFPID().GetLibNickname(), f.GetFPID().GetLibItemName())
    for p in f.Pads():
        net = p.GetNetname()
        if not signal_net(net) or net in rf:
            continue
        if RF_FP.search(fp):
            rf[net] = '%s.%s of %s' % (f.GetReference(), p.GetNumber(), fp.split(':')[-1])
        elif rf_name(p.GetPinFunction() or '', p.GetNumber()):
            rf[net] = '%s.%s %s' % (f.GetReference(), p.GetNumber(), p.GetPinFunction())
for net in npads:
    if net not in rf and signal_net(net) and '/' in net and rf_name(net.split('/')[-1]):
        rf[net] = 'net label %s' % net.split('/')[-1]

# on through the matching network: across two-pad passives, and through small RF parts (diplexers, baluns,
# filters: at most four pads, two or more of them RF already) to their other signal pads
todo = list(rf)
while todo:
    net = todo.pop()
    for p in npads[net]:
        f = p.GetParentFootprint()
        pads = [q for q in f.Pads() if q.GetNetCode() or q.IsOnCopperLayer()]
        if len(pads) == 2 and PASSIVE.match(f.GetReference()):
            nxt = [pads[1] if pads[0].GetNumber() == p.GetNumber() else pads[0]]
        elif len(pads) <= 4 and sum(q.GetNetname() in rf for q in pads) >= 2:
            nxt = pads
        else:
            continue
        for q in nxt:
            on = q.GetNetname()
            if signal_net(on) and on not in rf:
                rf[on] = 'through %s from %s' % (f.GetReference(), net)
                todo.append(on)

out, skipped = [], []
for net, why in sorted(rf.items()):
    layers = None
    for p in npads[net]:
        ls = {b.GetStandardLayerName(l) for l in p.GetLayerSet().CuStack()}
        if ls - {'F.Cu', 'B.Cu'} or len(ls) > 1:
            continue  # through-hole: on every layer
        layers = ls if layers is None else layers & ls
    if not layers:
        skipped.append((net, 'pads on no common outer layer'))
        continue
    layer = sorted(layers)[0]
    if layer not in m['layers']:
        skipped.append((net, 'layer %s not routed' % layer))
        continue
    out.append(dict(net=net, layer=layer, why=why, pads=['%s.%s' % (p.GetParentFootprint().GetReference(), p.GetNumber()) for p in npads[net]]))
m['rf_nets'] = out
(case / 'model.json').write_text(json.dumps(m, indent=2))
print(json.dumps(dict(rf_nets=[(x['net'], x['layer'], x['why']) for x in out], skipped=skipped)))
