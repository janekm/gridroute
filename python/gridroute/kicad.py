"""KiCad inputs for gridroute.board: the XML netlist, and footprint geometry dumped with KiCad's own Python.

    from gridroute.kicad import load_netlist, load_footprints
    comps, pin_net = load_netlist('project.xml')            # kicad-cli sch export netlist --format kicadxml
    footprints = load_footprints('footprints.json')

The footprint dump needs pcbnew, so it runs under KiCad's Python (the layout itself does not):

    /Applications/KiCad.app/Contents/Frameworks/Python.framework/Versions/Current/bin/python3 \\
        -m gridroute.kicad footprints project.xml footprints.json [--lib NAME=DIR ...]

Stock libraries are found in KiCad's footprint directory; project libraries are given with --lib.
"""
import json
import os
import re
import sys
import xml.etree.ElementTree as ET

STOCK = '/Applications/KiCad.app/Contents/SharedSupport/footprints'


def load_netlist(path):
    """KiCad XML netlist -> (comps, pin_net): comps {ref: {value, footprint, sheet}}, pin_net {(ref, pin): net}.
    KiCad's generated names 'Net-(...)' / 'unconnected-(...)' have '/' replaced by '{slash}' (as in pcbnew)."""
    root = ET.parse(path).getroot()
    comps = {}
    for c in root.find('components'):
        sp = c.find('sheetpath')
        comps[c.get('ref')] = dict(value=c.findtext('value'), footprint=c.findtext('footprint'),
                                   sheet=sp.get('names').strip('/') if sp is not None else '')
    pin_net = {}
    for n in root.find('nets'):
        name = n.get('name')
        m = re.match(r'^(unconnected|Net)-\((.*)\)$', name)
        if m:
            name = '%s-(%s)' % (m.group(1), m.group(2).replace('/', '{slash}'))
        for nd in n.findall('node'):
            pin_net[(nd.get('ref'), nd.get('pin'))] = name
    return comps, pin_net


def load_footprints(path):
    with open(path) as f:
        return json.load(f)


def dump_footprints(netlist, out, libs=None):
    """Write the pad and courtyard geometry of every footprint in the netlist to `out` (KiCad python only):
    {lib:name: {pads: [...], courtyard: [x0, y0, x1, y1]}}, mm relative to the footprint origin at 0 degrees, front
    side. libs: {library nickname: .pretty directory} for project libraries."""
    import pcbnew
    shapes = {pcbnew.PAD_SHAPE_RECT: 'rect', pcbnew.PAD_SHAPE_ROUNDRECT: 'roundrect', pcbnew.PAD_SHAPE_CIRCLE: 'circle',
              pcbnew.PAD_SHAPE_OVAL: 'oval'}
    mm = lambda v: round(pcbnew.ToMM(v), 4)
    root = ET.parse(netlist).getroot()
    names = sorted({c.findtext('footprint') for c in root.find('components') if c.findtext('footprint')})
    res = {}
    for full in names:
        lib, name = full.split(':')
        fp = pcbnew.FootprintLoad((libs or {}).get(lib) or os.path.join(STOCK, lib + '.pretty'), name)
        pads = []
        for p in fp.Pads():
            ls = p.GetLayerSet()
            tht = p.GetAttribute() in (pcbnew.PAD_ATTRIB_PTH, pcbnew.PAD_ATTRIB_NPTH)
            if not tht and not (ls.Contains(pcbnew.F_Cu) or ls.Contains(pcbnew.B_Cu)):
                continue                  # paste-only apertures (e.g. split exposed-pad stencil openings)
            sz = p.GetSize(pcbnew.F_Cu)
            drill = p.GetDrillSize()
            shape = p.GetShape(pcbnew.F_Cu)
            pads.append(dict(num=p.GetNumber(), x=mm(p.GetPosition().x), y=mm(p.GetPosition().y), w=mm(sz.x),
                             h=mm(sz.y), rot=round(p.GetOrientation().AsDegrees(), 2), shape=shapes.get(shape, 'rect'),
                             layers='all' if tht else ('F' if ls.Contains(pcbnew.F_Cu) else 'B'),
                             npth=p.GetAttribute() == pcbnew.PAD_ATTRIB_NPTH,
                             drill=mm(max(drill.x, drill.y)) if tht else 0.0,
                             rr=round(mm(p.GetRoundRectCornerRadius(pcbnew.F_Cu)), 4)
                             if shape == pcbnew.PAD_SHAPE_ROUNDRECT else 0.0))
        cy = fp.GetCourtyard(pcbnew.F_CrtYd)
        bb = cy.BBox() if cy.OutlineCount() else fp.GetBoundingBox(False)
        res[full] = dict(pads=pads, courtyard=[mm(bb.GetLeft()), mm(bb.GetTop()), mm(bb.GetRight()), mm(bb.GetBottom())])
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, 'w') as f:
        json.dump(res, f, indent=1)
    return len(res)


if __name__ == '__main__':
    a = sys.argv[1:]
    if len(a) < 3 or a[0] != 'footprints':
        sys.exit(__doc__)
    libs = dict(x.split('=', 1) for x in a[3:] if x != '--lib')
    print('wrote', a[2], dump_footprints(a[1], a[2], libs), 'footprints')
