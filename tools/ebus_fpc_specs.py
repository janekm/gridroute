"""fieldsolve specs for the eBUS_IF inline adapter's FPC antenna (KiCad coordinates, z up from the PCB underside).

The FPC is a stand-in for the Molex 146153 class part named in the enclosure README (34.9 x 9.0 x 0.1 mm), whose
copper pattern is not published: an off-centre-fed dipole, a 13 x 8 mm shield-side patch at the left and a
1.5 mm radiator arm folded round the right end, its return length `ret` the tuning knob. It sits in the cover's
pocket (centred x 20.0, y 16.4), copper at z 9.00 on a 0.1 mm polyimide carrier (8.95..9.05) under the ceiling.

    python3 ebus_fpc_specs.py free NAME RET [ff] [dry]            FPC on a 1.05 mm PC/ABS plate in free space
    python3 ebus_fpc_specs.py inst NAME RET [ff] [nocable] [dry]  in the enclosure over the board
    ... [yl1439] [air]: the FORCE YL1439-BT-V1 instead of the stand-in (RET ignored); "air" drops the lid plate

Run from a directory holding `gbr/` (the board's Gerbers, exports/fab/gerbers at p-hw 81d5260c) and, for
"free", `dummy_gbr/` (a 0.2 mm board of air with one 0.02 mm pad, so that only the FPC and plate are
modelled). ENC points at the enclosure's build/ directory (asm-base.stl, asm-cover.stl). RET 8.5 tunes the
stand-in to about 2.44 GHz on the plate. See docs/ebus-fpc-antenna-findings.md.
"""
import json, os, sys

ENC = os.path.expanduser(os.environ.get('ENC', '~/src/p-hw/.claude/worktrees/ebus-if-inline/Revision_7/eBUS_IF_inline/enclosure-moulded/build'))
X0, X1, Y0, Y1 = 20.0 - 34.9 / 2, 20.0 + 34.9 / 2, 16.4 - 4.5, 16.4 + 4.5  # FPC outline
ZC, ZS0, ZS1 = 9.0, 8.95, 9.05  # copper, substrate
PCABS = {'eps': 2.9, 'tand': 0.008}  # FR PC/ABS at 2.4 GHz (estimate; datasheets give 1 MHz values)
FEED_Y = Y0 + 0.75  # feed on the edge nearer J3 (KiCad y 4.48)
GX = X0 + 13.0  # ground patch right edge
DECOUPLING = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19',
              'C20', 'C21', 'C29', 'C30']


def fpc(ret):
    m, w, gap = 0.5, 1.5, 0.5
    ax0 = GX + gap
    rects = [
        {'rect': [X0 + m, Y0 + m, GX, Y1 - m], 'net': 'FPC_GND'},
        {'rect': [ax0, Y0 + m - 0.0, X1 - m, Y0 + m + w], 'net': 'ANT'},           # arm along the lower edge
        {'rect': [X1 - m - w, Y0 + m, X1 - m, Y1 - m], 'net': 'ANT'},              # fold up the right end
    ]
    if ret > 0:
        rects.append({'rect': [X1 - m - ret, Y1 - m - w, X1 - m, Y1 - m], 'net': 'ANT'})  # return along the top
    return rects


# FORCE (ShenZhen Powerful Technology) YL1439-BT-V1: 36 x 10 mm printed inverted-F, copper traced from the vendor's
# AutoCAD drawing (Downloads/YL1439-BT天线-FPC.pdf, 2022). Frame: mm from the FPC's top-left corner, y down. A
# 1.25 mm arm (y 4.80-6.05) runs from its open end at the left (with a tab down) to the feed pad; a 0.5 mm bridge
# over a diagonal slot shorts it to the cable block at the right, which carries a 1.4 mm ground strip back to
# x 14.84. Centre conductor on the feed pad (28.0-29.6), braid on the pad at x 31.0-32.6.
YL1439 = [
    ('rect', [0.22, 4.80, 28.01, 6.05]),                    # arm
    ('rect', [0.22, 6.05, 1.48, 7.79]),                     # tab at the open end
    ('rect', [28.01, 4.80, 29.58, 6.92]),                   # feed pad
    ('poly', [(29.33, 4.80), (29.89, 4.39), (31.0, 4.39), (31.0, 4.91), (30.87, 4.91), (29.58, 5.87), (29.58, 4.80)]),  # short bridge
    ('poly', [(31.0, 4.39), (35.83, 4.39), (35.83, 7.93), (14.84, 7.93), (14.84, 6.53), (27.53, 6.53), (27.53, 7.27),
              (30.08, 7.27), (30.08, 6.68), (31.0, 5.61)]),  # cable block and ground strip
]
YL_W, YL_H = 36.0, 10.0
YL_FEED = ((29.3, 6.4), (30.75, 6.4))   # feed pad -> block across the slot
YL_BRAID = (31.8, 6.4)


def yl1439(ox, oy):
    """The YL1439 copper with its top-left corner at KiCad (ox, oy)."""
    out = []
    for kind, g in YL1439:
        if kind == 'rect':
            out.append(('rect', {'rect': [ox + g[0], oy + g[1], ox + g[2], oy + g[3]], 'net': 'FPC'}))
        else:
            out.append(('poly', {'poly': [[ox + x, oy + y] for x, y in g], 'net': 'FPC'}))
    return out


def port():
    return {'name': 'FPC', 'layer': 'FPC', 'pos': [GX + 0.75, FEED_Y], 'neg': [GX - 0.25, FEED_Y], 'width': 1.0, 'r': 50}


def common(name, ret, ff):
    s = {
        'coords': 'kicad', 'cell': 0.1, 'air': 15.0, 'air_z': 15.0, 'fmin': 1e8, 'fmax': 8e9, 'nf': 400,
        'end_db': -60, 'max_steps': 200000, 'max_cell': 1.5, 'f_loss': 2.44e9, 'body_cell': 0.4,
        'out': name, 'ports': [port()], 'excite': 'all',
        'layers': [{'name': 'FPC', 'z': ZC, 'rects': fpc(ret)}],
        'bodies': [{'name': 'fpc_carrier', 'box': [X0, Y0, ZS0, X1, Y1, ZS1], 'eps': 3.4, 'tand': 0.005}],
        'roi': [X0 - 1.0, Y0 - 1.0, X1 + 1.0, Y1 + 1.0],
    }
    if ff:
        s['farfield'] = {'freqs': [round(f * 1e9) for f in [2.3, 2.35, 2.4, 2.42, 2.44, 2.46, 2.48, 2.5, 2.55, 2.6]], 'step_deg': 5, 'gap': 3, 'tol': 0.01}
    return s


def free(name, ret, ff):
    s = common(name, ret, ff)
    s['gerbers'] = 'dummy_gbr'
    s['dielectrics'] = [{'t': 0.8, 'er': 1.0, 'tand': 0.0}]
    s['preset'] = 'antenna'
    s['bodies'].append({'name': 'lid_plate', 'box': [X0 - 2.0, Y0 - 2.5, ZS1, X1 + 2.0, Y1 + 2.5, ZS1 + 1.05], **PCABS})
    return s


def inst(name, ret, ff, cable=True):
    s = common(name, ret, ff)
    s['gerbers'] = 'gbr'
    s['preset'] = 'antenna'
    s['bodies'] = [
        {'name': 'base', 'stl': os.path.join(ENC, 'asm-base.stl'), 'frame': 'y_up', **PCABS},
        {'name': 'cover', 'stl': os.path.join(ENC, 'asm-cover.stl'), 'frame': 'y_up', **PCABS},
    ] + s['bodies'] + [
        {'name': f'screw {h}', 'cylinder': [x, y, 1.1, -2.8, 8.3], 'metal': True} for h, x, y in [('H1', 19.7, 4.0), ('H2', 15.8, 24.0), ('H3', 45.0, 4.3)]
    ]
    # decoupling at RF: every capacitor of 47 pF or more as a short (the R7 A.4 BOM fits all of them); the
    # 12 pF crystal loads and the 1-1.2 pF RF match stay open. Unmodelled chips' signal pins get 3 pF.
    s['elements'] = [{'kind': 'short', 'pos': f'{c}.1', 'neg': f'{c}.2', 'layer': 'F.Cu'} for c in DECOUPLING]
    s['pin_loads'] = [{'ref': u, 'c': 3e-12} for u in ('U4', 'U2', 'U5', 'U6', 'U1')]
    s['coarse_copper_ok'] = [[-2.0, -2.0, 52.0, 30.0]]
    # the bus lines run to J1 and on into the eBUS cable (not modelled); TX_B and CC_S are transistor nodes
    s['open_nets_ok'] = ['/eBUS interface/EBUS_B', '/eBUS interface/TX_B', '/eBUS interface/EBUS_A', '/Power/CC_S']  # board copper away from the FPC carries no antenna current of its own
    if cable:
        # the cable's shield (1.13 mm coax, outer surface only): down from the FPC's shield pad to z 3.0, along to
        # the mated u.FL plug on J3, and the plug's shell down to J3's ground pads
        gx, gy = GX - 0.6, FEED_Y
        jx, jy = 24.227, 3.72
        zc = 3.0
        s['layers'].append({'name': 'CABLE', 'z': zc, 'net': 'SHIELD', 'polygons': [strip(gx, gy, jx, jy + 0.6, 1.13)],
                            'rects': [[jx - 1.6, jy - 0.6, jx + 1.6, jy + 1.6]]})
        s['bodies'] += [
            {'name': 'shield riser', 'cylinder': [gx, gy, 0.56, zc, ZC], 'metal': True},
            {'name': 'plug shell L', 'cylinder': [22.757 + 0.3, 3.72, 0.5, 0.69, zc], 'metal': True},
            {'name': 'plug shell R', 'cylinder': [25.697 - 0.3, 3.72, 0.5, 0.69, zc], 'metal': True},
        ]
    return s


def strip(xa, ya, xb, yb, w):
    import math
    dx, dy = xb - xa, yb - ya
    L = math.hypot(dx, dy)
    nx, ny = -dy / L * w / 2, dx / L * w / 2
    return [[xa + nx, ya + ny], [xb + nx, yb + ny], [xb - nx, yb - ny], [xa - nx, ya - ny]]


def use_yl1439(s, installed, cable):
    """Swap the stand-in for the YL1439-BT-V1, centred on the pocket (x 2.0-38.0, y 11.4-21.4: 0.5 x 0.4 mm larger
    than the pocket drawn for the Molex part), arm towards the board's left. Its 1.13 mm cable leaves along the
    FPC's axis past the right end; installed, the shield is modelled lying under the FPC (z 8.6) to x 38.5, down to
    z 3.0 there and along to J3."""
    ox, oy = 20.0 - YL_W / 2, 16.4 - YL_H / 2
    shapes = yl1439(ox, oy)
    s['layers'][0] = {'name': 'FPC', 'z': ZC, 'rects': [g for k, g in shapes if k == 'rect'], 'polygons': [g for k, g in shapes if k == 'poly']}
    (fa, fb) = YL_FEED
    s['ports'] = [{'name': 'FPC', 'layer': 'FPC', 'pos': [ox + fa[0], oy + fa[1]], 'neg': [ox + fb[0], oy + fb[1]], 'width': 0.5, 'r': 50}]
    for b in s['bodies']:
        if b['name'] == 'fpc_carrier':
            b['box'] = [ox, oy, ZS0, ox + YL_W, oy + YL_H, ZS1]
        if b['name'] == 'lid_plate':
            b['box'] = [ox - 2.0, oy - 2.5, ZS1, ox + YL_W + 2.0, oy + YL_H + 2.5, ZS1 + 1.05]
    s['roi'] = [ox - 1.0, oy - 1.0, ox + YL_W + 1.0, oy + YL_H + 1.0]
    # drop the stand-in's cable model
    s['layers'] = [l for l in s['layers'] if l['name'] not in ('CABLE', 'CABLE_TOP')]
    s['bodies'] = [b for b in s['bodies'] if b['name'] not in ('shield riser', 'plug shell L', 'plug shell R')]
    if installed and cable:
        bx, by = ox + YL_BRAID[0], oy + YL_BRAID[1]
        xe = ox + YL_W + 0.5
        jx, jy = 24.227, 3.72
        zt, zc = 8.6, 3.0
        s['layers'] += [
            {'name': 'CABLE_TOP', 'z': zt, 'net': 'SHIELD', 'rects': [[bx - 0.6, by - 0.565, xe + 0.565, by + 0.565]]},
            {'name': 'CABLE', 'z': zc, 'net': 'SHIELD', 'polygons': [strip(xe, by, jx, jy + 0.6, 1.13)], 'rects': [[jx - 1.6, jy - 0.6, jx + 1.6, jy + 1.6]]},
        ]
        s['bodies'] += [
            {'name': 'braid solder', 'cylinder': [bx, by, 0.5, zt, ZC], 'metal': True},
            {'name': 'cable bend', 'cylinder': [xe, by, 0.56, zc, zt], 'metal': True},
            {'name': 'plug shell L', 'cylinder': [22.757 + 0.3, 3.72, 0.5, 0.5, zc], 'metal': True},
            {'name': 'plug shell R', 'cylinder': [25.697 - 0.3, 3.72, 0.5, 0.5, zc], 'metal': True},
        ]
    return s


if __name__ == '__main__':
    kind, name, ret = sys.argv[1], sys.argv[2], float(sys.argv[3])
    ff = 'ff' in sys.argv[4:]
    cable = 'nocable' not in sys.argv[4:]
    s = free(name, ret, ff) if kind == 'free' else inst(name, ret, ff, cable=cable)
    if 'yl1439' in sys.argv[4:]:
        s = use_yl1439(s, kind == 'inst', cable)
        if kind == 'free' and 'air' in sys.argv[4:]:
            s['bodies'] = [b for b in s['bodies'] if b['name'] != 'lid_plate']
    if 'dry' in sys.argv[4:]:
        s['max_steps'] = 20
        s['scene'] = True
    json.dump(s, open(name + '.json', 'w'), indent=1)
    print('wrote', name + '.json')
