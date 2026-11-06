"""Far-field radiation pattern over the board in Blender.

    blender --background --python tools/blender_farfield.py -- PATTERN.json \\
        --board board_top.png --outline X0 Y0 X1 Y1 [--thickness 0.67] --centre X Y Z \\
        [--top DB] [--range 30] [--radius 30] [--render out.png] [--save out.blend] [--samples 32]

PATTERN.json is written by `fieldsolve cosim` with "farfield": true: theta/phi in degrees (theta from +z,
phi from +x, the solver's frame = Gerber orientation, mm, z up) and the EIRP per watt incident at the
driven port in dB on that grid.

The pattern is a surface around --centre (e.g. the chip's RF pin): its radius and colour both map
dB in [top - range, top] to [0, 1] (radius 0 below the floor). Pass the same --top and --range for
patterns that should be compared, so equal sizes and colours mean equal radiated levels. The board is a
slab with --board (a top-copper image covering --outline) on its top face.
"""
import json
import math
import pathlib
import sys

import bpy

argv = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else []
if not argv:
    raise SystemExit(__doc__)


def opt(k, default=None, n=1):
    if k not in argv:
        return default
    i = argv.index(k)
    return argv[i + 1] if n == 1 else argv[i + 1:i + 1 + n]


pat = json.loads(pathlib.Path(argv[0]).read_text())
theta = [math.radians(t) for t in pat['theta_deg']]
phi = [math.radians(p) for p in pat['phi_deg']]
g = pat['eirp_per_w_db']
nt, npz = len(theta), len(phi)
top = float(opt('--top', max(g)))
rng = float(opt('--range', 20.0))
R0 = float(opt('--radius', 25.0))
cx, cy, cz = (float(v) for v in opt('--centre', ['0', '0', '0'], 3))
x0, y0, x1, y1 = (float(v) for v in opt('--outline', ['0', '0', '10', '10'], 4))
th = float(opt('--thickness', 0.67))

# fresh scene
bpy.ops.wm.read_factory_settings(use_empty=True)
scene = bpy.context.scene

# --- board slab with the top-copper image on its top face
mesh = bpy.data.meshes.new('board')
verts = [(x0, y0, 0), (x1, y0, 0), (x1, y1, 0), (x0, y1, 0), (x0, y0, th), (x1, y0, th), (x1, y1, th), (x0, y1, th)]
faces = [(4, 5, 6, 7), (0, 3, 2, 1), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
mesh.from_pydata(verts, [], faces)
mesh.update()
uv = mesh.uv_layers.new(name='uv')
for poly in mesh.polygons:
    for li in poly.loop_indices:
        v = mesh.vertices[mesh.loops[li].vertex_index].co
        uv.data[li].uv = ((v.x - x0) / (x1 - x0), (v.y - y0) / (y1 - y0))
board = bpy.data.objects.new('board', mesh)
scene.collection.objects.link(board)
mb = bpy.data.materials.new('board')
mb.use_nodes = True
nt_ = mb.node_tree
bsdf = nt_.nodes['Principled BSDF']
bsdf.inputs['Roughness'].default_value = 0.5
img_path = opt('--board')
if img_path:
    tex = nt_.nodes.new('ShaderNodeTexImage')
    tex.image = bpy.data.images.load(str(pathlib.Path(img_path).resolve()))
    nt_.links.new(tex.outputs['Color'], bsdf.inputs['Base Color'])
else:
    bsdf.inputs['Base Color'].default_value = (0.05, 0.25, 0.1, 1)
board.data.materials.append(mb)

# --- pattern surface: one vertex per (theta, phi), radius and colour from dB
pm = bpy.data.meshes.new('pattern')
pv, pf, val = [], [], []
for it in range(nt):
    for ip in range(npz):
        x = (g[it * npz + ip] - (top - rng)) / rng
        x = min(max(x, 0.0), 1.0)
        r = R0 * x
        st, ct = math.sin(theta[it]), math.cos(theta[it])
        pv.append((cx + r * st * math.cos(phi[ip]), cy + r * st * math.sin(phi[ip]), cz + r * ct))
        val.append(x)
for it in range(nt - 1):
    for ip in range(npz):
        a, b = it * npz + ip, it * npz + (ip + 1) % npz
        pf.append((a, b, b + npz, a + npz))
pm.from_pydata(pv, [], pf)
pm.update()
attr = pm.attributes.new('level', 'FLOAT', 'POINT')
attr.data.foreach_set('value', val)
for poly in pm.polygons:
    poly.use_smooth = True
pobj = bpy.data.objects.new('pattern', pm)
scene.collection.objects.link(pobj)
mp = bpy.data.materials.new('pattern')
mp.use_nodes = True
n = mp.node_tree.nodes
l = mp.node_tree.links
pb = n['Principled BSDF']
at = n.new('ShaderNodeAttribute')
at.attribute_name = 'level'
ramp = n.new('ShaderNodeValToRGB')
cr = ramp.color_ramp
stops = [(0.0, (0.0, 0.0, 0.016)), (0.25, (0.34, 0.06, 0.43)), (0.5, (0.74, 0.22, 0.33)), (0.75, (0.98, 0.56, 0.04)), (1.0, (0.99, 1.0, 0.64))]
cr.elements[0].position, cr.elements[0].color = stops[0][0], (*stops[0][1], 1)
cr.elements[1].position, cr.elements[1].color = stops[-1][0], (*stops[-1][1], 1)
for pos, col in stops[1:-1]:
    e = cr.elements.new(pos)
    e.color = (*col, 1)
l.new(at.outputs['Fac'], ramp.inputs['Fac'])
l.new(ramp.outputs['Color'], pb.inputs['Base Color'])
l.new(ramp.outputs['Color'], pb.inputs['Emission Color'])
pb.inputs['Emission Strength'].default_value = 0.6
pb.inputs['Alpha'].default_value = 0.5
pb.inputs['Roughness'].default_value = 0.35
try:
    mp.surface_render_method = 'BLENDED'
except AttributeError:
    mp.blend_method = 'BLEND'
pobj.data.materials.append(mp)
# a fine dark wireframe over the surface so the lobe shape reads through the transparency
wire = bpy.data.objects.new('pattern_wire', pm)
scene.collection.objects.link(wire)
mod = wire.modifiers.new('wire', 'WIREFRAME')
mod.thickness = 0.06
mod.use_replace = True
mw = bpy.data.materials.new('wire')
mw.use_nodes = True
mw.node_tree.nodes['Principled BSDF'].inputs['Base Color'].default_value = (0.02, 0.02, 0.02, 1)
mw.node_tree.nodes['Principled BSDF'].inputs['Alpha'].default_value = 0.6
try:
    mw.surface_render_method = 'BLENDED'
except AttributeError:
    mw.blend_method = 'BLEND'
wire.data = pm.copy()
wire.data.materials.clear()
wire.data.materials.append(mw)

# --- light, world, camera
sun = bpy.data.objects.new('sun', bpy.data.lights.new('sun', 'SUN'))
sun.data.energy = 3.0
sun.rotation_euler = (math.radians(40), math.radians(10), math.radians(30))
scene.collection.objects.link(sun)
world = bpy.data.worlds.new('world')
world.use_nodes = True
world.node_tree.nodes['Background'].inputs['Color'].default_value = (0.02, 0.02, 0.025, 1)
world.node_tree.nodes['Background'].inputs['Strength'].default_value = 1.0
scene.world = world
cam = bpy.data.objects.new('cam', bpy.data.cameras.new('cam'))
scene.collection.objects.link(cam)
scene.camera = cam
bx, by = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
target = ((bx + cx) / 2, (by + cy) / 2, R0 * 0.25)
dist = max(x1 - x0, y1 - y0, 2 * R0) * 1.9
# view from the side of --view: 'rf' looks from the board end nearest the centre (the RF end), 'far' from the other end
side = 1.0 if (opt('--view', 'rf') == 'rf') == (cy > by) else -1.0
cam.location = (target[0] + dist * 0.55, target[1] + side * dist * 0.75, target[2] + dist * 0.62)
d = [target[i] - cam.location[i] for i in range(3)]
cam.rotation_euler = (math.atan2(math.hypot(d[0], d[1]), -d[2]), 0, math.atan2(d[1], d[0]) - math.pi / 2)
cam.data.lens = 40
cam.data.clip_end = 10000

out = opt('--render')
if out:
    for eng in ('BLENDER_EEVEE', 'BLENDER_EEVEE_NEXT', 'CYCLES'):
        try:
            scene.render.engine = eng
            break
        except TypeError:
            continue
    if scene.render.engine == 'CYCLES':
        scene.cycles.samples = int(opt('--samples', 32))
    scene.render.resolution_x, scene.render.resolution_y = 1400, 1100
    scene.render.filepath = str(pathlib.Path(out).resolve())
    bpy.ops.render.render(write_still=True)
save = opt('--save')
if save:
    bpy.ops.file.pack_all()  # keep the board texture inside the .blend
    bpy.ops.wm.save_as_mainfile(filepath=str(pathlib.Path(save).resolve()))
