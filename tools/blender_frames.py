"""Animate fieldsolve time-domain surface-current frames in Blender.

    /Applications/Blender.app/Contents/MacOS/Blender -b --python tools/blender_frames.py -- RUN_DIR \
        [--explode 3.0] [--gamma 0.5] [--every 1] [--render out.mp4] [--save out.blend] [--png-only]

RUN_DIR is a `fieldsolve` output directory with frames/ (from a spec's "frames": {"layers": [...],
"every_ps": 2}) and scene.json. For each layer, every frame's |J| (sheet current, A/m) becomes a colour-mapped
PNG (fixed scale per layer over the whole run, `gamma` compresses the range). They are written to
frames/png/<layer>_####.png, which ffmpeg can turn into a 2-D movie directly. The script then builds one
plane per layer with that image sequence:
- planes are separated vertically by `explode` mm so that both stay visible;
- copper outlines from scene.json are drawn as thin white lines (clear areas such as plane slots too);
- the camera is in perspective.

--render writes an MP4 (Eevee, PNG frames then ffmpeg). Scene units are mm, Gerber orientation.
"""
import json
import math
import pathlib
import shutil
import subprocess
import sys

import bpy
import numpy as np

argv = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else []
run = pathlib.Path(argv[0]).resolve()


def opt(k, d=None):
    return argv[argv.index(k) + 1] if k in argv else d


explode = float(opt('--explode', 3.0))
gamma = float(opt('--gamma', 0.5))
every = int(opt('--every', 1))
meta = json.loads((run / 'frames' / 'frames.json').read_text())
scene_j = json.loads((run / 'scene.json').read_text())
nx, ny = meta['shape']
X, Y = np.array(meta['x']), np.array(meta['y'])
pml = meta.get('pml', [0] * 6)
i0, i1 = pml[0], nx - 1 - pml[1]
j0, j1 = pml[2], ny - 1 - pml[3]
frames = meta['frames'][::every]
fidx = list(range(0, len(meta['frames']), every))

# uniform resampling of the interior (nearest node)
h = min(np.diff(X[i0:i1 + 1]).min(), np.diff(Y[j0:j1 + 1]).min())
ux = np.arange(X[i0], X[i1] + 1e-9, h)
uy = np.arange(Y[j0], Y[j1] + 1e-9, h)
ix = np.clip(np.searchsorted(X, ux), 0, nx - 1)
iy = np.clip(np.searchsorted(Y, uy), 0, ny - 1)


def cmap(v, alpha=False):
    """v in 0..1 -> RGBA, dark blue -> purple -> orange -> pale yellow; alpha follows v for see-through
    upper layers."""
    stops = np.array([[0.02, 0.02, 0.08], [0.25, 0.05, 0.5], [0.85, 0.25, 0.25], [1.0, 0.7, 0.15], [1.0, 1.0, 0.85]])
    pos = np.linspace(0, 1, len(stops))
    out = np.empty(v.shape + (4,), dtype=np.float32)
    for c in range(3):
        out[..., c] = np.interp(v, pos, stops[:, c])
    out[..., 3] = np.clip(1.6 * v, 0, 1) if alpha else 1.0
    return out


png_dir = run / 'frames' / 'png'
png_dir.mkdir(exist_ok=True)
seqs = {}
bottom = min(meta['layers'], key=lambda L: L['z'])['tag']
for L in meta['layers']:
    tag, jmax = L['tag'], max(L['j_max'], 1e-30)
    first = None
    for n, fi in enumerate(fidx):
        raw = np.fromfile(run / 'frames' / ('J_%s_%04d.f32' % (tag, fi)), dtype='<f4').reshape(2, ny, nx)
        mag = np.sqrt(raw[0] ** 2 + raw[1] ** 2)[np.ix_(iy, ix)]
        v = np.clip(mag / jmax, 0, 1) ** gamma
        img = bpy.data.images.new('tmp', width=len(ux), height=len(uy), alpha=True, float_buffer=False)
        img.pixels.foreach_set(cmap(v, alpha=tag != bottom).ravel())
        path = png_dir / ('%s_%04d.png' % (tag, n + 1))
        img.filepath_raw = str(path)
        img.file_format = 'PNG'
        img.save()
        bpy.data.images.remove(img)
        first = first or path
    seqs[tag] = (first, len(fidx), L)
print('blender_frames: wrote %d frames per layer to %s' % (len(fidx), png_dir))
if '--png-only' in argv:
    raise SystemExit(0)

for o in list(bpy.data.objects):
    bpy.data.objects.remove(o, do_unlink=True)
layers_sorted = sorted(seqs.values(), key=lambda s: s[2]['z'])
zpos = {s[2]['tag']: k * explode for k, s in enumerate(layers_sorted)}
for tag, (first, nfr, L) in seqs.items():
    img = bpy.data.images.load(str(first))
    img.source = 'SEQUENCE'
    m = bpy.data.materials.new('J ' + L['name'])
    m.use_nodes = True
    nt = m.node_tree
    for nd in list(nt.nodes):
        nt.nodes.remove(nd)
    out = nt.nodes.new('ShaderNodeOutputMaterial')
    tex = nt.nodes.new('ShaderNodeTexImage')
    tex.image = img
    tex.image_user.frame_duration = nfr
    tex.image_user.frame_start = 1
    tex.image_user.use_auto_refresh = True
    img.alpha_mode = 'STRAIGHT'
    em = nt.nodes.new('ShaderNodeEmission')
    em.inputs['Strength'].default_value = 1.0
    tr = nt.nodes.new('ShaderNodeBsdfTransparent')
    mix = nt.nodes.new('ShaderNodeMixShader')
    nt.links.new(tex.outputs['Color'], em.inputs['Color'])
    nt.links.new(tex.outputs['Alpha'], mix.inputs['Fac'])
    nt.links.new(tr.outputs['BSDF'], mix.inputs[1])
    nt.links.new(em.outputs['Emission'], mix.inputs[2])
    nt.links.new(mix.outputs['Shader'], out.inputs['Surface'])
    for attr, val in [('surface_render_method', 'BLENDED'), ('blend_method', 'BLEND')]:
        try:
            setattr(m, attr, val)
        except Exception:
            pass
    bpy.ops.mesh.primitive_plane_add(size=1.0, location=((ux[0] + ux[-1]) / 2, (uy[0] + uy[-1]) / 2, zpos[tag]))
    pl = bpy.context.active_object
    pl.name = 'J ' + L['name']
    pl.scale = (ux[-1] - ux[0], uy[-1] - uy[0], 1)
    pl.data.materials.append(m)
    # copper outlines of this layer
    sh = next((s for s in scene_j['sheets'] if s['name'] == L['name']), None)
    if sh:
        verts, edges = [], []
        for p in sh['prims']:
            base = len(verts)
            pts = p['pts']
            verts += [(x, y, zpos[tag] + 0.01) for x, y in pts]
            edges += [(base + k, base + (k + 1) % len(pts)) for k in range(len(pts))]
        me = bpy.data.meshes.new('outline ' + L['name'])
        me.from_pydata(verts, edges, [])
        ob = bpy.data.objects.new('outline ' + L['name'], me)
        bpy.context.scene.collection.objects.link(ob)
        mod = ob.modifiers.new('wire', 'NODES')
        ng = bpy.data.node_groups.new('outline_tube', 'GeometryNodeTree')
        ng.interface.new_socket('Geometry', in_out='INPUT', socket_type='NodeSocketGeometry')
        ng.interface.new_socket('Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')
        gi, go = ng.nodes.new('NodeGroupInput'), ng.nodes.new('NodeGroupOutput')
        m2c = ng.nodes.new('GeometryNodeMeshToCurve')
        prof = ng.nodes.new('GeometryNodeCurvePrimitiveCircle')
        prof.inputs['Radius'].default_value = 0.04
        prof.inputs['Resolution'].default_value = 4
        c2m = ng.nodes.new('GeometryNodeCurveToMesh')
        sm = ng.nodes.new('GeometryNodeSetMaterial')
        wm = bpy.data.materials.get('outline') or bpy.data.materials.new('outline')
        wm.use_nodes = True
        bsdf = wm.node_tree.nodes.get('Principled BSDF')
        if bsdf:
            bsdf.inputs['Emission Color'].default_value = (0.9, 0.9, 0.9, 1)
            bsdf.inputs['Emission Strength'].default_value = 1.0
        sm.inputs['Material'].default_value = wm
        ng.links.new(gi.outputs[0], m2c.inputs['Mesh'])
        ng.links.new(m2c.outputs['Curve'], c2m.inputs['Curve'])
        ng.links.new(prof.outputs['Curve'], c2m.inputs['Profile Curve'])
        ng.links.new(c2m.outputs['Mesh'], sm.inputs['Geometry'])
        ng.links.new(sm.outputs['Geometry'], go.inputs[0])
        mod.node_group = ng

sc = bpy.context.scene
cx, cy = (ux[0] + ux[-1]) / 2, (uy[0] + uy[-1]) / 2
span = max(ux[-1] - ux[0], uy[-1] - uy[0])
bpy.ops.object.camera_add(location=(cx - 0.3 * span, cy - 1.05 * span, 0.85 * span + explode))
cam = bpy.context.active_object
cam.data.lens = 40
direction = np.array([cx, cy, explode / 2]) - np.array(cam.location)
cam.rotation_euler = (math.atan2(math.hypot(direction[0], direction[1]), -direction[2]), 0, math.atan2(direction[1], direction[0]) - math.pi / 2)
cam.data.clip_end = 20 * span
sc.camera = cam
sc.world = sc.world or bpy.data.worlds.new('w')
sc.world.color = (0.0, 0.0, 0.0)
sc.frame_start, sc.frame_end = 1, len(fidx)
# frame time stamp
bpy.ops.object.text_add(location=(ux[0] + 0.05 * span, uy[-1] - 0.08 * span, explode * (len(seqs) - 1) + 0.02))
txt = bpy.context.active_object
txt.data.size = span / 25
tm = bpy.data.materials.new('label')
tm.use_nodes = True
tb = tm.node_tree.nodes.get('Principled BSDF')
if tb:
    tb.inputs['Emission Color'].default_value = (1, 1, 1, 1)
    tb.inputs['Emission Strength'].default_value = 1.0
txt.data.materials.append(tm)
times = [f['t_ps'] for f in frames]


def stamp(scene, *_):
    k = min(max(scene.frame_current - 1, 0), len(times) - 1)
    txt.data.body = '|J| on %s   t = %5.0f ps' % (' / '.join(s[2]['name'] for s in reversed(layers_sorted)), times[k])


bpy.app.handlers.frame_change_pre.append(stamp)
stamp(sc)

if opt('--save'):
    bpy.ops.wm.save_as_mainfile(filepath=str(pathlib.Path(opt('--save')).resolve()))
if opt('--render'):
    mp4 = pathlib.Path(opt('--render')).resolve()
    tmp = run / 'frames' / 'render'
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir()
    sc.render.engine = 'BLENDER_EEVEE' if 'BLENDER_EEVEE' in {e.identifier for e in bpy.types.RenderSettings.bl_rna.properties['engine'].enum_items} else 'BLENDER_EEVEE_NEXT'
    sc.render.resolution_x, sc.render.resolution_y = 1280, 720
    sc.render.image_settings.file_format = 'PNG'
    sc.render.filepath = str(tmp / 'f_')
    bpy.ops.render.render(animation=True)
    ff = shutil.which('ffmpeg')
    if ff:
        subprocess.run([ff, '-y', '-loglevel', 'error', '-framerate', '24', '-i', str(tmp / 'f_%04d.png'), '-pix_fmt', 'yuv420p', '-c:v', 'libx264', str(mp4)], check=True)
        print('blender_frames: wrote', mp4)
    else:
        print('blender_frames: ffmpeg not found; PNG frames are in', tmp)
