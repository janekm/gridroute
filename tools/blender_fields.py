"""Load fieldsolve field snapshots into Blender.

    /Applications/Blender.app/Contents/MacOS/Blender --python tools/blender_fields.py -- RUN_DIR/fields \
        [--freq 1.0000GHz] [--field E|H] [--layers F.Cu,In1.Cu] [--no-copper] [--no-volume] [--top]
        [--lines lines_H_<f>.json] [--lines-mode volume|tubes|both] [--fog-cut 0.5] [--frame-all] [--save out.blend] [--render out.png [--engine CYCLES|BLENDER_EEVEE] [--samples 64]]

Builds, in millimetre scene units (Gerber orientation, z up):
- one object per copper layer from copper_<layer>.obj (dim copper colour),
- one emissive plane per copper layer showing |J_s| (surface current density, log scale over 40 dB),
- the |E| (or |H|) magnitude as an OpenVDB volume (resampled to the finest cell; log scale, 50 dB).

- with --lines (from `fieldsolve lines`): field lines as an emissive fog volume (Geometry Nodes: points ->
  volume, thicker where the field is stronger) and/or glowing tubes; with several phases, one object per
  phase, shown on successive frames (an animation of the field over one period).

Uses only what Blender ships (numpy, openvdb).
"""
import json
import math
import pathlib
import sys

import bpy
import numpy as np

argv = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else []
if not argv:
    raise SystemExit('usage: blender --python blender_fields.py -- RUN_DIR/fields [--freq TAG] [--field E|H]')
d = pathlib.Path(argv[0])


def opt(k, default=None):
    return argv[argv.index(k) + 1] if k in argv else default


meta = json.loads((d / 'meta.json').read_text())
tag = opt('--freq', meta['snapshots'][0]['tag'])
field = opt('--field', 'E')
nx, ny, nz = meta['shape']
X, Y, Z = (np.array(meta[a]) for a in 'xyz')

for o in list(bpy.data.objects):
    bpy.data.objects.remove(o, do_unlink=True)


def emissive_ramp_material(name, image=None, attribute=None, volume=False, alpha=False):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    nt = m.node_tree
    for n in list(nt.nodes):
        nt.nodes.remove(n)
    out = nt.nodes.new('ShaderNodeOutputMaterial')
    ramp = nt.nodes.new('ShaderNodeValToRGB')
    cr = ramp.color_ramp
    cr.elements[0].position, cr.elements[0].color = 0.0, (0.0, 0.0, 0.05, 1)
    cr.elements[1].position, cr.elements[1].color = 1.0, (1.0, 1.0, 0.85, 1)
    for pos, col in [(0.35, (0.25, 0.0, 0.55, 1)), (0.6, (0.9, 0.2, 0.1, 1)), (0.8, (1.0, 0.7, 0.1, 1))]:
        e = cr.elements.new(pos)
        e.color = col
    if volume:
        attr = nt.nodes.new('ShaderNodeAttribute')
        attr.attribute_name = attribute
        vol = nt.nodes.new('ShaderNodeVolumePrincipled')
        nt.links.new(attr.outputs['Fac'], ramp.inputs['Fac'])
        nt.links.new(ramp.outputs['Color'], vol.inputs['Emission Color'])
        nt.links.new(attr.outputs['Fac'], vol.inputs['Emission Strength'])
        nt.links.new(attr.outputs['Fac'], vol.inputs['Density'])
        vol.inputs['Color'].default_value = (0, 0, 0, 1)
        nt.links.new(vol.outputs['Volume'], out.inputs['Volume'])
        return m
    tex = nt.nodes.new('ShaderNodeTexImage')
    tex.image = image
    tex.interpolation = 'Closest'
    em = nt.nodes.new('ShaderNodeEmission')
    tr = nt.nodes.new('ShaderNodeBsdfTransparent')
    mix = nt.nodes.new('ShaderNodeMixShader')
    nt.links.new(tex.outputs['Color'], ramp.inputs['Fac'])
    nt.links.new(ramp.outputs['Color'], em.inputs['Color'])
    em.inputs['Strength'].default_value = 2.0
    nt.links.new(tex.outputs['Color'], mix.inputs['Fac'])
    nt.links.new(tr.outputs['BSDF'], mix.inputs[1])
    nt.links.new(em.outputs['Emission'], mix.inputs[2])
    nt.links.new(mix.outputs['Shader'], out.inputs['Surface'])
    m.blend_method = 'BLEND' if hasattr(m, 'blend_method') else None
    return m


def uniform_axis(lines, h):
    n = max(2, int(round((lines[-1] - lines[0]) / h)) + 1)
    return np.linspace(lines[0], lines[-1], n)


def resample_idx(lines, u):
    # nearest original node for each uniform coordinate
    j = np.clip(np.searchsorted(lines, u), 1, len(lines) - 1)
    return np.where(np.abs(lines[j - 1] - u) <= np.abs(lines[j] - u), j - 1, j)


def log_norm(a, db):
    m = float(a.max()) or 1.0
    v = 20.0 * np.log10(np.maximum(a / m, 1e-12)) if db > 0 else a / m
    return np.clip((v + db) / db, 0.0, 1.0).astype(np.float32)


copper = bpy.data.materials.new('copper')
copper.diffuse_color = (0.72, 0.45, 0.2, 0.35)
h_xy = min(np.diff(X).min(), np.diff(Y).min())
ux, uy = uniform_axis(X, h_xy), uniform_axis(Y, h_xy)
ix, iy = resample_idx(X, ux), resample_idx(Y, uy)

only = set(opt('--layers').split(',')) if opt('--layers') else None
for L in meta['layers']:
    if only and L['name'] not in only:
        continue
    obj_path = d / ('copper_%s.obj' % L['tag'])
    if obj_path.exists() and '--no-copper' not in argv:
        bpy.ops.wm.obj_import(filepath=str(obj_path), forward_axis='Y', up_axis='Z')
        o = bpy.context.selected_objects[0]
        o.name = 'copper ' + L['name']
        o.data.materials.append(copper)
    jp = d / ('J_%s_%s.f32' % (L['tag'], tag))
    if not jp.exists():
        continue
    J = np.fromfile(jp, dtype='<f4').reshape(ny, nx)
    Ju = log_norm(J[np.ix_(iy, ix)], 40.0)
    img = bpy.data.images.new('J ' + L['name'], width=len(ux), height=len(uy), float_buffer=True)
    px = np.repeat(Ju.reshape(-1, 1), 4, axis=1)
    px[:, 3] = 1.0
    img.pixels.foreach_set(px.ravel())
    img.pack()
    bpy.ops.mesh.primitive_plane_add(size=1.0, location=((ux[0] + ux[-1]) / 2, (uy[0] + uy[-1]) / 2, L['z'] + 0.002))
    pl = bpy.context.active_object
    pl.name = 'J ' + L['name']
    pl.scale = (ux[-1] - ux[0], uy[-1] - uy[0], 1)
    pl.data.materials.append(emissive_ramp_material('J ' + L['name'], image=img))

# volume
vol_path = d / ('%s_abs_%s.f32' % (field, tag))
if vol_path.exists() and '--no-volume' not in argv:
    import openvdb as vdb
    A = np.fromfile(vol_path, dtype='<f4').reshape(nz, ny, nx)
    h_z = max(np.diff(Z).min(), h_xy / 2)
    uz = uniform_axis(Z, h_z)
    iz = resample_idx(Z, uz)
    Au = log_norm(A[np.ix_(iz, iy, ix)], 50.0)
    grid = vdb.FloatGrid()
    grid.copyFromArray(np.ascontiguousarray(Au.transpose(2, 1, 0)))  # [i, j, k]
    sx, sy, sz = ux[1] - ux[0], uy[1] - uy[0], uz[1] - uz[0]
    grid.transform = vdb.createLinearTransform([[sx, 0, 0, 0], [0, sy, 0, 0], [0, 0, sz, 0], [ux[0], uy[0], uz[0], 1]])
    grid.name = 'density'
    vp = d / ('%s_abs_%s.vdb' % (field, tag))
    vdb.write(str(vp), grids=[grid])
    bpy.ops.object.volume_import(filepath=str(vp))
    v = bpy.context.active_object
    v.name = '|%s| %s' % (field, tag)
    v.data.materials.append(emissive_ramp_material('|%s|' % field, attribute='density', volume=True))


def lines_material(name, volume):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    nt = m.node_tree
    for n in list(nt.nodes):
        nt.nodes.remove(n)
    out = nt.nodes.new('ShaderNodeOutputMaterial')
    if volume:
        vol = nt.nodes.new('ShaderNodeVolumePrincipled')
        vol.inputs['Color'].default_value = (0, 0, 0, 1)
        vol.inputs['Emission Color'].default_value = (0.35, 0.65, 1.0, 1)
        # density and glow from the grid itself (a constant would fill the volume's bounds uniformly)
        attr = nt.nodes.new('ShaderNodeAttribute')
        attr.attribute_name = 'density'
        dm = nt.nodes.new('ShaderNodeMath')
        dm.operation = 'MULTIPLY'
        dm.inputs[1].default_value = 2.0
        em = nt.nodes.new('ShaderNodeMath')
        em.operation = 'MULTIPLY'
        em.inputs[1].default_value = 20.0
        nt.links.new(attr.outputs['Fac'], dm.inputs[0])
        nt.links.new(attr.outputs['Fac'], em.inputs[0])
        nt.links.new(dm.outputs[0], vol.inputs['Density'])
        nt.links.new(em.outputs[0], vol.inputs['Emission Strength'])
        nt.links.new(vol.outputs['Volume'], out.inputs['Volume'])
    else:
        attr = nt.nodes.new('ShaderNodeAttribute')
        attr.attribute_name = 'mag'
        ramp = nt.nodes.new('ShaderNodeValToRGB')
        cr = ramp.color_ramp
        cr.elements[0].color = (0.05, 0.1, 0.6, 1)
        cr.elements[1].color = (0.9, 1.0, 1.0, 1)
        em = nt.nodes.new('ShaderNodeEmission')
        nt.links.new(attr.outputs['Fac'], ramp.inputs['Fac'])
        nt.links.new(ramp.outputs['Color'], em.inputs['Color'])
        math_ = nt.nodes.new('ShaderNodeMath')
        math_.operation = 'MULTIPLY'
        math_.inputs[1].default_value = 4.0
        nt.links.new(attr.outputs['Fac'], math_.inputs[0])
        nt.links.new(math_.outputs[0], em.inputs['Strength'])
        nt.links.new(em.outputs['Emission'], out.inputs['Surface'])
    return m


def lines_nodes(mode, radius, voxel, fog_cut=0.5):
    ng = bpy.data.node_groups.new('field_lines_' + mode, 'GeometryNodeTree')
    ng.interface.new_socket('Geometry', in_out='INPUT', socket_type='NodeSocketGeometry')
    ng.interface.new_socket('Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')
    N = ng.nodes
    L = ng.links
    gi, go = N.new('NodeGroupInput'), N.new('NodeGroupOutput')
    mag = N.new('GeometryNodeInputNamedAttribute')
    mag.data_type = 'FLOAT'
    mag.inputs['Name'].default_value = 'mag'
    rad = N.new('ShaderNodeMath')
    rad.operation = 'MULTIPLY_ADD'  # radius * (0.3 + mag)
    rad.inputs[1].default_value = radius
    rad.inputs[2].default_value = 0.3 * radius
    L.new(mag.outputs['Attribute'], rad.inputs[0])
    join = N.new('GeometryNodeJoinGeometry')
    if mode in ('volume', 'both'):
        # only the stronger part of the field becomes fog (mag is 0..1 over floor_db)
        cmp_ = N.new('FunctionNodeCompare')
        cmp_.data_type = 'FLOAT'
        cmp_.operation = 'LESS_THAN'
        cmp_.inputs['B'].default_value = fog_cut
        L.new(mag.outputs['Attribute'], cmp_.inputs['A'])
        dele = N.new('GeometryNodeDeleteGeometry')
        dele.domain = 'POINT'
        L.new(gi.outputs[0], dele.inputs['Geometry'])
        L.new(cmp_.outputs['Result'], dele.inputs['Selection'])
        pts = N.new('GeometryNodeMeshToPoints')
        pts.mode = 'VERTICES'
        L.new(dele.outputs['Geometry'], pts.inputs['Mesh'])
        L.new(rad.outputs[0], pts.inputs['Radius'])
        p2v = N.new('GeometryNodePointsToVolume')
        # Blender 4.x: a node property; 5.x: a menu input socket
        if hasattr(p2v, 'resolution_mode'):
            p2v.resolution_mode = 'VOXEL_SIZE'
        elif 'Resolution Mode' in p2v.inputs:
            p2v.inputs['Resolution Mode'].default_value = 'Size'
        p2v.inputs['Voxel Size'].default_value = voxel
        p2v.inputs['Density'].default_value = 1.0
        L.new(pts.outputs['Points'], p2v.inputs['Points'])
        L.new(rad.outputs[0], p2v.inputs['Radius'])
        sm = N.new('GeometryNodeSetMaterial')
        sm.inputs['Material'].default_value = lines_material('field lines (volume)', True)
        L.new(p2v.outputs['Volume'], sm.inputs['Geometry'])
        L.new(sm.outputs['Geometry'], join.inputs['Geometry'])
    if mode in ('tubes', 'both'):
        m2c = N.new('GeometryNodeMeshToCurve')
        L.new(gi.outputs[0], m2c.inputs['Mesh'])
        prof = N.new('GeometryNodeCurvePrimitiveCircle')
        prof.inputs['Resolution'].default_value = 6
        prof.inputs['Radius'].default_value = 0.6 * radius
        c2m = N.new('GeometryNodeCurveToMesh')
        L.new(m2c.outputs['Curve'], c2m.inputs['Curve'])
        L.new(prof.outputs['Curve'], c2m.inputs['Profile Curve'])
        sm = N.new('GeometryNodeSetMaterial')
        sm.inputs['Material'].default_value = lines_material('field lines (tubes)', False)
        L.new(c2m.outputs['Mesh'], sm.inputs['Geometry'])
        L.new(sm.outputs['Geometry'], join.inputs['Geometry'])
    L.new(join.outputs['Geometry'], go.inputs[0])
    return ng


if opt('--lines'):
    lp = pathlib.Path(opt('--lines'))
    if not lp.is_absolute():
        lp = d / lp
    LD = json.loads(lp.read_text())
    mode = opt('--lines-mode', 'volume')
    floor_db = LD.get('floor_db', -40.0)
    radius = float(opt('--line-radius', 0.15 * h_xy))
    ng = lines_nodes(mode, radius, voxel=radius / 1.5, fog_cut=float(opt('--fog-cut', 0.5)))
    nph = len(LD['phases'])
    for pi, PH in enumerate(LD['phases']):
        verts, edges, mags = [], [], []
        for ln in PH['lines']:
            P = np.array(ln['p']).reshape(-1, 3)
            M = np.array(ln['m'])
            base = len(verts)
            verts.extend(P.tolist())
            mags.extend(np.clip((20 * np.log10(np.maximum(M / LD['max'], 1e-12)) - floor_db) / -floor_db, 0, 1).tolist())
            edges.extend([(base + k, base + k + 1) for k in range(len(P) - 1)])
        me = bpy.data.meshes.new('lines %d' % pi)
        me.from_pydata(verts, edges, [])
        a = me.attributes.new('mag', 'FLOAT', 'POINT')
        a.data.foreach_set('value', np.array(mags, dtype=np.float32))
        ob = bpy.data.objects.new('%s lines %s ph%d' % (LD['field'], LD['freq'], pi), me)
        bpy.context.scene.collection.objects.link(ob)
        mod = ob.modifiers.new('field lines', 'NODES')
        mod.node_group = ng
        if nph > 1:
            # visible on frame pi + 1 only
            for f, vis in [(pi, False), (pi + 1, True), (pi + 2, False)]:
                if f < 1:
                    continue
                ob.hide_render = ob.hide_viewport = not vis
                ob.keyframe_insert('hide_render', frame=f)
                ob.keyframe_insert('hide_viewport', frame=f)
    if nph > 1:
        bpy.context.scene.frame_start, bpy.context.scene.frame_end = 1, nph
        bpy.context.scene.frame_set(1)

# camera and light for a quick look
cx, cy = (X[0] + X[-1]) / 2, (Y[0] + Y[-1]) / 2
span = max(X[-1] - X[0], Y[-1] - Y[0])
if opt('--lines') and '--frame-all' not in argv:
    # frame the field lines (where the action is) rather than the whole domain with its PML
    allp = np.concatenate([np.array(ln['p']).reshape(-1, 3) for PH in LD['phases'] for ln in PH['lines']])
    lo_, hi_ = allp.min(axis=0), allp.max(axis=0)
    cx, cy = (lo_[0] + hi_[0]) / 2, (lo_[1] + hi_[1]) / 2
    span = 1.6 * max(hi_[0] - lo_[0], hi_[1] - lo_[1])
if '--top' in argv:
    bpy.ops.object.camera_add(location=(cx, cy, Z[-1] + span), rotation=(0, 0, 0))
    cam = bpy.context.active_object
    cam.data.type = 'ORTHO'
    cam.data.ortho_scale = span
else:
    bpy.ops.object.camera_add(location=(cx, cy - 0.9 * span, 0.9 * span), rotation=(math.radians(45), 0, 0))
    cam = bpy.context.active_object
cam.data.clip_end = 10 * span
bpy.context.scene.camera = cam
bpy.context.scene.world = bpy.context.scene.world or bpy.data.worlds.new('w')
bpy.context.scene.world.color = (0.01, 0.01, 0.015)

if opt('--save'):
    bpy.ops.wm.save_as_mainfile(filepath=str(pathlib.Path(opt('--save')).resolve()))
# Eevee's volume froxels must span the scene, not the default 0.1..100 units, or mm-sized fog turns blocky
sc = bpy.context.scene
if hasattr(sc, 'eevee'):
    for k, v in [('volumetric_start', 0.2 * span), ('volumetric_end', 3.0 * span), ('volumetric_tile_size', '2'), ('volumetric_samples', 128)]:
        try:
            setattr(sc.eevee, k, v)
        except Exception:
            pass
if opt('--render'):
    engine = opt('--engine', 'CYCLES')
    sc.render.engine = engine
    if engine == 'CYCLES':
        sc.cycles.samples = int(opt('--samples', 64))
        sc.cycles.volume_step_rate = 0.25
        try:
            prefs = bpy.context.preferences.addons['cycles'].preferences
            prefs.compute_device_type = 'METAL'
            prefs.get_devices()
            for dv in prefs.devices:
                dv.use = True
            sc.cycles.device = 'GPU'
        except Exception:
            pass
    sc.render.resolution_x, sc.render.resolution_y = 1280, 900
    sc.render.filepath = str(pathlib.Path(opt('--render')).resolve())
    bpy.ops.render.render(write_still=True)
