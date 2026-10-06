"""Python (numpy + ctypes) binding of the gridroute kernels.

    import gridroute as gr
    blk = gr.dilate(occ, layers=[0, 3], excl=[nid], win=(i0, j0, i1, j1), span=gr.disc_span(r_cells))
    isl = gr.island(core, thru, nid, win, seeds)
    path = gr.astar(blk, vok, src_cells, tgt, mcost, lay_ok, vcost, turn, hmul)

The shared library is found through $GRIDROUTE_LIB, else ../../target/release/ next to this package; it is
(re)built with cargo when missing or older than the Rust sources (set GRIDROUTE_NO_BUILD=1 to disable).
`available()` is False when it cannot be loaded, so callers can keep a pure-numpy fallback.
"""
import ctypes
import glob
import math
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CRATE = os.path.normpath(os.path.join(HERE, '..', '..'))
_EXT = {'darwin': 'dylib', 'win32': 'dll'}.get(os.sys.platform, 'so')
_LIBNAME = ('gridroute.dll' if _EXT == 'dll' else 'libgridroute.' + _EXT)

BACKENDS = {'auto': 0, 'cpu': 1, 'metal': 2, 'naive': 3}
BACKEND_NAMES = {v: k for k, v in BACKENDS.items()}

_u8p = ctypes.POINTER(ctypes.c_uint8)
_i16p = ctypes.POINTER(ctypes.c_int16)
_i32p = ctypes.POINTER(ctypes.c_int32)
_f32p = ctypes.POINTER(ctypes.c_float)
_i64p = ctypes.POINTER(ctypes.c_int64)


def _stale(lib):
    srcs = glob.glob(os.path.join(CRATE, 'src', '*.rs')) + glob.glob(os.path.join(CRATE, 'shaders', '*')) + \
        [os.path.join(CRATE, 'Cargo.toml')]
    return not os.path.exists(lib) or os.path.getmtime(lib) < max(os.path.getmtime(s) for s in srcs)


def _load():
    lib = os.environ.get('GRIDROUTE_LIB')
    if not lib:
        lib = os.path.join(CRATE, 'target', 'release', _LIBNAME)
        if os.path.exists(os.path.join(CRATE, 'Cargo.toml')) and _stale(lib) and not os.environ.get('GRIDROUTE_NO_BUILD'):
            try:
                subprocess.run(['cargo', 'build', '--release', '--quiet'], cwd=CRATE, check=True)
            except (OSError, subprocess.CalledProcessError):
                pass
    try:
        so = ctypes.CDLL(lib)
    except OSError:
        return None
    so.gr_version.restype = ctypes.c_uint32
    so.gr_metal_available.restype = ctypes.c_int
    so.gr_device_name.argtypes = [ctypes.c_char_p, ctypes.c_size_t]
    so.gr_device_name.restype = ctypes.c_size_t
    so.gr_set_metal_threshold.argtypes = [ctypes.c_int64]
    so.gr_set_threads.argtypes = [ctypes.c_int]
    so.gr_set_field_factor.argtypes = [ctypes.c_int]
    so.gr_dilate.argtypes = [_i16p, ctypes.c_int, ctypes.c_int, ctypes.c_int, _i32p, ctypes.c_int, _i16p, ctypes.c_int,
                             _u8p, _u8p] + [ctypes.c_int] * 4 + [_i32p, ctypes.c_int, ctypes.c_int, ctypes.c_int, _u8p]
    so.gr_dilate.restype = ctypes.c_int
    so.gr_island.argtypes = [_i16p, ctypes.c_int, ctypes.c_int, ctypes.c_int, _i16p, ctypes.c_int16] + \
        [ctypes.c_int] * 4 + [_i32p, ctypes.c_int, _u8p]
    so.gr_island.restype = ctypes.c_int64
    so.gr_island_crop.argtypes = so.gr_island.argtypes[:-1] + [ctypes.c_int] * 4 + [_u8p]
    so.gr_island_crop.restype = ctypes.c_int64
    so.gr_astar.argtypes = [ctypes.c_int] * 3 + [_u8p] * 3 + [_i32p, ctypes.c_int, _f32p, _u8p] + \
        [ctypes.c_float] * 3 + [ctypes.c_int] * 4 + [ctypes.c_int64, _i32p, ctypes.c_int, _i64p]
    so.gr_astar.restype = ctypes.c_int
    so.gr_astar_hybrid.argtypes = so.gr_astar.argtypes[:-1] + [ctypes.c_int64, ctypes.c_float, ctypes.POINTER(ctypes.c_double), _u8p]
    so.gr_astar_hybrid_cost.argtypes = so.gr_astar_hybrid.argtypes + [_f32p]
    so.gr_astar_hybrid_cost.restype = ctypes.c_int
    so.gr_hash_tiles.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_void_p), _i32p, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                 _i32p, ctypes.POINTER(ctypes.c_uint64)]
    so.gr_astar_hybrid.restype = ctypes.c_int
    so.gr_hash_windows.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_int64)] + \
        [ctypes.c_int64] * 4 + [ctypes.c_uint64, ctypes.POINTER(ctypes.c_uint64)]
    so.gr_hash_buffers.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_int64),
                                   ctypes.c_uint64, ctypes.POINTER(ctypes.c_uint64)]
    so.gr_mask_scan.argtypes = [_u8p, _u8p, ctypes.c_int64, ctypes.c_int, ctypes.c_int, _i32p, ctypes.c_int64, _i32p]
    so.gr_mask_scan.restype = ctypes.c_int64
    so.gr_tile_hashes.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_int64),
                                  ctypes.c_int, ctypes.c_int, _i32p, ctypes.POINTER(ctypes.c_uint64)]
    so.gr_paint_capsule.argtypes = [_i16p] + [ctypes.c_int] * 4 + [ctypes.c_double] * 6 + \
        [ctypes.c_int16, ctypes.c_int, ctypes.c_int16, _i32p]
    so.gr_paint_capsule.restype = ctypes.c_int64
    so.gr_label_cells.argtypes = [_i16p, ctypes.c_int, ctypes.c_int, ctypes.c_int, _i16p, ctypes.c_int16, _i32p,
                                  ctypes.c_int, _i32p]
    so.gr_astar_ref.argtypes = so.gr_astar.argtypes
    so.gr_astar_ref.restype = ctypes.c_int
    return so


_LIB = _load()
stats = {'dilate': {}, 'astar_exp': 0, 'guided': 0, 'field_s': 0.0, 'field_nopath': 0}


def available():
    return _LIB is not None


def metal_available():
    return bool(_LIB and _LIB.gr_metal_available())


def device_name():
    if not _LIB:
        return ''
    buf = ctypes.create_string_buffer(256)
    _LIB.gr_device_name(buf, 256)
    return buf.value.decode()


def set_metal_threshold(work):
    """BACKEND 'auto' runs a dilation on the GPU when h * w * element rows * layers >= work (None: never)."""
    _LIB.gr_set_metal_threshold(-1 if work is None else int(work))


def set_field_factor(f):
    """Hybrid search field resolution: f x f search cells per field cell (1 = exact field)."""
    _LIB.gr_set_field_factor(int(f))


def set_threads(n):
    return bool(_LIB.gr_set_threads(int(n)))


def disc_span(r, strict=True, eps=1e-9):
    """Row half-widths of a disc of radius r (cells) for dilate(). strict: a cell at exactly r is outside
    (row skipped when di^2 >= r^2 - eps, half-width floor(sqrt(r^2 - di^2) - eps)); else inside."""
    rr = r ** 2
    rc = int(math.ceil(r))
    span = np.full(2 * rc + 1, -1, dtype=np.int32)
    for di in range(-rc, rc + 1):
        if strict:
            if di * di >= rr - eps:
                continue
            span[di + rc] = int(math.floor(math.sqrt(rr - di * di) - eps))
        else:
            if di * di > rr:
                continue
            span[di + rc] = int(math.floor(math.sqrt(rr - di * di)))
    return span


def _ptr(a, t):
    return a.ctypes.data_as(t)


def dilate(src, layers, win, span, excl=(), extra=None, extra_on=None, reduce_or=False, backend='auto'):
    """Cells of window win = (i0, j0, i1, j1) (inclusive) within the element `span` of an obstacle.

    src: int16 [nl, ny, nx] (or [ny, nx]); obstacle = value != 0 and not in excl, or extra (bool/uint8 [ny, nx])
    set on a layer whose extra_on[l] is true (default: all layers). Returns bool [len(layers), h, w], or [h, w]
    when reduce_or (OR over the layers)."""
    if src.ndim == 2:
        src = src[None]
    assert src.dtype == np.int16 and src.flags.c_contiguous
    nl, ny, nx = src.shape
    i0, j0, i1, j1 = win
    h, w = i1 - i0 + 1, j1 - j0 + 1
    lay = np.ascontiguousarray(layers, dtype=np.int32)
    ex = np.ascontiguousarray(excl, dtype=np.int16)
    sp = np.ascontiguousarray(span, dtype=np.int32)
    planes = 1 if reduce_or else len(lay)
    out = np.empty((planes, h, w), dtype=np.uint8)
    xp = eo = None
    if extra is not None:
        xp = np.ascontiguousarray(extra).view(np.uint8) if extra.dtype == np.bool_ else np.ascontiguousarray(extra, dtype=np.uint8)
        eo = np.ones(nl, dtype=np.uint8) if extra_on is None else np.ascontiguousarray(extra_on, dtype=np.uint8)
    rc = _LIB.gr_dilate(_ptr(src, _i16p), nl, ny, nx, _ptr(lay, _i32p), len(lay), _ptr(ex, _i16p) if len(ex) else None,
                        len(ex), _ptr(xp, _u8p) if xp is not None else None, _ptr(eo, _u8p) if eo is not None else None,
                        i0, j0, h, w, _ptr(sp, _i32p), len(sp), int(bool(reduce_or)), BACKENDS[backend],
                        _ptr(out, _u8p))
    if rc < 0:
        raise ValueError('gr_dilate failed: %d' % rc)
    b = BACKEND_NAMES[rc]
    stats['dilate'][b] = stats['dilate'].get(b, 0) + 1
    out = out.view(np.bool_)
    return out[0] if reduce_or else out


def island(core, thru, nid, win, seeds):
    """bool [nl, h, w]: window cells of net nid connected to seeds ((layer, y, x) relative to the window) through
    4-neighbour copper (core int16 [nl, ny, nx]) and holes (thru int16 [ny, nx])."""
    assert core.dtype == np.int16 and thru.dtype == np.int16 and core.flags.c_contiguous and thru.flags.c_contiguous
    nl, ny, nx = core.shape
    i0, j0, i1, j1 = win
    h, w = i1 - i0 + 1, j1 - j0 + 1
    sd = np.ascontiguousarray(np.asarray(seeds, dtype=np.int32).reshape(-1, 3))
    out = np.empty((nl, h, w), dtype=np.uint8)
    n = _LIB.gr_island(_ptr(core, _i16p), nl, ny, nx, _ptr(thru, _i16p), int(nid), i0, j0, h, w,
                       _ptr(sd, _i32p), len(sd), _ptr(out, _u8p))
    if n < 0:
        raise ValueError('gr_island failed: %d' % n)
    return out.view(np.bool_)


def paint_capsule(arr, layer, a, b, r, pitch, nid, keep=None):
    """Paint segment a-b of radius r (board units, grid pitch `pitch`) into int16 arr [nl, ny, nx] (or [ny, nx]):
    empty cells take nid; with keep (a value), cells of another value become keep. Returns the painted window
    (i0, j0, i1, j1). Bit-identical to the numpy capsule mask (cell point within r + 1e-9)."""
    a3 = arr if arr.ndim == 3 else arr[None]
    assert a3.dtype == np.int16 and a3.flags.c_contiguous
    nl, ny, nx = a3.shape
    w = (ctypes.c_int32 * 4)()
    _LIB.gr_paint_capsule(a3.ctypes.data_as(_i16p), nl, ny, nx, 0 if arr.ndim == 2 else int(layer), float(a[0]),
                          float(a[1]), float(b[0]), float(b[1]), float(r), float(pitch), int(nid),
                          0 if keep is None else 1, 0 if keep is None else int(keep), w)
    return tuple(w)


def mask_scan(a, b=None, bbox=False):
    """Flat indices (int32, in order) of a & ~b for bool arrays of the same shape (b optional); with bbox=True also
    the (x0, y0, x1, y1) box over the last two axes (None if empty). Native and parallel."""
    a8 = np.ascontiguousarray(a).view(np.uint8)
    b8 = None if b is None else np.ascontiguousarray(b).view(np.uint8)
    h, w = (a.shape[-2], a.shape[-1]) if a.ndim >= 2 else (1, a.size)
    cap = 1 << 16
    while True:
        out = np.empty(cap, dtype=np.int32)
        bb = np.empty(4, dtype=np.int32)
        n = _LIB.gr_mask_scan(_ptr(a8, _u8p), None if b8 is None else _ptr(b8, _u8p), a8.size, h, w, _ptr(out, _i32p),
                              cap, _ptr(bb, _i32p))
        if n <= cap:
            break
        cap = int(n)
    idx = out[:n]
    if bbox:
        return idx, (None if n == 0 else tuple(int(v) for v in bb))
    return idx


def label_cells(core, thru, nid, cells):
    """Copper component labels of net nid at (layer, y, x) cells over the whole grid: equal labels = connected,
    -1 = not the net's copper. One native pass over the net's copper reachable from the cells."""
    nl, ny, nx = core.shape
    c = np.ascontiguousarray(np.asarray(cells, dtype=np.int32).reshape(-1, 3))
    out = np.empty(len(c), dtype=np.int32)
    _LIB.gr_label_cells(_ptr(core, _i16p), nl, ny, nx, _ptr(thru, _i16p), int(nid), _ptr(c, _i32p), len(c), _ptr(out, _i32p))
    return out


def island_crop(core, thru, nid, win, seeds, crop):
    """island() traced over win but returned only for the crop window (i0, j0, i1, j1) inside it: bool [nl, h, w]."""
    nl, ny, nx = core.shape
    i0, j0, i1, j1 = win
    c0, d0, c1, d1 = crop
    sd = np.ascontiguousarray(np.asarray(seeds, dtype=np.int32).reshape(-1, 3))
    out = np.empty((nl, c1 - c0 + 1, d1 - d0 + 1), dtype=np.uint8)
    n = _LIB.gr_island_crop(_ptr(core, _i16p), nl, ny, nx, _ptr(thru, _i16p), int(nid), i0, j0, i1 - i0 + 1, j1 - j0 + 1,
                            _ptr(sd, _i32p), len(sd), c0, d0, c1 - c0 + 1, d1 - d0 + 1, _ptr(out, _u8p))
    if n < 0:
        raise ValueError('gr_island_crop failed: %d' % n)
    return out.view(np.bool_)


def astar(blk, vok, src, tgt, mcost, lay_ok, vcost, turn, hmul, max_exp, tbox=None, cap=1 << 22, reference=False,
          budget=None, weight=1.2, read_set=False, cost=None):
    """Multi-layer A* (see src/astar.rs). budget: hybrid search (gr_astar_hybrid): after `budget` expansions an
    open search switches to the GPU cost-to-target field and field-guided A* with `weight`. blk, tgt: bool [nl, h, w]; vok: bool [h, w]; src: int array of flat state
    indices; mcost: float32 [nl, 8]; lay_ok: bool [nl]. tbox: target bounding box (x0, y0, x1, y1), default from tgt.
    cost: optional float32 [nl, h, w], an extra cost (>= 0) for entering each state (always a hybrid search).
    Returns the flat state indices of the path (source first) or None."""
    nl, h, w = blk.shape
    b8 = np.ascontiguousarray(blk).view(np.uint8)
    v8 = np.ascontiguousarray(vok).view(np.uint8)
    t8 = np.ascontiguousarray(tgt).view(np.uint8)
    s32 = np.ascontiguousarray(src, dtype=np.int32)
    mc = np.ascontiguousarray(mcost, dtype=np.float32).ravel()
    lok = np.ascontiguousarray(lay_ok, dtype=np.uint8)
    if tbox is None:
        tbox = mask_scan(tgt, bbox=True)[1]
    out = np.empty(cap, dtype=np.int32)
    nexp = ctypes.c_int64(0)
    args = (nl, h, w, _ptr(b8, _u8p), _ptr(v8, _u8p), _ptr(t8, _u8p), _ptr(s32, _i32p), len(s32),
            _ptr(mc, _f32p), _ptr(lok, _u8p), float(vcost), float(turn), float(hmul), *tbox, int(max_exp), _ptr(out, _i32p), cap)
    if cost is not None and not reference and budget is None:
        budget = max_exp
    if budget is not None and not reference:
        info = (ctypes.c_double * 3)()
        touched = np.zeros(((h + 31) // 32, (w + 31) // 32), dtype=np.uint8) if read_set else None
        if cost is not None:
            c32 = np.ascontiguousarray(cost, dtype=np.float32)
            assert c32.shape == blk.shape
            n = _LIB.gr_astar_hybrid_cost(*args, int(budget), float(weight), info,
                                          _ptr(touched, _u8p) if read_set else None, _ptr(c32, _f32p))
        else:
            n = _LIB.gr_astar_hybrid(*args, int(budget), float(weight), info,
                                     _ptr(touched, _u8p) if read_set else None)
        stats['astar_exp'] += int(info[2])
        stats['guided'] += info[0] == 2
        stats['field_nopath'] += info[0] == 3
        stats['field_s'] += info[1]
    else:
        fn = _LIB.gr_astar_ref if reference else _LIB.gr_astar
        n = fn(*args, ctypes.byref(nexp))
        stats['astar_exp'] += nexp.value
    if n == -2:
        return astar(blk, vok, src, tgt, mcost, lay_ok, vcost, turn, hmul, max_exp, tbox, cap * 8, reference, budget, weight,
                     read_set, cost)
    path = out[:n].copy() if n >= 0 else None
    if read_set and budget is not None and not reference:
        return path, np.argwhere(touched).astype(np.int32)
    return path


def hash_windows(arrays, win, seed=0):
    """128-bit content hash (hex string) of window win = (i0, j0, i1, j1) (inclusive, clipped) of each C-contiguous
    array ([ny, nx] or [nl, ny, nx]); chained from seed. A few GB/s: cheap next to rebuilding masks or searching."""
    ptrs = (ctypes.c_void_p * len(arrays))()
    dims = (ctypes.c_int64 * (4 * len(arrays)))()
    for k, a in enumerate(arrays):
        assert a.flags.c_contiguous
        shp = a.shape if a.ndim == 3 else (1,) + a.shape
        ptrs[k] = a.ctypes.data
        dims[4 * k:4 * k + 4] = [shp[0], shp[1], shp[2], a.itemsize]
    out = (ctypes.c_uint64 * 2)()
    i0, j0, i1, j1 = win
    _LIB.gr_hash_windows(len(arrays), ptrs, dims, i0, j0, i1, j1, seed, out)
    return '%016x%016x' % (out[0], out[1])


def tile_hashes(arrays, tile, tiles):
    """uint64 [len(tiles), 2]: 128-bit hash of each (ty, tx) tile (tile x tile cells) across C-contiguous arrays
    ([ny, nx] or [nl, ny, nx], same ny, nx). Parallel; for incremental content hashing of a grid."""
    tl = np.ascontiguousarray(tiles, dtype=np.int32).reshape(-1, 2)
    ptrs = (ctypes.c_void_p * len(arrays))()
    dims = (ctypes.c_int64 * (4 * len(arrays)))()
    for k, a in enumerate(arrays):
        assert a.flags.c_contiguous
        shp = a.shape if a.ndim == 3 else (1,) + a.shape
        ptrs[k] = a.ctypes.data
        dims[4 * k:4 * k + 4] = [shp[0], shp[1], shp[2], a.itemsize]
    out = np.empty((len(tl), 2), dtype=np.uint64)
    _LIB.gr_tile_hashes(len(arrays), ptrs, dims, int(tile), len(tl), _ptr(tl, _i32p),
                        out.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)))
    return out


def hash_arrays(arrays, seed=0):
    """128-bit content hash (hex) of whole arrays (made contiguous; dtype and shape included), parallel for large ones."""
    arrs = [np.ascontiguousarray(a) for a in arrays]
    head = repr([(a.dtype.str, a.shape) for a in arrs]).encode()
    bufs = [np.frombuffer(head, dtype=np.uint8)] + [a.reshape(-1).view(np.uint8) for a in arrs]
    ptrs = (ctypes.c_void_p * len(bufs))(*[b.ctypes.data for b in bufs])
    lens = (ctypes.c_int64 * len(bufs))(*[b.nbytes for b in bufs])
    out = (ctypes.c_uint64 * 2)()
    _LIB.gr_hash_buffers(len(bufs), ptrs, lens, seed, out)
    return '%016x%016x' % (out[0], out[1])


def hash_tiles(arrays, tiles):
    """128-bit hash (hex) of the listed 32x32 tiles ((ty, tx) rows) of byte/bool arrays ([h, w] or [nl, h, w], all
    the same h x w): the fingerprint of a search's read set."""
    arrs = [np.ascontiguousarray(a if a.ndim == 3 else a[None]).view(np.uint8) for a in arrays]
    h, w = arrs[0].shape[1:]
    ptrs = (ctypes.c_void_p * len(arrs))(*[a.ctypes.data for a in arrs])
    nls = np.array([a.shape[0] for a in arrs], dtype=np.int32)
    tl = np.ascontiguousarray(tiles, dtype=np.int32).reshape(-1, 2)
    out = (ctypes.c_uint64 * 2)()
    _LIB.gr_hash_tiles(len(arrs), ptrs, _ptr(nls, _i32p), h, w, len(tl), _ptr(tl, _i32p), out)
    return '%016x%016x' % (out[0], out[1])


class SearchCache:
    """Persistent memo of search results keyed by content hashes (see hash_windows / hash_arrays): a re-run of a
    deterministic router only repeats the searches whose inputs changed. Stored with pickle at `path` (loaded on
    creation, written by save()). Entries not used for `keep` runs are dropped, so it follows the design."""

    def __init__(self, path, keep=4):
        import pickle
        self.path, self.keep, self.hits, self.misses = path, keep, 0, 0
        try:
            with open(path, 'rb') as f:
                self.gen, self.data = pickle.load(f)
        except (OSError, EOFError, pickle.UnpicklingError, AttributeError, ValueError, TypeError):
            self.gen, self.data = 0, {}
        self.gen += 1

    def get(self, key, valid=None):
        """The value for key; for a list value with `valid`, the first item for which valid(item) holds."""
        e = self.data.get(key)
        v = None
        if e is not None:
            e[0] = self.gen
            v = e[1] if valid is None else next((x for x in e[1] if valid(x)), None)
        if v is None:
            self.misses += 1
        else:
            self.hits += 1
        return v

    def peek(self, key):
        e = self.data.get(key)
        return None if e is None else e[1]

    def put(self, key, value):
        self.data[key] = [self.gen, value]

    def save(self):
        import pickle
        os.makedirs(os.path.dirname(self.path) or '.', exist_ok=True)
        live = {k: e for k, e in self.data.items() if e[0] > self.gen - self.keep}
        tmp = self.path + '.tmp'
        with open(tmp, 'wb') as f:
            pickle.dump((self.gen, live), f, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, self.path)
