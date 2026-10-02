// Cost-to-target field on the GPU (see src/field.rs for the definition): chaotic tiled Bellman-Ford.
//
// The window is cut into TS x TS tiles; a threadgroup owns one tile (all layers) and relaxes it in threadgroup
// memory until it stops changing, reading a one-cell halo of neighbour tiles from device memory. A tile runs in a
// round only if it or a neighbour changed in the previous round, so the work follows the wavefront. Values only
// ever decrease and every value written is the cost of a real path, so races between neighbouring tiles in the
// same round are harmless; the result is the exact fixpoint once a round changes nothing.
#include <metal_stdlib>
using namespace metal;

constant constexpr int TS = 32;
constant constexpr int TP = TS + 2;
constant constexpr int NLMAX = 6;
constant constexpr float INF = INFINITY;
constant int DI[8] = {0, 0, 1, -1, 1, 1, -1, -1};
constant int DJ[8] = {1, -1, 0, 0, 1, -1, 1, -1};

struct FieldParams {
    uint nl, h, w;
    uint tiles_x, tiles_y;
    uint lay_ok;      // bit mask
    float vcost;
    uint maxit;       // local iterations per round
    uint round;
};

kernel void field_init(device float *F [[buffer(0)]],
                       device const uchar *blk [[buffer(1)]],
                       device const uchar *tgt [[buffer(2)]],
                       constant uint &n [[buffer(3)]],
                       uint id [[thread_position_in_grid]])
{
    if (id >= n) return;
    F[id] = (tgt[id] && !blk[id]) ? 0.0f : INF;
}

// Active tiles of a round (this tile or a neighbour changed in the previous round; round 0: all) appended to
// `list`; args[0] counts them and is the indirect dispatch's threadgroup count (reset by field_reset).
kernel void field_reset(device atomic_uint *args [[buffer(0)]], uint id [[thread_position_in_grid]])
{
    if (id == 0) atomic_store_explicit(&args[0], 0u, memory_order_relaxed);
}

kernel void field_collect(constant FieldParams &p [[buffer(0)]],
                          device const uint *last_change [[buffer(1)]],
                          device atomic_uint *args [[buffer(2)]],
                          device uint *list [[buffer(3)]],
                          uint id [[thread_position_in_grid]])
{
    uint nt = p.tiles_x * p.tiles_y;
    if (id >= nt) return;
    int tx = int(id % p.tiles_x), ty = int(id / p.tiles_x);
    bool act = p.round == 0;
    for (int dy = -1; dy <= 1 && !act; dy++)
        for (int dx = -1; dx <= 1; dx++) {
            int ax = tx + dx, ay = ty + dy;
            if (ax >= 0 && ay >= 0 && ax < int(p.tiles_x) && ay < int(p.tiles_y) &&
                last_change[ay * p.tiles_x + ax] >= p.round)
                act = true;
        }
    if (act) list[atomic_fetch_add_explicit(&args[0], 1u, memory_order_relaxed)] = id;
}

kernel void field_round(device float *F [[buffer(0)]],
                        device const uchar *blk [[buffer(1)]],
                        device const uchar *vok [[buffer(2)]],
                        constant float *mcost [[buffer(3)]],
                        constant FieldParams &p [[buffer(4)]],
                        device uint *last_change [[buffer(5)]],
                        device atomic_uint *changed [[buffer(6)]],
                        device const uint *list [[buffer(7)]],
                        uint2 tg [[threadgroup_position_in_grid]],
                        uint2 lt [[thread_position_in_threadgroup]],
                        uint li [[thread_index_in_threadgroup]])
{
    threadgroup float sh[NLMAX][TP][TP];
    threadgroup uchar sblk[TP][TP];       // bit l: blocked on layer l (out of window: all set)
    threadgroup int tflag;

    uint tile = list[tg.x];
    int tx = int(tile % p.tiles_x), ty = int(tile / p.tiles_x);
    int H = int(p.h), W = int(p.w), NL = int(p.nl);
    ulong HW = ulong(H) * ulong(W);
    int y0 = ty * TS - 1, x0 = tx * TS - 1;
    for (int k = int(li); k < TP * TP; k += TS * TS) {
        int yy = k / TP, xx = k % TP;
        int y = y0 + yy, x = x0 + xx;
        bool in = y >= 0 && x >= 0 && y < H && x < W;
        uchar m = 0;
        for (int l = 0; l < NLMAX; l++) {
            float v = INF;
            if (l < NL) {
                if (in) {
                    ulong s = ulong(l) * HW + ulong(y) * ulong(W) + ulong(x);
                    v = F[s];
                    if (blk[s]) m |= uchar(1 << l);
                } else {
                    m |= uchar(1 << l);
                }
            }
            sh[l][yy][xx] = v;
        }
        sblk[yy][xx] = m;
    }
    if (li == 0) tflag = 0;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    int ly = int(lt.y) + 1, lx = int(lt.x) + 1;
    int y = y0 + ly, x = x0 + lx;
    bool mine = y < H && x < W;
    // allowed moves per layer (bounds, target free, diagonal corners free) and via targets
    uchar own = sblk[ly][lx];
    uchar mv[NLMAX];
    bool via = mine && vok[ulong(y) * ulong(W) + ulong(x)] != 0;
    for (int l = 0; l < NLMAX; l++) {
        uchar a = 0;
        if (l < NL && mine && !(own >> l & 1)) {
            for (int d = 0; d < 8; d++) {
                int yy = ly + DI[d], xx = lx + DJ[d];
                if (sblk[yy][xx] >> l & 1) continue;
                if (d >= 4 && ((sblk[ly][xx] >> l & 1) || (sblk[yy][lx] >> l & 1))) continue;
                a |= uchar(1 << d);
            }
        }
        mv[l] = a;
    }
    bool any = false;
    for (uint it = 0; it < p.maxit; it++) {
        float nv[NLMAX];
        bool ch = false;
        for (int l = 0; l < NL; l++) {
            float v = sh[l][ly][lx];
            nv[l] = v;
            if (!mine || (own >> l & 1)) continue;
            uchar a = mv[l];
            for (int d = 0; d < 8; d++)
                if (a >> d & 1) v = min(v, mcost[l * 8 + d] + sh[l][ly + DI[d]][lx + DJ[d]]);
            if (via)
                for (int l2 = 0; l2 < NL; l2++)
                    if (l2 != l && (p.lay_ok >> l2 & 1) && !(own >> l2 & 1)) v = min(v, p.vcost + sh[l2][ly][lx]);
            if (v < nv[l]) { nv[l] = v; ch = true; }
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
        if (ch) {
            for (int l = 0; l < NL; l++) sh[l][ly][lx] = nv[l];
            tflag = 1;
            any = true;
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
        int f = tflag;
        threadgroup_barrier(mem_flags::mem_threadgroup);
        if (f == 0) break;
        if (li == 0) tflag = 0;
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    threadgroup int tany;
    if (li == 0) tany = 0;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (any) {
        tany = 1;
        for (int l = 0; l < NL; l++)
            F[ulong(l) * HW + ulong(y) * ulong(W) + ulong(x)] = sh[l][ly][lx];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (li == 0 && tany) {
        last_change[ty * p.tiles_x + tx] = p.round + 1;
        atomic_fetch_add_explicit(&changed[p.round], 1u, memory_order_relaxed);
    }
}
