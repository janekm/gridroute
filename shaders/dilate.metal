// Obstacle dilation on the GPU; the same contract as src/dilate.rs (see there).
//
// Pass 1 (row_distance): one thread per (output column, touched source row) scans +-wmax columns of every layer of
// the group for the nearest obstacle and stores the distance (255 = none within wmax).
// Pass 2 (combine): one thread per output cell ORs `hd[row + k] <= span[k]` over the element's rows.
#include <metal_stdlib>
using namespace metal;

struct Params {
    uint nx, ny;            // source grid
    uint i0, j0, h, w;      // output window
    int r0;                 // first touched source row (i0 - rc, may be negative)
    uint nrows;             // touched rows (h + 2 rc)
    uint wmax;
    uint nspan;
    uint nlayers;           // layers in this group
    uint nexcl;
    uint has_extra;
    uint base_layer;        // layer index of the first plane in the source buffer
    uint base_row;          // source row at the start of the source buffer (per plane)
    uint plane_rows;        // rows per plane in the source buffer
    uint src_off;           // element offset of the first plane in src
    uint extra_off;         // byte offset of the first row in extra
    uint hd_off, out_off;   // byte offsets of this group's planes in hd / out
    uint layer[16];
    uint extra_on[16];
    short excl[8];
};

inline bool obstacle(short v, constant Params &p) {
    if (v == 0) return false;
    for (uint e = 0; e < p.nexcl; e++)
        if (v == p.excl[e]) return false;
    return true;
}

kernel void row_distance(device const short *src [[buffer(0)]],
                         device const uchar *extra [[buffer(1)]],
                         constant Params &p [[buffer(2)]],
                         device uchar *hd [[buffer(3)]],
                         uint2 gid [[thread_position_in_grid]])
{
    uint x = gid.x, k = gid.y;
    if (x >= p.w || k >= p.nrows) return;
    int i = p.r0 + int(k);
    uchar best = 255;
    if (i >= 0 && i < int(p.ny)) {
        int jc = int(p.j0 + x);
        int lo = max(jc - int(p.wmax), 0), hi = min(jc + int(p.wmax), int(p.nx) - 1);
        uint prow = uint(i) - p.base_row;
        for (uint n = 0; n < p.nlayers; n++) {
            uint l = p.layer[n];
            device const short *row = src + p.src_off + (ulong(l - p.base_layer) * p.plane_rows + prow) * p.nx;
            device const uchar *erow = extra + p.extra_off + ulong(prow) * p.nx;
            bool ex = p.has_extra && p.extra_on[n];
            for (int j = lo; j <= hi; j++) {
                uint d = uint(abs(j - jc));
                if (d >= best) continue;
                if (obstacle(row[j], p) || (ex && erow[j] != 0)) best = uchar(d);
            }
        }
    }
    hd[p.hd_off + k * p.w + x] = best;
}

kernel void combine(device const uchar *hd [[buffer(0)]],
                    constant Params &p [[buffer(1)]],
                    constant int *span [[buffer(2)]],
                    device uchar *out [[buffer(3)]],
                    uint2 gid [[thread_position_in_grid]])
{
    uint x = gid.x, y = gid.y;
    if (x >= p.w || y >= p.h) return;
    uchar hit = 0;
    for (uint k = 0; k < p.nspan; k++) {
        int s = span[k];
        if (s >= 0 && int(hd[p.hd_off + (y + k) * p.w + x]) <= s) { hit = 1; break; }
    }
    out[p.out_off + y * p.w + x] = hit;
}
