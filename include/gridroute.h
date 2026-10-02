/* gridroute: kernels for grid-based PCB routers (Rust, multi-threaded CPU + Metal GPU), C ABI.
 *
 * All grids are dense row-major arrays owned by the caller; the library keeps nothing between calls except GPU
 * pipelines and per-thread search buffers. Build: cargo build --release (target/release/libgridroute.{dylib,so}).
 */
#ifndef GRIDROUTE_H
#define GRIDROUTE_H
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum { GR_AUTO = 0, GR_CPU = 1, GR_METAL = 2, GR_NAIVE = 3 };

uint32_t gr_version(void);
int gr_metal_available(void);
size_t gr_device_name(char *buf, size_t cap);
/* GR_AUTO dilations run on the GPU from this work size (h * w * element rows * layers); negative = never. */
void gr_set_metal_threshold(int64_t work);
/* CPU worker threads; call before the first kernel (0 = one per core). Returns 1 on success. */
int gr_set_threads(int n);

/* Obstacle dilation: out = 1 where a source obstacle lies within the structuring element of the cell.
 * src [nl][ny][nx] int16 (0 = free); obstacle = value != 0 and not in excl[nexcl], or extra[ny][nx] != 0 on layers
 * with extra_on[l] set (extra may be NULL). Window rows i0..i0+h-1, cols j0..j0+w-1. span[nspan], nspan = 2 rc + 1:
 * half-width (columns) of the element on row offset di = k - rc, -1 = row unused. reduce_or != 0: one output plane
 * (OR over the layers), else one per layer. out: planes * h * w bytes. Returns the backend used, < 0 on error. */
int gr_dilate(const int16_t *src, int nl, int ny, int nx, const int32_t *layers, int nlayers, const int16_t *excl,
              int nexcl, const uint8_t *extra, const uint8_t *extra_on, int i0, int j0, int h, int w,
              const int32_t *span, int nspan, int reduce_or, int backend, uint8_t *out);

/* Copper island of net nid connected to the seeds: 4-neighbour copper in core [nl][ny][nx], layer jumps where
 * thru [ny][nx] == nid. seeds: nseeds x (layer, y, x) window-relative. out: nl * h * w bytes. Returns cells. */
int64_t gr_island(const int16_t *core, int nl, int ny, int nx, const int16_t *thru, int16_t nid, int i0, int j0,
                  int h, int w, const int32_t *seeds, int nseeds, uint8_t *out);

/* Multi-layer A*: states L*h*w + i*w + j. blk, tgt [nl][h][w]; vok [h][w] (via allowed); mcost [nl][8] cost per
 * move direction (E, W, S, N, SE, SW, NE, NW); lay_ok [nl]; vcost per via; turn per direction change; hmul times the
 * octile distance to the target box (tx0, ty0)-(tx1, ty1) is the heuristic. Writes the path (source first) to out.
 * Returns its length, -1 if none was found within max_exp expansions, -2 if out_cap is too small.
 * nexp (may be NULL) receives the number of expansions. */
int gr_astar(int nl, int h, int w, const uint8_t *blk, const uint8_t *vok, const uint8_t *tgt, const int32_t *src,
             int nsrc, const float *mcost, const uint8_t *lay_ok, float vcost, float turn, float hmul, int tx0,
             int ty0, int tx1, int ty1, int64_t max_exp, int32_t *out, int out_cap, int64_t *nexp);

/* gr_astar with the unoptimised reference search (same results; for tests and benchmarks). */
int gr_astar_ref(int nl, int h, int w, const uint8_t *blk, const uint8_t *vok, const uint8_t *tgt, const int32_t *src,
                 int nsrc, const float *mcost, const uint8_t *lay_ok, float vcost, float turn, float hmul, int tx0,
                 int ty0, int tx1, int ty1, int64_t max_exp, int32_t *out, int out_cap, int64_t *nexp);

/* Cost-to-target field: out[nl*h*w] = cheapest cost from each state to a target over the search graph without turn
 * penalties (+inf if unreachable). GR_CPU: exact Dijkstra; GR_METAL / GR_AUTO: tiled GPU relaxation (same values).
 * Returns the backend used; rounds (may be NULL) receives the GPU round count. */
int gr_field(int nl, int h, int w, const uint8_t *blk, const uint8_t *vok, const uint8_t *tgt, const float *mcost,
             const uint8_t *lay_ok, float vcost, int backend, float *out, int64_t *rounds);

/* Hybrid search: gr_astar up to `budget` expansions; a search still open then switches to the field (GPU) and
 * field-guided A* (f = g + weight * field) within max_exp. Easy searches return exactly gr_astar's path. info
 * (may be NULL): [phase 1 octile / 2 guided / 3 no path, field seconds, expansions]. touched (may be NULL):
 * [ceil(h/32)][ceil(w/32)] bytes set for the 32x32 tiles the search read (its read set, for memoisation). */
int gr_astar_hybrid(int nl, int h, int w, const uint8_t *blk, const uint8_t *vok, const uint8_t *tgt,
                    const int32_t *src, int nsrc, const float *mcost, const uint8_t *lay_ok, float vcost, float turn,
                    float hmul, int tx0, int ty0, int tx1, int ty1, int64_t max_exp, int32_t *out, int out_cap,
                    int64_t budget, float weight, double *info, uint8_t *touched);

/* Field-guided A* only (experimental; computes the field itself). */
int gr_astar_guided(int nl, int h, int w, const uint8_t *blk, const uint8_t *vok, const uint8_t *tgt,
                    const int32_t *src, int nsrc, const float *mcost, const uint8_t *lay_ok, float vcost, float turn,
                    float weight, int64_t max_exp, int32_t *out, int out_cap, int64_t *nexp, double *field_secs);

/* 128-bit content hashes (out[2]) for memoising searches: a window of several [nl][ny][nx] arrays (dims: nl, ny, nx,
 * element size per array); whole buffers (parallel); listed 32x32 tiles ((ty, tx) pairs) of [nl][h][w] byte arrays. */
void gr_hash_windows(int narr, const uint8_t *const *arrays, const int64_t *dims, int64_t i0, int64_t j0, int64_t i1,
                     int64_t j1, uint64_t seed, uint64_t *out);
void gr_hash_buffers(int n, const uint8_t *const *ptrs, const int64_t *lens, uint64_t seed, uint64_t *out);
void gr_hash_tiles(int narr, const uint8_t *const *arrays, const int32_t *nls, int h, int w, int ntiles,
                   const int32_t *tiles, uint64_t *out);

#ifdef __cplusplus
}
#endif
#endif
