// fieldsolve FDTD kernels: the same updates as src/bin/fieldsolve/fdtd.rs (CPU reference), one thread per
// cell. Index (k * ny + j) * nx + i for every component.
#include <metal_stdlib>
using namespace metal;

struct P {
    uint nx, ny, nz;
    float kmu;          // dt / mu0
    uint axis;          // CPML axis for pml kernels
    uint slab_len;      // CPML slab length for that axis
    uint step;          // probe row / excitation step
    float vs;           // source voltage this step
    uint nsrc, nprobe;
    float w[12];        // DFT weights (re, im) per component
};

kernel void update_h(device float* hx [[buffer(0)]], device float* hy [[buffer(1)]], device float* hz [[buffer(2)]],
                     device const float* ex [[buffer(3)]], device const float* ey [[buffer(4)]], device const float* ez [[buffer(5)]],
                     device const float* inv_d [[buffer(6)]], constant P& p [[buffer(7)]],
                     uint3 g [[thread_position_in_grid]]) {
    uint i = g.x, j = g.y, k = g.z;
    if (i >= p.nx || j >= p.ny || k >= p.nz) return;
    uint sy = p.nx, sz = p.nx * p.ny;
    uint id = (k * p.ny + j) * p.nx + i;
    device const float* idx = inv_d;
    device const float* idy = inv_d + p.nx;
    device const float* idz = inv_d + p.nx + p.ny;
    bool in_ = i + 1 < p.nx, jn = j + 1 < p.ny, kn = k + 1 < p.nz;
    if (jn && kn) hx[id] -= p.kmu * ((ez[id + sy] - ez[id]) * idy[j] - (ey[id + sz] - ey[id]) * idz[k]);
    if (in_ && kn) hy[id] -= p.kmu * ((ex[id + sz] - ex[id]) * idz[k] - (ez[id + 1] - ez[id]) * idx[i]);
    if (in_ && jn) hz[id] -= p.kmu * ((ey[id + 1] - ey[id]) * idx[i] - (ex[id + sy] - ex[id]) * idy[j]);
}

kernel void update_e(device float* ex [[buffer(0)]], device float* ey [[buffer(1)]], device float* ez [[buffer(2)]],
                     device const float* hx [[buffer(3)]], device const float* hy [[buffer(4)]], device const float* hz [[buffer(5)]],
                     device const float* inv_dd [[buffer(6)]], constant P& p [[buffer(7)]],
                     device const float* cax [[buffer(8)]], device const float* cay [[buffer(9)]], device const float* caz [[buffer(10)]],
                     device const float* cbx [[buffer(11)]], device const float* cby [[buffer(12)]], device const float* cbz [[buffer(13)]],
                     uint3 g [[thread_position_in_grid]]) {
    uint i = g.x, j = g.y, k = g.z;
    if (i >= p.nx || j >= p.ny || k >= p.nz) return;
    uint sy = p.nx, sz = p.nx * p.ny;
    uint id = (k * p.ny + j) * p.nx + i;
    device const float* idx = inv_dd;
    device const float* idy = inv_dd + p.nx;
    device const float* idz = inv_dd + p.nx + p.ny;
    float b;
    b = cbx[id];
    if (b != 0.0f && j > 0 && k > 0) ex[id] = cax[id] * ex[id] + b * ((hz[id] - hz[id - sy]) * idy[j] - (hy[id] - hy[id - sz]) * idz[k]);
    b = cby[id];
    if (b != 0.0f && i > 0 && k > 0) ey[id] = cay[id] * ey[id] + b * ((hx[id] - hx[id - sz]) * idz[k] - (hz[id] - hz[id - 1]) * idx[i]);
    b = cbz[id];
    if (b != 0.0f && i > 0 && j > 0) ez[id] = caz[id] * ez[id] + b * ((hy[id] - hy[id - 1]) * idx[i] - (hx[id] - hx[id - sy]) * idy[j]);
}

// slab cell -> (field index, psi index, position along the axis); grid dims: axis 0 (L, ny, nz), 1 (nx, L, nz), 2 (nx, ny, L)
static inline bool slab_cell(constant P& p, device const uint* slab, uint3 g, thread uint& id, thread uint& q, thread uint& pos, thread uint3& c) {
    uint L = p.slab_len;
    if (p.axis == 0) {
        if (g.x >= L || g.y >= p.ny || g.z >= p.nz) return false;
        pos = slab[g.x];
        c = uint3(pos, g.y, g.z);
        q = (g.z * p.ny + g.y) * L + g.x;
    } else if (p.axis == 1) {
        if (g.x >= p.nx || g.y >= L || g.z >= p.nz) return false;
        pos = slab[g.y];
        c = uint3(g.x, pos, g.z);
        q = (g.z * L + g.y) * p.nx + g.x;
    } else {
        if (g.x >= p.nx || g.y >= p.ny || g.z >= L) return false;
        pos = slab[g.z];
        c = uint3(g.x, g.y, pos);
        q = (g.z * p.ny + g.y) * p.nx + g.x;
    }
    id = (c.z * p.ny + c.y) * p.nx + c.x;
    return true;
}

static inline uint comp(uint3 c, uint a) { return a == 0 ? c.x : (a == 1 ? c.y : c.z); }

// H_b += kmu * psi(dE_c/da) ; H_c -= kmu * psi(dE_b/da)   (b = a+1, c = a+2)
kernel void pml_h(device float* f0 [[buffer(0)]], device float* f1 [[buffer(1)]], device float* f2 [[buffer(2)]],
                  device const float* e0 [[buffer(3)]], device const float* e1 [[buffer(4)]], device const float* e2 [[buffer(5)]],
                  device const float* inv_d [[buffer(6)]], constant P& p [[buffer(7)]],
                  device const uint* slab [[buffer(8)]], device const float* bh [[buffer(9)]], device const float* ch [[buffer(10)]],
                  device float* psi0 [[buffer(11)]], device float* psi1 [[buffer(12)]],
                  uint3 g [[thread_position_in_grid]]) {
    uint id, q, pos; uint3 c;
    if (!slab_cell(p, slab, g, id, q, pos, c)) return;
    uint a = p.axis, b = (a + 1) % 3, cc = (a + 2) % 3;
    uint n[3] = {p.nx, p.ny, p.nz};
    if (pos + 1 >= n[a]) return;
    uint st[3] = {1, p.nx, p.nx * p.ny};
    device float* h[3] = {f0, f1, f2};
    device const float* e[3] = {e0, e1, e2};
    uint off = a == 0 ? 0 : (a == 1 ? p.nx : p.nx + p.ny);
    float invd = inv_d[off + pos];
    uint s = (a == 0) ? g.x : (a == 1 ? g.y : g.z);
    float B = bh[s], C = ch[s];
    if (comp(c, cc) + 1 < n[cc]) {
        float der = (e[cc][id + st[a]] - e[cc][id]) * invd;
        float v = B * psi0[q] + C * der;
        psi0[q] = v;
        h[b][id] += p.kmu * v;
    }
    if (comp(c, b) + 1 < n[b]) {
        float der = (e[b][id + st[a]] - e[b][id]) * invd;
        float v = B * psi1[q] + C * der;
        psi1[q] = v;
        h[cc][id] -= p.kmu * v;
    }
}

// E_b -= cb * psi(dH_c/da) ; E_c += cb * psi(dH_b/da)
kernel void pml_e(device float* f0 [[buffer(0)]], device float* f1 [[buffer(1)]], device float* f2 [[buffer(2)]],
                  device const float* h0 [[buffer(3)]], device const float* h1 [[buffer(4)]], device const float* h2 [[buffer(5)]],
                  device const float* inv_dd [[buffer(6)]], constant P& p [[buffer(7)]],
                  device const uint* slab [[buffer(8)]], device const float* be [[buffer(9)]], device const float* ce [[buffer(10)]],
                  device float* psi0 [[buffer(11)]], device float* psi1 [[buffer(12)]],
                  device const float* cb0 [[buffer(13)]], device const float* cb1 [[buffer(14)]], device const float* cb2 [[buffer(15)]],
                  uint3 g [[thread_position_in_grid]]) {
    uint id, q, pos; uint3 c;
    if (!slab_cell(p, slab, g, id, q, pos, c)) return;
    if (pos == 0) return;
    uint a = p.axis, b = (a + 1) % 3, cc = (a + 2) % 3;
    uint st[3] = {1, p.nx, p.nx * p.ny};
    device float* e[3] = {f0, f1, f2};
    device const float* h[3] = {h0, h1, h2};
    device const float* cb[3] = {cb0, cb1, cb2};
    uint off = a == 0 ? 0 : (a == 1 ? p.nx : p.nx + p.ny);
    float invdd = inv_dd[off + pos];
    uint s = (a == 0) ? g.x : (a == 1 ? g.y : g.z);
    float B = be[s], C = ce[s];
    float k1 = cb[b][id];
    if (k1 != 0.0f) {
        float der = (h[cc][id] - h[cc][id - st[a]]) * invdd;
        float v = B * psi0[q] + C * der;
        psi0[q] = v;
        e[b][id] -= k1 * v;
    }
    float k2 = cb[cc][id];
    if (k2 != 0.0f) {
        float der = (h[b][id] - h[b][id - st[a]]) * invdd;
        float v = B * psi1[q] + C * der;
        psi1[q] = v;
        e[cc][id] += k2 * v;
    }
}

kernel void add_src(device float* ex [[buffer(0)]], device float* ey [[buffer(1)]], device float* ez [[buffer(2)]],
                    device const uint* sc [[buffer(3)]], device const uint* si [[buffer(4)]], device const float* sg [[buffer(5)]],
                    constant P& p [[buffer(7)]], uint t [[thread_position_in_grid]]) {
    if (t >= p.nsrc) return;
    device float* e[3] = {ex, ey, ez};
    e[sc[t]][si[t]] += sg[t] * p.vs;
}

// one thread per probe; terms [start[t], start[t+1]); comp 0..2 E, 3..5 H
kernel void probe(device const float* ex [[buffer(0)]], device const float* ey [[buffer(1)]], device const float* ez [[buffer(2)]],
                  device const float* hx [[buffer(3)]], device const float* hy [[buffer(4)]], device const float* hz [[buffer(5)]],
                  device const uint* start [[buffer(6)]], constant P& p [[buffer(7)]],
                  device const uint* tc [[buffer(8)]], device const uint* ti [[buffer(9)]], device const float* tw [[buffer(10)]],
                  device float* out [[buffer(11)]], uint t [[thread_position_in_grid]]) {
    if (t >= p.nprobe) return;
    device const float* f[6] = {ex, ey, ez, hx, hy, hz};
    float s = 0.0f;
    for (uint k = start[t]; k < start[t + 1]; k++) s += tw[k] * f[tc[k]][ti[k]];
    out[p.step * p.nprobe + t] = s;
}

// running DFT of one component: acc[2 id] += v * w_re, acc[2 id + 1] += v * w_im  (p.axis = component)
kernel void dft(device const float* v [[buffer(0)]], device float* acc [[buffer(1)]], constant P& p [[buffer(7)]],
                uint t [[thread_position_in_grid]]) {
    uint n = p.nx * p.ny * p.nz;
    if (t >= n) return;
    float x = v[t];
    acc[2 * t] += x * p.w[2 * p.axis];
    acc[2 * t + 1] += x * p.w[2 * p.axis + 1];
}

// per-row sums of E^2 + eta0^2 H^2
kernel void energy(device const float* ex [[buffer(0)]], device const float* ey [[buffer(1)]], device const float* ez [[buffer(2)]],
                   device const float* hx [[buffer(3)]], device const float* hy [[buffer(4)]], device const float* hz [[buffer(5)]],
                   device float* out [[buffer(6)]], constant P& p [[buffer(7)]], uint r [[thread_position_in_grid]]) {
    if (r >= p.ny * p.nz) return;
    const float eta2 = 376.730313668f * 376.730313668f;
    float s = 0.0f;
    uint base = r * p.nx;
    for (uint i = 0; i < p.nx; i++) {
        uint id = base + i;
        s += ex[id] * ex[id] + ey[id] * ey[id] + ez[id] * ez[id] + eta2 * (hx[id] * hx[id] + hy[id] * hy[id] + hz[id] * hz[id]);
    }
    out[r] = s;
}
