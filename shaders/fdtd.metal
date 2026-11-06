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
    uint ng;            // gathered samples for dft_gather
    uint pl0, pl1, pl2; // CPML slab lengths per axis (fused kernels)
    uint kd0, kd1, npole; // Debye planes and pole count (fused E kernel)
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

// E update with per-edge material indices (ushort) into a table of (ca, cb): 6 bytes per cell instead of
// 24 for the coefficients; identical arithmetic to update_e
kernel void update_e_idx(device float* ex [[buffer(0)]], device float* ey [[buffer(1)]], device float* ez [[buffer(2)]],
                         device const float* hx [[buffer(3)]], device const float* hy [[buffer(4)]], device const float* hz [[buffer(5)]],
                         device const float* inv_dd [[buffer(6)]], constant P& p [[buffer(7)]],
                         device const ushort* mx [[buffer(8)]], device const ushort* my [[buffer(9)]], device const ushort* mz [[buffer(10)]],
                         device const float2* tab [[buffer(11)]],
                         uint3 g [[thread_position_in_grid]]) {
    uint i = g.x, j = g.y, k = g.z;
    if (i >= p.nx || j >= p.ny || k >= p.nz) return;
    uint sy = p.nx, sz = p.nx * p.ny;
    uint id = (k * p.ny + j) * p.nx + i;
    device const float* idx = inv_dd;
    device const float* idy = inv_dd + p.nx;
    device const float* idz = inv_dd + p.nx + p.ny;
    float2 c;
    c = tab[mx[id]];
    if (c.y != 0.0f && j > 0 && k > 0) ex[id] = c.x * ex[id] + c.y * ((hz[id] - hz[id - sy]) * idy[j] - (hy[id] - hy[id - sz]) * idz[k]);
    c = tab[my[id]];
    if (c.y != 0.0f && i > 0 && k > 0) ey[id] = c.x * ey[id] + c.y * ((hx[id] - hx[id - sz]) * idz[k] - (hz[id] - hz[id - 1]) * idx[i]);
    c = tab[mz[id]];
    if (c.y != 0.0f && i > 0 && j > 0) ez[id] = c.x * ez[id] + c.y * ((hy[id] - hy[id - 1]) * idx[i] - (hx[id] - hx[id - sy]) * idy[j]);
}

// Fused H update: main curl update plus the CPML corrections of every axis whose slab contains the
// cell, accumulated in registers and stored once (the separate pml_h passes re-read and re-write H).
// Per axis a (b = a+1, c = a+2): psi0 <- dE_c/da, H_b += kmu psi0;  psi1 <- dE_b/da, H_c -= kmu psi1.
kernel void update_h_f(device float* hx [[buffer(0)]], device float* hy [[buffer(1)]], device float* hz [[buffer(2)]],
                       device const float* ex [[buffer(3)]], device const float* ey [[buffer(4)]], device const float* ez [[buffer(5)]],
                       device const float* inv_d [[buffer(6)]], constant P& p [[buffer(7)]],
                       device const uint* map0 [[buffer(8)]], device const float* bh0 [[buffer(9)]], device const float* ch0 [[buffer(10)]],
                       device float* ps00 [[buffer(11)]], device float* ps01 [[buffer(12)]],
                       device const uint* map1 [[buffer(13)]], device const float* bh1 [[buffer(14)]], device const float* ch1 [[buffer(15)]],
                       device float* ps10 [[buffer(16)]], device float* ps11 [[buffer(17)]],
                       device const uint* map2 [[buffer(18)]], device const float* bh2 [[buffer(19)]], device const float* ch2 [[buffer(20)]],
                       device float* ps20 [[buffer(21)]], device float* ps21 [[buffer(22)]],
                       uint3 g [[thread_position_in_grid]]) {
    uint i = g.x, j = g.y, k = g.z;
    if (i >= p.nx || j >= p.ny || k >= p.nz) return;
    uint sy = p.nx, sz = p.nx * p.ny;
    uint id = (k * p.ny + j) * p.nx + i;
    device const float* idx = inv_d;
    device const float* idy = inv_d + p.nx;
    device const float* idz = inv_d + p.nx + p.ny;
    bool in_ = i + 1 < p.nx, jn = j + 1 < p.ny, kn = k + 1 < p.nz;
    float vx = hx[id], vy = hy[id], vz = hz[id];
    float ex0 = ex[id], ey0 = ey[id], ez0 = ez[id];
    float exy = jn ? ex[id + sy] : 0.0f, exz = kn ? ex[id + sz] : 0.0f;
    float eyx = in_ ? ey[id + 1] : 0.0f, eyz = kn ? ey[id + sz] : 0.0f;
    float ezx = in_ ? ez[id + 1] : 0.0f, ezy = jn ? ez[id + sy] : 0.0f;
    if (jn && kn) vx -= p.kmu * ((ezy - ez0) * idy[j] - (eyz - ey0) * idz[k]);
    if (in_ && kn) vy -= p.kmu * ((exz - ex0) * idz[k] - (ezx - ez0) * idx[i]);
    if (in_ && jn) vz -= p.kmu * ((eyx - ey0) * idx[i] - (exy - ex0) * idy[j]);
    uint s = map0[i];
    if (s != 0xffffffffu && in_) {
        float B = bh0[s], C = ch0[s], d = idx[i];
        uint q = (k * p.ny + j) * p.pl0 + s;
        if (kn) { float v = B * ps00[q] + C * (ezx - ez0) * d; ps00[q] = v; vy += p.kmu * v; }
        if (jn) { float v = B * ps01[q] + C * (eyx - ey0) * d; ps01[q] = v; vz -= p.kmu * v; }
    }
    s = map1[j];
    if (s != 0xffffffffu && jn) {
        float B = bh1[s], C = ch1[s], d = idy[j];
        uint q = (k * p.pl1 + s) * p.nx + i;
        if (in_) { float v = B * ps10[q] + C * (exy - ex0) * d; ps10[q] = v; vz += p.kmu * v; }
        if (kn) { float v = B * ps11[q] + C * (ezy - ez0) * d; ps11[q] = v; vx -= p.kmu * v; }
    }
    s = map2[k];
    if (s != 0xffffffffu && kn) {
        float B = bh2[s], C = ch2[s], d = idz[k];
        uint q = (s * p.ny + j) * p.nx + i;
        if (jn) { float v = B * ps20[q] + C * (eyz - ey0) * d; ps20[q] = v; vx += p.kmu * v; }
        if (in_) { float v = B * ps21[q] + C * (exz - ex0) * d; ps21[q] = v; vy -= p.kmu * v; }
    }
    hx[id] = vx; hy[id] = vy; hz[id] = vz;
}

// Fused E update: material-index curl update plus the CPML corrections of every axis whose slab contains
// the cell. Per axis a: psi0 <- dH_c/da, E_b -= cb_b psi0;  psi1 <- dH_b/da, E_c += cb_c psi1.
kernel void update_e_f(device float* ex [[buffer(0)]], device float* ey [[buffer(1)]], device float* ez [[buffer(2)]],
                       device const float* hx [[buffer(3)]], device const float* hy [[buffer(4)]], device const float* hz [[buffer(5)]],
                       device const float* inv_dd [[buffer(6)]], constant P& p [[buffer(7)]],
                       device const ushort* mx [[buffer(8)]], device const ushort* my [[buffer(9)]], device const ushort* mz [[buffer(10)]],
                       device const float2* tab [[buffer(11)]],
                       device const uint* map0 [[buffer(12)]], device const float* be0 [[buffer(13)]], device const float* ce0 [[buffer(14)]],
                       device float* ps00 [[buffer(15)]], device float* ps01 [[buffer(16)]],
                       device const uint* map1 [[buffer(17)]], device const float* be1 [[buffer(18)]], device const float* ce1 [[buffer(19)]],
                       device float* ps10 [[buffer(20)]], device float* ps11 [[buffer(21)]],
                       device const uint* map2 [[buffer(22)]], device const float* be2 [[buffer(23)]], device const float* ce2 [[buffer(24)]],
                       device float* ps20 [[buffer(25)]], device float* ps21 [[buffer(26)]],
                       device float* dj [[buffer(27)]], device const float* db [[buffer(28)]], device const float2* kh [[buffer(29)]],
                       uint3 g [[thread_position_in_grid]]) {
    uint i = g.x, j = g.y, k = g.z;
    if (i >= p.nx || j >= p.ny || k >= p.nz) return;
    uint sy = p.nx, sz = p.nx * p.ny;
    uint id = (k * p.ny + j) * p.nx + i;
    device const float* idx = inv_dd;
    device const float* idy = inv_dd + p.nx;
    device const float* idz = inv_dd + p.nx + p.ny;
    float2 cx = tab[mx[id]], cy = tab[my[id]], cz = tab[mz[id]];
    if (cx.y == 0.0f && cy.y == 0.0f && cz.y == 0.0f) return; // metal or wall on all three: nothing changes
    float hx0 = hx[id], hy0 = hy[id], hz0 = hz[id];
    float hxy = j > 0 ? hx[id - sy] : 0.0f, hxz = k > 0 ? hx[id - sz] : 0.0f;
    float hyx = i > 0 ? hy[id - 1] : 0.0f, hyz = k > 0 ? hy[id - sz] : 0.0f;
    float hzx = i > 0 ? hz[id - 1] : 0.0f, hzy = j > 0 ? hz[id - sy] : 0.0f;
    float vx = ex[id], vy = ey[id], vz = ez[id];
    float e0x = vx, e0y = vy, e0z = vz;
    if (cx.y != 0.0f && j > 0 && k > 0) vx = cx.x * vx + cx.y * ((hz0 - hzy) * idy[j] - (hy0 - hyz) * idz[k]);
    if (cy.y != 0.0f && i > 0 && k > 0) vy = cy.x * vy + cy.y * ((hx0 - hxz) * idz[k] - (hz0 - hzx) * idx[i]);
    if (cz.y != 0.0f && i > 0 && j > 0) vz = cz.x * vz + cz.y * ((hy0 - hyx) * idx[i] - (hx0 - hxy) * idy[j]);
    uint s = map0[i];
    if (s != 0xffffffffu && i > 0) {
        float B = be0[s], C = ce0[s], d = idx[i];
        uint q = (k * p.ny + j) * p.pl0 + s;
        if (cy.y != 0.0f) { float v = B * ps00[q] + C * (hz0 - hzx) * d; ps00[q] = v; vy -= cy.y * v; }
        if (cz.y != 0.0f) { float v = B * ps01[q] + C * (hy0 - hyx) * d; ps01[q] = v; vz += cz.y * v; }
    }
    s = map1[j];
    if (s != 0xffffffffu && j > 0) {
        float B = be1[s], C = ce1[s], d = idy[j];
        uint q = (k * p.pl1 + s) * p.nx + i;
        if (cz.y != 0.0f) { float v = B * ps10[q] + C * (hx0 - hxy) * d; ps10[q] = v; vz -= cz.y * v; }
        if (cx.y != 0.0f) { float v = B * ps11[q] + C * (hz0 - hzy) * d; ps11[q] = v; vx += cx.y * v; }
    }
    s = map2[k];
    if (s != 0xffffffffu && k > 0) {
        float B = be2[s], C = ce2[s], d = idz[k];
        uint q = (s * p.ny + j) * p.nx + i;
        if (cx.y != 0.0f) { float v = B * ps20[q] + C * (hy0 - hyz) * d; ps20[q] = v; vx -= cx.y * v; }
        if (cy.y != 0.0f) { float v = B * ps21[q] + C * (hx0 - hxz) * d; ps21[q] = v; vy += cy.y * v; }
    }
    // Debye polarisation currents (disp.rs): E -= cb sum hk_p J_p, then J_p = k_p J_p + b_p (E_new - E_old)
    if (p.npole > 0 && k >= p.kd0 && k <= p.kd1) {
        uint np = p.npole;
        uint nd = p.nx * p.ny * (p.kd1 - p.kd0 + 1);
        uint qd = ((k - p.kd0) * p.ny + j) * p.nx + i;
        float cbs[3] = {cx.y, cy.y, cz.y};
        float vv[3] = {vx, vy, vz};
        float e0[3] = {e0x, e0y, e0z};
        for (uint c = 0; c < 3; c++) {
            if (cbs[c] == 0.0f) continue;
            uint base = (c * nd + qd) * np;
            float s = 0.0f;
            bool any = false;
            for (uint pp = 0; pp < np; pp++) { s += kh[pp].y * dj[base + pp]; any = any || db[base + pp] != 0.0f; }
            if (!any) continue;
            float e1 = vv[c] - cbs[c] * s;
            float de = e1 - e0[c];
            for (uint pp = 0; pp < np; pp++) dj[base + pp] = kh[pp].x * dj[base + pp] + db[base + pp] * de;
            vv[c] = e1;
        }
        vx = vv[0]; vy = vv[1]; vz = vv[2];
    }
    ex[id] = vx; ey[id] = vy; ez[id] = vz;
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

// Copper with a surface impedance (cond.rs): one thread per sheet edge; E = Z_s J_s with J_s the curl of H
// across the sheet times the dual length, Z_s a Foster chain R_0 + sum R_k s tau_k / (1 + s tau_k) whose
// section currents decay exactly over a step. p.axis = sections per edge, p.slab_len = edge count.
kernel void cond_e(device float* ex [[buffer(0)]], device float* ey [[buffer(1)]],
                   device const float* hx [[buffer(3)]], device const float* hy [[buffer(4)]], device const float* hz [[buffer(5)]],
                   device const float* inv_dd [[buffer(6)]], constant P& p [[buffer(7)]],
                   device const uint* ids [[buffer(8)]], device const float* fac [[buffer(9)]], device const uint* lay [[buffer(10)]],
                   device float* ci [[buffer(11)]], device const float* tab [[buffer(12)]], device const float* dec [[buffer(13)]],
                   uint t [[thread_position_in_grid]]) {
    if (t >= p.slab_len) return;
    uint w = ids[t];
    uint a = w >> 30, id = w & 0x3fffffffu;
    uint sy = p.nx, sz = p.nx * p.ny;
    uint i = id % p.nx, j = (id / p.nx) % p.ny, k = id / sz;
    device const float* idx = inv_dd;
    device const float* idy = inv_dd + p.nx;
    device const float* idz = inv_dd + p.nx + p.ny;
    float curl = (a == 0) ? (hz[id] - hz[id - sy]) * idy[j] - (hy[id] - hy[id - sz]) * idz[k]
                          : (hx[id] - hx[id - sz]) * idz[k] - (hz[id] - hz[id - 1]) * idx[i];
    float js = curl * fac[t];
    uint nb = p.axis;
    device const float* r = tab + lay[t] * (nb + 1);
    device float* ik = ci + t * nb;
    float e = r[0] * js;
    for (uint q = 0; q < nb; q++) {
        float v = js + (ik[q] - js) * dec[q];
        ik[q] = v;
        e += r[q + 1] * (js - v);
    }
    if (a == 0) ex[id] = e; else ey[id] = e;
}

// The fused E kernel has already updated the Debye currents from E_new - E_old; the source's share of the
// E change is added here to E and, on dispersive edges, to the currents (J_p += b_p dE), which equals
// injecting the source before the current update (as the CPU reference does).
kernel void add_src(device float* ex [[buffer(0)]], device float* ey [[buffer(1)]], device float* ez [[buffer(2)]],
                    device const uint* sc [[buffer(3)]], device const uint* si [[buffer(4)]], device const float* sg [[buffer(5)]],
                    constant P& p [[buffer(7)]], device float* dj [[buffer(8)]], device const float* db [[buffer(9)]],
                    uint t [[thread_position_in_grid]]) {
    if (t >= p.nsrc) return;
    device float* e[3] = {ex, ey, ez};
    uint c = sc[t], id = si[t];
    float de = sg[t] * p.vs;
    e[c][id] += de;
    if (p.npole > 0) {
        uint plane = p.nx * p.ny;
        uint k = id / plane;
        if (k >= p.kd0 && k <= p.kd1) {
            uint nd = plane * (p.kd1 - p.kd0 + 1);
            uint base = (c * nd + (id - p.kd0 * plane)) * p.npole;
            for (uint pp = 0; pp < p.npole; pp++) dj[base + pp] += db[base + pp] * de;
        }
    }
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

// running DFT of gathered samples (indices idx into one component array)
kernel void dft_gather(device const float* v [[buffer(0)]], device float* acc [[buffer(1)]], device const uint* idx [[buffer(2)]],
                       constant P& p [[buffer(7)]], uint t [[thread_position_in_grid]]) {
    if (t >= p.ng) return;
    float x = v[idx[t]];
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
