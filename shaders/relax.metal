// Elastic relaxation of routed copper (see src/relax.rs, which has the same arithmetic on the CPU): one thread
// per node; a tension step, then Jacobi projections of the constraints of the node's elements.
#include <metal_stdlib>
using namespace metal;

struct RelaxParams {
    uint iters, sub;
    float alpha_wire, alpha_via, beta, dmax, tol;
    uint n;
};

static inline float2 unit2(float2 a) {
    float l = length(a);
    return l > 1e-9f ? a / l : float2(0.0f);
}

// closest points of segments p1-q1 and p2-q2 (either may be a point): parameters on each
static inline float2 seg_seg(float2 p1, float2 q1, float2 p2, float2 q2) {
    float2 d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
    float a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
    const float EPS = 1e-12f;
    if (a <= EPS && e <= EPS) return float2(0.0f);
    if (a <= EPS) return float2(0.0f, clamp(f / e, 0.0f, 1.0f));
    float c = dot(d1, r);
    if (e <= EPS) return float2(clamp(-c / a, 0.0f, 1.0f), 0.0f);
    float b = dot(d1, d2);
    float den = a * e - b * b;
    float s = den > EPS ? clamp((b * f - c * e) / den, 0.0f, 1.0f) : 0.0f;
    float t = (b * s + f) / e;
    if (t < 0.0f) { t = 0.0f; s = clamp(-c / a, 0.0f, 1.0f); }
    else if (t > 1.0f) { t = 1.0f; s = clamp((b - c) / a, 0.0f, 1.0f); }
    return float2(s, t);
}

kernel void relax_tension(device const float2 *p [[buffer(0)]],
                          device float2 *out [[buffer(1)]],
                          device const uint *kind [[buffer(2)]],
                          device const uint *tn_off [[buffer(3)]],
                          device const uint *tn [[buffer(4)]],
                          constant RelaxParams &prm [[buffer(5)]],
                          device float *scale [[buffer(6)]],
                          device const float *viol [[buffer(7)]],
                          uint i [[thread_position_in_grid]])
{
    if (i >= prm.n) return;
    float2 x = p[i];
    float sc = viol[i] > prm.tol ? scale[i] * 0.5f : min(1.0f, scale[i] * 1.25f);
    scale[i] = sc;
    uint a = tn_off[i], b = tn_off[i + 1];
    float2 d;
    if (kind[i] == 1 && b - a == 2) {
        float2 u = p[tn[a]], w = p[tn[a + 1]];
        float2 g = unit2(u - x) + unit2(w - x);
        float2 t = unit2(w - u);
        float2 lap = 0.5f * (u + w) - x;
        d = g * prm.alpha_wire + t * (prm.beta * dot(lap, t));
    } else if (kind[i] == 2) {
        float2 g = float2(0.0f);
        for (uint k = a; k < b; k++) g += unit2(p[tn[k]] - x);
        d = g * prm.alpha_via;
    } else {
        out[i] = x;
        return;
    }
    d *= sc;
    float l = length(d);
    out[i] = x + (l > prm.dmax ? d * (prm.dmax / l) : d);
}

struct Model {
    device const uint *kind;
    device const uint *el_off;
    device const uint *el_node;
};

static inline bool is_fixed(Model m, uint f) {
    for (uint k = m.el_off[f]; k < m.el_off[f + 1]; k++)
        if (m.kind[m.el_node[k]] != 0) return false;
    return true;
}

// correction for node i (at xi, the rest at p) from the constraint keeping element e `need` from f, and the gap
// deficit; false when the gap holds
static inline bool contact(Model m, device const float2 *p, uint i, float2 xi, uint e, uint f, float need,
                           thread float2 &corr, thread float &pen)
{
    uint a0 = m.el_node[m.el_off[e]], a1 = m.el_node[m.el_off[e + 1] - 1];
    float2 pa = a0 == i ? xi : p[a0], qa = a1 == i ? xi : p[a1];
    float ma0 = m.kind[a0] != 0 ? 1.0f : 0.0f, ma1 = m.kind[a1] != 0 ? 1.0f : 0.0f;
    uint f0 = m.el_off[f], f1 = m.el_off[f + 1];
    uint nf = f1 - f0;
    uint edges = nf <= 2 ? 1 : nf;
    float bd2 = INFINITY, bs = 0.0f, bt = 0.0f;
    float2 bca = float2(0.0f), bcb = float2(0.0f);
    uint bb0 = 0, bb1 = 0;
    for (uint j = 0; j < edges; j++) {
        uint b0 = nf <= 2 ? m.el_node[f0] : m.el_node[f0 + j];
        uint b1 = nf <= 2 ? m.el_node[f1 - 1] : m.el_node[f0 + (j + 1) % nf];
        float2 pb = b0 == i ? xi : p[b0], qb = b1 == i ? xi : p[b1];
        float2 st = seg_seg(pa, qa, pb, qb);
        float2 ca = pa + (qa - pa) * st.x;
        float2 cb = pb + (qb - pb) * st.y;
        float2 dv = ca - cb;
        float d2 = dot(dv, dv);
        if (d2 < bd2) { bd2 = d2; bs = st.x; bt = st.y; bca = ca; bcb = cb; bb0 = b0; bb1 = b1; }
    }
    if (bd2 >= need * need) return false;
    float d = sqrt(bd2);
    if (d < 1e-7f) return false;
    float2 nrm = (bca - bcb) / d;
    float la = a0 == a1 ? 1.0f : (i == a0 ? 1.0f - bs : bs);
    float w = a0 == a1 ? ma0 : ma0 * (1.0f - bs) * (1.0f - bs) + ma1 * bs * bs;
    if (nf <= 2) {
        float mb0 = m.kind[bb0] != 0 ? 1.0f : 0.0f, mb1 = m.kind[bb1] != 0 ? 1.0f : 0.0f;
        w += bb0 == bb1 ? mb0 : mb0 * (1.0f - bt) * (1.0f - bt) + mb1 * bt * bt;
    }
    if (w < 1e-9f) return false;
    // the exact share for a contact near a pinned end is a long lever (penetration / s): capped, so a node
    // never leaps
    corr = nrm * ((need - d) * min(la / w, 2.0f));
    pen = need - d;
    return true;
}

// moving neighbours share a Jacobi average; fixed obstacles then apply in full from where that leaves the node,
// so no number of pushing neighbours can carry a node through a boundary
kernel void relax_project(device const float2 *p [[buffer(0)]],
                          device float2 *out [[buffer(1)]],
                          device const uint *kind [[buffer(2)]],
                          device const float4 *leash [[buffer(3)]],
                          device const uint *el_off [[buffer(4)]],
                          device const uint *el_node [[buffer(5)]],
                          device const uint *ne_off [[buffer(6)]],
                          device const uint *ne [[buffer(7)]],
                          device const uint *c_off [[buffer(8)]],
                          device const uint *c_el [[buffer(9)]],
                          device const float *c_need [[buffer(10)]],
                          constant RelaxParams &prm [[buffer(11)]],
                          device float *viol [[buffer(12)]],
                          uint i [[thread_position_in_grid]])
{
    if (i >= prm.n) return;
    float2 x = p[i];
    if (kind[i] == 0) { out[i] = x; viol[i] = 0.0f; return; }
    Model m = { kind, el_off, el_node };
    float2 acc = float2(0.0f);
    uint n = 0;
    float worst = 0.0f;
    float2 corr;
    float pen;
    for (uint k = ne_off[i]; k < ne_off[i + 1]; k++) {
        uint e = ne[k];
        for (uint c = c_off[e]; c < c_off[e + 1]; c++) {
            uint f = c_el[c];
            if (is_fixed(m, f)) continue;
            if (contact(m, p, i, x, e, f, c_need[c], corr, pen)) { acc += corr; n++; worst = max(worst, pen); }
        }
    }
    float2 y = n > 0 ? x + acc / float(n) : x;
    float2 fix = float2(0.0f);
    for (uint k = ne_off[i]; k < ne_off[i + 1]; k++) {
        uint e = ne[k];
        for (uint c = c_off[e]; c < c_off[e + 1]; c++) {
            uint f = c_el[c];
            if (!is_fixed(m, f)) continue;
            if (contact(m, p, i, y, e, f, c_need[c], corr, pen)) { fix += corr; worst = max(worst, pen); }
        }
    }
    y += fix;
    float4 lz = leash[i];
    if (lz.z >= 0.0f) {
        float2 o = y - lz.xy;
        float l = length(o);
        if (l > lz.z) y = lz.xy + o * (lz.z / l);
    }
    out[i] = y;
    viol[i] = worst;
}
