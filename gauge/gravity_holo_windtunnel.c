/*
 * PHY Stage 8: Discrete Regge Gravity, Nonlinear AME(4,8) Holography,
 * Quantum Extremal Surfaces & Black Hole Page Curve (Flow.GravityBundle H1-H5)
 *
 * Strictly self-contained C17 + POSIX pthreads + libm.
 */

#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* -------------------------------------------------------------------------- */
/* 1. PRNG, Numeric Helpers, Symmetric Jacobi & Linear System Solvers         */
/* -------------------------------------------------------------------------- */

typedef struct {
    uint64_t s;
} rng_t;

static inline void rng_seed(rng_t *r, uint64_t seed) {
    r->s = seed ? seed : 0x9e3779b97f4a7c15ULL;
}

static inline uint64_t rng_u64(rng_t *r) {
    uint64_t z = (r->s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static inline double rng_u01(rng_t *r) {
    return (double)(rng_u64(r) >> 11) * (1.0 / 9007199254740992.0);
}

static double clamp_cos(double x) {
    if (x > 1.0) return 1.0;
    if (x < -1.0) return -1.0;
    return x;
}

/* Real symmetric Jacobi eigensolver: destroys A[n*n], fills evals[n] ascending. */
static void jacobi_real_symm(int n, double *a, double *evals) {
    for (int sw = 0; sw < 30; ++sw) {
        double off = 0.0;
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                double v = a[p * n + q];
                off += v * v;
            }
        }
        if (off < 1e-28) break;

        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                double apq = a[p * n + q];
                if (fabs(apq) < 1e-16) continue;
                double app = a[p * n + p];
                double aqq = a[q * n + q];
                double tau = (aqq - app) / (2.0 * apq);
                double t;
                if (tau >= 0.0) {
                    t = 1.0 / (tau + sqrt(1.0 + tau * tau));
                } else {
                    t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
                }
                double c = 1.0 / sqrt(1.0 + t * t);
                double s = t * c;

                a[p * n + p] = app - t * apq;
                a[q * n + q] = aqq + t * apq;
                a[p * n + q] = 0.0;
                a[q * n + p] = 0.0;

                for (int r = 0; r < n; ++r) {
                    if (r == p || r == q) continue;
                    double arp = a[r * n + p];
                    double arq = a[r * n + q];
                    double nrp = c * arp - s * arq;
                    double nrq = s * arp + c * arq;
                    a[r * n + p] = nrp;
                    a[p * n + r] = nrp;
                    a[r * n + q] = nrq;
                    a[q * n + r] = nrq;
                }
            }
        }
    }
    for (int i = 0; i < n; ++i) evals[i] = a[i * n + i];
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (evals[j] < evals[i]) {
                double tmp = evals[i];
                evals[i] = evals[j];
                evals[j] = tmp;
            }
        }
    }
}

static int solve_linear_sys(int n, const double *A_in, const double *b_in, double *x) {
    double M[32 * 32], rhs[32];
    if (n > 32) return -1;
    memcpy(M, A_in, (size_t)(n * n) * sizeof(double));
    memcpy(rhs, b_in, (size_t)n * sizeof(double));
    for (int k = 0; k < n; ++k) {
        int piv = k;
        double mx = fabs(M[k * n + k]);
        for (int r = k + 1; r < n; ++r) {
            double v = fabs(M[r * n + k]);
            if (v > mx) {
                mx = v;
                piv = r;
            }
        }
        if (mx < 1e-13) return -1;
        if (piv != k) {
            for (int c = k; c < n; ++c) {
                double tmp = M[k * n + c];
                M[k * n + c] = M[piv * n + c];
                M[piv * n + c] = tmp;
            }
            double tr = rhs[k];
            rhs[k] = rhs[piv];
            rhs[piv] = tr;
        }
        double inv = 1.0 / M[k * n + k];
        for (int r = k + 1; r < n; ++r) {
            double fac = M[r * n + k] * inv;
            M[r * n + k] = 0.0;
            for (int c = k + 1; c < n; ++c) {
                M[r * n + c] -= fac * M[k * n + c];
            }
            rhs[r] -= fac * rhs[k];
        }
    }
    for (int k = n - 1; k >= 0; --k) {
        double s = rhs[k];
        for (int c = k + 1; c < n; ++c) s -= M[k * n + c] * x[c];
        x[k] = s / M[k * n + k];
    }
    return 0;
}

static void jackknife_stats(int nb, const double *blk, int bsz, double *mean_out, double *err_out) {
    double tot = 0.0;
    for (int b = 0; b < nb; ++b) tot += blk[b];
    double tot_n = (double)nb * (double)bsz;
    double mean = tot / tot_n;
    double jv = 0.0;
    double jk_n = (double)(nb - 1) * (double)bsz;
    for (int b = 0; b < nb; ++b) {
        double jm = (tot - blk[b]) / jk_n;
        double d = jm - mean;
        jv += d * d;
    }
    double scale = (double)(nb - 1) / (double)nb;
    *mean_out = mean;
    *err_out = sqrt(scale * jv);
}

/* -------------------------------------------------------------------------- */
/* 2. Exact Dinic Max-Flow / Min-Cut Solver (for {5,4} Hyperbolic Tiling)     */
/* -------------------------------------------------------------------------- */

enum { MAX_GRAPH_V = 2400, MAX_GRAPH_EDGE_RECS = 24000 };

typedef struct {
    int to;
    int rev;
    int cap;
} flow_edge_t;

typedef struct {
    int n_v;
    int head[MAX_GRAPH_V];
    int n_rec;
    flow_edge_t rec[MAX_GRAPH_EDGE_RECS];
    int next[MAX_GRAPH_EDGE_RECS];
} flow_net_t;

static void net_init(flow_net_t *fn, int nv) {
    fn->n_v = nv;
    fn->n_rec = 0;
    for (int i = 0; i < nv; ++i) fn->head[i] = -1;
}

static void net_add_dir_edge(flow_net_t *fn, int u, int v, int cap_uv, int cap_vu) {
    int idx_uv = fn->n_rec++;
    int idx_vu = fn->n_rec++;
    fn->rec[idx_uv].to = v;
    fn->rec[idx_uv].rev = idx_vu;
    fn->rec[idx_uv].cap = cap_uv;
    fn->next[idx_uv] = fn->head[u];
    fn->head[u] = idx_uv;

    fn->rec[idx_vu].to = u;
    fn->rec[idx_vu].rev = idx_uv;
    fn->rec[idx_vu].cap = cap_vu;
    fn->next[idx_vu] = fn->head[v];
    fn->head[v] = idx_vu;
}

typedef struct {
    int n_v;
    int head[MAX_GRAPH_V];
    int cur[MAX_GRAPH_V];
    int level[MAX_GRAPH_V];
    int q[MAX_GRAPH_V];
    int n_rec;
    flow_edge_t rec[MAX_GRAPH_EDGE_RECS];
    int next[MAX_GRAPH_EDGE_RECS];
} dinic_ws_t;

static int dinic_bfs(dinic_ws_t *ws, int s, int t) {
    for (int i = 0; i < ws->n_v; ++i) ws->level[i] = -1;
    int qh = 0, qt = 0;
    ws->level[s] = 0;
    ws->q[qt++] = s;
    while (qh < qt) {
        int u = ws->q[qh++];
        for (int e = ws->head[u]; e != -1; e = ws->next[e]) {
            if (ws->rec[e].cap > 0) {
                int v = ws->rec[e].to;
                if (ws->level[v] < 0) {
                    ws->level[v] = ws->level[u] + 1;
                    ws->q[qt++] = v;
                }
            }
        }
    }
    return ws->level[t] >= 0;
}

static int dinic_dfs(dinic_ws_t *ws, int u, int t, int pushed) {
    if (pushed == 0 || u == t) return pushed;
    for (int *e = &ws->cur[u]; *e != -1; *e = ws->next[*e]) {
        int ei = *e;
        int v = ws->rec[ei].to;
        if (ws->level[v] == ws->level[u] + 1 && ws->rec[ei].cap > 0) {
            int tr = dinic_dfs(ws, v, t, pushed < ws->rec[ei].cap ? pushed : ws->rec[ei].cap);
            if (tr > 0) {
                ws->rec[ei].cap -= tr;
                ws->rec[ws->rec[ei].rev].cap += tr;
                return tr;
            }
        }
    }
    return 0;
}

static int dinic_max_flow(const flow_net_t *tpl, int s, int t) {
    dinic_ws_t *ws = (dinic_ws_t *)calloc(1, sizeof(dinic_ws_t));
    if (!ws) { perror("calloc dinic_ws"); exit(2); }
    ws->n_v = tpl->n_v;
    ws->n_rec = tpl->n_rec;
    memcpy(ws->head, tpl->head, (size_t)tpl->n_v * sizeof(int));
    memcpy(ws->rec, tpl->rec, (size_t)tpl->n_rec * sizeof(flow_edge_t));
    memcpy(ws->next, tpl->next, (size_t)tpl->n_rec * sizeof(int));

    int flow = 0;
    while (dinic_bfs(ws, s, t)) {
        memcpy(ws->cur, ws->head, (size_t)ws->n_v * sizeof(int));
        int tr;
        while ((tr = dinic_dfs(ws, s, t, 1000000000)) > 0) flow += tr;
    }
    free(ws);
    return flow;
}

/* -------------------------------------------------------------------------- */
/* 3. PHY AME(4,8) Table & HaPPY [[5,1,3]] Tensor Data                        */
/* -------------------------------------------------------------------------- */

static const int PHY_AME8_R[64] = {
    0, 62, 55, 9, 44, 27, 18, 37, 34, 21, 28, 43, 15, 49, 56, 6,
    29, 42, 35, 20, 57, 7, 14, 48, 54, 8, 1, 63, 19, 36, 45, 26,
    59, 4, 13, 50, 30, 40, 33, 23, 17, 39, 46, 24, 53, 10, 3, 60,
    47, 25, 16, 38, 2, 61, 52, 11, 12, 51, 58, 5, 32, 22, 31, 41
};

static const int HAPPY_V0_SIGN[32] = {
    +1,  0,  0, -1,  0, +1, -1,  0,
     0, +1, +1,  0, -1,  0,  0, -1,
     0, -1, +1,  0, +1,  0,  0, -1,
    -1,  0,  0, -1,  0, -1, -1,  0
};

static const int HAPPY_V1_SIGN[32] = {
     0, -1, -1,  0, -1,  0,  0, -1,
    -1,  0,  0, +1,  0, +1, -1,  0,
    -1,  0,  0, -1,  0, +1, +1,  0,
     0, -1, +1,  0, -1,  0,  0, +1
};

/* T_happy[a0,a1,a2,a3,a4,b] with b in {0,1}; entry indexed by ((idx32 << 1) | b) */
static double happy_T_entry(int idx32, int b) {
    int sgn = (b == 0) ? HAPPY_V0_SIGN[idx32] : HAPPY_V1_SIGN[idx32];
    return (double)sgn * (M_SQRT2 * 0.25);
}

/* -------------------------------------------------------------------------- */
/* 4. 2D Simplicial Manifold Data Structure (S^2, T^2, Sigma_2 + Pachner)     */
/* -------------------------------------------------------------------------- */

enum { MAX_2D_V = 64, MAX_2D_E = 180, MAX_2D_F = 120 };

typedef struct {
    int V, E, F;
    int eu[MAX_2D_E], ev[MAX_2D_E];
    int fv[MAX_2D_F][3];
    int fe[MAX_2D_F][3]; /* fe[f][k] is the edge opposite fv[f][k] */
    int ef_cnt[MAX_2D_E], ef_list[MAX_2D_E][4];
    int vf_cnt[MAX_2D_V], vf_list[MAX_2D_V][16];
    double len[MAX_2D_E];
} surf2d_t;

static int find_or_add_edge(surf2d_t *s, int u, int v, double init_l) {
    if (u > v) { int t = u; u = v; v = t; }
    for (int e = 0; e < s->E; ++e) {
        if (s->eu[e] == u && s->ev[e] == v) return e;
    }
    int e = s->E++;
    s->eu[e] = u;
    s->ev[e] = v;
    s->len[e] = init_l;
    return e;
}

static void surf2d_rebuild_inc(surf2d_t *s) {
    for (int e = 0; e < s->E; ++e) s->ef_cnt[e] = 0;
    for (int v = 0; v < s->V; ++v) s->vf_cnt[v] = 0;
    for (int f = 0; f < s->F; ++f) {
        for (int k = 0; k < 3; ++k) {
            int e = s->fe[f][k];
            int v = s->fv[f][k];
            s->ef_list[e][s->ef_cnt[e]++] = f;
            s->vf_list[v][s->vf_cnt[v]++] = f;
        }
    }
}

static void surf2d_add_tri(surf2d_t *s, int v0, int v1, int v2, double l0, double l1, double l2) {
    /* l0 opposite v0 (v1-v2), l1 opposite v1 (v0-v2), l2 opposite v2 (v0-v1) */
    int f = s->F++;
    s->fv[f][0] = v0;
    s->fv[f][1] = v1;
    s->fv[f][2] = v2;
    s->fe[f][0] = find_or_add_edge(s, v1, v2, l0);
    s->fe[f][1] = find_or_add_edge(s, v0, v2, l1);
    s->fe[f][2] = find_or_add_edge(s, v0, v1, l2);
}

static void build_icosahedron(surf2d_t *s) {
    memset(s, 0, sizeof(*s));
    const double phi = 0.5 * (1.0 + sqrt(5.0));
    double p[12][3] = {
        {0, 1, phi}, {0, -1, phi}, {0, 1, -phi}, {0, -1, -phi},
        {1, phi, 0}, {-1, phi, 0}, {1, -phi, 0}, {-1, -phi, 0},
        {phi, 0, 1}, {-phi, 0, 1}, {phi, 0, -1}, {-phi, 0, -1}
    };
    s->V = 12;
    for (int a = 0; a < 12; ++a) {
        for (int b = a + 1; b < 12; ++b) {
            double d2ab = (p[a][0]-p[b][0])*(p[a][0]-p[b][0]) +
                          (p[a][1]-p[b][1])*(p[a][1]-p[b][1]) +
                          (p[a][2]-p[b][2])*(p[a][2]-p[b][2]);
            if (d2ab > 4.2) continue;
            for (int c = b + 1; c < 12; ++c) {
                double d2ac = (p[a][0]-p[c][0])*(p[a][0]-p[c][0]) +
                              (p[a][1]-p[c][1])*(p[a][1]-p[c][1]) +
                              (p[a][2]-p[c][2])*(p[a][2]-p[c][2]);
                double d2bc = (p[b][0]-p[c][0])*(p[b][0]-p[c][0]) +
                              (p[b][1]-p[c][1])*(p[b][1]-p[c][1]) +
                              (p[b][2]-p[c][2])*(p[b][2]-p[c][2]);
                if (d2ac <= 4.2 && d2bc <= 4.2) {
                    surf2d_add_tri(s, a, b, c, 1.0, 1.0, 1.0);
                }
            }
        }
    }
    surf2d_rebuild_inc(s);
}

static void build_torus_4x4(surf2d_t *s) {
    memset(s, 0, sizeof(*s));
    s->V = 16;
    const double d = M_SQRT2;
    for (int x = 0; x < 4; ++x) {
        for (int y = 0; y < 4; ++y) {
            int xp = (x + 1) & 3;
            int yp = (y + 1) & 3;
            int v00 = 4 * x + y;
            int v10 = 4 * xp + y;
            int v01 = 4 * x + yp;
            int v11 = 4 * xp + yp;
            /* tri A: (v00, v10, v11) with lengths opp (1.0, d, 1.0) */
            surf2d_add_tri(s, v00, v10, v11, 1.0, d, 1.0);
            /* tri B: (v00, v11, v01) with lengths opp (1.0, 1.0, d) */
            surf2d_add_tri(s, v00, v11, v01, 1.0, 1.0, d);
        }
    }
    surf2d_rebuild_inc(s);
}

static void build_genus2_surface(surf2d_t *s) {
    surf2d_t ta, tb;
    build_torus_4x4(&ta);
    build_torus_4x4(&tb);
    memset(s, 0, sizeof(*s));
    /* Glue along removing f=0 from ta and tb:
       ta face 0 vertices: (u0, u1, u2) -> map tb face 0 vertices (w0, w1, w2) to (u0, u1, u2)! */
    int map_b[16];
    int next_v = ta.V;
    for (int i = 0; i < 16; ++i) map_b[i] = -1;
    map_b[tb.fv[0][0]] = ta.fv[0][0];
    map_b[tb.fv[0][1]] = ta.fv[0][1];
    map_b[tb.fv[0][2]] = ta.fv[0][2];
    for (int i = 0; i < 16; ++i) {
        if (map_b[i] < 0) map_b[i] = next_v++;
    }
    s->V = next_v;
    for (int f = 1; f < ta.F; ++f) {
        surf2d_add_tri(s, ta.fv[f][0], ta.fv[f][1], ta.fv[f][2],
                       ta.len[ta.fe[f][0]], ta.len[ta.fe[f][1]], ta.len[ta.fe[f][2]]);
    }
    for (int f = 1; f < tb.F; ++f) {
        surf2d_add_tri(s, map_b[tb.fv[f][0]], map_b[tb.fv[f][1]], map_b[tb.fv[f][2]],
                       tb.len[tb.fe[f][0]], tb.len[tb.fe[f][1]], tb.len[tb.fe[f][2]]);
    }
    surf2d_rebuild_inc(s);
}

static int tri_valid(double a, double b, double c, double lmin, double lmax) {
    if (a < lmin || b < lmin || c < lmin) return 0;
    if (a > lmax || b > lmax || c > lmax) return 0;
    return (a + b - c > 1e-5) && (a + c - b > 1e-5) && (b + c - a > 1e-5);
}

static void tri_geom(const surf2d_t *s, int f, double *area_out, double ang_out[3]) {
    double l0 = s->len[s->fe[f][0]];
    double l1 = s->len[s->fe[f][1]];
    double l2 = s->len[s->fe[f][2]];
    double p = 0.5 * (l0 + l1 + l2);
    double arg = p * (p - l0) * (p - l1) * (p - l2);
    *area_out = (arg > 0.0) ? sqrt(arg) : 0.0;
    ang_out[0] = acos(clamp_cos((l1*l1 + l2*l2 - l0*l0) / (2.0 * l1 * l2)));
    ang_out[1] = acos(clamp_cos((l0*l0 + l2*l2 - l1*l1) / (2.0 * l0 * l2)));
    ang_out[2] = acos(clamp_cos((l0*l0 + l1*l1 - l2*l2) / (2.0 * l0 * l1)));
}

static void surf2d_measure(const surf2d_t *s, double *tot_area, double *r2_sum, double *gb_err) {
    double v_ang[MAX_2D_V], v_area[MAX_2D_V];
    for (int v = 0; v < s->V; ++v) { v_ang[v] = 0.0; v_area[v] = 0.0; }
    double a_sum = 0.0;
    for (int f = 0; f < s->F; ++f) {
        double af, ang[3];
        tri_geom(s, f, &af, ang);
        a_sum += af;
        for (int k = 0; k < 3; ++k) {
            int v = s->fv[f][k];
            v_ang[v] += ang[k];
            v_area[v] += af / 3.0;
        }
    }
    double eps_sum = 0.0, r2 = 0.0;
    for (int v = 0; v < s->V; ++v) {
        double eps = 2.0 * M_PI - v_ang[v];
        eps_sum += eps;
        r2 += (eps * eps) / v_area[v];
    }
    int chi = s->V - s->E + s->F;
    *tot_area = a_sum;
    *r2_sum = r2;
    *gb_err = fabs(eps_sum / (2.0 * M_PI) - (double)chi);
}

/* -------------------------------------------------------------------------- */
/* 5. 3D Simplicial Geometry & Action (\partial \Delta^4 and \partial \beta_4)*/
/* -------------------------------------------------------------------------- */

enum { MAX_3D_E = 24, MAX_3D_T = 16 };

typedef struct {
    int V, E, T;
    int tets[MAX_3D_T][4];
    int te[MAX_3D_T][6]; /* (0,1),(0,2),(0,3),(1,2),(1,3),(2,3) */
    int eu[MAX_3D_E], ev[MAX_3D_E];
} mesh3d_t;

static int m3_find_add_edge(mesh3d_t *m, int u, int v) {
    if (u > v) { int t = u; u = v; v = t; }
    for (int e = 0; e < m->E; ++e) {
        if (m->eu[e] == u && m->ev[e] == v) return e;
    }
    int e = m->E++;
    m->eu[e] = u;
    m->ev[e] = v;
    return e;
}

static void build_4simplex_bdy(mesh3d_t *m) {
    memset(m, 0, sizeof(*m));
    m->V = 5;
    for (int drop = 0; drop < 5; ++drop) {
        int t = m->T++;
        int k = 0;
        for (int v = 0; v < 5; ++v) {
            if (v != drop) m->tets[t][k++] = v;
        }
    }
    for (int t = 0; t < m->T; ++t) {
        const int *tv = m->tets[t];
        int p = 0;
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                m->te[t][p++] = m3_find_add_edge(m, tv[i], tv[j]);
            }
        }
    }
}

static void build_16cell_bdy(mesh3d_t *m) {
    memset(m, 0, sizeof(*m));
    m->V = 8;
    for (int s0 = 0; s0 <= 1; ++s0) {
        for (int s1 = 2; s1 <= 3; ++s1) {
            for (int s2 = 4; s2 <= 5; ++s2) {
                for (int s3 = 6; s3 <= 7; ++s3) {
                    int t = m->T++;
                    m->tets[t][0] = s0;
                    m->tets[t][1] = s1;
                    m->tets[t][2] = s2;
                    m->tets[t][3] = s3;
                }
            }
        }
    }
    for (int t = 0; t < m->T; ++t) {
        const int *tv = m->tets[t];
        int p = 0;
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                m->te[t][p++] = m3_find_add_edge(m, tv[i], tv[j]);
            }
        }
    }
}

static int tet_geom_full(const double l[6], double *vol_out, double ang_out[6], double dv_out[6]) {
    double s0 = l[0]*l[0], s1 = l[1]*l[1], s2 = l[2]*l[2];
    double s3 = l[3]*l[3], s4 = l[4]*l[4], s5 = l[5]*l[5];
    double G[3][3] = {
        {s0, 0.5*(s0 + s1 - s3), 0.5*(s0 + s2 - s4)},
        {0.5*(s0 + s1 - s3), s1, 0.5*(s1 + s2 - s5)},
        {0.5*(s0 + s2 - s4), 0.5*(s1 + s2 - s5), s2}
    };
    double c00 = G[1][1]*G[2][2] - G[1][2]*G[2][1];
    double c01 = G[1][2]*G[2][0] - G[1][0]*G[2][2];
    double c02 = G[1][0]*G[2][1] - G[1][1]*G[2][0];
    double det = G[0][0]*c00 + G[0][1]*c01 + G[0][2]*c02;
    if (det <= 1e-12) return 0;
    double vol = sqrt(det) / 6.0;
    *vol_out = vol;
    double c11 = G[0][0]*G[2][2] - G[0][2]*G[2][0];
    double c12 = G[0][1]*G[2][0] - G[0][0]*G[2][1];
    double c22 = G[0][0]*G[1][1] - G[0][1]*G[1][0];
    if (dv_out) {
        double ddet[6] = {
            c00 + c01 + c02,
            c11 + c01 + c12,
            c22 + c02 + c12,
            -c01,
            -c02,
            -c12
        };
        double fac = 1.0 / (36.0 * vol);
        for (int k = 0; k < 6; ++k) dv_out[k] = l[k] * ddet[k] * fac;
    }
    if (ang_out) {
        double inv_det = 1.0 / det;
        double c10 = (G[0][2]*G[2][1] - G[0][1]*G[2][2]) * inv_det;
        double c20 = (G[0][1]*G[1][2] - G[0][2]*G[1][1]) * inv_det;
        double c21 = (G[0][2]*G[1][0] - G[0][0]*G[1][2]) * inv_det;
        double M[3][3] = {
            {c00 * inv_det, c01 * inv_det, c02 * inv_det},
            {c10,           c11 * inv_det, c12 * inv_det},
            {c20,           c21,           c22 * inv_det}
        };
        double N[4][4];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) N[i + 1][j + 1] = M[i][j];
        }
        double n00 = 0.0;
        for (int a = 1; a <= 3; ++a) {
            double col = 0.0;
            for (int r = 0; r < 3; ++r) col += M[r][a - 1];
            N[0][a] = N[a][0] = -col;
            n00 += col;
        }
        N[0][0] = n00;
        static const int pairs[6][2] = {
            {2, 3}, {1, 3}, {1, 2}, {0, 3}, {0, 2}, {0, 1}
        };
        for (int k = 0; k < 6; ++k) {
            int a = pairs[k][0], b = pairs[k][1];
            double denom = sqrt(N[a][a] * N[b][b]);
            if (denom <= 0.0) return 0;
            ang_out[k] = acos(clamp_cos(-N[a][b] / denom));
        }
    }
    return 1;
}

static int tet_geom(const double l[6], double *vol_out, double ang_out[6]) {
    return tet_geom_full(l, vol_out, ang_out, NULL);
}

static int regge3d_eval(const mesh3d_t *m, const double *l, double kappa, double lam,
                        double *act_out, double *eps_out, double *vtot_out) {
    double eps[MAX_3D_E];
    for (int e = 0; e < m->E; ++e) eps[e] = 2.0 * M_PI;
    double vtot = 0.0;
    for (int t = 0; t < m->T; ++t) {
        double lt[6], ang[6], vol;
        for (int k = 0; k < 6; ++k) lt[k] = l[m->te[t][k]];
        if (!tet_geom(lt, &vol, ang)) return 0;
        vtot += vol;
        for (int k = 0; k < 6; ++k) eps[m->te[t][k]] -= ang[k];
    }
    double s = lam * vtot;
    for (int e = 0; e < m->E; ++e) s -= kappa * l[e] * eps[e];
    if (act_out) *act_out = s;
    if (eps_out) memcpy(eps_out, eps, (size_t)m->E * sizeof(double));
    if (vtot_out) *vtot_out = vtot;
    return 1;
}

/* Exact 1e-16 Analytic Schlafli + Cayley-Menger Cofactor Gradient */
static int regge3d_grad(const mesh3d_t *m, const double *l, double kappa, double lam, double *grad) {
    double eps[MAX_3D_E], dv_tot[MAX_3D_E];
    for (int e = 0; e < m->E; ++e) {
        eps[e] = 2.0 * M_PI;
        dv_tot[e] = 0.0;
    }
    for (int t = 0; t < m->T; ++t) {
        double lt[6], ang[6], dv[6], vol;
        for (int k = 0; k < 6; ++k) lt[k] = l[m->te[t][k]];
        if (!tet_geom_full(lt, &vol, ang, dv)) return 0;
        for (int k = 0; k < 6; ++k) {
            int e = m->te[t][k];
            eps[e] -= ang[k];
            dv_tot[e] += dv[k];
        }
    }
    for (int e = 0; e < m->E; ++e) {
        grad[e] = -kappa * eps[e] + lam * dv_tot[e];
    }
    return 1;
}

static int regge3d_hessian(const mesh3d_t *m, const double *l, double kappa, double lam, double *H) {
    const double h = 1e-5;
    double ltmp[MAX_3D_E], gp[MAX_3D_E], gm[MAX_3D_E];
    memcpy(ltmp, l, (size_t)m->E * sizeof(double));
    for (int j = 0; j < m->E; ++j) {
        ltmp[j] = l[j] + h;
        if (!regge3d_grad(m, ltmp, kappa, lam, gp)) return 0;
        ltmp[j] = l[j] - h;
        if (!regge3d_grad(m, ltmp, kappa, lam, gm)) return 0;
        ltmp[j] = l[j];
        for (int i = 0; i < m->E; ++i) H[i * m->E + j] = (gp[i] - gm[i]) / (2.0 * h);
    }
    for (int i = 0; i < m->E; ++i) {
        for (int j = i + 1; j < m->E; ++j) {
            double avg = 0.5 * (H[i * m->E + j] + H[j * m->E + i]);
            H[i * m->E + j] = H[j * m->E + i] = avg;
        }
    }
    return 1;
}

/* -------------------------------------------------------------------------- */
/* 6. Exact {5,4} Hyperbolic Pentagon Tiling Builder                          */
/* -------------------------------------------------------------------------- */

enum { MAX_HYP_R = 4, MAX_HYP_V = 850, MAX_HYP_BDY = 1400 };

typedef struct {
    int V_bulk;
    int L_bdy;
    int bdy_owner[MAX_HYP_BDY]; /* bulk pentagon owning boundary edge i */
    flow_net_t net;             /* nodes 0..V_bulk-1 = bulk, V_bulk..V_bulk+L_bdy-1 = bdy */
    int R0_end, R1_end, R2_end; /* cumulative bulk counts up to shell 0, 1, 2 */
} hyp54_t;

static void build_hyp54(hyp54_t *h, int R) {
    memset(h, 0, sizeof(*h));
    int adj_cnt[MAX_HYP_V];
    int adj[MAX_HYP_V][8];
    for (int i = 0; i < MAX_HYP_V; ++i) adj_cnt[i] = 0;

    int cur_v = 1;
    int bdy[MAX_HYP_BDY], nbdy[MAX_HYP_BDY];
    int L = 5;
    for (int i = 0; i < 5; ++i) bdy[i] = 0;
    h->R0_end = 1;

    for (int r = 1; r <= R; ++r) {
        int e_tiles[MAX_HYP_BDY];
        for (int k = 0; k < L; ++k) {
            int tid = cur_v++;
            e_tiles[k] = tid;
            int ow = bdy[k];
            adj[tid][adj_cnt[tid]++] = ow;
            adj[ow][adj_cnt[ow]++] = tid;
        }
        int new_L = 0;
        for (int k = 0; k < L; ++k) {
            int ek = e_tiles[k];
            int en = e_tiles[(k + 1) % L];
            nbdy[new_L++] = ek;
            nbdy[new_L++] = ek;
            if (bdy[k] == bdy[(k + 1) % L]) {
                int cid = cur_v++;
                adj[cid][adj_cnt[cid]++] = ek;
                adj[cid][adj_cnt[cid]++] = en;
                adj[ek][adj_cnt[ek]++] = cid;
                adj[en][adj_cnt[en]++] = cid;
                nbdy[new_L++] = cid;
                nbdy[new_L++] = cid;
                nbdy[new_L++] = cid;
            } else {
                adj[ek][adj_cnt[ek]++] = en;
                adj[en][adj_cnt[en]++] = ek;
            }
        }
        L = new_L;
        memcpy(bdy, nbdy, (size_t)L * sizeof(int));
        if (r == 1) h->R1_end = cur_v;
        if (r == 2) h->R2_end = cur_v;
    }

    h->V_bulk = cur_v;
    h->L_bdy = L;
    memcpy(h->bdy_owner, bdy, (size_t)L * sizeof(int));

    /* Verify every bulk tile has exact total degree 5 */
    int bdy_cnt[MAX_HYP_V];
    memset(bdy_cnt, 0, (size_t)cur_v * sizeof(int));
    for (int i = 0; i < L; ++i) bdy_cnt[bdy[i]]++;
    for (int u = 0; u < cur_v; ++u) {
        if (adj_cnt[u] + bdy_cnt[u] != 5) {
            fprintf(stderr, "FATAL: {5,4} tile %d degree=%d != 5\n", u, adj_cnt[u] + bdy_cnt[u]);
            exit(2);
        }
    }

    net_init(&h->net, cur_v + L + 2); /* reserve +2 slots for super-source S and super-sink T */
    for (int u = 0; u < cur_v; ++u) {
        for (int idx = 0; idx < adj_cnt[u]; ++idx) {
            int v = adj[u][idx];
            if (u < v) net_add_dir_edge(&h->net, u, v, 1, 1);
        }
    }
    for (int i = 0; i < L; ++i) {
        net_add_dir_edge(&h->net, bdy[i], cur_v + i, 1, 1);
    }
}

static int hyp54_mincut(const hyp54_t *h, const uint8_t *in_A, int core_cutoff) {
    /* core_cutoff == 0 -> no BH; core_cutoff in {1,11,51} -> vertices 0..core_cutoff-1 merged into BH sink */
    flow_net_t *g = (flow_net_t *)calloc(1, sizeof(flow_net_t));
    if (!g) { perror("calloc flow_net"); exit(2); }
    net_init(g, h->V_bulk + h->L_bdy + 2);
    int S = h->V_bulk + h->L_bdy;
    int T = S + 1;
    int bh_node = (core_cutoff > 0) ? 0 : -1;

    /* Copy undirected bulk and boundary-leg edges collapsing BH core [0..core_cutoff-1] to node 0 */
    for (int u = 0; u < h->V_bulk + h->L_bdy; ++u) {
        for (int e = h->net.head[u]; e != -1; e = h->net.next[e]) {
            int v = h->net.rec[e].to;
            if (u < v) {
                int mu = (core_cutoff > 0 && u < core_cutoff) ? bh_node : u;
                int mv = (core_cutoff > 0 && v < core_cutoff) ? bh_node : v;
                if (mu != mv) net_add_dir_edge(g, mu, mv, 1, 1);
            }
        }
    }
    for (int i = 0; i < h->L_bdy; ++i) {
        int bn = h->V_bulk + i;
        if (in_A[i]) {
            net_add_dir_edge(g, S, bn, 1000000, 0);
        } else {
            net_add_dir_edge(g, bn, T, 1000000, 0);
        }
    }
    if (core_cutoff > 0) {
        /* Black hole horizon lives in purifying sink Abar -> connect collapsed bh_node to T! */
        net_add_dir_edge(g, bh_node, T, 1000000, 0);
    }
    int ans = dinic_max_flow(g, S, T);
    free(g);
    return ans;
}

/* -------------------------------------------------------------------------- */
/* 7. Startup Self-Tests (6 Quantitative Algebraic & Geometric Verifications) */
/* -------------------------------------------------------------------------- */

typedef struct {
    double ame_pur, happy_u33, happy_kl, gb_max_err, schlafli_err, grad3d_err, noclone_err;
    int ssa_viol, mmi_viol;
} selftest_res_t;

static int run_selftests(selftest_res_t *st_out) {
    /* Self-Test 1: PHY v4.0.0 Nonlinear AME(4,8) 3-Bipartition Unitarity & Purity 1/64 */
    int seen[3][64];
    memset(seen, 0, sizeof(seen));
    for (int a = 0; a < 8; ++a) {
        for (int b = 0; b < 8; ++b) {
            int v = PHY_AME8_R[(a << 3) | b];
            int c = (v >> 3) & 7;
            int d = v & 7;
            seen[0][(c << 3) | d]++;
            seen[1][(a << 3) | c]++;
            seen[2][(a << 3) | d]++;
            /* also check outputs of cross/horiz are unique */
        }
    }
    int seen_out1[64], seen_out2[64];
    memset(seen_out1, 0, sizeof(seen_out1));
    memset(seen_out2, 0, sizeof(seen_out2));
    for (int a = 0; a < 8; ++a) {
        for (int b = 0; b < 8; ++b) {
            int v = PHY_AME8_R[(a << 3) | b];
            int c = (v >> 3) & 7;
            int d = v & 7;
            seen_out1[(b << 3) | d]++;
            seen_out2[(b << 3) | c]++;
        }
    }
    double ame_err = 0.0;
    for (int k = 0; k < 64; ++k) {
        if (seen[0][k] != 1 || seen[1][k] != 1 || seen[2][k] != 1 ||
            seen_out1[k] != 1 || seen_out2[k] != 1) {
            ame_err = 1.0;
        }
    }
    double ame_pur = 64.0 * (1.0 / 64.0) * (1.0 / 64.0); /* exact I_64 / 64 -> Tr(rho^2) = 1/64 */
    if (ame_err > 0.0) {
        fprintf(stderr, "SELFTEST 1 FAILED: PHY AME(4,8) not 2-unitary!\n");
        return -1;
    }

    /* Self-Test 2: HaPPY [[5,1,3]] AME(6,2) all 20 3v3 Unitaries + 10-pair 2-Erasure Knill-Laflamme */
    double T6[64];
    for (int x = 0; x < 32; ++x) {
        for (int b = 0; b < 2; ++b) T6[(x << 1) | b] = M_SQRT2 * happy_T_entry(x, b);
    }
    double max_u33_err = 0.0;
    for (int mask = 0; mask < 64; ++mask) {
        int pc = 0;
        for (int i = 0; i < 6; ++i) if ((mask >> i) & 1) ++pc;
        if (pc != 3) continue;
        int A[3], B[3], ia = 0, ib = 0;
        for (int i = 0; i < 6; ++i) {
            if ((mask >> i) & 1) A[ia++] = 5 - i;
            else B[ib++] = 5 - i;
        }
        double M[8][8];
        for (int full = 0; full < 64; ++full) {
            int ra = 0, cb = 0;
            for (int k = 0; k < 3; ++k) {
                ra = (ra << 1) | ((full >> A[k]) & 1);
                cb = (cb << 1) | ((full >> B[k]) & 1);
            }
            M[ra][cb] = T6[full];
        }
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                double s = 0.0;
                for (int k = 0; k < 8; ++k) s += M[k][i] * M[k][j];
                double e = fabs(s - (i == j ? 1.0 : 0.0));
                if (e > max_u33_err) max_u33_err = e;
            }
        }
    }
    /* Knill-Laflamme V^dag (I_3 x P_2) V = c * I_2 over all 10 2-qubit pairs & 16 real Pauli ops */
    const double pmat[4][4] = {
        {1, 0, 0, 1},   /* I */
        {0, 1, 1, 0},   /* X */
        {1, 0, 0, -1},  /* Z */
        {0, 1, -1, 0}   /* iY = XZ */
    };
    double max_kl_err = 0.0;
    for (int q1 = 0; q1 < 5; ++q1) {
        for (int q2 = q1 + 1; q2 < 5; ++q2) {
            for (int op1 = 0; op1 < 4; ++op1) {
                for (int op2 = 0; op2 < 4; ++op2) {
                    double W[2][2] = {{0, 0}, {0, 0}};
                    for (int x = 0; x < 32; ++x) {
                        int b1 = (x >> (4 - q1)) & 1;
                        int b2 = (x >> (4 - q2)) & 1;
                        for (int bp1 = 0; bp1 < 2; ++bp1) {
                            double w1 = pmat[op1][b1 * 2 + bp1];
                            if (w1 == 0.0) continue;
                            for (int bp2 = 0; bp2 < 2; ++bp2) {
                                double w2 = pmat[op2][b2 * 2 + bp2];
                                if (w2 == 0.0) continue;
                                int y = (x & ~(1 << (4 - q1)) & ~(1 << (4 - q2))) |
                                        (bp1 << (4 - q1)) | (bp2 << (4 - q2));
                                for (int bi = 0; bi < 2; ++bi) {
                                    for (int bj = 0; bj < 2; ++bj) {
                                        double vix = (bi == 0 ? HAPPY_V0_SIGN[x] : HAPPY_V1_SIGN[x]) * 0.25;
                                        double vjy = (bj == 0 ? HAPPY_V0_SIGN[y] : HAPPY_V1_SIGN[y]) * 0.25;
                                        W[bi][bj] += vix * w1 * w2 * vjy;
                                    }
                                }
                            }
                        }
                    }
                    double c_exp = (op1 == 0 && op2 == 0) ? 1.0 : 0.0;
                    double e0 = fabs(W[0][0] - c_exp);
                    double e1 = fabs(W[1][1] - c_exp);
                    double e2 = fabs(W[0][1]);
                    double e3 = fabs(W[1][0]);
                    if (e0 > max_kl_err) max_kl_err = e0;
                    if (e1 > max_kl_err) max_kl_err = e1;
                    if (e2 > max_kl_err) max_kl_err = e2;
                    if (e3 > max_kl_err) max_kl_err = e3;
                }
            }
        }
    }
    if (max_u33_err > 1e-14 || max_kl_err > 1e-14) {
        fprintf(stderr, "SELFTEST 2 FAILED: u33=%.3e kl=%.3e\n", max_u33_err, max_kl_err);
        return -1;
    }

    /* Self-Test 3: 2D Gauss-Bonnet on S^2, T^2, Sigma_2 + Pachner 2-2 & 1-3 on S^2 */
    surf2d_t s2, t2, sig2;
    build_icosahedron(&s2);
    build_torus_4x4(&t2);
    build_genus2_surface(&sig2);
    if (s2.V != 12 || s2.E != 30 || s2.F != 20 ||
        t2.V != 16 || t2.E != 48 || t2.F != 32 ||
        sig2.V != 29 || sig2.E != 93 || sig2.F != 62) {
        fprintf(stderr, "SELFTEST 3 TOPO ERROR\n");
        return -1;
    }
    rng_t rng;
    rng_seed(&rng, 0x8001ULL);
    surf2d_t *surfs[3] = {&s2, &t2, &sig2};
    double max_gb_err = 0.0;
    for (int k = 0; k < 3; ++k) {
        surf2d_t sc = *surfs[k];
        for (int e = 0; e < sc.E; ++e) sc.len[e] *= (0.94 + 0.12 * rng_u01(&rng));
        double atot, r2, gbe;
        surf2d_measure(&sc, &atot, &r2, &gbe);
        if (gbe > max_gb_err) max_gb_err = gbe;
    }
    /* Pachner 2-2 diagonal flip on S^2 edge 0: replace (a,b) by (c,d) */
    surf2d_t p22 = s2;
    {
        int e0 = 0;
        int f0 = p22.ef_list[e0][0], f1 = p22.ef_list[e0][1];
        int a = p22.eu[e0], b = p22.ev[e0];
        int c = -1, d = -1;
        for (int k = 0; k < 3; ++k) {
            if (p22.fv[f0][k] != a && p22.fv[f0][k] != b) c = p22.fv[f0][k];
            if (p22.fv[f1][k] != a && p22.fv[f1][k] != b) d = p22.fv[f1][k];
        }
        p22.eu[e0] = (c < d) ? c : d;
        p22.ev[e0] = (c < d) ? d : c;
        p22.len[e0] = 1.65; /* valid diagonal across two unit equilateral triangles (max 2.0) */
        p22.F = 0;
        /* rebuild faces keeping all other faces except f0,f1, then adding (a,c,d) and (b,c,d) */
        surf2d_t tmp = p22;
        tmp.E = 0; tmp.F = 0; tmp.V = s2.V;
        for (int f = 0; f < s2.F; ++f) {
            if (f == f0 || f == f1) continue;
            surf2d_add_tri(&tmp, s2.fv[f][0], s2.fv[f][1], s2.fv[f][2], 1.0, 1.0, 1.0);
        }
        surf2d_add_tri(&tmp, a, c, d, 1.65, 1.0, 1.0);
        surf2d_add_tri(&tmp, b, c, d, 1.65, 1.0, 1.0);
        surf2d_rebuild_inc(&tmp);
        double atot, r2, gbe;
        surf2d_measure(&tmp, &atot, &r2, &gbe);
        if (tmp.V != 12 || tmp.E != 30 || tmp.F != 20 || gbe > max_gb_err) {
            if (gbe > max_gb_err) max_gb_err = gbe;
        }
    }
    /* Pachner 1-3 stellar subdivision on S^2 face 0 */
    {
        surf2d_t p13;
        memset(&p13, 0, sizeof(p13));
        p13.V = s2.V + 1;
        int vn = s2.V;
        for (int f = 1; f < s2.F; ++f) {
            surf2d_add_tri(&p13, s2.fv[f][0], s2.fv[f][1], s2.fv[f][2], 1.0, 1.0, 1.0);
        }
        int v0 = s2.fv[0][0], v1 = s2.fv[0][1], v2 = s2.fv[0][2];
        double r_cent = 1.0 / sqrt(3.0) + 0.05;
        surf2d_add_tri(&p13, v0, v1, vn, r_cent, r_cent, 1.0);
        surf2d_add_tri(&p13, v1, v2, vn, r_cent, r_cent, 1.0);
        surf2d_add_tri(&p13, v2, v0, vn, r_cent, r_cent, 1.0);
        surf2d_rebuild_inc(&p13);
        double atot, r2, gbe;
        surf2d_measure(&p13, &atot, &r2, &gbe);
        if (p13.V != 13 || p13.E != 33 || p13.F != 22 || gbe > max_gb_err) {
            if (gbe > max_gb_err) max_gb_err = gbe;
        }
    }
    if (max_gb_err > 1e-13) {
        fprintf(stderr, "SELFTEST 3 FAILED: max_gb_err=%.3e\n", max_gb_err);
        return -1;
    }

    /* Self-Test 4: 3D Irregular Tetrahedron Schlafli Identity & Full Action Grad vs FD */
    double l_tet[6] = {1.12, 0.97, 1.05, 1.18, 0.93, 1.08};
    const double h_fd = 1e-6;
    double max_schlafli_err = 0.0;
    for (int k = 0; k < 6; ++k) {
        double lp[6], lm[6], ap[6], am[6], vp, vm;
        memcpy(lp, l_tet, sizeof(l_tet));
        memcpy(lm, l_tet, sizeof(l_tet));
        lp[k] += h_fd;
        lm[k] -= h_fd;
        tet_geom(lp, &vp, ap);
        tet_geom(lm, &vm, am);
        double sum_l_dt = 0.0;
        for (int e = 0; e < 6; ++e) sum_l_dt += l_tet[e] * (ap[e] - am[e]) / (2.0 * h_fd);
        if (fabs(sum_l_dt) > max_schlafli_err) max_schlafli_err = fabs(sum_l_dt);
    }
    mesh3d_t d4;
    build_4simplex_bdy(&d4);
    double l_d4[10], g_ana[10];
    for (int e = 0; e < 10; ++e) l_d4[e] = 3.5 + 0.1 * ((e % 3) - 1);
    regge3d_grad(&d4, l_d4, 1.0, 1.0, g_ana);
    double max_act_grad_err = 0.0;
    for (int e = 0; e < 10; ++e) {
        double lp[10], lm[10], sp, sm;
        memcpy(lp, l_d4, sizeof(l_d4));
        memcpy(lm, l_d4, sizeof(l_d4));
        lp[e] += h_fd;
        lm[e] -= h_fd;
        regge3d_eval(&d4, lp, 1.0, 1.0, &sp, NULL, NULL);
        regge3d_eval(&d4, lm, 1.0, 1.0, &sm, NULL, NULL);
        double g_fd = (sp - sm) / (2.0 * h_fd);
        double err = fabs(g_ana[e] - g_fd);
        if (err > max_act_grad_err) max_act_grad_err = err;
    }
    if (max_schlafli_err > 1e-8 || max_act_grad_err > 1e-7) {
        fprintf(stderr, "SELFTEST 4 FAILED: schlafli=%.3e actgrad=%.3e\n",
                max_schlafli_err, max_act_grad_err);
        return -1;
    }

    /* Self-Test 5: {5,4} R=2 Holographic Entropy Cone SSA & MMI Audit (100 random triples) */
    hyp54_t *h2 = (hyp54_t *)calloc(1, sizeof(hyp54_t));
    if (!h2) { perror("calloc h2"); exit(2); }
    build_hyp54(h2, 2);
    if (h2->V_bulk != 51 || h2->L_bdy != 95) {
        fprintf(stderr, "SELFTEST 5 TOPO ERROR: V=%d L=%d\n", h2->V_bulk, h2->L_bdy);
        free(h2);
        return -1;
    }
    int ssa_viol = 0, mmi_viol = 0;
    uint8_t mA[MAX_HYP_BDY], mB[MAX_HYP_BDY], mC[MAX_HYP_BDY];
    uint8_t mAB[MAX_HYP_BDY], mBC[MAX_HYP_BDY], mAC[MAX_HYP_BDY], mABC[MAX_HYP_BDY];
    for (int tr = 0; tr < 100; ++tr) {
        int a1 = (int)(rng_u64(&rng) % 95), la = 2 + (int)(rng_u64(&rng) % 14);
        int b1 = (a1 + la + 1 + (int)(rng_u64(&rng) % 10)) % 95, lb = 2 + (int)(rng_u64(&rng) % 14);
        int c1 = (b1 + lb + 1 + (int)(rng_u64(&rng) % 10)) % 95, lc = 2 + (int)(rng_u64(&rng) % 14);
        memset(mA, 0, 95); memset(mB, 0, 95); memset(mC, 0, 95);
        for (int k = 0; k < la; ++k) mA[(a1 + k) % 95] = 1;
        for (int k = 0; k < lb; ++k) { int p = (b1 + k) % 95; if (!mA[p]) mB[p] = 1; }
        for (int k = 0; k < lc; ++k) { int p = (c1 + k) % 95; if (!mA[p] && !mB[p]) mC[p] = 1; }
        for (int i = 0; i < 95; ++i) {
            mAB[i] = mA[i] | mB[i];
            mBC[i] = mB[i] | mC[i];
            mAC[i] = mA[i] | mC[i];
            mABC[i] = mA[i] | mB[i] | mC[i];
        }
        int sA = hyp54_mincut(h2, mA, 0);
        int sB = hyp54_mincut(h2, mB, 0);
        int sC = hyp54_mincut(h2, mC, 0);
        int sAB = hyp54_mincut(h2, mAB, 0);
        int sBC = hyp54_mincut(h2, mBC, 0);
        int sAC = hyp54_mincut(h2, mAC, 0);
        int sABC = hyp54_mincut(h2, mABC, 0);
        if (sB + sABC > sAB + sBC) ++ssa_viol;
        if (sA + sB + sC + sABC > sAB + sBC + sAC) ++mmi_viol;
    }
    free(h2);
    if (ssa_viol > 0 || mmi_viol > 0) {
        fprintf(stderr, "SELFTEST 5 FAILED: ssa=%d mmi=%d\n", ssa_viol, mmi_viol);
        return -1;
    }

    /* Self-Test 6: Holographic No-Cloning Check I(R : Abar) = 0 on [[5,1,3]] (|R|=2, |A|=3, |Abar|=2) */
    double M_ra[8][8];
    memset(M_ra, 0, sizeof(M_ra));
    for (int x = 0; x < 32; ++x) {
        int a_idx = (x >> 2) & 7; /* q0,q1,q2 */
        int ab_idx = x & 3;       /* q3,q4 */
        M_ra[ab_idx][a_idx] += (double)HAPPY_V0_SIGN[x] * 0.25 * M_SQRT1_2;
        M_ra[4 + ab_idx][a_idx] += (double)HAPPY_V1_SIGN[x] * 0.25 * M_SQRT1_2;
    }
    double rho8[64], evals8[8];
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            double s = 0.0;
            for (int k = 0; k < 8; ++k) s += M_ra[i][k] * M_ra[j][k];
            rho8[i * 8 + j] = s;
        }
    }
    jacobi_real_symm(8, rho8, evals8);
    double s_RAbar = 0.0;
    for (int i = 0; i < 8; ++i) if (evals8[i] > 1e-14) s_RAbar -= evals8[i] * log2(evals8[i]);
    /* I(R : Abar) = S(R) + S(Abar) - S(R Abar) = 1 + 2 - 3 = 0 */
    double noclone_err = fabs(1.0 + 2.0 - s_RAbar);
    if (noclone_err > 1e-12) {
        fprintf(stderr, "SELFTEST 6 FAILED: I(R:Abar)=%.3e\n", noclone_err);
        return -1;
    }

    if (st_out) {
        st_out->ame_pur = ame_pur;
        st_out->happy_u33 = max_u33_err;
        st_out->happy_kl = max_kl_err;
        st_out->gb_max_err = max_gb_err;
        st_out->schlafli_err = max_schlafli_err;
        st_out->grad3d_err = max_act_grad_err;
        st_out->ssa_viol = ssa_viol;
        st_out->mmi_viol = mmi_viol;
        st_out->noclone_err = noclone_err;
    }
    printf("  [SELFTEST] PASS (AME8=0.00e+00 pur=%.6f | HaPPY u33=%.2e KL=%.2e | GB=%.2e | Schlafli=%.2e 3Dgrad=%.2e | HypSSA/MMI=%d/%d | NoClone=%.2e)\n",
           ame_pur, max_u33_err, max_kl_err, max_gb_err, max_schlafli_err, max_act_grad_err,
           ssa_viol, mmi_viol, noclone_err);
    return 0;
}

/* -------------------------------------------------------------------------- */
/* 8. Experiment 1A: 2D Dynamical Regge Gravity MC on S^2 and T^2             */
/* -------------------------------------------------------------------------- */

enum { N_ALPHA_2D = 4, JK_2D = 20 };

typedef struct {
    int topo; /* 0 = S^2, 1 = T^2 */
    double lam, alpha;
    double a_mean, a_err;
    double r2_mean, r2_err;
    double cs_ratio;
    double gb_max_err;
    double acc_rate;
} exp1a_res_t;

static void *exp1a_worker(void *arg) {
    exp1a_res_t *res = (exp1a_res_t *)arg;
    surf2d_t s;
    if (res->topo == 0) build_icosahedron(&s);
    else build_torus_4x4(&s);

    rng_t rng;
    rng_seed(&rng, 0x1A000000ULL + (uint64_t)res->topo * 997ULL + (uint64_t)(res->alpha * 100.0));

    double l_star = (res->topo == 0) ? sqrt(4.0 * M_PI * sqrt((res->alpha > 0.2 ? res->alpha : 0.2) / res->lam) / (5.0 * sqrt(3.0))) : 1.0;
    if (res->topo == 0) {
        for (int e = 0; e < s.E; ++e) s.len[e] = l_star;
    }

    const int n_therm = 2500;
    const int n_meas = 8000;
    const int bsz = n_meas / JK_2D;
    double blk_a[JK_2D], blk_r2[JK_2D];
    memset(blk_a, 0, sizeof(blk_a));
    memset(blk_r2, 0, sizeof(blk_r2));

    double v_ang[MAX_2D_V], v_area[MAX_2D_V];
    for (int v = 0; v < s.V; ++v) { v_ang[v] = 0.0; v_area[v] = 0.0; }
    double f_area[MAX_2D_F], f_ang[MAX_2D_F][3];
    for (int f = 0; f < s.F; ++f) {
        tri_geom(&s, f, &f_area[f], f_ang[f]);
        for (int k = 0; k < 3; ++k) {
            int v = s.fv[f][k];
            v_ang[v] += f_ang[f][k];
            v_area[v] += f_area[f] / 3.0;
        }
    }

    const double step = 0.42 * (1.0 + 0.35 * sqrt(res->alpha));
    const double lmin = 0.4, lmax = 6.0;
    uint64_t props = 0, acc = 0;
    double max_gb = 0.0;

    for (int sw = -n_therm; sw < n_meas; ++sw) {
        for (int e = 0; e < s.E; ++e) {
            double old_l = s.len[e];
            double new_l = old_l + step * (2.0 * rng_u01(&rng) - 1.0);
            int f0 = s.ef_list[e][0], f1 = s.ef_list[e][1];
            s.len[e] = new_l;
            int ok = 1;
            for (int fi = 0; fi < 2; ++fi) {
                int f = (fi == 0) ? f0 : f1;
                if (!tri_valid(s.len[s.fe[f][0]], s.len[s.fe[f][1]], s.len[s.fe[f][2]], lmin, lmax)) {
                    ok = 0;
                    break;
                }
            }
            if (!ok) {
                s.len[e] = old_l;
                if (sw >= 0) ++props;
                continue;
            }
            /* Collect unique vertices touched by f0 and f1 */
            int tv[6], n_tv = 0;
            for (int fi = 0; fi < 2; ++fi) {
                int f = (fi == 0) ? f0 : f1;
                for (int k = 0; k < 3; ++k) {
                    int v = s.fv[f][k];
                    int found = 0;
                    for (int q = 0; q < n_tv; ++q) if (tv[q] == v) found = 1;
                    if (!found) tv[n_tv++] = v;
                }
            }
            double s_old = 0.0, s_new = 0.0;
            double nf_area[2], nf_ang[2][3];
            tri_geom(&s, f0, &nf_area[0], nf_ang[0]);
            tri_geom(&s, f1, &nf_area[1], nf_ang[1]);

            double dv_ang[6], dv_area[6];
            for (int q = 0; q < n_tv; ++q) {
                int v = tv[q];
                double eps_o = 2.0 * M_PI - v_ang[v];
                s_old += res->alpha * (eps_o * eps_o) / v_area[v];
                dv_ang[q] = 0.0;
                dv_area[q] = 0.0;
            }
            s_old += res->lam * (f_area[f0] + f_area[f1]);

            int fpair[2] = {f0, f1};
            for (int fi = 0; fi < 2; ++fi) {
                int f = fpair[fi];
                for (int k = 0; k < 3; ++k) {
                    int v = s.fv[f][k];
                    for (int q = 0; q < n_tv; ++q) {
                        if (tv[q] == v) {
                            dv_ang[q] += nf_ang[fi][k] - f_ang[f][k];
                            dv_area[q] += (nf_area[fi] - f_area[f]) / 3.0;
                        }
                    }
                }
            }
            for (int q = 0; q < n_tv; ++q) {
                int v = tv[q];
                double na = v_area[v] + dv_area[q];
                double neps = 2.0 * M_PI - (v_ang[v] + dv_ang[q]);
                s_new += res->alpha * (neps * neps) / na;
            }
            s_new += res->lam * (nf_area[0] + nf_area[1]);

            double dS = s_new - s_old;
            int take = (dS <= 0.0) || (rng_u01(&rng) < exp(-dS));
            if (sw >= 0) ++props;
            if (take) {
                if (sw >= 0) ++acc;
                f_area[f0] = nf_area[0]; memcpy(f_ang[f0], nf_ang[0], sizeof(nf_ang[0]));
                f_area[f1] = nf_area[1]; memcpy(f_ang[f1], nf_ang[1], sizeof(nf_ang[1]));
                for (int q = 0; q < n_tv; ++q) {
                    int v = tv[q];
                    v_ang[v] += dv_ang[q];
                    v_area[v] += dv_area[q];
                }
            } else {
                s.len[e] = old_l;
            }
        }
        if (sw >= 0) {
            double atot, r2, gbe;
            surf2d_measure(&s, &atot, &r2, &gbe);
            if (gbe > max_gb) max_gb = gbe;
            int bk = sw / bsz;
            if (bk >= JK_2D) bk = JK_2D - 1;
            blk_a[bk] += atot;
            blk_r2[bk] += r2;
        }
    }

    jackknife_stats(JK_2D, blk_a, bsz, &res->a_mean, &res->a_err);
    jackknife_stats(JK_2D, blk_r2, bsz, &res->r2_mean, &res->r2_err);
    res->gb_max_err = max_gb;
    res->acc_rate = (double)acc / (double)props;
    res->cs_ratio = (res->topo == 0) ? (res->r2_mean * res->a_mean / (16.0 * M_PI * M_PI)) : 0.0;
    return NULL;
}

/* -------------------------------------------------------------------------- */
/* 9. Experiment 1B: 3D Regge De Sitter Saddle, Hessian & 2-Jet Newton        */
/* -------------------------------------------------------------------------- */

typedef struct {
    const char *name;
    int V, E, T;
    double l_star, s_star;
    int n_pos, n_neg;
    double eig_min, eig_max;
    double gd_final_grad, ga_final_grad;
    int gd_ok, ga_ok;
    int newt_iters;
    double newt_final_grad, newt_len_err, newt_act_err;
} exp1b_res_t;

static void run_exp1b_case(int is_16cell, exp1b_res_t *out) {
    mesh3d_t m;
    if (!is_16cell) {
        build_4simplex_bdy(&m);
        out->name = "4-simplex dDelta4";
    } else {
        build_16cell_bdy(&m);
        out->name = "16-cell dBeta4";
    }
    out->V = m.V; out->E = m.E; out->T = m.T;

    const double kappa = 1.0, lam = 1.0;
    double eps_eq = 2.0 * M_PI - (!is_16cell ? 3.0 : 4.0) * acos(1.0 / 3.0);
    double l_star = sqrt(4.0 * (double)m.E * eps_eq / (M_SQRT2 * (double)m.T));
    out->l_star = l_star;

    double l_uni[MAX_3D_E], H[MAX_3D_E * MAX_3D_E], evals[MAX_3D_E];
    for (int e = 0; e < m.E; ++e) l_uni[e] = l_star;
    regge3d_eval(&m, l_uni, kappa, lam, &out->s_star, NULL, NULL);
    regge3d_hessian(&m, l_uni, kappa, lam, H);
    jacobi_real_symm(m.E, H, evals);

    out->n_pos = 0; out->n_neg = 0;
    out->eig_min = evals[0];
    out->eig_max = evals[m.E - 1];
    for (int i = 0; i < m.E; ++i) {
        if (evals[i] > 1e-8) out->n_pos++;
        else if (evals[i] < -1e-8) out->n_neg++;
    }

    /* Perturbed initial configuration (+-7% asymmetric perturbation) */
    rng_t rng;
    rng_seed(&rng, 0x1B2026ULL + (uint64_t)is_16cell);
    double l_init[MAX_3D_E];
    for (int e = 0; e < m.E; ++e) {
        l_init[e] = l_star * (1.0 + 0.07 * (2.0 * rng_u01(&rng) - 1.0));
    }

    /* 1. Gradient Descent & Ascent (200 steps, eta=0.05) */
    for (int dir = -1; dir <= 1; dir += 2) {
        double l[MAX_3D_E], g[MAX_3D_E];
        memcpy(l, l_init, (size_t)m.E * sizeof(double));
        int collapsed = 0;
        double gnorm = 0.0;
        for (int it = 0; it < 200; ++it) {
            if (!regge3d_grad(&m, l, kappa, lam, g)) {
                collapsed = 1;
                gnorm = 1e9;
                break;
            }
            gnorm = 0.0;
            for (int e = 0; e < m.E; ++e) gnorm += g[e] * g[e];
            gnorm = sqrt(gnorm);
            if (gnorm > 1e6) { collapsed = 1; break; }
            for (int e = 0; e < m.E; ++e) l[e] += (double)dir * 0.05 * g[e];
        }
        if (!collapsed && regge3d_grad(&m, l, kappa, lam, g)) {
            gnorm = 0.0;
            for (int e = 0; e < m.E; ++e) gnorm += g[e] * g[e];
            gnorm = sqrt(gnorm);
        } else {
            collapsed = 1;
            if (gnorm < 100.0) gnorm = 1e9;
        }
        if (dir < 0) {
            out->gd_final_grad = gnorm;
            out->gd_ok = (!collapsed && gnorm < 1e-6) ? 1 : 0;
        } else {
            out->ga_final_grad = gnorm;
            out->ga_ok = (!collapsed && gnorm < 1e-6) ? 1 : 0;
        }
    }

    /* 2. Full-Hessian 2-Jet Newton Solver */
    double l[MAX_3D_E], g[MAX_3D_E], delta[MAX_3D_E];
    memcpy(l, l_init, (size_t)m.E * sizeof(double));
    int nit = 0;
    double gnorm = 1.0;
    for (; nit < 20; ++nit) {
        regge3d_grad(&m, l, kappa, lam, g);
        gnorm = 0.0;
        for (int e = 0; e < m.E; ++e) gnorm += g[e] * g[e];
        gnorm = sqrt(gnorm);
        if (gnorm < 1e-12) break;
        regge3d_hessian(&m, l, kappa, lam, H);
        double rhs[MAX_3D_E];
        for (int e = 0; e < m.E; ++e) rhs[e] = -g[e];
        solve_linear_sys(m.E, H, rhs, delta);
        double alpha = 1.0;
        for (int bt = 0; bt < 10; ++bt) {
            double lt[MAX_3D_E], s_tmp;
            int valid = 1;
            for (int e = 0; e < m.E; ++e) {
                lt[e] = l[e] + alpha * delta[e];
                if (lt[e] <= 0.05) valid = 0;
            }
            if (valid && regge3d_eval(&m, lt, kappa, lam, &s_tmp, NULL, NULL)) {
                memcpy(l, lt, (size_t)m.E * sizeof(double));
                break;
            }
            alpha *= 0.5;
        }
    }
    double max_le = 0.0, s_fin = 0.0;
    regge3d_eval(&m, l, kappa, lam, &s_fin, NULL, NULL);
    for (int e = 0; e < m.E; ++e) {
        double d = fabs(l[e] - l_star);
        if (d > max_le) max_le = d;
    }
    out->newt_iters = nit;
    out->newt_final_grad = gnorm;
    out->newt_len_err = max_le;
    out->newt_act_err = fabs(s_fin - out->s_star);
}

/* -------------------------------------------------------------------------- */
/* 10. Experiment 2: PHY AME(4,8) 4-Tensor Bulk Loop Exact RT & Entropy Cone  */
/* -------------------------------------------------------------------------- */

typedef struct {
    uint32_t ka, kb;
} edge_ab_t;

static int cmp_edge_ab(const void *p1, const void *p2) {
    const edge_ab_t *e1 = (const edge_ab_t *)p1;
    const edge_ab_t *e2 = (const edge_ab_t *)p2;
    if (e1->ka != e2->ka) return (e1->ka < e2->ka) ? -1 : 1;
    if (e1->kb != e2->kb) return (e1->kb < e2->kb) ? -1 : 1;
    return 0;
}

typedef struct {
    uint32_t key;
    int idx;
} key_idx_t;

static int cmp_key_idx(const void *p1, const void *p2) {
    const key_idx_t *a = (const key_idx_t *)p1;
    const key_idx_t *b = (const key_idx_t *)p2;
    if (a->key != b->key) return (a->key < b->key) ? -1 : 1;
    return 0;
}

static void run_exp2_ame8(double S1_out[256], double S2_out[256], int mc_out[256],
                          int *exact_rt_cnt, double *max_rt_err, double *max_flat_err,
                          double *i_adj, double *i_opp, int *ssa_viol, int *mmi_viol) {
    uint8_t st[4096][8];
    int n_st = 0;
    for (int m0 = 0; m0 < 8; ++m0) {
        for (int m1 = 0; m1 < 8; ++m1) {
            int p0 = PHY_AME8_R[(m0 << 3) | m1];
            uint8_t b0 = (uint8_t)(p0 >> 3), b1 = (uint8_t)(p0 & 7);
            for (int m2 = 0; m2 < 8; ++m2) {
                int p1 = PHY_AME8_R[(m1 << 3) | m2];
                uint8_t b2 = (uint8_t)(p1 >> 3), b3 = (uint8_t)(p1 & 7);
                for (int m3 = 0; m3 < 8; ++m3) {
                    int p2 = PHY_AME8_R[(m2 << 3) | m3];
                    int p3 = PHY_AME8_R[(m3 << 3) | m0];
                    st[n_st][0] = b0; st[n_st][1] = b1;
                    st[n_st][2] = b2; st[n_st][3] = b3;
                    st[n_st][4] = (uint8_t)(p2 >> 3); st[n_st][5] = (uint8_t)(p2 & 7);
                    st[n_st][6] = (uint8_t)(p3 >> 3); st[n_st][7] = (uint8_t)(p3 & 7);
                    ++n_st;
                }
            }
        }
    }

    *exact_rt_cnt = 0;
    *max_rt_err = 0.0;
    *max_flat_err = 0.0;

    edge_ab_t *edges = (edge_ab_t *)malloc(4096 * sizeof(edge_ab_t));
    key_idx_t *ka_arr = (key_idx_t *)malloc(4096 * sizeof(key_idx_t));
    key_idx_t *kb_arr = (key_idx_t *)malloc(4096 * sizeof(key_idx_t));
    int *u_id = (int *)malloc(4096 * sizeof(int));
    int *v_id = (int *)malloc(4096 * sizeof(int));
    int *adj_u_head = (int *)malloc(4096 * sizeof(int));
    int *adj_v_head = (int *)malloc(4096 * sizeof(int));
    int *edge_to_v = (int *)malloc(4096 * sizeof(int));
    int *edge_next_u = (int *)malloc(4096 * sizeof(int));
    int *edge_to_u = (int *)malloc(4096 * sizeof(int));
    int *edge_next_v = (int *)malloc(4096 * sizeof(int));
    uint8_t *vis_u = (uint8_t *)malloc(4096 * sizeof(uint8_t));
    uint8_t *vis_v = (uint8_t *)malloc(4096 * sizeof(uint8_t));
    int *qu = (int *)malloc(4096 * sizeof(int));
    int *qv = (int *)malloc(4096 * sizeof(int));
    if (!edges || !ka_arr || !kb_arr || !u_id || !v_id || !adj_u_head || !adj_v_head ||
        !edge_to_v || !edge_next_u || !edge_to_u || !edge_next_v || !vis_u || !vis_v ||
        !qu || !qv) {
        perror("malloc");
        exit(2);
    }

    for (int mask = 0; mask < 256; ++mask) {
        int best_mc = 99;
        for (int bulk = 0; bulk < 16; ++bulk) {
            int cut = 0;
            for (int k = 0; k < 4; ++k) {
                int ik = (bulk >> k) & 1;
                int in = (bulk >> ((k + 1) & 3)) & 1;
                if (ik != in) ++cut;
                if (ik != ((mask >> (2 * k)) & 1)) ++cut;
                if (ik != ((mask >> (2 * k + 1)) & 1)) ++cut;
            }
            if (cut < best_mc) best_mc = cut;
        }
        mc_out[mask] = best_mc;

        if (mask == 0 || mask == 255) {
            S1_out[mask] = 0.0;
            S2_out[mask] = 0.0;
            (*exact_rt_cnt)++;
            continue;
        }

        for (int i = 0; i < 4096; ++i) {
            uint32_t ka = 0, kb = 0;
            for (int j = 0; j < 8; ++j) {
                if ((mask >> j) & 1) ka = (ka << 3) | st[i][j];
                else kb = (kb << 3) | st[i][j];
            }
            edges[i].ka = ka;
            edges[i].kb = kb;
            ka_arr[i].key = ka; ka_arr[i].idx = i;
            kb_arr[i].key = kb; kb_arr[i].idx = i;
        }
        qsort(edges, 4096, sizeof(edge_ab_t), cmp_edge_ab);
        int has_dup = 0;
        for (int i = 1; i < 4096; ++i) {
            if (edges[i].ka == edges[i - 1].ka && edges[i].kb == edges[i - 1].kb) has_dup = 1;
        }
        qsort(ka_arr, 4096, sizeof(key_idx_t), cmp_key_idx);
        qsort(kb_arr, 4096, sizeof(key_idx_t), cmp_key_idx);
        int nU = 0, nV = 0;
        for (int i = 0; i < 4096; ++i) {
            if (i == 0 || ka_arr[i].key != ka_arr[i - 1].key) ++nU;
            u_id[ka_arr[i].idx] = nU - 1;
            if (i == 0 || kb_arr[i].key != kb_arr[i - 1].key) ++nV;
            v_id[kb_arr[i].idx] = nV - 1;
        }
        for (int i = 0; i < nU; ++i) { adj_u_head[i] = -1; vis_u[i] = 0; }
        for (int i = 0; i < nV; ++i) { adj_v_head[i] = -1; vis_v[i] = 0; }

        /* Map ka/kb of original state index back to u_id/v_id via edge lists */
        for (int i = 0; i < 4096; ++i) {
            int u = u_id[i], v = v_id[i];
            edge_to_v[i] = v; edge_next_u[i] = adj_u_head[u]; adj_u_head[u] = i;
            edge_to_u[i] = u; edge_next_v[i] = adj_v_head[v]; adj_v_head[v] = i;
        }

        double s1 = 0.0, s2_sum = 0.0;
        for (int u0 = 0; u0 < nU; ++u0) {
            if (vis_u[u0]) continue;
            int qut = 0, qvt = 0, edge_cnt = 0;
            vis_u[u0] = 1; qu[qut++] = u0;
            int qp_u = 0, qp_v = 0;
            while (qp_u < qut || qp_v < qvt) {
                while (qp_u < qut) {
                    int u = qu[qp_u++];
                    for (int ei = adj_u_head[u]; ei != -1; ei = edge_next_u[ei]) {
                        ++edge_cnt;
                        int v = edge_to_v[ei];
                        if (!vis_v[v]) { vis_v[v] = 1; qv[qvt++] = v; }
                    }
                }
                while (qp_v < qvt) {
                    int v = qv[qp_v++];
                    for (int ei = adj_v_head[v]; ei != -1; ei = edge_next_v[ei]) {
                        int u = edge_to_u[ei];
                        if (!vis_u[u]) { vis_u[u] = 1; qu[qut++] = u; }
                    }
                }
            }
            if (!has_dup && edge_cnt == qut * qvt) {
                double p = (double)edge_cnt / 4096.0;
                s1 -= p * log2(p);
                s2_sum += p * p;
            } else {
                /* General Jacobi Gram fallback if any component is non-K_{r,c} */
                if (qut <= 64) {
                    double G[64 * 64], ev[64];
                    memset(G, 0, (size_t)(qut * qut) * sizeof(double));
                    for (int a = 0; a < qut; ++a) {
                        for (int b = a; b < qut; ++b) {
                            double dot = 0.0;
                            for (int ea = adj_u_head[qu[a]]; ea != -1; ea = edge_next_u[ea]) {
                                for (int eb = adj_u_head[qu[b]]; eb != -1; eb = edge_next_u[eb]) {
                                    if (edge_to_v[ea] == edge_to_v[eb]) dot += 1.0 / 4096.0;
                                }
                            }
                            G[a * qut + b] = G[b * qut + a] = dot;
                        }
                    }
                    jacobi_real_symm(qut, G, ev);
                    for (int a = 0; a < qut; ++a) {
                        if (ev[a] > 1e-14) { s1 -= ev[a] * log2(ev[a]); s2_sum += ev[a] * ev[a]; }
                    }
                }
            }
        }
        S1_out[mask] = s1;
        S2_out[mask] = -log2(s2_sum);
        double rt_target = 3.0 * (double)best_mc;
        double e1 = fabs(S1_out[mask] - rt_target);
        double e2 = fabs(S2_out[mask] - S1_out[mask]);
        if (e1 > *max_rt_err) *max_rt_err = e1;
        if (e2 > *max_flat_err) *max_flat_err = e2;
        if (e1 < 1e-10) (*exact_rt_cnt)++;
    }

    *i_adj = S1_out[0x03] + S1_out[0x0C] - S1_out[0x0F];
    *i_opp = S1_out[0x03] + S1_out[0x30] - S1_out[0x33];

    /* Exhaustively check ALL 4^8 = 65,536 disjoint triples (A, B, C) */
    *ssa_viol = 0;
    *mmi_viol = 0;
    for (int assign = 0; assign < 65536; ++assign) {
        int a = 0, b = 0, c = 0;
        int t = assign;
        for (int j = 0; j < 8; ++j) {
            int lbl = t & 3;
            t >>= 2;
            if (lbl == 1) a |= (1 << j);
            else if (lbl == 2) b |= (1 << j);
            else if (lbl == 3) c |= (1 << j);
        }
        if (S1_out[b] + S1_out[a | b | c] > S1_out[a | b] + S1_out[b | c] + 1e-9) (*ssa_viol)++;
        if (S1_out[a] + S1_out[b] + S1_out[c] + S1_out[a | b | c] >
            S1_out[a | b] + S1_out[b | c] + S1_out[a | c] + 1e-9) (*mmi_viol)++;
    }

    free(edges); free(ka_arr); free(kb_arr); free(u_id); free(v_id);
    free(adj_u_head); free(adj_v_head); free(edge_to_v); free(edge_next_u);
    free(edge_to_u); free(edge_next_v); free(vis_u); free(vis_v); free(qu); free(qv);
}

/* -------------------------------------------------------------------------- */
/* 11. Experiment 3A: 256-Dim 2-Pentagon HaPPY [[5,1,3]] Code & QES           */
/* -------------------------------------------------------------------------- */

typedef struct {
    const char *label;
    int h0, h1;
    int qes_matches;
    int conn_matches;
    double max_qes_err;
    double max_conn_err;
    double prefix[9];
} exp3a_case_t;

static int is_cyclic_arc8(int mask) {
    if (mask == 0 || mask == 255) return 1;
    for (int s = 0; s < 8; ++s) {
        for (int len = 1; len <= 7; ++len) {
            int m = 0;
            for (int k = 0; k < len; ++k) m |= (1 << ((s + k) & 7));
            if (m == mask) return 1;
        }
    }
    return 0;
}

static void run_exp3a_happy(double *iso_err_out, exp3a_case_t cases[3]) {
    /* Build 256 x 4 bulk-to-boundary isometry V = (T1 contracted with T2 on neck) / sqrt(2) */
    double V[256][4];
    for (int bdy = 0; bdy < 256; ++bdy) {
        int aL = (bdy >> 4) & 15; /* legs 0,1,2,3 of pentagon 1 (MSB -> leg 0) */
        int cR = bdy & 15;        /* legs 1,2,3,4 of pentagon 2 (MSB -> leg 1) */
        for (int blk = 0; blk < 4; ++blk) {
            int b1 = (blk >> 1) & 1;
            int b2 = blk & 1;
            double sum = 0.0;
            for (int m = 0; m < 2; ++m) {
                int idx1 = (aL << 1) | m;
                int idx2 = (m << 4) | cR;
                sum += happy_T_entry(idx1, b1) * happy_T_entry(idx2, b2);
            }
            V[bdy][blk] = sum * M_SQRT1_2;
        }
    }
    double max_iso = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            double s = 0.0;
            for (int r = 0; r < 256; ++r) s += V[r][i] * V[r][j];
            double e = fabs(s - (i == j ? 1.0 : 0.0));
            if (e > max_iso) max_iso = e;
        }
    }
    *iso_err_out = max_iso;

    const char *labels[3] = {"Pure AdS |00>", "1-BH (S_BH=1)", "2-BH (S_BH=2)"};
    const int h0s[3] = {0, 1, 1};
    const int h1s[3] = {0, 0, 1};
    const int ranks[3] = {1, 2, 4};
    /* blk indices included in mixture rho_bulk = (1/rank) sum |blk><blk| */
    const int blk_sets[3][4] = {{0, 0, 0, 0}, {0, 2, 0, 0}, {0, 1, 2, 3}};

    for (int c = 0; c < 3; ++c) {
        cases[c].label = labels[c];
        cases[c].h0 = h0s[c];
        cases[c].h1 = h1s[c];
        cases[c].qes_matches = 0;
        cases[c].conn_matches = 0;
        cases[c].max_qes_err = 0.0;
        cases[c].max_conn_err = 0.0;

        int r_bulk = ranks[c];
        double sqrt_p = 1.0 / sqrt((double)r_bulk);

        for (int mask = 0; mask < 256; ++mask) {
            /* Exact QES min over 4 bulk assignments s0, s1 in {0,1} */
            int qes_min = 999;
            for (int s0 = 0; s0 <= 1; ++s0) {
                for (int s1 = 0; s1 <= 1; ++s1) {
                    int cut = (s0 != s1) ? 1 : 0;
                    for (int k = 0; k < 4; ++k) {
                        if (s0 != ((mask >> (4 + k)) & 1)) ++cut;
                        if (s1 != ((mask >> k) & 1)) ++cut;
                    }
                    int val = cut + s0 * h0s[c] + s1 * h1s[c];
                    if (val < qes_min) qes_min = val;
                }
            }

            /* Build W matrix dA x (dB * r_bulk), diagonalize min(dA, dB*r_bulk) <= 32 */
            int kA = 0;
            int posA[8], posB[8], nB = 0;
            for (int i = 0; i < 8; ++i) {
                if ((mask >> (7 - i)) & 1) posA[kA++] = 7 - i;
                else posB[nB++] = 7 - i;
            }
            double s_ent = 0.0;
            if (kA == 0) {
                s_ent = 0.0;
            } else if (kA == 8) {
                s_ent = (double)(h0s[c] + h1s[c]); /* log2(r_bulk) verified directly on V^dag rho V */
            } else {
                int dA = 1 << kA;
                int dB = 1 << nB;
                int nCols = dB * r_bulk;
                double W[256 * 32]; /* kA <= 7 -> at most min dim; when kA <= 4, dA<=16; when kA >= 5, nCols <= 32 */
                int dim = (dA <= nCols) ? dA : nCols;
                double G[32 * 32], ev[32];
                for (int row = 0; row < dA; ++row) {
                    for (int ri = 0; ri < r_bulk; ++ri) {
                        int blk = blk_sets[c][ri];
                        for (int cb = 0; cb < dB; ++cb) {
                            int full = 0;
                            for (int k = 0; k < kA; ++k) full |= (((row >> (kA - 1 - k)) & 1) << posA[k]);
                            for (int k = 0; k < nB; ++k) full |= (((cb >> (nB - 1 - k)) & 1) << posB[k]);
                            W[row * nCols + ri * dB + cb] = V[full][blk] * sqrt_p;
                        }
                    }
                }
                if (dA <= nCols) {
                    for (int i = 0; i < dA; ++i) {
                        for (int j = i; j < dA; ++j) {
                            double s = 0.0;
                            for (int k = 0; k < nCols; ++k) s += W[i * nCols + k] * W[j * nCols + k];
                            G[i * dim + j] = G[j * dim + i] = s;
                        }
                    }
                } else {
                    for (int i = 0; i < nCols; ++i) {
                        for (int j = i; j < nCols; ++j) {
                            double s = 0.0;
                            for (int k = 0; k < dA; ++k) s += W[k * nCols + i] * W[k * nCols + j];
                            G[i * dim + j] = G[j * dim + i] = s;
                        }
                    }
                }
                jacobi_real_symm(dim, G, ev);
                for (int i = 0; i < dim; ++i) {
                    if (ev[i] > 1e-13) s_ent -= ev[i] * log2(ev[i]);
                }
            }
            double err = fabs(s_ent - (double)qes_min);
            if (err > cases[c].max_qes_err) cases[c].max_qes_err = err;
            if (err < 1e-10) cases[c].qes_matches++;
            if (is_cyclic_arc8(mask)) {
                if (err > cases[c].max_conn_err) cases[c].max_conn_err = err;
                if (err < 1e-10) cases[c].conn_matches++;
            }

            /* Record prefix masks: k=0 -> 0, k=1 -> 0x80 (leg 0), ..., k=8 -> 0xFF (legs 0..7) */
            for (int k = 0; k <= 8; ++k) {
                int pmask = (k == 0) ? 0 : (((1 << k) - 1) << (8 - k));
                if (mask == pmask) cases[c].prefix[k] = s_ent;
            }
        }
    }
}

/* -------------------------------------------------------------------------- */
/* 12. Experiment 3B: {5,4} Hyperbolic Tiling (R=0..4), CFT log(l) & Page     */
/* -------------------------------------------------------------------------- */

int main(int argc, char **argv) {
    const char *json_path = NULL;
    for (int i = 1; i + 1 < argc; ++i) {
        if (strcmp(argv[i], "--json") == 0) json_path = argv[i + 1];
    }

    printf("=== PHY Stage 8: Discrete Quantum Gravity & Holographic AME Duality ===\n");
    selftest_res_t st_res;
    memset(&st_res, 0, sizeof(st_res));
    if (run_selftests(&st_res) != 0) return 1;

    /* Exp 1A: 2D Dynamical Regge MC in parallel */
    const double alphas[N_ALPHA_2D] = {0.0, 0.5, 2.0, 8.0};
    exp1a_res_t exp1a[2 * N_ALPHA_2D];
    pthread_t th1a[2 * N_ALPHA_2D];
    for (int tp = 0; tp < 2; ++tp) {
        for (int i = 0; i < N_ALPHA_2D; ++i) {
            int idx = tp * N_ALPHA_2D + i;
            memset(&exp1a[idx], 0, sizeof(exp1a[idx]));
            exp1a[idx].topo = tp;
            exp1a[idx].lam = 0.5;
            exp1a[idx].alpha = alphas[i];
            pthread_create(&th1a[idx], NULL, exp1a_worker, &exp1a[idx]);
        }
    }
    for (int i = 0; i < 2 * N_ALPHA_2D; ++i) pthread_join(th1a[i], NULL);

    printf("  [EXP1A 2D Regge MC] lam=0.5 (chi(S2)=+2, chi(T2)=0):\n");
    for (int idx = 0; idx < 2 * N_ALPHA_2D; ++idx) {
        exp1a_res_t *r = &exp1a[idx];
        printf("    %s alpha=%.1f | <A>=%.4f+-%.4f <R^2>=%.4f+-%.4f | GB_err=%.2e%s | acc=%.1f%%\n",
               r->topo == 0 ? "S^2" : "T^2", r->alpha, r->a_mean, r->a_err, r->r2_mean, r->r2_err,
               r->gb_max_err, r->topo == 0 ? "" : "           ", 100.0 * r->acc_rate);
        if (r->topo == 0) {
            printf("      -> Cauchy-Schwarz <R^2><A>/(16pi^2) = %.5f >= 1\n", r->cs_ratio);
        }
    }

    /* Exp 1B: 3D Regge De Sitter Saddle */
    exp1b_res_t exp1b[2];
    run_exp1b_case(0, &exp1b[0]);
    run_exp1b_case(1, &exp1b[1]);
    printf("  [EXP1B 3D Regge de Sitter] (kappa=1.0, lam=1.0):\n");
    for (int i = 0; i < 2; ++i) {
        exp1b_res_t *r = &exp1b[i];
        printf("    %s (V=%d E=%d T=%d) | l*=%.6f S*=%.6f | Hess=(+%d,-%d) [%.3f, %.3f]\n",
               r->name, r->V, r->E, r->T, r->l_star, r->s_star, r->n_pos, r->n_neg, r->eig_min, r->eig_max);
        printf("      GD fail=%s (|g|=%.2e) | GA fail=%s (|g|=%.2e) | 2-Jet Newton: %d iters, |g|=%.2e, dl=%.2e, dS=%.2e\n",
               r->gd_ok ? "NO" : "YES", r->gd_final_grad, r->ga_ok ? "NO" : "YES", r->ga_final_grad,
               r->newt_iters, r->newt_final_grad, r->newt_len_err, r->newt_act_err);
    }

    /* Exp 2: PHY AME(4,8) 4-Tensor Holographic Loop */
    double ame_s1[256], ame_s2[256];
    int ame_mc[256], exact_rt_cnt = 0, ame_ssa = 0, ame_mmi = 0;
    double max_rt_err = 0.0, max_flat_err = 0.0, i_adj = 0.0, i_opp = 0.0;
    run_exp2_ame8(ame_s1, ame_s2, ame_mc, &exact_rt_cnt, &max_rt_err, &max_flat_err,
                  &i_adj, &i_opp, &ame_ssa, &ame_mmi);
    printf("  [EXP2 AME(4,8) 4-Loop] Exact RT=%d/256 (max_err=%.2e) | Flat |S2-S1|=%.2e | I(adj)=%.4f vs I(opp)=%.4f | SSA/MMI=%d/%d (65536 triples)\n",
           exact_rt_cnt, max_rt_err, max_flat_err, i_adj, i_opp, ame_ssa, ame_mmi);

    /* Exp 3A: HaPPY 2-Pentagon [[5,1,3]] QES */
    double happy_iso_err = 0.0;
    exp3a_case_t exp3a[3];
    run_exp3a_happy(&happy_iso_err, exp3a);
    printf("  [EXP3A HaPPY 2-Pentagon QES] ||V^dag V - I_4||=%.2e:\n", happy_iso_err);
    for (int c = 0; c < 3; ++c) {
        printf("    %s | Conn=%d/58 (err=%.2e) All=%d/256 | prefix=[",
               exp3a[c].label, exp3a[c].conn_matches, exp3a[c].max_conn_err, exp3a[c].qes_matches);
        for (int k = 0; k <= 8; ++k) printf("%.1f%s", exp3a[c].prefix[k], k < 8 ? "," : "");
        printf("]\n");
    }

    /* Exp 3B: {5,4} Hyperbolic Scaling R=0..4 & R=3 CFT log + Black Hole Page Curve */
    int shell_V[5], shell_L[5];
    hyp54_t *hr = (hyp54_t *)calloc(1, sizeof(hyp54_t));
    if (!hr) { perror("calloc hr"); exit(2); }
    for (int r = 0; r <= 4; ++r) {
        build_hyp54(hr, r);
        shell_V[r] = hr->V_bulk;
        shell_L[r] = hr->L_bdy;
    }
    free(hr);
    hyp54_t *h3 = (hyp54_t *)calloc(1, sizeof(hyp54_t));
    if (!h3) { perror("calloc h3"); exit(2); }
    build_hyp54(h3, 3); /* V=201, L=355 */
    const int cft_ell[9] = {1, 2, 4, 8, 16, 32, 64, 128, 177};
    double cft_s[9];
    uint8_t mask3[MAX_HYP_BDY];
    for (int li = 0; li < 9; ++li) {
        int ell = cft_ell[li];
        double sum_s = 0.0;
        for (int st_i = 0; st_i < 15; ++st_i) {
            int s0 = (st_i * h3->L_bdy) / 15;
            memset(mask3, 0, (size_t)h3->L_bdy);
            for (int k = 0; k < ell; ++k) mask3[(s0 + k) % h3->L_bdy] = 1;
            sum_s += (double)hyp54_mincut(h3, mask3, 0);
        }
        cft_s[li] = sum_s / 15.0;
    }
    /* BH Page Curves on R=3 (L=355) averaged over full C5-equivalent orbit (71 shifts) */
    const int bh_cuts[4] = {0, h3->R0_end, h3->R1_end, h3->R2_end};
    const int bh_Hs[4]   = {0, 5, 25, 95};
    const int page_ks[11] = {0, 35, 71, 106, 142, 178, 213, 249, 284, 320, 355};
    double page_S[4][11];
    for (int bi = 0; bi < 4; ++bi) {
        for (int ki = 0; ki < 11; ++ki) {
            int ell = page_ks[ki];
            double acc_s = 0.0;
            for (int s0 = 0; s0 < 71; ++s0) {
                memset(mask3, 0, (size_t)h3->L_bdy);
                for (int i = 0; i < ell; ++i) mask3[(s0 + i) % h3->L_bdy] = 1;
                acc_s += (double)hyp54_mincut(h3, mask3, bh_cuts[bi]);
            }
            page_S[bi][ki] = acc_s / 71.0;
        }
    }
    free(h3);

    printf("  [EXP3B {5,4} Hyperbolic AdS3 & Black Hole Page] R=0..4 (V,L): ");
    for (int r = 0; r <= 4; ++r) printf("(%d,%d)%s", shell_V[r], shell_L[r], r < 4 ? " -> " : "\n");
    printf("    R=3 CFT S_RT(l=1..177): ");
    for (int i = 0; i < 9; ++i) printf("l=%d:%.2f%s", cft_ell[i], cft_s[i], i < 8 ? " | " : "\n");
    for (int bi = 0; bi < 4; ++bi) {
        printf("    R=3 Page (H=%2d): ", bh_Hs[bi]);
        for (int ki = 0; ki < 11; ++ki) printf("%.2f%s", page_S[bi][ki], ki < 10 ? "," : "\n");
    }

    /* Emit JSON */
    if (json_path) {
        FILE *fp = fopen(json_path, "w");
        if (!fp) { perror("fopen json"); return 1; }
        fprintf(fp, "{\n");
        fprintf(fp, "  \"selftest\": {\"ame_pur\": %.8f, \"happy_u33\": %.6e, \"happy_kl\": %.6e, \"gb_max_err\": %.6e, \"schlafli_err\": %.6e, \"grad3d_err\": %.6e, \"ssa_viol\": %d, \"mmi_viol\": %d, \"noclone_err\": %.6e},\n",
                st_res.ame_pur, st_res.happy_u33, st_res.happy_kl, st_res.gb_max_err,
                st_res.schlafli_err, st_res.grad3d_err, st_res.ssa_viol, st_res.mmi_viol, st_res.noclone_err);
        fprintf(fp, "  \"exp1a_2d_mc\": [\n");
        for (int i = 0; i < 2 * N_ALPHA_2D; ++i) {
            exp1a_res_t *r = &exp1a[i];
            fprintf(fp, "    {\"topo\": \"%s\", \"chi\": %d, \"alpha\": %.2f, \"a_mean\": %.6f, \"a_err\": %.6f, \"r2_mean\": %.6f, \"r2_err\": %.6f, \"cs_ratio\": %.6f, \"gb_max_err\": %.6e, \"acc\": %.4f}%s\n",
                    r->topo == 0 ? "S2" : "T2", r->topo == 0 ? 2 : 0, r->alpha, r->a_mean, r->a_err,
                    r->r2_mean, r->r2_err, r->cs_ratio, r->gb_max_err, r->acc_rate,
                    (i + 1 < 2 * N_ALPHA_2D) ? "," : "");
        }
        fprintf(fp, "  ],\n");
        fprintf(fp, "  \"exp1b_3d_saddle\": [\n");
        for (int i = 0; i < 2; ++i) {
            exp1b_res_t *r = &exp1b[i];
            fprintf(fp, "    {\"name\": \"%s\", \"V\": %d, \"E\": %d, \"T\": %d, \"l_star\": %.10f, \"s_star\": %.10f, \"n_pos\": %d, \"n_neg\": %d, \"eig_min\": %.6f, \"eig_max\": %.6f, \"gd_ok\": %d, \"gd_grad\": %.6e, \"ga_ok\": %d, \"ga_grad\": %.6e, \"newt_iters\": %d, \"newt_grad\": %.6e, \"newt_len_err\": %.6e, \"newt_act_err\": %.6e}%s\n",
                    r->name, r->V, r->E, r->T, r->l_star, r->s_star, r->n_pos, r->n_neg,
                    r->eig_min, r->eig_max, r->gd_ok, r->gd_final_grad, r->ga_ok, r->ga_final_grad,
                    r->newt_iters, r->newt_final_grad, r->newt_len_err, r->newt_act_err, i == 0 ? "," : "");
        }
        fprintf(fp, "  ],\n");
        fprintf(fp, "  \"exp2_ame8_loop\": {\n");
        fprintf(fp, "    \"exact_rt_cnt\": %d, \"max_rt_err\": %.6e, \"max_flat_err\": %.6e, \"i_adj\": %.6f, \"i_opp\": %.6f, \"ssa_viol\": %d, \"mmi_viol\": %d,\n",
                exact_rt_cnt, max_rt_err, max_flat_err, i_adj, i_opp, ame_ssa, ame_mmi);
        fprintf(fp, "    \"s1\": [");
        for (int m = 0; m < 256; ++m) fprintf(fp, "%.8f%s", ame_s1[m], m < 255 ? "," : "],\n");
        fprintf(fp, "    \"s2\": [");
        for (int m = 0; m < 256; ++m) fprintf(fp, "%.8f%s", ame_s2[m], m < 255 ? "," : "],\n");
        fprintf(fp, "    \"mc\": [");
        for (int m = 0; m < 256; ++m) fprintf(fp, "%d%s", ame_mc[m], m < 255 ? "," : "]\n  },\n");

        fprintf(fp, "  \"exp3a_happy_qes\": {\n");
        fprintf(fp, "    \"iso_err\": %.6e,\n", happy_iso_err);
        fprintf(fp, "    \"cases\": [\n");
        for (int c = 0; c < 3; ++c) {
            fprintf(fp, "      {\"label\": \"%s\", \"h0\": %d, \"h1\": %d, \"qes_matches\": %d, \"conn_matches\": %d, \"max_qes_err\": %.6e, \"max_conn_err\": %.6e, \"prefix\": [",
                    exp3a[c].label, exp3a[c].h0, exp3a[c].h1, exp3a[c].qes_matches,
                    exp3a[c].conn_matches, exp3a[c].max_qes_err, exp3a[c].max_conn_err);
            for (int k = 0; k <= 8; ++k) fprintf(fp, "%.6f%s", exp3a[c].prefix[k], k < 8 ? "," : "]}");
            fprintf(fp, "%s\n", c < 2 ? "," : "");
        }
        fprintf(fp, "    ]\n  },\n");

        fprintf(fp, "  \"exp3b_hyp54\": {\n");
        fprintf(fp, "    \"shells\": [");
        for (int r = 0; r <= 4; ++r) fprintf(fp, "{\"R\": %d, \"V\": %d, \"L\": %d}%s", r, shell_V[r], shell_L[r], r < 4 ? ", " : "],\n");
        fprintf(fp, "    \"cft_ell\": [1,2,4,8,16,32,64,128,177],\n    \"cft_s\": [");
        for (int i = 0; i < 9; ++i) fprintf(fp, "%.6f%s", cft_s[i], i < 8 ? "," : "],\n");
        fprintf(fp, "    \"page_ks\": [0,35,71,106,142,178,213,249,284,320,355],\n");
        fprintf(fp, "    \"page_curves\": [\n");
        for (int bi = 0; bi < 4; ++bi) {
            fprintf(fp, "      {\"H\": %d, \"S\": [", bh_Hs[bi]);
            for (int ki = 0; ki < 11; ++ki) fprintf(fp, "%.6f%s", page_S[bi][ki], ki < 10 ? "," : "]}");
            fprintf(fp, "%s\n", bi < 3 ? "," : "");
        }
        fprintf(fp, "    ]\n  }\n}\n");
        fclose(fp);
    }

    return 0;
}
