/*
 * PHY Stage 9: Unified Physics — SM Anomaly Quantization, GUT SU(5)/SO(10) RG Unification,
 * 4D SU(2)-Higgs Electroweak Fradkin-Shenker Phase Transition & Coupled 2D Regge-Maxwell Monopoles
 * (Flow.UnifiedBundle Theorems U1-U5)
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
/* 1. PRNG, Quaternion SU(2) Algebra & Helpers                                */
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

static inline double rng_pos_u01(rng_t *r) {
    double u = rng_u01(r);
    return (u > 1e-16) ? u : 1e-16;
}

typedef struct {
    double w, x, y, z;
} quat_t;

static inline quat_t q_mk(double w, double x, double y, double z) {
    quat_t q = {w, x, y, z};
    return q;
}

static inline quat_t q_conj(quat_t a) {
    return q_mk(a.w, -a.x, -a.y, -a.z);
}

static inline quat_t q_mul(quat_t a, quat_t b) {
    return q_mk(
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w
    );
}

static inline quat_t q_add(quat_t a, quat_t b) {
    return q_mk(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline quat_t q_scale(double s, quat_t a) {
    return q_mk(s * a.w, s * a.x, s * a.y, s * a.z);
}

static inline double q_norm2(quat_t a) {
    return a.w * a.w + a.x * a.x + a.y * a.y + a.z * a.z;
}

static inline double q_dot(quat_t a, quat_t b) {
    return a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline quat_t q_norm(quat_t a) {
    double inv = 1.0 / sqrt(q_norm2(a));
    return q_scale(inv, a);
}

static quat_t q_rand_su2(rng_t *rng) {
    for (;;) {
        double u1 = 2.0 * rng_u01(rng) - 1.0;
        double u2 = 2.0 * rng_u01(rng) - 1.0;
        double s1 = u1 * u1 + u2 * u2;
        if (s1 >= 1.0) continue;
        double u3 = 2.0 * rng_u01(rng) - 1.0;
        double u4 = 2.0 * rng_u01(rng) - 1.0;
        double s2 = u3 * u3 + u4 * u4;
        if (s2 >= 1.0) continue;
        double sq = sqrt((1.0 - s1) / (s2 > 1e-16 ? s2 : 1e-16));
        return q_mk(u1, u2, u3 * sq, u4 * sq);
    }
}

/* Exact Kennedy-Pendleton heatbath: sample U ~ exp(q_dot(U, S_tot)) dU */
static quat_t kp_heatbath(rng_t *rng, quat_t S_tot) {
    double k = sqrt(q_norm2(S_tot));
    if (k < 1e-12) return q_rand_su2(rng);
    quat_t v = q_scale(1.0 / k, S_tot);
    double a0 = 0.0;
    for (;;) {
        double r1 = rng_pos_u01(rng);
        double r2 = rng_u01(rng);
        double r3 = rng_pos_u01(rng);
        double u4 = rng_u01(rng);
        double c = cos(2.0 * M_PI * r2);
        double lam2 = -(log(r1) + c * c * log(r3)) / (2.0 * k);
        if (u4 * u4 <= 1.0 - lam2) {
            a0 = 1.0 - 2.0 * lam2;
            break;
        }
    }
    double rad = sqrt((1.0 - a0 * a0 > 0.0) ? (1.0 - a0 * a0) : 0.0);
    double rx, ry, rz;
    for (;;) {
        rx = 2.0 * rng_u01(rng) - 1.0;
        ry = 2.0 * rng_u01(rng) - 1.0;
        rz = 2.0 * rng_u01(rng) - 1.0;
        double n2 = rx * rx + ry * ry + rz * rz;
        if (n2 > 1e-12 && n2 <= 1.0) {
            double sc = rad / sqrt(n2);
            rx *= sc; ry *= sc; rz *= sc;
            break;
        }
    }
    /* Wait: S_tot weight is exp(q_dot(U, S_tot)) = exp(k * q_dot(U, v)) = exp(k * Re(U * v_conj))! */
    /*Wait: Re(A * v) = q_dot(A, conj(v)), so for Re(W) = a0 with weight exp(k*Re(W)), W = U * conj(v) => U = W * v!*/
    quat_t W = q_mk(a0, rx, ry, rz);
    return q_norm(q_mul(W, v));
}

static double clamp_cos(double x) {
    if (x > 1.0) return 1.0;
    if (x < -1.0) return -1.0;
    return x;
}

static double wrap_angle(double a) {
    while (a > M_PI) a -= 2.0 * M_PI;
    while (a <= -M_PI) a += 2.0 * M_PI;
    return a;
}

static void jackknife_stats(int nb, const double *blk, int bsz, double *mean_out, double *err_out) {
    double tot = 0.0;
    for (int b = 0; b < nb; ++b) tot += blk[b];
    double mean = tot / ((double)nb * (double)bsz);
    double jv = 0.0;
    double jk_n = (double)(nb - 1) * (double)bsz;
    for (int b = 0; b < nb; ++b) {
        double d = (tot - blk[b]) / jk_n - mean;
        jv += d * d;
    }
    *mean_out = mean;
    *err_out = sqrt(((double)(nb - 1) / (double)nb) * jv);
}

/* -------------------------------------------------------------------------- */
/* 2. Standard Model / MSSM 1-Loop & 2-Loop RGE Integrator                    */
/* -------------------------------------------------------------------------- */

static void rg_2loop_rhs(const double b1[3], const double B2[3][3], const double y[3], double rhs[3]) {
    double inv_8pi2 = 1.0 / (8.0 * M_PI * M_PI);
    double inv_2pi  = 1.0 / (2.0 * M_PI);
    for (int i = 0; i < 3; ++i) {
        double corr = 0.0;
        for (int j = 0; j < 3; ++j) corr += B2[i][j] / y[j];
        rhs[i] = b1[i] * inv_2pi - corr * inv_8pi2;
    }
}

static void rg_rk4(const double y0[3], const double b1[3], const double B2[3][3],
                   double t_span, int steps, double y_out[3]) {
    double y[3] = {y0[0], y0[1], y0[2]};
    double dt = t_span / (double)steps;
    for (int s = 0; s < steps; ++s) {
        double k1[3], k2[3], k3[3], k4[3], yt[3];
        rg_2loop_rhs(b1, B2, y, k1);
        for (int i = 0; i < 3; ++i) yt[i] = y[i] + 0.5 * dt * k1[i];
        rg_2loop_rhs(b1, B2, yt, k2);
        for (int i = 0; i < 3; ++i) yt[i] = y[i] + 0.5 * dt * k2[i];
        rg_2loop_rhs(b1, B2, yt, k3);
        for (int i = 0; i < 3; ++i) yt[i] = y[i] + dt * k3[i];
        rg_2loop_rhs(b1, B2, yt, k4);
        for (int i = 0; i < 3; ++i) y[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    memcpy(y_out, y, 3 * sizeof(double));
}

/* -------------------------------------------------------------------------- */
/* 3. 4D Lattice (6^4) SU(2)-Higgs Topology & Helpers                         */
/* -------------------------------------------------------------------------- */

enum { L4 = 6, V4 = L4 * L4 * L4 * L4, JK_EW = 20 };

static inline int idx4(int t, int x, int y, int z) {
    return (((t * L4 + x) * L4 + y) * L4 + z);
}

static inline int shift4_idx(int s, int dir, int step) {
    int z = s % L4; int tmp = s / L4;
    int y = tmp % L4; tmp /= L4;
    int x = tmp % L4;
    int t = tmp / L4;
    if (dir == 0) t = (t + step + L4) % L4;
    else if (dir == 1) x = (x + step + L4) % L4;
    else if (dir == 2) y = (y + step + L4) % L4;
    else z = (z + step + L4) % L4;
    return idx4(t, x, y, z);
}

static quat_t su2_gauge_staple(const quat_t U[V4][4], int s, int mu) {
    quat_t st = q_mk(0, 0, 0, 0);
    int spm = shift4_idx(s, mu, +1);
    for (int nu = 0; nu < 4; ++nu) {
        if (nu == mu) continue;
        int spn = shift4_idx(s, nu, +1);
        /* Forward staple: U_nu(s) * U_mu(s+nu) * U_nu†(s+mu) -> U * st_fwd† has P_{mu,nu}? Wait: */
        /* Wilson plaq P = U_mu(s) * U_nu(s+mu) * U_mu†(s+nu) * U_nu†(s); Trace = qDot(U_mu(s), V)! */
        /* Here P = U_mu(s) * conj(V) => V = U_nu(s) * U_mu(s+nu) * U_nu†(s+mu)! */
        quat_t fwd = q_mul(q_mul(U[s][nu], U[spn][mu]), q_conj(U[spm][nu]));
        st = q_add(st, fwd);
        /* Backward staple: U_nu†(s-nu) * U_mu(s-nu) * U_nu(s+mu-nu) */
        int smn = shift4_idx(s, nu, -1);
        int spmn = shift4_idx(spm, nu, -1);
        quat_t bwd = q_mul(q_mul(q_conj(U[smn][nu]), U[smn][mu]), U[spmn][nu]);
        st = q_add(st, bwd);
    }
    return st;
}

static quat_t su2_higgs_action_eval(const quat_t U[V4][4], const quat_t h[V4],
                                    double beta, double kappa, double lam, double *act_out) {
    double sp = 0.0, sh = 0.0, sv = 0.0;
    for (int s = 0; s < V4; ++s) {
        double r2 = q_norm2(h[s]);
        double d1 = r2 - 1.0;
        sv += r2 + lam * d1 * d1;
        for (int mu = 0; mu < 4; ++mu) {
            int spm = shift4_idx(s, mu, +1);
            sh += q_dot(h[s], q_mul(U[s][mu], h[spm]));
            for (int nu = mu + 1; nu < 4; ++nu) {
                int spn = shift4_idx(s, nu, +1);
                quat_t p = q_mul(q_mul(U[s][mu], U[spm][nu]),
                                 q_mul(q_conj(U[spn][mu]), q_conj(U[s][nu])));
                sp += p.w;
            }
        }
    }
    *act_out = -beta * sp - 2.0 * kappa * sh + sv;
    return q_mk(sp / (6.0 * (double)V4), sh / (4.0 * (double)V4), sv / (double)V4, 0.0);
}

/* -------------------------------------------------------------------------- */
/* 4. Oriented Icosahedron S^2 for Regge-Maxwell + Monopole Backreaction       */
/* -------------------------------------------------------------------------- */

enum { ICO_V = 12, ICO_E = 30, ICO_F = 20 };

typedef struct {
    double coords[ICO_V][3];
    int eu[ICO_E], ev[ICO_E];
    int fv[ICO_F][3];
    int fe[ICO_F][3];  /* opposite vertex edge indices e0=(v1,v2), e1=(v0,v2), e2=(v0,v1) */
    int bf[ICO_F][3];  /* oriented boundary edges e_k in CCW order (v0->v1, v1->v2, v2->v0) */
    int bsgn[ICO_F][3];/* +1 if eu->ev matches CCW orientation, -1 if opposite */
    int ef_cnt[ICO_E], ef_list[ICO_E][2];
    double len[ICO_E];
    double th[ICO_E];
} ico_s2_t;

static void build_oriented_icosahedron(ico_s2_t *m) {
    memset(m, 0, sizeof(*m));
    const double phi = 0.5 * (1.0 + sqrt(5.0));
    double p[ICO_V][3] = {
        {0, 1, phi}, {0, -1, phi}, {0, 1, -phi}, {0, -1, -phi},
        {1, phi, 0}, {-1, phi, 0}, {1, -phi, 0}, {-1, -phi, 0},
        {phi, 0, 1}, {-phi, 0, 1}, {phi, 0, -1}, {-phi, 0, -1}
    };
    memcpy(m->coords, p, sizeof(p));

    int e_idx = 0;
    for (int a = 0; a < ICO_V; ++a) {
        for (int b = a + 1; b < ICO_V; ++b) {
            double d2 = 0.0;
            for (int k = 0; k < 3; ++k) d2 += (p[a][k] - p[b][k]) * (p[a][k] - p[b][k]);
            if (d2 <= 4.2) {
                m->eu[e_idx] = a;
                m->ev[e_idx] = b;
                m->len[e_idx] = 1.0;
                ++e_idx;
            }
        }
    }
    /* Hack cleanly: reuse exact field names without relying on m->E before setting */
    int f_idx = 0;
    for (int a = 0; a < ICO_V; ++a) {
        for (int b = a + 1; b < ICO_V; ++b) {
            double dab = 0.0;
            for (int k = 0; k < 3; ++k) dab += (p[a][k]-p[b][k])*(p[a][k]-p[b][k]);
            if (dab > 4.2) continue;
            for (int c = b + 1; c < ICO_V; ++c) {
                double dac = 0.0, dbc = 0.0;
                for (int k = 0; k < 3; ++k) {
                    dac += (p[a][k]-p[c][k])*(p[a][k]-p[c][k]);
                    dbc += (p[b][k]-p[c][k])*(p[b][k]-p[c][k]);
                }
                if (dac <= 4.2 && dbc <= 4.2) {
                    int v0 = a, v1 = b, v2 = c;
                    /* Orient CCW outward using dot((p1-p0)x(p2-p0), p0+p1+p2) > 0 */
                    double u[3], v[3], n[3], cent[3];
                    for (int k = 0; k < 3; ++k) {
                        u[k] = p[v1][k] - p[v0][k];
                        v[k] = p[v2][k] - p[v0][k];
                        cent[k] = p[v0][k] + p[v1][k] + p[v2][k];
                    }
                    n[0] = u[1]*v[2] - u[2]*v[1];
                    n[1] = u[2]*v[0] - u[0]*v[2];
                    n[2] = u[0]*v[1] - u[1]*v[0];
                    if (n[0]*cent[0] + n[1]*cent[1] + n[2]*cent[2] < 0.0) {
                        int tmp = v1; v1 = v2; v2 = tmp;
                    }
                    m->fv[f_idx][0] = v0;
                    m->fv[f_idx][1] = v1;
                    m->fv[f_idx][2] = v2;
                    ++f_idx;
                }
            }
        }
    }
    /* Match edge indices & orientations */
    for (int f = 0; f < ICO_F; ++f) {
        int v[3] = {m->fv[f][0], m->fv[f][1], m->fv[f][2]};
        /* opposite edges fe[f][0]=(v1,v2), fe[f][1]=(v0,v2), fe[f][2]=(v0,v1) */
        int opp[3][2] = {{v[1], v[2]}, {v[0], v[2]}, {v[0], v[1]}};
        for (int k = 0; k < 3; ++k) {
            int x = opp[k][0], y = opp[k][1];
            if (x > y) { int t = x; x = y; y = t; }
            for (int e = 0; e < ICO_E; ++e) {
                if (m->eu[e] == x && m->ev[e] == y) {
                    m->fe[f][k] = e;
                    m->ef_list[e][m->ef_cnt[e]++] = f;
                }
            }
        }
        /* CCW directed boundary edges: (v0->v1), (v1->v2), (v2->v0) */
        int bpair[3][2] = {{v[0], v[1]}, {v[1], v[2]}, {v[2], v[0]}};
        for (int k = 0; k < 3; ++k) {
            int p0 = bpair[k][0], p1 = bpair[k][1];
            int lo = p0 < p1 ? p0 : p1;
            int hi = p0 < p1 ? p1 : p0;
            for (int e = 0; e < ICO_E; ++e) {
                if (m->eu[e] == lo && m->ev[e] == hi) {
                    m->bf[f][k] = e;
                    m->bsgn[f][k] = (p0 < p1) ? +1 : -1;
                }
            }
        }
    }
}

static int ico_tri_valid(double a, double b, double c, double lmin, double lmax) {
    if (a < lmin || b < lmin || c < lmin || a > lmax || b > lmax || c > lmax) return 0;
    return (a + b - c > 1e-5) && (a + c - b > 1e-5) && (b + c - a > 1e-5);
}

static void ico_tri_geom(const ico_s2_t *m, int f, double *area_out, double ang[3]) {
    double l0 = m->len[m->fe[f][0]];
    double l1 = m->len[m->fe[f][1]];
    double l2 = m->len[m->fe[f][2]];
    double p = 0.5 * (l0 + l1 + l2);
    double arg = p * (p - l0) * (p - l1) * (p - l2);
    *area_out = (arg > 0.0) ? sqrt(arg) : 0.0;
    ang[0] = acos(clamp_cos((l1*l1 + l2*l2 - l0*l0) / (2.0 * l1 * l2)));
    ang[1] = acos(clamp_cos((l0*l0 + l2*l2 - l1*l1) / (2.0 * l0 * l2)));
    ang[2] = acos(clamp_cos((l0*l0 + l1*l1 - l2*l2) / (2.0 * l0 * l1)));
}

/* -------------------------------------------------------------------------- */
/* 5. Startup Self-Tests (6 Exact Algebraic, Uniqueness & Orientation Checks) */
/* -------------------------------------------------------------------------- */

typedef struct {
    int anom_all_zero;
    int uniq_total, uniq_bad;
    double gut_sin2w_err;
    double so4_norm_err;
    double ew_mgam_err, ew_rho_err;
    double su2_higgs_gi_err;
    double rg_rk4_err;
    double ico_bianchi_err, ico_q1_err, ico_q2_err;
} selftest_res_t;

static int run_selftests(selftest_res_t *st) {
    memset(st, 0, sizeof(*st));

    /* 1. SM Anomaly Cancellation & Brute-Force Hypercharge Uniqueness Audit */
    int a3 = 2*1 + (-4) + 2;
    int a2 = 3*1 + (-3);
    int ag = 6*1 + 3*(-4) + 3*2 + 2*(-3) + 6;
    int ac = 6*1 + 3*(-64) + 3*8 + 2*(-27) + 216;
    st->anom_all_zero = (a3 == 0 && a2 == 0 && ag == 0 && ac == 0 && ((3 + 1) % 2 == 0)) ? 1 : 0;

    for (int yQ = 1; yQ <= 3; ++yQ) {
        for (int yu = -12; yu <= 12; ++yu) {
            int yd = -(2 * yQ + yu);
            if (yd < -12 || yd > 12) continue;
            int yL = -3 * yQ;
            int ye = -(6 * yQ + 3 * yu + 3 * yd + 2 * yL);
            if (ye < -24 || ye > 24) continue;
            long long cub = 6LL*yQ*yQ*yQ + 3LL*yu*yu*yu + 3LL*yd*yd*yd + 2LL*yL*yL*yL + (long long)ye*ye*ye;
            if (cub == 0) {
                st->uniq_total++;
                int ok_u = ((yu == -4*yQ && yd == 2*yQ) || (yu == 2*yQ && yd == -4*yQ)) && (ye == 6*yQ);
                if (!ok_u) st->uniq_bad++;
            }
        }
    }
    if (!st->anom_all_zero || st->uniq_bad != 0 || st->uniq_total != 6) {
        fprintf(stderr, "SELFTEST 1 FAILED: anom=%d tot=%d bad=%d\n",
                st->anom_all_zero, st->uniq_total, st->uniq_bad);
        return -1;
    }

    /* 2. SU(5) sin^2(theta_W)=3/8 & SO(4) ~= SU(2)_L x SU(2)_R quaternion vierbein norm */
    const double q5[5] = {2, 2, 2, -6, 0};
    const double t3[5] = {0, 0, 0, -3, 3};
    double n_t = 0.0, d_q = 0.0;
    for (int i = 0; i < 5; ++i) { n_t += t3[i] * t3[i]; d_q += q5[i] * q5[i]; }
    st->gut_sin2w_err = fabs(n_t / d_q - 3.0 / 8.0);

    rng_t rng;
    rng_seed(&rng, 0x9001ULL);
    double max_so4 = 0.0;
    for (int i = 0; i < 100; ++i) {
        quat_t qL = q_rand_su2(&rng), qR = q_rand_su2(&rng);
        quat_t x = q_mk(2.0*rng_u01(&rng)-1.0, 2.0*rng_u01(&rng)-1.0,
                        2.0*rng_u01(&rng)-1.0, 2.0*rng_u01(&rng)-1.0);
        quat_t xp = q_mul(q_mul(qL, x), q_conj(qR));
        double e = fabs(q_norm2(xp) - q_norm2(x));
        if (e > max_so4) max_so4 = e;
    }
    st->so4_norm_err = max_so4;

    /* 3. Electroweak Neutral Gauge Mass Matrix & Custodial rho = 1 */
    const double g = 0.6517, gp = 0.3574;
    double g2 = g * g, gp2 = gp * gp, tr = g2 + gp2;
    double det = g2 * gp2 - (-g * gp) * (-g * gp);
    double lam_min = 0.5 * (tr - sqrt(tr * tr - 4.0 * det));
    double lam_max = 0.5 * (tr + sqrt(tr * tr - 4.0 * det));
    st->ew_mgam_err = fabs(lam_min);
    double mw2 = g2, cos2w = g2 / lam_max;
    st->ew_rho_err = fabs(mw2 / (lam_max * cos2w) - 1.0);

    /* 4. 4D 6^4 SU(2)-Higgs Exact Local Gauge Invariance */
    quat_t (*U)[4] = calloc(V4, sizeof(*U));
    quat_t (*Up)[4] = calloc(V4, sizeof(*Up));
    quat_t *h = calloc(V4, sizeof(*h));
    quat_t *hp = calloc(V4, sizeof(*hp));
    quat_t *gs = calloc(V4, sizeof(*gs));
    for (int s = 0; s < V4; ++s) {
        gs[s] = q_rand_su2(&rng);
        h[s] = q_mk(0.7+0.3*rng_u01(&rng), rng_u01(&rng)-0.5, rng_u01(&rng)-0.5, rng_u01(&rng)-0.5);
        for (int mu = 0; mu < 4; ++mu) U[s][mu] = q_rand_su2(&rng);
    }
    for (int s = 0; s < V4; ++s) {
        hp[s] = q_mul(gs[s], h[s]);
        for (int mu = 0; mu < 4; ++mu) {
            int spm = shift4_idx(s, mu, +1);
            Up[s][mu] = q_mul(q_mul(gs[s], U[s][mu]), q_conj(gs[spm]));
        }
    }
    double act0 = 0.0, act1 = 0.0;
    su2_higgs_action_eval(U, h, 2.3, 0.30, 0.5, &act0);
    su2_higgs_action_eval(Up, hp, 2.3, 0.30, 0.5, &act1);
    st->su2_higgs_gi_err = fabs(act1 - act0);
    free(U); free(Up); free(h); free(hp); free(gs);

    /* 5. 1-Loop RK4 vs Exact Closed-Form RG Agreement */
    const double b1_sm[3] = {-41.0 / 10.0, 19.0 / 6.0, 7.0};
    const double B0[3][3] = {{0}};
    double a0[3] = {59.0197, 29.5848, 8.4746}, a_rk[3];
    double t16 = log(1e16 / 91.1876);
    rg_rk4(a0, b1_sm, B0, t16, 2000, a_rk);
    double max_rg = 0.0;
    for (int i = 0; i < 3; ++i) {
        double ex = a0[i] + (b1_sm[i] / (2.0 * M_PI)) * t16;
        double e = fabs(a_rk[i] - ex);
        if (e > max_rg) max_rg = e;
    }
    st->rg_rk4_err = max_rg;

    /* 6. Oriented Icosahedron S^2 Exact Edge Orientation, Bianchi & Dirac Monopoles */
    ico_s2_t ico;
    build_oriented_icosahedron(&ico);
    int esum[ICO_E];
    memset(esum, 0, sizeof(esum));
    for (int f = 0; f < ICO_F; ++f) {
        for (int k = 0; k < 3; ++k) esum[ico.bf[f][k]] += ico.bsgn[f][k];
    }
    for (int e = 0; e < ICO_E; ++e) {
        if (esum[e] != 0) {
            fprintf(stderr, "SELFTEST 6 ORIENTATION FAILED at edge %d (sum=%d)\n", e, esum[e]);
            return -1;
        }
        ico.th[e] = 0.4 * (2.0 * rng_u01(&rng) - 1.0);
    }
    double sum_raw = 0.0, sum_q1 = 0.0, sum_q2 = 0.0;
    for (int f = 0; f < ICO_F; ++f) {
        double dth = 0.0;
        for (int k = 0; k < 3; ++k) dth += (double)ico.bsgn[f][k] * ico.th[ico.bf[f][k]];
        sum_raw += dth;
        sum_q1 += wrap_angle(dth + 2.0 * M_PI / (double)ICO_F);
        sum_q2 += wrap_angle(dth + 4.0 * M_PI / (double)ICO_F);
    }
    st->ico_bianchi_err = fabs(sum_raw);
    st->ico_q1_err = fabs(sum_q1 / (2.0 * M_PI) - 1.0);
    st->ico_q2_err = fabs(sum_q2 / (2.0 * M_PI) - 2.0);

    if (st->gut_sin2w_err > 1e-14 || st->so4_norm_err > 1e-13 ||
        st->ew_mgam_err > 1e-13 || st->ew_rho_err > 1e-13 ||
        st->su2_higgs_gi_err > 1e-10 || st->rg_rk4_err > 1e-11 ||
        st->ico_bianchi_err > 1e-13 || st->ico_q1_err > 1e-13 || st->ico_q2_err > 1e-13) {
        fprintf(stderr, "SELFTEST RESIDUAL FAILED\n");
        return -1;
    }

    printf("  [SELFTEST] PASS (Anom=4/4 uniq=%d bad=%d | SU5=3/8 err=%.2e | SO4=%.2e | EW m_gam=%.2e rho_err=%.2e | SU2HiggsGI=%.2e | RG=%.2e | S2 Bianchi=%.2e Q1=%.2e Q2=%.2e)\n",
           st->uniq_total, st->uniq_bad, st->gut_sin2w_err, st->so4_norm_err,
           st->ew_mgam_err, st->ew_rho_err, st->su2_higgs_gi_err, st->rg_rk4_err,
           st->ico_bianchi_err, st->ico_q1_err, st->ico_q2_err);
    return 0;
}

/* -------------------------------------------------------------------------- */
/* 6. Experiment 1: 4D SU(2)-Higgs Electroweak Transition (6^4 Lattice)        */
/* -------------------------------------------------------------------------- */

enum { N_KAPPA = 5 };

typedef struct {
    double beta, kappa, lam;
    double p_mean, p_err;
    double l_mean, l_err;
    double rho2_mean, rho2_err;
    double chi22_mean, chi22_err;
    double higgs_acc;
} exp1_ew_res_t;

static quat_t su2_link_staple_sm(const quat_t U[V4][4], int s, int mu) {
    quat_t st = q_mk(0, 0, 0, 0);
    int spm = shift4_idx(s, mu, +1);
    for (int nu = 0; nu < 4; ++nu) {
        if (nu == mu) continue;
        int spn = shift4_idx(s, nu, +1);
        quat_t fwd = q_mul(q_mul(U[s][nu], U[spn][mu]), q_conj(U[spm][nu]));
        st = q_add(st, fwd);
        int smn = shift4_idx(s, nu, -1);
        int spmn = shift4_idx(spm, nu, -1);
        quat_t bwd = q_mul(q_mul(q_conj(U[smn][nu]), U[smn][mu]), U[spmn][nu]);
        st = q_add(st, bwd);
    }
    return st;
}

static quat_t su2_wilson_loop(const quat_t U[V4][4], int s, int mu, int nu, int l1, int l2) {
    quat_t p = q_mk(1, 0, 0, 0);
    int cur = s;
    for (int i = 0; i < l1; ++i) { p = q_mul(p, U[cur][mu]); cur = shift4_idx(cur, mu, +1); }
    for (int i = 0; i < l2; ++i) { p = q_mul(p, U[cur][nu]); cur = shift4_idx(cur, nu, +1); }
    for (int i = 0; i < l1; ++i) { cur = shift4_idx(cur, mu, -1); p = q_mul(p, q_conj(U[cur][mu])); }
    for (int i = 0; i < l2; ++i) { cur = shift4_idx(cur, nu, -1); p = q_mul(p, q_conj(U[cur][nu])); }
    return p;
}

static void ape_smear_4d(const quat_t U_in[V4][4], quat_t U_out[V4][4], double eps) {
    for (int s = 0; s < V4; ++s) {
        for (int mu = 0; mu < 4; ++mu) {
            quat_t st = su2_link_staple_sm(U_in, s, mu);
            quat_t comb = q_add(U_in[s][mu], q_scale(eps, st));
            U_out[s][mu] = q_norm(comb);
        }
    }
}

static void *exp1_ew_worker(void *arg) {
    exp1_ew_res_t *res = (exp1_ew_res_t *)arg;
    quat_t (*U)[4]  = calloc(V4, sizeof(*U));
    quat_t (*Us1)[4] = calloc(V4, sizeof(*Us1));
    quat_t (*Us2)[4] = calloc(V4, sizeof(*Us2));
    quat_t *h       = calloc(V4, sizeof(*h));

    rng_t rng;
    rng_seed(&rng, 0xE1000000ULL + (uint64_t)(res->kappa * 1000.0));

    for (int s = 0; s < V4; ++s) {
        h[s] = q_mk(1.0, 0.0, 0.0, 0.0);
        for (int mu = 0; mu < 4; ++mu) U[s][mu] = q_mk(1.0, 0.0, 0.0, 0.0);
    }

    const int n_therm = 250;
    const int n_meas = 500;
    const int bsz = n_meas / JK_EW;
    double blk_p[JK_EW], blk_l[JK_EW], blk_r2[JK_EW], blk_w11[JK_EW], blk_w21[JK_EW], blk_w22[JK_EW];
    memset(blk_p, 0, sizeof(blk_p));
    memset(blk_l, 0, sizeof(blk_l));
    memset(blk_r2, 0, sizeof(blk_r2));
    memset(blk_w11, 0, sizeof(blk_w11));
    memset(blk_w21, 0, sizeof(blk_w21));
    memset(blk_w22, 0, sizeof(blk_w22));

    uint64_t h_prop = 0, h_acc = 0;
    const double h_step = 0.45;

    for (int sw = -n_therm; sw < n_meas; ++sw) {
        /* 1. Exact Kennedy-Pendleton SU(2) heatbath + 1 Adler OR sweep over gauge links */
        for (int s = 0; s < V4; ++s) {
            for (int mu = 0; mu < 4; ++mu) {
                int spm = shift4_idx(s, mu, +1);
                quat_t st_g = su2_gauge_staple(U, s, mu);
                quat_t st_h = q_mul(h[s], q_conj(h[spm]));
                quat_t st_tot = q_add(q_scale(res->beta, st_g), q_scale(2.0 * res->kappa, st_h));
                U[s][mu] = kp_heatbath(&rng, st_tot);
                /* Over-relaxation reflection across normalized staple */
                double n2 = q_norm2(st_tot);
                if (n2 > 1e-12) {
                    quat_t v = q_scale(1.0 / sqrt(n2), st_tot);
                    double c = q_dot(U[s][mu], v);
                    quat_t uor = q_add(q_scale(2.0 * c, v), q_scale(-1.0, U[s][mu]));
                    U[s][mu] = q_norm(uor);
                }
            }
        }
        /* 2. Cartesian 4-vector Higgs Metropolis (5 hits per site) */
        for (int s = 0; s < V4; ++s) {
            quat_t V_h = q_mk(0, 0, 0, 0);
            for (int mu = 0; mu < 4; ++mu) {
                int spm = shift4_idx(s, mu, +1);
                int smm = shift4_idx(s, mu, -1);
                V_h = q_add(V_h, q_mul(U[s][mu], h[spm]));
                V_h = q_add(V_h, q_mul(q_conj(U[smm][mu]), h[smm]));
            }
            double r2 = q_norm2(h[s]);
            double s_cur = r2 + res->lam * (r2 - 1.0) * (r2 - 1.0) - 2.0 * res->kappa * q_dot(h[s], V_h);
            for (int hit = 0; hit < 5; ++hit) {
                quat_t cand = q_mk(
                    h[s].w + h_step * (2.0 * rng_u01(&rng) - 1.0),
                    h[s].x + h_step * (2.0 * rng_u01(&rng) - 1.0),
                    h[s].y + h_step * (2.0 * rng_u01(&rng) - 1.0),
                    h[s].z + h_step * (2.0 * rng_u01(&rng) - 1.0)
                );
                double cr2 = q_norm2(cand);
                double s_new = cr2 + res->lam * (cr2 - 1.0) * (cr2 - 1.0) - 2.0 * res->kappa * q_dot(cand, V_h);
                double ds = s_new - s_cur;
                if (sw >= 0) ++h_prop;
                if (ds <= 0.0 || rng_pos_u01(&rng) < exp(-ds)) {
                    h[s] = cand;
                    s_cur = s_new;
                    if (sw >= 0) ++h_acc;
                }
            }
        }

        if (sw >= 0) {
            double act_dum = 0.0;
            quat_t obs = su2_higgs_action_eval(U, h, res->beta, res->kappa, res->lam, &act_dum);
            /* Normalized Higgs link order L_H = <qDot(h_x, U h_{x+mu}) / (|h_x| |h_{x+mu}|)> */
            double lh = 0.0, r2sum = 0.0;
            for (int s = 0; s < V4; ++s) {
                double nr0 = sqrt(q_norm2(h[s]));
                r2sum += nr0 * nr0;
                for (int mu = 0; mu < 4; ++mu) {
                    int spm = shift4_idx(s, mu, +1);
                    double nr1 = sqrt(q_norm2(h[spm]));
                    lh += q_dot(h[s], q_mul(U[s][mu], h[spm])) / (nr0 * nr1);
                }
            }
            lh /= (4.0 * (double)V4);
            r2sum /= (double)V4;

            /* 2-step 4D APE smeared Wilson loops for clean Creutz chi(2,2) */
            ape_smear_4d(U, Us1, 0.45);
            ape_smear_4d(Us1, Us2, 0.45);
            double w11 = 0.0, w21 = 0.0, w22 = 0.0;
            for (int s = 0; s < V4; ++s) {
                for (int mu = 0; mu < 4; ++mu) {
                    for (int nu = mu + 1; nu < 4; ++nu) {
                        w11 += su2_wilson_loop(Us2, s, mu, nu, 1, 1).w;
                        w21 += su2_wilson_loop(Us2, s, mu, nu, 2, 1).w;
                        w22 += su2_wilson_loop(Us2, s, mu, nu, 2, 2).w;
                    }
                }
            }
            double np = 6.0 * (double)V4;
            int bk = sw / bsz;
            if (bk >= JK_EW) bk = JK_EW - 1;
            blk_p[bk]   += obs.w;
            blk_l[bk]   += lh;
            blk_r2[bk]  += r2sum;
            blk_w11[bk] += w11 / np;
            blk_w21[bk] += w21 / np;
            blk_w22[bk] += w22 / np;
        }
    }

    jackknife_stats(JK_EW, blk_p, bsz, &res->p_mean, &res->p_err);
    jackknife_stats(JK_EW, blk_l, bsz, &res->l_mean, &res->l_err);
    jackknife_stats(JK_EW, blk_r2, bsz, &res->rho2_mean, &res->rho2_err);

    double t11 = 0.0, t21 = 0.0, t22 = 0.0;
    for (int b = 0; b < JK_EW; ++b) { t11 += blk_w11[b]; t21 += blk_w21[b]; t22 += blk_w22[b]; }
    double c_full = -log((t22 * t11) / (t21 * t21));
    double jv = 0.0;
    for (int b = 0; b < JK_EW; ++b) {
        double j11 = t11 - blk_w11[b], j21 = t21 - blk_w21[b], j22 = t22 - blk_w22[b];
        double d = -log((j22 * j11) / (j21 * j21)) - c_full;
        jv += d * d;
    }
    res->chi22_mean = c_full;
    res->chi22_err  = sqrt(((double)(JK_EW - 1) / (double)JK_EW) * jv);
    res->higgs_acc  = (double)h_acc / (double)h_prop;

    free(U); free(Us1); free(Us2); free(h);
    return NULL;
}

/* -------------------------------------------------------------------------- */
/* 7. Experiment 3: Coupled 2D Regge-Maxwell + Monopole Backreaction (S^2)    */
/* -------------------------------------------------------------------------- */

enum { N_EM_CASES = 4, JK_EM = 20 };

typedef struct {
    int q_mono;
    double lam, alpha, beta_em;
    double a_mean, a_err;
    double r2_mean, r2_err;
    double eem_mean, eem_err;
    double gb_max_err, bianchi_max_err;
    double acc_l, acc_th;
} exp3_em_res_t;

static void ico_measure_all(const ico_s2_t *m, int q_mono, double *atot, double *r2,
                            double *eem, double *gb_err, double *bianchi_err) {
    double v_ang[ICO_V], v_area[ICO_V];
    for (int v = 0; v < ICO_V; ++v) { v_ang[v] = 0.0; v_area[v] = 0.0; }
    double asum = 0.0, emsum = 0.0, qsum = 0.0;
    const double qbg = (2.0 * M_PI * (double)q_mono) / (double)ICO_F;

    for (int f = 0; f < ICO_F; ++f) {
        double af, ang[3];
        ico_tri_geom(m, f, &af, ang);
        asum += af;
        for (int k = 0; k < 3; ++k) {
            int v = m->fv[f][k];
            v_ang[v] += ang[k];
            v_area[v] += af / 3.0;
        }
        double dth = qbg;
        for (int k = 0; k < 3; ++k) dth += (double)m->bsgn[f][k] * m->th[m->bf[f][k]];
        double thf = wrap_angle(dth);
        qsum += thf;
        emsum += (1.0 - cos(thf)) / af;
    }
    double esum = 0.0, r2s = 0.0;
    for (int v = 0; v < ICO_V; ++v) {
        double eps = 2.0 * M_PI - v_ang[v];
        esum += eps;
        r2s += (eps * eps) / v_area[v];
    }
    *atot = asum;
    *r2 = r2s;
    *eem = emsum;
    *gb_err = fabs(esum / (2.0 * M_PI) - 2.0);
    *bianchi_err = fabs(qsum / (2.0 * M_PI) - (double)q_mono);
}

static void *exp3_em_worker(void *arg) {
    exp3_em_res_t *res = (exp3_em_res_t *)arg;
    ico_s2_t m;
    build_oriented_icosahedron(&m);
    const double qbg = (2.0 * M_PI * (double)res->q_mono) / (double)ICO_F;
    /* Start near classical Regge equilibrium */
    double l0 = sqrt(4.0 * M_PI * sqrt(res->alpha / res->lam) / (5.0 * sqrt(3.0)));
    for (int e = 0; e < ICO_E; ++e) { m.len[e] = l0; m.th[e] = 0.0; }

    rng_t rng;
    rng_seed(&rng, 0x3E000000ULL + (uint64_t)(res->q_mono + 1) * 1009ULL +
                   (res->beta_em > 0.0 ? 77ULL : 0ULL));

    const int n_therm = 2500;
    const int n_meas = 8000;
    const int bsz = n_meas / JK_EM;
    double blk_a[JK_EM], blk_r2[JK_EM], blk_em[JK_EM];
    memset(blk_a, 0, sizeof(blk_a));
    memset(blk_r2, 0, sizeof(blk_r2));
    memset(blk_em, 0, sizeof(blk_em));

    double v_ang[ICO_V], v_area[ICO_V], f_area[ICO_F], f_ang[ICO_F][3], f_th[ICO_F];
    for (int v = 0; v < ICO_V; ++v) { v_ang[v] = 0.0; v_area[v] = 0.0; }
    for (int f = 0; f < ICO_F; ++f) {
        ico_tri_geom(&m, f, &f_area[f], f_ang[f]);
        for (int k = 0; k < 3; ++k) {
            v_ang[m.fv[f][k]] += f_ang[f][k];
            v_area[m.fv[f][k]] += f_area[f] / 3.0;
        }
        f_th[f] = qbg;
    }

    const double step_l = 0.45, step_th = 0.40, lmin = 0.4, lmax = 8.0;
    uint64_t lp = 0, la = 0, tp = 0, ta = 0;
    double max_gb = 0.0, max_bi = 0.0;

    for (int sw = -n_therm; sw < n_meas; ++sw) {
        /* 1. Regge edge-length Metropolis coupled to electromagnetic Maxwell pressure */
        for (int e = 0; e < ICO_E; ++e) {
            double old_l = m.len[e];
            double new_l = old_l + step_l * (2.0 * rng_u01(&rng) - 1.0);
            int f0 = m.ef_list[e][0], f1 = m.ef_list[e][1];
            m.len[e] = new_l;
            if (!ico_tri_valid(m.len[m.fe[f0][0]], m.len[m.fe[f0][1]], m.len[m.fe[f0][2]], lmin, lmax) ||
                !ico_tri_valid(m.len[m.fe[f1][0]], m.len[m.fe[f1][1]], m.len[m.fe[f1][2]], lmin, lmax)) {
                m.len[e] = old_l;
                if (sw >= 0) ++lp;
                continue;
            }
            int tv[6], n_tv = 0;
            int fp[2] = {f0, f1};
            for (int fi = 0; fi < 2; ++fi) {
                for (int k = 0; k < 3; ++k) {
                    int v = m.fv[fp[fi]][k], seen = 0;
                    for (int q = 0; q < n_tv; ++q) if (tv[q] == v) seen = 1;
                    if (!seen) tv[n_tv++] = v;
                }
            }
            double nf_area[2], nf_ang[2][3];
            ico_tri_geom(&m, f0, &nf_area[0], nf_ang[0]);
            ico_tri_geom(&m, f1, &nf_area[1], nf_ang[1]);

            double s_old = res->lam * (f_area[f0] + f_area[f1]) +
                           res->beta_em * ((1.0 - cos(f_th[f0])) / f_area[f0] + (1.0 - cos(f_th[f1])) / f_area[f1]);
            double s_new = res->lam * (nf_area[0] + nf_area[1]) +
                           res->beta_em * ((1.0 - cos(f_th[f0])) / nf_area[0] + (1.0 - cos(f_th[f1])) / nf_area[1]);

            double dv_ang[6], dv_area[6];
            for (int q = 0; q < n_tv; ++q) {
                int v = tv[q];
                double eps_o = 2.0 * M_PI - v_ang[v];
                s_old += res->alpha * (eps_o * eps_o) / v_area[v];
                dv_ang[q] = 0.0; dv_area[q] = 0.0;
            }
            for (int fi = 0; fi < 2; ++fi) {
                int f = fp[fi];
                for (int k = 0; k < 3; ++k) {
                    int v = m.fv[f][k];
                    for (int q = 0; q < n_tv; ++q) if (tv[q] == v) {
                        dv_ang[q] += nf_ang[fi][k] - f_ang[f][k];
                        dv_area[q] += (nf_area[fi] - f_area[f]) / 3.0;
                    }
                }
            }
            for (int q = 0; q < n_tv; ++q) {
                int v = tv[q];
                double eps_n = 2.0 * M_PI - (v_ang[v] + dv_ang[q]);
                s_new += res->alpha * (eps_n * eps_n) / (v_area[v] + dv_area[q]);
            }
            double ds = s_new - s_old;
            if (sw >= 0) ++lp;
            if (ds <= 0.0 || rng_u01(&rng) < exp(-ds)) {
                if (sw >= 0) ++la;
                f_area[f0] = nf_area[0]; memcpy(f_ang[f0], nf_ang[0], sizeof(nf_ang[0]));
                f_area[f1] = nf_area[1]; memcpy(f_ang[f1], nf_ang[1], sizeof(nf_ang[1]));
                for (int q = 0; q < n_tv; ++q) {
                    v_ang[tv[q]] += dv_ang[q];
                    v_area[tv[q]] += dv_area[q];
                }
            } else {
                m.len[e] = old_l;
            }
        }

        /* 2. U(1) gauge link Metropolis (when beta_em > 0) within exact topological monopole sector Q */
        if (res->beta_em > 0.0) {
            for (int e = 0; e < ICO_E; ++e) {
                double old_t = m.th[e];
                double new_t = wrap_angle(old_t + step_th * (2.0 * rng_u01(&rng) - 1.0));
                int f0 = m.ef_list[e][0], f1 = m.ef_list[e][1];
                int s0 = 0, s1 = 0;
                for (int k = 0; k < 3; ++k) {
                    if (m.bf[f0][k] == e) s0 = m.bsgn[f0][k];
                    if (m.bf[f1][k] == e) s1 = m.bsgn[f1][k];
                }
                double nt0 = f_th[f0] + (double)s0 * (new_t - old_t);
                double nt1 = f_th[f1] + (double)s1 * (new_t - old_t);
                /* Admissibility corridor |theta_f| < pi - 0.1 preserves exact Chern monopole class Q */
                if (fabs(nt0) > M_PI - 0.1 || fabs(nt1) > M_PI - 0.1) {
                    if (sw >= 0) ++tp;
                    continue;
                }
                double s_old = res->beta_em * ((1.0 - cos(f_th[f0])) / f_area[f0] + (1.0 - cos(f_th[f1])) / f_area[f1]);
                double s_new = res->beta_em * ((1.0 - cos(nt0)) / f_area[f0] + (1.0 - cos(nt1)) / f_area[f1]);
                double ds = s_new - s_old;
                if (sw >= 0) ++tp;
                if (ds <= 0.0 || rng_u01(&rng) < exp(-ds)) {
                    if (sw >= 0) ++ta;
                    m.th[e] = new_t;
                    f_th[f0] = nt0;
                    f_th[f1] = nt1;
                }
            }
        }

        if (sw >= 0) {
            double atot, r2, eem, gbe, bie;
            ico_measure_all(&m, res->q_mono, &atot, &r2, &eem, &gbe, &bie);
            if (gbe > max_gb) max_gb = gbe;
            if (bie > max_bi) max_bi = bie;
            int bk = sw / bsz;
            if (bk >= JK_EM) bk = JK_EM - 1;
            blk_a[bk]  += atot;
            blk_r2[bk] += r2;
            blk_em[bk] += eem;
        }
    }

    jackknife_stats(JK_EM, blk_a, bsz, &res->a_mean, &res->a_err);
    jackknife_stats(JK_EM, blk_r2, bsz, &res->r2_mean, &res->r2_err);
    jackknife_stats(JK_EM, blk_em, bsz, &res->eem_mean, &res->eem_err);
    res->gb_max_err = max_gb;
    res->bianchi_max_err = max_bi;
    res->acc_l = (double)la / (double)lp;
    res->acc_th = (tp > 0) ? ((double)ta / (double)tp) : 0.0;
    return NULL;
}

/* -------------------------------------------------------------------------- */
/* 8. Main Driver                                                             */
/* -------------------------------------------------------------------------- */

int main(int argc, char **argv) {
    const char *json_path = NULL;
    for (int i = 1; i + 1 < argc; ++i) {
        if (strcmp(argv[i], "--json") == 0) json_path = argv[i + 1];
    }

    printf("=== PHY Stage 9: Complete Physics — SM Unification, 4D SU(2)-Higgs & Regge-Maxwell ===\n");
    selftest_res_t st;
    if (run_selftests(&st) != 0) return 1;

    /* Exp 1: 4D SU(2)-Higgs Electroweak Fradkin-Shenker Sweep across 5 parallel pthreads */
    const double kappas[N_KAPPA] = {0.10, 0.20, 0.27, 0.34, 0.45};
    exp1_ew_res_t exp1[N_KAPPA];
    pthread_t th1[N_KAPPA];
    for (int i = 0; i < N_KAPPA; ++i) {
        memset(&exp1[i], 0, sizeof(exp1[i]));
        exp1[i].beta = 2.3;
        exp1[i].kappa = kappas[i];
        exp1[i].lam = 0.5;
        pthread_create(&th1[i], NULL, exp1_ew_worker, &exp1[i]);
    }
    for (int i = 0; i < N_KAPPA; ++i) pthread_join(th1[i], NULL);

    printf("  [EXP1 6^4 SU(2)-Higgs EW Transition] beta=2.3, lam=0.5:\n");
    for (int i = 0; i < N_KAPPA; ++i) {
        exp1_ew_res_t *r = &exp1[i];
        printf("    kappa=%.2f | <P>=%.5f+-%.5f | <L_H>=%.5f+-%.5f | <rho^2>=%.4f+-%.4f | chi22=%.4f+-%.4f | h_acc=%.1f%%\n",
               r->kappa, r->p_mean, r->p_err, r->l_mean, r->l_err,
               r->rho2_mean, r->rho2_err, r->chi22_mean, r->chi22_err, 100.0 * r->higgs_acc);
    }

    /* Exp 2: SM vs MSSM 1-Loop & 2-Loop RG Gauge Unification */
    const double mz = 91.1876, m_susy = 1000.0, a_em_inv = 127.951, sw2 = 0.23122;
    const double a0[3] = {0.6 * a_em_inv * (1.0 - sw2), a_em_inv * sw2, 1.0 / 0.1180};
    const double b1_sm[3]   = {-41.0 / 10.0, 19.0 / 6.0, 7.0};
    const double b1_mssm[3] = {-33.0 / 5.0, -1.0, 3.0};
    const double B2_sm[3][3] = {
        {199.0/50.0, 27.0/10.0, 44.0/5.0},
        {9.0/10.0,   35.0/6.0,  12.0},
        {11.0/10.0,  9.0/2.0,   -26.0}
    };
    const double B2_mssm[3][3] = {
        {199.0/25.0, 27.0/5.0,  88.0/5.0},
        {9.0/5.0,    25.0,      24.0},
        {11.0/5.0,   9.0,       14.0}
    };

    /* Pairwise 1-loop crossings */
    double t12_sm = 2.0 * M_PI * (a0[1] - a0[0]) / (b1_sm[0] - b1_sm[1]);
    double t13_sm = 2.0 * M_PI * (a0[2] - a0[0]) / (b1_sm[0] - b1_sm[2]);
    double t23_sm = 2.0 * M_PI * (a0[2] - a0[1]) / (b1_sm[1] - b1_sm[2]);
    double lg12_sm = log10(mz * exp(t12_sm)), lg13_sm = log10(mz * exp(t13_sm)), lg23_sm = log10(mz * exp(t23_sm));

    double a_susy1[3], t_su = log(m_susy / mz);
    for (int i = 0; i < 3; ++i) a_susy1[i] = a0[i] + (b1_sm[i] / (2.0 * M_PI)) * t_su;
    double lg12_ms = log10(m_susy * exp(2.0 * M_PI * (a_susy1[1] - a_susy1[0]) / (b1_mssm[0] - b1_mssm[1])));
    double lg13_ms = log10(m_susy * exp(2.0 * M_PI * (a_susy1[2] - a_susy1[0]) / (b1_mssm[0] - b1_mssm[2])));
    double lg23_ms = log10(m_susy * exp(2.0 * M_PI * (a_susy1[2] - a_susy1[1]) / (b1_mssm[1] - b1_mssm[2])));

    /* 2-Loop scans */
    double sm_min_sp = 1e9, sm_best_lg = 0.0, sm_best_a[3] = {0};
    for (int k = 0; k <= 400; ++k) {
        double lg = 13.0 + 0.01 * (double)k;
        double y[3];
        rg_rk4(a0, b1_sm, B2_sm, log(pow(10.0, lg) / mz), 500, y);
        double mn = y[0] < y[1] ? y[0] : y[1]; if (y[2] < mn) mn = y[2];
        double mx = y[0] > y[1] ? y[0] : y[1]; if (y[2] > mx) mx = y[2];
        if (mx - mn < sm_min_sp) { sm_min_sp = mx - mn; sm_best_lg = lg; memcpy(sm_best_a, y, sizeof(y)); }
    }
    double a_susy2[3];
    rg_rk4(a0, b1_sm, B2_sm, t_su, 400, a_susy2);
    double ms_min_sp = 1e9, ms_best_lg = 0.0, ms_best_a[3] = {0};
    for (int k = 0; k <= 200; ++k) {
        double lg = 15.0 + 0.01 * (double)k;
        double y[3];
        rg_rk4(a_susy2, b1_mssm, B2_mssm, log(pow(10.0, lg) / m_susy), 500, y);
        double mn = y[0] < y[1] ? y[0] : y[1]; if (y[2] < mn) mn = y[2];
        double mx = y[0] > y[1] ? y[0] : y[1]; if (y[2] > mx) mx = y[2];
        if (mx - mn < ms_min_sp) { ms_min_sp = mx - mn; ms_best_lg = lg; memcpy(ms_best_a, y, sizeof(y)); }
    }
    double ir_sin2w_mssm = (3.0 * a0[1]) / (5.0 * a0[0] + 3.0 * a0[1]); /* at M_Z */
    double gut_sin2w = (3.0 * ms_best_a[1]) / (5.0 * ms_best_a[0] + 3.0 * ms_best_a[1]);

    printf("  [EXP2 SM vs MSSM RG Unification] M_Z=91.19 GeV (sin2W(MZ)=%.5f -> GUT=%.5f):\n",
           ir_sin2w_mssm, gut_sin2w);
    printf("    SM   1-loop crossings: (12=%.2f, 13=%.2f, 23=%.2f) | 2-loop min spread=%.3f at log10(mu)=%.2f\n",
           lg12_sm, lg13_sm, lg23_sm, sm_min_sp, sm_best_lg);
    printf("    MSSM 1-loop crossings: (12=%.2f, 13=%.2f, 23=%.2f) | 2-loop min spread=%.3f at log10(mu)=%.2f\n",
           lg12_ms, lg13_ms, lg23_ms, ms_min_sp, ms_best_lg);

    /* Exp 3: Coupled 2D Regge-Maxwell + Monopole Backreaction (4 parallel pthreads) */
    const char *em_names[N_EM_CASES] = {"Pure Grav (bEM=0, Q=0)", "Dyn Regge-Max (bEM=6, Q=0)", "Monopole Q=1 (bEM=6)", "Monopole Q=2 (bEM=6)"};
    const int qs[N_EM_CASES]       = {0, 0, 1, 2};
    const double bems[N_EM_CASES]  = {0.0, 6.0, 6.0, 6.0};
    exp3_em_res_t exp3[N_EM_CASES];
    pthread_t th3[N_EM_CASES];
    for (int i = 0; i < N_EM_CASES; ++i) {
        memset(&exp3[i], 0, sizeof(exp3[i]));
        exp3[i].q_mono = qs[i];
        exp3[i].lam = 0.5;
        exp3[i].alpha = 2.0;
        exp3[i].beta_em = bems[i];
        pthread_create(&th3[i], NULL, exp3_em_worker, &exp3[i]);
    }
    for (int i = 0; i < N_EM_CASES; ++i) pthread_join(th3[i], NULL);

    printf("  [EXP3 Coupled 2D Regge-Maxwell on S^2] lam=0.5, alpha=2.0:\n");
    for (int i = 0; i < N_EM_CASES; ++i) {
        exp3_em_res_t *r = &exp3[i];
        printf("    %-25s | <A>=%.4f+-%.4f <R2>=%.4f+-%.4f <E_EM>=%.4f | GB=%.2e Bi=%.2e | acc(l/t)=%.0f%%/%.0f%%\n",
               em_names[i], r->a_mean, r->a_err, r->r2_mean, r->r2_err, r->eem_mean,
               r->gb_max_err, r->bianchi_max_err, 100.0 * r->acc_l, 100.0 * r->acc_th);
    }

    /* Emit JSON */
    if (json_path) {
        FILE *fp = fopen(json_path, "w");
        if (!fp) { perror("fopen json"); return 1; }
        fprintf(fp, "{\n");
        fprintf(fp, "  \"selftest\": {\"anom_ok\": %d, \"uniq_total\": %d, \"uniq_bad\": %d, \"gut_sin2w_err\": %.6e, \"so4_err\": %.6e, \"ew_mgam_err\": %.6e, \"ew_rho_err\": %.6e, \"su2h_gi_err\": %.6e, \"rg_rk4_err\": %.6e, \"ico_bianchi_err\": %.6e, \"ico_q1_err\": %.6e, \"ico_q2_err\": %.6e},\n",
                st.anom_all_zero, st.uniq_total, st.uniq_bad, st.gut_sin2w_err, st.so4_norm_err,
                st.ew_mgam_err, st.ew_rho_err, st.su2_higgs_gi_err, st.rg_rk4_err,
                st.ico_bianchi_err, st.ico_q1_err, st.ico_q2_err);
        fprintf(fp, "  \"exp1_ew_higgs\": [\n");
        for (int i = 0; i < N_KAPPA; ++i) {
            exp1_ew_res_t *r = &exp1[i];
            fprintf(fp, "    {\"kappa\": %.2f, \"p_mean\": %.6f, \"p_err\": %.6f, \"l_mean\": %.6f, \"l_err\": %.6f, \"rho2_mean\": %.6f, \"rho2_err\": %.6f, \"chi22\": %.6f, \"chi22_err\": %.6f, \"h_acc\": %.4f}%s\n",
                    r->kappa, r->p_mean, r->p_err, r->l_mean, r->l_err, r->rho2_mean, r->rho2_err,
                    r->chi22_mean, r->chi22_err, r->higgs_acc, i + 1 < N_KAPPA ? "," : "");
        }
        fprintf(fp, "  ],\n");
        fprintf(fp, "  \"exp2_rg_unif\": {\n");
        fprintf(fp, "    \"sm_cross_1loop\": [%.6f, %.6f, %.6f], \"mssm_cross_1loop\": [%.6f, %.6f, %.6f],\n",
                lg12_sm, lg13_sm, lg23_sm, lg12_ms, lg13_ms, lg23_ms);
        fprintf(fp, "    \"sm_2loop\": {\"min_spread\": %.6f, \"log10_mu\": %.4f, \"alpha_inv\": [%.6f, %.6f, %.6f]},\n",
                sm_min_sp, sm_best_lg, sm_best_a[0], sm_best_a[1], sm_best_a[2]);
        fprintf(fp, "    \"mssm_2loop\": {\"min_spread\": %.6f, \"log10_mu\": %.4f, \"alpha_inv\": [%.6f, %.6f, %.6f]},\n",
                ms_min_sp, ms_best_lg, ms_best_a[0], ms_best_a[1], ms_best_a[2]);
        fprintf(fp, "    \"sin2w_ir\": %.6f, \"sin2w_gut\": %.6f\n  },\n", ir_sin2w_mssm, gut_sin2w);

        fprintf(fp, "  \"exp3_regge_maxwell\": [\n");
        for (int i = 0; i < N_EM_CASES; ++i) {
            exp3_em_res_t *r = &exp3[i];
            fprintf(fp, "    {\"name\": \"%s\", \"q_mono\": %d, \"beta_em\": %.2f, \"a_mean\": %.6f, \"a_err\": %.6f, \"r2_mean\": %.6f, \"r2_err\": %.6f, \"eem_mean\": %.6f, \"eem_err\": %.6f, \"gb_err\": %.6e, \"bi_err\": %.6e}%s\n",
                    em_names[i], r->q_mono, r->beta_em, r->a_mean, r->a_err, r->r2_mean, r->r2_err,
                    r->eem_mean, r->eem_err, r->gb_max_err, r->bianchi_max_err, i + 1 < N_EM_CASES ? "," : "");
        }
        fprintf(fp, "  ]\n}\n");
        fclose(fp);
    }
    return 0;
}
