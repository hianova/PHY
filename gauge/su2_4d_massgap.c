/* su2_4d_massgap.c — 4D Pure SU(2) Glueball 0++ Mass Gap & Direct OS Reflection Positivity Wind Tunnel
 * Tests Constructive QFT Theorems Q1–Q5 from `Flow.ConstructiveQFT` across beta in {2.30, 2.40, 2.50}:
 *   - Strictly time-slice-local 3D APE smearing (`smear_levels = {0, 2, 5, 10}`, temporal links untouched -> exact OS positivity)
 *   - 4 x 4 GEVP variational extraction of 0++ scalar glueball mass gap + 20-block jackknife cosh effective masses
 *   - Direct lattice check of C(0), C(1) positive eigenvalues (Theorem Q1) & Hankel log-convexity (Theorem Q2C)
 *   - 3D-smeared Wilson loops W_sm(R, T) & Luscher-subtracted static-quark string tension a*sqrt(sigma) -> m_{0++}/sqrt(sigma)
 */
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    NT_MG        = 12,
    NS_MG        = 10,
    N_BETAS      = 3,
    N_THERM      = 500,
    N_SWEEPS     = 6000,
    MEAS_STRIDE  = 2,
    N_MEAS       = N_SWEEPS / MEAS_STRIDE, /* 3000 measurements */
    S_BASIS      = 4,
    DT_MAX       = 4,
    W_MAX        = 4,
    JK_BLOCKS    = 20
};

static const int K_SMEAR_LEVELS[S_BASIS] = {0, 2, 5, 10};
static const double K_APE_ALPHA          = 0.5;
static const double K_PI                 = 3.14159265358979323846264338327950288;
static const double K_TWO_PI             = 6.28318530717958647692528676655900577;

/* --- Deterministic xoshiro256** PRNG --- */
typedef struct {
    uint64_t s[4];
} rng_t;

static inline uint64_t rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t splitmix64(uint64_t *state) {
    uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
    z          = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z          = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static void rng_seed(rng_t *r, uint64_t seed) {
    uint64_t st = seed;
    for (int i = 0; i < 4; ++i) {
        r->s[i] = splitmix64(&st);
    }
}

static inline uint64_t rng_u64(rng_t *r) {
    const uint64_t res = rotl64(r->s[1] * 5ULL, 7) * 9ULL;
    const uint64_t t   = r->s[1] << 17;
    r->s[2] ^= r->s[0];
    r->s[3] ^= r->s[1];
    r->s[1] ^= r->s[2];
    r->s[0] ^= r->s[3];
    r->s[2] ^= t;
    r->s[3] = rotl64(r->s[3], 45);
    return res;
}

static inline double rng_u01(rng_t *r) {
    return ((double)(rng_u64(r) >> 11)) * 0x1.0p-53;
}

static inline double rng_pos_u01(rng_t *r) {
    return ((double)(rng_u64(r) >> 11) + 0.5) * 0x1.0p-53;
}

/* --- SU(2) Quaternions matching Flow.YangMillsBundle --- */
typedef struct {
    double w, x, y, z;
} quat_t;

static inline quat_t q_id(void) {
    quat_t q = {1.0, 0.0, 0.0, 0.0};
    return q;
}

static inline quat_t q_conj(quat_t a) {
    quat_t r = {a.w, -a.x, -a.y, -a.z};
    return r;
}

static inline quat_t q_mul(quat_t a, quat_t b) {
    quat_t r = {
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w
    };
    return r;
}

static inline double q_inner(quat_t a, quat_t b) {
    return a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline double q_norm2(quat_t a) {
    return a.w * a.w + a.x * a.x + a.y * a.y + a.z * a.z;
}

static inline void q_normalize(quat_t *a) {
    double n = sqrt(q_norm2(*a));
    if (n > 1e-15) {
        double inv = 1.0 / n;
        a->w *= inv;
        a->x *= inv;
        a->y *= inv;
        a->z *= inv;
    } else {
        *a = q_id();
    }
}

/* --- Heap-Allocated 4D SU(2) Lattice --- */
typedef struct {
    int nt;
    int ns;
    int vol;
    quat_t *u;
} su2_lat_t;

static inline int next_i(int i, int n) {
    return (i + 1 < n) ? (i + 1) : 0;
}

static inline int prev_i(int i, int n) {
    return (i > 0) ? (i - 1) : (n - 1);
}

static inline size_t site_idx(const su2_lat_t *lat, int t, int x, int y, int z) {
    return (size_t)(((t * lat->ns + x) * lat->ns + y) * lat->ns + z);
}

static inline size_t link_idx(const su2_lat_t *lat, int t, int x, int y, int z, int mu) {
    return 4U * site_idx(lat, t, x, y, z) + (size_t)mu;
}

static void lat_alloc(su2_lat_t *lat, int nt, int ns) {
    lat->nt  = nt;
    lat->ns  = ns;
    lat->vol = nt * ns * ns * ns;
    lat->u   = (quat_t *)malloc(4U * (size_t)lat->vol * sizeof(quat_t));
    if (lat->u == NULL) {
        fprintf(stderr, "ERROR: OOM allocating %dx%d^3 SU(2) lattice\n", nt, ns);
        exit(1);
    }
}

static void lat_free(su2_lat_t *lat) {
    free(lat->u);
    lat->u = NULL;
}

static void lat_copy(su2_lat_t *dst, const su2_lat_t *src) {
    memcpy(dst->u, src->u, 4U * (size_t)src->vol * sizeof(quat_t));
}

static void lat_init_cold(su2_lat_t *lat) {
    for (int i = 0; i < 4 * lat->vol; ++i) {
        lat->u[i] = q_id();
    }
}

static inline void step_coord(const su2_lat_t *lat, const int c[4], int d, int dir, int out[4]) {
    out[0]  = c[0];
    out[1]  = c[1];
    out[2]  = c[2];
    out[3]  = c[3];
    int lim = (d == 0) ? lat->nt : lat->ns;
    out[d]  = (dir > 0) ? next_i(c[d], lim) : prev_i(c[d], lim);
}

static inline quat_t get_link(const su2_lat_t *lat, const int c[4], int mu) {
    return lat->u[link_idx(lat, c[0], c[1], c[2], c[3], mu)];
}

static quat_t plaq_quat(const su2_lat_t *lat, const int c[4], int mu, int nu) {
    int c_mu[4], c_nu[4];
    step_coord(lat, c, mu, +1, c_mu);
    step_coord(lat, c, nu, +1, c_nu);
    quat_t u1 = get_link(lat, c, mu);
    quat_t u2 = get_link(lat, c_mu, nu);
    quat_t u3 = get_link(lat, c_nu, mu);
    quat_t u4 = get_link(lat, c, nu);
    return q_mul(q_mul(u1, u2), q_mul(q_conj(u3), q_conj(u4)));
}

static quat_t compute_staple_4d(const su2_lat_t *lat, int t, int x, int y, int z, int mu) {
    const int c[4] = {t, x, y, z};
    int c_mu[4];
    step_coord(lat, c, mu, +1, c_mu);
    quat_t sum = {0.0, 0.0, 0.0, 0.0};
    for (int nu = 0; nu < 4; ++nu) {
        if (nu == mu) {
            continue;
        }
        int c_nu[4], c_mnu[4], c_mnumu[4];
        step_coord(lat, c, nu, +1, c_nu);
        step_coord(lat, c, nu, -1, c_mnu);
        step_coord(lat, c_mnu, mu, +1, c_mnumu);

        quat_t sf = q_mul(q_mul(get_link(lat, c, nu), get_link(lat, c_nu, mu)),
                          q_conj(get_link(lat, c_mu, nu)));
        quat_t sb = q_mul(q_mul(q_conj(get_link(lat, c_mnu, nu)), get_link(lat, c_mnu, mu)),
                          get_link(lat, c_mnumu, nu));

        sum.w += sf.w + sb.w;
        sum.x += sf.x + sb.x;
        sum.y += sf.y + sb.y;
        sum.z += sf.z + sb.z;
    }
    return sum;
}

static double cfg_total_action(const su2_lat_t *lat, double beta) {
    double s = 0.0;
    for (int t = 0; t < lat->nt; ++t) {
        for (int x = 0; x < lat->ns; ++x) {
            for (int y = 0; y < lat->ns; ++y) {
                for (int z = 0; z < lat->ns; ++z) {
                    const int c[4] = {t, x, y, z};
                    for (int mu = 0; mu < 4; ++mu) {
                        for (int nu = mu + 1; nu < 4; ++nu) {
                            s += 1.0 - plaq_quat(lat, c, mu, nu).w;
                        }
                    }
                }
            }
        }
    }
    return beta * s;
}

static double cfg_mean_plaq(const su2_lat_t *lat) {
    double sum_p = 0.0;
    for (int t = 0; t < lat->nt; ++t) {
        for (int x = 0; x < lat->ns; ++x) {
            for (int y = 0; y < lat->ns; ++y) {
                for (int z = 0; z < lat->ns; ++z) {
                    const int c[4] = {t, x, y, z};
                    for (int mu = 0; mu < 4; ++mu) {
                        for (int nu = mu + 1; nu < 4; ++nu) {
                            sum_p += plaq_quat(lat, c, mu, nu).w;
                        }
                    }
                }
            }
        }
    }
    return sum_p / (double)(6 * lat->vol);
}

static double cfg_abs_polyakov(const su2_lat_t *lat) {
    double sum_l    = 0.0;
    int space_sites = lat->ns * lat->ns * lat->ns;
    for (int x = 0; x < lat->ns; ++x) {
        for (int y = 0; y < lat->ns; ++y) {
            for (int z = 0; z < lat->ns; ++z) {
                quat_t p = q_id();
                for (int t = 0; t < lat->nt; ++t) {
                    p = q_mul(p, lat->u[link_idx(lat, t, x, y, z, 0)]);
                }
                sum_l += p.w;
            }
        }
    }
    return fabs(sum_l / (double)space_sites);
}

/* --- Creutz-Kennedy-Pendleton Exact SU(2) Heatbath & Checkerboard Over-Relaxation --- */
static double sample_su2_a0(double alpha, rng_t *rng) {
    if (alpha <= 1e-14) {
        for (;;) {
            double a0 = 2.0 * rng_u01(rng) - 1.0;
            if (rng_pos_u01(rng) < sqrt(fmax(0.0, 1.0 - a0 * a0))) {
                return a0;
            }
        }
    }
    if (alpha < 0.5) {
        for (;;) {
            double a0 = 2.0 * rng_u01(rng) - 1.0;
            double b  = sqrt(fmax(0.0, 1.0 - a0 * a0)) * exp(alpha * (a0 - 1.0));
            if (rng_pos_u01(rng) < b) {
                return a0;
            }
        }
    }
    for (;;) {
        double r1   = rng_pos_u01(rng);
        double r2   = rng_u01(rng);
        double r3   = rng_pos_u01(rng);
        double r4   = rng_u01(rng);
        double c    = cos(K_TWO_PI * r2);
        double lam2 = -(log(r1) + c * c * log(r3)) / (2.0 * alpha);
        if (lam2 > 1.0) {
            continue;
        }
        if (r4 * r4 <= 1.0 - lam2) {
            return 1.0 - 2.0 * lam2;
        }
    }
}

static quat_t sample_su2_boltzmann(double alpha, rng_t *rng) {
    double a0  = sample_su2_a0(alpha, rng);
    double z   = 2.0 * rng_u01(rng) - 1.0;
    double phi = K_TWO_PI * rng_u01(rng);
    double rho = sqrt(fmax(0.0, 1.0 - z * z));
    double rad = sqrt(fmax(0.0, 1.0 - a0 * a0));
    quat_t q   = {a0, rad * rho * cos(phi), rad * rho * sin(phi), rad * z};
    return q;
}

static void sweep_su2_hb_or(su2_lat_t *lat, double beta, rng_t *rng) {
    for (int par = 0; par < 2; ++par) {
        for (int t = 0; t < lat->nt; ++t) {
            for (int x = 0; x < lat->ns; ++x) {
                for (int y = 0; y < lat->ns; ++y) {
                    for (int z = 0; z < lat->ns; ++z) {
                        if (((t + x + y + z) & 1) != par) {
                            continue;
                        }
                        for (int mu = 0; mu < 4; ++mu) {
                            quat_t v  = compute_staple_4d(lat, t, x, y, z, mu);
                            double k  = sqrt(q_norm2(v));
                            quat_t xq = sample_su2_boltzmann(beta * k, rng);
                            size_t id = link_idx(lat, t, x, y, z, mu);
                            if (k > 1e-14) {
                                double inv = 1.0 / k;
                                quat_t vu  = {v.w * inv, v.x * inv, v.y * inv, v.z * inv};
                                quat_t u1  = q_mul(xq, vu);
                                q_normalize(&u1);
                                lat->u[id] = u1;
                            } else {
                                lat->u[id] = xq;
                            }
                        }
                    }
                }
            }
        }
    }

    for (int or_sw = 0; or_sw < 2; ++or_sw) {
        for (int par = 0; par < 2; ++par) {
            for (int t = 0; t < lat->nt; ++t) {
                for (int x = 0; x < lat->ns; ++x) {
                    for (int y = 0; y < lat->ns; ++y) {
                        for (int z = 0; z < lat->ns; ++z) {
                            if (((t + x + y + z) & 1) != par) {
                                continue;
                            }
                            for (int mu = 0; mu < 4; ++mu) {
                                quat_t v  = compute_staple_4d(lat, t, x, y, z, mu);
                                double v2 = q_norm2(v);
                                if (v2 > 1e-24) {
                                    size_t id = link_idx(lat, t, x, y, z, mu);
                                    quat_t u  = lat->u[id];
                                    double sc = (2.0 * q_inner(u, v)) / v2;
                                    quat_t u1 = {sc * v.w - u.w, sc * v.x - u.x,
                                                 sc * v.y - u.y, sc * v.z - u.z};
                                    q_normalize(&u1);
                                    lat->u[id] = u1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/* --- Strictly Time-Slice-Local 3D APE Smearing (Temporal links mu=0 UNTOUCHED -> Exact OS Positivity) --- */
static void smear_3d_ape_step(su2_lat_t *dst, const su2_lat_t *src, double c_ape) {
    for (int t = 0; t < src->nt; ++t) {
        for (int x = 0; x < src->ns; ++x) {
            for (int y = 0; y < src->ns; ++y) {
                for (int z = 0; z < src->ns; ++z) {
                    const int c[4] = {t, x, y, z};
                    /* Preserve temporal link untouched */
                    dst->u[link_idx(dst, t, x, y, z, 0)] = src->u[link_idx(src, t, x, y, z, 0)];

                    for (int mu = 1; mu <= 3; ++mu) {
                        int c_mu[4];
                        step_coord(src, c, mu, +1, c_mu);
                        quat_t st = {0.0, 0.0, 0.0, 0.0};
                        for (int nu = 1; nu <= 3; ++nu) {
                            if (nu == mu) {
                                continue;
                            }
                            int c_nu[4], c_mnu[4], c_mnumu[4];
                            step_coord(src, c, nu, +1, c_nu);
                            step_coord(src, c, nu, -1, c_mnu);
                            step_coord(src, c_mnu, mu, +1, c_mnumu);

                            quat_t sf = q_mul(q_mul(get_link(src, c, nu), get_link(src, c_nu, mu)),
                                              q_conj(get_link(src, c_mu, nu)));
                            quat_t sb = q_mul(q_mul(q_conj(get_link(src, c_mnu, nu)),
                                                    get_link(src, c_mnu, mu)),
                                              get_link(src, c_mnumu, nu));
                            st.w += sf.w + sb.w;
                            st.x += sf.x + sb.x;
                            st.y += sf.y + sb.y;
                            st.z += sf.z + sb.z;
                        }
                        quat_t u0  = get_link(src, c, mu);
                        quat_t out = {
                            u0.w + c_ape * st.w,
                            u0.x + c_ape * st.x,
                            u0.y + c_ape * st.y,
                            u0.z + c_ape * st.z
                        };
                        q_normalize(&out);
                        dst->u[link_idx(dst, t, x, y, z, mu)] = out;
                    }
                }
            }
        }
    }
}

/* 0++ scalar glueball time-slice operator averaged over all 3 spatial planes at time slice t */
static void measure_0pp_timeslices(const su2_lat_t *lat, double out_o[NT_MG]) {
    const double norm = 1.0 / (double)(3 * lat->ns * lat->ns * lat->ns);
    for (int t = 0; t < lat->nt; ++t) {
        double acc = 0.0;
        for (int x = 0; x < lat->ns; ++x) {
            for (int y = 0; y < lat->ns; ++y) {
                for (int z = 0; z < lat->ns; ++z) {
                    const int c[4] = {t, x, y, z};
                    acc += plaq_quat(lat, c, 1, 2).w;
                    acc += plaq_quat(lat, c, 1, 3).w;
                    acc += plaq_quat(lat, c, 2, 3).w;
                }
            }
        }
        out_o[t] = acc * norm;
    }
}

/* Cached incremental Wilson line workspace for 10x faster exact W(R, T), 1 <= R,T <= W_MAX */
typedef struct {
    quat_t *sp_line[3][W_MAX]; /* [dir-1][R-1][site] */
    quat_t *tm_line[W_MAX];    /* [T-1][site] */
} wloop_ws_t;

static void wloop_ws_alloc(wloop_ws_t *ws, int vol) {
    for (int d = 0; d < 3; ++d) {
        for (int r = 0; r < W_MAX; ++r) {
            ws->sp_line[d][r] = (quat_t *)malloc((size_t)vol * sizeof(quat_t));
            if (!ws->sp_line[d][r]) {
                fprintf(stderr, "ERROR: OOM in wloop_ws_alloc (sp)\n");
                exit(1);
            }
        }
    }
    for (int t = 0; t < W_MAX; ++t) {
        ws->tm_line[t] = (quat_t *)malloc((size_t)vol * sizeof(quat_t));
        if (!ws->tm_line[t]) {
            fprintf(stderr, "ERROR: OOM in wloop_ws_alloc (tm)\n");
            exit(1);
        }
    }
}

static void wloop_ws_free(wloop_ws_t *ws) {
    for (int d = 0; d < 3; ++d) {
        for (int r = 0; r < W_MAX; ++r) {
            free(ws->sp_line[d][r]);
        }
    }
    for (int t = 0; t < W_MAX; ++t) {
        free(ws->tm_line[t]);
    }
}

static void measure_smeared_wloops(const su2_lat_t *lat, wloop_ws_t *ws,
                                   double out_w[W_MAX + 1][W_MAX + 1]) {
    const int nt = lat->nt;
    const int ns = lat->ns;
    for (int t = 0; t < nt; ++t) {
        for (int x = 0; x < ns; ++x) {
            for (int y = 0; y < ns; ++y) {
                for (int z = 0; z < ns; ++z) {
                    size_t sid     = site_idx(lat, t, x, y, z);
                    for (int d = 0; d < 3; ++d) {
                        int mu       = d + 1;
                        int cur[4]   = {t, x, y, z};
                        quat_t prod  = q_id();
                        for (int r = 0; r < W_MAX; ++r) {
                            prod               = q_mul(prod, get_link(lat, cur, mu));
                            ws->sp_line[d][r][sid] = prod;
                            int nxt[4];
                            step_coord(lat, cur, mu, +1, nxt);
                            memcpy(cur, nxt, sizeof(cur));
                        }
                    }
                    int tcur[4] = {t, x, y, z};
                    quat_t tp   = q_id();
                    for (int tl = 0; tl < W_MAX; ++tl) {
                        tp               = q_mul(tp, get_link(lat, tcur, 0));
                        ws->tm_line[tl][sid] = tp;
                        tcur[0]          = next_i(tcur[0], nt);
                    }
                }
            }
        }
    }

    for (int r = 1; r <= W_MAX; ++r) {
        for (int tl = 1; tl <= W_MAX; ++tl) {
            out_w[r][tl] = 0.0;
        }
    }

    const double norm = 1.0 / (double)(3 * lat->vol);
    for (int t = 0; t < nt; ++t) {
        for (int x = 0; x < ns; ++x) {
            for (int y = 0; y < ns; ++y) {
                for (int z = 0; z < ns; ++z) {
                    size_t s0 = site_idx(lat, t, x, y, z);
                    int tp    = (t + 0) % nt;
                    for (int tl = 1; tl <= W_MAX; ++tl) {
                        tp           = (t + tl) % nt;
                        size_t s_top = site_idx(lat, tp, x, y, z);
                        quat_t t_left = ws->tm_line[tl - 1][s0];

                        /* d = 0 (x-dir), d = 1 (y-dir), d = 2 (z-dir) */
                        for (int r = 1; r <= W_MAX; ++r) {
                            size_t sr_x = site_idx(lat, t, (x + r) % ns, y, z);
                            size_t sr_y = site_idx(lat, t, x, (y + r) % ns, z);
                            size_t sr_z = site_idx(lat, t, x, y, (z + r) % ns);

                            out_w[r][tl] += q_inner(
                                q_mul(ws->sp_line[0][r - 1][s0], ws->tm_line[tl - 1][sr_x]),
                                q_mul(t_left, ws->sp_line[0][r - 1][s_top]));
                            out_w[r][tl] += q_inner(
                                q_mul(ws->sp_line[1][r - 1][s0], ws->tm_line[tl - 1][sr_y]),
                                q_mul(t_left, ws->sp_line[1][r - 1][s_top]));
                            out_w[r][tl] += q_inner(
                                q_mul(ws->sp_line[2][r - 1][s0], ws->tm_line[tl - 1][sr_z]),
                                q_mul(t_left, ws->sp_line[2][r - 1][s_top]));
                        }
                    }
                }
            }
        }
    }
    for (int r = 1; r <= W_MAX; ++r) {
        for (int tl = 1; tl <= W_MAX; ++tl) {
            out_w[r][tl] *= norm;
        }
    }
}

/* --- Exact 4x4 Jacobi Symmetric Eigensolver & Variational GEVP --- */
static void jacobi_sym_4x4(const double a_in[S_BASIS][S_BASIS], double evals[S_BASIS],
                           double evecs[S_BASIS][S_BASIS]) {
    double a[S_BASIS][S_BASIS];
    for (int i = 0; i < S_BASIS; ++i) {
        for (int j = 0; j < S_BASIS; ++j) {
            a[i][j]     = 0.5 * (a_in[i][j] + a_in[j][i]);
            evecs[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int it = 0; it < 64; ++it) {
        int p = 0, q = 1;
        double mx = 0.0;
        for (int i = 0; i < S_BASIS; ++i) {
            for (int j = i + 1; j < S_BASIS; ++j) {
                double v = fabs(a[i][j]);
                if (v > mx) {
                    mx = v;
                    p  = i;
                    q  = j;
                }
            }
        }
        if (mx < 1e-15) {
            break;
        }
        double app = a[p][p], aqq = a[q][q], apq = a[p][q];
        double tau = (aqq - app) / (2.0 * apq);
        double t   = (tau >= 0.0) ? (1.0 / (tau + sqrt(1.0 + tau * tau)))
                                  : (-1.0 / (-tau + sqrt(1.0 + tau * tau)));
        double c   = 1.0 / sqrt(1.0 + t * t);
        double s   = t * c;

        a[p][p] = app - t * apq;
        a[q][q] = aqq + t * apq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;
        for (int r = 0; r < S_BASIS; ++r) {
            if (r == p || r == q) {
                continue;
            }
            double arp = a[r][p], arq = a[r][q];
            a[r][p] = c * arp - s * arq;
            a[p][r] = a[r][p];
            a[r][q] = s * arp + c * arq;
            a[q][r] = a[r][q];
        }
        for (int r = 0; r < S_BASIS; ++r) {
            double vp = evecs[r][p], vq = evecs[r][q];
            evecs[r][p] = c * vp - s * vq;
            evecs[r][q] = s * vp + c * vq;
        }
    }

    for (int i = 0; i < S_BASIS; ++i) {
        evals[i] = a[i][i];
    }
    /* Sort strictly descending by eigenvalue */
    for (int i = 0; i < S_BASIS - 1; ++i) {
        int best = i;
        for (int j = i + 1; j < S_BASIS; ++j) {
            if (evals[j] > evals[best]) {
                best = j;
            }
        }
        if (best != i) {
            double tmp  = evals[i];
            evals[i]    = evals[best];
            evals[best] = tmp;
            for (int r = 0; r < S_BASIS; ++r) {
                double ev      = evecs[r][i];
                evecs[r][i]    = evecs[r][best];
                evecs[r][best] = ev;
            }
        }
    }
}

/* Solve C1 u = lambda C0 u (keeping positive subspace of C0), return largest lambda_0 & unit vector u_opt */
static double solve_gevp_4x4(const double c1[S_BASIS][S_BASIS], const double c0[S_BASIS][S_BASIS],
                             double u_opt[S_BASIS]) {
    double eval0[S_BASIS], evec0[S_BASIS][S_BASIS];
    jacobi_sym_4x4(c0, eval0, evec0);

    double dmax = fmax(fabs(eval0[0]), 1e-30);
    double Wmat[S_BASIS][S_BASIS];
    memset(Wmat, 0, sizeof(Wmat));
    for (int k = 0; k < S_BASIS; ++k) {
        if (eval0[k] > 1e-13 * dmax) {
            double inv_rt = 1.0 / sqrt(eval0[k]);
            for (int i = 0; i < S_BASIS; ++i) {
                Wmat[i][k] = evec0[i][k] * inv_rt;
            }
        }
    }

    double atilde[S_BASIS][S_BASIS];
    for (int i = 0; i < S_BASIS; ++i) {
        for (int j = 0; j < S_BASIS; ++j) {
            double s = 0.0;
            for (int p = 0; p < S_BASIS; ++p) {
                for (int q = 0; q < S_BASIS; ++q) {
                    s += Wmat[p][i] * c1[p][q] * Wmat[q][j];
                }
            }
            atilde[i][j] = s;
        }
    }

    double gevals[S_BASIS], gevecs[S_BASIS][S_BASIS];
    jacobi_sym_4x4(atilde, gevals, gevecs);

    double u2 = 0.0;
    for (int i = 0; i < S_BASIS; ++i) {
        double s = 0.0;
        for (int k = 0; k < S_BASIS; ++k) {
            s += Wmat[i][k] * gevecs[k][0];
        }
        u_opt[i] = s;
        u2 += s * s;
    }
    if (u2 > 1e-30) {
        double inv = 1.0 / sqrt(u2);
        for (int i = 0; i < S_BASIS; ++i) {
            u_opt[i] *= inv;
        }
    } else {
        u_opt[0] = 1.0;
        for (int i = 1; i < S_BASIS; ++i) {
            u_opt[i] = 0.0;
        }
    }
    return gevals[0];
}

static int extract_cosh_mass(double c_t, double c_tp1, int t, int nt, double *out_m) {
    if (c_t <= 0.0 || c_tp1 <= 0.0) {
        *out_m = 0.0;
        return 0;
    }
    double r = c_tp1 / c_t;
    if (r >= 1.0) {
        *out_m = 0.0;
        return 0;
    }
    double a   = (double)(nt / 2 - t);
    double b   = (double)(nt / 2 - (t + 1));
    double lo  = 0.0;
    double hi  = 8.0;
    if (cosh(b * hi) / cosh(a * hi) >= r) {
        *out_m = -log(r);
        return 1;
    }
    for (int it = 0; it < 64; ++it) {
        double mid = 0.5 * (lo + hi);
        double val = cosh(b * mid) / cosh(a * mid);
        if (val > r) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    *out_m = 0.5 * (lo + hi);
    return 1;
}

/* --- Startup Mathematical & Physics Self-Tests --- */
static void run_startup_selftests(void) {
    su2_lat_t g, gp, sm0, sm1;
    lat_alloc(&g, NT_MG, NS_MG);
    lat_alloc(&gp, NT_MG, NS_MG);
    lat_alloc(&sm0, NT_MG, NS_MG);
    lat_alloc(&sm1, NT_MG, NS_MG);

    rng_t rng;
    rng_seed(&rng, 0x50204A55ULL);

    for (int i = 0; i < 4 * g.vol; ++i) {
        g.u[i] = sample_su2_boltzmann(0.0, &rng);
    }

    /* 1. 3D APE smearing + O_a(t) + W_sm(R,T) exact 4D local SU(2) gauge covariance */
    quat_t *gmap = (quat_t *)malloc((size_t)g.vol * sizeof(quat_t));
    if (!gmap) {
        fprintf(stderr, "ERROR: OOM in selftest 1\n");
        exit(1);
    }
    for (int i = 0; i < g.vol; ++i) {
        gmap[i] = sample_su2_boltzmann(0.0, &rng);
    }
    for (int t = 0; t < g.nt; ++t) {
        for (int x = 0; x < g.ns; ++x) {
            for (int y = 0; y < g.ns; ++y) {
                for (int z = 0; z < g.ns; ++z) {
                    const int c[4] = {t, x, y, z};
                    quat_t gs      = gmap[site_idx(&g, t, x, y, z)];
                    for (int mu = 0; mu < 4; ++mu) {
                        int cp[4];
                        step_coord(&g, c, mu, +1, cp);
                        quat_t gsp = gmap[site_idx(&g, cp[0], cp[1], cp[2], cp[3])];
                        gp.u[link_idx(&gp, t, x, y, z, mu)] =
                            q_mul(q_mul(gs, g.u[link_idx(&g, t, x, y, z, mu)]), q_conj(gsp));
                    }
                }
            }
        }
    }
    free(gmap);

    lat_copy(&sm0, &g);
    lat_copy(&sm1, &gp);
    su2_lat_t tmp_sm;
    lat_alloc(&tmp_sm, NT_MG, NS_MG);
    for (int s = 0; s < 5; ++s) {
        smear_3d_ape_step(&tmp_sm, &sm0, K_APE_ALPHA);
        lat_copy(&sm0, &tmp_sm);
        smear_3d_ape_step(&tmp_sm, &sm1, K_APE_ALPHA);
        lat_copy(&sm1, &tmp_sm);
    }
    double o0[NT_MG], o1[NT_MG];
    measure_0pp_timeslices(&sm0, o0);
    measure_0pp_timeslices(&sm1, o1);
    double max_do = 0.0;
    for (int t = 0; t < NT_MG; ++t) {
        double d = fabs(o0[t] - o1[t]);
        if (d > max_do) {
            max_do = d;
        }
    }
    wloop_ws_t ws;
    wloop_ws_alloc(&ws, g.vol);
    double w0[W_MAX + 1][W_MAX + 1], w1[W_MAX + 1][W_MAX + 1];
    measure_smeared_wloops(&sm0, &ws, w0);
    measure_smeared_wloops(&sm1, &ws, w1);
    double max_dw = 0.0;
    for (int r = 1; r <= W_MAX; ++r) {
        for (int tl = 1; tl <= W_MAX; ++tl) {
            double d = fabs(w0[r][tl] - w1[r][tl]);
            if (d > max_dw) {
                max_dw = d;
            }
        }
    }
    if (max_do > 1e-10 || max_dw > 1e-10) {
        fprintf(stderr, "SELFTEST FAIL 1: 3D APE gauge invariance dO=%.3e dW=%.3e\n", max_do, max_dw);
        exit(1);
    }

    /* 2. Exact 4x4 Jacobi & GEVP recovery on known synthetic spectral matrix */
    const double masses[S_BASIS] = {0.6, 1.5, 2.3, 3.1};
    const double vmat[S_BASIS][S_BASIS] = {
        {1.0,  0.4, -0.2,  0.1},
        {0.5, -0.9,  0.3, -0.2},
        {0.2,  0.5,  0.8, -0.3},
        {0.1, -0.2,  0.4,  0.9}
    };
    double syn_c0[S_BASIS][S_BASIS], syn_c1[S_BASIS][S_BASIS];
    for (int a = 0; a < S_BASIS; ++a) {
        for (int b = 0; b < S_BASIS; ++b) {
            double s0 = 0.0, s1 = 0.0;
            for (int k = 0; k < S_BASIS; ++k) {
                s0 += vmat[k][a] * vmat[k][b];
                s1 += vmat[k][a] * vmat[k][b] * exp(-masses[k]);
            }
            syn_c0[a][b] = s0;
            syn_c1[a][b] = s1;
        }
    }
    double u_syn[S_BASIS];
    double lam_syn = solve_gevp_4x4(syn_c1, syn_c0, u_syn);
    double m_rec   = -log(lam_syn);
    if (fabs(m_rec - 0.6) > 1e-10) {
        fprintf(stderr, "SELFTEST FAIL 2: GEVP synthetic m0=%.12f != 0.6\n", m_rec);
        exit(1);
    }

    /* 3. CKP heatbath + 2x OR microcanonical action check on 12x10^3 */
    lat_copy(&gp, &g);
    double a_before = cfg_total_action(&gp, 2.4);
    for (int or_sw = 0; or_sw < 2; ++or_sw) {
        for (int par = 0; par < 2; ++par) {
            for (int t = 0; t < gp.nt; ++t) {
                for (int x = 0; x < gp.ns; ++x) {
                    for (int y = 0; y < gp.ns; ++y) {
                        for (int z = 0; z < gp.ns; ++z) {
                            if (((t + x + y + z) & 1) != par) {
                                continue;
                            }
                            for (int mu = 0; mu < 4; ++mu) {
                                quat_t v  = compute_staple_4d(&gp, t, x, y, z, mu);
                                double v2 = q_norm2(v);
                                if (v2 > 1e-24) {
                                    size_t id = link_idx(&gp, t, x, y, z, mu);
                                    quat_t u  = gp.u[id];
                                    double sc = (2.0 * q_inner(u, v)) / v2;
                                    quat_t u1 = {sc * v.w - u.w, sc * v.x - u.x,
                                                 sc * v.y - u.y, sc * v.z - u.z};
                                    q_normalize(&u1);
                                    gp.u[id] = u1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    double a_after = cfg_total_action(&gp, 2.4);
    if (fabs(a_after - a_before) > 1e-8) {
        fprintf(stderr, "SELFTEST FAIL 3: OR action drift %.3e\n", a_after - a_before);
        exit(1);
    }

    /* 4. OS Time-Slice Locality: randomizing ALL temporal links mu=0 leaves 3D APE O_a(t) 0.0 unchanged */
    lat_copy(&gp, &g);
    for (int t = 0; t < gp.nt; ++t) {
        for (int x = 0; x < gp.ns; ++x) {
            for (int y = 0; y < gp.ns; ++y) {
                for (int z = 0; z < gp.ns; ++z) {
                    gp.u[link_idx(&gp, t, x, y, z, 0)] = sample_su2_boltzmann(0.0, &rng);
                }
            }
        }
    }
    lat_copy(&sm0, &g);
    lat_copy(&sm1, &gp);
    for (int s = 0; s < 5; ++s) {
        smear_3d_ape_step(&tmp_sm, &sm0, K_APE_ALPHA);
        lat_copy(&sm0, &tmp_sm);
        smear_3d_ape_step(&tmp_sm, &sm1, K_APE_ALPHA);
        lat_copy(&sm1, &tmp_sm);
    }
    measure_0pp_timeslices(&sm0, o0);
    measure_0pp_timeslices(&sm1, o1);
    double max_dt0 = 0.0;
    for (int t = 0; t < NT_MG; ++t) {
        double d = fabs(o0[t] - o1[t]);
        if (d > max_dt0) {
            max_dt0 = d;
        }
    }
    if (max_dt0 > 1e-15) {
        fprintf(stderr, "SELFTEST FAIL 4: O_a(t) leaked temporal link dependence %.3e\n", max_dt0);
        exit(1);
    }

    printf("  [SELFTEST] PASS (3D APE gauge <%.2e, GEVP syn m0=%.6f exact, OR dS=%.2e, OS time-slice local dO=%.1e)\n",
           fmax(max_do, max_dw), m_rec, fabs(a_after - a_before), max_dt0);

    wloop_ws_free(&ws);
    lat_free(&tmp_sm);
    lat_free(&g);
    lat_free(&gp);
    lat_free(&sm0);
    lat_free(&sm1);
}

/* --- Simulation & 20-Block Jackknife Spectral Pipeline --- */
typedef struct {
    int idx;
    double beta;
    double mean_p;
    double err_p;
    double mean_abs_l;
    double err_abs_l;
    double eig_c0[S_BASIS];
    double eig_c1[S_BASIS];
    double c_raw[DT_MAX + 1];
    double err_c_raw[DT_MAX + 1];
    double c_gevp[DT_MAX + 1];
    double err_c_gevp[DT_MAX + 1];
    double hankel_ratio;
    double m_raw_01;
    double err_m_raw_01;
    double m_gevp_01;
    double err_m_gevp_01;
    double m_gevp_12;
    double err_m_gevp_12;
    double m_gevp_23;
    double err_m_gevp_23;
    int valid_m_gevp_23;
    double am_0pp;
    double err_am_0pp;
    double w_sm[W_MAX + 1][W_MAX + 1];
    double v_pot[W_MAX + 1];
    double a2_sigma;
    double a_sqrt_sigma;
    double err_a_sqrt_sigma;
    double chi_sm22;
    double err_chi_sm22;
    double chi_sm33;
    double err_chi_sm33;
    double ratio_m_over_sig;
    double err_ratio_m_over_sig;
} mg_res_t;

static void compute_simple_jk_err(const double *vals, int n, double *out_mean, double *out_err) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += vals[i];
    }
    double mean = sum / (double)n;
    int bsz     = n / JK_BLOCKS;
    double jk[JK_BLOCKS];
    double jk_m = 0.0;
    for (int b = 0; b < JK_BLOCKS; ++b) {
        double bs = 0.0;
        int lo = b * bsz, hi = (b == JK_BLOCKS - 1) ? n : ((b + 1) * bsz);
        int cnt = hi - lo;
        for (int i = lo; i < hi; ++i) {
            bs += vals[i];
        }
        jk[b] = (sum - bs) / (double)(n - cnt);
        jk_m += jk[b];
    }
    jk_m /= (double)JK_BLOCKS;
    double var = 0.0;
    for (int b = 0; b < JK_BLOCKS; ++b) {
        double d = jk[b] - jk_m;
        var += d * d;
    }
    var *= (double)(JK_BLOCKS - 1) / (double)JK_BLOCKS;
    *out_mean = mean;
    *out_err  = sqrt(var);
}

static void build_conn_corr_matrix(const double o_mean[S_BASIS],
                                     const double prod_mean[DT_MAX + 1][S_BASIS][S_BASIS],
                                     double c_mat[DT_MAX + 1][S_BASIS][S_BASIS]) {
    for (int dt = 0; dt <= DT_MAX; ++dt) {
        for (int a = 0; a < S_BASIS; ++a) {
            for (int b = 0; b < S_BASIS; ++b) {
                double sym_p = 0.5 * (prod_mean[dt][a][b] + prod_mean[dt][b][a]);
                c_mat[dt][a][b] = sym_p - o_mean[a] * o_mean[b];
            }
        }
    }
}

static void *massgap_worker(void *arg) {
    mg_res_t *res = (mg_res_t *)arg;
    su2_lat_t lat, sm_a, sm_b;
    lat_alloc(&lat, NT_MG, NS_MG);
    lat_alloc(&sm_a, NT_MG, NS_MG);
    lat_alloc(&sm_b, NT_MG, NS_MG);
    lat_init_cold(&lat);

    wloop_ws_t ws;
    wloop_ws_alloc(&ws, lat.vol);

    rng_t rng;
    rng_seed(&rng, 0xA4D00000ULL + (uint64_t)res->idx * 7777ULL);

    double *p_h    = (double *)malloc((size_t)N_MEAS * sizeof(double));
    double *absl_h = (double *)malloc((size_t)N_MEAS * sizeof(double));
    double (*o_h)[S_BASIS] = malloc((size_t)N_MEAS * sizeof(*o_h));
    double (*prod_h)[DT_MAX + 1][S_BASIS][S_BASIS] = malloc((size_t)N_MEAS * sizeof(*prod_h));
    double (*w_h)[W_MAX + 1][W_MAX + 1] = malloc((size_t)N_MEAS * sizeof(*w_h));
    if (!p_h || !absl_h || !o_h || !prod_h || !w_h) {
        fprintf(stderr, "ERROR: OOM in massgap_worker\n");
        exit(1);
    }

    for (int t = 0; t < N_THERM; ++t) {
        sweep_su2_hb_or(&lat, res->beta, &rng);
    }

    for (int m = 0; m < N_MEAS; ++m) {
        for (int s = 0; s < MEAS_STRIDE; ++s) {
            sweep_su2_hb_or(&lat, res->beta, &rng);
        }
        p_h[m]    = cfg_mean_plaq(&lat);
        absl_h[m] = cfg_abs_polyakov(&lat);

        double o_slices[S_BASIS][NT_MG];
        lat_copy(&sm_a, &lat);
        int cur_level = 0;
        for (int a = 0; a < S_BASIS; ++a) {
            int target = K_SMEAR_LEVELS[a];
            while (cur_level < target) {
                smear_3d_ape_step(&sm_b, &sm_a, K_APE_ALPHA);
                lat_copy(&sm_a, &sm_b);
                ++cur_level;
            }
            measure_0pp_timeslices(&sm_a, o_slices[a]);
            if (target == 5) {
                measure_smeared_wloops(&sm_a, &ws, w_h[m]);
            }
        }

        for (int a = 0; a < S_BASIS; ++a) {
            double so = 0.0;
            for (int t = 0; t < NT_MG; ++t) {
                so += o_slices[a][t];
            }
            o_h[m][a] = so / (double)NT_MG;
        }
        for (int dt = 0; dt <= DT_MAX; ++dt) {
            for (int a = 0; a < S_BASIS; ++a) {
                for (int b = 0; b < S_BASIS; ++b) {
                    double sp = 0.0;
                    for (int t = 0; t < NT_MG; ++t) {
                        int tp = (t + dt) % NT_MG;
                        sp += o_slices[a][t] * o_slices[b][tp];
                    }
                    prod_h[m][dt][a][b] = sp / (double)NT_MG;
                }
            }
        }
    }

    compute_simple_jk_err(p_h, N_MEAS, &res->mean_p, &res->err_p);
    compute_simple_jk_err(absl_h, N_MEAS, &res->mean_abs_l, &res->err_abs_l);

    /* Full-sample means + 20-block jackknife spectral analysis */
    double tot_o[S_BASIS];
    double tot_prod[DT_MAX + 1][S_BASIS][S_BASIS];
    double tot_w[W_MAX + 1][W_MAX + 1];
    memset(tot_o, 0, sizeof(tot_o));
    memset(tot_prod, 0, sizeof(tot_prod));
    memset(tot_w, 0, sizeof(tot_w));

    for (int m = 0; m < N_MEAS; ++m) {
        for (int a = 0; a < S_BASIS; ++a) {
            tot_o[a] += o_h[m][a];
        }
        for (int dt = 0; dt <= DT_MAX; ++dt) {
            for (int a = 0; a < S_BASIS; ++a) {
                for (int b = 0; b < S_BASIS; ++b) {
                    tot_prod[dt][a][b] += prod_h[m][dt][a][b];
                }
            }
        }
        for (int r = 1; r <= W_MAX; ++r) {
            for (int tl = 1; tl <= W_MAX; ++tl) {
                tot_w[r][tl] += w_h[m][r][tl];
            }
        }
    }

    double full_o[S_BASIS], full_prod[DT_MAX + 1][S_BASIS][S_BASIS];
    for (int a = 0; a < S_BASIS; ++a) {
        full_o[a] = tot_o[a] / (double)N_MEAS;
    }
    for (int dt = 0; dt <= DT_MAX; ++dt) {
        for (int a = 0; a < S_BASIS; ++a) {
            for (int b = 0; b < S_BASIS; ++b) {
                full_prod[dt][a][b] = tot_prod[dt][a][b] / (double)N_MEAS;
            }
        }
    }
    for (int r = 1; r <= W_MAX; ++r) {
        for (int tl = 1; tl <= W_MAX; ++tl) {
            res->w_sm[r][tl] = tot_w[r][tl] / (double)N_MEAS;
        }
    }

    double full_c[DT_MAX + 1][S_BASIS][S_BASIS];
    build_conn_corr_matrix(full_o, full_prod, full_c);

    double dummy_v[S_BASIS][S_BASIS];
    jacobi_sym_4x4(full_c[0], res->eig_c0, dummy_v);
    jacobi_sym_4x4(full_c[1], res->eig_c1, dummy_v);

    double u_full[S_BASIS];
    solve_gevp_4x4(full_c[1], full_c[0], u_full);

    for (int dt = 0; dt <= DT_MAX; ++dt) {
        res->c_raw[dt] = full_c[dt][0][0];
        double cg      = 0.0;
        for (int a = 0; a < S_BASIS; ++a) {
            for (int b = 0; b < S_BASIS; ++b) {
                cg += u_full[a] * full_c[dt][a][b] * u_full[b];
            }
        }
        res->c_gevp[dt] = cg;
    }
    res->hankel_ratio = (res->c_gevp[0] * res->c_gevp[2]) / (res->c_gevp[1] * res->c_gevp[1]);

    extract_cosh_mass(res->c_raw[0], res->c_raw[1], 0, NT_MG, &res->m_raw_01);
    extract_cosh_mass(res->c_gevp[0], res->c_gevp[1], 0, NT_MG, &res->m_gevp_01);
    extract_cosh_mass(res->c_gevp[1], res->c_gevp[2], 1, NT_MG, &res->m_gevp_12);
    res->valid_m_gevp_23 = extract_cosh_mass(res->c_gevp[2], res->c_gevp[3], 2, NT_MG, &res->m_gevp_23);
    res->am_0pp          = res->m_gevp_12;

    for (int r = 1; r <= W_MAX; ++r) {
        res->v_pot[r] = log(res->w_sm[r][2] / res->w_sm[r][3]);
    }
    res->a2_sigma     = res->v_pot[3] - res->v_pot[2] - (K_PI / 72.0);
    res->a_sqrt_sigma = (res->a2_sigma > 0.0) ? sqrt(res->a2_sigma) : 0.0;
    res->chi_sm22     = -log((res->w_sm[2][2] * res->w_sm[1][1]) / (res->w_sm[1][2] * res->w_sm[1][2]));
    res->chi_sm33     = -log((res->w_sm[3][3] * res->w_sm[2][2]) / (res->w_sm[2][3] * res->w_sm[2][3]));
    res->ratio_m_over_sig = (res->a_sqrt_sigma > 1e-12) ? (res->am_0pp / res->a_sqrt_sigma) : 0.0;

    /* Leave-one-out 20-block jackknife errors */
    int bsz = N_MEAS / JK_BLOCKS;
    double jk_craw[JK_BLOCKS][DT_MAX + 1];
    double jk_cgevp[JK_BLOCKS][DT_MAX + 1];
    double jk_mr01[JK_BLOCKS], jk_mg01[JK_BLOCKS], jk_mg12[JK_BLOCKS], jk_mg23[JK_BLOCKS];
    double jk_sig[JK_BLOCKS], jk_c22[JK_BLOCKS], jk_c33[JK_BLOCKS], jk_rat[JK_BLOCKS];
    int all_v23 = res->valid_m_gevp_23;

    for (int b = 0; b < JK_BLOCKS; ++b) {
        int lo = b * bsz, hi = (b == JK_BLOCKS - 1) ? N_MEAS : ((b + 1) * bsz);
        int rem = N_MEAS - (hi - lo);
        double bo[S_BASIS], bprod[DT_MAX + 1][S_BASIS][S_BASIS], bw[W_MAX + 1][W_MAX + 1];
        memset(bo, 0, sizeof(bo));
        memset(bprod, 0, sizeof(bprod));
        memset(bw, 0, sizeof(bw));
        for (int m = lo; m < hi; ++m) {
            for (int a = 0; a < S_BASIS; ++a) {
                bo[a] += o_h[m][a];
            }
            for (int dt = 0; dt <= DT_MAX; ++dt) {
                for (int a = 0; a < S_BASIS; ++a) {
                    for (int bb = 0; bb < S_BASIS; ++bb) {
                        bprod[dt][a][bb] += prod_h[m][dt][a][bb];
                    }
                }
            }
            for (int r = 1; r <= W_MAX; ++r) {
                for (int tl = 1; tl <= W_MAX; ++tl) {
                    bw[r][tl] += w_h[m][r][tl];
                }
            }
        }
        double jko[S_BASIS], jkprod[DT_MAX + 1][S_BASIS][S_BASIS], jkw[W_MAX + 1][W_MAX + 1];
        for (int a = 0; a < S_BASIS; ++a) {
            jko[a] = (tot_o[a] - bo[a]) / (double)rem;
        }
        for (int dt = 0; dt <= DT_MAX; ++dt) {
            for (int a = 0; a < S_BASIS; ++a) {
                for (int bb = 0; bb < S_BASIS; ++bb) {
                    jkprod[dt][a][bb] = (tot_prod[dt][a][bb] - bprod[dt][a][bb]) / (double)rem;
                }
            }
        }
        for (int r = 1; r <= W_MAX; ++r) {
            for (int tl = 1; tl <= W_MAX; ++tl) {
                jkw[r][tl] = (tot_w[r][tl] - bw[r][tl]) / (double)rem;
            }
        }
        double jkc[DT_MAX + 1][S_BASIS][S_BASIS];
        build_conn_corr_matrix(jko, jkprod, jkc);
        double u_b[S_BASIS];
        solve_gevp_4x4(jkc[1], jkc[0], u_b);
        for (int dt = 0; dt <= DT_MAX; ++dt) {
            jk_craw[b][dt] = jkc[dt][0][0];
            double cg      = 0.0;
            for (int a = 0; a < S_BASIS; ++a) {
                for (int bb = 0; bb < S_BASIS; ++bb) {
                    cg += u_b[a] * jkc[dt][a][bb] * u_b[bb];
                }
            }
            jk_cgevp[b][dt] = cg;
        }
        extract_cosh_mass(jk_craw[b][0], jk_craw[b][1], 0, NT_MG, &jk_mr01[b]);
        extract_cosh_mass(jk_cgevp[b][0], jk_cgevp[b][1], 0, NT_MG, &jk_mg01[b]);
        extract_cosh_mass(jk_cgevp[b][1], jk_cgevp[b][2], 1, NT_MG, &jk_mg12[b]);
        int v23 = extract_cosh_mass(jk_cgevp[b][2], jk_cgevp[b][3], 2, NT_MG, &jk_mg23[b]);
        if (!v23) {
            all_v23 = 0;
        }

        double v2  = log(jkw[2][2] / jkw[2][3]);
        double v3  = log(jkw[3][2] / jkw[3][3]);
        double a2s = v3 - v2 - (K_PI / 72.0);
        jk_sig[b]  = (a2s > 0.0) ? sqrt(a2s) : 0.0;
        jk_c22[b]  = -log((jkw[2][2] * jkw[1][1]) / (jkw[1][2] * jkw[1][2]));
        jk_c33[b]  = -log((jkw[3][3] * jkw[2][2]) / (jkw[2][3] * jkw[2][3]));
        jk_rat[b]  = (jk_sig[b] > 1e-12) ? (jk_mg12[b] / jk_sig[b]) : 0.0;
    }

    const double jk_fac = (double)(JK_BLOCKS - 1) / (double)JK_BLOCKS;
    for (int dt = 0; dt <= DT_MAX; ++dt) {
        double mr = 0.0, mg = 0.0;
        for (int b = 0; b < JK_BLOCKS; ++b) {
            mr += jk_craw[b][dt];
            mg += jk_cgevp[b][dt];
        }
        mr /= (double)JK_BLOCKS;
        mg /= (double)JK_BLOCKS;
        double vr = 0.0, vg = 0.0;
        for (int b = 0; b < JK_BLOCKS; ++b) {
            vr += (jk_craw[b][dt] - mr) * (jk_craw[b][dt] - mr);
            vg += (jk_cgevp[b][dt] - mg) * (jk_cgevp[b][dt] - mg);
        }
        res->err_c_raw[dt]  = sqrt(jk_fac * vr);
        res->err_c_gevp[dt] = sqrt(jk_fac * vg);
    }

    #define JK_ERR_1D(arr, out_field) \
        do { \
            double _m = 0.0, _v = 0.0; \
            for (int _b = 0; _b < JK_BLOCKS; ++_b) _m += (arr)[_b]; \
            _m /= (double)JK_BLOCKS; \
            for (int _b = 0; _b < JK_BLOCKS; ++_b) { \
                double _d = (arr)[_b] - _m; \
                _v += _d * _d; \
            } \
            (out_field) = sqrt(jk_fac * _v); \
        } while (0)

    JK_ERR_1D(jk_mr01, res->err_m_raw_01);
    JK_ERR_1D(jk_mg01, res->err_m_gevp_01);
    JK_ERR_1D(jk_mg12, res->err_m_gevp_12);
    JK_ERR_1D(jk_sig, res->err_a_sqrt_sigma);
    JK_ERR_1D(jk_c22, res->err_chi_sm22);
    JK_ERR_1D(jk_c33, res->err_chi_sm33);
    JK_ERR_1D(jk_rat, res->err_ratio_m_over_sig);
    res->err_am_0pp = res->err_m_gevp_12;
    if (all_v23) {
        JK_ERR_1D(jk_mg23, res->err_m_gevp_23);
    } else {
        res->valid_m_gevp_23 = 0;
        res->err_m_gevp_23   = 0.0;
    }
    #undef JK_ERR_1D

    free(p_h);
    free(absl_h);
    free(o_h);
    free(prod_h);
    free(w_h);
    wloop_ws_free(&ws);
    lat_free(&lat);
    lat_free(&sm_a);
    lat_free(&sm_b);
    return NULL;
}

int main(int argc, char **argv) {
    const char *json_path = NULL;
    for (int i = 1; i + 1 < argc; ++i) {
        if (strcmp(argv[i], "--json") == 0) {
            json_path = argv[i + 1];
        }
    }

    printf("=== 4D Pure SU(2) Glueball 0++ Mass Gap & OS Positivity Wind-Tunnel (%dx%d^3, N_meas=%d, 3 pthreads) ===\n",
           NT_MG, NS_MG, N_MEAS);
    run_startup_selftests();

    const double betas[N_BETAS] = {2.30, 2.40, 2.50};
    mg_res_t jobs[N_BETAS];
    pthread_t tids[N_BETAS];

    for (int b = 0; b < N_BETAS; ++b) {
        memset(&jobs[b], 0, sizeof(jobs[b]));
        jobs[b].idx  = b;
        jobs[b].beta = betas[b];
        if (pthread_create(&tids[b], NULL, massgap_worker, &jobs[b]) != 0) {
            fprintf(stderr, "ERROR: pthread_create failed\n");
            return 1;
        }
    }
    for (int b = 0; b < N_BETAS; ++b) {
        pthread_join(tids[b], NULL);
    }

    for (int b = 0; b < N_BETAS; ++b) {
        const mg_res_t *r = &jobs[b];
        printf("\n  --- beta = %.2f ---\n", r->beta);
        printf("    <P>=%.5f+-%.5f | <|L|>=%.5f+-%.5f | a*sqrt(sigma)=%.4f+-%.4f | chi_sm22=%.4f | chi_sm33=%.4f\n",
               r->mean_p, r->err_p, r->mean_abs_l, r->err_abs_l,
               r->a_sqrt_sigma, r->err_a_sqrt_sigma, r->chi_sm22, r->chi_sm33);
        printf("    OS eigs C(0)=[%.3e, %.3e, %.3e, %.3e] | C(1)=[%.3e, %.3e, %.3e, %.3e] | Hankel=%.4f\n",
               r->eig_c0[0], r->eig_c0[1], r->eig_c0[2], r->eig_c0[3],
               r->eig_c1[0], r->eig_c1[1], r->eig_c1[2], r->eig_c1[3], r->hankel_ratio);
        printf("    C_gevp(0..3)=[%.4e, %.4e, %.4e, %.4e]\n",
               r->c_gevp[0], r->c_gevp[1], r->c_gevp[2], r->c_gevp[3]);
        printf("    m_raw(0->1)=%.3f+-%.3f | m_gevp(0->1)=%.3f+-%.3f | m_gevp(1->2)=%.3f+-%.3f",
               r->m_raw_01, r->err_m_raw_01, r->m_gevp_01, r->err_m_gevp_01, r->am_0pp, r->err_am_0pp);
        if (r->valid_m_gevp_23) {
            printf(" | m_gevp(2->3)=%.3f+-%.3f", r->m_gevp_23, r->err_m_gevp_23);
        }
        printf(" | m_0++/sqrt(sigma)=%.2f+-%.2f\n", r->ratio_m_over_sig, r->err_ratio_m_over_sig);
    }

    if (json_path != NULL) {
        FILE *fp = fopen(json_path, "w");
        if (fp == NULL) {
            fprintf(stderr, "ERROR: cannot write %s\n", json_path);
            return 1;
        }
        fprintf(fp, "{\n");
        fprintf(fp, "  \"nt\": %d,\n", NT_MG);
        fprintf(fp, "  \"ns\": %d,\n", NS_MG);
        fprintf(fp, "  \"n_therm\": %d,\n", N_THERM);
        fprintf(fp, "  \"n_sweeps\": %d,\n", N_SWEEPS);
        fprintf(fp, "  \"n_meas\": %d,\n", N_MEAS);
        fprintf(fp, "  \"smear_levels\": [0, 2, 5, 10],\n");
        fprintf(fp, "  \"runs\": [\n");
        for (int b = 0; b < N_BETAS; ++b) {
            const mg_res_t *r = &jobs[b];
            fprintf(fp, "    {\n");
            fprintf(fp, "      \"beta\": %.2f,\n", r->beta);
            fprintf(fp, "      \"mean_p\": %.8f,\n", r->mean_p);
            fprintf(fp, "      \"err_p\": %.8f,\n", r->err_p);
            fprintf(fp, "      \"mean_abs_l\": %.8f,\n", r->mean_abs_l);
            fprintf(fp, "      \"err_abs_l\": %.8f,\n", r->err_abs_l);
            fprintf(fp, "      \"eig_c0\": [%.8e, %.8e, %.8e, %.8e],\n",
                    r->eig_c0[0], r->eig_c0[1], r->eig_c0[2], r->eig_c0[3]);
            fprintf(fp, "      \"eig_c1\": [%.8e, %.8e, %.8e, %.8e],\n",
                    r->eig_c1[0], r->eig_c1[1], r->eig_c1[2], r->eig_c1[3]);
            fprintf(fp, "      \"c_raw\": [%.8e, %.8e, %.8e, %.8e, %.8e],\n",
                    r->c_raw[0], r->c_raw[1], r->c_raw[2], r->c_raw[3], r->c_raw[4]);
            fprintf(fp, "      \"c_gevp\": [%.8e, %.8e, %.8e, %.8e, %.8e],\n",
                    r->c_gevp[0], r->c_gevp[1], r->c_gevp[2], r->c_gevp[3], r->c_gevp[4]);
            fprintf(fp, "      \"err_c_gevp\": [%.8e, %.8e, %.8e, %.8e, %.8e],\n",
                    r->err_c_gevp[0], r->err_c_gevp[1], r->err_c_gevp[2], r->err_c_gevp[3], r->err_c_gevp[4]);
            fprintf(fp, "      \"hankel_ratio\": %.6f,\n", r->hankel_ratio);
            fprintf(fp, "      \"m_raw_01\": %.6f,\n", r->m_raw_01);
            fprintf(fp, "      \"err_m_raw_01\": %.6f,\n", r->err_m_raw_01);
            fprintf(fp, "      \"m_gevp_01\": %.6f,\n", r->m_gevp_01);
            fprintf(fp, "      \"err_m_gevp_01\": %.6f,\n", r->err_m_gevp_01);
            fprintf(fp, "      \"m_gevp_12\": %.6f,\n", r->m_gevp_12);
            fprintf(fp, "      \"err_m_gevp_12\": %.6f,\n", r->err_m_gevp_12);
            fprintf(fp, "      \"m_gevp_23\": %.6f,\n", r->m_gevp_23);
            fprintf(fp, "      \"err_m_gevp_23\": %.6f,\n", r->err_m_gevp_23);
            fprintf(fp, "      \"valid_m_gevp_23\": %s,\n", r->valid_m_gevp_23 ? "true" : "false");
            fprintf(fp, "      \"am_0pp\": %.6f,\n", r->am_0pp);
            fprintf(fp, "      \"err_am_0pp\": %.6f,\n", r->err_am_0pp);
            fprintf(fp, "      \"a_sqrt_sigma\": %.6f,\n", r->a_sqrt_sigma);
            fprintf(fp, "      \"err_a_sqrt_sigma\": %.6f,\n", r->err_a_sqrt_sigma);
            fprintf(fp, "      \"chi_sm22\": %.6f,\n", r->chi_sm22);
            fprintf(fp, "      \"err_chi_sm22\": %.6f,\n", r->err_chi_sm22);
            fprintf(fp, "      \"chi_sm33\": %.6f,\n", r->chi_sm33);
            fprintf(fp, "      \"err_chi_sm33\": %.6f,\n", r->err_chi_sm33);
            fprintf(fp, "      \"ratio_m_over_sig\": %.6f,\n", r->ratio_m_over_sig);
            fprintf(fp, "      \"err_ratio_m_over_sig\": %.6f\n", r->err_ratio_m_over_sig);
            fprintf(fp, "    }%s\n", (b + 1 < N_BETAS) ? "," : "");
        }
        fprintf(fp, "  ]\n");
        fprintf(fp, "}\n");
        fclose(fp);
    }

    return 0;
}
