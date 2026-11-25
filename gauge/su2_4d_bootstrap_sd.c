/* su2_4d_bootstrap_sd.c — 4D Pure SU(2) Exact Schwinger-Dyson, OS Reflection Positivity & 7-Staple Gram Probe
 * Tests `Flow.LatticeBootstrap` Theorems B2 & B3 directly on 4D pure SU(2) Monte Carlo ensembles:
 *   - Exact 4D Schwinger-Dyson loop identity: 3*u_P == beta * [(1 - u_{P^2}) + (w_{1x2} - u_copl2) + 4*(w_chair - u_orth2)]
 *   - Exact quaternionic field-strength sum-of-squares: (1 - u_{P^2}) + (w_{1x2} - u_copl2) + 4*(w_chair - u_orth2) == |sum_s F_s|^2 / 6
 *   - Exact resultant staple norm-squared identity: R_l^2 == 6 + 6*w_{1x2} + 24*w_chair
 *   - Strict 4D OS positivity: w_{1x2} > u_copl2 >= u_P^2 and w_chair > u_orth2 >= u_P^2
 *   - Strict 7-link backtracking Gram positivity: 1 + w_{1x2} + 4*w_chair - 6*u_P^2 > 0
 */
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    L_BS        = 8,
    N_BETAS     = 5,
    N_THERM     = 400,
    N_SWEEPS    = 2000,
    MEAS_STRIDE = 2,
    N_MEAS      = N_SWEEPS / MEAS_STRIDE, /* 1000 measurements */
    JK_BLOCKS   = 20,
    N_OBS       = 11
};

static const double K_BETAS[N_BETAS] = {1.50, 2.00, 2.30, 2.40, 2.50};
static const double K_TWO_PI         = 6.28318530717958647692528676655900577;

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

/* --- 4D SU(2) Periodic Hypercubic Lattice --- */
typedef struct {
    int L;
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
    return (size_t)(((t * lat->L + x) * lat->L + y) * lat->L + z);
}

static inline size_t link_idx(const su2_lat_t *lat, int t, int x, int y, int z, int mu) {
    return 4U * site_idx(lat, t, x, y, z) + (size_t)mu;
}

static void lat_alloc(su2_lat_t *lat, int L) {
    lat->L   = L;
    lat->vol = L * L * L * L;
    lat->u   = (quat_t *)malloc(4U * (size_t)lat->vol * sizeof(quat_t));
    if (lat->u == NULL) {
        fprintf(stderr, "ERROR: OOM allocating %d^4 lattice\n", L);
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
    out[0] = c[0];
    out[1] = c[1];
    out[2] = c[2];
    out[3] = c[3];
    out[d] = (dir > 0) ? next_i(c[d], lat->L) : prev_i(c[d], lat->L);
}

static inline quat_t get_link(const su2_lat_t *lat, const int c[4], int mu) {
    return lat->u[link_idx(lat, c[0], c[1], c[2], c[3], mu)];
}

static void compute_6_staples(const su2_lat_t *lat, const int c[4], int mu, quat_t staples[6]) {
    int c_mu[4];
    step_coord(lat, c, mu, +1, c_mu);
    int k = 0;
    for (int nu = 0; nu < 4; ++nu) {
        if (nu == mu) {
            continue;
        }
        int c_nu[4], c_mnu[4], c_mnumu[4];
        step_coord(lat, c, nu, +1, c_nu);
        step_coord(lat, c, nu, -1, c_mnu);
        step_coord(lat, c_mnu, mu, +1, c_mnumu);

        staples[2 * k + 0] = q_mul(q_mul(get_link(lat, c, nu), get_link(lat, c_nu, mu)),
                                   q_conj(get_link(lat, c_mu, nu)));
        staples[2 * k + 1] = q_mul(q_mul(q_conj(get_link(lat, c_mnu, nu)), get_link(lat, c_mnu, mu)),
                                   get_link(lat, c_mnumu, nu));
        ++k;
    }
}

static quat_t compute_staple_4d(const su2_lat_t *lat, int t, int x, int y, int z, int mu) {
    const int c[4] = {t, x, y, z};
    quat_t v6[6];
    compute_6_staples(lat, c, mu, v6);
    quat_t sum = {0.0, 0.0, 0.0, 0.0};
    for (int s = 0; s < 6; ++s) {
        sum.w += v6[s].w;
        sum.x += v6[s].x;
        sum.y += v6[s].y;
        sum.z += v6[s].z;
    }
    return sum;
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

static double cfg_total_action(const su2_lat_t *lat, double beta) {
    double s = 0.0;
    for (int t = 0; t < lat->L; ++t) {
        for (int x = 0; x < lat->L; ++x) {
            for (int y = 0; y < lat->L; ++y) {
                for (int z = 0; z < lat->L; ++z) {
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

/* --- CKP SU(2) Heatbath + 2x Adler Checkerboard OR --- */
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

static void sweep_su2_or_only(su2_lat_t *lat) {
    for (int par = 0; par < 2; ++par) {
        for (int t = 0; t < lat->L; ++t) {
            for (int x = 0; x < lat->L; ++x) {
                for (int y = 0; y < lat->L; ++y) {
                    for (int z = 0; z < lat->L; ++z) {
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

static void sweep_su2_hb_or(su2_lat_t *lat, double beta, rng_t *rng) {
    for (int par = 0; par < 2; ++par) {
        for (int t = 0; t < lat->L; ++t) {
            for (int x = 0; x < lat->L; ++x) {
                for (int y = 0; y < lat->L; ++y) {
                    for (int z = 0; z < lat->L; ++z) {
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
    sweep_su2_or_only(lat);
    sweep_su2_or_only(lat);
}

/* --- Link-wise 6-Staple SD & Gram Measurement Kernel --- */
typedef struct {
    double v[N_OBS];
    /* 0: u_P
     * 1: u_P2
     * 2: w_1x2
     * 3: u_copl2
     * 4: w_chair
     * 5: u_orth2
     * 6: sum_y_sq_over6  (|sum_s F_s|^2 / 6)
     * 7: r_sq            (|sum_s V_s|^2)
     * 8: sd_res          (3*u_P - beta*((1 - u_P2) + (w_1x2 - u_copl2) + 4*(w_chair - u_orth2)))
     * 9: gram_diff       (((1 - u_P2) + (w_1x2 - u_copl2) + 4*(w_chair - u_orth2)) - sum_y_sq_over6)
     * 10: abs_l          (<|L|>)
     */
} sd_obs_t;

static sd_obs_t measure_sd_observables(const su2_lat_t *lat, double beta) {
    double sum_up       = 0.0;
    double sum_up2      = 0.0;
    double sum_w12      = 0.0;
    double sum_uco2     = 0.0;
    double sum_wch      = 0.0;
    double sum_uor2     = 0.0;
    double sum_ysq_o6   = 0.0;
    double sum_rsq      = 0.0;

    const double n_links = (double)(4 * lat->vol);

    for (int t = 0; t < lat->L; ++t) {
        for (int x = 0; x < lat->L; ++x) {
            for (int y = 0; y < lat->L; ++y) {
                for (int z = 0; z < lat->L; ++z) {
                    const int c[4] = {t, x, y, z};
                    for (int mu = 0; mu < 4; ++mu) {
                        quat_t ul = get_link(lat, c, mu);
                        quat_t v6[6];
                        compute_6_staples(lat, c, mu, v6);

                        quat_t q6[6];
                        double sx = 0.0, sy = 0.0, sz = 0.0;
                        double rw = 0.0, rx = 0.0, ry = 0.0, rz = 0.0;
                        for (int s = 0; s < 6; ++s) {
                            q6[s] = q_mul(ul, q_conj(v6[s]));
                            sum_up  += q6[s].w;
                            sum_up2 += q6[s].w * q6[s].w;
                            sx += q6[s].x;
                            sy += q6[s].y;
                            sz += q6[s].z;
                            rw += v6[s].w;
                            rx += v6[s].x;
                            ry += v6[s].y;
                            rz += v6[s].z;
                        }
                        sum_ysq_o6 += (sx * sx + sy * sy + sz * sz) / 6.0;
                        sum_rsq    += rw * rw + rx * rx + ry * ry + rz * rz;

                        /* 3 coplanar opposite pairs (0,1), (2,3), (4,5) */
                        for (int k = 0; k < 3; ++k) {
                            int s = 2 * k;
                            int p = 2 * k + 1;
                            sum_w12  += q_inner(v6[s], v6[p]);
                            sum_uco2 += q6[s].w * q6[p].w;
                        }

                        /* 12 unordered orthogonal pairs (2*k1+a, 2*k2+b) for k1 < k2 */
                        for (int k1 = 0; k1 < 3; ++k1) {
                            for (int k2 = k1 + 1; k2 < 3; ++k2) {
                                for (int a = 0; a < 2; ++a) {
                                    for (int b = 0; b < 2; ++b) {
                                        int s = 2 * k1 + a;
                                        int p = 2 * k2 + b;
                                        sum_wch  += q_inner(v6[s], v6[p]);
                                        sum_uor2 += q6[s].w * q6[p].w;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    double sum_l   = 0.0;
    int space_vol  = lat->L * lat->L * lat->L;
    for (int x = 0; x < lat->L; ++x) {
        for (int y = 0; y < lat->L; ++y) {
            for (int z = 0; z < lat->L; ++z) {
                quat_t pl = q_id();
                for (int t = 0; t < lat->L; ++t) {
                    pl = q_mul(pl, lat->u[link_idx(lat, t, x, y, z, 0)]);
                }
                sum_l += pl.w;
            }
        }
    }

    sd_obs_t out;
    double up    = sum_up / (6.0 * n_links);
    double up2   = sum_up2 / (6.0 * n_links);
    double w12   = sum_w12 / (3.0 * n_links);
    double uco2  = sum_uco2 / (3.0 * n_links);
    double wch   = sum_wch / (12.0 * n_links);
    double uor2  = sum_uor2 / (12.0 * n_links);
    double ysq6  = sum_ysq_o6 / n_links;
    double rsq   = sum_rsq / n_links;
    double sd_op = (1.0 - up2) + (w12 - uco2) + 4.0 * (wch - uor2);

    out.v[0]  = up;
    out.v[1]  = up2;
    out.v[2]  = w12;
    out.v[3]  = uco2;
    out.v[4]  = wch;
    out.v[5]  = uor2;
    out.v[6]  = ysq6;
    out.v[7]  = rsq;
    out.v[8]  = 3.0 * up - beta * sd_op;
    out.v[9]  = sd_op - ysq6;
    out.v[10] = fabs(sum_l / (double)space_vol);
    return out;
}

/* --- Exact Local Gauge & Quaternion Gram Startup Self-Tests --- */
static void run_selftests(void) {
    rng_t rng;
    rng_seed(&rng, 0xB00757A9504DULL);
    su2_lat_t lat, gauged;
    lat_alloc(&lat, 4);
    lat_alloc(&gauged, 4);
    lat_init_cold(&lat);
    for (int sw = 0; sw < 15; ++sw) {
        sweep_su2_hb_or(&lat, 2.30, &rng);
    }

    sd_obs_t o1 = measure_sd_observables(&lat, 2.30);

    /* 1. Random local SU(2) gauge transformation g(x) */
    quat_t *g = (quat_t *)malloc((size_t)lat.vol * sizeof(quat_t));
    if (!g) { fprintf(stderr, "ERROR: OOM allocating gauge field\n"); exit(1); }
    for (int i = 0; i < lat.vol; ++i) {
        g[i] = sample_su2_boltzmann(0.0, &rng);
    }
    for (int t = 0; t < lat.L; ++t) {
        for (int x = 0; x < lat.L; ++x) {
            for (int y = 0; y < lat.L; ++y) {
                for (int z = 0; z < lat.L; ++z) {
                    const int c[4] = {t, x, y, z};
                    size_t sid     = site_idx(&lat, t, x, y, z);
                    for (int mu = 0; mu < 4; ++mu) {
                        int c_mu[4];
                        step_coord(&lat, c, mu, +1, c_mu);
                        size_t sid_mu = site_idx(&lat, c_mu[0], c_mu[1], c_mu[2], c_mu[3]);
                        quat_t u      = get_link(&lat, c, mu);
                        quat_t ug     = q_mul(q_mul(g[sid], u), q_conj(g[sid_mu]));
                        q_normalize(&ug);
                        gauged.u[link_idx(&gauged, t, x, y, z, mu)] = ug;
                    }
                }
            }
        }
    }
    sd_obs_t o2    = measure_sd_observables(&gauged, 2.30);
    double max_dg  = 0.0;
    for (int k = 0; k < N_OBS; ++k) {
        double d = fabs(o1.v[k] - o2.v[k]);
        if (d > max_dg) {
            max_dg = d;
        }
    }

    /* 2. Exact algebraic Gram & resultant norm identities */
    double d_gram = fabs(o1.v[9]);
    double d_rsq  = fabs(o1.v[7] - (6.0 + 6.0 * o1.v[2] + 24.0 * o1.v[4]));

    /* 3. Microcanonical OR check */
    lat_copy(&gauged, &lat);
    double s_before = cfg_total_action(&gauged, 2.30);
    sweep_su2_or_only(&gauged);
    double s_after = cfg_total_action(&gauged, 2.30);
    double d_or    = fabs(s_after - s_before);

    free(g);
    lat_free(&lat);
    lat_free(&gauged);

    if (max_dg > 1e-13 || d_gram > 1e-14 || d_rsq > 1e-13 || d_or > 1e-9) {
        fprintf(stderr, "[SELFTEST] FAIL max_dg=%.3e d_gram=%.3e d_rsq=%.3e d_or=%.3e\n",
                max_dg, d_gram, d_rsq, d_or);
        exit(1);
    }
    printf("  [SELFTEST] PASS (gauge <%.2e, Gram id=%.2e, R^2 id=%.2e, OR dS=%.2e)\n",
           max_dg, d_gram, d_rsq, d_or);
    fflush(stdout);
}

/* --- Per-Beta Multi-Threaded Simulation --- */
typedef struct {
    int idx;
    double beta;
    double mean[N_OBS];
    double err[N_OBS];
} bs_res_t;

static void *run_beta_worker(void *arg) {
    bs_res_t *res = (bs_res_t *)arg;
    rng_t rng;
    rng_seed(&rng, 0x5D8B0070001ULL + (uint64_t)(res->idx + 1) * 0x9E3779B97F4A7C15ULL);

    su2_lat_t lat;
    lat_alloc(&lat, L_BS);
    lat_init_cold(&lat);

    for (int sw = 0; sw < N_THERM; ++sw) {
        sweep_su2_hb_or(&lat, res->beta, &rng);
    }

    double blk_sum[JK_BLOCKS][N_OBS];
    memset(blk_sum, 0, sizeof(blk_sum));
    double tot_sum[N_OBS];
    memset(tot_sum, 0, sizeof(tot_sum));
    const int blk_len = N_MEAS / JK_BLOCKS;

    for (int meas = 0; meas < N_MEAS; ++meas) {
        for (int s = 0; s < MEAS_STRIDE; ++s) {
            sweep_su2_hb_or(&lat, res->beta, &rng);
        }
        sd_obs_t ob = measure_sd_observables(&lat, res->beta);
        int b       = meas / blk_len;
        if (b >= JK_BLOCKS) {
            b = JK_BLOCKS - 1;
        }
        for (int k = 0; k < N_OBS; ++k) {
            blk_sum[b][k] += ob.v[k];
            tot_sum[k] += ob.v[k];
        }
    }

    const double nb_n = (double)(N_MEAS - blk_len);
    const double jf   = (double)(JK_BLOCKS - 1) / (double)JK_BLOCKS;
    for (int k = 0; k < N_OBS; ++k) {
        res->mean[k] = tot_sum[k] / (double)N_MEAS;
        double ssq   = 0.0;
        for (int b = 0; b < JK_BLOCKS; ++b) {
            double jk_val = (tot_sum[k] - blk_sum[b][k]) / nb_n;
            double diff   = jk_val - res->mean[k];
            ssq += diff * diff;
        }
        res->err[k] = sqrt(jf * ssq);
    }

    lat_free(&lat);
    return NULL;
}

int main(int argc, char **argv) {
    const char *json_path = NULL;
    for (int i = 1; i + 1 < argc; ++i) {
        if (strcmp(argv[i], "--json") == 0) {
            json_path = argv[i + 1];
        }
    }

    printf("=== 4D Pure SU(2) Schwinger-Dyson, OS & 7-Staple Gram Probe (%d^4, N_meas=%d, %d pthreads) ===\n",
           L_BS, N_MEAS, N_BETAS);
    run_selftests();

    bs_res_t results[N_BETAS];
    pthread_t tids[N_BETAS];
    for (int i = 0; i < N_BETAS; ++i) {
        memset(&results[i], 0, sizeof(results[i]));
        results[i].idx  = i;
        results[i].beta = K_BETAS[i];
        pthread_create(&tids[i], NULL, run_beta_worker, &results[i]);
    }
    for (int i = 0; i < N_BETAS; ++i) {
        pthread_join(tids[i], NULL);
    }

    for (int i = 0; i < N_BETAS; ++i) {
        const bs_res_t *r = &results[i];
        double gram7      = 1.0 + r->mean[2] + 4.0 * r->mean[4] - 6.0 * r->mean[0] * r->mean[0];
        printf("  beta=%.2f | u_P=%.5f+-%.5f | u_P2=%.5f | w_1x2=%.5f (co2=%.5f) | w_chair=%.5f (or2=%.5f) "
               "| sd_res=%+.5f+-%.5f | Gram7=%.5f | <|L|>=%.4f\n",
               r->beta, r->mean[0], r->err[0], r->mean[1], r->mean[2], r->mean[3],
               r->mean[4], r->mean[5], r->mean[8], r->err[8], gram7, r->mean[10]);
    }

    if (json_path != NULL) {
        FILE *fp = fopen(json_path, "w");
        if (fp == NULL) {
            fprintf(stderr, "ERROR: failed to open %s for writing\n", json_path);
            return 1;
        }
        fprintf(fp, "{\n");
        fprintf(fp, "  \"L\": %d,\n", L_BS);
        fprintf(fp, "  \"n_therm\": %d,\n", N_THERM);
        fprintf(fp, "  \"n_sweeps\": %d,\n", N_SWEEPS);
        fprintf(fp, "  \"n_meas\": %d,\n", N_MEAS);
        fprintf(fp, "  \"runs\": [\n");
        for (int i = 0; i < N_BETAS; ++i) {
            const bs_res_t *r = &results[i];
            fprintf(fp, "    {\n");
            fprintf(fp, "      \"beta\": %.2f,\n", r->beta);
            fprintf(fp, "      \"u_p\": %.14f,\n", r->mean[0]);
            fprintf(fp, "      \"err_u_p\": %.8f,\n", r->err[0]);
            fprintf(fp, "      \"u_p2\": %.14f,\n", r->mean[1]);
            fprintf(fp, "      \"err_u_p2\": %.8f,\n", r->err[1]);
            fprintf(fp, "      \"w_1x2\": %.14f,\n", r->mean[2]);
            fprintf(fp, "      \"err_w_1x2\": %.8f,\n", r->err[2]);
            fprintf(fp, "      \"u_copl2\": %.14f,\n", r->mean[3]);
            fprintf(fp, "      \"err_u_copl2\": %.8f,\n", r->err[3]);
            fprintf(fp, "      \"w_chair\": %.14f,\n", r->mean[4]);
            fprintf(fp, "      \"err_w_chair\": %.8f,\n", r->err[4]);
            fprintf(fp, "      \"u_orth2\": %.14f,\n", r->mean[5]);
            fprintf(fp, "      \"err_u_orth2\": %.8f,\n", r->err[5]);
            fprintf(fp, "      \"sum_y_sq_over6\": %.14f,\n", r->mean[6]);
            fprintf(fp, "      \"r_sq\": %.14f,\n", r->mean[7]);
            fprintf(fp, "      \"sd_res\": %.8f,\n", r->mean[8]);
            fprintf(fp, "      \"err_sd_res\": %.8f,\n", r->err[8]);
            fprintf(fp, "      \"gram_diff\": %.12e,\n", r->mean[9]);
            fprintf(fp, "      \"abs_l\": %.8f\n", r->mean[10]);
            fprintf(fp, "    }%s\n", (i + 1 < N_BETAS) ? "," : "");
        }
        fprintf(fp, "  ]\n");
        fprintf(fp, "}\n");
        fclose(fp);
    }

    return 0;
}
