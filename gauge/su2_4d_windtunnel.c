/* su2_4d_windtunnel.c — 4D Pure SU(2) Yang-Mills Wind Tunnel (Theorems Y1–Y5)
 * Tests Theorems Y1–Y5 from `Flow.YangMillsBundle`:
 *   - 6 Startup Mathematical & Physics Self-Tests (4D gauge invariance, Z_2 center flip, OR exactness,
 *     fast 4D defect-wall @swap, Creutz-Kennedy-Pendleton exactness vs I_2/I_1, BPST Q_cool = +-1 + roundtrip)
 *   - Experiment 1 (4 x 12^3, beta in {2.15, 2.25, 2.30, 2.35, 2.45}):
 *       Thermal Z_2 deconfinement & symbiotic_z2_center ergodicity restoration
 *   - Experiment 2 (10^4, beta in {2.25, 2.35, 2.45, 2.52}, eq_work = 2,400 replica-sweeps):
 *       Confinement scaling chi(2,2), chi(3,3), topological susceptibility chi_t^(1/4)/sqrt(chi22),
 *       Corollary 2B on 4D BPST instantons (bare_bpst_inst), and 4D defect-wall PTBC + Fisher 2-Jet ladder
 */
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    /* Experiment 1: Thermal Z_2 Deconfinement */
    NT_THERM     = 4,
    NS_THERM     = 12,
    N_BETAS_Z2   = 5,
    N_THERM_Z2   = 400,
    N_MEAS_Z2    = 1600,

    /* Experiment 2: Confined 4D SU(2) Confinement, Topology & Defect-Wall PTBC */
    L_CONF       = 10,
    L_DEFECT_4D  = 2,
    R_REPLICAS   = 8,
    N_BETAS_CONF = 4,
    N_EQ_MEAS    = 2400,
    N_THERM_1R   = 400,
    N_THERM_PT   = 200,
    N_MEAS_PT    = N_EQ_MEAS / R_REPLICAS, /* 300 -> exact equal work 8 * 300 = 2400 */
    COOL_STEPS   = 12,
    JK_BLOCKS    = 20
};

static const double K_PI     = 3.14159265358979323846264338327950288;
static const double K_TWO_PI = 6.28318530717958647692528676655900577;

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

static inline int rng_below(rng_t *r, int bound) {
    return (int)(((uint64_t)(uint32_t)rng_u64(r) * (uint64_t)(uint32_t)bound) >> 32);
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

static inline quat_t q_neg(quat_t a) {
    quat_t r = {-a.w, -a.x, -a.y, -a.z};
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

static inline quat_t q_exp_im(double vx, double vy, double vz) {
    double n = sqrt(vx * vx + vy * vy + vz * vz);
    if (n < 1e-14) {
        return q_id();
    }
    double s = sin(n) / n;
    quat_t r = {cos(n), vx * s, vy * s, vz * s};
    return r;
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
        fprintf(stderr, "ERROR: OOM allocating 4D SU(2) lattice %dx%d^3\n", nt, ns);
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

/* Step a coordinate array c[4] by +1 or -1 in direction d */
static inline void step_coord(const su2_lat_t *lat, const int c[4], int d, int dir, int out[4]) {
    out[0]     = c[0];
    out[1]     = c[1];
    out[2]     = c[2];
    out[3]     = c[3];
    int lim    = (d == 0) ? lat->nt : lat->ns;
    out[d]     = (dir > 0) ? next_i(c[d], lim) : prev_i(c[d], lim);
}

static inline quat_t get_link(const su2_lat_t *lat, const int c[4], int mu) {
    return lat->u[link_idx(lat, c[0], c[1], c[2], c[3], mu)];
}

static inline int is_defect_site(int l_def, int x, int y, int z) {
    return (x >= 0 && x < l_def) && (y >= 0 && y < l_def) && (z >= 0 && z < l_def);
}

static inline int is_defect_plaq(const su2_lat_t *lat, int l_def, const int c[4], int mu_lo, int nu_hi) {
    if (l_def <= 0 || mu_lo != 0 || c[0] != (lat->nt / 2)) {
        return 0;
    }
    if (is_defect_site(l_def, c[1], c[2], c[3])) {
        return 1;
    }
    int c2[4];
    step_coord(lat, c, nu_hi, +1, c2);
    return is_defect_site(l_def, c2[1], c2[2], c2[3]);
}

static inline double plaq_weight(const su2_lat_t *lat, int l_def, double c_r,
                                 const int c[4], int mu, int nu) {
    int mu_lo = (mu < nu) ? mu : nu;
    int nu_hi = (mu < nu) ? nu : mu;
    return is_defect_plaq(lat, l_def, c, mu_lo, nu_hi) ? c_r : 1.0;
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

/* Staple V = sum_{nu != mu} [ w_fwd * (S^+)^dag + w_bwd * (S^-)^dag ] such that Re(U * conj(V)) = sum w_P Re(P) */
static quat_t compute_staple(const su2_lat_t *lat, int t, int x, int y, int z, int mu,
                             double c_r, int l_def) {
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

        double wf  = plaq_weight(lat, l_def, c_r, c, mu, nu);
        quat_t sf  = q_mul(q_mul(get_link(lat, c, nu), get_link(lat, c_nu, mu)),
                           q_conj(get_link(lat, c_mu, nu)));

        double wb  = plaq_weight(lat, l_def, c_r, c_mnu, mu, nu);
        quat_t sb  = q_mul(q_mul(q_conj(get_link(lat, c_mnu, nu)), get_link(lat, c_mnu, mu)),
                           get_link(lat, c_mnumu, nu));

        sum.w += wf * sf.w + wb * sb.w;
        sum.x += wf * sf.x + wb * sb.x;
        sum.y += wf * sf.y + wb * sb.y;
        sum.z += wf * sf.z + wb * sb.z;
    }
    return sum;
}

static double cfg_total_action(const su2_lat_t *lat, double beta, double c_r, int l_def) {
    double s = 0.0;
    for (int t = 0; t < lat->nt; ++t) {
        for (int x = 0; x < lat->ns; ++x) {
            for (int y = 0; y < lat->ns; ++y) {
                for (int z = 0; z < lat->ns; ++z) {
                    const int c[4] = {t, x, y, z};
                    for (int mu = 0; mu < 4; ++mu) {
                        for (int nu = mu + 1; nu < 4; ++nu) {
                            double w = plaq_weight(lat, l_def, c_r, c, mu, nu);
                            quat_t p = plaq_quat(lat, c, mu, nu);
                            s += w * (1.0 - p.w);
                        }
                    }
                }
            }
        }
    }
    return beta * s;
}

static double cfg_defect_sum_re(const su2_lat_t *lat, int l_def, int *out_cnt) {
    double s = 0.0;
    int cnt  = 0;
    int t    = lat->nt / 2;
    for (int x = 0; x < lat->ns; ++x) {
        for (int y = 0; y < lat->ns; ++y) {
            for (int z = 0; z < lat->ns; ++z) {
                const int c[4] = {t, x, y, z};
                for (int nu = 1; nu < 4; ++nu) {
                    if (is_defect_plaq(lat, l_def, c, 0, nu)) {
                        s += plaq_quat(lat, c, 0, nu).w;
                        ++cnt;
                    }
                }
            }
        }
    }
    if (out_cnt != NULL) {
        *out_cnt = cnt;
    }
    return s;
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

static double cfg_polyakov_mean(const su2_lat_t *lat, double *out_abs) {
    double sum_l     = 0.0;
    int space_sites  = lat->ns * lat->ns * lat->ns;
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
    double l = sum_l / (double)space_sites;
    if (out_abs != NULL) {
        *out_abs = fabs(l);
    }
    return l;
}

static void cfg_flip_z2_slice(su2_lat_t *lat, int t0) {
    for (int x = 0; x < lat->ns; ++x) {
        for (int y = 0; y < lat->ns; ++y) {
            for (int z = 0; z < lat->ns; ++z) {
                size_t id  = link_idx(lat, t0, x, y, z, 0);
                lat->u[id] = q_neg(lat->u[id]);
            }
        }
    }
}

/* Planar Wilson loop <1/2 Tr W(R, T)> averaged over all 6 planes, both orientations (when R != T), and all sites */
static double cfg_wilson_loop(const su2_lat_t *lat, int r_len, int t_len) {
    double sum_w = 0.0;
    int planes   = 0;
    int same     = (r_len == t_len);
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = mu + 1; nu < 4; ++nu) {
            for (int ori = 0; ori < (same ? 1 : 2); ++ori) {
                int d1 = (ori == 0) ? mu : nu;
                int d2 = (ori == 0) ? nu : mu;
                for (int t = 0; t < lat->nt; ++t) {
                    for (int x = 0; x < lat->ns; ++x) {
                        for (int y = 0; y < lat->ns; ++y) {
                            for (int z = 0; z < lat->ns; ++z) {
                                int c[4] = {t, x, y, z};
                                int c2[4];
                                quat_t w = q_id();
                                for (int i = 0; i < r_len; ++i) {
                                    w = q_mul(w, get_link(lat, c, d1));
                                    step_coord(lat, c, d1, +1, c2);
                                    memcpy(c, c2, sizeof(c));
                                }
                                for (int j = 0; j < t_len; ++j) {
                                    w = q_mul(w, get_link(lat, c, d2));
                                    step_coord(lat, c, d2, +1, c2);
                                    memcpy(c, c2, sizeof(c));
                                }
                                for (int i = 0; i < r_len; ++i) {
                                    step_coord(lat, c, d1, -1, c2);
                                    memcpy(c, c2, sizeof(c));
                                    w = q_mul(w, q_conj(get_link(lat, c, d1)));
                                }
                                for (int j = 0; j < t_len; ++j) {
                                    step_coord(lat, c, d2, -1, c2);
                                    memcpy(c, c2, sizeof(c));
                                    w = q_mul(w, q_conj(get_link(lat, c, d2)));
                                }
                                sum_w += w.w;
                            }
                        }
                    }
                }
                ++planes;
            }
        }
    }
    return sum_w / (double)(planes * lat->vol);
}

/* 4D gauge-covariant Clover Topological Charge Q = (1 / (2*pi^2)) sum_x (f01.f23 - f02.f13 + f03.f12) */
static void clover_im_plaq(const su2_lat_t *lat, const int c[4], int mu, int nu, double out_f[3]) {
    int cp_mu[4], cp_nu[4], cm_mu[4], cm_nu[4], cm_mu_p_nu[4], cm_mu_m_nu[4], cp_mu_m_nu[4];
    step_coord(lat, c, mu, +1, cp_mu);
    step_coord(lat, c, nu, +1, cp_nu);
    step_coord(lat, c, mu, -1, cm_mu);
    step_coord(lat, c, nu, -1, cm_nu);
    step_coord(lat, cm_mu, nu, +1, cm_mu_p_nu);
    step_coord(lat, cm_mu, nu, -1, cm_mu_m_nu);
    step_coord(lat, cp_mu, nu, -1, cp_mu_m_nu);

    /* Leaf 1 (+mu, +nu): U_mu(x) * U_nu(x+mu) * conj(U_mu(x+nu)) * conj(U_nu(x)) */
    quat_t l1 = q_mul(q_mul(get_link(lat, c, mu), get_link(lat, cp_mu, nu)),
                      q_mul(q_conj(get_link(lat, cp_nu, mu)), q_conj(get_link(lat, c, nu))));

    /* Leaf 2 (+nu, -mu): U_nu(x) * conj(U_mu(x-mu+nu)) * conj(U_nu(x-mu)) * U_mu(x-mu) */
    quat_t l2 = q_mul(q_mul(get_link(lat, c, nu), q_conj(get_link(lat, cm_mu_p_nu, mu))),
                      q_mul(q_conj(get_link(lat, cm_mu, nu)), get_link(lat, cm_mu, mu)));

    /* Leaf 3 (-mu, -nu): conj(U_mu(x-mu)) * conj(U_nu(x-mu-nu)) * U_mu(x-mu-nu) * U_nu(x-nu) */
    quat_t l3 = q_mul(q_mul(q_conj(get_link(lat, cm_mu, mu)), q_conj(get_link(lat, cm_mu_m_nu, nu))),
                      q_mul(get_link(lat, cm_mu_m_nu, mu), get_link(lat, cm_nu, nu)));

    /* Leaf 4 (-nu, +mu): conj(U_nu(x-nu)) * U_mu(x-nu) * U_nu(x+mu-nu) * conj(U_mu(x)) */
    quat_t l4 = q_mul(q_mul(q_conj(get_link(lat, cm_nu, nu)), get_link(lat, cm_nu, mu)),
                      q_mul(get_link(lat, cp_mu_m_nu, nu), q_conj(get_link(lat, c, mu))));

    out_f[0] = 0.25 * (l1.x + l2.x + l3.x + l4.x);
    out_f[1] = 0.25 * (l1.y + l2.y + l3.y + l4.y);
    out_f[2] = 0.25 * (l1.z + l2.z + l3.z + l4.z);
}

static double cfg_topo_charge_clover(const su2_lat_t *lat) {
    double sum_q      = 0.0;
    const double norm = 1.0 / (2.0 * K_PI * K_PI);
    for (int t = 0; t < lat->nt; ++t) {
        for (int x = 0; x < lat->ns; ++x) {
            for (int y = 0; y < lat->ns; ++y) {
                for (int z = 0; z < lat->ns; ++z) {
                    const int c[4] = {t, x, y, z};
                    double f01[3], f23[3], f02[3], f13[3], f03[3], f12[3];
                    clover_im_plaq(lat, c, 0, 1, f01);
                    clover_im_plaq(lat, c, 2, 3, f23);
                    clover_im_plaq(lat, c, 0, 2, f02);
                    clover_im_plaq(lat, c, 1, 3, f13);
                    clover_im_plaq(lat, c, 0, 3, f03);
                    clover_im_plaq(lat, c, 1, 2, f12);
                    for (int k = 0; k < 3; ++k) {
                        sum_q += f01[k] * f23[k] - f02[k] * f13[k] + f03[k] * f12[k];
                    }
                }
            }
        }
    }
    return sum_q * norm;
}

static void cfg_cool_sweeps(const su2_lat_t *src, su2_lat_t *wrk, int steps) {
    const double eps = 0.35; /* damped APE/gradient-flow cooling prevents O(a^2) instanton dislocation */
    lat_copy(wrk, src);
    for (int sw = 0; sw < steps; ++sw) {
        for (int par = 0; par < 2; ++par) {
            for (int t = 0; t < wrk->nt; ++t) {
                for (int x = 0; x < wrk->ns; ++x) {
                    for (int y = 0; y < wrk->ns; ++y) {
                        for (int z = 0; z < wrk->ns; ++z) {
                            if (((t + x + y + z) & 1) != par) {
                                continue;
                            }
                            for (int mu = 0; mu < 4; ++mu) {
                                quat_t v  = compute_staple(wrk, t, x, y, z, mu, 1.0, 0);
                                double n2 = q_norm2(v);
                                if (n2 > 1e-24) {
                                    q_normalize(&v);
                                    size_t id = link_idx(wrk, t, x, y, z, mu);
                                    quat_t u0 = wrk->u[id];
                                    quat_t u1 = {
                                        (1.0 - eps) * u0.w + eps * v.w,
                                        (1.0 - eps) * u0.x + eps * v.x,
                                        (1.0 - eps) * u0.y + eps * v.y,
                                        (1.0 - eps) * u0.z + eps * v.z
                                    };
                                    q_normalize(&u1);
                                    wrk->u[id] = u1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

static long long cfg_topo_charge_cooled(const su2_lat_t *src, su2_lat_t *wrk, double *out_raw) {
    cfg_cool_sweeps(src, wrk, COOL_STEPS);
    double q = cfg_topo_charge_clover(wrk);
    if (out_raw != NULL) {
        *out_raw = q;
    }
    return llround(q);
}

/* --- Creutz-Kennedy-Pendleton Exact SU(2) Heatbath & Microcanonical Over-Relaxation --- */
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

static void sweep_su2_hb_or(su2_lat_t *lat, double beta, double c_r, int l_def, rng_t *rng) {
    /* 1. Exact Checkerboard Creutz-Kennedy-Pendleton Heatbath */
    for (int par = 0; par < 2; ++par) {
        for (int t = 0; t < lat->nt; ++t) {
            for (int x = 0; x < lat->ns; ++x) {
                for (int y = 0; y < lat->ns; ++y) {
                    for (int z = 0; z < lat->ns; ++z) {
                        if (((t + x + y + z) & 1) != par) {
                            continue;
                        }
                        for (int mu = 0; mu < 4; ++mu) {
                            quat_t v   = compute_staple(lat, t, x, y, z, mu, c_r, l_def);
                            double k   = sqrt(q_norm2(v));
                            quat_t xq  = sample_su2_boltzmann(beta * k, rng);
                            size_t id  = link_idx(lat, t, x, y, z, mu);
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

    /* 2. Two Microcanonical Checkerboard Over-Relaxation Sweeps (Theorem Y3: dS = 0, |U|^2 = 1) */
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
                                quat_t v   = compute_staple(lat, t, x, y, z, mu, c_r, l_def);
                                double v2  = q_norm2(v);
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

static void shift_torus_4d_random(su2_lat_t *lat, su2_lat_t *tmp, rng_t *rng) {
    int st = rng_below(rng, lat->nt);
    int sx = rng_below(rng, lat->ns);
    int sy = rng_below(rng, lat->ns);
    int sz = rng_below(rng, lat->ns);
    if (st == 0 && sx == 0 && sy == 0 && sz == 0) {
        return;
    }
    lat_copy(tmp, lat);
    for (int t = 0; t < lat->nt; ++t) {
        int nt_i = (t + st) % lat->nt;
        for (int x = 0; x < lat->ns; ++x) {
            int nx_i = (x + sx) % lat->ns;
            for (int y = 0; y < lat->ns; ++y) {
                int ny_i = (y + sy) % lat->ns;
                for (int z = 0; z < lat->ns; ++z) {
                    int nz_i = (z + sz) % lat->ns;
                    for (int mu = 0; mu < 4; ++mu) {
                        lat->u[link_idx(lat, nt_i, nx_i, ny_i, nz_i, mu)] =
                            tmp->u[link_idx(tmp, t, x, y, z, mu)];
                    }
                }
            }
        }
    }
}

/* --- Singular-Gauge Smoothstep 4D BPST Instanton Twist (Exact Unit-Jacobian Involution) --- */
static inline quat_t g_bpst_at(const double c[4], int sign) {
    double r2  = c[0] * c[0] + c[1] * c[1] + c[2] * c[2] + c[3] * c[3];
    double inv = 1.0 / sqrt(r2);
    double s   = (double)sign;
    quat_t r   = {c[0] * inv, -s * c[1] * inv, -s * c[2] * inv, -s * c[3] * inv};
    return r;
}

static inline int torus_delta(int d, int lim) {
    if (d >= lim / 2) {
        d -= lim;
    }
    if (d < -lim / 2) {
        d += lim;
    }
    return d;
}

static void apply_bpst_instanton(su2_lat_t *dst, const su2_lat_t *src,
                                 int t0, int x0, int y0, int z0,
                                 double rho, int sign, int inv_dir) {
    const double r_in  = 3.2;
    const double r_out = 4.6;
    const quat_t basis[4] = {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 1.0}
    };
    const int centers[4] = {t0, x0, y0, z0};

    lat_copy(dst, src);
    for (int t = 0; t < src->nt; ++t) {
        for (int x = 0; x < src->ns; ++x) {
            for (int y = 0; y < src->ns; ++y) {
                for (int z = 0; z < src->ns; ++z) {
                    const int idx4[4] = {t, x, y, z};
                    double c0[4];
                    for (int k = 0; k < 4; ++k) {
                        int lim = (k == 0) ? src->nt : src->ns;
                        c0[k]   = (double)torus_delta(idx4[k] - centers[k], lim) - 0.5;
                    }
                    quat_t qc0 = {c0[0], c0[1], c0[2], c0[3]};
                    quat_t g0  = g_bpst_at(c0, sign);

                    for (int mu = 0; mu < 4; ++mu) {
                        double c1[4] = {c0[0], c0[1], c0[2], c0[3]};
                        c1[mu] += 1.0;
                        double r2_mid = 0.0;
                        for (int k = 0; k < 4; ++k) {
                            double mk = 0.5 * (c0[k] + c1[k]);
                            r2_mid += mk * mk;
                        }
                        double rm = sqrt(r2_mid);
                        if (rm >= r_out) {
                            continue;
                        }

                        double perp2 = rho * rho;
                        for (int k = 0; k < 4; ++k) {
                            if (k != mu) {
                                perp2 += c0[k] * c0[k];
                            }
                        }
                        double s   = sqrt(perp2);
                        double ang = (atan(c1[mu] / s) - atan(c0[mu] / s)) / s;
                        quat_t q   = (sign > 0) ? q_mul(qc0, q_conj(basis[mu]))
                                                : q_mul(q_conj(qc0), basis[mu]);
                        quat_t u_reg  = q_exp_im(ang * q.x, ang * q.y, ang * q.z);
                        quat_t g1     = g_bpst_at(c1, sign);
                        quat_t u_sing = q_mul(q_mul(g0, u_reg), q_conj(g1));

                        double tap = 1.0;
                        if (rm > r_in) {
                            double u = 1.0 - (rm - r_in) / (r_out - r_in);
                            tap      = u * u * (3.0 - 2.0 * u);
                        }
                        double w_clamped = (u_sing.w > 1.0) ? 1.0 : ((u_sing.w < -1.0) ? -1.0 : u_sing.w);
                        double th        = acos(w_clamped) * tap;
                        double vn        = sqrt(u_sing.x * u_sing.x + u_sing.y * u_sing.y + u_sing.z * u_sing.z);
                        quat_t u_inst    = q_id();
                        if (vn > 1e-12) {
                            double sc  = sin(th) / vn;
                            u_inst.w   = cos(th);
                            u_inst.x   = u_sing.x * sc;
                            u_inst.y   = u_sing.y * sc;
                            u_inst.z   = u_sing.z * sc;
                        }
                        if (inv_dir < 0) {
                            u_inst = q_conj(u_inst);
                        }
                        size_t id  = link_idx(dst, t, x, y, z, mu);
                        quat_t res = q_mul(u_inst, src->u[id]);
                        q_normalize(&res);
                        dst->u[id] = res;
                    }
                }
            }
        }
    }
}

static int propose_bare_bpst_instanton(su2_lat_t *lat, su2_lat_t *prop, double beta, rng_t *rng) {
    int t0      = rng_below(rng, lat->nt);
    int x0      = rng_below(rng, lat->ns);
    int y0      = rng_below(rng, lat->ns);
    int z0      = rng_below(rng, lat->ns);
    double rho  = 1.3 + 0.9 * rng_u01(rng);
    uint64_t rw = rng_u64(rng);
    int sign    = (rw & 1ULL) ? 1 : -1;
    int inv_dir = (rw & 2ULL) ? 1 : -1;

    double s_old = cfg_total_action(lat, beta, 1.0, 0);
    apply_bpst_instanton(prop, lat, t0, x0, y0, z0, rho, sign, inv_dir);
    double s_new = cfg_total_action(prop, beta, 1.0, 0);
    double dS    = s_new - s_old;
    if (dS <= 0.0 || rng_pos_u01(rng) < exp(-dS)) {
        lat_copy(lat, prop);
        return 1;
    }
    return 0;
}

/* --- Startup Mathematical & Physics Self-Tests --- */
static void run_startup_selftests(void) {
    su2_lat_t g, gp, rot, vac, rt, wrk;
    lat_alloc(&g, L_CONF, L_CONF);
    lat_alloc(&gp, L_CONF, L_CONF);
    lat_alloc(&rot, L_CONF, L_CONF);
    lat_alloc(&vac, L_CONF, L_CONF);
    lat_alloc(&rt, L_CONF, L_CONF);
    lat_alloc(&wrk, L_CONF, L_CONF);

    rng_t rng;
    rng_seed(&rng, 0x5024D2026ULL);

    /* Create generic hot SU(2) lattice on 10^4 */
    for (int i = 0; i < 4 * g.vol; ++i) {
        g.u[i] = sample_su2_boltzmann(0.0, &rng);
    }

    /* 1. Random 4D local SU(2) gauge transformation invariance */
    quat_t *gauge_map = (quat_t *)malloc((size_t)g.vol * sizeof(quat_t));
    if (gauge_map == NULL) {
        fprintf(stderr, "ERROR: OOM in selftest gauge_map\n");
        exit(1);
    }
    for (int i = 0; i < g.vol; ++i) {
        gauge_map[i] = sample_su2_boltzmann(0.0, &rng);
    }
    for (int t = 0; t < g.nt; ++t) {
        for (int x = 0; x < g.ns; ++x) {
            for (int y = 0; y < g.ns; ++y) {
                for (int z = 0; z < g.ns; ++z) {
                    const int c[4] = {t, x, y, z};
                    quat_t gs      = gauge_map[site_idx(&g, t, x, y, z)];
                    for (int mu = 0; mu < 4; ++mu) {
                        int cp[4];
                        step_coord(&g, c, mu, +1, cp);
                        quat_t gsp = gauge_map[site_idx(&g, cp[0], cp[1], cp[2], cp[3])];
                        gp.u[link_idx(&gp, t, x, y, z, mu)] =
                            q_mul(q_mul(gs, g.u[link_idx(&g, t, x, y, z, mu)]), q_conj(gsp));
                    }
                }
            }
        }
    }
    free(gauge_map);

    double abs_l0 = 0.0, abs_lp = 0.0;
    double a0 = cfg_total_action(&g, 2.4, 0.43, L_DEFECT_4D);
    double ap = cfg_total_action(&gp, 2.4, 0.43, L_DEFECT_4D);
    double p0 = cfg_mean_plaq(&g), pp = cfg_mean_plaq(&gp);
    double l0 = cfg_polyakov_mean(&g, &abs_l0), lp = cfg_polyakov_mean(&gp, &abs_lp);
    double q0 = cfg_topo_charge_clover(&g), qp = cfg_topo_charge_clover(&gp);
    if (fabs(a0 - ap) > 1e-8 || fabs(p0 - pp) > 1e-11 ||
        fabs(l0 - lp) > 1e-11 || fabs(q0 - qp) > 1e-9) {
        fprintf(stderr, "SELFTEST FAIL 1: gauge invariance mismatch dS=%.3e dP=%.3e dL=%.3e dQ=%.3e\n",
                a0 - ap, p0 - pp, l0 - lp, q0 - qp);
        exit(1);
    }

    /* 2. Z_2 center flip exact symmetry on t = 0 */
    lat_copy(&gp, &g);
    cfg_flip_z2_slice(&gp, 0);
    double a_z2 = cfg_total_action(&gp, 2.4, 1.0, 0);
    double a_g1 = cfg_total_action(&g, 2.4, 1.0, 0);
    double l_z2 = cfg_polyakov_mean(&gp, NULL);
    if (fabs(a_z2 - a_g1) > 1e-9 || fabs(l_z2 + l0) > 1e-12) {
        fprintf(stderr, "SELFTEST FAIL 2: Z_2 flip mismatch dS=%.3e L_sum=%.3e\n", a_z2 - a_g1, l_z2 + l0);
        exit(1);
    }

    /* 3. Microcanonical checkerboard Over-Relaxation exactness (dS = 0, |q|^2 = 1) */
    lat_copy(&gp, &g);
    for (int par = 0; par < 2; ++par) {
        for (int t = 0; t < gp.nt; ++t) {
            for (int x = 0; x < gp.ns; ++x) {
                for (int y = 0; y < gp.ns; ++y) {
                    for (int z = 0; z < gp.ns; ++z) {
                        if (((t + x + y + z) & 1) != par) {
                            continue;
                        }
                        for (int mu = 0; mu < 4; ++mu) {
                            quat_t v  = compute_staple(&gp, t, x, y, z, mu, 0.43, L_DEFECT_4D);
                            double v2 = q_norm2(v);
                            if (v2 > 1e-24) {
                                size_t id  = link_idx(&gp, t, x, y, z, mu);
                                quat_t u   = gp.u[id];
                                double sc  = (2.0 * q_inner(u, v)) / v2;
                                quat_t u1  = {sc * v.w - u.w, sc * v.x - u.x,
                                              sc * v.y - u.y, sc * v.z - u.z};
                                double n2  = q_norm2(u1);
                                if (fabs(n2 - 1.0) > 1e-12) {
                                    fprintf(stderr, "SELFTEST FAIL 3: OR norm drift %.3e\n", n2 - 1.0);
                                    exit(1);
                                }
                                gp.u[id] = u1;
                            }
                        }
                    }
                }
            }
        }
    }
    double a_or = cfg_total_action(&gp, 2.4, 0.43, L_DEFECT_4D);
    if (fabs(a_or - a0) > 1e-8) {
        fprintf(stderr, "SELFTEST FAIL 3: OR dS=%.3e != 0\n", a_or - a0);
        exit(1);
    }

    /* 4. Fast 4D defect-wall @swap dS exactness (36 defect plaquettes vs full 60,000 plaquettes) */
    for (int i = 0; i < 4 * rot.vol; ++i) {
        rot.u[i] = sample_su2_boltzmann(0.0, &rng);
    }
    double c1 = 0.72, c2 = 0.31;
    double sa_1 = cfg_total_action(&g, 2.4, c1, L_DEFECT_4D);
    double sb_2 = cfg_total_action(&rot, 2.4, c2, L_DEFECT_4D);
    double sa_2 = cfg_total_action(&g, 2.4, c2, L_DEFECT_4D);
    double sb_1 = cfg_total_action(&rot, 2.4, c1, L_DEFECT_4D);
    int def_cntA = 0, def_cntB = 0;
    double sd_A = cfg_defect_sum_re(&g, L_DEFECT_4D, &def_cntA);
    double sd_B = cfg_defect_sum_re(&rot, L_DEFECT_4D, &def_cntB);
    double ds_exact = (sa_2 + sb_1) - (sa_1 + sb_2);
    double ds_fast  = 2.4 * (c1 - c2) * (sd_A - sd_B);
    if (def_cntA != 36 || def_cntB != 36 || fabs(ds_exact - ds_fast) > 1e-9) {
        fprintf(stderr, "SELFTEST FAIL 4: defect swap cnt=(%d,%d) dS %.3e vs %.3e\n",
                def_cntA, def_cntB, ds_exact, ds_fast);
        exit(1);
    }

    /* 5. Creutz-Kennedy-Pendleton sampler exactness at alpha=4.0 and alpha=0.2 vs I_2/I_1 */
    const double test_alphas[2] = {4.0, 0.2};
    for (int ia = 0; ia < 2; ++ia) {
        double al    = test_alphas[ia];
        double s1 = 0.0, s2 = 0.0;
        const int n_kp = 100000;
        for (int i = 0; i < n_kp; ++i) {
            double a = sample_su2_a0(al, &rng);
            s1 += a;
            s2 += a * a;
        }
        double mean_a = s1 / (double)n_kp;
        double err_a  = sqrt((s2 / (double)n_kp - mean_a * mean_a) / (double)n_kp);
        /* Analytical I_1(al) and I_2(al) power series */
        double i1 = 0.5 * al, i2 = 0.125 * al * al;
        double t1 = i1, t2 = i2;
        for (int k = 1; k < 40; ++k) {
            double half_al2 = 0.25 * al * al;
            t1 *= half_al2 / (double)(k * (k + 1));
            t2 *= half_al2 / (double)(k * (k + 2));
            i1 += t1;
            i2 += t2;
        }
        double exact_ratio = i2 / i1;
        if (fabs(mean_a - exact_ratio) > 4.0 * err_a) {
            fprintf(stderr, "SELFTEST FAIL 5: KP alpha=%.1f mean=%.6f exact=%.6f\n",
                    al, mean_a, exact_ratio);
            exit(1);
        }
    }

    /* 6. 4D singular-gauge BPST instanton on cold 10^4 vacuum: Q_cool=+-1, S/(2pi^2*beta) in [0.98, 1.05], rt <1e-12 */
    lat_init_cold(&vac);
    for (int is = 0; is < 2; ++is) {
        int sign = (is == 0) ? 1 : -1;
        apply_bpst_instanton(&gp, &vac, 5, 5, 5, 5, 1.6, sign, 1);
        apply_bpst_instanton(&rt, &gp, 5, 5, 5, 5, 1.6, sign, -1);
        double max_rt_err = 0.0;
        for (int i = 0; i < 4 * gp.vol; ++i) {
            double n2 = q_norm2(gp.u[i]);
            if (fabs(n2 - 1.0) > 1e-12) {
                fprintf(stderr, "SELFTEST FAIL 6: BPST link norm %.3e\n", n2);
                exit(1);
            }
            double dw = fabs(rt.u[i].w - vac.u[i].w);
            double dx = fabs(rt.u[i].x - vac.u[i].x);
            double dy = fabs(rt.u[i].y - vac.u[i].y);
            double dz = fabs(rt.u[i].z - vac.u[i].z);
            if (dw > max_rt_err) max_rt_err = dw;
            if (dx > max_rt_err) max_rt_err = dx;
            if (dy > max_rt_err) max_rt_err = dy;
            if (dz > max_rt_err) max_rt_err = dz;
        }
        double s_ratio   = cfg_total_action(&gp, 2.4, 1.0, 0) / (2.0 * K_PI * K_PI * 2.4);
        double q_cool_r  = 0.0;
        long long q_cool = cfg_topo_charge_cooled(&gp, &wrk, &q_cool_r);
        if (max_rt_err > 1e-12 || s_ratio < 0.98 || s_ratio > 1.05 || q_cool != sign) {
            fprintf(stderr, "SELFTEST FAIL 6: BPST sign=%+d S/2pi2b=%.4f Qc=%lld (%.3f) rt=%.3e\n",
                    sign, s_ratio, q_cool, q_cool_r, max_rt_err);
            exit(1);
        }
        if (sign == 1) {
            printf("  [SELFTEST] PASS (4D SU(2) gauge <1e-8, Z2 exact <1e-12, OR <1e-9, "
                   "36-plaq @swap <1e-9, CKP exact I2/I1, BPST Q=%+lld [%.3f] S/2pi2b=%.4f rt=%.2e)\n",
                   q_cool, q_cool_r, s_ratio, max_rt_err);
        }
    }

    lat_free(&g);
    lat_free(&gp);
    lat_free(&rot);
    lat_free(&vac);
    lat_free(&rt);
    lat_free(&wrk);
}

/* --- Statistics Helpers: Madras-Sokal + 20-Block Creutz Ratio Jackknife --- */
static void compute_madras_sokal(const double *x, int n, double *out_mean,
                                 double *out_var, double *out_tau_int,
                                 double *out_err, int *out_frozen) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += x[i];
    }
    double mean = sum / (double)n;
    double c0   = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = x[i] - mean;
        c0 += d * d;
    }
    c0 /= (double)n;
    *out_mean = mean;
    *out_var  = c0;

    if (c0 <= 1e-15) {
        *out_frozen  = 1;
        *out_tau_int = -1.0;
        *out_err     = 0.0;
        return;
    }

    *out_frozen = 0;
    double tau  = 0.5;
    int max_w   = n / 2;
    for (int w = 1; w <= max_w; ++w) {
        double ct = 0.0;
        int lim   = n - w;
        for (int i = 0; i < lim; ++i) {
            ct += (x[i] - mean) * (x[i + w] - mean);
        }
        ct /= (double)lim;
        tau += ct / c0;
        double tau_floor = (tau > 0.5) ? tau : 0.5;
        if ((double)w >= 6.0 * tau_floor) {
            break;
        }
    }
    if (tau < 0.5) {
        tau = 0.5;
    }
    *out_tau_int = tau;
    *out_err     = sqrt((2.0 * tau * c0) / (double)n);
}

static void compute_creutz_jackknife(const double *waa, const double *wab, const double *wbb,
                                     int n, double *out_chi, double *out_err) {
    double sa = 0.0, sab = 0.0, sb = 0.0;
    for (int i = 0; i < n; ++i) {
        sa += waa[i];
        sab += wab[i];
        sb += wbb[i];
    }
    double ratio    = (sb * sa) / (sab * sab);
    double chi_full = (ratio > 1e-12) ? (-log(ratio)) : 0.0;

    int bsz = n / JK_BLOCKS;
    double jk[JK_BLOCKS];
    double jk_mean = 0.0;
    for (int b = 0; b < JK_BLOCKS; ++b) {
        double bsa = 0.0, bsab = 0.0, bsb = 0.0;
        int lo = b * bsz;
        int hi = (b == JK_BLOCKS - 1) ? n : ((b + 1) * bsz);
        for (int i = lo; i < hi; ++i) {
            bsa += waa[i];
            bsab += wab[i];
            bsb += wbb[i];
        }
        double r_jk = ((sb - bsb) * (sa - bsa)) / ((sab - bsab) * (sab - bsab));
        jk[b]       = (r_jk > 1e-12) ? (-log(r_jk)) : 0.0;
        jk_mean += jk[b];
    }
    jk_mean /= (double)JK_BLOCKS;
    double jk_var = 0.0;
    for (int b = 0; b < JK_BLOCKS; ++b) {
        double d = jk[b] - jk_mean;
        jk_var += d * d;
    }
    jk_var *= (double)(JK_BLOCKS - 1) / (double)JK_BLOCKS;
    *out_chi = chi_full;
    *out_err = sqrt(jk_var);
}

/* Fisher 2-Jet geodesic ladder redistribution from empirical swap probabilities */
static void adapt_fisher_ladder(double c[R_REPLICAS], const double avg_p[R_REPLICAS - 1]) {
    double ds[R_REPLICAS - 1];
    double total_s = 0.0;
    for (int r = 0; r < R_REPLICAS - 1; ++r) {
        double p = avg_p[r];
        if (p < 0.01) {
            p = 0.01;
        }
        if (p > 0.99) {
            p = 0.99;
        }
        ds[r] = sqrt(-log(p));
        total_s += ds[r];
    }
    if (total_s <= 1e-12) {
        return;
    }
    double cum_c[R_REPLICAS];
    double cum_s[R_REPLICAS];
    cum_s[0] = 0.0;
    for (int r = 0; r < R_REPLICAS; ++r) {
        cum_c[r] = c[r];
        if (r > 0) {
            cum_s[r] = cum_s[r - 1] + ds[r - 1];
        }
    }

    double c_new[R_REPLICAS];
    c_new[0]              = 1.0;
    c_new[R_REPLICAS - 1] = 0.0;
    for (int r = 1; r < R_REPLICAS - 1; ++r) {
        double s_tgt = ((double)r / (double)(R_REPLICAS - 1)) * total_s;
        int k        = 0;
        while (k + 1 < R_REPLICAS - 1 && cum_s[k + 1] < s_tgt) {
            ++k;
        }
        double denom  = cum_s[k + 1] - cum_s[k];
        double t      = (denom > 1e-12) ? ((s_tgt - cum_s[k]) / denom) : 0.5;
        double interp = cum_c[k] + t * (cum_c[k + 1] - cum_c[k]);
        c_new[r]      = 0.5 * c[r] + 0.5 * interp;
    }
    for (int r = 1; r < R_REPLICAS - 1; ++r) {
        if (c_new[r] > c_new[r - 1] - 0.01) {
            c_new[r] = c_new[r - 1] - 0.01;
        }
        double min_floor = (double)(R_REPLICAS - 1 - r) * 0.01;
        if (c_new[r] < min_floor) {
            c_new[r] = min_floor;
        }
    }
    for (int r = 0; r < R_REPLICAS; ++r) {
        c[r] = c_new[r];
    }
}

/* --- Experiment 1: 4D Thermal Z_2 Deconfinement (4 x 12^3) --- */
typedef struct {
    const char *solver;
    double beta;
    double mean_p;
    double err_p;
    double mean_abs_l;
    double err_abs_l;
    double mean_l;
    double err_l;
    double tau_l;
    int frozen_l;
    long long z2_sign_flips;
} z2_res_t;

typedef struct {
    int idx;
    double beta;
    z2_res_t res[2];
} z2_job_t;

static z2_res_t run_z2_solver(int use_center_flip, const char *name, double beta, uint64_t seed) {
    su2_lat_t lat;
    lat_alloc(&lat, NT_THERM, NS_THERM);
    lat_init_cold(&lat);
    rng_t rng;
    rng_seed(&rng, seed);

    double *p_h    = (double *)malloc((size_t)N_MEAS_Z2 * sizeof(double));
    double *l_h    = (double *)malloc((size_t)N_MEAS_Z2 * sizeof(double));
    double *absl_h = (double *)malloc((size_t)N_MEAS_Z2 * sizeof(double));
    if (!p_h || !l_h || !absl_h) {
        fprintf(stderr, "ERROR: OOM in run_z2_solver\n");
        exit(1);
    }

    for (int t = 0; t < N_THERM_Z2; ++t) {
        sweep_su2_hb_or(&lat, beta, 1.0, 0, &rng);
        if (use_center_flip && (rng_u64(&rng) & 1ULL)) {
            cfg_flip_z2_slice(&lat, 0);
        }
    }

    long long flips = 0;
    double prev_l   = 0.0;
    for (int t = 0; t < N_MEAS_Z2; ++t) {
        sweep_su2_hb_or(&lat, beta, 1.0, 0, &rng);
        if (use_center_flip && (rng_u64(&rng) & 1ULL)) {
            cfg_flip_z2_slice(&lat, 0);
        }
        p_h[t] = cfg_mean_plaq(&lat);
        double al = 0.0;
        double lv = cfg_polyakov_mean(&lat, &al);
        l_h[t]    = lv;
        absl_h[t] = al;
        if (t > 0 && ((prev_l >= 0.0 && lv < 0.0) || (prev_l < 0.0 && lv >= 0.0))) {
            ++flips;
        }
        prev_l = lv;
    }

    double pm, pv, pt, pe, alm, alv, alt, ale, lm, lv, lt, le;
    int pf = 0, alf = 0, lf = 0;
    compute_madras_sokal(p_h, N_MEAS_Z2, &pm, &pv, &pt, &pe, &pf);
    compute_madras_sokal(absl_h, N_MEAS_Z2, &alm, &alv, &alt, &ale, &alf);
    compute_madras_sokal(l_h, N_MEAS_Z2, &lm, &lv, &lt, &le, &lf);

    free(p_h);
    free(l_h);
    free(absl_h);
    lat_free(&lat);

    z2_res_t res;
    res.solver        = name;
    res.beta          = beta;
    res.mean_p        = pm;
    res.err_p         = pe;
    res.mean_abs_l    = alm;
    res.err_abs_l     = ale;
    res.mean_l        = lm;
    res.err_l         = le;
    res.tau_l         = lt;
    res.frozen_l      = lf;
    res.z2_sign_flips = flips;
    return res;
}

static void *z2_worker(void *arg) {
    z2_job_t *job = (z2_job_t *)arg;
    uint64_t base = 0x50211100ULL + (uint64_t)job->idx * 100ULL;
    job->res[0]   = run_z2_solver(0, "local_hb_or", job->beta, base + 1ULL);
    job->res[1]   = run_z2_solver(1, "symbiotic_z2_center", job->beta, base + 2ULL);
    return NULL;
}

/* --- Experiment 2: 4D Confined Hypercube Confinement, BPST Corollary 2B & Defect-Wall PTBC (10^4) --- */
typedef enum {
    SU2_LOCAL_HB_OR    = 0,
    SU2_BARE_BPST_INST = 1,
    SU2_SYMBIOTIC_PTBC = 2
} su2_solver_t;

typedef struct {
    const char *solver;
    double beta;
    int n_replicas;
    int work_sweeps;
    double mean_p;
    double err_p;
    double mean_abs_l;
    double err_abs_l;
    double w11;
    double w12;
    double w22;
    double w23;
    double w33;
    double chi22;
    double err_chi22;
    double chi33;
    double err_chi33;
    double chi_t_cool;
    double err_chi_t_cool;
    double chi_quarter_over_sigma;
    double tau_work_q_cool;
    int frozen_q_cool;
    long long trans_q_cool;
    double bare_inst_acc;
    double min_swap_acc;
    double mean_swap_acc;
    double c_ladder[R_REPLICAS];
} su2_conf_res_t;

typedef struct {
    int idx;
    double beta;
    su2_conf_res_t res[3];
} su2_conf_job_t;

static su2_conf_res_t run_conf_single(su2_solver_t kind, const char *name, double beta, uint64_t seed) {
    su2_lat_t lat, prop, wrk;
    lat_alloc(&lat, L_CONF, L_CONF);
    lat_alloc(&prop, L_CONF, L_CONF);
    lat_alloc(&wrk, L_CONF, L_CONF);
    lat_init_cold(&lat);
    rng_t rng;
    rng_seed(&rng, seed);

    const int n_samp = N_EQ_MEAS / R_REPLICAS; /* 300 exact, spaced by R_REPLICAS = 8 sweeps */
    double *p_h    = (double *)malloc((size_t)n_samp * sizeof(double));
    double *absl_h = (double *)malloc((size_t)n_samp * sizeof(double));
    double *w11_h  = (double *)malloc((size_t)n_samp * sizeof(double));
    double *w12_h  = (double *)malloc((size_t)n_samp * sizeof(double));
    double *w22_h  = (double *)malloc((size_t)n_samp * sizeof(double));
    double *w23_h  = (double *)malloc((size_t)n_samp * sizeof(double));
    double *w33_h  = (double *)malloc((size_t)n_samp * sizeof(double));
    double *qc_h   = (double *)malloc((size_t)n_samp * sizeof(double));
    double *qc2_h  = (double *)malloc((size_t)n_samp * sizeof(double));
    if (!p_h || !absl_h || !w11_h || !w12_h || !w22_h || !w23_h || !w33_h || !qc_h || !qc2_h) {
        fprintf(stderr, "ERROR: OOM in run_conf_single\n");
        exit(1);
    }

    for (int t = 0; t < N_THERM_1R; ++t) {
        sweep_su2_hb_or(&lat, beta, 1.0, 0, &rng);
        if (kind == SU2_BARE_BPST_INST) {
            propose_bare_bpst_instanton(&lat, &prop, beta, &rng);
        }
    }

    long long inst_acc = 0;
    long long trans_qc = 0;
    long long prev_qc  = 0;
    int si             = 0;
    const double vol4  = (double)(L_CONF * L_CONF * L_CONF * L_CONF);
    for (int t = 0; t < N_EQ_MEAS; ++t) {
        sweep_su2_hb_or(&lat, beta, 1.0, 0, &rng);
        if (kind == SU2_BARE_BPST_INST) {
            inst_acc += propose_bare_bpst_instanton(&lat, &prop, beta, &rng);
        }
        if ((t % R_REPLICAS) == (R_REPLICAS - 1)) {
            double p_val  = cfg_mean_plaq(&lat);
            double al_val = 0.0;
            cfg_polyakov_mean(&lat, &al_val);
            p_h[si]    = p_val;
            absl_h[si] = al_val;
            w11_h[si]  = p_val;
            w12_h[si]  = cfg_wilson_loop(&lat, 1, 2);
            w22_h[si]  = cfg_wilson_loop(&lat, 2, 2);
            w23_h[si]  = cfg_wilson_loop(&lat, 2, 3);
            w33_h[si]  = cfg_wilson_loop(&lat, 3, 3);
            long long qc = cfg_topo_charge_cooled(&lat, &wrk, NULL);
            qc_h[si]     = (double)qc;
            qc2_h[si]    = (double)(qc * qc) / vol4;
            if (si > 0 && qc != prev_qc) {
                ++trans_qc;
            }
            prev_qc = qc;
            ++si;
        }
    }

    double pm, pv, pt, pe, alm, alv, alt, ale;
    double w11m, w12m, w22m, w23m, w33m, dummy;
    double qcm, qcv, qctau, qcerr, qc2m, qc2v, qc2tau, qc2err;
    double chi22, echi22, chi33, echi33;
    int pf = 0, alf = 0, df = 0, qcf = 0, qc2f = 0;

    compute_madras_sokal(p_h, n_samp, &pm, &pv, &pt, &pe, &pf);
    compute_madras_sokal(absl_h, n_samp, &alm, &alv, &alt, &ale, &alf);
    compute_madras_sokal(w11_h, n_samp, &w11m, &dummy, &dummy, &dummy, &df);
    compute_madras_sokal(w12_h, n_samp, &w12m, &dummy, &dummy, &dummy, &df);
    compute_madras_sokal(w22_h, n_samp, &w22m, &dummy, &dummy, &dummy, &df);
    compute_madras_sokal(w23_h, n_samp, &w23m, &dummy, &dummy, &dummy, &df);
    compute_madras_sokal(w33_h, n_samp, &w33m, &dummy, &dummy, &dummy, &df);
    compute_madras_sokal(qc_h, n_samp, &qcm, &qcv, &qctau, &qcerr, &qcf);
    compute_madras_sokal(qc2_h, n_samp, &qc2m, &qc2v, &qc2tau, &qc2err, &qc2f);
    compute_creutz_jackknife(w11_h, w12_h, w22_h, n_samp, &chi22, &echi22);
    compute_creutz_jackknife(w22_h, w23_h, w33_h, n_samp, &chi33, &echi33);

    free(p_h);
    free(absl_h);
    free(w11_h);
    free(w12_h);
    free(w22_h);
    free(w23_h);
    free(w33_h);
    free(qc_h);
    free(qc2_h);
    lat_free(&lat);
    lat_free(&prop);
    lat_free(&wrk);

    su2_conf_res_t res;
    memset(&res, 0, sizeof(res));
    res.solver                 = name;
    res.beta                   = beta;
    res.n_replicas             = 1;
    res.work_sweeps            = N_EQ_MEAS;
    res.mean_p                 = pm;
    res.err_p                  = pe;
    res.mean_abs_l             = alm;
    res.err_abs_l              = ale;
    res.w11                    = w11m;
    res.w12                    = w12m;
    res.w22                    = w22m;
    res.w23                    = w23m;
    res.w33                    = w33m;
    res.chi22                  = chi22;
    res.err_chi22              = echi22;
    res.chi33                  = chi33;
    res.err_chi33              = echi33;
    res.chi_t_cool             = qc2m;
    res.err_chi_t_cool         = qc2err;
    res.chi_quarter_over_sigma = (chi22 > 1e-12 && qc2m > 0.0) ? (pow(qc2m, 0.25) / sqrt(chi22)) : 0.0;
    res.tau_work_q_cool        = qcf ? -1.0 : ((double)R_REPLICAS * qctau);
    res.frozen_q_cool          = qcf;
    res.trans_q_cool           = trans_qc;
    res.bare_inst_acc          = (kind == SU2_BARE_BPST_INST) ? ((double)inst_acc / (double)N_EQ_MEAS) : 0.0;
    res.c_ladder[0]            = 1.0;
    return res;
}

static su2_conf_res_t run_conf_ptbc(double beta, uint64_t seed) {
    su2_lat_t rep_buf[R_REPLICAS];
    su2_lat_t *rep[R_REPLICAS];
    for (int r = 0; r < R_REPLICAS; ++r) {
        lat_alloc(&rep_buf[r], L_CONF, L_CONF);
        lat_init_cold(&rep_buf[r]);
        rep[r] = &rep_buf[r];
    }
    su2_lat_t tmp, wrk;
    lat_alloc(&tmp, L_CONF, L_CONF);
    lat_alloc(&wrk, L_CONF, L_CONF);

    rng_t rng;
    rng_seed(&rng, seed);

    double c[R_REPLICAS];
    for (int r = 0; r < R_REPLICAS; ++r) {
        c[r] = 1.0 - (double)r / (double)(R_REPLICAS - 1);
    }

    double *p_h    = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *absl_h = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *w11_h  = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *w12_h  = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *w22_h  = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *w23_h  = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *w33_h  = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *qc_h   = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *qc2_h  = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    if (!p_h || !absl_h || !w11_h || !w12_h || !w22_h || !w23_h || !w33_h || !qc_h || !qc2_h) {
        fprintf(stderr, "ERROR: OOM in run_conf_ptbc\n");
        exit(1);
    }

    double swap_prob_sum[R_REPLICAS - 1];
    long long swap_prob_cnt[R_REPLICAS - 1];
    memset(swap_prob_sum, 0, sizeof(swap_prob_sum));
    memset(swap_prob_cnt, 0, sizeof(swap_prob_cnt));

    /* Thermalization: adapt Fisher 2-Jet ladder during first half, then FREEZE cleanly! */
    for (int t = 0; t < N_THERM_PT; ++t) {
        for (int r = 0; r < R_REPLICAS; ++r) {
            sweep_su2_hb_or(rep[r], beta, c[r], L_DEFECT_4D, &rng);
        }
        for (int pass = 0; pass < 2; ++pass) {
            int parity = (t + pass) & 1;
            for (int r = parity; r + 1 < R_REPLICAS; r += 2) {
                double sd_a = cfg_defect_sum_re(rep[r], L_DEFECT_4D, NULL);
                double sd_b = cfg_defect_sum_re(rep[r + 1], L_DEFECT_4D, NULL);
                double dS   = beta * (c[r] - c[r + 1]) * (sd_a - sd_b);
                double p    = (dS <= 0.0) ? 1.0 : exp(-dS);
                if (t < N_THERM_PT / 2) {
                    swap_prob_sum[r] += p;
                    ++swap_prob_cnt[r];
                }
                if (dS <= 0.0 || rng_pos_u01(&rng) < p) {
                    su2_lat_t *sw = rep[r];
                    rep[r]        = rep[r + 1];
                    rep[r + 1]    = sw;
                }
            }
        }
        shift_torus_4d_random(rep[0], &tmp, &rng);

        if (t < N_THERM_PT / 2 && (t % 25) == 24) {
            double avg_p[R_REPLICAS - 1];
            for (int r = 0; r < R_REPLICAS - 1; ++r) {
                avg_p[r]         = (swap_prob_cnt[r] > 0) ? (swap_prob_sum[r] / (double)swap_prob_cnt[r]) : 0.5;
                swap_prob_sum[r] = 0.0;
                swap_prob_cnt[r] = 0;
            }
            adapt_fisher_ladder(c, avg_p);
        }
    }

    /* Frozen measurement phase (exact product-weight detailed balance!) */
    long long swap_acc[R_REPLICAS - 1];
    long long swap_try[R_REPLICAS - 1];
    memset(swap_acc, 0, sizeof(swap_acc));
    memset(swap_try, 0, sizeof(swap_try));
    long long trans_qc = 0, prev_qc = 0;
    const double vol4 = (double)(L_CONF * L_CONF * L_CONF * L_CONF);

    for (int t = 0; t < N_MEAS_PT; ++t) {
        for (int r = 0; r < R_REPLICAS; ++r) {
            sweep_su2_hb_or(rep[r], beta, c[r], L_DEFECT_4D, &rng);
        }
        for (int pass = 0; pass < 2; ++pass) {
            int parity = (t + pass) & 1;
            for (int r = parity; r + 1 < R_REPLICAS; r += 2) {
                double sd_a = cfg_defect_sum_re(rep[r], L_DEFECT_4D, NULL);
                double sd_b = cfg_defect_sum_re(rep[r + 1], L_DEFECT_4D, NULL);
                double dS   = beta * (c[r] - c[r + 1]) * (sd_a - sd_b);
                ++swap_try[r];
                if (dS <= 0.0 || rng_pos_u01(&rng) < exp(-dS)) {
                    su2_lat_t *sw = rep[r];
                    rep[r]        = rep[r + 1];
                    rep[r + 1]    = sw;
                    ++swap_acc[r];
                }
            }
        }
        shift_torus_4d_random(rep[0], &tmp, &rng);

        double p_val  = cfg_mean_plaq(rep[0]);
        double al_val = 0.0;
        cfg_polyakov_mean(rep[0], &al_val);
        p_h[t]    = p_val;
        absl_h[t] = al_val;
        w11_h[t]  = p_val;
        w12_h[t]  = cfg_wilson_loop(rep[0], 1, 2);
        w22_h[t]  = cfg_wilson_loop(rep[0], 2, 2);
        w23_h[t]  = cfg_wilson_loop(rep[0], 2, 3);
        w33_h[t]  = cfg_wilson_loop(rep[0], 3, 3);
        long long qc = cfg_topo_charge_cooled(rep[0], &wrk, NULL);
        qc_h[t]      = (double)qc;
        qc2_h[t]     = (double)(qc * qc) / vol4;
        if (t > 0 && qc != prev_qc) {
            ++trans_qc;
        }
        prev_qc = qc;
    }

    double min_sw = 1.0, sum_sw = 0.0;
    for (int r = 0; r < R_REPLICAS - 1; ++r) {
        double rate = (swap_try[r] > 0) ? ((double)swap_acc[r] / (double)swap_try[r]) : 0.0;
        if (rate < min_sw) {
            min_sw = rate;
        }
        sum_sw += rate;
    }

    double pm, pv, pt, pe, alm, alv, alt, ale;
    double w11m, w12m, w22m, w23m, w33m, dummy;
    double qcm, qcv, qctau, qcerr, qc2m, qc2v, qc2tau, qc2err;
    double chi22, echi22, chi33, echi33;
    int pf = 0, alf = 0, df = 0, qcf = 0, qc2f = 0;

    compute_madras_sokal(p_h, N_MEAS_PT, &pm, &pv, &pt, &pe, &pf);
    compute_madras_sokal(absl_h, N_MEAS_PT, &alm, &alv, &alt, &ale, &alf);
    compute_madras_sokal(w11_h, N_MEAS_PT, &w11m, &dummy, &dummy, &dummy, &df);
    compute_madras_sokal(w12_h, N_MEAS_PT, &w12m, &dummy, &dummy, &dummy, &df);
    compute_madras_sokal(w22_h, N_MEAS_PT, &w22m, &dummy, &dummy, &dummy, &df);
    compute_madras_sokal(w23_h, N_MEAS_PT, &w23m, &dummy, &dummy, &dummy, &df);
    compute_madras_sokal(w33_h, N_MEAS_PT, &w33m, &dummy, &dummy, &dummy, &df);
    compute_madras_sokal(qc_h, N_MEAS_PT, &qcm, &qcv, &qctau, &qcerr, &qcf);
    compute_madras_sokal(qc2_h, N_MEAS_PT, &qc2m, &qc2v, &qc2tau, &qc2err, &qc2f);
    compute_creutz_jackknife(w11_h, w12_h, w22_h, N_MEAS_PT, &chi22, &echi22);
    compute_creutz_jackknife(w22_h, w23_h, w33_h, N_MEAS_PT, &chi33, &echi33);

    free(p_h);
    free(absl_h);
    free(w11_h);
    free(w12_h);
    free(w22_h);
    free(w23_h);
    free(w33_h);
    free(qc_h);
    free(qc2_h);
    for (int r = 0; r < R_REPLICAS; ++r) {
        lat_free(&rep_buf[r]);
    }
    lat_free(&tmp);
    lat_free(&wrk);

    su2_conf_res_t res;
    memset(&res, 0, sizeof(res));
    res.solver                 = "symbiotic_ptbc_4d";
    res.beta                   = beta;
    res.n_replicas             = R_REPLICAS;
    res.work_sweeps            = R_REPLICAS * N_MEAS_PT;
    res.mean_p                 = pm;
    res.err_p                  = pe;
    res.mean_abs_l             = alm;
    res.err_abs_l              = ale;
    res.w11                    = w11m;
    res.w12                    = w12m;
    res.w22                    = w22m;
    res.w23                    = w23m;
    res.w33                    = w33m;
    res.chi22                  = chi22;
    res.err_chi22              = echi22;
    res.chi33                  = chi33;
    res.err_chi33              = echi33;
    res.chi_t_cool             = qc2m;
    res.err_chi_t_cool         = qc2err;
    res.chi_quarter_over_sigma = (chi22 > 1e-12 && qc2m > 0.0) ? (pow(qc2m, 0.25) / sqrt(chi22)) : 0.0;
    res.tau_work_q_cool        = qcf ? -1.0 : ((double)R_REPLICAS * qctau);
    res.frozen_q_cool          = qcf;
    res.trans_q_cool           = trans_qc;
    res.min_swap_acc           = min_sw;
    res.mean_swap_acc          = sum_sw / (double)(R_REPLICAS - 1);
    for (int r = 0; r < R_REPLICAS; ++r) {
        res.c_ladder[r] = c[r];
    }
    return res;
}

static void *conf_worker(void *arg) {
    su2_conf_job_t *job = (su2_conf_job_t *)arg;
    uint64_t base       = 0x50222200ULL + (uint64_t)job->idx * 1000ULL;
    job->res[0]         = run_conf_single(SU2_LOCAL_HB_OR, "local_hb_or", job->beta, base + 1ULL);
    job->res[1]         = run_conf_single(SU2_BARE_BPST_INST, "bare_bpst_inst", job->beta, base + 2ULL);
    job->res[2]         = run_conf_ptbc(job->beta, base + 3ULL);
    return NULL;
}

int main(int argc, char **argv) {
    const char *json_path = NULL;
    for (int i = 1; i + 1 < argc; ++i) {
        if (strcmp(argv[i], "--json") == 0) {
            json_path = argv[i + 1];
        }
    }

    printf("=== 4D Pure SU(2) Yang-Mills Wind Tunnel (Exp1: %dx%d^3 | Exp2: %d^4 eq_work=%d) ===\n",
           NT_THERM, NS_THERM, L_CONF, N_EQ_MEAS);
    run_startup_selftests();

    /* Run Experiment 1: 4 x 12^3 Z_2 Deconfinement */
    const double betas_therm[N_BETAS_Z2] = {2.15, 2.25, 2.30, 2.35, 2.45};
    z2_job_t z2_jobs[N_BETAS_Z2];
    pthread_t z2_tids[N_BETAS_Z2];
    for (int b = 0; b < N_BETAS_Z2; ++b) {
        z2_jobs[b].idx  = b;
        z2_jobs[b].beta = betas_therm[b];
        if (pthread_create(&z2_tids[b], NULL, z2_worker, &z2_jobs[b]) != 0) {
            fprintf(stderr, "ERROR: pthread_create failed (z2)\n");
            return 1;
        }
    }
    for (int b = 0; b < N_BETAS_Z2; ++b) {
        pthread_join(z2_tids[b], NULL);
    }

    printf("\n--- [Experiment 1] Thermal 4 x 12^3 Z_2 Deconfinement & Vacuum Unfreezing ---\n");
    for (int b = 0; b < N_BETAS_Z2; ++b) {
        printf("  beta=%.2f:\n", z2_jobs[b].beta);
        for (int s = 0; s < 2; ++s) {
            const z2_res_t *r = &z2_jobs[b].res[s];
            printf("    %-20s | <P>=%.5f+-%.5f | <|L|>=%.5f+-%.5f | <L>=%+.5f+-%.5f | tau_L=%.2f | flips=%lld\n",
                   r->solver, r->mean_p, r->err_p, r->mean_abs_l, r->err_abs_l,
                   r->mean_l, r->err_l, r->tau_l, r->z2_sign_flips);
        }
    }

    /* Run Experiment 2: 10^4 Confinement, BPST Corollary 2B & Defect-Wall PTBC */
    const double betas_conf[N_BETAS_CONF] = {2.25, 2.35, 2.45, 2.52};
    su2_conf_job_t conf_jobs[N_BETAS_CONF];
    pthread_t conf_tids[N_BETAS_CONF];
    for (int b = 0; b < N_BETAS_CONF; ++b) {
        conf_jobs[b].idx  = b;
        conf_jobs[b].beta = betas_conf[b];
        if (pthread_create(&conf_tids[b], NULL, conf_worker, &conf_jobs[b]) != 0) {
            fprintf(stderr, "ERROR: pthread_create failed (conf)\n");
            return 1;
        }
    }
    for (int b = 0; b < N_BETAS_CONF; ++b) {
        pthread_join(conf_tids[b], NULL);
    }

    printf("\n--- [Experiment 2] Confined 10^4 Scaling, BPST Corollary 2B & Defect-Wall PTBC ---\n");
    for (int b = 0; b < N_BETAS_CONF; ++b) {
        printf("  beta=%.2f:\n", conf_jobs[b].beta);
        for (int s = 0; s < 3; ++s) {
            const su2_conf_res_t *r = &conf_jobs[b].res[s];
            printf("    %-20s | <P>=%.5f+-%.5f | <|L|>=%.4f | chi22=%.4f+-%.4f | chi33=%.4f+-%.4f | "
                   "chi_c=%.6f+-%.6f (Q14/s=%.3f) | tauWc=%.1f (tc=%lld) | bare=%.4f | sw_min=%.3f\n",
                   r->solver, r->mean_p, r->err_p, r->mean_abs_l, r->chi22, r->err_chi22,
                   r->chi33, r->err_chi33, r->chi_t_cool, r->err_chi_t_cool,
                   r->chi_quarter_over_sigma, r->tau_work_q_cool, r->trans_q_cool,
                   r->bare_inst_acc, r->min_swap_acc);
        }
    }

    if (json_path != NULL) {
        FILE *fp = fopen(json_path, "w");
        if (fp == NULL) {
            fprintf(stderr, "ERROR: cannot write %s\n", json_path);
            return 1;
        }
        fprintf(fp, "{\n");
        fprintf(fp, "  \"nt_therm\": %d,\n", NT_THERM);
        fprintf(fp, "  \"ns_therm\": %d,\n", NS_THERM);
        fprintf(fp, "  \"l_conf\": %d,\n", L_CONF);
        fprintf(fp, "  \"l_defect\": %d,\n", L_DEFECT_4D);
        fprintf(fp, "  \"r_replicas\": %d,\n", R_REPLICAS);
        fprintf(fp, "  \"eq_work_sweeps\": %d,\n", N_EQ_MEAS);
        fprintf(fp, "  \"exp1_thermal\": [\n");
        int z2_total = N_BETAS_Z2 * 2, z2_idx = 0;
        for (int b = 0; b < N_BETAS_Z2; ++b) {
            for (int s = 0; s < 2; ++s) {
                const z2_res_t *r = &z2_jobs[b].res[s];
                ++z2_idx;
                fprintf(fp, "    {\n");
                fprintf(fp, "      \"beta\": %.2f,\n", r->beta);
                fprintf(fp, "      \"solver\": \"%s\",\n", r->solver);
                fprintf(fp, "      \"mean_p\": %.8f,\n", r->mean_p);
                fprintf(fp, "      \"err_p\": %.8f,\n", r->err_p);
                fprintf(fp, "      \"mean_abs_l\": %.8f,\n", r->mean_abs_l);
                fprintf(fp, "      \"err_abs_l\": %.8f,\n", r->err_abs_l);
                fprintf(fp, "      \"mean_l\": %.8f,\n", r->mean_l);
                fprintf(fp, "      \"err_l\": %.8f,\n", r->err_l);
                fprintf(fp, "      \"tau_l\": %.4f,\n", r->tau_l);
                fprintf(fp, "      \"frozen_l\": %s,\n", r->frozen_l ? "true" : "false");
                fprintf(fp, "      \"z2_sign_flips\": %lld\n", r->z2_sign_flips);
                fprintf(fp, "    }%s\n", (z2_idx < z2_total) ? "," : "");
            }
        }
        fprintf(fp, "  ],\n");
        fprintf(fp, "  \"exp2_confined\": [\n");
        int c_total = N_BETAS_CONF * 3, c_idx = 0;
        for (int b = 0; b < N_BETAS_CONF; ++b) {
            for (int s = 0; s < 3; ++s) {
                const su2_conf_res_t *r = &conf_jobs[b].res[s];
                ++c_idx;
                fprintf(fp, "    {\n");
                fprintf(fp, "      \"beta\": %.2f,\n", r->beta);
                fprintf(fp, "      \"solver\": \"%s\",\n", r->solver);
                fprintf(fp, "      \"n_replicas\": %d,\n", r->n_replicas);
                fprintf(fp, "      \"work_sweeps\": %d,\n", r->work_sweeps);
                fprintf(fp, "      \"mean_p\": %.8f,\n", r->mean_p);
                fprintf(fp, "      \"err_p\": %.8f,\n", r->err_p);
                fprintf(fp, "      \"mean_abs_l\": %.8f,\n", r->mean_abs_l);
                fprintf(fp, "      \"err_abs_l\": %.8f,\n", r->err_abs_l);
                fprintf(fp, "      \"w11\": %.8f,\n", r->w11);
                fprintf(fp, "      \"w12\": %.8f,\n", r->w12);
                fprintf(fp, "      \"w22\": %.8f,\n", r->w22);
                fprintf(fp, "      \"w23\": %.8f,\n", r->w23);
                fprintf(fp, "      \"w33\": %.8f,\n", r->w33);
                fprintf(fp, "      \"chi22\": %.6f,\n", r->chi22);
                fprintf(fp, "      \"err_chi22\": %.6f,\n", r->err_chi22);
                fprintf(fp, "      \"chi33\": %.6f,\n", r->chi33);
                fprintf(fp, "      \"err_chi33\": %.6f,\n", r->err_chi33);
                fprintf(fp, "      \"chi_t_cool\": %.8f,\n", r->chi_t_cool);
                fprintf(fp, "      \"err_chi_t_cool\": %.8f,\n", r->err_chi_t_cool);
                fprintf(fp, "      \"chi_quarter_over_sigma\": %.6f,\n", r->chi_quarter_over_sigma);
                fprintf(fp, "      \"tau_work_q_cool\": %.4f,\n", r->tau_work_q_cool);
                fprintf(fp, "      \"frozen_q_cool\": %s,\n", r->frozen_q_cool ? "true" : "false");
                fprintf(fp, "      \"trans_q_cool\": %lld,\n", r->trans_q_cool);
                fprintf(fp, "      \"bare_inst_acc\": %.6f,\n", r->bare_inst_acc);
                fprintf(fp, "      \"min_swap_acc\": %.6f,\n", r->min_swap_acc);
                fprintf(fp, "      \"mean_swap_acc\": %.6f,\n", r->mean_swap_acc);
                fprintf(fp, "      \"c_ladder\": [");
                for (int k = 0; k < r->n_replicas; ++k) {
                    fprintf(fp, "%.4f%s", r->c_ladder[k], (k + 1 < r->n_replicas) ? ", " : "");
                }
                fprintf(fp, "]\n");
                fprintf(fp, "    }%s\n", (c_idx < c_total) ? "," : "");
            }
        }
        fprintf(fp, "  ]\n");
        fprintf(fp, "}\n");
        fclose(fp);
    }

    return 0;
}
