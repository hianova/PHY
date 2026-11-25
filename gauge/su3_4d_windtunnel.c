/* su3_4d_windtunnel.c — 4D Pure SU(3) Color Yang-Mills Wind Tunnel (Theorems C1–C5)
 * Tests Theorems C1–C5 from `Flow.ColorGaugeBundle`:
 *   - 6 Startup Mathematical & Physics Self-Tests:
 *       1. Infinite-Z^4 SU(3) gauge invariance (|U†U - I| < 1e-13, |det U - 1| < 1e-13, dP < 1e-13)
 *       2. Eisenstein Z_3 center twist invariance (dS < 1e-10, L' = omega*L, adjoint L_8 invariant)
 *       3. Cabibbo-Marinari 3-subgroup Over-Relaxation microcanonical exactness (dS < 1e-10)
 *       4. Single-link Cabibbo-Marinari Heatbath (beta in {3.0, 6.0}, 200k samples vs exact Weyl Haar)
 *       5. Localized 3D defect-wall @swap action identity on 4^4 (dS_wall == dS_full < 1e-10)
 *       6. Gauge-covariant 3D APE smearing & 4D APE cooling + embedded SU(2)->SU(3) BPST Q_cool = +-1
 *   - Experiment 1 (4 x 8^3, beta in {5.55, 5.65, 5.69, 5.73, 6.15}):
 *       Thermal Z_3 first-order deconfinement, symbiotic_z3_center ergodicity restoration,
 *       and cold-vs-hot first-order hysteresis / bimodality check at beta = 5.69
 *   - Experiment 2 (8^4, beta in {5.70, 5.80, 6.00}, 3D APE smeared):
 *       Confinement scaling, fundamental Creutz chi_3(2,2), chi_3(3,3) vs adjoint chi_8(2,2) Casimir ratio
 *   - Experiment 3 (8^4, beta = 5.90, eq_work = 800 sweeps):
 *       4D SU(3) topological susceptibility chi_t & 4D defect-wall @swap PTBC (R = 8, Fisher 2-Jet ladder)
 */
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    /* Experiment 1: Thermal Z_3 First-Order Deconfinement */
    NT_Z3       = 4,
    NS_Z3       = 8,
    N_BETAS_Z3  = 5,
    N_THERM_Z3  = 180,
    N_MEAS_Z3   = 750,

    /* Experiment 2: Confined 4D SU(3) Fundamental vs Adjoint Casimir Scaling */
    L_SCALE     = 8,
    N_BETAS_SC  = 3,
    N_THERM_SC  = 180,
    N_MEAS_SC   = 280,
    SW_PER_SC   = 2,
    W_MAX_F     = 4,
    W_MAX_A     = 3,
    APE3D_STEPS = 4,

    /* Experiment 3: 4D SU(3) Topology & Defect-Wall PTBC */
    L_TOP       = 8,
    L_DEF_SU3   = 2,
    R_REP       = 8,
    N_MEAS_1R   = 80,
    SW_PER_1R   = 10,
    N_THERM_1R  = 140,
    N_THERM_PT  = 160,
    N_MEAS_PT   = 50,
    SW_PER_PT   = 2, /* 8 * 50 * 2 = 800 sweeps == 80 * 10 exact equal work */
    COOL_STEPS  = 12,
    JK_BLOCKS   = 20
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

/* --- Complex & SU(3) / SU(2) Algebra matching Flow.ColorGaugeBundle --- */
typedef struct {
    double re, im;
} cpx_t;

typedef struct {
    cpx_t e[3][3];
} su3_t;

typedef struct {
    double w, x, y, z;
} quat_t;

static inline cpx_t c_mk(double re, double im) {
    cpx_t c = {re, im};
    return c;
}

static inline cpx_t c_add(cpx_t a, cpx_t b) {
    return c_mk(a.re + b.re, a.im + b.im);
}

static inline cpx_t c_sub(cpx_t a, cpx_t b) {
    return c_mk(a.re - b.re, a.im - b.im);
}

static inline cpx_t c_conj(cpx_t a) {
    return c_mk(a.re, -a.im);
}

static inline cpx_t c_mul(cpx_t a, cpx_t b) {
    return c_mk(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}

static inline double c_abs2(cpx_t a) {
    return a.re * a.re + a.im * a.im;
}

static inline double c_abs(cpx_t a) {
    return sqrt(c_abs2(a));
}

static inline quat_t q_id(void) {
    quat_t q = {1.0, 0.0, 0.0, 0.0};
    return q;
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

static su3_t su3_id(void) {
    su3_t m;
    memset(&m, 0, sizeof(m));
    m.e[0][0].re = 1.0;
    m.e[1][1].re = 1.0;
    m.e[2][2].re = 1.0;
    return m;
}

static su3_t su3_zero(void) {
    su3_t m;
    memset(&m, 0, sizeof(m));
    return m;
}

static su3_t su3_mul(const su3_t *a, const su3_t *b) {
    su3_t r;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double sr = 0.0, si = 0.0;
            for (int k = 0; k < 3; ++k) {
                sr += a->e[i][k].re * b->e[k][j].re - a->e[i][k].im * b->e[k][j].im;
                si += a->e[i][k].re * b->e[k][j].im + a->e[i][k].im * b->e[k][j].re;
            }
            r.e[i][j].re = sr;
            r.e[i][j].im = si;
        }
    }
    return r;
}

static su3_t su3_mul_dag(const su3_t *a, const su3_t *b) {
    su3_t r;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double sr = 0.0, si = 0.0;
            for (int k = 0; k < 3; ++k) {
                /* a[i][k] * conj(b[j][k]) */
                sr += a->e[i][k].re * b->e[j][k].re + a->e[i][k].im * b->e[j][k].im;
                si += a->e[i][k].im * b->e[j][k].re - a->e[i][k].re * b->e[j][k].im;
            }
            r.e[i][j].re = sr;
            r.e[i][j].im = si;
        }
    }
    return r;
}

static su3_t su3_dag_mul(const su3_t *a, const su3_t *b) {
    su3_t r;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double sr = 0.0, si = 0.0;
            for (int k = 0; k < 3; ++k) {
                /* conj(a[k][i]) * b[k][j] */
                sr += a->e[k][i].re * b->e[k][j].re + a->e[k][i].im * b->e[k][j].im;
                si += a->e[k][i].re * b->e[k][j].im - a->e[k][i].im * b->e[k][j].re;
            }
            r.e[i][j].re = sr;
            r.e[i][j].im = si;
        }
    }
    return r;
}

static su3_t su3_dag(const su3_t *a) {
    su3_t r;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            r.e[i][j] = c_conj(a->e[j][i]);
        }
    }
    return r;
}

static inline cpx_t su3_tr(const su3_t *a) {
    return c_add(c_add(a->e[0][0], a->e[1][1]), a->e[2][2]);
}

static inline double su3_re_tr(const su3_t *a) {
    return a->e[0][0].re + a->e[1][1].re + a->e[2][2].re;
}

static cpx_t su3_det(const su3_t *m) {
    cpx_t c0 = c_sub(c_mul(m->e[1][1], m->e[2][2]), c_mul(m->e[1][2], m->e[2][1]));
    cpx_t c1 = c_sub(c_mul(m->e[1][0], m->e[2][2]), c_mul(m->e[1][2], m->e[2][0]));
    cpx_t c2 = c_sub(c_mul(m->e[1][0], m->e[2][1]), c_mul(m->e[1][1], m->e[2][0]));
    return c_add(c_sub(c_mul(m->e[0][0], c0), c_mul(m->e[0][1], c1)), c_mul(m->e[0][2], c2));
}

static double su3_unitarity_err(const su3_t *m) {
    su3_t p    = su3_dag_mul(m, m);
    double err = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double tre = (i == j) ? 1.0 : 0.0;
            double d   = hypot(p.e[i][j].re - tre, p.e[i][j].im);
            if (d > err) {
                err = d;
            }
        }
    }
    return err;
}

static void su3_project_gs(su3_t *u) {
    /* Row 0 normalize */
    double n0 = 0.0;
    for (int j = 0; j < 3; ++j) {
        n0 += c_abs2(u->e[0][j]);
    }
    if (n0 < 1e-28) {
        *u = su3_id();
        return;
    }
    double inv0 = 1.0 / sqrt(n0);
    for (int j = 0; j < 3; ++j) {
        u->e[0][j].re *= inv0;
        u->e[0][j].im *= inv0;
    }

    /* Row 1 <- row1 - (row1 . row0*) row0 */
    cpx_t dot = {0.0, 0.0};
    for (int j = 0; j < 3; ++j) {
        dot = c_add(dot, c_mul(u->e[1][j], c_conj(u->e[0][j])));
    }
    for (int j = 0; j < 3; ++j) {
        u->e[1][j] = c_sub(u->e[1][j], c_mul(dot, u->e[0][j]));
    }
    double n1 = 0.0;
    for (int j = 0; j < 3; ++j) {
        n1 += c_abs2(u->e[1][j]);
    }
    if (n1 < 1e-28) {
        *u = su3_id();
        return;
    }
    double inv1 = 1.0 / sqrt(n1);
    for (int j = 0; j < 3; ++j) {
        u->e[1][j].re *= inv1;
        u->e[1][j].im *= inv1;
    }

    /* Row 2 = conj(row0 x row1) -> guarantees U†U = I AND det(U) = +1 */
    u->e[2][0] = c_conj(c_sub(c_mul(u->e[0][1], u->e[1][2]), c_mul(u->e[0][2], u->e[1][1])));
    u->e[2][1] = c_conj(c_sub(c_mul(u->e[0][2], u->e[1][0]), c_mul(u->e[0][0], u->e[1][2])));
    u->e[2][2] = c_conj(c_sub(c_mul(u->e[0][0], u->e[1][1]), c_mul(u->e[0][1], u->e[1][0])));
}

/* Subblock projection matching `subblockQuat` & `embedQuat` in `Flow.ColorGaugeBundle` */
static inline quat_t su3_extract_subblock(const su3_t *w, int a, int b) {
    quat_t v = {
        w->e[a][a].re + w->e[b][b].re,
        -(w->e[a][b].im + w->e[b][a].im),
        w->e[b][a].re - w->e[a][b].re,
        w->e[b][b].im - w->e[a][a].im
    };
    return v;
}

static inline void su3_left_mul_subblock(su3_t *m, int a, int b, quat_t r) {
    cpx_t r00 = c_mk(r.w, r.z);
    cpx_t r01 = c_mk(r.y, r.x);
    cpx_t r10 = c_mk(-r.y, r.x);
    cpx_t r11 = c_mk(r.w, -r.z);
    for (int j = 0; j < 3; ++j) {
        cpx_t ma     = m->e[a][j];
        cpx_t mb     = m->e[b][j];
        m->e[a][j]   = c_add(c_mul(r00, ma), c_mul(r01, mb));
        m->e[b][j]   = c_add(c_mul(r10, ma), c_mul(r11, mb));
    }
}

/* Gauge-covariant SU(3) projection: maximizes Re Tr(U * target†) from covariant U_init */
static su3_t su3_project_max_tr(const su3_t *u_init, const su3_t *target, int sweeps) {
    const int sub_a[3] = {0, 0, 1};
    const int sub_b[3] = {1, 2, 2};
    su3_t u            = *u_init;
    su3_t w            = su3_mul_dag(&u, target);
    for (int sw = 0; sw < sweeps; ++sw) {
        double max_sin2 = 0.0;
        for (int s = 0; s < 3; ++s) {
            quat_t v  = su3_extract_subblock(&w, sub_a[s], sub_b[s]);
            double k2 = q_norm2(v);
            if (k2 > 1e-28) {
                double sin2 = (v.x * v.x + v.y * v.y + v.z * v.z) / k2;
                if (sin2 > max_sin2) {
                    max_sin2 = sin2;
                }
                double inv = 1.0 / sqrt(k2);
                quat_t r   = {v.w * inv, v.x * inv, v.y * inv, v.z * inv};
                su3_left_mul_subblock(&u, sub_a[s], sub_b[s], r);
                su3_left_mul_subblock(&w, sub_a[s], sub_b[s], r);
            }
        }
        if (max_sin2 < 1e-24) {
            break;
        }
    }
    su3_project_gs(&u);
    return u;
}

/* --- 4D SU(3) Lattice --- */
typedef struct {
    int nt;
    int ns;
    int vol;
    su3_t *u;
} su3_lat_t;

static inline int next_i(int i, int n) {
    return (i + 1 < n) ? (i + 1) : 0;
}

static inline int prev_i(int i, int n) {
    return (i > 0) ? (i - 1) : (n - 1);
}

static inline size_t site_idx(const su3_lat_t *lat, int t, int x, int y, int z) {
    return (size_t)(((t * lat->ns + x) * lat->ns + y) * lat->ns + z);
}

static inline size_t link_idx(const su3_lat_t *lat, int t, int x, int y, int z, int mu) {
    return 4U * site_idx(lat, t, x, y, z) + (size_t)mu;
}

static void lat_alloc(su3_lat_t *lat, int nt, int ns) {
    lat->nt  = nt;
    lat->ns  = ns;
    lat->vol = nt * ns * ns * ns;
    lat->u   = (su3_t *)malloc(4U * (size_t)lat->vol * sizeof(su3_t));
    if (lat->u == NULL) {
        fprintf(stderr, "ERROR: OOM allocating 4D SU(3) lattice %dx%d^3\n", nt, ns);
        exit(1);
    }
}

static void lat_free(su3_lat_t *lat) {
    free(lat->u);
    lat->u = NULL;
}

static void lat_copy(su3_lat_t *dst, const su3_lat_t *src) {
    memcpy(dst->u, src->u, 4U * (size_t)src->vol * sizeof(su3_t));
}

static void lat_init_cold(su3_lat_t *lat) {
    for (int i = 0; i < 4 * lat->vol; ++i) {
        lat->u[i] = su3_id();
    }
}

static inline void step_coord(const su3_lat_t *lat, const int c[4], int d, int dir, int out[4]) {
    out[0]  = c[0];
    out[1]  = c[1];
    out[2]  = c[2];
    out[3]  = c[3];
    int lim = (d == 0) ? lat->nt : lat->ns;
    out[d]  = (dir > 0) ? next_i(c[d], lim) : prev_i(c[d], lim);
}

static inline su3_t get_link(const su3_lat_t *lat, const int c[4], int mu) {
    return lat->u[link_idx(lat, c[0], c[1], c[2], c[3], mu)];
}

/* --- Creutz-Kennedy-Pendleton SU(2) Sampler & Haar Generator --- */
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

static su3_t random_su3_haar(rng_t *rng) {
    su3_t u            = su3_id();
    const int sub_a[3] = {0, 0, 1};
    const int sub_b[3] = {1, 2, 2};
    for (int s = 0; s < 3; ++s) {
        quat_t q = sample_su2_boltzmann(0.0, rng);
        su3_left_mul_subblock(&u, sub_a[s], sub_b[s], q);
    }
    su3_project_gs(&u);
    return u;
}

static void lat_init_hot(su3_lat_t *lat, rng_t *rng) {
    for (int i = 0; i < 4 * lat->vol; ++i) {
        lat->u[i] = random_su3_haar(rng);
    }
}

/* --- 4D Defect-Wall Mask & Plaquette / Staple Computations --- */
static inline int is_defect_site(int l_def, int x, int y, int z) {
    return (x >= 0 && x < l_def) && (y >= 0 && y < l_def) && (z >= 0 && z < l_def);
}

static inline int is_defect_plaq(const su3_lat_t *lat, int l_def, const int c[4], int mu_lo, int nu_hi) {
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

static inline double plaq_weight(const su3_lat_t *lat, int l_def, double c_r,
                                 const int c[4], int mu, int nu) {
    int mu_lo = (mu < nu) ? mu : nu;
    int nu_hi = (mu < nu) ? nu : mu;
    return is_defect_plaq(lat, l_def, c, mu_lo, nu_hi) ? c_r : 1.0;
}

static su3_t plaq_su3(const su3_lat_t *lat, const int c[4], int mu, int nu) {
    int c_mu[4], c_nu[4];
    step_coord(lat, c, mu, +1, c_mu);
    step_coord(lat, c, nu, +1, c_nu);
    su3_t u1  = get_link(lat, c, mu);
    su3_t u2  = get_link(lat, c_mu, nu);
    su3_t u3  = get_link(lat, c_nu, mu);
    su3_t u4  = get_link(lat, c, nu);
    su3_t p12 = su3_mul(&u1, &u2);
    su3_t u3d = su3_dag(&u3);
    su3_t p34 = su3_mul_dag(&u3d, &u4);
    return su3_mul(&p12, &p34);
}

/* Staple V = sum_{nu != mu} [ w_f * (U_nu(x+mu) U_mu†(x+nu) U_nu†(x))† + ... ]
 * such that Re Tr(U V†) = sum w_P Re Tr(U_P) */
static su3_t compute_staple_4d(const su3_lat_t *lat, int t, int x, int y, int z, int mu,
                               double c_r, int l_def) {
    const int c[4] = {t, x, y, z};
    int c_mu[4];
    step_coord(lat, c, mu, +1, c_mu);

    su3_t sum = su3_zero();
    for (int nu = 0; nu < 4; ++nu) {
        if (nu == mu) {
            continue;
        }
        int c_nu[4], c_mnu[4], c_mnumu[4];
        step_coord(lat, c, nu, +1, c_nu);
        step_coord(lat, c, nu, -1, c_mnu);
        step_coord(lat, c_mnu, mu, +1, c_mnumu);

        su3_t unu   = get_link(lat, c, nu);
        su3_t umu_n = get_link(lat, c_nu, mu);
        su3_t u_mnu = get_link(lat, c_mu, nu);
        su3_t t1    = su3_mul(&unu, &umu_n);
        su3_t sf    = su3_mul_dag(&t1, &u_mnu);
        double wf   = plaq_weight(lat, l_def, c_r, c, mu, nu);

        su3_t d_nu  = get_link(lat, c_mnu, nu);
        su3_t d_mu  = get_link(lat, c_mnu, mu);
        su3_t d_end = get_link(lat, c_mnumu, nu);
        su3_t t2    = su3_dag_mul(&d_nu, &d_mu);
        su3_t sb    = su3_mul(&t2, &d_end);
        double wb   = plaq_weight(lat, l_def, c_r, c_mnu, mu, nu);

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                sum.e[i][j].re += wf * sf.e[i][j].re + wb * sb.e[i][j].re;
                sum.e[i][j].im += wf * sf.e[i][j].im + wb * sb.e[i][j].im;
            }
        }
    }
    return sum;
}

static double cfg_total_action(const su3_lat_t *lat, double beta, double c_r, int l_def) {
    double s = 0.0;
    for (int t = 0; t < lat->nt; ++t) {
        for (int x = 0; x < lat->ns; ++x) {
            for (int y = 0; y < lat->ns; ++y) {
                for (int z = 0; z < lat->ns; ++z) {
                    const int c[4] = {t, x, y, z};
                    for (int mu = 0; mu < 4; ++mu) {
                        for (int nu = mu + 1; nu < 4; ++nu) {
                            double w = plaq_weight(lat, l_def, c_r, c, mu, nu);
                            su3_t p  = plaq_su3(lat, c, mu, nu);
                            s += w * (1.0 - su3_re_tr(&p) / 3.0);
                        }
                    }
                }
            }
        }
    }
    return beta * s;
}

static double cfg_defect_sum_re3(const su3_lat_t *lat, int l_def, int *out_cnt) {
    double s = 0.0;
    int cnt  = 0;
    int t    = lat->nt / 2;
    for (int x = 0; x < lat->ns; ++x) {
        for (int y = 0; y < lat->ns; ++y) {
            for (int z = 0; z < lat->ns; ++z) {
                const int c[4] = {t, x, y, z};
                for (int nu = 1; nu < 4; ++nu) {
                    if (is_defect_plaq(lat, l_def, c, 0, nu)) {
                        su3_t p = plaq_su3(lat, c, 0, nu);
                        s += su3_re_tr(&p) / 3.0;
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

static double cfg_mean_plaq(const su3_lat_t *lat) {
    double sum_p = 0.0;
    for (int t = 0; t < lat->nt; ++t) {
        for (int x = 0; x < lat->ns; ++x) {
            for (int y = 0; y < lat->ns; ++y) {
                for (int z = 0; z < lat->ns; ++z) {
                    const int c[4] = {t, x, y, z};
                    for (int mu = 0; mu < 4; ++mu) {
                        for (int nu = mu + 1; nu < 4; ++nu) {
                            su3_t p = plaq_su3(lat, c, mu, nu);
                            sum_p += su3_re_tr(&p) / 3.0;
                        }
                    }
                }
            }
        }
    }
    return sum_p / (double)(6 * lat->vol);
}

/* Normalized complex Polyakov loop <1/3 Tr L>, adjoint <1/8 (|Tr L|^2 - 1)> */
static cpx_t cfg_polyakov_mean(const su3_lat_t *lat, double *out_abs_mean, double *out_adj_mean) {
    double sum_re   = 0.0;
    double sum_im   = 0.0;
    double sum_adj  = 0.0;
    int space_sites = lat->ns * lat->ns * lat->ns;
    for (int x = 0; x < lat->ns; ++x) {
        for (int y = 0; y < lat->ns; ++y) {
            for (int z = 0; z < lat->ns; ++z) {
                su3_t p = su3_id();
                for (int t = 0; t < lat->nt; ++t) {
                    su3_t ut = lat->u[link_idx(lat, t, x, y, z, 0)];
                    p        = su3_mul(&p, &ut);
                }
                cpx_t tr  = su3_tr(&p);
                sum_re   += tr.re / 3.0;
                sum_im   += tr.im / 3.0;
                sum_adj  += (c_abs2(tr) - 1.0) / 8.0;
            }
        }
    }
    cpx_t l = c_mk(sum_re / (double)space_sites, sum_im / (double)space_sites);
    if (out_abs_mean != NULL) {
        *out_abs_mean = c_abs(l);
    }
    if (out_adj_mean != NULL) {
        *out_adj_mean = sum_adj / (double)space_sites;
    }
    return l;
}

static void cfg_twist_z3_slice(su3_lat_t *lat, int t0, int k_pow) {
    int km = ((k_pow % 3) + 3) % 3;
    if (km == 0) {
        return;
    }
    double ang = (double)km * (K_TWO_PI / 3.0);
    cpx_t z    = c_mk(cos(ang), sin(ang));
    for (int x = 0; x < lat->ns; ++x) {
        for (int y = 0; y < lat->ns; ++y) {
            for (int zc = 0; zc < lat->ns; ++zc) {
                size_t id = link_idx(lat, t0, x, y, zc, 0);
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        lat->u[id].e[i][j] = c_mul(z, lat->u[id].e[i][j]);
                    }
                }
            }
        }
    }
}

/* --- Cabibbo-Marinari 3-Subgroup Heatbath & Microcanonical Over-Relaxation --- */
static void cm_hb_update_link(su3_t *u_link, const su3_t *staple, double beta, rng_t *rng) {
    const int sub_a[3] = {0, 0, 1};
    const int sub_b[3] = {1, 2, 2};
    su3_t u            = *u_link;
    su3_t w            = su3_mul_dag(&u, staple);
    for (int s = 0; s < 3; ++s) {
        quat_t v  = su3_extract_subblock(&w, sub_a[s], sub_b[s]);
        double k  = sqrt(q_norm2(v));
        quat_t xq = sample_su2_boltzmann((beta / 3.0) * k, rng);
        quat_t r;
        if (k > 1e-14) {
            double inv = 1.0 / k;
            quat_t vh  = {v.w * inv, v.x * inv, v.y * inv, v.z * inv};
            r          = q_mul(xq, vh);
            q_normalize(&r);
        } else {
            r = xq;
        }
        su3_left_mul_subblock(&u, sub_a[s], sub_b[s], r);
        su3_left_mul_subblock(&w, sub_a[s], sub_b[s], r);
    }
    su3_project_gs(&u);
    *u_link = u;
}

static void cm_or_update_link(su3_t *u_link, const su3_t *staple) {
    const int sub_a[3] = {0, 0, 1};
    const int sub_b[3] = {1, 2, 2};
    su3_t u            = *u_link;
    su3_t w            = su3_mul_dag(&u, staple);
    for (int s = 0; s < 3; ++s) {
        quat_t v  = su3_extract_subblock(&w, sub_a[s], sub_b[s]);
        double v2 = q_norm2(v);
        if (v2 > 1e-24) {
            double sc = (2.0 * v.w) / v2;
            quat_t r  = {sc * v.w - 1.0, sc * v.x, sc * v.y, sc * v.z};
            su3_left_mul_subblock(&u, sub_a[s], sub_b[s], r);
            su3_left_mul_subblock(&w, sub_a[s], sub_b[s], r);
        }
    }
    su3_project_gs(&u);
    *u_link = u;
}

static void sweep_su3_hb_or(su3_lat_t *lat, double beta, double c_r, int l_def, rng_t *rng) {
    for (int par = 0; par < 2; ++par) {
        for (int t = 0; t < lat->nt; ++t) {
            for (int x = 0; x < lat->ns; ++x) {
                for (int y = 0; y < lat->ns; ++y) {
                    for (int z = 0; z < lat->ns; ++z) {
                        if (((t + x + y + z) & 1) != par) {
                            continue;
                        }
                        for (int mu = 0; mu < 4; ++mu) {
                            su3_t v   = compute_staple_4d(lat, t, x, y, z, mu, c_r, l_def);
                            size_t id = link_idx(lat, t, x, y, z, mu);
                            cm_hb_update_link(&lat->u[id], &v, beta, rng);
                        }
                    }
                }
            }
        }
    }

    for (int or_sw = 0; or_sw < 1; ++or_sw) {
        for (int par = 0; par < 2; ++par) {
            for (int t = 0; t < lat->nt; ++t) {
                for (int x = 0; x < lat->ns; ++x) {
                    for (int y = 0; y < lat->ns; ++y) {
                        for (int z = 0; z < lat->ns; ++z) {
                            if (((t + x + y + z) & 1) != par) {
                                continue;
                            }
                            for (int mu = 0; mu < 4; ++mu) {
                                su3_t v   = compute_staple_4d(lat, t, x, y, z, mu, c_r, l_def);
                                size_t id = link_idx(lat, t, x, y, z, mu);
                                cm_or_update_link(&lat->u[id], &v);
                            }
                        }
                    }
                }
            }
        }
    }
}

static void shift_torus_4d_random(su3_lat_t *lat, su3_lat_t *tmp, rng_t *rng) {
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

/* --- 3D Time-Slice-Local APE Smearing & Cached Wilson Loop Measurer --- */
static void smear_3d_ape_step(su3_lat_t *dst, const su3_lat_t *src, double c_ape) {
    for (int t = 0; t < src->nt; ++t) {
        for (int x = 0; x < src->ns; ++x) {
            for (int y = 0; y < src->ns; ++y) {
                for (int z = 0; z < src->ns; ++z) {
                    const int c[4] = {t, x, y, z};
                    /* Preserve temporal link untouched for exact OS positivity */
                    dst->u[link_idx(dst, t, x, y, z, 0)] = src->u[link_idx(src, t, x, y, z, 0)];

                    for (int mu = 1; mu <= 3; ++mu) {
                        int c_mu[4];
                        step_coord(src, c, mu, +1, c_mu);
                        su3_t u0  = get_link(src, c, mu);
                        su3_t sum = su3_zero();
                        for (int nu = 1; nu <= 3; ++nu) {
                            if (nu == mu) {
                                continue;
                            }
                            int c_nu[4], c_mnu[4], c_mnumu[4];
                            step_coord(src, c, nu, +1, c_nu);
                            step_coord(src, c, nu, -1, c_mnu);
                            step_coord(src, c_mnu, mu, +1, c_mnumu);

                            su3_t unu   = get_link(src, c, nu);
                            su3_t umu_n = get_link(src, c_nu, mu);
                            su3_t u_mnu = get_link(src, c_mu, nu);
                            su3_t t1    = su3_mul(&unu, &umu_n);
                            su3_t sf    = su3_mul_dag(&t1, &u_mnu);

                            su3_t d_nu  = get_link(src, c_mnu, nu);
                            su3_t d_mu  = get_link(src, c_mnu, mu);
                            su3_t d_end = get_link(src, c_mnumu, nu);
                            su3_t t2    = su3_dag_mul(&d_nu, &d_mu);
                            su3_t sb    = su3_mul(&t2, &d_end);

                            for (int i = 0; i < 3; ++i) {
                                for (int j = 0; j < 3; ++j) {
                                    sum.e[i][j].re += sf.e[i][j].re + sb.e[i][j].re;
                                    sum.e[i][j].im += sf.e[i][j].im + sb.e[i][j].im;
                                }
                            }
                        }
                        su3_t target;
                        for (int i = 0; i < 3; ++i) {
                            for (int j = 0; j < 3; ++j) {
                                target.e[i][j].re = (1.0 - c_ape) * u0.e[i][j].re + (c_ape / 4.0) * sum.e[i][j].re;
                                target.e[i][j].im = (1.0 - c_ape) * u0.e[i][j].im + (c_ape / 4.0) * sum.e[i][j].im;
                            }
                        }
                        dst->u[link_idx(dst, t, x, y, z, mu)] = su3_project_max_tr(&u0, &target, 12);
                    }
                }
            }
        }
    }
}

typedef struct {
    su3_t *sp_line[3][W_MAX_F];
    su3_t *tm_line[W_MAX_F];
} wloop_ws_t;

static void wloop_ws_alloc(wloop_ws_t *ws, int vol) {
    for (int d = 0; d < 3; ++d) {
        for (int r = 0; r < W_MAX_F; ++r) {
            ws->sp_line[d][r] = (su3_t *)malloc((size_t)vol * sizeof(su3_t));
            if (ws->sp_line[d][r] == NULL) {
                fprintf(stderr, "ERROR: OOM wloop_ws spatial\n");
                exit(1);
            }
        }
    }
    for (int r = 0; r < W_MAX_F; ++r) {
        ws->tm_line[r] = (su3_t *)malloc((size_t)vol * sizeof(su3_t));
        if (ws->tm_line[r] == NULL) {
            fprintf(stderr, "ERROR: OOM wloop_ws temporal\n");
            exit(1);
        }
    }
}

static void wloop_ws_free(wloop_ws_t *ws) {
    for (int d = 0; d < 3; ++d) {
        for (int r = 0; r < W_MAX_F; ++r) {
            free(ws->sp_line[d][r]);
        }
    }
    for (int r = 0; r < W_MAX_F; ++r) {
        free(ws->tm_line[r]);
    }
}

static void measure_wilson_loops(const su3_lat_t *lat, wloop_ws_t *ws,
                                 double wf[W_MAX_F][W_MAX_F],
                                 double wa[W_MAX_A][W_MAX_A]) {
    for (int t = 0; t < lat->nt; ++t) {
        for (int x = 0; x < lat->ns; ++x) {
            for (int y = 0; y < lat->ns; ++y) {
                for (int z = 0; z < lat->ns; ++z) {
                    size_t sid = site_idx(lat, t, x, y, z);
                    for (int d = 0; d < 3; ++d) {
                        ws->sp_line[d][0][sid] = lat->u[4U * sid + (size_t)(d + 1)];
                    }
                    ws->tm_line[0][sid] = lat->u[4U * sid];
                }
            }
        }
    }
    for (int len = 1; len < W_MAX_F; ++len) {
        for (int t = 0; t < lat->nt; ++t) {
            for (int x = 0; x < lat->ns; ++x) {
                for (int y = 0; y < lat->ns; ++y) {
                    for (int z = 0; z < lat->ns; ++z) {
                        size_t sid = site_idx(lat, t, x, y, z);
                        int dt     = (t + len) % lat->nt;
                        size_t n_t = site_idx(lat, dt, x, y, z);
                        ws->tm_line[len][sid] = su3_mul(&ws->tm_line[len - 1][sid], &ws->tm_line[0][n_t]);
                        for (int d = 0; d < 3; ++d) {
                            int c_s[4] = {t, x, y, z};
                            c_s[d + 1] = (c_s[d + 1] + len) % lat->ns;
                            size_t n_s = site_idx(lat, c_s[0], c_s[1], c_s[2], c_s[3]);
                            ws->sp_line[d][len][sid] =
                                su3_mul(&ws->sp_line[d][len - 1][sid], &ws->sp_line[d][0][n_s]);
                        }
                    }
                }
            }
        }
    }

    for (int r = 0; r < W_MAX_F; ++r) {
        for (int tl = 0; tl < W_MAX_F; ++tl) {
            double acc_f = 0.0;
            double acc_a = 0.0;
            int do_adj   = (r < W_MAX_A && tl < W_MAX_A);
            for (int t = 0; t < lat->nt; ++t) {
                int ts = (t + tl + 1) % lat->nt;
                for (int x = 0; x < lat->ns; ++x) {
                    for (int y = 0; y < lat->ns; ++y) {
                        for (int z = 0; z < lat->ns; ++z) {
                            size_t s0 = site_idx(lat, t, x, y, z);
                            size_t st = site_idx(lat, ts, x, y, z);
                            for (int d = 0; d < 3; ++d) {
                                int cs[4] = {t, x, y, z};
                                cs[d + 1] = (cs[d + 1] + r + 1) % lat->ns;
                                size_t sr = site_idx(lat, cs[0], cs[1], cs[2], cs[3]);

                                su3_t w1  = su3_mul(&ws->sp_line[d][r][s0], &ws->tm_line[tl][sr]);
                                su3_t w2  = su3_mul(&ws->tm_line[tl][s0], &ws->sp_line[d][r][st]);
                                su3_t lp  = su3_mul_dag(&w1, &w2);
                                cpx_t tr  = su3_tr(&lp);
                                acc_f    += tr.re / 3.0;
                                if (do_adj) {
                                    acc_a += (c_abs2(tr) - 1.0) / 8.0;
                                }
                            }
                        }
                    }
                }
            }
            double denom = (double)(3 * lat->vol);
            wf[r][tl]    = acc_f / denom;
            if (do_adj) {
                wa[r][tl] = acc_a / denom;
            }
        }
    }
}

/* --- 4D Gauge-Covariant Clover Topological Charge & Damped APE Cooling --- */
static su3_t clover_field_su3(const su3_lat_t *lat, const int c[4], int mu, int nu) {
    int cp_mu[4], cp_nu[4], cm_mu[4], cm_nu[4], cm_mu_p_nu[4], cm_mu_m_nu[4], cp_mu_m_nu[4];
    step_coord(lat, c, mu, +1, cp_mu);
    step_coord(lat, c, nu, +1, cp_nu);
    step_coord(lat, c, mu, -1, cm_mu);
    step_coord(lat, c, nu, -1, cm_nu);
    step_coord(lat, cm_mu, nu, +1, cm_mu_p_nu);
    step_coord(lat, cm_mu, nu, -1, cm_mu_m_nu);
    step_coord(lat, cp_mu, nu, -1, cp_mu_m_nu);

    su3_t u1, u2, u3, u4, t1, t2;

    /* Leaf 1 (+mu, +nu) */
    u1       = get_link(lat, c, mu);
    u2       = get_link(lat, cp_mu, nu);
    u3       = get_link(lat, cp_nu, mu);
    u4       = get_link(lat, c, nu);
    t1       = su3_mul(&u1, &u2);
    t2       = su3_mul(&u4, &u3);
    su3_t l1 = su3_mul_dag(&t1, &t2);

    /* Leaf 2 (+nu, -mu) */
    u1       = get_link(lat, c, nu);
    u2       = get_link(lat, cm_mu_p_nu, mu);
    u3       = get_link(lat, cm_mu, nu);
    u4       = get_link(lat, cm_mu, mu);
    t1       = su3_mul_dag(&u1, &u2);
    t2       = su3_dag_mul(&u3, &u4);
    su3_t l2 = su3_mul(&t1, &t2);

    /* Leaf 3 (-mu, -nu) */
    u1       = get_link(lat, cm_mu, mu);
    u2       = get_link(lat, cm_mu_m_nu, nu);
    u3       = get_link(lat, cm_mu_m_nu, mu);
    u4       = get_link(lat, cm_nu, nu);
    t1       = su3_mul(&u2, &u1);
    t2       = su3_mul(&u3, &u4);
    su3_t l3 = su3_dag_mul(&t1, &t2);

    /* Leaf 4 (-nu, +mu) */
    u1       = get_link(lat, cm_nu, nu);
    u2       = get_link(lat, cm_nu, mu);
    u3       = get_link(lat, cp_mu_m_nu, nu);
    u4       = get_link(lat, c, mu);
    t1       = su3_dag_mul(&u1, &u2);
    t2       = su3_mul_dag(&u3, &u4);
    su3_t l4 = su3_mul(&t1, &t2);

    /* Traceless anti-Hermitian part C = 1/8 sum (L_k - L_k†) - Tr(...)/3 * I */
    su3_t f;
    cpx_t tr_acc = {0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            cpx_t s = c_add(c_add(l1.e[i][j], l2.e[i][j]), c_add(l3.e[i][j], l4.e[i][j]));
            cpx_t d = c_add(c_add(c_conj(l1.e[j][i]), c_conj(l2.e[j][i])),
                            c_add(c_conj(l3.e[j][i]), c_conj(l4.e[j][i])));
            f.e[i][j].re = 0.125 * (s.re - d.re);
            f.e[i][j].im = 0.125 * (s.im - d.im);
        }
        tr_acc = c_add(tr_acc, f.e[i][i]);
    }
    cpx_t tr3 = c_mk(tr_acc.re / 3.0, tr_acc.im / 3.0);
    for (int i = 0; i < 3; ++i) {
        f.e[i][i] = c_sub(f.e[i][i], tr3);
    }
    return f;
}

static double cfg_topo_charge_clover(const su3_lat_t *lat) {
    double sum_q      = 0.0;
    const double norm = -1.0 / (4.0 * K_PI * K_PI);
    for (int t = 0; t < lat->nt; ++t) {
        for (int x = 0; x < lat->ns; ++x) {
            for (int y = 0; y < lat->ns; ++y) {
                for (int z = 0; z < lat->ns; ++z) {
                    const int c[4] = {t, x, y, z};
                    su3_t f01 = clover_field_su3(lat, c, 0, 1);
                    su3_t f23 = clover_field_su3(lat, c, 2, 3);
                    su3_t f02 = clover_field_su3(lat, c, 0, 2);
                    su3_t f13 = clover_field_su3(lat, c, 1, 3);
                    su3_t f03 = clover_field_su3(lat, c, 0, 3);
                    su3_t f12 = clover_field_su3(lat, c, 1, 2);

                    su3_t p1 = su3_mul(&f01, &f23);
                    su3_t p2 = su3_mul(&f02, &f13);
                    su3_t p3 = su3_mul(&f03, &f12);
                    sum_q += su3_re_tr(&p1) - su3_re_tr(&p2) + su3_re_tr(&p3);
                }
            }
        }
    }
    return sum_q * norm;
}

static void cfg_cool_sweeps(const su3_lat_t *src, su3_lat_t *wrk, int steps) {
    const double eps = 0.28;
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
                                su3_t v   = compute_staple_4d(wrk, t, x, y, z, mu, 1.0, 0);
                                size_t id = link_idx(wrk, t, x, y, z, mu);
                                su3_t u0  = wrk->u[id];
                                su3_t target;
                                for (int i = 0; i < 3; ++i) {
                                    for (int j = 0; j < 3; ++j) {
                                        target.e[i][j].re =
                                            (1.0 - eps) * u0.e[i][j].re + (eps / 6.0) * v.e[i][j].re;
                                        target.e[i][j].im =
                                            (1.0 - eps) * u0.e[i][j].im + (eps / 6.0) * v.e[i][j].im;
                                    }
                                }
                                wrk->u[id] = su3_project_max_tr(&u0, &target, 12);
                            }
                        }
                    }
                }
            }
        }
    }
}

static long long cfg_topo_charge_cooled(const su3_lat_t *src, su3_lat_t *wrk, double *out_raw) {
    cfg_cool_sweeps(src, wrk, COOL_STEPS);
    double q = cfg_topo_charge_clover(wrk);
    if (out_raw != NULL) {
        *out_raw = q;
    }
    return llround(q);
}

/* Embedded SU(2) -> SU(3) Smoothstep Singular-Gauge BPST Instanton for Clover Self-Test */
static inline quat_t q_conj(quat_t a) {
    quat_t r = {a.w, -a.x, -a.y, -a.z};
    return r;
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

static void apply_embedded_bpst_su3(su3_lat_t *dst, const su3_lat_t *src,
                                    int t0, int x0, int y0, int z0, double rho, int sign) {
    const double r_in      = 2.8;
    const double r_out     = 3.9;
    const quat_t basis[4]  = {
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
                        double s      = sqrt(perp2);
                        double ang    = (atan(c1[mu] / s) - atan(c0[mu] / s)) / s;
                        quat_t q      = (sign > 0) ? q_mul(qc0, q_conj(basis[mu]))
                                                   : q_mul(q_conj(qc0), basis[mu]);
                        quat_t u_reg  = q_exp_im(ang * q.x, ang * q.y, ang * q.z);
                        quat_t g1     = g_bpst_at(c1, sign);
                        quat_t u_sing = q_mul(q_mul(g0, u_reg), q_conj(g1));
                        double tap    = 1.0;
                        if (rm > r_in) {
                            double u = 1.0 - (rm - r_in) / (r_out - r_in);
                            tap      = u * u * (3.0 - 2.0 * u);
                        }
                        double w_cl = (u_sing.w > 1.0) ? 1.0 : ((u_sing.w < -1.0) ? -1.0 : u_sing.w);
                        double th   = acos(w_cl) * tap;
                        double vn   = sqrt(u_sing.x * u_sing.x + u_sing.y * u_sing.y + u_sing.z * u_sing.z);
                        quat_t u2   = q_id();
                        if (vn > 1e-12) {
                            double sc = sin(th) / vn;
                            u2.w      = cos(th);
                            u2.x      = u_sing.x * sc;
                            u2.y      = u_sing.y * sc;
                            u2.z      = u_sing.z * sc;
                        }
                        size_t id = link_idx(dst, t, x, y, z, mu);
                        su3_left_mul_subblock(&dst->u[id], 0, 1, q_conj(u2));
                    }
                }
            }
        }
    }
}

/* --- Statistics Helpers: Madras-Sokal, Jackknife & Fisher 2-Jet Ladder --- */
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

static void compute_creutz_jackknife(const double *waa, const double *wab,
                                     const double *wba, const double *wbb,
                                     int n, double *out_chi, double *out_err) {
    double saa = 0.0, sab = 0.0, sba = 0.0, sbb = 0.0;
    for (int i = 0; i < n; ++i) {
        saa += waa[i];
        sab += wab[i];
        sba += wba[i];
        sbb += wbb[i];
    }
    double den_full = sab * sba;
    double ratio    = (den_full > 1e-18) ? ((sbb * saa) / den_full) : 0.0;
    double chi_full = (ratio > 1e-12) ? (-log(ratio)) : 0.0;

    int bsz = n / JK_BLOCKS;
    double jk[JK_BLOCKS];
    double jk_mean = 0.0;
    for (int b = 0; b < JK_BLOCKS; ++b) {
        double bsaa = 0.0, bsab = 0.0, bsba = 0.0, bsbb = 0.0;
        int lo = b * bsz;
        int hi = (b == JK_BLOCKS - 1) ? n : ((b + 1) * bsz);
        for (int i = lo; i < hi; ++i) {
            bsaa += waa[i];
            bsab += wab[i];
            bsba += wba[i];
            bsbb += wbb[i];
        }
        double d_jk = (sab - bsab) * (sba - bsba);
        double r_jk = (d_jk > 1e-18) ? (((sbb - bsbb) * (saa - bsaa)) / d_jk) : 0.0;
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

static void adapt_fisher_ladder(double c[R_REP], const double avg_p[R_REP - 1]) {
    double ds[R_REP - 1];
    double total_s = 0.0;
    for (int r = 0; r < R_REP - 1; ++r) {
        double p = avg_p[r];
        if (p < 0.01) p = 0.01;
        if (p > 0.99) p = 0.99;
        ds[r] = sqrt(-log(p));
        total_s += ds[r];
    }
    if (total_s <= 1e-12) {
        return;
    }
    double cum_c[R_REP], cum_s[R_REP];
    cum_s[0] = 0.0;
    for (int r = 0; r < R_REP; ++r) {
        cum_c[r] = c[r];
        if (r > 0) {
            cum_s[r] = cum_s[r - 1] + ds[r - 1];
        }
    }
    double c_new[R_REP];
    c_new[0]         = 1.0;
    c_new[R_REP - 1] = 0.0;
    for (int r = 1; r < R_REP - 1; ++r) {
        double s_tgt = ((double)r / (double)(R_REP - 1)) * total_s;
        int k        = 0;
        while (k + 1 < R_REP - 1 && cum_s[k + 1] < s_tgt) {
            ++k;
        }
        double denom  = cum_s[k + 1] - cum_s[k];
        double t      = (denom > 1e-12) ? ((s_tgt - cum_s[k]) / denom) : 0.5;
        double interp = cum_c[k] + t * (cum_c[k + 1] - cum_c[k]);
        c_new[r]      = 0.5 * c[r] + 0.5 * interp;
    }
    for (int r = 1; r < R_REP - 1; ++r) {
        if (c_new[r] > c_new[r - 1] - 0.01) {
            c_new[r] = c_new[r - 1] - 0.01;
        }
        double floor_min = 0.01 * (double)(R_REP - 1 - r);
        if (c_new[r] < floor_min) {
            c_new[r] = floor_min;
        }
    }
    for (int r = 0; r < R_REP; ++r) {
        c[r] = c_new[r];
    }
}

/* --- 6 Startup Mathematical & Physics Self-Tests --- */
typedef struct {
    double max_unit_err;
    double max_det_err;
    double max_gauge_dP;
    double z3_dS;
    double z3_rot_err;
    double z3_adj_err;
    double or_dS;
    double or_inv_err;
    double hb_b3_mean;
    double hb_b3_err;
    double hb_b6_mean;
    double hb_b6_err;
    double swap_dS_err;
    double ape3d_gauge_err;
    double cool_gauge_err;
    double bpst_pos_q;
    double bpst_neg_q;
    long long bpst_pos_round;
    long long bpst_neg_round;
} selftest_res_t;

static void run_startup_selftests(selftest_res_t *out) {
    su3_lat_t g, gp, wrk, sm0, sm1;
    lat_alloc(&g, 4, 4);
    lat_alloc(&gp, 4, 4);
    lat_alloc(&wrk, 4, 4);
    lat_alloc(&sm0, 4, 4);
    lat_alloc(&sm1, 4, 4);

    rng_t rng;
    rng_seed(&rng, 0x5034D2026ULL);
    lat_init_hot(&g, &rng);

    /* 1. Random local SU(3) gauge invariance + unitarity + det = 1 */
    double max_u = 0.0, max_det = 0.0;
    su3_t *gmap = (su3_t *)malloc((size_t)g.vol * sizeof(su3_t));
    if (gmap == NULL) {
        fprintf(stderr, "ERROR: OOM selftest gmap\n");
        exit(1);
    }
    for (int i = 0; i < g.vol; ++i) {
        gmap[i] = random_su3_haar(&rng);
    }
    for (int i = 0; i < 4 * g.vol; ++i) {
        double ue = su3_unitarity_err(&g.u[i]);
        cpx_t dt  = su3_det(&g.u[i]);
        double de = hypot(dt.re - 1.0, dt.im);
        if (ue > max_u) max_u = ue;
        if (de > max_det) max_det = de;
    }
    for (int t = 0; t < g.nt; ++t) {
        for (int x = 0; x < g.ns; ++x) {
            for (int y = 0; y < g.ns; ++y) {
                for (int z = 0; z < g.ns; ++z) {
                    const int c[4] = {t, x, y, z};
                    su3_t gs       = gmap[site_idx(&g, t, x, y, z)];
                    for (int mu = 0; mu < 4; ++mu) {
                        int cp[4];
                        step_coord(&g, c, mu, +1, cp);
                        su3_t gsp = gmap[site_idx(&g, cp[0], cp[1], cp[2], cp[3])];
                        su3_t u0  = get_link(&g, c, mu);
                        su3_t gu  = su3_mul(&gs, &u0);
                        gp.u[link_idx(&gp, t, x, y, z, mu)] = su3_mul_dag(&gu, &gsp);
                    }
                }
            }
        }
    }
    free(gmap);

    double a0 = cfg_total_action(&g, 6.0, 0.43, L_DEF_SU3);
    double ap = cfg_total_action(&gp, 6.0, 0.43, L_DEF_SU3);
    double p0 = cfg_mean_plaq(&g), pp = cfg_mean_plaq(&gp);
    if (max_u > 1e-13 || max_det > 1e-13 || fabs(p0 - pp) > 1e-13 || fabs(a0 - ap) > 1e-10) {
        fprintf(stderr, "SELFTEST 1 FAIL: u=%.3e det=%.3e dP=%.3e dS=%.3e\n",
                max_u, max_det, fabs(p0 - pp), fabs(a0 - ap));
        exit(1);
    }

    /* 2. Z_3 center twist invariance (omega = e^{2*pi*i/3}) */
    lat_copy(&wrk, &g);
    double l0_abs = 0.0, l0_adj = 0.0;
    cpx_t l0 = cfg_polyakov_mean(&wrk, &l0_abs, &l0_adj);
    double s_orig = cfg_total_action(&wrk, 6.0, 1.0, 0);
    cfg_twist_z3_slice(&wrk, 1, 1);
    double l1_abs = 0.0, l1_adj = 0.0;
    cpx_t l1      = cfg_polyakov_mean(&wrk, &l1_abs, &l1_adj);
    double s_z3   = cfg_total_action(&wrk, 6.0, 1.0, 0);
    cpx_t omega   = c_mk(cos(K_TWO_PI / 3.0), sin(K_TWO_PI / 3.0));
    cpx_t om_l0   = c_mul(omega, l0);
    double rot_er = hypot(l1.re - om_l0.re, l1.im - om_l0.im);
    double adj_er = fabs(l1_adj - l0_adj);
    double z3_dS  = fabs(s_z3 - s_orig);
    if (z3_dS > 1e-10 || rot_er > 1e-13 || adj_er > 1e-13) {
        fprintf(stderr, "SELFTEST 2 FAIL: z3 dS=%.3e rot=%.3e adj=%.3e\n", z3_dS, rot_er, adj_er);
        exit(1);
    }

    /* 3. Cabibbo-Marinari Over-Relaxation exactness (dS = 0 and involution OR^2 on 1-link) */
    lat_copy(&wrk, &g);
    for (int par = 0; par < 2; ++par) {
        for (int t = 0; t < wrk.nt; ++t) {
            for (int x = 0; x < wrk.ns; ++x) {
                for (int y = 0; y < wrk.ns; ++y) {
                    for (int z = 0; z < wrk.ns; ++z) {
                        if (((t + x + y + z) & 1) != par) continue;
                        for (int mu = 0; mu < 4; ++mu) {
                            su3_t v   = compute_staple_4d(&wrk, t, x, y, z, mu, 0.43, L_DEF_SU3);
                            size_t id = link_idx(&wrk, t, x, y, z, mu);
                            cm_or_update_link(&wrk.u[id], &v);
                        }
                    }
                }
            }
        }
    }
    double s_or  = cfg_total_action(&wrk, 6.0, 0.43, L_DEF_SU3);
    double or_dS = fabs(s_or - a0);
    /* Single-subgroup OR involution test */
    su3_t u_s0 = random_su3_haar(&rng);
    su3_t v_st = compute_staple_4d(&g, 0, 0, 0, 0, 0, 1.0, 0);
    su3_t u_s2 = u_s0;
    for (int rep = 0; rep < 2; ++rep) {
        su3_t w  = su3_mul_dag(&u_s2, &v_st);
        quat_t v = su3_extract_subblock(&w, 0, 1);
        double sc = (2.0 * v.w) / q_norm2(v);
        quat_t r  = {sc * v.w - 1.0, sc * v.x, sc * v.y, sc * v.z};
        su3_left_mul_subblock(&u_s2, 0, 1, r);
    }
    double inv_err = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double d = hypot(u_s2.e[i][j].re - u_s0.e[i][j].re, u_s2.e[i][j].im - u_s0.e[i][j].im);
            if (d > inv_err) inv_err = d;
        }
    }
    if (or_dS > 1e-10 || inv_err > 1e-13) {
        fprintf(stderr, "SELFTEST 3 FAIL: or_dS=%.3e inv=%.3e\n", or_dS, inv_err);
        exit(1);
    }

    /* 4. Single-link Cabibbo-Marinari Heatbath at beta in {3.0, 6.0} vs Weyl Haar integral */
    const int n_hb_s      = 200000;
    const double betas[2] = {3.0, 6.0};
    double hb_m[2], hb_e[2];
    for (int bi = 0; bi < 2; ++bi) {
        su3_t u1    = su3_id();
        su3_t st_id = su3_id();
        for (int w = 0; w < 1000; ++w) {
            cm_hb_update_link(&u1, &st_id, betas[bi], &rng);
        }
        const int nb  = 200;
        const int bsz = n_hb_s / nb;
        double bmeans[200];
        double tot = 0.0;
        for (int b = 0; b < nb; ++b) {
            double bs = 0.0;
            for (int k = 0; k < bsz; ++k) {
                cm_hb_update_link(&u1, &st_id, betas[bi], &rng);
                bs += su3_re_tr(&u1) / 3.0;
            }
            bmeans[b] = bs / (double)bsz;
            tot += bmeans[b];
        }
        double mn  = tot / (double)nb;
        double bv  = 0.0;
        for (int b = 0; b < nb; ++b) {
            double d = bmeans[b] - mn;
            bv += d * d;
        }
        hb_m[bi] = mn;
        hb_e[bi] = sqrt(bv / (double)(nb * (nb - 1)));
    }

    /* 5. Defect-wall @swap identity on 4^4 */
    lat_init_hot(&wrk, &rng);
    double c1 = 0.72, c2 = 0.31, beta_sw = 6.0;
    double sa1 = cfg_total_action(&g, beta_sw, c1, L_DEF_SU3);
    double sb2 = cfg_total_action(&wrk, beta_sw, c2, L_DEF_SU3);
    double sa2 = cfg_total_action(&g, beta_sw, c2, L_DEF_SU3);
    double sb1 = cfg_total_action(&wrk, beta_sw, c1, L_DEF_SU3);
    int cnt_a = 0, cnt_b = 0;
    double da = cfg_defect_sum_re3(&g, L_DEF_SU3, &cnt_a);
    double db = cfg_defect_sum_re3(&wrk, L_DEF_SU3, &cnt_b);
    double dS_full = (sa2 + sb1) - (sa1 + sb2);
    double dS_fast = beta_sw * (c1 - c2) * (da - db);
    double sw_err  = fabs(dS_full - dS_fast);
    if (cnt_a != 36 || cnt_a != cnt_b || sw_err > 1e-10) {
        fprintf(stderr, "SELFTEST 5 FAIL: cnt=%d/%d sw_err=%.3e\n", cnt_a, cnt_b, sw_err);
        exit(1);
    }

    /* 6. Gauge-covariant 3D APE smearing, 4D cooling, and 8^4 embedded SU(2)->SU(3) BPST Q_cool = +-1 */
    smear_3d_ape_step(&sm0, &g, 0.5);
    smear_3d_ape_step(&sm1, &gp, 0.5);
    double ape3d_err = fabs(cfg_mean_plaq(&sm0) - cfg_mean_plaq(&sm1));
    cfg_cool_sweeps(&g, &sm0, 2);
    cfg_cool_sweeps(&gp, &sm1, 2);
    double cool_err = fabs(cfg_topo_charge_clover(&sm0) - cfg_topo_charge_clover(&sm1));

    su3_lat_t bvac, bins, bwrk;
    lat_alloc(&bvac, 8, 8);
    lat_alloc(&bins, 8, 8);
    lat_alloc(&bwrk, 8, 8);
    lat_init_cold(&bvac);

    apply_embedded_bpst_su3(&bins, &bvac, 4, 4, 4, 4, 1.85, +1);
    double qp_raw    = 0.0;
    long long qp_rnd = cfg_topo_charge_cooled(&bins, &bwrk, &qp_raw);

    apply_embedded_bpst_su3(&bins, &bvac, 4, 4, 4, 4, 1.85, -1);
    double qm_raw    = 0.0;
    long long qm_rnd = cfg_topo_charge_cooled(&bins, &bwrk, &qm_raw);

    lat_free(&bvac);
    lat_free(&bins);
    lat_free(&bwrk);
    lat_free(&g);
    lat_free(&gp);
    lat_free(&wrk);
    lat_free(&sm0);
    lat_free(&sm1);

    if (ape3d_err > 1e-11 || cool_err > 1e-10 || qp_rnd != 1 || qm_rnd != -1) {
        fprintf(stderr, "SELFTEST 6 FAIL: ape=%.3e cool=%.3e Q+=%lld(%.3f) Q-=%lld(%.3f)\n",
                ape3d_err, cool_err, qp_rnd, qp_raw, qm_rnd, qm_raw);
        exit(1);
    }

    out->max_unit_err    = max_u;
    out->max_det_err     = max_det;
    out->max_gauge_dP    = fabs(p0 - pp);
    out->z3_dS           = z3_dS;
    out->z3_rot_err      = rot_er;
    out->z3_adj_err      = adj_er;
    out->or_dS           = or_dS;
    out->or_inv_err      = inv_err;
    out->hb_b3_mean      = hb_m[0];
    out->hb_b3_err       = hb_e[0];
    out->hb_b6_mean      = hb_m[1];
    out->hb_b6_err       = hb_e[1];
    out->swap_dS_err     = sw_err;
    out->ape3d_gauge_err = ape3d_err;
    out->cool_gauge_err  = cool_err;
    out->bpst_pos_q      = qp_raw;
    out->bpst_neg_q      = qm_raw;
    out->bpst_pos_round  = qp_rnd;
    out->bpst_neg_round  = qm_rnd;

    printf("  [SELFTEST] PASS (unit=%.2e det=%.2e gauge=%.2e Z3_dS=%.2e rot=%.2e OR_inv=%.2e "
           "CM_b3=%.5f CM_b6=%.5f @swap=%.2e ape3D=%.2e BPST Q=%+lld[%.3f]/%+lld[%.3f])\n",
           max_u, max_det, out->max_gauge_dP, z3_dS, rot_er, inv_err,
           hb_m[0], hb_m[1], sw_err, ape3d_err, qp_rnd, qp_raw, qm_rnd, qm_raw);
    fflush(stdout);
}

/* --- Experiment 1: 4 x 8^3 Thermal Z_3 First-Order Deconfinement --- */
typedef struct {
    int idx;
    double beta;
    uint64_t seed;
    double p_loc, p_loc_err, abs_l_loc, abs_l_loc_err, mag_mean_l_loc;
    double re_l_loc, re_l_loc_err, im_l_loc, tau_rel_loc;
    double cos3_loc, adj_l_loc;
    double sec_loc[3];
    int hops_loc, frozen_loc;

    double p_sym, p_sym_err, abs_l_sym, abs_l_sym_err, mag_mean_l_sym;
    double re_l_sym, re_l_sym_err, im_l_sym, tau_rel_sym;
    double cos3_sym, adj_l_sym;
    double sec_sym[3];
    int hops_sym;
} exp1_res_t;

static int classify_z3_sector(cpx_t l) {
    double a = atan2(l.im, l.re);
    if (a >= -K_PI / 3.0 && a < K_PI / 3.0) return 0;
    if (a >= K_PI / 3.0) return 1;
    return 2;
}

static void *exp1_worker(void *arg) {
    exp1_res_t *res = (exp1_res_t *)arg;
    su3_lat_t lat;
    lat_alloc(&lat, NT_Z3, NS_Z3);
    rng_t rng;

    double *ps   = (double *)malloc((size_t)N_MEAS_Z3 * sizeof(double));
    double *abls = (double *)malloc((size_t)N_MEAS_Z3 * sizeof(double));
    double *rels = (double *)malloc((size_t)N_MEAS_Z3 * sizeof(double));
    double *imls = (double *)malloc((size_t)N_MEAS_Z3 * sizeof(double));
    double *c3s  = (double *)malloc((size_t)N_MEAS_Z3 * sizeof(double));
    double *adjs = (double *)malloc((size_t)N_MEAS_Z3 * sizeof(double));
    if (!ps || !abls || !rels || !imls || !c3s || !adjs) {
        fprintf(stderr, "ERROR: OOM in exp1\n");
        exit(1);
    }

    for (int mode = 0; mode < 2; ++mode) {
        lat_init_cold(&lat);
        rng_seed(&rng, res->seed + (uint64_t)mode * 999983ULL);
        for (int sw = 0; sw < N_THERM_Z3; ++sw) {
            sweep_su3_hb_or(&lat, res->beta, 1.0, 0, &rng);
            if (mode == 1) {
                int k = 1 + rng_below(&rng, 2); /* {1, 2} non-identity Z_3 twist */
                cfg_twist_z3_slice(&lat, 0, k);
            }
        }
        int hops     = 0;
        int prev_sec = -1;
        int scnt[3]  = {0, 0, 0};
        for (int sw = 0; sw < N_MEAS_Z3; ++sw) {
            sweep_su3_hb_or(&lat, res->beta, 1.0, 0, &rng);
            if (mode == 1) {
                int k = rng_below(&rng, 3); /* uniform Z_3 twist on {0, 1, 2} */
                cfg_twist_z3_slice(&lat, 0, k);
            }
            ps[sw]        = cfg_mean_plaq(&lat);
            double al = 0.0, adj = 0.0;
            cpx_t l   = cfg_polyakov_mean(&lat, &al, &adj);
            abls[sw]  = al;
            rels[sw]  = l.re;
            imls[sw]  = l.im;
            adjs[sw]  = adj;
            double ph = atan2(l.im, l.re);
            c3s[sw]   = cos(3.0 * ph);
            int sec   = classify_z3_sector(l);
            ++scnt[sec];
            if (prev_sec >= 0 && sec != prev_sec) {
                ++hops;
            }
            prev_sec = sec;
        }

        double vp, mp, ep, tp, val, mal, eal, tal, vrel, mrel, erel, trel;
        double vim, mim, eim, tim, vc3, mc3, ec3, tc3, vadj, madj, eadj, tadj;
        int fp = 0, fal = 0, frel = 0, fim = 0, fc3 = 0, fadj = 0;
        compute_madras_sokal(ps, N_MEAS_Z3, &mp, &vp, &tp, &ep, &fp);
        compute_madras_sokal(abls, N_MEAS_Z3, &mal, &val, &tal, &eal, &fal);
        compute_madras_sokal(rels, N_MEAS_Z3, &mrel, &vrel, &trel, &erel, &frel);
        compute_madras_sokal(imls, N_MEAS_Z3, &mim, &vim, &tim, &eim, &fim);
        compute_madras_sokal(c3s, N_MEAS_Z3, &mc3, &vc3, &tc3, &ec3, &fc3);
        compute_madras_sokal(adjs, N_MEAS_Z3, &madj, &vadj, &tadj, &eadj, &fadj);

        if (mode == 0) {
            res->p_loc          = mp;
            res->p_loc_err      = ep;
            res->abs_l_loc      = mal;
            res->abs_l_loc_err  = eal;
            res->re_l_loc       = mrel;
            res->re_l_loc_err   = erel;
            res->im_l_loc       = mim;
            res->mag_mean_l_loc = hypot(mrel, mim);
            res->tau_rel_loc    = trel;
            res->cos3_loc       = mc3;
            res->adj_l_loc      = madj;
            res->hops_loc       = hops;
            res->frozen_loc     = frel;
            for (int k = 0; k < 3; ++k) res->sec_loc[k] = (double)scnt[k] / (double)N_MEAS_Z3;
        } else {
            res->p_sym          = mp;
            res->p_sym_err      = ep;
            res->abs_l_sym      = mal;
            res->abs_l_sym_err  = eal;
            res->re_l_sym       = mrel;
            res->re_l_sym_err   = erel;
            res->im_l_sym       = mim;
            res->mag_mean_l_sym = hypot(mrel, mim);
            res->tau_rel_sym    = trel;
            res->cos3_sym       = mc3;
            res->adj_l_sym      = madj;
            res->hops_sym       = hops;
            for (int k = 0; k < 3; ++k) res->sec_sym[k] = (double)scnt[k] / (double)N_MEAS_Z3;
        }
    }

    free(ps);
    free(abls);
    free(rels);
    free(imls);
    free(c3s);
    free(adjs);
    lat_free(&lat);
    return NULL;
}

typedef struct {
    double cold_p, cold_p_err, cold_abs_l, cold_abs_l_err;
    double hot_p, hot_p_err, hot_abs_l, hot_abs_l_err;
    double mixed_frac_lo, mixed_frac_hi;
} exp1_hyst_t;

static void run_exp1_hysteresis_569(exp1_hyst_t *out) {
    su3_lat_t lat;
    lat_alloc(&lat, NT_Z3, NS_Z3);
    rng_t rng;
    const double beta = 5.69;
    enum { N_HYST_SEEDS = 4 };

    /* Pre-tunneling metastability at beta = 5.69 (4 seeds x [25 relax + 100 meas sweeps]) */
    for (int st = 0; st < 2; ++st) {
        double seed_p[N_HYST_SEEDS];
        double seed_l[N_HYST_SEEDS];
        for (int sd = 0; sd < N_HYST_SEEDS; ++sd) {
            rng_seed(&rng, 0x11457000ULL + (uint64_t)st * 7717ULL + (uint64_t)sd * 104729ULL);
            if (st == 0) {
                lat_init_cold(&lat);
            } else {
                lat_init_hot(&lat, &rng);
            }
            for (int sw = 0; sw < 25; ++sw) {
                sweep_su3_hb_or(&lat, beta, 1.0, 0, &rng);
            }
            double sum_p = 0.0, sum_l = 0.0;
            for (int sw = 0; sw < 100; ++sw) {
                sweep_su3_hb_or(&lat, beta, 1.0, 0, &rng);
                sum_p += cfg_mean_plaq(&lat);
                double al = 0.0;
                cfg_polyakov_mean(&lat, &al, NULL);
                sum_l += al;
            }
            seed_p[sd] = sum_p / 100.0;
            seed_l[sd] = sum_l / 100.0;
        }
        double mp = 0.0, ml = 0.0;
        for (int sd = 0; sd < N_HYST_SEEDS; ++sd) {
            mp += seed_p[sd];
            ml += seed_l[sd];
        }
        mp /= (double)N_HYST_SEEDS;
        ml /= (double)N_HYST_SEEDS;
        double vp = 0.0, vl = 0.0;
        for (int sd = 0; sd < N_HYST_SEEDS; ++sd) {
            double dp = seed_p[sd] - mp;
            double dl = seed_l[sd] - ml;
            vp += dp * dp;
            vl += dl * dl;
        }
        double ep = sqrt(vp / ((double)(N_HYST_SEEDS - 1) * (double)N_HYST_SEEDS));
        double el = sqrt(vl / ((double)(N_HYST_SEEDS - 1) * (double)N_HYST_SEEDS));
        if (st == 0) {
            out->cold_p         = mp;
            out->cold_p_err     = ep;
            out->cold_abs_l     = ml;
            out->cold_abs_l_err = el;
        } else {
            out->hot_p         = mp;
            out->hot_p_err     = ep;
            out->hot_abs_l     = ml;
            out->hot_abs_l_err = el;
        }
    }

    /* Long 1,200-sweep coexistence run at beta = 5.69 to sample both confined and deconfined phases */
    const int n_coex = 1200;
    double *l_long   = (double *)malloc((size_t)n_coex * sizeof(double));
    if (!l_long) exit(1);
    lat_init_cold(&lat);
    rng_seed(&rng, 0xC0E569ULL);
    for (int sw = 0; sw < 200; ++sw) {
        sweep_su3_hb_or(&lat, beta, 1.0, 0, &rng);
    }
    int c_lo = 0, c_hi = 0;
    for (int sw = 0; sw < n_coex; ++sw) {
        sweep_su3_hb_or(&lat, beta, 1.0, 0, &rng);
        double al = 0.0;
        cfg_polyakov_mean(&lat, &al, NULL);
        l_long[sw] = al;
        if (al < 0.09) ++c_lo;
        if (al > 0.14) ++c_hi;
    }
    out->mixed_frac_lo = (double)c_lo / (double)n_coex;
    out->mixed_frac_hi = (double)c_hi / (double)n_coex;

    free(l_long);
    lat_free(&lat);
}

static void *exp1_hyst_worker(void *arg) {
    run_exp1_hysteresis_569((exp1_hyst_t *)arg);
    return NULL;
}

/* --- Experiment 2: 8^4 Confined Scaling & Casimir Ratio chi_8 / chi_3 --- */
typedef struct {
    int idx;
    double beta;
    uint64_t seed;
    double p_mean, p_err, abs_l;
    double wf_mean[W_MAX_F][W_MAX_F];
    double wa_mean[W_MAX_A][W_MAX_A];
    double chi3_22, chi3_22_err;
    double chi3_33, chi3_33_err;
    double chi8_22, chi8_22_err;
    double casimir_ratio, casimir_ratio_err;
} exp2_res_t;

static void *exp2_worker(void *arg) {
    exp2_res_t *res = (exp2_res_t *)arg;
    su3_lat_t lat, sm_a, sm_b;
    lat_alloc(&lat, L_SCALE, L_SCALE);
    lat_alloc(&sm_a, L_SCALE, L_SCALE);
    lat_alloc(&sm_b, L_SCALE, L_SCALE);
    wloop_ws_t ws;
    wloop_ws_alloc(&ws, lat.vol);

    rng_t rng;
    rng_seed(&rng, res->seed);
    lat_init_cold(&lat);
    for (int sw = 0; sw < N_THERM_SC; ++sw) {
        sweep_su3_hb_or(&lat, res->beta, 1.0, 0, &rng);
    }

    double *ps    = (double *)malloc((size_t)N_MEAS_SC * sizeof(double));
    double *ls    = (double *)malloc((size_t)N_MEAS_SC * sizeof(double));
    if (!ps || !ls) { fprintf(stderr, "ERROR: OOM exp2_worker\n"); exit(1); }
    double *wf_ts[W_MAX_F][W_MAX_F];
    double *wa_ts[W_MAX_A][W_MAX_A];
    for (int r = 0; r < W_MAX_F; ++r) {
        for (int t = 0; t < W_MAX_F; ++t) {
            wf_ts[r][t] = (double *)malloc((size_t)N_MEAS_SC * sizeof(double));
            if (!wf_ts[r][t]) { fprintf(stderr, "ERROR: OOM wf_ts\n"); exit(1); }
        }
    }
    for (int r = 0; r < W_MAX_A; ++r) {
        for (int t = 0; t < W_MAX_A; ++t) {
            wa_ts[r][t] = (double *)malloc((size_t)N_MEAS_SC * sizeof(double));
            if (!wa_ts[r][t]) { fprintf(stderr, "ERROR: OOM wa_ts\n"); exit(1); }
        }
    }

    for (int m = 0; m < N_MEAS_SC; ++m) {
        for (int s = 0; s < SW_PER_SC; ++s) {
            sweep_su3_hb_or(&lat, res->beta, 1.0, 0, &rng);
        }
        ps[m]     = cfg_mean_plaq(&lat);
        double al = 0.0;
        cfg_polyakov_mean(&lat, &al, NULL);
        ls[m] = al;

        lat_copy(&sm_a, &lat);
        for (int ap = 0; ap < APE3D_STEPS; ++ap) {
            smear_3d_ape_step(&sm_b, &sm_a, 0.5);
            lat_copy(&sm_a, &sm_b);
        }
        double wf[W_MAX_F][W_MAX_F];
        double wa[W_MAX_A][W_MAX_A];
        measure_wilson_loops(&sm_a, &ws, wf, wa);
        for (int r = 0; r < W_MAX_F; ++r) {
            for (int t = 0; t < W_MAX_F; ++t) {
                wf_ts[r][t][m] = wf[r][t];
            }
        }
        for (int r = 0; r < W_MAX_A; ++r) {
            for (int t = 0; t < W_MAX_A; ++t) {
                wa_ts[r][t][m] = wa[r][t];
            }
        }
    }

    double vp, tp, vl, el, tl;
    int fp = 0, fl = 0;
    compute_madras_sokal(ps, N_MEAS_SC, &res->p_mean, &vp, &tp, &res->p_err, &fp);
    compute_madras_sokal(ls, N_MEAS_SC, &res->abs_l, &vl, &tl, &el, &fl);
    for (int r = 0; r < W_MAX_F; ++r) {
        for (int t = 0; t < W_MAX_F; ++t) {
            double s = 0.0;
            for (int m = 0; m < N_MEAS_SC; ++m) s += wf_ts[r][t][m];
            res->wf_mean[r][t] = s / (double)N_MEAS_SC;
        }
    }
    for (int r = 0; r < W_MAX_A; ++r) {
        for (int t = 0; t < W_MAX_A; ++t) {
            double s = 0.0;
            for (int m = 0; m < N_MEAS_SC; ++m) s += wa_ts[r][t][m];
            res->wa_mean[r][t] = s / (double)N_MEAS_SC;
        }
    }

    compute_creutz_jackknife(wf_ts[0][0], wf_ts[0][1], wf_ts[1][0], wf_ts[1][1],
                             N_MEAS_SC, &res->chi3_22, &res->chi3_22_err);
    compute_creutz_jackknife(wf_ts[1][1], wf_ts[1][2], wf_ts[2][1], wf_ts[2][2],
                             N_MEAS_SC, &res->chi3_33, &res->chi3_33_err);
    compute_creutz_jackknife(wa_ts[0][0], wa_ts[0][1], wa_ts[1][0], wa_ts[1][1],
                             N_MEAS_SC, &res->chi8_22, &res->chi8_22_err);

    /* Direct correlated jackknife on the Casimir ratio chi8(2,2) / chi3(2,2) */
    res->casimir_ratio = (res->chi3_22 > 1e-12) ? (res->chi8_22 / res->chi3_22) : 0.0;
    int bsz            = N_MEAS_SC / JK_BLOCKS;
    double rk[JK_BLOCKS];
    double rk_mean = 0.0;
    double sum_f00 = 0.0, sum_f01 = 0.0, sum_f10 = 0.0, sum_f11 = 0.0;
    double sum_a00 = 0.0, sum_a01 = 0.0, sum_a10 = 0.0, sum_a11 = 0.0;
    for (int i = 0; i < N_MEAS_SC; ++i) {
        sum_f00 += wf_ts[0][0][i]; sum_f01 += wf_ts[0][1][i];
        sum_f10 += wf_ts[1][0][i]; sum_f11 += wf_ts[1][1][i];
        sum_a00 += wa_ts[0][0][i]; sum_a01 += wa_ts[0][1][i];
        sum_a10 += wa_ts[1][0][i]; sum_a11 += wa_ts[1][1][i];
    }
    for (int b = 0; b < JK_BLOCKS; ++b) {
        double bf00 = 0.0, bf01 = 0.0, bf10 = 0.0, bf11 = 0.0;
        double ba00 = 0.0, ba01 = 0.0, ba10 = 0.0, ba11 = 0.0;
        int lo = b * bsz;
        int hi = (b == JK_BLOCKS - 1) ? N_MEAS_SC : ((b + 1) * bsz);
        for (int i = lo; i < hi; ++i) {
            bf00 += wf_ts[0][0][i]; bf01 += wf_ts[0][1][i];
            bf10 += wf_ts[1][0][i]; bf11 += wf_ts[1][1][i];
            ba00 += wa_ts[0][0][i]; ba01 += wa_ts[0][1][i];
            ba10 += wa_ts[1][0][i]; ba11 += wa_ts[1][1][i];
        }
        double rf = ((sum_f11 - bf11) * (sum_f00 - bf00)) / ((sum_f01 - bf01) * (sum_f10 - bf10));
        double ra = ((sum_a11 - ba11) * (sum_a00 - ba00)) / ((sum_a01 - ba01) * (sum_a10 - ba10));
        double cf = (rf > 1e-12) ? (-log(rf)) : 0.0;
        double ca = (ra > 1e-12) ? (-log(ra)) : 0.0;
        rk[b]     = (cf > 1e-12) ? (ca / cf) : 0.0;
        rk_mean += rk[b];
    }
    rk_mean /= (double)JK_BLOCKS;
    double rk_var = 0.0;
    for (int b = 0; b < JK_BLOCKS; ++b) {
        double d = rk[b] - rk_mean;
        rk_var += d * d;
    }
    rk_var *= (double)(JK_BLOCKS - 1) / (double)JK_BLOCKS;
    res->casimir_ratio_err = sqrt(rk_var);

    free(ps);
    free(ls);
    for (int r = 0; r < W_MAX_F; ++r) {
        for (int t = 0; t < W_MAX_F; ++t) free(wf_ts[r][t]);
    }
    for (int r = 0; r < W_MAX_A; ++r) {
        for (int t = 0; t < W_MAX_A; ++t) free(wa_ts[r][t]);
    }
    wloop_ws_free(&ws);
    lat_free(&sm_a);
    lat_free(&sm_b);
    lat_free(&lat);
    return NULL;
}

/* --- Experiment 3: 8^4 beta=6.00 Topology & Defect-Wall PTBC --- */
typedef struct {
    double p_loc, p_loc_err, abs_l_loc;
    double chi_t_loc, chi_t_loc_err, q_mean_loc, q2_mean_loc;
    double tau_q_loc, tau_q_work_loc;
    int tun_loc, frozen_loc;

    double p_pt, p_pt_err, abs_l_pt;
    double chi_t_pt, chi_t_pt_err, q_mean_pt, q2_mean_pt;
    double tau_q_pt, tau_q_work_pt;
    int tun_pt, frozen_pt;
    double sw_min, sw_mean, c_ladder[R_REP];
} exp3_res_t;

static void *exp3_local_worker(void *arg) {
    exp3_res_t *res   = (exp3_res_t *)arg;
    const double beta = 5.90;
    su3_lat_t wrk, lat;
    lat_alloc(&wrk, L_TOP, L_TOP);
    lat_alloc(&lat, L_TOP, L_TOP);
    lat_init_cold(&lat);
    rng_t rng;
    rng_seed(&rng, 0x33600001ULL);
    for (int sw = 0; sw < N_THERM_1R; ++sw) {
        sweep_su3_hb_or(&lat, beta, 1.0, 0, &rng);
    }
    double *ps = (double *)malloc((size_t)N_MEAS_1R * sizeof(double));
    double *ls = (double *)malloc((size_t)N_MEAS_1R * sizeof(double));
    double *qs = (double *)malloc((size_t)N_MEAS_1R * sizeof(double));
    double *q2 = (double *)malloc((size_t)N_MEAS_1R * sizeof(double));
    if (!ps || !ls || !qs || !q2) exit(1);

    int tun_loc       = 0;
    long long prev_ql = 0;
    for (int m = 0; m < N_MEAS_1R; ++m) {
        for (int s = 0; s < SW_PER_1R; ++s) {
            sweep_su3_hb_or(&lat, beta, 1.0, 0, &rng);
        }
        ps[m]     = cfg_mean_plaq(&lat);
        double al = 0.0;
        cfg_polyakov_mean(&lat, &al, NULL);
        ls[m]        = al;
        double qraw  = 0.0;
        long long ql = cfg_topo_charge_cooled(&lat, &wrk, &qraw);
        qs[m]        = (double)ql;
        q2[m]        = (double)(ql * ql) / (double)lat.vol;
        if (m > 0 && ql != prev_ql) {
            ++tun_loc;
        }
        prev_ql = ql;
    }

    double vp, tp, vl, el, tl, vq, eq, vq2, mq2, tq2, eq2;
    int fp = 0, fl = 0, fq = 0, fq2 = 0;
    compute_madras_sokal(ps, N_MEAS_1R, &res->p_loc, &vp, &tp, &res->p_loc_err, &fp);
    compute_madras_sokal(ls, N_MEAS_1R, &res->abs_l_loc, &vl, &tl, &el, &fl);
    compute_madras_sokal(qs, N_MEAS_1R, &res->q_mean_loc, &vq, &res->tau_q_loc, &eq, &fq);
    compute_madras_sokal(q2, N_MEAS_1R, &mq2, &vq2, &tq2, &eq2, &fq2);
    res->chi_t_loc      = mq2;
    res->chi_t_loc_err  = eq2;
    res->q2_mean_loc    = mq2 * (double)lat.vol;
    res->tau_q_work_loc = fq ? -1.0 : (res->tau_q_loc * (double)SW_PER_1R);
    res->tun_loc        = tun_loc;
    res->frozen_loc     = fq;
    free(ps);
    free(ls);
    free(qs);
    free(q2);
    lat_free(&lat);
    lat_free(&wrk);
    return NULL;
}

typedef struct {
    su3_lat_t *reps;
    const double *c;
    rng_t *rrng;
    double beta;
    int r_lo;
    int r_hi;
} ptbc_half_arg_t;

static void *ptbc_half_sweep_worker(void *arg) {
    ptbc_half_arg_t *a = (ptbc_half_arg_t *)arg;
    for (int r = a->r_lo; r < a->r_hi; ++r) {
        sweep_su3_hb_or(&a->reps[r], a->beta, a->c[r], L_DEF_SU3, &a->rrng[r]);
    }
    return NULL;
}

static void *exp3_ptbc_worker(void *arg) {
    exp3_res_t *res   = (exp3_res_t *)arg;
    const double beta = 5.90;
    su3_lat_t wrk;
    lat_alloc(&wrk, L_TOP, L_TOP);

    su3_lat_t reps[R_REP];
    double c[R_REP];
    rng_t rrng[R_REP];
    for (int r = 0; r < R_REP; ++r) {
        lat_alloc(&reps[r], L_TOP, L_TOP);
        lat_init_cold(&reps[r]);
        rng_seed(&rrng[r], 0x77760000ULL + (uint64_t)r * 1337ULL);
        c[r] = 1.0 - (double)r / (double)(R_REP - 1);
    }
    rng_t sw_rng;
    rng_seed(&sw_rng, 0x99960000ULL);
    double sw_acc[R_REP - 1];
    long sw_try[R_REP - 1];
    memset(sw_acc, 0, sizeof(sw_acc));
    memset(sw_try, 0, sizeof(sw_try));

    for (int sw = 0; sw < N_THERM_PT; ++sw) {
        ptbc_half_arg_t h_hi = {reps, c, rrng, beta, R_REP / 2, R_REP};
        pthread_t th_half;
        pthread_create(&th_half, NULL, ptbc_half_sweep_worker, &h_hi);
        for (int r = 0; r < R_REP / 2; ++r) {
            sweep_su3_hb_or(&reps[r], beta, c[r], L_DEF_SU3, &rrng[r]);
        }
        pthread_join(th_half, NULL);

        shift_torus_4d_random(&reps[R_REP - 1], &wrk, &rrng[R_REP - 1]);
        int eo = sw & 1;
        for (int r = eo; r < R_REP - 1; r += 2) {
            double da = cfg_defect_sum_re3(&reps[r], L_DEF_SU3, NULL);
            double db = cfg_defect_sum_re3(&reps[r + 1], L_DEF_SU3, NULL);
            double dS = beta * (c[r] - c[r + 1]) * (da - db);
            double p  = (dS <= 0.0) ? 1.0 : exp(-dS);
            sw_acc[r] += p;
            ++sw_try[r];
            if (dS <= 0.0 || rng_pos_u01(&sw_rng) < p) {
                su3_t *tmp    = reps[r].u;
                reps[r].u     = reps[r + 1].u;
                reps[r + 1].u = tmp;
            }
        }
        if ((sw + 1) % 30 == 0 && sw + 1 <= 120) {
            double avg_p[R_REP - 1];
            for (int r = 0; r < R_REP - 1; ++r) {
                avg_p[r] = (sw_try[r] > 0) ? (sw_acc[r] / (double)sw_try[r]) : 0.5;
            }
            adapt_fisher_ladder(c, avg_p);
            memset(sw_acc, 0, sizeof(sw_acc));
            memset(sw_try, 0, sizeof(sw_try));
        }
    }

    double *pp  = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *lp  = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *qp  = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *qp2 = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    if (!pp || !lp || !qp || !qp2) exit(1);

    memset(sw_acc, 0, sizeof(sw_acc));
    memset(sw_try, 0, sizeof(sw_try));
    int tun_pt        = 0;
    long long prev_ql = 0;
    for (int m = 0; m < N_MEAS_PT; ++m) {
        for (int s = 0; s < SW_PER_PT; ++s) {
            ptbc_half_arg_t h_hi = {reps, c, rrng, beta, R_REP / 2, R_REP};
            pthread_t th_half;
            pthread_create(&th_half, NULL, ptbc_half_sweep_worker, &h_hi);
            for (int r = 0; r < R_REP / 2; ++r) {
                sweep_su3_hb_or(&reps[r], beta, c[r], L_DEF_SU3, &rrng[r]);
            }
            pthread_join(th_half, NULL);

            shift_torus_4d_random(&reps[R_REP - 1], &wrk, &rrng[R_REP - 1]);
            int eo = (m * SW_PER_PT + s) & 1;
            for (int r = eo; r < R_REP - 1; r += 2) {
                double da = cfg_defect_sum_re3(&reps[r], L_DEF_SU3, NULL);
                double db = cfg_defect_sum_re3(&reps[r + 1], L_DEF_SU3, NULL);
                double dS = beta * (c[r] - c[r + 1]) * (da - db);
                double p  = (dS <= 0.0) ? 1.0 : exp(-dS);
                sw_acc[r] += p;
                ++sw_try[r];
                if (dS <= 0.0 || rng_pos_u01(&sw_rng) < p) {
                    su3_t *tmp    = reps[r].u;
                    reps[r].u     = reps[r + 1].u;
                    reps[r + 1].u = tmp;
                }
            }
        }
        pp[m]     = cfg_mean_plaq(&reps[0]);
        double al = 0.0;
        cfg_polyakov_mean(&reps[0], &al, NULL);
        lp[m]        = al;
        double qraw  = 0.0;
        long long ql = cfg_topo_charge_cooled(&reps[0], &wrk, &qraw);
        qp[m]        = (double)ql;
        qp2[m]       = (double)(ql * ql) / (double)reps[0].vol;
        if (m > 0 && ql != prev_ql) {
            ++tun_pt;
        }
        prev_ql = ql;
    }

    double vp, tp, vl, el, tl, vq, eq, vq2, mq2, tq2, eq2;
    int fp = 0, fl = 0, fq = 0, fq2 = 0;
    compute_madras_sokal(pp, N_MEAS_PT, &res->p_pt, &vp, &tp, &res->p_pt_err, &fp);
    compute_madras_sokal(lp, N_MEAS_PT, &res->abs_l_pt, &vl, &tl, &el, &fl);
    compute_madras_sokal(qp, N_MEAS_PT, &res->q_mean_pt, &vq, &res->tau_q_pt, &eq, &fq);
    compute_madras_sokal(qp2, N_MEAS_PT, &mq2, &vq2, &tq2, &eq2, &fq2);
    res->chi_t_pt      = mq2;
    res->chi_t_pt_err  = eq2;
    res->q2_mean_pt    = mq2 * (double)reps[0].vol;
    res->tau_q_work_pt = fq ? -1.0 : (res->tau_q_pt * (double)(R_REP * SW_PER_PT));
    res->tun_pt        = tun_pt;
    res->frozen_pt     = fq;

    double smin = 1.0, smean = 0.0;
    for (int r = 0; r < R_REP - 1; ++r) {
        double ar = sw_acc[r] / (double)sw_try[r];
        if (ar < smin) smin = ar;
        smean += ar;
    }
    res->sw_min  = smin;
    res->sw_mean = smean / (double)(R_REP - 1);
    for (int r = 0; r < R_REP; ++r) {
        res->c_ladder[r] = c[r];
    }

    free(pp);
    free(lp);
    free(qp);
    free(qp2);
    for (int r = 0; r < R_REP; ++r) lat_free(&reps[r]);
    lat_free(&wrk);
    return NULL;
}

int main(int argc, char **argv) {
    const char *json_path = "build/su3_4d_bench.json";
    for (int i = 1; i + 1 < argc; ++i) {
        if (strcmp(argv[i], "--json") == 0) {
            json_path = argv[i + 1];
        }
    }

    printf("=== PHY Stage 6: 4D Pure SU(3) Color Yang-Mills Wind Tunnel ===\n");
    selftest_res_t st;
    run_startup_selftests(&st);

    /* Experiment 1: 4 x 8^3 across beta_c(Nt=4) ~ 5.6925 in 5 parallel pthreads + hysteresis pthread */
    const double betas_z3[N_BETAS_Z3] = {5.55, 5.65, 5.69, 5.73, 6.15};
    exp1_res_t exp1[N_BETAS_Z3];
    pthread_t th1[N_BETAS_Z3];
    for (int i = 0; i < N_BETAS_Z3; ++i) {
        memset(&exp1[i], 0, sizeof(exp1[i]));
        exp1[i].idx  = i;
        exp1[i].beta = betas_z3[i];
        exp1[i].seed = 0x10003000ULL + (uint64_t)i * 4099ULL;
        pthread_create(&th1[i], NULL, exp1_worker, &exp1[i]);
    }
    exp1_hyst_t hyst;
    memset(&hyst, 0, sizeof(hyst));
    pthread_t th_hyst;
    pthread_create(&th_hyst, NULL, exp1_hyst_worker, &hyst);
    for (int i = 0; i < N_BETAS_Z3; ++i) {
        pthread_join(th1[i], NULL);
    }
    pthread_join(th_hyst, NULL);

    printf("  [EXP1 4x8^3] Z_3 Deconfinement & Symbiotic Twist Ergodicity:\n");
    for (int i = 0; i < N_BETAS_Z3; ++i) {
        exp1_res_t *e = &exp1[i];
        printf("    beta=%.2f | <P>=%.5f |<L>|=%.4f | loc: |<L>|_vec=%.4f tau=%.1f hops=%d | "
               "sym: |<L>|_vec=%.4f tau=%.2f hops=%d sec=(%.2f,%.2f,%.2f)\n",
               e->beta, e->p_sym, e->abs_l_sym, e->mag_mean_l_loc, e->tau_rel_loc, e->hops_loc,
               e->mag_mean_l_sym, e->tau_rel_sym, e->hops_sym,
               e->sec_sym[0], e->sec_sym[1], e->sec_sym[2]);
    }
    printf("    hyst@5.69 | cold:<P>=%.5f |<L>|=%.4f vs hot:<P>=%.5f |<L>|=%.4f | coex(lo=%.3f,hi=%.3f)\n",
           hyst.cold_p, hyst.cold_abs_l, hyst.hot_p, hyst.hot_abs_l,
           hyst.mixed_frac_lo, hyst.mixed_frac_hi);
    fflush(stdout);

    /* Experiments 2 & 3 concurrently (3 exp2_worker + exp3_local_worker + exp3_ptbc_worker = 5 threads) */
    const double betas_sc[N_BETAS_SC] = {5.70, 5.80, 6.00};
    exp2_res_t exp2[N_BETAS_SC];
    pthread_t th2[N_BETAS_SC];
    for (int i = 0; i < N_BETAS_SC; ++i) {
        memset(&exp2[i], 0, sizeof(exp2[i]));
        exp2[i].idx  = i;
        exp2[i].beta = betas_sc[i];
        exp2[i].seed = 0x20003000ULL + (uint64_t)i * 7919ULL;
        pthread_create(&th2[i], NULL, exp2_worker, &exp2[i]);
    }
    exp3_res_t exp3;
    memset(&exp3, 0, sizeof(exp3));
    pthread_t th3_loc, th3_pt;
    pthread_create(&th3_loc, NULL, exp3_local_worker, &exp3);
    pthread_create(&th3_pt, NULL, exp3_ptbc_worker, &exp3);

    for (int i = 0; i < N_BETAS_SC; ++i) {
        pthread_join(th2[i], NULL);
    }
    pthread_join(th3_loc, NULL);
    pthread_join(th3_pt, NULL);

    printf("  [EXP2 8^4] Fundamental vs Adjoint Casimir Scaling (C2(8)/C2(3) = 2.25):\n");
    for (int i = 0; i < N_BETAS_SC; ++i) {
        exp2_res_t *e = &exp2[i];
        printf("    beta=%.2f | <P>=%.5f |<L>|=%.4f | chi3(2,2)=%.4f+-%.4f chi3(3,3)=%.4f+-%.4f | "
               "chi8(2,2)=%.4f+-%.4f | chi8/chi3=%.3f+-%.3f\n",
               e->beta, e->p_mean, e->abs_l, e->chi3_22, e->chi3_22_err,
               e->chi3_33, e->chi3_33_err, e->chi8_22, e->chi8_22_err,
               e->casimir_ratio, e->casimir_ratio_err);
    }
    fflush(stdout);

    printf("  [EXP3 8^4@5.90] local <P>=%.5f chi_t=%.3e tau_w=%.1f tun=%d | "
           "PTBC <P>=%.5f chi_t=%.3e tau_w=%.1f tun=%d sw_min=%.1f%%\n",
           exp3.p_loc, exp3.chi_t_loc, exp3.tau_q_work_loc, exp3.tun_loc,
           exp3.p_pt, exp3.chi_t_pt, exp3.tau_q_work_pt, exp3.tun_pt, 100.0 * exp3.sw_min);
    fflush(stdout);

    /* Serialize JSON */
    FILE *fp = fopen(json_path, "w");
    if (fp == NULL) {
        fprintf(stderr, "ERROR: Failed to open %s\n", json_path);
        return 1;
    }
    fprintf(fp, "{\n");
    fprintf(fp, "  \"selftests\": {\n");
    fprintf(fp, "    \"max_unit_err\": %.14e,\n", st.max_unit_err);
    fprintf(fp, "    \"max_det_err\": %.14e,\n", st.max_det_err);
    fprintf(fp, "    \"max_gauge_dP\": %.14e,\n", st.max_gauge_dP);
    fprintf(fp, "    \"z3_dS\": %.14e,\n", st.z3_dS);
    fprintf(fp, "    \"z3_rot_err\": %.14e,\n", st.z3_rot_err);
    fprintf(fp, "    \"z3_adj_err\": %.14e,\n", st.z3_adj_err);
    fprintf(fp, "    \"or_dS\": %.14e,\n", st.or_dS);
    fprintf(fp, "    \"or_inv_err\": %.14e,\n", st.or_inv_err);
    fprintf(fp, "    \"cm_hb_beta3_mean\": %.14f,\n", st.hb_b3_mean);
    fprintf(fp, "    \"cm_hb_beta3_err\": %.14e,\n", st.hb_b3_err);
    fprintf(fp, "    \"cm_hb_beta6_mean\": %.14f,\n", st.hb_b6_mean);
    fprintf(fp, "    \"cm_hb_beta6_err\": %.14e,\n", st.hb_b6_err);
    fprintf(fp, "    \"swap_dS_err\": %.14e,\n", st.swap_dS_err);
    fprintf(fp, "    \"ape3d_gauge_err\": %.14e,\n", st.ape3d_gauge_err);
    fprintf(fp, "    \"cool_gauge_err\": %.14e,\n", st.cool_gauge_err);
    fprintf(fp, "    \"bpst_pos_q\": %.6f,\n", st.bpst_pos_q);
    fprintf(fp, "    \"bpst_neg_q\": %.6f,\n", st.bpst_neg_q);
    fprintf(fp, "    \"bpst_pos_round\": %lld,\n", st.bpst_pos_round);
    fprintf(fp, "    \"bpst_neg_round\": %lld\n", st.bpst_neg_round);
    fprintf(fp, "  },\n");

    fprintf(fp, "  \"exp1_z3\": [\n");
    for (int i = 0; i < N_BETAS_Z3; ++i) {
        exp1_res_t *e = &exp1[i];
        fprintf(fp, "    {\n");
        fprintf(fp, "      \"beta\": %.2f,\n", e->beta);
        fprintf(fp, "      \"p_loc\": %.8f, \"p_loc_err\": %.8f,\n", e->p_loc, e->p_loc_err);
        fprintf(fp, "      \"abs_l_loc\": %.8f, \"abs_l_loc_err\": %.8f,\n", e->abs_l_loc, e->abs_l_loc_err);
        fprintf(fp, "      \"mag_mean_l_loc\": %.8f, \"re_l_loc\": %.8f, \"im_l_loc\": %.8f,\n",
                e->mag_mean_l_loc, e->re_l_loc, e->im_l_loc);
        fprintf(fp, "      \"tau_rel_loc\": %.6f, \"cos3_loc\": %.6f, \"adj_l_loc\": %.6f,\n",
                e->tau_rel_loc, e->cos3_loc, e->adj_l_loc);
        fprintf(fp, "      \"sec_loc\": [%.6f, %.6f, %.6f], \"hops_loc\": %d,\n",
                e->sec_loc[0], e->sec_loc[1], e->sec_loc[2], e->hops_loc);
        fprintf(fp, "      \"p_sym\": %.8f, \"p_sym_err\": %.8f,\n", e->p_sym, e->p_sym_err);
        fprintf(fp, "      \"abs_l_sym\": %.8f, \"abs_l_sym_err\": %.8f,\n", e->abs_l_sym, e->abs_l_sym_err);
        fprintf(fp, "      \"mag_mean_l_sym\": %.8f, \"re_l_sym\": %.8f, \"im_l_sym\": %.8f,\n",
                e->mag_mean_l_sym, e->re_l_sym, e->im_l_sym);
        fprintf(fp, "      \"tau_rel_sym\": %.6f, \"cos3_sym\": %.6f, \"adj_l_sym\": %.6f,\n",
                e->tau_rel_sym, e->cos3_sym, e->adj_l_sym);
        fprintf(fp, "      \"sec_sym\": [%.6f, %.6f, %.6f], \"hops_sym\": %d\n",
                e->sec_sym[0], e->sec_sym[1], e->sec_sym[2], e->hops_sym);
        fprintf(fp, "    }%s\n", (i + 1 < N_BETAS_Z3) ? "," : "");
    }
    fprintf(fp, "  ],\n");

    fprintf(fp, "  \"exp1_hyst_569\": {\n");
    fprintf(fp, "    \"cold_p\": %.8f, \"cold_p_err\": %.8f,\n", hyst.cold_p, hyst.cold_p_err);
    fprintf(fp, "    \"cold_abs_l\": %.8f, \"cold_abs_l_err\": %.8f,\n", hyst.cold_abs_l, hyst.cold_abs_l_err);
    fprintf(fp, "    \"hot_p\": %.8f, \"hot_p_err\": %.8f,\n", hyst.hot_p, hyst.hot_p_err);
    fprintf(fp, "    \"hot_abs_l\": %.8f, \"hot_abs_l_err\": %.8f,\n", hyst.hot_abs_l, hyst.hot_abs_l_err);
    fprintf(fp, "    \"coex_frac_lo\": %.6f, \"coex_frac_hi\": %.6f\n",
            hyst.mixed_frac_lo, hyst.mixed_frac_hi);
    fprintf(fp, "  },\n");

    fprintf(fp, "  \"exp2_casimir\": [\n");
    for (int i = 0; i < N_BETAS_SC; ++i) {
        exp2_res_t *e = &exp2[i];
        fprintf(fp, "    {\n");
        fprintf(fp, "      \"beta\": %.2f, \"p_mean\": %.8f, \"p_err\": %.8f, \"abs_l\": %.8f,\n",
                e->beta, e->p_mean, e->p_err, e->abs_l);
        fprintf(fp, "      \"chi3_22\": %.8f, \"chi3_22_err\": %.8f,\n", e->chi3_22, e->chi3_22_err);
        fprintf(fp, "      \"chi3_33\": %.8f, \"chi3_33_err\": %.8f,\n", e->chi3_33, e->chi3_33_err);
        fprintf(fp, "      \"chi8_22\": %.8f, \"chi8_22_err\": %.8f,\n", e->chi8_22, e->chi8_22_err);
        fprintf(fp, "      \"casimir_ratio\": %.8f, \"casimir_ratio_err\": %.8f,\n",
                e->casimir_ratio, e->casimir_ratio_err);
        fprintf(fp, "      \"w3_11\": %.8f, \"w3_22\": %.8f, \"w3_33\": %.8f, \"w8_11\": %.8f, \"w8_22\": %.8f\n",
                e->wf_mean[0][0], e->wf_mean[1][1], e->wf_mean[2][2],
                e->wa_mean[0][0], e->wa_mean[1][1]);
        fprintf(fp, "    }%s\n", (i + 1 < N_BETAS_SC) ? "," : "");
    }
    fprintf(fp, "  ],\n");

    fprintf(fp, "  \"exp3_topology\": {\n");
    fprintf(fp, "    \"beta\": 5.90,\n");
    fprintf(fp, "    \"p_loc\": %.8f, \"p_loc_err\": %.8f, \"abs_l_loc\": %.8f,\n",
            exp3.p_loc, exp3.p_loc_err, exp3.abs_l_loc);
    fprintf(fp, "    \"chi_t_loc\": %.8e, \"chi_t_loc_err\": %.8e, \"q2_mean_loc\": %.8f,\n",
            exp3.chi_t_loc, exp3.chi_t_loc_err, exp3.q2_mean_loc);
    fprintf(fp, "    \"tau_q_work_loc\": %.6f, \"tun_loc\": %d,\n", exp3.tau_q_work_loc, exp3.tun_loc);
    fprintf(fp, "    \"p_pt\": %.8f, \"p_pt_err\": %.8f, \"abs_l_pt\": %.8f,\n",
            exp3.p_pt, exp3.p_pt_err, exp3.abs_l_pt);
    fprintf(fp, "    \"chi_t_pt\": %.8e, \"chi_t_pt_err\": %.8e, \"q2_mean_pt\": %.8f,\n",
            exp3.chi_t_pt, exp3.chi_t_pt_err, exp3.q2_mean_pt);
    fprintf(fp, "    \"tau_q_work_pt\": %.6f, \"tun_pt\": %d,\n", exp3.tau_q_work_pt, exp3.tun_pt);
    fprintf(fp, "    \"sw_min\": %.6f, \"sw_mean\": %.6f,\n", exp3.sw_min, exp3.sw_mean);
    fprintf(fp, "    \"c_ladder\": [");
    for (int r = 0; r < R_REP; ++r) {
        fprintf(fp, "%.4f%s", exp3.c_ladder[r], (r + 1 < R_REP) ? ", " : "");
    }
    fprintf(fp, "]\n");
    fprintf(fp, "  }\n");
    fprintf(fp, "}\n");
    fclose(fp);
    return 0;
}
