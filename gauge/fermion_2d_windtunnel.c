/* fermion_2d_windtunnel.c — 2D/4D Lattice Fermions, Atiyah-Singer Index & Sign Problem (Theorems F1–F5)
 * Tests Theorems F1–F5 from `Flow.FermionBundle`:
 *   - 6 Startup Mathematical & Algebraic Self-Tests:
 *       1. Exact n=6 (64-dim) Fock-space Jordan-Wigner CAR, Pauli exclusion & number idempotence
 *       2. 2D/4D Euclidean Clifford algebra + 2D Wilson gamma_5-Hermiticity & U(1) gauge covariance
 *       3. Nielsen-Ninomiya 2D free 4-corner doubling & Wilson lifting + 4D staggered pi-flux
 *       4. 6x6 Overlap Ginsparg-Wilson relation, normality & spectral circle (< 1e-11)
 *       5. 6x6 Staggered HMC exact fermion force vs. central finite difference (< 1e-6)
 *       6. Lefschetz quadratic saddle Im S(x + i mu/a) = 0 exact cancellation (< 1e-14)
 *   - Experiment 1 (8x8 U(1), Q in {-2,-1,0,+1,+2} smooth & gauge-fluctuated):
 *       Discrete Atiyah-Singer Index Theorem: index(D_ov) == n_+ - n_- == -Q_geo & chirality < 1e-10
 *   - Experiment 2 (6x6 Dynamical N_f = 2 Schwinger Model vs. Quenched):
 *       Exact LU-determinant Leapfrog HMC across beta in {1.5, 2.5, 4.0} & mass scan m in {0.05..0.50}
 *   - Experiment 3 (Finite-Density Sign Problem):
 *       3A: 1D complex U(1) chain (beta=1.0, mu=0.6, L in {8,16,32,64}): naive reweighting vs. positive dual
 *       3B: Non-Gaussian Lefschetz thimble (a=2.0, lam=0.05, mu in {0.5,1.5,2.5,3.5,4.0})
 */
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    JW_MODES       = 6,
    JW_DIM         = 1 << JW_MODES, /* 64 */
    L_IDX          = 8,
    N_IDX_SECTORS  = 5,             /* Q in {-2, -1, 0, +1, +2} */
    N_IDX_CFGS     = 10,            /* 5 smooth + 5 fluctuated */
    L_SCHW         = 6,
    V_SCHW         = L_SCHW * L_SCHW, /* 36 */
    HMC_STEPS      = 10,
    N_THERM_HMC    = 300,
    N_MEAS_HMC     = 800,
    N_BETAS_HMC    = 3,
    N_MASS_HMC     = 4,
    JK_BLOCKS      = 20,
    N_L_CHAIN      = 4,
    N_SAMP_3A      = 120000,
    N_MU_THIM      = 5,
    N_SAMP_3B      = 200000
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
    for (int i = 0; i < 4; ++i) r->s[i] = splitmix64(&st);
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
    return (double)(rng_u64(r) >> 11) * (1.0 / 9007199254740992.0);
}

static inline double rng_pos_u01(rng_t *r) {
    double u = rng_u01(r);
    return (u < 1e-16) ? 1e-16 : u;
}

static double rng_gauss(rng_t *r) {
    double u1 = rng_pos_u01(r);
    double u2 = rng_u01(r);
    return sqrt(-2.0 * log(u1)) * cos(K_TWO_PI * u2);
}

/* --- Complex Arithmetic --- */
typedef struct {
    double re;
    double im;
} cpx_t;

static inline cpx_t c_mk(double re, double im) { cpx_t z = {re, im}; return z; }
static inline cpx_t c_add(cpx_t a, cpx_t b) { return c_mk(a.re + b.re, a.im + b.im); }
static inline cpx_t c_sub(cpx_t a, cpx_t b) { return c_mk(a.re - b.re, a.im - b.im); }
static inline cpx_t c_mul(cpx_t a, cpx_t b) {
    return c_mk(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
static inline cpx_t c_conj(cpx_t a) { return c_mk(a.re, -a.im); }
static inline cpx_t c_scale(cpx_t a, double s) { return c_mk(a.re * s, a.im * s); }
static inline double c_abs2(cpx_t a) { return a.re * a.re + a.im * a.im; }
static inline cpx_t c_div(cpx_t a, cpx_t b) {
    double d = c_abs2(b);
    return c_mk((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
}
static inline cpx_t c_exp_i(double th) { return c_mk(cos(th), sin(th)); }

static double wrap_angle(double th) {
    double w = fmod(th + K_PI, K_TWO_PI);
    if (w < 0.0) w += K_TWO_PI;
    w -= K_PI;
    if (w <= -K_PI) w += K_TWO_PI;
    return w;
}

static double log_bessel_i0(double x) {
    double ax = fabs(x);
    if (ax < 3.75) {
        double t = ax / 3.75;
        double y = t * t;
        double s = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 +
                   y * (0.2659732 + y * (0.0360768 + y * 0.0045813)))));
        return log(s);
    } else {
        double y = 3.75 / ax;
        double s = 0.39894228 + y * (0.01328592 + y * (0.00225319 + y * (-0.00157565 +
                   y * (0.00916281 + y * (-0.02057706 + y * (0.02635537 +
                   y * (-0.01647633 + y * 0.00392377)))))));
        return ax + log(s / sqrt(ax));
    }
}

/* Exact Bessel I_n(x) for integer n >= 0 using Miller downward recurrence */
static double bessel_in(int n, double x) {
    if (n < 0) n = -n;
    if (fabs(x) < 1e-15) return (n == 0) ? 1.0 : 0.0;
    if (n == 0) return exp(log_bessel_i0(x));
    double i0_val = exp(log_bessel_i0(x));
    int m = n + 40 + (int)(2.0 * fabs(x));
    double bip = 0.0, bi = 1.0, b_res = 0.0, bi0 = 0.0;
    double scale = 2.0 / x;
    for (int j = m; j >= 1; --j) {
        double bim = bip + (double)j * scale * bi;
        bip = bi;
        bi  = bim;
        if (fabs(bi) > 1e100) {
            bip *= 1e-100;
            bi  *= 1e-100;
            b_res *= 1e-100;
            bi0   *= 1e-100;
        }
        if (j - 1 == n) b_res = bi;
        if (j - 1 == 0) bi0 = bi;
    }
    return b_res * (i0_val / bi0);
}

/* --- Cyclic Jacobi Complex Hermitian Eigensolver (V * diag(evals) * V† = A_init) --- */
static void jacobi_hermitian_eigen(int n, cpx_t *a, double *evals, cpx_t *evecs) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            evecs[i * n + j] = (i == j) ? c_mk(1.0, 0.0) : c_mk(0.0, 0.0);
        }
        a[i * n + i].im = 0.0;
    }

    const int max_sweeps = 25;
    for (int sw = 0; sw < max_sweeps; ++sw) {
        double off_norm = 0.0;
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                off_norm += c_abs2(a[p * n + q]);
            }
        }
        if (off_norm < 1e-26) break;

        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                cpx_t apq = a[p * n + q];
                double r2 = c_abs2(apq);
                if (r2 < 1e-32) continue;
                double r   = sqrt(r2);
                cpx_t zhat = c_scale(apq, 1.0 / r);
                double app = a[p * n + p].re;
                double aqq = a[q * n + q].re;
                double tau = (app - aqq) / (2.0 * r);
                double t;
                if (tau >= 0.0) {
                    t = 1.0 / (tau + sqrt(1.0 + tau * tau));
                } else {
                    t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
                }
                double c = 1.0 / sqrt(1.0 + t * t);
                double s = t * c;
                cpx_t szh  = c_scale(zhat, s);
                cpx_t szhc = c_scale(c_conj(zhat), s);

                /* Right update columns p, q of A: A' = A * U */
                for (int i = 0; i < n; ++i) {
                    cpx_t aip        = a[i * n + p];
                    cpx_t aiq        = a[i * n + q];
                    a[i * n + p]     = c_add(c_scale(aip, c), c_mul(szhc, aiq));
                    a[i * n + q]     = c_sub(c_scale(aiq, c), c_mul(szh, aip));
                }
                /* Left update rows p, q of A': A'' = U† * A' */
                for (int j = 0; j < n; ++j) {
                    cpx_t apj        = a[p * n + j];
                    cpx_t aqj        = a[q * n + j];
                    a[p * n + j]     = c_add(c_scale(apj, c), c_mul(szh, aqj));
                    a[q * n + j]     = c_sub(c_scale(aqj, c), c_mul(szhc, apj));
                }
                a[p * n + p].im = 0.0;
                a[q * n + q].im = 0.0;
                a[p * n + q]    = c_mk(0.0, 0.0);
                a[q * n + p]    = c_mk(0.0, 0.0);

                /* Update V <- V * U */
                for (int i = 0; i < n; ++i) {
                    cpx_t vip        = evecs[i * n + p];
                    cpx_t viq        = evecs[i * n + q];
                    evecs[i * n + p] = c_add(c_scale(vip, c), c_mul(szhc, viq));
                    evecs[i * n + q] = c_sub(c_scale(viq, c), c_mul(szh, vip));
                }
            }
        }
    }
    for (int i = 0; i < n; ++i) evals[i] = a[i * n + i].re;
}

/* --- Complex Matrix LU Inversion & log |det M| --- */
static double cpx_lu_logdet_inv(int n, const cpx_t *m_in, cpx_t *m_inv_out) {
    cpx_t *lu  = (cpx_t *)malloc((size_t)(n * n) * sizeof(cpx_t));
    int *perm  = (int *)malloc((size_t)n * sizeof(int));
    if (!lu || !perm) { fprintf(stderr, "OOM in LU\n"); exit(1); }
    memcpy(lu, m_in, (size_t)(n * n) * sizeof(cpx_t));
    for (int i = 0; i < n; ++i) perm[i] = i;

    double log_abs_det = 0.0;
    for (int k = 0; k < n; ++k) {
        int piv       = k;
        double max_a2 = c_abs2(lu[k * n + k]);
        for (int i = k + 1; i < n; ++i) {
            double v2 = c_abs2(lu[i * n + k]);
            if (v2 > max_a2) { max_a2 = v2; piv = i; }
        }
        if (max_a2 < 1e-28) {
            fprintf(stderr, "ERROR: singular matrix in LU\n");
            exit(1);
        }
        if (piv != k) {
            for (int j = 0; j < n; ++j) {
                cpx_t tmp      = lu[k * n + j];
                lu[k * n + j]  = lu[piv * n + j];
                lu[piv * n + j] = tmp;
            }
            int tp    = perm[k];
            perm[k]   = perm[piv];
            perm[piv] = tp;
        }
        cpx_t diag = lu[k * n + k];
        log_abs_det += 0.5 * log(c_abs2(diag));
        for (int i = k + 1; i < n; ++i) {
            cpx_t mult   = c_div(lu[i * n + k], diag);
            lu[i * n + k] = mult;
            for (int j = k + 1; j < n; ++j) {
                lu[i * n + j] = c_sub(lu[i * n + j], c_mul(mult, lu[k * n + j]));
            }
        }
    }

    if (m_inv_out != NULL) {
        cpx_t *col = (cpx_t *)malloc((size_t)n * sizeof(cpx_t));
        if (!col) exit(1);
        for (int rhs = 0; rhs < n; ++rhs) {
            for (int i = 0; i < n; ++i) {
                cpx_t s = (perm[i] == rhs) ? c_mk(1.0, 0.0) : c_mk(0.0, 0.0);
                for (int j = 0; j < i; ++j) {
                    s = c_sub(s, c_mul(lu[i * n + j], col[j]));
                }
                col[i] = s;
            }
            for (int i = n - 1; i >= 0; --i) {
                cpx_t s = col[i];
                for (int j = i + 1; j < n; ++j) {
                    s = c_sub(s, c_mul(lu[i * n + j], col[j]));
                }
                col[i] = c_div(s, lu[i * n + i]);
            }
            for (int i = 0; i < n; ++i) m_inv_out[i * n + rhs] = col[i];
        }
        free(col);
    }

    free(lu);
    free(perm);
    return log_abs_det;
}

/* --- Madras-Sokal Autocorrelation --- */
static void madras_sokal(const double *x, int n, double *out_mean, double *out_tau,
                         double *out_err, int *out_frozen) {
    double m = 0.0;
    for (int i = 0; i < n; ++i) m += x[i];
    m /= (double)n;
    double c0 = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = x[i] - m;
        c0 += d * d;
    }
    c0 /= (double)n;
    *out_mean = m;
    if (c0 < 1e-18) {
        *out_tau    = 0.5;
        *out_err    = 0.0;
        *out_frozen = 0;
        return;
    }
    double tau = 0.5;
    int max_t  = n / 4;
    for (int t = 1; t < max_t; ++t) {
        double ct = 0.0;
        for (int i = 0; i < n - t; ++i) ct += (x[i] - m) * (x[i + t] - m);
        ct /= (double)(n - t);
        double rho = ct / c0;
        if (rho <= 0.0 && (double)t > 6.0 * tau) break;
        if (rho > 0.0) tau += rho;
        if ((double)t > 6.0 * tau) break;
    }
    if (tau < 0.5) tau = 0.5;
    *out_tau    = tau;
    *out_err    = sqrt((2.0 * tau * c0) / (double)n);
    *out_frozen = 0;
}

/* --- 2D U(1) Gauge & Fermion Operators --- */
static inline int pbc(int x, int l) {
    int r = x % l;
    return (r < 0) ? (r + l) : r;
}

static void init_gamma2d(cpx_t g[2][2][2], cpx_t g5[2][2]) {
    memset(g, 0, sizeof(cpx_t) * 8);
    memset(g5, 0, sizeof(cpx_t) * 4);
    g[0][0][1] = c_mk(1.0, 0.0);  g[0][1][0] = c_mk(1.0, 0.0);
    g[1][0][1] = c_mk(0.0, -1.0); g[1][1][0] = c_mk(0.0, 1.0);
    g5[0][0]   = c_mk(1.0, 0.0);  g5[1][1]   = c_mk(-1.0, 0.0);
}

static void init_gamma4d(cpx_t g[4][4][4], cpx_t g5[4][4]) {
    memset(g, 0, sizeof(cpx_t) * 64);
    memset(g5, 0, sizeof(cpx_t) * 16);
    g[0][0][2] = c_mk(1.0, 0.0); g[0][1][3] = c_mk(1.0, 0.0);
    g[0][2][0] = c_mk(1.0, 0.0); g[0][3][1] = c_mk(1.0, 0.0);
    g[1][0][3] = c_mk(0.0, 1.0); g[1][1][2] = c_mk(0.0, 1.0);
    g[1][2][1] = c_mk(0.0, -1.0); g[1][3][0] = c_mk(0.0, -1.0);
    g[2][0][3] = c_mk(1.0, 0.0); g[2][1][2] = c_mk(-1.0, 0.0);
    g[2][2][1] = c_mk(-1.0, 0.0); g[2][3][0] = c_mk(1.0, 0.0);
    g[3][0][2] = c_mk(0.0, 1.0); g[3][1][3] = c_mk(0.0, -1.0);
    g[3][2][0] = c_mk(0.0, -1.0); g[3][3][1] = c_mk(0.0, 1.0);
    g5[0][0] = c_mk(1.0, 0.0);  g5[1][1] = c_mk(1.0, 0.0);
    g5[2][2] = c_mk(-1.0, 0.0); g5[3][3] = c_mk(-1.0, 0.0);
}

/* Build 2D Wilson-Dirac (r_w=1 or r_w=0 naive) matrix of dimension 2*L*L
 * t in [0..L-1] is mu=0; x in [0..L-1] is mu=1. */
static void build_wilson_2d(int l, const double *th, double m0, double r_w, int anti_t, cpx_t *d) {
    int vol = l * l;
    int dim = 2 * vol;
    memset(d, 0, (size_t)(dim * dim) * sizeof(cpx_t));

    cpx_t g[2][2][2], g5[2][2];
    init_gamma2d(g, g5);

    double diag_mass = m0 + 2.0 * r_w;
    for (int t = 0; t < l; ++t) {
        for (int x = 0; x < l; ++x) {
            int s0 = t * l + x;
            for (int a = 0; a < 2; ++a) {
                int idx = 2 * s0 + a;
                d[idx * dim + idx] = c_add(d[idx * dim + idx], c_mk(diag_mass, 0.0));
            }
            for (int mu = 0; mu < 2; ++mu) {
                int tn = (mu == 0) ? pbc(t + 1, l) : t;
                int xn = (mu == 1) ? pbc(x + 1, l) : x;
                int s1 = tn * l + xn;
                double bc = (anti_t && mu == 0 && t == l - 1) ? -1.0 : 1.0;
                cpx_t u   = c_scale(c_exp_i(th[2 * s0 + mu]), bc);
                cpx_t ud  = c_conj(u);
                for (int a = 0; a < 2; ++a) {
                    for (int b = 0; b < 2; ++b) {
                        cpx_t id_ab = (a == b) ? c_mk(r_w, 0.0) : c_mk(0.0, 0.0);
                        cpx_t km_ab = c_scale(c_sub(id_ab, g[mu][a][b]), -0.5);
                        cpx_t kp_ab = c_scale(c_add(id_ab, g[mu][a][b]), -0.5);
                        int r_f = 2 * s0 + a, c_f = 2 * s1 + b;
                        int r_b = 2 * s1 + a, c_b = 2 * s0 + b;
                        d[r_f * dim + c_f] = c_add(d[r_f * dim + c_f], c_mul(km_ab, u));
                        d[r_b * dim + c_b] = c_add(d[r_b * dim + c_b], c_mul(kp_ab, ud));
                    }
                }
            }
        }
    }
}

/* Construct Hermitian H_W = gamma_5 D_W, compute Overlap D_ov = I + gamma_5 sgn(H_W) */
static void construct_overlap_2d(int l, const double *th, double m0, cpx_t *d_ov, double *out_tr_hs) {
    int vol = l * l;
    int dim = 2 * vol;
    cpx_t *dw  = (cpx_t *)malloc((size_t)(dim * dim) * sizeof(cpx_t));
    cpx_t *hw  = (cpx_t *)malloc((size_t)(dim * dim) * sizeof(cpx_t));
    cpx_t *vec = (cpx_t *)malloc((size_t)(dim * dim) * sizeof(cpx_t));
    double *ev = (double *)malloc((size_t)dim * sizeof(double));
    if (!dw || !hw || !vec || !ev) exit(1);

    build_wilson_2d(l, th, m0, 1.0, 1, dw);
    for (int i = 0; i < dim; ++i) {
        double sg = ((i & 1) == 0) ? 1.0 : -1.0;
        for (int j = 0; j < dim; ++j) {
            hw[i * dim + j] = c_scale(dw[i * dim + j], sg);
        }
    }
    jacobi_hermitian_eigen(dim, hw, ev, vec);

    cpx_t *hs = (cpx_t *)calloc((size_t)(dim * dim), sizeof(cpx_t));
    if (!hs) exit(1);
    double tr_hs = 0.0;
    for (int k = 0; k < dim; ++k) {
        double sgn_k = (ev[k] >= 0.0) ? 1.0 : -1.0;
        tr_hs += sgn_k;
        for (int i = 0; i < dim; ++i) {
            cpx_t vik = c_scale(vec[i * dim + k], sgn_k);
            for (int j = 0; j < dim; ++j) {
                hs[i * dim + j] = c_add(hs[i * dim + j], c_mul(vik, c_conj(vec[j * dim + k])));
            }
        }
    }
    if (out_tr_hs != NULL) *out_tr_hs = tr_hs;

    for (int i = 0; i < dim; ++i) {
        double sg = ((i & 1) == 0) ? 1.0 : -1.0;
        for (int j = 0; j < dim; ++j) {
            cpx_t val        = c_scale(hs[i * dim + j], sg);
            d_ov[i * dim + j] = (i == j) ? c_add(c_mk(1.0, 0.0), val) : val;
        }
    }
    free(hs);
    free(dw);
    free(hw);
    free(vec);
    free(ev);
}

/* Build 2D harmonic instanton with exact winding Q_geo = round(sum wrap(theta_p)/(2pi)) */
static void build_harmonic_instanton_2d(int l, double q_tgt, double *th) {
    double f = (K_TWO_PI * q_tgt) / (double)(l * l);
    for (int t = 0; t < l; ++t) {
        for (int x = 0; x < l; ++x) {
            int s          = t * l + x;
            th[2 * s + 0]  = (t == l - 1) ? (-f * (double)l * (double)x) : 0.0;
            th[2 * s + 1]  = f * (double)t;
        }
    }
}

static long long u1_topo_charge_2d(int l, const double *th) {
    double sum = 0.0;
    for (int t = 0; t < l; ++t) {
        for (int x = 0; x < l; ++x) {
            int s0  = t * l + x;
            int st1 = pbc(t + 1, l) * l + x;
            int sx1 = t * l + pbc(x + 1, l);
            double tp = th[2 * s0 + 0] + th[2 * st1 + 1] - th[2 * sx1 + 0] - th[2 * s0 + 1];
            sum += wrap_angle(tp);
        }
    }
    return llround(sum / K_TWO_PI);
}

static double u1_mean_plaq_2d(int l, const double *th) {
    double sum = 0.0;
    for (int t = 0; t < l; ++t) {
        for (int x = 0; x < l; ++x) {
            int s0  = t * l + x;
            int st1 = pbc(t + 1, l) * l + x;
            int sx1 = t * l + pbc(x + 1, l);
            double tp = th[2 * s0 + 0] + th[2 * st1 + 1] - th[2 * sx1 + 0] - th[2 * s0 + 1];
            sum += cos(tp);
        }
    }
    return sum / (double)(l * l);
}

/* Measure Creutz chi(2,2) from Wilson loops W(1,1), W(1,2), W(2,1), W(2,2) */
static void u1_wilson_loops_2d(int l, const double *th, double *w11, double *w12, double *w21, double *w22) {
    cpx_t *u0 = (cpx_t *)malloc((size_t)(l * l) * sizeof(cpx_t));
    cpx_t *u1 = (cpx_t *)malloc((size_t)(l * l) * sizeof(cpx_t));
    if (!u0 || !u1) exit(1);
    for (int s = 0; s < l * l; ++s) {
        u0[s] = c_exp_i(th[2 * s + 0]);
        u1[s] = c_exp_i(th[2 * s + 1]);
    }
    double s11 = 0.0, s12 = 0.0, s21 = 0.0, s22 = 0.0;
    for (int t = 0; t < l; ++t) {
        for (int x = 0; x < l; ++x) {
            for (int r = 1; r <= 2; ++r) {
                for (int k = 1; k <= 2; ++k) {
                    cpx_t prod = c_mk(1.0, 0.0);
                    int ct = t, cx = x;
                    for (int i = 0; i < r; ++i) { prod = c_mul(prod, u0[ct * l + cx]); ct = pbc(ct + 1, l); }
                    for (int i = 0; i < k; ++i) { prod = c_mul(prod, u1[ct * l + cx]); cx = pbc(cx + 1, l); }
                    for (int i = 0; i < r; ++i) { ct = pbc(ct - 1, l); prod = c_mul(prod, c_conj(u0[ct * l + cx])); }
                    for (int i = 0; i < k; ++i) { cx = pbc(cx - 1, l); prod = c_mul(prod, c_conj(u1[ct * l + cx])); }
                    if (r == 1 && k == 1) s11 += prod.re;
                    if (r == 1 && k == 2) s12 += prod.re;
                    if (r == 2 && k == 1) s21 += prod.re;
                    if (r == 2 && k == 2) s22 += prod.re;
                }
            }
        }
    }
    double v = (double)(l * l);
    *w11 = s11 / v; *w12 = s12 / v; *w21 = s21 / v; *w22 = s22 / v;
    free(u0);
    free(u1);
}

/* Staggered Dirac operator D_stag(l, th, m) on V = l*l */
static void build_staggered_2d(int l, const double *th, double m, cpx_t *d) {
    int v = l * l;
    memset(d, 0, (size_t)(v * v) * sizeof(cpx_t));
    for (int t = 0; t < l; ++t) {
        for (int x = 0; x < l; ++x) {
            int s0 = t * l + x;
            d[s0 * v + s0] = c_mk(m, 0.0);
            for (int mu = 0; mu < 2; ++mu) {
                int tn    = (mu == 0) ? pbc(t + 1, l) : t;
                int xn    = (mu == 1) ? pbc(x + 1, l) : x;
                int s1    = tn * l + xn;
                double eta = (mu == 0) ? 1.0 : ((t & 1) ? -1.0 : 1.0);
                double bc  = (mu == 0 && t == l - 1) ? -1.0 : 1.0;
                cpx_t u    = c_scale(c_exp_i(th[2 * s0 + mu]), 0.5 * eta * bc);
                d[s0 * v + s1] = c_add(d[s0 * v + s1], u);
                d[s1 * v + s0] = c_sub(d[s1 * v + s0], c_conj(u));
            }
        }
    }
}

/* Exact action S = S_gauge + S_fermion (N_f = 2 staggered species: S_f = -2 log |det D_stag|)
 * Force = -dS/dtheta */
static double eval_schwinger_action_force(int l, const double *th, double beta, double m,
                                          int dynamical, double *force_out, double *out_pbp) {
    int v  = l * l;
    double sg = 0.0;
    double sin_p[L_SCHW * L_SCHW];
    for (int t = 0; t < l; ++t) {
        for (int x = 0; x < l; ++x) {
            int s0  = t * l + x;
            int st1 = pbc(t + 1, l) * l + x;
            int sx1 = t * l + pbc(x + 1, l);
            double tp = th[2 * s0 + 0] + th[2 * st1 + 1] - th[2 * sx1 + 0] - th[2 * s0 + 1];
            sg += beta * (1.0 - cos(tp));
            sin_p[s0] = sin(tp);
        }
    }

    if (force_out != NULL) {
        for (int t = 0; t < l; ++t) {
            for (int x = 0; x < l; ++x) {
                int s0  = t * l + x;
                int sxm = t * l + pbc(x - 1, l);
                int stm = pbc(t - 1, l) * l + x;
                force_out[2 * s0 + 0] = -beta * (sin_p[s0] - sin_p[sxm]);
                force_out[2 * s0 + 1] = -beta * (sin_p[stm] - sin_p[s0]);
            }
        }
    }

    double sf = 0.0;
    if (dynamical || out_pbp != NULL) {
        cpx_t *d    = (cpx_t *)malloc((size_t)(v * v) * sizeof(cpx_t));
        cpx_t *dinv = (cpx_t *)malloc((size_t)(v * v) * sizeof(cpx_t));
        if (!d || !dinv) exit(1);
        build_staggered_2d(l, th, m, d);
        int need_inv = (force_out != NULL || out_pbp != NULL);
        double ldet  = cpx_lu_logdet_inv(v, d, need_inv ? dinv : NULL);
        if (dynamical) sf = -2.0 * ldet;

        if (out_pbp != NULL) {
            double tr_re = 0.0;
            for (int i = 0; i < v; ++i) tr_re += dinv[i * v + i].re;
            *out_pbp = tr_re / (double)v;
        }
        if (dynamical && force_out != NULL) {
            for (int t = 0; t < l; ++t) {
                for (int x = 0; x < l; ++x) {
                    int s0 = t * l + x;
                    for (int mu = 0; mu < 2; ++mu) {
                        int tn     = (mu == 0) ? pbc(t + 1, l) : t;
                        int xn     = (mu == 1) ? pbc(x + 1, l) : x;
                        int s1     = tn * l + xn;
                        double eta = (mu == 0) ? 1.0 : ((t & 1) ? -1.0 : 1.0);
                        double bc  = (mu == 0 && t == l - 1) ? -1.0 : 1.0;
                        cpx_t iu  = c_mul(c_mk(0.0, 1.0), c_scale(c_exp_i(th[2 * s0 + mu]), 0.5 * eta * bc));
                        cpx_t iuc = c_mk(-iu.re, iu.im); /* +i * u* */
                        cpx_t tr  = c_add(c_mul(dinv[s1 * v + s0], iu), c_mul(dinv[s0 * v + s1], iuc));
                        force_out[2 * s0 + mu] += 2.0 * tr.re;
                    }
                }
            }
        }
        free(d);
        free(dinv);
    }
    return sg + sf;
}

/* --- Self-Test 1 helper: 64-dim Jordan-Wigner operators for n=6 modes --- */
static void jw_apply_annih(int k, const double *in, double *out) {
    int mask = (1 << k) - 1;
    int bit  = 1 << k;
    for (int s = 0; s < JW_DIM; ++s) {
        if ((s & bit) != 0) {
            out[s] = 0.0;
        } else {
            int par = __builtin_popcount((unsigned)(s & mask)) & 1;
            double sg = par ? -1.0 : 1.0;
            out[s]    = sg * in[s ^ bit];
        }
    }
}

static void jw_apply_create(int k, const double *in, double *out) {
    int mask = (1 << k) - 1;
    int bit  = 1 << k;
    for (int s = 0; s < JW_DIM; ++s) {
        if ((s & bit) == 0) {
            out[s] = 0.0;
        } else {
            int par = __builtin_popcount((unsigned)(s & mask)) & 1;
            double sg = par ? -1.0 : 1.0;
            out[s]    = sg * in[s ^ bit];
        }
    }
}

typedef struct {
    double max_jw_err;
    double max_clifford_err;
    double max_g5_herm_err;
    double max_gauge_cov_err;
    int nn_naive_zeros;
    int nn_wilson_zeros;
    double nn_wilson_doubler_min;
    int stag_pi_violations;
    double ov_gw_err;
    double ov_norm_err;
    double ov_circle_err;
    double hmc_force_err;
    double lefschetz_st_im_err;
} selftest_res_t;

static void run_startup_selftests(selftest_res_t *out) {
    /* 1. Jordan-Wigner CAR, Pauli exclusion, number idempotence across all 64 basis vectors (n=6) */
    double max_jw = 0.0;
    double e_s[JW_DIM], t1[JW_DIM], t2[JW_DIM], u1[JW_DIM], u2[JW_DIM];
    for (int b = 0; b < JW_DIM; ++b) {
        memset(e_s, 0, sizeof(e_s));
        e_s[b] = 1.0;
        for (int j = 0; j < JW_MODES; ++j) {
            for (int k = 0; k < JW_MODES; ++k) {
                jw_apply_annih(k, e_s, t1); jw_apply_annih(j, t1, t2);
                jw_apply_annih(j, e_s, u1); jw_apply_annih(k, u1, u2);
                for (int s = 0; s < JW_DIM; ++s) {
                    double v = fabs(t2[s] + u2[s]);
                    if (v > max_jw) max_jw = v;
                }
                jw_apply_create(k, e_s, t1); jw_apply_create(j, t1, t2);
                jw_apply_create(j, e_s, u1); jw_apply_create(k, u1, u2);
                for (int s = 0; s < JW_DIM; ++s) {
                    double v = fabs(t2[s] + u2[s]);
                    if (v > max_jw) max_jw = v;
                }
                jw_apply_create(k, e_s, t1); jw_apply_annih(j, t1, t2);
                jw_apply_annih(j, e_s, u1);  jw_apply_create(k, u1, u2);
                for (int s = 0; s < JW_DIM; ++s) {
                    double dlt = (j == k) ? e_s[s] : 0.0;
                    double v   = fabs(t2[s] + u2[s] - dlt);
                    if (v > max_jw) max_jw = v;
                }
            }
            jw_apply_annih(j, e_s, t1); jw_apply_create(j, t1, t2);
            jw_apply_annih(j, t2, u1);  jw_apply_create(j, u1, u2);
            for (int s = 0; s < JW_DIM; ++s) {
                double v = fabs(u2[s] - t2[s]);
                if (v > max_jw) max_jw = v;
            }
        }
    }

    /* 2. 2D/4D Clifford algebra & 2D Wilson gamma_5-Hermiticity + local U(1) gauge covariance */
    cpx_t g2[2][2][2], g5_2[2][2];
    cpx_t g4[4][4][4], g5_4[4][4];
    init_gamma2d(g2, g5_2);
    init_gamma4d(g4, g5_4);
    double max_cl = 0.0;
    for (int mu = 0; mu < 2; ++mu) {
        for (int nu = 0; nu < 2; ++nu) {
            for (int i = 0; i < 2; ++i) {
                double d_dag = sqrt(c_abs2(c_sub(c_conj(g2[mu][i][i]), g2[mu][i][i])));
                if (d_dag > max_cl) max_cl = d_dag;
                for (int j = 0; j < 2; ++j) {
                    cpx_t ac = c_mk(0.0, 0.0);
                    for (int k = 0; k < 2; ++k) {
                        ac = c_add(ac, c_add(c_mul(g2[mu][i][k], g2[nu][k][j]),
                                             c_mul(g2[nu][i][k], g2[mu][k][j])));
                    }
                    double tgt = (mu == nu && i == j) ? 2.0 : 0.0;
                    double err = sqrt(c_abs2(c_sub(ac, c_mk(tgt, 0.0))));
                    if (err > max_cl) max_cl = err;
                }
            }
        }
    }
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    cpx_t ac = c_mk(0.0, 0.0);
                    for (int k = 0; k < 4; ++k) {
                        ac = c_add(ac, c_add(c_mul(g4[mu][i][k], g4[nu][k][j]),
                                             c_mul(g4[nu][i][k], g4[mu][k][j])));
                    }
                    double tgt = (mu == nu && i == j) ? 2.0 : 0.0;
                    double err = sqrt(c_abs2(c_sub(ac, c_mk(tgt, 0.0))));
                    if (err > max_cl) max_cl = err;
                }
            }
        }
    }
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            cpx_t prod = c_mk(0.0, 0.0);
            for (int k1 = 0; k1 < 4; ++k1) {
                for (int k2 = 0; k2 < 4; ++k2) {
                    for (int k3 = 0; k3 < 4; ++k3) {
                        cpx_t t = c_mul(g4[0][i][k1], g4[1][k1][k2]);
                        t       = c_mul(t, g4[2][k2][k3]);
                        t       = c_mul(t, g4[3][k3][j]);
                        prod    = c_add(prod, t);
                    }
                }
            }
            double err = sqrt(c_abs2(c_sub(prod, g5_4[i][j])));
            if (err > max_cl) max_cl = err;
        }
    }

    const int l_st = L_SCHW;
    const int v_st = l_st * l_st;
    const int d_st = 2 * v_st;
    rng_t rng;
    rng_seed(&rng, 0x57A6E701ULL);
    double th[2 * L_SCHW * L_SCHW], th_g[2 * L_SCHW * L_SCHW], alpha[L_SCHW * L_SCHW];
    for (int i = 0; i < 2 * v_st; ++i) th[i] = (2.0 * rng_u01(&rng) - 1.0) * K_PI;
    for (int s = 0; s < v_st; ++s) alpha[s] = (2.0 * rng_u01(&rng) - 1.0) * K_PI;
    for (int t = 0; t < l_st; ++t) {
        for (int x = 0; x < l_st; ++x) {
            int s0 = t * l_st + x;
            int st = pbc(t + 1, l_st) * l_st + x;
            int sx = t * l_st + pbc(x + 1, l_st);
            th_g[2 * s0 + 0] = th[2 * s0 + 0] + alpha[s0] - alpha[st];
            th_g[2 * s0 + 1] = th[2 * s0 + 1] + alpha[s0] - alpha[sx];
        }
    }

    cpx_t *dw  = (cpx_t *)malloc((size_t)(d_st * d_st) * sizeof(cpx_t));
    cpx_t *dwg = (cpx_t *)malloc((size_t)(d_st * d_st) * sizeof(cpx_t));
    cpx_t *dov = (cpx_t *)malloc((size_t)(d_st * d_st) * sizeof(cpx_t));
    if (!dw || !dwg || !dov) exit(1);
    build_wilson_2d(l_st, th, -1.0, 1.0, 1, dw);
    build_wilson_2d(l_st, th_g, -1.0, 1.0, 1, dwg);

    double max_g5h = 0.0, max_gcov = 0.0;
    for (int i = 0; i < d_st; ++i) {
        double gi = ((i & 1) == 0) ? 1.0 : -1.0;
        cpx_t g_i = c_exp_i(alpha[i / 2]);
        for (int j = 0; j < d_st; ++j) {
            double gj  = ((j & 1) == 0) ? 1.0 : -1.0;
            cpx_t g_j  = c_exp_i(alpha[j / 2]);
            cpx_t lhs  = c_scale(dw[i * d_st + j], gi * gj);
            cpx_t rhs  = c_conj(dw[j * d_st + i]);
            double e5  = sqrt(c_abs2(c_sub(lhs, rhs)));
            if (e5 > max_g5h) max_g5h = e5;

            cpx_t cov  = c_mul(g_i, c_mul(dw[i * d_st + j], c_conj(g_j)));
            double ec  = sqrt(c_abs2(c_sub(dwg[i * d_st + j], cov)));
            if (ec > max_gcov) max_gcov = ec;
        }
    }

    /* 3. Nielsen-Ninomiya 8x8 doubling & Wilson lifting + 4^4 staggered pi-flux */
    const int l_nn = L_IDX;
    const int v_nn = l_nn * l_nn;
    const int d_nn = 2 * v_nn;
    double th_zero[2 * L_IDX * L_IDX];
    memset(th_zero, 0, sizeof(th_zero));
    cpx_t *d_naive = (cpx_t *)malloc((size_t)(d_nn * d_nn) * sizeof(cpx_t));
    cpx_t *d_wil   = (cpx_t *)malloc((size_t)(d_nn * d_nn) * sizeof(cpx_t));
    if (!d_naive || !d_wil) exit(1);
    build_wilson_2d(l_nn, th_zero, 0.0, 0.0, 0, d_naive);
    build_wilson_2d(l_nn, th_zero, 0.0, 1.0, 0, d_wil);
    int nn_z = 0, nw_z = 0;
    double nw_doub_min = 1e9;
    for (int kt = 0; kt < 2; ++kt) {
        for (int kx = 0; kx < 2; ++kx) {
            for (int sp = 0; sp < 2; ++sp) {
                double n_naive = 0.0, n_wil = 0.0;
                for (int t = 0; t < l_nn; ++t) {
                    for (int x = 0; x < l_nn; ++x) {
                        int s0     = t * l_nn + x;
                        double sg0 = ((kt ? (t & 1) : 0) ^ (kx ? (x & 1) : 0)) ? -1.0 : 1.0;
                        (void)sg0;
                        for (int a = 0; a < 2; ++a) {
                            int r     = 2 * s0 + a;
                            cpx_t yna = c_mk(0.0, 0.0), yw = c_mk(0.0, 0.0);
                            for (int t2 = 0; t2 < l_nn; ++t2) {
                                for (int x2 = 0; x2 < l_nn; ++x2) {
                                    double ph2 = ((kt ? (t2 & 1) : 0) ^ (kx ? (x2 & 1) : 0)) ? -1.0 : 1.0;
                                    int c2     = 2 * (t2 * l_nn + x2) + sp;
                                    yna = c_add(yna, c_scale(d_naive[r * d_nn + c2], ph2 / (double)l_nn));
                                    yw  = c_add(yw, c_scale(d_wil[r * d_nn + c2], ph2 / (double)l_nn));
                                }
                            }
                            n_naive += c_abs2(yna);
                            n_wil   += c_abs2(yw);
                        }
                    }
                }
                if (sqrt(n_naive) < 1e-12) ++nn_z;
                if (sqrt(n_wil) < 1e-12) {
                    ++nw_z;
                } else if (kt != 0 || kx != 0) {
                    double nm = sqrt(n_wil);
                    if (nm < nw_doub_min) nw_doub_min = nm;
                }
            }
        }
    }
    nn_z /= 2; /* Dirac modes (2 spinors each) */
    nw_z /= 2;
    free(d_naive);
    free(d_wil);

    int pi_viol = 0;
    for (int t = 0; t < 4; ++t) {
        for (int x = 0; x < 4; ++x) {
            for (int y = 0; y < 4; ++y) {
                for (int z = 0; z < 4; ++z) {
                    for (int mu = 0; mu < 4; ++mu) {
                        for (int nu = 0; nu < 4; ++nu) {
                            if (mu == nu) continue;
                            int c[4] = {t, x, y, z}, cm[4] = {t, x, y, z}, cn[4] = {t, x, y, z};
                            cm[mu]   = (cm[mu] + 1) & 3;
                            cn[nu]   = (cn[nu] + 1) & 3;
                            int eta_p_mu = (mu == 0) ? 1 : (mu == 1) ? ((c[0] & 1) ? -1 : 1) :
                                           (mu == 2) ? (((c[0] + c[1]) & 1) ? -1 : 1) :
                                                       (((c[0] + c[1] + c[2]) & 1) ? -1 : 1);
                            int eta_m_nu = (nu == 0) ? 1 : (nu == 1) ? ((cm[0] & 1) ? -1 : 1) :
                                           (nu == 2) ? (((cm[0] + cm[1]) & 1) ? -1 : 1) :
                                                       (((cm[0] + cm[1] + cm[2]) & 1) ? -1 : 1);
                            int eta_n_mu = (mu == 0) ? 1 : (mu == 1) ? ((cn[0] & 1) ? -1 : 1) :
                                           (mu == 2) ? (((cn[0] + cn[1]) & 1) ? -1 : 1) :
                                                       (((cn[0] + cn[1] + cn[2]) & 1) ? -1 : 1);
                            int eta_p_nu = (nu == 0) ? 1 : (nu == 1) ? ((c[0] & 1) ? -1 : 1) :
                                           (nu == 2) ? (((c[0] + c[1]) & 1) ? -1 : 1) :
                                                       (((c[0] + c[1] + c[2]) & 1) ? -1 : 1);
                            if (eta_p_mu * eta_m_nu * eta_n_mu * eta_p_nu != -1) ++pi_viol;
                        }
                    }
                }
            }
        }
    }

    /* 4. 6x6 Overlap Ginsparg-Wilson exactness, normality & circle spectrum */
    construct_overlap_2d(l_st, th, -1.0, dov, NULL);
    double gw_err = 0.0, norm_err = 0.0, circ_err = 0.0;
    for (int i = 0; i < d_st; ++i) {
        double gi = ((i & 1) == 0) ? 1.0 : -1.0;
        for (int j = 0; j < d_st; ++j) {
            double gj  = ((j & 1) == 0) ? 1.0 : -1.0;
            cpx_t g5d  = c_scale(dov[i * d_st + j], gi);
            cpx_t dg5  = c_scale(dov[i * d_st + j], gj);
            cpx_t dg5d = c_mk(0.0, 0.0);
            cpx_t dd_1 = c_mk(0.0, 0.0);
            cpx_t dd_2 = c_mk(0.0, 0.0);
            cpx_t m_ci = c_mk(0.0, 0.0);
            for (int k = 0; k < d_st; ++k) {
                double gk = ((k & 1) == 0) ? 1.0 : -1.0;
                dg5d = c_add(dg5d, c_mul(c_scale(dov[i * d_st + k], gk), dov[k * d_st + j]));
                dd_1 = c_add(dd_1, c_mul(c_conj(dov[k * d_st + i]), dov[k * d_st + j]));
                dd_2 = c_add(dd_2, c_mul(dov[i * d_st + k], c_conj(dov[j * d_st + k])));
                cpx_t a_ki = (k == i) ? c_sub(dov[k * d_st + i], c_mk(1.0, 0.0)) : dov[k * d_st + i];
                cpx_t a_kj = (k == j) ? c_sub(dov[k * d_st + j], c_mk(1.0, 0.0)) : dov[k * d_st + j];
                m_ci       = c_add(m_ci, c_mul(c_conj(a_ki), a_kj));
            }
            double egw = sqrt(c_abs2(c_sub(c_add(g5d, dg5), dg5d)));
            double en  = sqrt(c_abs2(c_sub(dd_1, dd_2)));
            double ec  = sqrt(c_abs2(c_sub(m_ci, (i == j) ? c_mk(1.0, 0.0) : c_mk(0.0, 0.0))));
            if (egw > gw_err) gw_err = egw;
            if (en > norm_err) norm_err = en;
            if (ec > circ_err) circ_err = ec;
        }
    }
    free(dw);
    free(dwg);
    free(dov);

    /* 5. Exact staggered HMC fermion force vs central FD on 6x6 */
    double f_ana[2 * L_SCHW * L_SCHW], th_p[2 * L_SCHW * L_SCHW], th_m[2 * L_SCHW * L_SCHW];
    (void)eval_schwinger_action_force(l_st, th, 2.0, 0.15, 1, f_ana, NULL);
    const double eps_fd = 1e-6;
    double max_ferr     = 0.0;
    for (int l = 0; l < 2 * v_st; ++l) {
        memcpy(th_p, th, (size_t)(2 * v_st) * sizeof(double));
        memcpy(th_m, th, (size_t)(2 * v_st) * sizeof(double));
        th_p[l] += eps_fd;
        th_m[l] -= eps_fd;
        double sp = eval_schwinger_action_force(l_st, th_p, 2.0, 0.15, 1, NULL, NULL);
        double sm = eval_schwinger_action_force(l_st, th_m, 2.0, 0.15, 1, NULL, NULL);
        double f_fd = -(sp - sm) / (2.0 * eps_fd);
        double d    = fabs(f_ana[l] - f_fd);
        if (d > max_ferr) max_ferr = d;
    }

    /* 6. Lefschetz quadratic saddle imaginary-action cancellation */
    double l_im_err = 0.0;
    for (int k = -20; k <= 20; ++k) {
        double a = 2.0, mu = 3.5, x = 0.25 * (double)k, y = mu / a;
        cpx_t z = c_mk(x, y);
        cpx_t s = c_sub(c_scale(c_mul(z, z), 0.5 * a), c_mul(c_mk(0.0, mu), z));
        if (fabs(s.im) > l_im_err) l_im_err = fabs(s.im);
    }

    if (max_jw > 1e-14 || max_cl > 1e-14 || max_g5h > 1e-13 || max_gcov > 1e-13 ||
        nn_z != 4 || nw_z != 1 || nw_doub_min < 1.9 || pi_viol != 0 ||
        gw_err > 1e-11 || norm_err > 1e-11 || circ_err > 1e-11 ||
        max_ferr > 1e-6 || l_im_err > 1e-14) {
        fprintf(stderr, "SELFTEST FAIL: jw=%.2e cl=%.2e g5h=%.2e cov=%.2e nn=%d->%d(%.2f) pi=%d "
                "gw=%.2e norm=%.2e circ=%.2e f=%.2e lef=%.2e\n",
                max_jw, max_cl, max_g5h, max_gcov, nn_z, nw_z, nw_doub_min, pi_viol,
                gw_err, norm_err, circ_err, max_ferr, l_im_err);
        exit(1);
    }

    out->max_jw_err           = max_jw;
    out->max_clifford_err     = max_cl;
    out->max_g5_herm_err      = max_g5h;
    out->max_gauge_cov_err    = max_gcov;
    out->nn_naive_zeros       = nn_z;
    out->nn_wilson_zeros      = nw_z;
    out->nn_wilson_doubler_min = nw_doub_min;
    out->stag_pi_violations   = pi_viol;
    out->ov_gw_err            = gw_err;
    out->ov_norm_err          = norm_err;
    out->ov_circle_err        = circ_err;
    out->hmc_force_err        = max_ferr;
    out->lefschetz_st_im_err  = l_im_err;

    printf("  [SELFTEST] PASS (JW64=%.2e Cliff=%.2e g5H=%.2e U1cov=%.2e NN=%d->%d piFlux=%d "
           "GW=%.2e HMC_F=%.2e LefIm=%.2e)\n",
           max_jw, max_cl, max_g5h, max_gcov, nn_z, nw_z, pi_viol, gw_err, max_ferr, l_im_err);
    fflush(stdout);
}

/* --- Experiment 1: 8x8 Discrete Atiyah-Singer Index Theorem --- */
typedef struct {
    int q_target;
    int noisy;
    long long q_geo;
    int n_plus;
    int n_minus;
    long long index_zero_mode;
    long long index_trace;
    double max_nonzero_chiral;
    double max_eig_err;
} idx_cfg_t;

static void run_exp1_index(idx_cfg_t res[N_IDX_CFGS]) {
    const int l = L_IDX, v = l * l, dim = 2 * v;
    const int q_list[N_IDX_SECTORS] = {-2, -1, 0, 1, 2};
    cpx_t *dov   = (cpx_t *)malloc((size_t)(dim * dim) * sizeof(cpx_t));
    cpx_t *hcomb = (cpx_t *)malloc((size_t)(dim * dim) * sizeof(cpx_t));
    cpx_t *vec   = (cpx_t *)malloc((size_t)(dim * dim) * sizeof(cpx_t));
    double *ev   = (double *)malloc((size_t)dim * sizeof(double));
    if (!dov || !hcomb || !vec || !ev) exit(1);

    rng_t rng;
    rng_seed(&rng, 0x1ADE4700ULL);
    double th[2 * L_IDX * L_IDX];

    for (int qi = 0; qi < N_IDX_SECTORS; ++qi) {
        int qt = q_list[qi];
        for (int ns = 0; ns < 2; ++ns) {
            int idx = 2 * qi + ns;
            memset(&res[idx], 0, sizeof(res[idx]));
            res[idx].q_target = qt;
            res[idx].noisy    = ns;
            build_harmonic_instanton_2d(l, (double)qt, th);
            if (ns) {
                for (int i = 0; i < 2 * v; ++i) th[i] += 0.12 * (2.0 * rng_u01(&rng) - 1.0);
            }
            res[idx].q_geo = u1_topo_charge_2d(l, th);

            double tr_hs = 0.0;
            construct_overlap_2d(l, th, -1.0, dov, &tr_hs);
            res[idx].index_trace = -llround(0.5 * tr_hs);

            /* Normal matrix D_ov = A + i B with Hermitian A = 0.5(D+D†), B = 1/(2i)(D-D†).
             * Simultaneous eigenvectors satisfy (A + phi_inv * B) v = ev * v. */
            const double alpha_ir = 0.6180339887498949;
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < dim; ++j) {
                    cpx_t d_ij = dov[i * dim + j];
                    cpx_t d_ji = c_conj(dov[j * dim + i]);
                    cpx_t a_ij = c_scale(c_add(d_ij, d_ji), 0.5);
                    cpx_t b_ij = c_mul(c_mk(0.0, -0.5), c_sub(d_ij, d_ji));
                    hcomb[i * dim + j] = c_add(a_ij, c_scale(b_ij, alpha_ir));
                }
            }
            jacobi_hermitian_eigen(dim, hcomb, ev, vec);

            double nplus_sum = 0.0, nminus_sum = 0.0;
            int nzero_tot = 0;
            double max_nz_chi = 0.0, max_eig_e = 0.0;
            for (int k = 0; k < dim; ++k) {
                /* Rayleigh quotient lambda_k = v_k† D_ov v_k and eigenvector residual */
                cpx_t lam = c_mk(0.0, 0.0);
                double g5_exp = 0.0;
                cpx_t dv_vec[2 * L_IDX * L_IDX];
                for (int i = 0; i < dim; ++i) {
                    cpx_t acc = c_mk(0.0, 0.0);
                    for (int j = 0; j < dim; ++j) {
                        acc = c_add(acc, c_mul(dov[i * dim + j], vec[j * dim + k]));
                    }
                    dv_vec[i] = acc;
                    lam       = c_add(lam, c_mul(c_conj(vec[i * dim + k]), acc));
                    double gi = ((i & 1) == 0) ? 1.0 : -1.0;
                    g5_exp   += gi * c_abs2(vec[i * dim + k]);
                }
                double res_e = 0.0;
                for (int i = 0; i < dim; ++i) {
                    res_e += c_abs2(c_sub(dv_vec[i], c_mul(lam, vec[i * dim + k])));
                }
                res_e = sqrt(res_e);
                if (res_e > max_eig_e) max_eig_e = res_e;

                double slam = sqrt(c_abs2(lam));
                if (slam < 1e-6) {
                    ++nzero_tot;
                    nplus_sum  += 0.5 * (1.0 + g5_exp);
                    nminus_sum += 0.5 * (1.0 - g5_exp);
                } else if (slam > 0.05 && slam < 1.95) {
                    if (fabs(g5_exp) > max_nz_chi) max_nz_chi = fabs(g5_exp);
                }
            }
            (void)nzero_tot;
            res[idx].n_plus             = (int)llround(nplus_sum);
            res[idx].n_minus            = (int)llround(nminus_sum);
            res[idx].index_zero_mode    = (long long)(res[idx].n_plus - res[idx].n_minus);
            res[idx].max_nonzero_chiral = max_nz_chi;
            res[idx].max_eig_err        = max_eig_e;
        }
    }
    free(dov);
    free(hcomb);
    free(vec);
    free(ev);
}

/* --- Experiment 2: 6x6 2D Schwinger Model Exact-Determinant HMC vs. Quenched --- */
typedef struct {
    double beta;
    double mass;
    int dynamical;
    uint64_t seed;
    double p_mean, p_err;
    double chi22, chi22_err;
    double pbp_mean, pbp_err;
    double acc_rate;
    double w11, w12, w21, w22;
} hmc_run_t;

static void *hmc_worker(void *arg) {
    hmc_run_t *r = (hmc_run_t *)arg;
    const int nl = 2 * V_SCHW;
    double th[2 * V_SCHW], th_tr[2 * V_SCHW], pi_mom[2 * V_SCHW], force[2 * V_SCHW];
    rng_t rng;
    rng_seed(&rng, r->seed);
    for (int i = 0; i < nl; ++i) th[i] = 0.0;

    double *p_ts   = (double *)malloc((size_t)N_MEAS_HMC * sizeof(double));
    double *pbp_ts = (double *)malloc((size_t)N_MEAS_HMC * sizeof(double));
    double *w11_ts = (double *)malloc((size_t)N_MEAS_HMC * sizeof(double));
    double *w12_ts = (double *)malloc((size_t)N_MEAS_HMC * sizeof(double));
    double *w21_ts = (double *)malloc((size_t)N_MEAS_HMC * sizeof(double));
    double *w22_ts = (double *)malloc((size_t)N_MEAS_HMC * sizeof(double));
    if (!p_ts || !pbp_ts || !w11_ts || !w12_ts || !w21_ts || !w22_ts) exit(1);

    const double dt = 0.06;
    int total_traj  = N_THERM_HMC + N_MEAS_HMC;
    int acc_meas    = 0;

    for (int tr = 0; tr < total_traj; ++tr) {
        memcpy(th_tr, th, (size_t)nl * sizeof(double));
        double k0 = 0.0;
        for (int i = 0; i < nl; ++i) {
            pi_mom[i] = rng_gauss(&rng);
            k0 += 0.5 * pi_mom[i] * pi_mom[i];
        }
        double s0 = eval_schwinger_action_force(L_SCHW, th, r->beta, r->mass, r->dynamical, force, NULL);
        double h0 = s0 + k0;

        for (int i = 0; i < nl; ++i) pi_mom[i] += 0.5 * dt * force[i];
        for (int st = 0; st < HMC_STEPS; ++st) {
            for (int i = 0; i < nl; ++i) th_tr[i] += dt * pi_mom[i];
            (void)eval_schwinger_action_force(L_SCHW, th_tr, r->beta, r->mass, r->dynamical, force, NULL);
            double mul = (st + 1 < HMC_STEPS) ? dt : (0.5 * dt);
            for (int i = 0; i < nl; ++i) pi_mom[i] += mul * force[i];
        }
        double k1 = 0.0;
        for (int i = 0; i < nl; ++i) k1 += 0.5 * pi_mom[i] * pi_mom[i];
        double s1 = eval_schwinger_action_force(L_SCHW, th_tr, r->beta, r->mass, r->dynamical, NULL, NULL);
        double dH = (s1 + k1) - h0;
        if (dH <= 0.0 || rng_pos_u01(&rng) < exp(-dH)) {
            memcpy(th, th_tr, (size_t)nl * sizeof(double));
            if (tr >= N_THERM_HMC) ++acc_meas;
        }
        for (int i = 0; i < nl; ++i) th[i] = wrap_angle(th[i]);

        if (tr >= N_THERM_HMC) {
            int m = tr - N_THERM_HMC;
            p_ts[m] = u1_mean_plaq_2d(L_SCHW, th);
            double pbp = 0.0;
            (void)eval_schwinger_action_force(L_SCHW, th, r->beta, r->mass, 0, NULL, &pbp);
            pbp_ts[m] = pbp;
            u1_wilson_loops_2d(L_SCHW, th, &w11_ts[m], &w12_ts[m], &w21_ts[m], &w22_ts[m]);
        }
    }

    double tp = 0.0, tpbp = 0.0;
    int fp = 0, fpbp = 0;
    madras_sokal(p_ts, N_MEAS_HMC, &r->p_mean, &tp, &r->p_err, &fp);
    madras_sokal(pbp_ts, N_MEAS_HMC, &r->pbp_mean, &tpbp, &r->pbp_err, &fpbp);
    r->acc_rate = (double)acc_meas / (double)N_MEAS_HMC;

    double s11 = 0.0, s12 = 0.0, s21 = 0.0, s22 = 0.0;
    for (int i = 0; i < N_MEAS_HMC; ++i) {
        s11 += w11_ts[i]; s12 += w12_ts[i]; s21 += w21_ts[i]; s22 += w22_ts[i];
    }
    r->w11 = s11 / (double)N_MEAS_HMC;
    r->w12 = s12 / (double)N_MEAS_HMC;
    r->w21 = s21 / (double)N_MEAS_HMC;
    r->w22 = s22 / (double)N_MEAS_HMC;

    double rat_full = (r->w22 * r->w11) / (r->w12 * r->w21);
    r->chi22        = (rat_full > 1e-12) ? (-log(rat_full)) : 0.0;

    int bsz = N_MEAS_HMC / JK_BLOCKS;
    double jk_c[JK_BLOCKS];
    double jk_m = 0.0;
    for (int b = 0; b < JK_BLOCKS; ++b) {
        double b11 = 0.0, b12 = 0.0, b21 = 0.0, b22 = 0.0;
        int lo = b * bsz, hi = (b == JK_BLOCKS - 1) ? N_MEAS_HMC : ((b + 1) * bsz);
        for (int i = lo; i < hi; ++i) {
            b11 += w11_ts[i]; b12 += w12_ts[i]; b21 += w21_ts[i]; b22 += w22_ts[i];
        }
        double rat = ((s22 - b22) * (s11 - b11)) / ((s12 - b12) * (s21 - b21));
        jk_c[b]    = (rat > 1e-12) ? (-log(rat)) : 0.0;
        jk_m += jk_c[b];
    }
    jk_m /= (double)JK_BLOCKS;
    double jk_v = 0.0;
    for (int b = 0; b < JK_BLOCKS; ++b) {
        double d = jk_c[b] - jk_m;
        jk_v += d * d;
    }
    r->chi22_err = sqrt(jk_v * (double)(JK_BLOCKS - 1) / (double)JK_BLOCKS);

    free(p_ts);
    free(pbp_ts);
    free(w11_ts);
    free(w12_ts);
    free(w21_ts);
    free(w22_ts);
    return NULL;
}

/* --- Experiment 3A: 1D Finite-Density U(1) Sign Problem vs. Exact Positive Dual --- */
typedef struct {
    int L;
    double beta, mu;
    double pq_sgn, pq_sgn_err;
    double n_rw, n_rw_err;
    double dual_sgn;
    double n_dual, n_dual_err;
} exp3a_res_t;

static void *exp3a_worker(void *arg) {
    exp3a_res_t *res = (exp3a_res_t *)arg;
    const int L      = res->L;
    const double b   = res->beta;
    const double mu  = res->mu;
    const double b_ch = b * cosh(mu);
    const double b_sh = b * sinh(mu);
    rng_t rng;
    rng_seed(&rng, 0x3A000000ULL + (uint64_t)L * 919ULL);

    /* 1. Positive Worldline Dual: p(k) propto I_k(b)^L * exp(k * mu * L) > 0 */
    enum { K_MAX = 28, K_CNT = 2 * K_MAX + 1 };
    double logw[K_CNT], cdf[K_CNT];
    double max_lw = -1e300;
    for (int k = -K_MAX; k <= K_MAX; ++k) {
        double iv = bessel_in(k, b);
        double lw = (iv > 1e-300) ? ((double)L * log(iv) + (double)k * mu * (double)L) : -1e300;
        logw[k + K_MAX] = lw;
        if (lw > max_lw) max_lw = lw;
    }
    double z_sum = 0.0;
    for (int i = 0; i < K_CNT; ++i) {
        z_sum += exp(logw[i] - max_lw);
        cdf[i] = z_sum;
    }
    for (int i = 0; i < K_CNT; ++i) cdf[i] /= z_sum;

    double sum_k = 0.0, sum_k2 = 0.0;
    for (int s = 0; s < N_SAMP_3A; ++s) {
        double u = rng_u01(&rng);
        int idx  = 0;
        while (idx + 1 < K_CNT && cdf[idx] < u) ++idx;
        double kv = (double)(idx - K_MAX);
        sum_k  += kv;
        sum_k2 += kv * kv;
    }
    double mk       = sum_k / (double)N_SAMP_3A;
    double vk       = (sum_k2 / (double)N_SAMP_3A) - mk * mk;
    res->dual_sgn   = 1.0;
    res->n_dual     = mk;
    res->n_dual_err = sqrt(((vk > 0.0) ? vk : 0.0) / (double)N_SAMP_3A);

    /* 2. Naive Phase-Quenched Metropolis on Closed Ring theta[0..L-1] */
    enum { N_PQ = 60000 };
    double th[64];
    for (int i = 0; i < L; ++i) th[i] = 0.0;
    for (int sw = 0; sw < 1000; ++sw) {
        for (int t = 0; t < L; ++t) {
            int tp = (t + 1 < L) ? (t + 1) : 0;
            int tm = (t > 0) ? (t - 1) : (L - 1);
            double cand = wrap_angle(th[t] + 1.5 * (2.0 * rng_u01(&rng) - 1.0));
            double ds = b_ch * (cos(th[tp] - th[t]) + cos(th[t] - th[tm]) -
                                cos(th[tp] - cand) - cos(cand - th[tm]));
            if (ds <= 0.0 || rng_pos_u01(&rng) < exp(-ds)) th[t] = cand;
        }
    }

    double blk_sgn[JK_BLOCKS], blk_num[JK_BLOCKS];
    memset(blk_sgn, 0, sizeof(blk_sgn));
    memset(blk_num, 0, sizeof(blk_num));
    int bsz = N_PQ / JK_BLOCKS;
    for (int s = 0; s < N_PQ; ++s) {
        for (int sw = 0; sw < 2; ++sw) {
            for (int t = 0; t < L; ++t) {
                int tp = (t + 1 < L) ? (t + 1) : 0;
                int tm = (t > 0) ? (t - 1) : (L - 1);
                double cand = wrap_angle(th[t] + 1.5 * (2.0 * rng_u01(&rng) - 1.0));
                double ds = b_ch * (cos(th[tp] - th[t]) + cos(th[t] - th[tm]) -
                                    cos(th[tp] - cand) - cos(cand - th[tm]));
                if (ds <= 0.0 || rng_pos_u01(&rng) < exp(-ds)) th[t] = cand;
            }
        }
        double phi = 0.0, nr = 0.0, ni = 0.0;
        for (int t = 0; t < L; ++t) {
            int tp    = (t + 1 < L) ? (t + 1) : 0;
            double dth = th[tp] - th[t];
            phi += b_sh * sin(dth);
            nr  += b * cos(dth) * sinh(mu);
            ni  += b * sin(dth) * cosh(mu);
        }
        nr /= (double)L;
        ni /= (double)L;
        double cp = cos(phi), sp = sin(phi);
        int bk = s / bsz;
        if (bk >= JK_BLOCKS) bk = JK_BLOCKS - 1;
        blk_sgn[bk] += cp;
        blk_num[bk] += nr * cp - ni * sp;
    }

    double tot_sgn = 0.0, tot_num = 0.0;
    for (int b_i = 0; b_i < JK_BLOCKS; ++b_i) {
        tot_sgn += blk_sgn[b_i];
        tot_num += blk_num[b_i];
    }
    res->pq_sgn = tot_sgn / (double)N_PQ;
    res->n_rw   = tot_num / tot_sgn;

    double j_s[JK_BLOCKS], j_n[JK_BLOCKS], ms = 0.0, mn = 0.0;
    for (int b_i = 0; b_i < JK_BLOCKS; ++b_i) {
        j_s[b_i] = (tot_sgn - blk_sgn[b_i]) / (double)(N_PQ - bsz);
        j_n[b_i] = (tot_num - blk_num[b_i]) / (tot_sgn - blk_sgn[b_i]);
        ms += j_s[b_i];
        mn += j_n[b_i];
    }
    ms /= (double)JK_BLOCKS;
    mn /= (double)JK_BLOCKS;
    double vs = 0.0, vn = 0.0;
    for (int b_i = 0; b_i < JK_BLOCKS; ++b_i) {
        double ds_v = j_s[b_i] - ms;
        double dn_v = j_n[b_i] - mn;
        vs += ds_v * ds_v;
        vn += dn_v * dn_v;
    }
    res->pq_sgn_err = sqrt(vs * (double)(JK_BLOCKS - 1) / (double)JK_BLOCKS);
    res->n_rw_err   = sqrt(vn * (double)(JK_BLOCKS - 1) / (double)JK_BLOCKS);
    return NULL;
}

/* --- Experiment 3B: Non-Gaussian Complex Action & Lefschetz Thimble Contour --- */
typedef struct {
    double mu, y0;
    double naive_sgn, naive_re_z2, naive_re_z2_err;
    double th_sgn, th_re_z2, th_re_z2_err, th_im_z, th_im_z_err;
} exp3b_res_t;

static void run_exp3b_thimble(exp3b_res_t out[N_MU_THIM]) {
    const double a = 2.0, lam = 0.05;
    const double mus[N_MU_THIM] = {0.5, 1.5, 2.5, 3.5, 4.0};
    rng_t rng;
    rng_seed(&rng, 0x3B000000ULL);
    int bsz = N_SAMP_3B / JK_BLOCKS;

    for (int mi = 0; mi < N_MU_THIM; ++mi) {
        double mu = mus[mi];
        memset(&out[mi], 0, sizeof(out[mi]));
        out[mi].mu = mu;

        double lo = 0.0, hi = sqrt(a / (3.0 * lam));
        for (int it = 0; it < 60; ++it) {
            double mid = 0.5 * (lo + hi);
            double fv  = a * mid - lam * mid * mid * mid;
            if (fv < mu) lo = mid; else hi = mid;
        }
        double y0 = 0.5 * (lo + hi);
        out[mi].y0 = y0;

        /* Naive real-axis importance sampling x ~ N(0, 1/a) */
        double sig_nv = 1.0 / sqrt(a);
        double bn_w[JK_BLOCKS], bn_re[JK_BLOCKS], bn_im[JK_BLOCKS], bn_o[JK_BLOCKS];
        memset(bn_w, 0, sizeof(bn_w)); memset(bn_re, 0, sizeof(bn_re));
        memset(bn_im, 0, sizeof(bn_im)); memset(bn_o, 0, sizeof(bn_o));
        for (int s = 0; s < N_SAMP_3B; ++s) {
            double x  = sig_nv * rng_gauss(&rng);
            double w  = exp(-0.25 * lam * x * x * x * x);
            double c  = cos(mu * x), sn = sin(mu * x);
            int bk    = s / bsz;
            if (bk >= JK_BLOCKS) bk = JK_BLOCKS - 1;
            bn_w[bk]  += w;
            bn_re[bk] += w * c;
            bn_im[bk] += w * sn;
            bn_o[bk]  += w * c * (x * x);
        }
        double tw = 0.0, tre = 0.0, tim = 0.0, to = 0.0;
        for (int b_i = 0; b_i < JK_BLOCKS; ++b_i) {
            tw += bn_w[b_i]; tre += bn_re[b_i]; tim += bn_im[b_i]; to += bn_o[b_i];
        }
        out[mi].naive_sgn   = hypot(tre, tim) / tw;
        out[mi].naive_re_z2 = to / tre;
        double mn_z2 = 0.0, j_nz2[JK_BLOCKS];
        for (int b_i = 0; b_i < JK_BLOCKS; ++b_i) {
            j_nz2[b_i] = (to - bn_o[b_i]) / (tre - bn_re[b_i]);
            mn_z2 += j_nz2[b_i];
        }
        mn_z2 /= (double)JK_BLOCKS;
        double vn_z2 = 0.0;
        for (int b_i = 0; b_i < JK_BLOCKS; ++b_i) {
            double d = j_nz2[b_i] - mn_z2;
            vn_z2 += d * d;
        }
        out[mi].naive_re_z2_err = sqrt(vn_z2 * (double)(JK_BLOCKS - 1) / (double)JK_BLOCKS);

        /* Lefschetz shifted saddle contour z = x + i y0, x ~ N(0, 1/(a - 3*lam*y0^2)) */
        double a_eff  = a - 3.0 * lam * y0 * y0;
        double sig_th = 1.0 / sqrt(a_eff);
        double bt_w[JK_BLOCKS], bt_re[JK_BLOCKS], bt_im[JK_BLOCKS], bt_z2[JK_BLOCKS], bt_iz[JK_BLOCKS];
        memset(bt_w, 0, sizeof(bt_w)); memset(bt_re, 0, sizeof(bt_re));
        memset(bt_im, 0, sizeof(bt_im)); memset(bt_z2, 0, sizeof(bt_z2)); memset(bt_iz, 0, sizeof(bt_iz));

        for (int s = 0; s < N_SAMP_3B; ++s) {
            double x    = sig_th * rng_gauss(&rng);
            double w    = exp(-0.25 * lam * x * x * x * x);
            double phi  = -lam * y0 * x * x * x;
            double c    = cos(phi), sn = sin(phi);
            double z2r  = x * x - y0 * y0;
            double z2i  = 2.0 * x * y0;
            int bk      = s / bsz;
            if (bk >= JK_BLOCKS) bk = JK_BLOCKS - 1;
            bt_w[bk]   += w;
            bt_re[bk]  += w * c;
            bt_im[bk]  += w * sn;
            bt_z2[bk]  += w * (c * z2r - sn * z2i);
            bt_iz[bk]  += w * (c * y0 + sn * x);
        }
        double twt = 0.0, tret = 0.0, timt = 0.0, tz2 = 0.0, tiz = 0.0;
        for (int b_i = 0; b_i < JK_BLOCKS; ++b_i) {
            twt += bt_w[b_i]; tret += bt_re[b_i]; timt += bt_im[b_i];
            tz2 += bt_z2[b_i]; tiz += bt_iz[b_i];
        }
        out[mi].th_sgn   = hypot(tret, timt) / twt;
        out[mi].th_re_z2 = tz2 / tret;
        out[mi].th_im_z  = tiz / tret;

        double j_z2[JK_BLOCKS], j_iz[JK_BLOCKS], mz2 = 0.0, miz = 0.0;
        for (int b_i = 0; b_i < JK_BLOCKS; ++b_i) {
            j_z2[b_i] = (tz2 - bt_z2[b_i]) / (tret - bt_re[b_i]);
            j_iz[b_i] = (tiz - bt_iz[b_i]) / (tret - bt_re[b_i]);
            mz2 += j_z2[b_i];
            miz += j_iz[b_i];
        }
        mz2 /= (double)JK_BLOCKS;
        miz /= (double)JK_BLOCKS;
        double vz2 = 0.0, viz = 0.0;
        for (int b_i = 0; b_i < JK_BLOCKS; ++b_i) {
            double d1 = j_z2[b_i] - mz2;
            double d2 = j_iz[b_i] - miz;
            vz2 += d1 * d1;
            viz += d2 * d2;
        }
        out[mi].th_re_z2_err = sqrt(vz2 * (double)(JK_BLOCKS - 1) / (double)JK_BLOCKS);
        out[mi].th_im_z_err  = sqrt(viz * (double)(JK_BLOCKS - 1) / (double)JK_BLOCKS);
    }
}

int main(int argc, char **argv) {
    const char *json_path = "build/fermion_2d_bench.json";
    for (int i = 1; i + 1 < argc; ++i) {
        if (strcmp(argv[i], "--json") == 0) json_path = argv[i + 1];
    }

    printf("=== PHY Stage 7: 2D/4D Lattice Fermions, Atiyah-Singer Index & Sign Problem ===\n");
    selftest_res_t st;
    run_startup_selftests(&st);

    /* Experiment 1: 8x8 Discrete Atiyah-Singer Index Theorem */
    idx_cfg_t exp1[N_IDX_CFGS];
    run_exp1_index(exp1);
    printf("  [EXP1 8x8] Discrete Atiyah-Singer Index Theorem (Overlap D_ov, m0=-1.0):\n");
    for (int i = 0; i < N_IDX_CFGS; ++i) {
        idx_cfg_t *r = &exp1[i];
        printf("    %s Q_tgt=%+d | Q_geo=%+lld n_+=%d n_-=%d | idx(n+-n-)=%+lld tr_idx=%+lld | max_nz_chi=%.2e eig_err=%.2e\n",
               r->noisy ? "noisy " : "smooth", r->q_target, r->q_geo, r->n_plus, r->n_minus,
               r->index_zero_mode, r->index_trace, r->max_nonzero_chiral, r->max_eig_err);
    }
    fflush(stdout);

    /* Experiment 2: 6x6 Dynamical N_f=2 Schwinger Model vs Quenched + Mass Scan (10 parallel pthreads) */
    const double betas_hmc[N_BETAS_HMC] = {1.5, 2.5, 4.0};
    const double masses_hmc[N_MASS_HMC] = {0.05, 0.15, 0.30, 0.50};
    hmc_run_t runs_q[N_BETAS_HMC], runs_d[N_BETAS_HMC], runs_m[N_MASS_HMC];
    pthread_t th_q[N_BETAS_HMC], th_d[N_BETAS_HMC], th_m[N_MASS_HMC];

    for (int i = 0; i < N_BETAS_HMC; ++i) {
        memset(&runs_q[i], 0, sizeof(runs_q[i]));
        runs_q[i].beta = betas_hmc[i]; runs_q[i].mass = 0.15; runs_q[i].dynamical = 0;
        runs_q[i].seed = 0x20007000ULL + (uint64_t)i * 101ULL;
        pthread_create(&th_q[i], NULL, hmc_worker, &runs_q[i]);

        memset(&runs_d[i], 0, sizeof(runs_d[i]));
        runs_d[i].beta = betas_hmc[i]; runs_d[i].mass = 0.15; runs_d[i].dynamical = 1;
        runs_d[i].seed = 0x21007000ULL + (uint64_t)i * 101ULL;
        pthread_create(&th_d[i], NULL, hmc_worker, &runs_d[i]);
    }
    for (int i = 0; i < N_MASS_HMC; ++i) {
        memset(&runs_m[i], 0, sizeof(runs_m[i]));
        runs_m[i].beta = 2.0; runs_m[i].mass = masses_hmc[i]; runs_m[i].dynamical = 1;
        runs_m[i].seed = 0x22007000ULL + (uint64_t)i * 307ULL;
        pthread_create(&th_m[i], NULL, hmc_worker, &runs_m[i]);
    }

    /* Experiment 3A in 4 pthreads concurrently */
    const int l_list[N_L_CHAIN] = {8, 16, 32, 64};
    exp3a_res_t exp3a[N_L_CHAIN];
    pthread_t th_3a[N_L_CHAIN];
    for (int i = 0; i < N_L_CHAIN; ++i) {
        memset(&exp3a[i], 0, sizeof(exp3a[i]));
        exp3a[i].L = l_list[i]; exp3a[i].beta = 1.0; exp3a[i].mu = 0.6;
        pthread_create(&th_3a[i], NULL, exp3a_worker, &exp3a[i]);
    }

    exp3b_res_t exp3b[N_MU_THIM];
    run_exp3b_thimble(exp3b);

    for (int i = 0; i < N_BETAS_HMC; ++i) {
        pthread_join(th_q[i], NULL);
        pthread_join(th_d[i], NULL);
    }
    for (int i = 0; i < N_MASS_HMC; ++i) pthread_join(th_m[i], NULL);
    for (int i = 0; i < N_L_CHAIN; ++i) pthread_join(th_3a[i], NULL);

    printf("  [EXP2 6x6] Dynamical N_f=2 Staggered Schwinger HMC vs. Quenched (m=0.15):\n");
    for (int i = 0; i < N_BETAS_HMC; ++i) {
        printf("    beta=%.1f | Quenched: <P>=%.5f chi22=%.4f <pbp>=%.5f (acc=%.1f%%) | "
               "Dyn Nf=2: <P>=%.5f chi22=%.4f <pbp>=%.5f (acc=%.1f%%)\n",
               betas_hmc[i], runs_q[i].p_mean, runs_q[i].chi22, runs_q[i].pbp_mean, 100.0 * runs_q[i].acc_rate,
               runs_d[i].p_mean, runs_d[i].chi22, runs_d[i].pbp_mean, 100.0 * runs_d[i].acc_rate);
    }
    printf("    Chiral condensate @ beta=2.0: ");
    for (int i = 0; i < N_MASS_HMC; ++i) {
        printf("m=%.2f->%.5f%s", runs_m[i].mass, runs_m[i].pbp_mean, (i + 1 < N_MASS_HMC) ? " | " : "\n");
    }
    printf("  [EXP3A 1D U(1) Sign] beta=1.0 mu=0.6 (Naive PQ vs. Positive Worldline Dual):\n");
    for (int i = 0; i < N_L_CHAIN; ++i) {
        printf("    L=%2d | Naive: <sgn>=%.5f+-%.5f <n>=%+.4f+-%.4f | Dual: <sgn>=%.6f <n>=%+.5f+-%.5f\n",
               exp3a[i].L, exp3a[i].pq_sgn, exp3a[i].pq_sgn_err, exp3a[i].n_rw, exp3a[i].n_rw_err,
               exp3a[i].dual_sgn, exp3a[i].n_dual, exp3a[i].n_dual_err);
    }
    printf("  [EXP3B Lefschetz Thimble] a=2.0 lam=0.05:\n");
    for (int i = 0; i < N_MU_THIM; ++i) {
        printf("    mu=%.1f y0=%.4f | Naive <sgn>=%.5f <z^2>=%.4f+-%.4f | "
               "Thimble <sgn>=%.5f <z^2>=%.4f+-%.4f <Im z>=%.4f+-%.4f\n",
               exp3b[i].mu, exp3b[i].y0, exp3b[i].naive_sgn, exp3b[i].naive_re_z2, exp3b[i].naive_re_z2_err,
               exp3b[i].th_sgn, exp3b[i].th_re_z2, exp3b[i].th_re_z2_err, exp3b[i].th_im_z, exp3b[i].th_im_z_err);
    }
    fflush(stdout);

    /* JSON Serialization */
    FILE *fp = fopen(json_path, "w");
    if (!fp) { fprintf(stderr, "ERROR: Failed to open %s\n", json_path); return 1; }
    fprintf(fp, "{\n");
    fprintf(fp, "  \"selftests\": {\n");
    fprintf(fp, "    \"max_jw_err\": %.14e,\n", st.max_jw_err);
    fprintf(fp, "    \"max_clifford_err\": %.14e,\n", st.max_clifford_err);
    fprintf(fp, "    \"max_g5_herm_err\": %.14e,\n", st.max_g5_herm_err);
    fprintf(fp, "    \"max_gauge_cov_err\": %.14e,\n", st.max_gauge_cov_err);
    fprintf(fp, "    \"nn_naive_zeros\": %d,\n", st.nn_naive_zeros);
    fprintf(fp, "    \"nn_wilson_zeros\": %d,\n", st.nn_wilson_zeros);
    fprintf(fp, "    \"nn_wilson_doubler_min\": %.6f,\n", st.nn_wilson_doubler_min);
    fprintf(fp, "    \"stag_pi_violations\": %d,\n", st.stag_pi_violations);
    fprintf(fp, "    \"ov_gw_err\": %.14e,\n", st.ov_gw_err);
    fprintf(fp, "    \"ov_norm_err\": %.14e,\n", st.ov_norm_err);
    fprintf(fp, "    \"ov_circle_err\": %.14e,\n", st.ov_circle_err);
    fprintf(fp, "    \"hmc_force_err\": %.14e,\n", st.hmc_force_err);
    fprintf(fp, "    \"lefschetz_st_im_err\": %.14e\n", st.lefschetz_st_im_err);
    fprintf(fp, "  },\n");

    fprintf(fp, "  \"exp1_index\": [\n");
    for (int i = 0; i < N_IDX_CFGS; ++i) {
        idx_cfg_t *r = &exp1[i];
        fprintf(fp, "    {\"q_target\": %d, \"noisy\": %d, \"q_geo\": %lld, \"n_plus\": %d, \"n_minus\": %d, "
                    "\"index_zero_mode\": %lld, \"index_trace\": %lld, \"max_nonzero_chiral\": %.6e, \"max_eig_err\": %.6e}%s\n",
                r->q_target, r->noisy, r->q_geo, r->n_plus, r->n_minus, r->index_zero_mode, r->index_trace,
                r->max_nonzero_chiral, r->max_eig_err, (i + 1 < N_IDX_CFGS) ? "," : "");
    }
    fprintf(fp, "  ],\n");

    fprintf(fp, "  \"exp2_schwinger\": [\n");
    for (int i = 0; i < N_BETAS_HMC; ++i) {
        fprintf(fp, "    {\n");
        fprintf(fp, "      \"beta\": %.2f,\n", betas_hmc[i]);
        fprintf(fp, "      \"q_p\": %.8f, \"q_p_err\": %.8f, \"q_chi22\": %.8f, \"q_chi22_err\": %.8f, \"q_pbp\": %.8f, \"q_pbp_err\": %.8f, \"q_acc\": %.6f,\n",
                runs_q[i].p_mean, runs_q[i].p_err, runs_q[i].chi22, runs_q[i].chi22_err, runs_q[i].pbp_mean, runs_q[i].pbp_err, runs_q[i].acc_rate);
        fprintf(fp, "      \"d_p\": %.8f, \"d_p_err\": %.8f, \"d_chi22\": %.8f, \"d_chi22_err\": %.8f, \"d_pbp\": %.8f, \"d_pbp_err\": %.8f, \"d_acc\": %.6f\n",
                runs_d[i].p_mean, runs_d[i].p_err, runs_d[i].chi22, runs_d[i].chi22_err, runs_d[i].pbp_mean, runs_d[i].pbp_err, runs_d[i].acc_rate);
        fprintf(fp, "    }%s\n", (i + 1 < N_BETAS_HMC) ? "," : "");
    }
    fprintf(fp, "  ],\n");

    fprintf(fp, "  \"exp2_mass_scan\": [\n");
    for (int i = 0; i < N_MASS_HMC; ++i) {
        fprintf(fp, "    {\"mass\": %.2f, \"p_mean\": %.8f, \"pbp_mean\": %.8f, \"pbp_err\": %.8f, \"chi22\": %.8f, \"acc\": %.6f}%s\n",
                runs_m[i].mass, runs_m[i].p_mean, runs_m[i].pbp_mean, runs_m[i].pbp_err, runs_m[i].chi22, runs_m[i].acc_rate,
                (i + 1 < N_MASS_HMC) ? "," : "");
    }
    fprintf(fp, "  ],\n");

    fprintf(fp, "  \"exp3a_chain\": [\n");
    for (int i = 0; i < N_L_CHAIN; ++i) {
        fprintf(fp, "    {\"L\": %d, \"pq_sgn\": %.8f, \"pq_sgn_err\": %.8f, \"n_rw\": %.8f, \"n_rw_err\": %.8f, "
                    "\"dual_sgn\": %.8f, \"n_dual\": %.8f, \"n_dual_err\": %.8f}%s\n",
                exp3a[i].L, exp3a[i].pq_sgn, exp3a[i].pq_sgn_err, exp3a[i].n_rw, exp3a[i].n_rw_err,
                exp3a[i].dual_sgn, exp3a[i].n_dual, exp3a[i].n_dual_err, (i + 1 < N_L_CHAIN) ? "," : "");
    }
    fprintf(fp, "  ],\n");

    fprintf(fp, "  \"exp3b_thimble\": [\n");
    for (int i = 0; i < N_MU_THIM; ++i) {
        fprintf(fp, "    {\"mu\": %.2f, \"y0\": %.8f, \"naive_sgn\": %.8f, \"naive_re_z2\": %.8f, \"naive_re_z2_err\": %.8f, "
                    "\"th_sgn\": %.8f, \"th_re_z2\": %.8f, \"th_re_z2_err\": %.8f, \"th_im_z\": %.8f, \"th_im_z_err\": %.8f}%s\n",
                exp3b[i].mu, exp3b[i].y0, exp3b[i].naive_sgn, exp3b[i].naive_re_z2, exp3b[i].naive_re_z2_err,
                exp3b[i].th_sgn, exp3b[i].th_re_z2, exp3b[i].th_re_z2_err, exp3b[i].th_im_z, exp3b[i].th_im_z_err,
                (i + 1 < N_MU_THIM) ? "," : "");
    }
    fprintf(fp, "  ]\n");
    fprintf(fp, "}\n");
    fclose(fp);
    return 0;
}

