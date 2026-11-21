/* u1_2d_windtunnel.c — 2D Compact U(1) Lattice Gauge Theory Topological Freezing Wind Tunnel */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    L_SIZE   = 16,
    V_SITES  = L_SIZE * L_SIZE,
    N_BETAS  = 5,
    N_THERM  = 2000,
    N_MEAS   = 50000
};

static const double K_PI     = 3.14159265358979323846264338327950288;
static const double K_TWO_PI = 6.28318530717958647692528676655900577;

/* --- Deterministic xoshiro256** PRNG (self-contained, zero external deps) --- */
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

static inline uint32_t rng_below(rng_t *r, uint32_t bound) {
    return (uint32_t)(((uint64_t)(uint32_t)rng_u64(r) * (uint64_t)bound) >> 32);
}

/* --- Compact U(1) Torus Geometry --- */
typedef struct {
    double ax[L_SIZE][L_SIZE];
    double ay[L_SIZE][L_SIZE];
} gauge2d_t;

static inline double wrap_pi(double x) {
    x = fmod(x + K_PI, K_TWO_PI);
    if (x < 0.0) {
        x += K_TWO_PI;
    }
    x -= K_PI;
    if (x <= -K_PI) {
        x += K_TWO_PI;
    }
    return x;
}

static inline int idx_next(int i) {
    return (i + 1 < L_SIZE) ? (i + 1) : 0;
}

static inline int idx_prev(int i) {
    return (i > 0) ? (i - 1) : (L_SIZE - 1);
}

static inline double plaq_theta(const gauge2d_t *g, int x, int y) {
    int xp = idx_next(x);
    int yp = idx_next(y);
    return wrap_pi(g->ax[x][y] + g->ay[xp][y] - g->ax[x][yp] - g->ay[x][y]);
}

static void gauge_zero(gauge2d_t *g) {
    memset(g, 0, sizeof(*g));
}

static void gauge_observables(const gauge2d_t *g, double beta,
                              double *out_action, double *out_mean_plaq, long long *out_q) {
    double s_cos   = 0.0;
    double s_theta = 0.0;
    for (int x = 0; x < L_SIZE; ++x) {
        for (int y = 0; y < L_SIZE; ++y) {
            double th = plaq_theta(g, x, y);
            s_cos += cos(th);
            s_theta += th;
        }
    }
    *out_mean_plaq = s_cos / (double)V_SITES;
    *out_action    = beta * ((double)V_SITES - s_cos);
    *out_q         = llround(s_theta / K_TWO_PI);
}

static void init_harmonic_instanton(gauge2d_t *h, int sign) {
    const double delta0 = (double)sign * (K_TWO_PI / (double)V_SITES);
    for (int x = 0; x < L_SIZE; ++x) {
        for (int y = 0; y < L_SIZE; ++y) {
            h->ax[x][y] = wrap_pi(-(double)y * delta0);
            h->ay[x][y] = (y + 1 == L_SIZE) ? wrap_pi((double)L_SIZE * (double)x * delta0) : 0.0;
        }
    }
}

/* --- Local Metropolis (fixed step_rad during measurement -> exact detailed balance) --- */
static void sweep_local_metropolis(gauge2d_t *g, double beta, double step_rad,
                                   rng_t *rng, double *out_acc_rate) {
    long long tries = 0;
    long long acc   = 0;

    for (int x = 0; x < L_SIZE; ++x) {
        int xp = idx_next(x);
        int xm = idx_prev(x);
        for (int y = 0; y < L_SIZE; ++y) {
            int yp = idx_next(y);
            int ym = idx_prev(y);

            /* Update ax[x][y]: affects plaq(x, y) and plaq(x, ym) */
            {
                double old_l = g->ax[x][y];
                double r_up  = g->ay[xp][y] - g->ax[x][yp] - g->ay[x][y];
                double r_dn  = g->ax[x][ym] + g->ay[xp][ym] - g->ay[x][ym];
                double th_up_old = old_l + r_up;
                double th_dn_old = r_dn - old_l;
                double d_link    = (2.0 * rng_u01(rng) - 1.0) * step_rad;
                double new_l     = wrap_pi(old_l + d_link);
                double th_up_new = new_l + r_up;
                double th_dn_new = r_dn - new_l;
                double dS = beta * ((cos(th_up_old) + cos(th_dn_old)) -
                                    (cos(th_up_new) + cos(th_dn_new)));
                ++tries;
                if (dS <= 0.0 || rng_pos_u01(rng) < exp(-dS)) {
                    g->ax[x][y] = new_l;
                    ++acc;
                }
            }

            /* Update ay[x][y]: affects plaq(x, y) and plaq(xm, y) */
            {
                double old_l = g->ay[x][y];
                double r_rt  = g->ax[x][y] + g->ay[xp][y] - g->ax[x][yp];
                double r_lf  = g->ax[xm][y] - g->ax[xm][yp] - g->ay[xm][y];
                double th_rt_old = r_rt - old_l;
                double th_lf_old = old_l + r_lf;
                double d_link    = (2.0 * rng_u01(rng) - 1.0) * step_rad;
                double new_l     = wrap_pi(old_l + d_link);
                double th_rt_new = r_rt - new_l;
                double th_lf_new = new_l + r_lf;
                double dS = beta * ((cos(th_rt_old) + cos(th_lf_old)) -
                                    (cos(th_rt_new) + cos(th_lf_new)));
                ++tries;
                if (dS <= 0.0 || rng_pos_u01(rng) < exp(-dS)) {
                    g->ay[x][y] = new_l;
                    ++acc;
                }
            }
        }
    }
    if (out_acc_rate != NULL) {
        *out_acc_rate = (double)acc / (double)tries;
    }
}

/* --- Microcanonical Over-Relaxation (exact dS = 0 reflection across staple phase) --- */
static void sweep_overrelaxation(gauge2d_t *g) {
    for (int x = 0; x < L_SIZE; ++x) {
        int xp = idx_next(x);
        int xm = idx_prev(x);
        for (int y = 0; y < L_SIZE; ++y) {
            int yp = idx_next(y);
            int ym = idx_prev(y);

            /* ax[x][y] */
            {
                double s_up = g->ax[x][yp] + g->ay[x][y] - g->ay[xp][y];
                double s_dn = g->ax[x][ym] + g->ay[xp][ym] - g->ay[x][ym];
                double rx   = cos(s_up) + cos(s_dn);
                double ry   = sin(s_up) + sin(s_dn);
                if (rx * rx + ry * ry > 1e-24) {
                    double phi  = atan2(ry, rx);
                    g->ax[x][y] = wrap_pi(2.0 * phi - g->ax[x][y]);
                }
            }

            /* ay[x][y] */
            {
                double s_rt = g->ax[x][y] + g->ay[xp][y] - g->ax[x][yp];
                double s_lf = g->ax[xm][yp] + g->ay[xm][y] - g->ax[xm][y];
                double rx   = cos(s_rt) + cos(s_lf);
                double ry   = sin(s_rt) + sin(s_lf);
                if (rx * rx + ry * ry > 1e-24) {
                    double phi  = atan2(ry, rx);
                    g->ay[x][y] = wrap_pi(2.0 * phi - g->ay[x][y]);
                }
            }
        }
    }
}

/* --- Instanton Proposals (Symmetric, det J = 1, Exact Metropolis Acceptance) --- */
static int propose_bare_core_instanton(gauge2d_t *g, double beta, rng_t *rng) {
    int x0       = (int)rng_below(rng, L_SIZE);
    int y0       = (int)rng_below(rng, L_SIZE);
    double sign  = (rng_u64(rng) & 1ULL) ? 1.0 : -1.0;
    int x1       = idx_next(x0);
    int y1       = idx_next(y0);

    double s_old = 0.0;
    double s_new = 0.0;
    const double d_half = sign * (0.5 * K_PI);
    const double d_pi   = sign * K_PI;

    /* Core links modified: ax(x1,y1) -= d_pi, ay(x1,y0) += d_half, ay(x1,y1) += d_half.
       Affected plaquettes: 4 core plaquettes + up to 2 neighbor plaquettes whose raw delta = 0
       modulo 2pi! Actually ONLY raw plaquettes (x0,y0),(x1,y0),(x0,y1),(x1,y1) change modulo 2pi. */
    int xs[2] = {x0, x1};
    int ys[2] = {y0, y1};
    for (int ix = 0; ix < 2; ++ix) {
        for (int iy = 0; iy < 2; ++iy) {
            double th = plaq_theta(g, xs[ix], ys[iy]);
            s_old += 1.0 - cos(th);
        }
    }

    g->ax[x1][y1] = wrap_pi(g->ax[x1][y1] - d_pi);
    g->ay[x1][y0] = wrap_pi(g->ay[x1][y0] + d_half);
    g->ay[x1][y1] = wrap_pi(g->ay[x1][y1] + d_half);

    for (int ix = 0; ix < 2; ++ix) {
        for (int iy = 0; iy < 2; ++iy) {
            double th = plaq_theta(g, xs[ix], ys[iy]);
            s_new += 1.0 - cos(th);
        }
    }

    double dS = beta * (s_new - s_old);
    if (dS <= 0.0 || rng_pos_u01(rng) < exp(-dS)) {
        return 1;
    }

    /* Reject: exact undo */
    g->ax[x1][y1] = wrap_pi(g->ax[x1][y1] + d_pi);
    g->ay[x1][y0] = wrap_pi(g->ay[x1][y0] - d_half);
    g->ay[x1][y1] = wrap_pi(g->ay[x1][y1] - d_half);
    return 0;
}

static int propose_harmonic_instanton(gauge2d_t *g, const gauge2d_t *harm_pos,
                                      double beta, rng_t *rng) {
    double sign = (rng_u64(rng) & 1ULL) ? 1.0 : -1.0;
    double d_sum_cos = 0.0;
    const double delta0 = sign * (K_TWO_PI / (double)V_SITES);

    for (int x = 0; x < L_SIZE; ++x) {
        for (int y = 0; y < L_SIZE; ++y) {
            double th_old = plaq_theta(g, x, y);
            double th_new = th_old + delta0;
            d_sum_cos += cos(th_old) - cos(th_new);
        }
    }
    double dS = beta * d_sum_cos;
    if (dS <= 0.0 || rng_pos_u01(rng) < exp(-dS)) {
        for (int x = 0; x < L_SIZE; ++x) {
            for (int y = 0; y < L_SIZE; ++y) {
                g->ax[x][y] = wrap_pi(g->ax[x][y] + sign * harm_pos->ax[x][y]);
                g->ay[x][y] = wrap_pi(g->ay[x][y] + sign * harm_pos->ay[x][y]);
            }
        }
        return 1;
    }
    return 0;
}

/* --- Built-In Startup Invariance & Physics Self-Tests --- */
static void run_startup_selftests(void) {
    gauge2d_t g, gp, h;
    double lam[L_SIZE][L_SIZE];
    rng_t rng;
    rng_seed(&rng, 0xF10E2026ULL);

    /* 1. Local U(1) gauge invariance of plaquettes, action, and topological charge */
    for (int x = 0; x < L_SIZE; ++x) {
        for (int y = 0; y < L_SIZE; ++y) {
            g.ax[x][y]  = (2.0 * rng_u01(&rng) - 1.0) * K_PI;
            g.ay[x][y]  = (2.0 * rng_u01(&rng) - 1.0) * K_PI;
            lam[x][y]   = (2.0 * rng_u01(&rng) - 1.0) * 5.0 * K_PI;
        }
    }
    for (int x = 0; x < L_SIZE; ++x) {
        int xp = idx_next(x);
        for (int y = 0; y < L_SIZE; ++y) {
            int yp       = idx_next(y);
            gp.ax[x][y]  = wrap_pi(g.ax[x][y] + lam[xp][y] - lam[x][y]);
            gp.ay[x][y]  = wrap_pi(g.ay[x][y] + lam[x][yp] - lam[x][y]);
        }
    }
    double s0, p0, s1, p1;
    long long q0, q1;
    gauge_observables(&g, 6.0, &s0, &p0, &q0);
    gauge_observables(&gp, 6.0, &s1, &p1, &q1);
    for (int x = 0; x < L_SIZE; ++x) {
        for (int y = 0; y < L_SIZE; ++y) {
            double d = fabs(plaq_theta(&g, x, y) - plaq_theta(&gp, x, y));
            if (d > 1e-11) {
                fprintf(stderr, "SELFTEST FAIL: gauge plaq diff %.3e at (%d,%d)\n", d, x, y);
                exit(1);
            }
        }
    }
    if (fabs(s0 - s1) > 1e-10 || q0 != q1) {
        fprintf(stderr, "SELFTEST FAIL: gauge action/Q mismatch\n");
        exit(1);
    }

    /* 2. Over-relaxation microcanonical exactness (dS = 0) */
    sweep_overrelaxation(&gp);
    double s_or, p_or;
    long long q_or;
    gauge_observables(&gp, 6.0, &s_or, &p_or, &q_or);
    if (fabs(s1 - s_or) > 1e-9) {
        fprintf(stderr, "SELFTEST FAIL: overrelaxation dS=%.3e != 0\n", s1 - s_or);
        exit(1);
    }

    /* 3. Harmonic 2-Jet instanton plaquette uniformity, Q=+1, exact action */
    init_harmonic_instanton(&h, 1);
    const double delta0 = K_TWO_PI / (double)V_SITES;
    for (int x = 0; x < L_SIZE; ++x) {
        for (int y = 0; y < L_SIZE; ++y) {
            double d = fabs(plaq_theta(&h, x, y) - delta0);
            if (d > 1e-12) {
                fprintf(stderr, "SELFTEST FAIL: harmonic plaq mismatch %.3e at (%d,%d)\n", d, x, y);
                exit(1);
            }
        }
    }
    double sh, ph;
    long long qh;
    gauge_observables(&h, 10.0, &sh, &ph, &qh);
    double sh_exp = 10.0 * (double)V_SITES * (1.0 - cos(delta0));
    if (qh != 1 || fabs(sh - sh_exp) > 1e-10) {
        fprintf(stderr, "SELFTEST FAIL: harmonic Q=%lld S=%.8f exp=%.8f\n", qh, sh, sh_exp);
        exit(1);
    }

    /* 4. Bare 3-link core instanton check: 4 core plaquettes = +pi/2, 252 zero, Q=+1, S = 4*beta */
    gauge_zero(&g);
    g.ax[1][1] = -K_PI;
    g.ay[1][0] = 0.5 * K_PI;
    g.ay[1][1] = 0.5 * K_PI;
    for (int x = 0; x < L_SIZE; ++x) {
        for (int y = 0; y < L_SIZE; ++y) {
            double th  = plaq_theta(&g, x, y);
            double exp = ((x <= 1) && (y <= 1)) ? (0.5 * K_PI) : 0.0;
            if (fabs(th - exp) > 1e-12) {
                fprintf(stderr, "SELFTEST FAIL: bare core plaq (%d,%d) got %.6f exp %.6f\n",
                        x, y, th, exp);
                exit(1);
            }
        }
    }
    gauge_observables(&g, 10.0, &s0, &p0, &q0);
    if (q0 != 1 || fabs(s0 - 40.0) > 1e-10) {
        fprintf(stderr, "SELFTEST FAIL: bare core Q=%lld S=%.6f\n", q0, s0);
        exit(1);
    }
    printf("  [SELFTEST] PASS (gauge invariance <1e-11, OR microcanonical <1e-9, harmonic Q=+1, bare Q=+1)\n");
}

/* --- Madras-Sokal Automatic-Windowed Integrated Autocorrelation Time --- */
static void compute_madras_sokal(const double *x, int n, double *out_mean,
                                 double *out_var, double *out_tau_int,
                                 double *out_err, int *out_frozen) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += x[i];
    }
    double mean = sum / (double)n;

    double c0 = 0.0;
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

    *out_frozen   = 0;
    double tau    = 0.5;
    int max_w     = n / 2;
    for (int w = 1; w <= max_w; ++w) {
        double ct = 0.0;
        int lim   = n - w;
        for (int i = 0; i < lim; ++i) {
            ct += (x[i] - mean) * (x[i + w] - mean);
        }
        ct /= (double)lim;
        double rho = ct / c0;
        tau += rho;
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

typedef enum {
    SOLVER_LOCAL_METRO_OR = 0,
    SOLVER_BARE_INSTANTON = 1,
    SOLVER_SYMBIOTIC_FLOW = 2
} solver_kind_t;

typedef struct {
    const char *name;
    double beta;
    double step_rad;
    double mean_plaq;
    double err_plaq;
    double tau_int_plaq;
    double chi_t;
    double err_chi_t;
    double tau_int_q;
    double tau_int_q2;
    int frozen_q;
    long long q_transitions;
    double q_trans_rate;
    double inst_accept_rate;
} bench_res_t;

static bench_res_t run_chain(solver_kind_t kind, const char *name, double beta,
                             const gauge2d_t *harm_pos, uint64_t seed) {
    gauge2d_t g;
    gauge_zero(&g);
    rng_t rng;
    rng_seed(&rng, seed);

    double *p_hist  = (double *)malloc((size_t)N_MEAS * sizeof(double));
    double *q_hist  = (double *)malloc((size_t)N_MEAS * sizeof(double));
    double *q2_hist = (double *)malloc((size_t)N_MEAS * sizeof(double));
    if (p_hist == NULL || q_hist == NULL || q2_hist == NULL) {
        fprintf(stderr, "ERROR: OOM allocating histories\n");
        exit(1);
    }

    double step_rad = 1.6 / sqrt(beta);
    if (step_rad > K_PI) {
        step_rad = K_PI;
    }

    /* Thermalization: adapt step_rad slightly during first half, then FREEZE before measurement */
    for (int t = 0; t < N_THERM; ++t) {
        double acc = 0.0;
        sweep_local_metropolis(&g, beta, step_rad, &rng, &acc);
        sweep_overrelaxation(&g);
        if (kind == SOLVER_BARE_INSTANTON) {
            propose_bare_core_instanton(&g, beta, &rng);
        } else if (kind == SOLVER_SYMBIOTIC_FLOW) {
            propose_harmonic_instanton(&g, harm_pos, beta, &rng);
        }
        if (t < N_THERM / 2 && (t % 25) == 24) {
            if (acc > 0.58) {
                step_rad *= 1.08;
                if (step_rad > K_PI) {
                    step_rad = K_PI;
                }
            } else if (acc < 0.42) {
                step_rad *= 0.92;
            }
        }
    }

    /* Frozen step_rad measurement phase -> exact detailed balance */
    long long inst_acc   = 0;
    long long q_trans    = 0;
    long long prev_q     = 0;
    for (int t = 0; t < N_MEAS; ++t) {
        sweep_local_metropolis(&g, beta, step_rad, &rng, NULL);
        sweep_overrelaxation(&g);
        if (kind == SOLVER_BARE_INSTANTON) {
            inst_acc += propose_bare_core_instanton(&g, beta, &rng);
        } else if (kind == SOLVER_SYMBIOTIC_FLOW) {
            inst_acc += propose_harmonic_instanton(&g, harm_pos, beta, &rng);
        }
        double act, mp;
        long long q;
        gauge_observables(&g, beta, &act, &mp, &q);
        p_hist[t]  = mp;
        q_hist[t]  = (double)q;
        q2_hist[t] = (double)(q * q) / (double)V_SITES;
        if (t > 0 && q != prev_q) {
            ++q_trans;
        }
        prev_q = q;
    }

    double p_mean, p_var, p_tau, p_err;
    double q_mean, q_var, q_tau, q_err;
    double q2_mean, q2_var, q2_tau, q2_err;
    int p_frz, q_frz, q2_frz;

    compute_madras_sokal(p_hist, N_MEAS, &p_mean, &p_var, &p_tau, &p_err, &p_frz);
    compute_madras_sokal(q_hist, N_MEAS, &q_mean, &q_var, &q_tau, &q_err, &q_frz);
    compute_madras_sokal(q2_hist, N_MEAS, &q2_mean, &q2_var, &q2_tau, &q2_err, &q2_frz);

    free(p_hist);
    free(q_hist);
    free(q2_hist);

    bench_res_t res;
    res.name             = name;
    res.beta             = beta;
    res.step_rad         = step_rad;
    res.mean_plaq        = p_mean;
    res.err_plaq         = p_err;
    res.tau_int_plaq     = p_tau;
    res.chi_t            = q2_mean;
    res.err_chi_t        = q2_err;
    res.tau_int_q        = q_tau;
    res.tau_int_q2       = q2_tau;
    res.frozen_q         = q_frz;
    res.q_transitions    = q_trans;
    res.q_trans_rate     = (double)q_trans / (double)(N_MEAS - 1);
    res.inst_accept_rate = (kind == SOLVER_LOCAL_METRO_OR) ? 0.0 : ((double)inst_acc / (double)N_MEAS);
    return res;
}

int main(int argc, char **argv) {
    const char *json_path = NULL;
    for (int i = 1; i + 1 < argc; ++i) {
        if (strcmp(argv[i], "--json") == 0) {
            json_path = argv[i + 1];
        }
    }

    printf("=== 2D Compact U(1) Topological Freezing Wind-Tunnel (L=%d, V=%d, N_meas=%d) ===\n",
           L_SIZE, V_SITES, N_MEAS);
    run_startup_selftests();

    gauge2d_t harm_pos;
    init_harmonic_instanton(&harm_pos, 1);

    const double betas[N_BETAS]      = {2.0, 4.0, 6.0, 8.0, 10.0};
    const solver_kind_t kinds[3]     = {SOLVER_LOCAL_METRO_OR, SOLVER_BARE_INSTANTON, SOLVER_SYMBIOTIC_FLOW};
    const char *names[3]             = {"local_metro_or", "bare_instanton", "symbiotic_flow"};
    bench_res_t results[N_BETAS][3];

    for (int b = 0; b < N_BETAS; ++b) {
        printf("\n  --- beta = %.1f ---\n", betas[b]);
        for (int s = 0; s < 3; ++s) {
            uint64_t seed = 0x202611210000ULL + (uint64_t)b * 100ULL + (uint64_t)s;
            results[b][s] = run_chain(kinds[s], names[s], betas[b], &harm_pos, seed);
            const bench_res_t *r = &results[b][s];
            if (r->frozen_q) {
                printf("    %-16s | <P>=%.6f+-%.6f | chi_t=%.6f (FROZEN) | tau_Q=INF   | dQ/sweep=%.5f | inst_acc=%.4f\n",
                       r->name, r->mean_plaq, r->err_plaq, r->chi_t, r->q_trans_rate, r->inst_accept_rate);
            } else {
                printf("    %-16s | <P>=%.6f+-%.6f | chi_t=%.6f+-%.6f | tau_Q=%-6.2f | dQ/sweep=%.5f | inst_acc=%.4f\n",
                       r->name, r->mean_plaq, r->err_plaq, r->chi_t, r->err_chi_t,
                       r->tau_int_q, r->q_trans_rate, r->inst_accept_rate);
            }
        }
    }

    if (json_path != NULL) {
        FILE *fp = fopen(json_path, "w");
        if (fp == NULL) {
            fprintf(stderr, "ERROR: cannot write %s\n", json_path);
            return 1;
        }
        fprintf(fp, "{\n");
        fprintf(fp, "  \"L\": %d,\n", L_SIZE);
        fprintf(fp, "  \"V\": %d,\n", V_SITES);
        fprintf(fp, "  \"n_therm\": %d,\n", N_THERM);
        fprintf(fp, "  \"n_meas\": %d,\n", N_MEAS);
        fprintf(fp, "  \"runs\": [\n");
        int total = N_BETAS * 3;
        int idx   = 0;
        for (int b = 0; b < N_BETAS; ++b) {
            for (int s = 0; s < 3; ++s) {
                const bench_res_t *r = &results[b][s];
                ++idx;
                fprintf(fp, "    {\n");
                fprintf(fp, "      \"beta\": %.1f,\n", r->beta);
                fprintf(fp, "      \"solver\": \"%s\",\n", r->name);
                fprintf(fp, "      \"step_rad\": %.6f,\n", r->step_rad);
                fprintf(fp, "      \"mean_plaq\": %.10f,\n", r->mean_plaq);
                fprintf(fp, "      \"err_plaq\": %.10f,\n", r->err_plaq);
                fprintf(fp, "      \"tau_int_plaq\": %.4f,\n", r->tau_int_plaq);
                fprintf(fp, "      \"chi_t\": %.10f,\n", r->chi_t);
                fprintf(fp, "      \"err_chi_t\": %.10f,\n", r->err_chi_t);
                fprintf(fp, "      \"tau_int_q\": %.4f,\n", r->tau_int_q);
                fprintf(fp, "      \"tau_int_q2\": %.4f,\n", r->tau_int_q2);
                fprintf(fp, "      \"frozen_q\": %s,\n", r->frozen_q ? "true" : "false");
                fprintf(fp, "      \"q_transitions\": %lld,\n", r->q_transitions);
                fprintf(fp, "      \"q_trans_rate\": %.8f,\n", r->q_trans_rate);
                fprintf(fp, "      \"inst_accept_rate\": %.8f\n", r->inst_accept_rate);
                fprintf(fp, "    }%s\n", (idx < total) ? "," : "");
            }
        }
        fprintf(fp, "  ]\n");
        fprintf(fp, "}\n");
        fclose(fp);
    }

    return 0;
}
