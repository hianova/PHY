/* cpn_2d_windtunnel.c — 2D CP^{20} (N=21) Nonlinear Sigma Model Equal-Work Wind Tunnel
 * Tests Theorems P1–P4 from `Flow.ProjectiveBundle` across beta in {0.60, 0.65, 0.70, 0.73}:
 *   1. local_or_metro      (R=1,  30,000 measurement sweeps)
 *   2. bare_bp_inst        (R=1,  30,000 measurement sweeps + 1 BP SU(2) instanton proposal/sweep)
 *   3. symbiotic_ptbc_2jet (R=12 replicas, 2,500 PTBC sweeps = 30,000 equal replica-sweeps, Fisher 2-Jet ladder)
 */
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    L_SIZE       = 32,
    V_SITES      = L_SIZE * L_SIZE,
    N_COMP       = 21,
    L_DEFECT     = 8,
    R_REPLICAS   = 12,
    N_BETAS      = 4,
    N_EQ_MEAS    = 30000,
    N_THERM_1R   = 2400,
    N_THERM_PT   = 800,
    N_MEAS_PT    = N_EQ_MEAS / R_REPLICAS, /* 2500 -> exact equal work 12*2500 = 30000 */
    COOL_STEPS   = 10,
    JK_BLOCKS    = 25
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

/* --- Lattice Structure & Memory Layout --- */
typedef struct {
    double re[L_SIZE][L_SIZE][N_COMP];
    double im[L_SIZE][L_SIZE][N_COMP];
    double ux_re[L_SIZE][L_SIZE];
    double ux_im[L_SIZE][L_SIZE];
    double uy_re[L_SIZE][L_SIZE];
    double uy_im[L_SIZE][L_SIZE];
} cpn_cfg_t;

static inline int idx_next(int i) {
    return (i + 1 < L_SIZE) ? (i + 1) : 0;
}

static inline int idx_prev(int i) {
    return (i > 0) ? (i - 1) : (L_SIZE - 1);
}

/* Defect slit: L_DEFECT vertical links (x, y_slit) -> (x, y_slit + 1) for x in [0, L_DEFECT) */
static inline int is_defect_y_link(int x, int y) {
    return (y == (L_SIZE / 2)) && (x >= 0 && x < L_DEFECT);
}

static void cfg_init_ordered(cpn_cfg_t *g) {
    memset(g, 0, sizeof(*g));
    for (int x = 0; x < L_SIZE; ++x) {
        for (int y = 0; y < L_SIZE; ++y) {
            g->re[x][y][0] = 1.0;
            g->ux_re[x][y] = 1.0;
            g->uy_re[x][y] = 1.0;
        }
    }
}

static void normalize_site(double *re, double *im) {
    double n2 = 0.0;
    for (int a = 0; a < N_COMP; ++a) {
        n2 += re[a] * re[a] + im[a] * im[a];
    }
    double inv = 1.0 / sqrt(n2);
    for (int a = 0; a < N_COMP; ++a) {
        re[a] *= inv;
        im[a] *= inv;
    }
}

static inline void inner_dot(const double *are, const double *aim,
                             const double *bre, const double *bim,
                             double *out_re, double *out_im) {
    double sr = 0.0;
    double si = 0.0;
    for (int a = 0; a < N_COMP; ++a) {
        sr += are[a] * bre[a] + aim[a] * bim[a];
        si += are[a] * bim[a] - aim[a] * bre[a];
    }
    *out_re = sr;
    *out_im = si;
}

static inline double hop_re(const cpn_cfg_t *g, int x1, int y1,
                            double ur, double ui, int x2, int y2) {
    double dr, di;
    inner_dot(g->re[x1][y1], g->im[x1][y1], g->re[x2][y2], g->im[x2][y2], &dr, &di);
    return ur * dr - ui * di;
}

static double cfg_defect_sum_re(const cpn_cfg_t *g) {
    double s    = 0.0;
    const int y = L_SIZE / 2;
    const int yp = idx_next(y);
    for (int x = 0; x < L_DEFECT; ++x) {
        s += hop_re(g, x, y, g->uy_re[x][y], g->uy_im[x][y], x, yp);
    }
    return s;
}

static double cfg_total_action(const cpn_cfg_t *g, double beta, double c_r) {
    double sum_hop = 0.0;
    for (int x = 0; x < L_SIZE; ++x) {
        int xp = idx_next(x);
        for (int y = 0; y < L_SIZE; ++y) {
            int yp = idx_next(y);
            sum_hop += hop_re(g, x, y, g->ux_re[x][y], g->ux_im[x][y], xp, y);
            double cy = is_defect_y_link(x, y) ? c_r : 1.0;
            sum_hop += cy * hop_re(g, x, y, g->uy_re[x][y], g->uy_im[x][y], x, yp);
        }
    }
    return -2.0 * (double)N_COMP * beta * sum_hop;
}

static double cfg_mean_link_energy(const cpn_cfg_t *g) {
    double s = 0.0;
    for (int x = 0; x < L_SIZE; ++x) {
        int xp = idx_next(x);
        for (int y = 0; y < L_SIZE; ++y) {
            int yp = idx_next(y);
            s += hop_re(g, x, y, g->ux_re[x][y], g->ux_im[x][y], xp, y);
            s += hop_re(g, x, y, g->uy_re[x][y], g->uy_im[x][y], x, yp);
        }
    }
    return s / (double)(2 * V_SITES);
}

/* Compute Pancharatnam-Berry triangle phase arg( (z1† z2)(z2† z3)(z3† z1) ) */
static inline double triangle_phase(const double *z1r, const double *z1i,
                                    const double *z2r, const double *z2i,
                                    const double *z3r, const double *z3i) {
    double a1r, a1i, a2r, a2i, a3r, a3i;
    inner_dot(z1r, z1i, z2r, z2i, &a1r, &a1i);
    inner_dot(z2r, z2i, z3r, z3i, &a2r, &a2i);
    inner_dot(z3r, z3i, z1r, z1i, &a3r, &a3i);

    double pr = a1r * a2r - a1i * a2i;
    double pi = a1r * a2i + a1i * a2r;
    double wr = pr * a3r - pi * a3i;
    double wi = pr * a3i + pi * a3r;
    return atan2(wi, wr);
}

static long long cfg_topo_charge_raw(const cpn_cfg_t *g, double *out_sum_rad) {
    double sum_rad = 0.0;
    for (int x = 0; x < L_SIZE; ++x) {
        int xp = idx_next(x);
        for (int y = 0; y < L_SIZE; ++y) {
            int yp = idx_next(y);
            sum_rad += triangle_phase(g->re[x][y], g->im[x][y],
                                      g->re[xp][y], g->im[xp][y],
                                      g->re[xp][yp], g->im[xp][yp]);
            sum_rad += triangle_phase(g->re[x][y], g->im[x][y],
                                      g->re[xp][yp], g->im[xp][yp],
                                      g->re[x][yp], g->im[x][yp]);
        }
    }
    if (out_sum_rad != NULL) {
        *out_sum_rad = sum_rad;
    }
    return llround(sum_rad / K_TWO_PI);
}

/* Deterministic local cooling on a scratch copy (strictly periodic c=1) */
static long long cfg_topo_charge_cooled(const cpn_cfg_t *src, cpn_cfg_t *wrk) {
    memcpy(wrk, src, sizeof(*wrk));
    for (int sweep = 0; sweep < COOL_STEPS; ++sweep) {
        /* 1. Align links deterministically to U_mu(x) = z_x z_{x+mu}† / |.| (conjugate of z_x† z_{x+mu}) */
        for (int x = 0; x < L_SIZE; ++x) {
            int xp = idx_next(x);
            for (int y = 0; y < L_SIZE; ++y) {
                int yp = idx_next(y);
                double dr, di, mag;
                inner_dot(wrk->re[x][y], wrk->im[x][y], wrk->re[xp][y], wrk->im[xp][y], &dr, &di);
                mag = hypot(dr, di);
                if (mag > 1e-12) {
                    wrk->ux_re[x][y] = dr / mag;
                    wrk->ux_im[x][y] = -di / mag;
                }
                inner_dot(wrk->re[x][y], wrk->im[x][y], wrk->re[x][yp], wrk->im[x][yp], &dr, &di);
                mag = hypot(dr, di);
                if (mag > 1e-12) {
                    wrk->uy_re[x][y] = dr / mag;
                    wrk->uy_im[x][y] = -di / mag;
                }
            }
        }
        /* 2. Align each site z_x directly to its normalized covariant staple F_x */
        for (int x = 0; x < L_SIZE; ++x) {
            int xp = idx_next(x);
            int xm = idx_prev(x);
            for (int y = 0; y < L_SIZE; ++y) {
                int yp = idx_next(y);
                int ym = idx_prev(y);
                double fr[N_COMP];
                double fi[N_COMP];
                double n2 = 0.0;
                for (int a = 0; a < N_COMP; ++a) {
                    double r = 0.0;
                    double i = 0.0;
                    /* +x: U_x(x,y) * z(xp,y) */
                    r += wrk->ux_re[x][y] * wrk->re[xp][y][a] - wrk->ux_im[x][y] * wrk->im[xp][y][a];
                    i += wrk->ux_re[x][y] * wrk->im[xp][y][a] + wrk->ux_im[x][y] * wrk->re[xp][y][a];
                    /* -x: conj U_x(xm,y) * z(xm,y) */
                    r += wrk->ux_re[xm][y] * wrk->re[xm][y][a] + wrk->ux_im[xm][y] * wrk->im[xm][y][a];
                    i += wrk->ux_re[xm][y] * wrk->im[xm][y][a] - wrk->ux_im[xm][y] * wrk->re[xm][y][a];
                    /* +y: U_y(x,y) * z(x,yp) */
                    r += wrk->uy_re[x][y] * wrk->re[x][yp][a] - wrk->uy_im[x][y] * wrk->im[x][yp][a];
                    i += wrk->uy_re[x][y] * wrk->im[x][yp][a] + wrk->uy_im[x][y] * wrk->re[x][yp][a];
                    /* -y: conj U_y(x,ym) * z(x,ym) */
                    r += wrk->uy_re[x][ym] * wrk->re[x][ym][a] + wrk->uy_im[x][ym] * wrk->im[x][ym][a];
                    i += wrk->uy_re[x][ym] * wrk->im[x][ym][a] - wrk->uy_im[x][ym] * wrk->re[x][ym][a];
                    fr[a] = r;
                    fi[a] = i;
                    n2 += r * r + i * i;
                }
                if (n2 > 1e-24) {
                    double inv = 1.0 / sqrt(n2);
                    for (int a = 0; a < N_COMP; ++a) {
                        wrk->re[x][y][a] = fr[a] * inv;
                        wrk->im[x][y][a] = fi[a] * inv;
                    }
                }
            }
        }
    }
    return cfg_topo_charge_raw(wrk, NULL);
}

/* Gauge-invariant projector Fourier power Gtilde(0) and 0.5*(Gtilde(p1,0) + Gtilde(0,p1)) */
static void cfg_projector_fourier(const cpn_cfg_t *g, double *out_g0, double *out_gp1) {
    double pr0[N_COMP][N_COMP];
    double pi0[N_COMP][N_COMP];
    double prx_r[N_COMP][N_COMP], prx_i[N_COMP][N_COMP];
    double pry_r[N_COMP][N_COMP], pry_i[N_COMP][N_COMP];

    memset(pr0, 0, sizeof(pr0));
    memset(pi0, 0, sizeof(pi0));
    memset(prx_r, 0, sizeof(prx_r));
    memset(prx_i, 0, sizeof(prx_i));
    memset(pry_r, 0, sizeof(pry_r));
    memset(pry_i, 0, sizeof(pry_i));

    double cx[L_SIZE], sx[L_SIZE];
    for (int k = 0; k < L_SIZE; ++k) {
        double ang = K_TWO_PI * (double)k / (double)L_SIZE;
        cx[k]      = cos(ang);
        sx[k]      = sin(ang);
    }

    for (int x = 0; x < L_SIZE; ++x) {
        for (int y = 0; y < L_SIZE; ++y) {
            const double *zr = g->re[x][y];
            const double *zi = g->im[x][y];
            double cpx = cx[x], spx = sx[x];
            double cpy = cx[y], spy = sx[y];
            for (int a = 0; a < N_COMP; ++a) {
                double ar = zr[a], ai = zi[a];
                for (int b = 0; b < N_COMP; ++b) {
                    /* P_{ab} = z_a * conj(z_b) - delta_{ab}/N */
                    double er = ar * zr[b] + ai * zi[b];
                    double ei = ai * zr[b] - ar * zi[b];
                    if (a == b) {
                        er -= 1.0 / (double)N_COMP;
                    }
                    pr0[a][b] += er;
                    pi0[a][b] += ei;
                    prx_r[a][b] += cpx * er - spx * ei;
                    prx_i[a][b] += cpx * ei + spx * er;
                    pry_r[a][b] += cpy * er - spy * ei;
                    pry_i[a][b] += cpy * ei + spy * er;
                }
            }
        }
    }

    double g0 = 0.0, gx = 0.0, gy = 0.0;
    for (int a = 0; a < N_COMP; ++a) {
        for (int b = 0; b < N_COMP; ++b) {
            g0 += pr0[a][b] * pr0[a][b] + pi0[a][b] * pi0[a][b];
            gx += prx_r[a][b] * prx_r[a][b] + prx_i[a][b] * prx_i[a][b];
            gy += pry_r[a][b] * pry_r[a][b] + pry_i[a][b] * pry_i[a][b];
        }
    }
    *out_g0  = g0 / (double)V_SITES;
    *out_gp1 = 0.5 * (gx + gy) / (double)V_SITES;
}

/* --- Site Neighborhood Covariant Staple --- */
static void compute_site_staple(const cpn_cfg_t *g, int x, int y, double c_r,
                                double *fr, double *fi) {
    int xp     = idx_next(x);
    int xm     = idx_prev(x);
    int yp     = idx_next(y);
    int ym     = idx_prev(y);
    double cyp = is_defect_y_link(x, y) ? c_r : 1.0;
    double cym = is_defect_y_link(x, ym) ? c_r : 1.0;

    for (int a = 0; a < N_COMP; ++a) {
        double r = 0.0;
        double i = 0.0;
        r += g->ux_re[x][y] * g->re[xp][y][a] - g->ux_im[x][y] * g->im[xp][y][a];
        i += g->ux_re[x][y] * g->im[xp][y][a] + g->ux_im[x][y] * g->re[xp][y][a];

        r += g->ux_re[xm][y] * g->re[xm][y][a] + g->ux_im[xm][y] * g->im[xm][y][a];
        i += g->ux_re[xm][y] * g->im[xm][y][a] - g->ux_im[xm][y] * g->re[xm][y][a];

        r += cyp * (g->uy_re[x][y] * g->re[x][yp][a] - g->uy_im[x][y] * g->im[x][yp][a]);
        i += cyp * (g->uy_re[x][y] * g->im[x][yp][a] + g->uy_im[x][y] * g->re[x][yp][a]);

        r += cym * (g->uy_re[x][ym] * g->re[x][ym][a] + g->uy_im[x][ym] * g->im[x][ym][a]);
        i += cym * (g->uy_re[x][ym] * g->im[x][ym][a] - g->uy_im[x][ym] * g->re[x][ym][a]);

        fr[a] = r;
        fi[a] = i;
    }
}

/* Best-Fisher (1979) exact von Mises sampler for p(psi) propto exp(kappa * cos(psi)) */
static double sample_von_mises(double kappa, rng_t *rng) {
    if (kappa < 1e-10) {
        return (2.0 * rng_u01(rng) - 1.0) * K_PI;
    }
    double tau = 1.0 + sqrt(1.0 + 4.0 * kappa * kappa);
    double rho = (tau - sqrt(2.0 * tau)) / (2.0 * kappa);
    double r   = (1.0 + rho * rho) / (2.0 * rho);
    for (;;) {
        double u1 = rng_u01(rng);
        double z  = cos(K_PI * u1);
        double f  = (1.0 + r * z) / (r + z);
        double c  = kappa * (r - f);
        double u2 = rng_pos_u01(rng);
        if (c * (2.0 - c) - u2 > 0.0 || log(c / u2) + 1.0 - c >= 0.0) {
            double u3  = rng_u01(rng);
            if (f > 1.0) {
                f = 1.0;
            }
            if (f < -1.0) {
                f = -1.0;
            }
            double psi = acos(f);
            return (u3 < 0.5) ? -psi : psi;
        }
    }
}

/* Local Sweep: 1 SU(N) Givens+phase Metropolis + 2 Over-Relaxation + 1 von Mises U(1) heatbath */
static void sweep_cpn_local(cpn_cfg_t *g, double beta, double c_r,
                            double step_eps, rng_t *rng, double *out_site_acc) {
    const double coeff  = 2.0 * (double)N_COMP * beta;
    long long tries     = 0;
    long long acc       = 0;

    /* 1. Symmetric Givens SU(2) + component phase-twist Metropolis per site */
    for (int x = 0; x < L_SIZE; ++x) {
        for (int y = 0; y < L_SIZE; ++y) {
            double fr[N_COMP], fi[N_COMP];
            compute_site_staple(g, x, y, c_r, fr, fi);
            double *zr = g->re[x][y];
            double *zi = g->im[x][y];
            double cur_re, cur_im;
            inner_dot(zr, zi, fr, fi, &cur_re, &cur_im);

            /* 6 random SU(2) Givens plane rotations */
            for (int k = 0; k < 6; ++k) {
                int a = rng_below(rng, N_COMP);
                int b = rng_below(rng, N_COMP - 1);
                if (b >= a) {
                    ++b;
                }
                double th    = (2.0 * rng_u01(rng) - 1.0) * step_eps;
                double phi   = (2.0 * rng_u01(rng) - 1.0) * K_PI;
                double ct    = cos(th);
                double st    = sin(th);
                double cp    = cos(phi);
                double sp    = sin(phi);
                /* b_r + i b_i = st * exp(i*phi) */
                double br = st * cp;
                double bi = st * sp;

                double za_r = zr[a], za_i = zi[a];
                double zb_r = zr[b], zb_i = zi[b];
                /* za' = ct * za - conj(beta) * zb ; zb' = beta * za + ct * zb */
                double nza_r = ct * za_r - (br * zb_r + bi * zb_i);
                double nza_i = ct * za_i - (br * zb_i - bi * zb_r);
                double nzb_r = (br * za_r - bi * za_i) + ct * zb_r;
                double nzb_i = (br * za_i + bi * za_r) + ct * zb_i;

                double old_ab = (za_r * fr[a] + za_i * fi[a]) + (zb_r * fr[b] + zb_i * fi[b]);
                double new_ab = (nza_r * fr[a] + nza_i * fi[a]) + (nzb_r * fr[b] + nzb_i * fi[b]);
                double dS     = coeff * (old_ab - new_ab);
                ++tries;
                if (dS <= 0.0 || rng_pos_u01(rng) < exp(-dS)) {
                    zr[a] = nza_r;
                    zi[a] = nza_i;
                    zr[b] = nzb_r;
                    zi[b] = nzb_i;
                    cur_re += (new_ab - old_ab);
                    ++acc;
                }
            }

            /* 2 random component phase twists */
            for (int k = 0; k < 2; ++k) {
                int a     = rng_below(rng, N_COMP);
                double ph = (2.0 * rng_u01(rng) - 1.0) * step_eps;
                double cp = cos(ph);
                double sp = sin(ph);
                double oar = zr[a], oai = zi[a];
                double nar = oar * cp - oai * sp;
                double nai = oar * sp + oai * cp;
                double old_a = oar * fr[a] + oai * fi[a];
                double new_a = nar * fr[a] + nai * fi[a];
                double dS    = coeff * (old_a - new_a);
                ++tries;
                if (dS <= 0.0 || rng_pos_u01(rng) < exp(-dS)) {
                    zr[a] = nar;
                    zi[a] = nai;
                    ++acc;
                }
            }
            normalize_site(zr, zi);
        }
    }

    /* 2. Two microcanonical S^{2N-1} Over-Relaxation checkerboard sweeps (exact dS = 0) */
    for (int or_sw = 0; or_sw < 2; ++or_sw) {
        for (int par = 0; par < 2; ++par) {
            for (int x = 0; x < L_SIZE; ++x) {
                for (int y = 0; y < L_SIZE; ++y) {
                    if (((x + y) & 1) != par) {
                        continue;
                    }
                    double fr[N_COMP], fi[N_COMP];
                    compute_site_staple(g, x, y, c_r, fr, fi);
                    double fn2 = 0.0;
                    double zf  = 0.0;
                    for (int a = 0; a < N_COMP; ++a) {
                        fn2 += fr[a] * fr[a] + fi[a] * fi[a];
                        zf += g->re[x][y][a] * fr[a] + g->im[x][y][a] * fi[a];
                    }
                    if (fn2 > 1e-24) {
                        double sc = 2.0 * zf / fn2;
                        for (int a = 0; a < N_COMP; ++a) {
                            g->re[x][y][a] = sc * fr[a] - g->re[x][y][a];
                            g->im[x][y][a] = sc * fi[a] - g->im[x][y][a];
                        }
                        normalize_site(g->re[x][y], g->im[x][y]);
                    }
                }
            }
        }
    }

    /* 3. Exact Best-Fisher von Mises heatbath on all auxiliary U(1) link phases */
    for (int x = 0; x < L_SIZE; ++x) {
        int xp = idx_next(x);
        for (int y = 0; y < L_SIZE; ++y) {
            int yp = idx_next(y);
            double dr, di;
            inner_dot(g->re[x][y], g->im[x][y], g->re[xp][y], g->im[xp][y], &dr, &di);
            double magx = hypot(dr, di);
            double psix = sample_von_mises(coeff * magx, rng);
            if (magx > 1e-14) {
                double phiz     = atan2(di, dr);
                double th       = psix - phiz;
                g->ux_re[x][y]  = cos(th);
                g->ux_im[x][y]  = sin(th);
            } else {
                double th      = (2.0 * rng_u01(rng) - 1.0) * K_PI;
                g->ux_re[x][y] = cos(th);
                g->ux_im[x][y] = sin(th);
            }

            inner_dot(g->re[x][y], g->im[x][y], g->re[x][yp], g->im[x][yp], &dr, &di);
            double cy   = is_defect_y_link(x, y) ? c_r : 1.0;
            double magy = hypot(dr, di);
            double psiy = sample_von_mises(coeff * cy * magy, rng);
            if (magy > 1e-14) {
                double phiz    = atan2(di, dr);
                double th      = psiy - phiz;
                g->uy_re[x][y] = cos(th);
                g->uy_im[x][y] = sin(th);
            } else {
                double th      = (2.0 * rng_u01(rng) - 1.0) * K_PI;
                g->uy_re[x][y] = cos(th);
                g->uy_im[x][y] = sin(th);
            }
        }
    }

    if (out_site_acc != NULL) {
        *out_site_acc = (double)acc / (double)tries;
    }
}

/* Uniform torus translation on replica r=0 (exact volume-preserving symmetry of periodic torus) */
static void shift_torus_random(cpn_cfg_t *g, cpn_cfg_t *tmp, rng_t *rng) {
    int sx = rng_below(rng, L_SIZE);
    int sy = rng_below(rng, L_SIZE);
    if (sx == 0 && sy == 0) {
        return;
    }
    memcpy(tmp, g, sizeof(*tmp));
    for (int x = 0; x < L_SIZE; ++x) {
        int nx = (x + sx) % L_SIZE;
        for (int y = 0; y < L_SIZE; ++y) {
            int ny = (y + sy) % L_SIZE;
            memcpy(g->re[nx][ny], tmp->re[x][y], N_COMP * sizeof(double));
            memcpy(g->im[nx][ny], tmp->im[x][y], N_COMP * sizeof(double));
            g->ux_re[nx][ny] = tmp->ux_re[x][y];
            g->ux_im[nx][ny] = tmp->ux_im[x][y];
            g->uy_re[nx][ny] = tmp->uy_re[x][y];
            g->uy_im[nx][ny] = tmp->uy_im[x][y];
        }
    }
}

/* Belavin-Polyakov SU(2) subset SU(N) instanton twist (smooth tapered core inside R_c) */
static void apply_bp_su2_twist(cpn_cfg_t *dst, const cpn_cfg_t *src,
                               int x0, int y0, double rho, int a, int b, int sign, int inv_dir) {
    const double rc  = 10.0;
    const double rc2 = rc * rc;
    memcpy(dst, src, sizeof(*dst));

    double site_c[L_SIZE][L_SIZE];
    double site_sr[L_SIZE][L_SIZE];
    double site_si[L_SIZE][L_SIZE];

    for (int x = 0; x < L_SIZE; ++x) {
        int dx_i = x - x0;
        if (dx_i >= L_SIZE / 2) {
            dx_i -= L_SIZE;
        }
        if (dx_i < -L_SIZE / 2) {
            dx_i += L_SIZE;
        }
        double dx = (double)dx_i - 0.5;
        for (int y = 0; y < L_SIZE; ++y) {
            int dy_i = y - y0;
            if (dy_i >= L_SIZE / 2) {
                dy_i -= L_SIZE;
            }
            if (dy_i < -L_SIZE / 2) {
                dy_i += L_SIZE;
            }
            double dy = (double)dy_i - 0.5;
            double r2 = dx * dx + dy * dy;
            double r  = sqrt(r2);
            double c  = 1.0, sr = 0.0, si = 0.0;
            double phi = -atan2(dy, dx) * (double)sign;
            if (r2 < rc2) {
                double tap   = 1.0 - r2 / rc2;
                double theta = 2.0 * atan((rho / r) * tap * tap);
                c            = cos(0.5 * theta);
                double s     = sin(0.5 * theta);
                sr           = s * cos(phi);
                si           = s * sin(phi);
            }
            site_c[x][y]  = c;
            site_sr[x][y] = sr;
            site_si[x][y] = si;

            double zar    = src->re[x][y][a], zai = src->im[x][y][a];
            double zbr    = src->re[x][y][b], zbi = src->im[x][y][b];
            double sr_rot = (inv_dir > 0) ? sr : -sr;
            double si_rot = (inv_dir > 0) ? si : -si;
            /* R = [[c, -conj(s)], [s, c]] (inv_dir=+1) or R† = [[c, conj(s)], [-s, c]] (inv_dir=-1) */
            dst->re[x][y][a] = c * zar - (sr_rot * zbr + si_rot * zbi);
            dst->im[x][y][a] = c * zai - (sr_rot * zbi - si_rot * zbr);
            dst->re[x][y][b] = (sr_rot * zar - si_rot * zai) + c * zbr;
            dst->im[x][y][b] = (sr_rot * zai + si_rot * zar) + c * zbi;
        }
    }

    /* Covariant Berry link transport factor U' = U * p (inv_dir=+1) or U * conj(p) (inv_dir=-1) */
    for (int x = 0; x < L_SIZE; ++x) {
        int xp = idx_next(x);
        for (int y = 0; y < L_SIZE; ++y) {
            int yp = idx_next(y);
            double gr, gi, mag;
            gr  = site_c[x][y] * site_c[xp][y] + (site_sr[x][y] * site_sr[xp][y] + site_si[x][y] * site_si[xp][y]);
            gi  = -(site_sr[x][y] * site_si[xp][y] - site_si[x][y] * site_sr[xp][y]);
            mag = hypot(gr, gi);
            if (mag > 1e-12) {
                double pr = gr / mag;
                double pi = (inv_dir > 0) ? (gi / mag) : (-gi / mag);
                double ur = src->ux_re[x][y], ui = src->ux_im[x][y];
                dst->ux_re[x][y] = ur * pr - ui * pi;
                dst->ux_im[x][y] = ur * pi + ui * pr;
            }
            gr  = site_c[x][y] * site_c[x][yp] + (site_sr[x][y] * site_sr[x][yp] + site_si[x][y] * site_si[x][yp]);
            gi  = -(site_sr[x][y] * site_si[x][yp] - site_si[x][y] * site_sr[x][yp]);
            mag = hypot(gr, gi);
            if (mag > 1e-12) {
                double pr = gr / mag;
                double pi = (inv_dir > 0) ? (gi / mag) : (-gi / mag);
                double ur = src->uy_re[x][y], ui = src->uy_im[x][y];
                dst->uy_re[x][y] = ur * pr - ui * pi;
                dst->uy_im[x][y] = ur * pi + ui * pr;
            }
        }
    }
}

static int propose_bare_bp_instanton(cpn_cfg_t *g, cpn_cfg_t *prop, double beta, rng_t *rng) {
    int x0     = rng_below(rng, L_SIZE);
    int y0     = rng_below(rng, L_SIZE);
    double rho = 1.5 + 2.5 * rng_u01(rng);
    int a      = rng_below(rng, N_COMP);
    int b      = rng_below(rng, N_COMP - 1);
    if (b >= a) {
        ++b;
    }
    uint64_t rw = rng_u64(rng);
    int sign    = (rw & 1ULL) ? 1 : -1;
    int inv_dir = (rw & 2ULL) ? 1 : -1;

    double s_old = cfg_total_action(g, beta, 1.0);
    apply_bp_su2_twist(prop, g, x0, y0, rho, a, b, sign, inv_dir);
    double s_new = cfg_total_action(prop, beta, 1.0);
    double dS    = s_new - s_old;
    if (dS <= 0.0 || rng_pos_u01(rng) < exp(-dS)) {
        memcpy(g, prop, sizeof(*g));
        return 1;
    }
    return 0;
}

/* --- Startup Mathematical & Physics Self-Tests --- */
static void run_startup_selftests(void) {
    cpn_cfg_t *g   = (cpn_cfg_t *)malloc(sizeof(cpn_cfg_t));
    cpn_cfg_t *gp  = (cpn_cfg_t *)malloc(sizeof(cpn_cfg_t));
    cpn_cfg_t *rot = (cpn_cfg_t *)malloc(sizeof(cpn_cfg_t));
    cpn_cfg_t *vac = (cpn_cfg_t *)malloc(sizeof(cpn_cfg_t));
    cpn_cfg_t *rt  = (cpn_cfg_t *)malloc(sizeof(cpn_cfg_t));
    if (!g || !gp || !rot || !vac || !rt) {
        fprintf(stderr, "ERROR: OOM in run_startup_selftests\n");
        exit(1);
    }
    rng_t rng;
    rng_seed(&rng, 0xC9092026ULL);

    /* Create generic normalized random CP^{N-1} config */
    for (int x = 0; x < L_SIZE; ++x) {
        for (int y = 0; y < L_SIZE; ++y) {
            for (int a = 0; a < N_COMP; ++a) {
                g->re[x][y][a] = 2.0 * rng_u01(&rng) - 1.0;
                g->im[x][y][a] = 2.0 * rng_u01(&rng) - 1.0;
            }
            normalize_site(g->re[x][y], g->im[x][y]);
            double tx = (2.0 * rng_u01(&rng) - 1.0) * K_PI;
            double ty = (2.0 * rng_u01(&rng) - 1.0) * K_PI;
            g->ux_re[x][y] = cos(tx);
            g->ux_im[x][y] = sin(tx);
            g->uy_re[x][y] = cos(ty);
            g->uy_im[x][y] = sin(ty);
        }
    }

    /* 1. Random local U(1) gauge transform invariance + exact integer Q_geom */
    double lam[L_SIZE][L_SIZE];
    for (int x = 0; x < L_SIZE; ++x) {
        for (int y = 0; y < L_SIZE; ++y) {
            lam[x][y] = (2.0 * rng_u01(&rng) - 1.0) * K_PI;
        }
    }
    for (int x = 0; x < L_SIZE; ++x) {
        int xp = idx_next(x);
        for (int y = 0; y < L_SIZE; ++y) {
            int yp    = idx_next(y);
            double c  = cos(lam[x][y]), s = sin(lam[x][y]);
            for (int a = 0; a < N_COMP; ++a) {
                gp->re[x][y][a] = c * g->re[x][y][a] - s * g->im[x][y][a];
                gp->im[x][y][a] = c * g->im[x][y][a] + s * g->re[x][y][a];
            }
            double dx = wrap_pi(atan2(g->ux_im[x][y], g->ux_re[x][y]) + lam[x][y] - lam[xp][y]);
            double dy = wrap_pi(atan2(g->uy_im[x][y], g->uy_re[x][y]) + lam[x][y] - lam[x][yp]);
            gp->ux_re[x][y] = cos(dx);
            gp->ux_im[x][y] = sin(dx);
            gp->uy_re[x][y] = cos(dy);
            gp->uy_im[x][y] = sin(dy);
        }
    }

    /* Also verify global SU(N) rotation invariance via Givens mixing */
    memcpy(rot, g, sizeof(*rot));
    double gt = 0.73, cg = cos(gt), sg = sin(gt);
    for (int x = 0; x < L_SIZE; ++x) {
        for (int y = 0; y < L_SIZE; ++y) {
            double ar = rot->re[x][y][1], ai = rot->im[x][y][1];
            double br = rot->re[x][y][5], bi = rot->im[x][y][5];
            rot->re[x][y][1] = cg * ar - sg * br;
            rot->im[x][y][1] = cg * ai - sg * bi;
            rot->re[x][y][5] = sg * ar + cg * br;
            rot->im[x][y][5] = sg * ai + cg * bi;
        }
    }

    double sum_r0, sum_rp, sum_rr;
    long long q0 = cfg_topo_charge_raw(g, &sum_r0);
    long long qp = cfg_topo_charge_raw(gp, &sum_rp);
    long long qr = cfg_topo_charge_raw(rot, &sum_rr);
    if (fabs(sum_r0 / K_TWO_PI - (double)q0) > 1e-10) {
        fprintf(stderr, "SELFTEST FAIL: Q_geom not integer! %.12f\n", sum_r0 / K_TWO_PI);
        exit(1);
    }
    double a0 = cfg_total_action(g, 0.8, 0.37);
    double ap = cfg_total_action(gp, 0.8, 0.37);
    double ar = cfg_total_action(rot, 0.8, 0.37);
    double g0_a, gp1_a, g0_b, gp1_b, g0_c, gp1_c;
    cfg_projector_fourier(g, &g0_a, &gp1_a);
    cfg_projector_fourier(gp, &g0_b, &gp1_b);
    cfg_projector_fourier(rot, &g0_c, &gp1_c);
    if (q0 != qp || q0 != qr || fabs(a0 - ap) > 1e-9 || fabs(a0 - ar) > 1e-9 ||
        fabs(g0_a - g0_b) > 1e-10 || fabs(g0_a - g0_c) > 1e-10 || fabs(gp1_a - gp1_c) > 1e-10) {
        fprintf(stderr, "SELFTEST FAIL: gauge/SU(N) invariance mismatch\n");
        exit(1);
    }

    /* 2. Microcanonical Over-Relaxation exactness (dS = 0, |z| = 1) */
    memcpy(gp, g, sizeof(*gp));
    for (int par = 0; par < 2; ++par) {
        for (int x = 0; x < L_SIZE; ++x) {
            for (int y = 0; y < L_SIZE; ++y) {
                if (((x + y) & 1) != par) {
                    continue;
                }
                double fr[N_COMP], fi[N_COMP];
                compute_site_staple(gp, x, y, 0.37, fr, fi);
                double fn2 = 0.0, zf = 0.0;
                for (int a = 0; a < N_COMP; ++a) {
                    fn2 += fr[a] * fr[a] + fi[a] * fi[a];
                    zf += gp->re[x][y][a] * fr[a] + gp->im[x][y][a] * fi[a];
                }
                double sc = 2.0 * zf / fn2;
                for (int a = 0; a < N_COMP; ++a) {
                    gp->re[x][y][a] = sc * fr[a] - gp->re[x][y][a];
                    gp->im[x][y][a] = sc * fi[a] - gp->im[x][y][a];
                }
                double n2 = 0.0;
                for (int a = 0; a < N_COMP; ++a) {
                    n2 += gp->re[x][y][a] * gp->re[x][y][a] + gp->im[x][y][a] * gp->im[x][y][a];
                }
                if (fabs(n2 - 1.0) > 1e-12) {
                    fprintf(stderr, "SELFTEST FAIL: OR norm drift %.3e\n", n2 - 1.0);
                    exit(1);
                }
            }
        }
    }
    double a_or = cfg_total_action(gp, 0.8, 0.37);
    if (fabs(a_or - a0) > 1e-9) {
        fprintf(stderr, "SELFTEST FAIL: OR dS=%.3e != 0\n", a_or - a0);
        exit(1);
    }

    /* 3. Fast defect @swap deltaS exactness */
    double c1 = 0.65, c2 = 0.35;
    double sa_1 = cfg_total_action(g, 0.8, c1);
    double sb_2 = cfg_total_action(rot, 0.8, c2);
    double sa_2 = cfg_total_action(g, 0.8, c2);
    double sb_1 = cfg_total_action(rot, 0.8, c1);
    double sd_1 = cfg_defect_sum_re(g);
    double sd_2 = cfg_defect_sum_re(rot);
    double ds_exact = (sa_2 + sb_1) - (sa_1 + sb_2);
    double ds_fast  = 2.0 * (double)N_COMP * 0.8 * (c1 - c2) * (sd_1 - sd_2);
    if (fabs(ds_exact - ds_fast) > 1e-9) {
        fprintf(stderr, "SELFTEST FAIL: swap dS mismatch %.3e vs %.3e\n", ds_exact, ds_fast);
        exit(1);
    }

    /* 4. Best-Fisher von Mises kappa=3.0 check vs analytical I1(3)/I0(3) */
    double sum_cos_vm = 0.0, sum_cos2_vm = 0.0;
    const int n_vm = 100000;
    for (int i = 0; i < n_vm; ++i) {
        double cv = cos(sample_von_mises(3.0, &rng));
        sum_cos_vm += cv;
        sum_cos2_vm += cv * cv;
    }
    double m_vm   = sum_cos_vm / (double)n_vm;
    double var_vm = sum_cos2_vm / (double)n_vm - m_vm * m_vm;
    double err_vm = sqrt(var_vm / (double)n_vm);
    /* Exact I_1(3.0) / I_0(3.0) from power series: */
    double i0 = 1.0, i1 = 1.5, t0 = 1.0, t1 = 1.5;
    for (int k = 1; k < 40; ++k) {
        t0 *= (0.25 * 9.0) / (double)(k * k);
        t1 *= (0.25 * 9.0) / (double)(k * (k + 1));
        i0 += t0;
        i1 += t1;
    }
    double vm_exact = i1 / i0;
    if (fabs(m_vm - vm_exact) > 4.0 * err_vm) {
        fprintf(stderr, "SELFTEST FAIL: von Mises mean %.6f vs exact %.6f\n", m_vm, vm_exact);
        exit(1);
    }

    /* 5. Belavin-Polyakov SU(2) instanton on vacuum: |z|=1, Q_geom=+1, Bogomol'nyi bound S >= 2*pi*N*beta + exact roundtrip */
    cfg_init_ordered(vac);
    apply_bp_su2_twist(gp, vac, L_SIZE / 2, L_SIZE / 2, 2.5, 0, 1, 1, 1);
    apply_bp_su2_twist(rt, gp, L_SIZE / 2, L_SIZE / 2, 2.5, 0, 1, 1, -1);
    long long q_bp = cfg_topo_charge_raw(gp, NULL);
    double s_bp    = cfg_total_action(gp, 0.8, 1.0) - cfg_total_action(vac, 0.8, 1.0);
    double bog_bnd = K_TWO_PI * (double)N_COMP * 0.8;
    double max_rt_err = 0.0;
    for (int x = 0; x < L_SIZE; ++x) {
        for (int y = 0; y < L_SIZE; ++y) {
            double n2 = 0.0;
            for (int a = 0; a < N_COMP; ++a) {
                n2 += gp->re[x][y][a] * gp->re[x][y][a] + gp->im[x][y][a] * gp->im[x][y][a];
                double er = fabs(rt->re[x][y][a] - vac->re[x][y][a]);
                double ei = fabs(rt->im[x][y][a] - vac->im[x][y][a]);
                if (er > max_rt_err) {
                    max_rt_err = er;
                }
                if (ei > max_rt_err) {
                    max_rt_err = ei;
                }
            }
            double exr = fabs(rt->ux_re[x][y] - vac->ux_re[x][y]);
            double exi = fabs(rt->ux_im[x][y] - vac->ux_im[x][y]);
            double eyr = fabs(rt->uy_re[x][y] - vac->uy_re[x][y]);
            double eyi = fabs(rt->uy_im[x][y] - vac->uy_im[x][y]);
            if (exr > max_rt_err) max_rt_err = exr;
            if (exi > max_rt_err) max_rt_err = exi;
            if (eyr > max_rt_err) max_rt_err = eyr;
            if (eyi > max_rt_err) max_rt_err = eyi;
            if (fabs(n2 - 1.0) > 1e-12) {
                fprintf(stderr, "SELFTEST FAIL: BP norm %.3e\n", n2);
                exit(1);
            }
        }
    }
    if (q_bp != 1 || s_bp + 1e-6 < bog_bnd || max_rt_err > 1e-12) {
        fprintf(stderr, "SELFTEST FAIL: BP Q=%lld S_exc=%.4f bogomolnyi=%.4f rt=%.3e\n",
                q_bp, s_bp, bog_bnd, max_rt_err);
        exit(1);
    }
    printf("  [SELFTEST] PASS (gauge/SU(%d) <1e-9, OR <1e-9, fast swap <1e-9, vM exact, BP Q=+1 S=%.2f >= %.2f, rt=%.2e)\n",
           N_COMP, s_bp, bog_bnd, max_rt_err);

    free(g);
    free(gp);
    free(rot);
    free(vac);
    free(rt);
}

/* --- Madras-Sokal Autocorrelation + Jackknife Correlation Length --- */
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

static void compute_xi_jackknife(const double *g0_h, const double *gp1_h, int n,
                                 double *out_xi, double *out_err_xi) {
    double sum_g0 = 0.0, sum_gp1 = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_g0 += g0_h[i];
        sum_gp1 += gp1_h[i];
    }
    const double inv_2sin = 1.0 / (2.0 * sin(K_PI / (double)L_SIZE));
    double r_all          = sum_g0 / sum_gp1;
    double xi_full        = (r_all > 1.0) ? (inv_2sin * sqrt(r_all - 1.0)) : 0.0;

    int bsz        = n / JK_BLOCKS;
    double jk[JK_BLOCKS];
    double jk_mean = 0.0;
    for (int b = 0; b < JK_BLOCKS; ++b) {
        double bg0 = 0.0, bgp1 = 0.0;
        int lo = b * bsz;
        int hi = (b == JK_BLOCKS - 1) ? n : ((b + 1) * bsz);
        for (int i = lo; i < hi; ++i) {
            bg0 += g0_h[i];
            bgp1 += gp1_h[i];
        }
        double r_jk = (sum_g0 - bg0) / (sum_gp1 - bgp1);
        jk[b]       = (r_jk > 1.0) ? (inv_2sin * sqrt(r_jk - 1.0)) : 0.0;
        jk_mean += jk[b];
    }
    jk_mean /= (double)JK_BLOCKS;
    double jk_var = 0.0;
    for (int b = 0; b < JK_BLOCKS; ++b) {
        double d = jk[b] - jk_mean;
        jk_var += d * d;
    }
    jk_var *= (double)(JK_BLOCKS - 1) / (double)JK_BLOCKS;
    *out_xi     = xi_full;
    *out_err_xi = sqrt(jk_var);
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
    c_new[0]               = 1.0;
    c_new[R_REPLICAS - 1]  = 0.0;
    for (int r = 1; r < R_REPLICAS - 1; ++r) {
        double s_tgt = ((double)r / (double)(R_REPLICAS - 1)) * total_s;
        int k        = 0;
        while (k + 1 < R_REPLICAS - 1 && cum_s[k + 1] < s_tgt) {
            ++k;
        }
        double denom = cum_s[k + 1] - cum_s[k];
        double t     = (denom > 1e-12) ? ((s_tgt - cum_s[k]) / denom) : 0.5;
        double interp = cum_c[k] + t * (cum_c[k + 1] - cum_c[k]);
        c_new[r]     = 0.5 * c[r] + 0.5 * interp; /* damped geodesic update */
    }
    /* Enforce strict monotonicity & boundary gaps */
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

typedef enum {
    CPN_LOCAL_OR_METRO   = 0,
    CPN_BARE_BP_INST     = 1,
    CPN_SYMBIOTIC_PTBC   = 2
} cpn_solver_t;

typedef struct {
    const char *name;
    double beta;
    int n_replicas;
    int work_sweeps;
    double mean_e;
    double err_e;
    double xi_g;
    double err_xi_g;
    double chi_t_geom;
    double err_chi_t_geom;
    double tau_work_q_geom;
    int frozen_q_geom;
    long long trans_q_geom;
    double chi_t_cool;
    double err_chi_t_cool;
    double tau_work_q_cool;
    int frozen_q_cool;
    long long trans_q_cool;
    double chi_xi2_cool;
    double bare_inst_acc;
    double min_swap_acc;
    double mean_swap_acc;
    double c_ladder[R_REPLICAS];
} cpn_res_t;

static cpn_res_t run_single_replica_solver(cpn_solver_t kind, const char *name,
                                           double beta, uint64_t seed) {
    cpn_cfg_t *g    = (cpn_cfg_t *)malloc(sizeof(cpn_cfg_t));
    cpn_cfg_t *prop = (cpn_cfg_t *)malloc(sizeof(cpn_cfg_t));
    cpn_cfg_t *wrk  = (cpn_cfg_t *)malloc(sizeof(cpn_cfg_t));
    if (!g || !prop || !wrk) {
        fprintf(stderr, "ERROR: OOM allocating cfg in run_single_replica_solver\n");
        exit(1);
    }
    cfg_init_ordered(g);
    rng_t rng;
    rng_seed(&rng, seed);

    const int n_cool_samples = N_EQ_MEAS / R_REPLICAS; /* spaced by R_REPLICAS sweeps = exact equal work */
    double *e_h     = (double *)malloc((size_t)N_EQ_MEAS * sizeof(double));
    double *g0_h    = (double *)malloc((size_t)n_cool_samples * sizeof(double));
    double *gp1_h   = (double *)malloc((size_t)n_cool_samples * sizeof(double));
    double *qg_h    = (double *)malloc((size_t)N_EQ_MEAS * sizeof(double));
    double *qg2_h   = (double *)malloc((size_t)N_EQ_MEAS * sizeof(double));
    double *qc_h    = (double *)malloc((size_t)n_cool_samples * sizeof(double));
    double *qc2_h   = (double *)malloc((size_t)n_cool_samples * sizeof(double));
    if (!e_h || !g0_h || !gp1_h || !qg_h || !qg2_h || !qc_h || !qc2_h) {
        fprintf(stderr, "ERROR: OOM in run_single_replica_solver\n");
        exit(1);
    }

    double step_eps = 0.8;
    for (int t = 0; t < N_THERM_1R; ++t) {
        double acc = 0.0;
        sweep_cpn_local(g, beta, 1.0, step_eps, &rng, &acc);
        if (kind == CPN_BARE_BP_INST) {
            propose_bare_bp_instanton(g, prop, beta, &rng);
        }
        if (t < N_THERM_1R / 2 && (t % 25) == 24) {
            if (acc > 0.56) {
                step_eps *= 1.08;
                if (step_eps > K_PI) {
                    step_eps = K_PI;
                }
            } else if (acc < 0.42) {
                step_eps *= 0.92;
            }
        }
    }

    /* Frozen step_eps measurement phase */
    long long inst_acc   = 0;
    long long trans_qg   = 0;
    long long trans_qc   = 0;
    long long prev_qg    = 0;
    long long prev_qc    = 0;
    int cool_idx         = 0;
    for (int t = 0; t < N_EQ_MEAS; ++t) {
        sweep_cpn_local(g, beta, 1.0, step_eps, &rng, NULL);
        if (kind == CPN_BARE_BP_INST) {
            inst_acc += propose_bare_bp_instanton(g, prop, beta, &rng);
        }
        e_h[t]        = cfg_mean_link_energy(g);
        long long qg  = cfg_topo_charge_raw(g, NULL);
        qg_h[t]       = (double)qg;
        qg2_h[t]      = (double)(qg * qg) / (double)V_SITES;
        if (t > 0 && qg != prev_qg) {
            ++trans_qg;
        }
        prev_qg = qg;

        if ((t % R_REPLICAS) == (R_REPLICAS - 1)) {
            cfg_projector_fourier(g, &g0_h[cool_idx], &gp1_h[cool_idx]);
            long long qc    = cfg_topo_charge_cooled(g, wrk);
            qc_h[cool_idx]  = (double)qc;
            qc2_h[cool_idx] = (double)(qc * qc) / (double)V_SITES;
            if (cool_idx > 0 && qc != prev_qc) {
                ++trans_qc;
            }
            prev_qc = qc;
            ++cool_idx;
        }
    }

    double em, ev, etau, eerr;
    double qgm, qgv, qgtau, qgerr;
    double qg2m, qg2v, qg2tau, qg2err;
    double qcm, qcv, qctau, qcerr;
    double qc2m, qc2v, qc2tau, qc2err;
    double xi, exi;
    int efrz, qgfrz, qg2frz, qcfrz, qc2frz;

    compute_madras_sokal(e_h, N_EQ_MEAS, &em, &ev, &etau, &eerr, &efrz);
    compute_madras_sokal(qg_h, N_EQ_MEAS, &qgm, &qgv, &qgtau, &qgerr, &qgfrz);
    compute_madras_sokal(qg2_h, N_EQ_MEAS, &qg2m, &qg2v, &qg2tau, &qg2err, &qg2frz);
    compute_madras_sokal(qc_h, n_cool_samples, &qcm, &qcv, &qctau, &qcerr, &qcfrz);
    compute_madras_sokal(qc2_h, n_cool_samples, &qc2m, &qc2v, &qc2tau, &qc2err, &qc2frz);
    compute_xi_jackknife(g0_h, gp1_h, n_cool_samples, &xi, &exi);

    free(g);
    free(prop);
    free(wrk);
    free(e_h);
    free(g0_h);
    free(gp1_h);
    free(qg_h);
    free(qg2_h);
    free(qc_h);
    free(qc2_h);

    cpn_res_t res;
    memset(&res, 0, sizeof(res));
    res.name            = name;
    res.beta            = beta;
    res.n_replicas      = 1;
    res.work_sweeps     = N_EQ_MEAS;
    res.mean_e          = em;
    res.err_e           = eerr;
    res.xi_g            = xi;
    res.err_xi_g        = exi;
    res.chi_t_geom      = qg2m;
    res.err_chi_t_geom  = qg2err;
    res.tau_work_q_geom = qgfrz ? -1.0 : qgtau;
    res.frozen_q_geom   = qgfrz;
    res.trans_q_geom    = trans_qg;
    res.chi_t_cool      = qc2m;
    res.err_chi_t_cool  = qc2err;
    res.tau_work_q_cool = qcfrz ? -1.0 : ((double)R_REPLICAS * qctau);
    res.frozen_q_cool   = qcfrz;
    res.trans_q_cool    = trans_qc;
    res.chi_xi2_cool    = qc2m * xi * xi;
    res.bare_inst_acc   = (kind == CPN_BARE_BP_INST) ? ((double)inst_acc / (double)N_EQ_MEAS) : 0.0;
    res.c_ladder[0]     = 1.0;
    return res;
}

static cpn_res_t run_ptbc_solver(double beta, uint64_t seed) {
    cpn_cfg_t *rep[R_REPLICAS];
    for (int r = 0; r < R_REPLICAS; ++r) {
        rep[r] = (cpn_cfg_t *)malloc(sizeof(cpn_cfg_t));
        if (!rep[r]) {
            fprintf(stderr, "ERROR: OOM allocating replica %d\n", r);
            exit(1);
        }
        cfg_init_ordered(rep[r]);
    }
    cpn_cfg_t *tmp = (cpn_cfg_t *)malloc(sizeof(cpn_cfg_t));
    cpn_cfg_t *wrk = (cpn_cfg_t *)malloc(sizeof(cpn_cfg_t));
    if (!tmp || !wrk) {
        fprintf(stderr, "ERROR: OOM allocating scratch in run_ptbc_solver\n");
        exit(1);
    }
    rng_t rng;
    rng_seed(&rng, seed);

    double c[R_REPLICAS];
    double step_eps[R_REPLICAS];
    for (int r = 0; r < R_REPLICAS; ++r) {
        c[r]        = 1.0 - (double)r / (double)(R_REPLICAS - 1);
        step_eps[r] = 0.8;
    }

    double *e_h   = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *g0_h  = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *gp1_h = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *qg_h  = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *qg2_h = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *qc_h  = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    double *qc2_h = (double *)malloc((size_t)N_MEAS_PT * sizeof(double));
    if (!e_h || !g0_h || !gp1_h || !qg_h || !qg2_h || !qc_h || !qc2_h) {
        fprintf(stderr, "ERROR: OOM in run_ptbc_solver\n");
        exit(1);
    }

    double swap_prob_sum[R_REPLICAS - 1];
    long long swap_prob_cnt[R_REPLICAS - 1];
    memset(swap_prob_sum, 0, sizeof(swap_prob_sum));
    memset(swap_prob_cnt, 0, sizeof(swap_prob_cnt));

    const double coeff = 2.0 * (double)N_COMP * beta;

    /* Thermalization: adapt step_eps & Fisher 2-Jet ladder during first half, then FREEZE cleanly! */
    for (int t = 0; t < N_THERM_PT; ++t) {
        for (int r = 0; r < R_REPLICAS; ++r) {
            double acc = 0.0;
            sweep_cpn_local(rep[r], beta, c[r], step_eps[r], &rng, &acc);
            if (t < N_THERM_PT / 2 && (t % 25) == 24) {
                if (acc > 0.56) {
                    step_eps[r] *= 1.08;
                    if (step_eps[r] > K_PI) {
                        step_eps[r] = K_PI;
                    }
                } else if (acc < 0.42) {
                    step_eps[r] *= 0.92;
                }
            }
        }
        /* Symmetric even+odd adjacent @swap rounds (Theorem P3B) via zero-copy pointer swap */
        for (int pass = 0; pass < 2; ++pass) {
            int parity = (t + pass) & 1;
            for (int r = parity; r + 1 < R_REPLICAS; r += 2) {
                double sd_a = cfg_defect_sum_re(rep[r]);
                double sd_b = cfg_defect_sum_re(rep[r + 1]);
                double dS   = coeff * (c[r] - c[r + 1]) * (sd_a - sd_b);
                double p    = (dS <= 0.0) ? 1.0 : exp(-dS);
                if (t < N_THERM_PT / 2) {
                    swap_prob_sum[r] += p;
                    ++swap_prob_cnt[r];
                }
                if (dS <= 0.0 || rng_pos_u01(&rng) < p) {
                    cpn_cfg_t *sw = rep[r];
                    rep[r]        = rep[r + 1];
                    rep[r + 1]    = sw;
                }
            }
        }
        shift_torus_random(rep[0], tmp, &rng);

        /* Fisher 2-Jet geodesic ladder update every 50 sweeps during first half */
        if (t < N_THERM_PT / 2 && (t % 50) == 49) {
            double avg_p[R_REPLICAS - 1];
            for (int r = 0; r < R_REPLICAS - 1; ++r) {
                avg_p[r]         = (swap_prob_cnt[r] > 0) ? (swap_prob_sum[r] / (double)swap_prob_cnt[r]) : 0.5;
                swap_prob_sum[r] = 0.0;
                swap_prob_cnt[r] = 0;
            }
            adapt_fisher_ladder(c, avg_p);
        }
    }

    /* Frozen measurement phase (exact detailed balance!) */
    long long swap_acc[R_REPLICAS - 1];
    long long swap_try[R_REPLICAS - 1];
    memset(swap_acc, 0, sizeof(swap_acc));
    memset(swap_try, 0, sizeof(swap_try));
    long long trans_qg = 0, trans_qc = 0;
    long long prev_qg  = 0, prev_qc  = 0;

    for (int t = 0; t < N_MEAS_PT; ++t) {
        for (int r = 0; r < R_REPLICAS; ++r) {
            sweep_cpn_local(rep[r], beta, c[r], step_eps[r], &rng, NULL);
        }
        for (int pass = 0; pass < 2; ++pass) {
            int parity = (t + pass) & 1;
            for (int r = parity; r + 1 < R_REPLICAS; r += 2) {
                double sd_a = cfg_defect_sum_re(rep[r]);
                double sd_b = cfg_defect_sum_re(rep[r + 1]);
                double dS   = coeff * (c[r] - c[r + 1]) * (sd_a - sd_b);
                ++swap_try[r];
                if (dS <= 0.0 || rng_pos_u01(&rng) < exp(-dS)) {
                    cpn_cfg_t *sw = rep[r];
                    rep[r]        = rep[r + 1];
                    rep[r + 1]    = sw;
                    ++swap_acc[r];
                }
            }
        }
        shift_torus_random(rep[0], tmp, &rng);

        e_h[t] = cfg_mean_link_energy(rep[0]);
        cfg_projector_fourier(rep[0], &g0_h[t], &gp1_h[t]);
        long long qg = cfg_topo_charge_raw(rep[0], NULL);
        long long qc = cfg_topo_charge_cooled(rep[0], wrk);
        qg_h[t]      = (double)qg;
        qg2_h[t]     = (double)(qg * qg) / (double)V_SITES;
        qc_h[t]      = (double)qc;
        qc2_h[t]     = (double)(qc * qc) / (double)V_SITES;
        if (t > 0 && qg != prev_qg) {
            ++trans_qg;
        }
        if (t > 0 && qc != prev_qc) {
            ++trans_qc;
        }
        prev_qg = qg;
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

    double em, ev, etau, eerr;
    double qgm, qgv, qgtau, qgerr;
    double qg2m, qg2v, qg2tau, qg2err;
    double qcm, qcv, qctau, qcerr;
    double qc2m, qc2v, qc2tau, qc2err;
    double xi, exi;
    int efrz, qgfrz, qg2frz, qcfrz, qc2frz;

    compute_madras_sokal(e_h, N_MEAS_PT, &em, &ev, &etau, &eerr, &efrz);
    compute_madras_sokal(qg_h, N_MEAS_PT, &qgm, &qgv, &qgtau, &qgerr, &qgfrz);
    compute_madras_sokal(qg2_h, N_MEAS_PT, &qg2m, &qg2v, &qg2tau, &qg2err, &qg2frz);
    compute_madras_sokal(qc_h, N_MEAS_PT, &qcm, &qcv, &qctau, &qcerr, &qcfrz);
    compute_madras_sokal(qc2_h, N_MEAS_PT, &qc2m, &qc2v, &qc2tau, &qc2err, &qc2frz);
    compute_xi_jackknife(g0_h, gp1_h, N_MEAS_PT, &xi, &exi);

    for (int r = 0; r < R_REPLICAS; ++r) {
        free(rep[r]);
    }
    free(tmp);
    free(wrk);
    free(e_h);
    free(g0_h);
    free(gp1_h);
    free(qg_h);
    free(qg2_h);
    free(qc_h);
    free(qc2_h);

    cpn_res_t res;
    memset(&res, 0, sizeof(res));
    res.name            = "symbiotic_ptbc_2jet";
    res.beta            = beta;
    res.n_replicas      = R_REPLICAS;
    res.work_sweeps     = R_REPLICAS * N_MEAS_PT;
    res.mean_e          = em;
    res.err_e           = eerr;
    res.xi_g            = xi;
    res.err_xi_g        = exi;
    res.chi_t_geom      = qg2m;
    res.err_chi_t_geom  = qg2err;
    res.tau_work_q_geom = qgfrz ? -1.0 : ((double)R_REPLICAS * qgtau);
    res.frozen_q_geom   = qgfrz;
    res.trans_q_geom    = trans_qg;
    res.chi_t_cool      = qc2m;
    res.err_chi_t_cool  = qc2err;
    res.tau_work_q_cool = qcfrz ? -1.0 : ((double)R_REPLICAS * qctau);
    res.frozen_q_cool   = qcfrz;
    res.trans_q_cool    = trans_qc;
    res.chi_xi2_cool    = qc2m * xi * xi;
    res.min_swap_acc    = min_sw;
    res.mean_swap_acc   = sum_sw / (double)(R_REPLICAS - 1);
    for (int r = 0; r < R_REPLICAS; ++r) {
        res.c_ladder[r] = c[r];
    }
    return res;
}

typedef struct {
    int b_idx;
    double beta;
    cpn_res_t res[3];
} beta_job_t;

static void *beta_worker(void *arg) {
    beta_job_t *job = (beta_job_t *)arg;
    uint64_t base   = 0xC90920260000ULL + (uint64_t)job->b_idx * 1000ULL;
    job->res[0]     = run_single_replica_solver(CPN_LOCAL_OR_METRO, "local_or_metro", job->beta, base + 1ULL);
    job->res[1]     = run_single_replica_solver(CPN_BARE_BP_INST, "bare_bp_inst", job->beta, base + 2ULL);
    job->res[2]     = run_ptbc_solver(job->beta, base + 3ULL);
    return NULL;
}

int main(int argc, char **argv) {
    const char *json_path = NULL;
    for (int i = 1; i + 1 < argc; ++i) {
        if (strcmp(argv[i], "--json") == 0) {
            json_path = argv[i + 1];
        }
    }

    printf("=== 2D CP^%d (N=%d) Equal-Work Wind-Tunnel (L=%d, V=%d, eq_work=%d sweeps/beta) ===\n",
           N_COMP - 1, N_COMP, L_SIZE, V_SITES, N_EQ_MEAS);
    run_startup_selftests();

    const double betas[N_BETAS] = {0.60, 0.65, 0.70, 0.73};
    beta_job_t jobs[N_BETAS];
    pthread_t tids[N_BETAS];

    for (int b = 0; b < N_BETAS; ++b) {
        jobs[b].b_idx = b;
        jobs[b].beta  = betas[b];
        if (pthread_create(&tids[b], NULL, beta_worker, &jobs[b]) != 0) {
            fprintf(stderr, "ERROR: pthread_create failed\n");
            return 1;
        }
    }
    for (int b = 0; b < N_BETAS; ++b) {
        pthread_join(tids[b], NULL);
    }

    for (int b = 0; b < N_BETAS; ++b) {
        printf("\n  --- beta = %.2f ---\n", jobs[b].beta);
        for (int s = 0; s < 3; ++s) {
            const cpn_res_t *r = &jobs[b].res[s];
            printf("    %-20s | E=%.5f+-%.5f | xi=%.3f+-%.3f | chi_c=%.6f+-%.6f | chi*xi2=%.4f | "
                   "tauWc=%s | tc=%lld | acc_sw=%.3f (min=%.3f) | bare_acc=%.4f\n",
                   r->name, r->mean_e, r->err_e, r->xi_g, r->err_xi_g,
                   r->chi_t_cool, r->err_chi_t_cool, r->chi_xi2_cool,
                   r->frozen_q_cool ? "FROZEN" : "", r->trans_q_cool,
                   r->mean_swap_acc, r->min_swap_acc, r->bare_inst_acc);
            if (!r->frozen_q_cool) {
                printf("      -> tau_work_q_cool=%.2f, tau_work_q_geom=%.2f (trans_geom=%lld)\n",
                       r->tau_work_q_cool, r->tau_work_q_geom, r->trans_q_geom);
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
        fprintf(fp, "  \"N\": %d,\n", N_COMP);
        fprintf(fp, "  \"L\": %d,\n", L_SIZE);
        fprintf(fp, "  \"V\": %d,\n", V_SITES);
        fprintf(fp, "  \"L_defect\": %d,\n", L_DEFECT);
        fprintf(fp, "  \"R_replicas\": %d,\n", R_REPLICAS);
        fprintf(fp, "  \"eq_work_sweeps\": %d,\n", N_EQ_MEAS);
        fprintf(fp, "  \"runs\": [\n");
        int idx   = 0;
        int total = N_BETAS * 3;
        for (int b = 0; b < N_BETAS; ++b) {
            for (int s = 0; s < 3; ++s) {
                const cpn_res_t *r = &jobs[b].res[s];
                ++idx;
                fprintf(fp, "    {\n");
                fprintf(fp, "      \"beta\": %.2f,\n", r->beta);
                fprintf(fp, "      \"solver\": \"%s\",\n", r->name);
                fprintf(fp, "      \"n_replicas\": %d,\n", r->n_replicas);
                fprintf(fp, "      \"work_sweeps\": %d,\n", r->work_sweeps);
                fprintf(fp, "      \"mean_e\": %.8f,\n", r->mean_e);
                fprintf(fp, "      \"err_e\": %.8f,\n", r->err_e);
                fprintf(fp, "      \"xi_g\": %.6f,\n", r->xi_g);
                fprintf(fp, "      \"err_xi_g\": %.6f,\n", r->err_xi_g);
                fprintf(fp, "      \"chi_t_geom\": %.8f,\n", r->chi_t_geom);
                fprintf(fp, "      \"err_chi_t_geom\": %.8f,\n", r->err_chi_t_geom);
                fprintf(fp, "      \"tau_work_q_geom\": %.4f,\n", r->tau_work_q_geom);
                fprintf(fp, "      \"frozen_q_geom\": %s,\n", r->frozen_q_geom ? "true" : "false");
                fprintf(fp, "      \"trans_q_geom\": %lld,\n", r->trans_q_geom);
                fprintf(fp, "      \"chi_t_cool\": %.8f,\n", r->chi_t_cool);
                fprintf(fp, "      \"err_chi_t_cool\": %.8f,\n", r->err_chi_t_cool);
                fprintf(fp, "      \"tau_work_q_cool\": %.4f,\n", r->tau_work_q_cool);
                fprintf(fp, "      \"frozen_q_cool\": %s,\n", r->frozen_q_cool ? "true" : "false");
                fprintf(fp, "      \"trans_q_cool\": %lld,\n", r->trans_q_cool);
                fprintf(fp, "      \"chi_xi2_cool\": %.6f,\n", r->chi_xi2_cool);
                fprintf(fp, "      \"bare_inst_acc\": %.6f,\n", r->bare_inst_acc);
                fprintf(fp, "      \"min_swap_acc\": %.6f,\n", r->min_swap_acc);
                fprintf(fp, "      \"mean_swap_acc\": %.6f,\n", r->mean_swap_acc);
                fprintf(fp, "      \"c_ladder\": [");
                for (int k = 0; k < r->n_replicas; ++k) {
                    fprintf(fp, "%.4f%s", r->c_ladder[k], (k + 1 < r->n_replicas) ? ", " : "");
                }
                fprintf(fp, "]\n");
                fprintf(fp, "    }%s\n", (idx < total) ? "," : "");
            }
        }
        fprintf(fp, "  ]\n");
        fprintf(fp, "}\n");
        fclose(fp);
    }

    return 0;
}
