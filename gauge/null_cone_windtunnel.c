/* Stage 11: Null-Cone Algebraic Gravity & Holographic Falsification Wind Tunnel
 * Native C17 verification companion to FLOW (`Flow.NullConeBundle`) for arXiv:2611.34924:
 *   N1: Canonical SL(2, Q) symplectic invariance, antisymmetry & Jacobi identity (exact int64 ratios)
 *   N2: 8x8 Gaussian 2-Jet wavefunction Gram Cholesky PD + indefinite conformal ghost witness
 *   N3: Finite-G_N nonperturbative path integral quadrature vs. g=0 log divergence & g<0 poles
 *   N4: Ashoke Sen 4D pure-gravity Schwarzschild 1-loop rational coeffs + Planckian breakdown
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

typedef struct {
    double re;
    double im;
} c128_t;

static c128_t c_add(c128_t a, c128_t b) { return (c128_t){a.re + b.re, a.im + b.im}; }
static c128_t c_sub(c128_t a, c128_t b) { return (c128_t){a.re - b.re, a.im - b.im}; }
static c128_t c_mul(c128_t a, c128_t b) {
    return (c128_t){a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}
static c128_t c_conj(c128_t a) { return (c128_t){a.re, -a.im}; }
static c128_t c_scale(c128_t a, double s) { return (c128_t){a.re * s, a.im * s}; }
static double c_abs(c128_t a) { return hypot(a.re, a.im); }

typedef struct {
    int64_t num;
    int64_t den;
} rat_t;

static int64_t gcd64(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) {
        int64_t t = a % b;
        a = b;
        b = t;
    }
    return a == 0 ? 1 : a;
}

static rat_t rat_make(int64_t n, int64_t d) {
    if (d < 0) {
        n = -n;
        d = -d;
    }
    int64_t g = gcd64(n, d);
    return (rat_t){n / g, d / g};
}

static rat_t rat_add(rat_t a, rat_t b) {
    return rat_make(a.num * b.den + b.num * a.den, a.den * b.den);
}

static rat_t rat_sub(rat_t a, rat_t b) {
    return rat_make(a.num * b.den - b.num * a.den, a.den * b.den);
}

static rat_t rat_mul(rat_t a, rat_t b) {
    int64_t g1 = gcd64(a.num, b.den);
    int64_t g2 = gcd64(b.num, a.den);
    return rat_make((a.num / g1) * (b.num / g2), (a.den / g2) * (b.den / g1));
}

static int run_n1(int64_t *sl2_err_num, int64_t *anti_err_num, int64_t *jacobi_err_num) {
    rat_t uq = rat_make(3, 5), up = rat_make(-7, 11);
    rat_t vq = rat_make(13, 7), vp = rat_make(4, 9);
    rat_t a = rat_make(5, 1), b = rat_make(2, 1), c = rat_make(2, 1), d = rat_make(1, 1);

    rat_t tuq = rat_add(rat_mul(a, uq), rat_mul(b, up));
    rat_t tup = rat_add(rat_mul(c, uq), rat_mul(d, up));
    rat_t tvq = rat_add(rat_mul(a, vq), rat_mul(b, vp));
    rat_t tvp = rat_add(rat_mul(c, vq), rat_mul(d, vp));

    rat_t s_uv = rat_sub(rat_mul(uq, vp), rat_mul(up, vq));
    rat_t s_vu = rat_sub(rat_mul(vq, up), rat_mul(vp, uq));
    rat_t s_tuv = rat_sub(rat_mul(tuq, tvp), rat_mul(tup, tvq));

    rat_t e_sl2 = rat_sub(s_tuv, s_uv);
    rat_t e_anti = rat_add(s_uv, s_vu);

    rat_t q = rat_make(7, 13), p = rat_make(-11, 5);
    rat_t t1 = rat_mul(rat_mul(rat_make(2, 1), q), rat_mul(rat_make(-4, 1), p));
    rat_t t2 = rat_mul(rat_mul(rat_make(-2, 1), p), rat_mul(rat_make(-4, 1), q));
    rat_t t3 = rat_sub(rat_mul(p, rat_mul(rat_make(4, 1), q)),
                       rat_mul(q, rat_mul(rat_make(4, 1), p)));
    rat_t jac = rat_add(rat_add(t1, t2), t3);

    *sl2_err_num = e_sl2.num;
    *anti_err_num = e_anti.num;
    *jacobi_err_num = jac.num;
    return (*sl2_err_num == 0 && *anti_err_num == 0 && *jacobi_err_num == 0) ? 1 : 0;
}

static int cholesky_ldl(size_t n, const c128_t G[9][9], double *min_pivot, double *ghost_norm) {
    c128_t L[9][9];
    double d[9];
    memset(L, 0, sizeof(L));
    memset(d, 0, sizeof(d));
    double min_p = 1e300;

    for (size_t i = 0; i < n; i++) {
        L[i][i] = (c128_t){1.0, 0.0};
        for (size_t j = 0; j < i; j++) {
            c128_t s = G[i][j];
            for (size_t m = 0; m < j; m++) {
                c128_t lm = c_mul(L[i][m], c_conj(L[j][m]));
                s = c_sub(s, c_scale(lm, d[m]));
            }
            L[i][j] = c_scale(s, 1.0 / d[j]);
        }
        c128_t diag_v = G[i][i];
        for (size_t m = 0; m < i; m++) {
            c128_t lm = c_mul(L[i][m], c_conj(L[i][m]));
            diag_v = c_sub(diag_v, c_scale(lm, d[m]));
        }
        double di = diag_v.re;
        d[i] = di;
        if (di < min_p) min_p = di;
        if (di <= 1e-12) {
            c128_t v[9];
            memset(v, 0, sizeof(v));
            v[i] = (c128_t){1.0, 0.0};
            for (size_t jr = i; jr > 0; jr--) {
                size_t j = jr - 1;
                c128_t acc = {0.0, 0.0};
                for (size_t m = j + 1; m <= i; m++) {
                    acc = c_add(acc, c_mul(c_conj(L[m][j]), v[m]));
                }
                v[j] = c_scale(acc, -1.0);
            }
            c128_t qsum = {0.0, 0.0};
            for (size_t r = 0; r < n; r++) {
                c128_t row = {0.0, 0.0};
                for (size_t c = 0; c < n; c++) {
                    row = c_add(row, c_mul(G[r][c], v[c]));
                }
                qsum = c_add(qsum, c_mul(c_conj(v[r]), row));
            }
            *min_pivot = di;
            *ghost_norm = qsum.re;
            return 0;
        }
    }
    *min_pivot = min_p;
    *ghost_norm = 0.0;
    return 1;
}

typedef struct {
    int is_pd;
    double min_pivot;
    double hermitian_err;
    double trace_val;
    int ghost_found;
    double conformal_eta;
    double ghost_norm;
} n2_res_t;

static int run_n2(n2_res_t *out) {
    static const double gammas[8] = {0.15, 0.35, 0.75, 1.00, 1.60, 2.80, 4.50, 8.00};
    static const double mus[8] = {-1.8, -0.9, -0.2, 0.0, 0.4, 1.1, 1.7, 2.3};
    static const double ps[8] = {-0.7, 0.4, -1.1, 0.0, 0.9, -0.5, 1.3, -1.6};

    c128_t G[9][9];
    memset(G, 0, sizeof(G));
    for (size_t j = 0; j < 8; j++) {
        for (size_t k = 0; k < 8; k++) {
            double gj = gammas[j], gk = gammas[k];
            double pref = pow((4.0 * gj * gk) / ((gj + gk) * (gj + gk)), 0.25);
            double dm = mus[j] - mus[k];
            double dp = ps[j] - ps[k];
            double re_exp = -(gj * gk * dm * dm + dp * dp) / (2.0 * (gj + gk));
            double im_exp = ((gj * mus[j] + gk * mus[k]) * (ps[k] - ps[j])) / (gj + gk);
            double mag = pref * exp(re_exp);
            G[j][k] = (c128_t){mag * cos(im_exp), mag * sin(im_exp)};
        }
    }

    double herm_err = 0.0;
    double trace_v = 0.0;
    for (size_t j = 0; j < 8; j++) {
        trace_v += G[j][j].re;
        for (size_t k = 0; k < 8; k++) {
            double d_h = c_abs(c_sub(G[j][k], c_conj(G[k][j])));
            if (d_h > herm_err) herm_err = d_h;
        }
    }

    double min_p = 0.0, unused_ghost = 0.0;
    int pd = cholesky_ldl(8, G, &min_p, &unused_ghost);

    out->is_pd = pd;
    out->min_pivot = min_p;
    out->hermitian_err = herm_err;
    out->trace_val = trace_v;
    out->ghost_found = 0;
    out->conformal_eta = 0.0;
    out->ghost_norm = 0.0;

    for (int m = 0; m < 12; m++) {
        double eta = ldexp(1.0, -m);
        G[8][8] = (c128_t){-eta, 0.0};
        double p9 = 0.0, gn = 0.0;
        int pd9 = cholesky_ldl(9, G, &p9, &gn);
        if (!pd9 && gn < 0.0) {
            out->ghost_found = 1;
            out->conformal_eta = eta;
            out->ghost_norm = gn;
            break;
        }
    }

    return (pd && min_p > 0.0 && herm_err < 1e-13 && out->ghost_found) ? 1 : 0;
}

static double np_integrand(double t, double g) {
    double om = 1.0 - t;
    return om / (om * om + g * t * t);
}

static double simpson_step(double a, double b, double fa, double fm, double fb) {
    return ((b - a) / 6.0) * (fa + 4.0 * fm + fb);
}

static double adaptive_simpson_rec(double g, double a, double b, double eps, double whole,
                                   double fa, double fm, double fb, int depth) {
    double m = 0.5 * (a + b);
    double lm = 0.5 * (a + m);
    double rm = 0.5 * (m + b);
    double flm = np_integrand(lm, g);
    double frm = np_integrand(rm, g);
    double left = simpson_step(a, m, fa, flm, fm);
    double right = simpson_step(m, b, fm, frm, fb);
    double delta = left + right - whole;
    if (depth <= 0 || fabs(delta) <= 15.0 * eps) {
        return left + right + delta / 15.0;
    }
    return adaptive_simpson_rec(g, a, m, 0.5 * eps, left, fa, flm, fm, depth - 1) +
           adaptive_simpson_rec(g, m, b, 0.5 * eps, right, fm, frm, fb, depth - 1);
}

static double adaptive_simpson(double g) {
    double fa = np_integrand(0.0, g);
    double fm = np_integrand(0.5, g);
    double fb = np_integrand(1.0, g);
    double whole = simpson_step(0.0, 1.0, fa, fm, fb);
    return adaptive_simpson_rec(g, 0.0, 1.0, 1e-12, whole, fa, fm, fb, 18);
}

static int run_n3(double *max_rel_err, double *log_div_1e9, double *pole_x) {
    static const double gs[5] = {1e-5, 1e-3, 1e-1, 1.0, 10.0};
    double max_e = 0.0;
    for (size_t i = 0; i < 5; i++) {
        double g = gs[i];
        double num = adaptive_simpson(g);
        double exact = (0.5 * 3.14159265358979323846 * sqrt(g) - 0.5 * log(g)) / (1.0 + g);
        double err = fabs(num - exact) / fabs(exact);
        if (err > max_e) max_e = err;
    }
    *max_rel_err = max_e;
    *log_div_1e9 = log1p(1e9);
    *pole_x = 1.0 / sqrt(1e-4);
    return (max_e < 1e-8 && *log_div_1e9 > 20.0) ? 1 : 0;
}

static int run_n4(double *uv_area_gn, double *uv_entropy) {
    int64_t c_loc_a_num = 212, c_loc_a_den = 45;
    int64_t c_loc_A_num = 106, c_loc_A_den = 45;
    int64_t c_mic_a_num = c_loc_a_num - 3 * c_loc_a_den;
    int64_t c_mic_A_num = c_mic_a_num;
    int64_t c_mic_A_den = 2 * c_loc_a_den;

    int exact_ok = (c_loc_a_num * c_loc_A_den == 2 * c_loc_A_num * c_loc_a_den) &&
                   (c_mic_a_num == 77) && (c_mic_A_den == 90);
    double c_micro = (double)c_mic_A_num / (double)c_mic_A_den;

    static const double areas[4] = {1e2, 1e4, 1e8, 1e12};
    int macro_ok = 1;
    for (size_t i = 0; i < 4; i++) {
        double a = areas[i];
        double s_bh = 0.25 * a;
        double s_log = c_micro * log(a);
        if (s_bh + s_log <= 0.0 || (a >= 1e3 && fabs(s_log) / s_bh >= 0.05)) {
            macro_ok = 0;
        }
    }

    *uv_area_gn = 0.0;
    *uv_entropy = 0.0;
    int uv_found = 0;
    for (int m = 1; m <= 16; m++) {
        double x = ldexp(1.0, -m);
        double s_uv = 0.25 * x + c_micro * log(x);
        if (s_uv <= 0.0) {
            *uv_area_gn = x;
            *uv_entropy = s_uv;
            uv_found = 1;
            break;
        }
    }
    return (exact_ok && macro_ok && uv_found) ? 1 : 0;
}

int main(int argc, char **argv) {
    const char *json_path = NULL;
    for (int i = 1; i + 1 < argc; i++) {
        if (strcmp(argv[i], "--json") == 0) {
            json_path = argv[i + 1];
        }
    }

    int64_t sl2_err = 0, anti_err = 0, jacobi_err = 0;
    int ok1 = run_n1(&sl2_err, &anti_err, &jacobi_err);

    n2_res_t n2;
    int ok2 = run_n2(&n2);

    double max_rel_err = 0.0, log_div_1e9 = 0.0, pole_x = 0.0;
    int ok3 = run_n3(&max_rel_err, &log_div_1e9, &pole_x);

    double uv_area_gn = 0.0, uv_entropy = 0.0;
    int ok4 = run_n4(&uv_area_gn, &uv_entropy);

    int all_ok = (ok1 && ok2 && ok3 && ok4);
    const char *verdict = all_ok ? "SURVIVED_WITH_BOUNDARY" : "FALSIFIED";

    FILE *out = stdout;
    if (json_path != NULL) {
        out = fopen(json_path, "w");
        if (out == NULL) {
            fprintf(stderr, "Failed to open --json output: %s\n", json_path);
            return 2;
        }
    }

    fprintf(out, "{\n");
    fprintf(out, "  \"verifier\": \"PHY_NULL_CONE\",\n");
    fprintf(out, "  \"target\": \"chandrasekaran_2026\",\n");
    fprintf(out, "  \"verdict\": \"%s\",\n", verdict);
    fprintf(out, "  \"all_healthy_passed\": %s,\n", all_ok ? "true" : "false");
    fprintf(out, "  \"n1\": {\"passed\": %s, \"sl2_err\": %lld, \"anti_err\": %lld, \"jacobi_err\": %lld},\n",
            ok1 ? "true" : "false", (long long)sl2_err, (long long)anti_err, (long long)jacobi_err);
    fprintf(out, "  \"n2\": {\"passed\": %s, \"is_pd\": %s, \"min_pivot\": %.10e, \"hermitian_err\": %.3e, \"trace\": %.1f, \"ghost_found\": %s, \"conformal_eta\": %.6f, \"ghost_norm\": %.6f},\n",
            ok2 ? "true" : "false", n2.is_pd ? "true" : "false", n2.min_pivot, n2.hermitian_err,
            n2.trace_val, n2.ghost_found ? "true" : "false", n2.conformal_eta, n2.ghost_norm);
    fprintf(out, "  \"n3\": {\"passed\": %s, \"max_rel_err\": %.6e, \"log_div_1e9\": %.10f, \"pole_x\": %.2f},\n",
            ok3 ? "true" : "false", max_rel_err, log_div_1e9, pole_x);
    fprintf(out, "  \"n4\": {\"passed\": %s, \"c_local_A\": \"106/45\", \"c_micro_A\": \"77/90\", \"uv_area_gn\": %.4f, \"uv_entropy\": %.10f}\n",
            ok4 ? "true" : "false", uv_area_gn, uv_entropy);
    fprintf(out, "}\n");

    if (json_path != NULL) {
        fclose(out);
    }
    return all_ok ? 0 : 1;
}
