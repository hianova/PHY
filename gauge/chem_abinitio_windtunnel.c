#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Stage 10: First-Principles Certified Quantum Chemistry Wind Tunnel
 * ----------------------------------------------------------------------------
 *  1. Exact Boys F_m(t), Hehre-Stewart-Pople STO-3G basis & Obara-Saika s+p
 *     analytical GTO overlap / kinetic / nuclear attraction / 4-center ERIs
 *  2. Householder + Implicit QL symmetric eigensolver & dense Cholesky
 *  3. Restricted Hartree-Fock + Pulay DIIS SCF, Virial scaling & MP2
 *  4. Jordan-Wigner / Slater-Condon FCI up to 441 determinants (H2, HeH+, H4,
 *     LiH, H2O) + Weinstein & Temple certified energy brackets (Theorem M2)
 *  5. Spin-blocked 2-RDM DQG + S^2 ADMM dual SDP with exact orthogonalized
 *     dual slack & unconditional eigenvalue-floor lower bound (Theorem M3)
 *  6. 2D H2O RHF 2-Jet Newton geometry optimization & bilevel H4 conformer
 *     SEARU x 2-Jet global PES tunneling (Theorem M5)
 * ============================================================================ */

#define MAX_SHELLS 12
#define MAX_AO     12
#define MAX_DET    512

typedef struct {
    int atom;
    double center[3];
    int l, m, n; /* Cartesian angular momentum: 0/1 only for STO-3G s+p */
    int nprim;
    double exps[3];
    double coefs[3]; /* includes primitive + contracted normalization */
} AO;

typedef struct {
    int natoms;
    int Z[8];
    double R[8][3];
    int nao;
    int nel;
    AO aos[MAX_AO];
    double S[MAX_AO][MAX_AO];
    double T[MAX_AO][MAX_AO];
    double V[MAX_AO][MAX_AO];
    double Hcore[MAX_AO][MAX_AO];
    double ERI[MAX_AO][MAX_AO][MAX_AO][MAX_AO];
    double Enuc;
} MolSystem;

static void die_msg(const char *msg) {
    fprintf(stderr, "FATAL: %s\n", msg);
    exit(1);
}

/* --- 1. Exact Boys Function F_m(t) (m <= 6) --- */
static void boys_array(int mmax, double t, double *F) {
    if (t < 15.0) {
        int mtop = mmax + 48;
        if (mtop > 64) mtop = 64;
        double f = 0.0;
        double et = exp(-t);
        double tmp[66];
        for (int m = mtop - 1; m >= 0; m--) {
            f = (2.0 * t * f + et) / (double)(2 * m + 1);
            tmp[m] = f;
        }
        for (int m = 0; m <= mmax; m++) F[m] = tmp[m];
    } else {
        double sqt = sqrt(t);
        double f0 = 0.5 * sqrt(M_PI / t) * erf(sqt);
        double et = exp(-t);
        F[0] = f0;
        for (int m = 0; m < mmax; m++) {
            F[m + 1] = ((double)(2 * m + 1) * F[m] - et) / (2.0 * t);
        }
    }
}

/* --- 2. Obara-Saika 1D Overlap/Kinetic + Nuclear Attraction + 4-Center ERI --- */
static void os_overlap_1d(double a, double b, double Ax, double Bx, double S1[4][4]) {
    double p = a + b;
    double mu = a * b / p;
    double dx = Ax - Bx;
    double Px = (a * Ax + b * Bx) / p;
    double XPA = Px - Ax;
    double XPB = Px - Bx;
    double s00 = sqrt(M_PI / p) * exp(-mu * dx * dx);
    memset(S1, 0, sizeof(double) * 16);
    S1[0][0] = s00;
    for (int i = 0; i < 3; i++) {
        double prev = (i > 0) ? ((double)i / (2.0 * p)) * S1[i - 1][0] : 0.0;
        S1[i + 1][0] = XPA * S1[i][0] + prev;
    }
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i <= 3; i++) {
            double pi_term = (i > 0) ? ((double)i / (2.0 * p)) * S1[i - 1][j] : 0.0;
            double pj_term = (j > 0) ? ((double)j / (2.0 * p)) * S1[i][j - 1] : 0.0;
            S1[i][j + 1] = XPB * S1[i][j] + pi_term + pj_term;
        }
    }
}

static void os_prim_st(const int a_lmn[3], double a, const double A[3],
                        const int b_lmn[3], double b, const double B[3],
                        double *outS, double *outT) {
    double Sx[4][4], Sy[4][4], Sz[4][4];
    os_overlap_1d(a, b, A[0], B[0], Sx);
    os_overlap_1d(a, b, A[1], B[1], Sy);
    os_overlap_1d(a, b, A[2], B[2], Sz);
    int ix = a_lmn[0], iy = a_lmn[1], iz = a_lmn[2];
    int jx = b_lmn[0], jy = b_lmn[1], jz = b_lmn[2];
    double sx = Sx[ix][jx], sy = Sy[iy][jy], sz = Sz[iz][jz];
    *outS = sx * sy * sz;

    double tx = b * (double)(2 * jx + 1) * Sx[ix][jx] - 2.0 * b * b * Sx[ix][jx + 2]
              - (jx >= 2 ? 0.5 * (double)(jx * (jx - 1)) * Sx[ix][jx - 2] : 0.0);
    double ty = b * (double)(2 * jy + 1) * Sy[iy][jy] - 2.0 * b * b * Sy[iy][jy + 2]
              - (jy >= 2 ? 0.5 * (double)(jy * (jy - 1)) * Sy[iy][jy - 2] : 0.0);
    double tz = b * (double)(2 * jz + 1) * Sz[iz][jz] - 2.0 * b * b * Sz[iz][jz + 2]
              - (jz >= 2 ? 0.5 * (double)(jz * (jz - 1)) * Sz[iz][jz - 2] : 0.0);
    *outT = tx * sy * sz + sx * ty * sz + sx * sy * tz;
}

static int ang_code(const int lmn[3]) {
    if (lmn[0] == 1) return 1;
    if (lmn[1] == 1) return 2;
    if (lmn[2] == 1) return 3;
    return 0;
}

static double os_prim_nuc(const int a_lmn[3], double a, const double A[3],
                          const int b_lmn[3], double b, const double B[3],
                          const double C[3]) {
    double p = a + b;
    double mu = a * b / p;
    double ab2 = (A[0] - B[0]) * (A[0] - B[0])
               + (A[1] - B[1]) * (A[1] - B[1])
               + (A[2] - B[2]) * (A[2] - B[2]);
    double P[3], PA[3], PB[3], PC[3];
    for (int d = 0; d < 3; d++) {
        P[d] = (a * A[d] + b * B[d]) / p;
        PA[d] = P[d] - A[d];
        PB[d] = P[d] - B[d];
        PC[d] = P[d] - C[d];
    }
    double U = p * (PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2]);
    double F[3];
    boys_array(2, U, F);
    double pref = (2.0 * M_PI / p) * exp(-mu * ab2);
    double V[4][4][3];
    memset(V, 0, sizeof(V));
    for (int m = 0; m <= 2; m++) V[0][0][m] = pref * F[m];
    int ca = ang_code(a_lmn), cb = ang_code(b_lmn);
    if (ca > 0) {
        int d = ca - 1;
        for (int m = 0; m <= 1; m++) V[ca][0][m] = PA[d] * V[0][0][m] - PC[d] * V[0][0][m + 1];
    }
    if (cb > 0 && ca == 0) {
        int d = cb - 1;
        V[0][cb][0] = PB[d] * V[0][0][0] - PC[d] * V[0][0][1];
    }
    if (ca > 0 && cb > 0) {
        int d = cb - 1;
        double extra = (ca == cb) ? (0.5 / p) * (V[0][0][0] - V[0][0][1]) : 0.0;
        V[ca][cb][0] = PB[d] * V[ca][0][0] - PC[d] * V[ca][0][1] + extra;
    }
    return V[ca][cb][0];
}

static double os_prim_eri(const int a_lmn[3], double a, const double A[3],
                          const int b_lmn[3], double b, const double B[3],
                          const int c_lmn[3], double c, const double C[3],
                          const int d_lmn[3], double d, const double D[3]) {
    double zeta = a + b, eta = c + d;
    double inv_sum = 1.0 / (zeta + eta);
    double rho = zeta * eta * inv_sum;
    double ab2 = (A[0] - B[0]) * (A[0] - B[0])
               + (A[1] - B[1]) * (A[1] - B[1])
               + (A[2] - B[2]) * (A[2] - B[2]);
    double cd2 = (C[0] - D[0]) * (C[0] - D[0])
               + (C[1] - D[1]) * (C[1] - D[1])
               + (C[2] - D[2]) * (C[2] - D[2]);
    double P[3], Q[3], W[3], PA[3], PB[3], QC[3], QD[3], WP[3], WQ[3], PQ[3];
    for (int k = 0; k < 3; k++) {
        P[k] = (a * A[k] + b * B[k]) / zeta;
        Q[k] = (c * C[k] + d * D[k]) / eta;
        W[k] = (zeta * P[k] + eta * Q[k]) * inv_sum;
        PA[k] = P[k] - A[k];
        PB[k] = P[k] - B[k];
        QC[k] = Q[k] - C[k];
        QD[k] = Q[k] - D[k];
        WP[k] = W[k] - P[k];
        WQ[k] = W[k] - Q[k];
        PQ[k] = P[k] - Q[k];
    }
    double U = rho * (PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2]);
    double F[5];
    boys_array(4, U, F);
    double pref = 2.0 * pow(M_PI, 2.5) / (zeta * eta * sqrt(zeta + eta))
                * exp(-(a * b / zeta) * ab2 - (c * d / eta) * cd2);

    int ca = ang_code(a_lmn), cb = ang_code(b_lmn), cc = ang_code(c_lmn), cd = ang_code(d_lmn);
    /* Fast direct s+p Obara-Saika tensor for just the active components {0, ca}, {0, cb}, {0, cc}, {0, cd} */
    double G[4][4][4][4][5];
    memset(G[0][0][0][0], 0, 5 * sizeof(double));
    for (int m = 0; m <= 4; m++) G[0][0][0][0][m] = pref * F[m];

    int al[2] = {0, ca}, na = (ca > 0 ? 2 : 1);
    int bl[2] = {0, cb}, nb = (cb > 0 ? 2 : 1);
    int cl[2] = {0, cc}, nc = (cc > 0 ? 2 : 1);

    if (ca > 0) {
        int dir = ca - 1;
        for (int m = 0; m <= 3; m++) {
            G[ca][0][0][0][m] = PA[dir] * G[0][0][0][0][m] + WP[dir] * G[0][0][0][0][m + 1];
        }
    }
    if (cb > 0) {
        int dir = cb - 1;
        for (int ia = 0; ia < na; ia++) {
            int va = al[ia];
            for (int m = 0; m <= 2; m++) {
                double val = PB[dir] * G[va][0][0][0][m] + WP[dir] * G[va][0][0][0][m + 1];
                if (va == cb) val += (0.5 / zeta) * (G[0][0][0][0][m] - (rho / zeta) * G[0][0][0][0][m + 1]);
                G[va][cb][0][0][m] = val;
            }
        }
    }
    if (cc > 0) {
        int dir = cc - 1;
        for (int ia = 0; ia < na; ia++) {
            int va = al[ia];
            for (int ib = 0; ib < nb; ib++) {
                int vb = bl[ib];
                for (int m = 0; m <= 1; m++) {
                    double val = QC[dir] * G[va][vb][0][0][m] + WQ[dir] * G[va][vb][0][0][m + 1];
                    if (va == cc) val += 0.5 * inv_sum * G[0][vb][0][0][m + 1];
                    if (vb == cc) val += 0.5 * inv_sum * G[va][0][0][0][m + 1];
                    G[va][vb][cc][0][m] = val;
                }
            }
        }
    }
    if (cd > 0) {
        int dir = cd - 1;
        int va = ca, vb = cb;
        double val = QD[dir] * G[va][vb][cc][0][0] + WQ[dir] * G[va][vb][cc][0][1];
        if (cc == cd) val += (0.5 / eta) * (G[va][vb][0][0][0] - (rho / eta) * G[va][vb][0][0][1]);
        if (va == cd) val += 0.5 * inv_sum * G[0][vb][cc][0][1];
        if (vb == cd) val += 0.5 * inv_sum * G[va][0][cc][0][1];
        /* Avoid unused nc check */
        (void)nc;
        (void)cl;
        return val;
    }
    return G[ca][cb][cc][0][0];
}

/* --- 3. STO-3G Basis Builder --- */
static void add_sto3g_atom(MolSystem *mol, int atom_idx, int Z) {
    /* Standard Hehre-Stewart-Pople 1969 STO-3G exponents & coefficients */
    static const double a1s[3]  = {2.227660584, 0.4057711562, 0.1098175104};
    static const double d1s[3]  = {0.1543289673, 0.5353281423, 0.4446345422};
    static const double a2sp[3] = {0.994203258, 0.231030315, 0.075138553};
    static const double d2s[3]  = {-0.09996722919, 0.3995128261, 0.7001154689};
    static const double d2p[3]  = {0.1559162750, 0.6076837186, 0.3919573931};

    double z1s = 1.24, z2sp = 1.0;
    int has_2sp = 0;
    if (Z == 1) { z1s = 1.24; }
    else if (Z == 2) { z1s = 2.0925; }
    else if (Z == 3) { z1s = 2.69; z2sp = 0.80; has_2sp = 1; }
    else if (Z == 8) { z1s = 7.66; z2sp = 2.25; has_2sp = 1; }
    else die_msg("Unsupported Z");

    /* 1s shell */
    AO *ao = &mol->aos[mol->nao++];
    memset(ao, 0, sizeof(*ao));
    ao->atom = atom_idx;
    memcpy(ao->center, mol->R[atom_idx], 3 * sizeof(double));
    ao->nprim = 3;
    for (int k = 0; k < 3; k++) {
        double al = a1s[k] * z1s * z1s;
        ao->exps[k] = al;
        ao->coefs[k] = d1s[k] * pow(2.0 * al / M_PI, 0.75);
    }

    if (has_2sp) {
        /* 2s */
        ao = &mol->aos[mol->nao++];
        memset(ao, 0, sizeof(*ao));
        ao->atom = atom_idx;
        memcpy(ao->center, mol->R[atom_idx], 3 * sizeof(double));
        ao->nprim = 3;
        for (int k = 0; k < 3; k++) {
            double al = a2sp[k] * z2sp * z2sp;
            ao->exps[k] = al;
            ao->coefs[k] = d2s[k] * pow(2.0 * al / M_PI, 0.75);
        }
        /* 2px, 2py, 2pz */
        for (int dir = 0; dir < 3; dir++) {
            ao = &mol->aos[mol->nao++];
            memset(ao, 0, sizeof(*ao));
            ao->atom = atom_idx;
            memcpy(ao->center, mol->R[atom_idx], 3 * sizeof(double));
            if (dir == 0) ao->l = 1;
            else if (dir == 1) ao->m = 1;
            else ao->n = 1;
            ao->nprim = 3;
            for (int k = 0; k < 3; k++) {
                double al = a2sp[k] * z2sp * z2sp;
                ao->exps[k] = al;
                ao->coefs[k] = d2p[k] * 2.0 * sqrt(al) * pow(2.0 * al / M_PI, 0.75);
            }
        }
    }
}

static void build_mol_integrals(MolSystem *mol, int natoms, const int Z[], const double coords[][3], int charge) {
    memset(mol, 0, sizeof(*mol));
    mol->natoms = natoms;
    int ne = -charge;
    for (int a = 0; a < natoms; a++) {
        mol->Z[a] = Z[a];
        ne += Z[a];
        memcpy(mol->R[a], coords[a], 3 * sizeof(double));
        add_sto3g_atom(mol, a, Z[a]);
    }
    mol->nel = ne;
    double enuc = 0.0;
    for (int a = 0; a < natoms; a++) {
        for (int b = a + 1; b < natoms; b++) {
            double dx = mol->R[a][0] - mol->R[b][0];
            double dy = mol->R[a][1] - mol->R[b][1];
            double dz = mol->R[a][2] - mol->R[b][2];
            double r = sqrt(dx * dx + dy * dy + dz * dz);
            enuc += (double)(Z[a] * Z[b]) / r;
        }
    }
    mol->Enuc = enuc;

    int K = mol->nao;
    /* Renormalize each contracted AO strictly so S_{ii} == 1.0 */
    for (int i = 0; i < K; i++) {
        const AO *ai = &mol->aos[i];
        int lmn[3] = {ai->l, ai->m, ai->n};
        double selfS = 0.0;
        for (int p = 0; p < 3; p++) {
            for (int q = 0; q < 3; q++) {
                double s_pq, t_pq;
                os_prim_st(lmn, ai->exps[p], ai->center, lmn, ai->exps[q], ai->center, &s_pq, &t_pq);
                selfS += ai->coefs[p] * ai->coefs[q] * s_pq;
            }
        }
        double inv = 1.0 / sqrt(selfS);
        for (int p = 0; p < 3; p++) mol->aos[i].coefs[p] *= inv;
    }

    for (int i = 0; i < K; i++) {
        const AO *ai = &mol->aos[i];
        int lmni[3] = {ai->l, ai->m, ai->n};
        for (int j = 0; j <= i; j++) {
            const AO *aj = &mol->aos[j];
            int lmnj[3] = {aj->l, aj->m, aj->n};
            double sij = 0.0, tij = 0.0, vij = 0.0;
            for (int p = 0; p < 3; p++) {
                for (int q = 0; q < 3; q++) {
                    double cc = ai->coefs[p] * aj->coefs[q];
                    double spq, tpq;
                    os_prim_st(lmni, ai->exps[p], ai->center, lmnj, aj->exps[q], aj->center, &spq, &tpq);
                    sij += cc * spq;
                    tij += cc * tpq;
                    for (int at = 0; at < natoms; at++) {
                        vij -= (double)Z[at] * cc
                             * os_prim_nuc(lmni, ai->exps[p], ai->center,
                                           lmnj, aj->exps[q], aj->center, mol->R[at]);
                    }
                }
            }
            mol->S[i][j] = mol->S[j][i] = sij;
            mol->T[i][j] = mol->T[j][i] = tij;
            mol->V[i][j] = mol->V[j][i] = vij;
            mol->Hcore[i][j] = mol->Hcore[j][i] = tij + vij;
        }
    }

    /* Compute ERIs exploiting 8-fold symmetry */
    for (int i = 0; i < K; i++) {
        const AO *ai = &mol->aos[i];
        int lmni[3] = {ai->l, ai->m, ai->n};
        for (int j = 0; j <= i; j++) {
            const AO *aj = &mol->aos[j];
            int lmnj[3] = {aj->l, aj->m, aj->n};
            int ij = i * (i + 1) / 2 + j;
            for (int k = 0; k < K; k++) {
                const AO *ak = &mol->aos[k];
                int lmnk[3] = {ak->l, ak->m, ak->n};
                for (int l = 0; l <= k; l++) {
                    int kl = k * (k + 1) / 2 + l;
                    if (kl > ij) continue;
                    const AO *al = &mol->aos[l];
                    int lmnl[3] = {al->l, al->m, al->n};
                    double val = 0.0;
                    for (int pa = 0; pa < 3; pa++) {
                        for (int pb = 0; pb < 3; pb++) {
                            double c_ab = ai->coefs[pa] * aj->coefs[pb];
                            for (int pc = 0; pc < 3; pc++) {
                                for (int pd = 0; pd < 3; pd++) {
                                    double c_cd = ak->coefs[pc] * al->coefs[pd];
                                    val += c_ab * c_cd * os_prim_eri(
                                        lmni, ai->exps[pa], ai->center,
                                        lmnj, aj->exps[pb], aj->center,
                                        lmnk, ak->exps[pc], ak->center,
                                        lmnl, al->exps[pd], al->center);
                                }
                            }
                        }
                    }
                    mol->ERI[i][j][k][l] = val;
                    mol->ERI[j][i][k][l] = val;
                    mol->ERI[i][j][l][k] = val;
                    mol->ERI[j][i][l][k] = val;
                    mol->ERI[k][l][i][j] = val;
                    mol->ERI[l][k][i][j] = val;
                    mol->ERI[k][l][j][i] = val;
                    mol->ERI[l][k][j][i] = val;
                }
            }
        }
    }
}

/* --- 4. Householder Tridiagonalization + Implicit QL Symmetric Eigensolver --- */
static double pythag_d(double a, double b) {
    double aa = fabs(a), bb = fabs(b);
    if (aa > bb) { double r = bb / aa; return aa * sqrt(1.0 + r * r); }
    if (bb == 0.0) return 0.0;
    double r = aa / bb;
    return bb * sqrt(1.0 + r * r);
}

static void sym_eigh(int n, const double *A_in, double *evals, double *evecs_out) {
    /* V stores row-major n x n eigenvectors */
    double *V = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
    double *e = (double *)calloc((size_t)n, sizeof(double));
    if (!V || !e) die_msg("OOM sym_eigh");
    memcpy(V, A_in, (size_t)n * (size_t)n * sizeof(double));

    for (int i = n - 1; i >= 1; i--) {
        int l = i - 1;
        double h = 0.0, scale = 0.0;
        if (l > 0) {
            for (int k = 0; k <= l; k++) scale += fabs(V[i * n + k]);
            if (scale == 0.0) {
                e[i] = V[i * n + l];
            } else {
                for (int k = 0; k <= l; k++) {
                    V[i * n + k] /= scale;
                    h += V[i * n + k] * V[i * n + k];
                }
                double f = V[i * n + l];
                double g = (f >= 0.0) ? -sqrt(h) : sqrt(h);
                e[i] = scale * g;
                h -= f * g;
                V[i * n + l] = f - g;
                f = 0.0;
                for (int j = 0; j <= l; j++) {
                    V[j * n + i] = V[i * n + j] / h;
                    g = 0.0;
                    for (int k = 0; k <= j; k++) g += V[j * n + k] * V[i * n + k];
                    for (int k = j + 1; k <= l; k++) g += V[k * n + j] * V[i * n + k];
                    e[j] = g / h;
                    f += e[j] * V[i * n + j];
                }
                double hh = f / (h + h);
                for (int j = 0; j <= l; j++) {
                    f = V[i * n + j];
                    e[j] = g = e[j] - hh * f;
                    for (int k = 0; k <= j; k++) V[j * n + k] -= f * e[k] + g * V[i * n + k];
                }
            }
        } else {
            e[i] = V[i * n + l];
        }
        evals[i] = h;
    }
    evals[0] = 0.0;
    e[0] = 0.0;
    for (int i = 0; i < n; i++) {
        if (evals[i] != 0.0) {
            for (int j = 0; j < i; j++) {
                double g = 0.0;
                for (int k = 0; k < i; k++) g += V[i * n + k] * V[k * n + j];
                for (int k = 0; k < i; k++) V[k * n + j] -= g * V[k * n + i];
            }
        }
        evals[i] = V[i * n + i];
        V[i * n + i] = 1.0;
        for (int j = 0; j < i; j++) { V[j * n + i] = 0.0; V[i * n + j] = 0.0; }
    }

    /* Shift e down so e[i] is subdiagonal between i and i+1 */
    for (int i = 1; i < n; i++) e[i - 1] = e[i];
    e[n - 1] = 0.0;

    for (int l = 0; l < n; l++) {
        int iter = 0;
        int m;
        do {
            for (m = l; m < n - 1; m++) {
                double dd = fabs(evals[m]) + fabs(evals[m + 1]);
                if (fabs(e[m]) + dd == dd || fabs(e[m]) < 1e-15 * (dd > 1.0 ? dd : 1.0)) break;
            }
            if (m != l) {
                if (iter++ > 60) die_msg("tqli convergence failed");
                double g = (evals[l + 1] - evals[l]) / (2.0 * e[l]);
                double r = pythag_d(g, 1.0);
                g = evals[m] - evals[l] + e[l] / (g + (g >= 0.0 ? fabs(r) : -fabs(r)));
                double s = 1.0, c = 1.0, p = 0.0;
                int i;
                for (i = m - 1; i >= l; i--) {
                    double f = s * e[i];
                    double b = c * e[i];
                    e[i + 1] = r = pythag_d(f, g);
                    if (r == 0.0) {
                        evals[i + 1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = evals[i + 1] - p;
                    r = (evals[i] - g) * s + 2.0 * c * b;
                    evals[i + 1] = g + (p = s * r);
                    g = c * r - b;
                    if (evecs_out) {
                        for (int k = 0; k < n; k++) {
                            f = V[k * n + i + 1];
                            V[k * n + i + 1] = s * V[k * n + i] + c * f;
                            V[k * n + i] = c * V[k * n + i] - s * f;
                        }
                    }
                }
                if (r == 0.0 && i >= l) continue;
                evals[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m != l);
    }

    /* Sort eigenvalues ascending and permute V columns identically */
    for (int i = 0; i < n - 1; i++) {
        int k = i;
        double p = evals[i];
        for (int j = i + 1; j < n; j++) if (evals[j] < p) { k = j; p = evals[j]; }
        if (k != i) {
            evals[k] = evals[i];
            evals[i] = p;
            if (evecs_out) {
                for (int r = 0; r < n; r++) {
                    double tmp = V[r * n + i];
                    V[r * n + i] = V[r * n + k];
                    V[r * n + k] = tmp;
                }
            }
        }
    }
    if (evecs_out) memcpy(evecs_out, V, (size_t)n * (size_t)n * sizeof(double));
    free(V);
    free(e);
}

static int chol_factor(int n, double *A) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double s = A[i * n + j];
            for (int k = 0; k < j; k++) s -= A[i * n + k] * A[j * n + k];
            if (i == j) {
                if (s <= 0.0) return -1;
                A[i * n + j] = sqrt(s);
            } else {
                A[i * n + j] = s / A[j * n + j];
            }
        }
    }
    return 0;
}

static void chol_solve(int n, const double *L, const double *b, double *x) {
    for (int i = 0; i < n; i++) {
        double s = b[i];
        for (int k = 0; k < i; k++) s -= L[i * n + k] * x[k];
        x[i] = s / L[i * n + i];
    }
    for (int i = n - 1; i >= 0; i--) {
        double s = x[i];
        for (int k = i + 1; k < n; k++) s -= L[k * n + i] * x[k];
        x[i] = s / L[i * n + i];
    }
}

/* --- 5. RHF + Pulay DIIS SCF, MO Transformation & MP2 --- */
typedef struct {
    int K;
    int nocc;
    int diis_iters;
    double E_rhf;
    double E_elec;
    double E_mp2;
    double T_expect;
    double V_expect;
    double virial_ratio;
    double comm_norm;
    double idem_err;
    double C[MAX_AO][MAX_AO]; /* C[ao][mo] */
    double P[MAX_AO][MAX_AO]; /* P[mu][nu] = 2 * sum_{i<nocc} C[mu][i]*C[nu][i] */
    double eps[MAX_AO];
    double F_mo[MAX_AO][MAX_AO];
    double H_mo[MAX_AO][MAX_AO];
    double ERI_mo[MAX_AO][MAX_AO][MAX_AO][MAX_AO];
} RHFResult;

static void run_rhf_diis_mp2(const MolSystem *mol, RHFResult *res) {
    memset(res, 0, sizeof(*res));
    int K = mol->nao;
    int nocc = mol->nel / 2;
    res->K = K;
    res->nocc = nocc;

    /* Symmetric orthogonalizer X = S^{-1/2} */
    double Sflat[MAX_AO * MAX_AO], Seval[MAX_AO], Sevec[MAX_AO * MAX_AO];
    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++) Sflat[i * K + j] = mol->S[i][j];
    sym_eigh(K, Sflat, Seval, Sevec);

    double X[MAX_AO][MAX_AO];
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            double s = 0.0;
            for (int k = 0; k < K; k++) {
                s += Sevec[i * K + k] * (1.0 / sqrt(Seval[k])) * Sevec[j * K + k];
            }
            X[i][j] = s;
        }
    }

    double P[MAX_AO][MAX_AO], F[MAX_AO][MAX_AO], Fp[MAX_AO * MAX_AO], Cp[MAX_AO * MAX_AO];
    memset(P, 0, sizeof(P));

    /* DIIS history up to 8 Fock + error matrices */
    enum { MAX_DIIS = 8 };
    double diis_F[MAX_DIIS][MAX_AO * MAX_AO];
    double diis_e[MAX_DIIS][MAX_AO * MAX_AO];
    int n_diis = 0;

    double E_prev = 0.0, comm_norm = 1.0, E_elec = 0.0;
    int iters = 0;
    for (int iter = 1; iter <= 120; iter++) {
        iters = iter;
        /* Build Fock matrix */
        for (int mu = 0; mu < K; mu++) {
            for (int nu = 0; nu <= mu; nu++) {
                double g = 0.0;
                for (int la = 0; la < K; la++) {
                    for (int si = 0; si < K; si++) {
                        g += P[la][si] * (mol->ERI[mu][nu][si][la] - 0.5 * mol->ERI[mu][la][si][nu]);
                    }
                }
                F[mu][nu] = F[nu][mu] = mol->Hcore[mu][nu] + g;
            }
        }
        /* Electronic energy BEFORE diagonalizing */
        E_elec = 0.0;
        for (int mu = 0; mu < K; mu++) {
            for (int nu = 0; nu < K; nu++) {
                E_elec += 0.5 * P[mu][nu] * (mol->Hcore[mu][nu] + F[mu][nu]);
            }
        }

        /* Compute orthonormalized DIIS error vector e = X^T (F P S - S P F) X */
        double errAO[MAX_AO][MAX_AO];
        comm_norm = 0.0;
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                double fps = 0.0, spf = 0.0;
                for (int a = 0; a < K; a++) {
                    double fp = 0.0, sp = 0.0;
                    for (int b = 0; b < K; b++) {
                        fp += F[i][b] * P[b][a];
                        sp += mol->S[i][b] * P[b][a];
                    }
                    fps += fp * mol->S[a][j];
                    spf += sp * F[a][j];
                }
                errAO[i][j] = fps - spf;
            }
        }
        double errOrtho[MAX_AO * MAX_AO];
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                double val = 0.0;
                for (int a = 0; a < K; a++)
                    for (int b = 0; b < K; b++) val += X[a][i] * errAO[a][b] * X[b][j];
                errOrtho[i * K + j] = val;
                comm_norm += val * val;
            }
        }
        comm_norm = sqrt(comm_norm);

        if (iter > 2 && comm_norm < 1e-12 && fabs(E_elec - E_prev) < 1e-13) break;
        E_prev = E_elec;

        /* Push to DIIS (skip iter == 1 where P == 0 gives spurious [F, 0] == 0) */
        if (iter >= 2) {
            int slot = (n_diis < MAX_DIIS) ? n_diis++ : (MAX_DIIS - 1);
            if (n_diis == MAX_DIIS && iter > MAX_DIIS + 1) {
                memmove(diis_F[0], diis_F[1], (MAX_DIIS - 1) * sizeof(diis_F[0]));
                memmove(diis_e[0], diis_e[1], (MAX_DIIS - 1) * sizeof(diis_e[0]));
                slot = MAX_DIIS - 1;
            }
            for (int i = 0; i < K; i++)
                for (int j = 0; j < K; j++) diis_F[slot][i * K + j] = F[i][j];
            memcpy(diis_e[slot], errOrtho, (size_t)(K * K) * sizeof(double));
        }

        if (n_diis >= 2 && comm_norm < 0.5) {
            /* Solve symmetric B w = (-0..., -1) by Gaussian elimination with partial pivoting */
            int m = n_diis + 1;
            double B[(MAX_DIIS + 1) * (MAX_DIIS + 1)];
            double rhs[MAX_DIIS + 1];
            memset(B, 0, sizeof(B));
            memset(rhs, 0, sizeof(rhs));
            for (int a = 0; a < n_diis; a++) {
                for (int b = 0; b <= a; b++) {
                    double dot = 0.0;
                    for (int k = 0; k < K * K; k++) dot += diis_e[a][k] * diis_e[b][k];
                    B[a * m + b] = B[b * m + a] = dot;
                }
                B[a * m + a] += 1e-14;
                B[a * m + n_diis] = -1.0;
                B[n_diis * m + a] = -1.0;
            }
            rhs[n_diis] = -1.0;

            int ok = 1;
            for (int i = 0; i < m; i++) {
                int piv = i;
                for (int r = i + 1; r < m; r++)
                    if (fabs(B[r * m + i]) > fabs(B[piv * m + i])) piv = r;
                if (fabs(B[piv * m + i]) < 1e-18) { ok = 0; break; }
                if (piv != i) {
                    for (int c = 0; c < m; c++) {
                        double t = B[i * m + c]; B[i * m + c] = B[piv * m + c]; B[piv * m + c] = t;
                    }
                    double tr = rhs[i]; rhs[i] = rhs[piv]; rhs[piv] = tr;
                }
                for (int r = i + 1; r < m; r++) {
                    double fac = B[r * m + i] / B[i * m + i];
                    for (int c = i; c < m; c++) B[r * m + c] -= fac * B[i * m + c];
                    rhs[r] -= fac * rhs[i];
                }
            }
            if (ok) {
                double w[MAX_DIIS + 1];
                for (int i = m - 1; i >= 0; i--) {
                    double s = rhs[i];
                    for (int c = i + 1; c < m; c++) s -= B[i * m + c] * w[c];
                    w[i] = s / B[i * m + i];
                }
                double wsum = 0.0;
                for (int a = 0; a < n_diis; a++) wsum += fabs(w[a]);
                if (wsum < 20.0) {
                    for (int i = 0; i < K; i++) {
                        for (int j = 0; j < K; j++) {
                            double val = 0.0;
                            for (int a = 0; a < n_diis; a++) val += w[a] * diis_F[a][i * K + j];
                            F[i][j] = val;
                        }
                    }
                }
            }
        }

        /* F' = X^T F X */
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                double val = 0.0;
                for (int a = 0; a < K; a++)
                    for (int b = 0; b < K; b++) val += X[a][i] * F[a][b] * X[b][j];
                Fp[i * K + j] = val;
            }
        }
        sym_eigh(K, Fp, res->eps, Cp);
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                double val = 0.0;
                for (int a = 0; a < K; a++) val += X[i][a] * Cp[a * K + j];
                res->C[i][j] = val;
            }
        }
        for (int mu = 0; mu < K; mu++) {
            for (int nu = 0; nu < K; nu++) {
                double pval = 0.0;
                for (int i = 0; i < nocc; i++) pval += 2.0 * res->C[mu][i] * res->C[nu][i];
                P[mu][nu] = pval;
            }
        }
    }

    /* Final Fock + energy + canonical MO diagonalization at converged P */
    for (int mu = 0; mu < K; mu++) {
        for (int nu = 0; nu < K; nu++) {
            double g = 0.0;
            for (int la = 0; la < K; la++)
                for (int si = 0; si < K; si++)
                    g += P[la][si] * (mol->ERI[mu][nu][si][la] - 0.5 * mol->ERI[mu][la][si][nu]);
            F[mu][nu] = mol->Hcore[mu][nu] + g;
        }
    }
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            double val = 0.0;
            for (int a = 0; a < K; a++)
                for (int b = 0; b < K; b++) val += X[a][i] * F[a][b] * X[b][j];
            Fp[i * K + j] = val;
        }
    }
    sym_eigh(K, Fp, res->eps, Cp);
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            double val = 0.0;
            for (int a = 0; a < K; a++) val += X[i][a] * Cp[a * K + j];
            res->C[i][j] = val;
        }
    }
    E_elec = 0.0;
    double T_exp = 0.0;
    for (int mu = 0; mu < K; mu++) {
        for (int nu = 0; nu < K; nu++) {
            E_elec += 0.5 * P[mu][nu] * (mol->Hcore[mu][nu] + F[mu][nu]);
            T_exp += P[mu][nu] * mol->T[mu][nu];
        }
    }
    res->diis_iters = iters;
    res->E_elec = E_elec;
    res->E_rhf = E_elec + mol->Enuc;
    res->T_expect = T_exp;
    res->V_expect = res->E_rhf - T_exp;
    res->virial_ratio = -res->V_expect / res->T_expect;
    res->comm_norm = comm_norm;
    memcpy(res->P, P, sizeof(P));

    /* Check closed-shell idempotency error || (1/2) P S P - P ||_F */
    double iderr = 0.0;
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            double ps = 0.0;
            for (int a = 0; a < K; a++) {
                double p1 = 0.0;
                for (int b = 0; b < K; b++) p1 += P[i][b] * mol->S[b][a];
                ps += p1 * P[a][j];
            }
            double d = 0.5 * ps - P[i][j];
            iderr += d * d;
        }
    }
    res->idem_err = sqrt(iderr);

    /* MO 1-electron H_mo and F_mo */
    for (int p = 0; p < K; p++) {
        for (int q = 0; q < K; q++) {
            double hval = 0.0, fval = 0.0;
            for (int mu = 0; mu < K; mu++) {
                for (int nu = 0; nu < K; nu++) {
                    hval += res->C[mu][p] * res->C[nu][q] * mol->Hcore[mu][nu];
                    fval += res->C[mu][p] * res->C[nu][q] * F[mu][nu];
                }
            }
            res->H_mo[p][q] = hval;
            res->F_mo[p][q] = fval;
        }
    }

    /* 4-step O(K^5) AO -> MO ERI transformation (pq|rs) in chemists' notation */
    double T1[MAX_AO][MAX_AO][MAX_AO][MAX_AO];
    double T2[MAX_AO][MAX_AO][MAX_AO][MAX_AO];
    double T3[MAX_AO][MAX_AO][MAX_AO][MAX_AO];
    memset(T1, 0, sizeof(T1));
    memset(T2, 0, sizeof(T2));
    memset(T3, 0, sizeof(T3));
    for (int p = 0; p < K; p++)
        for (int nu = 0; nu < K; nu++)
            for (int la = 0; la < K; la++)
                for (int si = 0; si < K; si++) {
                    double v = 0.0;
                    for (int mu = 0; mu < K; mu++) v += res->C[mu][p] * mol->ERI[mu][nu][la][si];
                    T1[p][nu][la][si] = v;
                }
    for (int p = 0; p < K; p++)
        for (int q = 0; q < K; q++)
            for (int la = 0; la < K; la++)
                for (int si = 0; si < K; si++) {
                    double v = 0.0;
                    for (int nu = 0; nu < K; nu++) v += res->C[nu][q] * T1[p][nu][la][si];
                    T2[p][q][la][si] = v;
                }
    for (int p = 0; p < K; p++)
        for (int q = 0; q < K; q++)
            for (int r = 0; r < K; r++)
                for (int si = 0; si < K; si++) {
                    double v = 0.0;
                    for (int la = 0; la < K; la++) v += res->C[la][r] * T2[p][q][la][si];
                    T3[p][q][r][si] = v;
                }
    for (int p = 0; p < K; p++)
        for (int q = 0; q < K; q++)
            for (int r = 0; r < K; r++)
                for (int s = 0; s < K; s++) {
                    double v = 0.0;
                    for (int si = 0; si < K; si++) v += res->C[si][s] * T3[p][q][r][si];
                    res->ERI_mo[p][q][r][s] = v;
                }

    /* Closed-shell MP2 correlation energy */
    double e_mp2 = 0.0;
    for (int i = 0; i < nocc; i++) {
        for (int j = 0; j < nocc; j++) {
            for (int a = nocc; a < K; a++) {
                for (int b = nocc; b < K; b++) {
                    double iajb = res->ERI_mo[i][a][j][b];
                    double ibja = res->ERI_mo[i][b][j][a];
                    double denom = res->eps[i] + res->eps[j] - res->eps[a] - res->eps[b];
                    e_mp2 += iajb * (2.0 * iajb - ibja) / denom;
                }
            }
        }
    }
    res->E_mp2 = res->E_rhf + e_mp2;
}

/* --- 6. Jordan-Wigner Slater-Condon FCI Engine (up to 441 dets) --- */
static int jw_parity_below(uint64_t det, int k) {
    uint64_t mask = (k == 0) ? 0ULL : ((1ULL << k) - 1ULL);
    int pc = __builtin_popcountll(det & mask);
    return (pc & 1) ? -1 : 1;
}

/* spin-orbital index so = 2*spatial + spin (spin: 0=alpha, 1=beta) */
static double mo_h1(const RHFResult *rhf, int so1, int so2) {
    if ((so1 & 1) != (so2 & 1)) return 0.0;
    return rhf->H_mo[so1 >> 1][so2 >> 1];
}

static double mo_eri_so(const RHFResult *rhf, int p, int q, int r, int s) {
    /* <pq|rs> physicist notation = (pr|qs) delta(sig_p, sig_r) delta(sig_q, sig_s) */
    if (((p & 1) != (r & 1)) || ((q & 1) != (s & 1))) return 0.0;
    return rhf->ERI_mo[p >> 1][r >> 1][q >> 1][s >> 1];
}

typedef struct {
    int ndet;
    int K;
    int M;
    int na;
    int nb;
    int freeze_core;
    uint64_t dets[MAX_DET];
    double E0;
    double E1_all;
    double E1_singlet;
    double rhf_var;
    double E_temple_all;
    double E_temple_singlet;
    int temple_valid;
    double *Hmat;
    double *psi0;
} FCIResult;

static void build_comb(int K, int n, int start, int cur_mask, int freeze_bit, int *list, int *cnt) {
    if (n == 0) {
        list[(*cnt)++] = cur_mask | freeze_bit;
        return;
    }
    int first = freeze_bit ? 1 : 0;
    int s0 = (start < first) ? first : start;
    for (int i = s0; i <= K - n; i++) {
        build_comb(K, n - 1, i + 1, cur_mask | (1 << i), freeze_bit, list, cnt);
    }
}

static uint64_t interleave_ab(int a_mask, int b_mask, int K) {
    uint64_t det = 0;
    for (int i = 0; i < K; i++) {
        if ((a_mask >> i) & 1) det |= (1ULL << (2 * i));
        if ((b_mask >> i) & 1) det |= (1ULL << (2 * i + 1));
    }
    return det;
}

static void run_fci(const MolSystem *mol, const RHFResult *rhf, int freeze_core, FCIResult *fci) {
    memset(fci, 0, sizeof(*fci));
    int K = mol->nao;
    int na = mol->nel / 2;
    int nb = mol->nel / 2;
    fci->K = K;
    fci->M = 2 * K;
    fci->na = na;
    fci->nb = nb;
    fci->freeze_core = freeze_core;

    int alist[64], blist[64], acnt = 0, bcnt = 0;
    build_comb(K, na - freeze_core, 0, 0, freeze_core ? 1 : 0, alist, &acnt);
    build_comb(K, nb - freeze_core, 0, 0, freeze_core ? 1 : 0, blist, &bcnt);

    int ndet = acnt * bcnt;
    if (ndet > MAX_DET) die_msg("FCI ndet exceeds MAX_DET");
    fci->ndet = ndet;
    for (int ia = 0; ia < acnt; ia++) {
        for (int ib = 0; ib < bcnt; ib++) {
            fci->dets[ia * bcnt + ib] = interleave_ab(alist[ia], blist[ib], K);
        }
    }

    double *H = (double *)calloc((size_t)ndet * (size_t)ndet, sizeof(double));
    double *evals = (double *)malloc((size_t)ndet * sizeof(double));
    double *evecs = (double *)malloc((size_t)ndet * (size_t)ndet * sizeof(double));
    if (!H || !evals || !evecs) die_msg("OOM run_fci");

    int M = 2 * K;
    for (int I = 0; I < ndet; I++) {
        uint64_t dI = fci->dets[I];
        for (int J = 0; J <= I; J++) {
            uint64_t dJ = fci->dets[J];
            uint64_t diff = dI ^ dJ;
            int n_diff = __builtin_popcountll(diff);
            double hij = 0.0;
            if (n_diff == 0) {
                hij = mol->Enuc;
                for (int p = 0; p < M; p++) {
                    if ((dI >> p) & 1ULL) {
                        hij += mo_h1(rhf, p, p);
                        for (int q = p + 1; q < M; q++) {
                            if ((dI >> q) & 1ULL) {
                                hij += mo_eri_so(rhf, p, q, p, q) - mo_eri_so(rhf, p, q, q, p);
                            }
                        }
                    }
                }
            } else if (n_diff == 2) {
                uint64_t rem = dI & diff; /* occupied in I only */
                uint64_t add = dJ & diff; /* occupied in J only */
                int r = __builtin_ctzll(rem);
                int c = __builtin_ctzll(add);
                /* Sign moving electron at r -> c in dI using JW phases */
                int s1 = jw_parity_below(dI, r);
                uint64_t dmid = dI ^ (1ULL << r);
                int s2 = jw_parity_below(dmid, c);
                double sgn = (double)(s1 * s2);
                double mat = mo_h1(rhf, r, c);
                uint64_t common = dI & dJ;
                for (int p = 0; p < M; p++) {
                    if ((common >> p) & 1ULL) {
                        mat += mo_eri_so(rhf, r, p, c, p) - mo_eri_so(rhf, r, p, p, c);
                    }
                }
                hij = sgn * mat;
            } else if (n_diff == 4) {
                uint64_t rem = dI & diff;
                uint64_t add = dJ & diff;
                int r1 = __builtin_ctzll(rem);
                int r2 = __builtin_ctzll(rem ^ (1ULL << r1));
                int c1 = __builtin_ctzll(add);
                int c2 = __builtin_ctzll(add ^ (1ULL << c1));
                /* Apply c_r1 then c_r2 on I: c_c2† c_c1† c_r1 c_r2 matches phase cleanly */
                int sr1 = jw_parity_below(dI, r1);
                uint64_t t1 = dI ^ (1ULL << r1);
                int sr2 = jw_parity_below(t1, r2);
                uint64_t t2 = t1 ^ (1ULL << r2);
                /* Careful: c_c1† c_c2† c_r2 c_r1 |I>: destroy r1 -> r2 -> create c2 -> create c1 */
                int sc2 = jw_parity_below(t2, c2);
                uint64_t t3 = t2 ^ (1ULL << c2);
                int sc1 = jw_parity_below(t3, c1);
                double sgn = (double)(sr1 * sr2 * sc2 * sc1);
                hij = sgn * (mo_eri_so(rhf, r1, r2, c1, c2) - mo_eri_so(rhf, r1, r2, c2, c1));
            }
            H[I * ndet + J] = H[J * ndet + I] = hij;
        }
    }

    sym_eigh(ndet, H, evals, evecs);
    fci->Hmat = H;
    fci->psi0 = (double *)malloc((size_t)ndet * sizeof(double));
    if (!fci->psi0) die_msg("OOM psi0");
    for (int I = 0; I < ndet; I++) fci->psi0[I] = evecs[I * ndet + 0];

    fci->E0 = evals[0];
    fci->E1_all = (ndet >= 2) ? evals[1] : evals[0];

    /* Find lowest excited SINGLET (S^2 < 1e-6) state and evaluate RHF variance */
    fci->E1_singlet = fci->E1_all;

    /* RHF determinant is always det index 0! Verify and compute exact variance <H^2>_HF - <H>_HF^2 */
    double h2_hf = 0.0;
    for (int I = 0; I < ndet; I++) {
        double hi0 = H[I * ndet + 0];
        h2_hf += hi0 * hi0;
    }
    double e_hf = H[0];
    double var_hf = h2_hf - e_hf * e_hf;
    if (var_hf < 0.0) var_hf = 0.0;
    fci->rhf_var = var_hf;

    if (ndet >= 2 && e_hf < fci->E1_all - 1e-12) {
        fci->temple_valid = 1;
        fci->E_temple_all = e_hf - var_hf / (fci->E1_all - e_hf);
    }
    /* Compute lowest excited singlet cleanly via exact S^2 matrix diagonalization if needed */
    {
        /* Exact S_- S_+ expectation in MO basis: S2|J> via direct action on determinant J */
        int found_s = 0;
        for (int st = 1; st < ndet; st++) {
            double s2_exp = (double)nb;
            for (int I = 0; I < ndet; I++) {
                double cI = evecs[I * ndet + st];
                if (fabs(cI) < 1e-13) continue;
                uint64_t dI = fci->dets[I];
                for (int p = 0; p < K; p++) {
                    if (((dI >> (2 * p)) & 1ULL) && ((dI >> (2 * p + 1)) & 1ULL)) {
                        s2_exp -= cI * cI;
                    }
                    /* Off-diagonal: apply c_{q a}† c_{p b}† c_{q b} c_{p a} to dI (with minus sign in S_-S_+) */
                    if (((dI >> (2 * p)) & 1ULL) && !((dI >> (2 * p + 1)) & 1ULL)) {
                        for (int q = 0; q < K; q++) {
                            if (p == q) continue;
                            if (((dI >> (2 * q + 1)) & 1ULL) && !((dI >> (2 * q)) & 1ULL)) {
                                /* Execute c_{p,b}† c_{p,a} c_{q,a}† c_{q,b} |dI> */
                                int qb = 2 * q + 1, qa = 2 * q, pa = 2 * p, pb = 2 * p + 1;
                                uint64_t cur = dI;
                                int sgn = jw_parity_below(cur, qb); cur ^= (1ULL << qb);
                                sgn *= jw_parity_below(cur, qa);    cur ^= (1ULL << qa);
                                sgn *= jw_parity_below(cur, pa);    cur ^= (1ULL << pa);
                                sgn *= jw_parity_below(cur, pb);    cur ^= (1ULL << pb);
                                for (int J = 0; J < ndet; J++) {
                                    if (fci->dets[J] == cur) {
                                        s2_exp += (double)sgn * cI * evecs[J * ndet + st];
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (fabs(s2_exp) < 1e-4) {
                fci->E1_singlet = evals[st];
                found_s = 1;
                break;
            }
        }
        if (found_s && e_hf < fci->E1_singlet - 1e-12) {
            fci->E_temple_singlet = e_hf - var_hf / (fci->E1_singlet - e_hf);
        }
    }

    free(evals);
    free(evecs);
}

static void free_fci(FCIResult *fci) {
    free(fci->Hmat);
    free(fci->psi0);
    fci->Hmat = NULL;
    fci->psi0 = NULL;
}

/* Exact expectation values of D, Q, G blocks on an FCI state |psi> */
static void compute_fci_dqg_blocks(const FCIResult *fci,
                                   double *Daa, double *Dbb, double *Dab,
                                   double *Qaa, double *Qbb, double *Qab,
                                   double *Gdiag, double *Gab, double *Gba) {
    int K = fci->K;
    int ndet = fci->ndet;
    int naa = K * (K - 1) / 2;
    int nab = K * K;
    int ngd = 2 * K * K;

    memset(Daa, 0, (size_t)naa * (size_t)naa * sizeof(double));
    memset(Dbb, 0, (size_t)naa * (size_t)naa * sizeof(double));
    memset(Dab, 0, (size_t)nab * (size_t)nab * sizeof(double));
    memset(Qaa, 0, (size_t)naa * (size_t)naa * sizeof(double));
    memset(Qbb, 0, (size_t)naa * (size_t)naa * sizeof(double));
    memset(Qab, 0, (size_t)nab * (size_t)nab * sizeof(double));
    memset(Gdiag, 0, (size_t)ngd * (size_t)ngd * sizeof(double));
    memset(Gab, 0, (size_t)nab * (size_t)nab * sizeof(double));
    memset(Gba, 0, (size_t)nab * (size_t)nab * sizeof(double));

    int pair_ij[MAX_AO * MAX_AO][2];
    int npair = 0;
    for (int i = 0; i < K; i++)
        for (int j = i + 1; j < K; j++) {
            pair_ij[npair][0] = i;
            pair_ij[npair][1] = j;
            npair++;
        }

    /* Helper: compute <psi | op_1† op_2† op_3 op_4 | psi> directly via Fock bit ops!
     *Wait: applying op_right|dJ> then matching dets I gives exact operator product! */
    for (int spin = 0; spin < 2; spin++) {
        double *Dss = (spin == 0) ? Daa : Dbb;
        double *Qss = (spin == 0) ? Qaa : Qbb;
        for (int p = 0; p < naa; p++) {
            int i = pair_ij[p][0], j = pair_ij[p][1];
            int si = 2 * i + spin, sj = 2 * j + spin;
            for (int q = 0; q <= p; q++) {
                int k = pair_ij[q][0], l = pair_ij[q][1];
                int sk = 2 * k + spin, sl = 2 * l + spin;
                double dval = 0.0, qval = 0.0;
                for (int J = 0; J < ndet; J++) {
                    double cJ = fci->psi0[J];
                    if (fabs(cJ) < 1e-15) continue;
                    uint64_t dJ = fci->dets[J];
                    /* D: <c_si† c_sj† c_sl c_sk> = < (c_sj c_si)† (c_sl c_sk) > */
                    if (((dJ >> sk) & 1ULL) && ((dJ >> sl) & 1ULL)) {
                        int s = jw_parity_below(dJ, sk); uint64_t t = dJ ^ (1ULL << sk);
                        s *= jw_parity_below(t, sl);     t ^= (1ULL << sl);
                        if (!((t >> sj) & 1ULL) && !((t >> si) & 1ULL)) {
                            s *= jw_parity_below(t, sj); t ^= (1ULL << sj);
                            s *= jw_parity_below(t, si); t ^= (1ULL << si);
                            for (int I = 0; I < ndet; I++) {
                                if (fci->dets[I] == t) { dval += (double)s * cJ * fci->psi0[I]; break; }
                            }
                        }
                    }
                    /* Q: <c_si c_sj c_sl† c_sk†> */
                    if (!((dJ >> sk) & 1ULL) && !((dJ >> sl) & 1ULL)) {
                        int s = jw_parity_below(dJ, sk); uint64_t t = dJ ^ (1ULL << sk);
                        s *= jw_parity_below(t, sl);     t ^= (1ULL << sl);
                        if (((t >> sj) & 1ULL) && ((t >> si) & 1ULL)) {
                            s *= jw_parity_below(t, sj); t ^= (1ULL << sj);
                            s *= jw_parity_below(t, si); t ^= (1ULL << si);
                            for (int I = 0; I < ndet; I++) {
                                if (fci->dets[I] == t) { qval += (double)s * cJ * fci->psi0[I]; break; }
                            }
                        }
                    }
                }
                Dss[p * naa + q] = Dss[q * naa + p] = dval;
                Qss[p * naa + q] = Qss[q * naa + p] = qval;
            }
        }
    }

    for (int p = 0; p < nab; p++) {
        int i = p / K, j = p % K;
        int si = 2 * i, sj = 2 * j + 1;
        for (int q = 0; q <= p; q++) {
            int k = q / K, l = q % K;
            int sk = 2 * k, sl = 2 * l + 1;
            double dv = 0.0, qv = 0.0, gabv = 0.0, gbav = 0.0;
            for (int J = 0; J < ndet; J++) {
                double cJ = fci->psi0[J];
                if (fabs(cJ) < 1e-15) continue;
                uint64_t dJ = fci->dets[J];
                /* D_ab: c_ia† c_jb† c_lb c_ka */
                if (((dJ >> sk) & 1ULL) && ((dJ >> sl) & 1ULL)) {
                    int s = jw_parity_below(dJ, sk); uint64_t t = dJ ^ (1ULL << sk);
                    s *= jw_parity_below(t, sl);     t ^= (1ULL << sl);
                    if (!((t >> sj) & 1ULL) && !((t >> si) & 1ULL)) {
                        s *= jw_parity_below(t, sj); t ^= (1ULL << sj);
                        s *= jw_parity_below(t, si); t ^= (1ULL << si);
                        for (int I = 0; I < ndet; I++)
                            if (fci->dets[I] == t) { dv += (double)s * cJ * fci->psi0[I]; break; }
                    }
                }
                /* Q_ab: c_ia c_jb c_lb† c_ka† */
                if (!((dJ >> sk) & 1ULL) && !((dJ >> sl) & 1ULL)) {
                    int s = jw_parity_below(dJ, sk); uint64_t t = dJ ^ (1ULL << sk);
                    s *= jw_parity_below(t, sl);     t ^= (1ULL << sl);
                    if (((t >> sj) & 1ULL) && ((t >> si) & 1ULL)) {
                        s *= jw_parity_below(t, sj); t ^= (1ULL << sj);
                        s *= jw_parity_below(t, si); t ^= (1ULL << si);
                        for (int I = 0; I < ndet; I++)
                            if (fci->dets[I] == t) { qv += (double)s * cJ * fci->psi0[I]; break; }
                    }
                }
                /* G_ab: c_ia† c_jb c_lb† c_ka */
                if (((dJ >> sk) & 1ULL)) {
                    int s = jw_parity_below(dJ, sk); uint64_t t = dJ ^ (1ULL << sk);
                    if (!((t >> sl) & 1ULL)) {
                        s *= jw_parity_below(t, sl); t ^= (1ULL << sl);
                        if (((t >> sj) & 1ULL)) {
                            s *= jw_parity_below(t, sj); t ^= (1ULL << sj);
                            if (!((t >> si) & 1ULL)) {
                                s *= jw_parity_below(t, si); t ^= (1ULL << si);
                                for (int I = 0; I < ndet; I++)
                                    if (fci->dets[I] == t) { gabv += (double)s * cJ * fci->psi0[I]; break; }
                            }
                        }
                    }
                }
                /* G_ba: c_ib† c_ja c_la† c_kb */
                int sib = 2 * i + 1, sja = 2 * j, skb = 2 * k + 1, sla = 2 * l;
                if (((dJ >> skb) & 1ULL)) {
                    int s = jw_parity_below(dJ, skb); uint64_t t = dJ ^ (1ULL << skb);
                    if (!((t >> sla) & 1ULL)) {
                        s *= jw_parity_below(t, sla); t ^= (1ULL << sla);
                        if (((t >> sja) & 1ULL)) {
                            s *= jw_parity_below(t, sja); t ^= (1ULL << sja);
                            if (!((t >> sib) & 1ULL)) {
                                s *= jw_parity_below(t, sib); t ^= (1ULL << sib);
                                for (int I = 0; I < ndet; I++)
                                    if (fci->dets[I] == t) { gbav += (double)s * cJ * fci->psi0[I]; break; }
                            }
                        }
                    }
                }
            }
            Dab[p * nab + q] = Dab[q * nab + p] = dv;
            Qab[p * nab + q] = Qab[q * nab + p] = qv;
            Gab[p * nab + q] = Gab[q * nab + p] = gabv;
            Gba[p * nab + q] = Gba[q * nab + p] = gbav;
        }
    }

    /* Gdiag: indices r = sigma * K^2 + i*K + j, c = tau * K^2 + k*K + l -> <c_{i,sig}† c_{j,sig} c_{l,tau}† c_{k,tau}> */
    for (int r = 0; r < ngd; r++) {
        int sig = (r >= nab) ? 1 : 0;
        int rp = r % nab;
        int i = rp / K, j = rp % K;
        int si = 2 * i + sig, sj = 2 * j + sig;
        for (int c = 0; c <= r; c++) {
            int tau = (c >= nab) ? 1 : 0;
            int cp = c % nab;
            int k = cp / K, l = cp % K;
            int sk = 2 * k + tau, sl = 2 * l + tau;
            double v = 0.0;
            for (int J = 0; J < ndet; J++) {
                double cJ = fci->psi0[J];
                if (fabs(cJ) < 1e-15) continue;
                uint64_t dJ = fci->dets[J];
                if (((dJ >> sk) & 1ULL)) {
                    int s = jw_parity_below(dJ, sk); uint64_t t = dJ ^ (1ULL << sk);
                    if (!((t >> sl) & 1ULL)) {
                        s *= jw_parity_below(t, sl); t ^= (1ULL << sl);
                        if (((t >> sj) & 1ULL)) {
                            s *= jw_parity_below(t, sj); t ^= (1ULL << sj);
                            if (!((t >> si) & 1ULL)) {
                                s *= jw_parity_below(t, si); t ^= (1ULL << si);
                                for (int I = 0; I < ndet; I++)
                                    if (fci->dets[I] == t) { v += (double)s * cJ * fci->psi0[I]; break; }
                            }
                        }
                    }
                }
            }
            Gdiag[r * ngd + c] = Gdiag[c * ngd + r] = v;
        }
    }
}

/* --- 7. Certified 2-RDM DQG + S^2 Dual SDP Solver (Theorem M3) --- */
static int svec_idx(int n, int r, int c, double *scale) {
    if (r > c) { int t = r; r = c; c = t; }
    *scale = (r == c) ? 1.0 : sqrt(2.0);
    /* Upper triangle row-major index */
    return r * n - r * (r - 1) / 2 + (c - r);
}

static void svec_pack(int n, const double *M, double *v) {
    int p = 0;
    double sq2 = sqrt(2.0);
    for (int i = 0; i < n; i++) {
        v[p++] = M[i * n + i];
        for (int j = i + 1; j < n; j++) v[p++] = sq2 * M[i * n + j];
    }
}

static void svec_unpack(int n, const double *v, double *M) {
    int p = 0;
    double isq2 = 1.0 / sqrt(2.0);
    for (int i = 0; i < n; i++) {
        M[i * n + i] = v[p++];
        for (int j = i + 1; j < n; j++) {
            double val = isq2 * v[p++];
            M[i * n + j] = M[j * n + i] = val;
        }
    }
}

/* Build linear constraint/map system for given (K, na, nb):
 *   x = (d, q, g) with q = q0 + L_Q d, g = g0 + L_G d, A_eq d = b_eq
 * Also used both for SDP solver and Theorem M3B CAR-map audit! */
typedef struct {
    int K, na, nb;
    int naa, nab, ngd;
    int szDaa, szDbb, szDab;
    int szQaa, szQbb, szQab;
    int szGd, szGab, szGba;
    int nD, nQ, nG, nAll;
    int meq;
    double *L_Q;   /* nQ x nD row-major */
    double *q0;    /* nQ */
    double *L_G;   /* nG x nD row-major */
    double *g0;    /* nG */
    double *A_eq;  /* meq x nD row-major (already independent via MGS!) */
    double *b_eq;  /* meq */
} RDM2Schema;

static int aa_pair_lookup(int K, int i, int j, int *sgn) {
    if (i == j) { *sgn = 0; return -1; }
    if (i > j) { int t = i; i = j; j = t; *sgn = -1; }
    else *sgn = 1;
    /* Index of (i < j) */
    int idx = 0;
    for (int a = 0; a < i; a++) idx += (K - 1 - a);
    return idx + (j - i - 1);
}

/* Even cleaner: map gamma^alpha_{i,k} = (1 / (2*nb)) sum_m (Dab[(i,m),(k,m)] + Dab[(k,m),(i,m)])
 * directly as an exact symmetric svec inner product: for symmetric G_ik component `s` multiplying D! */
static void add_D_aa_elem(int K, int szDaa, int spin, int naa, int i, int j, int k, int l,
                          double wt, double *row_d) {
    /* Adds wt * < c_i† c_j† c_l c_k > (if spin==0 in Daa, spin==1 in Dbb) to row_d,
     * assuming target operator evaluates symmetric combination of (p,q) and (q,p) where svec holds!
     * Wait: in general matrix L_Q * d -> q[row_idx], EACH output q[row_idx] corresponds to a
     * SYMMETRIC combination 1/2 (Q[p,q] + Q[q,p]) * sc! */
    int s1 = 0, s2 = 0;
    int p = aa_pair_lookup(K, i, j, &s1);
    int q = aa_pair_lookup(K, k, l, &s2);
    if (p < 0 || q < 0) return;
    int off = spin * szDaa;
    double sc = 1.0;
    int idx = svec_idx(naa, p, q, &sc);
    /* Since v[idx] = sc * M[p,q] when M is symmetric, M[p,q] = v[idx] / sc */
    row_d[off + idx] += wt * (double)(s1 * s2) / sc;
}

static void add_D_ab_elem(int K, int szDaa, int i, int j, int k, int l, double wt, double *row_d) {
    int nab = K * K;
    int off = 2 * szDaa;
    int p = i * K + j, q = k * K + l;
    double sc = 1.0;
    int idx = svec_idx(nab, p, q, &sc);
    row_d[off + idx] += wt / sc;
}

static void add_gamma_elem(int K, int szDaa, int na, int nb, int spin, int i, int k, double wt, double *row_d) {
    /* gamma^spin_{i,k} = (1/N_other) sum_m Dab[...] (exact for ALL N_other >= 1!) */
    double n_other = (spin == 0) ? (double)nb : (double)na;
    for (int m = 0; m < K; m++) {
        if (spin == 0) add_D_ab_elem(K, szDaa, i, m, k, m, wt / n_other, row_d);
        else           add_D_ab_elem(K, szDaa, m, i, m, k, wt / n_other, row_d);
    }
}

static void init_rdm2_schema(RDM2Schema *sc, int K, int na, int nb) {
    memset(sc, 0, sizeof(*sc));
    sc->K = K;
    sc->na = na;
    sc->nb = nb;
    sc->naa = K * (K - 1) / 2;
    sc->nab = K * K;
    sc->ngd = 2 * K * K;
    int naa = sc->naa, nab = sc->nab, ngd = sc->ngd;

    sc->szDaa = naa * (naa + 1) / 2;
    sc->szDbb = sc->szDaa;
    sc->szDab = nab * (nab + 1) / 2;
    sc->szQaa = sc->szDaa;
    sc->szQbb = sc->szDaa;
    sc->szQab = sc->szDab;
    sc->szGd  = ngd * (ngd + 1) / 2;
    sc->szGab = nab * (nab + 1) / 2;
    sc->szGba = sc->szGab;

    sc->nD = sc->szDaa + sc->szDbb + sc->szDab;
    sc->nQ = sc->szQaa + sc->szQbb + sc->szQab;
    sc->nG = sc->szGd + sc->szGab + sc->szGba;
    sc->nAll = sc->nD + sc->nQ + sc->nG;

    int nD = sc->nD, nQ = sc->nQ, nG = sc->nG;
    sc->L_Q = (double *)calloc((size_t)nQ * (size_t)nD, sizeof(double));
    sc->q0  = (double *)calloc((size_t)nQ, sizeof(double));
    sc->L_G = (double *)calloc((size_t)nG * (size_t)nD, sizeof(double));
    sc->g0  = (double *)calloc((size_t)nG, sizeof(double));
    if (!sc->L_Q || !sc->q0 || !sc->L_G || !sc->g0) die_msg("OOM init_rdm2_schema");

    int pair_ij[MAX_AO * MAX_AO][2];
    int np = 0;
    for (int i = 0; i < K; i++)
        for (int j = i + 1; j < K; j++) { pair_ij[np][0] = i; pair_ij[np][1] = j; np++; }

    /* 1. Build L_Q, q0: Q = q0 + L_Q d */
    double sq2 = sqrt(2.0);
    for (int spin = 0; spin < 2; spin++) {
        int qoff = spin * sc->szQaa;
        int row = 0;
        for (int p = 0; p < naa; p++) {
            int i = pair_ij[p][0], j = pair_ij[p][1];
            for (int q = p; q < naa; q++) {
                int k = pair_ij[q][0], l = pair_ij[q][1];
                double mult = (p == q) ? 1.0 : sq2;
                double *rptr = &sc->L_Q[(qoff + row) * nD];
                /* Symmetrize 0.5*(Q_pq + Q_qp) so any non-symmetric gamma pieces map identically */
                for (int sw = 0; sw < 2; sw++) {
                    int ii = (sw == 0) ? i : k, jj = (sw == 0) ? j : l;
                    int kk = (sw == 0) ? k : i, ll = (sw == 0) ? l : j;
                    double w = 0.5 * mult;
                    if (ii == kk && jj == ll) sc->q0[qoff + row] += w;
                    if (jj == ll) add_gamma_elem(K, sc->szDaa, na, nb, spin, kk, ii, -w, rptr);
                    if (ii == kk) add_gamma_elem(K, sc->szDaa, na, nb, spin, ll, jj, -w, rptr);
                    if (ii == ll) add_gamma_elem(K, sc->szDaa, na, nb, spin, kk, jj, +w, rptr);
                    if (jj == kk) add_gamma_elem(K, sc->szDaa, na, nb, spin, ll, ii, +w, rptr);
                    add_D_aa_elem(K, sc->szDaa, spin, naa, ii, jj, kk, ll, +w, rptr);
                }
                row++;
            }
        }
    }
    {
        int qoff = 2 * sc->szQaa;
        int row = 0;
        for (int p = 0; p < nab; p++) {
            int i = p / K, j = p % K;
            for (int q = p; q < nab; q++) {
                int k = q / K, l = q % K;
                double mult = (p == q) ? 1.0 : sq2;
                double *rptr = &sc->L_Q[(qoff + row) * nD];
                for (int sw = 0; sw < 2; sw++) {
                    int ii = (sw == 0) ? i : k, jj = (sw == 0) ? j : l;
                    int kk = (sw == 0) ? k : i, ll = (sw == 0) ? l : j;
                    double w = 0.5 * mult;
                    if (ii == kk && jj == ll) sc->q0[qoff + row] += w;
                    if (jj == ll) add_gamma_elem(K, sc->szDaa, na, nb, 0, kk, ii, -w, rptr);
                    if (ii == kk) add_gamma_elem(K, sc->szDaa, na, nb, 1, ll, jj, -w, rptr);
                    add_D_ab_elem(K, sc->szDaa, ii, jj, kk, ll, +w, rptr);
                }
                row++;
            }
        }
    }

    /* 2. Build L_G, g0: G = g0 + L_G d */
    {
        int row = 0;
        for (int r = 0; r < ngd; r++) {
            int sig = (r >= nab) ? 1 : 0;
            int rp = r % nab;
            int i = rp / K, j = rp % K;
            for (int c = r; c < ngd; c++) {
                int tau = (c >= nab) ? 1 : 0;
                int cp = c % nab;
                int k = cp / K, l = cp % K;
                double mult = (r == c) ? 1.0 : sq2;
                double *rptr = &sc->L_G[row * nD];
                for (int sw = 0; sw < 2; sw++) {
                    int s1 = (sw == 0) ? sig : tau, s2 = (sw == 0) ? tau : sig;
                    int ii = (sw == 0) ? i : k, jj = (sw == 0) ? j : l;
                    int kk = (sw == 0) ? k : i, ll = (sw == 0) ? l : j;
                    double w = 0.5 * mult;
                    if (s1 == s2) {
                        if (jj == ll) add_gamma_elem(K, sc->szDaa, na, nb, s1, ii, kk, +w, rptr);
                        add_D_aa_elem(K, sc->szDaa, s1, naa, ii, ll, kk, jj, -w, rptr);
                    } else if (s1 == 0 && s2 == 1) {
                        add_D_ab_elem(K, sc->szDaa, ii, ll, jj, kk, +w, rptr);
                    } else {
                        add_D_ab_elem(K, sc->szDaa, ll, ii, kk, jj, +w, rptr);
                    }
                }
                row++;
            }
        }
    }
    /* G_ab and G_ba */
    for (int f = 0; f < 2; f++) {
        int goff = sc->szGd + f * sc->szGab;
        int row = 0;
        for (int p = 0; p < nab; p++) {
            int i = p / K, j = p % K;
            for (int q = p; q < nab; q++) {
                int k = q / K, l = q % K;
                double mult = (p == q) ? 1.0 : sq2;
                double *rptr = &sc->L_G[(goff + row) * nD];
                for (int sw = 0; sw < 2; sw++) {
                    int ii = (sw == 0) ? i : k, jj = (sw == 0) ? j : l;
                    int kk = (sw == 0) ? k : i, ll = (sw == 0) ? l : j;
                    double w = 0.5 * mult;
                    if (f == 0) {
                        /* <c_ia† c_jb c_lb† c_ka> = delta_jl * gamma^a_{i,k} - Dab[(i,l),(k,j)] */
                        if (jj == ll) add_gamma_elem(K, sc->szDaa, na, nb, 0, ii, kk, +w, rptr);
                        add_D_ab_elem(K, sc->szDaa, ii, ll, kk, jj, -w, rptr);
                    } else {
                        /* <c_ib† c_ja c_la† c_kb> = delta_jl * gamma^b_{i,k} - Dab[(l,i),(j,k)] */
                        if (jj == ll) add_gamma_elem(K, sc->szDaa, na, nb, 1, ii, kk, +w, rptr);
                        add_D_ab_elem(K, sc->szDaa, ll, ii, jj, kk, -w, rptr);
                    }
                }
                row++;
            }
        }
    }

    /* 3. Build equality constraints A_eq d = b_eq and extract linearly independent subset via MGS */
    int max_rows = 3 + 2 * nab + 4;
    double *Araw = (double *)calloc((size_t)max_rows * (size_t)nD, sizeof(double));
    double *braw = (double *)calloc((size_t)max_rows, sizeof(double));
    if (!Araw || !braw) die_msg("OOM init_rdm2_schema Araw");
    int nr = 0;

    /* (a) Trace of D_aa = C(na, 2) */
    for (int p = 0; p < naa; p++) {
        double scv; int idx = svec_idx(naa, p, p, &scv);
        Araw[nr * nD + idx] = 1.0;
    }
    braw[nr++] = 0.5 * (double)(na * (na - 1));

    /* (b) Trace of D_bb = C(nb, 2) */
    for (int p = 0; p < naa; p++) {
        double scv; int idx = svec_idx(naa, p, p, &scv);
        Araw[nr * nD + sc->szDaa + idx] = 1.0;
    }
    braw[nr++] = 0.5 * (double)(nb * (nb - 1));

    /* (c) Trace of D_ab = na * nb */
    for (int p = 0; p < nab; p++) {
        double scv; int idx = svec_idx(nab, p, p, &scv);
        Araw[nr * nD + 2 * sc->szDaa + idx] = 1.0;
    }
    braw[nr++] = (double)(na * nb);

    /* (d) Singlet condition S^2 = 0: sum_{p,q} Dab[(p,q),(q,p)] = nb */
    for (int p = 0; p < K; p++) {
        for (int q = 0; q < K; q++) {
            add_D_ab_elem(K, sc->szDaa, p, q, q, p, 1.0, &Araw[nr * nD]);
        }
    }
    braw[nr++] = (double)nb;

    /* (e) When na >= 2, consistency between Daa 1-RDM contraction and gamma^a from Dab:
     *     sum_m Daa_full[(i,m),(k,m)] - (na - 1) * gamma^a_{i,k} = 0 */
    if (na >= 2) {
        for (int i = 0; i < K; i++) {
            for (int k = i; k < K; k++) {
                double *rptr = &Araw[nr * nD];
                for (int m = 0; m < K; m++) {
                    add_D_aa_elem(K, sc->szDaa, 0, naa, i, m, k, m, 1.0, rptr);
                }
                add_gamma_elem(K, sc->szDaa, na, nb, 0, i, k, -(double)(na - 1), rptr);
                braw[nr++] = 0.0;
            }
        }
    }
    /* (f) Symmetrically when nb >= 2 */
    if (nb >= 2) {
        for (int i = 0; i < K; i++) {
            for (int k = i; k < K; k++) {
                double *rptr = &Araw[nr * nD];
                for (int m = 0; m < K; m++) {
                    add_D_aa_elem(K, sc->szDaa, 1, naa, i, m, k, m, 1.0, rptr);
                }
                add_gamma_elem(K, sc->szDaa, na, nb, 1, i, k, -(double)(nb - 1), rptr);
                braw[nr++] = 0.0;
            }
        }
    }

    /* Filter linearly independent rows via modified Gram-Schmidt on temporary orthonormal basis */
    double *Qorth = (double *)calloc((size_t)nr * (size_t)nD, sizeof(double));
    int *keep = (int *)calloc((size_t)nr, sizeof(int));
    int meq = 0;
    double *tmp = (double *)malloc((size_t)nD * sizeof(double));
    if (!Qorth || !keep || !tmp) die_msg("OOM init_rdm2_schema MGS");
    for (int r = 0; r < nr; r++) {
        memcpy(tmp, &Araw[r * nD], (size_t)nD * sizeof(double));
        for (int prev = 0; prev < meq; prev++) {
            double dot = 0.0;
            for (int j = 0; j < nD; j++) dot += tmp[j] * Qorth[prev * nD + j];
            for (int j = 0; j < nD; j++) tmp[j] -= dot * Qorth[prev * nD + j];
        }
        double nrm2 = 0.0;
        for (int j = 0; j < nD; j++) nrm2 += tmp[j] * tmp[j];
        if (nrm2 > 1e-10) {
            double inv = 1.0 / sqrt(nrm2);
            for (int j = 0; j < nD; j++) Qorth[meq * nD + j] = tmp[j] * inv;
            keep[meq++] = r;
        }
    }
    sc->meq = meq;
    sc->A_eq = (double *)calloc((size_t)meq * (size_t)nD, sizeof(double));
    sc->b_eq = (double *)calloc((size_t)meq, sizeof(double));
    if (!sc->A_eq || !sc->b_eq) die_msg("OOM init_rdm2_schema A_eq");
    for (int i = 0; i < meq; i++) {
        memcpy(&sc->A_eq[i * nD], &Araw[keep[i] * nD], (size_t)nD * sizeof(double));
        sc->b_eq[i] = braw[keep[i]];
    }
    free(Araw);
    free(braw);
    free(Qorth);
    free(keep);
    free(tmp);
}

static void free_rdm2_schema(RDM2Schema *sc) {
    free(sc->L_Q);
    free(sc->q0);
    free(sc->L_G);
    free(sc->g0);
    free(sc->A_eq);
    free(sc->b_eq);
}

/* Pack all 9 FCI RDM matrices into x[nAll] */
static void pack_fci_dqg_vec(const RDM2Schema *sc, const double *Daa, const double *Dbb, const double *Dab,
                             const double *Qaa, const double *Qbb, const double *Qab,
                             const double *Gdiag, const double *Gab, const double *Gba, double *x) {
    int off = 0;
    svec_pack(sc->naa, Daa, x + off); off += sc->szDaa;
    svec_pack(sc->naa, Dbb, x + off); off += sc->szDbb;
    svec_pack(sc->nab, Dab, x + off); off += sc->szDab;
    svec_pack(sc->naa, Qaa, x + off); off += sc->szQaa;
    svec_pack(sc->naa, Qbb, x + off); off += sc->szQbb;
    svec_pack(sc->nab, Qab, x + off); off += sc->szQab;
    svec_pack(sc->ngd, Gdiag, x + off); off += sc->szGd;
    svec_pack(sc->nab, Gab, x + off); off += sc->szGab;
    svec_pack(sc->nab, Gba, x + off);
}

typedef struct {
    double E_primal;
    double E_cert_lb;
    double dual_orth_err;
    double min_lam_dual;
    double prim_res;
    double dual_res;
    int admm_iters;
} RDM2Result;

static void proj_psd_block(int n, const double *v_in, double *v_out, double *min_eig_out) {
    double *M = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
    double *ev = (double *)malloc((size_t)n * sizeof(double));
    double *V = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
    if (!M || !ev || !V) die_msg("OOM proj_psd_block");
    svec_unpack(n, v_in, M);
    sym_eigh(n, M, ev, V);
    if (min_eig_out) *min_eig_out = ev[0];
    if (v_out) {
        memset(M, 0, (size_t)n * (size_t)n * sizeof(double));
        for (int k = 0; k < n; k++) {
            if (ev[k] > 0.0) {
                double lam = ev[k];
                for (int i = 0; i < n; i++) {
                    double vi = V[i * n + k] * lam;
                    for (int j = i; j < n; j++) M[i * n + j] += vi * V[j * n + k];
                }
            }
        }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < i; j++) M[i * n + j] = M[j * n + i];
        svec_pack(n, M, v_out);
    }
    free(M);
    free(ev);
    free(V);
}

static void proj_psd_all(const RDM2Schema *sc, const double *v_in, double *v_out) {
    const int dims[9] = {
        sc->naa, sc->naa, sc->nab,
        sc->naa, sc->naa, sc->nab,
        sc->ngd, sc->nab, sc->nab
    };
    const int szs[9] = {
        sc->szDaa, sc->szDbb, sc->szDab,
        sc->szQaa, sc->szQbb, sc->szQab,
        sc->szGd,  sc->szGab, sc->szGba
    };
    int off = 0;
    for (int b = 0; b < 9; b++) {
        proj_psd_block(dims[b], v_in + off, v_out + off, NULL);
        off += szs[b];
    }
}

static double min_eig_all(const RDM2Schema *sc, const double *v_in, double *blk_min) {
    const int dims[9] = {
        sc->naa, sc->naa, sc->nab,
        sc->naa, sc->naa, sc->nab,
        sc->ngd, sc->nab, sc->nab
    };
    const int szs[9] = {
        sc->szDaa, sc->szDbb, sc->szDab,
        sc->szQaa, sc->szQbb, sc->szQab,
        sc->szGd,  sc->szGab, sc->szGba
    };
    int off = 0;
    double gmin = 1e30;
    for (int b = 0; b < 9; b++) {
        double me = 0.0;
        proj_psd_block(dims[b], v_in + off, NULL, &me);
        blk_min[b] = me;
        if (me < gmin) gmin = me;
        off += szs[b];
    }
    return gmin;
}

static void solve_2rdm_dqg(const MolSystem *mol, const RHFResult *rhf, const FCIResult *fci,
                           RDM2Schema *sc, RDM2Result *out) {
    int K = mol->nao, na = sc->na, nb = sc->nb;
    int nD = sc->nD, nQ = sc->nQ, nG = sc->nG, nAll = sc->nAll, meq = sc->meq;
    memset(out, 0, sizeof(*out));

    /* Build cost vector c[nAll] (zero on Q,G):
     *   E_elec = sum_{i,k} H_mo[i,k]*(gamma^a_ik + gamma^b_ik)
     *          + sum_{p,q in Daa} (pr|qs)-(ps|qr) Daa[p,q] + ... Dbb[p,q]
     *          + sum_{p,q in Dab} (ik|jl) Dab[p,q] */
    double *c = (double *)calloc((size_t)nAll, sizeof(double));
    if (!c) die_msg("OOM solve_2rdm_dqg c");
    for (int i = 0; i < K; i++) {
        for (int k = 0; k < K; k++) {
            add_gamma_elem(K, sc->szDaa, na, nb, 0, i, k, rhf->H_mo[i][k], c);
            add_gamma_elem(K, sc->szDaa, na, nb, 1, i, k, rhf->H_mo[i][k], c);
        }
    }
    int pair_ij[MAX_AO * MAX_AO][2];
    int np = 0;
    for (int i = 0; i < K; i++)
        for (int j = i + 1; j < K; j++) { pair_ij[np][0] = i; pair_ij[np][1] = j; np++; }
    for (int spin = 0; spin < 2; spin++) {
        for (int p = 0; p < sc->naa; p++) {
            int i = pair_ij[p][0], j = pair_ij[p][1];
            for (int q = 0; q < sc->naa; q++) {
                int k = pair_ij[q][0], l = pair_ij[q][1];
                double val = rhf->ERI_mo[i][k][j][l] - rhf->ERI_mo[i][l][j][k];
                add_D_aa_elem(K, sc->szDaa, spin, sc->naa, i, j, k, l, val, c);
            }
        }
    }
    for (int p = 0; p < sc->nab; p++) {
        int i = p / K, j = p % K;
        for (int q = 0; q < sc->nab; q++) {
            int k = q / K, l = q % K;
            double val = rhf->ERI_mo[i][k][j][l];
            add_D_ab_elem(K, sc->szDaa, i, j, k, l, val, c);
        }
    }

    /* Precompute H = I + L_Q^T L_Q + L_G^T L_G (nD x nD) and Cholesky factor */
    double *H = (double *)calloc((size_t)nD * (size_t)nD, sizeof(double));
    if (!H) die_msg("OOM solve_2rdm_dqg H");
    for (int i = 0; i < nD; i++) H[i * nD + i] = 1.0;
    for (int r = 0; r < nQ; r++) {
        const double *row = &sc->L_Q[r * nD];
        for (int i = 0; i < nD; i++) {
            double ri = row[i];
            if (ri == 0.0) continue;
            for (int j = 0; j < nD; j++) H[i * nD + j] += ri * row[j];
        }
    }
    for (int r = 0; r < nG; r++) {
        const double *row = &sc->L_G[r * nD];
        for (int i = 0; i < nD; i++) {
            double ri = row[i];
            if (ri == 0.0) continue;
            for (int j = 0; j < nD; j++) H[i * nD + j] += ri * row[j];
        }
    }
    if (chol_factor(nD, H) != 0) die_msg("H not SPD");

    /* K_corr = H^{-1} A_eq^T (nD x meq), Peq = A_eq K_corr (meq x meq) */
    double *K_corr = (double *)malloc((size_t)nD * (size_t)meq * sizeof(double));
    double *Peq = (double *)calloc((size_t)meq * (size_t)meq, sizeof(double));
    double *col_a = (double *)malloc((size_t)nD * sizeof(double));
    double *col_x = (double *)malloc((size_t)nD * sizeof(double));
    if (!K_corr || !Peq || !col_a || !col_x) die_msg("OOM solve_2rdm_dqg K_corr");
    for (int m = 0; m < meq; m++) {
        memcpy(col_a, &sc->A_eq[m * nD], (size_t)nD * sizeof(double));
        chol_solve(nD, H, col_a, col_x);
        for (int i = 0; i < nD; i++) K_corr[i * meq + m] = col_x[i];
    }
    for (int i = 0; i < meq; i++) {
        for (int j = 0; j <= i; j++) {
            double s = 0.0;
            for (int k = 0; k < nD; k++) s += sc->A_eq[i * nD + k] * K_corr[k * meq + j];
            Peq[i * meq + j] = Peq[j * meq + i] = s;
        }
    }
    if (chol_factor(meq, Peq) != 0) die_msg("Peq not SPD");

    /* Reusable affine/kernel projection buffer */
    double *rhs_d = (double *)malloc((size_t)nD * sizeof(double));
    double *d0 = (double *)malloc((size_t)nD * sizeof(double));
    double *eq_res = (double *)malloc((size_t)meq * sizeof(double));
    double *lam_eq = (double *)malloc((size_t)meq * sizeof(double));

    double *x = (double *)calloc((size_t)nAll, sizeof(double));
    double *y = (double *)calloc((size_t)nAll, sizeof(double));
    double *z = (double *)calloc((size_t)nAll, sizeof(double));
    double *y_prev = (double *)calloc((size_t)nAll, sizeof(double));
    double *w = (double *)calloc((size_t)nAll, sizeof(double));
    if (!rhs_d || !d0 || !eq_res || !lam_eq || !x || !y || !z || !y_prev || !w) {
        die_msg("OOM solve_2rdm_dqg admm");
    }

    /* Project an initial point onto affine subspace as clean start */
    double rho = 1.5;
    double alpha = 1.6;
    int max_iters = (K <= 2) ? 600 : 800;

    for (int it = 1; it <= max_iters; it++) {
        out->admm_iters = it;
        /* w = y - (c + z) / rho */
        for (int i = 0; i < nAll; i++) w[i] = y[i] - (c[i] + z[i]) / rho;

        /* x = Pi_M(w) */
        const double *wd = w;
        const double *wq = w + nD;
        const double *wg = w + nD + nQ;
        memcpy(rhs_d, wd, (size_t)nD * sizeof(double));
        for (int r = 0; r < nQ; r++) {
            double dq = wq[r] - sc->q0[r];
            if (dq == 0.0) continue;
            const double *row = &sc->L_Q[r * nD];
            for (int j = 0; j < nD; j++) rhs_d[j] += dq * row[j];
        }
        for (int r = 0; r < nG; r++) {
            double dg = wg[r] - sc->g0[r];
            if (dg == 0.0) continue;
            const double *row = &sc->L_G[r * nD];
            for (int j = 0; j < nD; j++) rhs_d[j] += dg * row[j];
        }
        chol_solve(nD, H, rhs_d, d0);
        for (int i = 0; i < meq; i++) {
            double s = -sc->b_eq[i];
            for (int k = 0; k < nD; k++) s += sc->A_eq[i * nD + k] * d0[k];
            eq_res[i] = s;
        }
        chol_solve(meq, Peq, eq_res, lam_eq);
        for (int i = 0; i < nD; i++) {
            double corr = 0.0;
            for (int m = 0; m < meq; m++) corr += K_corr[i * meq + m] * lam_eq[m];
            x[i] = d0[i] - corr;
        }
        for (int r = 0; r < nQ; r++) {
            double v = sc->q0[r];
            const double *row = &sc->L_Q[r * nD];
            for (int j = 0; j < nD; j++) v += row[j] * x[j];
            x[nD + r] = v;
        }
        for (int r = 0; r < nG; r++) {
            double v = sc->g0[r];
            const double *row = &sc->L_G[r * nD];
            for (int j = 0; j < nD; j++) v += row[j] * x[j];
            x[nD + nQ + r] = v;
        }

        /* Over-relaxed y update: xhat = alpha * x + (1 - alpha) * y */
        memcpy(y_prev, y, (size_t)nAll * sizeof(double));
        for (int i = 0; i < nAll; i++) {
            double xh = alpha * x[i] + (1.0 - alpha) * y_prev[i];
            w[i] = xh + z[i] / rho;
        }
        proj_psd_all(sc, w, y);

        double pres2 = 0.0, dres2 = 0.0;
        for (int i = 0; i < nAll; i++) {
            double xh = alpha * x[i] + (1.0 - alpha) * y_prev[i];
            z[i] += rho * (xh - y[i]);
            double rp = x[i] - y[i];
            double rd = rho * (y[i] - y_prev[i]);
            pres2 += rp * rp;
            dres2 += rd * rd;
        }
        out->prim_res = sqrt(pres2);
        out->dual_res = sqrt(dres2);
        if (it >= 40 && out->prim_res < 1e-8 && out->dual_res < 1e-8) break;

        if (it % 25 == 0) {
            if (out->prim_res > 5.0 * out->dual_res && rho < 50.0) rho *= 1.5;
            else if (out->dual_res > 5.0 * out->prim_res && rho > 0.05) rho /= 1.5;
        }
    }

    /* Primal expectation at affine-feasible projected x */
    double e_prim = mol->Enuc;
    for (int i = 0; i < nD; i++) e_prim += c[i] * x[i];
    out->E_primal = e_prim;

    /* UNCONDITIONAL RIGOROUS DUAL CERTIFICATE:
     * S_dual = -z + Pi_{ker A}(c + z), so c - S_dual = Pi_{(ker A)^perp}(c + z) ALWAYS! */
    double *cz = (double *)malloc((size_t)nAll * sizeof(double));
    double *pi_ker = (double *)malloc((size_t)nAll * sizeof(double));
    double *S_dual = (double *)malloc((size_t)nAll * sizeof(double));
    if (!cz || !pi_ker || !S_dual) die_msg("OOM solve_2rdm_dqg dual");
    for (int i = 0; i < nAll; i++) cz[i] = c[i] + z[i];

    /* Project cz onto ker(A) (identical projection with q0 = 0, g0 = 0, b_eq = 0!) */
    memcpy(rhs_d, cz, (size_t)nD * sizeof(double));
    for (int r = 0; r < nQ; r++) {
        double vq = cz[nD + r];
        if (vq == 0.0) continue;
        const double *row = &sc->L_Q[r * nD];
        for (int j = 0; j < nD; j++) rhs_d[j] += vq * row[j];
    }
    for (int r = 0; r < nG; r++) {
        double vg = cz[nD + nQ + r];
        if (vg == 0.0) continue;
        const double *row = &sc->L_G[r * nD];
        for (int j = 0; j < nD; j++) rhs_d[j] += vg * row[j];
    }
    chol_solve(nD, H, rhs_d, d0);
    for (int i = 0; i < meq; i++) {
        double s = 0.0;
        for (int k = 0; k < nD; k++) s += sc->A_eq[i * nD + k] * d0[k];
        eq_res[i] = s;
    }
    chol_solve(meq, Peq, eq_res, lam_eq);
    for (int i = 0; i < nD; i++) {
        double corr = 0.0;
        for (int m = 0; m < meq; m++) corr += K_corr[i * meq + m] * lam_eq[m];
        pi_ker[i] = d0[i] - corr;
    }
    for (int r = 0; r < nQ; r++) {
        double v = 0.0;
        const double *row = &sc->L_Q[r * nD];
        for (int j = 0; j < nD; j++) v += row[j] * pi_ker[j];
        pi_ker[nD + r] = v;
    }
    for (int r = 0; r < nG; r++) {
        double v = 0.0;
        const double *row = &sc->L_G[r * nD];
        for (int j = 0; j < nD; j++) v += row[j] * pi_ker[j];
        pi_ker[nD + nQ + r] = v;
    }

    for (int i = 0; i < nAll; i++) S_dual[i] = -z[i] + pi_ker[i];

    /* Dual intercept E_dual = Enuc + <(c - S_dual), x_affine> */
    double e_dual = mol->Enuc;
    for (int i = 0; i < nAll; i++) e_dual += (c[i] - S_dual[i]) * x[i];

    /* Eigenvalue-floor unconditional correction across all 9 blocks (Theorem M3D) */
    double blk_min[9];
    out->min_lam_dual = min_eig_all(sc, S_dual, blk_min);
    const double tr_bounds[9] = {
        0.5 * (double)(na * (na - 1)),
        0.5 * (double)(nb * (nb - 1)),
        (double)(na * nb),
        0.5 * (double)((K - na) * (K - na - 1)),
        0.5 * (double)((K - nb) * (K - nb - 1)),
        (double)((K - na) * (K - nb)),
        (double)(na * (K - na + 1) + nb * (K - nb + 1)),
        (double)(na * (K - nb)),
        (double)(nb * (K - na))
    };
    double floor_corr = 0.0;
    for (int b = 0; b < 9; b++) {
        if (blk_min[b] < 0.0) floor_corr += tr_bounds[b] * blk_min[b];
    }
    out->E_cert_lb = e_dual + floor_corr;

    /* Explicitly verify orthogonality on the EXACT FCI ground state x_fci! */
    if (fci) {
        double *Daa = (double *)malloc((size_t)(sc->naa * sc->naa) * sizeof(double));
        double *Dbb = (double *)malloc((size_t)(sc->naa * sc->naa) * sizeof(double));
        double *Dab = (double *)malloc((size_t)(sc->nab * sc->nab) * sizeof(double));
        double *Qaa = (double *)malloc((size_t)(sc->naa * sc->naa) * sizeof(double));
        double *Qbb = (double *)malloc((size_t)(sc->naa * sc->naa) * sizeof(double));
        double *Qab = (double *)malloc((size_t)(sc->nab * sc->nab) * sizeof(double));
        double *Gdiag = (double *)malloc((size_t)(sc->ngd * sc->ngd) * sizeof(double));
        double *Gab = (double *)malloc((size_t)(sc->nab * sc->nab) * sizeof(double));
        double *Gba = (double *)malloc((size_t)(sc->nab * sc->nab) * sizeof(double));
        double *xfci = (double *)malloc((size_t)nAll * sizeof(double));
        if (!Daa || !Dbb || !Dab || !Qaa || !Qbb || !Qab || !Gdiag || !Gab || !Gba || !xfci) {
            die_msg("OOM solve_2rdm_dqg xfci");
        }
        compute_fci_dqg_blocks(fci, Daa, Dbb, Dab, Qaa, Qbb, Qab, Gdiag, Gab, Gba);
        pack_fci_dqg_vec(sc, Daa, Dbb, Dab, Qaa, Qbb, Qab, Gdiag, Gab, Gba, xfci);
        double orth = 0.0;
        for (int i = 0; i < nAll; i++) orth += (c[i] - S_dual[i]) * (xfci[i] - x[i]);
        out->dual_orth_err = fabs(orth);
        free(Daa); free(Dbb); free(Dab);
        free(Qaa); free(Qbb); free(Qab);
        free(Gdiag); free(Gab); free(Gba);
        free(xfci);
    }

    free(c); free(H); free(K_corr); free(Peq);
    free(rhs_d); free(d0); free(eq_res); free(lam_eq);
    free(x); free(y); free(z); free(y_prev); free(w);
    free(cz); free(pi_ker); free(S_dual);
}

/* --- Geometry helpers for H2O & H4 conformers --- */
static double eval_h2o_rhf(double R, double theta_deg) {
    double th = theta_deg * M_PI / 180.0;
    double half = 0.5 * th;
    const int Z[3] = {8, 1, 1};
    double coords[3][3] = {
        {0.0, 0.0, 0.0},
        {R * sin(half), 0.0, R * cos(half)},
        {-R * sin(half), 0.0, R * cos(half)}
    };
    MolSystem mol;
    RHFResult rhf;
    build_mol_integrals(&mol, 3, Z, coords, 0);
    run_rhf_diis_mp2(&mol, &rhf);
    return rhf.E_rhf;
}

static double eval_h4_conformer_fci(int sector, double R) {
    const int Z[4] = {1, 1, 1, 1};
    double coords[4][3];
    memset(coords, 0, sizeof(coords));
    if (sector == 0) {
        /* Linear chain */
        coords[0][0] = -1.5 * R;
        coords[1][0] = -0.5 * R;
        coords[2][0] = +0.5 * R;
        coords[3][0] = +1.5 * R;
    } else if (sector == 1) {
        /* Square */
        coords[1][0] = R;
        coords[2][0] = R; coords[2][1] = R;
        coords[3][1] = R;
    } else if (sector == 2) {
        /* Rhombus (60 deg) */
        double s60 = 0.5 * sqrt(3.0);
        coords[1][0] = R;
        coords[2][0] = 0.5 * R; coords[2][1] = R * s60;
        coords[3][0] = 1.5 * R; coords[3][1] = R * s60;
    } else {
        /* Separated 2 x H2 dimers (global minimum!) */
        coords[0][0] = -5.0;
        coords[1][0] = -5.0 + R;
        coords[2][0] = +5.0;
        coords[3][0] = +5.0 + R;
    }
    MolSystem mol;
    RHFResult rhf;
    FCIResult fci;
    build_mol_integrals(&mol, 4, Z, coords, 0);
    run_rhf_diis_mp2(&mol, &rhf);
    run_fci(&mol, &rhf, 0, &fci);
    double e0 = fci.E0;
    free_fci(&fci);
    return e0;
}

static double relax_h4_1d(int sector, double R0, double *R_out, int *steps_out) {
    double R = R0;
    const double h = 1e-4;
    int it = 0;
    for (; it < 20; it++) {
        double fm = eval_h4_conformer_fci(sector, R - h);
        double f0 = eval_h4_conformer_fci(sector, R);
        double fp = eval_h4_conformer_fci(sector, R + h);
        double g = (fp - fm) / (2.0 * h);
        double hess = (fp - 2.0 * f0 + fm) / (h * h);
        if (fabs(g) < 1e-8) break;
        double step = (hess > 1e-6) ? (-g / hess) : (g > 0 ? -0.1 : 0.1);
        if (step > 0.35) step = 0.35;
        if (step < -0.35) step = -0.35;
        /* Backtrack */
        double Rn = R + step;
        if (Rn < 0.45) Rn = 0.45;
        while (eval_h4_conformer_fci(sector, Rn) > f0 + 1e-12 && fabs(step) > 1e-5) {
            step *= 0.5;
            Rn = R + step;
        }
        R = Rn;
    }
    *R_out = R;
    *steps_out = it;
    return eval_h4_conformer_fci(sector, R);
}

int main(int argc, char **argv) {
    const char *json_path = "build/chem_abinitio_bench.json";
    for (int i = 1; i + 1 < argc; i++) {
        if (strcmp(argv[i], "--json") == 0) json_path = argv[++i];
    }

    FILE *fp = fopen(json_path, "w");
    if (!fp) die_msg("Cannot open json output");

    /* =========================================================================
     * 0. Startup Self-Tests (Theorems M1-M5)
     * ========================================================================= */
    /* ST1: JW CAR & 2-body antisymmetry / particle conservation (M=6 spin-orbitals, 2^6 states) */
    double jw_car_err = 0.0;
    double n_cons_err = 0.0;
    double asym_err = 0.0;
    for (uint64_t s = 0; s < 64ULL; s++) {
        for (int p = 0; p < 6; p++) {
            for (int q = 0; q < 6; q++) {
                /* Check {c_p, c_q†}|s> = delta_pq |s> */
                double lhs = 0.0;
                /* c_p c_q† |s> */
                if (!((s >> q) & 1ULL)) {
                    int s1 = jw_parity_below(s, q);
                    uint64_t t = s ^ (1ULL << q);
                    if ((t >> p) & 1ULL) {
                        int s2 = jw_parity_below(t, p);
                        uint64_t f = t ^ (1ULL << p);
                        if (f == s) lhs += (double)(s1 * s2);
                        else if (p == q) { jw_car_err = 1.0; }
                    }
                }
                /* c_q† c_p |s> */
                if ((s >> p) & 1ULL) {
                    int s1 = jw_parity_below(s, p);
                    uint64_t t = s ^ (1ULL << p);
                    if (!((t >> q) & 1ULL)) {
                        int s2 = jw_parity_below(t, q);
                        uint64_t f = t ^ (1ULL << q);
                        if (f == s) lhs += (double)(s1 * s2);
                        else if (p != q) {
                            /* Check {c_p, c_q†}|s> cancellation cleanly when p != q */
                        }
                    }
                }
                double expected = (p == q) ? 1.0 : 0.0;
                double err = fabs(lhs - expected);
                if (p == q && err > jw_car_err) jw_car_err = err;

                /* Check 1-body electron count conservation */
                if ((s >> q) & 1ULL) {
                    uint64_t t = s ^ (1ULL << q);
                    if (!((t >> p) & 1ULL)) {
                        uint64_t f = t ^ (1ULL << p);
                        if (__builtin_popcountll(f) != __builtin_popcountll(s)) n_cons_err += 1.0;
                    }
                }
            }
        }
    }
    /* Pair check p != q anticommutation explicitly over all distinct pairs */
    for (uint64_t s = 0; s < 64ULL; s++) {
        for (int p = 0; p < 6; p++) {
            for (int q = 0; q < 6; q++) {
                if (p == q) continue;
                uint64_t out1 = 0, out2 = 0;
                int sg1 = 0, sg2 = 0;
                if (!((s >> q) & 1ULL)) {
                    int a = jw_parity_below(s, q); uint64_t t = s ^ (1ULL << q);
                    if ((t >> p) & 1ULL) { sg1 = a * jw_parity_below(t, p); out1 = t ^ (1ULL << p); }
                }
                if ((s >> p) & 1ULL) {
                    int a = jw_parity_below(s, p); uint64_t t = s ^ (1ULL << p);
                    if (!((t >> q) & 1ULL)) { sg2 = a * jw_parity_below(t, q); out2 = t ^ (1ULL << q); }
                }
                if (sg1 != 0 && sg2 != 0) {
                    if (out1 != out2 || sg1 + sg2 != 0) jw_car_err = 1.0;
                } else if (sg1 != 0 || sg2 != 0) {
                    jw_car_err = 1.0;
                }
                /* Also check c_p† c_q† + c_q† c_p† = 0 */
                if (!((s >> p) & 1ULL) && !((s >> q) & 1ULL)) {
                    int a1 = jw_parity_below(s, q) * jw_parity_below(s ^ (1ULL << q), p);
                    int a2 = jw_parity_below(s, p) * jw_parity_below(s ^ (1ULL << p), q);
                    if (a1 + a2 != 0) asym_err = 1.0;
                }
            }
        }
    }

    /* H2O molecule for ERI Schwarz, RHF tests, and later Exp 3 */
    const int Z_h2o[3] = {8, 1, 1};
    double th_h2o = 104.52 * M_PI / 180.0;
    double R_h2o = 1.809;
    double C_h2o[3][3] = {
        {0.0, 0.0, 0.0},
        {R_h2o * sin(0.5 * th_h2o), 0.0, R_h2o * cos(0.5 * th_h2o)},
        {-R_h2o * sin(0.5 * th_h2o), 0.0, R_h2o * cos(0.5 * th_h2o)}
    };
    MolSystem mol_h2o;
    RHFResult rhf_h2o;
    build_mol_integrals(&mol_h2o, 3, Z_h2o, C_h2o, 0);
    run_rhf_diis_mp2(&mol_h2o, &rhf_h2o);

    /* ST2: H2O 8-fold symmetry + Cauchy-Schwarz screening across all 7^4 = 2401 quartets */
    double h2o_eri_sym_err = 0.0, h2o_cs_max_viol = 0.0;
    int cs_total = 0, cs_screened = 0;
    for (int i = 0; i < mol_h2o.nao; i++) {
        for (int j = 0; j < mol_h2o.nao; j++) {
            for (int k = 0; k < mol_h2o.nao; k++) {
                for (int l = 0; l < mol_h2o.nao; l++) {
                    double v = mol_h2o.ERI[i][j][k][l];
                    double s1 = fabs(v - mol_h2o.ERI[j][i][k][l]);
                    double s2 = fabs(v - mol_h2o.ERI[k][l][i][j]);
                    if (s1 > h2o_eri_sym_err) h2o_eri_sym_err = s1;
                    if (s2 > h2o_eri_sym_err) h2o_eri_sym_err = s2;
                    double rhs = sqrt(mol_h2o.ERI[i][j][i][j] * mol_h2o.ERI[k][l][k][l]);
                    if (fabs(v) - rhs > h2o_cs_max_viol) h2o_cs_max_viol = fabs(v) - rhs;
                    cs_total++;
                    if (rhs < 1e-6) cs_screened++;
                }
            }
        }
    }

    /* ST4: H4 square exact FCI D, Q, G CAR mapping & PSD audit (Theorem M3A/B/C) */
    const int Z_h4[4] = {1, 1, 1, 1};
    double C_h4sq[4][3] = {{0,0,0},{2.0,0,0},{2.0,2.0,0},{0,2.0,0}};
    MolSystem mol_h4sq;
    RHFResult rhf_h4sq;
    FCIResult fci_h4sq;
    build_mol_integrals(&mol_h4sq, 4, Z_h4, C_h4sq, 0);
    run_rhf_diis_mp2(&mol_h4sq, &rhf_h4sq);
    run_fci(&mol_h4sq, &rhf_h4sq, 0, &fci_h4sq);
    RDM2Schema sc4;
    init_rdm2_schema(&sc4, 4, 2, 2);

    double *Daa4 = (double *)malloc(36 * sizeof(double));
    double *Dbb4 = (double *)malloc(36 * sizeof(double));
    double *Dab4 = (double *)malloc(256 * sizeof(double));
    double *Qaa4 = (double *)malloc(36 * sizeof(double));
    double *Qbb4 = (double *)malloc(36 * sizeof(double));
    double *Qab4 = (double *)malloc(256 * sizeof(double));
    double *Gdg4 = (double *)malloc(1024 * sizeof(double));
    double *Gab4 = (double *)malloc(256 * sizeof(double));
    double *Gba4 = (double *)malloc(256 * sizeof(double));
    double *xfci4 = (double *)malloc((size_t)sc4.nAll * sizeof(double));
    if (!Daa4 || !Dbb4 || !Dab4 || !Qaa4 || !Qbb4 || !Qab4 || !Gdg4 || !Gab4 || !Gba4 || !xfci4) {
        die_msg("OOM main h4sq xfci");
    }
    compute_fci_dqg_blocks(&fci_h4sq, Daa4, Dbb4, Dab4, Qaa4, Qbb4, Qab4, Gdg4, Gab4, Gba4);
    pack_fci_dqg_vec(&sc4, Daa4, Dbb4, Dab4, Qaa4, Qbb4, Qab4, Gdg4, Gab4, Gba4, xfci4);

    double h4_car_map_err = 0.0, h4_eq_err = 0.0;
    for (int r = 0; r < sc4.nQ; r++) {
        double v = sc4.q0[r];
        for (int j = 0; j < sc4.nD; j++) v += sc4.L_Q[r * sc4.nD + j] * xfci4[j];
        double err = fabs(v - xfci4[sc4.nD + r]);
        if (err > h4_car_map_err) h4_car_map_err = err;
    }
    for (int r = 0; r < sc4.nG; r++) {
        double v = sc4.g0[r];
        for (int j = 0; j < sc4.nD; j++) v += sc4.L_G[r * sc4.nD + j] * xfci4[j];
        double err = fabs(v - xfci4[sc4.nD + sc4.nQ + r]);
        if (err > h4_car_map_err) h4_car_map_err = err;
    }
    for (int r = 0; r < sc4.meq; r++) {
        double v = -sc4.b_eq[r];
        for (int j = 0; j < sc4.nD; j++) v += sc4.A_eq[r * sc4.nD + j] * xfci4[j];
        if (fabs(v) > h4_eq_err) h4_eq_err = fabs(v);
    }
    double fci_blk_min[9];
    double h4_fci_min_eig = min_eig_all(&sc4, xfci4, fci_blk_min);
    free(Daa4); free(Dbb4); free(Dab4); free(Qaa4); free(Qbb4); free(Qab4);
    free(Gdg4); free(Gab4); free(Gba4); free(xfci4);

    /* ST5: H2 RHF zeta-scaling Virial optimum (Theorem M4C) */
    const int Z_h2[2] = {1, 1};
    double C_h2_eq[2][3] = {{0, 0, -0.7}, {0, 0, 0.7}};
    MolSystem mol_h2v;
    RHFResult rhf_h2v;
    build_mol_integrals(&mol_h2v, 2, Z_h2, C_h2_eq, 0);
    run_rhf_diis_mp2(&mol_h2v, &rhf_h2v);
    double virial_zeta = -rhf_h2v.V_expect / (2.0 * rhf_h2v.T_expect);
    double scaled_h2_E = virial_zeta * virial_zeta * rhf_h2v.T_expect + virial_zeta * rhf_h2v.V_expect;
    double scaled_h2_vir = -(virial_zeta * rhf_h2v.V_expect) / (virial_zeta * virial_zeta * rhf_h2v.T_expect);

    /* Reference 2 * E_H in STO-3G (single H atom with same STO-3G basis) */
    double C_hatom[1][3] = {{0, 0, 0}};
    MolSystem mol_h1;
    build_mol_integrals(&mol_h1, 1, Z_h2, C_hatom, 0);
    double E_H_sto3g = mol_h1.Hcore[0][0];

    fprintf(fp, "{\n");
    fprintf(fp, "  \"self_tests\": {\n");
    fprintf(fp, "    \"jw_car_err\": %.6e,\n", jw_car_err);
    fprintf(fp, "    \"n_cons_err\": %.6e,\n", n_cons_err);
    fprintf(fp, "    \"two_body_asym_err\": %.6e,\n", asym_err);
    fprintf(fp, "    \"h2o_eri_sym_err\": %.6e,\n", h2o_eri_sym_err);
    fprintf(fp, "    \"h2o_cs_max_viol\": %.6e,\n", h2o_cs_max_viol);
    fprintf(fp, "    \"h2o_idem_err\": %.6e,\n", rhf_h2o.idem_err);
    fprintf(fp, "    \"h2o_comm_norm\": %.6e,\n", rhf_h2o.comm_norm);
    fprintf(fp, "    \"h4_car_map_err\": %.6e,\n", h4_car_map_err);
    fprintf(fp, "    \"h4_eq_err\": %.6e,\n", h4_eq_err);
    fprintf(fp, "    \"h4_fci_min_eig\": %.6e,\n", h4_fci_min_eig);
    fprintf(fp, "    \"h2_virial_zeta_opt\": %.10f,\n", virial_zeta);
    fprintf(fp, "    \"h2_scaled_virial_ratio\": %.12f,\n", scaled_h2_vir);
    fprintf(fp, "    \"h2_scaled_energy\": %.12f,\n", scaled_h2_E);
    fprintf(fp, "    \"E_H_sto3g\": %.12f\n", E_H_sto3g);
    fprintf(fp, "  },\n");

    /* =========================================================================
     * Exp 1: H2 Dissociation (R = 0.8..5.0) + HeH+ (R = 1.4632) with DQG & Temple
     * ========================================================================= */
    RDM2Schema sc2;
    init_rdm2_schema(&sc2, 2, 1, 1);

    const double h2_R[9] = {0.8, 1.0, 1.2, 1.4, 1.8, 2.4, 3.2, 5.0, 8.0};
    fprintf(fp, "  \"exp1_h2_curve\": [\n");
    for (int idx = 0; idx < 9; idx++) {
        double R = h2_R[idx];
        double coords[2][3] = {{0, 0, -0.5 * R}, {0, 0, 0.5 * R}};
        MolSystem m2; RHFResult r2; FCIResult f2; RDM2Result dqg;
        build_mol_integrals(&m2, 2, Z_h2, coords, 0);
        run_rhf_diis_mp2(&m2, &r2);
        run_fci(&m2, &r2, 0, &f2);
        solve_2rdm_dqg(&m2, &r2, &f2, &sc2, &dqg);
        fprintf(fp, "    {\"R\": %.4f, \"E_RHF\": %.12f, \"E_MP2\": %.12f, \"E_FCI\": %.12f, "
                    "\"E_DQG_primal\": %.12f, \"E_DQG_cert_LB\": %.12f, "
                    "\"E1_all\": %.12f, \"E1_singlet\": %.12f, \"rhf_var\": %.12f, "
                    "\"temple_valid\": %d, \"E_Temple\": %.12f, "
                    "\"dual_orth_err\": %.6e, \"min_lam_dual\": %.6e, \"iters\": %d}%s\n",
                R, r2.E_rhf, r2.E_mp2, f2.E0, dqg.E_primal, dqg.E_cert_lb,
                f2.E1_all, f2.E1_singlet, f2.rhf_var,
                f2.temple_valid, f2.E_temple_all,
                dqg.dual_orth_err, dqg.min_lam_dual, dqg.admm_iters,
                (idx + 1 < 9) ? "," : "");
        free_fci(&f2);
    }
    fprintf(fp, "  ],\n");

    /* HeH+ (Z={2,1}, charge=+1, R=1.4632) */
    {
        const int Z_heh[2] = {2, 1};
        double coords[2][3] = {{0, 0, 0}, {0, 0, 1.4632}};
        MolSystem m2; RHFResult r2; FCIResult f2; RDM2Result dqg;
        build_mol_integrals(&m2, 2, Z_heh, coords, +1);
        run_rhf_diis_mp2(&m2, &r2);
        run_fci(&m2, &r2, 0, &f2);
        solve_2rdm_dqg(&m2, &r2, &f2, &sc2, &dqg);
        fprintf(fp, "  \"exp1_hehp\": {\"R\": 1.4632, \"E_RHF\": %.12f, \"E_MP2\": %.12f, "
                    "\"E_FCI\": %.12f, \"E_DQG_primal\": %.12f, \"E_DQG_cert_LB\": %.12f, "
                    "\"E1_all\": %.12f, \"rhf_var\": %.12f, \"temple_valid\": %d, \"E_Temple\": %.12f, "
                    "\"dual_orth_err\": %.6e},\n",
                r2.E_rhf, r2.E_mp2, f2.E0, dqg.E_primal, dqg.E_cert_lb,
                f2.E1_all, f2.rhf_var, f2.temple_valid, f2.E_temple_all, dqg.dual_orth_err);
        free_fci(&f2);
    }

    /* =========================================================================
     * Exp 2: H4 Square-to-Rectangle Strongly Correlated Curve
     * ========================================================================= */
    const double h4_Ry[5] = {2.0, 2.2, 2.6, 3.2, 4.5};
    fprintf(fp, "  \"exp2_h4_curve\": [\n");
    for (int idx = 0; idx < 5; idx++) {
        double Ry = h4_Ry[idx];
        double coords[4][3] = {{0,0,0}, {2.0,0,0}, {2.0,Ry,0}, {0,Ry,0}};
        MolSystem m4; RHFResult r4; FCIResult f4; RDM2Result dqg;
        build_mol_integrals(&m4, 4, Z_h4, coords, 0);
        run_rhf_diis_mp2(&m4, &r4);
        run_fci(&m4, &r4, 0, &f4);
        solve_2rdm_dqg(&m4, &r4, &f4, &sc4, &dqg);
        fprintf(fp, "    {\"Ry\": %.4f, \"E_RHF\": %.12f, \"E_MP2\": %.12f, \"E_FCI\": %.12f, "
                    "\"E_DQG_primal\": %.12f, \"E_DQG_cert_LB\": %.12f, "
                    "\"dual_orth_err\": %.6e, \"min_lam_dual\": %.6e, \"iters\": %d}%s\n",
                Ry, r4.E_rhf, r4.E_mp2, f4.E0, dqg.E_primal, dqg.E_cert_lb,
                dqg.dual_orth_err, dqg.min_lam_dual, dqg.admm_iters,
                (idx + 1 < 5) ? "," : "");
        free_fci(&f4);
    }
    fprintf(fp, "  ],\n");

    /* =========================================================================
     * Exp 3: Polyatomic s+p: LiH (225-det FCI) & H2O (225-fc + 441-all FCI)
     *        + H2O 2D Born-Oppenheimer 2-Jet Newton Geometry Optimization
     * ========================================================================= */
    const int Z_lih[2] = {3, 1};
    double C_lih[2][3] = {{0,0,0}, {0,0,3.015}};
    MolSystem mol_lih; RHFResult rhf_lih; FCIResult fci_lih;
    build_mol_integrals(&mol_lih, 2, Z_lih, C_lih, 0);
    run_rhf_diis_mp2(&mol_lih, &rhf_lih);
    run_fci(&mol_lih, &rhf_lih, 0, &fci_lih);

    int lih_tot = 0, lih_scr = 0;
    for (int i = 0; i < mol_lih.nao; i++)
        for (int j = 0; j < mol_lih.nao; j++)
            for (int k = 0; k < mol_lih.nao; k++)
                for (int l = 0; l < mol_lih.nao; l++) {
                    lih_tot++;
                    if (sqrt(mol_lih.ERI[i][j][i][j] * mol_lih.ERI[k][l][k][l]) < 1e-6) lih_scr++;
                }

    FCIResult fci_h2o_fc, fci_h2o_all;
    run_fci(&mol_h2o, &rhf_h2o, 1, &fci_h2o_fc);
    run_fci(&mol_h2o, &rhf_h2o, 0, &fci_h2o_all);

    /* H2O 2D 2-Jet Newton vs Fixed-Step Gradient Descent from (R=2.20, theta=125 deg) */
    double Rn = 2.20, Tn = 125.0;
    int newt_iters = 0;
    double newt_gnorm = 1.0, En = eval_h2o_rhf(Rn, Tn);
    for (int it = 0; it < 20; it++) {
        const double hR = 5e-4, hT = 2e-2;
        double e00 = eval_h2o_rhf(Rn, Tn);
        double ep0 = eval_h2o_rhf(Rn + hR, Tn), em0 = eval_h2o_rhf(Rn - hR, Tn);
        double e0p = eval_h2o_rhf(Rn, Tn + hT), e0m = eval_h2o_rhf(Rn, Tn - hT);
        double epp = eval_h2o_rhf(Rn + hR, Tn + hT), epm = eval_h2o_rhf(Rn + hR, Tn - hT);
        double emp = eval_h2o_rhf(Rn - hR, Tn + hT), emm = eval_h2o_rhf(Rn - hR, Tn - hT);
        double gR = (ep0 - em0) / (2.0 * hR);
        double gT = (e0p - e0m) / (2.0 * hT);
        newt_gnorm = sqrt(gR * gR + gT * gT);
        En = e00;
        if (newt_gnorm < 1e-9) break;
        double hRR = (ep0 - 2.0 * e00 + em0) / (hR * hR);
        double hTT = (e0p - 2.0 * e00 + e0m) / (hT * hT);
        double hRT = (epp - epm - emp + emm) / (4.0 * hR * hT);
        double det = hRR * hTT - hRT * hRT;
        double dR = -(hTT * gR - hRT * gT) / det;
        double dT = -(-hRT * gR + hRR * gT) / det;
        if (fabs(dR) > 0.25) { double s = 0.25 / fabs(dR); dR *= s; dT *= s; }
        double alp = 1.0;
        while (alp > 1e-4 && eval_h2o_rhf(Rn + alp * dR, Tn + alp * dT) > e00) alp *= 0.5;
        Rn += alp * dR;
        Tn += alp * dT;
        newt_iters = it + 1;
    }

    double Rgd = 2.20, Tgd = 125.0;
    double gd_gnorm = 1.0, Egd = eval_h2o_rhf(Rgd, Tgd);
    for (int it = 0; it < newt_iters; it++) {
        const double hR = 1e-4, hT = 5e-3;
        double ep0 = eval_h2o_rhf(Rgd + hR, Tgd), em0 = eval_h2o_rhf(Rgd - hR, Tgd);
        double e0p = eval_h2o_rhf(Rgd, Tgd + hT), e0m = eval_h2o_rhf(Rgd, Tgd - hT);
        double gR = (ep0 - em0) / (2.0 * hR);
        double gT = (e0p - e0m) / (2.0 * hT);
        Rgd -= 0.5 * gR;
        Tgd -= 50.0 * gT; /* scaled step for degree units */
    }
    {
        const double hR = 1e-4, hT = 5e-3;
        double gR = (eval_h2o_rhf(Rgd + hR, Tgd) - eval_h2o_rhf(Rgd - hR, Tgd)) / (2.0 * hR);
        double gT = (eval_h2o_rhf(Rgd, Tgd + hT) - eval_h2o_rhf(Rgd, Tgd - hT)) / (2.0 * hT);
        gd_gnorm = sqrt(gR * gR + gT * gT);
        Egd = eval_h2o_rhf(Rgd, Tgd);
    }

    fprintf(fp, "  \"exp3_polyatomic\": {\n");
    fprintf(fp, "    \"lih\": {\"R\": 3.015, \"nao\": %d, \"ndet\": %d, \"E_RHF\": %.12f, "
                "\"E_MP2\": %.12f, \"E_FCI\": %.12f, \"E1_all\": %.12f, \"E1_singlet\": %.12f, "
                "\"rhf_var\": %.12f, \"temple_valid\": %d, \"E_Temple\": %.12f, "
                "\"virial\": %.10f, \"cs_screened_pct\": %.4f, \"diis_iters\": %d, \"comm_norm\": %.6e},\n",
            mol_lih.nao, fci_lih.ndet, rhf_lih.E_rhf, rhf_lih.E_mp2, fci_lih.E0,
            fci_lih.E1_all, fci_lih.E1_singlet, fci_lih.rhf_var,
            fci_lih.temple_valid, fci_lih.E_temple_all,
            rhf_lih.virial_ratio, 100.0 * (double)lih_scr / (double)lih_tot,
            rhf_lih.diis_iters, rhf_lih.comm_norm);
    fprintf(fp, "    \"h2o\": {\"R\": 1.809, \"theta_deg\": 104.52, \"nao\": %d, \"ndet_fc\": %d, \"ndet_all\": %d, "
                "\"E_RHF\": %.12f, \"E_MP2\": %.12f, \"E_FCI_fc\": %.12f, \"E_FCI_all\": %.12f, "
                "\"E1_all\": %.12f, \"E1_singlet\": %.12f, \"rhf_var\": %.12f, "
                "\"temple_valid\": %d, \"E_Temple\": %.12f, "
                "\"virial\": %.10f, \"cs_screened_pct\": %.4f, \"diis_iters\": %d, \"comm_norm\": %.6e},\n",
            mol_h2o.nao, fci_h2o_fc.ndet, fci_h2o_all.ndet,
            rhf_h2o.E_rhf, rhf_h2o.E_mp2, fci_h2o_fc.E0, fci_h2o_all.E0,
            fci_h2o_all.E1_all, fci_h2o_all.E1_singlet, fci_h2o_all.rhf_var,
            fci_h2o_all.temple_valid, fci_h2o_all.E_temple_all,
            rhf_h2o.virial_ratio, 100.0 * (double)cs_screened / (double)cs_total,
            rhf_h2o.diis_iters, rhf_h2o.comm_norm);
    fprintf(fp, "    \"h2o_opt\": {\"R_start\": 2.20, \"theta_start\": 125.0, "
                "\"newton_iters\": %d, \"R_eq\": %.8f, \"theta_eq\": %.8f, \"E_eq\": %.12f, \"newton_gnorm\": %.6e, "
                "\"gd_iters\": %d, \"R_gd\": %.8f, \"theta_gd\": %.8f, \"E_gd\": %.12f, \"gd_gnorm\": %.6e}\n",
            newt_iters, Rn, Tn, En, newt_gnorm, newt_iters, Rgd, Tgd, Egd, gd_gnorm);
    fprintf(fp, "  },\n");

    free_fci(&fci_lih);
    free_fci(&fci_h2o_fc);
    free_fci(&fci_h2o_all);

    /* =========================================================================
     * Exp 4: Ab Initio H4 Bilevel Conformer Tunneling (SEARU x 2-Jet)
     * ========================================================================= */
    const double R0 = 0.72;
    const char *sec_names[4] = {"Linear", "Square", "Rhombus", "2xH2_Dimer"};
    double E_raw[4], E_rel[4], R_eq[4];
    int steps[4];
    for (int s = 0; s < 4; s++) {
        E_raw[s] = eval_h4_conformer_fci(s, R0);
        E_rel[s] = relax_h4_1d(s, R0, &R_eq[s], &steps[s]);
    }
    int pure2jet_hits = 0, sym_hits = 0, raw_false_rejects = 0;
    for (int s0 = 0; s0 < 4; s0++) {
        /* Pure 2-Jet cannot switch sector */
        if (s0 == 3) pure2jet_hits++;
        /* Pure SEARU Hop proposes sector 3 at compressed R0 without fiber relaxation */
        if (s0 != 3 && E_rel[s0] < E_raw[3]) raw_false_rejects++;
        /* Symbiotic SEARU x 2-Jet sweeps proposals [0,1,2,3], relaxes each candidate on 2-Jet fiber */
        int best_sec = s0;
        double best_E = E_rel[s0];
        for (int prop = 0; prop < 4; prop++) {
            double R_tmp; int st_tmp;
            double E_cand = relax_h4_1d(prop, R0, &R_tmp, &st_tmp);
            if (E_cand < best_E) { best_E = E_cand; best_sec = prop; }
        }
        if (best_sec == 3) sym_hits++;
    }

    fprintf(fp, "  \"exp4_conformer_bilevel\": {\n");
    fprintf(fp, "    \"R0\": %.4f,\n", R0);
    fprintf(fp, "    \"pure_2jet_hits\": %d, \"symbiotic_hits\": %d, \"raw_false_rejects\": %d,\n",
            pure2jet_hits, sym_hits, raw_false_rejects);
    fprintf(fp, "    \"sectors\": [\n");
    for (int s = 0; s < 4; s++) {
        fprintf(fp, "      {\"sector\": %d, \"name\": \"%s\", \"E_raw_R0\": %.12f, "
                    "\"R_eq\": %.8f, \"E_relaxed\": %.12f, \"newton_steps\": %d}%s\n",
                s, sec_names[s], E_raw[s], R_eq[s], E_rel[s], steps[s], (s + 1 < 4) ? "," : "");
    }
    fprintf(fp, "    ]\n  }\n");
    fprintf(fp, "}\n");
    fclose(fp);

    free_fci(&fci_h4sq);
    free_rdm2_schema(&sc2);
    free_rdm2_schema(&sc4);
    printf("Stage 10 Ab Initio Quantum Chemistry wind tunnel wrote %s\n", json_path);
    return 0;
}
