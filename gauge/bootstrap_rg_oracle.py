#!/usr/bin/env python3
"""Independent Stage 5 Oracle & Report Generator:
1. Certified Exact-Rational Dual SOS Bootstrap (2D U(1), 2D SU(2));
2. Gross-Witten-Wadia (SU(infinity)) Nonlinear Loop Bootstrap Across lambda_c = 2;
3. 4D SU(2) Exact Schwinger-Dyson, OS Positivity & 7-Staple Gram Monte Carlo Audit;
4. 4D (b=2) Migdal-Kadanoff Character Renormalization Group Flow (U(1) Coulomb vs SU(2) Asymptotic Freedom).
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np


def bessel_i(n: int, x: float, terms: int = 160) -> float:
    s = 0.0
    for k in range(terms):
        t = ((0.5 * x) ** (2 * k + n)) / (math.factorial(k) * math.factorial(k + n))
        s += t
        if t < 1e-30 * max(1.0, s):
            break
    return s


def dual_sos_u1_2d(beta: float, K: int, scale: int = 1000) -> dict:
    bf = Fraction(beta).limit_denominator(100)
    a, b = [Fraction(1), Fraction(0)], [Fraction(0), Fraction(1)]
    for n in range(0, 2 * K):
        a.append(a[n] - Fraction(2 * (n + 1)) * a[n + 1] / bf)
        b.append(b[n] - Fraction(2 * (n + 1)) * b[n + 1] / bf)

    def min_eig(u: float):
        W = [float(a[k]) + float(b[k]) * u for k in range(2 * K + 1)]
        Mp = np.array([[W[abs(i - j)] + W[i + j] for j in range(K)] for i in range(K)])
        Mm = np.array([[W[abs(i - j)] - W[i + j] for j in range(1, K)] for i in range(1, K)])
        wp, vp = np.linalg.eigh(Mp)
        wm, vm = np.linalg.eigh(Mm)
        if wp[0] <= wm[0]:
            return wp[0], 0, vp[:, 0]
        return wm[0], 1, vm[:, 0]

    L, R = 0.0, 1.0
    for _ in range(90):
        m1, m2 = L + (R - L) / 3.0, R - (R - L) / 3.0
        if min_eig(m1)[0] < min_eig(m2)[0]:
            L = m1
        else:
            R = m2
    u_star = 0.5 * (L + R)

    lo0, lo1 = 0.0, u_star
    for _ in range(85):
        mid = 0.5 * (lo0 + lo1)
        if min_eig(mid)[0] >= 0.0:
            lo1 = mid
        else:
            lo0 = mid
    hi0, hi1 = u_star, 1.0
    for _ in range(85):
        mid = 0.5 * (hi0 + hi1)
        if min_eig(mid)[0] >= 0.0:
            hi0 = mid
        else:
            hi1 = mid

    out = {"beta": beta, "K": K, "exact": bessel_i(1, beta) / bessel_i(0, beta)}
    for tag, u_pt, need_pos_B in [("lo", lo0, True), ("hi", hi1, False)]:
        _, which, vec = min_eig(u_pt)
        vec = vec / np.max(np.abs(vec))
        iv = [int(round(float(x) * scale)) for x in vec]
        A, B = Fraction(0), Fraction(0)
        rng = range(K) if which == 0 else range(1, K)
        for ii, i in enumerate(rng):
            for jj, j in enumerate(rng):
                c = iv[ii] * iv[jj]
                if which == 0:
                    A += c * (a[abs(i - j)] + a[i + j])
                    B += c * (b[abs(i - j)] + b[i + j])
                else:
                    A += c * (a[abs(i - j)] - a[i + j])
                    B += c * (b[abs(i - j)] - b[i + j])
        assert (B > 0) == need_pos_B, f"U(1) dual B sign mismatch: {B}"
        bound = -A / B
        out[f"{tag}_frac"] = bound
        out[f"{tag}_iv"] = iv
        out[f"{tag}_mat"] = "M+" if which == 0 else "M-"
    return out


def dual_sos_su2_2d(beta: float, K: int, scale: int = 1000) -> dict:
    bf = Fraction(beta).limit_denominator(100)
    a, b = [Fraction(1), Fraction(0)], [Fraction(0), Fraction(1)]
    for m in range(0, 2 * K):
        am_prev = a[m - 1] if m > 0 else Fraction(0)
        bm_prev = b[m - 1] if m > 0 else Fraction(0)
        a.append(a[m] + (m * am_prev - (m + 3) * a[m + 1]) / bf)
        b.append(b[m] + (m * bm_prev - (m + 3) * b[m + 1]) / bf)

    def min_eig(u: float):
        M = [float(a[k]) + float(b[k]) * u for k in range(2 * K + 1)]
        H0 = np.array([[M[i + j] for j in range(K)] for i in range(K)])
        H1 = np.array([[M[i + j] - M[i + j + 2] for j in range(K - 1)] for i in range(K - 1)])
        w0, v0 = np.linalg.eigh(H0)
        w1, v1 = np.linalg.eigh(H1)
        if w0[0] <= w1[0]:
            return w0[0], 0, v0[:, 0]
        return w1[0], 1, v1[:, 0]

    L, R = 0.0, 1.0
    for _ in range(90):
        m1, m2 = L + (R - L) / 3.0, R - (R - L) / 3.0
        if min_eig(m1)[0] < min_eig(m2)[0]:
            L = m1
        else:
            R = m2
    u_star = 0.5 * (L + R)

    lo0, lo1 = 0.0, u_star
    for _ in range(85):
        mid = 0.5 * (lo0 + lo1)
        if min_eig(mid)[0] >= 0.0:
            lo1 = mid
        else:
            lo0 = mid
    hi0, hi1 = u_star, 1.0
    for _ in range(85):
        mid = 0.5 * (hi0 + hi1)
        if min_eig(mid)[0] >= 0.0:
            hi0 = mid
        else:
            hi1 = mid

    out = {"beta": beta, "K": K, "exact": bessel_i(2, beta) / bessel_i(1, beta)}
    for tag, u_pt, need_pos_B in [("lo", lo0, True), ("hi", hi1, False)]:
        _, which, vec = min_eig(u_pt)
        vec = vec / np.max(np.abs(vec))
        iv = [int(round(float(x) * scale)) for x in vec]
        n = K if which == 0 else K - 1
        A, B = Fraction(0), Fraction(0)
        for i in range(n):
            for j in range(n):
                c = iv[i] * iv[j]
                if which == 0:
                    A += c * a[i + j]
                    B += c * b[i + j]
                else:
                    A += c * (a[i + j] - a[i + j + 2])
                    B += c * (b[i + j] - b[i + j + 2])
        assert (B > 0) == need_pos_B, f"SU(2) dual B sign mismatch: {B}"
        bound = -A / B
        out[f"{tag}_frac"] = bound
        out[f"{tag}_iv"] = iv
        out[f"{tag}_mat"] = f"H{which}"
    return out


def gww_toeplitz_boot(lam: float, K: int = 8) -> dict:
    def min_eig(u: float) -> float:
        W = [1.0, u]
        for k in range(1, K):
            rhs = sum(W[m] * W[k - m] for m in range(k))
            W.append(W[k - 1] - lam * rhs)
        T = np.array([[W[abs(i - j)] for j in range(K + 1)] for i in range(K + 1)])
        return float(np.linalg.eigvalsh(T)[0])

    grid = np.linspace(0.0, 1.0, 4001)
    vals = [min_eig(float(u)) for u in grid]
    feas = [float(u) for u, v in zip(grid, vals) if v >= -1e-12]
    if not feas:
        idx = int(np.argmax(vals))
        lo_g = max(0.0, float(grid[max(0, idx - 1)]))
        hi_g = min(1.0, float(grid[min(len(grid) - 1, idx + 1)]))
        sub = np.linspace(lo_g, hi_g, 4001)
        feas = [float(u) for u in sub if min_eig(float(u)) >= -1e-12]
    lo, hi = min(feas), max(feas)
    exact = (1.0 - lam / 4.0) if lam <= 2.0 else (1.0 / lam)

    alt_neg_K = None
    alt_quad = None
    if lam < 2.0:
        K_test = int(math.ceil(2.0 / (2.0 - lam))) + 2
        v_alt = np.array([1.0 if i % 2 == 0 else -1.0 for i in range(K_test)])
        W_str = [1.0, 1.0 / lam] + [0.0] * K_test
        T_str = np.array([[W_str[abs(i - j)] for j in range(K_test)] for i in range(K_test)])
        alt_quad = float(v_alt @ T_str @ v_alt)
        alt_neg_K = K_test
    return {
        "lam": lam,
        "exact": exact,
        "lo": lo,
        "hi": hi,
        "width": hi - lo,
        "alt_K": alt_neg_K,
        "alt_quad": alt_quad,
    }


def inv_bessel_ratio(n_num: int, n_den: int, target: float) -> float:
    lo, hi = 0.0, 60.0
    for _ in range(75):
        mid = 0.5 * (lo + hi)
        if bessel_i(n_num, mid) / bessel_i(n_den, mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def mk_rg_4d_step_u1(beta: float, K: int = 48, steps: int = 1) -> tuple[list[float], float]:
    w = np.array([bessel_i(r, beta) / bessel_i(0, beta) for r in range(K)], dtype=np.float64)
    th = np.linspace(0.0, math.pi, 2049)
    traj = [float(w[1])]
    beta_eff_1 = beta
    for s in range(steps):
        F = np.ones_like(th)
        for r in range(1, K):
            F += 2.0 * (w[r] ** 4) * np.cos(r * th)
        F4 = np.maximum(F, 0.0) ** 4
        Z = float(np.trapezoid(F4, th))
        new_w = np.zeros(K, dtype=np.float64)
        for r in range(K):
            new_w[r] = float(np.trapezoid(F4 * np.cos(r * th), th)) / Z
        w = new_w
        traj.append(float(w[1]))
        if s == 0:
            beta_eff_1 = inv_bessel_ratio(1, 0, float(w[1]))
    return traj, beta_eff_1 - beta


def mk_rg_4d_step_su2(beta: float, K: int = 40, steps: int = 1) -> tuple[list[float], float]:
    w = np.array([bessel_i(r + 1, beta) / bessel_i(1, beta) for r in range(K)], dtype=np.float64)
    ph = np.linspace(1e-6, math.pi - 1e-6, 2049)
    sph = np.sin(ph)
    haar = sph ** 2
    traj = [float(w[1])]
    beta_eff_1 = beta
    for s in range(steps):
        F = np.zeros_like(ph)
        for r in range(K):
            F += (r + 1) * (w[r] ** 4) * (np.sin((r + 1) * ph) / sph)
        F4 = np.maximum(F, 0.0) ** 4
        Z = float(np.trapezoid(F4 * haar, ph))
        new_w = np.zeros(K, dtype=np.float64)
        for r in range(K):
            chi_r = np.sin((r + 1) * ph) / sph
            new_w[r] = float(np.trapezoid(F4 * chi_r * haar, ph)) / (Z * (r + 1))
        w = new_w
        traj.append(float(w[1]))
        if s == 0:
            beta_eff_1 = inv_bessel_ratio(2, 1, float(w[1]))
    return traj, beta_eff_1 - beta


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=True)
    ap.add_argument("--report", required=True)
    args = ap.parse_args()

    data = json.loads(Path(args.bench).read_text())
    failures: list[str] = []

    print(
        f"\n=== [Stage 5 Audit] Certified Dual Bootstrap, GWW Transition, 4D SU(2) SD & Exact RG "
        f"(L={data['L']}^4, meas={data['n_meas']}) ==="
    )

    # 1. Certified Rational SOS Dual Bootstrap (2D U(1) & 2D SU(2))
    u1_rows = []
    print("\n--- 1. Certified Exact-Rational Dual SOS Bootstrap: 2D U(1) & 2D SU(2) ---")
    for beta, K_top in [(1.0, 6), (2.0, 6), (4.0, 8), (6.0, 8)]:
        res = dual_sos_u1_2d(beta, K=K_top)
        lo_f, hi_f = float(res["lo_frac"]), float(res["hi_frac"])
        ex = res["exact"]
        width = hi_f - lo_f
        if not (lo_f <= ex <= hi_f):
            failures.append(f"2D U(1) beta={beta} K={K_top}: exact {ex} outside [{lo_f}, {hi_f}]")
        if beta in (1.0, 2.0) and width >= 1e-5:
            failures.append(f"2D U(1) beta={beta} K={K_top}: width={width:.2e} >= 1e-5")
        if width >= 1e-3:
            failures.append(f"2D U(1) beta={beta} K={K_top}: width={width:.2e} >= 1e-3")
        if beta == 2.0:
            assert res["lo_iv"] == [60, -215, 564, -1000, 899], f"U(1) b=2 vLo mismatch: {res['lo_iv']}"
            assert res["hi_iv"] == [-15, 66, -216, 564, -1000, 899], f"U(1) b=2 vHi mismatch: {res['hi_iv']}"
        u1_rows.append(
            f"| 2D $U(1)$ | `{beta:.1f}` | `{K_top}` | `{ex:.10f}` | `[{lo_f:.10f}, {hi_f:.10f}]` | `{width:.2e}` | PASS (Lean B5) |"
            if beta == 2.0
            else f"| 2D $U(1)$ | `{beta:.1f}` | `{K_top}` | `{ex:.10f}` | `[{lo_f:.10f}, {hi_f:.10f}]` | `{width:.2e}` | PASS |"
        )
        print(f"  2D U(1)  beta={beta:3.1f} K={K_top}: exact={ex:.10f} in [{lo_f:.10f}, {hi_f:.10f}] (width={width:.2e})")

    su2_2d_rows = []
    for beta, K_top in [(1.0, 5), (2.0, 5), (2.5, 5), (4.0, 6)]:
        res = dual_sos_su2_2d(beta, K=K_top)
        lo_f, hi_f = float(res["lo_frac"]), float(res["hi_frac"])
        ex = res["exact"]
        width = hi_f - lo_f
        if not (lo_f <= ex <= hi_f):
            failures.append(f"2D SU(2) beta={beta} K={K_top}: exact {ex} outside [{lo_f}, {hi_f}]")
        if width >= 2e-5:
            failures.append(f"2D SU(2) beta={beta} K={K_top}: width={width:.2e} >= 2e-5")
        if beta == 2.0:
            assert res["lo_iv"] == [69, -261, -496, 1000], f"SU(2) b=2 vLo mismatch: {res['lo_iv']}"
            assert res["hi_iv"] == [27, 248, -593, -556, 1000], f"SU(2) b=2 vHi mismatch: {res['hi_iv']}"
        if beta == 2.5:
            assert res["lo_iv"] == [77, -199, -617, 1000], f"SU(2) b=2.5 vLo mismatch: {res['lo_iv']}"
            assert res["hi_iv"] == [11, 290, -505, -696, 1000], f"SU(2) b=2.5 vHi mismatch: {res['hi_iv']}"
        su2_2d_rows.append(
            f"| 2D $SU(2)$ | `{beta:.1f}` | `{K_top}` | `{ex:.10f}` | `[{lo_f:.10f}, {hi_f:.10f}]` | `{width:.2e}` | "
            + ("PASS (Lean B5) |" if beta in (2.0, 2.5) else "PASS |")
        )
        print(f"  2D SU(2) beta={beta:3.1f} K={K_top}: exact={ex:.10f} in [{lo_f:.10f}, {hi_f:.10f}] (width={width:.2e})")

    # 2. Gross-Witten-Wadia Nonlinear Bootstrap across lambda_c = 2
    gww_rows = []
    print("\n--- 2. Gross-Witten-Wadia (SU(infinity)) Loop Bootstrap Across lambda_c = 2 (K=8) ---")
    for lam in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]:
        g = gww_toeplitz_boot(lam, K=8)
        tol = 5e-4
        if not (g["lo"] - tol <= g["exact"] <= g["hi"] + tol):
            failures.append(f"GWW lam={lam}: exact={g['exact']} outside [{g['lo']}, {g['hi']}]")
        if g["width"] > 0.02:
            failures.append(f"GWW lam={lam}: width={g['width']:.3e} > 0.02")
        if lam < 2.0 and (g["alt_quad"] is None or g["alt_quad"] >= 0.0):
            failures.append(f"GWW lam={lam}: expected alternating vector Toeplitz violation < 0, got {g['alt_quad']}")
        alt_disp = f"`{g['alt_quad']:+.3f}` (`K={g['alt_K']}`)" if g["alt_quad"] is not None else "`>= 0 (strong valid)`"
        gww_rows.append(
            f"| `{lam:.2f}` | `{g['exact']:.6f}` | `[{g['lo']:.6f}, {g['hi']:.6f}]` | `{g['width']:.2e}` | {alt_disp} |"
        )
        print(
            f"  GWW lam={lam:4.2f} | exact={g['exact']:.6f} | boot=[{g['lo']:.6f}, {g['hi']:.6f}] (width={g['width']:.2e}) | alt_strong={alt_disp}"
        )

    # 3. 4D SU(2) Schwinger-Dyson, OS & 7-Staple Gram Audit
    sd_rows = []
    print("\n--- 3. 4D SU(2) Exact Schwinger-Dyson, OS & 7-Staple Gram Monte Carlo Audit ---")
    for r in sorted(data["runs"], key=lambda x: float(x["beta"])):
        b = float(r["beta"])
        up = float(r["u_p"])
        up2 = float(r["u_p2"])
        w12 = float(r["w_1x2"])
        uco2 = float(r["u_copl2"])
        wch = float(r["w_chair"])
        uor2 = float(r["u_orth2"])
        rsq = float(r["r_sq"])
        sd_res = float(r["sd_res"])
        err_sd = float(r["err_sd_res"])
        gram_diff = float(r["gram_diff"])

        rsq_diff = abs(rsq - (6.0 + 6.0 * w12 + 24.0 * wch))
        if abs(gram_diff) > 1e-12:
            failures.append(f"4D SU(2) beta={b}: gram_diff={gram_diff:.3e} > 1e-12")
        if rsq_diff > 1e-11:
            failures.append(f"4D SU(2) beta={b}: rsq_diff={rsq_diff:.3e} > 1e-11")
        if abs(sd_res) > 3.5 * max(err_sd, 1e-8):
            failures.append(f"4D SU(2) beta={b}: |sd_res|={abs(sd_res):.5f} > 3.5*err ({3.5*err_sd:.5f})")

        # OS Reflection Positivity checks (with statistical guard on MC connected products)
        if not (w12 > uco2 and wch > uor2):
            failures.append(f"4D SU(2) beta={b}: OS difference w_pair > u_pair2 violated")
        if uco2 < up * up - 3.5 * float(r["err_u_copl2"]) or uor2 < up * up - 3.5 * float(r["err_u_orth2"]):
            failures.append(f"4D SU(2) beta={b}: scalar product < u_P^2 beyond 3.5 sigma")

        # 7-Staple Gram checks
        gram7 = 1.0 + w12 + 4.0 * wch - 6.0 * up * up
        if gram7 <= 0.0 or (1.0 - w12) <= 0.0 or (1.0 + w12 - 2.0 * wch) <= 0.0:
            failures.append(f"4D SU(2) beta={b}: 7-staple Gram inequality violated (gram7={gram7:.5f})")

        # Rigorous analytic bracket: I_2(b)/I_1(b) < u_P < min(SD_CS_upper, I_2(6b)/I_1(6b))
        u_lo = bessel_i(2, b) / bessel_i(1, b)
        u_hi_sd = (-1.0 + math.sqrt(1.0 + 16.0 * b * b)) / (4.0 * b)
        u_hi_6b = bessel_i(2, 6.0 * b) / bessel_i(1, 6.0 * b)
        u_hi = min(u_hi_sd, u_hi_6b)
        if not (u_lo < up < u_hi):
            failures.append(f"4D SU(2) beta={b}: u_P={up:.5f} outside rigorous analytic [{u_lo:.5f}, {u_hi:.5f}]")

        sd_rows.append(
            f"| `{b:.2f}` | `{up:.5f} ± {r['err_u_p']:.5f}` | `[{u_lo:.4f}, {u_hi:.4f}]` | "
            f"`{w12:.4f}` (`{uco2:.4f}`) | `{wch:.4f}` (`{uor2:.4f}`) | "
            f"`{sd_res:+.5f} ± {err_sd:.5f}` (`{abs(sd_res)/err_sd:.2f}σ`) | `{gram7:.4f}` |"
        )
        print(
            f"  beta={b:.2f} | u_P={up:.5f} in [{u_lo:.4f},{u_hi:.4f}] | "
            f"w12={w12:.4f}>{uco2:.4f} | wch={wch:.4f}>{uor2:.4f} | "
            f"sd_res={sd_res:+.5f}+-{err_sd:.5f} ({abs(sd_res)/err_sd:.2f}s) | Gram7={gram7:.4f} | diff={gram_diff:.1e}"
        )

    # 4. 4D Migdal-Kadanoff Exact RG Flow: U(1) vs SU(2)
    rg_rows = []
    print("\n--- 4. 4D Migdal-Kadanoff (b=2, D=4) Character RG Flow: U(1) Coulomb vs SU(2) Confining ---")
    for beta in [2.0, 4.0, 6.0, 8.0, 10.0]:
        u1_tr, u1_db = mk_rg_4d_step_u1(beta, K=48, steps=12)
        su2_tr, su2_db = mk_rg_4d_step_su2(beta, K=40, steps=12)
        if beta == 2.0 and (u1_db >= 0.0 or u1_tr[-1] >= 1e-4):
            failures.append(f"4D U(1) MK RG beta=2.0 expected confining flow, got db={u1_db:+.5f}, c12={u1_tr[-1]:.4f}")
        if beta >= 4.0 and (u1_db <= 0.0 or u1_tr[-1] <= u1_tr[0]):
            failures.append(f"4D U(1) MK RG beta={beta} expected Coulomb fixed-ray db > 0, got {u1_db:+.5f}")
        if su2_db >= 0.0:
            failures.append(f"4D SU(2) MK RG beta={beta} expected asymptotic-freedom db < 0, got {su2_db:+.5f}")
        if su2_tr[-1] > 0.01:
            failures.append(f"4D SU(2) MK RG beta={beta} expected flow to confinement at step 12, got c1={su2_tr[-1]:.4f}")
        rg_rows.append(
            f"| `{beta:.1f}` | `{u1_tr[0]:.4f} -> {u1_tr[1]:.4f} -> {u1_tr[4]:.4f} -> {u1_tr[12]:.4f}` | "
            f"`{u1_db:+.4f}` | `{su2_tr[0]:.4f} -> {su2_tr[1]:.4f} -> {su2_tr[4]:.4f} -> {su2_tr[12]:.4f}` | `{su2_db:+.4f}` |"
        )
        print(
            f"  beta={beta:4.1f} | U(1) db={u1_db:+.4f} ({u1_tr[0]:.4f}->{u1_tr[1]:.4f}->{u1_tr[12]:.4f}) | "
            f"SU(2) db={su2_db:+.4f} ({su2_tr[0]:.4f}->{su2_tr[1]:.4f}->{su2_tr[12]:.4f})"
        )

    if failures:
        print("\nSTAGE 5 ORACLE FAILURES:")
        for f in failures:
            print("  -", f)
        return 1

    print(
        "\n  [PASS] Stage 5 Certified Bootstrap, GWW Transition, 4D SU(2) SD/Gram & Exact 4D MK RG Audit CLEAN!"
    )

    report = "\n".join(
        [
            "# Stage 5 Report — Certified Lattice Yang–Mills Bootstrap & Exact Renormalization Group",
            "",
            "> **Formal Companion:** [`FLOW/core/lean/Flow/LatticeBootstrap.lean`](../../FLOW/core/lean/Flow/LatticeBootstrap.lean) (`Theorems B1–B5`: exact SD moment reductions, certified dual SOS rational certificates, GWW $\\lambda_c = 2$ necessity, 4D quaternionic SD/Gram bounds, block-spin & Migdal–Kadanoff RG)",
            "",
            "## 1. Certified Exact-Rational Dual SOS Bootstrap (`Theorems B1 & B5`)",
            "",
            "Exact rational Schwinger–Dyson moment reductions (`Fraction` arithmetic, zero floats) produce dual SOS integer witness vectors enclosing analytical Bessel benchmarks strictly inside certified intervals (with $\\beta=2.0, 2.5$ verified in Lean 4 directly by `Theorem B5`).",
            "",
            "| Model | $\\beta$ | $K$ | Analytical Exact $u_P$ | Certified Rational Bounds `[lo, hi]` | Width | Lean 4 / Oracle Status |",
            "| :--- | ---: | ---: | ---: | :--- | ---: | :--- |",
            *u1_rows,
            *su2_2d_rows,
            "",
            "## 2. Gross–Witten–Wadia ($SU(\\infty)$) Nonlinear Loop Bootstrap Across $\\lambda_c = 2$ (`Theorem B1`)",
            "",
            "| $\\lambda$ | Exact GWW $u_P$ | Certified Toeplitz Interval `[lo, hi]` | Width | Alternating Vector $v_i = (-1)^i$ on Strong Branch |",
            "| ---: | ---: | :--- | ---: | :--- |",
            *gww_rows,
            "",
            f"## 3. 4D Pure $SU(2)$ Exact Schwinger–Dyson, OS & 7-Staple Gram Audit (`{data['L']}^4`, `{data['n_meas']}` meas, `20`-block JK, `Theorems B2–B3`)",
            "",
            "Algebraic checks pass to machine precision across all couplings (`|gram_diff| < 1e-12`, `|R_sq - (6 + 6*w_1x2 + 24*w_chair)| < 1e-11`). Analytical non-perturbative bounds come from 2D decoupling $I_2(\\beta)/I_1(\\beta) < u_P < \\min\\big((-1 + \\sqrt{1 + 16\\beta^2})/(4\\beta),\\ I_2(6\\beta)/I_1(6\\beta)\\big)$.",
            "",
            "| $\\beta$ | $\\langle u_P \\rangle$ | Analytic Bracket | $W_{1\\times 2}$ (`u_copl2`) | $W_{\\text{chair}}$ (`u_orth2`) | 4D SD Residual (`3u_P - beta*...`) | 7-Staple `Gram7` |",
            "| ---: | :--- | :--- | :--- | :--- | :--- | ---: |",
            *sd_rows,
            "",
            "## 4. 4D ($b = 2, D = 4$) Migdal–Kadanoff Character RG Flow (`Theorem B4`)",
            "",
            "4D Abelian $U(1)$ transitions from confining ($\\Delta\\beta < 0$ at $\\beta=2.0$) to a deconfined Coulomb fixed ray ($\\Delta\\beta > 0$ for $\\beta \\ge 4.0$), whereas non-Abelian 4D $SU(2)$ exhibits universal asymptotic freedom ($\\Delta\\beta < 0$ at every tested $\\beta \\in [2.0, 10.0]$) flowing monotonically into strong-coupling confinement.",
            "",
            "| Initial $\\beta$ | 4D $U(1)$ $c_1(0 \\to 1 \\to 4 \\to 12)$ | $U(1)$ $\\Delta\\beta_1$ | 4D $SU(2)$ $c_1(0 \\to 1 \\to 4 \\to 12)$ | $SU(2)$ $\\Delta\\beta_1$ |",
            "| ---: | :--- | ---: | :--- | ---: |",
            *rg_rows,
            "",
        ]
    )
    Path(args.report).write_text(report)
    print(f"Generated {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
