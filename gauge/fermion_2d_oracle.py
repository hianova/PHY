#!/usr/bin/env python3
"""Independent Stage 7 audit & report generator for 2D/4D Lattice Fermions (Theorems F1–F5)."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import numpy as np


def free_staggered_condensate_2d(l: int, m: float) -> float:
    """Exact 2D LxL antiperiodic-t periodic-x free staggered chiral condensate (1/V) Re Tr D^{-1}."""
    kt = np.array([(2.0 * n + 1.0) * math.pi / l for n in range(l)])
    kx = np.array([(2.0 * n) * math.pi / l for n in range(l)])
    s2 = np.sin(kt)[:, None] ** 2 + np.sin(kx)[None, :] ** 2
    return float(np.mean(m / (m * m + s2)))


def exact_bessel_chain_1d(beta: float, mu: float, l: int, n_quad: int = 256) -> tuple[float, float]:
    """Exact (Z/Z_pq average sign, exact <n>, Z/Z_0) for 1D U(1) closed ring via quadrature Fourier modes."""
    th = np.linspace(-math.pi, math.pi, n_quad, endpoint=False)
    k = np.arange(-28, 29)[:, None]
    ik = np.mean(np.exp(beta * np.cos(th))[None, :] * np.cos(k * th[None, :]), axis=1)
    log_w = np.where(ik > 1e-300, l * np.log(np.maximum(ik, 1e-300)) + k[:, 0] * mu * l, -1e9)
    max_lw = float(np.max(log_w))
    w = np.exp(log_w - max_lw)
    n_exact = float(np.sum(k[:, 0] * w) / np.sum(w))
    log_z = max_lw + math.log(float(np.sum(w)))

    b_pq = beta * math.cosh(mu)
    ik_pq = np.mean(np.exp(b_pq * np.cos(th))[None, :] * np.cos(k * th[None, :]), axis=1)
    log_wpq = np.where(ik_pq > 1e-300, l * np.log(np.maximum(ik_pq, 1e-300)), -1e9)
    max_lpq = float(np.max(log_wpq))
    log_zpq = max_lpq + math.log(float(np.sum(np.exp(log_wpq - max_lpq))))
    exact_sgn = math.exp(log_z - log_zpq)
    return exact_sgn, n_exact


def exact_thimble_quad_1d(a: float, lam: float, mu: float, grid_n: int = 2001) -> tuple[float, float, float]:
    """Exact 2001-point quadrature for S(x) = 0.5 a x^2 + 0.25 lam x^4 - i mu x."""
    x = np.linspace(-10.0, 10.0, grid_n)
    sr = 0.5 * a * x * x + 0.25 * lam * (x ** 4)
    sr_min = float(np.min(sr))
    w = np.exp(-(sr - sr_min))
    c = np.cos(mu * x)
    s = np.sin(mu * x)
    z_re = float(np.sum(w * c))
    z_abs = float(np.sum(w))
    naive_sgn = abs(z_re) / z_abs
    re_z2 = float(np.sum(w * c * (x * x)) / z_re)
    im_z = float(np.sum(w * s * x) / z_re)
    return naive_sgn, re_z2, im_z


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=True)
    ap.add_argument("--report", required=True)
    args = ap.parse_args()

    data = json.loads(Path(args.bench).read_text())
    st = data["selftests"]
    exp1 = data["exp1_index"]
    exp2 = sorted(data["exp2_schwinger"], key=lambda r: float(r["beta"]))
    mass_scan = sorted(data["exp2_mass_scan"], key=lambda r: float(r["mass"]))
    exp3a = sorted(data["exp3a_chain"], key=lambda r: int(r["L"]))
    exp3b = sorted(data["exp3b_thimble"], key=lambda r: float(r["mu"]))

    failures: list[str] = []
    print("\n=== [Stage 7 Audit] 2D/4D Lattice Fermions, Atiyah-Singer Index & Sign Problem (Theorems F1–F5) ===")

    # 0. Self-Tests
    if st["max_jw_err"] >= 1e-14:
        failures.append(f"SelfTest JW CAR err {st['max_jw_err']:.3e} >= 1e-14")
    if st["max_clifford_err"] >= 1e-14:
        failures.append(f"SelfTest Clifford err {st['max_clifford_err']:.3e} >= 1e-14")
    if st["max_g5_herm_err"] >= 1e-13 or st["max_gauge_cov_err"] >= 1e-13:
        failures.append("SelfTest Wilson g5H / gauge cov tolerance exceeded")
    if st["nn_naive_zeros"] != 4 or st["nn_wilson_zeros"] != 1 or st["nn_wilson_doubler_min"] < 1.9:
        failures.append("SelfTest Nielsen-Ninomiya doubling/lifting mismatch")
    if st["stag_pi_violations"] != 0:
        failures.append(f"SelfTest 4D staggered pi-flux violations={st['stag_pi_violations']}")
    if st["ov_gw_err"] >= 1e-11 or st["ov_norm_err"] >= 1e-11 or st["ov_circle_err"] >= 1e-11:
        failures.append("SelfTest Overlap GW/norm/circle tolerance exceeded")
    if st["hmc_force_err"] >= 1e-6:
        failures.append(f"SelfTest HMC force err {st['hmc_force_err']:.3e} >= 1e-6")
    if st["lefschetz_st_im_err"] >= 1e-14:
        failures.append(f"SelfTest Lefschetz saddle Im S {st['lefschetz_st_im_err']:.3e} >= 1e-14")

    # 1. Experiment 1: Atiyah-Singer Index Theorem
    print("  --- Exp 1: 8x8 Discrete Atiyah-Singer Index Theorem (Overlap D_ov) ---")
    md_exp1: list[str] = []
    orient = None
    for r in exp1:
        qg = int(r["q_geo"])
        iz = int(r["index_zero_mode"])
        it = int(r["index_trace"])
        if qg != int(r["q_target"]):
            failures.append(f"Exp1 Q_geo={qg} != q_target={r['q_target']}")
        if qg != 0:
            s = 1 if iz == qg else (-1 if iz == -qg else 0)
            if s == 0:
                failures.append(f"Exp1 Atiyah-Singer mismatch: Q_geo={qg}, n+-n-={iz}")
            elif orient is None:
                orient = s
            elif orient != s:
                failures.append(f"Exp1 inconsistent index orientation at Q={qg}: expected {orient * qg}, got {iz}")
        elif iz != 0:
            failures.append(f"Exp1 Q=0 nonzero index {iz}")
        if r["max_nonzero_chiral"] >= 1e-10:
            failures.append(f"Exp1 non-zero-mode chirality {r['max_nonzero_chiral']:.3e} >= 1e-10")
        if r["max_eig_err"] >= 1e-8:
            failures.append(f"Exp1 eigenvector residual {r['max_eig_err']:.3e} >= 1e-8")
        tag = "noisy" if r["noisy"] else "smooth"
        print(
            f"    {tag:6s} Q_geo={qg:+d} | n_+={r['n_plus']} n_-={r['n_minus']} | "
            f"idx={iz:+d} tr_idx={it:+d} | max_nz_chi={r['max_nonzero_chiral']:.2e}"
        )
        md_exp1.append(
            f"| `{tag}` | `{qg:+d}` | `({r['n_plus']}, {r['n_minus']})` | **`{iz:+d}`** (`tr={it:+d}`) | `{r['max_nonzero_chiral']:.2e}` |"
        )

    # 2. Experiment 2: 6x6 Dynamical Schwinger Model vs Quenched
    print("  --- Exp 2: 6x6 Dynamical N_f=2 Staggered Schwinger HMC vs. Quenched ---")
    md_exp2: list[str] = []
    for r in exp2:
        b = float(r["beta"])
        if r["q_acc"] < 0.65 or r["d_acc"] < 0.65:
            failures.append(f"Exp2 beta={b}: HMC acceptance < 0.65")
        if r["d_p"] <= r["q_p"]:
            failures.append(f"Exp2 beta={b}: <P>_dyn={r['d_p']:.5f} <= <P>_quench={r['q_p']:.5f}")
        if r["d_chi22"] >= r["q_chi22"]:
            failures.append(f"Exp2 beta={b}: chi_dyn={r['d_chi22']:.4f} >= chi_quench={r['q_chi22']:.4f}")
        print(
            f"    beta={b:.1f} | Q: <P>={r['q_p']:.5f}±{r['q_p_err']:.5f} chi22={r['q_chi22']:.4f}±{r['q_chi22_err']:.4f} "
            f"<pbp>={r['q_pbp']:.5f} (acc={100*r['q_acc']:.1f}%) | "
            f"D: <P>={r['d_p']:.5f}±{r['d_p_err']:.5f} chi22={r['d_chi22']:.4f}±{r['d_chi22_err']:.4f} "
            f"<pbp>={r['d_pbp']:.5f} (acc={100*r['d_acc']:.1f}%)"
        )
        md_exp2.append(
            f"| `{b:.1f}` | `{r['q_p']:.5f} ± {r['q_p_err']:.5f}` | **`{r['d_p']:.5f} ± {r['d_p_err']:.5f}`** | "
            f"`{r['q_chi22']:.4f} ± {r['q_chi22_err']:.4f}` -> **`{r['d_chi22']:.4f} ± {r['d_chi22_err']:.4f}`** | "
            f"`{r['d_pbp']:.5f} ± {r['d_pbp_err']:.5f}` (`acc={100*r['d_acc']:.1f}%`) |"
        )

    md_mass: list[str] = []
    for r in mass_scan:
        m = float(r["mass"])
        free_pbp = free_staggered_condensate_2d(6, m)
        if r["pbp_mean"] <= 0.0:
            failures.append(f"Exp2 mass={m}: non-positive condensate")
        md_mass.append(
            f"| `{m:.2f}` | `{r['p_mean']:.5f}` | **`{r['pbp_mean']:.5f} ± {r['pbp_err']:.5f}`** | `{free_pbp:.5f}` | `{r['chi22']:.4f}` |"
        )

    # 3A. 1D Finite-Density Sign Problem vs Worldline Dual
    print("  --- Exp 3A: 1D Finite-Density U(1) Sign Problem vs. Positive Worldline Dual ---")
    md_3a: list[str] = []
    prev_sgn = 2.0
    for r in exp3a:
        L = int(r["L"])
        ex_sgn, ex_n = exact_bessel_chain_1d(1.0, 0.6, L)
        z_sgn = abs(r["pq_sgn"] - ex_sgn) / max(r["pq_sgn_err"], 1e-6)
        z_dual = abs(r["n_dual"] - ex_n) / max(r["n_dual_err"], 1e-6)
        if r["pq_sgn"] >= prev_sgn:
            failures.append(f"Exp3A L={L}: naive sign not strictly decreasing")
        prev_sgn = r["pq_sgn"]
        if L == 64 and r["pq_sgn"] >= 0.05:
            failures.append(f"Exp3A L=64: expected severe naive sign < 0.05, got {r['pq_sgn']:.4f}")
        if abs(r["dual_sgn"] - 1.0) > 1e-12:
            failures.append(f"Exp3A L={L}: dual sign {r['dual_sgn']} != 1.0")
        if z_sgn > 4.0:
            failures.append(f"Exp3A L={L}: naive pq <sgn> {r['pq_sgn']:.5f} vs exact {ex_sgn:.5f} ({z_sgn:.2f}σ)")
        if z_dual > 3.5:
            failures.append(f"Exp3A L={L}: dual <n> {r['n_dual']:.5f} vs exact {ex_n:.5f} ({z_dual:.2f}σ)")
        print(
            f"    L={L:2d} | Naive <sgn>={r['pq_sgn']:.5f} (exact {ex_sgn:.5f}) <n>_rw={r['n_rw']:+.4f}±{r['n_rw_err']:.4f} | "
            f"Dual <sgn>={r['dual_sgn']:.6f} <n>={r['n_dual']:+.5f}±{r['n_dual_err']:.5f} (exact {ex_n:+.5f}, {z_dual:.2f}σ)"
        )
        md_3a.append(
            f"| `{L}` | `{r['pq_sgn']:.5f} ± {r['pq_sgn_err']:.5f}` (`exact={ex_sgn:.5f}`) | `{r['n_rw']:+.4f} ± {r['n_rw_err']:.4f}` | "
            f"**`1.000000`** | **`{r['n_dual']:+.5f} ± {r['n_dual_err']:.5f}`** (`exact={ex_n:+.5f}`, `{z_dual:.2f}σ`) |"
        )

    # 3B. Non-Gaussian Lefschetz Thimble
    print("  --- Exp 3B: Non-Gaussian Lefschetz Thimble ---")
    md_3b: list[str] = []
    prev_ns = 2.0
    for r in exp3b:
        mu = float(r["mu"])
        ex_nsgn, ex_z2, ex_iz = exact_thimble_quad_1d(2.0, 0.05, mu)
        z_z2 = abs(r["th_re_z2"] - ex_z2) / max(r["th_re_z2_err"], 1e-6)
        z_iz = abs(r["th_im_z"] - ex_iz) / max(r["th_im_z_err"], 1e-6)
        if r["naive_sgn"] >= prev_ns:
            failures.append(f"Exp3B mu={mu}: naive sign not strictly decreasing")
        prev_ns = r["naive_sgn"]
        if mu >= 3.9 and r["naive_sgn"] >= 0.05:
            failures.append(f"Exp3B mu={mu}: expected severe naive sign < 0.05, got {r['naive_sgn']:.4f}")
        if r["th_sgn"] <= 0.85:
            failures.append(f"Exp3B mu={mu}: Lefschetz thimble sign {r['th_sgn']:.4f} <= 0.85")
        if z_z2 > 3.5 or z_iz > 3.5:
            failures.append(f"Exp3B mu={mu}: thimble mismatch z2={z_z2:.2f}σ iz={z_iz:.2f}σ")
        print(
            f"    mu={mu:.1f} y0={r['y0']:.4f} | Naive <sgn>={r['naive_sgn']:.5f} (exact {ex_nsgn:.5f}) "
            f"<z^2>={r['naive_re_z2']:+.3f}±{r['naive_re_z2_err']:.3f} | "
            f"Thimble <sgn>={r['th_sgn']:.5f} <z^2>={r['th_re_z2']:+.4f}±{r['th_re_z2_err']:.4f} (ex {ex_z2:+.4f}, {z_z2:.2f}σ) "
            f"<Im z>={r['th_im_z']:+.4f}±{r['th_im_z_err']:.4f} (ex {ex_iz:+.4f}, {z_iz:.2f}σ)"
        )
        md_3b.append(
            f"| `{mu:.1f}` (`y0={r['y0']:.4f}`) | `{r['naive_sgn']:.5f}` (`<z^2>={r['naive_re_z2']:+.3f}±{r['naive_re_z2_err']:.3f}`) | "
            f"**`{r['th_sgn']:.5f}`** | **`{r['th_re_z2']:+.4f} ± {r['th_re_z2_err']:.4f}`** (`ex={ex_z2:+.4f}`, `{z_z2:.2f}σ`) | "
            f"**`{r['th_im_z']:+.4f} ± {r['th_im_z_err']:.4f}`** (`ex={ex_iz:+.4f}`, `{z_iz:.2f}σ`) |"
        )

    if failures:
        print("\nSTAGE7 ORACLE FAILURES:")
        for f in failures:
            print("  FAIL:", f)
        return 1

    print("\nALL STAGE 7 LATTICE FERMION VERIFICATION CHECKS PASSED.")

    orient_str = "+Q_geo" if (orient == 1) else "-Q_geo"
    report_lines = [
        "# Stage 7 Report: 2D/4D Lattice Fermions, Atiyah–Singer Index & Sign Problem (`Theorems F1–F5`)",
        "",
        "Generated from [`build/fermion_2d_bench.json`](build/fermion_2d_bench.json) via [`fermion_2d_windtunnel.c`](fermion_2d_windtunnel.c) and audited against [`fermion_2d_oracle.py`](fermion_2d_oracle.py).",
        "",
        "## 1. Startup Mathematical & Algebraic Self-Tests (`Flow.FermionBundle`)",
        "",
        "| Check | Theorem Link | Measured Residual / Value |",
        "| :--- | :--- | :--- |",
        f"| 64-Dim ($n=6$) Jordan–Wigner CAR, Pauli & $n_k^2$ | `Theorems F1A–F1F` | `max_err = {st['max_jw_err']:.2e} < 1e-14` |",
        f"| 2D & 4D Euclidean Clifford Algebra + Wilson `γ_5H` & `U(1)` Cov | `Theorems F2A–F2C` | `Cliff = {st['max_clifford_err']:.2e}`, `γ_5H = {st['max_g5_herm_err']:.2e}`, `cov = {st['max_gauge_cov_err']:.2e}` |",
        f"| Nielsen–Ninomiya 2D Doubling / Wilson Lifting + 4D Staggered $\\pi$-Flux | `Theorems F3A–F3D` | Naive zeros `{st['nn_naive_zeros']}` -> Wilson zeros `{st['nn_wilson_zeros']}` (`m_min = {st['nn_wilson_doubler_min']:.3f}`), 4D `pi` violations `{st['stag_pi_violations']}` |",
        f"| $6\\times 6$ Overlap GW Relation, Normality & Unit Circle | `Theorems F4A–F4C` | `GW = {st['ov_gw_err']:.2e}`, `norm = {st['ov_norm_err']:.2e}`, `circle = {st['ov_circle_err']:.2e}` |",
        f"| $6\\times 6$ Staggered HMC Exact Fermion Force vs. Central FD | `Theorem F2D` | `max|F_ana - F_FD| = {st['hmc_force_err']:.2e} < 1e-6` |",
        f"| Lefschetz Quadratic Saddle Cancellation | `Theorem F5E` | `max|Im S| = {st['lefschetz_st_im_err']:.2e} < 1e-14` |",
        "",
        f"## 2. Experiment 1 — $8\\times 8$ Discrete Atiyah–Singer Index Theorem (`Theorem F4`, `n_+ - n_- = {orient_str}`)",
        "",
        "| Background | `Q_geo` | Zero Modes `(n_+, n_-)` | `index(D_ov) = n_+ - n_-` | Max Non-Zero Chirality `|⟨v, γ_5 v⟩|` |",
        "| :--- | :---: | :---: | :---: | :---: |",
        *md_exp1,
        "",
        "## 3. Experiment 2 — $6\\times 6$ Dynamical $N_f = 2$ Staggered Schwinger Model vs. Quenched (`m = 0.15`)",
        "",
        "| $\\beta$ | Quenched `<P>` | Dynamical $N_f=2$ `<P>` | Screening $\\chi(2,2)$ (`Quench -> Dyn`) | Dynamical $\\langle \\bar\\psi\\psi \\rangle$ (`acc`) |",
        "| ---: | :---: | :---: | :---: | :---: |",
        *md_exp2,
        "",
        "### Chiral Condensate Mass Scan at $\\beta = 2.0$ ($6\\times 6$ $N_f=2$ vs. Free Field)",
        "",
        "| Mass $m$ | `<P>` | Dynamical $\\langle \\bar\\psi\\psi \\rangle$ | Free Staggered $\\langle \\bar\\psi\\psi \\rangle_{\\text{free}}$ | $\\chi(2,2)$ |",
        "| ---: | :---: | :---: | :---: | :---: |",
        *md_mass,
        "",
        "## 4. Experiment 3A — 1D Finite-Density $U(1)$ Sign Problem vs. Positive Worldline Dual (`β=1.0, μ=0.6`)",
        "",
        "| Chain $L$ | Naive PQ `<sgn>` | Naive Reweighted `<n>` | Dual `<sgn>` | Positive Dual `<n>` (vs. Exact Bessel) |",
        "| ---: | :---: | :---: | :---: | :--- |",
        *md_3a,
        "",
        "## 5. Experiment 3B — Non-Gaussian Complex Action & Lefschetz Thimble (`a=2.0, λ=0.05`)",
        "",
        "| $\\mu$ (`y_0`) | Naive Real-Axis `<sgn>` & `<z^2>` | Thimble `<sgn>` | Thimble `<z^2>` (vs. Exact Quad) | Thimble `<Im z>` (vs. Exact Quad) |",
        "| :--- | :--- | :---: | :--- | :--- |",
        *md_3b,
        "",
    ]
    Path(args.report).write_text("\n".join(report_lines))
    print(f"Wrote Stage 7 report: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
