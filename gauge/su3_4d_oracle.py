#!/usr/bin/env python3
"""Independent Stage 6 audit & report generator for 4D Pure SU(3) Color Yang-Mills."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import numpy as np


def weyl_haar_su3_plaq_exact(beta: float, grid_n: int = 256) -> float:
    """Exact 2D Weyl Haar integral for single-link SU(3) with S = -(beta/3) Re Tr U."""
    t = np.linspace(-math.pi, math.pi, grid_n, endpoint=False)
    t1, t2 = np.meshgrid(t, t, indexing="ij")
    t3 = -(t1 + t2)
    re_tr = np.cos(t1) + np.cos(t2) + np.cos(t3)
    vdm2 = (
        (2.0 - 2.0 * np.cos(t1 - t2))
        * (2.0 - 2.0 * np.cos(t1 - t3))
        * (2.0 - 2.0 * np.cos(t2 - t3))
    )
    w = vdm2 * np.exp((beta / 3.0) * (re_tr - np.max(re_tr)))
    return float(np.sum(w * (re_tr / 3.0)) / np.sum(w))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=True)
    ap.add_argument("--report", required=True)
    args = ap.parse_args()

    data = json.loads(Path(args.bench).read_text())
    st = data["selftests"]
    exp1 = sorted(data["exp1_z3"], key=lambda r: float(r["beta"]))
    hyst = data["exp1_hyst_569"]
    exp2 = sorted(data["exp2_casimir"], key=lambda r: float(r["beta"]))
    exp3 = data["exp3_topology"]

    failures: list[str] = []

    print("\n=== [Stage 6 Audit] 4D Pure SU(3) Color Yang-Mills (Theorems C1–C5) ===")

    # 0. Self-tests & exact 2D Weyl Haar single-link oracle check
    if st["max_unit_err"] >= 1e-13:
        failures.append(f"SelfTest unitarity err {st['max_unit_err']:.3e} >= 1e-13")
    if st["max_det_err"] >= 1e-13:
        failures.append(f"SelfTest det err {st['max_det_err']:.3e} >= 1e-13")
    if st["max_gauge_dP"] >= 1e-13:
        failures.append(f"SelfTest gauge dP {st['max_gauge_dP']:.3e} >= 1e-13")
    if st["z3_dS"] >= 1e-10 or st["z3_rot_err"] >= 1e-13 or st["z3_adj_err"] >= 1e-13:
        failures.append("SelfTest Z3 twist tolerance exceeded")
    if st["or_dS"] >= 1e-10 or st["or_inv_err"] >= 1e-13:
        failures.append("SelfTest CM OR exactness/involution tolerance exceeded")
    if st["swap_dS_err"] >= 1e-10:
        failures.append(f"SelfTest @swap dS err {st['swap_dS_err']:.3e} >= 1e-10")
    if st["ape3d_gauge_err"] >= 1e-11 or st["cool_gauge_err"] >= 1e-10:
        failures.append("SelfTest 3D APE / 4D cooling gauge covariance tolerance exceeded")
    if int(st["bpst_pos_round"]) != 1 or int(st["bpst_neg_round"]) != -1:
        failures.append(
            f"SelfTest embedded BPST Q_cool failed: +{st['bpst_pos_round']}, {st['bpst_neg_round']}"
        )

    w3_exact = weyl_haar_su3_plaq_exact(3.0)
    w6_exact = weyl_haar_su3_plaq_exact(6.0)
    z3_sig = abs(st["cm_hb_beta3_mean"] - w3_exact) / max(st["cm_hb_beta3_err"], 1e-9)
    z6_sig = abs(st["cm_hb_beta6_mean"] - w6_exact) / max(st["cm_hb_beta6_err"], 1e-9)
    if z3_sig > 3.5:
        failures.append(
            f"SelfTest CM single-link beta=3.0 mismatch: {st['cm_hb_beta3_mean']:.6f} vs Weyl {w3_exact:.6f} ({z3_sig:.2f}s)"
        )
    if z6_sig > 3.5:
        failures.append(
            f"SelfTest CM single-link beta=6.0 mismatch: {st['cm_hb_beta6_mean']:.6f} vs Weyl {w6_exact:.6f} ({z6_sig:.2f}s)"
        )
    print(
        f"  [SelfTests] Weyl Haar single-link: beta=3.0 CM={st['cm_hb_beta3_mean']:.6f}±{st['cm_hb_beta3_err']:.6f} "
        f"vs exact {w3_exact:.6f} ({z3_sig:.2f}σ) | beta=6.0 CM={st['cm_hb_beta6_mean']:.6f}±{st['cm_hb_beta6_err']:.6f} "
        f"vs exact {w6_exact:.6f} ({z6_sig:.2f}σ)"
    )

    # 1. Experiment 1: 4 x 8^3 Thermal Z_3 First-Order Deconfinement
    print("  --- Exp 1: 4 x 8^3 Thermal Z_3 First-Order Deconfinement ---")
    md_exp1_rows: list[str] = []
    for r in exp1:
        b = float(r["beta"])
        dp_sig = abs(r["p_sym"] - r["p_loc"]) / max(math.hypot(r["p_sym_err"], r["p_loc_err"]), 1e-6)
        dal_sig = abs(r["abs_l_sym"] - r["abs_l_loc"]) / max(
            math.hypot(r["abs_l_sym_err"], r["abs_l_loc_err"]), 1e-5
        )
        if b != 5.69 and dp_sig > 4.5:
            failures.append(f"Exp1 beta={b}: <P> bias between loc and sym ({dp_sig:.2f}σ)")
        if b != 5.69 and dal_sig > 4.5:
            failures.append(f"Exp1 beta={b}: <|L|> bias between loc and sym ({dal_sig:.2f}σ)")

        if b <= 5.56 and r["abs_l_sym"] >= 0.12:
            failures.append(f"Exp1 beta={b}: expected confined <|L|><0.12, got {r['abs_l_sym']:.4f}")
        if b >= 6.14:
            if r["abs_l_sym"] <= 0.25:
                failures.append(f"Exp1 beta={b}: expected deconfined <|L|>>0.25, got {r['abs_l_sym']:.4f}")
            if r["cos3_sym"] <= 0.45:
                failures.append(f"Exp1 beta={b}: expected Z_3 locking <cos(3 arg L)>>0.45, got {r['cos3_sym']:.4f}")
            if r["hops_loc"] >= 8:
                failures.append(f"Exp1 beta={b}: local_cm_hb_or should trap in one Z_3 sector, hops={r['hops_loc']}")
            if r["hops_sym"] <= 300:
                failures.append(f"Exp1 beta={b}: symbiotic_z3_center hops={r['hops_sym']} <= 300")
            if r["tau_rel_sym"] >= 2.0:
                failures.append(f"Exp1 beta={b}: symbiotic_z3_center tau_ReL={r['tau_rel_sym']:.2f} >= 2.0")
            if r["mag_mean_l_sym"] >= 0.06:
                failures.append(
                    f"Exp1 beta={b}: symbiotic_z3_center |<L>|={r['mag_mean_l_sym']:.4f} failed to restore center"
                )
            for si, occ in enumerate(r["sec_sym"]):
                if abs(occ - 1.0 / 3.0) > 0.08:
                    failures.append(f"Exp1 beta={b}: sector {si} occupancy {occ:.3f} far from 1/3")

        print(
            f"    beta={b:.2f} | <P>={r['p_sym']:.5f} <|L|>={r['abs_l_sym']:.4f} cos3={r['cos3_sym']:+.3f} "
            f"<L_8>={r['adj_l_sym']:+.4f} | loc: |<L>|={r['mag_mean_l_loc']:.4f} hops={r['hops_loc']} | "
            f"sym: |<L>|={r['mag_mean_l_sym']:.4f} tau={r['tau_rel_sym']:.2f} hops={r['hops_sym']}"
        )
        md_exp1_rows.append(
            f"| `{b:.2f}` | `{r['p_sym']:.5f} ± {r['p_sym_err']:.5f}` (`loc={r['p_loc']:.5f}`) | "
            f"`{r['abs_l_sym']:.4f} ± {r['abs_l_sym_err']:.4f}` (`cos3={r['cos3_sym']:+.3f}`, `<L_8>={r['adj_l_sym']:+.4f}`) | "
            f"`{r['mag_mean_l_loc']:.4f}` (`hops={r['hops_loc']}`, `sec={r['sec_loc'][0]:.2f}/{r['sec_loc'][1]:.2f}/{r['sec_loc'][2]:.2f}`) | "
            f"**`{r['mag_mean_l_sym']:.4f}`** (`τ={r['tau_rel_sym']:.2f}`, **`{r['hops_sym']}`** hops, `sec={r['sec_sym'][0]:.2f}/{r['sec_sym'][1]:.2f}/{r['sec_sym'][2]:.2f}`) |"
        )

    # First-order hysteresis check at beta = 5.69
    if hyst["cold_abs_l"] <= hyst["hot_abs_l"] + 0.02:
        failures.append(
            f"Exp1 hysteresis@5.69 weak: cold <|L|>={hyst['cold_abs_l']:.4f} vs hot {hyst['hot_abs_l']:.4f}"
        )
    if hyst["coex_frac_lo"] < 0.05 or hyst["coex_frac_hi"] < 0.05:
        failures.append(
            f"Exp1 coexistence@5.69 missing two-state overlap: lo={hyst['coex_frac_lo']:.3f}, hi={hyst['coex_frac_hi']:.3f}"
        )
    print(
        f"    hyst@5.69 | cold:<P>={hyst['cold_p']:.5f} <|L|>={hyst['cold_abs_l']:.4f} vs "
        f"hot:<P>={hyst['hot_p']:.5f} <|L|>={hyst['hot_abs_l']:.4f} | "
        f"coex(<0.09:{hyst['coex_frac_lo']:.3f}, >0.14:{hyst['coex_frac_hi']:.3f})"
    )

    # 2. Experiment 2: 8^4 Scaling & Casimir Scaling C2(8)/C2(3) = 9/4 = 2.25
    print("  --- Exp 2: 8^4 Fundamental vs Adjoint Casimir Scaling ---")
    md_exp2_rows: list[str] = []
    prev_p = 0.0
    prev_chi3 = 1e9
    for r in exp2:
        b = float(r["beta"])
        if r["p_mean"] <= prev_p:
            failures.append(f"Exp2 beta={b}: <P> not strictly increasing")
        prev_p = r["p_mean"]
        if r["chi3_22"] <= 0.0 or r["chi3_22"] >= prev_chi3:
            failures.append(f"Exp2 beta={b}: chi3(2,2)={r['chi3_22']:.4f} failed monotonic drop")
        prev_chi3 = r["chi3_22"]
        if r["chi3_33"] <= 0.0 or r["chi8_22"] <= 0.0:
            failures.append(f"Exp2 beta={b}: positive string tension violated")
        if r["abs_l"] >= 0.05:
            failures.append(f"Exp2 beta={b}: 8^4 escaped confined phase (<|L|>={r['abs_l']:.4f})")
        if not (1.70 <= r["casimir_ratio"] <= 2.60):
            failures.append(f"Exp2 beta={b}: Casimir ratio {r['casimir_ratio']:.3f} outside [1.70, 2.60]")
        print(
            f"    beta={b:.2f} | <P>={r['p_mean']:.5f} <|L|>={r['abs_l']:.4f} | "
            f"chi3(2,2)={r['chi3_22']:.4f}±{r['chi3_22_err']:.4f} chi3(3,3)={r['chi3_33']:.4f}±{r['chi3_33_err']:.4f} | "
            f"chi8(2,2)={r['chi8_22']:.4f}±{r['chi8_22_err']:.4f} | ratio={r['casimir_ratio']:.3f}±{r['casimir_ratio_err']:.3f}"
        )
        md_exp2_rows.append(
            f"| `{b:.2f}` | `{r['p_mean']:.5f} ± {r['p_err']:.5f}` (`{r['abs_l']:.4f}`) | "
            f"`{r['chi3_22']:.4f} ± {r['chi3_22_err']:.4f}` / `{r['chi3_33']:.4f} ± {r['chi3_33_err']:.4f}` | "
            f"`{r['chi8_22']:.4f} ± {r['chi8_22_err']:.4f}` | **`{r['casimir_ratio']:.3f} ± {r['casimir_ratio_err']:.3f}`** (vs `9/4=2.25`) |"
        )

    # 3. Experiment 3: 8^4 beta=5.90 Topology & Defect-Wall PTBC
    print("  --- Exp 3: 8^4 beta=5.90 Topology & Defect-Wall @swap PTBC ---")
    dp_exp3 = abs(exp3["p_pt"] - exp3["p_loc"]) / max(math.hypot(exp3["p_pt_err"], exp3["p_loc_err"]), 1e-6)
    dchi_exp3 = abs(exp3["chi_t_pt"] - exp3["chi_t_loc"]) / max(
        math.hypot(exp3["chi_t_pt_err"], exp3["chi_t_loc_err"]), 1e-8
    )
    if dp_exp3 > 4.0:
        failures.append(f"Exp3 <P> PTBC vs local bias {dp_exp3:.2f}σ (> 4.0)")
    if dchi_exp3 > 4.0:
        failures.append(f"Exp3 chi_t PTBC vs local bias {dchi_exp3:.2f}σ (> 4.0)")
    if exp3["chi_t_loc"] <= 0.0 or exp3["chi_t_pt"] <= 0.0:
        failures.append("Exp3 chi_t must be strictly positive")
    if exp3["sw_min"] < 0.30:
        failures.append(f"Exp3 PTBC min swap acceptance {exp3['sw_min']:.3f} < 0.30")
    if len(exp3["c_ladder"]) != 8 or exp3["c_ladder"][0] != 1.0 or exp3["c_ladder"][-1] != 0.0:
        failures.append("Exp3 c_ladder invalid endpoints")
    for k in range(len(exp3["c_ladder"]) - 1):
        if exp3["c_ladder"][k] <= exp3["c_ladder"][k + 1]:
            failures.append(f"Exp3 c_ladder not strictly decreasing at {k}")
    print(
        f"    local <P>={exp3['p_loc']:.5f}±{exp3['p_loc_err']:.5f} chi_t={exp3['chi_t_loc']:.3e}±{exp3['chi_t_loc_err']:.3e} "
        f"<Q^2>={exp3['q2_mean_loc']:.3f} tau_w={exp3['tau_q_work_loc']:.1f} tun={exp3['tun_loc']}\n"
        f"    PTBC  <P>={exp3['p_pt']:.5f}±{exp3['p_pt_err']:.5f} chi_t={exp3['chi_t_pt']:.3e}±{exp3['chi_t_pt_err']:.3e} "
        f"<Q^2>={exp3['q2_mean_pt']:.3f} tau_w={exp3['tau_q_work_pt']:.1f} tun={exp3['tun_pt']} "
        f"sw_min={100*exp3['sw_min']:.1f}% sw_mean={100*exp3['sw_mean']:.1f}%"
    )

    if failures:
        print("\nSTAGE6 ORACLE FAILURES:")
        for f in failures:
            print("  FAIL:", f)
        return 1

    print("\nALL STAGE 6 4D SU(3) VERIFICATION CHECKS PASSED.")

    ladder_str = ", ".join(f"`{v:.3f}`" for v in exp3["c_ladder"])
    report_lines = [
        "# Stage 6 Report: 4D Pure $SU(3)$ Color Yang–Mills Wind Tunnel (`Theorems C1–C5`)",
        "",
        "Generated from [`build/su3_4d_bench.json`](build/su3_4d_bench.json) via [`su3_4d_windtunnel.c`](su3_4d_windtunnel.c) and audited against [`su3_4d_oracle.py`](su3_4d_oracle.py).",
        "",
        "## 1. Startup Mathematical & Algebraic Self-Tests (`Flow.ColorGaugeBundle`)",
        "",
        "| Check | Theorem Link | Measured Residual / Comparison |",
        "| :--- | :--- | :--- |",
        f"| Unitarity & Determinant | `Theorems C2A–C2E` | `max‖U†U - I‖ = {st['max_unit_err']:.2e} < 1e-13`, `max|det U - 1| = {st['max_det_err']:.2e} < 1e-13` |",
        f"| 4D Local $SU(3)$ Gauge Invariance | `Theorems C1A–C1C` | `max|ΔP| = {st['max_gauge_dP']:.2e} < 1e-13` |",
        f"| $\\mathbb{{Z}}_3$ Slice Twist (`t0=1, ω`) | `Theorems C3A–C4A` | `|ΔS| = {st['z3_dS']:.2e}`, `|L' - ωL| = {st['z3_rot_err']:.2e}`, `|ΔL_8| = {st['z3_adj_err']:.2e}` |",
        f"| Cabibbo–Marinari OR Involution | `Theorems C2C–C2D` | `|ΔS_sweep| = {st['or_dS']:.2e}`, `‖OR^2(U)-U‖ = {st['or_inv_err']:.2e}` |",
        f"| Single-Link CM Heatbath vs Weyl Haar (`200k` each) | `Theorem C2C` | `β=3.0`: `{st['cm_hb_beta3_mean']:.6f}±{st['cm_hb_beta3_err']:.6f}` vs exact `{w3_exact:.6f}` (`{z3_sig:.2f}σ`); `β=6.0`: `{st['cm_hb_beta6_mean']:.6f}±{st['cm_hb_beta6_err']:.6f}` vs exact `{w6_exact:.6f}` (`{z6_sig:.2f}σ`) |",
        f"| 36-Plaquette 4D Defect-Wall `@swap` | `Theorem C5A` | `|ΔS_wall - ΔS_full| = {st['swap_dS_err']:.2e} < 1e-10` |",
        f"| 3D APE & 4D Cooling + Embedded BPST $Q$ | `Theorems C1C–C2E` | `ΔP_3D = {st['ape3d_gauge_err']:.2e}`, `ΔQ_cool = {st['cool_gauge_err']:.2e}`, BPST `Q = {st['bpst_pos_round']:+d} [{st['bpst_pos_q']:+.3f}] / {st['bpst_neg_round']:+d} [{st['bpst_neg_q']:+.3f}]` |",
        "",
        "## 2. Experiment 1 — $4 \\times 8^3$ $\\mathbb{Z}_3$ First-Order Deconfinement & Symbiotic Unfreezing",
        "",
        "| $\\beta$ | `<P>` (`Sym` vs `Loc`) | `<\|L\|>`, `<cos(3 arg L)>` & `<L_8>` | `local_cm_hb_or` (`|<L>|`, `hops`, `sec`) | `symbiotic_z3_center` (`|<L>|`, `τ_ReL`, `hops`, `sec`) |",
        "| ---: | :--- | :--- | :--- | :--- |",
        *md_exp1_rows,
        "",
        f"**First-Order Metastability / Hysteresis at $\\beta = 5.69$:** Cold start gives `<P> = {hyst['cold_p']:.5f} ± {hyst['cold_p_err']:.5f}`, `<|L|> = {hyst['cold_abs_l']:.4f} ± {hyst['cold_abs_l_err']:.4f}` vs. hot start `<P> = {hyst['hot_p']:.5f} ± {hyst['hot_p_err']:.5f}`, `<|L|> = {hyst['hot_abs_l']:.4f} ± {hyst['hot_abs_l_err']:.4f}`; in the `1,200`-sweep coexistence run, both confined (`|L|<0.09`: `{hyst['coex_frac_lo']:.3f}`) and deconfined (`|L|>0.14`: `{hyst['coex_frac_hi']:.3f}`) states are visited repeatedly.",
        "",
        "## 3. Experiment 2 — $8^4$ Confined Scaling & Fundamental vs. Adjoint Casimir Ratio (`Theorem C4C`)",
        "",
        "| $\\beta$ | `<P>` (`<\|L\|>`) | Fundamental $\\chi_3(2,2)$ / $\\chi_3(3,3)$ | Adjoint $\\chi_8(2,2)$ | Correlated Casimir Ratio $\\chi_8(2,2)/\\chi_3(2,2)$ |",
        "| ---: | :--- | :---: | :---: | :---: |",
        *md_exp2_rows,
        "",
        "## 4. Experiment 3 — $8^4$ Topological Susceptibility & 4D Defect-Wall PTBC (`β = 5.90`, `eq_work = 800` sweeps)",
        "",
        f"- **Local `local_cm_hb_or` (`80` meas × `10` sweeps):** `<P> = {exp3['p_loc']:.5f} ± {exp3['p_loc_err']:.5f}`, `<|L|> = {exp3['abs_l_loc']:.4f}`, `<Q^2> = {exp3['q2_mean_loc']:.3f}`, `χ_t = {exp3['chi_t_loc']:.3e} ± {exp3['chi_t_loc_err']:.3e}`, `τ_int^work(Q) = {exp3['tau_q_work_loc']:.1f}`, `tunnelings = {exp3['tun_loc']}`.",
        f"- **Symbiotic `symbiotic_ptbc_4d` (`R = 8`, `50` meas × `2` sweeps, `L_d = 2` [`36` defect plaq], Fisher `2-Jet` ladder [{ladder_str}]):** `<P> = {exp3['p_pt']:.5f} ± {exp3['p_pt_err']:.5f}`, `<|L|> = {exp3['abs_l_pt']:.4f}`, `<Q^2> = {exp3['q2_mean_pt']:.3f}`, `χ_t = {exp3['chi_t_pt']:.3e} ± {exp3['chi_t_pt_err']:.3e}`, `τ_int^work(Q) = {exp3['tau_q_work_pt']:.1f}`, `tunnelings = {exp3['tun_pt']}`, `sw_min = {100.0 * exp3['sw_min']:.1f}%`, `sw_mean = {100.0 * exp3['sw_mean']:.1f}%`.",
        f"- **Unbiased Agreement:** `Δ<P> = {dp_exp3:.2f}σ`, `Δχ_t = {dchi_exp3:.2f}σ`.",
        "",
    ]
    Path(args.report).write_text("\n".join(report_lines))
    print(f"Wrote Stage 6 report: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
