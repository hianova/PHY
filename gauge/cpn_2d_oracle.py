#!/usr/bin/env python3
"""Independent Stage 2 audit & report generator for 2D CP^{N-1} Nonlinear Sigma Model."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=True)
    ap.add_argument("--report", required=True)
    args = ap.parse_args()

    data = json.loads(Path(args.bench).read_text())
    n_comp = int(data["N"])
    L = int(data["L"])
    vol = int(data["V"])
    l_def = int(data["L_defect"])
    n_rep = int(data["R_replicas"])
    eq_work = int(data["eq_work_sweeps"])

    by_beta: dict[float, dict[str, dict]] = {}
    for row in data["runs"]:
        by_beta.setdefault(float(row["beta"]), {})[str(row["solver"])] = row

    betas = sorted(by_beta.keys())
    if not betas:
        print("STAGE2 ORACLE FAIL: empty benchmark data")
        return 1

    chi_xi2_asymp = 1.0 / (2.0 * math.pi * float(n_comp))
    chi_xi2_lo = 0.4 * chi_xi2_asymp
    chi_xi2_hi = 3.0 * chi_xi2_asymp
    beta_lo = betas[0]
    beta_hi = betas[-1]
    failures: list[str] = []

    print(
        f"\n=== [Stage 2 Audit] 2D CP^{n_comp - 1} "
        f"(N={n_comp}, L={L}, L_d={l_def}, R={n_rep}, eq_work={eq_work}) ==="
    )
    prev_xi = 0.0
    md_obs_rows: list[str] = []
    md_tau_rows: list[str] = []

    for b in betas:
        sol = by_beta[b]
        loc = sol["local_or_metro"]
        bare = sol["bare_bp_inst"]
        sym = sol["symbiotic_ptbc_2jet"]

        de_sig = abs(sym["mean_e"] - loc["mean_e"]) / max(
            math.hypot(sym["err_e"], loc["err_e"]), 1e-6
        )
        dxi_sig = abs(sym["xi_g"] - loc["xi_g"]) / max(
            math.hypot(sym["err_xi_g"], loc["err_xi_g"]), 1e-4
        )
        dchi_sig = abs(sym["chi_t_cool"] - loc["chi_t_cool"]) / max(
            math.hypot(sym["err_chi_t_cool"], loc["err_chi_t_cool"]), 1e-6
        )

        l_over_xi = L / sym["xi_g"] if sym["xi_g"] > 0.0 else 0.0
        loc_tau_disp = "FROZEN" if loc["frozen_q_cool"] else f"{loc['tau_work_q_cool']:.1f}"
        bare_tau_disp = "FROZEN" if bare["frozen_q_cool"] else f"{bare['tau_work_q_cool']:.1f}"
        sym_tau_disp = "FROZEN" if sym["frozen_q_cool"] else f"{sym['tau_work_q_cool']:.1f}"
        if sym["frozen_q_cool"] or sym["tau_work_q_cool"] <= 0.0:
            speedup_disp = "N/A"
        elif loc["frozen_q_cool"]:
            speedup_disp = "inf (FROZEN)"
        else:
            speedup_disp = f"{loc['tau_work_q_cool'] / sym['tau_work_q_cool']:.2f}x"

        print(
            f"  beta={b:.2f} | Sym(E={sym['mean_e']:.5f}+-{sym['err_e']:.5f}, "
            f"xi={sym['xi_g']:.3f}+-{sym['err_xi_g']:.3f}, L/xi={l_over_xi:.2f}, "
            f"chi_c={sym['chi_t_cool']:.6f}+-{sym['err_chi_t_cool']:.6f}, "
            f"chi*xi2={sym['chi_xi2_cool']:.4f}, tauWc={sym_tau_disp}, tc={sym['trans_q_cool']}, "
            f"sw_min={sym['min_swap_acc']*100:.1f}%) | "
            f"Loc(xi={loc['xi_g']:.3f}, chi_c={loc['chi_t_cool']:.6f}, "
            f"tauWc={loc_tau_disp}, tc={loc['trans_q_cool']}) | "
            f"Bare(acc={bare['bare_inst_acc']*100:.4f}%, tauWc={bare_tau_disp}) | "
            f"speedup={speedup_disp}"
        )

        # 1. Zero-bias agreement at lowest beta against independent local_or_metro
        if abs(b - beta_lo) < 1e-6:
            if loc["frozen_q_cool"]:
                failures.append(f"beta={b}: baseline local_or_metro froze at reference anchor beta")
            if de_sig > 4.0:
                failures.append(f"beta={b}: energy deviation {de_sig:.2f}s (> 4.0)")
            if dxi_sig > 4.0:
                failures.append(f"beta={b}: xi_G deviation {dxi_sig:.2f}s (> 4.0)")
            if dchi_sig > 4.0:
                failures.append(f"beta={b}: chi_t_cool deviation {dchi_sig:.2f}s (> 4.0)")

        # 2. Corollary 2B: bare BP instanton acceptance < 1% for beta >= 0.65
        if b >= 0.65 - 1e-9 and bare["bare_inst_acc"] >= 0.01:
            failures.append(f"beta={b}: bare_bp_inst acc={bare['bare_inst_acc']} >= 0.01")

        # 3. Fisher 2-Jet ladder balance: min_swap_acc > 15%
        if sym["min_swap_acc"] <= 0.15:
            failures.append(f"beta={b}: min_swap_acc={sym['min_swap_acc']:.4f} <= 0.15")

        # 4. Strict monotonic asymptotic freedom scaling xi_G(beta), L/xi_G >= 5.0, and O(1/(2*pi*N)) chi*xi^2
        if sym["xi_g"] <= prev_xi:
            failures.append(f"beta={b}: xi_G={sym['xi_g']:.3f} did not increase above {prev_xi:.3f}")
        prev_xi = sym["xi_g"]
        if l_over_xi < 5.0:
            failures.append(f"beta={b}: finite-volume squeeze L/xi_G={l_over_xi:.2f} < 5.0")
        if not (chi_xi2_lo <= sym["chi_xi2_cool"] <= chi_xi2_hi):
            failures.append(
                f"beta={b}: unphysical chi_t*xi_G^2={sym['chi_xi2_cool']:.4f} "
                f"(expected in [{chi_xi2_lo:.4f}, {chi_xi2_hi:.4f}], asymp={chi_xi2_asymp:.4f})"
            )

        # 5. Strict equal-work unfreezing at highest beta
        if abs(b - beta_hi) < 1e-6:
            if sym["frozen_q_cool"] or sym["trans_q_cool"] < 20:
                failures.append(f"beta={b}: symbiotic_ptbc_2jet froze (tc={sym['trans_q_cool']})")
            if not loc["frozen_q_cool"]:
                if loc["tau_work_q_cool"] <= 3.0 * sym["tau_work_q_cool"]:
                    failures.append(
                        f"beta={b}: local_tau_work={loc['tau_work_q_cool']:.1f} "
                        f"<= 3 * sym_tau_work={sym['tau_work_q_cool']:.1f}"
                    )

        loc_tau_md = "**FROZEN**" if loc["frozen_q_cool"] else f"`{loc['tau_work_q_cool']:.1f}`"
        bare_tau_md = "**FROZEN**" if bare["frozen_q_cool"] else f"`{bare['tau_work_q_cool']:.1f}`"
        sym_tau_md = f"**`{sym['tau_work_q_cool']:.1f}`**"
        speedup_md = (
            "**`∞ (FROZEN)`**"
            if loc["frozen_q_cool"]
            else f"**`{loc['tau_work_q_cool'] / sym['tau_work_q_cool']:.2f}×`**"
        )
        ladder_str = ", ".join(f"{v:.2f}" for v in sym["c_ladder"])

        md_obs_rows.append(
            f"| `{b:.2f}` | `{sym['mean_e']:.5f} ± {sym['err_e']:.5f}` (`loc={loc['mean_e']:.5f}`) | "
            f"`{sym['xi_g']:.3f} ± {sym['err_xi_g']:.3f}` | `{l_over_xi:.2f}` | "
            f"`{sym['chi_t_cool']:.6f} ± {sym['err_chi_t_cool']:.6f}` | "
            f"`{sym['chi_xi2_cool']:.4f}` | `{loc['chi_t_cool']:.6f} ± {loc['err_chi_t_cool']:.6f}` |"
        )
        md_tau_rows.append(
            f"| `{b:.2f}` | {loc_tau_md} (`{loc['trans_q_cool']}` trans) | "
            f"{bare_tau_md} (`acc={bare['bare_inst_acc']*100:.4f}%`) | "
            f"{sym_tau_md} (`{sym['trans_q_cool']}` trans, `swap min={sym['min_swap_acc']*100:.1f}%`, "
            f"`mean={sym['mean_swap_acc']*100:.1f}%`) | {speedup_md} | `{ladder_str}` |"
        )

    if failures:
        for msg in failures:
            print(f"STAGE2 ORACLE FAIL: {msg}")
        return 1

    print(
        f"  [PASS] Stage 2 2D CP^{n_comp - 1} Equal-Work Oracle Audit CLEAN "
        f"(zero bias at beta={beta_lo:.2f}, L/xi>=5.0, Fisher ladder >=15% min swap, >=3x unfreezing at beta={beta_hi:.2f}!)"
    )

    md = [
        f"# Stage 2 Report — 2D $CP^{{{n_comp - 1}}}$ ($N = {n_comp}$) Nonlinear $\\sigma$ Model Equal-Work Topological Unfreezing",
        "",
        (
            f"> **Lattice:** `L = {L}` (`V = {vol}` sites) Periodic Torus | "
            f"**Target:** $CP^{{{n_comp - 1}}}$ ($N = {n_comp}$) with Exact Auxiliary $U(1)$ "
            f"(`Best–Fisher` von Mises HB + `2× OR` + `6× SU(2)` Givens) | "
            f"**Equal Work Budget:** `{eq_work:,}` single-lattice sweeps/$\\beta$ "
            f"(`R = {n_rep}`, `L_d = {l_def}` defect slit) | "
            f"**Formal Companion:** [`FLOW/core/lean/Flow/ProjectiveBundle.lean`](../../FLOW/core/lean/Flow/ProjectiveBundle.lean)"
        ),
        "",
        f"## 1. Continuum Observables, Finite-Volume Ratio & Asymptotic Freedom Scaling ($L = {L}$)",
        "",
        (
            f"Leading large-$N$ asymptotics predict $\\chi_t \\xi_G^2 \\approx \\frac{{1}}{{2\\pi N}} "
            f"= {chi_xi2_asymp:.5f}$ for $N = {n_comp}$ (with small positive $1/N$ corrections). "
            f"Literature cross-check against Hasenbusch (2017, Table II for $N=21$ periodic BC): "
            f"$\\beta = 0.65$ ($L=32$, exact match) gives $\\xi_{{2\\text{{nd}}}} = 2.6993(7)$, $\\chi_t = 0.001170(9)$; "
            f"$\\beta = 0.70$ ($L=44$) gives $\\xi_{{2\\text{{nd}}}} = 3.744(2)$, $\\chi_t = 0.000570(16)$."
        ),
        "",
        "| $\\beta$ | Mean Link Energy $E$ (`Sym` vs `Loc`) | Corr. Length $\\xi_G$ | $L / \\xi_G$ | Cooled $\\chi_t$ (`Symbiotic`) | Continuum $\\chi_t \\xi_G^2$ | Cooled $\\chi_t$ (`local_or_metro`) |",
        "| ---: | :--- | ---: | ---: | :--- | ---: | :--- |",
        *md_obs_rows,
        "",
        f"## 2. Work-Normalized $\\tau_{{\\text{{int}}}}^{{\\text{{work}}}}(Q_{{\\text{{cool}}}})$ (Equal `{eq_work:,}` Replica-Sweeps) & Fisher `2-Jet` Ladder",
        "",
        f"| $\\beta$ | `local_or_metro` ($R=1$) | `bare_bp_inst` ($R=1$) | `symbiotic_ptbc_2jet` ($R={n_rep}$) | Equal-Work Speedup | Frozen Fisher `2-Jet` Ladder $c_0 \\dots c_{{{n_rep - 1}}}$ |",
        "| ---: | :--- | :--- | :--- | :---: | :--- |",
        *md_tau_rows,
        "",
        "## 3. Structural Takeaways Linked to Lean 4 (`Flow.ProjectiveBundle`)",
        "",
        f"1. **Hopf $U(1)$-Bundle & Exact Bogomol'nyi Bound (`Theorems P1–P2`)**: Auxiliary $U(1)$ heatbath + $SU(2) \\subset SU({n_comp})$ Givens updates preserve exact local gauge invariance and unit norm $\\|z_x\\| = 1$ for all $N$, while Berg–Lüscher geometric charge $Q_{{\\text{{geom}}}}$ bounds lattice action strictly from below.",
        f"2. **Corollary 2B on Non-Abelian Moduli (`bare_bp_inst`)**: Direct Belavin–Polyakov $SU(2)$ instanton proposals face ultraviolet component mismatch across thermal backgrounds, collapsing instanton acceptance below `1%` at $\\beta \\ge 0.65$.",
        f"3. **Defect-Slit Bianchi Opening & Fisher `2-Jet` Ladder (`Theorems P3–P4`)**: Opening an $L_d = {l_def}$ defect slit (`Theorem P3A`) with exact product-weight `@swap` detailed balance (`Theorem P3B`) and a Cauchy–Schwarz-optimal uniform Fisher arc ladder (`Theorem P4`) maintains balanced swap rates (`> 15%` min pair acceptance across all `{n_rep - 1}` links) while keeping $L/\\xi_G \\ge 5.0$ across $\\beta \\in \\{{{beta_lo:.2f} \\dots {beta_hi:.2f}\\}}$ and breaking $CP^{{{n_comp - 1}}}$ topological freezing at equal compute cost.",
        "",
    ]
    Path(args.report).write_text("\n".join(md))
    print(f"Generated {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
