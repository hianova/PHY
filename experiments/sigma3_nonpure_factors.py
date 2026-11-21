#!/usr/bin/env python3
"""
PHY Experiment: Why 3-Sigma Is "Not Pure Enough" — Winner's Curse, Look-Elsewhere,
SPAM Readout Asymmetry, Temporal Correlation, and Heavy-Tailed Glitches on [[5, 1, 3]]
=======================================================================================

Empirical complement to `scale_mismatch_2jet.py`:
  - EXP A1 (Direction A): Post-selection Winner's Curse overshoot & de-biasing
    (Naive 3σ vs Brute-Force 5σ vs Soft-Threshold 3σ vs Split-Sample 3σ)
  - EXP A2 (Direction A): Sequential Look-Elsewhere Effect (200 quiet epochs + 100 drift epochs)
    (Single-step 3σ/5σ vs Anytime-Valid Normal-Mixture e-Process vs Page-CUSUM)
  - EXP B1 (Direction B): Asymmetric SPAM readout corruption (exact PS cancellation)
    + axis-dependent SPAM falsifier wall
  - EXP B2 (Direction B): Temporally correlated AR(1) noise (FWER blowup vs causal whitening / LRV)
  - EXP B3 (A+B Combined): Realistic non-pure stress test (asymmetric SPAM + AR(1) + 3% glitches)
    showing Robust Huber-Clipped + Whitened e-Process + Median/Soft Parameter-Shift 2-Jet
"""

from __future__ import annotations

import math
import sys
import time
from pathlib import Path
from typing import List, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scale_mismatch_2jet import (  # noqa: E402
    _lap_var,
    apply_global_depolarizing_5q,
    apply_qec,
    apply_unitary,
    avg_logical_infidelity,
    build_513_code,
    zz01_unitary,
)


def eval_zz_infidelity_exact(code: dict, phi: float, p_incoh: float = 0.0) -> float:
    rhos = apply_unitary(code, code["rhos"], zz01_unitary(phi))
    if p_incoh > 0.0:
        rhos = apply_global_depolarizing_5q(code, rhos, p_incoh)
    rho_rec, _ = apply_qec(code, rhos)
    return avg_logical_infidelity(code, rho_rec)


class FastZZ513Evaluator:
    """Exact sin^2(phi/2) reduced evaluator verified to < 1e-12 against full 32x32 CPTP."""

    def __init__(self, code: dict) -> None:
        self.code = code
        self.coeffs: dict[float, Tuple[float, float]] = {}
        for p_incoh in (0.0, 1e-3):
            c0 = eval_zz_infidelity_exact(code, 0.0, p_incoh=p_incoh)
            c1 = eval_zz_infidelity_exact(code, math.pi, p_incoh=p_incoh)
            # Check exact agreement against full 32x32 density matrix across 10 arbitrary angles
            for phi in (0.001, 0.005, 0.01, 0.03, 0.05, 0.1, 0.25, 0.7, 1.4, 2.3):
                exact = eval_zz_infidelity_exact(code, phi, p_incoh=p_incoh)
                pred = c0 + (c1 - c0) * (math.sin(0.5 * phi) ** 2)
                assert abs(exact - pred) < 1e-12, f"Mismatch at phi={phi}, p={p_incoh}: {exact} vs {pred}"
            self.coeffs[p_incoh] = (c0, c1)

    def inf(self, residual_phi: float, p_incoh: float = 0.0) -> float:
        if p_incoh in self.coeffs:
            c0, c1 = self.coeffs[p_incoh]
            return c0 + (c1 - c0) * (math.sin(0.5 * residual_phi) ** 2)
        return eval_zz_infidelity_exact(self.code, residual_phi, p_incoh=p_incoh)


def probe_signal_nonpure(
    delta_true: float,
    theta_cal: float,
    baseline: float = 1e-4,
    e01: float = 0.0,
    e10: float = 0.0,
    additive_bias: float = 0.0,
    n_shots: int | None = None,
    rng: np.random.Generator | None = None,
) -> float:
    diff = delta_true - theta_cal
    p_ideal = baseline + (1.0 - 2.0 * baseline) * (math.sin(0.5 * diff) ** 2)
    p_meas = e01 + (1.0 - e01 - e10) * p_ideal + additive_bias
    p_meas = min(1.0, max(0.0, p_meas))
    if n_shots is not None and rng is not None:
        return float(rng.binomial(n_shots, p_meas) / n_shots)
    return p_meas


def probe_ps_nonpure(
    delta_true: float,
    theta: float,
    baseline: float = 1e-4,
    e01: float = 0.0,
    e10: float = 0.0,
    delta_e_axis: float = 0.0,
    bias_m: float = 0.0,
    bias_0: float = 0.0,
    bias_p: float = 0.0,
    n_shots: int | None = None,
    rng: np.random.Generator | None = None,
) -> Tuple[float, float, float, float, float, float, float]:
    h_ps = 0.5 * math.pi
    fm = probe_signal_nonpure(
        delta_true, theta - h_ps, baseline, e01, e10, -0.5 * delta_e_axis + bias_m, n_shots, rng
    )
    f0 = probe_signal_nonpure(delta_true, theta, baseline, e01, e10, bias_0, n_shots, rng)
    fp = probe_signal_nonpure(
        delta_true, theta + h_ps, baseline, e01, e10, 0.5 * delta_e_axis + bias_p, n_shots, rng
    )
    j0 = f0
    j1 = 0.5 * (fp - fm)
    j2 = 0.5 * (fp + fm - 2.0 * f0)
    if n_shots is None:
        return j0, j1, j2, 0.0, 0.0, float("inf"), float("inf")

    vm, v0, vp = _lap_var(fm, n_shots), _lap_var(f0, n_shots), _lap_var(fp, n_shots)
    se_j1 = 0.5 * math.sqrt(vp + vm)
    se_j2 = 0.5 * math.sqrt(vp + 4.0 * v0 + vm)
    snr_j1 = abs(j1) / se_j1 if se_j1 > 0.0 else 0.0
    snr_j2 = (j2 / se_j2) if se_j2 > 0.0 else 0.0
    return j0, j1, j2, se_j1, se_j2, snr_j1, snr_j2


def clip_dt(step: float, dt_max: float = 0.8) -> float:
    return max(-dt_max, min(dt_max, step))


def _logsumexp(xs: List[float]) -> float:
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs))


# ---------------------------------------------------------------------------
# EXP A1: Winner's Curse Post-Selection Overshoot & De-biasing
# ---------------------------------------------------------------------------


def run_exp_a1(ev: FastZZ513Evaluator, rng: np.random.Generator) -> None:
    print("\n=== EXP A1: Winner's Curse Post-Selection Overshoot (N_shots=10,000 total or split, 400 MC trials) ===")
    n_trials = 400
    n_shots = 10_000
    dt_max = 0.8

    for phi in (0.0, 0.01, 0.03):
        base_inf = ev.inf(phi)
        print(f"  -- True coherent ZZ drift phi={phi:.3f} rad (Pure QEC baseline inf={base_inf:.4e}) --")
        hdr = "    Controller                  | Mean Inf   | p95 Inf    | Hurt/400 | Act/400 | Mean(theta|act)"
        print(hdr)

        policies = [
            "Pure QEC",
            "Naive 3σ Schmitt (JS)",
            "Brute-Force 5σ Schmitt",
            "Soft-Threshold 3σ (z=2.0)",
            "Split-Sample 3σ (5k+5k)",
            "Split-Sample 3σ (10k+5k)",
        ]
        res = {p: [] for p in policies}
        acts = {p: [] for p in policies}

        for _ in range(n_trials):
            # 1. Pure QEC
            res["Pure QEC"].append(base_inf)

            # Single 10k PS probe shared conceptually by 1-sample methods
            j0, j1, j2, se1, _, snr1, snr2 = probe_ps_nonpure(phi, 0.0, n_shots=n_shots, rng=rng)

            # 2. Naive 3σ (James-Stein shrinkage + noisy Lyapunov check)
            if snr1 >= 3.0 and j2 > 1e-9 and snr2 >= 2.0:
                step = (-j1 / j2) * max(0.0, 1.0 - 1.0 / (snr1 * snr1))
                step = clip_dt(step, dt_max)
                f_cand = probe_signal_nonpure(phi, step, n_shots=n_shots, rng=rng)
                if f_cand <= j0 and abs(step) > 0.0:
                    acts["Naive 3σ Schmitt (JS)"].append(step)
                    res["Naive 3σ Schmitt (JS)"].append(ev.inf(phi - step))
                else:
                    res["Naive 3σ Schmitt (JS)"].append(base_inf)
            else:
                res["Naive 3σ Schmitt (JS)"].append(base_inf)

            # 3. Brute-Force 5σ
            if snr1 >= 5.0 and j2 > 1e-9 and snr2 >= 2.0:
                step = (-j1 / j2) * max(0.0, 1.0 - 1.0 / (snr1 * snr1))
                step = clip_dt(step, dt_max)
                acts["Brute-Force 5σ Schmitt"].append(step)
                res["Brute-Force 5σ Schmitt"].append(ev.inf(phi - step))
            else:
                res["Brute-Force 5σ Schmitt"].append(base_inf)

            # 4. Soft-Threshold 3σ (subtracts post-selection excess bias ~ 2.0 * SE)
            if snr1 >= 3.0 and j2 > 1e-9 and snr2 >= 2.0:
                j1_deb = math.copysign(max(0.0, abs(j1) - 2.0 * se1), j1)
                step = clip_dt(-j1_deb / j2, dt_max)
                if abs(step) > 0.0:
                    acts["Soft-Threshold 3σ (z=2.0)"].append(step)
                    res["Soft-Threshold 3σ (z=2.0)"].append(ev.inf(phi - step))
                else:
                    res["Soft-Threshold 3σ (z=2.0)"].append(base_inf)
            else:
                res["Soft-Threshold 3σ (z=2.0)"].append(base_inf)

            # 5. Split-Sample 3σ (5k detect + 5k independent estimate -> same 10k total budget!)
            _, _, _, _, _, snr_d5, _ = probe_ps_nonpure(phi, 0.0, n_shots=5_000, rng=rng)
            if snr_d5 >= 3.0:
                _, j1_e, j2_e, _, _, snr1_e, snr2_e = probe_ps_nonpure(phi, 0.0, n_shots=5_000, rng=rng)
                if j2_e > 1e-9 and snr2_e >= 2.0 and snr1_e > 0.0:
                    step = clip_dt((-j1_e / j2_e) * max(0.0, 1.0 - 1.0 / (snr1_e * snr1_e)), dt_max)
                    if abs(step) > 0.0:
                        acts["Split-Sample 3σ (5k+5k)"].append(step)
                        res["Split-Sample 3σ (5k+5k)"].append(ev.inf(phi - step))
                    else:
                        res["Split-Sample 3σ (5k+5k)"].append(base_inf)
                else:
                    res["Split-Sample 3σ (5k+5k)"].append(base_inf)
            else:
                res["Split-Sample 3σ (5k+5k)"].append(base_inf)

            # 6. Split-Sample 3σ (10k detect + 5k independent verify/estimate)
            if snr1 >= 3.0:
                _, j1_e, j2_e, _, _, snr1_e, snr2_e = probe_ps_nonpure(phi, 0.0, n_shots=5_000, rng=rng)
                if j2_e > 1e-9 and snr2_e >= 2.0 and snr1_e > 0.0:
                    step = clip_dt((-j1_e / j2_e) * max(0.0, 1.0 - 1.0 / (snr1_e * snr1_e)), dt_max)
                    if abs(step) > 0.0:
                        acts["Split-Sample 3σ (10k+5k)"].append(step)
                        res["Split-Sample 3σ (10k+5k)"].append(ev.inf(phi - step))
                    else:
                        res["Split-Sample 3σ (10k+5k)"].append(base_inf)
                else:
                    res["Split-Sample 3σ (10k+5k)"].append(base_inf)
            else:
                res["Split-Sample 3σ (10k+5k)"].append(base_inf)

        for p in policies:
            arr = np.array(res[p])
            hurt = int(np.sum(arr > base_inf * (1.0 + 1e-9)))
            n_act = len(acts[p])
            mean_th = f"{np.mean(acts[p]):+.4f}" if n_act > 0 else "   N/A "
            p95 = float(np.percentile(arr, 95))
            print(f"    {p:<27s} | {np.mean(arr):.4e} | {p95:.4e} |   {hurt:3d}    |   {n_act:3d}   | {mean_th}")


# ---------------------------------------------------------------------------
# EXP A2: Look-Elsewhere Effect & Weak Persistent Drift (e-Process vs 3σ/5σ)
# ---------------------------------------------------------------------------


def _mixture_e_value(log_pos: List[float], log_neg: List[float]) -> float:
    all_logs = log_pos + log_neg
    return math.exp(_logsumexp(all_logs) - math.log(len(all_logs)))


def _update_mixture_logs(log_pos: List[float], log_neg: List[float], lambdas: Tuple[float, ...], z: float) -> float:
    for i, lam in enumerate(lambdas):
        half_lam2 = 0.5 * lam * lam
        log_pos[i] += lam * z - half_lam2
        log_neg[i] += -lam * z - half_lam2
    return _mixture_e_value(log_pos, log_neg)


def run_exp_a2(ev: FastZZ513Evaluator, rng: np.random.Generator) -> None:
    print("\n=== EXP A2: Look-Elsewhere Effect (200 Quiet Epochs) & Weak Drift phi=0.01 (100 Epochs, N=2000, 100 runs) ===")
    n_runs = 100
    t_quiet = 200
    t_drift = 100
    n_shots = 2_000
    lambdas = (0.25, 0.5, 1.0, 2.0)

    names = [
        "Single-Step 3σ",
        "Single-Step 5σ",
        "Anytime e-Process (α=0.01)",
        "Page-CUSUM (k=0.5, h=6.5)",
    ]
    null_false_runs = {n: 0 for n in names}
    null_acts = {n: 0 for n in names}
    drift_det_runs = {n: 0 for n in names}
    drift_delays: dict[str, List[int]] = {n: [] for n in names}
    drift_final_infs: dict[str, List[float]] = {n: [] for n in names}

    # Phase 1: Pure Null phi = 0 for 200 epochs
    for _ in range(n_runs):
        lp = [0.0] * len(lambdas)
        ln = [0.0] * len(lambdas)
        sp, sn = 0.0, 0.0
        seen = {n: False for n in names}
        for _t in range(t_quiet):
            _, j1, _, _, _, snr1, _ = probe_ps_nonpure(0.0, 0.0, n_shots=n_shots, rng=rng)
            z = math.copysign(snr1, -j1)
            if abs(z) >= 3.0:
                seen["Single-Step 3σ"] = True
                null_acts["Single-Step 3σ"] += 1
            if abs(z) >= 5.0:
                seen["Single-Step 5σ"] = True
                null_acts["Single-Step 5σ"] += 1
            e_val = _update_mixture_logs(lp, ln, lambdas, z)
            if e_val >= 100.0:
                seen["Anytime e-Process (α=0.01)"] = True
                null_acts["Anytime e-Process (α=0.01)"] += 1
                lp = [0.0] * len(lambdas)
                ln = [0.0] * len(lambdas)
            sp = max(0.0, sp + z - 0.5)
            sn = max(0.0, sn - z - 0.5)
            if max(sp, sn) >= 6.5:
                seen["Page-CUSUM (k=0.5, h=6.5)"] = True
                null_acts["Page-CUSUM (k=0.5, h=6.5)"] += 1
                sp, sn = 0.0, 0.0
        for n in names:
            if seen[n]:
                null_false_runs[n] += 1

    # Phase 2: Weak persistent drift phi = 0.01 rad for 100 epochs (closed loop!)
    phi_drift = 0.01
    for _ in range(n_runs):
        thetas = {n: 0.0 for n in names}
        det_epoch: dict[str, int | None] = {n: None for n in names}

        lp = [0.0] * len(lambdas)
        ln = [0.0] * len(lambdas)
        win_e_j1, win_e_j2, win_e_var1 = [], [], []

        sp, sn = 0.0, 0.0
        win_c_j1, win_c_j2 = [], []

        for t in range(t_drift):
            # Shared epoch noise RNG sequence per controller independently evaluated at its theta
            # 3sigma
            _, j1, j2, _, _, snr1, snr2 = probe_ps_nonpure(phi_drift, thetas["Single-Step 3σ"], n_shots=n_shots, rng=rng)
            if snr1 >= 3.0 and j2 > 1e-9 and snr2 >= 2.0:
                if det_epoch["Single-Step 3σ"] is None:
                    det_epoch["Single-Step 3σ"] = t + 1
                thetas["Single-Step 3σ"] += clip_dt((-j1 / j2) * max(0.0, 1.0 - 1.0 / (snr1 * snr1)))

            # 5sigma
            _, j1, j2, _, _, snr1, snr2 = probe_ps_nonpure(phi_drift, thetas["Single-Step 5σ"], n_shots=n_shots, rng=rng)
            if snr1 >= 5.0 and j2 > 1e-9 and snr2 >= 2.0:
                if det_epoch["Single-Step 5σ"] is None:
                    det_epoch["Single-Step 5σ"] = t + 1
                thetas["Single-Step 5σ"] += clip_dt((-j1 / j2) * max(0.0, 1.0 - 1.0 / (snr1 * snr1)))

            # Anytime e-Process (window-averaged coherent correction upon stopping)
            _, j1, j2, se1, _, snr1, _ = probe_ps_nonpure(
                phi_drift, thetas["Anytime e-Process (α=0.01)"], n_shots=n_shots, rng=rng
            )
            z_e = (-j1 / se1) if se1 > 0.0 else 0.0
            win_e_j1.append(j1)
            win_e_j2.append(j2)
            win_e_var1.append(se1 * se1)
            if _update_mixture_logs(lp, ln, lambdas, z_e) >= 100.0:
                if det_epoch["Anytime e-Process (α=0.01)"] is None:
                    det_epoch["Anytime e-Process (α=0.01)"] = t + 1
                m_j1 = float(np.mean(win_e_j1))
                m_j2 = float(np.mean(win_e_j2))
                se_mean = math.sqrt(float(np.sum(win_e_var1))) / len(win_e_var1)
                snr_mean = abs(m_j1) / se_mean if se_mean > 0.0 else 0.0
                if m_j2 > 1e-9:
                    shrink = max(0.0, 1.0 - 1.0 / max(snr_mean * snr_mean, 1.0))
                    thetas["Anytime e-Process (α=0.01)"] += clip_dt((-m_j1 / m_j2) * shrink)
                lp = [0.0] * len(lambdas)
                ln = [0.0] * len(lambdas)
                win_e_j1, win_e_j2, win_e_var1 = [], [], []

            # Page-CUSUM
            _, j1, j2, se1, _, _, _ = probe_ps_nonpure(
                phi_drift, thetas["Page-CUSUM (k=0.5, h=6.5)"], n_shots=n_shots, rng=rng
            )
            z_c = (-j1 / se1) if se1 > 0.0 else 0.0
            win_c_j1.append(j1)
            win_c_j2.append(j2)
            sp = max(0.0, sp + z_c - 0.5)
            sn = max(0.0, sn - z_c - 0.5)
            if max(sp, sn) >= 6.5:
                if det_epoch["Page-CUSUM (k=0.5, h=6.5)"] is None:
                    det_epoch["Page-CUSUM (k=0.5, h=6.5)"] = t + 1
                m_j1 = float(np.mean(win_c_j1[-10:]))
                m_j2 = float(np.mean(win_c_j2[-10:]))
                if m_j2 > 1e-9:
                    thetas["Page-CUSUM (k=0.5, h=6.5)"] += clip_dt(-0.8 * m_j1 / m_j2)
                sp, sn = 0.0, 0.0
                win_c_j1, win_c_j2 = [], []

        for n in names:
            if det_epoch[n] is not None:
                drift_det_runs[n] += 1
                drift_delays[n].append(det_epoch[n])
            drift_final_infs[n].append(ev.inf(phi_drift - thetas[n]))

    base_drift = ev.inf(phi_drift)
    print(f"    (Weak drift Pure QEC baseline inf = {base_drift:.4e})")
    print("    Controller                 | Null FWER | Null Act/run | Drift Det% | Mean Delay | Final Drift Inf")
    for n in names:
        fwer = 100.0 * null_false_runs[n] / n_runs
        na = null_acts[n] / n_runs
        dr = 100.0 * drift_det_runs[n] / n_runs
        md = f"{np.mean(drift_delays[n]):5.1f} ep" if drift_delays[n] else "   N/A  "
        fi = float(np.mean(drift_final_infs[n]))
        print(f"    {n:<26s} |   {fwer:5.1f}%  |    {na:6.2f}    |   {dr:5.1f}%   |  {md}  | {fi:.4e}")


# ---------------------------------------------------------------------------
# EXP B1: Asymmetric SPAM Cancellation & Axis-Dependent SPAM Falsifier
# ---------------------------------------------------------------------------


def run_exp_b1(ev: FastZZ513Evaluator, rng: np.random.Generator) -> None:
    print("\n=== EXP B1: Asymmetric SPAM Readout (e01=0.02, e10=0.08, η=0.90) & Axis-Dependent SPAM Falsifier ===")
    e01, e10 = 0.02, 0.08
    n_trials = 100
    n_shots = 10_000

    print("  [B1a] Symmetric Asymmetric-SPAM invariance test:")
    print("    phi   | Naive j0>1e-3 trig% | PS Mean(theta_hat) | Residual Bias | [[5,1,3]] Inf (vs Base)")
    for phi in (0.0, 0.05, 0.15):
        naive_trig = 0
        thetas = []
        infs = []
        base = ev.inf(phi)
        for _ in range(n_trials):
            f0 = probe_signal_nonpure(phi, 0.0, e01=e01, e10=e10, n_shots=n_shots, rng=rng)
            if f0 > 1e-3:
                naive_trig += 1
            _, j1, j2, _, _, snr1, snr2 = probe_ps_nonpure(phi, 0.0, e01=e01, e10=e10, n_shots=n_shots, rng=rng)
            th = clip_dt(-j1 / j2) if (snr1 >= 3.0 and j2 > 1e-9 and snr2 >= 2.0) else 0.0
            thetas.append(th)
            infs.append(ev.inf(phi - th))
        m_th = float(np.mean(thetas))
        print(
            f"    {phi:5.2f} |       {100.0 * naive_trig / n_trials:5.1f}%        |"
            f"     {m_th:+8.5f}       |  {abs(m_th - phi):.2e}     | {np.mean(infs):.3e} (base {base:.3e})"
        )

    print("  [B1b Falsifier Wall] Axis-Dependent SPAM Δe_axis breaks Parameter-Shift invariance (N=20,000, phi=0.05):")
    phi = 0.05
    eta = 1.0 - e01 - e10
    print("    Δe_axis | Pred Bias ≈ -Δe/η | PS Mean(theta_hat) | Observed Bias | PS+QEC Inf (Pure QEC=4.166e-4)")
    for de in (0.0, 0.002, 0.010, 0.030):
        thetas = []
        infs = []
        for _ in range(n_trials):
            _, j1, j2, _, _, _, _ = probe_ps_nonpure(
                phi, 0.0, e01=e01, e10=e10, delta_e_axis=de, n_shots=20_000, rng=rng
            )
            th = clip_dt(-j1 / j2) if j2 > 1e-9 else 0.0
            thetas.append(th)
            infs.append(ev.inf(phi - th))
        m_th = float(np.mean(thetas))
        pred_bias = -de / eta
        print(
            f"    {de:7.3f} |      {pred_bias:+8.5f}      |      {m_th:+8.5f}       |"
            f"   {m_th - phi:+8.5f}    | {np.mean(infs):.4e}"
        )


# ---------------------------------------------------------------------------
# EXP B2: AR(1) Temporally Correlated Noise Null FWER Stress Test
# ---------------------------------------------------------------------------


def run_exp_b2(rng: np.random.Generator) -> None:
    print("\n=== EXP B2: AR(1) Correlated Noise (Exact N(0,1) Marginal, T=200 Quiet Epochs, 200 Runs) ===")
    n_runs = 200
    t_len = 200
    lambdas = (0.25, 0.5, 1.0, 2.0)

    ctrls = [
        "Single-Step 3σ",
        "Naive i.i.d. e-Process",
        "Whitened e-Process (oracle)",
        "Whitened e-Process (online ρ̂)",
    ]
    print("    rho  | " + " | ".join(f"{c}" for c in ctrls) + "  (Null FWER% / Act)")
    for rho in (0.0, 0.5, 0.85):
        fwer_cnt = {c: 0 for c in ctrls}
        act_cnt = {c: 0 for c in ctrls}
        innov_std = math.sqrt(max(1e-9, 1.0 - rho * rho))
        for _ in range(n_runs):
            seen = {c: False for c in ctrls}
            z = float(rng.normal(0.0, 1.0))

            lp_iid, ln_iid = [0.0] * len(lambdas), [0.0] * len(lambdas)
            lp_orc, ln_orc = [0.0] * len(lambdas), [0.0] * len(lambdas)
            lp_onl, ln_onl = [0.0] * len(lambdas), [0.0] * len(lambdas)

            z_prev = z
            sum_z2 = z * z + 6.0  # small Bayesian prior prevents early division instability
            sum_z_lag = 0.0

            for t in range(t_len):
                if t > 0:
                    z = rho * z + innov_std * float(rng.normal(0.0, 1.0))

                # 1. Single-Step 3σ
                if abs(z) >= 3.0:
                    seen["Single-Step 3σ"] = True
                    act_cnt["Single-Step 3σ"] += 1

                # 2. Naive i.i.d. e-Process
                if _update_mixture_logs(lp_iid, ln_iid, lambdas, z) >= 100.0:
                    seen["Naive i.i.d. e-Process"] = True
                    act_cnt["Naive i.i.d. e-Process"] += 1
                    lp_iid, ln_iid = [0.0] * len(lambdas), [0.0] * len(lambdas)

                # 3. Oracle AR(1)-whitened e-Process
                if t == 0:
                    w_orc = z
                else:
                    w_orc = (z - rho * z_prev) / innov_std
                if _update_mixture_logs(lp_orc, ln_orc, lambdas, w_orc) >= 100.0:
                    seen["Whitened e-Process (oracle)"] = True
                    act_cnt["Whitened e-Process (oracle)"] += 1
                    lp_orc, ln_orc = [0.0] * len(lambdas), [0.0] * len(lambdas)

                # 4. Online causal ρ̂ whitened e-Process (warm-up LRV inflated for first 12 epochs)
                if t == 0:
                    w_onl = z
                else:
                    rho_hat = max(0.0, min(0.88, sum_z_lag / sum_z2))
                    w_onl = (z - rho_hat * z_prev) / math.sqrt(max(1e-6, 1.0 - rho_hat * rho_hat))
                    if t < 12:
                        w_onl *= 0.7  # conservative warm-up damping while ρ̂ calibrates
                    sum_z_lag += z * z_prev
                    sum_z2 += z * z
                if _update_mixture_logs(lp_onl, ln_onl, lambdas, w_onl) >= 100.0:
                    seen["Whitened e-Process (online ρ̂)"] = True
                    act_cnt["Whitened e-Process (online ρ̂)"] += 1
                    lp_onl, ln_onl = [0.0] * len(lambdas), [0.0] * len(lambdas)

                z_prev = z

            for c in ctrls:
                if seen[c]:
                    fwer_cnt[c] += 1

        cells = [
            f"{100.0 * fwer_cnt[c] / n_runs:5.1f}% ({act_cnt[c] / n_runs:4.2f})" for c in ctrls
        ]
        print(f"    {rho:4.2f} | " + "     |    ".join(cells))


# ---------------------------------------------------------------------------
# EXP B3: Combined "Non-Pure" Stress Test (Asym SPAM + AR(1) + 3% Glitches)
# ---------------------------------------------------------------------------


def run_exp_b3(ev: FastZZ513Evaluator, rng: np.random.Generator) -> None:
    print(
        "\n=== EXP B3: Combined Non-Pure Stress Test (60 Epochs x 40 Runs, Asym SPAM + AR(1) ρ=0.6 + 3% ±0.25 Glitches) ==="
    )
    n_runs = 40
    t_len = 60
    n_shots = 2_000
    p_incoh = 1e-3
    e01, e10 = 0.02, 0.08
    rho_corr = 0.6
    corr_std = 0.008  # realistic ~1σ correlated background modulation
    glitch_prob = 0.03
    lambdas = (0.25, 0.5, 1.0, 2.0)

    def phi_schedule(t: int) -> float:
        if t < 20:
            return 0.002
        if t < 40:
            return 0.002 + (0.14 - 0.002) * (t - 20) / 19.0
        return 0.26

    ctrls = [
        "Pure QEC",
        "Single-Step 3σ Schmitt",
        "Single-Step 5σ Schmitt",
        "Unclipped e-Process",
        "Robust Huber-LRV e-Process + Median PS",
    ]
    quiet_infs = {c: [] for c in ctrls}
    drift_infs = {c: [] for c in ctrls}
    all_infs = {c: [] for c in ctrls}
    quiet_acts = {c: 0 for c in ctrls}
    tot_acts = {c: 0 for c in ctrls}

    innov_std = math.sqrt(1.0 - rho_corr * rho_corr)
    lrv_eff = 1.35  # combines i.i.d. shot noise + correlated AR(1) modulation
    # Note: 100% AR(1) would use sqrt((1+rho)/(1-rho))=2.0; combined variance gives moderate ~1.35 factor

    for _ in range(n_runs):
        thetas = {c: 0.0 for c in ctrls}
        sch_act = {"Single-Step 3σ Schmitt": False, "Single-Step 5σ Schmitt": False}

        lp_unc, ln_unc = [0.0] * len(lambdas), [0.0] * len(lambdas)
        win_unc_j1, win_unc_j2 = [], []

        lp_rob, ln_rob = [0.0] * len(lambdas), [0.0] * len(lambdas)
        win_rob_j1, win_rob_j2, win_rob_se1 = [], [], []

        u_ar = float(rng.normal(0.0, corr_std))
        for t in range(t_len):
            phi_t = phi_schedule(t)
            if t > 0:
                u_ar = rho_corr * u_ar + innov_std * float(rng.normal(0.0, corr_std))

            # Single-epoch heavy-tailed readout / cosmic-ray glitch on fp or fm
            g_m, g_p = 0.0, 0.0
            if float(rng.random()) < glitch_prob:
                g_sign = 0.25 if float(rng.random()) < 0.5 else -0.25
                if float(rng.random()) < 0.5:
                    g_p += g_sign
                else:
                    g_m += g_sign

            bias_m = -0.5 * u_ar + g_m
            bias_p = 0.5 * u_ar + g_p

            # Pure QEC
            inf_q = ev.inf(phi_t, p_incoh=p_incoh)
            all_infs["Pure QEC"].append(inf_q)
            (quiet_infs if t < 20 else drift_infs)["Pure QEC"].append(inf_q)

            # Evaluate each feedback controller cleanly with independent shot noise on identical physical environment
            for c_name, z_ent in (("Single-Step 3σ Schmitt", 3.0), ("Single-Step 5σ Schmitt", 5.0)):
                _, j1, j2, _, _, snr1, snr2 = probe_ps_nonpure(
                    phi_t, thetas[c_name], e01=e01, e10=e10, bias_m=bias_m, bias_p=bias_p, n_shots=n_shots, rng=rng
                )
                act = sch_act[c_name]
                if (not act and snr1 >= z_ent) or (act and snr1 >= 1.5):
                    sch_act[c_name] = True
                    if j2 > 1e-9 and snr2 >= 2.0:
                        step = clip_dt((-j1 / j2) * max(0.0, 1.0 - 1.0 / max(snr1 * snr1, 1.0)))
                        thetas[c_name] += step
                        tot_acts[c_name] += 1
                        if t < 20:
                            quiet_acts[c_name] += 1
                else:
                    sch_act[c_name] = False
                inf_c = ev.inf(phi_t - thetas[c_name], p_incoh=p_incoh)
                all_infs[c_name].append(inf_c)
                (quiet_infs if t < 20 else drift_infs)[c_name].append(inf_c)

            # Unclipped i.i.d. e-Process
            _, j1, j2, se1, _, _, _ = probe_ps_nonpure(
                phi_t, thetas["Unclipped e-Process"], e01=e01, e10=e10, bias_m=bias_m, bias_p=bias_p, n_shots=n_shots, rng=rng
            )
            z_u = (-j1 / se1) if se1 > 0.0 else 0.0
            win_unc_j1.append(j1)
            win_unc_j2.append(j2)
            if _update_mixture_logs(lp_unc, ln_unc, lambdas, z_u) >= 100.0:
                m_j1 = float(np.mean(win_unc_j1[-6:]))
                m_j2 = float(np.mean(win_unc_j2[-6:]))
                if m_j2 > 1e-9:
                    thetas["Unclipped e-Process"] += clip_dt(-m_j1 / m_j2)
                    tot_acts["Unclipped e-Process"] += 1
                    if t < 20:
                        quiet_acts["Unclipped e-Process"] += 1
                lp_unc, ln_unc = [0.0] * len(lambdas), [0.0] * len(lambdas)
                win_unc_j1, win_unc_j2 = [], []
            inf_u = ev.inf(phi_t - thetas["Unclipped e-Process"], p_incoh=p_incoh)
            all_infs["Unclipped e-Process"].append(inf_u)
            (quiet_infs if t < 20 else drift_infs)["Unclipped e-Process"].append(inf_u)

            # Robust Huber-Clipped LRV e-Process + Median/Soft PS 2-Jet
            _, j1, j2, se1, _, _, _ = probe_ps_nonpure(
                phi_t,
                thetas["Robust Huber-LRV e-Process + Median PS"],
                e01=e01,
                e10=e10,
                bias_m=bias_m,
                bias_p=bias_p,
                n_shots=n_shots,
                rng=rng,
            )
            z_raw = (-j1 / se1) if se1 > 0.0 else 0.0
            z_huber = max(-3.0, min(3.0, z_raw)) / lrv_eff
            win_rob_j1.append(j1)
            win_rob_j2.append(j2)
            win_rob_se1.append(se1)
            if _update_mixture_logs(lp_rob, ln_rob, lambdas, z_huber) >= 100.0:
                tail_j1 = win_rob_j1[-6:]
                tail_j2 = win_rob_j2[-6:]
                med_j1 = float(np.median(tail_j1))
                med_j2 = float(np.median(tail_j2))
                se_med = lrv_eff * float(np.median(win_rob_se1[-6:])) / math.sqrt(len(tail_j1))
                j1_deb = math.copysign(max(0.0, abs(med_j1) - 1.5 * se_med), med_j1)
                if med_j2 > 1e-9 and abs(j1_deb) > 0.0:
                    thetas["Robust Huber-LRV e-Process + Median PS"] += clip_dt(-j1_deb / med_j2)
                    tot_acts["Robust Huber-LRV e-Process + Median PS"] += 1
                    if t < 20:
                        quiet_acts["Robust Huber-LRV e-Process + Median PS"] += 1
                lp_rob, ln_rob = [0.0] * len(lambdas), [0.0] * len(lambdas)
                win_rob_j1, win_rob_j2, win_rob_se1 = [], [], []
            inf_r = ev.inf(phi_t - thetas["Robust Huber-LRV e-Process + Median PS"], p_incoh=p_incoh)
            all_infs["Robust Huber-LRV e-Process + Median PS"].append(inf_r)
            (quiet_infs if t < 20 else drift_infs)["Robust Huber-LRV e-Process + Median PS"].append(inf_r)

    print("    Controller                             | Quiet(t<20) | Drift+TLS(t>=20) | Overall Mean | Quiet FalseAct | Total Act/60")
    for c in ctrls:
        print(
            f"    {c:<38s} |  {np.mean(quiet_infs[c]):.3e}  |    {np.mean(drift_infs[c]):.3e}     |"
            f"  {np.mean(all_infs[c]):.3e}   |     {quiet_acts[c] / n_runs:5.2f}      |   {tot_acts[c] / n_runs:5.2f} / 60"
        )


def main() -> int:
    t0 = time.perf_counter()
    code = build_513_code()
    ev = FastZZ513Evaluator(code)
    print(
        "[SELF-CHECK PASS] [[5,1,3]] exact CPTP vs closed-form sin^2(phi/2) residual < 1e-12 across 20 checks (p_incoh=0, 1e-3)"
    )

    rng = np.random.default_rng(20261121)
    run_exp_a1(ev, rng)
    run_exp_a2(ev, rng)
    run_exp_b1(ev, rng)
    run_exp_b2(rng)
    run_exp_b3(ev, rng)

    print(f"\nCompleted non-pure 3σ benchmark in {time.perf_counter() - t0:.2f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
