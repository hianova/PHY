#!/usr/bin/env python3
"""
PHY Experiment: Scale-Mismatch Quantitative Benchmark on [[5, 1, 3]] Perfect Code
==================================================================================
Tests the core hypothesis connecting FLOW/docs, PHY, and QUA:
  - Microscopic coherent unitary drift lives on a continuous geometric manifold
    (controllable via FLOW 2-Jet curvature calibration).
  - Macroscopic QEC lives on discrete stabilizer syndrome sectors (vulnerable to
    unrelaxed continuous strain & weight-2 correlated mis-decoding).

Includes exact 32x32 density-matrix verification + explicit falsification tests
(incoherent decoherence floor, shot-noise over-calibration, inflection breakdown).
"""

from __future__ import annotations

import math
import time
from typing import Dict, List, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# 1. Exact [[5, 1, 3]] Perfect Code Construction & Self-Consistency Assertions
# ---------------------------------------------------------------------------

PAULIS: Dict[str, np.ndarray] = {
    "I": np.array([[1, 0], [0, 1]], dtype=np.complex128),
    "X": np.array([[0, 1], [1, 0]], dtype=np.complex128),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
    "Z": np.array([[1, 0], [0, -1]], dtype=np.complex128),
}


def kron_op(pauli_str: str) -> np.ndarray:
    op = PAULIS[pauli_str[0]]
    for ch in pauli_str[1:]:
        op = np.kron(op, PAULIS[ch])
    return op


def comm_sign(a: str, b: str) -> int:
    """Return 0 if Pauli strings commute, 1 if they anticommute."""
    anti = 0
    for ca, cb in zip(a, b):
        if ca != "I" and cb != "I" and ca != cb:
            anti ^= 1
    return anti


def build_513_code():
    I32 = np.eye(32, dtype=np.complex128)
    stab_strs = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]
    stabs = [kron_op(s) for s in stab_strs]
    xl = kron_op("XXXXX")
    zl = kron_op("ZZZZZ")

    # Verify stabilizer & logical algebra exactly
    for i, si in enumerate(stabs):
        assert np.allclose(si @ si, I32)
        assert np.allclose(si @ xl, xl @ si)
        assert np.allclose(si @ zl, zl @ si)
        for j, sj in enumerate(stabs):
            assert np.allclose(si @ sj, sj @ si)
    assert np.allclose(xl @ zl + zl @ xl, 0.0)

    # All 16 syndrome projectors Pi_s (s in 0..15)
    projs = []
    for s in range(16):
        p = I32.copy()
        for k in range(4):
            sgn = -1.0 if ((s >> k) & 1) else 1.0
            p = p @ (0.5 * (I32 + sgn * stabs[k]))
        projs.append(p)
    projs_arr = np.stack(projs, axis=0)
    assert np.isclose(np.trace(projs_arr[0]).real, 2.0)
    assert np.allclose(projs_arr.sum(axis=0), I32)

    # Exact unique weight-1 recovery lookup R_s for syndromes 1..15
    rec_ops = [None] * 16
    rec_ops[0] = I32
    seen_syndromes = set()
    for q in range(5):
        for p_ch in ("X", "Y", "Z"):
            estr = "".join(p_ch if idx == q else "I" for idx in range(5))
            syn = 0
            for k, s_str in enumerate(stab_strs):
                if comm_sign(estr, s_str):
                    syn |= 1 << k
            assert syn != 0 and syn not in seen_syndromes
            seen_syndromes.add(syn)
            rec_ops[syn] = kron_op(estr)
    assert seen_syndromes == set(range(1, 16))
    rec_arr = np.stack(rec_ops, axis=0)  # (16, 32, 32)
    kraus = np.einsum("kij,kjl->kil", rec_arr, projs_arr)  # (16, 32, 32)

    # Logical computational basis |0_L>, |1_L> inside P_C
    ket00000 = np.zeros(32, dtype=np.complex128)
    ket00000[0] = 1.0
    ket0_l = projs_arr[0] @ ket00000
    ket0_l /= np.linalg.norm(ket0_l)
    ket1_l = xl @ ket0_l
    assert np.allclose(projs_arr[0] @ ket0_l, ket0_l)
    assert np.allclose(zl @ ket0_l, ket0_l)
    assert np.allclose(zl @ ket1_l, -ket1_l)
    assert np.isclose(np.vdot(ket0_l, ket1_l), 0.0)

    # 6 cardinal Bloch states (exact 2-design for average state fidelity)
    s2 = 1.0 / math.sqrt(2.0)
    kets = [
        ket0_l,
        ket1_l,
        s2 * (ket0_l + ket1_l),
        s2 * (ket0_l - ket1_l),
        s2 * (ket0_l + 1j * ket1_l),
        s2 * (ket0_l - 1j * ket1_l),
    ]
    rhos = np.stack([np.outer(psi, psi.conj()) for psi in kets], axis=0)  # (6, 32, 32)
    return {
        "I32": I32,
        "stab_strs": stab_strs,
        "stabs": stabs,
        "xl": xl,
        "zl": zl,
        "projs": projs_arr,
        "kraus": kraus,
        "kets": np.stack(kets, axis=0),
        "rhos": rhos,
    }


def apply_unitary(code: dict, rhos: np.ndarray, u: np.ndarray) -> np.ndarray:
    return u @ rhos @ u.conj().T


def apply_qec(code: dict, rhos: np.ndarray) -> Tuple[np.ndarray, float]:
    """Apply syndrome measurement + recovery; return (rho_out, syndrome_trigger_prob)."""
    pi0 = code["projs"][0]
    p_zero = float(np.mean(np.real(np.einsum("ij,bji->b", pi0, rhos))))
    p_trig = max(0.0, min(1.0, 1.0 - p_zero))
    k = code["kraus"]
    # k: (16, 32, 32), rhos: (6, 32, 32) -> (16, 6, 32, 32) summed over 16 Kraus terms
    out = np.sum(k[:, None, :, :] @ rhos[None, :, :, :] @ k.conj().transpose(0, 2, 1)[:, None, :, :], axis=0)
    return out, p_trig


def apply_global_depolarizing_5q(code: dict, rhos: np.ndarray, p_incoh: float) -> np.ndarray:
    """Apply independent single-qubit depolarizing channel with error prob p_incoh per qubit."""
    if p_incoh <= 0.0:
        return rhos
    cur = rhos
    p_pauli = p_incoh / 3.0
    p_keep = 1.0 - p_incoh
    for q in range(5):
        ops = [
            (p_keep, code["I32"]),
            (p_pauli, kron_op("".join("X" if i == q else "I" for i in range(5)))),
            (p_pauli, kron_op("".join("Y" if i == q else "I" for i in range(5)))),
            (p_pauli, kron_op("".join("Z" if i == q else "I" for i in range(5)))),
        ]
        nxt = np.zeros_like(cur)
        for prob, op in ops:
            nxt += prob * (op @ cur @ op.conj().T)
        cur = nxt
    return cur


def avg_logical_infidelity(code: dict, rhos: np.ndarray) -> float:
    """Exact 1 - average Bloch state fidelity over all 6 cardinal states."""
    kets = code["kets"]
    fids = np.real(np.einsum("bi,bij,bj->b", kets.conj(), rhos, kets))
    return float(max(0.0, 1.0 - np.mean(fids)))


def rx_layer_unitary(theta: float) -> np.ndarray:
    c = math.cos(0.5 * theta)
    s = -1j * math.sin(0.5 * theta)
    u1 = np.array([[c, s], [s, c]], dtype=np.complex128)
    u = u1
    for _ in range(4):
        u = np.kron(u, u1)
    return u


def zz01_unitary(phi: float) -> np.ndarray:
    c = math.cos(0.5 * phi)
    s = -1j * math.sin(0.5 * phi)
    zz = PAULIS["Z"]
    z0z1 = np.kron(zz, np.kron(zz, np.eye(8, dtype=np.complex128)))
    return c * np.eye(32, dtype=np.complex128) + s * z0z1


# ---------------------------------------------------------------------------
# 2. Microscopic FLOW 2-Jet Calibration Operators (vs Naive Newton)
# ---------------------------------------------------------------------------


def probe_error_signal(
    delta_true: float,
    theta_cal: float,
    baseline: float = 1e-4,
    n_shots: int | None = None,
    rng: np.random.Generator | None = None,
) -> float:
    """Microscopic gate-echo / Ramsey error probability: baseline + (1 - 2*baseline)*sin^2((delta-theta)/2)."""
    diff = delta_true - theta_cal
    p = baseline + (1.0 - 2.0 * baseline) * (math.sin(0.5 * diff) ** 2)
    p = min(1.0, max(0.0, p))
    if n_shots is not None and rng is not None:
        return float(rng.binomial(n_shots, p) / n_shots)
    return p


def _lap_var(p_hat: float, n_shots: int) -> float:
    k = p_hat * n_shots
    p_t = (k + 0.5) / (n_shots + 1.0)
    return p_t * (1.0 - p_t) / (n_shots + 1.0)


def probe_jets_fd_with_se(
    delta_true: float,
    theta: float,
    baseline: float = 1e-4,
    h: float = 0.05,
    n_shots: int | None = None,
    rng: np.random.Generator | None = None,
) -> Tuple[float, float, float, float, float]:
    """Classical narrow finite-difference 3-point probe (demonstrates 1/h^2 shot-noise blowup)."""
    fm = probe_error_signal(delta_true, theta - h, baseline, n_shots, rng)
    f0 = probe_error_signal(delta_true, theta, baseline, n_shots, rng)
    fp = probe_error_signal(delta_true, theta + h, baseline, n_shots, rng)
    j0 = f0
    j1 = (fp - fm) / (2.0 * h)
    j2 = (fp - 2.0 * f0 + fm) / (h * h)
    if n_shots is None:
        return j0, j1, j2, float("inf"), float("inf")

    vm, v0, vp = _lap_var(fm, n_shots), _lap_var(f0, n_shots), _lap_var(fp, n_shots)
    se_j1 = math.sqrt(vp + vm) / (2.0 * h)
    se_j2 = math.sqrt(vp + 4.0 * v0 + vm) / (h * h)
    snr_j1 = abs(j1) / se_j1 if se_j1 > 0.0 else 0.0
    snr_j2 = (j2 / se_j2) if se_j2 > 0.0 else 0.0
    return j0, j1, j2, snr_j1, snr_j2


def probe_jets_ps_with_se(
    delta_true: float,
    theta: float,
    baseline: float = 1e-4,
    n_shots: int | None = None,
    rng: np.random.Generator | None = None,
) -> Tuple[float, float, float, float, float]:
    """Exact Quantum Parameter-Shift (h=pi/2) 2-Jet probe: zero h^2 bias & ~300x lower j2 variance."""
    h_ps = 0.5 * math.pi
    fm = probe_error_signal(delta_true, theta - h_ps, baseline, n_shots, rng)
    f0 = probe_error_signal(delta_true, theta, baseline, n_shots, rng)
    fp = probe_error_signal(delta_true, theta + h_ps, baseline, n_shots, rng)
    j0 = f0
    j1 = 0.5 * (fp - fm)
    j2 = 0.5 * (fp + fm - 2.0 * f0)
    if n_shots is None:
        return j0, j1, j2, float("inf"), float("inf")

    vm, v0, vp = _lap_var(fm, n_shots), _lap_var(f0, n_shots), _lap_var(fp, n_shots)
    se_j1 = 0.5 * math.sqrt(vp + vm)
    se_j2 = 0.5 * math.sqrt(vp + 4.0 * v0 + vm)
    snr_j1 = abs(j1) / se_j1 if se_j1 > 0.0 else 0.0
    snr_j2 = (j2 / se_j2) if se_j2 > 0.0 else 0.0
    return j0, j1, j2, snr_j1, snr_j2


def _flow_vertex_step(
    j0: float,
    j1: float,
    j2: float,
    snr_j2: float,
    dt_max: float,
    n_shots: int | None,
    z_curv: float,
    vertex_mode: bool = True,
) -> float:
    if j2 > 1e-9 and (n_shots is None or snr_j2 >= z_curv):
        if vertex_mode:
            return -j1 / j2
        disc = j1 * j1 - 2.0 * j2 * j0
        if disc < 0.0:
            return -j1 / j2
        q_cit = j1 + math.copysign(math.sqrt(disc), j1 if j1 != 0.0 else 1.0)
        return (-2.0 * j0 / q_cit) if abs(q_cit) > 1e-12 else (-j1 / j2)
    return -math.copysign(min(dt_max, abs(j1) * 2.0 + 0.4), j1 if abs(j1) > 1e-12 else -1.0)


def run_2jet_schmitt_step(
    delta_true: float,
    theta: float,
    active: bool,
    dt_max: float = 0.8,
    baseline: float = 1e-4,
    h: float = 0.05,
    use_ps: bool = True,
    n_shots: int | None = None,
    rng: np.random.Generator | None = None,
    z_exit: float = 1.5,
    z_enter: float = 3.0,
    z_curv: float = 2.0,
    shrinkage: bool = True,
) -> Tuple[float, bool, bool]:
    """Single step of FLOW SNR-Schmitt guarded vertex 2-Jet calibration -> (theta_new, active_new, actuated)."""
    if use_ps:
        j0, j1, j2, snr_j1, snr_j2 = probe_jets_ps_with_se(
            delta_true, theta, baseline=baseline, n_shots=n_shots, rng=rng
        )
    else:
        j0, j1, j2, snr_j1, snr_j2 = probe_jets_fd_with_se(
            delta_true, theta, baseline=baseline, h=h, n_shots=n_shots, rng=rng
        )
    if not active:
        if snr_j1 < z_enter:
            return theta, False, False
        active_next = True
    else:
        if snr_j1 < z_exit:
            return theta, False, False
        active_next = True

    step = _flow_vertex_step(
        j0, j1, j2, snr_j2, dt_max=dt_max, n_shots=n_shots, z_curv=z_curv, vertex_mode=True
    )

    if n_shots is not None and shrinkage and snr_j1 > 0.0 and math.isfinite(snr_j1):
        shrink = max(0.0, 1.0 - 1.0 / (snr_j1 * snr_j1))
        step *= shrink

    step = max(-dt_max, min(dt_max, step))
    if abs(step) < 1e-15:
        return theta, active_next, False

    theta_cand = theta + step
    f_cand = probe_error_signal(delta_true, theta_cand, baseline=baseline, n_shots=n_shots, rng=rng)
    if f_cand <= j0:
        return theta_cand, active_next, True
    return theta, active_next, False


def run_2jet_calibration(
    delta_true: float,
    steps: int = 1,
    guarded: bool = True,
    dt_max: float = 0.8,
    baseline: float = 1e-4,
    h: float = 0.05,
    use_ps: bool = True,
    n_shots: int | None = None,
    rng: np.random.Generator | None = None,
    snr_schmitt: bool = False,
    z_exit: float = 1.5,
    z_enter: float = 3.0,
    z_curv: float = 2.0,
    shrinkage: bool = True,
) -> float:
    """
    Reproduces FLOW/core/flow.h 2-Jet mechanics:
      - Un-guarded (Naive Newton): theta <- theta - j1 / j2 directly
      - Guarded (FLOW vertex + Trust-Region + Lyapunov)
      - Guarded + SNR Schmitt (adds FLOW [b_exit, b_enter] SNR hysteresis + shrinkage)
    """
    theta = 0.0
    active = False
    for _ in range(steps):
        if snr_schmitt:
            theta, active, _ = run_2jet_schmitt_step(
                delta_true,
                theta,
                active,
                dt_max=dt_max,
                baseline=baseline,
                h=h,
                use_ps=use_ps,
                n_shots=n_shots,
                rng=rng,
                z_exit=z_exit,
                z_enter=z_enter,
                z_curv=z_curv,
                shrinkage=shrinkage,
            )
            continue

        if use_ps:
            j0, j1, j2, snr_j1, snr_j2 = probe_jets_ps_with_se(
                delta_true, theta, baseline=baseline, n_shots=n_shots, rng=rng
            )
        else:
            j0, j1, j2, snr_j1, snr_j2 = probe_jets_fd_with_se(
                delta_true, theta, baseline=baseline, h=h, n_shots=n_shots, rng=rng
            )
        _ = snr_j1
        if not guarded:
            step = 0.0 if abs(j2) < 1e-12 else (-j1 / j2)
            theta += step
            continue

        step = _flow_vertex_step(
            j0, j1, j2, snr_j2, dt_max=dt_max, n_shots=n_shots, z_curv=z_curv, vertex_mode=True
        )
        step = max(-dt_max, min(dt_max, step))
        theta_cand = theta + step
        f_cand = probe_error_signal(delta_true, theta_cand, baseline=baseline, n_shots=n_shots, rng=rng)
        if abs(f_cand) <= abs(j0) or (n_shots is not None and abs(f_cand - j0) < 3.0 / math.sqrt(n_shots)):
            if abs(f_cand) <= abs(j0):
                theta = theta_cand
            elif n_shots is not None and f_cand <= j0 + 2.0 * math.sqrt(max(j0 * (1.0 - j0), 1e-6) / n_shots):
                theta = theta_cand
    return theta


def fit_loglog_slope(xs: List[float], ys: List[float]) -> float:
    pts = [(math.log10(x), math.log10(y)) for x, y in zip(xs, ys) if x > 0 and y > 1e-14]
    if len(pts) < 2:
        return float("nan")
    lx = np.array([p[0] for p in pts])
    ly = np.array([p[1] for p in pts])
    slope, _ = np.polyfit(lx, ly, 1)
    return float(slope)


# ---------------------------------------------------------------------------
# 3. Four Quantitative Experiments (+ Falsification Regimes)
# ---------------------------------------------------------------------------


def main() -> int:
    t0 = time.perf_counter()
    code = build_513_code()

    # Exact verification that all 15 weight-1 errors are corrected to fidelity 1.0
    for q in range(5):
        for ch in ("X", "Y", "Z"):
            op = kron_op("".join(ch if i == q else "I" for i in range(5)))
            rho_rec, ptrig = apply_qec(code, apply_unitary(code, code["rhos"], op))
            assert abs(ptrig - 1.0) < 1e-12
            assert avg_logical_infidelity(code, rho_rec) < 1e-12
    print("[SELF-CHECK PASS] [[5,1,3]] exact CPTP stabilizer recovery: 15/15 weight-1 Pauli errors infidelity < 1e-12")

    rng = np.random.default_rng(2026)

    # ------------------------------------------------------------------
    # EXP 1: Weight-1 Coherent Over-Rotation Rx(delta)^5 (1 round & 20 rounds)
    # ------------------------------------------------------------------
    print("\n=== EXP 1: Microscopic Rx(delta)^5 Over-Rotation (1 Round & T=20 Rounds Unmitigated Accumulation) ===")
    deltas = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
    inf_qec_1, inf_2j1_1, inf_2j2_1, inf_qec_20, inf_2j1_20 = [], [], [], [], []
    print("  delta | QEC_1R (trig)      | 2Jet-1step + QEC  | Noisy2J(1e4)+QEC | Schmitt2J(1e4)+QEC | 20R_QEC    | 20R_2J1+QEC")
    for d in deltas:
        u_err = rx_layer_unitary(d)
        rho_qec, trig_q = apply_qec(code, apply_unitary(code, code["rhos"], u_err))
        lq1 = avg_logical_infidelity(code, rho_qec)

        th1 = run_2jet_calibration(d, steps=1, guarded=True)
        rho_2j1, trig_2j1 = apply_qec(code, apply_unitary(code, code["rhos"], rx_layer_unitary(d - th1)))
        l2j1 = avg_logical_infidelity(code, rho_2j1)

        th2 = run_2jet_calibration(d, steps=2, guarded=True)
        rho_2j2, _ = apply_qec(code, apply_unitary(code, code["rhos"], rx_layer_unitary(d - th2)))
        l2j2 = avg_logical_infidelity(code, rho_2j2)

        # Noisy finite-shot 2-Jet vs SNR Schmitt (N_shots = 10,000, mean over 40 MC trials)
        mc_vals, mc_sch = [], []
        for _ in range(40):
            th_n = run_2jet_calibration(d, steps=1, guarded=True, n_shots=10_000, rng=rng)
            r_n, _ = apply_qec(code, apply_unitary(code, code["rhos"], rx_layer_unitary(d - th_n)))
            mc_vals.append(avg_logical_infidelity(code, r_n))

            th_s = run_2jet_calibration(d, steps=1, guarded=True, n_shots=10_000, rng=rng, snr_schmitt=True)
            r_s, _ = apply_qec(code, apply_unitary(code, code["rhos"], rx_layer_unitary(d - th_s)))
            mc_sch.append(avg_logical_infidelity(code, r_s))
        l_noisy = float(np.mean(mc_vals))
        l_schmitt = float(np.mean(mc_sch))

        # 20 rounds of interleaved rotation + QEC (no re-calibration vs calibrated once)
        r20_q, r20_2j = code["rhos"], code["rhos"]
        u_res = rx_layer_unitary(d - th1)
        for _ in range(20):
            r20_q, _ = apply_qec(code, apply_unitary(code, r20_q, u_err))
            r20_2j, _ = apply_qec(code, apply_unitary(code, r20_2j, u_res))
        l20_q = avg_logical_infidelity(code, r20_q)
        l20_2j = avg_logical_infidelity(code, r20_2j)

        inf_qec_1.append(lq1)
        inf_2j1_1.append(l2j1)
        inf_2j2_1.append(l2j2)
        inf_qec_20.append(l20_q)
        inf_2j1_20.append(l20_2j)
        print(
            f"  {d:5.3f} | {lq1:9.3e} ({trig_q:6.4f}) | {l2j1:9.3e} ({trig_2j1:5.1e}) | "
            f"{l_noisy:9.3e}        | {l_schmitt:9.3e}          | {l20_q:9.3e}  | {l20_2j:9.3e}"
        )
    print(
        f"  -> Empirical scaling exponent (1-round): Pure QEC slope ~ {fit_loglog_slope(deltas, inf_qec_1):.2f}, "
        f"1-step 2-Jet slope ~ {fit_loglog_slope(deltas[3:], inf_2j1_1[3:]):.2f}"
    )

    # ------------------------------------------------------------------
    # EXP 2: Weight-2 Coherent ZZ Crosstalk exp(-i phi/2 Z0 Z1)
    # ------------------------------------------------------------------
    print("\n=== EXP 2: Weight-2 Correlated ZZ Crosstalk exp(-i*phi/2 * Z0Z1) -- Mis-Decoding Catastrophe ===")
    zz_qec, zz_2j1 = [], []
    print("  phi   | QEC_1R (trig / misdecode!) | 2Jet-1step + QEC   | Improvement factor")
    for p in deltas:
        rho_q, t_q = apply_qec(code, apply_unitary(code, code["rhos"], zz01_unitary(p)))
        lq = avg_logical_infidelity(code, rho_q)

        th1 = run_2jet_calibration(p, steps=1, guarded=True)
        rho_2j, _ = apply_qec(code, apply_unitary(code, code["rhos"], zz01_unitary(p - th1)))
        l2j = avg_logical_infidelity(code, rho_2j)

        zz_qec.append(lq)
        zz_2j1.append(l2j)
        ratio_str = f"{lq / max(l2j, 1e-18):.2e}x"
        print(f"  {p:5.3f} | {lq:10.3e} ({t_q:8.6f})      | {l2j:10.3e}         | {ratio_str}")
    print(
        f"  -> ZZ scaling exponent: Pure QEC slope ~ {fit_loglog_slope(deltas, zz_qec):.2f} (O(phi^2)!), "
        f"1-step 2-Jet slope ~ {fit_loglog_slope(deltas[2:], zz_2j1[2:]):.2f}"
    )

    # ------------------------------------------------------------------
    # EXP 3: Falsification Boundary -- Mixed Coherent + Incoherent Depolarizing Floor
    # ------------------------------------------------------------------
    print("\n=== EXP 3: Physical Wall / Falsification -- Pure Incoherent vs Coherent+Depolarizing Crossover ===")
    print("  [Falsifier 3A] Pure incoherent (delta=0, ZZ=0): 2-Jet unitary calibration cannot fix stochastic decoherence:")
    for pin in (1e-3, 1e-2, 3e-2):
        rq, _ = apply_qec(code, apply_global_depolarizing_5q(code, code["rhos"], pin))
        lq = avg_logical_infidelity(code, rq)
        # Apply a zero-signal 2-Jet probe run (ideal zero drift) -> zero unitary effect
        assert abs(run_2jet_calibration(0.0, steps=1, guarded=True)) < 1e-12
        print(f"    p_incoh={pin:5.3f} -> Pure QEC inf={lq:.4e}, 2-Jet gain=1.00x (strict ceiling!)")

    print("  [Crossover 3B] Fixed incoherent floor p_incoh=1e-3 + varying ZZ coherent phi:")
    p_floor = 1e-3
    for phi in (0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2):
        r_pure_in = apply_global_depolarizing_5q(code, apply_unitary(code, code["rhos"], zz01_unitary(phi)), p_floor)
        r_pure_qec, _ = apply_qec(code, r_pure_in)
        lp = avg_logical_infidelity(code, r_pure_qec)

        th1 = run_2jet_calibration(phi, steps=1, guarded=True)
        r_2j_in = apply_global_depolarizing_5q(code, apply_unitary(code, code["rhos"], zz01_unitary(phi - th1)), p_floor)
        r_2j_qec, _ = apply_qec(code, r_2j_in)
        l2 = avg_logical_infidelity(code, r_2j_qec)
        print(f"    phi={phi:5.3f}, p_incoh=0.001 | Pure QEC={lp:.4e} | 2-Jet+QEC={l2:.4e} | gain={lp/l2:6.2f}x")

    # ------------------------------------------------------------------
    # EXP 4: Curvature Breakdown (Large Drift > pi/2) & Shot-Noise Over-Calibration Falsifier
    # ------------------------------------------------------------------
    print("\n=== EXP 4: Failure & Falsification Boundaries of 2-Jet Calibration ===")
    print("  [4A] Non-quadratic inflection breakdown (|delta| crossing pi/2 = 1.571, using 4 steps):")
    for d in (0.3, 0.8, 1.3, 1.6, 2.0, 2.5):
        th_naive = run_2jet_calibration(d, steps=4, guarded=False)
        th_flow = run_2jet_calibration(d, steps=4, guarded=True, dt_max=0.6)
        err_naive = abs((d - th_naive + math.pi) % (2 * math.pi) - math.pi)
        err_flow = abs((d - th_flow + math.pi) % (2 * math.pi) - math.pi)
        tag = "BLOWUP!" if err_naive > 0.1 else "ok"
        print(f"    delta={d:4.1f} | Naive 4-step residual={err_naive:8.4e} ({tag:7s}) | FLOW Guarded residual={err_flow:8.4e}")

    print("  [4B] Shot-noise over-calibration falsifier vs FLOW SNR-Schmitt Gatekeeper (ZZ crosstalk, 100 MC trials):")
    for phi_test in (0.01, 0.05):
        r_base, _ = apply_qec(code, apply_unitary(code, code["rhos"], zz01_unitary(phi_test)))
        inf_base = avg_logical_infidelity(code, r_base)
        print(f"    -- True coherent ZZ drift phi={phi_test:.2f} rad (Pure QEC baseline inf={inf_base:.4e}) --")
        print("    N_shots | Ungated mean (p95)       hurt  | Schmitt mean (p95)       hurt | act/100 | Schmitt vs Base")
        for n_sh in (100, 1000, 10_000, 100_000):
            infs_u, infs_s = [], []
            hurt_u, hurt_s, act_s = 0, 0, 0
            for _ in range(100):
                th_u = run_2jet_calibration(phi_test, steps=1, guarded=True, n_shots=n_sh, rng=rng, snr_schmitt=False)
                r_u, _ = apply_qec(code, apply_unitary(code, code["rhos"], zz01_unitary(phi_test - th_u)))
                vu = avg_logical_infidelity(code, r_u)
                infs_u.append(vu)
                if vu > inf_base + 1e-12:
                    hurt_u += 1

                th_s, _, did_act = run_2jet_schmitt_step(
                    phi_test, 0.0, False, n_shots=n_sh, rng=rng
                )
                if did_act:
                    act_s += 1
                r_s, _ = apply_qec(code, apply_unitary(code, code["rhos"], zz01_unitary(phi_test - th_s)))
                vs = avg_logical_infidelity(code, r_s)
                infs_s.append(vs)
                if vs > inf_base + 1e-12:
                    hurt_s += 1

            mu, p95u = float(np.mean(infs_u)), float(np.percentile(infs_u, 95))
            ms, p95s = float(np.mean(infs_s)), float(np.percentile(infs_s, 95))
            if ms > inf_base * 1.02:
                verdict = f"HURTS ({ms/inf_base:.2f}x)"
            elif abs(ms - inf_base) <= inf_base * 0.02:
                verdict = "SAFE SLEEP (1.00x)"
            else:
                verdict = f"WINS ({inf_base/ms:.2f}x)"
            print(
                f"    {n_sh:7d} | {mu:9.3e} ({p95u:9.3e}) {hurt_u:3d}  | "
                f"{ms:9.3e} ({p95s:9.3e}) {hurt_s:3d}  | {act_s:5d}   | {verdict}"
            )

    # ------------------------------------------------------------------
    # EXP 5: Dynamic 60-Epoch Coherent ZZ Drift + TLS Step-Jump Tracking
    # ------------------------------------------------------------------
    print("\n=== EXP 5: Dynamic 60-Epoch Coherent ZZ Drift + TLS Step-Jump (N_shots=2000/epoch, p_incoh=1e-3, 20 runs) ===")
    phi_traj = np.empty(60, dtype=np.float64)
    for t in range(60):
        if t < 15:
            phi_traj[t] = 0.003
        elif t < 35:
            phi_traj[t] = 0.003 + (0.15 - 0.003) * ((t - 14) / 20.0)
        elif t < 40:
            phi_traj[t] = 0.27
        else:
            phi_traj[t] = 0.27

    ctrl_names = (
        "Pure QEC",
        "Always-On 2Jet (PS)",
        "FLOW Schmitt (FD 0.05)",
        "Single-Thresh(2.0 PS)",
        "FLOW Schmitt(1.5,3.0 PS)",
    )
    quiet_infs = {k: [] for k in ctrl_names}
    active_infs = {k: [] for k in ctrl_names}
    total_infs = {k: [] for k in ctrl_names}
    acts_tot = {k: 0 for k in ctrl_names}
    chatter_tot = {k: 0 for k in ctrl_names}

    n_runs = 20
    n_shots_dyn = 2000
    # Precompute Pure QEC trajectory (deterministic)
    for t in range(60):
        rp, _ = apply_qec(
            code, apply_global_depolarizing_5q(code, apply_unitary(code, code["rhos"], zz01_unitary(float(phi_traj[t]))), 1e-3)
        )
        ip = avg_logical_infidelity(code, rp)
        total_infs["Pure QEC"].append(ip)
        (quiet_infs["Pure QEC"] if t < 15 else active_infs["Pure QEC"]).append(ip)

    for _ in range(n_runs):
        st_always = 0.0
        st_fd, act_fd = 0.0, False
        st_single, act_single = 0.0, False
        st_sch, act_sch = 0.0, False
        for t in range(60):
            pt = float(phi_traj[t])

            # 1) Always-On un-gated Parameter-Shift 2-Jet
            j0, j1, j2, _, snr_j2 = probe_jets_ps_with_se(pt, st_always, n_shots=n_shots_dyn, rng=rng)
            sp = _flow_vertex_step(
                j0, j1, j2, snr_j2, dt_max=0.8, n_shots=n_shots_dyn, z_curv=2.0, vertex_mode=True
            )
            sp = max(-0.8, min(0.8, sp))
            fc = probe_error_signal(pt, st_always + sp, n_shots=n_shots_dyn, rng=rng)
            if fc <= j0 + 2.0 * math.sqrt(max(j0 * (1.0 - j0), 1e-6) / n_shots_dyn):
                st_always += sp
            acts_tot["Always-On 2Jet (PS)"] += 1

            # 2) Classical narrow FD (h=0.05) FLOW Schmitt baseline
            prev_fd = act_fd
            st_fd, act_fd, did_fd = run_2jet_schmitt_step(
                pt, st_fd, act_fd, use_ps=False, h=0.05, n_shots=n_shots_dyn, rng=rng, z_exit=1.5, z_enter=3.0
            )
            if did_fd:
                acts_tot["FLOW Schmitt (FD 0.05)"] += 1
            if act_fd != prev_fd:
                chatter_tot["FLOW Schmitt (FD 0.05)"] += 1

            # 3) Single-threshold (no hysteresis) Parameter-Shift z_enter = z_exit = 2.0
            prev_s = act_single
            st_single, act_single, did_single = run_2jet_schmitt_step(
                pt, st_single, act_single, use_ps=True, n_shots=n_shots_dyn, rng=rng, z_exit=2.0, z_enter=2.0
            )
            if did_single:
                acts_tot["Single-Thresh(2.0 PS)"] += 1
            if act_single != prev_s:
                chatter_tot["Single-Thresh(2.0 PS)"] += 1

            # 4) FLOW SNR-Schmitt Parameter-Shift (hysteresis z_exit=1.5, z_enter=3.0)
            prev_f = act_sch
            st_sch, act_sch, did_sch = run_2jet_schmitt_step(
                pt, st_sch, act_sch, use_ps=True, n_shots=n_shots_dyn, rng=rng, z_exit=1.5, z_enter=3.0
            )
            if did_sch:
                acts_tot["FLOW Schmitt(1.5,3.0 PS)"] += 1
            if act_sch != prev_f:
                chatter_tot["FLOW Schmitt(1.5,3.0 PS)"] += 1

            for name, th_val in (
                ("Always-On 2Jet (PS)", st_always),
                ("FLOW Schmitt (FD 0.05)", st_fd),
                ("Single-Thresh(2.0 PS)", st_single),
                ("FLOW Schmitt(1.5,3.0 PS)", st_sch),
            ):
                rr, _ = apply_qec(
                    code,
                    apply_global_depolarizing_5q(
                        code, apply_unitary(code, code["rhos"], zz01_unitary(pt - th_val)), 1e-3
                    ),
                )
                val = avg_logical_infidelity(code, rr)
                total_infs[name].append(val)
                (quiet_infs[name] if t < 15 else active_infs[name]).append(val)

    print("  Controller               | Quiet(t<15) | Drift+TLS(t>=15) | Overall Mean | Actuations/run | Toggles/run")
    for name in ctrl_names:
        mq = float(np.mean(quiet_infs[name]))
        ma = float(np.mean(active_infs[name]))
        mt = float(np.mean(total_infs[name]))
        ap = acts_tot[name] / float(n_runs)
        cp = chatter_tot[name] / float(n_runs)
        print(f"  {name:24s} | {mq:10.3e}  | {ma:10.3e}       | {mt:10.3e}   | {ap:8.1f} / 60  | {cp:6.1f}")

    print(f"\nCompleted exact density-matrix benchmark in {time.perf_counter() - t0:.2f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
