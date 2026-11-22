#!/usr/bin/env python3
"""
PHY v4.0.0 — 51-Mutant Single-Fault Sabotage Gate (Gate 6).
Verifies 100.0% mutation-kill efficacy (51/51 killed, 0 compile errors) across 17 categories:
  1. hott_fiber_merkle (3)
  2. mersenne61_rref (3)
  3. noether_telescoping (3)
  4. separable_anf_guard (3)
  5. yangbaxter_soliton (3)
  6. searu_hunt_fir_lift (3)
  7. higher_order_q4 (3)
  8. yb8_noninvol_atlas (3)
  9. dual4_ame8 (3)
 10. fracton_hsf (3)
 11. biquandle8_knots (3)
 12. tetrahedron_4 (3)
 13. ame4_nogo (3)
 14. simplex4_4d (3)
 15. chiral6_cohomology (3)
 16. qscar_dtc (3)
 17. fracton_comp (3)
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import List


ROOT = Path(__file__).resolve().parent.parent
FLOW_TOOLS_DIR = (
    Path(os.environ.get("FLOW_DIR", str(ROOT.parent / "FLOW"))).resolve()
    / "sdk"
    / "tools"
)
if str(FLOW_TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(FLOW_TOOLS_DIR))

from flow_sabotage_runner import FlowMutant, run_multifile_sabotage_suite


def Mutant(
    mid: str,
    category: str,
    relpath: str,
    target: str,
    replacement: str,
    reason: str = "",
) -> FlowMutant:
    return FlowMutant(
        name=mid,
        category=category,
        target=target,
        replacement=replacement,
        relpath=relpath,
        reason=reason,
    )


MUTANTS: List[FlowMutant] = [
    # 1. hott_fiber_merkle (3)
    Mutant(
        "M01_fnv1a_prime_zero",
        "hott_fiber_merkle",
        "src/phy_ring0.c",
        "h *= UINT64_C(0x100000001B3);",
        "h *= UINT64_C(0);",
        "Collapses FNV-1a avalanche prime to 0",
    ),
    Mutant(
        "M02_rule_merkle_dropped",
        "hott_fiber_merkle",
        "src/phy_ring0.c",
        "mh = phy_fnv1a64_u64(mh, (uint64_t)(uint32_t)out->q[u]);\n    }\n    out->merkle = mh;",
        "mh = phy_fnv1a64_u64(mh, (uint64_t)(uint32_t)out->q[u]);\n    }\n    out->merkle = UINT64_C(0);",
        "Suppresses rule HoTT Merkle digest",
    ),
    Mutant(
        "M03_fiber_merkle_frozen",
        "hott_fiber_merkle",
        "src/phy_ring0.c",
        "fib.merkle = phy_fnv1a64_u64(fib.merkle, fib.state_lo ^ fib.state_hi ^ diff_mask);",
        "fib.merkle = rule->merkle;",
        "Freezes trajectory fiber Merkle at initial rule digest",
    ),
    # 2. mersenne61_rref (3)
    Mutant(
        "M04_m61_fold_shift_bug",
        "mersenne61_rref",
        "src/phy_ring0.c",
        "uint64_t mid = (uint64_t)(x >> 61U) & PHY_M61;",
        "uint64_t mid = (uint64_t)(x >> 60U) & PHY_M61;",
        "Corrupts 2^61 - 1 Mersenne folding exponent",
    ),
    Mutant(
        "M05_m61_inv_exponent",
        "mersenne61_rref",
        "src/phy_ring0.c",
        "uint64_t exp = PHY_M61 - UINT64_C(2);",
        "uint64_t exp = PHY_M61 - UINT64_C(3);",
        "Corrupts Fermat inversion exponent M61 - 2",
    ),
    Mutant(
        "M06_rref_pivot_unscaled",
        "mersenne61_rref",
        "src/phy_ring0.c",
        "uint64_t inv = phy_m61_inv(mrow[lead]);",
        "uint64_t inv = UINT64_C(1);",
        "Skips inverse pivot normalization during streaming M61 RREF",
    ),
    # 3. noether_telescoping (3)
    Mutant(
        "M07_step3_nibble_extract",
        "noether_telescoping",
        "src/phy_ring0.c",
        "return phy_pack2(phy_lo2(oddL), phy_hi2(oddR));",
        "return phy_pack2(phy_hi2(oddL), phy_lo2(oddR));",
        "Takes outer instead of inner 2-bit sites in brick-wall step3",
    ),
    Mutant(
        "M08_continuity_row_sign",
        "noether_telescoping",
        "src/phy_ring0.c",
        "row[u4 - 1U] += INT64_C(1);",
        "row[u4 - 1U] -= INT64_C(1);",
        "Flips sign of canonCurrent f0[a][b] term in continuity matrix",
    ),
    Mutant(
        "M09_canon_current_sign",
        "noether_telescoping",
        "src/phy_ring0.c",
        "int32_t val = hb + (q[f0[b][c]] - q[b]);",
        "int32_t val = hb - (q[f0[b][c]] - q[b]);",
        "Corrupts canonical Noether current j[b][c] subtraction sign",
    ),
    # 4. separable_anf_guard (3)
    Mutant(
        "M10_sep_half_subspace",
        "separable_anf_guard",
        "src/phy_ring0.c",
        "for (uint32_t s = 0U; s < 6U; ++s) {",
        "for (uint32_t s = 0U; s < 3U; ++s) {",
        "Drops right-site chiral separable subspace basis (only projects 3 of 6 dims)",
    ),
    Mutant(
        "M11_anf_deg_suppress",
        "separable_anf_guard",
        "src/phy_ring0.c",
        "if (d > max_deg) {",
        "if (d > max_deg && d < 3U) {",
        "Caps GF(2) ANF degree below 3",
    ),
    Mutant(
        "M12_factorizable_guard_bypass",
        "separable_anf_guard",
        "src/phy_ring0.c",
        "*out_r1 = b;\n                return false;",
        "*out_r1 = b;\n                return true;",
        "Returns factorizable=true early on left-output cross-dependence",
    ),
    # 5. yangbaxter_soliton (3)
    Mutant(
        "M13_yb_sense_inverted",
        "yangbaxter_soliton",
        "src/phy_ring0.c",
        "uint16_t rhs = (uint16_t)(((uint32_t)a0 << 4U) | (uint32_t)bc1);\n\n                if (lhs != rhs) {\n                    ++viol;\n                }",
        "uint16_t rhs = (uint16_t)(((uint32_t)a0 << 4U) | (uint32_t)bc1);\n\n                if (lhs == rhs) {\n                    ++viol;\n                }",
        "Inverts Yang-Baxter equality violation counter",
    ),
    Mutant(
        "M14_ring64_odd_shift",
        "yangbaxter_soliton",
        "src/phy_ring0.c",
        "uint32_t s0 = (2U * m + 1U) & 63U;",
        "uint32_t s0 = (2U * m + 2U) & 63U;",
        "Misaligns odd brick-wall layer in 64-site periodic Floquet update",
    ),
    Mutant(
        "M15_flow_jet_front_zero",
        "yangbaxter_soliton",
        "src/phy_ring0.c",
        "int64_t span_q16 = (int64_t)((uint64_t)span << 16U);",
        "int64_t span_q16 = INT64_C(0);",
        "Zeroes soliton light-cone span target fed to FLOW 2-Jet tracker",
    ),
    # 6. searu_hunt_fir_lift (3)
    Mutant(
        "M16_hunt_yb_inverted",
        "searu_hunt_fir_lift",
        "src/phy_ring0.c",
        "return (cand->flags & req_yb) == req_yb;",
        "return (cand->flags & req_yb) != req_yb;",
        "Accepts Yang-Baxter failures instead of exact solutions in phy_hunt",
    ),
    Mutant(
        "M17_fir_tt_input_flip",
        "searu_hunt_fir_lift",
        "src/phy_cert.h",
        "if (((lut[x] >> o) & 1U) != 0U) {",
        "if (((lut[x ^ 1U] >> o) & 1U) != 0U) {",
        "Flips LSB input index when building searu_mc_gf2 target truth table",
    ),
    Mutant(
        "M18_fir_concat_reverse",
        "searu_hunt_fir_lift",
        "src/phy_cert.h",
        "p->ops[out_op].ins[o] = obits[o];",
        "p->ops[out_op].ins[o] = obits[3U - o];",
        "Reverses bit order in FLOW_IR_OP_CONCAT 4-bit FIR output packing",
    ),
    # 7. higher_order_q4 (3)
    Mutant(
        "M19_q4_continuity_sign",
        "higher_order_q4",
        "src/phy_ring0.c",
        "phy_q4_add_term(ws->work_row, (fabc << 4U) | fbcd, INT64_C(1));",
        "phy_q4_add_term(ws->work_row, (fabc << 4U) | fbcd, INT64_C(-1));",
        "Flips sign of evolved pair term q4(f(a,b,c), f(b,c,d)) in 255-col Q4 RREF",
    ),
    Mutant(
        "M20_q4_triv_chiral_drop",
        "higher_order_q4",
        "src/phy_ring0.c",
        "for (uint8_t h1 = 1U; h1 < 4U; ++h1) {",
        "for (uint8_t h1 = 1U; h1 < 1U; ++h1) {",
        "Drops right-chiral pair correlations from 48-dim trivial Q4 quotient subspace",
    ),
    Mutant(
        "M21_butterfly_vel_zero",
        "higher_order_q4",
        "src/phy_ring0.c",
        "vel_sum += (int64_t)tr.vel_num_q16;",
        "vel_sum = INT64_C(0);",
        "Zeroes mean butterfly velocity v_B accumulation across the 8-site lattice defect sweep",
    ),
    # 8. yb8_noninvol_atlas (3)
    Mutant(
        "M22_yb8_check_inverted",
        "yb8_noninvol_atlas",
        "src/phy_ring0.c",
        "if (lhs8 != rhs8) {\n                    ++yb8_viol;\n                }",
        "if (lhs8 == rhs8) {\n                    ++yb8_viol;\n                }",
        "Inverts |X|=8 Yang-Baxter equality violation counter across 512 triples",
    ),
    Mutant(
        "M23_yb8_noninvol_disabled",
        "yb8_noninvol_atlas",
        "src/phy_ring0.c",
        "if (out->r[out->r[u]] != u) {\n            out->is_involutive = 0U;",
        "if (out->r[out->r[u]] != u) {\n            out->is_involutive = 1U;",
        "Suppresses R^2 != id non-involutivity witness detection in phy_analyze_yb8",
    ),
    Mutant(
        "M24_yb8_fir_concat_reverse",
        "yb8_noninvol_atlas",
        "src/phy_cert.h",
        "p->ops[yb_out_op].ins[o] = yb_obits[o];",
        "p->ops[yb_out_op].ins[o] = yb_obits[5U - o];",
        "Reverses bit order in FLOW_IR_OP_CONCAT 6-bit YB8 FIR output packing",
    ),
    # 9. dual4_ame8 (3)
    Mutant(
        "M25_dual4_bijective_only",
        "dual4_ame8",
        "src/phy_ring0.c",
        "out->is_dual_unitary = (out->is_bijective != 0U && nondeg) ? UINT8_C(1) : UINT8_C(0);",
        "out->is_dual_unitary = out->is_bijective;",
        "Treats any bijection on BitVec 4 as dual-unitary without space-time channel checks",
    ),
    Mutant(
        "M26_dual4_mix_m4_bypass",
        "dual4_ame8",
        "src/phy_ring0.c",
        "if (acc == UINT32_C(0)) {\n                return false;\n            }",
        "if (acc == UINT32_C(999999)) {\n                return false;\n            }",
        "Bypasses strict positivity check on M_+^4 and M_-^4 transfer matrices",
    ),
    Mutant(
        "M27_ame8_crossbij_forced",
        "dual4_ame8",
        "src/phy_ring0.c",
        "out->is_crossbij = cross ? UINT8_C(1) : UINT8_C(0);",
        "out->is_crossbij = UINT8_C(1);",
        "Forces 2-unitary cross-bijectivity true even on non-AME(4,8) gates",
    ),
    # 10. fracton_hsf (3)
    Mutant(
        "M28_fracton_laplacian_sign",
        "fracton_hsf",
        "src/phy_ring0.c",
        "if ((ob - b) != -(2 * phi) || (oc - c) != phi) {\n            lap_ok = false;\n        } else if (phi != 0 && !got_wit) {",
        "if ((ob - b) != -(phi) || (oc - c) != phi) {\n            lap_ok = false;\n        } else if (phi != 0 && !got_wit) {",
        "Corrupts central (-2) coefficient of 3-site discrete dipole-Laplacian law",
    ),
    Mutant(
        "M29_fracton_iso_span_relaxed",
        "fracton_hsf",
        "src/phy_ring0.c",
        "if (out->iso_span == 0U) {\n        fl |= PHY_FRACTON_FLAG_ISO_FROZEN;\n    }",
        "if (out->iso_span <= 10U) {\n        fl |= PHY_FRACTON_FLAG_ISO_FROZEN;\n    }",
        "Relaxes isolated-charge zero-spreading requirement on L=32 circuit probe",
    ),
    Mutant(
        "M30_fracton_circuit_noop",
        "fracton_hsf",
        "src/phy_ring0.c",
        "for (uint32_t layer = 0U; layer < 3U; ++layer) {",
        "for (uint32_t layer = 0U; layer < 0U; ++layer) {",
        "Short-circuits L=32 staggered 3-site fracton circuit update to no-op",
    ),
    # 11. biquandle8_knots (3)
    Mutant(
        "M31_bq8_r1_disabled",
        "biquandle8_knots",
        "src/phy_ring0.c",
        "if (out->base.r[phy_pack3(x, x)] != phy_pack3(x, x)) {",
        "if (out->base.r[phy_pack3(x, x)] > 64U) {",
        "Allows diagonal non-idempotence while still marking R1 biquandle flag",
    ),
    Mutant(
        "M32_bq8_irretract_forced",
        "biquandle8_knots",
        "src/phy_ring0.c",
        "out->is_irretract = (ret_cls == 8U) ? UINT8_C(1) : UINT8_C(0);",
        "out->is_irretract = UINT8_C(1);",
        "Forces irretractable flag regardless of congruence-class collapse",
    ),
    Mutant(
        "M33_bq8_fix_pow_clamped",
        "biquandle8_knots",
        "src/phy_ring0.c",
        "for (uint32_t s = 0U; s < k; ++s) {",
        "for (uint32_t s = 0U; s < (k ? 1U : 0U); ++s) {",
        "Clamps T(2,k) torus braid closure coloring iteration count to 1",
    ),
    # 12. tetrahedron_4 (3)
    Mutant(
        "M34_tetra2_nl_suppressed",
        "tetrahedron_4",
        "src/phy_ring0.c",
        "if (md >= 2U) {\n        *nl += 1U;",
        "if (md >= 99U) {\n        *nl += 1U;",
        "Suppresses nonlinear |X|=2 Zamolodchikov classification counter",
    ),
    Mutant(
        "M35_tetra4_check_inverted",
        "tetrahedron_4",
        "src/phy_ring0.c",
        "2U, 4U, 5U);\n        if (l != r) {\n            ++viol;\n        }",
        "2U, 4U, 5U);\n        if (l == r) {\n            ++viol;\n        }",
        "Inverts |X|=4 4096-state tetrahedron identity violation counter",
    ),
    Mutant(
        "M36_fir6_clean_flag_dropped",
        "tetrahedron_4",
        "src/phy_cert.h",
        "p->flags = FLOW_IR_FLAG_CLEAN;\n    p->num_inputs = 1U;\n    p->inputs[0] = 0U;\n    p->vals[0].width = 6U;",
        "p->flags = 0U;\n    p->num_inputs = 1U;\n    p->inputs[0] = 0U;\n    p->vals[0].width = 6U;",
        "Drops FLOW_IR_FLAG_CLEAN unitary admissibility flag from 6-bit quantum-linear .fir",
    ),
    # 13. ame4_nogo (3)
    Mutant(
        "M37_ame4_2unitary_row_bypass",
        "ame4_nogo",
        "src/phy_ring0.c",
        "if (rc_m != UINT32_C(0xF) || rd_m != UINT32_C(0xF) || cc_m != UINT32_C(0xF) ||\n            cd_m != UINT32_C(0xF)) {",
        "if (rc_m == UINT32_C(0) && rd_m == UINT32_C(0) && cc_m == UINT32_C(0) &&\n            cd_m == UINT32_C(0)) {",
        "Bypasses row/column orthogonality masks in phy_is_2unitary4",
    ),
    Mutant(
        "M38_ame4_ols_seen_relaxed",
        "ame4_nogo",
        "src/phy_ring0.c",
        "if (seen_mask != UINT32_C(0xFFFF)) {\n                continue;\n            }",
        "if (seen_mask == UINT32_C(0)) {\n                continue;\n            }",
        "Accepts non-orthogonal 4x4 Latin square pairs into OLS census",
    ),
    Mutant(
        "M39_ame4_deg_misclassified",
        "ame4_nogo",
        "src/phy_ring0.c",
        "if (deg == 1U) {\n                ++affine_ols;\n            } else {\n                ++nonlin_ols;\n            }",
        "if (deg == 2U) {\n                ++affine_ols;\n            } else {\n                ++nonlin_ols;\n            }",
        "Misclassifies GF(2)-affine OLS gates as nonlinear in AME(4,4) rigidity census",
    ),
    # 14. simplex4_4d (3)
    Mutant(
        "M40_s4_2_ccx_viol_inverted",
        "simplex4_4d",
        "src/phy_ring0.c",
        "if (l != r) {\n            if (!got && out_first_bad != NULL) {",
        "if (l == r) {\n            if (!got && out_first_bad != NULL) {",
        "Inverts |X|=2 1024-state 4-simplex equation violation check",
    ),
    Mutant(
        "M41_s4_4_check_inverted",
        "simplex4_4d",
        "src/phy_ring0.c",
        "r = phy_simplex4_apply4(s, r, 6U, 7U, 8U, 9U);\n        if (l != r) {\n            ++viol;\n        }",
        "r = phy_simplex4_apply4(s, r, 6U, 7U, 8U, 9U);\n        if (l == r) {\n            ++viol;\n        }",
        "Inverts |X|=4 1048576-state 4-simplex equation violation counter",
    ),
    Mutant(
        "M42_s4_4ch_coupled_forced",
        "simplex4_4d",
        "src/phy_ring0.c",
        "return (c0 && c1 && c2 && c3) ? UINT8_C(1) : UINT8_C(0);",
        "return (c0 && c1 && c2 && c3) ? UINT8_C(0) : UINT8_C(1);",
        "Inverts 4-channel coupling detector in phy_check_4ch_coupled4",
    ),
    # 15. chiral6_cohomology (3)
    Mutant(
        "M43_q6_cocycle_viol_inverted",
        "chiral6_cohomology",
        "src/phy_ring0.c",
        "UINT32_C(3);\n                if (lhs != rhs) {\n                    ++viol;\n                }",
        "UINT32_C(3);\n                if (lhs == rhs) {\n                    ++viol;\n                }",
        "Inverts mod-4 2-cocycle equality check over octahedral quandle Q_6",
    ),
    Mutant(
        "M44_q6_neg_braid_sign_flipped",
        "chiral6_cohomology",
        "src/phy_ring0.c",
        "s = (s + (UINT32_C(4) - ((uint32_t)phi[y] & UINT32_C(3)))) & UINT32_C(3);",
        "s = (s + ((uint32_t)phi[y] & UINT32_C(3))) & UINT32_C(3);",
        "Flips crossing sign of negative braid T(2,-3) state-sum weight",
    ),
    Mutant(
        "M45_q6_noncobound_forced",
        "chiral6_cohomology",
        "src/phy_ring0.c",
        "out->is_non_coboundary = (out->cocycle_viol == 0U && non_cob) ? UINT8_C(1) : UINT8_C(0);",
        "out->is_non_coboundary = (out->cocycle_viol == 0U || non_cob) ? UINT8_C(1) : UINT8_C(0);",
        "Relaxes non-coboundary guard on trivial zero cocycle via OR instead of AND",
    ),
    # 16. qscar_dtc (3)
    Mutant(
        "M46_qscar_dual_unitary_forced",
        "qscar_dtc",
        "src/phy_ring0.c",
        "if (out_dual != NULL) {\n        *out_dual = dual ? UINT8_C(1) : UINT8_C(0);\n    }",
        "if (out_dual != NULL) {\n        *out_dual = dual ? UINT8_C(1) : UINT8_C(1);\n    }",
        "Forces dual-unitary reshuffle check true even on identity matrix",
    ),
    Mutant(
        "M47_qscar_ctrl_inverted",
        "qscar_dtc",
        "src/phy_ring0.c",
        "bool keep = inv_ctrl ? (b0 == b3) : (b0 != b3);",
        "bool keep = inv_ctrl ? (b0 != b3) : (b0 == b3);",
        "Inverts projective scar projector kinetic constraint on Floquet circuit",
    ),
    Mutant(
        "M48_qscar_purity_conj_flipped",
        "qscar_dtc",
        "src/phy_ring0.c",
        "mi += im[idx1] * re[idx2] - re[idx1] * im[idx2];",
        "mi += im[idx1] * re[idx2] + re[idx1] * im[idx2];",
        "Corrupts complex conjugation sign in Z[i] reduced density matrix purity",
    ),
    # 17. fracton_comp (3)
    Mutant(
        "M49_fcomp_toff_swap_dropped",
        "fracton_comp",
        "src/phy_ring0.c",
        "t[25] = 34U;\n    t[34] = 25U;",
        "t[25] = 25U;\n    t[34] = 34U;",
        "Drops dipole-carrier Toffoli transposition between (1,2,1) and (2,0,2)",
    ),
    Mutant(
        "M50_fcomp_ccx_check_relaxed",
        "fracton_comp",
        "src/phy_ring0.c",
        "uint32_t exp_tgt = tgt ^ (ca & cc);",
        "uint32_t exp_tgt = tgt ^ (ca | cc);",
        "Replaces Toffoli AND control logic with OR in dipole-Toffoli verifier",
    ),
    Mutant(
        "M51_fcomp_hsf_strict_inverted",
        "fracton_comp",
        "src/phy_ring0.c",
        "out->scaling_strict = strict ? UINT8_C(1) : UINT8_C(0);",
        "out->scaling_strict = strict ? UINT8_C(0) : UINT8_C(1);",
        "Inverts L=6,7,8 strict HSF Krylov-density collapse flag",
    ),
]


def main() -> int:
    assert len(MUTANTS) == 51, f"Expected 51 mutants, got {len(MUTANTS)}"
    cat_counts: dict[str, int] = {}
    for m in MUTANTS:
        cat_counts[m.category] = cat_counts.get(m.category, 0) + 1
    assert len(cat_counts) == 17 and all(v == 3 for v in cat_counts.values()), cat_counts

    return run_multifile_sabotage_suite(
        repo_root=ROOT,
        ring0_rel="src/phy_ring0.c",
        cli_rel="src/phy_cli.c",
        extra_header_rels=("src/phy_cert.h",),
        mutants=MUTANTS,
        prefix="phy_sab_",
    )


if __name__ == "__main__":
    raise SystemExit(main())
