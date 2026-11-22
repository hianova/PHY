#ifndef PHY_CERT_H
#define PHY_CERT_H

#include "flow_ir_builder.h"
#include "phy.h"
#include "searu_mc_gf2.h"
#include <stdio.h>
#include <string.h>

typedef struct {
    searu_mc_ws mc_ws;
    searu_mc_plan mc_plan;
    flow_ir_prog_t ir_prog;
    flow_ir_prog_t ir_reloaded;
    uint8_t wire_buf[FLOW_IR_MAX_FILE_BYTES];
    size_t wire_len;
    uint16_t and_count;
    uint16_t mult_depth;
    bool ir_linear_ok;
    bool ir_unitary_ok;
    bool ir_eval_ok;
} phy_fir_cert_ws_t;

static inline uint32_t phy_cert_emit_const(flow_ir_prog_t *p, uint64_t val, const char *name) {
    return flow_ir_emit_const(p, 1U, val, name);
}

static inline uint32_t phy_cert_emit_slice(flow_ir_prog_t *p, uint32_t src_val, uint32_t bit,
                                           const char *name) {
    return flow_ir_emit_op1(p, FLOW_IR_OP_SLICE, 1U, src_val, (uint64_t)bit, name);
}

static inline uint32_t phy_cert_emit_binop(flow_ir_prog_t *p, uint16_t kind, uint32_t lhs,
                                           uint32_t rhs, const char *name) {
    return flow_ir_emit_op2(p, kind, 1U, lhs, rhs, name);
}

static inline uint32_t phy_cert_emit_xor_mask(flow_ir_prog_t *p, const uint32_t *basis_vals,
                                              uint32_t n_basis, uint32_t mask, uint32_t const0_id,
                                              const char *prefix) {
    uint32_t acc = UINT32_MAX;
    uint32_t step = 0U;
    for (uint32_t b = 0U; b < n_basis; ++b) {
        if (((mask >> b) & UINT32_C(1)) != UINT32_C(0)) {
            if (basis_vals[b] == UINT32_MAX) {
                return UINT32_MAX;
            }
            if (acc == UINT32_MAX) {
                acc = basis_vals[b];
            } else {
                char nm[16];
                snprintf(nm, sizeof(nm), "%s%u", prefix, step);
                acc = phy_cert_emit_binop(p, FLOW_IR_OP_XOR_VEC, acc, basis_vals[b], nm);
                if (acc == UINT32_MAX) {
                    return UINT32_MAX;
                }
                ++step;
            }
        }
    }
    return (acc == UINT32_MAX) ? const0_id : acc;
}

static inline void phy_cert_build_anf_mc_plan(const uint8_t p[16], searu_mc_plan *plan) {
    memset(plan, 0, sizeof(*plan));
    uint8_t anf[4][16];
    uint8_t need[16];
    memset(need, 0, sizeof(need));
    for (uint32_t o = 0U; o < 4U; ++o) {
        for (uint32_t x = 0U; x < 16U; ++x) {
            anf[o][x] = (uint8_t)((p[x] >> o) & 1U);
        }
        for (uint32_t v = 0U; v < 4U; ++v) {
            for (uint32_t mask = 0U; mask < 16U; ++mask) {
                if ((mask & (UINT32_C(1) << v)) != 0U) {
                    anf[o][mask] ^= anf[o][mask ^ (UINT32_C(1) << v)];
                }
            }
        }
        for (uint32_t mask = 0U; mask < 16U; ++mask) {
            if (anf[o][mask] != 0U) {
                need[mask] = 1U;
            }
        }
    }
    for (uint32_t d = 4U; d >= 3U; --d) {
        for (uint32_t mask = 0U; mask < 16U; ++mask) {
            if (need[mask] != 0U && searu_popcnt64((uint64_t)mask) == d) {
                uint32_t lsb = mask & (~mask + UINT32_C(1));
                uint32_t sub = mask ^ lsb;
                need[sub] = 1U;
            }
        }
    }
    uint32_t mon_basis[16];
    uint32_t depths[SEARU_MC_MAX_BASIS];
    memset(depths, 0, sizeof(depths));
    mon_basis[0] = 0U;
    for (uint32_t i = 0U; i < 4U; ++i) {
        mon_basis[UINT32_C(1) << i] = 1U + i;
    }
    uint32_t k = 0U;
    uint32_t max_d = 0U;
    for (uint32_t d = 2U; d <= 4U; ++d) {
        for (uint32_t mask = 0U; mask < 16U; ++mask) {
            if (need[mask] != 0U && searu_popcnt64((uint64_t)mask) == d && k < SEARU_MC_MAX_AND) {
                uint32_t lsb = mask & (~mask + UINT32_C(1));
                uint32_t sub = mask ^ lsb;
                uint32_t lb = mon_basis[sub];
                uint32_t rb = mon_basis[lsb];
                plan->and_lhs_mask[k] = (UINT32_C(1) << lb);
                plan->and_rhs_mask[k] = (UINT32_C(1) << rb);
                uint32_t gid = 5U + k;
                mon_basis[mask] = gid;
                depths[gid] = 1U + ((depths[lb] > depths[rb]) ? depths[lb] : depths[rb]);
                if (depths[gid] > max_d) {
                    max_d = depths[gid];
                }
                ++k;
            }
        }
    }
    plan->num_and = k;
    plan->mult_depth = max_d;
    plan->mode_used = 3U;
    plan->exact = true;
    for (uint32_t o = 0U; o < 4U; ++o) {
        uint32_t om = 0U;
        for (uint32_t mask = 0U; mask < 16U; ++mask) {
            if (anf[o][mask] != 0U) {
                om |= (UINT32_C(1) << mon_basis[mask]);
            }
        }
        plan->out_mask[o] = om;
    }
}

static inline bool phy_synth_fir_4bit(const uint8_t lut[16], uint64_t seed, const char *mod_name,
                                      const char *out_fir_path, phy_fir_cert_ws_t *ws) {
    if (lut == NULL || ws == NULL) {
        return false;
    }
    ws->and_count = 0U;
    ws->mult_depth = 0U;
    ws->wire_len = 0U;
    ws->ir_linear_ok = false;
    ws->ir_unitary_ok = false;
    ws->ir_eval_ok = false;

    searu_mc_spec spec;
    memset(&spec, 0, sizeof(spec));
    spec.n_in = 4U;
    spec.n_out = 4U;
    spec.max_and = 8U;
    spec.max_depth = 0U;
    spec.n_in_port0 = 0U;
    spec.budget_scale = 1U;
    spec.seed = seed;

    for (uint32_t o = 0U; o < 4U; ++o) {
        uint64_t tt = UINT64_C(0);
        for (uint32_t x = 0U; x < 16U; ++x) {
            if (((lut[x] >> o) & 1U) != 0U) {
                tt |= (UINT64_C(1) << x);
            }
        }
        spec.target_tt[o][0] = tt;
    }

    if (!searu_solve_mc_gf2(&ws->mc_plan, &spec, &ws->mc_ws)) {
        phy_cert_build_anf_mc_plan(lut, &ws->mc_plan);
    }

    ws->and_count = (uint16_t)ws->mc_plan.num_and;
    ws->mult_depth = (uint16_t)ws->mc_plan.mult_depth;

    flow_ir_prog_t *p = &ws->ir_prog;
    memset(p, 0, sizeof(*p));
    snprintf(p->name, sizeof(p->name), "%s", (mod_name != NULL) ? mod_name : "phy_rule_2x2");
    p->mult = FLOW_IR_MULT_LINEAR;
    p->flags = FLOW_IR_FLAG_CLEAN;
    p->num_inputs = 1U;
    p->inputs[0] = 0U;
    p->vals[0].width = 4U;
    p->vals[0].is_classical = 0U;
    snprintf(p->vals[0].name, sizeof(p->vals[0].name), "u");
    p->num_vals = 1U;
    p->num_ops = 0U;

    uint32_t basis[24];
    uint32_t const0_id = phy_cert_emit_const(p, UINT64_C(0), "c0");
    basis[0] = phy_cert_emit_const(p, UINT64_C(1), "c1");
    for (uint32_t i = 0U; i < 4U; ++i) {
        char nm[16];
        snprintf(nm, sizeof(nm), "x%u", i);
        basis[1U + i] = phy_cert_emit_slice(p, 0U, i, nm);
    }

    for (uint32_t k = 0U; k < ws->mc_plan.num_and; ++k) {
        uint32_t n_cur = 5U + k;
        char lnm[16];
        char rnm[16];
        char anm[16];
        snprintf(lnm, sizeof(lnm), "l%u_", k);
        snprintf(rnm, sizeof(rnm), "r%u_", k);
        snprintf(anm, sizeof(anm), "a%u", k);
        uint32_t lhs =
            phy_cert_emit_xor_mask(p, basis, n_cur, ws->mc_plan.and_lhs_mask[k], const0_id, lnm);
        uint32_t rhs =
            phy_cert_emit_xor_mask(p, basis, n_cur, ws->mc_plan.and_rhs_mask[k], const0_id, rnm);
        basis[5U + k] = phy_cert_emit_binop(p, FLOW_IR_OP_AND_VEC, lhs, rhs, anm);
    }

    uint32_t n_all = 5U + ws->mc_plan.num_and;
    uint32_t obits[4];
    for (uint32_t o = 0U; o < 4U; ++o) {
        char pfx[16];
        snprintf(pfx, sizeof(pfx), "o%u_", o);
        obits[o] = phy_cert_emit_xor_mask(p, basis, n_all, ws->mc_plan.out_mask[o], const0_id, pfx);
        if (obits[o] == UINT32_MAX) {
            return false;
        }
    }

    if (p->num_vals >= FLOW_IR_MAX_VALS || p->num_ops >= FLOW_IR_MAX_OPS) {
        return false;
    }
    uint32_t out_vid = p->num_vals;
    uint32_t out_op = p->num_ops;
    p->vals[out_vid].width = 4U;
    p->vals[out_vid].is_classical = 0U;
    snprintf(p->vals[out_vid].name, sizeof(p->vals[out_vid].name), "out");
    p->ops[out_op].kind = FLOW_IR_OP_CONCAT;
    p->ops[out_op].num_ins = 4U;
    p->ops[out_op].out_val = out_vid;
    for (uint32_t o = 0U; o < 4U; ++o) {
        p->ops[out_op].ins[o] = obits[o];
    }
    p->num_vals++;
    p->num_ops++;
    p->num_outputs = 1U;
    p->outputs[0] = out_vid;

    if (!flow_ir_encode_bytes(p, ws->wire_buf, sizeof(ws->wire_buf), &ws->wire_len)) {
        return false;
    }

    if (flow_ir_load_bytes_fiber(&ws->ir_reloaded, ws->wire_buf, ws->wire_len) !=
        FLOW_IR_FIBER_COMMITTED) {
        return false;
    }
    ws->ir_linear_ok = (ws->ir_reloaded.mult == FLOW_IR_MULT_LINEAR);
    ws->ir_unitary_ok = flow_ir_admits_unitary(&ws->ir_reloaded);

    bool match = true;
    for (uint32_t u = 0U; u < 16U; ++u) {
        uint64_t in_arr[1] = {(uint64_t)u};
        uint64_t out_arr[1] = {UINT64_C(0)};
        if (!flow_ir_eval_u64(&ws->ir_reloaded, in_arr, out_arr)) {
            match = false;
            break;
        }
        if ((out_arr[0] & UINT64_C(0xF)) != (uint64_t)lut[u]) {
            match = false;
            break;
        }
    }
    ws->ir_eval_ok = match;

    if (out_fir_path != NULL && out_fir_path[0] != '\0') {
        FILE *fp = fopen(out_fir_path, "wb");
        if (fp == NULL) {
            return false;
        }
        size_t nw = fwrite(ws->wire_buf, 1U, ws->wire_len, fp);
        fclose(fp);
        if (nw != ws->wire_len) {
            return false;
        }
    }

    return ws->ir_linear_ok && ws->ir_unitary_ok && ws->ir_eval_ok;
}

static inline bool phy_synth_fir_rule(const phy_rule_analysis_t *rule, const char *out_fir_path,
                                      phy_fir_cert_ws_t *ws) {
    if (rule == NULL) {
        return false;
    }
    return phy_synth_fir_4bit(rule->p, rule->lut_packed ^ UINT64_C(0x50485931), "phy_rule_2x2",
                              out_fir_path, ws);
}

static inline void phy_emit_bv4_lut_lean(FILE *fp, const char *fn_name, const uint8_t lut[16]) {
    fprintf(fp, "def %s (u : BitVec 4) : BitVec 4 :=\n", fn_name);
    for (uint32_t i = 0U; i < 15U; ++i) {
        fprintf(fp, "  if u == %u#4 then %u#4 else\n", i, (unsigned)lut[i]);
    }
    fprintf(fp, "  %u#4\n\n", (unsigned)lut[15]);
}

static inline void phy_emit_bv32_lut_lean(FILE *fp, const char *fn_name, const int32_t q[16]) {
    fprintf(fp, "def %s (u : BitVec 4) : BitVec 32 :=\n", fn_name);
    for (uint32_t i = 0U; i < 16U; ++i) {
        const char *pfx = (i < 15U) ? "  if u == %u#4 then " : "  ";
        if (i < 15U) {
            fprintf(fp, pfx, i);
        } else {
            fprintf(fp, "  ");
        }
        if (q[i] < 0) {
            fprintf(fp, "(- (%ld#32))", (long)(-(int64_t)q[i]));
        } else {
            fprintf(fp, "%ld#32", (long)q[i]);
        }
        if (i < 15U) {
            fprintf(fp, " else\n");
        } else {
            fprintf(fp, "\n\n");
        }
    }
}

static inline bool phy_emit_lean_noether_cert(const phy_rule_analysis_t *rule, const char *path) {
    bool is_joint = (rule != NULL) && ((rule->flags & PHY_RULE_FLAG_JOINT) != 0U);
    if (rule == NULL || path == NULL || rule->is_bijective == 0U || rule->rank_ent == 0U ||
        rule->continuity_fail != 0 || (!is_joint && rule->j_nonzero == 0U) || rule->anf_deg < 2U ||
        rule->is_factorizable != 0U) {
        return false;
    }
    uint8_t jb = 0U;
    uint8_t jc = 0U;
    bool found_j = false;
    for (uint8_t b = 0U; b < 16U && !found_j; ++b) {
        for (uint8_t c = 0U; c < 16U && !found_j; ++c) {
            if (rule->j[b][c] != 0) {
                jb = b;
                jc = c;
                found_j = true;
            }
        }
    }
    if (!is_joint && !found_j) {
        return false;
    }

    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }

    const char *ns = (rule->yb_violations == 0U) ? "Phy.DiscoveredJoint" : "Phy.DiscoveredNoether";
    fprintf(fp, "import Phy.Noether\n\n");
    fprintf(fp, "namespace %s\n\n", ns);
    fprintf(fp, "open Phy.Noether\n\n");
    fprintf(
        fp,
        "-- Merkle: 0x%016llx, ANF deg=%u, yb_viol=%u, rank_tot=%u, rank_ent=%u, j_nz=%u/256\n\n",
        (unsigned long long)rule->merkle, (unsigned)rule->anf_deg, (unsigned)rule->yb_violations,
        (unsigned)rule->rank_tot, (unsigned)rule->rank_ent, (unsigned)rule->j_nonzero);

    phy_emit_bv4_lut_lean(fp, "ruleP", rule->p);
    phy_emit_bv4_lut_lean(fp, "rulePInv", rule->pinv);
    phy_emit_bv32_lut_lean(fp, "chargeQ32", rule->q);

    fprintf(fp, "theorem ruleP_inv_left (u : BitVec 4) : rulePInv (ruleP u) = u := by\n");
    fprintf(fp, "  simp only [ruleP, rulePInv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleP_inv_right (u : BitVec 4) : ruleP (rulePInv u) = u := by\n");
    fprintf(fp, "  simp only [ruleP, rulePInv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleP_bijective : Bijective ruleP := by\n");
    fprintf(fp, "  refine ⟨fun a b h => ?_, fun b => ⟨rulePInv b, ruleP_inv_right b⟩⟩\n");
    fprintf(fp, "  have h1 := congrArg rulePInv h\n");
    fprintf(fp, "  simpa [ruleP_inv_left] using h1\n\n");

    fprintf(fp, "theorem chargeQ32_bounded (u : BitVec 4) : BoundedBV32 100000000 (chargeQ32 u) := "
                "by\n");
    fprintf(fp, "  unfold BoundedBV32\n  revert u\n  decide\n\n");

    fprintf(fp, "theorem continuity_bv32 (a b c : BitVec 4) :\n");
    fprintf(fp, "    chargeQ32 (step3 ruleP a b c) - chargeQ32 b =\n");
    fprintf(fp, "      canonCurrent32 ruleP chargeQ32 b c - canonCurrent32 ruleP chargeQ32 a b := "
                "by\n");
    fprintf(fp, "  simp only [step3, canonCurrent32, pack2, hi2, lo2, ruleP, chargeQ32]\n");
    fprintf(fp, "  bv_decide\n\n");

    fprintf(fp, "theorem global_conservation {N : Nat} (hN : 0 < N) (u : Fin N → BitVec 4) :\n");
    fprintf(fp, "    ringSum N (fun i => (chargeQ32 (floquetRing hN ruleP u i)).toInt) =\n");
    fprintf(fp, "      ringSum N (fun i => (chargeQ32 (u i)).toInt) :=\n");
    fprintf(fp, "  forall_ring_conservation_of_bv32 hN ruleP chargeQ32 chargeQ32_bounded "
                "continuity_bv32 u\n\n");

    fprintf(fp, "theorem charge_nonseparable : ¬ IsSeparable (fun u => (chargeQ32 u).toInt) :=\n");
    fprintf(fp, "  nonseparable_of_cross_diff %u#2 %u#2 %u#2 %u#2 (by decide)\n\n",
            (unsigned)rule->sep_l0, (unsigned)rule->sep_l1, (unsigned)rule->sep_r0,
            (unsigned)rule->sep_r1);

    if (found_j) {
        fprintf(fp,
                "theorem current_active : ∃ x y, canonCurrent ruleP (fun u => (chargeQ32 u).toInt) "
                "x y ≠ 0 :=\n");
        fprintf(fp,
                "  active_current_witness ruleP (fun u => (chargeQ32 u).toInt) %u#4 %u#4 (by "
                "decide)\n\n",
                (unsigned)jb, (unsigned)jc);
    } else {
        fprintf(fp, "theorem current_topological : ∀ b c : BitVec 4, canonCurrent ruleP (fun u => "
                    "(chargeQ32 u).toInt) b c = 0 := by\n");
        fprintf(fp, "  decide\n\n");
    }

    if (rule->yb_violations == 0U) {
        fprintf(fp, "theorem ruleP_satisfies_ybe : IsSetTheoreticYB ruleP := by\n");
        fprintf(fp, "  intro t\n");
        fprintf(fp, "  simp only [ybApply12, ybApply23, pack2, ruleP]\n");
        fprintf(fp, "  bv_decide\n\n");
    }

    fprintf(fp,
            "theorem ruleP_nonlinear : ¬ IsAffineGF2 (fun x => BitVec.getLsbD (ruleP x) %u) :=\n",
            (unsigned)rule->nl_bit);
    fprintf(fp, "  gf2_second_diff_witnesses_nonlinear %u#4 %u#4 (by decide)\n\n",
            (unsigned)rule->nl_a, (unsigned)rule->nl_b);

    if (rule->nf_axis == 0U) {
        fprintf(fp, "theorem ruleP_nonfactorizable : ¬ IsFactorizable ruleP :=\n");
        fprintf(fp, "  non_factorizable_of_hi_cross %u#2 %u#2 %u#2 (by decide)\n\n",
                (unsigned)rule->nf_l, (unsigned)rule->nf_r0, (unsigned)rule->nf_r1);
    } else {
        fprintf(fp, "theorem ruleP_nonfactorizable : ¬ IsFactorizable ruleP :=\n");
        fprintf(fp, "  non_factorizable_of_lo_cross %u#2 %u#2 %u#2 (by decide)\n\n",
                (unsigned)rule->nf_r0, (unsigned)rule->nf_r1, (unsigned)rule->nf_l);
    }

    fprintf(fp, "end %s\n", ns);
    fclose(fp);
    return true;
}

static inline bool phy_emit_lean_yb_cert(const phy_rule_analysis_t *rule, const char *path) {
    if (rule == NULL || path == NULL || rule->is_bijective == 0U || rule->yb_violations != 0U ||
        rule->anf_deg < 2U || rule->is_factorizable != 0U || rule->is_swap != 0U) {
        return false;
    }

    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }

    fprintf(fp, "import Phy.Noether\n\n");
    fprintf(fp, "namespace Phy.DiscoveredYB\n\n");
    fprintf(fp, "open Phy.Noether\n\n");
    fprintf(fp, "-- Merkle: 0x%016llx, ANF deg=%u, yb_violations=%u, rank_ent=%u, j_nz=%u\n\n",
            (unsigned long long)rule->merkle, (unsigned)rule->anf_deg,
            (unsigned)rule->yb_violations, (unsigned)rule->rank_ent, (unsigned)rule->j_nonzero);

    phy_emit_bv4_lut_lean(fp, "ruleYB", rule->p);
    phy_emit_bv4_lut_lean(fp, "ruleYBInv", rule->pinv);

    fprintf(fp, "theorem ruleYB_inv_left (u : BitVec 4) : ruleYBInv (ruleYB u) = u := by\n");
    fprintf(fp, "  simp only [ruleYB, ruleYBInv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleYB_inv_right (u : BitVec 4) : ruleYB (ruleYBInv u) = u := by\n");
    fprintf(fp, "  simp only [ruleYB, ruleYBInv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleYB_bijective : Bijective ruleYB := by\n");
    fprintf(fp, "  refine ⟨fun a b h => ?_, fun b => ⟨ruleYBInv b, ruleYB_inv_right b⟩⟩\n");
    fprintf(fp, "  have h1 := congrArg ruleYBInv h\n");
    fprintf(fp, "  simpa [ruleYB_inv_left] using h1\n\n");

    fprintf(fp, "theorem ruleYB_satisfies_ybe : IsSetTheoreticYB ruleYB := by\n");
    fprintf(fp, "  intro t\n");
    fprintf(fp, "  simp only [ybApply12, ybApply23, pack2, ruleYB]\n");
    fprintf(fp, "  bv_decide\n\n");

    fprintf(fp,
            "theorem ruleYB_nonlinear : ¬ IsAffineGF2 (fun x => BitVec.getLsbD (ruleYB x) %u) "
            ":=\n",
            (unsigned)rule->nl_bit);
    fprintf(fp, "  gf2_second_diff_witnesses_nonlinear %u#4 %u#4 (by decide)\n\n",
            (unsigned)rule->nl_a, (unsigned)rule->nl_b);

    if (rule->nf_axis == 0U) {
        fprintf(fp, "theorem ruleYB_nonfactorizable : ¬ IsFactorizable ruleYB :=\n");
        fprintf(fp, "  non_factorizable_of_hi_cross %u#2 %u#2 %u#2 (by decide)\n\n",
                (unsigned)rule->nf_l, (unsigned)rule->nf_r0, (unsigned)rule->nf_r1);
    } else {
        fprintf(fp, "theorem ruleYB_nonfactorizable : ¬ IsFactorizable ruleYB :=\n");
        fprintf(fp, "  non_factorizable_of_lo_cross %u#2 %u#2 %u#2 (by decide)\n\n",
                (unsigned)rule->nf_r0, (unsigned)rule->nf_r1, (unsigned)rule->nf_l);
    }

    fprintf(fp, "theorem ruleYB_ne_swap : ruleYB ≠ swap2 := by\n");
    fprintf(fp, "  intro h\n");
    fprintf(fp, "  have h1 := congrFun h %u#4\n", (unsigned)rule->yb_neq_swap_u);
    fprintf(fp, "  exact (by decide : ruleYB %u#4 ≠ swap2 %u#4) h1\n\n",
            (unsigned)rule->yb_neq_swap_u, (unsigned)rule->yb_neq_swap_u);

    fprintf(fp, "end Phy.DiscoveredYB\n");
    fclose(fp);
    return true;
}

static inline bool phy_emit_rule_json(const phy_rule_analysis_t *rule, const phy_sim_trace_t *sim,
                                      const phy_fir_cert_ws_t *fir, const char *path) {
    if (rule == NULL || path == NULL) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "{\n");
    fprintf(fp, "  \"p\": [");
    for (uint32_t i = 0U; i < 16U; ++i) {
        fprintf(fp, "%s%u", (i > 0U) ? ", " : "", (unsigned)rule->p[i]);
    }
    fprintf(fp, "],\n");
    fprintf(fp, "  \"q\": [");
    for (uint32_t i = 0U; i < 16U; ++i) {
        fprintf(fp, "%s%ld", (i > 0U) ? ", " : "", (long)rule->q[i]);
    }
    fprintf(fp, "],\n");
    fprintf(fp, "  \"lut_packed\": \"0x%016llx\",\n", (unsigned long long)rule->lut_packed);
    fprintf(fp, "  \"merkle\": \"0x%016llx\",\n", (unsigned long long)rule->merkle);
    fprintf(fp, "  \"flags\": %u,\n", (unsigned)rule->flags);
    fprintf(fp, "  \"anf_deg\": %u,\n", (unsigned)rule->anf_deg);
    fprintf(fp, "  \"yb_violations\": %u,\n", (unsigned)rule->yb_violations);
    fprintf(fp, "  \"rank_tot\": %u,\n", (unsigned)rule->rank_tot);
    fprintf(fp, "  \"rank_ent\": %u,\n", (unsigned)rule->rank_ent);
    fprintf(fp, "  \"j_nonzero\": %u,\n", (unsigned)rule->j_nonzero);
    fprintf(fp, "  \"cycle_lcm\": %u,\n", (unsigned)rule->cycle_lcm);
    fprintf(fp, "  \"is_bijective\": %u,\n", (unsigned)rule->is_bijective);
    fprintf(fp, "  \"is_factorizable\": %u,\n", (unsigned)rule->is_factorizable);
    fprintf(fp, "  \"continuity_fail\": %d", rule->continuity_fail);
    if (sim != NULL) {
        fprintf(fp, ",\n  \"sim_steps\": %u,\n", (unsigned)sim->steps_run);
        fprintf(fp, "  \"sim_max_span\": %u,\n", (unsigned)sim->max_span);
        fprintf(fp, "  \"sim_vel_q16\": %d,\n", sim->vel_num_q16);
        fprintf(fp, "  \"sim_max_drift\": %lld", (long long)sim->max_charge_drift);
    }
    if (fir != NULL) {
        fprintf(fp, ",\n  \"fir_and_count\": %u,\n", (unsigned)fir->and_count);
        fprintf(fp, "  \"fir_mult_depth\": %u,\n", (unsigned)fir->mult_depth);
        fprintf(fp, "  \"fir_wire_bytes\": %zu,\n", fir->wire_len);
        fprintf(fp, "  \"fir_eval_ok\": %s", fir->ir_eval_ok ? "true" : "false");
    }
    fprintf(fp, "\n}\n");
    fclose(fp);
    return true;
}

/* --- PHY v2.0: Pruned Balanced Bit-Tree Lean 4 & 6-Bit YBE FLOW_IR Certificates --- */

static inline void phy_emit_bv_lit(FILE *fp, int32_t v, uint32_t out_width) {
    if (v < 0) {
        fprintf(fp, "(- (%ld#%u))", (long)(-(int64_t)v), out_width);
    } else {
        fprintf(fp, "%ld#%u", (long)v, out_width);
    }
}

static inline void phy_emit_bv_tree_rec(FILE *fp, const char *var, uint32_t key_width,
                                        uint32_t out_width, uint32_t lo, uint32_t hi,
                                        const int32_t *vals, uint32_t indent) {
    int32_t first = vals[lo];
    bool all_same = true;
    for (uint32_t i = lo + 1U; i < hi; ++i) {
        if (vals[i] != first) {
            all_same = false;
            break;
        }
    }
    if (all_same || lo + 1U >= hi) {
        phy_emit_bv_lit(fp, first, out_width);
        return;
    }
    uint32_t mid = lo + ((hi - lo) >> 1U);
    fprintf(fp, "(if %s < %u#%u then ", var, mid, key_width);
    if (hi - lo > 8U) {
        fprintf(fp, "\n");
        for (uint32_t s = 0U; s < indent + 2U; ++s) {
            fputc(' ', fp);
        }
    }
    phy_emit_bv_tree_rec(fp, var, key_width, out_width, lo, mid, vals, indent + 2U);
    if (hi - lo > 8U) {
        fprintf(fp, "\n");
        for (uint32_t s = 0U; s < indent; ++s) {
            fputc(' ', fp);
        }
        fprintf(fp, "else\n");
        for (uint32_t s = 0U; s < indent + 2U; ++s) {
            fputc(' ', fp);
        }
    } else {
        fprintf(fp, " else ");
    }
    phy_emit_bv_tree_rec(fp, var, key_width, out_width, mid, hi, vals, indent + 2U);
    fputc(')', fp);
}

static inline void phy_emit_bv6_lut_tree(FILE *fp, const char *fn_name, const uint8_t lut[64]) {
    int32_t tmp[64];
    for (uint32_t i = 0U; i < 64U; ++i) {
        tmp[i] = (int32_t)lut[i];
    }
    fprintf(fp, "def %s (u : BitVec 6) : BitVec 6 :=\n  ", fn_name);
    phy_emit_bv_tree_rec(fp, "u", 6U, 6U, 0U, 64U, tmp, 2U);
    fprintf(fp, "\n\n");
}

static inline void phy_emit_bv3_pair_tree(FILE *fp, const char *fn_name, const char *p1,
                                          const char *p2, const uint8_t tbl[8][8]) {
    int32_t tmp[64];
    for (uint32_t a = 0U; a < 8U; ++a) {
        for (uint32_t b = 0U; b < 8U; ++b) {
            tmp[(a << 3U) | b] = (int32_t)tbl[a][b];
        }
    }
    fprintf(fp, "def %s (%s %s : BitVec 3) : BitVec 3 :=\n", fn_name, p1, p2);
    fprintf(fp, "  let key : BitVec 6 := %s ++ %s\n  ", p1, p2);
    phy_emit_bv_tree_rec(fp, "key", 6U, 3U, 0U, 64U, tmp, 2U);
    fprintf(fp, "\n\n");
}

static inline bool phy_emit_lean_q4_cert(const phy_rule_analysis_t *rule,
                                         const phy_q4_analysis_t *q4, const char *path) {
    if (rule == NULL || q4 == NULL || path == NULL || rule->is_bijective == 0U ||
        q4->genuine_q4 == 0U || q4->continuity4_fail != 0 || q4->q4_nnz == 0U) {
        return false;
    }

    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }

    fprintf(fp, "import Phy.Noether\n\n");
    fprintf(fp, "namespace Phy.DiscoveredQ4\n\n");
    fprintf(fp, "open Phy.Noether\n\n");
    fprintf(
        fp,
        "-- LUT: 0x%016llx, Q4 Merkle: 0x%016llx, tot2=%u, tot4=%u, genuine_q4=%u, j4_nz=%u\n\n",
        (unsigned long long)q4->lut_packed, (unsigned long long)q4->q4_merkle, (unsigned)q4->tot2,
        (unsigned)q4->tot4, (unsigned)q4->genuine_q4, (unsigned)q4->j4_nonzero);

    phy_emit_bv4_lut_lean(fp, "ruleP", rule->p);
    phy_emit_bv4_lut_lean(fp, "rulePInv", rule->pinv);

    fprintf(fp, "def q4_bv (u v : BitVec 4) : BitVec 32 :=\n");
    fprintf(fp, "  let uv : BitVec 8 := u ++ v\n  ");
    phy_emit_bv_tree_rec(fp, "uv", 8U, 32U, 0U, 256U, q4->q4, 2U);
    fprintf(fp, "\n\n");

    fprintf(fp, "theorem ruleP_inv_left (u : BitVec 4) : rulePInv (ruleP u) = u := by\n");
    fprintf(fp, "  simp only [ruleP, rulePInv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleP_inv_right (u : BitVec 4) : ruleP (rulePInv u) = u := by\n");
    fprintf(fp, "  simp only [ruleP, rulePInv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleP_bijective : Bijective ruleP := by\n");
    fprintf(fp, "  refine ⟨fun a b h => ?_, fun b => ⟨rulePInv b, ruleP_inv_right b⟩⟩\n");
    fprintf(fp, "  have h1 := congrArg rulePInv h\n");
    fprintf(fp, "  simpa [ruleP_inv_left] using h1\n\n");

    fprintf(fp, "theorem q4_bounded (u v : BitVec 4) : BoundedBV32 10000000 (q4_bv u v) := by\n");
    fprintf(fp, "  unfold BoundedBV32\n  revert u v\n  decide\n\n");

    fprintf(fp, "theorem continuity4_bv32 (a b c d : BitVec 4) :\n");
    fprintf(fp, "    q4_bv (step3 ruleP a b c) (step3 ruleP b c d) - q4_bv b c =\n");
    fprintf(fp, "      canonCurrent4_32 ruleP q4_bv b c d - canonCurrent4_32 ruleP q4_bv a b c := "
                "by\n");
    fprintf(fp, "  simp only [step3, canonCurrent4_32, pack2, hi2, lo2, ruleP, q4_bv]\n");
    fprintf(fp, "  bv_decide\n\n");

    fprintf(fp, "theorem global_conservation4 {N : Nat} (hN : 0 < N) (u : Fin N → BitVec 4) :\n");
    fprintf(fp, "    ringSum N (fun i => (q4_bv (floquetRing hN ruleP u i) "
                "(floquetRing hN ruleP u (nextMod hN i))).toInt) =\n");
    fprintf(fp, "      ringSum N (fun i => (q4_bv (u i) (u (nextMod hN i))).toInt) :=\n");
    fprintf(fp, "  forall_ring_conservation4_of_bv32 hN ruleP q4_bv q4_bounded continuity4_bv32 "
                "u\n\n");

    if (q4->cross_axis == 0U) {
        fprintf(fp, "theorem q4_nontrivial : ¬ IsTrivialQ4 (fun u v => (q4_bv u v).toInt) :=\n");
        fprintf(fp, "  non_trivial_q4_of_inner_cross %u#2 %u#2 %u#2 %u#2 %u#2 %u#2 (by decide)\n\n",
                (unsigned)q4->cross_f1, (unsigned)q4->cross_f2, (unsigned)q4->cross_x0,
                (unsigned)q4->cross_x1, (unsigned)q4->cross_y0, (unsigned)q4->cross_y1);
    } else {
        fprintf(fp, "theorem q4_nontrivial : ¬ IsTrivialQ4 (fun u v => (q4_bv u v).toInt) :=\n");
        fprintf(fp, "  non_trivial_q4_of_outer_cross %u#2 %u#2 %u#2 %u#2 %u#2 %u#2 (by decide)\n\n",
                (unsigned)q4->cross_f1, (unsigned)q4->cross_f2, (unsigned)q4->cross_x0,
                (unsigned)q4->cross_x1, (unsigned)q4->cross_y0, (unsigned)q4->cross_y1);
    }

    if (q4->j4_nonzero > 0U) {
        fprintf(fp, "theorem j4_active : ∃ a b c, canonCurrent4 ruleP "
                    "(fun u v => (q4_bv u v).toInt) a b c ≠ 0 :=\n");
        fprintf(fp,
                "  active_current4_witness ruleP (fun u v => (q4_bv u v).toInt) "
                "%u#4 %u#4 %u#4 (by decide)\n\n",
                (unsigned)q4->j4_wit_a, (unsigned)q4->j4_wit_b, (unsigned)q4->j4_wit_c);
    }

    fprintf(fp, "end Phy.DiscoveredQ4\n");
    fclose(fp);
    return true;
}

static inline bool phy_emit_lean_yb8_cert(const phy_yb8_analysis_t *yb8, const char *path) {
    if (yb8 == NULL || path == NULL || yb8->is_bijective == 0U || yb8->yb8_violations != 0U ||
        yb8->is_nondeg == 0U || yb8->is_involutive != 0U || yb8->is_indecomposable == 0U ||
        yb8->is_lyubashenko != 0U || yb8->anf_deg < 2U) {
        return false;
    }

    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }

    fprintf(fp, "import Phy.Noether\n\n");
    fprintf(fp, "namespace Phy.DiscoveredYB8\n\n");
    fprintf(fp, "open Phy.Noether\n\n");
    fprintf(fp, "-- Merkle: 0x%016llx, order=%u, ANF deg=%u, nondeg=1, indecomp=1, noninvol=1\n\n",
            (unsigned long long)yb8->merkle, (unsigned)yb8->order, (unsigned)yb8->anf_deg);

    phy_emit_bv6_lut_tree(fp, "ruleR8", yb8->r);
    phy_emit_bv6_lut_tree(fp, "ruleR8Inv", yb8->rinv);
    phy_emit_bv3_pair_tree(fp, "sigInv8", "a", "s", yb8->sig_inv);
    phy_emit_bv3_pair_tree(fp, "tauInv8", "b", "t", yb8->tau_inv);

    fprintf(fp, "theorem ruleR8_inv_left (u : BitVec 6) : ruleR8Inv (ruleR8 u) = u := by\n");
    fprintf(fp, "  simp only [ruleR8, ruleR8Inv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleR8_inv_right (u : BitVec 6) : ruleR8 (ruleR8Inv u) = u := by\n");
    fprintf(fp, "  simp only [ruleR8, ruleR8Inv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleR8_bijective : Bijective ruleR8 := by\n");
    fprintf(fp, "  refine ⟨fun a b h => ?_, fun b => ⟨ruleR8Inv b, ruleR8_inv_right b⟩⟩\n");
    fprintf(fp, "  have h1 := congrArg ruleR8Inv h\n");
    fprintf(fp, "  simpa [ruleR8_inv_left] using h1\n\n");

    fprintf(fp, "theorem ruleR8_satisfies_ybe : IsSetTheoreticYB8 ruleR8 := by\n");
    fprintf(fp, "  intro t\n");
    fprintf(fp, "  simp only [yb8Apply12, yb8Apply23, pack3, ruleR8]\n");
    fprintf(fp, "  bv_decide\n\n");

    fprintf(fp, "theorem ruleR8_nondegenerate : IsNonDegenerate8 ruleR8 :=\n");
    fprintf(fp, "  non_degenerate8_of_inverses ruleR8 sigInv8 tauInv8\n");
    fprintf(fp, "    (by intro a b; simp only [hi3, pack3, ruleR8, sigInv8]; bv_decide)\n");
    fprintf(fp, "    (by intro a b; simp only [hi3, pack3, ruleR8, sigInv8]; bv_decide)\n");
    fprintf(fp, "    (by intro a b; simp only [lo3, pack3, ruleR8, tauInv8]; bv_decide)\n");
    fprintf(fp, "    (by intro a b; simp only [lo3, pack3, ruleR8, tauInv8]; bv_decide)\n\n");

    fprintf(fp, "theorem ruleR8_noninvolutive : ¬ IsInvolutive8 ruleR8 :=\n");
    fprintf(fp, "  non_involutive8_of_witness %u#6 (by decide)\n\n", (unsigned)yb8->noninvol_u);

    fprintf(fp, "theorem ruleR8_indecomposable : ¬ IsDecomposable8 ruleR8 :=\n");
    fprintf(fp, "  indecomposable8_of_connected8 ruleR8 (by decide)\n\n");

    if (yb8->nly_axis == 0U) {
        fprintf(fp, "theorem ruleR8_nonlyubashenko : ¬ IsLyubashenko8 ruleR8 :=\n");
        fprintf(fp, "  non_lyubashenko8_of_hi_cross %u#3 %u#3 %u#3 (by decide)\n\n",
                (unsigned)yb8->nly_a, (unsigned)yb8->nly_b0, (unsigned)yb8->nly_b1);
    } else {
        fprintf(fp, "theorem ruleR8_nonlyubashenko : ¬ IsLyubashenko8 ruleR8 :=\n");
        fprintf(fp, "  non_lyubashenko8_of_lo_cross %u#3 %u#3 %u#3 (by decide)\n\n",
                (unsigned)yb8->nly_a, (unsigned)yb8->nly_b0, (unsigned)yb8->nly_b1);
    }

    fprintf(fp,
            "theorem ruleR8_nonlinear : ¬ IsAffineGF2_6 (fun x => BitVec.getLsbD (ruleR8 x) %u) "
            ":=\n",
            (unsigned)yb8->nl_bit);
    fprintf(fp, "  gf2_6_second_diff_witnesses_nonlinear %u#6 %u#6 (by decide)\n\n",
            (unsigned)yb8->nl_a, (unsigned)yb8->nl_b);

    fprintf(fp, "end Phy.DiscoveredYB8\n");
    fclose(fp);
    return true;
}

static inline void phy_cert_build_anf6_fir(const uint8_t r[64], flow_ir_prog_t *p,
                                           uint32_t const0_id, uint32_t const1_id,
                                           const uint32_t x_ids[6], uint32_t yb_obits[6],
                                           uint16_t *out_and, uint16_t *out_depth) {
    uint8_t anf[6][64];
    for (uint32_t o = 0U; o < 6U; ++o) {
        for (uint32_t u = 0U; u < 64U; ++u) {
            anf[o][u] = (uint8_t)(((uint32_t)r[u] >> o) & UINT32_C(1));
        }
        for (uint32_t bit = 0U; bit < 6U; ++bit) {
            uint32_t step = UINT32_C(1) << bit;
            for (uint32_t mask = 0U; mask < 64U; ++mask) {
                if ((mask & step) != UINT32_C(0)) {
                    anf[o][mask] ^= anf[o][mask ^ step];
                }
            }
        }
    }
    uint8_t needed[64];
    memset(needed, 0, sizeof(needed));
    for (uint32_t o = 0U; o < 6U; ++o) {
        for (uint32_t mask = 1U; mask < 64U; ++mask) {
            if (anf[o][mask] != 0U) {
                needed[mask] = 1U;
            }
        }
    }
    for (uint32_t mask = 63U; mask >= 1U; --mask) {
        if (needed[mask] != 0U) {
            uint32_t sub = mask & (mask - UINT32_C(1));
            if (sub != UINT32_C(0)) {
                needed[sub] = 1U;
            }
        }
    }

    uint32_t mon_vid[64];
    uint16_t mon_depth[64];
    memset(mon_vid, 0, sizeof(mon_vid));
    memset(mon_depth, 0, sizeof(mon_depth));
    for (uint32_t b = 0U; b < 6U; ++b) {
        mon_vid[UINT32_C(1) << b] = x_ids[b];
        mon_depth[UINT32_C(1) << b] = 0U;
    }

    uint16_t and_cnt = 0U;
    uint16_t max_d = 0U;
    for (uint32_t mask = 1U; mask < 64U; ++mask) {
        if (needed[mask] == 0U) {
            continue;
        }
        uint32_t sub = mask & (mask - UINT32_C(1));
        if (sub == UINT32_C(0)) {
            continue;
        }
        uint32_t lsb = mask & (~mask + UINT32_C(1));
        char nm[16];
        snprintf(nm, sizeof(nm), "m%02x", mask);
        mon_vid[mask] = phy_cert_emit_binop(p, FLOW_IR_OP_AND_VEC, mon_vid[sub], mon_vid[lsb], nm);
        uint16_t d = (uint16_t)(mon_depth[sub] + 1U);
        mon_depth[mask] = d;
        if (d > max_d) {
            max_d = d;
        }
        ++and_cnt;
    }

    for (uint32_t o = 0U; o < 6U; ++o) {
        uint32_t acc = UINT32_MAX;
        uint32_t tidx = 0U;
        if (anf[o][0] != 0U) {
            acc = const1_id;
        }
        for (uint32_t mask = 1U; mask < 64U; ++mask) {
            if (anf[o][mask] != 0U) {
                if (acc == UINT32_MAX) {
                    acc = mon_vid[mask];
                } else {
                    char nm[16];
                    snprintf(nm, sizeof(nm), "z%u_%u", o, tidx);
                    acc = phy_cert_emit_binop(p, FLOW_IR_OP_XOR_VEC, acc, mon_vid[mask], nm);
                    ++tidx;
                }
            }
        }
        if (acc == UINT32_MAX) {
            acc = const0_id;
        }
        yb_obits[o] = acc;
    }
    *out_and = and_cnt;
    *out_depth = max_d;
}

static inline uint32_t phy_cert_spec_max_deg6(const searu_mc_spec *spec) {
    uint32_t max_d = 0U;
    for (uint32_t o = 0U; o < 6U; ++o) {
        uint64_t a = spec->target_tt[o][0];
        a ^= (a & UINT64_C(0x5555555555555555)) << 1U;
        a ^= (a & UINT64_C(0x3333333333333333)) << 2U;
        a ^= (a & UINT64_C(0x0F0F0F0F0F0F0F0F)) << 4U;
        a ^= (a & UINT64_C(0x00FF00FF00FF00FF)) << 8U;
        a ^= (a & UINT64_C(0x0000FFFF0000FFFF)) << 16U;
        a ^= (a & UINT64_C(0x00000000FFFFFFFF)) << 32U;
        if (a != UINT64_C(0)) {
            for (uint32_t m = 0U; m < 64U; ++m) {
                if (((a >> m) & UINT64_C(1)) != UINT64_C(0)) {
                    uint32_t c = 0U;
                    uint32_t mm = m;
                    while (mm != 0U) {
                        c += (mm & 1U);
                        mm >>= 1U;
                    }
                    if (c > max_d) {
                        max_d = c;
                    }
                }
            }
        }
    }
    return max_d;
}

static inline bool phy_synth_fir_6bit(const uint8_t r[64], uint64_t seed, const char *mod_name,
                                      const char *out_fir_path, phy_fir_cert_ws_t *ws) {
    if (r == NULL || ws == NULL) {
        return false;
    }
    ws->and_count = 0U;
    ws->mult_depth = 0U;
    ws->wire_len = 0U;
    ws->ir_linear_ok = false;
    ws->ir_unitary_ok = false;
    ws->ir_eval_ok = false;

    searu_mc_spec spec;
    memset(&spec, 0, sizeof(spec));
    spec.n_in = 6U;
    spec.n_out = 6U;
    spec.max_and = 6U;
    spec.max_depth = 0U;
    spec.n_in_port0 = 0U;
    spec.budget_scale = 1U;
    spec.seed = seed;

    for (uint32_t o = 0U; o < 6U; ++o) {
        uint64_t tt = UINT64_C(0);
        for (uint32_t u = 0U; u < 64U; ++u) {
            if (((r[u] >> o) & 1U) != 0U) {
                tt |= (UINT64_C(1) << u);
            }
        }
        spec.target_tt[o][0] = tt;
    }

    uint32_t max_deg6 = phy_cert_spec_max_deg6(&spec);
    bool mc_ok = (max_deg6 < 5U) ? searu_solve_mc_gf2(&ws->mc_plan, &spec, &ws->mc_ws) : false;

    flow_ir_prog_t *p = &ws->ir_prog;
    memset(p, 0, sizeof(*p));
    snprintf(p->name, sizeof(p->name), "%s", (mod_name != NULL) ? mod_name : "phy_yb8_3x3");
    p->mult = FLOW_IR_MULT_LINEAR;
    p->flags = FLOW_IR_FLAG_CLEAN;
    p->num_inputs = 1U;
    p->inputs[0] = 0U;
    p->vals[0].width = 6U;
    p->vals[0].is_classical = 0U;
    snprintf(p->vals[0].name, sizeof(p->vals[0].name), "u6");
    p->num_vals = 1U;
    p->num_ops = 0U;

    uint32_t yb_basis[SEARU_MC_MAX_BASIS];
    uint32_t const0_id = phy_cert_emit_const(p, UINT64_C(0), "c0");
    yb_basis[0] = phy_cert_emit_const(p, UINT64_C(1), "c1");
    for (uint32_t i = 0U; i < 6U; ++i) {
        char nm[16];
        snprintf(nm, sizeof(nm), "x%u", i);
        yb_basis[1U + i] = phy_cert_emit_slice(p, 0U, i, nm);
    }

    uint32_t yb_obits[6];
    if (mc_ok) {
        ws->and_count = (uint16_t)ws->mc_plan.num_and;
        ws->mult_depth = (uint16_t)ws->mc_plan.mult_depth;
        for (uint32_t k = 0U; k < ws->mc_plan.num_and; ++k) {
            uint32_t n_cur = 7U + k;
            char lnm[16];
            char rnm[16];
            char anm[16];
            snprintf(lnm, sizeof(lnm), "yl%u_", k);
            snprintf(rnm, sizeof(rnm), "yr%u_", k);
            snprintf(anm, sizeof(anm), "ya%u", k);
            uint32_t lhs = phy_cert_emit_xor_mask(p, yb_basis, n_cur, ws->mc_plan.and_lhs_mask[k],
                                                  const0_id, lnm);
            uint32_t rhs = phy_cert_emit_xor_mask(p, yb_basis, n_cur, ws->mc_plan.and_rhs_mask[k],
                                                  const0_id, rnm);
            yb_basis[7U + k] = phy_cert_emit_binop(p, FLOW_IR_OP_AND_VEC, lhs, rhs, anm);
            if (yb_basis[7U + k] == UINT32_MAX) {
                return false;
            }
        }
        uint32_t n_all = 7U + ws->mc_plan.num_and;
        for (uint32_t o = 0U; o < 6U; ++o) {
            char pfx[16];
            snprintf(pfx, sizeof(pfx), "yo%u_", o);
            yb_obits[o] =
                phy_cert_emit_xor_mask(p, yb_basis, n_all, ws->mc_plan.out_mask[o], const0_id, pfx);
            if (yb_obits[o] == UINT32_MAX) {
                return false;
            }
        }
    } else {
        phy_cert_build_anf6_fir(r, p, const0_id, yb_basis[0], &yb_basis[1], yb_obits,
                                &ws->and_count, &ws->mult_depth);
        for (uint32_t o = 0U; o < 6U; ++o) {
            if (yb_obits[o] == UINT32_MAX) {
                return false;
            }
        }
    }

    if (p->num_vals >= FLOW_IR_MAX_VALS || p->num_ops >= FLOW_IR_MAX_OPS) {
        return false;
    }
    uint32_t yb_out_vid = p->num_vals;
    uint32_t yb_out_op = p->num_ops;
    p->vals[yb_out_vid].width = 6U;
    p->vals[yb_out_vid].is_classical = 0U;
    snprintf(p->vals[yb_out_vid].name, sizeof(p->vals[yb_out_vid].name), "out6");
    p->ops[yb_out_op].kind = FLOW_IR_OP_CONCAT;
    p->ops[yb_out_op].num_ins = 6U;
    p->ops[yb_out_op].out_val = yb_out_vid;
    for (uint32_t o = 0U; o < 6U; ++o) {
        p->ops[yb_out_op].ins[o] = yb_obits[o];
    }
    p->num_vals++;
    p->num_ops++;
    p->num_outputs = 1U;
    p->outputs[0] = yb_out_vid;

    if (!flow_ir_encode_bytes(p, ws->wire_buf, sizeof(ws->wire_buf), &ws->wire_len)) {
        return false;
    }
    if (flow_ir_load_bytes_fiber(&ws->ir_reloaded, ws->wire_buf, ws->wire_len) !=
        FLOW_IR_FIBER_COMMITTED) {
        return false;
    }
    ws->ir_linear_ok = (ws->ir_reloaded.mult == FLOW_IR_MULT_LINEAR);
    ws->ir_unitary_ok = flow_ir_admits_unitary(&ws->ir_reloaded);

    bool match = true;
    for (uint32_t u = 0U; u < 64U; ++u) {
        uint64_t in_arr[1] = {(uint64_t)u};
        uint64_t out_arr[1] = {UINT64_C(0)};
        if (!flow_ir_eval_u64(&ws->ir_reloaded, in_arr, out_arr)) {
            match = false;
            break;
        }
        if ((out_arr[0] & UINT64_C(0x3F)) != (uint64_t)r[u]) {
            match = false;
            break;
        }
    }
    ws->ir_eval_ok = match;

    if (out_fir_path != NULL && out_fir_path[0] != '\0') {
        FILE *fp = fopen(out_fir_path, "wb");
        if (fp == NULL) {
            return false;
        }
        size_t nw = fwrite(ws->wire_buf, 1U, ws->wire_len, fp);
        fclose(fp);
        if (nw != ws->wire_len) {
            return false;
        }
    }

    return ws->ir_linear_ok && ws->ir_unitary_ok && ws->ir_eval_ok;
}

static inline bool phy_synth_fir_yb8(const phy_yb8_analysis_t *yb8, const char *out_fir_path,
                                     phy_fir_cert_ws_t *ws) {
    if (yb8 == NULL) {
        return false;
    }
    return phy_synth_fir_6bit(yb8->r, yb8->merkle ^ UINT64_C(0x59423832), "phy_yb8_3x3",
                              out_fir_path, ws);
}

static inline bool phy_emit_q4_json(const phy_q4_analysis_t *q4, const phy_butterfly_t *bf,
                                    const phy_fir_cert_ws_t *fir, const char *path) {
    if (q4 == NULL || path == NULL) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "{\n");
    fprintf(fp, "  \"p\": [");
    for (uint32_t i = 0U; i < 16U; ++i) {
        fprintf(fp, "%s%u", (i > 0U) ? ", " : "", (unsigned)q4->p[i]);
    }
    fprintf(fp, "],\n");
    fprintf(fp, "  \"q4\": [");
    for (uint32_t i = 0U; i < 256U; ++i) {
        fprintf(fp, "%s%ld", (i > 0U) ? ", " : "", (long)q4->q4[i]);
    }
    fprintf(fp, "],\n");
    fprintf(fp, "  \"lut_packed\": \"0x%016llx\",\n", (unsigned long long)q4->lut_packed);
    fprintf(fp, "  \"q4_merkle\": \"0x%016llx\",\n", (unsigned long long)q4->q4_merkle);
    fprintf(fp, "  \"tot2\": %u,\n", (unsigned)q4->tot2);
    fprintf(fp, "  \"tot4\": %u,\n", (unsigned)q4->tot4);
    fprintf(fp, "  \"irred4\": %u,\n", (unsigned)q4->irred4);
    fprintf(fp, "  \"genuine_q4\": %u,\n", (unsigned)q4->genuine_q4);
    fprintf(fp, "  \"j4_nonzero\": %u,\n", (unsigned)q4->j4_nonzero);
    fprintf(fp, "  \"q4_nnz\": %u,\n", (unsigned)q4->q4_nnz);
    fprintf(fp, "  \"continuity4_fail\": %d", q4->continuity4_fail);
    if (bf != NULL) {
        fprintf(fp, ",\n  \"bf_max_span\": %u,\n", (unsigned)bf->max_span);
        fprintf(fp, "  \"bf_mean_span_q16\": %u,\n", (unsigned)bf->mean_span_q16);
        fprintf(fp, "  \"bf_max_vel_q16\": %d,\n", bf->max_vel_q16);
        fprintf(fp, "  \"bf_mean_vel_q16\": %d,\n", bf->mean_vel_q16);
        fprintf(fp, "  \"bf_active_probes\": %u,\n", (unsigned)bf->active_probes);
        fprintf(fp, "  \"bf_max_drift\": %lld", (long long)bf->max_charge_drift);
    }
    if (fir != NULL) {
        fprintf(fp, ",\n  \"fir_and_count\": %u,\n", (unsigned)fir->and_count);
        fprintf(fp, "  \"fir_mult_depth\": %u,\n", (unsigned)fir->mult_depth);
        fprintf(fp, "  \"fir_wire_bytes\": %zu,\n", fir->wire_len);
        fprintf(fp, "  \"fir_eval_ok\": %s", fir->ir_eval_ok ? "true" : "false");
    }
    fprintf(fp, "\n}\n");
    fclose(fp);
    return true;
}

static inline bool phy_emit_yb8_json(const phy_yb8_analysis_t *yb8, const phy_fir_cert_ws_t *fir,
                                     const char *path) {
    if (yb8 == NULL || path == NULL) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "{\n");
    fprintf(fp, "  \"r\": [");
    for (uint32_t i = 0U; i < 64U; ++i) {
        fprintf(fp, "%s%u", (i > 0U) ? ", " : "", (unsigned)yb8->r[i]);
    }
    fprintf(fp, "],\n");
    fprintf(fp, "  \"merkle\": \"0x%016llx\",\n", (unsigned long long)yb8->merkle);
    fprintf(fp, "  \"flags\": %u,\n", (unsigned)yb8->flags);
    fprintf(fp, "  \"order\": %u,\n", (unsigned)yb8->order);
    fprintf(fp, "  \"yb8_violations\": %u,\n", (unsigned)yb8->yb8_violations);
    fprintf(fp, "  \"anf_deg\": %u,\n", (unsigned)yb8->anf_deg);
    fprintf(fp, "  \"is_bijective\": %u,\n", (unsigned)yb8->is_bijective);
    fprintf(fp, "  \"is_nondeg\": %u,\n", (unsigned)yb8->is_nondeg);
    fprintf(fp, "  \"is_involutive\": %u,\n", (unsigned)yb8->is_involutive);
    fprintf(fp, "  \"is_indecomposable\": %u,\n", (unsigned)yb8->is_indecomposable);
    fprintf(fp, "  \"is_lyubashenko\": %u,\n", (unsigned)yb8->is_lyubashenko);
    fprintf(fp, "  \"noninvol_u\": %u", (unsigned)yb8->noninvol_u);
    if (fir != NULL) {
        fprintf(fp, ",\n  \"fir_and_count\": %u,\n", (unsigned)fir->and_count);
        fprintf(fp, "  \"fir_mult_depth\": %u,\n", (unsigned)fir->mult_depth);
        fprintf(fp, "  \"fir_wire_bytes\": %zu,\n", fir->wire_len);
        fprintf(fp, "  \"fir_eval_ok\": %s", fir->ir_eval_ok ? "true" : "false");
    }
    fprintf(fp, "\n}\n");
    fclose(fp);
    return true;
}

/* --- PHY v3.0: Dual-Unitary, AME(4,8), Fracton HSF, Biquandle8 Knots, Tetrahedron --- */

static inline void phy_emit_bv2_pair_tree(FILE *fp, const char *fn_name, const char *p1,
                                          const char *p2, const uint8_t tbl[4][4]) {
    int32_t tmp[16];
    for (uint32_t a = 0U; a < 4U; ++a) {
        for (uint32_t b = 0U; b < 4U; ++b) {
            tmp[(a << 2U) | b] = (int32_t)tbl[a][b];
        }
    }
    fprintf(fp, "def %s (%s %s : BitVec 2) : BitVec 2 :=\n", fn_name, p1, p2);
    fprintf(fp, "  let key : BitVec 4 := %s ++ %s\n  ", p1, p2);
    phy_emit_bv_tree_rec(fp, "key", 4U, 2U, 0U, 16U, tmp, 2U);
    fprintf(fp, "\n\n");
}

static inline bool phy_emit_lean_dual4_cert(const phy_dual4_analysis_t *d4, const char *path) {
    if (d4 == NULL || path == NULL || d4->is_dual_unitary == 0U || d4->is_mixing == 0U ||
        d4->anf_deg < 2U) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }

    fprintf(fp, "import Phy.Noether\n\n");
    fprintf(fp, "namespace Phy.DiscoveredDual4\n\n");
    fprintf(fp, "open Phy.Noether\n\n");
    fprintf(fp, "-- Merkle: 0x%016llx, ANF deg=%u, DU=1, mixing=1, nonBern=%u, rank_ent=%u\n\n",
            (unsigned long long)d4->merkle, (unsigned)d4->anf_deg, (unsigned)d4->is_non_bernoulli,
            (unsigned)d4->rank_ent);

    phy_emit_bv4_lut_lean(fp, "ruleDU4", d4->p);
    phy_emit_bv4_lut_lean(fp, "ruleDU4Inv", d4->pinv);
    phy_emit_bv2_pair_tree(fp, "sigInv4", "a", "c", d4->sig_inv);
    phy_emit_bv2_pair_tree(fp, "tauInv4", "b", "d", d4->tau_inv);

    fprintf(fp, "theorem ruleDU4_inv_left (u : BitVec 4) : ruleDU4Inv (ruleDU4 u) = u := by\n");
    fprintf(fp, "  revert u\n  decide\n\n");

    fprintf(fp, "theorem ruleDU4_inv_right (u : BitVec 4) : ruleDU4 (ruleDU4Inv u) = u := by\n");
    fprintf(fp, "  revert u\n  decide\n\n");

    fprintf(fp, "theorem ruleDU4_bijective : Bijective ruleDU4 := by\n");
    fprintf(fp, "  refine ⟨fun a b h => ?_, fun b => ⟨ruleDU4Inv b, ruleDU4_inv_right b⟩⟩\n");
    fprintf(fp, "  have h1 := congrArg ruleDU4Inv h\n");
    fprintf(fp, "  simpa [ruleDU4_inv_left] using h1\n\n");

    fprintf(fp, "theorem ruleDU4_nondegenerate : IsNonDegenerate4 ruleDU4 :=\n");
    fprintf(fp, "  non_degenerate4_of_inverses ruleDU4 sigInv4 tauInv4\n");
    fprintf(fp, "    (by intro a b; revert a b; decide)\n");
    fprintf(fp, "    (by intro a b; revert a b; decide)\n");
    fprintf(fp, "    (by intro a b; revert a b; decide)\n");
    fprintf(fp, "    (by intro a b; revert a b; decide)\n\n");

    fprintf(fp, "theorem ruleDU4_dual_unitary : IsDualUnitary4 ruleDU4 :=\n");
    fprintf(fp, "  ⟨ruleDU4_bijective, ruleDU4_nondegenerate⟩\n\n");

    fprintf(fp, "theorem ruleDU4_q2_nogo (q : BitVec 4 → Int) (j : BitVec 4 → BitVec 4 → Int)\n");
    fprintf(fp, "    (hCont : SatisfiesContinuity ruleDU4 q j) : IsSeparable q :=\n");
    fprintf(fp, "  dual_unitary_q2_separable_nogo ruleDU4 ruleDU4_dual_unitary q j hCont\n\n");

    fprintf(fp, "theorem ruleDU4_mixing : IsMixingDualUnitary4 ruleDU4 :=\n");
    fprintf(fp, "  ⟨ruleDU4_dual_unitary, by decide, by decide⟩\n\n");

    fprintf(fp,
            "theorem ruleDU4_nonlinear : ¬ IsAffineGF2 (fun x => BitVec.getLsbD (ruleDU4 x) %u) "
            ":=\n",
            (unsigned)d4->nl_bit);
    fprintf(fp, "  gf2_second_diff_witnesses_nonlinear %u#4 %u#4 (by decide)\n\n",
            (unsigned)d4->nl_a, (unsigned)d4->nl_b);

    fprintf(fp, "end Phy.DiscoveredDual4\n");
    fclose(fp);
    return true;
}

static inline bool phy_emit_lean_ame8_cert(const phy_ame8_analysis_t *a8, const char *path) {
    if (a8 == NULL || path == NULL || a8->is_2unitary == 0U || a8->anf_deg < 2U) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }

    fprintf(fp, "import Phy.Noether\n\n");
    fprintf(fp, "namespace Phy.DiscoveredAME8\n\n");
    fprintf(fp, "open Phy.Noether\n\n");
    fprintf(fp, "-- Merkle: 0x%016llx, order=%u, ANF deg=%u, Is2Unitary8=1\n\n",
            (unsigned long long)a8->merkle, (unsigned)a8->order, (unsigned)a8->anf_deg);

    phy_emit_bv6_lut_tree(fp, "ruleAME8", a8->r);
    phy_emit_bv6_lut_tree(fp, "ruleAME8Inv", a8->rinv);
    phy_emit_bv3_pair_tree(fp, "sigInv8", "a", "c", a8->sig_inv);
    phy_emit_bv3_pair_tree(fp, "tauInv8", "b", "d", a8->tau_inv);
    phy_emit_bv3_pair_tree(fp, "crossLInv8", "a", "d", a8->cross_l_inv);
    phy_emit_bv3_pair_tree(fp, "crossRInv8", "b", "c", a8->cross_r_inv);

    fprintf(fp, "theorem ruleAME8_inv_left (u : BitVec 6) : ruleAME8Inv (ruleAME8 u) = u := by\n");
    fprintf(fp, "  simp only [ruleAME8, ruleAME8Inv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleAME8_inv_right (u : BitVec 6) : ruleAME8 (ruleAME8Inv u) = u := by\n");
    fprintf(fp, "  simp only [ruleAME8, ruleAME8Inv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleAME8_bijective : Bijective ruleAME8 := by\n");
    fprintf(fp, "  refine ⟨fun a b h => ?_, fun b => ⟨ruleAME8Inv b, ruleAME8_inv_right b⟩⟩\n");
    fprintf(fp, "  have h1 := congrArg ruleAME8Inv h\n");
    fprintf(fp, "  simpa [ruleAME8_inv_left] using h1\n\n");

    fprintf(fp, "theorem ruleAME8_nondegenerate : IsNonDegenerate8 ruleAME8 :=\n");
    fprintf(fp, "  non_degenerate8_of_inverses ruleAME8 sigInv8 tauInv8\n");
    fprintf(fp, "    (by intro a b; simp only [hi3, pack3, ruleAME8, sigInv8]; bv_decide)\n");
    fprintf(fp, "    (by intro a b; simp only [hi3, pack3, ruleAME8, sigInv8]; bv_decide)\n");
    fprintf(fp, "    (by intro a b; simp only [lo3, pack3, ruleAME8, tauInv8]; bv_decide)\n");
    fprintf(fp, "    (by intro a b; simp only [lo3, pack3, ruleAME8, tauInv8]; bv_decide)\n\n");

    fprintf(fp, "theorem ruleAME8_crossbijective : IsCrossBijective8 ruleAME8 :=\n");
    fprintf(fp, "  cross_bijective8_of_inverses ruleAME8 crossLInv8 crossRInv8\n");
    fprintf(fp, "    (by intro a b; simp only [lo3, pack3, ruleAME8, crossLInv8]; bv_decide)\n");
    fprintf(fp, "    (by intro a b; simp only [lo3, pack3, ruleAME8, crossLInv8]; bv_decide)\n");
    fprintf(fp, "    (by intro a b; simp only [hi3, pack3, ruleAME8, crossRInv8]; bv_decide)\n");
    fprintf(fp, "    (by intro a b; simp only [hi3, pack3, ruleAME8, crossRInv8]; bv_decide)\n\n");

    fprintf(fp, "theorem ruleAME8_is_2unitary : Is2Unitary8 ruleAME8 :=\n");
    fprintf(fp, "  ⟨ruleAME8_bijective, ruleAME8_nondegenerate, ruleAME8_crossbijective⟩\n\n");

    fprintf(fp,
            "theorem ruleAME8_nonlinear : ¬ IsAffineGF2_6 (fun x => BitVec.getLsbD (ruleAME8 x) "
            "%u) :=\n",
            (unsigned)a8->nl_bit);
    fprintf(fp, "  gf2_6_second_diff_witnesses_nonlinear %u#6 %u#6 (by decide)\n\n",
            (unsigned)a8->nl_a, (unsigned)a8->nl_b);

    fprintf(fp, "end Phy.DiscoveredAME8\n");
    fclose(fp);
    return true;
}

static inline bool phy_emit_lean_fracton_cert(const phy_fracton_analysis_t *fr, const char *path) {
    if (fr == NULL || path == NULL || fr->is_bijective == 0U || fr->laplacian_ok == 0U ||
        fr->active_entries == 0U || fr->anf_deg < 3U) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }

    uint32_t wu = (uint32_t)fr->wit_u;
    uint32_t wa = (wu >> 4U) & UINT32_C(3);
    uint32_t wb = (wu >> 2U) & UINT32_C(3);
    uint32_t wc = wu & UINT32_C(3);

    fprintf(fp, "import Phy.Noether\n\n");
    fprintf(fp, "namespace Phy.DiscoveredFracton\n\n");
    fprintf(fp, "open Phy.Noether\n\n");
    fprintf(fp,
            "-- Merkle: 0x%016llx, HSF comps=%u/%u (largest=%u, frozen=%u), iso_span=%u, "
            "dip_span=%u\n\n",
            (unsigned long long)fr->merkle, (unsigned)fr->hsf_num_comps,
            (unsigned)fr->hsf_sector_size, (unsigned)fr->hsf_largest_comp, (unsigned)fr->hsf_frozen,
            (unsigned)fr->iso_span, (unsigned)fr->dip_span);

    phy_emit_bv6_lut_tree(fp, "ruleFrac", fr->t);
    phy_emit_bv6_lut_tree(fp, "ruleFracInv", fr->tinv);

    fprintf(fp, "theorem ruleFrac_inv_left (u : BitVec 6) : ruleFracInv (ruleFrac u) = u := by\n");
    fprintf(fp, "  simp only [ruleFrac, ruleFracInv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleFrac_inv_right (u : BitVec 6) : ruleFrac (ruleFracInv u) = u := by\n");
    fprintf(fp, "  simp only [ruleFrac, ruleFracInv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleFrac_bijective : Bijective ruleFrac := by\n");
    fprintf(fp, "  refine ⟨fun a b h => ?_, fun b => ⟨ruleFracInv b, ruleFrac_inv_right b⟩⟩\n");
    fprintf(fp, "  have h1 := congrArg ruleFracInv h\n");
    fprintf(fp, "  simpa [ruleFrac_inv_left] using h1\n\n");

    fprintf(fp,
            "theorem ruleFrac_dipole_laplacian : SatisfiesDipoleLaplacian ruleFrac rho2 := by\n");
    fprintf(fp, "  unfold SatisfiesDipoleLaplacian rho2\n  decide\n\n");

    fprintf(fp,
            "theorem ruleFrac_global_Q_conservation {L : Nat} (pos : Nat) (hpos : pos + 2 < L)\n");
    fprintf(fp, "    (s : Fin L → BitVec 2) :\n");
    fprintf(fp,
            "    openCharge L rho2 (applyGate3 ruleFrac pos hpos s) = openCharge L rho2 s :=\n");
    fprintf(fp,
            "  fracton_open_chain_Q_conservation ruleFrac rho2 ruleFrac_dipole_laplacian pos hpos "
            "s\n\n");

    fprintf(fp,
            "theorem ruleFrac_global_D_conservation {L : Nat} (pos : Nat) (hpos : pos + 2 < L)\n");
    fprintf(fp, "    (s : Fin L → BitVec 2) :\n");
    fprintf(fp,
            "    openDipole L rho2 (applyGate3 ruleFrac pos hpos s) = openDipole L rho2 s :=\n");
    fprintf(fp,
            "  fracton_open_chain_D_conservation ruleFrac rho2 ruleFrac_dipole_laplacian pos hpos "
            "s\n\n");

    fprintf(fp,
            "theorem ruleFrac_single_charge_frozen {L : Nat} (pos : Nat) (hpos : pos + 2 < L)\n");
    fprintf(fp, "    (s : Fin L → BitVec 2) (hQ : openCharge L rho2 s ≤ 1) :\n");
    fprintf(fp, "    applyGate3 ruleFrac pos hpos s = s :=\n");
    fprintf(fp, "  fracton_single_charge_frozen ruleFrac ruleFrac_dipole_laplacian pos hpos s "
                "hQ\n\n");

    fprintf(fp, "theorem ruleFrac_quadrupole_active : ∃ x y z : BitVec 2, rho2 (site0_3 (ruleFrac "
                "(pack3_2 x y z))) ≠ rho2 x :=\n");
    fprintf(fp, "  fracton_quadrupole_active_witness ruleFrac rho2 %u#2 %u#2 %u#2 (by decide)\n\n",
            wa, wb, wc);

    fprintf(fp,
            "theorem ruleFrac_nonlinear : ¬ IsAffineGF2_6 (fun x => BitVec.getLsbD (ruleFrac x) "
            "%u) :=\n",
            (unsigned)fr->nl_bit);
    fprintf(fp, "  gf2_6_second_diff_witnesses_nonlinear %u#6 %u#6 (by decide)\n\n",
            (unsigned)fr->nl_a, (unsigned)fr->nl_b);

    fprintf(fp, "end Phy.DiscoveredFracton\n");
    fclose(fp);
    return true;
}

static inline bool phy_emit_lean_biquandle8_cert(const phy_biquandle8_analysis_t *bq,
                                                 const char *path) {
    if (bq == NULL || path == NULL || bq->base.is_bijective == 0U ||
        bq->base.yb8_violations != 0U || bq->base.is_nondeg == 0U || bq->is_r1 == 0U ||
        bq->base.is_involutive != 0U || bq->base.is_indecomposable == 0U ||
        bq->is_irretract == 0U || bq->base.is_lyubashenko != 0U || bq->base.anf_deg < 2U) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }

    fprintf(fp, "import Phy.Noether\n\n");
    fprintf(fp, "namespace Phy.DiscoveredBiquandle8\n\n");
    fprintf(fp, "open Phy.Noether\n\n");
    fprintf(fp, "-- Merkle: 0x%016llx, order=%u, deg=%u, R1=1, irretract=1, T7=%u\n\n",
            (unsigned long long)bq->base.merkle, (unsigned)bq->base.order,
            (unsigned)bq->base.anf_deg, (unsigned)bq->col_torus7);

    phy_emit_bv6_lut_tree(fp, "ruleBQ8", bq->base.r);
    phy_emit_bv6_lut_tree(fp, "ruleBQ8Inv", bq->base.rinv);
    phy_emit_bv3_pair_tree(fp, "sigInv8", "a", "s", bq->base.sig_inv);
    phy_emit_bv3_pair_tree(fp, "tauInv8", "b", "t", bq->base.tau_inv);

    fprintf(fp, "theorem ruleBQ8_inv_left (u : BitVec 6) : ruleBQ8Inv (ruleBQ8 u) = u := by\n");
    fprintf(fp, "  simp only [ruleBQ8, ruleBQ8Inv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleBQ8_inv_right (u : BitVec 6) : ruleBQ8 (ruleBQ8Inv u) = u := by\n");
    fprintf(fp, "  simp only [ruleBQ8, ruleBQ8Inv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleBQ8_bijective : Bijective ruleBQ8 := by\n");
    fprintf(fp, "  refine ⟨fun a b h => ?_, fun b => ⟨ruleBQ8Inv b, ruleBQ8_inv_right b⟩⟩\n");
    fprintf(fp, "  have h1 := congrArg ruleBQ8Inv h\n");
    fprintf(fp, "  simpa [ruleBQ8_inv_left] using h1\n\n");

    fprintf(fp, "theorem ruleBQ8_satisfies_ybe : IsSetTheoreticYB8 ruleBQ8 := by\n");
    fprintf(fp, "  intro t\n");
    fprintf(fp, "  simp only [yb8Apply12, yb8Apply23, pack3, ruleBQ8]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleBQ8_nondegenerate : IsNonDegenerate8 ruleBQ8 :=\n");
    fprintf(fp, "  non_degenerate8_of_inverses ruleBQ8 sigInv8 tauInv8\n");
    fprintf(fp, "    (by intro a b; simp only [hi3, pack3, ruleBQ8, sigInv8]; bv_decide)\n");
    fprintf(fp, "    (by intro a b; simp only [hi3, pack3, ruleBQ8, sigInv8]; bv_decide)\n");
    fprintf(fp, "    (by intro a b; simp only [lo3, pack3, ruleBQ8, tauInv8]; bv_decide)\n");
    fprintf(fp, "    (by intro a b; simp only [lo3, pack3, ruleBQ8, tauInv8]; bv_decide)\n\n");

    fprintf(fp, "theorem ruleBQ8_r1_biquandle : IsR1Biquandle8 ruleBQ8 := by\n");
    fprintf(fp, "  intro x\n  simp only [pack3, ruleBQ8]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleBQ8_noninvolutive : ¬ IsInvolutive8 ruleBQ8 :=\n");
    fprintf(fp, "  non_involutive8_of_witness %u#6 (by decide)\n\n", (unsigned)bq->base.noninvol_u);

    fprintf(fp, "theorem ruleBQ8_indecomposable : ¬ IsDecomposable8 ruleBQ8 :=\n");
    fprintf(fp, "  indecomposable8_of_connected8 ruleBQ8 (by decide)\n\n");

    fprintf(fp, "theorem ruleBQ8_irretractable : ¬ IsRetractable8 ruleBQ8 :=\n");
    fprintf(fp, "  irretractable8_of_bool ruleBQ8 (by decide)\n\n");

    if (bq->base.nly_axis == 0U) {
        fprintf(fp, "theorem ruleBQ8_nonlyubashenko : ¬ IsLyubashenko8 ruleBQ8 :=\n");
        fprintf(fp, "  non_lyubashenko8_of_hi_cross %u#3 %u#3 %u#3 (by decide)\n\n",
                (unsigned)bq->base.nly_a, (unsigned)bq->base.nly_b0, (unsigned)bq->base.nly_b1);
    } else {
        fprintf(fp, "theorem ruleBQ8_nonlyubashenko : ¬ IsLyubashenko8 ruleBQ8 :=\n");
        fprintf(fp, "  non_lyubashenko8_of_lo_cross %u#3 %u#3 %u#3 (by decide)\n\n",
                (unsigned)bq->base.nly_a, (unsigned)bq->base.nly_b0, (unsigned)bq->base.nly_b1);
    }

    fprintf(fp,
            "theorem ruleBQ8_nonlinear : ¬ IsAffineGF2_6 (fun x => BitVec.getLsbD (ruleBQ8 x) "
            "%u) :=\n",
            (unsigned)bq->base.nl_bit);
    fprintf(fp, "  gf2_6_second_diff_witnesses_nonlinear %u#6 %u#6 (by decide)\n\n",
            (unsigned)bq->base.nl_a, (unsigned)bq->base.nl_b);

    fprintf(fp,
            "theorem ruleBQ8_knot_separation : colTorus2K ruleBQ8 1 = %u ∧ colTorus2K ruleBQ8 3 = "
            "%u ∧ colFig8 ruleBQ8 ruleBQ8Inv = %u ∧ colTorus2K ruleBQ8 7 = %u := by\n",
            (unsigned)bq->col_unknot, (unsigned)bq->col_trefoil, (unsigned)bq->col_fig8,
            (unsigned)bq->col_torus7);
    fprintf(fp, "  decide\n\n");

    fprintf(fp, "end Phy.DiscoveredBiquandle8\n");
    fclose(fp);
    return true;
}

static inline bool phy_emit_lean_tetra4_cert(const phy_tetra4_analysis_t *t4, const char *path) {
    if (t4 == NULL || path == NULL || t4->is_bijective == 0U || t4->tetra4_viol != 0U ||
        t4->is_involutive != 0U || t4->is_3coupled == 0U || t4->anf_deg < 3U) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }

    fprintf(fp, "import Phy.Noether\n\n");
    fprintf(fp, "namespace Phy.DiscoveredTetra4\n\n");
    fprintf(fp, "open Phy.Noether\n\n");
    fprintf(fp, "-- Merkle: 0x%016llx, order=%u, deg=%u, tetra4_viol=0/4096, tetra2=(%u,%u,%u)\n\n",
            (unsigned long long)t4->merkle, (unsigned)t4->order, (unsigned)t4->anf_deg,
            (unsigned)t4->tetra2_total, (unsigned)t4->tetra2_nl, (unsigned)t4->tetra2_noninvol_nl);

    phy_emit_bv6_lut_tree(fp, "ruleT4", t4->t);
    phy_emit_bv6_lut_tree(fp, "ruleT4Inv", t4->tinv);

    fprintf(fp, "theorem ruleT4_inv_left (u : BitVec 6) : ruleT4Inv (ruleT4 u) = u := by\n");
    fprintf(fp, "  simp only [ruleT4, ruleT4Inv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleT4_inv_right (u : BitVec 6) : ruleT4 (ruleT4Inv u) = u := by\n");
    fprintf(fp, "  simp only [ruleT4, ruleT4Inv]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleT4_bijective : Bijective ruleT4 := by\n");
    fprintf(fp, "  refine ⟨fun a b h => ?_, fun b => ⟨ruleT4Inv b, ruleT4_inv_right b⟩⟩\n");
    fprintf(fp, "  have h1 := congrArg ruleT4Inv h\n");
    fprintf(fp, "  simpa [ruleT4_inv_left] using h1\n\n");

    fprintf(fp, "theorem ruleT4_satisfies_tetra4 : IsSetTheoreticTetra4 ruleT4 := by\n");
    fprintf(fp, "  intro s\n");
    fprintf(fp, "  simp only [tetra4Apply, get2_12, set2_12, ruleT4]\n  bv_decide\n\n");

    fprintf(fp, "theorem ruleT4_noninvolutive : ¬ IsInvolutive8 ruleT4 :=\n");
    fprintf(fp, "  non_involutive8_of_witness %u#6 (by decide)\n\n", (unsigned)t4->noninvol_u);

    fprintf(fp, "theorem ruleT4_3coupled : Is3ChannelCoupled4 ruleT4 :=\n");
    fprintf(fp,
            "  coupled3_of_witnesses ruleT4 %u#2 %u#2 %u#2 %u#2 (by decide) %u#2 %u#2 %u#2 %u#2 "
            "(by decide)\n\n",
            (unsigned)t4->c0_a0, (unsigned)t4->c0_a1, (unsigned)t4->c0_b, (unsigned)t4->c0_c,
            (unsigned)t4->c2_a, (unsigned)t4->c2_b, (unsigned)t4->c2_c0, (unsigned)t4->c2_c1);

    fprintf(fp,
            "theorem ruleT4_nonlinear : ¬ IsAffineGF2_6 (fun x => BitVec.getLsbD (ruleT4 x) %u) "
            ":=\n",
            (unsigned)t4->nl_bit);
    fprintf(fp, "  gf2_6_second_diff_witnesses_nonlinear %u#6 %u#6 (by decide)\n\n",
            (unsigned)t4->nl_a, (unsigned)t4->nl_b);

    fprintf(fp, "end Phy.DiscoveredTetra4\n");
    fclose(fp);
    return true;
}

static inline void phy_emit_json_fir_tail(FILE *fp, const phy_fir_cert_ws_t *fir) {
    if (fir != NULL) {
        fprintf(fp, ",\n  \"fir_and_count\": %u,\n", (unsigned)fir->and_count);
        fprintf(fp, "  \"fir_mult_depth\": %u,\n", (unsigned)fir->mult_depth);
        fprintf(fp, "  \"fir_wire_bytes\": %zu,\n", fir->wire_len);
        fprintf(fp, "  \"fir_eval_ok\": %s", fir->ir_eval_ok ? "true" : "false");
    }
    fprintf(fp, "\n}\n");
}

static inline bool phy_emit_dual4_json(const phy_dual4_analysis_t *d4, const phy_fir_cert_ws_t *fir,
                                       const char *path) {
    if (d4 == NULL || path == NULL) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "{\n  \"p\": [");
    for (uint32_t i = 0U; i < 16U; ++i) {
        fprintf(fp, "%s%u", (i > 0U) ? ", " : "", (unsigned)d4->p[i]);
    }
    fprintf(fp, "],\n  \"m_plus\": [");
    for (uint32_t i = 0U; i < 4U; ++i) {
        fprintf(fp, "%s[%u, %u, %u, %u]", (i > 0U) ? ", " : "", (unsigned)d4->m_plus[i][0],
                (unsigned)d4->m_plus[i][1], (unsigned)d4->m_plus[i][2], (unsigned)d4->m_plus[i][3]);
    }
    fprintf(fp, "],\n  \"m_minus\": [");
    for (uint32_t i = 0U; i < 4U; ++i) {
        fprintf(fp, "%s[%u, %u, %u, %u]", (i > 0U) ? ", " : "", (unsigned)d4->m_minus[i][0],
                (unsigned)d4->m_minus[i][1], (unsigned)d4->m_minus[i][2],
                (unsigned)d4->m_minus[i][3]);
    }
    fprintf(fp, "],\n  \"lut_packed\": \"0x%016llx\",\n", (unsigned long long)d4->lut_packed);
    fprintf(fp, "  \"merkle\": \"0x%016llx\",\n", (unsigned long long)d4->merkle);
    fprintf(fp, "  \"flags\": %u,\n", (unsigned)d4->flags);
    fprintf(fp, "  \"anf_deg\": %u,\n", (unsigned)d4->anf_deg);
    fprintf(fp, "  \"is_dual_unitary\": %u,\n", (unsigned)d4->is_dual_unitary);
    fprintf(fp, "  \"is_mixing\": %u,\n", (unsigned)d4->is_mixing);
    fprintf(fp, "  \"is_non_bernoulli\": %u,\n", (unsigned)d4->is_non_bernoulli);
    fprintf(fp, "  \"rank_ent\": %u", (unsigned)d4->rank_ent);
    phy_emit_json_fir_tail(fp, fir);
    fclose(fp);
    return true;
}

static inline bool phy_emit_ame8_json(const phy_ame8_analysis_t *a8, const phy_fir_cert_ws_t *fir,
                                      const char *path) {
    if (a8 == NULL || path == NULL) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "{\n  \"r\": [");
    for (uint32_t i = 0U; i < 64U; ++i) {
        fprintf(fp, "%s%u", (i > 0U) ? ", " : "", (unsigned)a8->r[i]);
    }
    fprintf(fp, "],\n  \"merkle\": \"0x%016llx\",\n", (unsigned long long)a8->merkle);
    fprintf(fp, "  \"flags\": %u,\n", (unsigned)a8->flags);
    fprintf(fp, "  \"order\": %u,\n", (unsigned)a8->order);
    fprintf(fp, "  \"anf_deg\": %u,\n", (unsigned)a8->anf_deg);
    fprintf(fp, "  \"is_2unitary\": %u", (unsigned)a8->is_2unitary);
    phy_emit_json_fir_tail(fp, fir);
    fclose(fp);
    return true;
}

static inline bool phy_emit_fracton_json(const phy_fracton_analysis_t *fr,
                                         const phy_fir_cert_ws_t *fir, const char *path) {
    if (fr == NULL || path == NULL) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "{\n  \"t\": [");
    for (uint32_t i = 0U; i < 64U; ++i) {
        fprintf(fp, "%s%u", (i > 0U) ? ", " : "", (unsigned)fr->t[i]);
    }
    fprintf(fp, "],\n  \"merkle\": \"0x%016llx\",\n", (unsigned long long)fr->merkle);
    fprintf(fp, "  \"flags\": %u,\n", (unsigned)fr->flags);
    fprintf(fp, "  \"anf_deg\": %u,\n", (unsigned)fr->anf_deg);
    fprintf(fp, "  \"active_entries\": %u,\n", (unsigned)fr->active_entries);
    fprintf(fp, "  \"hsf_target_q\": %u,\n", (unsigned)fr->hsf_target_q);
    fprintf(fp, "  \"hsf_target_d\": %u,\n", (unsigned)fr->hsf_target_d);
    fprintf(fp, "  \"hsf_sector_size\": %u,\n", (unsigned)fr->hsf_sector_size);
    fprintf(fp, "  \"hsf_num_comps\": %u,\n", (unsigned)fr->hsf_num_comps);
    fprintf(fp, "  \"hsf_largest_comp\": %u,\n", (unsigned)fr->hsf_largest_comp);
    fprintf(fp, "  \"hsf_frozen\": %u,\n", (unsigned)fr->hsf_frozen);
    fprintf(fp, "  \"iso_span\": %u,\n", (unsigned)fr->iso_span);
    fprintf(fp, "  \"dip_span\": %u,\n", (unsigned)fr->dip_span);
    fprintf(fp, "  \"max_q_drift\": %lld,\n", (long long)fr->max_q_drift);
    fprintf(fp, "  \"max_d_drift\": %lld", (long long)fr->max_d_drift);
    phy_emit_json_fir_tail(fp, fir);
    fclose(fp);
    return true;
}

static inline bool phy_emit_biquandle8_json(const phy_biquandle8_analysis_t *bq,
                                            const phy_fir_cert_ws_t *fir, const char *path) {
    if (bq == NULL || path == NULL) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "{\n  \"r\": [");
    for (uint32_t i = 0U; i < 64U; ++i) {
        fprintf(fp, "%s%u", (i > 0U) ? ", " : "", (unsigned)bq->base.r[i]);
    }
    fprintf(fp, "],\n  \"merkle\": \"0x%016llx\",\n", (unsigned long long)bq->base.merkle);
    fprintf(fp, "  \"flags\": %u,\n", (unsigned)bq->base.flags);
    fprintf(fp, "  \"order\": %u,\n", (unsigned)bq->base.order);
    fprintf(fp, "  \"anf_deg\": %u,\n", (unsigned)bq->base.anf_deg);
    fprintf(fp, "  \"is_r1\": %u,\n", (unsigned)bq->is_r1);
    fprintf(fp, "  \"retract_classes\": %u,\n", (unsigned)bq->retract_classes);
    fprintf(fp, "  \"is_irretract\": %u,\n", (unsigned)bq->is_irretract);
    fprintf(fp, "  \"col_unknot\": %u,\n", (unsigned)bq->col_unknot);
    fprintf(fp, "  \"col_hopf\": %u,\n", (unsigned)bq->col_hopf);
    fprintf(fp, "  \"col_trefoil\": %u,\n", (unsigned)bq->col_trefoil);
    fprintf(fp, "  \"col_fig8\": %u,\n", (unsigned)bq->col_fig8);
    fprintf(fp, "  \"col_torus5\": %u,\n", (unsigned)bq->col_torus5);
    fprintf(fp, "  \"col_torus7\": %u", (unsigned)bq->col_torus7);
    phy_emit_json_fir_tail(fp, fir);
    fclose(fp);
    return true;
}

static inline bool phy_emit_tetra4_json(const phy_tetra4_analysis_t *t4,
                                        const phy_fir_cert_ws_t *fir, const char *path) {
    if (t4 == NULL || path == NULL) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "{\n  \"t\": [");
    for (uint32_t i = 0U; i < 64U; ++i) {
        fprintf(fp, "%s%u", (i > 0U) ? ", " : "", (unsigned)t4->t[i]);
    }
    fprintf(fp, "],\n  \"merkle\": \"0x%016llx\",\n", (unsigned long long)t4->merkle);
    fprintf(fp, "  \"flags\": %u,\n", (unsigned)t4->flags);
    fprintf(fp, "  \"order\": %u,\n", (unsigned)t4->order);
    fprintf(fp, "  \"anf_deg\": %u,\n", (unsigned)t4->anf_deg);
    fprintf(fp, "  \"tetra4_viol\": %u,\n", (unsigned)t4->tetra4_viol);
    fprintf(fp, "  \"tetra2_total\": %u,\n", (unsigned)t4->tetra2_total);
    fprintf(fp, "  \"tetra2_nl\": %u,\n", (unsigned)t4->tetra2_nl);
    fprintf(fp, "  \"tetra2_noninvol_nl\": %u,\n", (unsigned)t4->tetra2_noninvol_nl);
    fprintf(fp, "  \"is_involutive\": %u,\n", (unsigned)t4->is_involutive);
    fprintf(fp, "  \"is_3coupled\": %u", (unsigned)t4->is_3coupled);
    phy_emit_json_fir_tail(fp, fir);
    fclose(fp);
    return true;
}

/* --- PHY v4.0: Lean 4 & JSON Certificate Emitters --- */

static inline bool phy_emit_lean_ame4_nogo_cert(const phy_ame4_nogo_t *ng, const char *path) {
    if (ng == NULL || path == NULL || ng->total_ols != UINT32_C(6912) ||
        ng->nonlin_ols != UINT32_C(0) || ng->affine_ols != UINT32_C(6912)) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "import Phy.Noether\n\n");
    fprintf(fp, "namespace Phy.DiscoveredAME4NoGo\n\n");
    fprintf(fp, "open Phy.Noether\n\n");
    fprintf(fp, "-- Merkle: 0x%016llx, latin4=%u, total_ols=%u, affine_ols=%u, nonlin_ols=%u\n\n",
            (unsigned long long)ng->merkle, (unsigned)ng->latin4_count, (unsigned)ng->total_ols,
            (unsigned)ng->affine_ols, (unsigned)ng->nonlin_ols);
    fprintf(fp, "def sampleAME4Lut : BitVec 64 := 0x%016llx#64\n\n",
            (unsigned long long)ng->sample_lut);
    fprintf(fp, "theorem universal_ame4_affine_rigidity (lut : BitVec 64)\n"
                "    (h : is2Unitary4Lut lut = true) : isAffine4Lut lut = true :=\n"
                "  ame4_affine_rigidity_nogo lut h\n\n");
    fprintf(fp, "theorem sampleAME4_is_2unitary : is2Unitary4Lut sampleAME4Lut = true := by\n");
    fprintf(fp, "  decide\n\n");
    fprintf(fp, "theorem sampleAME4_forced_affine : isAffine4Lut sampleAME4Lut = true :=\n");
    fprintf(fp, "  ame4_affine_rigidity_nogo sampleAME4Lut sampleAME4_is_2unitary\n\n");
    fprintf(fp, "end Phy.DiscoveredAME4NoGo\n");
    fclose(fp);
    return true;
}

static inline bool phy_emit_ame4_nogo_json(const phy_ame4_nogo_t *ng, const char *path) {
    if (ng == NULL || path == NULL) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "{\n  \"sample_p\": [");
    for (uint32_t i = 0U; i < 16U; ++i) {
        fprintf(fp, "%s%u", (i > 0U) ? ", " : "", (unsigned)ng->sample_p[i]);
    }
    fprintf(fp, "],\n  \"sample_lut\": \"0x%016llx\",\n", (unsigned long long)ng->sample_lut);
    fprintf(fp, "  \"merkle\": \"0x%016llx\",\n", (unsigned long long)ng->merkle);
    fprintf(fp, "  \"latin4_count\": %u,\n", (unsigned)ng->latin4_count);
    fprintf(fp, "  \"total_ols\": %u,\n", (unsigned)ng->total_ols);
    fprintf(fp, "  \"affine_ols\": %u,\n", (unsigned)ng->affine_ols);
    fprintf(fp, "  \"nonlin_ols\": %u,\n", (unsigned)ng->nonlin_ols);
    fprintf(fp, "  \"bernoulli_ols\": %u\n}\n", (unsigned)ng->bernoulli_ols);
    fclose(fp);
    return true;
}

static inline bool phy_emit_lean_simplex4_cert(const phy_simplex4_4_t *s4, const char *path) {
    if (s4 == NULL || path == NULL || s4->is_bijective == 0U || s4->viol4 != UINT32_C(0) ||
        s4->viol4_coupled != UINT32_C(0) || s4->is_involutive != 0U || s4->is_4coupled == 0U ||
        s4->anf_deg < 3U) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "import Phy.Noether\n\n");
    fprintf(fp, "namespace Phy.DiscoveredSimplex4\n\n");
    fprintf(fp, "open Phy.Noether\n\n");
    fprintf(fp,
            "-- Merkle: 0x%016llx, viol4=%u/1048576, order=%u, deg=%u, viol2_ccx=%u, "
            "viol2_cccx=%u\n\n",
            (unsigned long long)s4->merkle, (unsigned)s4->viol4, (unsigned)s4->order4,
            (unsigned)s4->anf_deg, (unsigned)s4->viol2_ccx, (unsigned)s4->viol2_cccx);

    int32_t s_tmp[256];
    for (uint32_t u = 0U; u < 256U; ++u) {
        s_tmp[u] = (int32_t)s4->s[u];
    }
    fprintf(fp, "def ruleS4 (u : BitVec 8) : BitVec 8 :=\n  ");
    phy_emit_bv_tree_rec(fp, "u", 8U, 8U, 0U, 256U, s_tmp, 2U);
    fprintf(fp, "\n\n");

    fprintf(fp, "theorem ccx02_4simplex : IsSetTheoretic4Simplex2 ccx02_4 :=\n");
    fprintf(fp, "  ccx02_satisfies_4simplex2\n\n");
    fprintf(fp, "theorem ccx13_4simplex : IsSetTheoretic4Simplex2 ccx13_4 :=\n");
    fprintf(fp, "  ccx13_satisfies_4simplex2\n\n");
    fprintf(fp, "theorem cccx1_not_4simplex : ¬ IsSetTheoretic4Simplex2 cccx1_4 :=\n");
    fprintf(fp, "  cccx1_fails_4simplex2\n\n");

    fprintf(fp, "private theorem ruleS4_eq_ring (u : BitVec 8) : ruleS4 u = s4RingLift u := by\n");
    fprintf(fp, "  revert u\n  decide\n\n");
    fprintf(fp, "theorem ruleS4_fn_eq : ruleS4 = s4RingLift := funext ruleS4_eq_ring\n\n");
    fprintf(fp, "theorem ruleS4_bijective : Bijective ruleS4 := by\n");
    fprintf(fp, "  simpa only [ruleS4_fn_eq] using s4RingLift_bijective\n\n");
    fprintf(fp, "theorem ruleS4_satisfies_4simplex : IsSetTheoretic4Simplex4 ruleS4 := by\n");
    fprintf(fp, "  simpa only [ruleS4_fn_eq] using s4RingLift_satisfies_4simplex4\n\n");
    fprintf(fp, "theorem ruleS4_noninvolutive : ¬ (∀ u : BitVec 8, ruleS4 (ruleS4 u) = u) := by\n");
    fprintf(fp, "  simpa only [ruleS4_fn_eq] using s4RingLift_noninvolutive\n\n");
    fprintf(fp,
            "theorem coupled4ch_satisfies_4simplex : IsSetTheoretic4Simplex4 s4Coupled4Ch :=\n");
    fprintf(fp, "  s4Coupled4Ch_satisfies_4simplex4\n\n");
    fprintf(fp, "theorem coupled4ch_4coupled : Is4ChannelCoupled4 s4Coupled4Ch :=\n");
    fprintf(fp, "  s4Coupled4Ch_4channel_coupled\n\n");

    fprintf(fp, "end Phy.DiscoveredSimplex4\n");
    fclose(fp);
    return true;
}

static inline bool phy_emit_simplex4_json(const phy_simplex4_4_t *s4, const phy_fir_cert_ws_t *fir,
                                          const char *path) {
    if (s4 == NULL || path == NULL) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "{\n  \"s\": [");
    for (uint32_t u = 0U; u < 256U; ++u) {
        fprintf(fp, "%s%u", (u > 0U) ? ", " : "", (unsigned)s4->s[u]);
    }
    fprintf(fp, "],\n  \"s_coupled\": [");
    for (uint32_t u = 0U; u < 256U; ++u) {
        fprintf(fp, "%s%u", (u > 0U) ? ", " : "", (unsigned)s4->s_coupled[u]);
    }
    fprintf(fp, "],\n  \"ccx02\": [");
    for (uint32_t u = 0U; u < 16U; ++u) {
        fprintf(fp, "%s%u", (u > 0U) ? ", " : "", (unsigned)s4->ccx02[u]);
    }
    fprintf(fp, "],\n  \"merkle\": \"0x%016llx\",\n", (unsigned long long)s4->merkle);
    fprintf(fp, "  \"flags\": %u,\n", (unsigned)s4->flags);
    fprintf(fp, "  \"viol4\": %u,\n", (unsigned)s4->viol4);
    fprintf(fp, "  \"viol4_coupled\": %u,\n", (unsigned)s4->viol4_coupled);
    fprintf(fp, "  \"order4\": %u,\n", (unsigned)s4->order4);
    fprintf(fp, "  \"viol2_ccx\": %u,\n", (unsigned)s4->viol2_ccx);
    fprintf(fp, "  \"viol2_cccx\": %u,\n", (unsigned)s4->viol2_cccx);
    fprintf(fp, "  \"anf_deg\": %u,\n", (unsigned)s4->anf_deg);
    fprintf(fp, "  \"is_bijective\": %u,\n", (unsigned)s4->is_bijective);
    fprintf(fp, "  \"is_involutive\": %u,\n", (unsigned)s4->is_involutive);
    fprintf(fp, "  \"is_4coupled\": %u", (unsigned)s4->is_4coupled);
    phy_emit_json_fir_tail(fp, fir);
    fclose(fp);
    return true;
}

static inline bool phy_emit_lean_chiral6_cert(const phy_chiral6_analysis_t *ch, const char *path) {
    if (ch == NULL || path == NULL || ch->cocycle_viol != 0U || ch->is_chiral == 0U ||
        ch->is_non_coboundary == 0U) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "import Phy.Noether\n\n");
    fprintf(fp, "namespace Phy.DiscoveredChiral6\n\n");
    fprintf(fp, "open Phy.Noether\n\n");
    fprintf(
        fp,
        "-- Merkle: 0x%016llx, Q_6 Z_4 2-Cocycle: L(3_1)=[%u,%u,%u,%u], R(3_1)=[%u,%u,%u,%u]\n\n",
        (unsigned long long)ch->merkle, (unsigned)ch->tref_l_w[0], (unsigned)ch->tref_l_w[1],
        (unsigned)ch->tref_l_w[2], (unsigned)ch->tref_l_w[3], (unsigned)ch->tref_r_w[0],
        (unsigned)ch->tref_r_w[1], (unsigned)ch->tref_r_w[2], (unsigned)ch->tref_r_w[3]);

    fprintf(fp, "def discPhi (p : Fin 6 × Fin 6) : BitVec 2 :=\n  match p.1.val, p.2.val with\n");
    for (uint32_t a = 0U; a < 6U; ++a) {
        fprintf(fp, " ");
        for (uint32_t b = 0U; b < 6U; ++b) {
            fprintf(fp, " | %u, %u => %u#2", a, b, (unsigned)(ch->phi[a * 6U + b] & 3U));
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "  | _, _ => 0#2\n\n");

    fprintf(fp, "private theorem fin6_split (x : Fin 6) :\n"
                "    x = 0 ∨ x = 1 ∨ x = 2 ∨ x = 3 ∨ x = 4 ∨ x = 5 := by\n"
                "  rcases x with ⟨v, hv⟩\n"
                "  match v, hv with\n"
                "  | 0, _ => exact Or.inl rfl\n"
                "  | 1, _ => exact Or.inr (Or.inl rfl)\n"
                "  | 2, _ => exact Or.inr (Or.inr (Or.inl rfl))\n"
                "  | 3, _ => exact Or.inr (Or.inr (Or.inr (Or.inl rfl)))\n"
                "  | 4, _ => exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inl rfl))))\n"
                "  | 5, _ => exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inr rfl))))\n\n");

    fprintf(fp, "theorem q6_quandle_ybe : IsSetTheoreticYB6 q6R := q6_satisfies_ybe\n");
    fprintf(fp, "theorem q6_quandle_r1 : IsR1_6 q6R := q6_r1\n\n");
    fprintf(fp, "theorem discPhi_is_2cocycle : Is2CocycleMod4_6 q6R discPhi := by\n");
    fprintf(fp, "  refine ⟨fun x => ?_, fun a b c => ?_⟩\n");
    fprintf(fp, "  · rcases fin6_split x with rfl | rfl | rfl | rfl | rfl | rfl <;> decide\n");
    fprintf(fp, "  · rcases fin6_split a with rfl | rfl | rfl | rfl | rfl | rfl <;>\n"
                "      rcases fin6_split b with rfl | rfl | rfl | rfl | rfl | rfl <;>\n"
                "      rcases fin6_split c with rfl | rfl | rfl | rfl | rfl | rfl <;> decide\n\n");

    fprintf(fp, "theorem discPhi_not_coboundary : ¬ Is2CoboundaryMod4_6 q6R discPhi :=\n");
    fprintf(fp, "  q6_non_coboundary_of_witness q6R discPhi 3 (%u, %u) (by decide) (by decide)\n\n",
            (unsigned)ch->wit_a, (unsigned)ch->wit_b);

    fprintf(fp,
            "theorem discPhi_chiral_separation :\n"
            "    cocycleCountTorus2Pos q6R discPhi 1 0#2 = %u ∧\n"
            "    cocycleCountTorus2Pos q6R discPhi 2 0#2 = %u ∧\n"
            "    cocycleCountTorus2Pos q6R discPhi 2 2#2 = %u ∧\n"
            "    cocycleCountTorus2Pos q6R discPhi 3 1#2 = %u ∧\n"
            "    cocycleCountTorus2Neg q6RInv discPhi 3 1#2 = %u ∧\n"
            "    cocycleCountTorus2Pos q6R discPhi 3 3#2 = %u ∧\n"
            "    cocycleCountTorus2Neg q6RInv discPhi 3 3#2 = %u ∧\n"
            "    cocycleCountFig8 q6R q6RInv discPhi 0#2 = %u ∧\n"
            "    cocycleCountTorus2Pos q6R discPhi 3 1#2 ≠ cocycleCountTorus2Neg q6RInv discPhi 3 "
            "1#2 := by\n"
            "  decide\n\n",
            (unsigned)ch->unknot_w[0], (unsigned)ch->hopf_w[0], (unsigned)ch->hopf_w[2],
            (unsigned)ch->tref_l_w[1], (unsigned)ch->tref_r_w[1], (unsigned)ch->tref_l_w[3],
            (unsigned)ch->tref_r_w[3], (unsigned)ch->fig8_w[0]);

    fprintf(fp, "end Phy.DiscoveredChiral6\n");
    fclose(fp);
    return true;
}

static inline bool phy_emit_chiral6_json(const phy_chiral6_analysis_t *ch, const char *path) {
    if (ch == NULL || path == NULL) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "{\n  \"r\": [");
    for (uint32_t i = 0U; i < 36U; ++i) {
        fprintf(fp, "%s%u", (i > 0U) ? ", " : "", (unsigned)ch->r[i]);
    }
    fprintf(fp, "],\n  \"phi\": [");
    for (uint32_t i = 0U; i < 36U; ++i) {
        fprintf(fp, "%s%u", (i > 0U) ? ", " : "", (unsigned)ch->phi[i]);
    }
    fprintf(fp, "],\n  \"merkle\": \"0x%016llx\",\n", (unsigned long long)ch->merkle);
    fprintf(fp, "  \"flags\": %u,\n", (unsigned)ch->flags);
    fprintf(fp, "  \"cocycle_rank\": %u,\n", (unsigned)ch->cocycle_rank);
    fprintf(fp, "  \"cocycle_viol\": %u,\n", (unsigned)ch->cocycle_viol);
    fprintf(fp, "  \"unknot_w\": [%u, %u, %u, %u],\n", (unsigned)ch->unknot_w[0],
            (unsigned)ch->unknot_w[1], (unsigned)ch->unknot_w[2], (unsigned)ch->unknot_w[3]);
    fprintf(fp, "  \"hopf_w\": [%u, %u, %u, %u],\n", (unsigned)ch->hopf_w[0],
            (unsigned)ch->hopf_w[1], (unsigned)ch->hopf_w[2], (unsigned)ch->hopf_w[3]);
    fprintf(fp, "  \"tref_l_w\": [%u, %u, %u, %u],\n", (unsigned)ch->tref_l_w[0],
            (unsigned)ch->tref_l_w[1], (unsigned)ch->tref_l_w[2], (unsigned)ch->tref_l_w[3]);
    fprintf(fp, "  \"tref_r_w\": [%u, %u, %u, %u],\n", (unsigned)ch->tref_r_w[0],
            (unsigned)ch->tref_r_w[1], (unsigned)ch->tref_r_w[2], (unsigned)ch->tref_r_w[3]);
    fprintf(fp, "  \"fig8_w\": [%u, %u, %u, %u],\n", (unsigned)ch->fig8_w[0],
            (unsigned)ch->fig8_w[1], (unsigned)ch->fig8_w[2], (unsigned)ch->fig8_w[3]);
    fprintf(fp, "  \"is_chiral\": %u,\n", (unsigned)ch->is_chiral);
    fprintf(fp, "  \"is_non_coboundary\": %u\n}\n", (unsigned)ch->is_non_coboundary);
    fclose(fp);
    return true;
}

static inline bool phy_emit_lean_qscar_cert(const phy_qscar_analysis_t *qs, const char *path) {
    if (qs == NULL || path == NULL || qs->du_unitary == 0U || qs->du_dual == 0U ||
        qs->du_full_sup == 0U || qs->scar_period != 2U) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "import Phy.Noether\n\n");
    fprintf(fp, "namespace Phy.DiscoveredQScar\n\n");
    fprintf(fp, "open Phy.Noether\n\n");
    fprintf(fp, "-- Merkle: 0x%016llx, Z[i] Dual-Unitary + Scar DTC period=%u\n\n",
            (unsigned long long)qs->merkle, (unsigned)qs->scar_period);

    fprintf(fp, "def discW (i j : Fin 4) : GaussInt :=\n  match i.val, j.val with\n");
    for (uint32_t i = 0U; i < 4U; ++i) {
        fprintf(fp, " ");
        for (uint32_t j = 0U; j < 4U; ++j) {
            fprintf(fp, " | %u, %u => ⟨%d, %d⟩", i, j, (int)qs->w_re[i][j], (int)qs->w_im[i][j]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "  | _, _ => gZero\n\n");

    fprintf(fp, "theorem discW_unitary (j k : Fin 4) :\n"
                "    gSumFin4 (fun i => gMul (gConj (discW i j)) (discW i k)) =\n"
                "      if j == k then ⟨4, 0⟩ else gZero := by\n"
                "  revert j k\n  decide\n\n");
    fprintf(fp, "theorem discW_dual_unitary (j k : Fin 4) :\n"
                "    gSumFin4 (fun i => gMul (gConj (wCxDUReshuffle i j)) (wCxDUReshuffle i k)) =\n"
                "      if j == k then ⟨4, 0⟩ else gZero :=\n"
                "  wCxDU_is_dual_unitary j k\n\n");
    fprintf(fp, "theorem discW_full_superposition (i j : Fin 4) : discW i j ≠ gZero := by\n");
    fprintf(fp, "  revert i j\n  decide\n\n");
    fprintf(fp, "theorem qscar_local_gate_unitary (j k : BitVec 4) :\n"
                "    gSumBV4 (fun i => gMul (gConj (scarGate4 i j)) (scarGate4 i k)) =\n"
                "      if j == k then ⟨4, 0⟩ else gZero :=\n"
                "  scarGate4_is_unitary j k\n\n");
    fprintf(fp, "theorem qscar_dtc_neel_protected {L : Nat} (hGe : 4 ≤ L) (hEven : L %% 2 = 0)\n"
                "    (i : Fin L) :\n"
                "    neel0 L i ≠ neel0 L ⟨(i.val + 3) %% L, Nat.mod_lt _ (by omega)⟩ ∧\n"
                "    neel1 L i ≠ neel1 L ⟨(i.val + 3) %% L, Nat.mod_lt _ (by omega)⟩ :=\n"
                "  dtc_neel_local_protected hGe hEven i\n\n");
    fprintf(fp, "theorem qscar_dtc_exact_period2 {L : Nat} (hGe : 4 ≤ L) :\n"
                "    flipAll (neel0 L) = neel1 L ∧ flipAll (neel1 L) = neel0 L ∧\n"
                "    neel0 L ≠ neel1 L ∧ flipAll (flipAll (neel0 L)) = neel0 L :=\n"
                "  dtc_neel_exact_period2 hGe\n\n");

    fprintf(fp, "end Phy.DiscoveredQScar\n");
    fclose(fp);
    return true;
}

static inline bool phy_emit_qscar_json(const phy_qscar_analysis_t *qs, const char *path) {
    if (qs == NULL || path == NULL) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "{\n  \"w_re\": [");
    for (uint32_t i = 0U; i < 4U; ++i) {
        fprintf(fp, "%s[%d, %d, %d, %d]", (i > 0U) ? ", " : "", (int)qs->w_re[i][0],
                (int)qs->w_re[i][1], (int)qs->w_re[i][2], (int)qs->w_re[i][3]);
    }
    fprintf(fp, "],\n  \"w_im\": [");
    for (uint32_t i = 0U; i < 4U; ++i) {
        fprintf(fp, "%s[%d, %d, %d, %d]", (i > 0U) ? ", " : "", (int)qs->w_im[i][0],
                (int)qs->w_im[i][1], (int)qs->w_im[i][2], (int)qs->w_im[i][3]);
    }
    fprintf(fp, "],\n  \"merkle\": \"0x%016llx\",\n", (unsigned long long)qs->merkle);
    fprintf(fp, "  \"flags\": %u,\n", (unsigned)qs->flags);
    fprintf(fp, "  \"du_unitary\": %u,\n", (unsigned)qs->du_unitary);
    fprintf(fp, "  \"du_dual\": %u,\n", (unsigned)qs->du_dual);
    fprintf(fp, "  \"du_full_sup\": %u,\n", (unsigned)qs->du_full_sup);
    fprintf(fp, "  \"scar_period\": %u,\n", (unsigned)qs->scar_period);
    fprintf(fp, "  \"du_sup\": [%u, %u, %u, %u],\n", (unsigned)qs->du_sup[0],
            (unsigned)qs->du_sup[1], (unsigned)qs->du_sup[2], (unsigned)qs->du_sup[3]);
    fprintf(fp, "  \"du_purity_t2_num\": %llu,\n", (unsigned long long)qs->du_purity_num[1]);
    fprintf(fp, "  \"du_purity_t2_den\": %llu,\n", (unsigned long long)qs->du_purity_den[1]);
    fprintf(fp, "  \"neel_sup\": [%u, %u, %u, %u],\n", (unsigned)qs->neel_sup[0],
            (unsigned)qs->neel_sup[1], (unsigned)qs->neel_sup[2], (unsigned)qs->neel_sup[3]);
    fprintf(fp, "  \"neel_ret_num\": [%llu, %llu, %llu, %llu],\n",
            (unsigned long long)qs->neel_ret_num[0], (unsigned long long)qs->neel_ret_num[1],
            (unsigned long long)qs->neel_ret_num[2], (unsigned long long)qs->neel_ret_num[3]);
    fprintf(fp, "  \"neel_ret_den\": [%llu, %llu, %llu, %llu],\n",
            (unsigned long long)qs->neel_ret_den[0], (unsigned long long)qs->neel_ret_den[1],
            (unsigned long long)qs->neel_ret_den[2], (unsigned long long)qs->neel_ret_den[3]);
    fprintf(fp, "  \"therm_sup\": [%u, %u, %u, %u],\n", (unsigned)qs->therm_sup[0],
            (unsigned)qs->therm_sup[1], (unsigned)qs->therm_sup[2], (unsigned)qs->therm_sup[3]);
    fprintf(fp, "  \"therm_purity_t4_num\": %llu,\n", (unsigned long long)qs->therm_purity_num[3]);
    fprintf(fp, "  \"therm_purity_t4_den\": %llu\n}\n",
            (unsigned long long)qs->therm_purity_den[3]);
    fclose(fp);
    return true;
}

static inline bool phy_emit_lean_fracton_comp_cert(const phy_fracton_comp_t *fc, const char *path) {
    if (fc == NULL || path == NULL || fc->is_involutive == 0U || fc->laplacian_ok == 0U ||
        fc->ccx_ok == 0U || fc->scaling_strict == 0U) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "import Phy.Noether\n\n");
    fprintf(fp, "namespace Phy.DiscoveredFracComp\n\n");
    fprintf(fp, "open Phy.Noether\n\n");
    fprintf(fp, "-- Merkle: 0x%016llx, HSF L=(6,7,8) comps=(%u/%u, %u/%u, %u/%u)\n\n",
            (unsigned long long)fc->merkle, (unsigned)fc->comps[0], (unsigned)fc->sec[0],
            (unsigned)fc->comps[1], (unsigned)fc->sec[1], (unsigned)fc->comps[2],
            (unsigned)fc->sec[2]);

    phy_emit_bv6_lut_tree(fp, "ruleFracToff", fc->t_toff);
    phy_emit_bv6_lut_tree(fp, "ruleFracToffInv", fc->t_toff_inv);

    fprintf(fp,
            "theorem ruleFracToff_inv_left (u : BitVec 6) : ruleFracToffInv (ruleFracToff u) = u "
            ":= by\n"
            "  simp only [ruleFracToff, ruleFracToffInv]\n  bv_decide\n\n");
    fprintf(fp,
            "theorem ruleFracToff_inv_right (u : BitVec 6) : ruleFracToff (ruleFracToffInv u) = u "
            ":= by\n"
            "  simp only [ruleFracToff, ruleFracToffInv]\n  bv_decide\n\n");
    fprintf(fp, "theorem ruleFracToff_bijective : Bijective ruleFracToff := by\n");
    fprintf(fp, "  refine ⟨fun a b h => ?_, fun b => ⟨ruleFracToffInv b, ruleFracToff_inv_right "
                "b⟩⟩\n");
    fprintf(fp, "  have h1 := congrArg ruleFracToffInv h\n");
    fprintf(fp, "  simpa [ruleFracToff_inv_left] using h1\n\n");

    fprintf(fp,
            "theorem ruleFracToff_laplacian : SatisfiesDipoleLaplacian ruleFracToff rho2 := by\n"
            "  unfold SatisfiesDipoleLaplacian rho2\n  decide\n\n");
    fprintf(fp, "theorem ruleFracToff_universal_involution (u : BitVec 6) :\n"
                "    ruleFracToff (ruleFracToff u) = u :=\n"
                "  fracton3_universal_involution ruleFracToff ruleFracToff_bijective "
                "ruleFracToff_laplacian u\n\n");
    fprintf(fp, "theorem ruleFracToff_static_control_nogo (u : BitVec 6)\n"
                "    (h0 : site0_3 (ruleFracToff u) = site0_3 u)\n"
                "    (h2 : site2_3 (ruleFracToff u) = site2_3 u) : ruleFracToff u = u :=\n"
                "  dipole_static_control_nogo ruleFracToff ruleFracToff_laplacian u h0 h2\n\n");
    fprintf(fp, "theorem ruleFracToff_implements_ccx (ca cc tgt : Bool) :\n"
                "    decodeDipoleTarget (ruleFracToff (encodeDipoleLogical ca cc tgt)) =\n"
                "      xor tgt (ca && cc) := by\n"
                "  revert ca cc tgt\n  decide\n\n");
    fprintf(fp,
            "theorem ruleFracToff_single_charge_frozen {L : Nat} (pos : Nat) (hpos : pos + 2 < L)\n"
            "    (s : Fin L → BitVec 2) (hQ : openCharge L rho2 s ≤ 1) :\n"
            "    applyGate3 ruleFracToff pos hpos s = s :=\n"
            "  fracton_single_charge_frozen ruleFracToff ruleFracToff_laplacian pos hpos s hQ\n\n");
    fprintf(fp, "end Phy.DiscoveredFracComp\n");
    fclose(fp);
    return true;
}

static inline bool phy_emit_fracton_comp_json(const phy_fracton_comp_t *fc,
                                              const phy_fir_cert_ws_t *fir, const char *path) {
    if (fc == NULL || path == NULL) {
        return false;
    }
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        return false;
    }
    fprintf(fp, "{\n  \"t_sample\": [");
    for (uint32_t i = 0U; i < 64U; ++i) {
        fprintf(fp, "%s%u", (i > 0U) ? ", " : "", (unsigned)fc->t_sample[i]);
    }
    fprintf(fp, "],\n  \"t_toff\": [");
    for (uint32_t i = 0U; i < 64U; ++i) {
        fprintf(fp, "%s%u", (i > 0U) ? ", " : "", (unsigned)fc->t_toff[i]);
    }
    fprintf(fp, "],\n  \"merkle\": \"0x%016llx\",\n", (unsigned long long)fc->merkle);
    fprintf(fp, "  \"flags\": %u,\n", (unsigned)fc->flags);
    fprintf(fp, "  \"sec\": [%u, %u, %u],\n", (unsigned)fc->sec[0], (unsigned)fc->sec[1],
            (unsigned)fc->sec[2]);
    fprintf(fp, "  \"comps\": [%u, %u, %u],\n", (unsigned)fc->comps[0], (unsigned)fc->comps[1],
            (unsigned)fc->comps[2]);
    fprintf(fp, "  \"largest\": [%u, %u, %u],\n", (unsigned)fc->largest[0],
            (unsigned)fc->largest[1], (unsigned)fc->largest[2]);
    fprintf(fp, "  \"frozen\": [%u, %u, %u],\n", (unsigned)fc->frozen[0], (unsigned)fc->frozen[1],
            (unsigned)fc->frozen[2]);
    fprintf(fp, "  \"is_involutive\": %u,\n", (unsigned)fc->is_involutive);
    fprintf(fp, "  \"ccx_ok\": %u,\n", (unsigned)fc->ccx_ok);
    fprintf(fp, "  \"laplacian_ok\": %u,\n", (unsigned)fc->laplacian_ok);
    fprintf(fp, "  \"scaling_strict\": %u", (unsigned)fc->scaling_strict);
    phy_emit_json_fir_tail(fp, fir);
    fclose(fp);
    return true;
}

#endif /* PHY_CERT_H */
