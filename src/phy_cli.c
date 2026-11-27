#include "phy.h"
#include "phy_cert.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static phy_arena_u g_arena;
static phy_arena_u g_arena_rt;
static phy_q4_ws_t g_q4_ws;
static phy_fracton_hsf_ws_t g_frac_ws;
static phy_fracton_scale_ws_t g_frac_scale_ws;
static phy_fir_cert_ws_t g_fir_ws;
static flow_ir_prog_t g_eval_prog;
static uint8_t g_eval_buf[FLOW_IR_MAX_FILE_BYTES];
static uint8_t g_trans_out_buf[FLOW_IR_MAX_FILE_BYTES];
static uint8_t g_arena_wire[PHY_ARENA_PACK_MAX_BYTES];

static uint64_t phy_rotl64_2bit_left(uint64_t *lo, uint64_t *hi) {
    /* Rotate 64 2-bit sites toward lower index by 1 site (site s <- old site (s+1)%64) */
    uint64_t old_lo = *lo;
    uint64_t old_hi = *hi;
    uint64_t lo_s0 = old_lo & UINT64_C(3);
    uint64_t hi_s0 = old_hi & UINT64_C(3);
    *lo = (old_lo >> 2U) | (hi_s0 << 62U);
    *hi = (old_hi >> 2U) | (lo_s0 << 62U);
    return old_lo;
}

static void phy_rotl64_2bit_right(uint64_t *lo, uint64_t *hi) {
    uint64_t old_lo = *lo;
    uint64_t old_hi = *hi;
    uint64_t lo_s31 = (old_lo >> 62U) & UINT64_C(3);
    uint64_t hi_s31 = (old_hi >> 62U) & UINT64_C(3);
    *lo = (old_lo << 2U) | hi_s31;
    *hi = (old_hi << 2U) | lo_s31;
}

static int phy_selftest(bool fail_fast) {
    int passed = 0;
    const int total = 26;

    /* T1: M61 arithmetic, Stein GCD, binary division, FNV-1a */
    bool ok1 = true;
    uint64_t vals[3] = {UINT64_C(3), UINT64_C(0x123456789ABCDEF), PHY_M61 - UINT64_C(1)};
    for (int i = 0; i < 3; ++i) {
        uint64_t inv = phy_m61_inv(vals[i]);
        if (phy_m61_mul(vals[i], inv) != UINT64_C(1)) {
            ok1 = false;
        }
    }
    if (phy_m61_add(PHY_M61 - UINT64_C(1), UINT64_C(2)) != UINT64_C(1)) {
        ok1 = false;
    }
    if (phy_m61_sub(UINT64_C(1), UINT64_C(2)) != (PHY_M61 - UINT64_C(1))) {
        ok1 = false;
    }
    if (phy_stein_gcd64(UINT64_C(2310), UINT64_C(130)) != UINT64_C(10) ||
        phy_stein_gcd64(UINT64_C(0), UINT64_C(42)) != UINT64_C(42)) {
        ok1 = false;
    }
    uint64_t rem = UINT64_C(0);
    uint64_t quo = phy_div_u64(UINT64_C(123456789012345), UINT64_C(987654321), &rem);
    if (quo * UINT64_C(987654321) + rem != UINT64_C(123456789012345) ||
        rem >= UINT64_C(987654321)) {
        ok1 = false;
    }
    if (phy_fnv1a64_u64(UINT64_C(1), UINT64_C(2)) == phy_fnv1a64_u64(UINT64_C(1), UINT64_C(3))) {
        ok1 = false;
    }
    printf("T1 (M61/Stein/Div/FNV): %s\n", ok1 ? "PASS" : "FAIL");
    if (ok1) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T2: ID & SWAP separable/YB baselines */
    bool ok2 = true;
    uint8_t p_id[16];
    uint8_t p_sw[16];
    for (uint8_t u = 0U; u < 16U; ++u) {
        p_id[u] = u;
        p_sw[u] = (uint8_t)(((u & 3U) << 2U) | ((u >> 2U) & 3U));
    }
    phy_rule_analysis_t ra_id;
    phy_rule_analysis_t ra_sw;
    phy_analyze_rule(p_id, 1ULL, &ra_id);
    phy_analyze_rule(p_sw, 2ULL, &ra_sw);
    if (ra_id.rank_tot != 15U || ra_id.j_nonzero != 0U || ra_id.rank_ent != 9U ||
        ra_id.yb_violations != 0U || ra_id.is_id != 1U ||
        (ra_id.flags & PHY_RULE_FLAG_YANGBOX) != 0U) {
        ok2 = false;
    }
    if (ra_sw.rank_tot != 6U || ra_sw.rank_ent != 0U || ra_sw.yb_violations != 0U ||
        ra_sw.anf_deg != 1U || ra_sw.is_swap != 1U || (ra_sw.flags & PHY_RULE_FLAG_YANGBOX) != 0U) {
        ok2 = false;
    }
    printf("T2 (ID/SWAP separable guard): %s (id_tot=%u sw_tot=%u)\n", ok2 ? "PASS" : "FAIL",
           (unsigned)ra_id.rank_tot, (unsigned)ra_sw.rank_tot);
    if (ok2) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T3: Reference nonlinear entangled active-current Noether universe */
    const uint8_t p_ref[16] = {0U, 4U, 8U,  7U, 1U,  5U,  9U,  11U,
                               2U, 6U, 10U, 3U, 14U, 12U, 13U, 15U};
    phy_rule_analysis_t ra_ref;
    phy_analyze_rule(p_ref, 3ULL, &ra_ref);
    bool ok3 = (ra_ref.is_bijective == 1U && ra_ref.anf_deg == 3U && ra_ref.is_factorizable == 0U &&
                ra_ref.rank_tot == 5U && ra_ref.rank_ent == 1U && ra_ref.j_nonzero >= 200U &&
                ra_ref.continuity_fail == 0);
    if (ok3) {
        /* Verify PT conjugate symmetry preserves rank_ent & continuity */
        uint8_t p_pt[16];
        for (uint8_t u = 0U; u < 16U; ++u) {
            uint8_t su = p_sw[u];
            p_pt[u] = p_sw[p_ref[su]];
        }
        phy_rule_analysis_t ra_pt;
        phy_analyze_rule(p_pt, 4ULL, &ra_pt);
        if (ra_pt.rank_ent != ra_ref.rank_ent || ra_pt.continuity_fail != 0) {
            ok3 = false;
        }
        for (uint8_t a = 0U; a < 16U && ok3; ++a) {
            for (uint8_t b = 0U; b < 16U && ok3; ++b) {
                for (uint8_t c = 0U; c < 16U; ++c) {
                    int64_t lhs =
                        (int64_t)ra_ref.q[phy_step3(ra_ref.p, a, b, c)] - (int64_t)ra_ref.q[b];
                    int64_t rhs = (int64_t)ra_ref.j[b][c] - (int64_t)ra_ref.j[a][b];
                    if (lhs != rhs) {
                        ok3 = false;
                        break;
                    }
                }
            }
        }
    }
    printf("T3 (Ref nonlinear entangled Noether): %s (deg=%u tot=%u ent=%u j_nz=%u)\n",
           ok3 ? "PASS" : "FAIL", (unsigned)ra_ref.anf_deg, (unsigned)ra_ref.rank_tot,
           (unsigned)ra_ref.rank_ent, (unsigned)ra_ref.j_nonzero);
    if (ok3) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T4: 64-site periodic Floquet exact reversibility + zero charge drift */
    bool ok4 = true;
    uint64_t init_lo = UINT64_C(0x123456789ABCDEF0);
    uint64_t init_hi = UINT64_C(0x0FEDCBA987654321);
    int64_t q_init = phy_ring_charge_64(ra_ref.q, init_lo, init_hi);
    uint64_t clo = init_lo;
    uint64_t chi = init_hi;
    for (int t = 0; t < 17; ++t) {
        uint64_t nlo = UINT64_C(0);
        uint64_t nhi = UINT64_C(0);
        phy_step_64ring(ra_ref.p, clo, chi, &nlo, &nhi);
        clo = nlo;
        chi = nhi;
        if (phy_ring_charge_64(ra_ref.q, clo, chi) != q_init) {
            ok4 = false;
        }
    }
    for (int t = 0; t < 17; ++t) {
        uint64_t nlo = UINT64_C(0);
        uint64_t nhi = UINT64_C(0);
        (void)phy_rotl64_2bit_left(&clo, &chi);
        phy_step_64ring(ra_ref.pinv, clo, chi, &nlo, &nhi);
        phy_rotl64_2bit_right(&nlo, &nhi);
        clo = nlo;
        chi = nhi;
    }
    if (clo != init_lo || chi != init_hi) {
        ok4 = false;
    }
    printf("T4 (64-site Floquet reversibility & q-drift): %s\n", ok4 ? "PASS" : "FAIL");
    if (ok4) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T5: FLOW 2-Jet soliton damage tracker */
    phy_arena_init(&g_arena, 5ULL);
    phy_sim_trace_t tr5;
    phy_sim_soliton_64(&g_arena, &ra_ref, init_lo, init_hi, 32U, 32U, 3U, &tr5);
    bool ok5 = (tr5.max_charge_drift == 0 && tr5.frames_recorded == 32U && tr5.max_span > 0U &&
                tr5.vel_num_q16 > 0 && tr5.final_merkle != UINT64_C(0) &&
                tr5.final_merkle != ra_ref.merkle && g_arena.map.sim_frames[31].step_idx == 31U &&
                g_arena.map.sim_frames[31].total_charge == q_init);
    printf("T5 (FLOW 2-Jet soliton tracker): %s (span=%u vel_q16=%d drift=%lld)\n",
           ok5 ? "PASS" : "FAIL", (unsigned)tr5.max_span, tr5.vel_num_q16,
           (long long)tr5.max_charge_drift);
    if (ok5) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T6: PHY_HUNT_ENTANGLED discovery */
    phy_rule_analysis_t best_ent;
    uint32_t cnt_ent = phy_hunt(&g_arena, PHY_HUNT_ENTANGLED, 2026ULL, 2000U, &best_ent);
    bool ok6 =
        (cnt_ent >= 1U && best_ent.is_bijective == 1U && best_ent.rank_ent >= 1U &&
         best_ent.j_nonzero > 0U && best_ent.anf_deg >= 2U && best_ent.is_factorizable == 0U &&
         best_ent.continuity_fail == 0 && g_arena.map.img_hdr.magic == PHY_MAGIC);
    printf("T6 (Hunt ENTANGLED): %s (found=%u deg=%u ent=%u j_nz=%u merkle=0x%016llx)\n",
           ok6 ? "PASS" : "FAIL", (unsigned)cnt_ent, (unsigned)best_ent.anf_deg,
           (unsigned)best_ent.rank_ent, (unsigned)best_ent.j_nonzero,
           (unsigned long long)best_ent.merkle);
    if (ok6) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T7: PHY_HUNT_YANGBOX discovery */
    phy_rule_analysis_t best_yb;
    uint32_t cnt_yb = phy_hunt(&g_arena, PHY_HUNT_YANGBOX, 2026ULL, 6000U, &best_yb);
    bool ok7 = (cnt_yb >= 1U && best_yb.is_bijective == 1U && best_yb.yb_violations == 0U &&
                best_yb.anf_deg >= 2U && best_yb.is_factorizable == 0U && best_yb.is_swap == 0U &&
                best_yb.is_id == 0U);
    printf("T7 (Hunt YANGBOX): %s (found=%u deg=%u yb_viol=%u lcm=%u)\n", ok7 ? "PASS" : "FAIL",
           (unsigned)cnt_yb, (unsigned)best_yb.anf_deg, (unsigned)best_yb.yb_violations,
           (unsigned)best_yb.cycle_lcm);
    if (ok7) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T8: PHY_HUNT_JOINT discovery */
    phy_rule_analysis_t best_jnt;
    uint32_t cnt_jnt = phy_hunt(&g_arena, PHY_HUNT_JOINT, 2026ULL, 25000U, &best_jnt);
    bool ok8 = (cnt_jnt >= 1U && best_jnt.is_bijective == 1U && best_jnt.yb_violations == 0U &&
                best_jnt.rank_ent >= 1U && best_jnt.anf_deg >= 2U &&
                best_jnt.is_factorizable == 0U && best_jnt.continuity_fail == 0);
    printf("T8 (Hunt JOINT): %s (found=%u deg=%u yb_viol=%u tot=%u ent=%u j_nz=%u)\n",
           ok8 ? "PASS" : "FAIL", (unsigned)cnt_jnt, (unsigned)best_jnt.anf_deg,
           (unsigned)best_jnt.yb_violations, (unsigned)best_jnt.rank_tot,
           (unsigned)best_jnt.rank_ent, (unsigned)best_jnt.j_nonzero);
    if (ok8) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T9: searu_solve_mc_gf2 + FLOW_IR_MULT_LINEAR .fir synthesis & 16/16 eval parity */
    bool ok9 = ok6 && ok7 && ok8;
    uint16_t ent_and = 0U;
    uint16_t yb_and = 0U;
    uint16_t jnt_and = 0U;
    if (ok9) {
        if (!phy_synth_fir_rule(&best_ent, NULL, &g_fir_ws)) {
            ok9 = false;
        } else {
            ent_and = g_fir_ws.and_count;
        }
        if (ok9 && !phy_synth_fir_rule(&best_yb, NULL, &g_fir_ws)) {
            ok9 = false;
        } else if (ok9) {
            yb_and = g_fir_ws.and_count;
        }
        if (ok9 && !phy_synth_fir_rule(&best_jnt, NULL, &g_fir_ws)) {
            ok9 = false;
        } else if (ok9) {
            jnt_and = g_fir_ws.and_count;
        }
    }
    printf("T9 (FLOW_IR_MULT_LINEAR .fir synthesis): %s (ent_and=%u yb_and=%u jnt_and=%u "
           "unitary=true)\n",
           ok9 ? "PASS" : "FAIL", (unsigned)ent_and, (unsigned)yb_and, (unsigned)jnt_and);
    if (ok9) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T10: Arena ABI pack/unpack round-trip, Merkle determinism, step3 vs 64-ring consistency */
    bool ok10 = true;
    phy_rule_analysis_t ra_again;
    phy_analyze_rule(p_ref, 3ULL, &ra_again);
    if (ra_again.merkle != ra_ref.merkle || ra_again.lut_packed != ra_ref.lut_packed ||
        ra_ref.merkle == UINT64_C(0) || ra_ref.merkle == ra_id.merkle) {
        ok10 = false;
    }
    flow_cell packed;
    phy_rule_cell_t rc_in;
    memset(&rc_in, 0, sizeof(rc_in));
    rc_in.lut_packed = ra_ref.lut_packed;
    rc_in.pinv_packed = ra_ref.pinv_packed;
    rc_in.merkle = ra_ref.merkle;
    rc_in.seed = ra_ref.seed;
    rc_in.flags = ra_ref.flags;
    rc_in.anf_deg = ra_ref.anf_deg;
    rc_in.yb_violations = ra_ref.yb_violations;
    rc_in.rank_tot = ra_ref.rank_tot;
    rc_in.rank_ent = ra_ref.rank_ent;
    rc_in.j_nonzero = ra_ref.j_nonzero;
    rc_in.cycle_lcm = ra_ref.cycle_lcm;
    PHY_PhyRuleCell_pack(packed, &rc_in);
    phy_rule_cell_t rc_out;
    PHY_PhyRuleCell_unpack(packed, &rc_out);
    if (rc_out.lut_packed != rc_in.lut_packed || rc_out.merkle != rc_in.merkle ||
        rc_out.j_nonzero != rc_in.j_nonzero || PHY_SCHEMA_HASH == UINT64_C(0)) {
        ok10 = false;
    }
    uint64_t rlo = UINT64_C(0x713A94D2C5501B83);
    uint64_t rhi = UINT64_C(0x4E8129A0B7C6F012);
    uint64_t slo = UINT64_C(0);
    uint64_t shi = UINT64_C(0);
    phy_step_64ring(ra_ref.p, rlo, rhi, &slo, &shi);
    for (uint32_t m = 0U; m < 32U; ++m) {
        uint32_t mp = (m + 31U) & 31U;
        uint32_t mn = (m + 1U) & 31U;
        uint64_t w_prev = (mp < 16U) ? rlo : rhi;
        uint64_t w_cur = (m < 16U) ? rlo : rhi;
        uint64_t w_next = (mn < 16U) ? rlo : rhi;
        uint8_t a = (uint8_t)((w_prev >> ((mp & 15U) * 4U)) & UINT64_C(0xF));
        uint8_t b = (uint8_t)((w_cur >> ((m & 15U) * 4U)) & UINT64_C(0xF));
        uint8_t c = (uint8_t)((w_next >> ((mn & 15U) * 4U)) & UINT64_C(0xF));
        /* Pack order: site 2m is bits 0..1 of pair, site 2m+1 is bits 2..3 of pair in word,
           whereas step3 takes (hi=left=2m, lo=right=2m+1), so swap nibbles when reading/writing! */
        uint8_t a_sw = (uint8_t)(((a & 3U) << 2U) | ((a >> 2U) & 3U));
        uint8_t b_sw = (uint8_t)(((b & 3U) << 2U) | ((b >> 2U) & 3U));
        uint8_t c_sw = (uint8_t)(((c & 3U) << 2U) | ((c >> 2U) & 3U));
        uint8_t expect_sw = phy_step3(ra_ref.p, a_sw, b_sw, c_sw);
        uint64_t w_out = (m < 16U) ? slo : shi;
        uint8_t got = (uint8_t)((w_out >> ((m & 15U) * 4U)) & UINT64_C(0xF));
        uint8_t got_sw = (uint8_t)(((got & 3U) << 2U) | ((got >> 2U) & 3U));
        if (got_sw != expect_sw) {
            ok10 = false;
        }
    }
    printf("T10 (Arena ABI / step3-64ring parity): %s\n", ok10 ? "PASS" : "FAIL");
    if (ok10) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T11: Q4 exact 4-site hierarchy & 64-ring zero-drift */
    phy_q4_analysis_t q4_sw;
    phy_q4_analysis_t q4_ref;
    phy_q4_analysis_t q4_yb;
    phy_analyze_q4(p_sw, &g_q4_ws, &q4_sw);
    phy_analyze_q4(p_ref, &g_q4_ws, &q4_ref);
    const uint8_t p_yb_q4[16] = {0U, 4U, 8U,  14U, 1U, 5U, 9U,  13U,
                                 2U, 6U, 10U, 12U, 3U, 7U, 11U, 15U};
    phy_analyze_q4(p_yb_q4, &g_q4_ws, &q4_yb);
    bool ok11 = (q4_sw.genuine_q4 == 0U && q4_sw.tot2 == 6U && q4_ref.genuine_q4 >= 1U &&
                 q4_ref.continuity4_fail == 0 && q4_yb.tot2 == 5U && q4_yb.genuine_q4 == 1U &&
                 q4_yb.j4_nonzero > 0U && q4_yb.continuity4_fail == 0);
    uint64_t q4_lo = UINT64_C(0x123456789ABCDEF0);
    uint64_t q4_hi = UINT64_C(0x0FEDCBA987654321);
    int64_t c0_q4 = phy_ring_charge4_64(q4_yb.q4, q4_lo, q4_hi);
    for (uint32_t st = 0U; st < 32U && ok11; ++st) {
        uint64_t nlo = UINT64_C(0);
        uint64_t nhi = UINT64_C(0);
        phy_step_64ring(p_yb_q4, q4_lo, q4_hi, &nlo, &nhi);
        q4_lo = nlo;
        q4_hi = nhi;
        if (phy_ring_charge4_64(q4_yb.q4, q4_lo, q4_hi) != c0_q4) {
            ok11 = false;
        }
    }
    printf("T11 (Q4 hierarchy / 64-ring drift=0): %s (sw_gen=%u ref_gen=%u yb_gen=%u j4_nz=%u)\n",
           ok11 ? "PASS" : "FAIL", (unsigned)q4_sw.genuine_q4, (unsigned)q4_ref.genuine_q4,
           (unsigned)q4_yb.genuine_q4, (unsigned)q4_yb.j4_nonzero);
    if (ok11) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T12: Multi-point butterfly velocity v_B */
    phy_butterfly_t bf;
    phy_measure_butterfly_64(&ra_ref, UINT64_C(0x123456789ABCDEF0), UINT64_C(0x0FEDCBA987654321),
                             24U, &bf);
    bool ok12 = (bf.active_probes == 6U && bf.max_span == 5U && bf.mean_vel_q16 > 0 &&
                 bf.max_vel_q16 >= bf.mean_vel_q16 && bf.max_charge_drift == 0);
    printf("T12 (Multi-point butterfly v_B): %s (probes=%u span=%u vmax_q16=%d drift=%lld)\n",
           ok12 ? "PASS" : "FAIL", (unsigned)bf.active_probes, (unsigned)bf.max_span,
           bf.max_vel_q16, (long long)bf.max_charge_drift);
    if (ok12) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T13: |X|=8 non-involutive non-degenerate indecomposable YBE atlas */
    uint8_t r_sw3[64];
    uint8_t r_bad[64];
    for (uint8_t u = 0U; u < 64U; ++u) {
        r_sw3[u] = (uint8_t)(((u & 7U) << 3U) | ((u >> 3U) & 7U));
        r_bad[u] = r_sw3[u];
    }
    uint8_t r_tmp = r_bad[0];
    r_bad[0] = r_bad[1];
    r_bad[1] = r_tmp;
    phy_yb8_analysis_t yb8;
    uint32_t cnt_yb8 = phy_hunt_yb8(&g_arena, 2026ULL, 4000U, &yb8);
    bool ok13 = (phy_yb8_check(r_sw3) == 0U && phy_yb8_check(r_bad) > 0U && cnt_yb8 >= 1U &&
                 yb8.is_bijective == 1U && yb8.yb8_violations == 0U && yb8.is_nondeg == 1U &&
                 yb8.is_involutive == 0U && yb8.is_indecomposable == 1U &&
                 yb8.is_lyubashenko == 0U && yb8.anf_deg >= 2U && yb8.order >= 4U);
    printf("T13 (|X|=8 non-involutive YBE): %s (found=%u order=%u deg=%u nondeg=%u indecomp=%u)\n",
           ok13 ? "PASS" : "FAIL", (unsigned)cnt_yb8, (unsigned)yb8.order, (unsigned)yb8.anf_deg,
           (unsigned)yb8.is_nondeg, (unsigned)yb8.is_indecomposable);
    if (ok13) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T14: |X|=8 6-bit FLOW_IR_MULT_LINEAR .fir synthesis & 64/64 eval parity */
    bool ok14 = ok13;
    if (ok14) {
        bool synth_ok = phy_synth_fir_yb8(&yb8, NULL, &g_fir_ws);
        if (!synth_ok || !g_fir_ws.ir_linear_ok || !g_fir_ws.ir_unitary_ok ||
            !g_fir_ws.ir_eval_ok || g_fir_ws.wire_len == 0U) {
            ok14 = false;
        }
    }
    printf("T14 (|X|=8 6-bit FLOW_IR .fir parity): %s (and=%u depth=%u bytes=%zu eval=%s)\n",
           ok14 ? "PASS" : "FAIL", (unsigned)g_fir_ws.and_count, (unsigned)g_fir_ws.mult_depth,
           g_fir_ws.wire_len, g_fir_ws.ir_eval_ok ? "64/64" : "FAIL");
    if (ok14) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T15: Nonlinear Ergodic & Mixing Dual-Unitary + Universal Q2 No-Go */
    phy_dual4_analysis_t du4;
    phy_dual4_analysis_t du4_id;
    phy_dual4_analysis_t du4_sw;
    uint8_t p_id4[16];
    uint8_t p_sw4[16];
    for (uint8_t u = 0U; u < 16U; ++u) {
        p_id4[u] = u;
        p_sw4[u] = (uint8_t)(((u & 3U) << 2U) | ((u >> 2U) & 3U));
    }
    phy_analyze_dual4(p_id4, &du4_id);
    phy_analyze_dual4(p_sw4, &du4_sw);
    uint32_t cnt_du4 = phy_hunt_dual4(&g_arena, 2026ULL, 2500U, &du4);
    bool ok15 =
        (cnt_du4 >= 1U && du4.is_bijective == 1U && du4.is_nondeg == 1U &&
         du4.is_dual_unitary == 1U && du4.is_mixing == 1U && du4.is_non_bernoulli == 1U &&
         du4.anf_deg >= 3U && du4.rank_ent == 0U && du4_id.is_bijective == 1U &&
         du4_id.is_dual_unitary == 0U && du4_sw.is_dual_unitary == 1U && du4_sw.is_mixing == 0U);
    if (ok15) {
        if (!phy_synth_fir_4bit(du4.p, du4.merkle ^ UINT64_C(0x44553431), "phy_dual4_2x2", NULL,
                                &g_fir_ws)) {
            ok15 = false;
        }
    }
    printf("T15 (Dual-Unitary mixing / Q2 No-Go): %s (found=%u deg=%u ent=%u mix=%u and=%u)\n",
           ok15 ? "PASS" : "FAIL", (unsigned)cnt_du4, (unsigned)du4.anf_deg, (unsigned)du4.rank_ent,
           (unsigned)du4.is_mixing, (unsigned)g_fir_ws.and_count);
    if (ok15) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T16: Nonlinear 2-Unitary Perfect Tensor AME(4,8) on BitVec 6 */
    phy_ame8_analysis_t ame8;
    phy_ame8_analysis_t ame_neg;
    phy_analyze_ame8(r_sw3, &ame_neg);
    uint32_t cnt_ame8 = phy_hunt_ame8(&g_arena, 2026ULL, 400U, &ame8);
    bool ok16 = (cnt_ame8 >= 1U && ame8.is_bijective == 1U && ame8.is_nondeg == 1U &&
                 ame8.is_crossbij == 1U && ame8.is_2unitary == 1U && ame8.anf_deg >= 2U &&
                 ame_neg.is_crossbij == 0U && ame_neg.is_2unitary == 0U);
    if (ok16) {
        if (!phy_synth_fir_6bit(ame8.r, ame8.merkle ^ UINT64_C(0x414D4538), "phy_ame8_3x3", NULL,
                                &g_fir_ws)) {
            ok16 = false;
        }
    }
    printf("T16 (2-Unitary Perfect Tensor AME(4,8)): %s (found=%u deg=%u order=%u and=%u)\n",
           ok16 ? "PASS" : "FAIL", (unsigned)cnt_ame8, (unsigned)ame8.anf_deg, (unsigned)ame8.order,
           (unsigned)g_fir_ws.and_count);
    if (ok16) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T17: 3-Site Fracton Dipole-Laplacian + Strong HSF + Single-Charge Immobility */
    uint8_t t_osc[64];
    for (uint8_t u = 0U; u < 64U; ++u) {
        t_osc[u] = u;
    }
    t_osc[16] = 8U;
    t_osc[8] = 16U;
    phy_fracton_analysis_t frac;
    phy_fracton_analysis_t frac_neg;
    phy_analyze_fracton(t_osc, &g_frac_ws, &frac_neg);
    uint32_t cnt_frac = phy_hunt_fracton(&g_arena, &g_frac_ws, 2026ULL, 250U, &frac);
    bool ok17 =
        (cnt_frac >= 1U && frac.is_bijective == 1U && frac.laplacian_ok == 1U &&
         frac.cons_qd_4096_ok == 1U && frac.max_q_drift == 0 && frac.max_d_drift == 0 &&
         frac.hsf_num_comps >= 8U && (uint32_t)frac.hsf_largest_comp * 2U < frac.hsf_sector_size &&
         frac.iso_span == 0U && frac.dip_span > 0U && frac.anf_deg >= 3U &&
         frac_neg.iso_span == 1U && (frac_neg.flags & PHY_FRACTON_FLAG_ISO_FROZEN) == 0U);
    if (ok17) {
        if (!phy_synth_fir_6bit(frac.t, frac.merkle ^ UINT64_C(0x46524143), "phy_fracton_3s", NULL,
                                &g_fir_ws)) {
            ok17 = false;
        }
    }
    printf("T17 (Fracton Laplacian / Strong HSF): %s (found=%u comps=%u/%u largest=%u iso=%u "
           "dip=%u)\n",
           ok17 ? "PASS" : "FAIL", (unsigned)cnt_frac, (unsigned)frac.hsf_num_comps,
           (unsigned)frac.hsf_sector_size, (unsigned)frac.hsf_largest_comp, (unsigned)frac.iso_span,
           (unsigned)frac.dip_span);
    if (ok17) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T18: Irretractable |X|=8 R1+R2+R3 Biquandle & Knot Coloring Separation */
    uint8_t r_nonr1[64];
    for (uint8_t u = 0U; u < 64U; ++u) {
        r_nonr1[u] = r_sw3[u];
    }
    r_nonr1[0] = 9U;
    r_nonr1[9] = 0U;
    phy_biquandle8_analysis_t bq8;
    phy_biquandle8_analysis_t bq_neg_r1;
    phy_biquandle8_analysis_t bq_neg_ret;
    phy_analyze_biquandle8(r_nonr1, &bq_neg_r1);
    phy_analyze_biquandle8(r_sw3, &bq_neg_ret);
    uint32_t cnt_bq8 = phy_hunt_biquandle8(&g_arena, 2026ULL, 400U, &bq8);
    bool ok18 = (cnt_bq8 >= 1U && bq8.is_r1 == 1U && bq8.is_irretract == 1U &&
                 bq8.base.is_indecomposable == 1U && bq8.base.is_involutive == 0U &&
                 bq8.base.yb8_violations == 0U && bq8.col_unknot == 8U && bq8.col_torus7 == 64U &&
                 bq8.base.anf_deg >= 3U && bq_neg_r1.is_r1 == 0U &&
                 bq_neg_ret.retract_classes == 1U && bq_neg_ret.is_irretract == 0U);
    if (ok18) {
        if (!phy_synth_fir_6bit(bq8.base.r, bq8.base.merkle ^ UINT64_C(0x42513831), "phy_bq8_3x3",
                                NULL, &g_fir_ws)) {
            ok18 = false;
        }
    }
    printf("T18 (Irretractable |X|=8 Biquandle Knots): %s (ret=%u R1=%u T7=%u and=%u)\n",
           ok18 ? "PASS" : "FAIL", (unsigned)bq8.retract_classes, (unsigned)bq8.is_r1,
           (unsigned)bq8.col_torus7, (unsigned)g_fir_ws.and_count);
    if (ok18) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T19: Exhaustive S_8 |X|=2 Tetrahedron + Nonlinear Non-Involutive |X|=4 Tetrahedron */
    uint16_t t2_tot = 0U;
    uint16_t t2_nl = 0U;
    uint16_t t2_ni = 0U;
    phy_classify_tetra2(&t2_tot, &t2_nl, &t2_ni);
    phy_tetra4_analysis_t t4;
    uint32_t cnt_t4 = phy_hunt_tetra4(&g_arena, 2026ULL, 64U, &t4);
    bool ok19 = (t2_tot == 26U && t2_nl == 2U && t2_ni == 0U && cnt_t4 >= 1U &&
                 t4.is_bijective == 1U && t4.tetra4_viol == 0U && t4.is_involutive == 0U &&
                 t4.is_3coupled == 1U && t4.anf_deg >= 3U);
    if (ok19) {
        if (!phy_synth_fir_6bit(t4.t, t4.merkle ^ UINT64_C(0x54455434), "phy_tetra4_3c", NULL,
                                &g_fir_ws)) {
            ok19 = false;
        }
    }
    printf("T19 (2+1D Tetrahedron |X|=2 & |X|=4): %s (S8=%u/%u/%u t4_viol=%u deg=%u and=%u)\n",
           ok19 ? "PASS" : "FAIL", (unsigned)t2_tot, (unsigned)t2_nl, (unsigned)t2_ni,
           (unsigned)t4.tetra4_viol, (unsigned)t4.anf_deg, (unsigned)g_fir_ws.and_count);
    if (ok19) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T20: v3 Schema ABI Pack/Unpack Round-Trip + Quantum-Linear .fir Flag Verification */
    bool ok20 = ok15 && ok17 && ok19;
    if (ok20) {
        flow_cell c_d4;
        phy_dual4_meta_t d4_in;
        phy_dual4_meta_t d4_out;
        memset(&d4_in, 0, sizeof(d4_in));
        d4_in.lut_packed = du4.lut_packed;
        d4_in.merkle = du4.merkle;
        d4_in.mplus_packed = du4.mplus_packed;
        d4_in.mminus_packed = du4.mminus_packed;
        d4_in.flags = du4.flags;
        d4_in.anf_deg = du4.anf_deg;
        PHY_PhyDual4MetaCell_pack(c_d4, &d4_in);
        PHY_PhyDual4MetaCell_unpack(c_d4, &d4_out);

        flow_cell c_fr;
        phy_fracton_meta_t fr_in;
        phy_fracton_meta_t fr_out;
        memset(&fr_in, 0, sizeof(fr_in));
        fr_in.merkle = frac.merkle;
        fr_in.hsf_num_comps = frac.hsf_num_comps;
        fr_in.hsf_largest_comp = frac.hsf_largest_comp;
        PHY_PhyFractonMetaCell_pack(c_fr, &fr_in);
        PHY_PhyFractonMetaCell_unpack(c_fr, &fr_out);

        flow_cell c_t4;
        phy_tetra4_meta_t t4_in;
        phy_tetra4_meta_t t4_out;
        memset(&t4_in, 0, sizeof(t4_in));
        t4_in.merkle = t4.merkle;
        t4_in.tetra2_total = t2_tot;
        t4_in.is_3coupled = t4.is_3coupled;
        PHY_PhyTetra4MetaCell_pack(c_t4, &t4_in);
        PHY_PhyTetra4MetaCell_unpack(c_t4, &t4_out);

        if (d4_out.merkle != du4.merkle || fr_out.hsf_num_comps != frac.hsf_num_comps ||
            t4_out.tetra2_total != 26U) {
            ok20 = false;
        }
        if ((g_fir_ws.ir_reloaded.flags & FLOW_IR_FLAG_CLEAN) == 0U ||
            g_fir_ws.ir_reloaded.num_vals == 0U ||
            g_fir_ws.ir_reloaded.vals[0].is_classical != 0U) {
            ok20 = false;
        }
    }
    printf("T20 (v3 Arena ABI / Quantum-Linear .fir Gate): %s\n", ok20 ? "PASS" : "FAIL");
    if (ok20) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T21: Universal AME(4,4) Exhaustive Census & Affine Rigidity No-Go */
    phy_ame4_nogo_t ng;
    phy_classify_ame4(&g_arena, &ng);
    uint8_t p_bad2u[16];
    for (uint8_t u = 0U; u < 16U; ++u) {
        p_bad2u[u] = ng.sample_p[u];
    }
    uint8_t tmp2u = p_bad2u[0];
    p_bad2u[0] = p_bad2u[4];
    p_bad2u[4] = tmp2u;
    bool ok21 =
        (ng.latin4_count == 576U && ng.total_ols == 6912U && ng.affine_ols == 6912U &&
         ng.nonlin_ols == 0U && ng.bernoulli_ols == 6912U && phy_is_2unitary4(ng.sample_p) &&
         !phy_is_2unitary4(p_id4) && !phy_is_2unitary4(p_bad2u));
    printf("T21 (AME(4,4) Rigidity No-Go): %s (ls=%u ols=%u affine=%u nl=%u)\n",
           ok21 ? "PASS" : "FAIL", (unsigned)ng.latin4_count, (unsigned)ng.total_ols,
           (unsigned)ng.affine_ols, (unsigned)ng.nonlin_ols);
    if (ok21) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T22: 3+1D 4-Simplex Equation (|X|=2 & |X|=4) + Quantum-Linear .fir */
    phy_simplex4_2_t s2;
    phy_analyze_simplex4_2(&s2);
    phy_simplex4_4_t s4;
    uint32_t cnt_s4 = phy_hunt_simplex4(&g_arena, 2026ULL, 16U, &s4);
    uint8_t s_bad4[256];
    for (uint32_t u = 0U; u < 256U; ++u) {
        s_bad4[u] = s4.s[u];
    }
    uint8_t ts4 = s_bad4[0];
    s_bad4[0] = s_bad4[1];
    s_bad4[1] = ts4;
    bool ok22 =
        (s2.viol2_ccx02 == 0U && s2.viol2_ccx13 == 0U && s2.viol2_cccx1 > 0U &&
         s2.nl_count_b_ctrl == 2U && cnt_s4 >= 1U && s4.is_bijective == 1U && s4.viol4 == 0U &&
         s4.viol4_coupled == 0U && s4.is_involutive == 0U && s4.is_4coupled == 1U &&
         s4.order4 == 4U && s4.anf_deg >= 3U && phy_simplex4_check4(s_bad4) > 0U);
    if (ok22) {
        if (!phy_synth_fir_4bit(s4.ccx02, s4.merkle ^ UINT64_C(0x53343242), "phy_simplex4_ccx",
                                NULL, &g_fir_ws)) {
            ok22 = false;
        }
    }
    printf("T22 (3+1D 4-Simplex |X|=2 & |X|=4): %s (found=%u deg=%u ord=%u viol4=%u and=%u)\n",
           ok22 ? "PASS" : "FAIL", (unsigned)cnt_s4, (unsigned)s4.anf_deg, (unsigned)s4.order4,
           (unsigned)s4.viol4, (unsigned)g_fir_ws.and_count);
    if (ok22) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T23: Octahedral Q_6 Z_4 Knot Cohomology & Trefoil Chirality */
    phy_chiral6_analysis_t ch;
    uint32_t cnt_ch = phy_solve_chiral6(&g_arena, &ch);
    uint8_t phi_zero[36];
    uint8_t phi_bad[36];
    memset(phi_zero, 0, sizeof(phi_zero));
    memset(phi_bad, 0, sizeof(phi_bad));
    phi_bad[1] = 1U;
    phy_chiral6_analysis_t ch_z;
    phy_chiral6_analysis_t ch_b;
    phy_analyze_chiral6_custom(ch.act, phi_zero, &ch_z);
    phy_analyze_chiral6_custom(ch.act, phi_bad, &ch_b);
    bool ok23 = (cnt_ch == 1U && ch.cocycle_viol == 0U && ch.is_chiral == 1U &&
                 ch.is_non_coboundary == 1U && ch.unknot_w[0] == 6U && ch.hopf_w[0] == 6U &&
                 ch.hopf_w[2] == 6U && ch.tref_l_w[1] == 24U && ch.tref_r_w[1] == 0U &&
                 ch.tref_l_w[3] == 0U && ch.tref_r_w[3] == 24U && ch.fig8_w[0] == 6U &&
                 ch_z.cocycle_viol == 0U && ch_z.is_chiral == 0U && ch_z.is_non_coboundary == 0U &&
                 ch_b.cocycle_viol > 0U);
    printf("T23 (Q_6 Z_4 Knot Cohomology): %s (3_1L=[%u,%u,%u,%u] 3_1R=[%u,%u,%u,%u])\n",
           ok23 ? "PASS" : "FAIL", (unsigned)ch.tref_l_w[0], (unsigned)ch.tref_l_w[1],
           (unsigned)ch.tref_l_w[2], (unsigned)ch.tref_l_w[3], (unsigned)ch.tref_r_w[0],
           (unsigned)ch.tref_r_w[1], (unsigned)ch.tref_r_w[2], (unsigned)ch.tref_r_w[3]);
    if (ok23) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T24: Exact Z[i] Quantum Dual-Unitary + Scar / Period-2 DTC */
    phy_qscar_analysis_t qs;
    uint32_t ok_qs_flag = phy_analyze_qscar(&g_arena, &qs);
    static const int8_t w_id2_re[4][4] = {{2, 0, 0, 0}, {0, 2, 0, 0}, {0, 0, 2, 0}, {0, 0, 0, 2}};
    static const int8_t w_zero_im[4][4] = {{0}};
    uint8_t bad_u = 0U;
    uint8_t bad_d = 0U;
    uint8_t bad_f = 0U;
    phy_verify_zi_du4(w_id2_re, w_zero_im, &bad_u, &bad_d, &bad_f);
    phy_qscar_analysis_t qs_inv;
    (void)phy_analyze_qscar_custom(&g_arena, qs.w_re, qs.w_im, true, &qs_inv);
    bool ok24 =
        (ok_qs_flag == 1U && qs.du_unitary == 1U && qs.du_dual == 1U && qs.du_full_sup == 1U &&
         qs.du_purity_num[1] == UINT64_C(1) && qs.du_purity_den[1] == UINT64_C(8) &&
         qs.scar_period == 2U && qs.neel_sup[0] == 1U && qs.neel_sup[1] == 1U &&
         qs.neel_ret_num[0] == 0U && qs.neel_ret_num[1] == 1U && qs.therm_sup[3] >= 32U &&
         qs.therm_purity_num[3] > 0U &&
         qs.therm_purity_num[3] * UINT64_C(3) < qs.therm_purity_den[3] && bad_u == 1U &&
         bad_d == 0U && bad_f == 0U && qs_inv.neel_sup[0] > 1U);
    printf("T24 (Z[i] Quantum DU & Scar DTC): %s (du_t2=%llu/%llu dtc=%u therm_t4=%llu/%llu)\n",
           ok24 ? "PASS" : "FAIL", (unsigned long long)qs.du_purity_num[1],
           (unsigned long long)qs.du_purity_den[1], (unsigned)qs.scar_period,
           (unsigned long long)qs.therm_purity_num[3], (unsigned long long)qs.therm_purity_den[3]);
    if (ok24) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T25: Fracton Involution, Dipole-Toffoli, HSF L=6,7,8 Scaling & v4 ABI Round-Trip */
    uint8_t t_toff[64];
    phy_build_dipole_toffoli(t_toff);
    phy_fracton_comp_t fc;
    uint32_t ok_fc_flag = phy_analyze_fracton_comp(&g_arena, frac.t, t_toff, &g_frac_scale_ws, &fc);
    bool ok25 =
        (ok17 && ok21 && ok22 && ok23 && ok24 && ok_fc_flag == 1U && fc.is_involutive == 1U &&
         fc.ccx_ok == 1U && fc.laplacian_ok == 1U && fc.scaling_strict == 1U && fc.sec[0] == 43U &&
         fc.comps[0] == 34U && fc.largest[0] == 3U && fc.sec[1] == 130U && fc.comps[1] == 84U &&
         fc.largest[1] == 6U && fc.sec[2] == 410U && fc.comps[2] == 264U && fc.largest[2] == 6U);
    if (ok25) {
        if (!phy_synth_fir_6bit(fc.t_toff, fc.merkle ^ UINT64_C(0x4643544F), "phy_dipole_toff",
                                NULL, &g_fir_ws)) {
            ok25 = false;
        }
        flow_cell c_a4;
        phy_ame4_nogo_meta_t a4_in;
        phy_ame4_nogo_meta_t a4_out;
        memset(&a4_in, 0, sizeof(a4_in));
        a4_in.sample_lut = ng.sample_lut;
        a4_in.total_ols = ng.total_ols;
        PHY_PhyAme4NogoMetaCell_pack(c_a4, &a4_in);
        PHY_PhyAme4NogoMetaCell_unpack(c_a4, &a4_out);

        flow_cell c_s4;
        phy_simplex4_meta_t s4_in;
        phy_simplex4_meta_t s4_out;
        memset(&s4_in, 0, sizeof(s4_in));
        s4_in.merkle = s4.merkle;
        s4_in.order4 = s4.order4;
        PHY_PhySimplex4MetaCell_pack(c_s4, &s4_in);
        PHY_PhySimplex4MetaCell_unpack(c_s4, &s4_out);

        flow_cell c_ch;
        phy_chiral6_meta_t ch_in;
        phy_chiral6_meta_t ch_out;
        memset(&ch_in, 0, sizeof(ch_in));
        ch_in.merkle = ch.merkle;
        ch_in.tref_l_w1 = ch.tref_l_w[1];
        PHY_PhyChiral6MetaCell_pack(c_ch, &ch_in);
        PHY_PhyChiral6MetaCell_unpack(c_ch, &ch_out);

        flow_cell c_qs;
        phy_qscar_meta_t qs_in;
        phy_qscar_meta_t qs_out;
        memset(&qs_in, 0, sizeof(qs_in));
        qs_in.merkle = qs.merkle;
        qs_in.scar_period = qs.scar_period;
        PHY_PhyQscarMetaCell_pack(c_qs, &qs_in);
        PHY_PhyQscarMetaCell_unpack(c_qs, &qs_out);

        flow_cell c_fc;
        phy_frac_comp_meta_t fc_in;
        phy_frac_comp_meta_t fc_out;
        memset(&fc_in, 0, sizeof(fc_in));
        fc_in.merkle = fc.merkle;
        fc_in.sec8 = fc.sec[2];
        PHY_PhyFracCompMetaCell_pack(c_fc, &fc_in);
        PHY_PhyFracCompMetaCell_unpack(c_fc, &fc_out);

        if (a4_out.total_ols != 6912U || s4_out.order4 != 4U || ch_out.tref_l_w1 != 24U ||
            qs_out.scar_period != 2U || fc_out.sec8 != 410U) {
            ok25 = false;
        }
    }
    printf(
        "T25 (Fracton Comp / HSF L=6,7,8 / v4 ABI): %s (comps=%u,%u,%u largest=%u,%u,%u and=%u)\n",
        ok25 ? "PASS" : "FAIL", (unsigned)fc.comps[0], (unsigned)fc.comps[1], (unsigned)fc.comps[2],
        (unsigned)fc.largest[0], (unsigned)fc.largest[1], (unsigned)fc.largest[2],
        (unsigned)g_fir_ws.and_count);
    if (ok25) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    /* T26: Canonical "FLWA" Sparse Arena & "FIR2"/"FIR1" Wire Pack/Unpack Bit-Exact Roundtrip */
    size_t flwa_len = 0U;
    bool ok26 = ok25 &&
                phy_arena_pack_wire(&g_arena, g_arena_wire, sizeof(g_arena_wire), &flwa_len) &&
                phy_wire_is_flwa(g_arena_wire, flwa_len) && flwa_len < PHY_ARENA_BYTES;
    if (ok26) {
        bool was_packed = false;
        memset(&g_arena_rt, 0xA5, sizeof(g_arena_rt));
        if (!phy_arena_unpack_auto(g_arena_wire, flwa_len, &g_arena_rt, &was_packed) ||
            !was_packed || memcmp(&g_arena, &g_arena_rt, PHY_ARENA_BYTES) != 0) {
            ok26 = false;
        }
        /* Transparent legacy raw 32KB arena auto-decode */
        memset(&g_arena_rt, 0x5A, sizeof(g_arena_rt));
        if (!phy_arena_unpack_auto((const uint8_t *)&g_arena, PHY_ARENA_BYTES, &g_arena_rt,
                                   &was_packed) ||
            was_packed || memcmp(&g_arena, &g_arena_rt, PHY_ARENA_BYTES) != 0) {
            ok26 = false;
        }
        /* Reject tampered FLWA trailer/body */
        g_arena_wire[flwa_len - 1U] ^= UINT8_C(0x01);
        if (phy_arena_unpack_auto(g_arena_wire, flwa_len, &g_arena_rt, &was_packed)) {
            ok26 = false;
        }
        g_arena_wire[flwa_len - 1U] ^= UINT8_C(0x01);
    }
    if (ok26) {
        size_t rt_len = 0U;
        if (!phy_wire_is_fir2(g_fir_ws.wire_buf, g_fir_ws.wire_len) ||
            !phy_wire_is_fir1(g_fir_ws.legacy_wire_buf, g_fir_ws.legacy_wire_len) ||
            g_fir_ws.wire_len >= g_fir_ws.legacy_wire_len) {
            ok26 = false;
        }
        /* Legacy FIR1 -> auto-decode -> canonical FIR2 bit-exact match */
        if (ok26 &&
            (flow_ir_load_bytes_fiber(&g_eval_prog, g_fir_ws.legacy_wire_buf,
                                      g_fir_ws.legacy_wire_len) != FLOW_IR_FIBER_COMMITTED ||
             !flow_ir_encode_packed_bytes(&g_eval_prog, g_trans_out_buf, sizeof(g_trans_out_buf),
                                          &rt_len) ||
             rt_len != g_fir_ws.wire_len ||
             memcmp(g_trans_out_buf, g_fir_ws.wire_buf, rt_len) != 0)) {
            ok26 = false;
        }
        /* Canonical FIR2 -> auto-decode -> legacy FIR1 bit-exact match */
        if (ok26 && (flow_ir_load_bytes_fiber(&g_eval_prog, g_fir_ws.wire_buf, g_fir_ws.wire_len) !=
                         FLOW_IR_FIBER_COMMITTED ||
                     !flow_ir_encode_bytes(&g_eval_prog, g_trans_out_buf, sizeof(g_trans_out_buf),
                                           &rt_len) ||
                     rt_len != g_fir_ws.legacy_wire_len ||
                     memcmp(g_trans_out_buf, g_fir_ws.legacy_wire_buf, rt_len) != 0)) {
            ok26 = false;
        }
        /* Reject tampered FIR2 MAC */
        if (ok26) {
            memcpy(g_eval_buf, g_fir_ws.wire_buf, g_fir_ws.wire_len);
            g_eval_buf[g_fir_ws.wire_len - 1U] ^= UINT8_C(0x01);
            if (flow_ir_load_bytes_fiber(&g_eval_prog, g_eval_buf, g_fir_ws.wire_len) ==
                FLOW_IR_FIBER_COMMITTED) {
                ok26 = false;
            }
        }
    }
    printf("T26 (FLWA Arena & FIR2/FIR1 Wire Pack): %s (arena=%zu->%zuB fir=%zu->%zuB)\n",
           ok26 ? "PASS" : "FAIL", (size_t)PHY_ARENA_BYTES, flwa_len, g_fir_ws.legacy_wire_len,
           g_fir_ws.wire_len);
    if (ok26) {
        ++passed;
    } else if (fail_fast) {
        return 1;
    }

    printf("SELFTEST: %d/%d PASS\n", passed, total);
    return (passed == total) ? 0 : 1;
}

static int cmd_hunt_exec(uint32_t mode, const char *mode_name, uint64_t seed, uint64_t trials,
                         const char *out_lean, const char *out_fir, const char *out_json) {

    if (mode == PHY_HUNT_Q4_HIERARCHY) {
        phy_rule_analysis_t best;
        phy_q4_analysis_t best_q4;
        uint32_t found = phy_hunt_q4(&g_arena, &g_q4_ws, seed, trials, &best, &best_q4);
        if (found == 0U) {
            fprintf(stderr, "hunt mode=q4 found 0 matching universes\n");
            return 2;
        }
        phy_butterfly_t bf;
        phy_measure_butterfly_64(&best, UINT64_C(0x123456789ABCDEF0), UINT64_C(0x0FEDCBA987654321),
                                 24U, &bf);
        if (!phy_synth_fir_rule(&best, out_fir, &g_fir_ws)) {
            fprintf(stderr, "FIR synthesis failed for q4 rule\n");
            return 3;
        }
        if (out_lean != NULL && !phy_emit_lean_q4_cert(&best, &best_q4, out_lean)) {
            fprintf(stderr, "Lean cert emission failed for %s\n", out_lean);
            return 4;
        }
        if (out_json != NULL && !phy_emit_q4_json(&best_q4, &bf, &g_fir_ws, out_json)) {
            fprintf(stderr, "JSON emission failed for %s\n", out_json);
            return 5;
        }
        printf("HUNT q4: found=%u deg=%u yb_viol=%u tot2=%u tot4=%u genuine_q4=%u j4_nz=%u "
               "nnz=%u bf_span=%u bf_vmax_q16=%d and=%u merkle=0x%016llx\n",
               (unsigned)found, (unsigned)best.anf_deg, (unsigned)best.yb_violations,
               (unsigned)best_q4.tot2, (unsigned)best_q4.tot4, (unsigned)best_q4.genuine_q4,
               (unsigned)best_q4.j4_nonzero, (unsigned)best_q4.q4_nnz, (unsigned)bf.max_span,
               bf.max_vel_q16, (unsigned)g_fir_ws.and_count, (unsigned long long)best_q4.q4_merkle);
        return 0;
    }

    if (mode == PHY_HUNT_YB8) {
        phy_yb8_analysis_t best_yb8;
        uint32_t found = phy_hunt_yb8(&g_arena, seed, trials, &best_yb8);
        if (found == 0U) {
            fprintf(stderr, "hunt mode=yb8 found 0 matching universes\n");
            return 2;
        }
        if (!phy_synth_fir_yb8(&best_yb8, out_fir, &g_fir_ws)) {
            fprintf(stderr, "FIR synthesis failed for yb8\n");
            return 3;
        }
        if (out_lean != NULL && !phy_emit_lean_yb8_cert(&best_yb8, out_lean)) {
            fprintf(stderr, "Lean cert emission failed for %s\n", out_lean);
            return 4;
        }
        if (out_json != NULL && !phy_emit_yb8_json(&best_yb8, &g_fir_ws, out_json)) {
            fprintf(stderr, "JSON emission failed for %s\n", out_json);
            return 5;
        }
        printf("HUNT yb8: found=%u order=%u deg=%u yb8_viol=%u nondeg=%u noninvol=%u indecomp=%u "
               "nonlyub=%u and=%u merkle=0x%016llx\n",
               (unsigned)found, (unsigned)best_yb8.order, (unsigned)best_yb8.anf_deg,
               (unsigned)best_yb8.yb8_violations, (unsigned)best_yb8.is_nondeg,
               (unsigned)(best_yb8.is_involutive == 0U ? 1U : 0U),
               (unsigned)best_yb8.is_indecomposable,
               (unsigned)(best_yb8.is_lyubashenko == 0U ? 1U : 0U), (unsigned)g_fir_ws.and_count,
               (unsigned long long)best_yb8.merkle);
        return 0;
    }

    if (mode == PHY_HUNT_DUAL4) {
        phy_dual4_analysis_t best;
        uint32_t found = phy_hunt_dual4(&g_arena, seed, trials, &best);
        if (found == 0U) {
            fprintf(stderr, "hunt mode=dual4 found 0 matching universes\n");
            return 2;
        }
        if (!phy_synth_fir_4bit(best.p, best.merkle ^ UINT64_C(0x44553431), "phy_dual4_2x2",
                                out_fir, &g_fir_ws)) {
            fprintf(stderr, "FIR synthesis failed for dual4\n");
            return 3;
        }
        if (out_lean != NULL && !phy_emit_lean_dual4_cert(&best, out_lean)) {
            fprintf(stderr, "Lean cert emission failed for %s\n", out_lean);
            return 4;
        }
        if (out_json != NULL && !phy_emit_dual4_json(&best, &g_fir_ws, out_json)) {
            fprintf(stderr, "JSON emission failed for %s\n", out_json);
            return 5;
        }
        printf("HUNT dual4: found=%u deg=%u du=1 mix=1 ent=%u and=%u merkle=0x%016llx\n",
               (unsigned)found, (unsigned)best.anf_deg, (unsigned)best.rank_ent,
               (unsigned)g_fir_ws.and_count, (unsigned long long)best.merkle);
        return 0;
    }

    if (mode == PHY_HUNT_AME8) {
        phy_ame8_analysis_t best;
        uint32_t found = phy_hunt_ame8(&g_arena, seed, trials, &best);
        if (found == 0U) {
            fprintf(stderr, "hunt mode=ame8 found 0 matching universes\n");
            return 2;
        }
        if (!phy_synth_fir_6bit(best.r, best.merkle ^ UINT64_C(0x414D4538), "phy_ame8_3x3", out_fir,
                                &g_fir_ws)) {
            fprintf(stderr, "FIR synthesis failed for ame8\n");
            return 3;
        }
        if (out_lean != NULL && !phy_emit_lean_ame8_cert(&best, out_lean)) {
            fprintf(stderr, "Lean cert emission failed for %s\n", out_lean);
            return 4;
        }
        if (out_json != NULL && !phy_emit_ame8_json(&best, &g_fir_ws, out_json)) {
            fprintf(stderr, "JSON emission failed for %s\n", out_json);
            return 5;
        }
        printf("HUNT ame8: found=%u order=%u deg=%u 2unitary=1 and=%u merkle=0x%016llx\n",
               (unsigned)found, (unsigned)best.order, (unsigned)best.anf_deg,
               (unsigned)g_fir_ws.and_count, (unsigned long long)best.merkle);
        return 0;
    }

    if (mode == PHY_HUNT_FRACTON) {
        phy_fracton_analysis_t best;
        uint32_t found = phy_hunt_fracton(&g_arena, &g_frac_ws, seed, trials, &best);
        if (found == 0U) {
            fprintf(stderr, "hunt mode=fracton found 0 matching universes\n");
            return 2;
        }
        if (!phy_synth_fir_6bit(best.t, best.merkle ^ UINT64_C(0x46524143), "phy_fracton_3s",
                                out_fir, &g_fir_ws)) {
            fprintf(stderr, "FIR synthesis failed for fracton\n");
            return 3;
        }
        if (out_lean != NULL && !phy_emit_lean_fracton_cert(&best, out_lean)) {
            fprintf(stderr, "Lean cert emission failed for %s\n", out_lean);
            return 4;
        }
        if (out_json != NULL && !phy_emit_fracton_json(&best, &g_fir_ws, out_json)) {
            fprintf(stderr, "JSON emission failed for %s\n", out_json);
            return 5;
        }
        printf("HUNT fracton: found=%u deg=%u comps=%u/%u largest=%u frozen=%u iso=%u dip=%u "
               "and=%u merkle=0x%016llx\n",
               (unsigned)found, (unsigned)best.anf_deg, (unsigned)best.hsf_num_comps,
               (unsigned)best.hsf_sector_size, (unsigned)best.hsf_largest_comp,
               (unsigned)best.hsf_frozen, (unsigned)best.iso_span, (unsigned)best.dip_span,
               (unsigned)g_fir_ws.and_count, (unsigned long long)best.merkle);
        return 0;
    }

    if (mode == PHY_HUNT_BIQUANDLE8) {
        phy_biquandle8_analysis_t best;
        uint32_t found = phy_hunt_biquandle8(&g_arena, seed, trials, &best);
        if (found == 0U) {
            fprintf(stderr, "hunt mode=biquandle8 found 0 matching universes\n");
            return 2;
        }
        if (!phy_synth_fir_6bit(best.base.r, best.base.merkle ^ UINT64_C(0x42513831), "phy_bq8_3x3",
                                out_fir, &g_fir_ws)) {
            fprintf(stderr, "FIR synthesis failed for biquandle8\n");
            return 3;
        }
        if (out_lean != NULL && !phy_emit_lean_biquandle8_cert(&best, out_lean)) {
            fprintf(stderr, "Lean cert emission failed for %s\n", out_lean);
            return 4;
        }
        if (out_json != NULL && !phy_emit_biquandle8_json(&best, &g_fir_ws, out_json)) {
            fprintf(stderr, "JSON emission failed for %s\n", out_json);
            return 5;
        }
        printf("HUNT biquandle8: found=%u order=%u deg=%u ret=%u unknot=%u tref=%u fig8=%u T7=%u "
               "and=%u merkle=0x%016llx\n",
               (unsigned)found, (unsigned)best.base.order, (unsigned)best.base.anf_deg,
               (unsigned)best.retract_classes, (unsigned)best.col_unknot,
               (unsigned)best.col_trefoil, (unsigned)best.col_fig8, (unsigned)best.col_torus7,
               (unsigned)g_fir_ws.and_count, (unsigned long long)best.base.merkle);
        return 0;
    }

    if (mode == PHY_HUNT_TETRA4) {
        phy_tetra4_analysis_t best;
        uint32_t found = phy_hunt_tetra4(&g_arena, seed, trials, &best);
        if (found == 0U) {
            fprintf(stderr, "hunt mode=tetra4 found 0 matching universes\n");
            return 2;
        }
        if (!phy_synth_fir_6bit(best.t, best.merkle ^ UINT64_C(0x54455434), "phy_tetra4_3c",
                                out_fir, &g_fir_ws)) {
            fprintf(stderr, "FIR synthesis failed for tetra4\n");
            return 3;
        }
        if (out_lean != NULL && !phy_emit_lean_tetra4_cert(&best, out_lean)) {
            fprintf(stderr, "Lean cert emission failed for %s\n", out_lean);
            return 4;
        }
        if (out_json != NULL && !phy_emit_tetra4_json(&best, &g_fir_ws, out_json)) {
            fprintf(stderr, "JSON emission failed for %s\n", out_json);
            return 5;
        }
        printf("HUNT tetra4: found=%u order=%u deg=%u viol=%u tetra2=(%u,%u,%u) and=%u "
               "merkle=0x%016llx\n",
               (unsigned)found, (unsigned)best.order, (unsigned)best.anf_deg,
               (unsigned)best.tetra4_viol, (unsigned)best.tetra2_total, (unsigned)best.tetra2_nl,
               (unsigned)best.tetra2_noninvol_nl, (unsigned)g_fir_ws.and_count,
               (unsigned long long)best.merkle);
        return 0;
    }

    if (mode == PHY_HUNT_AME4_NOGO) {
        phy_ame4_nogo_t ng;
        phy_classify_ame4(&g_arena, &ng);
        if (out_lean != NULL && !phy_emit_lean_ame4_nogo_cert(&ng, out_lean)) {
            fprintf(stderr, "Lean cert emission failed for %s\n", out_lean);
            return 4;
        }
        if (out_json != NULL && !phy_emit_ame4_nogo_json(&ng, out_json)) {
            fprintf(stderr, "JSON emission failed for %s\n", out_json);
            return 5;
        }
        printf("HUNT ame4-nogo: ls=%u ols=%u affine=%u nl=%u bern=%u merkle=0x%016llx\n",
               (unsigned)ng.latin4_count, (unsigned)ng.total_ols, (unsigned)ng.affine_ols,
               (unsigned)ng.nonlin_ols, (unsigned)ng.bernoulli_ols, (unsigned long long)ng.merkle);
        return 0;
    }

    if (mode == PHY_HUNT_SIMPLEX4) {
        phy_simplex4_4_t best;
        uint32_t found = phy_hunt_simplex4(&g_arena, seed, trials, &best);
        if (found == 0U) {
            fprintf(stderr, "hunt mode=simplex4 found 0 matching universes\n");
            return 2;
        }
        if (!phy_synth_fir_4bit(best.ccx02, best.merkle ^ UINT64_C(0x53343242), "phy_simplex4_ccx",
                                out_fir, &g_fir_ws)) {
            fprintf(stderr, "FIR synthesis failed for simplex4\n");
            return 3;
        }
        if (out_lean != NULL && !phy_emit_lean_simplex4_cert(&best, out_lean)) {
            fprintf(stderr, "Lean cert emission failed for %s\n", out_lean);
            return 4;
        }
        if (out_json != NULL && !phy_emit_simplex4_json(&best, &g_fir_ws, out_json)) {
            fprintf(stderr, "JSON emission failed for %s\n", out_json);
            return 5;
        }
        printf("HUNT simplex4: found=%u order=%u deg=%u viol4=%u coupled=%u viol2=(%u,%u) and=%u "
               "merkle=0x%016llx\n",
               (unsigned)found, (unsigned)best.order4, (unsigned)best.anf_deg, (unsigned)best.viol4,
               (unsigned)best.is_4coupled, (unsigned)best.viol2_ccx, (unsigned)best.viol2_cccx,
               (unsigned)g_fir_ws.and_count, (unsigned long long)best.merkle);
        return 0;
    }

    if (mode == PHY_HUNT_CHIRAL6) {
        phy_chiral6_analysis_t ch;
        uint32_t found = phy_solve_chiral6(&g_arena, &ch);
        if (found == 0U) {
            fprintf(stderr, "hunt mode=chiral6 failed\n");
            return 2;
        }
        if (out_lean != NULL && !phy_emit_lean_chiral6_cert(&ch, out_lean)) {
            fprintf(stderr, "Lean cert emission failed for %s\n", out_lean);
            return 4;
        }
        if (out_json != NULL && !phy_emit_chiral6_json(&ch, out_json)) {
            fprintf(stderr, "JSON emission failed for %s\n", out_json);
            return 5;
        }
        printf("HUNT chiral6: found=%u rank=%u chiral=%u noncobound=%u L3_1=%u R3_1=%u "
               "merkle=0x%016llx\n",
               (unsigned)found, (unsigned)ch.cocycle_rank, (unsigned)ch.is_chiral,
               (unsigned)ch.is_non_coboundary, (unsigned)ch.tref_l_w[1], (unsigned)ch.tref_r_w[1],
               (unsigned long long)ch.merkle);
        return 0;
    }

    if (mode == PHY_HUNT_QSCAR) {
        phy_qscar_analysis_t qs;
        uint32_t found = phy_analyze_qscar(&g_arena, &qs);
        if (found == 0U) {
            fprintf(stderr, "hunt mode=qscar failed\n");
            return 2;
        }
        if (out_lean != NULL && !phy_emit_lean_qscar_cert(&qs, out_lean)) {
            fprintf(stderr, "Lean cert emission failed for %s\n", out_lean);
            return 4;
        }
        if (out_json != NULL && !phy_emit_qscar_json(&qs, out_json)) {
            fprintf(stderr, "JSON emission failed for %s\n", out_json);
            return 5;
        }
        printf("HUNT qscar: found=%u du=(%u,%u,%u) t2=%llu/%llu dtc=%u therm=%llu/%llu "
               "merkle=0x%016llx\n",
               (unsigned)found, (unsigned)qs.du_unitary, (unsigned)qs.du_dual,
               (unsigned)qs.du_full_sup, (unsigned long long)qs.du_purity_num[1],
               (unsigned long long)qs.du_purity_den[1], (unsigned)qs.scar_period,
               (unsigned long long)qs.therm_purity_num[3],
               (unsigned long long)qs.therm_purity_den[3], (unsigned long long)qs.merkle);
        return 0;
    }

    if (mode == PHY_HUNT_FRACTON_COMP) {
        phy_fracton_analysis_t fr_best;
        uint32_t cnt_f = phy_hunt_fracton(&g_arena, &g_frac_ws, seed, trials, &fr_best);
        if (cnt_f == 0U) {
            fprintf(stderr, "hunt mode=fracton-comp found 0 fracton rules\n");
            return 2;
        }
        uint8_t t_toff[64];
        phy_build_dipole_toffoli(t_toff);
        phy_fracton_comp_t fc;
        uint32_t found =
            phy_analyze_fracton_comp(&g_arena, fr_best.t, t_toff, &g_frac_scale_ws, &fc);
        if (found == 0U) {
            fprintf(stderr, "hunt mode=fracton-comp failed\n");
            return 2;
        }
        if (!phy_synth_fir_6bit(fc.t_toff, fc.merkle ^ UINT64_C(0x4643544F), "phy_dipole_toff",
                                out_fir, &g_fir_ws)) {
            fprintf(stderr, "FIR synthesis failed for fracton-comp\n");
            return 3;
        }
        if (out_lean != NULL && !phy_emit_lean_fracton_comp_cert(&fc, out_lean)) {
            fprintf(stderr, "Lean cert emission failed for %s\n", out_lean);
            return 4;
        }
        if (out_json != NULL && !phy_emit_fracton_comp_json(&fc, &g_fir_ws, out_json)) {
            fprintf(stderr, "JSON emission failed for %s\n", out_json);
            return 5;
        }
        printf("HUNT fracton-comp: found=%u ccx=%u invol=%u strict=%u comps=(%u/%u,%u/%u,%u/%u) "
               "and=%u merkle=0x%016llx\n",
               (unsigned)found, (unsigned)fc.ccx_ok, (unsigned)fc.is_involutive,
               (unsigned)fc.scaling_strict, (unsigned)fc.comps[0], (unsigned)fc.sec[0],
               (unsigned)fc.comps[1], (unsigned)fc.sec[1], (unsigned)fc.comps[2],
               (unsigned)fc.sec[2], (unsigned)g_fir_ws.and_count, (unsigned long long)fc.merkle);
        return 0;
    }

    phy_rule_analysis_t best;
    uint32_t found = phy_hunt(&g_arena, mode, seed, trials, &best);
    if (found == 0U) {
        fprintf(stderr, "hunt mode=%s found 0 matching universes\n", mode_name);
        return 2;
    }

    phy_sim_trace_t sim;
    phy_sim_soliton_64(&g_arena, &best, UINT64_C(0x123456789ABCDEF0), UINT64_C(0x0FEDCBA987654321),
                       32U, 32U, 3U, &sim);

    if (!phy_synth_fir_rule(&best, out_fir, &g_fir_ws)) {
        fprintf(stderr, "FIR synthesis failed (and=%u depth=%u lin=%d unit=%d eval=%d len=%zu)\n",
                (unsigned)g_fir_ws.and_count, (unsigned)g_fir_ws.mult_depth,
                g_fir_ws.ir_linear_ok ? 1 : 0, g_fir_ws.ir_unitary_ok ? 1 : 0,
                g_fir_ws.ir_eval_ok ? 1 : 0, g_fir_ws.wire_len);
        return 3;
    }

    if (out_lean != NULL) {
        bool cert_ok = (mode == PHY_HUNT_YANGBOX) ? phy_emit_lean_yb_cert(&best, out_lean)
                                                  : phy_emit_lean_noether_cert(&best, out_lean);
        if (!cert_ok) {
            fprintf(stderr, "Lean cert emission failed for %s\n", out_lean);
            return 4;
        }
    }

    if (out_json != NULL) {
        if (!phy_emit_rule_json(&best, &sim, &g_fir_ws, out_json)) {
            fprintf(stderr, "JSON emission failed for %s\n", out_json);
            return 5;
        }
    }

    printf("HUNT %s: found=%u deg=%u yb_viol=%u rank_tot=%u rank_ent=%u j_nz=%u and=%u "
           "span=%u vel_q16=%d merkle=0x%016llx\n",
           mode_name, (unsigned)found, (unsigned)best.anf_deg, (unsigned)best.yb_violations,
           (unsigned)best.rank_tot, (unsigned)best.rank_ent, (unsigned)best.j_nonzero,
           (unsigned)g_fir_ws.and_count, (unsigned)sim.max_span, sim.vel_num_q16,
           (unsigned long long)best.merkle);
    return 0;
}

static int cmd_hunt(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: phy hunt "
                        "<entangled|yangbox|joint|q4|yb8|dual4|ame8|fracton|biquandle8|tetra4|"
                        "ame4-nogo|simplex4|chiral6|qscar|fracton-comp> "
                        "[--seed S] [--trials N] [--out-lean P] [--out-fir P] [--out-json P] "
                        "[--out-arena P]\n");
        return 1;
    }
    uint32_t mode = PHY_HUNT_ENTANGLED;
    uint64_t trials = UINT64_C(4000);
    if (strcmp(argv[2], "entangled") == 0) {
        mode = PHY_HUNT_ENTANGLED;
        trials = UINT64_C(4000);
    } else if (strcmp(argv[2], "yangbox") == 0) {
        mode = PHY_HUNT_YANGBOX;
        trials = UINT64_C(8000);
    } else if (strcmp(argv[2], "joint") == 0) {
        mode = PHY_HUNT_JOINT;
        trials = UINT64_C(8000);
    } else if (strcmp(argv[2], "q4") == 0) {
        mode = PHY_HUNT_Q4_HIERARCHY;
        trials = UINT64_C(12000);
    } else if (strcmp(argv[2], "yb8") == 0) {
        mode = PHY_HUNT_YB8;
        trials = UINT64_C(4000);
    } else if (strcmp(argv[2], "dual4") == 0) {
        mode = PHY_HUNT_DUAL4;
        trials = UINT64_C(2500);
    } else if (strcmp(argv[2], "ame8") == 0) {
        mode = PHY_HUNT_AME8;
        trials = UINT64_C(400);
    } else if (strcmp(argv[2], "fracton") == 0) {
        mode = PHY_HUNT_FRACTON;
        trials = UINT64_C(250);
    } else if (strcmp(argv[2], "biquandle8") == 0) {
        mode = PHY_HUNT_BIQUANDLE8;
        trials = UINT64_C(400);
    } else if (strcmp(argv[2], "tetra4") == 0) {
        mode = PHY_HUNT_TETRA4;
        trials = UINT64_C(64);
    } else if (strcmp(argv[2], "ame4-nogo") == 0) {
        mode = PHY_HUNT_AME4_NOGO;
        trials = UINT64_C(1);
    } else if (strcmp(argv[2], "simplex4") == 0) {
        mode = PHY_HUNT_SIMPLEX4;
        trials = UINT64_C(16);
    } else if (strcmp(argv[2], "chiral6") == 0) {
        mode = PHY_HUNT_CHIRAL6;
        trials = UINT64_C(1);
    } else if (strcmp(argv[2], "qscar") == 0) {
        mode = PHY_HUNT_QSCAR;
        trials = UINT64_C(1);
    } else if (strcmp(argv[2], "fracton-comp") == 0) {
        mode = PHY_HUNT_FRACTON_COMP;
        trials = UINT64_C(250);
    } else {
        fprintf(stderr, "Unknown hunt mode: %s\n", argv[2]);
        return 1;
    }

    uint64_t seed = UINT64_C(2026);
    const char *out_lean = NULL;
    const char *out_fir = NULL;
    const char *out_json = NULL;
    const char *out_arena = NULL;

    for (int i = 3; i < argc; ++i) {
        if (strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
            seed = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--trials") == 0 && i + 1 < argc) {
            trials = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--out-lean") == 0 && i + 1 < argc) {
            out_lean = argv[++i];
        } else if (strcmp(argv[i], "--out-fir") == 0 && i + 1 < argc) {
            out_fir = argv[++i];
        } else if (strcmp(argv[i], "--out-json") == 0 && i + 1 < argc) {
            out_json = argv[++i];
        } else if (strcmp(argv[i], "--out-arena") == 0 && i + 1 < argc) {
            out_arena = argv[++i];
        }
    }

    int rc = cmd_hunt_exec(mode, argv[2], seed, trials, out_lean, out_fir, out_json);
    if (rc == 0 && out_arena != NULL) {
        size_t flwa_len = 0U;
        if (!phy_arena_save_flwa(&g_arena, out_arena, g_arena_wire, sizeof(g_arena_wire),
                                 &flwa_len)) {
            fprintf(stderr, "FLWA arena emission failed for %s\n", out_arena);
            return 6;
        }
    }
    return rc;
}

static int cmd_inspect(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: phy inspect <p0,p1,...,p15> [--out-lean P] [--out-fir P]\n");
        return 1;
    }
    uint8_t p[16];
    const char *s = argv[2];
    for (int i = 0; i < 16; ++i) {
        char *endp = NULL;
        long v = strtol(s, &endp, 10);
        if (endp == s || v < 0 || v > 15) {
            fprintf(stderr, "Invalid permutation element at %d\n", i);
            return 1;
        }
        p[i] = (uint8_t)v;
        if (i < 15) {
            if (*endp != ',') {
                fprintf(stderr, "Expected comma at %d\n", i);
                return 1;
            }
            s = endp + 1;
        }
    }
    const char *out_lean = NULL;
    const char *out_fir = NULL;
    for (int i = 3; i < argc; ++i) {
        if (strcmp(argv[i], "--out-lean") == 0 && i + 1 < argc) {
            out_lean = argv[++i];
        } else if (strcmp(argv[i], "--out-fir") == 0 && i + 1 < argc) {
            out_fir = argv[++i];
        }
    }
    phy_rule_analysis_t ra;
    phy_analyze_rule(p, UINT64_C(0), &ra);
    if (!phy_synth_fir_rule(&ra, out_fir, &g_fir_ws)) {
        fprintf(stderr, "FIR synthesis failed\n");
        return 2;
    }
    if (out_lean != NULL) {
        bool ok = (ra.rank_ent > 0U && ra.j_nonzero > 0U)
                      ? phy_emit_lean_noether_cert(&ra, out_lean)
                      : phy_emit_lean_yb_cert(&ra, out_lean);
        if (!ok) {
            fprintf(stderr, "Lean cert emission failed\n");
            return 3;
        }
    }
    printf("INSPECT: bij=%u deg=%u fact=%u yb_viol=%u tot=%u ent=%u j_nz=%u and=%u "
           "merkle=0x%016llx\n",
           (unsigned)ra.is_bijective, (unsigned)ra.anf_deg, (unsigned)ra.is_factorizable,
           (unsigned)ra.yb_violations, (unsigned)ra.rank_tot, (unsigned)ra.rank_ent,
           (unsigned)ra.j_nonzero, (unsigned)g_fir_ws.and_count, (unsigned long long)ra.merkle);
    return 0;
}

static int cmd_eval_fir_bits(int argc, char **argv, uint32_t bits) {
    if (argc < 3) {
        fprintf(stderr, "Usage: phy %s <path.fir>\n", (bits == 6U) ? "eval-fir6" : "eval-fir");
        return 1;
    }
    size_t len = 0U;
    bool was_packed = false;
    if (!phy_fir_load_auto_file(argv[2], &g_eval_prog, g_eval_buf, sizeof(g_eval_buf), &len,
                                &was_packed)) {
        fprintf(stderr, "Invalid FIR file: %s\n", argv[2]);
        return 2;
    }
    if (g_eval_prog.mult != FLOW_IR_MULT_LINEAR || !flow_ir_admits_unitary(&g_eval_prog)) {
        fprintf(stderr, "FIR failed linear/unitary admission\n");
        return 3;
    }
    uint32_t count = UINT32_C(1) << bits;
    uint64_t mask = (UINT64_C(1) << bits) - UINT64_C(1);
    for (uint32_t u = 0U; u < count; ++u) {
        uint64_t in_a[1] = {(uint64_t)u};
        uint64_t out_a[1] = {UINT64_C(0)};
        if (!flow_ir_eval_u64(&g_eval_prog, in_a, out_a)) {
            fprintf(stderr, "FIR eval error at u=%u\n", u);
            return 4;
        }
        printf("%s%llu", (u > 0U) ? " " : "", (unsigned long long)(out_a[0] & mask));
    }
    printf("\n");
    return 0;
}

static int cmd_arena_pack(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: phy arena-pack <in.bin|in.flwa> <out.flwa>\n");
        return 1;
    }
    size_t in_len = 0U;
    bool was_packed = false;
    if (!phy_arena_load_auto_file(argv[2], &g_arena_rt, g_arena_wire, sizeof(g_arena_wire), &in_len,
                                  &was_packed)) {
        fprintf(stderr, "Invalid arena file: %s\n", argv[2]);
        return 2;
    }
    size_t out_len = 0U;
    if (!phy_arena_save_flwa(&g_arena_rt, argv[3], g_arena_wire, sizeof(g_arena_wire), &out_len)) {
        fprintf(stderr, "Failed to write FLWA arena: %s\n", argv[3]);
        return 3;
    }
    printf("ARENA_PACK: %s (%zu B, %s) -> %s (%zu B, FLWA)\n", argv[2], in_len,
           was_packed ? "FLWA" : "raw", argv[3], out_len);
    return 0;
}

static int cmd_arena_unpack(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: phy arena-unpack <in.flwa|in.bin> <out.bin>\n");
        return 1;
    }
    size_t in_len = 0U;
    bool was_packed = false;
    if (!phy_arena_load_auto_file(argv[2], &g_arena_rt, g_arena_wire, sizeof(g_arena_wire), &in_len,
                                  &was_packed)) {
        fprintf(stderr, "Invalid arena file: %s\n", argv[2]);
        return 2;
    }
    if (!flow_ir_write_bytes_file(argv[3], (const uint8_t *)&g_arena_rt, PHY_ARENA_BYTES)) {
        fprintf(stderr, "Failed to write raw arena: %s\n", argv[3]);
        return 3;
    }
    printf("ARENA_UNPACK: %s (%zu B, %s) -> %s (%zu B, raw)\n", argv[2], in_len,
           was_packed ? "FLWA" : "raw", argv[3], (size_t)PHY_ARENA_BYTES);
    return 0;
}

static int cmd_fir_transcode(int argc, char **argv, bool pack) {
    if (argc < 4) {
        fprintf(stderr, "Usage: phy %s <in.fir> <out.fir>\n", pack ? "fir-pack" : "fir-unpack");
        return 1;
    }
    size_t in_len = 0U;
    size_t out_len = 0U;
    bool ok =
        pack ? phy_fir_pack_file(argv[2], argv[3], &g_eval_prog, g_eval_buf, sizeof(g_eval_buf),
                                 g_trans_out_buf, sizeof(g_trans_out_buf), &in_len, &out_len)
             : phy_fir_unpack_file(argv[2], argv[3], &g_eval_prog, g_eval_buf, sizeof(g_eval_buf),
                                   g_trans_out_buf, sizeof(g_trans_out_buf), &in_len, &out_len);
    if (!ok) {
        fprintf(stderr, "FIR %s failed: %s -> %s\n", pack ? "pack" : "unpack", argv[2], argv[3]);
        return 2;
    }
    printf("%s: %s (%zu B) -> %s (%zu B, %s)\n", pack ? "FIR_PACK" : "FIR_UNPACK", argv[2], in_len,
           argv[3], out_len, pack ? "FIR2" : "FIR1");
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "PHY v4.0.0 — Exact Noether, Dual-Unitary, 4-Simplex, Knot Cohomology, "
                        "Z[i] Quantum Scar & Fracton Engine\n");
        fprintf(stderr, "Usage: phy <selftest|hunt|inspect|eval-fir|eval-fir6|fir-pack|"
                        "fir-unpack|arena-pack|arena-unpack|classify-tetra2|classify-ame4> ...\n");
        return 1;
    }
    if (strcmp(argv[1], "selftest") == 0) {
        bool ff = (argc >= 3 && strcmp(argv[2], "--fail-fast") == 0);
        return phy_selftest(ff);
    }
    if (strcmp(argv[1], "hunt") == 0) {
        return cmd_hunt(argc, argv);
    }
    if (strcmp(argv[1], "inspect") == 0) {
        return cmd_inspect(argc, argv);
    }
    if (strcmp(argv[1], "eval-fir") == 0) {
        return cmd_eval_fir_bits(argc, argv, 4U);
    }
    if (strcmp(argv[1], "eval-fir6") == 0) {
        return cmd_eval_fir_bits(argc, argv, 6U);
    }
    if (strcmp(argv[1], "fir-pack") == 0) {
        return cmd_fir_transcode(argc, argv, true);
    }
    if (strcmp(argv[1], "fir-unpack") == 0) {
        return cmd_fir_transcode(argc, argv, false);
    }
    if (strcmp(argv[1], "arena-pack") == 0) {
        return cmd_arena_pack(argc, argv);
    }
    if (strcmp(argv[1], "arena-unpack") == 0) {
        return cmd_arena_unpack(argc, argv);
    }
    if (strcmp(argv[1], "classify-tetra2") == 0) {
        uint16_t tot = 0U;
        uint16_t nl = 0U;
        uint16_t ni = 0U;
        phy_classify_tetra2(&tot, &nl, &ni);
        printf("TETRA2 S_8: total=%u nonlinear=%u noninvol_nl=%u\n", (unsigned)tot, (unsigned)nl,
               (unsigned)ni);
        return 0;
    }
    if (strcmp(argv[1], "classify-ame4") == 0) {
        phy_ame4_nogo_t ng;
        phy_classify_ame4(&g_arena, &ng);
        printf("AME4 OLS: latin4=%u total=%u affine=%u nonlinear=%u bernoulli=%u\n",
               (unsigned)ng.latin4_count, (unsigned)ng.total_ols, (unsigned)ng.affine_ols,
               (unsigned)ng.nonlin_ols, (unsigned)ng.bernoulli_ols);
        return 0;
    }
    fprintf(stderr, "Unknown command: %s\n", argv[1]);
    return 1;
}
