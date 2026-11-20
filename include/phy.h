#ifndef PHY_H
#define PHY_H

#include "flow.h"
#include "phy_schema_abi.h"
#include "searu.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define PHY_MAGIC                   UINT32_C(0x50485931) /* "PHY1" */
#define PHY_VERSION                 UINT32_C(0x00030000) /* v3.0.0 */
#define PHY_M61                     UINT64_C(0x1FFFFFFFFFFFFFFF)

#define PHY_RULE_FLAG_VALID         (UINT32_C(1) << 0U)
#define PHY_RULE_FLAG_BIJECTIVE     (UINT32_C(1) << 1U)
#define PHY_RULE_FLAG_ENTANGLED     (UINT32_C(1) << 2U)
#define PHY_RULE_FLAG_ACTIVE_J      (UINT32_C(1) << 3U)
#define PHY_RULE_FLAG_NONLINEAR     (UINT32_C(1) << 4U)
#define PHY_RULE_FLAG_NONFACT       (UINT32_C(1) << 5U)
#define PHY_RULE_FLAG_YANGBOX       (UINT32_C(1) << 6U)
#define PHY_RULE_FLAG_JOINT         (UINT32_C(1) << 7U)

#define PHY_YB8_FLAG_BIJECTIVE      (UINT32_C(1) << 0U)
#define PHY_YB8_FLAG_YBE            (UINT32_C(1) << 1U)
#define PHY_YB8_FLAG_NONDEG         (UINT32_C(1) << 2U)
#define PHY_YB8_FLAG_NONINVOL       (UINT32_C(1) << 3U)
#define PHY_YB8_FLAG_INDECOMP       (UINT32_C(1) << 4U)
#define PHY_YB8_FLAG_NONLYUB        (UINT32_C(1) << 5U)
#define PHY_YB8_FLAG_NONLINEAR      (UINT32_C(1) << 6U)

#define PHY_BQ8_FLAG_R1             (UINT32_C(1) << 7U)
#define PHY_BQ8_FLAG_IRRETRACT      (UINT32_C(1) << 8U)
#define PHY_BQ8_FLAG_KNOT_SEP       (UINT32_C(1) << 9U)

#define PHY_DUAL4_FLAG_BIJECTIVE    (UINT32_C(1) << 0U)
#define PHY_DUAL4_FLAG_NONDEG       (UINT32_C(1) << 1U)
#define PHY_DUAL4_FLAG_DUAL         (UINT32_C(1) << 2U)
#define PHY_DUAL4_FLAG_MIXING       (UINT32_C(1) << 3U)
#define PHY_DUAL4_FLAG_NONBERN      (UINT32_C(1) << 4U)
#define PHY_DUAL4_FLAG_NONLINEAR    (UINT32_C(1) << 5U)

#define PHY_AME8_FLAG_BIJECTIVE     (UINT32_C(1) << 0U)
#define PHY_AME8_FLAG_NONDEG        (UINT32_C(1) << 1U)
#define PHY_AME8_FLAG_CROSSBIJ      (UINT32_C(1) << 2U)
#define PHY_AME8_FLAG_2UNITARY      (UINT32_C(1) << 3U)
#define PHY_AME8_FLAG_NONLINEAR     (UINT32_C(1) << 4U)

#define PHY_FRACTON_FLAG_BIJECTIVE  (UINT32_C(1) << 0U)
#define PHY_FRACTON_FLAG_LAPLACIAN  (UINT32_C(1) << 1U)
#define PHY_FRACTON_FLAG_ACTIVE     (UINT32_C(1) << 2U)
#define PHY_FRACTON_FLAG_STRONG_HSF (UINT32_C(1) << 3U)
#define PHY_FRACTON_FLAG_ISO_FROZEN (UINT32_C(1) << 4U)
#define PHY_FRACTON_FLAG_DIP_MOBILE (UINT32_C(1) << 5U)
#define PHY_FRACTON_FLAG_NONLINEAR  (UINT32_C(1) << 6U)

#define PHY_TETRA4_FLAG_BIJECTIVE   (UINT32_C(1) << 0U)
#define PHY_TETRA4_FLAG_TETRA       (UINT32_C(1) << 1U)
#define PHY_TETRA4_FLAG_NONINVOL    (UINT32_C(1) << 2U)
#define PHY_TETRA4_FLAG_3COUPLED    (UINT32_C(1) << 3U)
#define PHY_TETRA4_FLAG_NONLINEAR   (UINT32_C(1) << 4U)

#define PHY_HUNT_ENTANGLED          UINT32_C(0)
#define PHY_HUNT_YANGBOX            UINT32_C(1)
#define PHY_HUNT_JOINT              UINT32_C(2)
#define PHY_HUNT_Q4_HIERARCHY       UINT32_C(3)
#define PHY_HUNT_YB8                UINT32_C(4)
#define PHY_HUNT_DUAL4              UINT32_C(5)
#define PHY_HUNT_AME8               UINT32_C(6)
#define PHY_HUNT_FRACTON            UINT32_C(7)
#define PHY_HUNT_BIQUANDLE8         UINT32_C(8)
#define PHY_HUNT_TETRA4             UINT32_C(9)

typedef struct {
    uint8_t p[16];
    uint8_t pinv[16];
    int32_t q[16];
    int32_t j[16][16];
    uint64_t lut_packed;
    uint64_t pinv_packed;
    uint64_t merkle;
    uint64_t seed;
    uint32_t flags;
    uint32_t cycle_lcm;
    uint16_t j_nonzero;
    uint8_t anf_deg;
    uint8_t yb_violations;
    uint8_t rank_tot;
    uint8_t rank_ent;
    uint8_t is_bijective;
    uint8_t is_swap;
    uint8_t is_id;
    uint8_t is_factorizable;
    uint8_t nl_bit;
    uint8_t nl_a;
    uint8_t nl_b;
    uint8_t sep_l0;
    uint8_t sep_l1;
    uint8_t sep_r0;
    uint8_t sep_r1;
    uint8_t nf_axis;
    uint8_t nf_l;
    uint8_t nf_r0;
    uint8_t nf_r1;
    uint8_t yb_neq_swap_u;
    int32_t continuity_fail;
} phy_rule_analysis_t;

typedef struct {
    uint64_t state_lo;
    uint64_t state_hi;
    uint64_t pert_lo;
    uint64_t pert_hi;
    uint64_t merkle;
    uint64_t flow_state[8];
    int64_t total_charge;
    uint32_t step_idx;
    uint16_t damage_count;
    uint8_t front_left;
    uint8_t front_right;
} phy_fiber_t;

typedef struct {
    uint32_t steps_run;
    uint32_t max_span;
    int32_t vel_num_q16;
    int64_t initial_charge;
    int64_t final_charge;
    int64_t max_charge_drift;
    uint64_t final_merkle;
    uint32_t frames_recorded;
} phy_sim_trace_t;

typedef struct {
    uint32_t max_span;
    uint32_t mean_span_q16;
    int32_t max_vel_q16;
    int32_t mean_vel_q16;
    uint32_t active_probes;
    int64_t max_charge_drift;
} phy_butterfly_t;

typedef struct {
    uint64_t rref[255U * 255U];
    uint64_t triv_ech[255U * 255U];
    uint64_t work_row[255];
    uint64_t null_vec[255];
    uint8_t has_piv[255];
    uint8_t triv_has[255];
    uint8_t triv_base_ready;
    uint8_t f3[16][16][16];
} phy_q4_ws_t;

typedef struct {
    uint8_t p[16];
    int32_t q4[256];
    uint64_t lut_packed;
    uint64_t q4_merkle;
    uint16_t tot2;
    uint16_t tot4;
    uint16_t irred4;
    uint16_t genuine_q4;
    uint16_t j4_nonzero;
    uint16_t q4_nnz;
    int32_t continuity4_fail;
    uint8_t cross_axis;
    uint8_t cross_f1;
    uint8_t cross_f2;
    uint8_t cross_x0;
    uint8_t cross_x1;
    uint8_t cross_y0;
    uint8_t cross_y1;
    uint8_t j4_wit_a;
    uint8_t j4_wit_b;
    uint8_t j4_wit_c;
} phy_q4_analysis_t;

typedef struct {
    uint8_t r[64];
    uint8_t rinv[64];
    uint8_t sig_inv[8][8];
    uint8_t tau_inv[8][8];
    uint64_t merkle;
    uint32_t flags;
    uint32_t order;
    uint16_t yb8_violations;
    uint8_t anf_deg;
    uint8_t is_bijective;
    uint8_t is_nondeg;
    uint8_t is_involutive;
    uint8_t is_indecomposable;
    uint8_t is_lyubashenko;
    uint8_t noninvol_u;
    uint8_t nl_bit;
    uint8_t nl_a;
    uint8_t nl_b;
    uint8_t nly_axis;
    uint8_t nly_a;
    uint8_t nly_b0;
    uint8_t nly_b1;
} phy_yb8_analysis_t;

typedef struct {
    uint8_t p[16];
    uint8_t pinv[16];
    uint8_t sig_inv[4][4];
    uint8_t tau_inv[4][4];
    uint8_t m_plus[4][4];
    uint8_t m_minus[4][4];
    uint64_t lut_packed;
    uint64_t merkle;
    uint64_t mplus_packed;
    uint64_t mminus_packed;
    uint32_t flags;
    uint8_t anf_deg;
    uint8_t is_bijective;
    uint8_t is_nondeg;
    uint8_t is_dual_unitary;
    uint8_t is_mixing;
    uint8_t is_non_bernoulli;
    uint8_t yb_violations;
    uint8_t rank_ent;
    uint8_t nl_bit;
    uint8_t nl_a;
    uint8_t nl_b;
} phy_dual4_analysis_t;

typedef struct {
    uint8_t r[64];
    uint8_t rinv[64];
    uint8_t sig_inv[8][8];
    uint8_t tau_inv[8][8];
    uint8_t cross_l_inv[8][8];
    uint8_t cross_r_inv[8][8];
    uint64_t merkle;
    uint32_t flags;
    uint32_t order;
    uint8_t anf_deg;
    uint8_t is_bijective;
    uint8_t is_nondeg;
    uint8_t is_crossbij;
    uint8_t is_2unitary;
    uint8_t nl_bit;
    uint8_t nl_a;
    uint8_t nl_b;
} phy_ame8_analysis_t;

typedef struct {
    uint8_t vis[4096];
    uint16_t queue[4096];
    uint16_t sec_states[128];
} phy_fracton_hsf_ws_t;

typedef struct {
    uint8_t t[64];
    uint8_t tinv[64];
    uint64_t merkle;
    uint32_t flags;
    uint32_t order;
    int64_t max_q_drift;
    int64_t max_d_drift;
    uint16_t active_entries;
    uint16_t hsf_sector_size;
    uint16_t hsf_num_comps;
    uint16_t hsf_largest_comp;
    uint16_t hsf_frozen;
    uint16_t iso_span;
    uint16_t dip_span;
    uint8_t anf_deg;
    uint8_t is_bijective;
    uint8_t laplacian_ok;
    uint8_t cons_qd_4096_ok;
    uint8_t wit_u;
    uint8_t hsf_target_q;
    uint8_t hsf_target_d;
    uint8_t nl_bit;
    uint8_t nl_a;
    uint8_t nl_b;
} phy_fracton_analysis_t;

typedef struct {
    phy_yb8_analysis_t base;
    uint16_t col_unknot;
    uint16_t col_hopf;
    uint16_t col_trefoil;
    uint16_t col_fig8;
    uint16_t col_torus5;
    uint16_t col_torus7;
    uint8_t retract_classes;
    uint8_t is_r1;
    uint8_t is_irretract;
} phy_biquandle8_analysis_t;

typedef struct {
    uint8_t t[64];
    uint8_t tinv[64];
    uint64_t merkle;
    uint32_t flags;
    uint32_t order;
    uint16_t tetra4_viol;
    uint16_t tetra2_total;
    uint16_t tetra2_nl;
    uint16_t tetra2_noninvol_nl;
    uint8_t anf_deg;
    uint8_t is_bijective;
    uint8_t is_involutive;
    uint8_t is_3coupled;
    uint8_t noninvol_u;
    uint8_t nl_bit;
    uint8_t nl_a;
    uint8_t nl_b;
    uint8_t c0_a0;
    uint8_t c0_a1;
    uint8_t c0_b;
    uint8_t c0_c;
    uint8_t c2_a;
    uint8_t c2_b;
    uint8_t c2_c0;
    uint8_t c2_c1;
} phy_tetra4_analysis_t;

uint64_t phy_fnv1a64_u64(uint64_t h, uint64_t v);
uint64_t phy_m61_add(uint64_t a, uint64_t b);
uint64_t phy_m61_sub(uint64_t a, uint64_t b);
uint64_t phy_m61_mul(uint64_t a, uint64_t b);
uint64_t phy_m61_inv(uint64_t a);
uint64_t phy_stein_gcd64(uint64_t a, uint64_t b);
uint64_t phy_div_u64(uint64_t num, uint64_t den, uint64_t *rem);

uint64_t phy_pack_lut(const uint8_t p[16]);
void phy_unpack_lut(uint64_t packed, uint8_t p[16]);
uint8_t phy_step3(const uint8_t p[16], uint8_t a, uint8_t b, uint8_t c);
uint8_t phy_anf_degree4(const uint8_t p[16], uint8_t *out_bit, uint8_t *out_a, uint8_t *out_b);
uint8_t phy_yb_check(const uint8_t p[16]);
uint32_t phy_perm_order(const uint8_t p[16]);

void phy_analyze_rule(const uint8_t p[16], uint64_t seed, phy_rule_analysis_t *out);

void phy_step_64ring(const uint8_t p[16], uint64_t lo, uint64_t hi, uint64_t *out_lo,
                     uint64_t *out_hi);
int64_t phy_ring_charge_64(const int32_t q[16], uint64_t lo, uint64_t hi);
void phy_sim_soliton_64(phy_arena_u *arena, const phy_rule_analysis_t *rule, uint64_t init_lo,
                        uint64_t init_hi, uint32_t steps, uint8_t defect_site, uint8_t defect_xor,
                        phy_sim_trace_t *out_trace);
void phy_measure_butterfly_64(const phy_rule_analysis_t *rule, uint64_t init_lo, uint64_t init_hi,
                              uint32_t steps, phy_butterfly_t *out);

void phy_arena_init(phy_arena_u *arena, uint64_t seed);
uint32_t phy_hunt(phy_arena_u *arena, uint32_t mode, uint64_t seed, uint64_t max_trials,
                  phy_rule_analysis_t *out_best);

/* --- PHY v2.0: Exact 4-Site (Q4) Hierarchy & |X|=8 Non-Involutive YBE Engine --- */

void phy_analyze_q4(const uint8_t p[16], phy_q4_ws_t *ws, phy_q4_analysis_t *out);
int64_t phy_ring_charge4_64(const int32_t q4[256], uint64_t lo, uint64_t hi);
uint32_t phy_hunt_q4(phy_arena_u *arena, phy_q4_ws_t *ws, uint64_t seed, uint64_t max_trials,
                     phy_rule_analysis_t *out_rule, phy_q4_analysis_t *out_q4);

uint16_t phy_yb8_check(const uint8_t r[64]);
void phy_analyze_yb8(const uint8_t r[64], phy_yb8_analysis_t *out);
uint32_t phy_hunt_yb8(phy_arena_u *arena, uint64_t seed, uint64_t max_trials,
                      phy_yb8_analysis_t *out_best);

/* --- PHY v3.0: Dual-Unitary, AME(4,8), Fracton HSF, Biquandle8 Knots, Tetrahedron --- */

void phy_analyze_dual4(const uint8_t p[16], phy_dual4_analysis_t *out);
uint32_t phy_hunt_dual4(phy_arena_u *arena, uint64_t seed, uint64_t max_trials,
                        phy_dual4_analysis_t *out_best);

void phy_analyze_ame8(const uint8_t r[64], phy_ame8_analysis_t *out);
uint32_t phy_hunt_ame8(phy_arena_u *arena, uint64_t seed, uint64_t max_trials,
                       phy_ame8_analysis_t *out_best);

void phy_analyze_fracton(const uint8_t t[64], phy_fracton_hsf_ws_t *ws,
                         phy_fracton_analysis_t *out);
uint32_t phy_hunt_fracton(phy_arena_u *arena, phy_fracton_hsf_ws_t *ws, uint64_t seed,
                          uint64_t max_trials, phy_fracton_analysis_t *out_best);

void phy_analyze_biquandle8(const uint8_t r[64], phy_biquandle8_analysis_t *out);
uint32_t phy_hunt_biquandle8(phy_arena_u *arena, uint64_t seed, uint64_t max_trials,
                             phy_biquandle8_analysis_t *out_best);

void phy_classify_tetra2(uint16_t *out_total, uint16_t *out_nl, uint16_t *out_noninvol_nl);
uint16_t phy_tetra4_check(const uint8_t t[64]);
void phy_analyze_tetra4(const uint8_t t[64], phy_tetra4_analysis_t *out);
uint32_t phy_hunt_tetra4(phy_arena_u *arena, uint64_t seed, uint64_t max_trials,
                         phy_tetra4_analysis_t *out_best);

#ifdef __cplusplus
}
#endif

#endif /* PHY_H */
