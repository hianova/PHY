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
#define PHY_VERSION                 UINT32_C(0x00040000) /* v4.0.0 */
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

#define PHY_S4_FLAG_BIJECTIVE       (UINT32_C(1) << 0U)
#define PHY_S4_FLAG_4SIMPLEX        (UINT32_C(1) << 1U)
#define PHY_S4_FLAG_NONINVOL        (UINT32_C(1) << 2U)
#define PHY_S4_FLAG_4COUPLED        (UINT32_C(1) << 3U)
#define PHY_S4_FLAG_NONLINEAR       (UINT32_C(1) << 4U)

#define PHY_CHIRAL6_FLAG_YBE        (UINT32_C(1) << 0U)
#define PHY_CHIRAL6_FLAG_R1         (UINT32_C(1) << 1U)
#define PHY_CHIRAL6_FLAG_COCYCLE    (UINT32_C(1) << 2U)
#define PHY_CHIRAL6_FLAG_CHIRAL     (UINT32_C(1) << 3U)
#define PHY_CHIRAL6_FLAG_NONCOBOUND (UINT32_C(1) << 4U)

#define PHY_QSCAR_FLAG_UNITARY      (UINT32_C(1) << 0U)
#define PHY_QSCAR_FLAG_DUAL         (UINT32_C(1) << 1U)
#define PHY_QSCAR_FLAG_FULL_SUP     (UINT32_C(1) << 2U)
#define PHY_QSCAR_FLAG_DTC_P2       (UINT32_C(1) << 3U)
#define PHY_QSCAR_FLAG_THERMAL      (UINT32_C(1) << 4U)

#define PHY_FCOMP_FLAG_BIJECTIVE    (UINT32_C(1) << 0U)
#define PHY_FCOMP_FLAG_INVOLUTIVE   (UINT32_C(1) << 1U)
#define PHY_FCOMP_FLAG_LAPLACIAN    (UINT32_C(1) << 2U)
#define PHY_FCOMP_FLAG_CCX          (UINT32_C(1) << 3U)
#define PHY_FCOMP_FLAG_SCALING      (UINT32_C(1) << 4U)

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
#define PHY_HUNT_SIMPLEX4           UINT32_C(10)
#define PHY_HUNT_CHIRAL6            UINT32_C(11)
#define PHY_HUNT_QSCAR              UINT32_C(12)
#define PHY_HUNT_FRACTON_COMP       UINT32_C(13)
#define PHY_HUNT_AME4_NOGO          UINT32_C(14)

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

/* --- PHY v4.0 Analysis Types --- */

typedef struct {
    uint8_t sample_p[16];
    uint64_t sample_lut;
    uint64_t merkle;
    uint32_t latin4_count;
    uint32_t total_ols;
    uint32_t affine_ols;
    uint32_t nonlin_ols;
    uint32_t bernoulli_ols;
} phy_ame4_nogo_t;

typedef struct {
    uint8_t ccx02[16];
    uint8_t ccx13[16];
    uint8_t cccx1[16];
    uint16_t viol2_ccx02;
    uint16_t viol2_ccx13;
    uint16_t viol2_cccx1;
    uint16_t cccx1_wit_s;
    uint8_t nl_count_b_ctrl;
} phy_simplex4_2_t;

typedef struct {
    uint8_t s[256];
    uint8_t sinv[256];
    uint8_t s_coupled[256];
    uint8_t ccx02[16];
    uint64_t merkle;
    uint32_t flags;
    uint32_t viol4;
    uint32_t viol4_coupled;
    uint32_t order4;
    uint16_t viol2_ccx;
    uint16_t viol2_cccx;
    uint8_t anf_deg;
    uint8_t is_bijective;
    uint8_t is_involutive;
    uint8_t is_4coupled;
    uint8_t noninvol_u;
} phy_simplex4_4_t;

typedef struct {
    uint8_t act[6][6];
    uint8_t r[36];
    uint8_t rinv[36];
    uint8_t phi[36];
    uint64_t merkle;
    uint32_t flags;
    uint16_t cocycle_viol;
    uint16_t cocycle_rank;
    uint16_t unknot_w[4];
    uint16_t hopf_w[4];
    uint16_t tref_l_w[4];
    uint16_t tref_r_w[4];
    uint16_t fig8_w[4];
    uint8_t is_chiral;
    uint8_t is_non_coboundary;
    uint8_t wit_a;
    uint8_t wit_b;
} phy_chiral6_analysis_t;

typedef struct {
    int8_t w_re[4][4];
    int8_t w_im[4][4];
    uint64_t merkle;
    uint64_t du_purity_num[4];
    uint64_t du_purity_den[4];
    uint64_t therm_purity_num[4];
    uint64_t therm_purity_den[4];
    uint64_t neel_ret_num[4];
    uint64_t neel_ret_den[4];
    uint32_t flags;
    uint16_t du_sup[4];
    uint16_t neel_sup[4];
    uint16_t therm_sup[4];
    uint16_t scar_period;
    uint8_t du_unitary;
    uint8_t du_dual;
    uint8_t du_full_sup;
    uint8_t scar_unitary;
} phy_qscar_analysis_t;

typedef struct {
    uint8_t vis[65536];
    uint32_t queue[2048];
    uint32_t sec_states[2048];
} phy_fracton_scale_ws_t;

typedef struct {
    uint8_t t_sample[64];
    uint8_t t_toff[64];
    uint8_t t_toff_inv[64];
    uint64_t merkle;
    uint32_t flags;
    uint16_t sec[3];
    uint16_t comps[3];
    uint16_t largest[3];
    uint16_t frozen[3];
    uint8_t is_involutive;
    uint8_t ccx_ok;
    uint8_t laplacian_ok;
    uint8_t scaling_strict;
} phy_fracton_comp_t;

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

/* --- PHY v4.0: Exact Rigidity, 3+1D 4-Simplex, Chiral Cohomology, Z[i] Scar, Fracton Univ --- */

bool phy_is_2unitary4(const uint8_t p[16]);
void phy_classify_ame4(phy_arena_u *arena, phy_ame4_nogo_t *out);

uint16_t phy_simplex4_check2(const uint8_t s[16], uint16_t *out_first_bad);
void phy_analyze_simplex4_2(phy_simplex4_2_t *out);
uint32_t phy_simplex4_check4(const uint8_t s[256]);
void phy_analyze_simplex4_4(const uint8_t s[256], phy_simplex4_4_t *out);
uint32_t phy_hunt_simplex4(phy_arena_u *arena, uint64_t seed, uint64_t max_trials,
                           phy_simplex4_4_t *out_best);

uint16_t phy_verify_q6_cocycle(const uint8_t r[36], const uint8_t phi[36]);
void phy_analyze_chiral6_custom(const uint8_t act[6][6], const uint8_t phi_in[36],
                                phy_chiral6_analysis_t *out);
uint32_t phy_solve_chiral6(phy_arena_u *arena, phy_chiral6_analysis_t *out);

void phy_verify_zi_du4(const int8_t re[4][4], const int8_t im[4][4], uint8_t *out_unitary,
                       uint8_t *out_dual, uint8_t *out_full_sup);
uint32_t phy_analyze_qscar_custom(phy_arena_u *arena, const int8_t re[4][4], const int8_t im[4][4],
                                  bool invert_scar_ctrl, phy_qscar_analysis_t *out);
uint32_t phy_analyze_qscar(phy_arena_u *arena, phy_qscar_analysis_t *out);

void phy_build_dipole_toffoli(uint8_t t[64]);
uint32_t phy_analyze_fracton_comp(phy_arena_u *arena, const uint8_t t_sample[64],
                                  const uint8_t t_toff[64], phy_fracton_scale_ws_t *ws,
                                  phy_fracton_comp_t *out);

#ifdef __cplusplus
}
#endif

#endif /* PHY_H */
