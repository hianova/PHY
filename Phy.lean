import Flow.Schema
import Flow.Morphism

namespace Phy

open Flow.Schema

def hdrLayout : CellLayout64 :=
  CellLayout64.stdHeader64 "PhyHdrCell" "phy_hdr_cell_t" "active_rules"
    ["merkle_root", "seed", "flags", "reserved"]

def huntStatsLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyHuntStatsCell" "phy_hunt_stats_t" "stats"
    [ ("trials_done", .u64, 1)
    , ("bijective_tested", .u64, 1)
    , ("entangled_found", .u64, 1)
    , ("yangbaxter_found", .u64, 1)
    , ("joint_found", .u64, 1)
    , ("rejected_separable", .u64, 1)
    , ("rejected_linear", .u64, 1)
    , ("best_rule_merkle", .u64, 1)
    ]

def ruleLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyRuleCell" "phy_rule_cell_t" "rule"
    [ ("lut_packed", .u64, 1)
    , ("pinv_packed", .u64, 1)
    , ("merkle", .u64, 1)
    , ("seed", .u64, 1)
    , ("flags", .u32, 1)
    , ("anf_deg", .u8, 1)
    , ("yb_violations", .u8, 1)
    , ("rank_tot", .u8, 1)
    , ("rank_ent", .u8, 1)
    , ("j_nonzero", .u16, 1)
    , ("mc_and_count", .u16, 1)
    , ("cycle_lcm", .u32, 1)
    , ("soliton_vel_num", .i32, 1)
    , ("soliton_max_span", .u32, 1)
    , ("pad0", .u64, 1)
    ]

def chargeLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyChargeCell" "phy_charge_cell_t" "charge"
    [ ("q", .i32, 16)
    ]

def currentLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyCurrentCell" "phy_current_cell_t" "current"
    [ ("j", .i32, 16)
    ]

def simFrameLayout : CellLayout64 :=
  CellLayout64.packed64 "PhySimFrameCell" "phy_sim_frame_t" "frame"
    [ ("state_lo", .u64, 1)
    , ("state_hi", .u64, 1)
    , ("pert_lo", .u64, 1)
    , ("pert_hi", .u64, 1)
    , ("diff_mask", .u64, 1)
    , ("total_charge", .i64, 1)
    , ("step_idx", .u32, 1)
    , ("damage_count", .u16, 1)
    , ("front_left", .u8, 1)
    , ("front_right", .u8, 1)
    , ("jet_width_lo", .u32, 1)
    , ("frame_merkle", .u32, 1)
    ]

def auxLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyAuxCell" "phy_aux_cell_t" "aux"
    [ ("w", .u64, 8)
    ]

def q4MetaLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyQ4MetaCell" "phy_q4_meta_t" "q4_meta"
    [ ("lut_packed", .u64, 1)
    , ("q4_merkle", .u64, 1)
    , ("tot2", .u16, 1)
    , ("tot4", .u16, 1)
    , ("irred4", .u16, 1)
    , ("genuine_q4", .u16, 1)
    , ("j4_nonzero", .u16, 1)
    , ("q4_nnz", .u16, 1)
    , ("continuity4_fail", .i32, 1)
    , ("pad", .u64, 4)
    ]

def q4DensityLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyQ4DensityCell" "phy_q4_density_t" "q4_dens"
    [ ("q", .i32, 16)
    ]

def yb8RuleLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyYb8RuleCell" "phy_yb8_rule_t" "yb8_rule"
    [ ("r", .u8, 64)
    ]

def yb8MetaLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyYb8MetaCell" "phy_yb8_meta_t" "yb8_meta"
    [ ("merkle", .u64, 1)
    , ("flags", .u32, 1)
    , ("order", .u32, 1)
    , ("yb8_violations", .u16, 1)
    , ("mc_and_count", .u16, 1)
    , ("anf_deg", .u8, 1)
    , ("is_bijective", .u8, 1)
    , ("is_nondeg", .u8, 1)
    , ("is_involutive", .u8, 1)
    , ("is_indecomposable", .u8, 1)
    , ("is_lyubashenko", .u8, 1)
    , ("noninvol_u", .u8, 1)
    , ("pad8", .u8, 1)
    , ("pad32", .u32, 1)
    , ("pad64", .u64, 4)
    ]

def dual4MetaLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyDual4MetaCell" "phy_dual4_meta_t" "dual4_meta"
    [ ("lut_packed", .u64, 1)
    , ("merkle", .u64, 1)
    , ("mplus_packed", .u64, 1)
    , ("mminus_packed", .u64, 1)
    , ("flags", .u32, 1)
    , ("anf_deg", .u8, 1)
    , ("is_dual_unitary", .u8, 1)
    , ("is_mixing", .u8, 1)
    , ("is_non_bernoulli", .u8, 1)
    , ("yb_violations", .u8, 1)
    , ("rank_ent", .u8, 1)
    , ("pad16", .u16, 1)
    , ("pad32", .u32, 1)
    , ("pad64", .u64, 2)
    ]

def ame8MetaLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyAme8MetaCell" "phy_ame8_meta_t" "ame8_meta"
    [ ("merkle", .u64, 1)
    , ("flags", .u32, 1)
    , ("order", .u32, 1)
    , ("anf_deg", .u8, 1)
    , ("is_bijective", .u8, 1)
    , ("is_nondeg", .u8, 1)
    , ("is_crossbij", .u8, 1)
    , ("is_2unitary", .u8, 1)
    , ("nl_bit", .u8, 1)
    , ("nl_a", .u8, 1)
    , ("nl_b", .u8, 1)
    , ("pad64", .u64, 5)
    ]

def fractonMetaLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyFractonMetaCell" "phy_fracton_meta_t" "fracton_meta"
    [ ("merkle", .u64, 1)
    , ("flags", .u32, 1)
    , ("active_entries", .u16, 1)
    , ("hsf_sector_size", .u16, 1)
    , ("hsf_num_comps", .u16, 1)
    , ("hsf_largest_comp", .u16, 1)
    , ("hsf_frozen", .u16, 1)
    , ("iso_span", .u16, 1)
    , ("dip_span", .u16, 1)
    , ("anf_deg", .u8, 1)
    , ("is_bijective", .u8, 1)
    , ("laplacian_ok", .u8, 1)
    , ("wit_u", .u8, 1)
    , ("hsf_target_q", .u8, 1)
    , ("hsf_target_d", .u8, 1)
    , ("pad64", .u64, 4)
    ]

def biquandle8MetaLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyBiquandle8MetaCell" "phy_biquandle8_meta_t" "bq8_meta"
    [ ("merkle", .u64, 1)
    , ("flags", .u32, 1)
    , ("order", .u32, 1)
    , ("col_unknot", .u16, 1)
    , ("col_hopf", .u16, 1)
    , ("col_trefoil", .u16, 1)
    , ("col_fig8", .u16, 1)
    , ("col_torus5", .u16, 1)
    , ("col_torus7", .u16, 1)
    , ("anf_deg", .u8, 1)
    , ("retract_classes", .u8, 1)
    , ("is_r1", .u8, 1)
    , ("is_irretract", .u8, 1)
    , ("pad64", .u64, 4)
    ]

def tetra4MetaLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyTetra4MetaCell" "phy_tetra4_meta_t" "tetra4_meta"
    [ ("merkle", .u64, 1)
    , ("flags", .u32, 1)
    , ("order", .u32, 1)
    , ("tetra4_viol", .u16, 1)
    , ("tetra2_total", .u16, 1)
    , ("tetra2_nl", .u16, 1)
    , ("tetra2_noninvol_nl", .u16, 1)
    , ("anf_deg", .u8, 1)
    , ("is_bijective", .u8, 1)
    , ("is_involutive", .u8, 1)
    , ("is_3coupled", .u8, 1)
    , ("noninvol_u", .u8, 1)
    , ("pad8", .u8, 1)
    , ("pad16", .u16, 1)
    , ("pad64", .u64, 4)
    ]

def ame4NogoMetaLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyAme4NogoMetaCell" "phy_ame4_nogo_meta_t" "ame4_nogo"
    [ ("sample_lut", .u64, 1)
    , ("merkle", .u64, 1)
    , ("latin4_count", .u32, 1)
    , ("total_ols", .u32, 1)
    , ("affine_ols", .u32, 1)
    , ("nonlin_ols", .u32, 1)
    , ("bernoulli_ols", .u32, 1)
    , ("pad32", .u32, 1)
    , ("pad64", .u64, 3)
    ]

def simplex4MetaLayout : CellLayout64 :=
  CellLayout64.packed64 "PhySimplex4MetaCell" "phy_simplex4_meta_t" "simplex4_meta"
    [ ("merkle", .u64, 1)
    , ("flags", .u32, 1)
    , ("viol4", .u32, 1)
    , ("order4", .u32, 1)
    , ("viol2_ccx", .u16, 1)
    , ("viol2_cccx", .u16, 1)
    , ("anf_deg", .u8, 1)
    , ("is_bijective", .u8, 1)
    , ("is_involutive", .u8, 1)
    , ("is_4coupled", .u8, 1)
    , ("pad32", .u32, 1)
    , ("pad64", .u64, 4)
    ]

def chiral6MetaLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyChiral6MetaCell" "phy_chiral6_meta_t" "chiral6_meta"
    [ ("merkle", .u64, 1)
    , ("flags", .u32, 1)
    , ("cocycle_rank", .u16, 1)
    , ("unknot_w0", .u16, 1)
    , ("hopf_w0", .u16, 1)
    , ("hopf_w2", .u16, 1)
    , ("tref_l_w1", .u16, 1)
    , ("tref_r_w3", .u16, 1)
    , ("fig8_w0", .u16, 1)
    , ("is_chiral", .u8, 1)
    , ("is_non_coboundary", .u8, 1)
    , ("pad32", .u32, 1)
    , ("pad64", .u64, 4)
    ]

def qscarMetaLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyQscarMetaCell" "phy_qscar_meta_t" "qscar_meta"
    [ ("merkle", .u64, 1)
    , ("du_purity_t2_num", .u64, 1)
    , ("du_purity_t2_den", .u64, 1)
    , ("therm_purity_t4_num", .u64, 1)
    , ("therm_purity_t4_den", .u64, 1)
    , ("flags", .u32, 1)
    , ("scar_period", .u16, 1)
    , ("therm_max_sup", .u16, 1)
    , ("du_unitary", .u8, 1)
    , ("du_dual", .u8, 1)
    , ("du_full_sup", .u8, 1)
    , ("pad8", .u8, 5)
    , ("pad64", .u64, 1)
    ]

def fracCompMetaLayout : CellLayout64 :=
  CellLayout64.packed64 "PhyFracCompMetaCell" "phy_frac_comp_meta_t" "frac_comp_meta"
    [ ("merkle", .u64, 1)
    , ("flags", .u32, 1)
    , ("sec6", .u16, 1)
    , ("sec7", .u16, 1)
    , ("sec8", .u16, 1)
    , ("comp6", .u16, 1)
    , ("comp7", .u16, 1)
    , ("comp8", .u16, 1)
    , ("largest6", .u16, 1)
    , ("largest7", .u16, 1)
    , ("largest8", .u16, 1)
    , ("frozen6", .u16, 1)
    , ("frozen7", .u16, 1)
    , ("frozen8", .u16, 1)
    , ("ccx_ok", .u8, 1)
    , ("laplacian_ok", .u8, 1)
    , ("is_involutive", .u8, 1)
    , ("scaling_strict", .u8, 1)
    , ("pad64", .u64, 3)
    ]

def phyArenaSchema : ArenaSpec :=
  { cPrefix := "PHY"
  , headerGuard := "PHY_SCHEMA_ABI_H"
  , totalCells := 512
  , cCellUnionName := "phy_cell_u"
  , cArenaMapName := "phy_arena_map"
  , cArenaUnionName := "phy_arena_u"
  , regions :=
    [ { name := "img_hdr"
      , baseSlot := 0
      , count := 1
      , cSlotMacro := "PHY_SLOT_HDR"
      , cCountMacro := "PHY_COUNT_HDR"
      , mapDecl := "phy_hdr_cell_t img_hdr;"
      }
    , { name := "hunt_stats"
      , baseSlot := 1
      , count := 1
      , cSlotMacro := "PHY_SLOT_STATS"
      , cCountMacro := "PHY_COUNT_STATS"
      , mapDecl := "phy_hunt_stats_t hunt_stats;"
      }
    , { name := "rules"
      , baseSlot := 2
      , count := 16
      , cSlotMacro := "PHY_SLOT_RULES"
      , cCountMacro := "PHY_COUNT_RULES"
      , mapDecl := "phy_rule_cell_t rules[16];"
      }
    , { name := "charges"
      , baseSlot := 18
      , count := 16
      , cSlotMacro := "PHY_SLOT_CHARGES"
      , cCountMacro := "PHY_COUNT_CHARGES"
      , mapDecl := "phy_charge_cell_t charges[16];"
      }
    , { name := "currents"
      , baseSlot := 34
      , count := 128
      , cSlotMacro := "PHY_SLOT_CURRENTS"
      , cCountMacro := "PHY_COUNT_CURRENTS"
      , mapDecl := "phy_current_cell_t currents[128];"
      }
    , { name := "sim_frames"
      , baseSlot := 162
      , count := 32
      , cSlotMacro := "PHY_SLOT_SIM_FRAMES"
      , cCountMacro := "PHY_COUNT_SIM_FRAMES"
      , mapDecl := "phy_sim_frame_t sim_frames[32];"
      }
    , { name := "flow_aux"
      , baseSlot := 194
      , count := 62
      , cSlotMacro := "PHY_SLOT_FLOW_AUX"
      , cCountMacro := "PHY_COUNT_FLOW_AUX"
      , mapDecl := "phy_aux_cell_t flow_aux[62];"
      }
    , { name := "q4_meta"
      , baseSlot := 256
      , count := 16
      , cSlotMacro := "PHY_SLOT_Q4_META"
      , cCountMacro := "PHY_COUNT_Q4_META"
      , mapDecl := "phy_q4_meta_t q4_meta[16];"
      }
    , { name := "q4_densities"
      , baseSlot := 272
      , count := 64
      , cSlotMacro := "PHY_SLOT_Q4_DENSITIES"
      , cCountMacro := "PHY_COUNT_Q4_DENSITIES"
      , mapDecl := "phy_q4_density_t q4_densities[64];"
      }
    , { name := "yb8_rules"
      , baseSlot := 336
      , count := 16
      , cSlotMacro := "PHY_SLOT_YB8_RULES"
      , cCountMacro := "PHY_COUNT_YB8_RULES"
      , mapDecl := "phy_yb8_rule_t yb8_rules[16];"
      }
    , { name := "yb8_meta"
      , baseSlot := 352
      , count := 16
      , cSlotMacro := "PHY_SLOT_YB8_META"
      , cCountMacro := "PHY_COUNT_YB8_META"
      , mapDecl := "phy_yb8_meta_t yb8_meta[16];"
      }
    , { name := "dual4_meta"
      , baseSlot := 368
      , count := 8
      , cSlotMacro := "PHY_SLOT_DUAL4_META"
      , cCountMacro := "PHY_COUNT_DUAL4_META"
      , mapDecl := "phy_dual4_meta_t dual4_meta[8];"
      }
    , { name := "ame8_meta"
      , baseSlot := 376
      , count := 8
      , cSlotMacro := "PHY_SLOT_AME8_META"
      , cCountMacro := "PHY_COUNT_AME8_META"
      , mapDecl := "phy_ame8_meta_t ame8_meta[8];"
      }
    , { name := "fracton_meta"
      , baseSlot := 384
      , count := 8
      , cSlotMacro := "PHY_SLOT_FRACTON_META"
      , cCountMacro := "PHY_COUNT_FRACTON_META"
      , mapDecl := "phy_fracton_meta_t fracton_meta[8];"
      }
    , { name := "bq8_meta"
      , baseSlot := 392
      , count := 8
      , cSlotMacro := "PHY_SLOT_BQ8_META"
      , cCountMacro := "PHY_COUNT_BQ8_META"
      , mapDecl := "phy_biquandle8_meta_t bq8_meta[8];"
      }
    , { name := "tetra4_meta"
      , baseSlot := 400
      , count := 8
      , cSlotMacro := "PHY_SLOT_TETRA4_META"
      , cCountMacro := "PHY_COUNT_TETRA4_META"
      , mapDecl := "phy_tetra4_meta_t tetra4_meta[8];"
      }
    , { name := "v3_rules6"
      , baseSlot := 408
      , count := 16
      , cSlotMacro := "PHY_SLOT_V3_RULES6"
      , cCountMacro := "PHY_COUNT_V3_RULES6"
      , mapDecl := "phy_yb8_rule_t v3_rules6[16];"
      }
    , { name := "ame4_nogo"
      , baseSlot := 424
      , count := 1
      , cSlotMacro := "PHY_SLOT_AME4_NOGO"
      , cCountMacro := "PHY_COUNT_AME4_NOGO"
      , mapDecl := "phy_ame4_nogo_meta_t ame4_nogo;"
      }
    , { name := "simplex4_meta"
      , baseSlot := 425
      , count := 8
      , cSlotMacro := "PHY_SLOT_SIMPLEX4_META"
      , cCountMacro := "PHY_COUNT_SIMPLEX4_META"
      , mapDecl := "phy_simplex4_meta_t simplex4_meta[8];"
      }
    , { name := "chiral6_meta"
      , baseSlot := 433
      , count := 4
      , cSlotMacro := "PHY_SLOT_CHIRAL6_META"
      , cCountMacro := "PHY_COUNT_CHIRAL6_META"
      , mapDecl := "phy_chiral6_meta_t chiral6_meta[4];"
      }
    , { name := "qscar_meta"
      , baseSlot := 437
      , count := 4
      , cSlotMacro := "PHY_SLOT_QSCAR_META"
      , cCountMacro := "PHY_COUNT_QSCAR_META"
      , mapDecl := "phy_qscar_meta_t qscar_meta[4];"
      }
    , { name := "frac_comp_meta"
      , baseSlot := 441
      , count := 4
      , cSlotMacro := "PHY_SLOT_FRAC_COMP_META"
      , cCountMacro := "PHY_COUNT_FRAC_COMP_META"
      , mapDecl := "phy_frac_comp_meta_t frac_comp_meta[4];"
      }
    , { name := "v4_aux"
      , baseSlot := 445
      , count := 67
      , cSlotMacro := "PHY_SLOT_V4_AUX"
      , cCountMacro := "PHY_COUNT_V4_AUX"
      , mapDecl := "phy_aux_cell_t v4_aux[67];"
      }
    ]
  , layouts :=
    [ hdrLayout
    , huntStatsLayout
    , ruleLayout
    , chargeLayout
    , currentLayout
    , simFrameLayout
    , auxLayout
    , q4MetaLayout
    , q4DensityLayout
    , yb8RuleLayout
    , yb8MetaLayout
    , dual4MetaLayout
    , ame8MetaLayout
    , fractonMetaLayout
    , biquandle8MetaLayout
    , tetra4MetaLayout
    , ame4NogoMetaLayout
    , simplex4MetaLayout
    , chiral6MetaLayout
    , qscarMetaLayout
    , fracCompMetaLayout
    ]
  }

def emitSchemaOnly (outPath : System.FilePath) : IO Unit := do
  let ok ← ArenaSpec.emitC17HeaderChecked outPath phyArenaSchema
  if !ok then
    throw <| IO.userError s!"Phy: ArenaSpec verification failed emitting {outPath}"

end Phy

def main (args : List String) : IO Unit := do
  match args with
  | ["emit-schema", outPath] => Phy.emitSchemaOnly outPath
  | _ => Phy.emitSchemaOnly "build/phy_schema_abi.h"
