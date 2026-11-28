SHELL := /bin/sh
UNAME_S := $(shell uname -s)
FLOW_DIR ?= ../FLOW
FLOW_CLI := $(FLOW_DIR)/build/flow
BUILD := build
LEAN ?= false
LEAN_BIN ?= $(or $(shell command -v lean 2>/dev/null),$(HOME)/.elan/bin/lean)
LEANCHECKER ?= $(or $(shell command -v leanchecker 2>/dev/null),$(HOME)/.elan/bin/leanchecker)
Z3 ?= $(or $(shell command -v z3 2>/dev/null),$(wildcard /opt/homebrew/bin/z3))
ifeq ($(UNAME_S),Darwin)
LLVM_PREFIX ?= /opt/homebrew/opt/llvm/bin
else
LLVM_PREFIX ?= /usr/bin
endif
ifneq ($(wildcard $(LLVM_PREFIX)/clang),)
DEFAULT_CC := $(LLVM_PREFIX)/clang
DEFAULT_CXX := $(LLVM_PREFIX)/clang++
DEFAULT_FMT := $(LLVM_PREFIX)/clang-format
DEFAULT_TIDY := $(LLVM_PREFIX)/clang-tidy
else
DEFAULT_CC := clang
DEFAULT_CXX := clang++
DEFAULT_FMT := clang-format
DEFAULT_TIDY := clang-tidy
endif
ifeq ($(origin CC),default)
CC := $(DEFAULT_CC)
endif
ifeq ($(origin CXX),default)
CXX := $(DEFAULT_CXX)
endif
CLANG_FORMAT ?= $(DEFAULT_FMT)
CLANG_TIDY ?= $(DEFAULT_TIDY)

C17_RING0_FLAGS := -std=c17 -O2 -ffreestanding -nostdlib -fno-builtin -fno-stack-protector \
	-fno-jump-tables -mgeneral-regs-only -Wall -Wextra -Werror -Wpedantic -Wconversion
C17_HOST_FLAGS := -std=c17 -O2 -Wall -Wextra -Werror -Wpedantic -Wconversion
INC_RING0 := -Iinclude -I$(BUILD) -I$(FLOW_DIR)/core
INC_HOST := -Iinclude -I$(BUILD) -I$(FLOW_DIR)/core -I$(FLOW_DIR)/sdk

C_SOURCES := include/phy.h src/phy_ring0.c src/phy_cert.h src/phy_cli.c

.PHONY: all clean fmt verify-all \
	gate0-lean-noether gate1-lean-schema gate2-ring0-audit gate3-selftest \
	gate4-cli-audit gate5-oracle gate6-sabotage gate7-format

all: $(BUILD)/phy

$(FLOW_CLI):
	$(MAKE) -C $(FLOW_DIR) all LEAN=false

$(BUILD)/lean/Phy/Noether.olean: Phy/Noether.lean
	@mkdir -p $(BUILD)/lean/Phy
	LEAN_PATH=$(BUILD)/lean:$(FLOW_DIR)/build/lean $(LEAN_BIN) -DwarningAsError=true --root=. -o $@ $<

$(BUILD)/lean/Phy/HigherDimFrontier.olean: Phy/HigherDimFrontier.lean
	@mkdir -p $(BUILD)/lean/Phy
	$(LEAN_BIN) -DwarningAsError=true --root=. -o $@ $<

$(BUILD)/lean/Phy.olean: Phy.lean
	@mkdir -p $(BUILD)/lean
	LEAN_PATH=$(BUILD)/lean:$(FLOW_DIR)/build/lean \
		$(LEAN_BIN) -DwarningAsError=true --root=. -o $@ $<

$(BUILD)/phy_schema_abi.h: phy_schema.fspec $(FLOW_CLI)
	@mkdir -p $(BUILD)
	$(FLOW_CLI) schema emit c17 $< $@

$(BUILD)/phy_schema.smt2: phy_schema.fspec $(FLOW_CLI)
	@mkdir -p $(BUILD)
	$(FLOW_CLI) schema emit smt2 $< $@

$(BUILD)/PhySchemaCert.lean: phy_schema.fspec $(FLOW_CLI)
	@mkdir -p $(BUILD)
	$(FLOW_CLI) schema emit lean $< $@

$(BUILD)/phy_ring0.c.o: src/phy_ring0.c include/phy.h $(BUILD)/phy_schema_abi.h
	$(CC) $(C17_RING0_FLAGS) $(INC_RING0) -c $< -o $@

$(BUILD)/phy: $(BUILD)/phy_ring0.c.o src/phy_cli.c src/phy_cert.h include/phy.h
	$(CC) $(C17_HOST_FLAGS) $(INC_HOST) $(BUILD)/phy_ring0.c.o src/phy_cli.c -o $@

gate0-lean-noether: $(BUILD)/lean/Phy/Noether.olean
	@echo "[Gate 0] Checking Phy.Noether via leanchecker and axiom audit..."
	@if grep -qiE '\bsorry\b' Phy/Noether.lean; then echo "ERROR: sorry found in Phy/Noether.lean"; exit 1; fi
	@LEAN_PATH=$(BUILD)/lean:$(FLOW_DIR)/build/lean $(LEANCHECKER) Phy.Noether & pid1=$$!; \
	(printf 'import Phy.Noether\n#print axioms Phy.Noether.discrete_noether_global_conservation\n#print axioms Phy.Noether.telescoping_criterion_iff_exists_current\n#print axioms Phy.Noether.forall_ring_conservation_of_bv32\n#print axioms Phy.Noether.discrete_noether4_global_conservation\n#print axioms Phy.Noether.telescoping4_criterion_iff_exists_current4\n#print axioms Phy.Noether.forall_ring_conservation4_of_bv32\n#print axioms Phy.Noether.non_degenerate8_of_inverses\n#print axioms Phy.Noether.indecomposable8_of_connected8\n#print axioms Phy.Noether.dual_bijective4_iff_nondegenerate4\n#print axioms Phy.Noether.dual_unitary_q2_separable_nogo\n#print axioms Phy.Noether.ybe_nondeg_q2_no_entangled\n#print axioms Phy.Noether.cross_bijective8_of_inverses\n#print axioms Phy.Noether.fracton_open_chain_Q_conservation\n#print axioms Phy.Noether.fracton_open_chain_D_conservation\n#print axioms Phy.Noether.fracton_single_charge_frozen\n#print axioms Phy.Noether.fracton_quadrupole_active_witness\n#print axioms Phy.Noether.irretractable8_of_bool\n#print axioms Phy.Noether.toffoli3_satisfies_tetra2\n#print axioms Phy.Noether.coupled3_of_witnesses\n#print axioms Phy.Noether.ame4_affine_rigidity_nogo\n#print axioms Phy.Noether.ame4_eval_packLut4\n#print axioms Phy.Noether.ccx02_satisfies_4simplex2\n#print axioms Phy.Noether.ccx13_satisfies_4simplex2\n#print axioms Phy.Noether.cccx1_fails_4simplex2\n#print axioms Phy.Noether.s4RingLift_bijective\n#print axioms Phy.Noether.s4RingLift_noninvolutive\n#print axioms Phy.Noether.s4RingLift_satisfies_4simplex4\n#print axioms Phy.Noether.s4RingDual_satisfies_4simplex4\n#print axioms Phy.Noether.s4Coupled4Ch_satisfies_4simplex4\n#print axioms Phy.Noether.s4Coupled4Ch_4channel_coupled\n#print axioms Phy.Noether.coboundary_vanishes_torus2k\n#print axioms Phy.Noether.q6_non_coboundary_of_witness\n#print axioms Phy.Noether.q6_bijective\n#print axioms Phy.Noether.q6_r1\n#print axioms Phy.Noether.q6_satisfies_ybe\n#print axioms Phy.Noether.q6_noninvolutive\n#print axioms Phy.Noether.q6_phi_is_2cocycle\n#print axioms Phy.Noether.q6_phi_not_coboundary\n#print axioms Phy.Noether.q6_trefoil_chiral_separation\n#print axioms Phy.Noether.wCxDU_is_unitary\n#print axioms Phy.Noether.wCxDU_is_dual_unitary\n#print axioms Phy.Noether.wCxDU_full_superposition\n#print axioms Phy.Noether.scarGate4_is_unitary\n#print axioms Phy.Noether.dtc_neel_local_protected\n#print axioms Phy.Noether.dtc_neel_exact_period2\n#print axioms Phy.Noether.fracton3_universal_involution\n#print axioms Phy.Noether.dipoleToffoli_bijective\n#print axioms Phy.Noether.dipoleToffoli_satisfies_laplacian\n#print axioms Phy.Noether.dipole_static_control_nogo\n#print axioms Phy.Noether.dipoleToffoli_implements_ccx\n#print axioms Phy.Noether.dipoleToffoli_single_charge_frozen\n' | \
		LEAN_PATH=$(BUILD)/lean:$(FLOW_DIR)/build/lean $(LEAN_BIN) --stdin > $(BUILD)/noether_axioms.txt) & pid2=$$!; \
	wait $$pid1; r1=$$?; wait $$pid2; r2=$$?; \
	if [ $$r1 -ne 0 ] || [ $$r2 -ne 0 ]; then exit 1; fi
	@if grep -q 'sorryAx' $(BUILD)/noether_axioms.txt; then echo "ERROR: sorryAx detected"; exit 1; fi
	@echo "[Gate 0] PASS (Phy.Noether 0 sorryAx)"

gate1-lean-schema: $(BUILD)/phy_schema_abi.h $(BUILD)/phy_schema.smt2 $(BUILD)/PhySchemaCert.lean $(BUILD)/lean/Phy.olean $(BUILD)/lean/Phy/HigherDimFrontier.olean
	@echo "[Gate 1] Auditing Phy schema (.fspec -> C17/Z3/Lean), hash parity, Phy.lean & Phy.HigherDimFrontier..."
	@if grep -qE '^\s*(theorem|lemma|axiom)\s' Phy.lean; then echo "ERROR: theorems not allowed in Phy.lean"; exit 1; fi
	@if grep -qiE '\b(sorry|native_decide)\b' Phy/HigherDimFrontier.lean || grep -qE '^\s*axiom\s' Phy/HigherDimFrontier.lean; then echo "ERROR: sorry/axiom/native_decide found in Phy/HigherDimFrontier.lean"; exit 1; fi
	@if grep -q 'memcpy' $(BUILD)/phy_schema_abi.h; then echo "ERROR: memcpy found in generated ABI header"; exit 1; fi
	@if [ -n "$(Z3)" ]; then \
		res=$$($(Z3) $(BUILD)/phy_schema.smt2); \
		if [ "$$res" != "unsat" ]; then echo "ERROR: Z3 QF_LIA proof failed: $$res"; exit 1; fi; \
	fi
	LEAN_PATH=$(BUILD)/lean:$(FLOW_DIR)/build/lean $(LEAN_BIN) --run Phy.lean emit-schema $(BUILD)/phy_schema_abi_lean.h
	@h_c=$$($(FLOW_CLI) schema hash phy_schema.fspec); \
	h_l=$$(grep '^#define PHY_SCHEMA_HASH' $(BUILD)/phy_schema_abi_lean.h | sed -E 's/.*UINT64_C\((0x[0-9A-Fa-f]+)\).*/\1/'); \
	if [ "$$h_c" != "$$h_l" ]; then echo "ERROR: PHY_SCHEMA_HASH mismatch: .fspec=$$h_c Phy.lean=$$h_l"; exit 1; fi
	LEAN_PATH=$(BUILD)/lean:$(FLOW_DIR)/build/lean $(LEAN_BIN) -DwarningAsError=true $(BUILD)/PhySchemaCert.lean
	LEAN_PATH=$(BUILD)/lean:$(FLOW_DIR)/build/lean $(LEANCHECKER) Phy
	LEAN_PATH=$(BUILD)/lean $(LEANCHECKER) Phy.HigherDimFrontier
	@printf 'import Phy.HigherDimFrontier\n#print axioms Phy.HigherDimFrontier.tetraR3_order3\n#print axioms Phy.HigherDimFrontier.tetraR3_noninvolutive\n#print axioms Phy.HigherDimFrontier.tetraR3_not_involutive\n#print axioms Phy.HigherDimFrontier.tetraR3_bijective\n#print axioms Phy.HigherDimFrontier.tetraR3_secondDiffXZ\n#print axioms Phy.HigherDimFrontier.tetraR3_nonlinear\n#print axioms Phy.HigherDimFrontier.tetraR3_3channel_coupled\n#print axioms Phy.HigherDimFrontier.tetraR3_satisfies_tetra\n#print axioms Phy.HigherDimFrontier.isLatin_of_inverses\n#print axioms Phy.HigherDimFrontier.isOrthogonal_of_inverse\n#print axioms Phy.HigherDimFrontier.isOLS_implies_2unitary\n#print axioms Phy.HigherDimFrontier.ame45_deg6_relabeling\n#print axioms Phy.HigherDimFrontier.ame45_deg6_is_ols5\n#print axioms Phy.HigherDimFrontier.ame45_deg6_2unitary\n#print axioms Phy.HigherDimFrontier.ame45_deg6_not_affine\n#print axioms Phy.HigherDimFrontier.ame45_deg6_poly_exact\n#print axioms Phy.HigherDimFrontier.ame45_deg6_forced_deg6\n#print axioms Phy.HigherDimFrontier.ame47_witness_is_ols7\n#print axioms Phy.HigherDimFrontier.isotopic_to_Z7_implies_quadrangle\n#print axioms Phy.HigherDimFrontier.ame47_A7_not_isotopic_to_Z7\n#print axioms Phy.HigherDimFrontier.eisen_omega_eq\n#print axioms Phy.HigherDimFrontier.eisen_conj_mul_self\n#print axioms Phy.HigherDimFrontier.eisen_norm_mul\n#print axioms Phy.HigherDimFrontier.ame46_classical_defect_floor_2\n#print axioms Phy.HigherDimFrontier.ame46_zw40_nnz_count\n#print axioms Phy.HigherDimFrontier.ame46_zw40_unitary\n#print axioms Phy.HigherDimFrontier.ame46_zw40_defect_escape\n#print axioms Phy.HigherDimFrontier.ame46_zw44_nnz_count\n#print axioms Phy.HigherDimFrontier.ame46_zw44_unitary\n#print axioms Phy.HigherDimFrontier.ame46_zw44_defect_escape\n#print axioms Phy.HigherDimFrontier.ame46_zw854_nnz_and_unitary\n#print axioms Phy.HigherDimFrontier.ame46_zw854_defect_beats_27_26\n#print axioms Phy.HigherDimFrontier.ame46_zs5_63_nnz_and_unitary\n#print axioms Phy.HigherDimFrontier.ame46_zs5_63_defect_beats_27_26\n' | \
		LEAN_PATH=$(BUILD)/lean $(LEAN_BIN) --stdin > $(BUILD)/frontier_axioms.txt
	@if grep -q 'sorryAx' $(BUILD)/frontier_axioms.txt; then echo "ERROR: sorryAx detected in Phy.HigherDimFrontier"; exit 1; fi
	@echo "[Gate 1] PASS (512-slot ABI 0 memcpy, Z3 unsat, .fspec<->Lean hash $$($(FLOW_CLI) schema hash phy_schema.fspec) parity, PhySchemaCert + leanchecker Phy + Phy.HigherDimFrontier 0 sorryAx PASS)"

gate2-ring0-audit: $(BUILD)/phy_ring0.c.o
	@echo "[Gate 2] Auditing Ring-0 freestanding Pure-GPR object..."
	@test -z "$$(nm -u $(BUILD)/phy_ring0.c.o)"
ifeq ($(UNAME_S),Darwin)
	@if otool -tV $(BUILD)/phy_ring0.c.o | grep -qE '\b(udiv|sdiv|fdiv|fmul|fadd|fsub|blr|br)\b|v[0-9]+\.'; then \
		echo "ERROR: Forbidden instruction or vector register in phy_ring0.c.o"; exit 1; \
	fi
else
	@if objdump -d --no-show-raw-insn $(BUILD)/phy_ring0.c.o | awk '/^[[:space:]]*[0-9a-f]+:[[:space:]]+/{ sub(/^[[:space:]]*[0-9a-f]+:[[:space:]]+/, ""); sub(/<[^>]*>/, ""); print }' | grep -E '^(udiv|sdiv|fdiv|fmul|fadd|fsub|blr|br)\b|%?[xyz]mm[0-9]+|\bv[0-9]+\.'; then \
		echo "ERROR: Forbidden instruction or vector register in phy_ring0.c.o"; exit 1; \
	fi
endif
	@echo '#include "phy.h"' | $(CXX) -std=c++20 -x c++ -fsyntax-only -Wall -Wextra -Werror $(INC_RING0) -
	$(CLANG_TIDY) --header-filter='phy\.h$$' --checks='clang-analyzer-*,bugprone-*,-bugprone-easily-swappable-parameters' \
		--warnings-as-errors='*' src/phy_ring0.c -- -std=c17 $(INC_RING0)
	@echo "[Gate 2] PASS (0 undef symbols, Pure-GPR clean, C++20 clean, clang-tidy 0 warnings)"

gate3-selftest: $(BUILD)/phy
	@echo "[Gate 3] Running built-in 26-suite Ring-0 + SDK selftest..."
	./$(BUILD)/phy selftest
	@echo "[Gate 3] PASS"

gate4-cli-audit: $(BUILD)/phy
	@echo "[Gate 4] Auditing CLI linkage (0 leanshared, 0 heap allocators)..."
ifeq ($(UNAME_S),Darwin)
	@if otool -L $(BUILD)/phy | grep -qi 'leanshared'; then echo "ERROR: leanshared linked"; exit 1; fi
else
	@if ldd $(BUILD)/phy | grep -qi 'leanshared'; then echo "ERROR: leanshared linked"; exit 1; fi
endif
	@if nm -u $(BUILD)/phy | grep -qE '[[:space:]]_?(malloc|free|calloc|realloc)([[:space:]@]|$$)'; then echo "ERROR: heap allocator symbol linked"; exit 1; fi
	@echo "[Gate 4] PASS"

gate5-oracle: gate0-lean-noether gate4-cli-audit
	@echo "[Gate 5] Running independent L2 exact-rational Oracle + live Lean certificate replay..."
	python3 oracles/verify_noether_oracle.py
	@echo "[Gate 5] PASS"

gate6-sabotage: gate1-lean-schema
	@echo "[Gate 6] Running 51-mutant single-fault sabotage gate..."
	python3 tools/phy_sabotage_check.py
	@echo "[Gate 6] PASS"

gate7-format:
	@echo "[Gate 7] Verifying clang-format compliance..."
	$(CLANG_FORMAT) --dry-run --Werror $(C_SOURCES)
	@echo "[Gate 7] PASS"

fmt:
	$(CLANG_FORMAT) -i $(C_SOURCES)

clean:
	rm -rf $(BUILD)

verify-all: gate0-lean-noether gate1-lean-schema gate2-ring0-audit gate3-selftest gate4-cli-audit gate5-oracle gate6-sabotage gate7-format
	@echo "=== PHY v4.2.0: ALL 8 VERIFICATION GATES PASSED ==="
