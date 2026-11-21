UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
LLVM_CLANG := $(wildcard /opt/homebrew/opt/llvm/bin/clang)
CC ?= $(if $(LLVM_CLANG),$(LLVM_CLANG),clang)
else
CC ?= clang
endif
PYTHON ?= python3

CFLAGS := -std=c17 -O3 -Wall -Wextra -Werror -Wpedantic -Wconversion
BUILD := build

# Stage 1: 2D Compact U(1)
BIN_U1     := $(BUILD)/u1_2d_windtunnel
BENCH_U1   := $(BUILD)/u1_2d_bench.json
REPORT_U1  := GAUGE_STAGE1_REPORT.md

# Stage 2: 2D CP^20 (N=21) Equal-Work Wind Tunnel
BIN_CPN    := $(BUILD)/cpn_2d_windtunnel
BENCH_CPN  := $(BUILD)/cpn_2d_bench.json
REPORT_CPN := GAUGE_STAGE2_CPN_REPORT.md

# Stage 3: 4D Pure SU(2) Yang-Mills Wind Tunnel
BIN_SU2    := $(BUILD)/su2_4d_windtunnel
BENCH_SU2  := $(BUILD)/su2_4d_bench.json
REPORT_SU2 := GAUGE_STAGE3_SU2_REPORT.md

# Stage 4: 4D Pure SU(2) Glueball 0++ Mass Gap & OS Positivity
BIN_MG     := $(BUILD)/su2_4d_massgap
BENCH_MG   := $(BUILD)/su2_4d_massgap.json
REPORT_MG  := GAUGE_STAGE4_MASSGAP_REPORT.md

# Stage 5: Certified Lattice Bootstrap & Exact Renormalization Group
BIN_BS     := $(BUILD)/su2_4d_bootstrap_sd
BENCH_BS   := $(BUILD)/su2_4d_bootstrap_sd.json
REPORT_BS  := GAUGE_STAGE5_BOOTSTRAP_RG_REPORT.md

# Stage 6: 4D Pure SU(3) Color Yang-Mills (Cabibbo-Marinari, Z_3 First-Order & Casimir Scaling)
BIN_SU3    := $(BUILD)/su3_4d_windtunnel
BENCH_SU3  := $(BUILD)/su3_4d_bench.json
REPORT_SU3 := GAUGE_STAGE6_SU3_REPORT.md

# Stage 7: 2D/4D Lattice Fermions, Chiral Symmetry, Atiyah-Singer Index & Sign Problem
BIN_FERM    := $(BUILD)/fermion_2d_windtunnel
BENCH_FERM  := $(BUILD)/fermion_2d_bench.json
REPORT_FERM := GAUGE_STAGE7_FERMION_REPORT.md

# Stage 8: Discrete Quantum Gravity & Holographic AME Duality
BIN_GRAV    := $(BUILD)/gravity_holo_windtunnel
BENCH_GRAV  := $(BUILD)/gravity_holo_bench.json
REPORT_GRAV := GAUGE_STAGE8_GRAVITY_HOLO_REPORT.md

# Stage 9: Complete Physics — SM Unification, 4D SU(2)-Higgs & Regge-Maxwell
BIN_UNIF    := $(BUILD)/unified_sm_windtunnel
BENCH_UNIF  := $(BUILD)/unified_sm_bench.json
REPORT_UNIF := GAUGE_STAGE9_UNIFIED_REPORT.md

# Stage 10: First-Principles Certified Quantum Chemistry (STO-3G FCI, 2-RDM DQG SDP & Conformer PES)
BIN_CHEM    := $(BUILD)/chem_abinitio_windtunnel
BENCH_CHEM  := $(BUILD)/chem_abinitio_bench.json
REPORT_CHEM := GAUGE_STAGE10_CHEM_REPORT.md

.PHONY: all stage1 stage2 stage3 stage4 stage5 stage6 stage7 stage8 stage9 stage10 compile run verify compile-cpn run-cpn verify-cpn compile-su2 run-su2 verify-su2 compile-mg run-mg verify-mg compile-bs run-bs verify-bs compile-su3 run-su3 verify-su3 compile-ferm run-ferm verify-ferm compile-grav run-grav verify-grav compile-unif run-unif verify-unif compile-chem run-chem verify-chem clean

all: stage1 stage2 stage3 stage4 stage5 stage6 stage7 stage8 stage9 stage10

stage1: compile run verify

stage2: compile-cpn run-cpn verify-cpn

stage3: compile-su2 run-su2 verify-su2

stage4: compile-mg run-mg verify-mg

stage5: compile-bs run-bs verify-bs

stage6: compile-su3 run-su3 verify-su3

stage7: compile-ferm run-ferm verify-ferm

stage8: compile-grav run-grav verify-grav

stage9: compile-unif run-unif verify-unif

stage10: compile-chem run-chem verify-chem

$(BUILD):
	@mkdir -p $(BUILD)

$(BIN_U1): u1_2d_windtunnel.c | $(BUILD)
	$(CC) $(CFLAGS) $< -lm -o $@

compile: $(BIN_U1)

run: $(BIN_U1)
	./$(BIN_U1) --json $(BENCH_U1)

verify: $(BENCH_U1) u1_2d_oracle.py
	$(PYTHON) u1_2d_oracle.py --bench $(BENCH_U1) --report $(REPORT_U1)

$(BIN_CPN): cpn_2d_windtunnel.c | $(BUILD)
	$(CC) $(CFLAGS) $< -lm -lpthread -o $@

compile-cpn: $(BIN_CPN)

run-cpn: $(BIN_CPN)
	./$(BIN_CPN) --json $(BENCH_CPN)

verify-cpn: $(BENCH_CPN) cpn_2d_oracle.py
	$(PYTHON) cpn_2d_oracle.py --bench $(BENCH_CPN) --report $(REPORT_CPN)

$(BIN_SU2): su2_4d_windtunnel.c | $(BUILD)
	$(CC) $(CFLAGS) $< -lm -lpthread -o $@

compile-su2: $(BIN_SU2)

run-su2: $(BIN_SU2)
	./$(BIN_SU2) --json $(BENCH_SU2)

verify-su2: $(BENCH_SU2) su2_4d_oracle.py
	$(PYTHON) su2_4d_oracle.py --bench $(BENCH_SU2) --report $(REPORT_SU2)

$(BIN_MG): su2_4d_massgap.c | $(BUILD)
	$(CC) $(CFLAGS) $< -lm -lpthread -o $@

compile-mg: $(BIN_MG)

run-mg: $(BIN_MG)
	./$(BIN_MG) --json $(BENCH_MG)

verify-mg: $(BENCH_MG) su2_4d_massgap_oracle.py
	$(PYTHON) su2_4d_massgap_oracle.py --bench $(BENCH_MG) --report $(REPORT_MG)

$(BIN_BS): su2_4d_bootstrap_sd.c | $(BUILD)
	$(CC) $(CFLAGS) $< -lm -lpthread -o $@

compile-bs: $(BIN_BS)

run-bs: $(BIN_BS)
	./$(BIN_BS) --json $(BENCH_BS)

verify-bs: $(BENCH_BS) bootstrap_rg_oracle.py
	$(PYTHON) bootstrap_rg_oracle.py --bench $(BENCH_BS) --report $(REPORT_BS)

$(BIN_SU3): su3_4d_windtunnel.c | $(BUILD)
	$(CC) $(CFLAGS) $< -lm -lpthread -o $@

compile-su3: $(BIN_SU3)

run-su3: $(BIN_SU3)
	./$(BIN_SU3) --json $(BENCH_SU3)

verify-su3: $(BENCH_SU3) su3_4d_oracle.py
	$(PYTHON) su3_4d_oracle.py --bench $(BENCH_SU3) --report $(REPORT_SU3)

$(BIN_FERM): fermion_2d_windtunnel.c | $(BUILD)
	$(CC) $(CFLAGS) $< -lm -lpthread -o $@

compile-ferm: $(BIN_FERM)

run-ferm: $(BIN_FERM)
	./$(BIN_FERM) --json $(BENCH_FERM)

verify-ferm: $(BENCH_FERM) fermion_2d_oracle.py
	$(PYTHON) fermion_2d_oracle.py --bench $(BENCH_FERM) --report $(REPORT_FERM)

$(BIN_GRAV): gravity_holo_windtunnel.c | $(BUILD)
	$(CC) $(CFLAGS) $< -lm -lpthread -o $@

compile-grav: $(BIN_GRAV)

run-grav: $(BIN_GRAV)
	./$(BIN_GRAV) --json $(BENCH_GRAV)

verify-grav: $(BENCH_GRAV) gravity_holo_oracle.py
	$(PYTHON) gravity_holo_oracle.py --bench $(BENCH_GRAV) --report $(REPORT_GRAV)

$(BIN_UNIF): unified_sm_windtunnel.c | $(BUILD)
	$(CC) $(CFLAGS) $< -lm -lpthread -o $@

compile-unif: $(BIN_UNIF)

run-unif: $(BIN_UNIF)
	./$(BIN_UNIF) --json $(BENCH_UNIF)

verify-unif: $(BENCH_UNIF) unified_sm_oracle.py
	$(PYTHON) unified_sm_oracle.py --bench $(BENCH_UNIF) --report $(REPORT_UNIF)

$(BIN_CHEM): chem_abinitio_windtunnel.c | $(BUILD)
	$(CC) $(CFLAGS) $< -lm -o $@

compile-chem: $(BIN_CHEM)

run-chem: $(BIN_CHEM)
	./$(BIN_CHEM) --json $(BENCH_CHEM)

verify-chem: $(BENCH_CHEM) chem_abinitio_oracle.py
	$(PYTHON) chem_abinitio_oracle.py --bench $(BENCH_CHEM) --report $(REPORT_CHEM)

clean:
	rm -rf $(BUILD)
