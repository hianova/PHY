# Stage 8 Report: Discrete Quantum Gravity & Holographic AME Duality (`Theorems H1–H5`)

Generated from [`build/gravity_holo_bench.json`](build/gravity_holo_bench.json) via [`gravity_holo_windtunnel.c`](gravity_holo_windtunnel.c) and audited against [`gravity_holo_oracle.py`](gravity_holo_oracle.py).

## 1. Startup Algebraic, Topological & Holographic Self-Tests (`Flow.GravityBundle`)

| Check | Theorem Link | Measured Residual / Value |
| :--- | :--- | :--- |
| `PHY v4.0.0` Nonlinear `AME(4,8)` 3-Bipartition Unitarity & Purity | `Theorem H2` | 3/3 bijections exact (`max‖U†U-I‖ = 0.00e+00`), `Tr(ρ²) = 0.015625 = 1/64` |
| HaPPY `[[5,1,3]]` `AME(6,2)` 20 `3v3` Unitaries & 10-Pair 2-Erasure KL | `Theorems H2–H3` | `max‖U†U-I_8‖ = 4.44e-16`, 160 Pauli checks `max‖V†EV-cI‖ = 0.00e+00` |
| 2D Gauss–Bonnet ($S^2, T^2, \Sigma_2$) + Pachner `2-2` & `1-3` | `Theorems H1A–H1D` | `max|Σϵ/(2π) - χ| = 8.88e-16` |
| 3D Schläfli Identity & $\partial\Delta^4$ Action Gradient vs. Central FD | `Theorem H1E` | `Schläfli = 2.13e-10`, `ActionGrad = 1.34e-08` |
| `{5,4}` `R=2` (`V=51, L=95`) Holographic Entropy Cone (100 Triples) | `Theorem H4` | `SSA violations = 0`, `MMI (I_3 <= 0) violations = 0` |
| Holographic No-Cloning Complementary Erasure Mutual Information | `Theorem H3C` | `I(R : Abar) = 0.00e+00` |

## 2. Experiment 1A — 2D Dynamical Regge Gravity Monte Carlo (`λ = 0.5`)

| Topology (`χ`) | $\alpha$ | $\langle A_{\text{tot}} \rangle$ | $\langle R^2 \rangle = \langle \sum \epsilon_v^2 / A_v \rangle$ | Cauchy–Schwarz Ratio | `max|Σϵ/(2π)-χ|` | `acc` |
| :--- | ---: | :---: | :---: | :---: | :---: | ---: |
| $S^2$ (`χ=+2`) | `0.0` | `29.6717 ± 1.0589` | **`24.4032 ± 0.7949`** | `4.58533` | `7.33e-15` | `73.5%` |
| $S^2$ (`χ=+2`) | `0.5` | `39.2477 ± 0.5794` | **`9.3200 ± 0.1273`** | `2.31639` | `9.99e-15` | `74.0%` |
| $S^2$ (`χ=+2`) | `2.0` | `48.5174 ± 0.5569` | **`5.0844 ± 0.0438`** | `1.56212` | `4.22e-15` | `69.2%` |
| $S^2$ (`χ=+2`) | `8.0` | `67.3781 ± 0.4149` | **`2.8554 ± 0.0161`** | `1.21832` | `4.44e-15` | `59.2%` |
| $T^2$ (`χ=0`) | `0.0` | `47.0818 ± 1.0238` | **`25.6547 ± 0.7803`** | `—` | `1.63e-14` | `73.5%` |
| $T^2$ (`χ=0`) | `0.5` | `54.0021 ± 0.7282` | **`7.8104 ± 0.1147`** | `—` | `9.19e-15` | `72.5%` |
| $T^2$ (`χ=0`) | `2.0` | `58.6632 ± 1.0208` | **`3.0010 ± 0.0306`** | `—` | `9.47e-15` | `63.8%` |
| $T^2$ (`χ=0`) | `8.0` | `59.0645 ± 1.3869` | **`0.8669 ± 0.0064`** | `—` | `1.07e-14` | `42.0%` |

## 3. Experiment 1B — 3D Regge De Sitter Saddle, Indefinite Hessian & `2-Jet` Newton (`κ=1, Λ=1`)

| 3-Sphere Triangulation | Exact De Sitter Saddle | Hessian Signature (`GHP`) | Gradient Descent & Ascent | `2-Jet` Full-Hessian Newton |
| :--- | :--- | :--- | :--- | :--- |
| `4-simplex dDelta4` (`V=5,E=10,T=5`) | `l*=3.827922` (`S*=-66.1033`) | `(+1, -9)` (`[-2.603, 1.353]`) | GD/GA diverge | **`4` iters** (`|g|=3.35e-15`, `dl=8.88e-16`) |
| `16-cell dBeta4` (`V=8,E=24,T=16`) | `l*=2.401504` (`S*=-52.2317`) | `(+3, -21)` (`[-3.268, 1.132]`) | GD/GA diverge | **`4` iters** (`|g|=5.03e-15`, `dl=1.78e-15`) |

## 4. Experiment 2 — `PHY AME(4,8)` 4-Tensor Bulk Loop (`16,777,216`-Dim State, 256 Cuts & 65,536 Triples)

- **Exact Ryu–Takayanagi Match:** **`256/256`** boundary subsets satisfy $S_1(A) = 3\,\text{minCut}(A)$ (`max_err = 0.00e+00`).
- **Flat Fixed-Area Entanglement Spectrum:** `max|S_2(A) - S_1(A)| = 0.00e+00` across all 256 bipartitions.
- **Ryu–Takayanagi Mutual Information Transition:** Adjacent tensor pair `I({0,1}:{2,3}) = 6.0000` bits vs. opposite pair `I({0,1}:{4,5}) = 0.0000` bits.
- **Exhaustive Entropy Cone Verification:** **`0` `SSA` violations** and **`0` `MMI` ($I_3 \le 0$) violations** across all **`65,536` ($4^8$)** disjoint boundary partitions `(A, B, C)`.

## 5. Experiment 3A — 256-Dim 2-Pentagon HaPPY `[[5,1,3]]` Code & Quantum Extremal Surface (`Theorem H3 & H5`)

Bulk-to-boundary isometry precision: `||V†V - I_4|| = 2.22e-16` (`numpy` independent oracle: `2.22e-16`). Fixed pure stabilizer `|00>` shows exact RT on `9/9` prefix & `56/58` contiguous intervals (`244/256` total, with `2` symmetric `2+2` neck-straddling stabilizer parity locks `S=3 < 4` per Pastawski–Yoshida–Harlow–Preskill); adding bulk horizon entropy (`1-BH`, `2-BH`) restores **`58/58`** and **`256/256`** exact QES saturation (`<= 8.88e-16`).

| Bulk State (`S_BH`) | Connected Arcs (`/58`) | All Subsets (`/256`) | Prefix Entropy Curve `S(k=0..8)` | `numpy` Cross-Check Err |
| :--- | :---: | :---: | :--- | :---: |
| `Pure AdS |00>` | **`56/58`** | **`244/256`** (`max_err=1.00e+00`) | `[0, 1, 2, 2, 1, 2, 2, 1, 0]` | `1.92e-15` |
| `1-BH (S_BH=1)` | **`58/58`** | **`256/256`** (`max_err=8.88e-16`) | `[0, 1, 2, 3, 2, 3, 3, 2, 1]` | `8.88e-16` |
| `2-BH (S_BH=2)` | **`58/58`** | **`256/256`** (`max_err=8.88e-16`) | `[0, 1, 2, 3, 2, 3, 4, 3, 2]` | `8.88e-16` |

## 6. Experiment 3B — `{{5,4}}` Hyperbolic $\mathrm{{AdS}}_3$ Tiling & Black Hole Page Curve (`Theorem H5`)

- **Hyperbolic Exponential Growth (`R = 0..4`):** `R=0: (V=1, L=5)` -> `R=1: (V=11, L=25)` -> `R=2: (V=51, L=95)` -> `R=3: (V=201, L=355)` -> `R=4: (V=761, L=1325)` (`V/L = 0.57434 -> 1/sqrt(3) = 0.57735`).
- **`R=3` (`V=201, L=355`) CFT Logarithmic Area Law:** $S_{RT}(\ell) \approx 1.479\log_2(\ell) + 0.802$ (`R^2 = 0.99394`, effective `c_log2 = 3*slope = 4.437` across $\ell \in \{1..177\}$).

| Horizon Sink `H` | Endpoint & Turnover (`S(0) -> peak -> S(L)`) | Page Curve `S(k)` at `k/L ∈ {0.0 .. 1.0}` |
| :--- | :--- | :--- |
| `H=0` | `0.00` -> peak **`11.86`** (`k/L=0.50`) -> **`0.00`** | `[0.00, 8.66, 10.59, 11.14, 11.65, 11.86, 11.65, 11.14, 10.59, 8.66, 0.00]` |
| `H=5` | `0.00` -> peak **`13.65`** (`k/L=0.80`) -> **`5.00`** | `[0.00, 8.66, 10.59, 11.14, 11.65, 12.15, 12.65, 13.15, 13.65, 13.49, 5.00]` |
| `H=25` | `0.00` -> peak **`28.73`** (`k/L=0.90`) -> **`25.00`** | `[0.00, 8.66, 11.20, 13.66, 16.20, 18.73, 21.20, 23.73, 26.20, 28.73, 25.00]` |
| `H=95` | `0.00` -> peak **`95.00`** (`k/L=1.00`) -> **`95.00`** | `[0.00, 12.25, 21.89, 31.25, 40.89, 50.52, 59.89, 69.52, 78.89, 88.52, 95.00]` |
