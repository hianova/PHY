# Stage 5 Report — Certified Lattice Yang–Mills Bootstrap & Exact Renormalization Group

> **Formal Companion:** [`FLOW/core/lean/Flow/LatticeBootstrap.lean`](../../FLOW/core/lean/Flow/LatticeBootstrap.lean) (`Theorems B1–B5`: exact SD moment reductions, certified dual SOS rational certificates, GWW $\lambda_c = 2$ necessity, 4D quaternionic SD/Gram bounds, block-spin & Migdal–Kadanoff RG)

## 1. Certified Exact-Rational Dual SOS Bootstrap (`Theorems B1 & B5`)

Exact rational Schwinger–Dyson moment reductions (`Fraction` arithmetic, zero floats) produce dual SOS integer witness vectors enclosing analytical Bessel benchmarks strictly inside certified intervals (with $\beta=2.0, 2.5$ verified in Lean 4 directly by `Theorem B5`).

| Model | $\beta$ | $K$ | Analytical Exact $u_P$ | Certified Rational Bounds `[lo, hi]` | Width | Lean 4 / Oracle Status |
| :--- | ---: | ---: | ---: | :--- | ---: | :--- |
| 2D $U(1)$ | `1.0` | `6` | `0.4463899659` | `[0.4463899626, 0.4463899692]` | `6.51e-09` | PASS |
| 2D $U(1)$ | `2.0` | `6` | `0.6977746580` | `[0.6977741898, 0.6977751261]` | `9.36e-07` | PASS (Lean B5) |
| 2D $U(1)$ | `4.0` | `8` | `0.8635226110` | `[0.8635226037, 0.8635226184]` | `1.47e-08` | PASS |
| 2D $U(1)$ | `6.0` | `8` | `0.9123593044` | `[0.9123592782, 0.9123593303]` | `5.21e-08` | PASS |
| 2D $SU(2)$ | `1.0` | `5` | `0.2401937239` | `[0.2401936990, 0.2401937402]` | `4.12e-08` | PASS |
| 2D $SU(2)$ | `2.0` | `5` | `0.4331274267` | `[0.4331258964, 0.4331283882]` | `2.49e-06` | PASS (Lean B5) |
| 2D $SU(2)$ | `2.5` | `5` | `0.5071951000` | `[0.5071906916, 0.5071977647]` | `7.07e-06` | PASS (Lean B5) |
| 2D $SU(2)$ | `4.0` | `6` | `0.6580472674` | `[0.6580466831, 0.6580476261]` | `9.43e-07` | PASS |

## 2. Gross–Witten–Wadia ($SU(\infty)$) Nonlinear Loop Bootstrap Across $\lambda_c = 2$ (`Theorem B1`)

| $\lambda$ | Exact GWW $u_P$ | Certified Toeplitz Interval `[lo, hi]` | Width | Alternating Vector $v_i = (-1)^i$ on Strong Branch |
| ---: | ---: | :--- | ---: | :--- |
| `0.50` | `0.875000` | `[0.874500, 0.875000]` | `5.00e-04` | `-8.000` (`K=4`) |
| `1.00` | `0.750000` | `[0.745500, 0.750000]` | `4.50e-03` | `-2.000` (`K=4`) |
| `1.50` | `0.625000` | `[0.612250, 0.625000]` | `1.28e-02` | `-0.667` (`K=6`) |
| `2.00` | `0.500000` | `[0.490000, 0.502250]` | `1.22e-02` | `>= 0 (strong valid)` |
| `2.50` | `0.400000` | `[0.398000, 0.401750]` | `3.75e-03` | `>= 0 (strong valid)` |
| `3.00` | `0.333333` | `[0.333000, 0.333750]` | `7.50e-04` | `>= 0 (strong valid)` |
| `4.00` | `0.250000` | `[0.250000, 0.250000]` | `0.00e+00` | `>= 0 (strong valid)` |

## 3. 4D Pure $SU(2)$ Exact Schwinger–Dyson, OS & 7-Staple Gram Audit (`8^4`, `1000` meas, `20`-block JK, `Theorems B2–B3`)

Algebraic checks pass to machine precision across all couplings (`|gram_diff| < 1e-12`, `|R_sq - (6 + 6*w_1x2 + 24*w_chair)| < 1e-11`). Analytical non-perturbative bounds come from 2D decoupling $I_2(\beta)/I_1(\beta) < u_P < \min\big((-1 + \sqrt{1 + 16\beta^2})/(4\beta),\ I_2(6\beta)/I_1(6\beta)\big)$.

| $\beta$ | $\langle u_P \rangle$ | Analytic Bracket | $W_{1\times 2}$ (`u_copl2`) | $W_{\text{chair}}$ (`u_orth2`) | 4D SD Residual (`3u_P - beta*...`) | 7-Staple `Gram7` |
| ---: | :--- | :--- | :--- | :--- | :--- | ---: |
| `1.50` | `0.36248 ± 0.00012` | `[0.3441, 0.8386]` | `0.1324` (`0.1316`) | `0.1446` (`0.1338`) | `-0.00029 ± 0.00034` (`0.87σ`) | `0.9224` |
| `2.00` | `0.50125 ± 0.00018` | `[0.4331, 0.8778]` | `0.2592` (`0.2520`) | `0.2889` (`0.2559`) | `-0.00027 ± 0.00048` (`0.57σ`) | `0.9071` |
| `2.30` | `0.60229 ± 0.00016` | `[0.4793, 0.8934]` | `0.3859` (`0.3638`) | `0.4216` (`0.3669`) | `+0.00015 ± 0.00036` (`0.42σ`) | `0.8960` |
| `2.40` | `0.63059 ± 0.00017` | `[0.4935, 0.8978]` | `0.4260` (`0.3987`) | `0.4610` (`0.4013`) | `+0.00095 ± 0.00038` (`2.52σ`) | `0.8843` |
| `2.50` | `0.65257 ± 0.00010` | `[0.5072, 0.9018]` | `0.4576` (`0.4268`) | `0.4916` (`0.4291`) | `+0.00059 ± 0.00048` (`1.23σ`) | `0.8688` |

## 4. 4D ($b = 2, D = 4$) Migdal–Kadanoff Character RG Flow (`Theorem B4`)

4D Abelian $U(1)$ transitions from confining ($\Delta\beta < 0$ at $\beta=2.0$) to a deconfined Coulomb fixed ray ($\Delta\beta > 0$ for $\beta \ge 4.0$), whereas non-Abelian 4D $SU(2)$ exhibits universal asymptotic freedom ($\Delta\beta < 0$ at every tested $\beta \in [2.0, 10.0]$) flowing monotonically into strong-coupling confinement.

| Initial $\beta$ | 4D $U(1)$ $c_1(0 \to 1 \to 4 \to 12)$ | $U(1)$ $\Delta\beta_1$ | 4D $SU(2)$ $c_1(0 \to 1 \to 4 \to 12)$ | $SU(2)$ $\Delta\beta_1$ |
| ---: | :--- | ---: | :--- | ---: |
| `2.0` | `0.6978 -> 0.6649 -> 0.1594 -> 0.0000` | `-0.1849` | `0.4331 -> 0.1384 -> 0.0000 -> 0.0000` | `-1.4390` |
| `4.0` | `0.8635 -> 0.8680 -> 0.8683 -> 0.8682` | `+0.1206` | `0.6580 -> 0.5472 -> 0.0000 -> 0.0000` | `-1.1825` |
| `6.0` | `0.9124 -> 0.9140 -> 0.9141 -> 0.9141` | `+0.1073` | `0.7627 -> 0.7195 -> 0.2815 -> 0.0000` | `-0.9974` |
| `8.0` | `0.9352 -> 0.9361 -> 0.9361 -> 0.9361` | `+0.1019` | `0.8192 -> 0.7963 -> 0.6410 -> 0.0000` | `-0.9418` |
| `10.0` | `0.9486 -> 0.9491 -> 0.9492 -> 0.9492` | `+0.0996` | `0.8542 -> 0.8400 -> 0.7602 -> 0.0000` | `-0.9164` |
