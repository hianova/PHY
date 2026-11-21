# Stage 2 Report — 2D $CP^{20}$ ($N = 21$) Nonlinear $\sigma$ Model Equal-Work Topological Unfreezing

> **Lattice:** `L = 32` (`V = 1024` sites) Periodic Torus | **Target:** $CP^{20}$ ($N = 21$) with Exact Auxiliary $U(1)$ (`Best–Fisher` von Mises HB + `2× OR` + `6× SU(2)` Givens) | **Equal Work Budget:** `30,000` single-lattice sweeps/$\beta$ (`R = 12`, `L_d = 8` defect slit) | **Formal Companion:** [`FLOW/core/lean/Flow/ProjectiveBundle.lean`](../../FLOW/core/lean/Flow/ProjectiveBundle.lean)

## 1. Continuum Observables, Finite-Volume Ratio & Asymptotic Freedom Scaling ($L = 32$)

Leading large-$N$ asymptotics predict $\chi_t \xi_G^2 \approx \frac{1}{2\pi N} = 0.00758$ for $N = 21$ (with small positive $1/N$ corrections). Literature cross-check against Hasenbusch (2017, Table II for $N=21$ periodic BC): $\beta = 0.65$ ($L=32$, exact match) gives $\xi_{2\text{nd}} = 2.6993(7)$, $\chi_t = 0.001170(9)$; $\beta = 0.70$ ($L=44$) gives $\xi_{2\text{nd}} = 3.744(2)$, $\chi_t = 0.000570(16)$.

| $\beta$ | Mean Link Energy $E$ (`Sym` vs `Loc`) | Corr. Length $\xi_G$ | $L / \xi_G$ | Cooled $\chi_t$ (`Symbiotic`) | Continuum $\chi_t \xi_G^2$ | Cooled $\chi_t$ (`local_or_metro`) |
| ---: | :--- | ---: | ---: | :--- | ---: | :--- |
| `0.60` | `0.56491 ± 0.00017` (`loc=0.56478`) | `1.941 ± 0.012` | `16.49` | `0.002352 ± 0.000154` | `0.0089` | `0.002372 ± 0.000140` |
| `0.65` | `0.60040 ± 0.00013` (`loc=0.60030`) | `2.715 ± 0.023` | `11.79` | `0.001396 ± 0.000205` | `0.0103` | `0.000921 ± 0.000102` |
| `0.70` | `0.63076 ± 0.00014` (`loc=0.63065`) | `3.997 ± 0.053` | `8.01` | `0.000452 ± 0.000041` | `0.0072` | `0.000641 ± 0.000249` |
| `0.73` | `0.64716 ± 0.00014` (`loc=0.64677`) | `5.193 ± 0.058` | `6.16` | `0.000231 ± 0.000040` | `0.0062` | `0.000603 ± 0.000140` |

## 2. Work-Normalized $\tau_{\text{int}}^{\text{work}}(Q_{\text{cool}})$ (Equal `30,000` Replica-Sweeps) & Fisher `2-Jet` Ladder

| $\beta$ | `local_or_metro` ($R=1$) | `bare_bp_inst` ($R=1$) | `symbiotic_ptbc_2jet` ($R=12$) | Equal-Work Speedup | Frozen Fisher `2-Jet` Ladder $c_0 \dots c_{11}$ |
| ---: | :--- | :--- | :--- | :---: | :--- |
| `0.60` | `53.7` (`1494` trans) | `52.3` (`acc=0.0000%`) | **`123.9`** (`1256` trans, `swap min=47.1%`, `mean=49.9%`) | **`0.43×`** | `1.00, 0.88, 0.77, 0.66, 0.57, 0.49, 0.41, 0.33, 0.25, 0.17, 0.09, 0.00` |
| `0.65` | `191.9` (`261` trans) | `274.6` (`acc=0.0000%`) | **`121.7`** (`912` trans, `swap min=41.1%`, `mean=47.4%`) | **`1.58×`** | `1.00, 0.87, 0.76, 0.65, 0.55, 0.46, 0.38, 0.31, 0.24, 0.16, 0.09, 0.00` |
| `0.70` | `1231.2` (`13` trans) | `2099.3` (`acc=0.0000%`) | **`135.6`** (`726` trans, `swap min=41.5%`, `mean=45.1%`) | **`9.08×`** | `1.00, 0.87, 0.74, 0.63, 0.53, 0.44, 0.36, 0.29, 0.22, 0.15, 0.08, 0.00` |
| `0.73` | `1301.5` (`5` trans) | `515.9` (`acc=0.0000%`) | **`261.6`** (`431` trans, `swap min=39.9%`, `mean=43.8%`) | **`4.97×`** | `1.00, 0.86, 0.73, 0.62, 0.51, 0.41, 0.34, 0.27, 0.21, 0.14, 0.07, 0.00` |

## 3. Structural Takeaways Linked to Lean 4 (`Flow.ProjectiveBundle`)

1. **Hopf $U(1)$-Bundle & Exact Bogomol'nyi Bound (`Theorems P1–P2`)**: Auxiliary $U(1)$ heatbath + $SU(2) \subset SU(21)$ Givens updates preserve exact local gauge invariance and unit norm $\|z_x\| = 1$ for all $N$, while Berg–Lüscher geometric charge $Q_{\text{geom}}$ bounds lattice action strictly from below.
2. **Corollary 2B on Non-Abelian Moduli (`bare_bp_inst`)**: Direct Belavin–Polyakov $SU(2)$ instanton proposals face ultraviolet component mismatch across thermal backgrounds, collapsing instanton acceptance below `1%` at $\beta \ge 0.65$.
3. **Defect-Slit Bianchi Opening & Fisher `2-Jet` Ladder (`Theorems P3–P4`)**: Opening an $L_d = 8$ defect slit (`Theorem P3A`) with exact product-weight `@swap` detailed balance (`Theorem P3B`) and a Cauchy–Schwarz-optimal uniform Fisher arc ladder (`Theorem P4`) maintains balanced swap rates (`> 15%` min pair acceptance across all `11` links) while keeping $L/\xi_G \ge 5.0$ across $\beta \in \{0.60 \dots 0.73\}$ and breaking $CP^{20}$ topological freezing at equal compute cost.
