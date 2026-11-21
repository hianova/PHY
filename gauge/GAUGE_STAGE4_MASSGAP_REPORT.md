# Stage 4 Report — 4D Pure $SU(2)$ Glueball $0^{++}$ Mass Gap & Direct Lattice Osterwalder–Schrader Positivity

> **Formal Companion:** [`FLOW/core/lean/Flow/ConstructiveQFT.lean`](../../FLOW/core/lean/Flow/ConstructiveQFT.lean) (`Theorems Q1–Q5`: OS reflection positivity, transfer Gram positivity, Hankel log-convexity, area law & variational GEVP bounds) | **Lattice:** `$12 \times 10^3$`, `6,000` sweeps (`3,000` meas = `36,000` source time-slices/$\beta$), 3D APE `smear_levels = [0, 2, 5, 10]` (temporal links untouched) | `20`-block jackknife

## 1. $0^{++}$ Scalar Glueball Mass Gap Spectrum & Continuum Ratio $m_{0^{++}} / \sqrt{\sigma}$

Static quark string tension $a\sqrt{\sigma}$ is extracted from $3\text{D}$-smeared (`level=5`) Wilson loops via Lüscher-subtracted potential $a^2 \sigma = V(3) - V(2) - \pi/72$ with $V(R) = \ln(W_{\text{sm}}(R,2)/W_{\text{sm}}(R,3))$. Published continuum $SU(2)$ benchmark: $m_{0^{++}}/\sqrt{\sigma} \approx 3.60 \pm 0.20$ (Teper 1987; Lucini & Teper 2001).

| $\beta$ | $\langle P \rangle$ (`<|L|>`) | $a\sqrt{\sigma}$ (`χ_sm(3,3)`) | $am_{\text{raw}}(0\to 1)$ | $am_{\text{gevp}}(0\to 1)$ | $am_{\text{gevp}}(1\to 2)$ ($am_{0^{++}}$) | $am_{\text{gevp}}(2\to 3)$ | $m_{0^{++}} / \sqrt{\sigma}$ |
| ---: | :--- | :--- | ---: | ---: | :--- | :--- | :---: |
| `2.30` | `0.60231 ± 0.00007` (`<|L|>=0.0131`) | `0.3791 ± 0.0011` (`χ33=-0.1014`) | `2.070 ± 0.038` | `1.394 ± 0.021` | **`1.231 ± 0.082`** | `0.870 ± 0.193` | **`3.25 ± 0.22`** |
| `2.40` | `0.63004 ± 0.00004` (`<|L|>=0.0142`) | `0.2787 ± 0.0008` (`χ33=-0.1761`) | `2.262 ± 0.047` | `1.178 ± 0.020` | **`0.956 ± 0.055`** | `0.975 ± 0.114` | **`3.43 ± 0.20`** |
| `2.50` | `0.65214 ± 0.00004` (`<|L|>=0.0180`) | `0.1984 ± 0.0009` (`χ33=-0.2068`) | `2.606 ± 0.078` | `1.039 ± 0.016` | **`0.709 ± 0.026`** | `0.592 ± 0.043` | **`3.57 ± 0.13`** |

## 2. Direct Verification of Osterwalder–Schrader Positivity (`Theorem Q1`) & Hankel Convexity (`Theorem Q2C`)

| $\beta$ | $\lambda_{\min}(C(0))$ | $C(1)$ Eigenvalues $\lambda_0 \ge \lambda_1 \ge \lambda_2 \ge \lambda_3$ | GEVP Correlator $C(0) > C(1) > C(2)$ | Hankel Ratio $C(0)C(2)/C(1)^2$ |
| ---: | ---: | :--- | :--- | ---: |
| `2.30` | `1.11e-07` (`> 0`) | `9.47e-06`, `5.57e-07`, `2.44e-08`, `6.34e-10` | `6.792e-06` → `1.685e-06` → `4.921e-07` | `1.1774` |
| `2.40` | `5.05e-08` (`> 0`) | `6.20e-06`, `4.43e-07`, `2.40e-08`, `1.16e-09` | `5.311e-07` → `1.634e-07` → `6.286e-08` | `1.2499` |
| `2.50` | `2.59e-08` (`> 0`) | `3.62e-06`, `3.09e-07`, `2.62e-08`, `5.96e-10` | `1.141e-07` → `4.036e-08` → `1.992e-08` | `1.3953` |

## 3. Structural Takeaways Linked to Lean 4 (`Flow.ConstructiveQFT`)

1. **Osterwalder–Schrader Time-Slice Positivity (`Theorem Q1`)**: Preserving temporal links untouched under 3D APE smearing yields positive $C(0)$ eigenvalues (`> 0`) and cleanly positive leading transfer-matrix eigenvalues in $C(1)$ across all couplings.
2. **Variational Rayleigh–Ritz Contraction (`Theorems Q2C & Q4`)**: Every multi-smearing GEVP strictly improves upon unsmeared raw operators ($am_{\text{gevp}}(0\to 1) < am_{\text{raw}}(0\to 1)$) and descends monotonically ($am_{\text{gevp}}(1\to 2) \le am_{\text{gevp}}(0\to 1)$) with convex Hankel decay $C(0)C(2)/C(1)^2 \gtrsim 1$.
3. **Asymptotic Continuum Mass Gap**: As lattice spacing $a(\beta)$ shrinks by `~2×` from $\beta = 2.30$ to $2.50$, both $a\sqrt{\sigma}$ and $am_{0^{++}}$ scale down together while leaving the physical dimensionless gap $m_{0^{++}}/\sqrt{\sigma}$ robustly positive near the Lucini–Teper continuum plateau.
