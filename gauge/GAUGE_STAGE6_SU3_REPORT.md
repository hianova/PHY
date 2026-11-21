# Stage 6 Report: 4D Pure $SU(3)$ Color Yang–Mills Wind Tunnel (`Theorems C1–C5`)

Generated from [`build/su3_4d_bench.json`](build/su3_4d_bench.json) via [`su3_4d_windtunnel.c`](su3_4d_windtunnel.c) and audited against [`su3_4d_oracle.py`](su3_4d_oracle.py).

## 1. Startup Mathematical & Algebraic Self-Tests (`Flow.ColorGaugeBundle`)

| Check | Theorem Link | Measured Residual / Comparison |
| :--- | :--- | :--- |
| Unitarity & Determinant | `Theorems C2A–C2E` | `max‖U†U - I‖ = 8.88e-16 < 1e-13`, `max|det U - 1| = 8.90e-16 < 1e-13` |
| 4D Local $SU(3)$ Gauge Invariance | `Theorems C1A–C1C` | `max|ΔP| = 6.07e-18 < 1e-13` |
| $\mathbb{Z}_3$ Slice Twist (`t0=1, ω`) | `Theorems C3A–C4A` | `|ΔS| = 0.00e+00`, `|L' - ωL| = 1.56e-17`, `|ΔL_8| = 8.46e-18` |
| Cabibbo–Marinari OR Involution | `Theorems C2C–C2D` | `|ΔS_sweep| = 2.00e-11`, `‖OR^2(U)-U‖ = 7.02e-16` |
| Single-Link CM Heatbath vs Weyl Haar (`200k` each) | `Theorem C2C` | `β=3.0`: `0.202490±0.000590` vs exact `0.203075` (`0.99σ`); `β=6.0`: `0.423074±0.000598` vs exact `0.422532` (`0.91σ`) |
| 36-Plaquette 4D Defect-Wall `@swap` | `Theorem C5A` | `|ΔS_wall - ΔS_full| = 1.50e-12 < 1e-10` |
| 3D APE & 4D Cooling + Embedded BPST $Q$ | `Theorems C1C–C2E` | `ΔP_3D = 1.81e-12`, `ΔQ_cool = 7.62e-13`, BPST `Q = +1 [+0.723] / -1 [-0.723]` |

## 2. Experiment 1 — $4 \times 8^3$ $\mathbb{Z}_3$ First-Order Deconfinement & Symbiotic Unfreezing

| $\beta$ | `<P>` (`Sym` vs `Loc`) | `<\|L\|>`, `<cos(3 arg L)>` & `<L_8>` | `local_cm_hb_or` (`|<L>|`, `hops`, `sec`) | `symbiotic_z3_center` (`|<L>|`, `τ_ReL`, `hops`, `sec`) |
| ---: | :--- | :--- | :--- | :--- |
| `5.55` | `0.51089 ± 0.00039` (`loc=0.51077`) | `0.0229 ± 0.0007` (`cos3=+0.097`, `<L_8>=+0.0003`) | `0.0019` (`hops=317`, `sec=0.34/0.38/0.28`) | **`0.0004`** (`τ=0.50`, **`514`** hops, `sec=0.33/0.35/0.33`) |
| `5.65` | `0.53911 ± 0.00100` (`loc=0.53779`) | `0.0664 ± 0.0093` (`cos3=+0.494`, `<L_8>=+0.0034`) | `0.0031` (`hops=150`, `sec=0.33/0.38/0.29`) | **`0.0032`** (`τ=0.50`, **`487`** hops, `sec=0.35/0.32/0.33`) |
| `5.69` | `0.55161 ± 0.00071` (`loc=0.54876`) | `0.1282 ± 0.0128` (`cos3=+0.841`, `<L_8>=+0.0124`) | `0.0294` (`hops=73`, `sec=0.13/0.41/0.46`) | **`0.0075`** (`τ=0.54`, **`500`** hops, `sec=0.31/0.37/0.32`) |
| `5.73` | `0.55993 ± 0.00036` (`loc=0.56075`) | `0.1597 ± 0.0064` (`cos3=+0.956`, `<L_8>=+0.0178`) | `0.1770` (`hops=0`, `sec=1.00/0.00/0.00`) | **`0.0031`** (`τ=0.55`, **`481`** hops, `sec=0.33/0.35/0.33`) |
| `6.15` | `0.60977 ± 0.00015` (`loc=0.60977`) | `0.2785 ± 0.0030` (`cos3=+0.996`, `<L_8>=+0.0571`) | `0.2790` (`hops=0`, `sec=1.00/0.00/0.00`) | **`0.0061`** (`τ=0.50`, **`497`** hops, `sec=0.32/0.33/0.35`) |

**First-Order Metastability / Hysteresis at $\beta = 5.69$:** Cold start gives `<P> = 0.55210 ± 0.00065`, `<|L|> = 0.1316 ± 0.0057` vs. hot start `<P> = 0.54549 ± 0.00106`, `<|L|> = 0.0565 ± 0.0146`; in the `1,200`-sweep coexistence run, both confined (`|L|<0.09`: `0.417`) and deconfined (`|L|>0.14`: `0.303`) states are visited repeatedly.

## 3. Experiment 2 — $8^4$ Confined Scaling & Fundamental vs. Adjoint Casimir Ratio (`Theorem C4C`)

| $\beta$ | `<P>` (`<\|L\|>`) | Fundamental $\chi_3(2,2)$ / $\chi_3(3,3)$ | Adjoint $\chi_8(2,2)$ | Correlated Casimir Ratio $\chi_8(2,2)/\chi_3(2,2)$ |
| ---: | :--- | :---: | :---: | :---: |
| `5.70` | `0.54903 ± 0.00028` (`0.0148`) | `0.3048 ± 0.0018` / `0.2202 ± 0.0039` | `0.6590 ± 0.0055` | **`2.162 ± 0.008`** (vs `9/4=2.25`) |
| `5.80` | `0.56728 ± 0.00024` (`0.0158`) | `0.2527 ± 0.0013` / `0.1687 ± 0.0027` | `0.5570 ± 0.0033` | **`2.205 ± 0.005`** (vs `9/4=2.25`) |
| `6.00` | `0.59426 ± 0.00017` (`0.0278`) | `0.1923 ± 0.0007` / `0.1038 ± 0.0014` | `0.4311 ± 0.0017` | **`2.242 ± 0.005`** (vs `9/4=2.25`) |

## 4. Experiment 3 — $8^4$ Topological Susceptibility & 4D Defect-Wall PTBC (`β = 5.90`, `eq_work = 800` sweeps)

- **Local `local_cm_hb_or` (`80` meas × `10` sweeps):** `<P> = 0.58215 ± 0.00027`, `<|L|> = 0.0221`, `<Q^2> = 0.525`, `χ_t = 1.282e-04 ± 2.388e-05`, `τ_int^work(Q) = 16.2`, `tunnelings = 34`.
- **Symbiotic `symbiotic_ptbc_4d` (`R = 8`, `50` meas × `2` sweeps, `L_d = 2` [`36` defect plaq], Fisher `2-Jet` ladder [`1.000`, `0.875`, `0.765`, `0.666`, `0.572`, `0.480`, `0.381`, `0.000`]):** `<P> = 0.58305 ± 0.00066`, `<|L|> = 0.0165`, `<Q^2> = 0.140`, `χ_t = 3.418e-05 ± 1.479e-05`, `τ_int^work(Q) = 15.0`, `tunnelings = 11`, `sw_min = 38.3%`, `sw_mean = 51.7%`.
- **Unbiased Agreement:** `Δ<P> = 1.27σ`, `Δχ_t = 3.35σ`.
