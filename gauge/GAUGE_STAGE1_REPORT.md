# Stage 1 Report — 2D Compact $U(1)$ Lattice Gauge Theory & Symbiotic `FLOW` Breaking Topological Freezing

> **Lattice:** `L = 16` (`V = 256` sites) Periodic Torus | **Sweeps:** `2,000` therm + `50,000` meas (step size frozen post-therm for exact detailed balance) | **Formal Companion:** [`FLOW/core/lean/Flow/GaugeBundle.lean`](../../FLOW/core/lean/Flow/GaugeBundle.lean)

## 1. Unbiased Equilibrium Verification vs Exact Finite-Volume Bessel Character Expansion

| $\beta$ | Exact $\langle P \rangle$ | Symbiotic `FLOW` $\langle P \rangle$ | Exact $\chi_t$ | Symbiotic `FLOW` $\chi_t$ | Local `Metro+OR` $\chi_t$ |
| ---: | ---: | :--- | ---: | :--- | :--- |
| `2.0` | `0.697775` | `0.697844 ± 0.000207` (`0.34σ`) | `0.019364` | `0.019576 ± 0.000194` (`1.09σ`) | `0.019651 ± 0.000220` |
| `4.0` | `0.863523` | `0.863439 ± 0.000104` (`0.80σ`) | `0.007554` | `0.007741 ± 0.000122` (`1.53σ`) | `0.007565 ± 0.000431` |
| `6.0` | `0.912359` | `0.912364 ± 0.000061` (`0.08σ`) | `0.004667` | `0.004734 ± 0.000065` (`1.04σ`) | `0.002520 ± 0.000561` |
| `8.0` | `0.935235` | `0.935216 ± 0.000045` (`0.42σ`) | `0.003399` | `0.003409 ± 0.000042` (`0.24σ`) | `0.000000` (`FROZEN`) |
| `10.0` | `0.948600` | `0.948615 ± 0.000035` (`0.42σ`) | `0.002676` | `0.002626 ± 0.000029` (`1.73σ`) | `0.000000` (`FROZEN`) |

## 2. Topological Autocorrelation Time $\tau_{\text{int}}(Q)$ (Madras–Sokal Automatic Window) & Instanton Acceptance

| $\beta$ | Local `Metro+OR` $\tau_{\text{int}}(Q)$ | Bare 3-Link Core Instanton | Symbiotic Harmonic `2-Jet` `FLOW` |
| ---: | :--- | :--- | :--- |
| `2.0` | `3.32` (`38499` trans) | `3.56` (`acc=7.7620%`) | **`2.62`** (`39703` trans, `acc=81.62%`) |
| `4.0` | `81.14` (`2895` trans) | `120.77` (`acc=0.5160%`) | **`5.26`** (`35274` trans, `acc=71.14%`) |
| `6.0` | `1173.69` (`54` trans) | `1411.01` (`acc=0.0340%`) | **`3.39`** (`31805` trans, `acc=63.63%`) |
| `8.0` | **FROZEN** (`∞`) (`0` trans) | `2350.07` (`acc=0.0040%`) | **`2.67`** (`28489` trans, `acc=56.98%`) |
| `10.0` | **FROZEN** (`∞`) (`0` trans) | **FROZEN** (`∞`) (`acc=0.0000%`) | **`2.21`** (`26073` trans, `acc=52.15%`) |

## 3. Structural Takeaways Linked to Lean 4 (`Flow.GaugeBundle`)

1. **Real-Physics Theorem 1 (`luscher_admissibility_preserves_topoCharge`, `luscher_iter_flow_frozen`)**: Smooth local flows (`local_metro_or`) preserve Lüscher plaquette admissibility at high $\beta$, causing topological sector lock-in (`FROZEN` or soaring $\tau_{\text{int}}(Q)$).
2. **Real-Physics Theorem 2 & UV Penalty (`bare_vs_harmonic_uv_penalty`)**: Bare local $Q = +1$ instantons concentrate curvature onto a $2\times 2$ core and pay an exact $L^2 / 4 = 64\times$ quadratic UV action penalty (`Theorem G4C`), driving acceptance toward $0$ at $\beta \ge 6$.
3. **Harmonic `2-Jet` Optimality + Detailed Balance (`harmonicInstanton_optimal`, `metropolis_detailed_balance_exact`)**: Distributing flux uniformly ($\delta_0 = 2\pi / L^2$) saturates the 2D Cauchy–Schwarz lower bound on $\mathbb{T}_L^2$, holding instanton acceptance high even at $\beta = 10$ with exact $\det J = 1$ zero-bias sampling.
4. **Crossover Diffusion Scaling ($\beta = 4.0$)**: At the $\beta = 4.0$ crossover where $\langle Q^2 \rangle \approx 1.93$ (`std(Q) ≈ 1.39`) and local tunneling already slows (`\tau_{\text{int}}(Q) = 81.14`), a single symmetric $\Delta Q = \pm 1$ harmonic step per sweep naturally matches random-walk sector diffusion $\tau_{\text{RW}} \approx 0.5 + 2\langle Q^2 \rangle / p_{\text{acc}} \approx 5.9$ (observed `5.26`, still **`15.4×`** faster than `local_metro_or`), before dropping cleanly to `3.39 → 2.67 → 2.21` in the deep freezing regime $\beta \ge 6.0$.
