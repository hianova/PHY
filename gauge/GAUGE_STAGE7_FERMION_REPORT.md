# Stage 7 Report: 2D/4D Lattice Fermions, Atiyah–Singer Index & Sign Problem (`Theorems F1–F5`)

Generated from [`build/fermion_2d_bench.json`](build/fermion_2d_bench.json) via [`fermion_2d_windtunnel.c`](fermion_2d_windtunnel.c) and audited against [`fermion_2d_oracle.py`](fermion_2d_oracle.py).

## 1. Startup Mathematical & Algebraic Self-Tests (`Flow.FermionBundle`)

| Check | Theorem Link | Measured Residual / Value |
| :--- | :--- | :--- |
| 64-Dim ($n=6$) Jordan–Wigner CAR, Pauli & $n_k^2$ | `Theorems F1A–F1F` | `max_err = 0.00e+00 < 1e-14` |
| 2D & 4D Euclidean Clifford Algebra + Wilson `γ_5H` & `U(1)` Cov | `Theorems F2A–F2C` | `Cliff = 0.00e+00`, `γ_5H = 0.00e+00`, `cov = 3.24e-16` |
| Nielsen–Ninomiya 2D Doubling / Wilson Lifting + 4D Staggered $\pi$-Flux | `Theorems F3A–F3D` | Naive zeros `4` -> Wilson zeros `1` (`m_min = 2.000`), 4D `pi` violations `0` |
| $6\times 6$ Overlap GW Relation, Normality & Unit Circle | `Theorems F4A–F4C` | `GW = 2.35e-14`, `norm = 4.04e-15`, `circle = 2.40e-14` |
| $6\times 6$ Staggered HMC Exact Fermion Force vs. Central FD | `Theorem F2D` | `max|F_ana - F_FD| = 1.13e-08 < 1e-6` |
| Lefschetz Quadratic Saddle Cancellation | `Theorem F5E` | `max|Im S| = 0.00e+00 < 1e-14` |

## 2. Experiment 1 — $8\times 8$ Discrete Atiyah–Singer Index Theorem (`Theorem F4`, `n_+ - n_- = -Q_geo`)

| Background | `Q_geo` | Zero Modes `(n_+, n_-)` | `index(D_ov) = n_+ - n_-` | Max Non-Zero Chirality `|⟨v, γ_5 v⟩|` |
| :--- | :---: | :---: | :---: | :---: |
| `smooth` | `-2` | `(2, 0)` | **`+2`** (`tr=+2`) | `2.93e-14` |
| `noisy` | `-2` | `(2, 0)` | **`+2`** (`tr=+2`) | `8.21e-15` |
| `smooth` | `-1` | `(1, 0)` | **`+1`** (`tr=+1`) | `1.36e-14` |
| `noisy` | `-1` | `(1, 0)` | **`+1`** (`tr=+1`) | `9.80e-15` |
| `smooth` | `+0` | `(0, 0)` | **`+0`** (`tr=+0`) | `2.37e-14` |
| `noisy` | `+0` | `(0, 0)` | **`+0`** (`tr=+0`) | `8.03e-15` |
| `smooth` | `+1` | `(0, 1)` | **`-1`** (`tr=-1`) | `1.08e-14` |
| `noisy` | `+1` | `(0, 1)` | **`-1`** (`tr=-1`) | `5.30e-15` |
| `smooth` | `+2` | `(0, 2)` | **`-2`** (`tr=-2`) | `8.54e-15` |
| `noisy` | `+2` | `(0, 2)` | **`-2`** (`tr=-2`) | `7.30e-15` |

## 3. Experiment 2 — $6\times 6$ Dynamical $N_f = 2$ Staggered Schwinger Model vs. Quenched (`m = 0.15`)

| $\beta$ | Quenched `<P>` | Dynamical $N_f=2$ `<P>` | Screening $\chi(2,2)$ (`Quench -> Dyn`) | Dynamical $\langle \bar\psi\psi \rangle$ (`acc`) |
| ---: | :---: | :---: | :---: | :---: |
| `1.5` | `0.60319 ± 0.00520` | **`0.70221 ± 0.00365`** | `0.5100 ± 0.0376` -> **`0.2321 ± 0.0148`** | `0.37544 ± 0.00744` (`acc=98.8%`) |
| `2.5` | `0.76642 ± 0.00359` | **`0.81573 ± 0.00246`** | `0.2631 ± 0.0150` -> **`0.1238 ± 0.0043`** | `0.27616 ± 0.00457` (`acc=98.5%`) |
| `4.0` | `0.86420 ± 0.00224` | **`0.88479 ± 0.00253`** | `0.1235 ± 0.0091` -> **`0.0806 ± 0.0041`** | `0.24855 ± 0.00279` (`acc=97.6%`) |

### Chiral Condensate Mass Scan at $\beta = 2.0$ ($6\times 6$ $N_f=2$ vs. Free Field)

| Mass $m$ | `<P>` | Dynamical $\langle \bar\psi\psi \rangle$ | Free Staggered $\langle \bar\psi\psi \rangle_{\text{free}}$ | $\chi(2,2)$ |
| ---: | :---: | :---: | :---: | :---: |
| `0.05` | `0.78123` | **`0.11046 ± 0.00525`** | `0.07805` | `0.1457` |
| `0.15` | `0.77049` | **`0.31996 ± 0.00817`** | `0.22263` | `0.1751` |
| `0.30` | `0.75177` | **`0.50231 ± 0.00623`** | `0.38522` | `0.2376` |
| `0.50` | `0.72361` | **`0.58154 ± 0.00147`** | `0.50000` | `0.3068` |

## 4. Experiment 3A — 1D Finite-Density $U(1)$ Sign Problem vs. Positive Worldline Dual (`β=1.0, μ=0.6`)

| Chain $L$ | Naive PQ `<sgn>` | Naive Reweighted `<n>` | Dual `<sgn>` | Positive Dual `<n>` (vs. Exact Bessel) |
| ---: | :---: | :---: | :---: | :--- |
| `8` | `0.58746 ± 0.00437` (`exact=0.58127`) | `+0.1657 ± 0.0039` | **`1.000000`** | **`+0.16134 ± 0.00106`** (`exact=+0.16116`, `0.17σ`) |
| `16` | `0.25413 ± 0.00452` (`exact=0.25099`) | `+0.0408 ± 0.0078` | **`1.000000`** | **`+0.03514 ± 0.00053`** (`exact=+0.03540`, `0.49σ`) |
| `32` | `0.06168 ± 0.00290` (`exact=0.05870`) | `+0.0060 ± 0.0169` | **`1.000000`** | **`+0.00114 ± 0.00010`** (`exact=+0.00135`, `2.09σ`) |
| `64` | `0.00260 ± 0.00313` (`exact=0.00344`) | `-0.2550 ± 0.7446` | **`1.000000`** | **`+0.00001 ± 0.00001`** (`exact=+0.00000`, `0.78σ`) |

## 5. Experiment 3B — Non-Gaussian Complex Action & Lefschetz Thimble (`a=2.0, λ=0.05`)

| $\mu$ (`y_0`) | Naive Real-Axis `<sgn>` & `<z^2>` | Thimble `<sgn>` | Thimble `<z^2>` (vs. Exact Quad) | Thimble `<Im z>` (vs. Exact Quad) |
| :--- | :--- | :---: | :--- | :--- |
| `0.5` (`y0=0.2504`) | `0.94140` (`<z^2>=+0.426±0.001`) | **`0.99987`** | **`+0.4264 ± 0.0014`** (`ex=+0.4263`, `0.05σ`) | **`+0.2417 ± 0.0001`** (`ex=+0.2418`, `0.47σ`) |
| `1.5` (`y0=0.7610`) | `0.57962` (`<z^2>=-0.037±0.002`) | **`0.99872`** | **`-0.0375 ± 0.0015`** (`ex=-0.0367`, `0.55σ`) | **`+0.7330 ± 0.0002`** (`ex=+0.7330`, `0.57σ`) |
| `2.5` (`y0=1.3056`) | `0.21466` (`<z^2>=-1.027±0.014`) | **`0.99525`** | **`-1.0225 ± 0.0019`** (`ex=-1.0244`, `0.99σ`) | **`+1.2494 ± 0.0003`** (`ex=+1.2497`, `1.03σ`) |
| `3.5` (`y0=1.9296`) | `0.04682` (`<z^2>=-2.705±0.110`) | **`0.98421`** | **`-2.6990 ± 0.0025`** (`ex=-2.7001`, `0.44σ`) | **`+1.8197 ± 0.0004`** (`ex=+1.8198`, `0.35σ`) |
| `4.0` (`y0=2.3069`) | `0.01670` (`<z^2>=-4.132±0.403`) | **`0.96745`** | **`-3.9065 ± 0.0036`** (`ex=-3.9061`, `0.09σ`) | **`+2.1404 ± 0.0006`** (`ex=+2.1406`, `0.46σ`) |
