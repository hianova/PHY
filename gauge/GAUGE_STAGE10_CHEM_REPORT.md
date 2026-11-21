# Stage 10 Report: First-Principles Certified Quantum Chemistry (`Flow.MolecularBundle`)

**Theorems Certified:** `Flow.MolecularBundle.M1–M5` (22 Lean 4 theorems, 0 `sorryAx`) + Independent McMurchie–Davidson/`numpy` Oracle (`< 1e-10` Ha)

## 1. Startup Self-Tests (`Theorems M1–M5`)
- **JW CAR (`M=6`, 64 states):** `0.00e+00`, Particle-Number Cons `0.00e+00`, 2-Body Antisym `0.00e+00`
- **H2O 8-Fold ERI Symmetry & Schwarz (`2,401` quartets):** `sym_err = 0.00e+00`, `cs_max_viol = 0.00e+00`
- **H2O RHF Idempotency & `[F, P]` Commutator:** `idem_err = 5.12e-15`, `comm_norm = 9.35e-13`
- **H4 FCI CAR `D->(Q,G)` Maps & PSD:** `car_map_err = 9.58e-16`, `eq_err = 1.60e-15`, `min_eig = -3.32e-16`
- **H2 Virial Scaling (`ζ = 0.96487944`):** `-V/T = 2.000000000000`, `2*E_H = -0.93316370 Ha`

## 2. Independent McMurchie–Davidson vs. C17 Obara–Saika Oracle
| System | $|\Delta E_{\text{RHF}}|$ (Ha) | $|\Delta E_{\text{MP2}}|$ (Ha) | $|\Delta E_{\text{FCI}}|$ (Ha) |
| :--- | :---: | :---: | :---: |
| H2 (R=1.4) | `4.47e-13` | `2.00e-13` | `3.31e-13` |
| H2 (R=5.0) | `3.27e-13` | `1.59e-13` | `1.40e-13` |
| HeH+ (R=1.4632) | `3.69e-14` | `9.33e-14` | `2.78e-13` |
| H4 (Ry=2.0 FCI) | `sym-brk` ($D_{4h}$) | `sym-brk` ($D_{4h}$) | `1.27e-13` |
| H4 (Ry=2.2) | `3.48e-13` | `2.67e-13` | `3.96e-13` |
| H4 (Ry=4.5) | `2.18e-13` | `3.81e-13` | `2.22e-13` |
| LiH (225-det FCI) | `2.42e-13` | `2.31e-14` | `1.34e-13` |
| H2O (441-det FCI) | `3.55e-13` | `5.97e-13` | `1.99e-13` |

## 3. Exp 1: H2 Dissociation & HeH+ Certified Two-Sided Brackets (Ha)
| System / $R$ ($a_0$) | $E_{\text{DQG cert-LB}}$ | $E_{\text{FCI}}$ | $E_{\text{RHF}}$ | $E_{\text{Temple}}$ | Gap $E_{\text{FCI}}-E_{\text{LB}}$ (mHa) | Orth Err |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| H2 `R=0.8` | `-0.95759853` | `-0.95759853` | `-0.94730793` | `-0.97145058` | `0.0000` | `2.21e-16` |
| H2 `R=1.0` | `-1.07896977` | `-1.07896977` | `-1.06599946` | `-1.09764839` | `0.0000` | `1.43e-16` |
| H2 `R=1.2` | `-1.12669883` | `-1.12669882` | `-1.11033388` | `-1.15231676` | `0.0000` | `2.34e-16` |
| H2 `R=1.4` | `-1.13727594` | `-1.13727594` | `-1.11671433` | `-1.17288470` | `0.0000` | `6.83e-16` |
| H2 `R=1.8` | `-1.11084598` | `-1.11084598` | `-1.07898221` | `-1.18311486` | `0.0000` | `4.12e-16` |
| H2 `R=2.4` | `-1.04148934` | `-1.04148934` | `-0.98269933` | `-1.32400763` | `0.0000` | `8.90e-17` |
| H2 `R=3.2` | `-0.97214738` | `-0.97214738` | `-0.85609692` | N/A ($E_{\text{HF}}\ge E_1$) | `0.0000` | `3.81e-16` |
| H2 `R=5.0` | `-0.93488936` | `-0.93488935` | `-0.68641593` | N/A ($E_{\text{HF}}\ge E_1$) | `0.0000` | `1.84e-16` |
| H2 `R=8.0` | `-0.93316639` | `-0.93316638` | `-0.61003890` | N/A ($E_{\text{HF}}\ge E_1$) | `0.0000` | `1.52e-16` |
| HeH+ `R=1.4632` | `-2.88070820` | `-2.88070820` | `-2.86065849` | `-2.91738612` | `0.0000` | `1.96e-16` |

## 4. Exp 2: H4 Square-to-Rectangle Strong Correlation (`36`-det FCI vs `DQG+S^2` Dual SDP)
| $R_y$ ($a_0$) | $E_{\text{DQG cert-LB}}$ | $E_{\text{FCI}}$ | $E_{\text{MP2}}$ | $E_{\text{RHF}}$ | $E_{\text{FCI}}-E_{\text{LB}}$ (mHa) | RHF Err (mHa) | Orth Err |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `2.0` | `-1.94044920` | `-1.93943162` | `-1.73600190` | `-1.70736010` | `1.0176` | `232.0715` | `1.57e-15` |
| `2.2` | `-1.97947801` | `-1.97770201` | `-1.91967648` | `-1.86021052` | `1.7760` | `117.4915` | `4.94e-15` |
| `2.6` | `-2.06484641` | `-2.06395840` | `-2.02290265` | `-1.96895873` | `0.8880` | `94.9997` | `1.17e-14` |
| `3.2` | `-2.13354313` | `-2.13324890` | `-2.09740406` | `-2.04811144` | `0.2942` | `85.1375` | `2.06e-15` |
| `4.5` | `-2.17256055` | `-2.17233524` | `-2.13874375` | `-2.09276362` | `0.2253` | `79.5716` | `1.29e-15` |

## 5. Exp 3: Polyatomic `s+p` Chemistry (`LiH` `225`-det & `H2O` `441`-det) & `H2O` `2-Jet` Newton Opt
- **LiH (`R=3.015`, `K=6`, `225`-det FCI):** `E_Temple = -8.07107425 <= E_FCI = -7.88239509 <= E_MP2 = -7.87487773 <= E_RHF = -7.86200941` (`DIIS=62 iters`, `comm=9.04e-13`, `Schwarz screened=0.0%`, `-V/T=1.988178`)
- **H2O (`R=1.809, θ=104.52°`, `K=7`, `441`-det FCI):** `E_Temple = -75.25299319 <= E_FCI(all) = -75.01242615 < E_FCI(fc) = -75.01234800 <= E_MP2 = -74.99843975 <= E_RHF = -74.96294062` (`DIIS=54 iters`, `comm=9.35e-13`, `Schwarz screened=0.0%`, `-V/T=2.005018`)
- **H2O 2D Geometry Optimization `(2.20, 125°) -> eq`:** `2-Jet Newton` converged in **`8` steps** to `R_eq=1.869713 a0, θ_eq=100.0269°, E=-74.96590179 Ha` (`|∇E|=4.26e-11`), vs Gradient Descent in `8` steps (`R=1.845492, θ=124.2651°, E=-74.94313906, |∇E|=1.81e-03`).

## 6. Exp 4: Ab Initio H4 Bilevel Conformer Tunneling (`SEARU × 2-Jet`, `R0 = 0.72`)
- **Global Ground Hits:** `Pure 2-Jet = 1/4` vs **`Symbiotic SEARU x 2-Jet = 4/4`** (`raw_false_rejects = 3/4`)
| Sector | Name | $E_{\text{raw}}(R_0)$ (Ha) | $R_{\text{eq}}$ ($a_0$) | $E_{\text{relaxed}}$ (Ha) |
| :---: | :--- | :---: | :---: | :---: |
| `0` | Linear | `-0.84650273` | `1.678481` | `-2.18050559` |
| `1` | Square | `0.67379151` | `2.413848` | `-1.97059898` |
| `2` | Rhombus | `0.64545243` | `2.139894` | `-1.92567605` |
| `3` | 2xH2_Dimer | `-1.74486719` | `1.388688` | `-2.27461092` |

ALL STAGE 10 FIRST-PRINCIPLES CERTIFIED QUANTUM CHEMISTRY CHECKS PASSED.
