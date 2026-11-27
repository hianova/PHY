namespace Phy.HigherDimFrontier

/-!
  PHY v4.1.0 — Higher-Dimensional Frontier Formalizations (Pure-`Init` Lean 4.36.0).
  Zero imports outside Lean 4 `Init`; pure kernel proofs only.

  Covers three referee-certified discovery frontiers:
  1. Exact `|X| = 3` (27-state) non-involutive nonlinear 2+1D Zamolodchikov
     tetrahedron family `R_{c,a}` across all 729 6-channel inputs.
  2. `AME(4,5)` degree-6 `F_5` algebraic-rigidity counterexample via single-channel
     relabeling and `d* = 7` non-group orthogonal Latin square onset `(A7, B7)`.
  3. `AME(4,6)` exact Eisenstein integer `Z[omega]` quantum escape below the
     classical `2`-defect floor (`D_total = 3/2 < 2` at 40 nnz and `27/26 < 2` at 44 nnz).
-/

def Bijective {α β : Type} (f : α → β) : Prop :=
  (∀ a b, f a = f b → a = b) ∧ (∀ b, ∃ a, f a = b)

/-! ## Section 1: Exact `|X| = 3` Non-Involutive Nonlinear 2+1D Tetrahedron Family -/

abbrev State3 := Fin 3 × Fin 3 × Fin 3
abbrev State6_3 := Fin 3 × Fin 3 × Fin 3 × Fin 3 × Fin 3 × Fin 3

def tetraR3 (c a : Fin 3) (s : State3) : State3 :=
  (s.1, s.2.1 + c * (s.1 - a) * (s.2.2 - a), s.2.2)

private theorem fin3_cases (x : Fin 3) : x = 0 ∨ x = 1 ∨ x = 2 := by
  rcases x with ⟨v, hv⟩
  match v, hv with
  | 0, _ => exact Or.inl rfl
  | 1, _ => exact Or.inr (Or.inl rfl)
  | 2, _ => exact Or.inr (Or.inr rfl)

/-- Exact order 3 (`R^3 = id`) for the entire `R_{c,a}` family over `F_3`. -/
theorem tetraR3_order3 (c a : Fin 3) (s : State3) :
    tetraR3 c a (tetraR3 c a (tetraR3 c a s)) = s := by
  rcases s with ⟨x, y, z⟩
  simp only [tetraR3]
  have h3 : ∀ y0 u : Fin 3, y0 + u + u + u = y0 := by decide
  rw [h3 y (c * (x - a) * (z - a))]

/-- Strict non-involutivity witness (`R^2 ≠ id`) for every `c ∈ {1, 2}` and `a ∈ {0, 1, 2}`. -/
theorem tetraR3_noninvolutive (c a : Fin 3) (hc : c ≠ 0) :
    ∃ s : State3, tetraR3 c a (tetraR3 c a s) ≠ s := by
  refine ⟨(a + 1, 0, a + 1), ?_⟩
  rcases fin3_cases c with rfl | rfl | rfl
  · exact False.elim (hc rfl)
  · rcases fin3_cases a with rfl | rfl | rfl <;> decide
  · rcases fin3_cases a with rfl | rfl | rfl <;> decide

/-- Universal non-involutivity (`¬ ∀ s, R^2(s) = s`) for every `c ≠ 0` and `a : Fin 3`. -/
theorem tetraR3_not_involutive (c a : Fin 3) (hc : c ≠ 0) :
    ¬ (∀ s : State3, tetraR3 c a (tetraR3 c a s) = s) := by
  intro h
  rcases tetraR3_noninvolutive c a hc with ⟨s, hs⟩
  exact hs (h s)

/-- Unconditional bijectivity of `R_{c,a}` via its explicit inverse `R_{c,a}^2`. -/
theorem tetraR3_bijective (c a : Fin 3) : Bijective (tetraR3 c a) := by
  refine ⟨fun s1 s2 h => ?_, fun s => ⟨tetraR3 c a (tetraR3 c a s), tetraR3_order3 c a s⟩⟩
  have h2 := congrArg (fun t => tetraR3 c a (tetraR3 c a t)) h
  simpa [tetraR3_order3 c a] using h2

def secondDiffXZ (f : State3 → Fin 3) (x y z : Fin 3) : Fin 3 :=
  f (x + 1, y, z + 1) - f (x + 1, y, z) - f (x, y, z + 1) + f (x, y, z)

def IsF3Affine3 (f : State3 → Fin 3) : Prop :=
  ∃ k0 kx ky kz : Fin 3, ∀ x y z : Fin 3, f (x, y, z) = k0 + kx * x + ky * y + kz * z

/-- Exact nonzero 2nd-order finite difference `Δ_x Δ_z (R_{c,a})_y = c`. -/
theorem tetraR3_secondDiffXZ (c a x y z : Fin 3) :
    secondDiffXZ (fun s => (tetraR3 c a s).2.1) x y z = c := by
  revert c a x y z
  decide

private theorem affine3_secondDiffXZ_zero (f : State3 → Fin 3) (hAff : IsF3Affine3 f) :
    secondDiffXZ f 0 0 0 = 0 := by
  rcases hAff with ⟨k0, kx, ky, kz, hf⟩
  simp only [secondDiffXZ, hf]
  have h0 : ∀ a0 ax ay az : Fin 3,
      (a0 + ax * (0 + 1) + ay * 0 + az * (0 + 1)) -
      (a0 + ax * (0 + 1) + ay * 0 + az * 0) -
      (a0 + ax * 0 + ay * 0 + az * (0 + 1)) +
      (a0 + ax * 0 + ay * 0 + az * 0) = 0 := by decide
  exact h0 k0 kx ky kz

/-- Strict `F_3` nonlinearity (`deg > 1`) of `R_{c,a}` for every `c ∈ {1, 2}` and `a ∈ {0, 1, 2}`. -/
theorem tetraR3_nonlinear (c a : Fin 3) (hc : c ≠ 0) :
    ¬ IsF3Affine3 (fun s => (tetraR3 c a s).2.1) := by
  intro hAff
  have h0 := affine3_secondDiffXZ_zero (fun s => (tetraR3 c a s).2.1) hAff
  rw [tetraR3_secondDiffXZ c a 0 0 0] at h0
  exact hc h0

def Is3ChannelCoupled3 (R : State3 → State3) : Prop :=
  (∃ x0 x1 y z : Fin 3, (R (x1, y, z)).1 ≠ (R (x0, y, z)).1) ∧
  (∃ x0 x1 y z : Fin 3, (R (x1, y, z)).2.1 ≠ (R (x0, y, z)).2.1) ∧
  (∃ x y0 y1 z : Fin 3, (R (x, y1, z)).2.1 ≠ (R (x, y0, z)).2.1) ∧
  (∃ x y z0 z1 : Fin 3, (R (x, y, z1)).2.1 ≠ (R (x, y, z0)).2.1) ∧
  (∃ x y z0 z1 : Fin 3, (R (x, y, z1)).2.2 ≠ (R (x, y, z0)).2.2)

/-- Genuine 3-channel coupling (all 5 active channel dependencies) for every `c ≠ 0`. -/
theorem tetraR3_3channel_coupled (c a : Fin 3) (hc : c ≠ 0) :
    Is3ChannelCoupled3 (tetraR3 c a) := by
  refine ⟨⟨0, 1, 0, 0, by revert c a; decide⟩,
          ⟨a, a + 1, 0, a + 1, ?_⟩,
          ⟨0, 0, 1, 0, ?_⟩,
          ⟨a + 1, 0, a, a + 1, ?_⟩,
          ⟨0, 0, 0, 1, by revert c a; decide⟩⟩
  · rcases fin3_cases c with rfl | rfl | rfl
    · exact False.elim (hc rfl)
    · rcases fin3_cases a with rfl | rfl | rfl <;> decide
    · rcases fin3_cases a with rfl | rfl | rfl <;> decide
  · rcases fin3_cases c with rfl | rfl | rfl <;>
      rcases fin3_cases a with rfl | rfl | rfl <;> decide
  · rcases fin3_cases c with rfl | rfl | rfl
    · exact False.elim (hc rfl)
    · rcases fin3_cases a with rfl | rfl | rfl <;> decide
    · rcases fin3_cases a with rfl | rfl | rfl <;> decide

def applyT012 (R : State3 → State3) (s : State6_3) : State6_3 :=
  let (x1, x2, x3, x4, x5, x6) := s
  let (y1, y2, y3) := R (x1, x2, x3)
  (y1, y2, y3, x4, x5, x6)

def applyT034 (R : State3 → State3) (s : State6_3) : State6_3 :=
  let (x1, x2, x3, x4, x5, x6) := s
  let (y1, y4, y5) := R (x1, x4, x5)
  (y1, x2, x3, y4, y5, x6)

def applyT135 (R : State3 → State3) (s : State6_3) : State6_3 :=
  let (x1, x2, x3, x4, x5, x6) := s
  let (y2, y4, y6) := R (x2, x4, x6)
  (x1, y2, x3, y4, x5, y6)

def applyT245 (R : State3 → State3) (s : State6_3) : State6_3 :=
  let (x1, x2, x3, x4, x5, x6) := s
  let (y3, y5, y6) := R (x3, x5, x6)
  (x1, x2, y3, x4, y5, y6)

def IsSetTheoreticTetra3 (R : State3 → State3) : Prop :=
  ∀ x1 x2 x3 x4 x5 x6 : Fin 3,
    applyT012 R (applyT034 R (applyT135 R (applyT245 R (x1, x2, x3, x4, x5, x6)))) =
      applyT245 R (applyT135 R (applyT034 R (applyT012 R (x1, x2, x3, x4, x5, x6))))

/-- Exact 6-channel Zamolodchikov Tetrahedron equation across all `729` states for all `(c, a)`. -/
theorem tetraR3_satisfies_tetra (c a : Fin 3) :
    IsSetTheoreticTetra3 (tetraR3 c a) := by
  intro x1 x2 x3 x4 x5 x6
  simp only [applyT012, applyT034, applyT135, applyT245, tetraR3]
  have hch4 : x4 + c * (x2 - a) * (x6 - a) + c * (x1 - a) * (x5 + c * (x3 - a) * (x6 - a) - a) =
      x4 + c * (x1 - a) * (x5 - a) + c * (x2 + c * (x1 - a) * (x3 - a) - a) * (x6 - a) := by
    revert c a x1 x2 x3 x4 x5 x6
    decide
  simp only [hch4]

/-! ## Section 2: `AME(4,5)` Degree-6 Counterexample & `d* = 7` Non-Group OLS Onset -/

def IsLatin (d : Nat) (M : Fin d → Fin d → Fin d) : Prop :=
  (∀ r c1 c2 : Fin d, M r c1 = M r c2 → c1 = c2) ∧
  (∀ r v : Fin d, ∃ c : Fin d, M r c = v) ∧
  (∀ c r1 r2 : Fin d, M r1 c = M r2 c → r1 = r2) ∧
  (∀ c v : Fin d, ∃ r : Fin d, M r c = v)

def IsOrthogonalPair (d : Nat) (A B : Fin d → Fin d → Fin d) : Prop :=
  (∀ r1 c1 r2 c2 : Fin d, A r1 c1 = A r2 c2 → B r1 c1 = B r2 c2 → r1 = r2 ∧ c1 = c2) ∧
  (∀ u v : Fin d, ∃ r c : Fin d, A r c = u ∧ B r c = v)

def IsOLS5 (A B : Fin 5 → Fin 5 → Fin 5) : Prop :=
  IsLatin 5 A ∧ IsLatin 5 B ∧ IsOrthogonalPair 5 A B

def IsOLS7 (A B : Fin 7 → Fin 7 → Fin 7) : Prop :=
  IsLatin 7 A ∧ IsLatin 7 B ∧ IsOrthogonalPair 7 A B

def Is2UnitaryPair (d : Nat) (A B : Fin d → Fin d → Fin d) : Prop :=
  Bijective (fun p : Fin d × Fin d => (A p.1 p.2, B p.1 p.2)) ∧
  Bijective (fun p : Fin d × Fin d => (p.1, A p.1 p.2)) ∧
  Bijective (fun p : Fin d × Fin d => (p.2, A p.1 p.2)) ∧
  Bijective (fun p : Fin d × Fin d => (p.1, B p.1 p.2)) ∧
  Bijective (fun p : Fin d × Fin d => (p.2, B p.1 p.2))

theorem isLatin_of_inverses {d : Nat} (M rowInv colInv : Fin d → Fin d → Fin d)
    (hRowL : ∀ r c, rowInv r (M r c) = c)
    (hRowR : ∀ r v, M r (rowInv r v) = v)
    (hColL : ∀ r c, colInv c (M r c) = r)
    (hColR : ∀ c v, M (colInv c v) c = v) : IsLatin d M := by
  refine ⟨fun r c1 c2 h => ?_, fun r v => ⟨rowInv r v, hRowR r v⟩,
          fun c r1 r2 h => ?_, fun c v => ⟨colInv c v, hColR c v⟩⟩
  · have h1 := congrArg (rowInv r) h
    simpa [hRowL] using h1
  · have h1 := congrArg (colInv c) h
    simpa [hColL] using h1

theorem isOrthogonal_of_inverse {d : Nat} (A B invR invC : Fin d → Fin d → Fin d)
    (hInvR : ∀ r c, invR (A r c) (B r c) = r)
    (hInvC : ∀ r c, invC (A r c) (B r c) = c)
    (hSurjA : ∀ u v, A (invR u v) (invC u v) = u)
    (hSurjB : ∀ u v, B (invR u v) (invC u v) = v) : IsOrthogonalPair d A B := by
  refine ⟨fun r1 c1 r2 c2 hA hB => ?_, fun u v => ⟨invR u v, invC u v, hSurjA u v, hSurjB u v⟩⟩
  constructor
  · have h1 := congrArg (fun x => invR x (B r1 c1)) hA
    rw [hInvR r1 c1, hB, hInvR r2 c2] at h1
    exact h1
  · have h2 := congrArg (fun x => invC x (B r1 c1)) hA
    rw [hInvC r1 c1, hB, hInvC r2 c2] at h2
    exact h2

/-- Every orthogonal Latin square pair `(A, B)` induces a 2-unitary permutation across all 2-to-2 partitions. -/
theorem isOLS_implies_2unitary {d : Nat} (A B : Fin d → Fin d → Fin d)
    (hLA : IsLatin d A) (hLB : IsLatin d B) (hOrth : IsOrthogonalPair d A B) :
    Is2UnitaryPair d A B := by
  refine ⟨⟨?_, ?_⟩, ⟨?_, ?_⟩, ⟨?_, ?_⟩, ⟨?_, ?_⟩, ⟨?_, ?_⟩⟩
  · intro p q h
    rcases Prod.ext_iff.mp h with ⟨hA, hB⟩
    rcases hOrth.1 p.1 p.2 q.1 q.2 hA hB with ⟨hr, hc⟩
    exact Prod.ext hr hc
  · intro ⟨u, v⟩
    rcases hOrth.2 u v with ⟨r, c, hr, hc⟩
    exact ⟨(r, c), Prod.ext hr hc⟩
  · intro ⟨p1, p2⟩ ⟨q1, q2⟩ h
    rcases Prod.ext_iff.mp h with ⟨h1, hA⟩
    dsimp at h1 hA
    subst h1
    exact Prod.ext rfl (hLA.1 p1 p2 q2 hA)
  · intro ⟨r, v⟩
    rcases hLA.2.1 r v with ⟨c, hc⟩
    exact ⟨(r, c), Prod.ext rfl hc⟩
  · intro ⟨p1, p2⟩ ⟨q1, q2⟩ h
    rcases Prod.ext_iff.mp h with ⟨h2, hA⟩
    dsimp at h2 hA
    subst h2
    exact Prod.ext (hLA.2.2.1 p2 p1 q1 hA) rfl
  · intro ⟨c, v⟩
    rcases hLA.2.2.2 c v with ⟨r, hr⟩
    exact ⟨(r, c), Prod.ext rfl hr⟩
  · intro ⟨p1, p2⟩ ⟨q1, q2⟩ h
    rcases Prod.ext_iff.mp h with ⟨h1, hB⟩
    dsimp at h1 hB
    subst h1
    exact Prod.ext rfl (hLB.1 p1 p2 q2 hB)
  · intro ⟨r, v⟩
    rcases hLB.2.1 r v with ⟨c, hc⟩
    exact ⟨(r, c), Prod.ext rfl hc⟩
  · intro ⟨p1, p2⟩ ⟨q1, q2⟩ h
    rcases Prod.ext_iff.mp h with ⟨h2, hB⟩
    dsimp at h2 hB
    subst h2
    exact Prod.ext (hLB.2.2.1 p2 p1 q1 hB) rfl
  · intro ⟨c, v⟩
    rcases hLB.2.2.2 c v with ⟨r, hr⟩
    exact ⟨(r, c), Prod.ext rfl hr⟩

def pi23 (x : Fin 5) : Fin 5 :=
  match x.val with
  | 0 => 0 | 1 => 1 | 2 => 2 | 3 => 4 | _ => 3

def A23_5 (r c : Fin 5) : Fin 5 :=
  match r.val, c.val with
  | 0, 0 => 0 | 0, 1 => 1 | 0, 2 => 2 | 0, 3 => 3 | 0, 4 => 4
  | 1, 0 => 1 | 1, 1 => 2 | 1, 2 => 4 | 1, 3 => 0 | 1, 4 => 3
  | 2, 0 => 2 | 2, 1 => 4 | 2, 2 => 3 | 2, 3 => 1 | 2, 4 => 0
  | 3, 0 => 3 | 3, 1 => 0 | 3, 2 => 1 | 3, 3 => 4 | 3, 4 => 2
  | _, 0 => 4 | _, 1 => 3 | _, 2 => 0 | _, 3 => 2 | _, _ => 1

def B23_5 (r c : Fin 5) : Fin 5 :=
  match r.val, c.val with
  | 0, 0 => 0 | 0, 1 => 1 | 0, 2 => 2 | 0, 3 => 3 | 0, 4 => 4
  | 1, 0 => 2 | 1, 1 => 4 | 1, 2 => 3 | 1, 3 => 1 | 1, 4 => 0
  | 2, 0 => 3 | 2, 1 => 0 | 2, 2 => 1 | 2, 3 => 4 | 2, 4 => 2
  | 3, 0 => 4 | 3, 1 => 3 | 3, 2 => 0 | 3, 3 => 2 | 3, 4 => 1
  | _, 0 => 1 | _, 1 => 2 | _, 2 => 4 | _, 3 => 0 | _, _ => 3

private def rowInvA23 (r v : Fin 5) : Fin 5 := pi23 (pi23 v - pi23 r)
private def colInvA23 (c v : Fin 5) : Fin 5 := pi23 (pi23 v - pi23 c)
private def rowInvB23 (r v : Fin 5) : Fin 5 := pi23 (pi23 v - 2 * pi23 r)
private def colInvB23 (c v : Fin 5) : Fin 5 := pi23 (3 * (pi23 v - pi23 c))
private def invR23 (u v : Fin 5) : Fin 5 := pi23 (pi23 v - pi23 u)
private def invC23 (u v : Fin 5) : Fin 5 := pi23 (2 * pi23 u - pi23 v)

/-- Exact single-channel relabeling formula for the degree-6 `AME(4,5)` counterexample `(A23_5, B23_5)`. -/
theorem ame45_deg6_relabeling (x y : Fin 5) :
    A23_5 x y = pi23 (pi23 x + pi23 y) ∧
    B23_5 x y = pi23 (2 * pi23 x + pi23 y) := by
  revert x y
  decide

/-- `(A23_5, B23_5)` is a valid 5×5 orthogonal Latin square pair (`IsOLS5`). -/
theorem ame45_deg6_is_ols5 : IsOLS5 A23_5 B23_5 := by
  refine ⟨isLatin_of_inverses A23_5 rowInvA23 colInvA23
            (by decide) (by decide) (by decide) (by decide),
          isLatin_of_inverses B23_5 rowInvB23 colInvB23
            (by decide) (by decide) (by decide) (by decide),
          isOrthogonal_of_inverse A23_5 B23_5 invR23 invC23
            (by decide) (by decide) (by decide) (by decide)⟩

/-- `(A23_5, B23_5)` defines a 2-unitary permutation on `Fin 5 × Fin 5`. -/
theorem ame45_deg6_2unitary : Is2UnitaryPair 5 A23_5 B23_5 :=
  have h := ame45_deg6_is_ols5
  isOLS_implies_2unitary A23_5 B23_5 h.1 h.2.1 h.2.2

def IsF5Affine2 (f : Fin 5 → Fin 5 → Fin 5) : Prop :=
  ∃ k0 kx ky : Fin 5, ∀ x y : Fin 5, f x y = k0 + kx * x + ky * y

/-- `A23_5` has `F_5` algebraic degree `> 1` (not `F_5`-affine). -/
theorem ame45_deg6_not_affine : ¬ IsF5Affine2 A23_5 := by
  rintro ⟨k0, kx, ky, hf⟩
  have h00 := hf 0 0
  have h10 := hf 1 0
  have h01 := hf 0 1
  have h22 := hf 2 2
  have hk0 : k0 = 0 := by
    have h0 : A23_5 0 0 = 0 := rfl
    rw [h0] at h00
    revert k0 kx ky h00
    decide
  subst hk0
  have hkx : kx = 1 := by
    have h1 : A23_5 1 0 = 1 := rfl
    rw [h1] at h10
    revert kx ky h10
    decide
  subst hkx
  have hky : ky = 1 := by
    have h2 : A23_5 0 1 = 1 := rfl
    rw [h2] at h01
    revert ky h01
    decide
  subst hky
  revert h22
  decide

def sumFin5 (f : Fin 5 → Fin 5) : Fin 5 :=
  f 0 + f 1 + f 2 + f 3 + f 4

def finPow5 (x n : Fin 5) : Fin 5 :=
  match n.val with
  | 0 => 1
  | 1 => x
  | 2 => x * x
  | 3 => x * x * x
  | _ => x * x * x * x

def evalBivarPoly5 (c : Fin 5 → Fin 5 → Fin 5) (x y : Fin 5) : Fin 5 :=
  sumFin5 fun i => sumFin5 fun j => c i j * finPow5 x i * finPow5 y j

def polyA23Coeff (i j : Fin 5) : Fin 5 :=
  match i.val, j.val with
  | 0, 1 => 1
  | 1, 0 => 1 | 1, 1 => 3 | 1, 2 => 3 | 1, 3 => 2
  | 2, 1 => 3 | 2, 2 => 1
  | 3, 1 => 2 | 3, 3 => 1
  | _, _ => 0

/-- Exact degree-6 bivariate `F_5` polynomial realization of `A23_5` (with `x^3 y^3` coefficient `1`). -/
theorem ame45_deg6_poly_exact (x y : Fin 5) :
    evalBivarPoly5 polyA23Coeff x y = A23_5 x y := by
  revert x y
  decide

private theorem fin5_mul_add (a b c : Fin 5) : a * (b + c) = a * b + a * c := by
  revert a b c; decide

private theorem fin5_add_mul (a b c : Fin 5) : (a + b) * c = a * c + b * c := by
  revert a b c; decide

private theorem fin5_mul_comm (a b : Fin 5) : a * b = b * a := by
  revert a b; decide

private theorem fin5_mul_assoc (a b c : Fin 5) : (a * b) * c = a * (b * c) := by
  revert a b c; decide

private theorem fin5_add_assoc (a b c : Fin 5) : (a + b) + c = a + (b + c) := by
  revert a b c; decide

private theorem fin5_add_comm (a b : Fin 5) : a + b = b + a := by
  revert a b; decide

private theorem sumFin5_const_mul (a : Fin 5) (f : Fin 5 → Fin 5) :
    sumFin5 (fun x => a * f x) = a * sumFin5 f := by
  simp only [sumFin5, fin5_mul_add]

private theorem sumFin5_mul_const (f : Fin 5 → Fin 5) (a : Fin 5) :
    sumFin5 (fun x => f x * a) = sumFin5 f * a := by
  simp only [sumFin5, fin5_add_mul]

private theorem sumFin5_add (f g : Fin 5 → Fin 5) :
    sumFin5 (fun x => f x + g x) = sumFin5 f + sumFin5 g := by
  have h4 : ∀ a b c d : Fin 5, (a + b) + (c + d) = (a + c) + (b + d) := by
    intro a b c d; revert a b c d; decide
  simp only [sumFin5]
  rw [h4 (f 0) (g 0) (f 1) (g 1)]
  rw [h4 (f 0 + f 1) (g 0 + g 1) (f 2) (g 2)]
  rw [h4 (f 0 + f 1 + f 2) (g 0 + g 1 + g 2) (f 3) (g 3)]
  rw [h4 (f 0 + f 1 + f 2 + f 3) (g 0 + g 1 + g 2 + g 3) (f 4) (g 4)]

private theorem sumFin5_comm (f : Fin 5 → Fin 5 → Fin 5) :
    sumFin5 (fun x => sumFin5 (fun y => f x y)) =
      sumFin5 (fun y => sumFin5 (fun x => f x y)) := by
  change sumFin5 (fun x => f x 0 + f x 1 + f x 2 + f x 3 + f x 4) =
    sumFin5 (fun x => f x 0) + sumFin5 (fun x => f x 1) + sumFin5 (fun x => f x 2) +
      sumFin5 (fun x => f x 3) + sumFin5 (fun x => f x 4)
  simp only [sumFin5_add]

private theorem fin5_pow_ortho (i : Fin 5) :
    sumFin5 (fun x => x * finPow5 x i) = if i == 3 then 4 else 0 := by
  revert i
  decide

private theorem fin5_extract_x_sum (g : Fin 5 → Fin 5) :
    sumFin5 (fun x => x * sumFin5 (fun i => g i * finPow5 x i)) = 4 * g 3 := by
  have h1 : ∀ x : Fin 5,
      x * sumFin5 (fun i => g i * finPow5 x i) = sumFin5 (fun i => g i * (x * finPow5 x i)) := by
    intro x
    rw [← sumFin5_const_mul]
    congr 1
    funext i
    rw [← fin5_mul_assoc, fin5_mul_comm x (g i), fin5_mul_assoc]
  have hz : ∀ a : Fin 5, a * (0 : Fin 5) = 0 := fun a => by revert a; decide
  have h2 : ∀ a : Fin 5, (0 + 0 + 0 + a + 0 : Fin 5) = a := by decide
  calc
    sumFin5 (fun x => x * sumFin5 (fun i => g i * finPow5 x i))
        = sumFin5 (fun x => sumFin5 (fun i => g i * (x * finPow5 x i))) := by
            simp only [h1]
    _ = sumFin5 (fun i => sumFin5 (fun x => g i * (x * finPow5 x i))) := sumFin5_comm _
    _ = sumFin5 (fun i => g i * sumFin5 (fun x => x * finPow5 x i)) := by
            congr 1; funext i; exact sumFin5_const_mul _ _
    _ = sumFin5 (fun i => g i * (if i == 3 then 4 else 0)) := by
            simp only [fin5_pow_ortho]
    _ = 4 * g 3 := by
            simp only [sumFin5]
            change (((g 0 * 0 + g 1 * 0) + g 2 * 0) + g 3 * 4) + g 4 * 0 = 4 * g 3
            simp only [hz]
            rw [h2 (g 3 * 4), fin5_mul_comm (g 3) 4]

/-- Unconditional lower bound: ANY `F_5` polynomial representing `A23_5` MUST have `c_{3,3} = 1 ≠ 0` (degree `≥ 6`). -/
theorem ame45_deg6_forced_deg6 (c : Fin 5 → Fin 5 → Fin 5)
    (hRep : ∀ x y : Fin 5, evalBivarPoly5 c x y = A23_5 x y) :
    c 3 3 = 1 := by
  let ext : (Fin 5 → Fin 5 → Fin 5) → Fin 5 :=
    fun f => sumFin5 fun x => sumFin5 fun y => x * y * f x y
  have hEq : ext (evalBivarPoly5 c) = ext A23_5 := by
    simp only [ext, hRep]
  have hA23 : ext A23_5 = 1 := by decide
  have hPoly : ext (evalBivarPoly5 c) = c 3 3 := by
    have hinner : ∀ x y : Fin 5,
        x * y * evalBivarPoly5 c x y =
          y * sumFin5 (fun j => (x * sumFin5 (fun i => c i j * finPow5 x i)) * finPow5 y j) := by
      intro x y
      unfold evalBivarPoly5
      rw [sumFin5_comm]
      have hmul : ∀ g : Fin 5 → Fin 5,
          x * y * sumFin5 (fun j => g j * finPow5 y j) =
            y * sumFin5 (fun j => (x * g j) * finPow5 y j) := by
        intro g
        rw [← sumFin5_const_mul, ← sumFin5_const_mul]
        congr 1
        funext j
        rw [fin5_mul_comm x y, fin5_mul_assoc, fin5_mul_assoc]
      have hr : (fun j => sumFin5 (fun i => c i j * finPow5 x i * finPow5 y j)) =
          (fun j => sumFin5 (fun i => c i j * finPow5 x i) * finPow5 y j) := by
        funext j
        exact sumFin5_mul_const (fun i => c i j * finPow5 x i) (finPow5 y j)
      rw [hr]
      exact hmul (fun j => sumFin5 (fun i => c i j * finPow5 x i))
    have h44 : ∀ u : Fin 5, 4 * (4 * u) = u := by decide
    calc
      ext (evalBivarPoly5 c)
          = sumFin5 (fun x => sumFin5 (fun y =>
              y * sumFin5 (fun j => (x * sumFin5 (fun i => c i j * finPow5 x i)) * finPow5 y j))) := by
              simp only [ext, hinner]
      _ = sumFin5 (fun x => 4 * (x * sumFin5 (fun i => c i 3 * finPow5 x i))) := by
              congr 1; funext x
              exact fin5_extract_x_sum (fun j => x * sumFin5 (fun i => c i j * finPow5 x i))
      _ = 4 * sumFin5 (fun x => x * sumFin5 (fun i => c i 3 * finPow5 x i)) :=
              sumFin5_const_mul 4 _
      _ = 4 * (4 * c 3 3) := by rw [fin5_extract_x_sum (fun i => c i 3)]
      _ = c 3 3 := h44 (c 3 3)
  rw [hPoly, hA23] at hEq
  exact hEq

def A7 (r c : Fin 7) : Fin 7 :=
  match r.val, c.val with
  | 0, 0 => 0 | 0, 1 => 1 | 0, 2 => 2 | 0, 3 => 3 | 0, 4 => 4 | 0, 5 => 5 | 0, 6 => 6
  | 1, 0 => 1 | 1, 1 => 0 | 1, 2 => 3 | 1, 3 => 2 | 1, 4 => 5 | 1, 5 => 6 | 1, 6 => 4
  | 2, 0 => 2 | 2, 1 => 3 | 2, 2 => 0 | 2, 3 => 1 | 2, 4 => 6 | 2, 5 => 4 | 2, 6 => 5
  | 3, 0 => 3 | 3, 1 => 4 | 3, 2 => 5 | 3, 3 => 6 | 3, 4 => 1 | 3, 5 => 2 | 3, 6 => 0
  | 4, 0 => 4 | 4, 1 => 2 | 4, 2 => 6 | 4, 3 => 5 | 4, 4 => 0 | 4, 5 => 3 | 4, 6 => 1
  | 5, 0 => 5 | 5, 1 => 6 | 5, 2 => 1 | 5, 3 => 4 | 5, 4 => 3 | 5, 5 => 0 | 5, 6 => 2
  | _, 0 => 6 | _, 1 => 5 | _, 2 => 4 | _, 3 => 0 | _, 4 => 2 | _, 5 => 1 | _, _ => 3

def B7 (r c : Fin 7) : Fin 7 :=
  match r.val, c.val with
  | 0, 0 => 0 | 0, 1 => 1 | 0, 2 => 2 | 0, 3 => 3 | 0, 4 => 4 | 0, 5 => 5 | 0, 6 => 6
  | 1, 0 => 2 | 1, 1 => 4 | 1, 2 => 1 | 1, 3 => 5 | 1, 4 => 6 | 1, 5 => 0 | 1, 6 => 3
  | 2, 0 => 3 | 2, 1 => 0 | 2, 2 => 5 | 2, 3 => 4 | 2, 4 => 2 | 2, 5 => 6 | 2, 6 => 1
  | 3, 0 => 6 | 3, 1 => 5 | 3, 2 => 3 | 3, 3 => 1 | 3, 4 => 0 | 3, 5 => 4 | 3, 6 => 2
  | 4, 0 => 1 | 4, 1 => 6 | 4, 2 => 4 | 4, 3 => 0 | 4, 4 => 3 | 4, 5 => 2 | 4, 6 => 5
  | 5, 0 => 4 | 5, 1 => 3 | 5, 2 => 6 | 5, 3 => 2 | 5, 4 => 5 | 5, 5 => 1 | 5, 6 => 0
  | _, 0 => 5 | _, 1 => 2 | _, 2 => 0 | _, 3 => 6 | _, 4 => 1 | _, 5 => 3 | _, _ => 4

private def rowInvA7 (r v : Fin 7) : Fin 7 :=
  if A7 r 0 == v then 0 else if A7 r 1 == v then 1 else if A7 r 2 == v then 2
  else if A7 r 3 == v then 3 else if A7 r 4 == v then 4 else if A7 r 5 == v then 5 else 6

private def colInvA7 (c v : Fin 7) : Fin 7 :=
  if A7 0 c == v then 0 else if A7 1 c == v then 1 else if A7 2 c == v then 2
  else if A7 3 c == v then 3 else if A7 4 c == v then 4 else if A7 5 c == v then 5 else 6

private def rowInvB7 (r v : Fin 7) : Fin 7 :=
  if B7 r 0 == v then 0 else if B7 r 1 == v then 1 else if B7 r 2 == v then 2
  else if B7 r 3 == v then 3 else if B7 r 4 == v then 4 else if B7 r 5 == v then 5 else 6

private def colInvB7 (c v : Fin 7) : Fin 7 :=
  if B7 0 c == v then 0 else if B7 1 c == v then 1 else if B7 2 c == v then 2
  else if B7 3 c == v then 3 else if B7 4 c == v then 4 else if B7 5 c == v then 5 else 6

private def invR7 (u v : Fin 7) : Fin 7 :=
  match u.val, v.val with
  | 0, 0 => 0 | 0, 1 => 5 | 0, 2 => 3 | 0, 3 => 4 | 0, 4 => 1 | 0, 5 => 2 | 0, 6 => 6
  | 1, 0 => 3 | 1, 1 => 0 | 1, 2 => 1 | 1, 3 => 6 | 1, 4 => 2 | 1, 5 => 4 | 1, 6 => 5
  | 2, 0 => 5 | 2, 1 => 6 | 2, 2 => 0 | 2, 3 => 2 | 2, 4 => 3 | 2, 5 => 1 | 2, 6 => 4
  | 3, 0 => 2 | 3, 1 => 1 | 3, 2 => 4 | 3, 3 => 0 | 3, 4 => 6 | 3, 5 => 5 | 3, 6 => 3
  | 4, 0 => 6 | 4, 1 => 4 | 4, 2 => 5 | 4, 3 => 1 | 4, 4 => 0 | 4, 5 => 3 | 4, 6 => 2
  | 5, 0 => 4 | 5, 1 => 2 | 5, 2 => 6 | 5, 3 => 3 | 5, 4 => 5 | 5, 5 => 0 | 5, 6 => 1
  | _, 0 => 1 | _, 1 => 3 | _, 2 => 2 | _, 3 => 5 | _, 4 => 4 | _, 5 => 6 | _, _ => 0

private def invC7 (u v : Fin 7) : Fin 7 :=
  match u.val, v.val with
  | 0, 0 => 0 | 0, 1 => 5 | 0, 2 => 6 | 0, 3 => 4 | 0, 4 => 1 | 0, 5 => 2 | 0, 6 => 3
  | 1, 0 => 4 | 1, 1 => 1 | 1, 2 => 0 | 1, 3 => 5 | 1, 4 => 3 | 1, 5 => 6 | 1, 6 => 2
  | 2, 0 => 6 | 2, 1 => 4 | 2, 2 => 2 | 2, 3 => 0 | 2, 4 => 5 | 2, 5 => 3 | 2, 6 => 1
  | 3, 0 => 1 | 3, 1 => 2 | 3, 2 => 5 | 3, 3 => 3 | 3, 4 => 6 | 3, 5 => 4 | 3, 6 => 0
  | 4, 0 => 2 | 4, 1 => 0 | 4, 2 => 3 | 4, 3 => 6 | 4, 4 => 4 | 4, 5 => 1 | 4, 6 => 5
  | 5, 0 => 3 | 5, 1 => 6 | 5, 2 => 1 | 5, 3 => 2 | 5, 4 => 0 | 5, 5 => 5 | 5, 6 => 4
  | _, 0 => 5 | _, 1 => 3 | _, 2 => 4 | _, 3 => 1 | _, 4 => 2 | _, 5 => 0 | _, _ => 6

/-- Certification that the `d* = 7` witness `(A7, B7)` is a valid 7×7 orthogonal Latin square pair (`IsOLS7`). -/
theorem ame47_witness_is_ols7 : IsOLS7 A7 B7 := by
  refine ⟨isLatin_of_inverses A7 rowInvA7 colInvA7
            (by decide) (by decide) (by decide) (by decide),
          isLatin_of_inverses B7 rowInvB7 colInvB7
            (by decide) (by decide) (by decide) (by decide),
          isOrthogonal_of_inverse A7 B7 invR7 invC7
            (by decide) (by decide) (by decide) (by decide)⟩

def IsIsotopicToZ7 (A : Fin 7 → Fin 7 → Fin 7) : Prop :=
  ∃ α β γ : Fin 7 → Fin 7, Bijective α ∧ Bijective β ∧ Bijective γ ∧
    ∀ x y : Fin 7, γ (A (α x) (β y)) = x + y

def SatisfiesQuadrangle7 (A : Fin 7 → Fin 7 → Fin 7) : Prop :=
  ∀ r1 r2 c1 c2 r1' c1' r2' c2' : Fin 7,
    A r1 c1 = A r1' c1' →
    A r1 c2 = A r1' c2' →
    A r2 c1 = A r2' c1' →
    A r2 c2 = A r2' c2'

private theorem fin7_add_comm (a b : Fin 7) : a + b = b + a := by
  revert a b; decide

private theorem fin7_add_assoc (a b c : Fin 7) : a + b + c = a + (b + c) := by
  revert a b c; decide

private theorem fin7_add_sub_cancel (a b : Fin 7) : a + b - a = b := by
  revert a b; decide

private theorem fin7_quad_sum (a b c d : Fin 7) :
    (a + b) + (c + d) - (a + d) = c + b := by
  rw [fin7_add_assoc, fin7_add_comm b (c + d), fin7_add_assoc c d b, fin7_add_comm d b]
  rw [← fin7_add_assoc c b d, fin7_add_comm (c + b) d, ← fin7_add_assoc a d (c + b)]
  exact fin7_add_sub_cancel (a + d) (c + b)

private theorem fin7_quad_alg (a1 a2 a3 a4 a1' a2' a3' a4' : Fin 7)
    (h1 : a1 = a1') (h2 : a2 = a2') (h3 : a3 = a3')
    (x1 x2 y1 y2 x1' x2' y1' y2' : Fin 7)
    (ha1 : a1 = x1 + y1) (ha2 : a2 = x1 + y2) (ha3 : a3 = x2 + y1) (ha4 : a4 = x2 + y2)
    (hb1 : a1' = x1' + y1') (hb2 : a2' = x1' + y2') (hb3 : a3' = x2' + y1') (hb4 : a4' = x2' + y2') :
    a4 = a4' := by
  have hL : a2 + a3 - a1 = a4 := by
    rw [ha1, ha2, ha3, ha4]
    exact fin7_quad_sum x1 y2 x2 y1
  have hR : a2' + a3' - a1' = a4' := by
    rw [hb1, hb2, hb3, hb4]
    exact fin7_quad_sum x1' y2' x2' y1'
  rw [← hL, ← hR, h1, h2, h3]

/-- Every 7×7 Latin square isotopic to `Z_7` satisfies the Reidemeister quadrangle criterion. -/
theorem isotopic_to_Z7_implies_quadrangle (A : Fin 7 → Fin 7 → Fin 7)
    (hIso : IsIsotopicToZ7 A) : SatisfiesQuadrangle7 A := by
  rcases hIso with ⟨α, β, γ, hα, hβ, hγ, hEq⟩
  let αInv : Fin 7 → Fin 7 := fun r => Classical.choose (hα.2 r)
  let βInv : Fin 7 → Fin 7 := fun c => Classical.choose (hβ.2 c)
  have hαInv : ∀ r, α (αInv r) = r := fun r => Classical.choose_spec (hα.2 r)
  have hβInv : ∀ c, β (βInv c) = c := fun c => Classical.choose_spec (hβ.2 c)
  have hLin : ∀ r c : Fin 7, γ (A r c) = αInv r + βInv c := by
    intro r c
    have h := hEq (αInv r) (βInv c)
    rwa [hαInv r, hβInv c] at h
  intro r1 r2 c1 c2 r1' c1' r2' c2' h1 h2 h3
  have hg1 : γ (A r1 c1) = γ (A r1' c1') := congrArg γ h1
  have hg2 : γ (A r1 c2) = γ (A r1' c2') := congrArg γ h2
  have hg3 : γ (A r2 c1) = γ (A r2' c1') := congrArg γ h3
  have hg4 : γ (A r2 c2) = γ (A r2' c2') :=
    fin7_quad_alg _ _ _ _ _ _ _ _ hg1 hg2 hg3
      (αInv r1) (αInv r2) (βInv c1) (βInv c2)
      (αInv r1') (αInv r2') (βInv c1') (βInv c2')
      (hLin r1 c1) (hLin r1 c2) (hLin r2 c1) (hLin r2 c2)
      (hLin r1' c1') (hLin r1' c2') (hLin r2' c1') (hLin r2' c2')
  exact hγ.1 _ _ hg4

/-- `A7` is strictly NOT isotopic to the cyclic group `Z_7` (explicit quadrangle counterexample). -/
theorem ame47_A7_not_isotopic_to_Z7 : ¬ IsIsotopicToZ7 A7 := by
  intro hIso
  have hQuad := isotopic_to_Z7_implies_quadrangle A7 hIso
  have h22 := hQuad 0 1 0 1 2 2 5 3 (by decide) (by decide) (by decide)
  revert h22
  decide

/-! ## Section 3: `AME(4,6)` Exact Eisenstein `Z[omega]` Quantum Escape Below Defect `2` -/

/-- Exact Eisenstein integer `re + om * ω`, where `ω^2 + ω + 1 = 0`. -/
structure EisenInt where
  re : Int
  om : Int
  deriving DecidableEq, Repr

def eZero : EisenInt := ⟨0, 0⟩
def eOne : EisenInt := ⟨1, 0⟩
def eOmega : EisenInt := ⟨0, 1⟩
def eAdd (x y : EisenInt) : EisenInt := ⟨x.re + y.re, x.om + y.om⟩
def eSub (x y : EisenInt) : EisenInt := ⟨x.re - y.re, x.om - y.om⟩
def eMul (x y : EisenInt) : EisenInt :=
  ⟨x.re * y.re - x.om * y.om, x.re * y.om + x.om * y.re - x.om * y.om⟩
def eConj (x : EisenInt) : EisenInt := ⟨x.re - x.om, -x.om⟩
def eNorm (x : EisenInt) : Int := x.re * x.re - x.re * x.om + x.om * x.om

/-- Fundamental Eisenstein cyclotomic identity `ω^2 + ω + 1 = 0`. -/
theorem eisen_omega_eq : eAdd (eAdd (eMul eOmega eOmega) eOmega) eOne = eZero := by
  decide

/-- Exact conjugation norm identity `z^* z = N(z) ∈ Z`. -/
theorem eisen_conj_mul_self (x : EisenInt) : eMul (eConj x) x = ⟨eNorm x, 0⟩ := by
  rcases x with ⟨a, b⟩
  simp only [eMul, eConj, eNorm]
  congr 1
  · rw [Int.sub_mul, Int.neg_mul, Int.sub_neg, Int.mul_comm b a]
  · rw [Int.sub_mul, Int.neg_mul, Int.neg_mul, Int.sub_neg, Int.mul_comm b a]
    omega

private theorem int_sq_eisen_id (a b c d : Int) :
    (a * c - b * d) * (a * c - b * d) - (a * c - b * d) * (a * d + b * c - b * d) +
      (a * d + b * c - b * d) * (a * d + b * c - b * d) =
      (a * a - a * b + b * b) * (c * c - c * d + d * d) := by
  have h1 : (a * c - b * d) * (a * c - b * d) =
      a * a * (c * c) - 2 * (a * b * (c * d)) + b * b * (d * d) := by
    rw [Int.sub_mul, Int.mul_sub, Int.mul_sub]
    have ha : a * c * (a * c) = a * a * (c * c) := by
      rw [Int.mul_assoc, Int.mul_left_comm c a c, ← Int.mul_assoc]
    have hb : b * d * (b * d) = b * b * (d * d) := by
      rw [Int.mul_assoc, Int.mul_left_comm d b d, ← Int.mul_assoc]
    have hc : a * c * (b * d) = a * b * (c * d) := by
      rw [Int.mul_assoc, Int.mul_left_comm c b d, ← Int.mul_assoc]
    have hd : b * d * (a * c) = a * b * (c * d) := by
      rw [Int.mul_comm (b * d) (a * c), hc]
    rw [ha, hb, hc, hd]
    omega
  have h2 : (a * c - b * d) * (a * d + b * c - b * d) =
      a * a * (c * d) + a * b * (c * c) - a * b * (c * d) -
      a * b * (d * d) - b * b * (c * d) + b * b * (d * d) := by
    rw [Int.sub_mul, Int.mul_sub, Int.mul_sub, Int.mul_add, Int.mul_add]
    have p1 : a * c * (a * d) = a * a * (c * d) := by
      rw [Int.mul_assoc, Int.mul_left_comm c a d, ← Int.mul_assoc]
    have p2 : a * c * (b * c) = a * b * (c * c) := by
      rw [Int.mul_assoc, Int.mul_left_comm c b c, ← Int.mul_assoc]
    have p3 : a * c * (b * d) = a * b * (c * d) := by
      rw [Int.mul_assoc, Int.mul_left_comm c b d, ← Int.mul_assoc]
    have p4 : b * d * (a * d) = a * b * (d * d) := by
      rw [Int.mul_assoc, Int.mul_left_comm d a d, ← Int.mul_assoc, Int.mul_comm b a]
    have p5 : b * d * (b * c) = b * b * (c * d) := by
      rw [Int.mul_assoc, Int.mul_left_comm d b c, Int.mul_comm d c, ← Int.mul_assoc]
    have p6 : b * d * (b * d) = b * b * (d * d) := by
      rw [Int.mul_assoc, Int.mul_left_comm d b d, ← Int.mul_assoc]
    rw [p1, p2, p3, p4, p5, p6]
    omega
  have h3 : (a * d + b * c - b * d) * (a * d + b * c - b * d) =
      a * a * (d * d) + b * b * (c * c) + b * b * (d * d) +
      2 * (a * b * (c * d)) - 2 * (a * b * (d * d)) - 2 * (b * b * (c * d)) := by
    simp only [Int.sub_mul, Int.mul_sub, Int.add_mul, Int.mul_add]
    have q1 : a * d * (a * d) = a * a * (d * d) := by
      rw [Int.mul_assoc, Int.mul_left_comm d a d, ← Int.mul_assoc]
    have q2 : a * d * (b * c) = a * b * (c * d) := by
      rw [Int.mul_assoc, Int.mul_left_comm d b c, Int.mul_comm d c, ← Int.mul_assoc]
    have q3 : b * c * (a * d) = a * b * (c * d) := by
      rw [Int.mul_comm (b * c) (a * d), q2]
    have q4 : b * c * (b * c) = b * b * (c * c) := by
      rw [Int.mul_assoc, Int.mul_left_comm c b c, ← Int.mul_assoc]
    have q6 : b * d * (a * d) = a * b * (d * d) := by
      rw [Int.mul_assoc, Int.mul_left_comm d a d, ← Int.mul_assoc, Int.mul_comm b a]
    have q7 : b * d * (b * c) = b * b * (c * d) := by
      rw [Int.mul_assoc, Int.mul_left_comm d b c, Int.mul_comm d c, ← Int.mul_assoc]
    have q8 : b * d * (b * d) = b * b * (d * d) := by
      rw [Int.mul_assoc, Int.mul_left_comm d b d, ← Int.mul_assoc]
    have q9 : (a * d) * (b * d) = a * b * (d * d) := by
      rw [Int.mul_assoc, Int.mul_left_comm d b d, ← Int.mul_assoc]
    have q10 : (b * c) * (b * d) = b * b * (c * d) := by
      rw [Int.mul_assoc, Int.mul_left_comm c b d, ← Int.mul_assoc]
    rw [q1, q2, q3, q4, q6, q7, q8, q9, q10]
    omega
  have h4 : (a * a - a * b + b * b) * (c * c - c * d + d * d) =
      a * a * (c * c) - a * a * (c * d) + a * a * (d * d) -
      (a * b * (c * c) - a * b * (c * d) + a * b * (d * d)) +
      (b * b * (c * c) - b * b * (c * d) + b * b * (d * d)) := by
    rw [Int.add_mul, Int.sub_mul, Int.mul_add, Int.mul_sub, Int.mul_add, Int.mul_sub, Int.mul_add, Int.mul_sub]
  rw [h1, h2, h3, h4]
  omega

/-- Multiplicativity of the Eisenstein norm `N(xy) = N(x)N(y)`. -/
theorem eisen_norm_mul (x y : EisenInt) : eNorm (eMul x y) = eNorm x * eNorm y := by
  rcases x with ⟨a, b⟩
  rcases y with ⟨c, d⟩
  simp only [eNorm, eMul]
  exact int_sq_eisen_id a b c d

/-- Exact quadratic extension `Z[ω][√13]` represented as `u + v * √13`. -/
structure EisenSqrt13 where
  u : EisenInt
  v : EisenInt
  deriving DecidableEq, Repr

def qZero : EisenSqrt13 := ⟨eZero, eZero⟩
def qAdd (x y : EisenSqrt13) : EisenSqrt13 := ⟨eAdd x.u y.u, eAdd x.v y.v⟩
def qSub (x y : EisenSqrt13) : EisenSqrt13 := ⟨eSub x.u y.u, eSub x.v y.v⟩
def qScale13 (x : EisenInt) : EisenInt := ⟨13 * x.re, 13 * x.om⟩
def qMul (x y : EisenSqrt13) : EisenSqrt13 :=
  ⟨eAdd (eMul x.u y.u) (qScale13 (eMul x.v y.v)), eAdd (eMul x.u y.v) (eMul x.v y.u)⟩
def qConj (x : EisenSqrt13) : EisenSqrt13 := ⟨eConj x.u, eConj x.v⟩

private theorem eAdd_zeroL (x : EisenInt) : eAdd eZero x = x := by
  rcases x with ⟨a, b⟩; simp [eAdd, eZero]

private theorem eAdd_zeroR (x : EisenInt) : eAdd x eZero = x := by
  rcases x with ⟨a, b⟩; simp [eAdd, eZero]

private theorem eAdd_assoc (x y z : EisenInt) : eAdd (eAdd x y) z = eAdd x (eAdd y z) := by
  rcases x with ⟨a, b⟩; rcases y with ⟨c, d⟩; rcases z with ⟨e, f⟩
  simp [eAdd, Int.add_assoc]

private theorem eAdd_comm (x y : EisenInt) : eAdd x y = eAdd y x := by
  rcases x with ⟨a, b⟩; rcases y with ⟨c, d⟩
  simp [eAdd, Int.add_comm]

private theorem eConj_zero : eConj eZero = eZero := rfl

private theorem eMul_zeroL (x : EisenInt) : eMul eZero x = eZero := by
  rcases x with ⟨a, b⟩; simp [eMul, eZero]

private theorem eMul_zeroR (x : EisenInt) : eMul x eZero = eZero := by
  rcases x with ⟨a, b⟩; simp [eMul, eZero]

private theorem qAdd_zeroL (x : EisenSqrt13) : qAdd qZero x = x := by
  rcases x with ⟨u, v⟩; simp [qAdd, qZero, eAdd_zeroL]

private theorem qAdd_zeroR (x : EisenSqrt13) : qAdd x qZero = x := by
  rcases x with ⟨u, v⟩; simp [qAdd, qZero, eAdd_zeroR]

private theorem qAdd_assoc (x y z : EisenSqrt13) : qAdd (qAdd x y) z = qAdd x (qAdd y z) := by
  rcases x with ⟨u1, v1⟩; rcases y with ⟨u2, v2⟩; rcases z with ⟨u3, v3⟩
  simp [qAdd, eAdd_assoc]

private theorem qAdd_comm (x y : EisenSqrt13) : qAdd x y = qAdd y x := by
  rcases x with ⟨u1, v1⟩; rcases y with ⟨u2, v2⟩
  simp [qAdd, eAdd_comm]

private theorem qConj_zero : qConj qZero = qZero := rfl

private theorem qMul_zeroL (x : EisenSqrt13) : qMul qZero x = qZero := by
  rcases x with ⟨u, v⟩
  simp only [qMul, qZero]
  have hu : eMul eZero u = eZero := eMul_zeroL u
  have hv : eMul eZero v = eZero := eMul_zeroL v
  rw [hu, hv]
  rfl

private theorem qMul_zeroR (x : EisenSqrt13) : qMul x qZero = qZero := by
  rcases x with ⟨u, v⟩
  simp only [qMul, qZero]
  have hu : eMul u eZero = eZero := eMul_zeroR u
  have hv : eMul v eZero = eZero := eMul_zeroR v
  rw [hu, hv]
  rfl

def allFin6 (p : Fin 6 → Bool) : Bool :=
  p 0 && p 1 && p 2 && p 3 && p 4 && p 5

def allFin36 (p : Fin 36 → Bool) : Bool :=
  allFin6 fun q => allFin6 fun r =>
    p ⟨q.val * 6 + r.val, by have := q.isLt; have := r.isLt; omega⟩

private theorem allFin6_forall {p : Fin 6 → Bool} (h : allFin6 p = true) :
    ∀ x : Fin 6, p x = true := by
  unfold allFin6 at h
  simp only [Bool.and_eq_true] at h
  rcases h with ⟨⟨⟨⟨⟨h0, h1⟩, h2⟩, h3⟩, h4⟩, h5⟩
  intro x
  rcases x with ⟨v, hv⟩
  match v, hv with
  | 0, _ => exact h0
  | 1, _ => exact h1
  | 2, _ => exact h2
  | 3, _ => exact h3
  | 4, _ => exact h4
  | 5, _ => exact h5

private theorem allFin36_forall {p : Fin 36 → Bool} (h : allFin36 p = true) :
    ∀ x : Fin 36, p x = true := by
  intro x
  have hq : allFin6 (fun r => p ⟨(x.val / 6) * 6 + r.val, by have := x.isLt; omega⟩) = true :=
    allFin6_forall h ⟨x.val / 6, by have := x.isLt; omega⟩
  have hr := allFin6_forall hq ⟨x.val % 6, Nat.mod_lt _ (by omega)⟩
  have heq : (⟨(x.val / 6) * 6 + x.val % 6, by have := x.isLt; omega⟩ : Fin 36) = x := by
    apply Fin.ext; dsimp; omega
  rwa [heq] at hr

def sumFin6 {α : Type} (add : α → α → α) (f : Fin 6 → α) : α :=
  add (add (add (f 0) (f 1)) (add (f 2) (f 3))) (add (f 4) (f 5))

def sumFin36 {α : Type} (add : α → α → α) (f : Fin 36 → α) : α :=
  sumFin6 add fun q => sumFin6 add fun r =>
    f ⟨q.val * 6 + r.val, by have := q.isLt; have := r.isLt; omega⟩

private theorem sumFin6_zero {α : Type} (zero : α) (add : α → α → α)
    (hadd0L : ∀ x, add zero x = x)
    (f : Fin 6 → α) (hz : ∀ r : Fin 6, f r = zero) :
    sumFin6 add f = zero := by
  simp only [sumFin6, hz, hadd0L]

private theorem fin6_cases (x : Fin 6) : x = 0 ∨ x = 1 ∨ x = 2 ∨ x = 3 ∨ x = 4 ∨ x = 5 := by
  rcases x with ⟨v, hv⟩
  match v, hv with
  | 0, _ => exact Or.inl rfl
  | 1, _ => exact Or.inr (Or.inl rfl)
  | 2, _ => exact Or.inr (Or.inr (Or.inl rfl))
  | 3, _ => exact Or.inr (Or.inr (Or.inr (Or.inl rfl)))
  | 4, _ => exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inl rfl))))
  | 5, _ => exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inr rfl))))

private theorem sumFin6_single {α : Type} (zero : α) (add : α → α → α)
    (hadd0L : ∀ x, add zero x = x) (hadd0R : ∀ x, add x zero = x)
    (f : Fin 6 → α) (r0 : Fin 6) (hz : ∀ r : Fin 6, r ≠ r0 → f r = zero) :
    sumFin6 add f = f r0 := by
  rcases fin6_cases r0 with rfl | rfl | rfl | rfl | rfl | rfl
  · simp [sumFin6, hz 1 (by decide), hz 2 (by decide), hz 3 (by decide), hz 4 (by decide), hz 5 (by decide), hadd0L, hadd0R]
  · simp [sumFin6, hz 0 (by decide), hz 2 (by decide), hz 3 (by decide), hz 4 (by decide), hz 5 (by decide), hadd0L, hadd0R]
  · simp [sumFin6, hz 0 (by decide), hz 1 (by decide), hz 3 (by decide), hz 4 (by decide), hz 5 (by decide), hadd0L, hadd0R]
  · simp [sumFin6, hz 0 (by decide), hz 1 (by decide), hz 2 (by decide), hz 4 (by decide), hz 5 (by decide), hadd0L, hadd0R]
  · simp [sumFin6, hz 0 (by decide), hz 1 (by decide), hz 2 (by decide), hz 3 (by decide), hz 5 (by decide), hadd0L, hadd0R]
  · simp [sumFin6, hz 0 (by decide), hz 1 (by decide), hz 2 (by decide), hz 3 (by decide), hz 4 (by decide), hadd0L]

private theorem sumFin6_add {α : Type} (add : α → α → α)
    (hassoc : ∀ x y z, add (add x y) z = add x (add y z))
    (hcomm : ∀ x y, add x y = add y x)
    (f g : Fin 6 → α) :
    sumFin6 add (fun r => add (f r) (g r)) = add (sumFin6 add f) (sumFin6 add g) := by
  have h4 : ∀ a b c d : α, add (add a b) (add c d) = add (add a c) (add b d) := by
    intro a b c d
    rw [hassoc, ← hassoc b c d, hcomm b c, hassoc c b d, ← hassoc]
  simp only [sumFin6, h4]

private theorem sumFin36_single {α : Type} (zero : α) (add : α → α → α)
    (hadd0L : ∀ x, add zero x = x) (hadd0R : ∀ x, add x zero = x)
    (f : Fin 36 → α) (a : Fin 36) (hz : ∀ k : Fin 36, k ≠ a → f k = zero) :
    sumFin36 add f = f a := by
  let q0 : Fin 6 := ⟨a.val / 6, by have := a.isLt; omega⟩
  let r0 : Fin 6 := ⟨a.val % 6, Nat.mod_lt _ (by omega)⟩
  unfold sumFin36
  rw [sumFin6_single zero add hadd0L hadd0R _ q0]
  · rw [sumFin6_single zero add hadd0L hadd0R _ r0]
    · have ha : (⟨q0.val * 6 + r0.val, by have := q0.isLt; have := r0.isLt; omega⟩ : Fin 36) = a := by
        apply Fin.ext; dsimp [q0, r0]; omega
      rw [ha]
    · intro r hr
      apply hz
      intro heq
      apply hr
      apply Fin.ext
      have := congrArg Fin.val heq
      dsimp [q0, r0] at *
      omega
  · intro q hq
    apply sumFin6_zero zero add hadd0L
    intro r
    apply hz
    intro heq
    apply hq
    apply Fin.ext
    have := congrArg Fin.val heq
    dsimp [q0] at *
    have := r.isLt
    omega

private theorem sumFin36_add {α : Type} (add : α → α → α)
    (hassoc : ∀ x y z, add (add x y) z = add x (add y z))
    (hcomm : ∀ x y, add x y = add y x)
    (f g : Fin 36 → α) :
    sumFin36 add (fun k => add (f k) (g k)) = add (sumFin36 add f) (sumFin36 add g) := by
  unfold sumFin36
  have h1 : ∀ q : Fin 6,
      sumFin6 add (fun r => add (f ⟨q.val * 6 + r.val, by have := q.isLt; have := r.isLt; omega⟩)
                                (g ⟨q.val * 6 + r.val, by have := q.isLt; have := r.isLt; omega⟩)) =
      add (sumFin6 add (fun r => f ⟨q.val * 6 + r.val, by have := q.isLt; have := r.isLt; omega⟩))
          (sumFin6 add (fun r => g ⟨q.val * 6 + r.val, by have := q.isLt; have := r.isLt; omega⟩)) :=
    fun q => sumFin6_add add hassoc hcomm _ _
  simp only [h1]
  exact sumFin6_add add hassoc hcomm _ _

private theorem sumFin36_eq_two_terms {α : Type} (zero : α) (add : α → α → α)
    (hadd0L : ∀ x, add zero x = x) (hadd0R : ∀ x, add x zero = x)
    (hassoc : ∀ x y z, add (add x y) z = add x (add y z))
    (hcomm : ∀ x y, add x y = add y x)
    (f : Fin 36 → α) (a b : Fin 36)
    (hsupp : ∀ k : Fin 36, k ≠ a → k ≠ b → f k = zero) :
    sumFin36 add f = if a == b then f a else add (f a) (f b) := by
  by_cases hab : a = b
  · subst hab
    simp only [beq_self_eq_true, ↓reduceIte]
    exact sumFin36_single zero add hadd0L hadd0R f a (fun k hk => hsupp k hk hk)
  · have hneB : (a == b) = false := beq_false_of_ne hab
    simp only [hneB, Bool.false_eq_true, ↓reduceIte]
    let fa : Fin 36 → α := fun k => if k == a then f k else zero
    let fb : Fin 36 → α := fun k => if k == b then f k else zero
    have hsplit : f = fun k => add (fa k) (fb k) := by
      funext k
      by_cases hka : k = a
      · subst hka
        have hkb : (k == b) = false := beq_false_of_ne hab
        dsimp [fa, fb]
        rw [beq_self_eq_true, hkb]
        exact (hadd0R (f k)).symm
      · by_cases hkb : k = b
        · subst hkb
          have hkaB : (k == a) = false := beq_false_of_ne hka
          dsimp [fa, fb]
          rw [hkaB, beq_self_eq_true]
          exact (hadd0L (f k)).symm
        · have hkaB : (k == a) = false := beq_false_of_ne hka
          have hkbB : (k == b) = false := beq_false_of_ne hkb
          dsimp [fa, fb]
          rw [hkaB, hkbB, hsupp k hka hkb]
          exact (hadd0L zero).symm
    have hsum : sumFin36 add f = add (sumFin36 add fa) (sumFin36 add fb) := by
      conv => lhs; rw [hsplit]
      exact sumFin36_add add hassoc hcomm fa fb
    have ha : sumFin36 add fa = f a := by
      rw [sumFin36_single zero add hadd0L hadd0R fa a]
      · simp [fa]
      · intro k hk
        dsimp [fa]
        rw [beq_false_of_ne hk]
        rfl
    have hb : sumFin36 add fb = f b := by
      rw [sumFin36_single zero add hadd0L hadd0R fb b]
      · simp [fb]
      · intro k hk
        dsimp [fb]
        rw [beq_false_of_ne hk]
        rfl
    rw [hsum, ha, hb]

def matDagMat {α : Type} (add mul : α → α → α) (conj : α → α)
    (M : Fin 36 → Fin 36 → α) (i j : Fin 36) : α :=
  sumFin36 add (fun k => mul (conj (M k i)) (M k j))

def sparseDotFast {α : Type} (zero : α) (add mul : α → α → α) (conj : α → α)
    (M : Fin 36 → Fin 36 → α) (r1 r2 : Fin 36 → Fin 36) (i j : Fin 36) : α :=
  if !(r1 i == r1 j || r1 i == r2 j || r2 i == r1 j || r2 i == r2 j) then zero
  else if r1 i == r2 i then mul (conj (M (r1 i) i)) (M (r1 i) j)
  else add (mul (conj (M (r1 i) i)) (M (r1 i) j)) (mul (conj (M (r2 i) i)) (M (r2 i) j))

private theorem matDagMat_eq_sparseDotFast {α : Type} (zero : α) (add mul : α → α → α) (conj : α → α)
    (hadd0L : ∀ x, add zero x = x) (hadd0R : ∀ x, add x zero = x)
    (hassoc : ∀ x y z, add (add x y) z = add x (add y z))
    (hcomm : ∀ x y, add x y = add y x)
    (hconj0 : conj zero = zero)
    (hmul0L : ∀ x, mul zero x = zero) (hmul0R : ∀ x, mul x zero = zero)
    (M : Fin 36 → Fin 36 → α) (r1 r2 : Fin 36 → Fin 36)
    (hsupp : ∀ k j : Fin 36, k ≠ r1 j → k ≠ r2 j → M k j = zero)
    (i j : Fin 36) :
    matDagMat add mul conj M i j = sparseDotFast zero add mul conj M r1 r2 i j := by
  have hf : ∀ k : Fin 36, k ≠ r1 i → k ≠ r2 i → mul (conj (M k i)) (M k j) = zero := by
    intro k hk1 hk2
    rw [hsupp k i hk1 hk2, hconj0, hmul0L]
  have h1 := sumFin36_eq_two_terms zero add hadd0L hadd0R hassoc hcomm
    (fun k => mul (conj (M k i)) (M k j)) (r1 i) (r2 i) hf
  change sumFin36 add (fun k => mul (conj (M k i)) (M k j)) =
    sparseDotFast zero add mul conj M r1 r2 i j at *
  rw [h1]
  unfold sparseDotFast
  by_cases hov : (r1 i == r1 j || r1 i == r2 j || r2 i == r1 j || r2 i == r2 j) = true
  · simp [hov]
  · have hne11 : r1 i ≠ r1 j := fun h => hov (by simp [h])
    have hne12 : r1 i ≠ r2 j := fun h => hov (by simp [h])
    have hne21 : r2 i ≠ r1 j := fun h => hov (by simp [h])
    have hne22 : r2 i ≠ r2 j := fun h => hov (by simp [h])
    have hz1 : M (r1 i) j = zero := hsupp (r1 i) j hne11 hne12
    have hz2 : M (r2 i) j = zero := hsupp (r2 i) j hne21 hne22
    simp [hov, hz1, hz2, hmul0R, hadd0R]

def reshuffleGamma {α : Type} (M : Fin 36 → Fin 36 → α) (r c : Fin 36) : α :=
  M ⟨(r.val / 6) * 6 + c.val % 6, by have := r.isLt; have := c.isLt; omega⟩
    ⟨(c.val / 6) * 6 + r.val % 6, by have := r.isLt; have := c.isLt; omega⟩

def reshuffleR {α : Type} (M : Fin 36 → Fin 36 → α) (r c : Fin 36) : α :=
  M ⟨(r.val % 6) * 6 + c.val % 6, by have := r.isLt; have := c.isLt; omega⟩
    ⟨(r.val / 6) * 6 + c.val / 6, by have := r.isLt; have := c.isLt; omega⟩

/-- Classical 6×6 witness `(A0_6, B0_6)` achieving the classical 2-collision defect floor. -/
def perm0_6 (r : Fin 36) : Fin 36 :=
  match r.val with
  | 0 => 0   | 1 => 7   | 2 => 14  | 3 => 21  | 4 => 28  | 5 => 35
  | 6 => 9   | 7 => 2   | 8 => 23  | 9 => 16  | 10 => 30 | 11 => 25
  | 12 => 29 | 13 => 34 | 14 => 1  | 15 => 6  | 16 => 15 | 17 => 20
  | 18 => 32 | 19 => 27 | 20 => 10 | 21 => 5  | 22 => 19 | 23 => 12
  | 24 => 13 | 25 => 18 | 26 => 33 | 27 => 26 | 28 => 4  | 29 => 11
  | 30 => 22 | 31 => 17 | 32 => 24 | 33 => 31 | 34 => 8  | _ => 3

def perm0Inv_6 (c : Fin 36) : Fin 36 :=
  match c.val with
  | 0 => 0   | 1 => 14  | 2 => 7   | 3 => 35  | 4 => 28  | 5 => 21
  | 6 => 15  | 7 => 1   | 8 => 34  | 9 => 6   | 10 => 20 | 11 => 29
  | 12 => 23 | 13 => 24 | 14 => 2  | 15 => 16 | 16 => 9  | 17 => 31
  | 18 => 25 | 19 => 22 | 20 => 17 | 21 => 3  | 22 => 30 | 23 => 8
  | 24 => 32 | 25 => 11 | 26 => 27 | 27 => 19 | 28 => 4  | 29 => 12
  | 30 => 10 | 31 => 33 | 32 => 18 | 33 => 26 | 34 => 13 | _ => 5

def U0_6 (r c : Fin 36) : EisenInt :=
  if perm0_6 r == c then eOne else eZero

/-- Scaled 40-nonzero `(1/2) Z[ω]` 2-block gate `U40_2` (entries in `Z[ω]`, `U = (1/2) U40_2`). -/
def U40_2 (r c : Fin 36) : EisenInt :=
  match r.val, c.val with
  | 4, 28 => ⟨1, -1⟩  | 4, 35 => ⟨-1, 0⟩
  | 5, 28 => ⟨1, 0⟩   | 5, 35 => ⟨2, 1⟩
  | 28, 4 => ⟨1, -1⟩  | 28, 11 => ⟨1, 0⟩
  | 29, 4 => ⟨-1, 0⟩  | 29, 11 => ⟨2, 1⟩
  | _, _ =>
    if (r.val == 4) || (r.val == 5) || (r.val == 28) || (r.val == 29) then eZero
    else if perm0_6 r == c then ⟨2, 0⟩ else eZero

/-- Scaled 44-nonzero `Z[ω]` 4-block gate `U44_2s13` (entries in `Z[ω][√13]`, `U = (1/(2√13)) U44_2s13`). -/
def U44_2s13 (r c : Fin 36) : EisenSqrt13 :=
  match r.val, c.val with
  | 4, 28  => ⟨eZero, ⟨1, -1⟩⟩  | 4, 35  => ⟨eZero, ⟨-1, 0⟩⟩
  | 5, 28  => ⟨eZero, ⟨1, 0⟩⟩   | 5, 35  => ⟨eZero, ⟨2, 1⟩⟩
  | 10, 30 => ⟨⟨4, -4⟩, eZero⟩  | 10, 25 => ⟨⟨-2, 0⟩, eZero⟩
  | 11, 30 => ⟨⟨2, 0⟩, eZero⟩   | 11, 25 => ⟨⟨8, 4⟩, eZero⟩
  | 28, 4  => ⟨eZero, ⟨1, -1⟩⟩  | 28, 11 => ⟨eZero, ⟨1, 0⟩⟩
  | 29, 4  => ⟨eZero, ⟨-1, 0⟩⟩  | 29, 11 => ⟨eZero, ⟨2, 1⟩⟩
  | 34, 8  => ⟨⟨4, -4⟩, eZero⟩  | 34, 3  => ⟨⟨2, 0⟩, eZero⟩
  | 35, 8  => ⟨⟨-2, 0⟩, eZero⟩  | 35, 3  => ⟨⟨8, 4⟩, eZero⟩
  | _, _ =>
    if (r.val == 4) || (r.val == 5) || (r.val == 10) || (r.val == 11) ||
       (r.val == 28) || (r.val == 29) || (r.val == 34) || (r.val == 35) then qZero
    else if perm0_6 r == c then ⟨eZero, ⟨2, 0⟩⟩ else qZero

private def suppU_r1 (j : Fin 36) : Fin 36 :=
  match j.val with
  | 0 => 0   | 1 => 14  | 2 => 7   | 3 => 34  | 4 => 28  | 5 => 21
  | 6 => 15  | 7 => 1   | 8 => 34  | 9 => 6   | 10 => 20 | 11 => 28
  | 12 => 23 | 13 => 24 | 14 => 2  | 15 => 16 | 16 => 9  | 17 => 31
  | 18 => 25 | 19 => 22 | 20 => 17 | 21 => 3  | 22 => 30 | 23 => 8
  | 24 => 32 | 25 => 10 | 26 => 27 | 27 => 19 | 28 => 4  | 29 => 12
  | 30 => 10 | 31 => 33 | 32 => 18 | 33 => 26 | 34 => 13 | _ => 4

private def suppU_r2 (j : Fin 36) : Fin 36 :=
  match j.val with
  | 0 => 0   | 1 => 14  | 2 => 7   | 3 => 35  | 4 => 29  | 5 => 21
  | 6 => 15  | 7 => 1   | 8 => 35  | 9 => 6   | 10 => 20 | 11 => 29
  | 12 => 23 | 13 => 24 | 14 => 2  | 15 => 16 | 16 => 9  | 17 => 31
  | 18 => 25 | 19 => 22 | 20 => 17 | 21 => 3  | 22 => 30 | 23 => 8
  | 24 => 32 | 25 => 11 | 26 => 27 | 27 => 19 | 28 => 5  | 29 => 12
  | 30 => 11 | 31 => 33 | 32 => 18 | 33 => 26 | 34 => 13 | _ => 5

private def suppUG_r1 (j : Fin 36) : Fin 36 :=
  match j.val with
  | 0 => 0   | 1 => 8   | 2 => 13  | 3 => 23  | 4 => 28  | 5 => 28
  | 6 => 9   | 7 => 1   | 8 => 22  | 9 => 12  | 10 => 29 | 11 => 29
  | 12 => 25 | 13 => 35 | 14 => 2  | 15 => 10 | 16 => 15 | 17 => 18
  | 18 => 34 | 19 => 24 | 20 => 11 | 21 => 3  | 22 => 19 | 23 => 14
  | 24 => 17 | 25 => 21 | 26 => 30 | 27 => 26 | 28 => 4  | 29 => 4
  | 30 => 20 | 31 => 16 | 32 => 27 | 33 => 31 | 34 => 5  | _ => 5

private def suppUG_r2 (j : Fin 36) : Fin 36 :=
  match j.val with
  | 0 => 0   | 1 => 8   | 2 => 13  | 3 => 23  | 4 => 33  | 5 => 33
  | 6 => 9   | 7 => 1   | 8 => 22  | 9 => 12  | 10 => 32 | 11 => 32
  | 12 => 25 | 13 => 35 | 14 => 2  | 15 => 10 | 16 => 15 | 17 => 18
  | 18 => 34 | 19 => 24 | 20 => 11 | 21 => 3  | 22 => 19 | 23 => 14
  | 24 => 17 | 25 => 21 | 26 => 30 | 27 => 26 | 28 => 7  | 29 => 7
  | 30 => 20 | 31 => 16 | 32 => 27 | 33 => 31 | 34 => 6  | _ => 6

private def suppUR_r1 (j : Fin 36) : Fin 36 :=
  match j.val with
  | 0 => 0   | 1 => 22  | 2 => 29  | 3 => 8   | 4 => 31  | 5 => 15
  | 6 => 16  | 7 => 6   | 8 => 2   | 9 => 35  | 10 => 21 | 11 => 25
  | 12 => 33 | 13 => 1  | 14 => 12 | 15 => 28 | 16 => 11 | 17 => 11
  | 18 => 7  | 19 => 27 | 20 => 34 | 21 => 18 | 22 => 5  | 23 => 5
  | 24 => 23 | 25 => 32 | 26 => 9  | 27 => 13 | 28 => 4  | 29 => 4
  | 30 => 26 | 31 => 17 | 32 => 19 | 33 => 3  | 34 => 10 | _ => 10

private def suppUR_r2 (j : Fin 36) : Fin 36 :=
  match j.val with
  | 0 => 0   | 1 => 22  | 2 => 29  | 3 => 8   | 4 => 31  | 5 => 31
  | 6 => 16  | 7 => 6   | 8 => 2   | 9 => 35  | 10 => 25 | 25 => 32
  | 11 => 25 | 12 => 33 | 13 => 1  | 14 => 12 | 15 => 28 | 16 => 11
  | 17 => 20 | 18 => 7  | 19 => 27 | 20 => 34 | 21 => 18 | 22 => 14
  | 23 => 5  | 24 => 23 | 26 => 9  | 27 => 13 | 28 => 24 | 29 => 24
  | 30 => 26 | 31 => 17 | 32 => 19 | 33 => 3  | 34 => 30 | _ => 30

def countNnz36 {α : Type} [DecidableEq α] (zero : α) (M : Fin 36 → Fin 36 → α) : Nat :=
  sumFin36 (· + ·) fun r => sumFin36 (· + ·) fun c => if M r c == zero then 0 else 1

def defectNumEisen (scaleSq : Int) (M : Fin 36 → Fin 36 → EisenInt) : Int :=
  sumFin36 (· + ·) fun i => sumFin36 (· + ·) fun j =>
    let d := matDagMat eAdd eMul eConj M i j
    let err := if i == j then eSub d ⟨scaleSq, 0⟩ else d
    eNorm err

private def defectNumEisenSparse (scaleSq : Int) (M : Fin 36 → Fin 36 → EisenInt)
    (r1 r2 : Fin 36 → Fin 36) : Int :=
  sumFin36 (· + ·) fun i => sumFin36 (· + ·) fun j =>
    let d := sparseDotFast eZero eAdd eMul eConj M r1 r2 i j
    let err := if i == j then eSub d ⟨scaleSq, 0⟩ else d
    eNorm err

def defectSumSqrt13 (scaleSq : Int) (M : Fin 36 → Fin 36 → EisenSqrt13) : EisenSqrt13 :=
  sumFin36 qAdd fun i => sumFin36 qAdd fun j =>
    let d := matDagMat qAdd qMul qConj M i j
    let err := if i == j then qSub d ⟨⟨scaleSq, 0⟩, eZero⟩ else d
    qMul (qConj err) err

private def defectSumSqrt13Sparse (scaleSq : Int) (M : Fin 36 → Fin 36 → EisenSqrt13)
    (r1 r2 : Fin 36 → Fin 36) : EisenSqrt13 :=
  sumFin36 qAdd fun i => sumFin36 qAdd fun j =>
    let d := sparseDotFast qZero qAdd qMul qConj M r1 r2 i j
    let err := if i == j then qSub d ⟨⟨scaleSq, 0⟩, eZero⟩ else d
    qMul (qConj err) err

private theorem defectNumEisen_eq_sparse (scaleSq : Int) (M : Fin 36 → Fin 36 → EisenInt)
    (r1 r2 : Fin 36 → Fin 36)
    (hsupp : ∀ k j : Fin 36, k ≠ r1 j → k ≠ r2 j → M k j = eZero) :
    defectNumEisen scaleSq M = defectNumEisenSparse scaleSq M r1 r2 := by
  simp only [defectNumEisen, defectNumEisenSparse]
  have h : ∀ i j : Fin 36,
      matDagMat eAdd eMul eConj M i j =
        sparseDotFast eZero eAdd eMul eConj M r1 r2 i j :=
    matDagMat_eq_sparseDotFast eZero eAdd eMul eConj
      eAdd_zeroL eAdd_zeroR eAdd_assoc eAdd_comm eConj_zero eMul_zeroL eMul_zeroR M r1 r2 hsupp
  simp only [h]

private theorem defectSumSqrt13_eq_sparse (scaleSq : Int) (M : Fin 36 → Fin 36 → EisenSqrt13)
    (r1 r2 : Fin 36 → Fin 36)
    (hsupp : ∀ k j : Fin 36, k ≠ r1 j → k ≠ r2 j → M k j = qZero) :
    defectSumSqrt13 scaleSq M = defectSumSqrt13Sparse scaleSq M r1 r2 := by
  simp only [defectSumSqrt13, defectSumSqrt13Sparse]
  have h : ∀ i j : Fin 36,
      matDagMat qAdd qMul qConj M i j =
        sparseDotFast qZero qAdd qMul qConj M r1 r2 i j :=
    matDagMat_eq_sparseDotFast qZero qAdd qMul qConj
      qAdd_zeroL qAdd_zeroR qAdd_assoc qAdd_comm qConj_zero qMul_zeroL qMul_zeroR M r1 r2 hsupp
  simp only [h]

private theorem supp_of_allFin36 {α : Type} [DecidableEq α] (zero : α)
    (M : Fin 36 → Fin 36 → α) (r1 r2 : Fin 36 → Fin 36)
    (h : allFin36 (fun k => allFin36 (fun j =>
      (k == r1 j) || (k == r2 j) || (M k j == zero))) = true) :
    ∀ k j : Fin 36, k ≠ r1 j → k ≠ r2 j → M k j = zero := by
  intro k j hk1 hk2
  have hkj := allFin36_forall (allFin36_forall h k) j
  have hb1 : (k == r1 j) = false := beq_false_of_ne hk1
  have hb2 : (k == r2 j) = false := beq_false_of_ne hk2
  rw [hb1, hb2] at hkj
  simp only [Bool.false_or, beq_iff_eq] at hkj
  exact hkj

private theorem u0_supp : ∀ k j : Fin 36, k ≠ perm0Inv_6 j → k ≠ perm0Inv_6 j → U0_6 k j = eZero :=
  supp_of_allFin36 eZero U0_6 perm0Inv_6 perm0Inv_6 (by decide)

private theorem ug0_supp : ∀ k j : Fin 36, k ≠ suppUG_r1 j → k ≠ suppUG_r2 j → reshuffleGamma U0_6 k j = eZero :=
  supp_of_allFin36 eZero (reshuffleGamma U0_6) suppUG_r1 suppUG_r2 (by decide)

private theorem ur0_supp : ∀ k j : Fin 36, k ≠ suppUR_r1 j → k ≠ suppUR_r2 j → reshuffleR U0_6 k j = eZero :=
  supp_of_allFin36 eZero (reshuffleR U0_6) suppUR_r1 suppUR_r2 (by decide)

private theorem u40_supp : ∀ k j : Fin 36, k ≠ suppU_r1 j → k ≠ suppU_r2 j → U40_2 k j = eZero :=
  supp_of_allFin36 eZero U40_2 suppU_r1 suppU_r2 (by decide)

private theorem ug40_supp : ∀ k j : Fin 36, k ≠ suppUG_r1 j → k ≠ suppUG_r2 j → reshuffleGamma U40_2 k j = eZero :=
  supp_of_allFin36 eZero (reshuffleGamma U40_2) suppUG_r1 suppUG_r2 (by decide)

private theorem ur40_supp : ∀ k j : Fin 36, k ≠ suppUR_r1 j → k ≠ suppUR_r2 j → reshuffleR U40_2 k j = eZero :=
  supp_of_allFin36 eZero (reshuffleR U40_2) suppUR_r1 suppUR_r2 (by decide)

private theorem u44_supp : ∀ k j : Fin 36, k ≠ suppU_r1 j → k ≠ suppU_r2 j → U44_2s13 k j = qZero :=
  supp_of_allFin36 qZero U44_2s13 suppU_r1 suppU_r2 (by decide)

private theorem ug44_supp : ∀ k j : Fin 36, k ≠ suppUG_r1 j → k ≠ suppUG_r2 j → reshuffleGamma U44_2s13 k j = qZero :=
  supp_of_allFin36 qZero (reshuffleGamma U44_2s13) suppUG_r1 suppUG_r2 (by decide)

private theorem ur44_supp : ∀ k j : Fin 36, k ≠ suppUR_r1 j → k ≠ suppUR_r2 j → reshuffleR U44_2s13 k j = qZero :=
  supp_of_allFin36 qZero (reshuffleR U44_2s13) suppUR_r1 suppUR_r2 (by decide)

/-- Classical 6×6 witness `U0_6` is unitary (`U^† U = I`) with `D_Gamma = 0`, `D_R = 2`, `D_total = 2`. -/
theorem ame46_classical_defect_floor_2 :
    (∀ i j : Fin 36, matDagMat eAdd eMul eConj U0_6 i j = if i == j then eOne else eZero) ∧
    defectNumEisen 1 (reshuffleGamma U0_6) = 0 ∧
    defectNumEisen 1 (reshuffleR U0_6) = 4 ∧
    (0 : Int) + 4 = 2 * 2 := by
  refine ⟨?_, ?_, ?_, rfl⟩
  · intro i j
    rw [matDagMat_eq_sparseDotFast eZero eAdd eMul eConj
      eAdd_zeroL eAdd_zeroR eAdd_assoc eAdd_comm eConj_zero eMul_zeroL eMul_zeroR U0_6 perm0Inv_6 perm0Inv_6 u0_supp]
    have h : allFin36 (fun i => allFin36 (fun j =>
        sparseDotFast eZero eAdd eMul eConj U0_6 perm0Inv_6 perm0Inv_6 i j ==
          (if i == j then eOne else eZero))) = true := by decide
    exact beq_iff_eq.mp (allFin36_forall (allFin36_forall h i) j)
  · rw [defectNumEisen_eq_sparse 1 (reshuffleGamma U0_6) suppUG_r1 suppUG_r2 ug0_supp]
    decide
  · rw [defectNumEisen_eq_sparse 1 (reshuffleR U0_6) suppUR_r1 suppUR_r2 ur0_supp]
    decide

/-- Exact sparsity (`nnz = 40`) of the `(1/2) Z[ω]` 2-block gate `U40_2`. -/
theorem ame46_zw40_nnz_count : countNnz36 eZero U40_2 = 40 := by
  decide

/-- Exact unitarity (`U40_2^† U40_2 = 4 I`) of the 40-nonzero `(1/2) Z[ω]` gate. -/
theorem ame46_zw40_unitary (i j : Fin 36) :
    matDagMat eAdd eMul eConj U40_2 i j = if i == j then ⟨4, 0⟩ else eZero := by
  rw [matDagMat_eq_sparseDotFast eZero eAdd eMul eConj
    eAdd_zeroL eAdd_zeroR eAdd_assoc eAdd_comm eConj_zero eMul_zeroL eMul_zeroR U40_2 suppU_r1 suppU_r2 u40_supp]
  have h : allFin36 (fun i => allFin36 (fun j =>
      sparseDotFast eZero eAdd eMul eConj U40_2 suppU_r1 suppU_r2 i j ==
        (if i == j then ⟨4, 0⟩ else eZero))) = true := by decide
  exact beq_iff_eq.mp (allFin36_forall (allFin36_forall h i) j)

/-- Quantum escape below classical defect `2`: `D_Gamma = 1`, `D_R = 1/2`, `D_total = 3/2 < 2` (scale `2 * 4^2 = 32`). -/
theorem ame46_zw40_defect_escape :
    defectNumEisen 4 (reshuffleGamma U40_2) = 32 ∧
    defectNumEisen 4 (reshuffleR U40_2) = 16 ∧
    (32 : Int) + 16 < 2 * 32 := by
  refine ⟨?_, ?_, by decide⟩
  · rw [defectNumEisen_eq_sparse 4 (reshuffleGamma U40_2) suppUG_r1 suppUG_r2 ug40_supp]
    decide
  · rw [defectNumEisen_eq_sparse 4 (reshuffleR U40_2) suppUR_r1 suppUR_r2 ur40_supp]
    decide

/-- Exact sparsity (`nnz = 44`) of the `Z[ω][√13]` 4-block gate `U44_2s13`. -/
theorem ame46_zw44_nnz_count : countNnz36 qZero U44_2s13 = 44 := by
  decide

/-- Exact unitarity (`U44_2s13^† U44_2s13 = 52 I`) of the 44-nonzero `Z[ω]` gate. -/
theorem ame46_zw44_unitary (i j : Fin 36) :
    matDagMat qAdd qMul qConj U44_2s13 i j =
      if i == j then ⟨⟨52, 0⟩, eZero⟩ else qZero := by
  rw [matDagMat_eq_sparseDotFast qZero qAdd qMul qConj
    qAdd_zeroL qAdd_zeroR qAdd_assoc qAdd_comm qConj_zero qMul_zeroL qMul_zeroR U44_2s13 suppU_r1 suppU_r2 u44_supp]
  have h : allFin36 (fun i => allFin36 (fun j =>
      sparseDotFast qZero qAdd qMul qConj U44_2s13 suppU_r1 suppU_r2 i j ==
        (if i == j then ⟨⟨52, 0⟩, eZero⟩ else qZero))) = true := by decide
  exact beq_iff_eq.mp (allFin36_forall (allFin36_forall h i) j)

/-- Exact `D_Gamma = 3/13`, `D_R = 21/26`, `D_total = 27/26 < 2` (scale `2 * 52^2 = 5408`) for `U44_2s13`. -/
theorem ame46_zw44_defect_escape :
    defectSumSqrt13 52 (reshuffleGamma U44_2s13) = ⟨⟨1248, 0⟩, eZero⟩ ∧
    defectSumSqrt13 52 (reshuffleR U44_2s13) = ⟨⟨4368, 0⟩, eZero⟩ ∧
    (1248 : Int) * 13 = 3 * 5408 ∧
    (4368 : Int) * 26 = 21 * 5408 ∧
    (1248 + 4368 : Int) * 26 = 27 * 5408 ∧
    (1248 + 4368 : Int) < 2 * 5408 := by
  refine ⟨?_, ?_, by decide, by decide, by decide, by decide⟩
  · rw [defectSumSqrt13_eq_sparse 52 (reshuffleGamma U44_2s13) suppUG_r1 suppUG_r2 ug44_supp]
    decide
  · rw [defectSumSqrt13_eq_sparse 52 (reshuffleR U44_2s13) suppUR_r1 suppUR_r2 ur44_supp]
    decide

end Phy.HigherDimFrontier
