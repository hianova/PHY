import Flow.FiniteSum
import Std.Tactic.BVDecide

namespace Phy.Noether

/-!
  PHY v1.0.0 — Zero-Ghost Discrete Noether, Set-Theoretic Yang–Baxter & HoTT Core.
  Pure Lean 4 theorem module (zero `sorryAx`, zero extra axioms).
-/

/-! ## Section 1: 2-Bit-per-Site Pair Algebra (`BitVec 4 ↔ BitVec 2 × BitVec 2`) -/

def hi2 (u : BitVec 4) : BitVec 2 :=
  BitVec.extractLsb' 2 2 u

def lo2 (u : BitVec 4) : BitVec 2 :=
  BitVec.extractLsb' 0 2 u

def pack2 (l r : BitVec 2) : BitVec 4 :=
  l ++ r

def swap2 (u : BitVec 4) : BitVec 4 :=
  pack2 (lo2 u) (hi2 u)

def Bijective {α β : Type} (f : α → β) : Prop :=
  (∀ a b, f a = f b → a = b) ∧ (∀ b, ∃ a, f a = b)

theorem Bijective.comp {α β γ : Type} {g : β → γ} {f : α → β}
    (hg : Bijective g) (hf : Bijective f) : Bijective (fun x => g (f x)) := by
  refine ⟨fun a b h => hf.1 a b (hg.1 _ _ h), fun c => ?_⟩
  rcases hg.2 c with ⟨b, hb⟩
  rcases hf.2 b with ⟨a, ha⟩
  exact ⟨a, by dsimp; rw [ha, hb]⟩

theorem pack2_hi2_lo2 (u : BitVec 4) : pack2 (hi2 u) (lo2 u) = u := by
  simp only [pack2, hi2, lo2]
  bv_decide

theorem hi2_pack2 (l r : BitVec 2) : hi2 (pack2 l r) = l := by
  simp only [pack2, hi2]
  bv_decide

theorem lo2_pack2 (l r : BitVec 2) : lo2 (pack2 l r) = r := by
  simp only [pack2, lo2]
  bv_decide

theorem swap2_involutive (u : BitVec 4) : swap2 (swap2 u) = u := by
  simp only [swap2, pack2, hi2, lo2]
  bv_decide

theorem swap2_bijective : Bijective swap2 := by
  refine ⟨fun a b h => ?_, fun b => ⟨swap2 b, swap2_involutive b⟩⟩
  have h1 := congrArg swap2 h
  simpa [swap2_involutive] using h1

/-! ## Section 2: Brick-Wall Even/Odd Staggered Layers & Exact Reversibility -/

abbrev prevMod {N : Nat} (hN : 0 < N) (i : Fin N) : Fin N :=
  Flow.FiniteSum.prevIdx hN i

abbrev nextMod {N : Nat} (hN : 0 < N) (i : Fin N) : Fin N :=
  Flow.FiniteSum.nextIdx hN i

theorem nextMod_prevMod {N : Nat} (hN : 0 < N) (i : Fin N) :
    nextMod hN (prevMod hN i) = i :=
  Flow.FiniteSum.nextIdx_prevIdx hN i

theorem prevMod_nextMod {N : Nat} (hN : 0 < N) (i : Fin N) :
    prevMod hN (nextMod hN i) = i :=
  Flow.FiniteSum.prevIdx_nextIdx hN i

def evenLayer {N : Nat} (P : BitVec 4 → BitVec 4) (u : Fin N → BitVec 4) :
    Fin N → BitVec 4 :=
  fun i => P (u i)

def shiftHalf {N : Nat} (next : Fin N → Fin N) (u : Fin N → BitVec 4) :
    Fin N → BitVec 4 :=
  fun i => pack2 (lo2 (u i)) (hi2 (u (next i)))

def unshiftHalf {N : Nat} (prev : Fin N → Fin N) (v : Fin N → BitVec 4) :
    Fin N → BitVec 4 :=
  fun i => pack2 (lo2 (v (prev i))) (hi2 (v i))

theorem shiftHalf_unshiftHalf {N : Nat} (hN : 0 < N) (v : Fin N → BitVec 4) :
    shiftHalf (nextMod hN) (unshiftHalf (prevMod hN) v) = v := by
  funext i
  simp only [shiftHalf, unshiftHalf, lo2_pack2, hi2_pack2, prevMod_nextMod, pack2_hi2_lo2]

theorem unshiftHalf_shiftHalf {N : Nat} (hN : 0 < N) (u : Fin N → BitVec 4) :
    unshiftHalf (prevMod hN) (shiftHalf (nextMod hN) u) = u := by
  funext i
  simp only [shiftHalf, unshiftHalf, lo2_pack2, hi2_pack2, nextMod_prevMod, pack2_hi2_lo2]

theorem shiftHalf_bijective {N : Nat} (hN : 0 < N) :
    Bijective (shiftHalf (nextMod hN)) := by
  refine ⟨fun u v h => ?_, fun v => ⟨unshiftHalf (prevMod hN) v, shiftHalf_unshiftHalf hN v⟩⟩
  have h1 := congrArg (unshiftHalf (prevMod hN)) h
  simpa [unshiftHalf_shiftHalf hN] using h1

theorem evenLayer_bijective {N : Nat} {P : BitVec 4 → BitVec 4}
    (hP : Bijective P) : Bijective (@evenLayer N P) := by
  refine ⟨fun u v h => ?_, fun v => ?_⟩
  · funext i
    exact hP.1 _ _ (congrFun h i)
  · let inv : BitVec 4 → BitVec 4 := fun y => Classical.choose (hP.2 y)
    have hinv : ∀ y, P (inv y) = y := fun y => Classical.choose_spec (hP.2 y)
    exact ⟨fun i => inv (v i), funext fun i => hinv (v i)⟩

def floquetRing {N : Nat} (hN : 0 < N) (P : BitVec 4 → BitVec 4)
    (u : Fin N → BitVec 4) : Fin N → BitVec 4 :=
  unshiftHalf (prevMod hN) (evenLayer P (shiftHalf (nextMod hN) (evenLayer P u)))

theorem floquetRing_bijective {N : Nat} (hN : 0 < N) {P : BitVec 4 → BitVec 4}
    (hP : Bijective P) : Bijective (floquetRing hN P) := by
  have hShift := shiftHalf_bijective hN
  have hUnshift : Bijective (unshiftHalf (prevMod hN)) := by
    refine ⟨fun u v h => ?_, fun u => ⟨shiftHalf (nextMod hN) u, unshiftHalf_shiftHalf hN u⟩⟩
    have h1 := congrArg (shiftHalf (nextMod hN)) h
    simpa [shiftHalf_unshiftHalf hN] using h1
  have hEven := @evenLayer_bijective N P hP
  exact hUnshift.comp (hEven.comp (hShift.comp hEven))

def step3 (P : BitVec 4 → BitVec 4) (a b c : BitVec 4) : BitVec 4 :=
  let pa := P a
  let pb := P b
  let pc := P c
  let oddL := P (pack2 (lo2 pa) (hi2 pb))
  let oddR := P (pack2 (lo2 pb) (hi2 pc))
  pack2 (lo2 oddL) (hi2 oddR)

theorem floquetRing_eq_step3 {N : Nat} (hN : 0 < N) (P : BitVec 4 → BitVec 4)
    (u : Fin N → BitVec 4) (i : Fin N) :
    floquetRing hN P u i = step3 P (u (prevMod hN i)) (u i) (u (nextMod hN i)) := by
  simp only [floquetRing, unshiftHalf, evenLayer, shiftHalf, step3, nextMod_prevMod]

theorem floquetRing_vacuum_invariant {N : Nat} (hN : 0 < N) (P : BitVec 4 → BitVec 4)
    (hP0 : P 0#4 = 0#4) (i : Fin N) :
    floquetRing hN P (fun _ => 0#4) i = 0#4 := by
  rw [floquetRing_eq_step3]
  have hhi : hi2 0#4 = 0#2 := by simp only [hi2]; bv_decide
  have hlo : lo2 0#4 = 0#2 := by simp only [lo2]; bv_decide
  have hp : pack2 0#2 0#2 = 0#4 := by simp only [pack2]; bv_decide
  simp only [step3, hP0, hhi, hlo, hp]

/-! ## Section 3: Exact Finite Ring Sums, Telescoping & Discrete Noether Conservation -/

abbrev sumPrefix {N : Nat} (f : Fin N → Int) : Nat → Int :=
  Flow.FiniteSum.sumPrefix f

abbrev ringSum (N : Nat) (f : Fin N → Int) : Int :=
  Flow.FiniteSum.sumFin N f

theorem sumPrefix_add {N : Nat} (f g : Fin N → Int) (k : Nat) :
    sumPrefix (fun i => f i + g i) k = sumPrefix f k + sumPrefix g k :=
  Flow.FiniteSum.sumPrefix_add f g k

theorem ringSum_add (N : Nat) (f g : Fin N → Int) :
    ringSum N (fun i => f i + g i) = ringSum N f + ringSum N g :=
  Flow.FiniteSum.sumFin_add N f g

theorem sumPrefix_sub {N : Nat} (f g : Fin N → Int) (k : Nat) :
    sumPrefix (fun i => f i - g i) k = sumPrefix f k - sumPrefix g k :=
  Flow.FiniteSum.sumPrefix_sub f g k

theorem ringSum_sub (N : Nat) (f g : Fin N → Int) :
    ringSum N (fun i => f i - g i) = ringSum N f - ringSum N g :=
  Flow.FiniteSum.sumFin_sub N f g

theorem ringSum_shift_next {N : Nat} (hN : 0 < N) (f : Fin N → Int) :
    ringSum N (fun i => f (nextMod hN i)) = ringSum N f :=
  Flow.FiniteSum.sumFin_shift_next hN f

theorem ringSum_telescoping {N : Nat} (hN : 0 < N) (j : Fin N → Int) :
    ringSum N (fun i => j (nextMod hN i) - j i) = 0 :=
  Flow.FiniteSum.sumFin_telescoping hN j

def SatisfiesContinuity (P : BitVec 4 → BitVec 4) (q : BitVec 4 → Int)
    (j : BitVec 4 → BitVec 4 → Int) : Prop :=
  ∀ a b c : BitVec 4, q (step3 P a b c) - q b = j b c - j a b

theorem discrete_noether_global_conservation {N : Nat} (hN : 0 < N)
    (P : BitVec 4 → BitVec 4) (q : BitVec 4 → Int) (j : BitVec 4 → BitVec 4 → Int)
    (hCont : SatisfiesContinuity P q j) (u : Fin N → BitVec 4) :
    ringSum N (fun i => q (floquetRing hN P u i)) = ringSum N (fun i => q (u i)) := by
  let edgeJ : Fin N → Int := fun i => j (u (prevMod hN i)) (u i)
  have hdiff : ringSum N (fun i => q (floquetRing hN P u i) - q (u i)) =
      ringSum N (fun i => edgeJ (nextMod hN i) - edgeJ i) := by
    congr 1
    funext i
    rw [floquetRing_eq_step3]
    have he := hCont (u (prevMod hN i)) (u i) (u (nextMod hN i))
    have hpn : prevMod hN (nextMod hN i) = i := prevMod_nextMod hN i
    simp only [edgeJ, hpn]
    exact he
  rw [ringSum_sub] at hdiff
  rw [ringSum_telescoping hN edgeJ] at hdiff
  omega

def canonCurrent (P : BitVec 4 → BitVec 4) (q : BitVec 4 → Int) (b c : BitVec 4) : Int :=
  (q (step3 P 0#4 0#4 b) - q 0#4) + (q (step3 P 0#4 b c) - q b)

theorem telescoping_criterion_iff_exists_current
    (P : BitVec 4 → BitVec 4) (q : BitVec 4 → Int) :
    (∃ j, SatisfiesContinuity P q j) ↔
      SatisfiesContinuity P q (canonCurrent P q) := by
  constructor
  · rintro ⟨j, hj⟩ a b c
    have h1 := hj 0#4 0#4 a
    have h2 := hj 0#4 a b
    have h3 := hj 0#4 0#4 b
    have h4 := hj 0#4 b c
    have h5 := hj a b c
    simp only [canonCurrent]
    omega
  · intro hc
    exact ⟨canonCurrent P q, hc⟩

theorem discrete_noether_n1_invariant (P : BitVec 4 → BitVec 4) (q : BitVec 4 → Int)
    (j : BitVec 4 → BitVec 4 → Int) (hCont : SatisfiesContinuity P q j)
    (b : BitVec 4) : q (step3 P b b b) = q b := by
  have h := hCont b b b
  omega

theorem discrete_noether_n2_invariant (P : BitVec 4 → BitVec 4) (q : BitVec 4 → Int)
    (j : BitVec 4 → BitVec 4 → Int) (hCont : SatisfiesContinuity P q j)
    (a b : BitVec 4) :
    q (step3 P b a b) + q (step3 P a b a) = q a + q b := by
  have h1 := hCont b a b
  have h2 := hCont a b a
  omega

theorem discrete_noether_n3_invariant (P : BitVec 4 → BitVec 4) (q : BitVec 4 → Int)
    (j : BitVec 4 → BitVec 4 → Int) (hCont : SatisfiesContinuity P q j)
    (a b c : BitVec 4) :
    q (step3 P c a b) + q (step3 P a b c) + q (step3 P b c a) = q a + q b + q c := by
  have h1 := hCont c a b
  have h2 := hCont a b c
  have h3 := hCont b c a
  omega

/-! ## Section 4: `BitVec 32` Verified `bv_decide` Bridge to Unbounded `Int` -/

def canonCurrent32 (P : BitVec 4 → BitVec 4) (q : BitVec 4 → BitVec 32)
    (b c : BitVec 4) : BitVec 32 :=
  (q (step3 P 0#4 0#4 b) - q 0#4) + (q (step3 P 0#4 b c) - q b)

def BoundedBV32 (bound : Int) (x : BitVec 32) : Prop :=
  -bound ≤ x.toInt ∧ x.toInt ≤ bound

private theorem bmod_pow32_eq_self {x : Int}
    (hlo : -(2 ^ 31 : Int) ≤ x) (hhi : x < 2 ^ 31) :
    x.bmod (2 ^ 32) = x := by
  unfold Int.bmod
  dsimp only
  have hn : ((2 ^ 32 : Nat) : Int) = 4294967296 := rfl
  rw [hn]
  by_cases hnonneg : 0 ≤ x
  · have hm : x % 4294967296 = x := Int.emod_eq_of_lt hnonneg (by omega)
    rw [hm]
    split <;> omega
  · have hxpos : 0 ≤ x + 4294967296 := by omega
    have hxlt : x + 4294967296 < 4294967296 := by omega
    have hm : x % 4294967296 = x + 4294967296 := by
      have : x % 4294967296 = (x + 4294967296) % 4294967296 := by omega
      rw [this]
      exact Int.emod_eq_of_lt hxpos hxlt
    rw [hm]
    split <;> omega

private theorem toInt_add_of_range {a b : BitVec 32}
    (hlo : -(2 ^ 31 : Int) ≤ a.toInt + b.toInt) (hhi : a.toInt + b.toInt < 2 ^ 31) :
    (a + b).toInt = a.toInt + b.toInt := by
  rw [BitVec.toInt_add]
  exact bmod_pow32_eq_self hlo hhi

private theorem toInt_sub_of_range {a b : BitVec 32}
    (hlo : -(2 ^ 31 : Int) ≤ a.toInt - b.toInt) (hhi : a.toInt - b.toInt < 2 ^ 31) :
    (a - b).toInt = a.toInt - b.toInt := by
  rw [BitVec.toInt_sub]
  exact bmod_pow32_eq_self hlo hhi

private theorem canonCurrent32_toInt (P : BitVec 4 → BitVec 4) (q : BitVec 4 → BitVec 32)
    (hb : ∀ u : BitVec 4, BoundedBV32 100000000 (q u)) (b c : BitVec 4) :
    (canonCurrent32 P q b c).toInt = canonCurrent P (fun u => (q u).toInt) b c ∧
      -(400000000 : Int) ≤ (canonCurrent32 P q b c).toInt ∧
      (canonCurrent32 P q b c).toInt ≤ 400000000 := by
  have h1 := hb (step3 P 0#4 0#4 b)
  have h2 := hb 0#4
  have h3 := hb (step3 P 0#4 b c)
  have h4 := hb b
  simp only [BoundedBV32] at h1 h2 h3 h4
  have hd1 : (q (step3 P 0#4 0#4 b) - q 0#4).toInt =
      (q (step3 P 0#4 0#4 b)).toInt - (q 0#4).toInt :=
    toInt_sub_of_range (by omega) (by omega)
  have hd2 : (q (step3 P 0#4 b c) - q b).toInt =
      (q (step3 P 0#4 b c)).toInt - (q b).toInt :=
    toInt_sub_of_range (by omega) (by omega)
  have hadd : (canonCurrent32 P q b c).toInt =
      (q (step3 P 0#4 0#4 b) - q 0#4).toInt + (q (step3 P 0#4 b c) - q b).toInt := by
    unfold canonCurrent32
    exact toInt_add_of_range (by omega) (by omega)
  refine ⟨?_, by omega, by omega⟩
  simp only [canonCurrent]
  omega

theorem bitvec32_continuity_lifts_to_int
    (P : BitVec 4 → BitVec 4) (q : BitVec 4 → BitVec 32)
    (hb : ∀ u : BitVec 4, BoundedBV32 100000000 (q u))
    (h32 : ∀ a b c : BitVec 4,
      q (step3 P a b c) - q b = canonCurrent32 P q b c - canonCurrent32 P q a b) :
    SatisfiesContinuity P (fun u => (q u).toInt) (canonCurrent P (fun u => (q u).toInt)) := by
  intro a b c
  have eq32 := congrArg BitVec.toInt (h32 a b c)
  have ha := hb (step3 P a b c)
  have hbv := hb b
  simp only [BoundedBV32] at ha hbv
  have hL : (q (step3 P a b c) - q b).toInt = (q (step3 P a b c)).toInt - (q b).toInt :=
    toInt_sub_of_range (by omega) (by omega)
  obtain ⟨hbcEq, hbcLo, hbcHi⟩ := canonCurrent32_toInt P q hb b c
  obtain ⟨habEq, habLo, habHi⟩ := canonCurrent32_toInt P q hb a b
  have hR : (canonCurrent32 P q b c - canonCurrent32 P q a b).toInt =
      (canonCurrent32 P q b c).toInt - (canonCurrent32 P q a b).toInt :=
    toInt_sub_of_range (by omega) (by omega)
  dsimp only
  omega

theorem forall_ring_conservation_of_bv32 {N : Nat} (hN : 0 < N)
    (P : BitVec 4 → BitVec 4) (q : BitVec 4 → BitVec 32)
    (hb : ∀ u : BitVec 4, BoundedBV32 100000000 (q u))
    (h32 : ∀ a b c : BitVec 4,
      q (step3 P a b c) - q b = canonCurrent32 P q b c - canonCurrent32 P q a b)
    (u : Fin N → BitVec 4) :
    ringSum N (fun i => (q (floquetRing hN P u i)).toInt) =
      ringSum N (fun i => (q (u i)).toInt) :=
  discrete_noether_global_conservation hN P (fun x => (q x).toInt)
    (canonCurrent P (fun x => (q x).toInt))
    (bitvec32_continuity_lifts_to_int P q hb h32) u

/-! ## Section 5: Non-Separability, Nonlinearity (GF(2) ANF Deg ≥ 2), & Non-Factorizability -/

def IsSeparable (q : BitVec 4 → Int) : Prop :=
  ∃ f g : BitVec 2 → Int, ∀ l r : BitVec 2, q (pack2 l r) = f l + g r

theorem nonseparable_of_cross_diff {q : BitVec 4 → Int}
    (l0 l1 r0 r1 : BitVec 2)
    (h : q (pack2 l1 r1) - q (pack2 l1 r0) - q (pack2 l0 r1) + q (pack2 l0 r0) ≠ 0) :
    ¬ IsSeparable q := by
  rintro ⟨f, g, hfg⟩
  rw [hfg l1 r1, hfg l1 r0, hfg l0 r1, hfg l0 r0] at h
  omega

theorem active_current_witness (P : BitVec 4 → BitVec 4) (q : BitVec 4 → Int)
    (b c : BitVec 4) (h : canonCurrent P q b c ≠ 0) :
    ∃ x y, canonCurrent P q x y ≠ 0 :=
  ⟨b, c, h⟩

def IsAffineGF2 (f : BitVec 4 → Bool) : Prop :=
  ∃ c c0 c1 c2 c3 : Bool, ∀ x : BitVec 4,
    f x = Bool.xor c
      (Bool.xor (and c0 (BitVec.getLsbD x 0))
        (Bool.xor (and c1 (BitVec.getLsbD x 1))
          (Bool.xor (and c2 (BitVec.getLsbD x 2)) (and c3 (BitVec.getLsbD x 3)))))

private theorem xor4_cancel (c c0 c1 c2 c3 a0 a1 a2 a3 b0 b1 b2 b3 : Bool) :
    Bool.xor
      (Bool.xor
        (Bool.xor c (Bool.xor (c0 && Bool.xor a0 b0)
          (Bool.xor (c1 && Bool.xor a1 b1)
            (Bool.xor (c2 && Bool.xor a2 b2) (c3 && Bool.xor a3 b3)))))
        (Bool.xor c (Bool.xor (c0 && a0)
          (Bool.xor (c1 && a1) (Bool.xor (c2 && a2) (c3 && a3))))))
      (Bool.xor
        (Bool.xor c (Bool.xor (c0 && b0)
          (Bool.xor (c1 && b1) (Bool.xor (c2 && b2) (c3 && b3)))))
        c) = false := by
  bv_decide

theorem gf2_second_diff_witnesses_nonlinear {f : BitVec 4 → Bool} (a b : BitVec 4)
    (h : Bool.xor (Bool.xor (f (a ^^^ b)) (f a)) (Bool.xor (f b) (f 0#4)) = true) :
    ¬ IsAffineGF2 f := by
  rintro ⟨c, c0, c1, c2, c3, hf⟩
  have h0 : BitVec.getLsbD 0#4 0 = false ∧ BitVec.getLsbD 0#4 1 = false ∧
      BitVec.getLsbD 0#4 2 = false ∧ BitVec.getLsbD 0#4 3 = false := by decide
  rw [hf (a ^^^ b), hf a, hf b, hf 0#4] at h
  simp only [BitVec.getLsbD_xor, h0.1, h0.2.1, h0.2.2.1, h0.2.2.2,
    Bool.and_false, Bool.xor_false] at h
  have hcancel := xor4_cancel c c0 c1 c2 c3
    (BitVec.getLsbD a 0) (BitVec.getLsbD a 1) (BitVec.getLsbD a 2) (BitVec.getLsbD a 3)
    (BitVec.getLsbD b 0) (BitVec.getLsbD b 1) (BitVec.getLsbD b 2) (BitVec.getLsbD b 3)
  rw [hcancel] at h
  exact Bool.false_ne_true h

def IsFactorizable (P : BitVec 4 → BitVec 4) : Prop :=
  ∃ fl fr : BitVec 2 → BitVec 2, ∀ l r : BitVec 2, P (pack2 l r) = pack2 (fl l) (fr r)

theorem non_factorizable_of_hi_cross {P : BitVec 4 → BitVec 4}
    (l r0 r1 : BitVec 2) (h : hi2 (P (pack2 l r1)) ≠ hi2 (P (pack2 l r0))) :
    ¬ IsFactorizable P := by
  rintro ⟨fl, fr, hflfr⟩
  rw [hflfr l r1, hflfr l r0, hi2_pack2, hi2_pack2] at h
  exact h rfl

theorem non_factorizable_of_lo_cross {P : BitVec 4 → BitVec 4}
    (l0 l1 r : BitVec 2) (h : lo2 (P (pack2 l1 r)) ≠ lo2 (P (pack2 l0 r))) :
    ¬ IsFactorizable P := by
  rintro ⟨fl, fr, hflfr⟩
  rw [hflfr l1 r, hflfr l0 r, lo2_pack2, lo2_pack2] at h
  exact h rfl

/-! ## Section 6: Set-Theoretic Yang–Baxter Equation, PT Involution & HoTT Fibers -/

def ybApply12 (R : BitVec 4 → BitVec 4) (t : BitVec 6) : BitVec 6 :=
  let a := BitVec.extractLsb' 4 2 t
  let b := BitVec.extractLsb' 2 2 t
  let c := BitVec.extractLsb' 0 2 t
  R (pack2 a b) ++ c

def ybApply23 (R : BitVec 4 → BitVec 4) (t : BitVec 6) : BitVec 6 :=
  let a := BitVec.extractLsb' 4 2 t
  let b := BitVec.extractLsb' 2 2 t
  let c := BitVec.extractLsb' 0 2 t
  a ++ R (pack2 b c)

def IsSetTheoreticYB (R : BitVec 4 → BitVec 4) : Prop :=
  ∀ t : BitVec 6,
    ybApply12 R (ybApply23 R (ybApply12 R t)) = ybApply23 R (ybApply12 R (ybApply23 R t))

theorem yb_id_is_solution : IsSetTheoreticYB id := by
  intro t
  simp only [ybApply12, ybApply23, pack2, id]
  bv_decide

theorem yb_swap_is_solution : IsSetTheoreticYB swap2 := by
  intro t
  simp only [ybApply12, ybApply23, swap2, hi2, lo2, pack2]
  bv_decide

def ptConjugate (P : BitVec 4 → BitVec 4) (u : BitVec 4) : BitVec 4 :=
  swap2 (P (swap2 u))

theorem pt_involution_self_inverse (P : BitVec 4 → BitVec 4) (u : BitVec 4) :
    ptConjugate (ptConjugate P) u = P u := by
  simp only [ptConjugate, swap2_involutive]

/-! ## Section 7: Division-Free $M_{61}$ Arithmetic & Binary Stein GCD Foundations -/

def M61 : Nat := 2 ^ 61 - 1

theorem m61_pow61_mod : 2 ^ 61 % M61 = 1 := by
  decide

theorem m61_fold_congruent (x : Nat) :
    ((x % 2 ^ 61) + (x / 2 ^ 61)) % M61 = x % M61 := by
  have hdiv : (x / 2 ^ 61) * 2 ^ 61 + (x % 2 ^ 61) = x := Nat.div_add_mod' x (2 ^ 61)
  symm
  calc x % M61
      = ((x / 2 ^ 61) * 2 ^ 61 + (x % 2 ^ 61)) % M61 := by rw [hdiv]
    _ = (((x / 2 ^ 61) * 2 ^ 61) % M61 + (x % 2 ^ 61) % M61) % M61 := Nat.add_mod _ _ _
    _ = ((((x / 2 ^ 61) % M61 * (2 ^ 61 % M61)) % M61) + (x % 2 ^ 61) % M61) % M61 := by
        rw [Nat.mul_mod]
    _ = (((x / 2 ^ 61) % M61 * 1) % M61 + (x % 2 ^ 61) % M61) % M61 := by
        rw [m61_pow61_mod]
    _ = ((x / 2 ^ 61) % M61 + (x % 2 ^ 61) % M61) % M61 := by
        rw [Nat.mul_one, Nat.mod_mod]
    _ = ((x / 2 ^ 61) + (x % 2 ^ 61)) % M61 := (Nat.add_mod _ _ _).symm
    _ = ((x % 2 ^ 61) + (x / 2 ^ 61)) % M61 := by rw [Nat.add_comm]

theorem stein_gcd_zero_right (a : Nat) : Nat.gcd a 0 = a :=
  Nat.gcd_zero_right a

theorem stein_gcd_even_left (a b : Nat) (hb : b % 2 = 1) :
    Nat.gcd (2 * a) b = Nat.gcd a b := by
  have hcop : Nat.Coprime 2 b := by
    rw [Nat.Coprime]
    have hd2 : Nat.gcd 2 b ∣ 2 := Nat.gcd_dvd_left 2 b
    have hdb : Nat.gcd 2 b ∣ b := Nat.gcd_dvd_right 2 b
    have hle : Nat.gcd 2 b ≤ 2 := Nat.le_of_dvd (by decide) hd2
    have hpos : 0 < Nat.gcd 2 b := Nat.gcd_pos_of_pos_left b (by decide)
    have hne2 : Nat.gcd 2 b ≠ 2 := by
      intro h2
      rw [h2] at hdb
      have hmod : b % 2 = 0 := Nat.mod_eq_zero_of_dvd hdb
      omega
    omega
  rw [Nat.Coprime.gcd_mul_left_cancel a hcop]

structure PhyFiber (α : Type) where
  state : α
  merkle : UInt64
  steps : Nat
  deriving DecidableEq, Repr

theorem phy_fiber_path_refl {α : Type} (x y : PhyFiber α) (h : x = y) :
    x.merkle = y.merkle := by
  rw [h]

theorem phy_fiber_path_symm {α : Type} {x y : PhyFiber α} (h : x = y) : y = x :=
  h.symm

theorem phy_fiber_path_trans {α : Type} {x y z : PhyFiber α}
    (h1 : x = y) (h2 : y = z) : x = z :=
  h1.trans h2

/-! ## Section 8: 4-Site (2-Block) Higher-Order Discrete Noether Conservation Laws ($Q_4$) -/

def SatisfiesContinuity4 (P : BitVec 4 → BitVec 4) (q4 : BitVec 4 → BitVec 4 → Int)
    (j4 : BitVec 4 → BitVec 4 → BitVec 4 → Int) : Prop :=
  ∀ a b c d : BitVec 4, q4 (step3 P a b c) (step3 P b c d) - q4 b c = j4 b c d - j4 a b c

def canonCurrent4 (P : BitVec 4 → BitVec 4) (q4 : BitVec 4 → BitVec 4 → Int)
    (a b c : BitVec 4) : Int :=
  (q4 (step3 P 0#4 0#4 0#4) (step3 P 0#4 0#4 a) - q4 0#4 0#4) +
    (q4 (step3 P 0#4 0#4 a) (step3 P 0#4 a b) - q4 0#4 a) +
    (q4 (step3 P 0#4 a b) (step3 P a b c) - q4 a b)

theorem telescoping4_criterion_iff_exists_current4
    (P : BitVec 4 → BitVec 4) (q4 : BitVec 4 → BitVec 4 → Int) :
    (∃ j4, SatisfiesContinuity4 P q4 j4) ↔
      SatisfiesContinuity4 P q4 (canonCurrent4 P q4) := by
  constructor
  · rintro ⟨j4, hj⟩ a b c d
    have ha1 := hj 0#4 0#4 0#4 a
    have ha2 := hj 0#4 0#4 a b
    have ha3 := hj 0#4 a b c
    have hb1 := hj 0#4 0#4 0#4 b
    have hb2 := hj 0#4 0#4 b c
    have hb3 := hj 0#4 b c d
    have habcd := hj a b c d
    simp only [canonCurrent4]
    omega
  · intro hc
    exact ⟨canonCurrent4 P q4, hc⟩

theorem shift_gauge_sum_zero {N : Nat} (hN : 0 < N) (g : BitVec 4 → Int)
    (u : Fin N → BitVec 4) :
    ringSum N (fun i => g (u (nextMod hN i)) - g (u i)) = 0 :=
  ringSum_telescoping hN (fun i => g (u i))

theorem discrete_noether4_global_conservation {N : Nat} (hN : 0 < N)
    (P : BitVec 4 → BitVec 4) (q4 : BitVec 4 → BitVec 4 → Int)
    (j4 : BitVec 4 → BitVec 4 → BitVec 4 → Int)
    (hCont4 : SatisfiesContinuity4 P q4 j4) (u : Fin N → BitVec 4) :
    ringSum N (fun i => q4 (floquetRing hN P u i) (floquetRing hN P u (nextMod hN i))) =
      ringSum N (fun i => q4 (u i) (u (nextMod hN i))) := by
  let edgeJ4 : Fin N → Int :=
    fun i => j4 (u (prevMod hN i)) (u i) (u (nextMod hN i))
  have hdiff :
      ringSum N
          (fun i => q4 (floquetRing hN P u i) (floquetRing hN P u (nextMod hN i)) -
            q4 (u i) (u (nextMod hN i))) =
        ringSum N (fun i => edgeJ4 (nextMod hN i) - edgeJ4 i) := by
    congr 1
    funext i
    have hpn : prevMod hN (nextMod hN i) = i := prevMod_nextMod hN i
    rw [floquetRing_eq_step3 hN P u i, floquetRing_eq_step3 hN P u (nextMod hN i), hpn]
    have he :=
      hCont4 (u (prevMod hN i)) (u i) (u (nextMod hN i)) (u (nextMod hN (nextMod hN i)))
    simpa [edgeJ4, hpn] using he
  rw [ringSum_sub, ringSum_telescoping hN edgeJ4] at hdiff
  omega

def canonCurrent4_32 (P : BitVec 4 → BitVec 4) (q4 : BitVec 4 → BitVec 4 → BitVec 32)
    (a b c : BitVec 4) : BitVec 32 :=
  (q4 (step3 P 0#4 0#4 0#4) (step3 P 0#4 0#4 a) - q4 0#4 0#4) +
    (q4 (step3 P 0#4 0#4 a) (step3 P 0#4 a b) - q4 0#4 a) +
    (q4 (step3 P 0#4 a b) (step3 P a b c) - q4 a b)

private theorem canonCurrent4_32_toInt (P : BitVec 4 → BitVec 4)
    (q4 : BitVec 4 → BitVec 4 → BitVec 32)
    (hb : ∀ u v : BitVec 4, BoundedBV32 10000000 (q4 u v)) (a b c : BitVec 4) :
    (canonCurrent4_32 P q4 a b c).toInt =
        canonCurrent4 P (fun u v => (q4 u v).toInt) a b c ∧
      -(60000000 : Int) ≤ (canonCurrent4_32 P q4 a b c).toInt ∧
      (canonCurrent4_32 P q4 a b c).toInt ≤ 60000000 := by
  have h1 := hb (step3 P 0#4 0#4 0#4) (step3 P 0#4 0#4 a)
  have h2 := hb 0#4 0#4
  have h3 := hb (step3 P 0#4 0#4 a) (step3 P 0#4 a b)
  have h4 := hb 0#4 a
  have h5 := hb (step3 P 0#4 a b) (step3 P a b c)
  have h6 := hb a b
  simp only [BoundedBV32] at h1 h2 h3 h4 h5 h6
  have hd1 : (q4 (step3 P 0#4 0#4 0#4) (step3 P 0#4 0#4 a) - q4 0#4 0#4).toInt =
      (q4 (step3 P 0#4 0#4 0#4) (step3 P 0#4 0#4 a)).toInt - (q4 0#4 0#4).toInt :=
    toInt_sub_of_range (by omega) (by omega)
  have hd2 : (q4 (step3 P 0#4 0#4 a) (step3 P 0#4 a b) - q4 0#4 a).toInt =
      (q4 (step3 P 0#4 0#4 a) (step3 P 0#4 a b)).toInt - (q4 0#4 a).toInt :=
    toInt_sub_of_range (by omega) (by omega)
  have hd3 : (q4 (step3 P 0#4 a b) (step3 P a b c) - q4 a b).toInt =
      (q4 (step3 P 0#4 a b) (step3 P a b c)).toInt - (q4 a b).toInt :=
    toInt_sub_of_range (by omega) (by omega)
  have had1 :
      ((q4 (step3 P 0#4 0#4 0#4) (step3 P 0#4 0#4 a) - q4 0#4 0#4) +
        (q4 (step3 P 0#4 0#4 a) (step3 P 0#4 a b) - q4 0#4 a)).toInt =
      (q4 (step3 P 0#4 0#4 0#4) (step3 P 0#4 0#4 a) - q4 0#4 0#4).toInt +
        (q4 (step3 P 0#4 0#4 a) (step3 P 0#4 a b) - q4 0#4 a).toInt :=
    toInt_add_of_range (by omega) (by omega)
  have had2 : (canonCurrent4_32 P q4 a b c).toInt =
      ((q4 (step3 P 0#4 0#4 0#4) (step3 P 0#4 0#4 a) - q4 0#4 0#4) +
        (q4 (step3 P 0#4 0#4 a) (step3 P 0#4 a b) - q4 0#4 a)).toInt +
        (q4 (step3 P 0#4 a b) (step3 P a b c) - q4 a b).toInt := by
    unfold canonCurrent4_32
    exact toInt_add_of_range (by omega) (by omega)
  refine ⟨?_, by omega, by omega⟩
  simp only [canonCurrent4]
  omega

theorem bitvec32_continuity4_lifts_to_int
    (P : BitVec 4 → BitVec 4) (q4 : BitVec 4 → BitVec 4 → BitVec 32)
    (hb : ∀ u v : BitVec 4, BoundedBV32 10000000 (q4 u v))
    (h32 : ∀ a b c d : BitVec 4,
      q4 (step3 P a b c) (step3 P b c d) - q4 b c =
        canonCurrent4_32 P q4 b c d - canonCurrent4_32 P q4 a b c) :
    SatisfiesContinuity4 P (fun u v => (q4 u v).toInt)
      (canonCurrent4 P (fun u v => (q4 u v).toInt)) := by
  intro a b c d
  have eq32 := congrArg BitVec.toInt (h32 a b c d)
  have ha := hb (step3 P a b c) (step3 P b c d)
  have hbc := hb b c
  simp only [BoundedBV32] at ha hbc
  have hL : (q4 (step3 P a b c) (step3 P b c d) - q4 b c).toInt =
      (q4 (step3 P a b c) (step3 P b c d)).toInt - (q4 b c).toInt :=
    toInt_sub_of_range (by omega) (by omega)
  obtain ⟨hbcdEq, hbcdLo, hbcdHi⟩ := canonCurrent4_32_toInt P q4 hb b c d
  obtain ⟨habcEq, habcLo, habcHi⟩ := canonCurrent4_32_toInt P q4 hb a b c
  have hR : (canonCurrent4_32 P q4 b c d - canonCurrent4_32 P q4 a b c).toInt =
      (canonCurrent4_32 P q4 b c d).toInt - (canonCurrent4_32 P q4 a b c).toInt :=
    toInt_sub_of_range (by omega) (by omega)
  dsimp only
  omega

theorem forall_ring_conservation4_of_bv32 {N : Nat} (hN : 0 < N)
    (P : BitVec 4 → BitVec 4) (q4 : BitVec 4 → BitVec 4 → BitVec 32)
    (hb : ∀ u v : BitVec 4, BoundedBV32 10000000 (q4 u v))
    (h32 : ∀ a b c d : BitVec 4,
      q4 (step3 P a b c) (step3 P b c d) - q4 b c =
        canonCurrent4_32 P q4 b c d - canonCurrent4_32 P q4 a b c)
    (u : Fin N → BitVec 4) :
    ringSum N (fun i => (q4 (floquetRing hN P u i) (floquetRing hN P u (nextMod hN i))).toInt) =
      ringSum N (fun i => (q4 (u i) (u (nextMod hN i))).toInt) :=
  discrete_noether4_global_conservation hN P (fun x y => (q4 x y).toInt)
    (canonCurrent4 P (fun x y => (q4 x y).toInt))
    (bitvec32_continuity4_lifts_to_int P q4 hb h32) u

def IsTrivialQ4 (q4 : BitVec 4 → BitVec 4 → Int) : Prop :=
  ∃ (g f : BitVec 4 → Int) (h k : BitVec 2 → BitVec 2 → Int), ∀ u v : BitVec 4,
    q4 u v = (g v - g u) + f u + h (hi2 u) (hi2 v) + k (lo2 u) (lo2 v)

theorem non_trivial_q4_of_inner_cross {q4 : BitVec 4 → BitVec 4 → Int}
    (aL bR r0 r1 l0 l1 : BitVec 2)
    (hc :
      q4 (pack2 aL r1) (pack2 l1 bR) - q4 (pack2 aL r1) (pack2 l0 bR) -
        q4 (pack2 aL r0) (pack2 l1 bR) + q4 (pack2 aL r0) (pack2 l0 bR) ≠ 0) :
    ¬ IsTrivialQ4 q4 := by
  rintro ⟨g, f, h, k, hq⟩
  rw [hq (pack2 aL r1) (pack2 l1 bR), hq (pack2 aL r1) (pack2 l0 bR),
      hq (pack2 aL r0) (pack2 l1 bR), hq (pack2 aL r0) (pack2 l0 bR)] at hc
  simp only [hi2_pack2, lo2_pack2] at hc
  omega

theorem non_trivial_q4_of_outer_cross {q4 : BitVec 4 → BitVec 4 → Int}
    (aR bL l0 l1 r0 r1 : BitVec 2)
    (hc :
      q4 (pack2 l1 aR) (pack2 bL r1) - q4 (pack2 l1 aR) (pack2 bL r0) -
        q4 (pack2 l0 aR) (pack2 bL r1) + q4 (pack2 l0 aR) (pack2 bL r0) ≠ 0) :
    ¬ IsTrivialQ4 q4 := by
  rintro ⟨g, f, h, k, hq⟩
  rw [hq (pack2 l1 aR) (pack2 bL r1), hq (pack2 l1 aR) (pack2 bL r0),
      hq (pack2 l0 aR) (pack2 bL r1), hq (pack2 l0 aR) (pack2 bL r0)] at hc
  simp only [hi2_pack2, lo2_pack2] at hc
  omega

theorem active_current4_witness (P : BitVec 4 → BitVec 4)
    (q4 : BitVec 4 → BitVec 4 → Int) (a b c : BitVec 4)
    (h : canonCurrent4 P q4 a b c ≠ 0) :
    ∃ x y z, canonCurrent4 P q4 x y z ≠ 0 :=
  ⟨a, b, c, h⟩

/-! ## Section 9: $|X| = 8$ Non-Involutive Non-Degenerate Indecomposable YBE Theory -/

def hi3 (u : BitVec 6) : BitVec 3 :=
  BitVec.extractLsb' 3 3 u

def lo3 (u : BitVec 6) : BitVec 3 :=
  BitVec.extractLsb' 0 3 u

def pack3 (l r : BitVec 3) : BitVec 6 :=
  l ++ r

def swap3 (u : BitVec 6) : BitVec 6 :=
  pack3 (lo3 u) (hi3 u)

theorem pack3_hi3_lo3 (u : BitVec 6) : pack3 (hi3 u) (lo3 u) = u := by
  simp only [pack3, hi3, lo3]
  bv_decide

theorem hi3_pack3 (l r : BitVec 3) : hi3 (pack3 l r) = l := by
  simp only [pack3, hi3]
  bv_decide

theorem lo3_pack3 (l r : BitVec 3) : lo3 (pack3 l r) = r := by
  simp only [pack3, lo3]
  bv_decide

theorem swap3_involutive (u : BitVec 6) : swap3 (swap3 u) = u := by
  simp only [swap3, pack3, hi3, lo3]
  bv_decide

def yb8Apply12 (R : BitVec 6 → BitVec 6) (t : BitVec 9) : BitVec 9 :=
  let a := BitVec.extractLsb' 6 3 t
  let b := BitVec.extractLsb' 3 3 t
  let c := BitVec.extractLsb' 0 3 t
  R (pack3 a b) ++ c

def yb8Apply23 (R : BitVec 6 → BitVec 6) (t : BitVec 9) : BitVec 9 :=
  let a := BitVec.extractLsb' 6 3 t
  let b := BitVec.extractLsb' 3 3 t
  let c := BitVec.extractLsb' 0 3 t
  a ++ R (pack3 b c)

def IsSetTheoreticYB8 (R : BitVec 6 → BitVec 6) : Prop :=
  ∀ t : BitVec 9,
    yb8Apply12 R (yb8Apply23 R (yb8Apply12 R t)) =
      yb8Apply23 R (yb8Apply12 R (yb8Apply23 R t))

theorem yb8_id_is_solution : IsSetTheoreticYB8 id := by
  intro t
  simp only [yb8Apply12, yb8Apply23, pack3, id]
  bv_decide

theorem yb8_swap_is_solution : IsSetTheoreticYB8 swap3 := by
  intro t
  simp only [yb8Apply12, yb8Apply23, swap3, hi3, lo3, pack3]
  bv_decide

def IsNonDegenerate8 (R : BitVec 6 → BitVec 6) : Prop :=
  (∀ a : BitVec 3, Bijective (fun b => hi3 (R (pack3 a b)))) ∧
    (∀ b : BitVec 3, Bijective (fun a => lo3 (R (pack3 a b))))

theorem non_degenerate8_of_inverses (R : BitVec 6 → BitVec 6)
    (sigInv tauInv : BitVec 3 → BitVec 3 → BitVec 3)
    (hsL : ∀ a b, sigInv a (hi3 (R (pack3 a b))) = b)
    (hsR : ∀ a b, hi3 (R (pack3 a (sigInv a b))) = b)
    (htL : ∀ a b, tauInv b (lo3 (R (pack3 a b))) = a)
    (htR : ∀ a b, lo3 (R (pack3 (tauInv b a) b)) = a) :
    IsNonDegenerate8 R := by
  refine ⟨fun a => ⟨fun x y hxy => ?_, fun y => ⟨sigInv a y, hsR a y⟩⟩,
          fun b => ⟨fun x y hxy => ?_, fun y => ⟨tauInv b y, htR y b⟩⟩⟩
  · have h := congrArg (sigInv a) hxy
    dsimp only at h
    rw [hsL a x, hsL a y] at h
    exact h
  · have h := congrArg (tauInv b) hxy
    dsimp only at h
    rw [htL x b, htL y b] at h
    exact h

def IsInvolutive8 (R : BitVec 6 → BitVec 6) : Prop :=
  ∀ u : BitVec 6, R (R u) = u

theorem non_involutive8_of_witness {R : BitVec 6 → BitVec 6} (u : BitVec 6)
    (h : R (R u) ≠ u) : ¬ IsInvolutive8 R :=
  fun hi => h (hi u)

def IsLyubashenko8 (R : BitVec 6 → BitVec 6) : Prop :=
  ∃ s t : BitVec 3 → BitVec 3, ∀ a b : BitVec 3, R (pack3 a b) = pack3 (s b) (t a)

theorem non_lyubashenko8_of_hi_cross {R : BitVec 6 → BitVec 6}
    (a0 a1 b : BitVec 3) (h : hi3 (R (pack3 a1 b)) ≠ hi3 (R (pack3 a0 b))) :
    ¬ IsLyubashenko8 R := by
  rintro ⟨s, t, hst⟩
  rw [hst a1 b, hst a0 b, hi3_pack3, hi3_pack3] at h
  exact h rfl

theorem non_lyubashenko8_of_lo_cross {R : BitVec 6 → BitVec 6}
    (a b0 b1 : BitVec 3) (h : lo3 (R (pack3 a b1)) ≠ lo3 (R (pack3 a b0))) :
    ¬ IsLyubashenko8 R := by
  rintro ⟨s, t, hst⟩
  rw [hst a b1, hst a b0, lo3_pack3, lo3_pack3] at h
  exact h rfl

def IsAffineGF2_6 (f : BitVec 6 → Bool) : Prop :=
  ∃ c c0 c1 c2 c3 c4 c5 : Bool, ∀ x : BitVec 6,
    f x = Bool.xor c
      (Bool.xor (and c0 (BitVec.getLsbD x 0))
        (Bool.xor (and c1 (BitVec.getLsbD x 1))
          (Bool.xor (and c2 (BitVec.getLsbD x 2))
            (Bool.xor (and c3 (BitVec.getLsbD x 3))
              (Bool.xor (and c4 (BitVec.getLsbD x 4)) (and c5 (BitVec.getLsbD x 5)))))))

private theorem xor6_cancel (c c0 c1 c2 c3 c4 c5 a0 a1 a2 a3 a4 a5 b0 b1 b2 b3 b4 b5 : Bool) :
    Bool.xor
      (Bool.xor
        (Bool.xor c (Bool.xor (c0 && Bool.xor a0 b0)
          (Bool.xor (c1 && Bool.xor a1 b1)
            (Bool.xor (c2 && Bool.xor a2 b2)
              (Bool.xor (c3 && Bool.xor a3 b3)
                (Bool.xor (c4 && Bool.xor a4 b4) (c5 && Bool.xor a5 b5)))))))
        (Bool.xor c (Bool.xor (c0 && a0)
          (Bool.xor (c1 && a1)
            (Bool.xor (c2 && a2)
              (Bool.xor (c3 && a3) (Bool.xor (c4 && a4) (c5 && a5))))))))
      (Bool.xor
        (Bool.xor c (Bool.xor (c0 && b0)
          (Bool.xor (c1 && b1)
            (Bool.xor (c2 && b2)
              (Bool.xor (c3 && b3) (Bool.xor (c4 && b4) (c5 && b5)))))))
        c) = false := by
  bv_decide

theorem gf2_6_second_diff_witnesses_nonlinear {f : BitVec 6 → Bool} (a b : BitVec 6)
    (h : Bool.xor (Bool.xor (f (a ^^^ b)) (f a)) (Bool.xor (f b) (f 0#6)) = true) :
    ¬ IsAffineGF2_6 f := by
  rintro ⟨c, c0, c1, c2, c3, c4, c5, hf⟩
  have h0 : BitVec.getLsbD 0#6 0 = false ∧ BitVec.getLsbD 0#6 1 = false ∧
      BitVec.getLsbD 0#6 2 = false ∧ BitVec.getLsbD 0#6 3 = false ∧
      BitVec.getLsbD 0#6 4 = false ∧ BitVec.getLsbD 0#6 5 = false := by decide
  rw [hf (a ^^^ b), hf a, hf b, hf 0#6] at h
  simp only [BitVec.getLsbD_xor, h0.1, h0.2.1, h0.2.2.1, h0.2.2.2.1, h0.2.2.2.2.1, h0.2.2.2.2.2,
    Bool.and_false, Bool.xor_false] at h
  have hcancel := xor6_cancel c c0 c1 c2 c3 c4 c5
    (BitVec.getLsbD a 0) (BitVec.getLsbD a 1) (BitVec.getLsbD a 2)
    (BitVec.getLsbD a 3) (BitVec.getLsbD a 4) (BitVec.getLsbD a 5)
    (BitVec.getLsbD b 0) (BitVec.getLsbD b 1) (BitVec.getLsbD b 2)
    (BitVec.getLsbD b 3) (BitVec.getLsbD b 4) (BitVec.getLsbD b 5)
  rw [hcancel] at h
  exact Bool.false_ne_true h

def IsDecomposable8 (R : BitVec 6 → BitVec 6) : Prop :=
  ∃ C : BitVec 3 → Bool, (∃ x, C x = true) ∧ (∃ y, C y = false) ∧
    (∀ a b : BitVec 3, C (hi3 (R (pack3 a b))) = C b ∧ C (lo3 (R (pack3 a b))) = C a)

def anyBV3 (f : BitVec 3 → Bool) : Bool :=
  f 0#3 || f 1#3 || f 2#3 || f 3#3 || f 4#3 || f 5#3 || f 6#3 || f 7#3

private theorem anyBV3_exists {f : BitVec 3 → Bool} (h : anyBV3 f = true) :
    ∃ x : BitVec 3, f x = true := by
  unfold anyBV3 at h
  simp only [Bool.or_eq_true] at h
  rcases h with (((((((h0 | h1) | h2) | h3) | h4) | h5) | h6) | h7)
  · exact ⟨0#3, h0⟩
  · exact ⟨1#3, h1⟩
  · exact ⟨2#3, h2⟩
  · exact ⟨3#3, h3⟩
  · exact ⟨4#3, h4⟩
  · exact ⟨5#3, h5⟩
  · exact ⟨6#3, h6⟩
  · exact ⟨7#3, h7⟩

def adj8 (R : BitVec 6 → BitVec 6) (u v : BitVec 3) : Bool :=
  anyBV3 (fun x =>
    (hi3 (R (pack3 x u)) == v) || (hi3 (R (pack3 x v)) == u) ||
    (lo3 (R (pack3 u x)) == v) || (lo3 (R (pack3 v x)) == u))

def reach8Step (R : BitVec 6 → BitVec 6) (reached : BitVec 3 → Bool) (v : BitVec 3) : Bool :=
  reached v || anyBV3 (fun u => reached u && adj8 R u v)

def reach8K (R : BitVec 6 → BitVec 6) : Nat → BitVec 3 → Bool
  | 0 => fun x => x == 0#3
  | k + 1 => reach8Step R (reach8K R k)

def reach8 (R : BitVec 6 → BitVec 6) : BitVec 3 → Bool :=
  reach8K R 7

private theorem adj8_preserves_const (R : BitVec 6 → BitVec 6) (C : BitVec 3 → Bool)
    (hInv : ∀ a b : BitVec 3, C (hi3 (R (pack3 a b))) = C b ∧ C (lo3 (R (pack3 a b))) = C a)
    {u v : BitVec 3} (hadj : adj8 R u v = true) : C v = C u := by
  obtain ⟨x, hx⟩ := anyBV3_exists hadj
  simp only [Bool.or_eq_true, beq_iff_eq] at hx
  rcases hx with ((h1 | h2) | h3) | h4
  · rw [← h1, (hInv x u).1]
  · have hm := (hInv x v).1
    rw [h2] at hm
    exact hm.symm
  · rw [← h3, (hInv u x).2]
  · have hm := (hInv v x).2
    rw [h4] at hm
    exact hm.symm

private theorem reach8K_const (R : BitVec 6 → BitVec 6) (C : BitVec 3 → Bool)
    (hInv : ∀ a b : BitVec 3, C (hi3 (R (pack3 a b))) = C b ∧ C (lo3 (R (pack3 a b))) = C a) :
    ∀ k v, reach8K R k v = true → C v = C 0#3
  | 0, v, hv => by
    simp only [reach8K, beq_iff_eq] at hv
    rw [hv]
  | k + 1, v, hv => by
    simp only [reach8K, reach8Step, Bool.or_eq_true] at hv
    rcases hv with hprev | hstep
    · exact reach8K_const R C hInv k v hprev
    · obtain ⟨u, hu⟩ := anyBV3_exists hstep
      simp only [Bool.and_eq_true] at hu
      rw [adj8_preserves_const R C hInv hu.2, reach8K_const R C hInv k u hu.1]

theorem indecomposable8_of_connected8 (R : BitVec 6 → BitVec 6)
    (hConn : ∀ x : BitVec 3, reach8 R x = true) : ¬ IsDecomposable8 R := by
  rintro ⟨C, ⟨x, hx⟩, ⟨y, hy⟩, hInv⟩
  have hcx := reach8K_const R C hInv 7 x (hConn x)
  have hcy := reach8K_const R C hInv 7 y (hConn y)
  rw [hx] at hcx
  rw [hy] at hcy
  have htf : true = false := hcx.trans hcy.symm
  exact Bool.false_ne_true htf.symm

/-! ## Section 10: Dual-Unitary Equivalence, Universal $Q_2$ No-Go Theorem, Mixing $M_\pm$, & $\text{AME}(4,8)$ -/

def IsNonDegenerate4 (P : BitVec 4 → BitVec 4) : Prop :=
  (∀ a : BitVec 2, Bijective (fun b => hi2 (P (pack2 a b)))) ∧
    (∀ b : BitVec 2, Bijective (fun a => lo2 (P (pack2 a b))))

def IsDualBijective4 (P : BitVec 4 → BitVec 4) : Prop :=
  Bijective (fun u : BitVec 4 => pack2 (hi2 u) (hi2 (P u))) ∧
    Bijective (fun u : BitVec 4 => pack2 (lo2 u) (lo2 (P u)))

theorem dual_bijective4_iff_nondegenerate4 (P : BitVec 4 → BitVec 4) :
    IsDualBijective4 P ↔ IsNonDegenerate4 P := by
  constructor
  · rintro ⟨hL, hR⟩
    refine ⟨fun a => ⟨fun b1 b2 h => ?_, fun c => ?_⟩,
            fun b => ⟨fun a1 a2 h => ?_, fun d => ?_⟩⟩
    · have hphi : pack2 (hi2 (pack2 a b1)) (hi2 (P (pack2 a b1))) =
          pack2 (hi2 (pack2 a b2)) (hi2 (P (pack2 a b2))) := by
        simp only [hi2_pack2, h]
      have heq := hL.1 _ _ hphi
      have hlo := congrArg lo2 heq
      simpa only [lo2_pack2] using hlo
    · obtain ⟨u, hu⟩ := hL.2 (pack2 a c)
      have hha : hi2 u = a := by
        have ht := congrArg hi2 hu
        simpa only [hi2_pack2] using ht
      have hhc : hi2 (P u) = c := by
        have ht := congrArg lo2 hu
        simpa only [lo2_pack2] using ht
      refine ⟨lo2 u, ?_⟩
      dsimp only
      have hpack : pack2 a (lo2 u) = u := by rw [← hha, pack2_hi2_lo2]
      rw [hpack, hhc]
    · have hpsi : pack2 (lo2 (pack2 a1 b)) (lo2 (P (pack2 a1 b))) =
          pack2 (lo2 (pack2 a2 b)) (lo2 (P (pack2 a2 b))) := by
        simp only [lo2_pack2, h]
      have heq := hR.1 _ _ hpsi
      have hhi := congrArg hi2 heq
      simpa only [hi2_pack2] using hhi
    · obtain ⟨u, hu⟩ := hR.2 (pack2 b d)
      have hhb : lo2 u = b := by
        have ht := congrArg hi2 hu
        simpa only [hi2_pack2] using ht
      have hhd : lo2 (P u) = d := by
        have ht := congrArg lo2 hu
        simpa only [lo2_pack2] using ht
      refine ⟨hi2 u, ?_⟩
      dsimp only
      have hpack : pack2 (hi2 u) b = u := by rw [← hhb, pack2_hi2_lo2]
      rw [hpack, hhd]
  · rintro ⟨hL, hR⟩
    refine ⟨⟨fun u v huv => ?_, fun w => ?_⟩, ⟨fun u v huv => ?_, fun w => ?_⟩⟩
    · have h1 : hi2 u = hi2 v := by
        have ht := congrArg hi2 huv
        simpa only [hi2_pack2] using ht
      have h2 : hi2 (P u) = hi2 (P v) := by
        have ht := congrArg lo2 huv
        simpa only [lo2_pack2] using ht
      have hpu : hi2 (P (pack2 (hi2 u) (lo2 u))) = hi2 (P (pack2 (hi2 u) (lo2 v))) := by
        rw [pack2_hi2_lo2, h1, pack2_hi2_lo2]
        exact h2
      have hlo := (hL (hi2 u)).1 _ _ hpu
      rw [← pack2_hi2_lo2 u, ← pack2_hi2_lo2 v, h1, hlo]
    · obtain ⟨b, hb⟩ := (hL (hi2 w)).2 (lo2 w)
      refine ⟨pack2 (hi2 w) b, ?_⟩
      simp only [hi2_pack2, hb, pack2_hi2_lo2]
    · have h1 : lo2 u = lo2 v := by
        have ht := congrArg hi2 huv
        simpa only [hi2_pack2] using ht
      have h2 : lo2 (P u) = lo2 (P v) := by
        have ht := congrArg lo2 huv
        simpa only [lo2_pack2] using ht
      have hpu : lo2 (P (pack2 (hi2 u) (lo2 u))) = lo2 (P (pack2 (hi2 v) (lo2 u))) := by
        rw [pack2_hi2_lo2, h1, pack2_hi2_lo2]
        exact h2
      have hhi := (hR (lo2 u)).1 _ _ hpu
      rw [← pack2_hi2_lo2 u, ← pack2_hi2_lo2 v, h1, hhi]
    · obtain ⟨a', ha'⟩ := (hR (hi2 w)).2 (lo2 w)
      refine ⟨pack2 a' (hi2 w), ?_⟩
      simp only [lo2_pack2, ha', pack2_hi2_lo2]

theorem non_degenerate4_of_inverses (P : BitVec 4 → BitVec 4)
    (sigInv4 tauInv4 : BitVec 2 → BitVec 2 → BitVec 2)
    (hsL : ∀ a b, sigInv4 a (hi2 (P (pack2 a b))) = b)
    (hsR : ∀ a b, hi2 (P (pack2 a (sigInv4 a b))) = b)
    (htL : ∀ a b, tauInv4 b (lo2 (P (pack2 a b))) = a)
    (htR : ∀ a b, lo2 (P (pack2 (tauInv4 b a) b)) = a) :
    IsNonDegenerate4 P := by
  refine ⟨fun a => ⟨fun x y hxy => ?_, fun y => ⟨sigInv4 a y, hsR a y⟩⟩,
          fun b => ⟨fun x y hxy => ?_, fun y => ⟨tauInv4 b y, htR y b⟩⟩⟩
  · have h := congrArg (sigInv4 a) hxy
    dsimp only at h
    rw [hsL a x, hsL a y] at h
    exact h
  · have h := congrArg (tauInv4 b) hxy
    dsimp only at h
    rw [htL x b, htL y b] at h
    exact h

def IsDualUnitary4 (P : BitVec 4 → BitVec 4) : Prop :=
  Bijective P ∧ IsNonDegenerate4 P

theorem dual_unitary_q2_separable_nogo
    (P : BitVec 4 → BitVec 4) (hDU : IsDualUnitary4 P)
    (q : BitVec 4 → Int) (j : BitVec 4 → BitVec 4 → Int)
    (hCont : SatisfiesContinuity P q j) : IsSeparable q := by
  let bL := hi2 (P 0#4)
  let bR := lo2 (P 0#4)
  let uOf : BitVec 2 → BitVec 2 := fun x => Classical.choose ((hDU.2.2 bL).2 x)
  have huOf : ∀ x, lo2 (P (pack2 (uOf x) bL)) = x := fun x => Classical.choose_spec ((hDU.2.2 bL).2 x)
  let axOf : BitVec 2 → BitVec 4 := fun x => Classical.choose (hDU.1.2 (pack2 0#2 (uOf x)))
  have haxOf : ∀ x, P (axOf x) = pack2 0#2 (uOf x) :=
    fun x => Classical.choose_spec (hDU.1.2 (pack2 0#2 (uOf x)))
  let vOf : BitVec 2 → BitVec 2 := fun y => Classical.choose ((hDU.2.1 bR).2 y)
  have hvOf : ∀ y, hi2 (P (pack2 bR (vOf y))) = y := fun y => Classical.choose_spec ((hDU.2.1 bR).2 y)
  let cyOf : BitVec 2 → BitVec 4 := fun y => Classical.choose (hDU.1.2 (pack2 (vOf y) 0#2))
  have hcyOf : ∀ y, P (cyOf y) = pack2 (vOf y) 0#2 :=
    fun y => Classical.choose_spec (hDU.1.2 (pack2 (vOf y) 0#2))
  refine ⟨fun x => q 0#4 - j (axOf x) 0#4, fun y => j 0#4 (cyOf y), fun x y => ?_⟩
  have hstep : step3 P (axOf x) 0#4 (cyOf y) = pack2 x y := by
    simp only [step3, haxOf, hcyOf, lo2_pack2, hi2_pack2, bL, bR, huOf, hvOf]
  have hc := hCont (axOf x) 0#4 (cyOf y)
  rw [hstep] at hc
  dsimp only
  omega

theorem ybe_nondeg_q2_no_entangled
    (P : BitVec 4 → BitVec 4) (hDU : IsDualUnitary4 P) (_hYB : IsSetTheoreticYB P) :
    ¬ ∃ q j, ¬ IsSeparable q ∧ SatisfiesContinuity P q j := by
  rintro ⟨q, j, hnq, hCont⟩
  exact hnq (dual_unitary_q2_separable_nogo P hDU q j hCont)

def countBV2 (p : BitVec 2 → Bool) : Nat :=
  (if p 0#2 then 1 else 0) + (if p 1#2 then 1 else 0) +
    (if p 2#2 then 1 else 0) + (if p 3#2 then 1 else 0)

def allBV2 (p : BitVec 2 → Bool) : Bool :=
  p 0#2 && p 1#2 && p 2#2 && p 3#2

def lcPlus4 (P : BitVec 4 → BitVec 4) (a d : BitVec 2) : Nat :=
  countBV2 (fun b => lo2 (P (pack2 a b)) == d)

def lcMinus4 (P : BitVec 4 → BitVec 4) (b c : BitVec 2) : Nat :=
  countBV2 (fun a => hi2 (P (pack2 a b)) == c)

def mat4Mul (A B : BitVec 2 → BitVec 2 → Nat) (i j : BitVec 2) : Nat :=
  A i 0#2 * B 0#2 j + A i 1#2 * B 1#2 j + A i 2#2 * B 2#2 j + A i 3#2 * B 3#2 j

def mat4Pow4AllPos (M : BitVec 2 → BitVec 2 → Nat) : Bool :=
  let m2 := mat4Mul M M
  let m4 := mat4Mul m2 m2
  allBV2 (fun i => allBV2 (fun j => decide (0 < m4 i j)))

def IsMixingDualUnitary4 (P : BitVec 4 → BitVec 4) : Prop :=
  IsDualUnitary4 P ∧ mat4Pow4AllPos (lcPlus4 P) = true ∧ mat4Pow4AllPos (lcMinus4 P) = true

def IsCrossBijective8 (R : BitVec 6 → BitVec 6) : Prop :=
  (∀ a : BitVec 3, Bijective (fun b => lo3 (R (pack3 a b)))) ∧
    (∀ b : BitVec 3, Bijective (fun a => hi3 (R (pack3 a b))))

theorem cross_bijective8_of_inverses (R : BitVec 6 → BitVec 6)
    (crossLInv crossRInv : BitVec 3 → BitVec 3 → BitVec 3)
    (hcL1 : ∀ a b, crossLInv a (lo3 (R (pack3 a b))) = b)
    (hcL2 : ∀ a b, lo3 (R (pack3 a (crossLInv a b))) = b)
    (hcR1 : ∀ a b, crossRInv b (hi3 (R (pack3 a b))) = a)
    (hcR2 : ∀ a b, hi3 (R (pack3 (crossRInv b a) b)) = a) :
    IsCrossBijective8 R := by
  refine ⟨fun a => ⟨fun x y hxy => ?_, fun y => ⟨crossLInv a y, hcL2 a y⟩⟩,
          fun b => ⟨fun x y hxy => ?_, fun y => ⟨crossRInv b y, hcR2 y b⟩⟩⟩
  · have h := congrArg (crossLInv a) hxy
    dsimp only at h
    rw [hcL1 a x, hcL1 a y] at h
    exact h
  · have h := congrArg (crossRInv b) hxy
    dsimp only at h
    rw [hcR1 x b, hcR1 y b] at h
    exact h

def Is2Unitary8 (R : BitVec 6 → BitVec 6) : Prop :=
  Bijective R ∧ IsNonDegenerate8 R ∧ IsCrossBijective8 R

/-! ## Section 11: 3-Site Fracton Dipole-Laplacian Noether Theorem & Single-Charge Immobility -/

def site0_3 (u : BitVec 6) : BitVec 2 :=
  BitVec.extractLsb' 4 2 u

def site1_3 (u : BitVec 6) : BitVec 2 :=
  BitVec.extractLsb' 2 2 u

def site2_3 (u : BitVec 6) : BitVec 2 :=
  BitVec.extractLsb' 0 2 u

def pack3_2 (a b c : BitVec 2) : BitVec 6 :=
  a ++ b ++ c

theorem pack3_2_sites (u : BitVec 6) : pack3_2 (site0_3 u) (site1_3 u) (site2_3 u) = u := by
  simp only [pack3_2, site0_3, site1_3, site2_3]
  bv_decide

theorem site0_3_pack3_2 (a b c : BitVec 2) : site0_3 (pack3_2 a b c) = a := by
  simp only [pack3_2, site0_3]
  bv_decide

theorem site1_3_pack3_2 (a b c : BitVec 2) : site1_3 (pack3_2 a b c) = b := by
  simp only [pack3_2, site1_3]
  bv_decide

theorem site2_3_pack3_2 (a b c : BitVec 2) : site2_3 (pack3_2 a b c) = c := by
  simp only [pack3_2, site2_3]
  bv_decide

def rho2 (x : BitVec 2) : Int :=
  (x.toNat : Int)

def SatisfiesDipoleLaplacian (T : BitVec 6 → BitVec 6) (rho : BitVec 2 → Int) : Prop :=
  ∀ a b c : BitVec 2,
    let out := T (pack3_2 a b c)
    let phi := rho (site0_3 out) - rho a
    rho (site1_3 out) - rho b = -(2 * phi) ∧ rho (site2_3 out) - rho c = phi

theorem dipole_laplacian_local_QD (T : BitVec 6 → BitVec 6) (rho : BitVec 2 → Int)
    (hL : SatisfiesDipoleLaplacian T rho) (a b c : BitVec 2) (k : Int) :
    let out := T (pack3_2 a b c)
    let a' := site0_3 out
    let b' := site1_3 out
    let c' := site2_3 out
    (rho a' + rho b' + rho c' = rho a + rho b + rho c) ∧
      (k * rho a' + (k + 1) * rho b' + (k + 2) * rho c' =
        k * rho a + (k + 1) * rho b + (k + 2) * rho c) := by
  have hL_abc := hL a b c
  dsimp only at hL_abc ⊢
  rcases hL_abc with ⟨h1, h2⟩
  have hQ : rho (site0_3 (T (pack3_2 a b c))) + rho (site1_3 (T (pack3_2 a b c))) +
      rho (site2_3 (T (pack3_2 a b c))) = rho a + rho b + rho c := by
    omega
  refine ⟨hQ, ?_⟩
  have hmul : k * rho (site0_3 (T (pack3_2 a b c))) + k * rho (site1_3 (T (pack3_2 a b c))) +
      k * rho (site2_3 (T (pack3_2 a b c))) =
      k * rho a + k * rho b + k * rho c := by
    calc k * rho (site0_3 (T (pack3_2 a b c))) + k * rho (site1_3 (T (pack3_2 a b c))) +
          k * rho (site2_3 (T (pack3_2 a b c)))
        = k * (rho (site0_3 (T (pack3_2 a b c))) + rho (site1_3 (T (pack3_2 a b c))) +
            rho (site2_3 (T (pack3_2 a b c)))) := by
          rw [Int.mul_add, Int.mul_add]
      _ = k * (rho a + rho b + rho c) := by rw [hQ]
      _ = k * rho a + k * rho b + k * rho c := by rw [Int.mul_add, Int.mul_add]
  have hr1 : (k + 1) * rho (site1_3 (T (pack3_2 a b c))) =
      k * rho (site1_3 (T (pack3_2 a b c))) + rho (site1_3 (T (pack3_2 a b c))) := by
    rw [Int.add_mul, Int.one_mul]
  have hr2 : (k + 2) * rho (site2_3 (T (pack3_2 a b c))) =
      k * rho (site2_3 (T (pack3_2 a b c))) + 2 * rho (site2_3 (T (pack3_2 a b c))) := by
    rw [Int.add_mul]
  have hl1 : (k + 1) * rho b = k * rho b + rho b := by rw [Int.add_mul, Int.one_mul]
  have hl2 : (k + 2) * rho c = k * rho c + 2 * rho c := by rw [Int.add_mul]
  rw [hr1, hr2, hl1, hl2]
  omega

def applyGate3 {L : Nat} (T : BitVec 6 → BitVec 6) (pos : Nat) (hpos : pos + 2 < L)
    (s : Fin L → BitVec 2) : Fin L → BitVec 2 :=
  let out := T (pack3_2 (s ⟨pos, by omega⟩) (s ⟨pos + 1, by omega⟩) (s ⟨pos + 2, by omega⟩))
  fun i =>
    if i.val = pos then site0_3 out
    else if i.val = pos + 1 then site1_3 out
    else if i.val = pos + 2 then site2_3 out
    else s i

def openCharge (L : Nat) (rho : BitVec 2 → Int) (s : Fin L → BitVec 2) : Int :=
  ringSum L (fun i => rho (s i))

def openDipole (L : Nat) (rho : BitVec 2 → Int) (s : Fin L → BitVec 2) : Int :=
  ringSum L (fun i => (i.val : Int) * rho (s i))

def pointDelta (L idx : Nat) (d : Int) (i : Fin L) : Int :=
  if i.val = idx then d else 0

private theorem sumPrefix_pointDelta (L idx : Nat) (d : Int) (hidx : idx < L) :
    ∀ k : Nat, k ≤ L →
      sumPrefix (pointDelta L idx d) k = if idx < k then d else 0
  | 0, _ => by simp [Flow.FiniteSum.sumPrefix]
  | k + 1, hk => by
    have hkL : k < L := by omega
    have ih := sumPrefix_pointDelta L idx d hidx k (by omega)
    simp only [Flow.FiniteSum.sumPrefix, hkL, ↓reduceDIte, pointDelta, ih]
    by_cases heq : k = idx
    · subst heq
      have h1 : ¬ k < k := Nat.lt_irrefl k
      have h2 : k < k + 1 := Nat.lt_succ_self k
      simp [h1, h2]
    · have h1 : (idx < k) ↔ (idx < k + 1) := by omega
      simp [heq, h1]

theorem ringSum_pointDelta {L idx : Nat} (hidx : idx < L) (d : Int) :
    ringSum L (pointDelta L idx d) = d := by
  change sumPrefix (pointDelta L idx d) L = d
  rw [sumPrefix_pointDelta L idx d hidx L (Nat.le_refl L)]
  simp [hidx]

theorem fracton_open_chain_Q_conservation {L : Nat} (T : BitVec 6 → BitVec 6)
    (rho : BitVec 2 → Int) (hL : SatisfiesDipoleLaplacian T rho)
    (pos : Nat) (hpos : pos + 2 < L) (s : Fin L → BitVec 2) :
    openCharge L rho (applyGate3 T pos hpos s) = openCharge L rho s := by
  let a := s ⟨pos, by omega⟩
  let b := s ⟨pos + 1, by omega⟩
  let c := s ⟨pos + 2, by omega⟩
  let out := T (pack3_2 a b c)
  let d0 := rho (site0_3 out) - rho a
  let d1 := rho (site1_3 out) - rho b
  let d2 := rho (site2_3 out) - rho c
  have hfn : (fun i => rho (applyGate3 T pos hpos s i) - rho (s i)) =
      (fun i => pointDelta L pos d0 i + pointDelta L (pos + 1) d1 i +
        pointDelta L (pos + 2) d2 i) := by
    funext i
    simp only [applyGate3, pointDelta]
    by_cases h0 : i.val = pos
    · have h01 : pos ≠ pos + 1 := by omega
      have h02 : pos ≠ pos + 2 := by omega
      have hidx : i = ⟨pos, by omega⟩ := Fin.ext h0
      simp only [h0, h01, h02, ↓reduceIte]
      rw [hidx]
      simp [d0, out, a, b, c]
    · by_cases h1 : i.val = pos + 1
      · have h10 : pos + 1 ≠ pos := by omega
        have h12 : pos + 1 ≠ pos + 2 := by omega
        have hidx : i = ⟨pos + 1, by omega⟩ := Fin.ext h1
        simp only [h1, h10, h12, ↓reduceIte]
        rw [hidx]
        simp [d1, out, a, b, c]
      · by_cases h2 : i.val = pos + 2
        · have h20 : pos + 2 ≠ pos := by omega
          have h21 : pos + 2 ≠ pos + 1 := by omega
          have hidx : i = ⟨pos + 2, by omega⟩ := Fin.ext h2
          simp only [h2, h20, h21, ↓reduceIte]
          rw [hidx]
          simp [d2, out, a, b, c]
        · simp [h0, h1, h2]
  have hsub : openCharge L rho (applyGate3 T pos hpos s) - openCharge L rho s =
      d0 + d1 + d2 := by
    unfold openCharge
    rw [← ringSum_sub, hfn, ringSum_add, ringSum_add,
        ringSum_pointDelta (by omega : pos < L) d0,
        ringSum_pointDelta (by omega : pos + 1 < L) d1,
        ringSum_pointDelta (by omega : pos + 2 < L) d2]
  have hloc := (dipole_laplacian_local_QD T rho hL a b c 0).1
  dsimp only [d0, d1, d2, out] at hsub hloc
  omega

theorem fracton_open_chain_D_conservation {L : Nat} (T : BitVec 6 → BitVec 6)
    (rho : BitVec 2 → Int) (hL : SatisfiesDipoleLaplacian T rho)
    (pos : Nat) (hpos : pos + 2 < L) (s : Fin L → BitVec 2) :
    openDipole L rho (applyGate3 T pos hpos s) = openDipole L rho s := by
  let a := s ⟨pos, by omega⟩
  let b := s ⟨pos + 1, by omega⟩
  let c := s ⟨pos + 2, by omega⟩
  let out := T (pack3_2 a b c)
  let e0 := (pos : Int) * rho (site0_3 out) - (pos : Int) * rho a
  let e1 := ((pos : Int) + 1) * rho (site1_3 out) - ((pos : Int) + 1) * rho b
  let e2 := ((pos : Int) + 2) * rho (site2_3 out) - ((pos : Int) + 2) * rho c
  have hfn : (fun i : Fin L => (i.val : Int) * rho (applyGate3 T pos hpos s i) -
        (i.val : Int) * rho (s i)) =
      (fun i => pointDelta L pos e0 i + pointDelta L (pos + 1) e1 i +
        pointDelta L (pos + 2) e2 i) := by
    funext i
    simp only [applyGate3, pointDelta]
    by_cases h0 : i.val = pos
    · have h01 : pos ≠ pos + 1 := by omega
      have h02 : pos ≠ pos + 2 := by omega
      have hidx : i = ⟨pos, by omega⟩ := Fin.ext h0
      simp only [h0, h01, h02, ↓reduceIte]
      rw [hidx]
      simp [e0, out, a, b, c]
    · by_cases h1 : i.val = pos + 1
      · have h10 : pos + 1 ≠ pos := by omega
        have h12 : pos + 1 ≠ pos + 2 := by omega
        have hidx : i = ⟨pos + 1, by omega⟩ := Fin.ext h1
        simp only [h1, h10, h12, ↓reduceIte]
        rw [hidx]
        simp [e1, out, a, b, c]
      · by_cases h2 : i.val = pos + 2
        · have h20 : pos + 2 ≠ pos := by omega
          have h21 : pos + 2 ≠ pos + 1 := by omega
          have hidx : i = ⟨pos + 2, by omega⟩ := Fin.ext h2
          simp only [h2, h20, h21, ↓reduceIte]
          rw [hidx]
          simp [e2, out, a, b, c]
        · simp [h0, h1, h2]
  have hsub : openDipole L rho (applyGate3 T pos hpos s) - openDipole L rho s =
      e0 + e1 + e2 := by
    unfold openDipole
    rw [← ringSum_sub, hfn, ringSum_add, ringSum_add,
        ringSum_pointDelta (by omega : pos < L) e0,
        ringSum_pointDelta (by omega : pos + 1 < L) e1,
        ringSum_pointDelta (by omega : pos + 2 < L) e2]
  have hloc := (dipole_laplacian_local_QD T rho hL a b c (pos : Int)).2
  dsimp only [e0, e1, e2, out] at hsub hloc
  omega

theorem dipole_laplacian_low_charge_fixed (T : BitVec 6 → BitVec 6)
    (hL : SatisfiesDipoleLaplacian T rho2) (a b c : BitVec 2)
    (hle : a.toNat + b.toNat + c.toNat ≤ 1) :
    T (pack3_2 a b c) = pack3_2 a b c := by
  have hL_abc := hL a b c
  dsimp only [rho2] at hL_abc
  rcases hL_abc with ⟨h1, h2⟩
  have ha : a.toNat < 4 := a.isLt
  have hb : b.toNat < 4 := b.isLt
  have hc : c.toNat < 4 := c.isLt
  have ha' : (site0_3 (T (pack3_2 a b c))).toNat < 4 := (site0_3 (T (pack3_2 a b c))).isLt
  have hb' : (site1_3 (T (pack3_2 a b c))).toNat < 4 := (site1_3 (T (pack3_2 a b c))).isLt
  have hc' : (site2_3 (T (pack3_2 a b c))).toNat < 4 := (site2_3 (T (pack3_2 a b c))).isLt
  have h0eq : (site0_3 (T (pack3_2 a b c))).toNat = a.toNat := by omega
  have h1eq : (site1_3 (T (pack3_2 a b c))).toNat = b.toNat := by omega
  have h2eq : (site2_3 (T (pack3_2 a b c))).toNat = c.toNat := by omega
  have hs0 : site0_3 (T (pack3_2 a b c)) = a := BitVec.eq_of_toNat_eq h0eq
  have hs1 : site1_3 (T (pack3_2 a b c)) = b := BitVec.eq_of_toNat_eq h1eq
  have hs2 : site2_3 (T (pack3_2 a b c)) = c := BitVec.eq_of_toNat_eq h2eq
  rw [← pack3_2_sites (T (pack3_2 a b c)), hs0, hs1, hs2]

private theorem sumPrefix_nonneg {L : Nat} (f : Fin L → Int) (hf : ∀ i, 0 ≤ f i) :
    ∀ k : Nat, 0 ≤ sumPrefix f k
  | 0 => by simp [Flow.FiniteSum.sumPrefix]
  | k + 1 => by
    have ih := sumPrefix_nonneg f hf k
    dsimp only [sumPrefix] at ih ⊢
    simp only [Flow.FiniteSum.sumPrefix]
    split
    · have := hf ⟨k, by assumption⟩
      omega
    · omega

private theorem sumPrefix_mono {L : Nat} (f : Fin L → Int) (hf : ∀ i, 0 ≤ f i)
    {m n : Nat} (hmn : m ≤ n) : sumPrefix f m ≤ sumPrefix f n := by
  induction n with
  | zero =>
    have hm : m = 0 := by omega
    subst hm
    exact Int.le_refl _
  | succ n ih =>
    by_cases hmeq : m = n + 1
    · subst hmeq
      exact Int.le_refl _
    · have hmlen : m ≤ n := by omega
      have h1 := ih hmlen
      dsimp only [sumPrefix] at h1 ⊢
      simp only [Flow.FiniteSum.sumPrefix]
      split
      · have := hf ⟨n, by assumption⟩
        omega
      · omega

theorem fracton_single_charge_frozen {L : Nat} (T : BitVec 6 → BitVec 6)
    (hL : SatisfiesDipoleLaplacian T rho2)
    (pos : Nat) (hpos : pos + 2 < L) (s : Fin L → BitVec 2)
    (hQ1 : openCharge L rho2 s ≤ 1) :
    applyGate3 T pos hpos s = s := by
  have hf : ∀ i : Fin L, 0 ≤ rho2 (s i) := fun i => by
    simp only [rho2]
    omega
  have hmono : sumPrefix (fun i => rho2 (s i)) (pos + 3) ≤ openCharge L rho2 s :=
    sumPrefix_mono _ hf (by omega)
  have hnonneg : 0 ≤ sumPrefix (fun i => rho2 (s i)) pos := sumPrefix_nonneg _ hf pos
  have hstep : sumPrefix (fun i => rho2 (s i)) (pos + 3) =
      sumPrefix (fun i => rho2 (s i)) pos +
        rho2 (s ⟨pos, by omega⟩) + rho2 (s ⟨pos + 1, by omega⟩) + rho2 (s ⟨pos + 2, by omega⟩) := by
    simp only [Flow.FiniteSum.sumPrefix, (by omega : pos < L), (by omega : pos + 1 < L),
               (by omega : pos + 2 < L), ↓reduceDIte]
  have hwin : (s ⟨pos, by omega⟩).toNat + (s ⟨pos + 1, by omega⟩).toNat +
      (s ⟨pos + 2, by omega⟩).toNat ≤ 1 := by
    simp only [rho2] at hstep hmono hnonneg
    omega
  have hfix := dipole_laplacian_low_charge_fixed T hL
    (s ⟨pos, by omega⟩) (s ⟨pos + 1, by omega⟩) (s ⟨pos + 2, by omega⟩) hwin
  funext i
  simp only [applyGate3, hfix, site0_3_pack3_2, site1_3_pack3_2, site2_3_pack3_2]
  by_cases h0 : i.val = pos
  · have hidx : (⟨pos, by omega⟩ : Fin L) = i := Fin.ext h0.symm
    simp only [h0, ↓reduceIte, hidx]
  · by_cases h1 : i.val = pos + 1
    · have h10 : pos + 1 ≠ pos := by omega
      have hidx : (⟨pos + 1, by omega⟩ : Fin L) = i := Fin.ext h1.symm
      simp only [h1, h10, ↓reduceIte, hidx]
    · by_cases h2 : i.val = pos + 2
      · have h20 : pos + 2 ≠ pos := by omega
        have h21 : pos + 2 ≠ pos + 1 := by omega
        have hidx : (⟨pos + 2, by omega⟩ : Fin L) = i := Fin.ext h2.symm
        simp only [h2, h20, h21, ↓reduceIte, hidx]
      · simp [h0, h1, h2]

theorem fracton_quadrupole_active_witness (T : BitVec 6 → BitVec 6) (rho : BitVec 2 → Int)
    (a b c : BitVec 2) (h : rho (site0_3 (T (pack3_2 a b c))) ≠ rho a) :
    ∃ x y z : BitVec 2, rho (site0_3 (T (pack3_2 x y z))) ≠ rho x :=
  ⟨a, b, c, h⟩

/-! ## Section 12: Irretractable $|X|=8$ R1+R2+R3 Biquandles & Topological Knot Invariants -/

def IsR1Biquandle8 (R : BitVec 6 → BitVec 6) : Prop :=
  ∀ x : BitVec 3, R (pack3 x x) = pack3 x x

def IsRetractable8 (R : BitVec 6 → BitVec 6) : Prop :=
  ∃ x y : BitVec 3, x ≠ y ∧
    (∀ b : BitVec 3, hi3 (R (pack3 x b)) = hi3 (R (pack3 y b))) ∧
    (∀ a : BitVec 3, lo3 (R (pack3 a x)) = lo3 (R (pack3 a y)))

def allBV3 (f : BitVec 3 → Bool) : Bool :=
  f 0#3 && f 1#3 && f 2#3 && f 3#3 && f 4#3 && f 5#3 && f 6#3 && f 7#3

private theorem bv3_cases (x : BitVec 3) :
    x = 0#3 ∨ x = 1#3 ∨ x = 2#3 ∨ x = 3#3 ∨ x = 4#3 ∨ x = 5#3 ∨ x = 6#3 ∨ x = 7#3 := by
  bv_decide

private theorem allBV3_forall {f : BitVec 3 → Bool} (h : allBV3 f = true) :
    ∀ x : BitVec 3, f x = true := by
  unfold allBV3 at h
  simp only [Bool.and_eq_true] at h
  rcases h with ⟨⟨⟨⟨⟨⟨⟨h0, h1⟩, h2⟩, h3⟩, h4⟩, h5⟩, h6⟩, h7⟩
  intro x
  rcases bv3_cases x with h | h | h | h | h | h | h | h <;> rw [h] <;> assumption

private theorem allBV3_of_forall {f : BitVec 3 → Bool} (h : ∀ x : BitVec 3, f x = true) :
    allBV3 f = true := by
  unfold allBV3
  simp [h 0#3, h 1#3, h 2#3, h 3#3, h 4#3, h 5#3, h 6#3, h 7#3]

def pairSameSigTau8 (R : BitVec 6 → BitVec 6) (x y : BitVec 3) : Bool :=
  allBV3 (fun b => hi3 (R (pack3 x b)) == hi3 (R (pack3 y b))) &&
    allBV3 (fun a => lo3 (R (pack3 a x)) == lo3 (R (pack3 a y)))

def irretractable8Bool (R : BitVec 6 → BitVec 6) : Bool :=
  allBV3 (fun x => allBV3 (fun y => (x == y) || !(pairSameSigTau8 R x y)))

theorem irretractable8_of_bool (R : BitVec 6 → BitVec 6)
    (h : irretractable8Bool R = true) : ¬ IsRetractable8 R := by
  rintro ⟨x, y, hne, hsig, htau⟩
  have hxy : ((x == y) || !(pairSameSigTau8 R x y)) = true :=
    allBV3_forall (allBV3_forall h x) y
  have hneB : (x == y) = false := beq_false_of_ne hne
  have hpair : pairSameSigTau8 R x y = true := by
    simp only [pairSameSigTau8, Bool.and_eq_true]
    exact ⟨allBV3_of_forall fun b => beq_iff_eq.mpr (hsig b),
           allBV3_of_forall fun a => beq_iff_eq.mpr (htau a)⟩
  rw [hneB, hpair] at hxy
  simp at hxy

def iterPow6 (R : BitVec 6 → BitVec 6) : Nat → BitVec 6 → BitVec 6
  | 0, u => u
  | k + 1, u => iterPow6 R k (R u)

def sumBV3 (f : BitVec 3 → Nat) : Nat :=
  f 0#3 + f 1#3 + f 2#3 + f 3#3 + f 4#3 + f 5#3 + f 6#3 + f 7#3

def colTorus2K (R : BitVec 6 → BitVec 6) (k : Nat) : Nat :=
  sumBV3 fun a => sumBV3 fun b =>
    let u : BitVec 6 := a ++ b
    if iterPow6 R k u == u then 1 else 0

def fig8Step (R Rinv : BitVec 6 → BitVec 6) (t : BitVec 9) : BitVec 9 :=
  yb8Apply23 Rinv (yb8Apply12 R t)

def colFig8 (R Rinv : BitVec 6 → BitVec 6) : Nat :=
  sumBV3 fun a => sumBV3 fun b => sumBV3 fun c =>
    let t : BitVec 9 := a ++ b ++ c
    if fig8Step R Rinv (fig8Step R Rinv t) == t then 1 else 0

/-! ## Section 13: 2+1D Zamolodchikov Tetrahedron Equation ($|X|=2$ & $|X|=4$) -/

def get1_6 (s : BitVec 6) (idx : Nat) : BitVec 1 :=
  BitVec.extractLsb' (5 - idx) 1 s

def set1_6 (s : BitVec 6) (idx : Nat) (v : BitVec 1) : BitVec 6 :=
  let sh := 5 - idx
  let clr : BitVec 6 := ~~~((1#6) <<< sh)
  (s &&& clr) ||| ((BitVec.zeroExtend 6 v) <<< sh)

def tetra2Apply (T : BitVec 3 → BitVec 3) (i j k : Nat) (s : BitVec 6) : BitVec 6 :=
  let u : BitVec 3 := get1_6 s i ++ get1_6 s j ++ get1_6 s k
  let o : BitVec 3 := T u
  let s1 := set1_6 s i (BitVec.extractLsb' 2 1 o)
  let s2 := set1_6 s1 j (BitVec.extractLsb' 1 1 o)
  set1_6 s2 k (BitVec.extractLsb' 0 1 o)

def IsSetTheoreticTetra2 (T : BitVec 3 → BitVec 3) : Prop :=
  ∀ s : BitVec 6,
    tetra2Apply T 0 1 2 (tetra2Apply T 0 3 4 (tetra2Apply T 1 3 5 (tetra2Apply T 2 4 5 s))) =
      tetra2Apply T 2 4 5 (tetra2Apply T 1 3 5 (tetra2Apply T 0 3 4 (tetra2Apply T 0 1 2 s)))

def toffoli3 (u : BitVec 3) : BitVec 3 :=
  u ^^^ (((BitVec.extractLsb' 2 1 u &&& BitVec.extractLsb' 0 1 u).zeroExtend 3) <<< 1)

theorem toffoli3_satisfies_tetra2 : IsSetTheoreticTetra2 toffoli3 := by
  intro s
  simp only [tetra2Apply, get1_6, set1_6, toffoli3]
  bv_decide

def get2_12 (s : BitVec 12) (idx : Nat) : BitVec 2 :=
  BitVec.extractLsb' (2 * (5 - idx)) 2 s

def set2_12 (s : BitVec 12) (idx : Nat) (v : BitVec 2) : BitVec 12 :=
  let sh := 2 * (5 - idx)
  let clr : BitVec 12 := ~~~((3#12) <<< sh)
  (s &&& clr) ||| ((BitVec.zeroExtend 12 v) <<< sh)

def tetra4Apply (T : BitVec 6 → BitVec 6) (i j k : Nat) (s : BitVec 12) : BitVec 12 :=
  let u : BitVec 6 := get2_12 s i ++ get2_12 s j ++ get2_12 s k
  let o : BitVec 6 := T u
  let s1 := set2_12 s i (BitVec.extractLsb' 4 2 o)
  let s2 := set2_12 s1 j (BitVec.extractLsb' 2 2 o)
  set2_12 s2 k (BitVec.extractLsb' 0 2 o)

def IsSetTheoreticTetra4 (T : BitVec 6 → BitVec 6) : Prop :=
  ∀ s : BitVec 12,
    tetra4Apply T 0 1 2 (tetra4Apply T 0 3 4 (tetra4Apply T 1 3 5 (tetra4Apply T 2 4 5 s))) =
      tetra4Apply T 2 4 5 (tetra4Apply T 1 3 5 (tetra4Apply T 0 3 4 (tetra4Apply T 0 1 2 s)))

def Is3ChannelCoupled4 (T : BitVec 6 → BitVec 6) : Prop :=
  (∃ a0 a1 b c : BitVec 2, site1_3 (T (pack3_2 a1 b c)) ≠ site1_3 (T (pack3_2 a0 b c))) ∧
    (∃ a b c0 c1 : BitVec 2, site1_3 (T (pack3_2 a b c1)) ≠ site1_3 (T (pack3_2 a b c0)))

theorem coupled3_of_witnesses (T : BitVec 6 → BitVec 6)
    (a0 a1 b0 c0 : BitVec 2)
    (hL : site1_3 (T (pack3_2 a1 b0 c0)) ≠ site1_3 (T (pack3_2 a0 b0 c0)))
    (a2 b1 c1 c2 : BitVec 2)
    (hR : site1_3 (T (pack3_2 a2 b1 c2)) ≠ site1_3 (T (pack3_2 a2 b1 c1))) :
    Is3ChannelCoupled4 T :=
  ⟨⟨a0, a1, b0, c0, hL⟩, ⟨a2, b1, c1, c2, hR⟩⟩

/-! ## Section 14: Universal AME(4,4) Affine Rigidity No-Go Theorem (`BitVec 4`) -/

def getNib4 (lut : BitVec 64) (i : Nat) : BitVec 4 :=
  BitVec.extractLsb' (i * 4) 4 lut

def allDistinct4BV2 (x0 x1 x2 x3 : BitVec 2) : Bool :=
  (x0 != x1) && (x0 != x2) && (x0 != x3) && (x1 != x2) && (x1 != x3) && (x2 != x3)

def allDistinct16Lut (lut : BitVec 64) : Bool :=
  (List.range 16).all fun i =>
    (List.range 16).all fun j =>
      (i == j) || (getNib4 lut i != getNib4 lut j)

def isLatinC4 (lut : BitVec 64) : Bool :=
  (List.range 4).all (fun a =>
    allDistinct4BV2 (hi2 (getNib4 lut (4 * a))) (hi2 (getNib4 lut (4 * a + 1)))
                    (hi2 (getNib4 lut (4 * a + 2))) (hi2 (getNib4 lut (4 * a + 3)))) &&
  (List.range 4).all (fun b =>
    allDistinct4BV2 (hi2 (getNib4 lut b)) (hi2 (getNib4 lut (4 + b)))
                    (hi2 (getNib4 lut (8 + b))) (hi2 (getNib4 lut (12 + b))))

def isLatinD4 (lut : BitVec 64) : Bool :=
  (List.range 4).all (fun a =>
    allDistinct4BV2 (lo2 (getNib4 lut (4 * a))) (lo2 (getNib4 lut (4 * a + 1)))
                    (lo2 (getNib4 lut (4 * a + 2))) (lo2 (getNib4 lut (4 * a + 3)))) &&
  (List.range 4).all (fun b =>
    allDistinct4BV2 (lo2 (getNib4 lut b)) (lo2 (getNib4 lut (4 + b)))
                    (lo2 (getNib4 lut (8 + b))) (lo2 (getNib4 lut (12 + b))))

def is2Unitary4Lut (lut : BitVec 64) : Bool :=
  allDistinct16Lut lut && isLatinC4 lut && isLatinD4 lut

def isAffine4Lut (lut : BitVec 64) : Bool :=
  let v0 := getNib4 lut 0
  (List.range 16).all fun x =>
    let exp :=
      v0 ^^^
      (if (x &&& 1) != 0 then getNib4 lut 1 ^^^ v0 else 0#4) ^^^
      (if (x &&& 2) != 0 then getNib4 lut 2 ^^^ v0 else 0#4) ^^^
      (if (x &&& 4) != 0 then getNib4 lut 4 ^^^ v0 else 0#4) ^^^
      (if (x &&& 8) != 0 then getNib4 lut 8 ^^^ v0 else 0#4)
    getNib4 lut x == exp

def packLut4 (P : BitVec 4 → BitVec 4) : BitVec 64 :=
  ((BitVec.zeroExtend 64 (P 0#4)) <<< 0) |||
  ((BitVec.zeroExtend 64 (P 1#4)) <<< 4) |||
  ((BitVec.zeroExtend 64 (P 2#4)) <<< 8) |||
  ((BitVec.zeroExtend 64 (P 3#4)) <<< 12) |||
  ((BitVec.zeroExtend 64 (P 4#4)) <<< 16) |||
  ((BitVec.zeroExtend 64 (P 5#4)) <<< 20) |||
  ((BitVec.zeroExtend 64 (P 6#4)) <<< 24) |||
  ((BitVec.zeroExtend 64 (P 7#4)) <<< 28) |||
  ((BitVec.zeroExtend 64 (P 8#4)) <<< 32) |||
  ((BitVec.zeroExtend 64 (P 9#4)) <<< 36) |||
  ((BitVec.zeroExtend 64 (P 10#4)) <<< 40) |||
  ((BitVec.zeroExtend 64 (P 11#4)) <<< 44) |||
  ((BitVec.zeroExtend 64 (P 12#4)) <<< 48) |||
  ((BitVec.zeroExtend 64 (P 13#4)) <<< 52) |||
  ((BitVec.zeroExtend 64 (P 14#4)) <<< 56) |||
  ((BitVec.zeroExtend 64 (P 15#4)) <<< 60)

theorem ame4_affine_rigidity_nogo (lut : BitVec 64)
    (h : is2Unitary4Lut lut = true) : isAffine4Lut lut = true := by
  simp only [is2Unitary4Lut, allDistinct16Lut, isLatinC4, isLatinD4, isAffine4Lut,
             allDistinct4BV2, getNib4, hi2, lo2, List.range, List.range.loop, List.all,
             Bool.and_true] at h ⊢
  bv_decide

theorem ame4_eval_packLut4 (P : BitVec 4 → BitVec 4)
    (h : is2Unitary4Lut (packLut4 P) = true) : isAffine4Lut (packLut4 P) = true :=
  ame4_affine_rigidity_nogo (packLut4 P) h

/-! ## Section 15: 3+1D 4-Simplex Equation (`|X|=2` on `BitVec 10` & `|X|=4` on `BitVec 20`) -/

def get1_10 (s : BitVec 10) (idx : Nat) : BitVec 1 :=
  BitVec.extractLsb' (9 - idx) 1 s

def set1_10 (s : BitVec 10) (idx : Nat) (v : BitVec 1) : BitVec 10 :=
  let sh := 9 - idx
  let clr : BitVec 10 := ~~~((1#10) <<< sh)
  (s &&& clr) ||| ((BitVec.zeroExtend 10 v) <<< sh)

def simplex4Apply2 (S : BitVec 4 → BitVec 4) (i j k l : Nat) (s : BitVec 10) : BitVec 10 :=
  let u : BitVec 4 := get1_10 s i ++ get1_10 s j ++ get1_10 s k ++ get1_10 s l
  let o : BitVec 4 := S u
  let s1 := set1_10 s i (BitVec.extractLsb' 3 1 o)
  let s2 := set1_10 s1 j (BitVec.extractLsb' 2 1 o)
  let s3 := set1_10 s2 k (BitVec.extractLsb' 1 1 o)
  set1_10 s3 l (BitVec.extractLsb' 0 1 o)

def IsSetTheoretic4Simplex2 (S : BitVec 4 → BitVec 4) : Prop :=
  ∀ s : BitVec 10,
    simplex4Apply2 S 0 1 3 6 (simplex4Apply2 S 0 2 4 7 (simplex4Apply2 S 1 2 5 8
      (simplex4Apply2 S 3 4 5 9 (simplex4Apply2 S 6 7 8 9 s)))) =
    simplex4Apply2 S 6 7 8 9 (simplex4Apply2 S 3 4 5 9 (simplex4Apply2 S 1 2 5 8
      (simplex4Apply2 S 0 2 4 7 (simplex4Apply2 S 0 1 3 6 s))))

def ccx02_4 (u : BitVec 4) : BitVec 4 :=
  u ^^^ (((BitVec.extractLsb' 3 1 u &&& BitVec.extractLsb' 1 1 u).zeroExtend 4) <<< 2)

def ccx13_4 (u : BitVec 4) : BitVec 4 :=
  u ^^^ (((BitVec.extractLsb' 2 1 u &&& BitVec.extractLsb' 0 1 u).zeroExtend 4) <<< 1)

def cccx1_4 (u : BitVec 4) : BitVec 4 :=
  u ^^^ (((BitVec.extractLsb' 3 1 u &&& BitVec.extractLsb' 1 1 u &&&
            BitVec.extractLsb' 0 1 u).zeroExtend 4) <<< 2)

theorem ccx02_satisfies_4simplex2 : IsSetTheoretic4Simplex2 ccx02_4 := by
  intro s
  simp only [simplex4Apply2, get1_10, set1_10, ccx02_4]
  bv_decide

theorem ccx13_satisfies_4simplex2 : IsSetTheoretic4Simplex2 ccx13_4 := by
  intro s
  simp only [simplex4Apply2, get1_10, set1_10, ccx13_4]
  bv_decide

theorem cccx1_witness_fails_4simplex2 (s : BitVec 10)
    (h : simplex4Apply2 cccx1_4 0 1 3 6 (simplex4Apply2 cccx1_4 0 2 4 7 (simplex4Apply2 cccx1_4 1 2 5 8
           (simplex4Apply2 cccx1_4 3 4 5 9 (simplex4Apply2 cccx1_4 6 7 8 9 s)))) ≠
         simplex4Apply2 cccx1_4 6 7 8 9 (simplex4Apply2 cccx1_4 3 4 5 9 (simplex4Apply2 cccx1_4 1 2 5 8
           (simplex4Apply2 cccx1_4 0 2 4 7 (simplex4Apply2 cccx1_4 0 1 3 6 s))))) :
    ¬ IsSetTheoretic4Simplex2 cccx1_4 :=
  fun hAll => h (hAll s)

theorem cccx1_fails_4simplex2 : ¬ IsSetTheoretic4Simplex2 cccx1_4 :=
  cccx1_witness_fails_4simplex2 0x3DF#10 (by decide)

def get2_20 (s : BitVec 20) (idx : Nat) : BitVec 2 :=
  BitVec.extractLsb' (2 * (9 - idx)) 2 s

def set2_20 (s : BitVec 20) (idx : Nat) (v : BitVec 2) : BitVec 20 :=
  let sh := 2 * (9 - idx)
  let clr : BitVec 20 := ~~~((3#20) <<< sh)
  (s &&& clr) ||| ((BitVec.zeroExtend 20 v) <<< sh)

def simplex4Apply4 (S : BitVec 8 → BitVec 8) (i j k l : Nat) (s : BitVec 20) : BitVec 20 :=
  let u : BitVec 8 := get2_20 s i ++ get2_20 s j ++ get2_20 s k ++ get2_20 s l
  let o : BitVec 8 := S u
  let s1 := set2_20 s i (BitVec.extractLsb' 6 2 o)
  let s2 := set2_20 s1 j (BitVec.extractLsb' 4 2 o)
  let s3 := set2_20 s2 k (BitVec.extractLsb' 2 2 o)
  set2_20 s3 l (BitVec.extractLsb' 0 2 o)

def IsSetTheoretic4Simplex4 (S : BitVec 8 → BitVec 8) : Prop :=
  ∀ s : BitVec 20,
    simplex4Apply4 S 0 1 3 6 (simplex4Apply4 S 0 2 4 7 (simplex4Apply4 S 1 2 5 8
      (simplex4Apply4 S 3 4 5 9 (simplex4Apply4 S 6 7 8 9 s)))) =
    simplex4Apply4 S 6 7 8 9 (simplex4Apply4 S 3 4 5 9 (simplex4Apply4 S 1 2 5 8
      (simplex4Apply4 S 0 2 4 7 (simplex4Apply4 S 0 1 3 6 s))))

def s4RingLift (u : BitVec 8) : BitVec 8 :=
  let a := BitVec.extractLsb' 6 2 u
  let b := BitVec.extractLsb' 4 2 u
  let c := BitVec.extractLsb' 2 2 u
  let d := BitVec.extractLsb' 0 2 u
  a ++ (b + a * c) ++ c ++ d

def s4RingLiftInv (u : BitVec 8) : BitVec 8 :=
  let a := BitVec.extractLsb' 6 2 u
  let b := BitVec.extractLsb' 4 2 u
  let c := BitVec.extractLsb' 2 2 u
  let d := BitVec.extractLsb' 0 2 u
  a ++ (b - a * c) ++ c ++ d

def s4RingDual (u : BitVec 8) : BitVec 8 :=
  let a := BitVec.extractLsb' 6 2 u
  let b := BitVec.extractLsb' 4 2 u
  let c := BitVec.extractLsb' 2 2 u
  let d := BitVec.extractLsb' 0 2 u
  a ++ b ++ (c + b * d) ++ d

def s4Coupled4Ch (u : BitVec 8) : BitVec 8 :=
  let a := BitVec.extractLsb' 6 2 u
  let b := BitVec.extractLsb' 4 2 u
  let c := BitVec.extractLsb' 2 2 u
  let d := BitVec.extractLsb' 0 2 u
  let b' := b ^^^ ((a &&& c) &&& 1#2)
  let c' := c ^^^ ((b &&& d) &&& 2#2)
  a ++ b' ++ c' ++ d

theorem s4RingLift_inv_left (u : BitVec 8) : s4RingLiftInv (s4RingLift u) = u := by
  simp only [s4RingLift, s4RingLiftInv]
  bv_decide

theorem s4RingLift_inv_right (u : BitVec 8) : s4RingLift (s4RingLiftInv u) = u := by
  simp only [s4RingLift, s4RingLiftInv]
  bv_decide

theorem s4RingLift_bijective : Bijective s4RingLift := by
  refine ⟨fun a b h => ?_, fun b => ⟨s4RingLiftInv b, s4RingLift_inv_right b⟩⟩
  have h1 := congrArg s4RingLiftInv h
  simpa [s4RingLift_inv_left] using h1

theorem s4RingLift_noninvolutive : ¬ (∀ u : BitVec 8, s4RingLift (s4RingLift u) = u) := by
  intro h
  have hw := h (1#2 ++ 0#2 ++ 1#2 ++ 0#2)
  revert hw
  decide

theorem s4RingLift_satisfies_4simplex4 : IsSetTheoretic4Simplex4 s4RingLift := by
  intro s
  simp only [simplex4Apply4, get2_20, set2_20, s4RingLift]
  bv_decide

theorem s4RingDual_satisfies_4simplex4 : IsSetTheoretic4Simplex4 s4RingDual := by
  intro s
  simp only [simplex4Apply4, get2_20, set2_20, s4RingDual]
  bv_decide

theorem s4Coupled4Ch_satisfies_4simplex4 : IsSetTheoretic4Simplex4 s4Coupled4Ch := by
  intro s
  simp only [simplex4Apply4, get2_20, set2_20, s4Coupled4Ch]
  bv_decide

def Is4ChannelCoupled4 (S : BitVec 8 → BitVec 8) : Prop :=
  (∃ a0 a1 b c d : BitVec 2,
      BitVec.extractLsb' 4 2 (S (a1 ++ b ++ c ++ d)) ≠
        BitVec.extractLsb' 4 2 (S (a0 ++ b ++ c ++ d))) ∧
  (∃ a b0 b1 c d : BitVec 2,
      BitVec.extractLsb' 2 2 (S (a ++ b1 ++ c ++ d)) ≠
        BitVec.extractLsb' 2 2 (S (a ++ b0 ++ c ++ d))) ∧
  (∃ a b c0 c1 d : BitVec 2,
      BitVec.extractLsb' 4 2 (S (a ++ b ++ c1 ++ d)) ≠
        BitVec.extractLsb' 4 2 (S (a ++ b ++ c0 ++ d))) ∧
  (∃ a b c d0 d1 : BitVec 2,
      BitVec.extractLsb' 2 2 (S (a ++ b ++ c ++ d1)) ≠
        BitVec.extractLsb' 2 2 (S (a ++ b ++ c ++ d0)))

theorem s4Coupled4Ch_4channel_coupled : Is4ChannelCoupled4 s4Coupled4Ch :=
  ⟨⟨0#2, 1#2, 0#2, 1#2, 0#2, by decide⟩,
   ⟨0#2, 0#2, 2#2, 0#2, 2#2, by decide⟩,
   ⟨1#2, 0#2, 0#2, 1#2, 0#2, by decide⟩,
   ⟨0#2, 2#2, 0#2, 0#2, 2#2, by decide⟩⟩

/-! ## Section 16: Knot Cohomology (`H^2_Q(Q_6, Z_4)`) & Trefoil Chirality -/

def sumFin6 (f : Fin 6 → Nat) : Nat :=
  f 0 + f 1 + f 2 + f 3 + f 4 + f 5

def allFin6 (p : Fin 6 → Bool) : Bool :=
  p 0 && p 1 && p 2 && p 3 && p 4 && p 5

def q6Act (b a : Fin 6) : Fin 6 :=
  match b.val, a.val with
  | 0, 0 => 0 | 0, 1 => 1 | 0, 2 => 3 | 0, 3 => 4 | 0, 4 => 5 | 0, 5 => 2
  | 1, 0 => 0 | 1, 1 => 1 | 1, 2 => 5 | 1, 3 => 2 | 1, 4 => 3 | 1, 5 => 4
  | 2, 0 => 5 | 2, 1 => 3 | 2, 2 => 2 | 2, 3 => 0 | 2, 4 => 4 | 2, 5 => 1
  | 3, 0 => 2 | 3, 1 => 4 | 3, 2 => 1 | 3, 3 => 3 | 3, 4 => 0 | 3, 5 => 5
  | 4, 0 => 3 | 4, 1 => 5 | 4, 2 => 2 | 4, 3 => 1 | 4, 4 => 4 | 4, 5 => 0
  | 5, 0 => 4 | 5, 1 => 2 | 5, 2 => 0 | 5, 3 => 3 | 5, 4 => 1 | 5, 5 => 5
  | _, _ => a

def q6R (p : Fin 6 × Fin 6) : Fin 6 × Fin 6 :=
  (p.2, q6Act p.2 p.1)

def q6RInv (p : Fin 6 × Fin 6) : Fin 6 × Fin 6 :=
  (q6Act p.1 (q6Act p.1 (q6Act p.1 p.2)), p.1)

def q6Phi (p : Fin 6 × Fin 6) : BitVec 2 :=
  match p.1.val, p.2.val with
  | 0, 0 => 0#2 | 0, 1 => 1#2 | 0, 2 => 0#2 | 0, 3 => 0#2 | 0, 4 => 0#2 | 0, 5 => 1#2
  | 1, 0 => 1#2 | 1, 1 => 0#2 | 1, 2 => 2#2 | 1, 3 => 2#2 | 1, 4 => 0#2 | 1, 5 => 3#2
  | 2, 0 => 0#2 | 2, 1 => 1#2 | 2, 2 => 0#2 | 2, 3 => 2#2 | 2, 4 => 1#2 | 2, 5 => 1#2
  | 3, 0 => 1#2 | 3, 1 => 1#2 | 3, 2 => 1#2 | 3, 3 => 0#2 | 3, 4 => 3#2 | 3, 5 => 1#2
  | 4, 0 => 3#2 | 4, 1 => 0#2 | 4, 2 => 1#2 | 4, 3 => 0#2 | 4, 4 => 0#2 | 4, 5 => 3#2
  | 5, 0 => 0#2 | 5, 1 => 2#2 | 5, 2 => 1#2 | 5, 3 => 1#2 | 5, 4 => 1#2 | 5, 5 => 0#2
  | _, _ => 0#2

def IsSetTheoreticYB6 (R : Fin 6 × Fin 6 → Fin 6 × Fin 6) : Prop :=
  ∀ a b c : Fin 6,
    let ab := R (a, b)
    let bc := R (ab.2, c)
    let ab2 := R (ab.1, bc.1)
    let bc0 := R (b, c)
    let ab0 := R (a, bc0.1)
    let bc1 := R (ab0.2, bc0.2)
    (ab2.1, ab2.2, bc.2) = (ab0.1, bc1.1, bc1.2)

def IsR1_6 (R : Fin 6 × Fin 6 → Fin 6 × Fin 6) : Prop :=
  ∀ x : Fin 6, R (x, x) = (x, x)

def Is2CocycleMod4_6 (R : Fin 6 × Fin 6 → Fin 6 × Fin 6) (phi : Fin 6 × Fin 6 → BitVec 2) :
    Prop :=
  (∀ x : Fin 6, phi (x, x) = 0#2) ∧
    ∀ a b c : Fin 6,
      let ab := R (a, b)
      let bc := R (ab.2, c)
      let bc0 := R (b, c)
      let ab0 := R (a, bc0.1)
      phi (a, b) + phi (ab.2, c) + phi (ab.1, bc.1) =
        phi (b, c) + phi (a, bc0.1) + phi (ab0.2, bc0.2)

def coboundary2 {α : Type} (R : α × α → α × α) (f : α → BitVec 2) (u : α × α) : BitVec 2 :=
  let v := R u
  (f v.1 + f v.2) - (f u.1 + f u.2)

def Is2CoboundaryMod4_6 (R : Fin 6 × Fin 6 → Fin 6 × Fin 6)
    (phi : Fin 6 × Fin 6 → BitVec 2) : Prop :=
  ∃ f : Fin 6 → BitVec 2, ∀ u : Fin 6 × Fin 6, phi u = coboundary2 R f u

def iterPair {α : Type} (R : α × α → α × α) : Nat → α × α → α × α
  | 0, u => u
  | k + 1, u => iterPair R k (R u)

def orbitCocycleSum {α : Type} (R : α × α → α × α) (phi : α × α → BitVec 2) :
    Nat → α × α → BitVec 2
  | 0, _ => 0#2
  | k + 1, u => phi u + orbitCocycleSum R phi k (R u)

private theorem bv2_tel (a b c : BitVec 2) : (b - a) + (c - b) = c - a := by
  bv_decide

private theorem bv2_sub_self (a : BitVec 2) : a - a = 0#2 := by
  bv_decide

theorem orbitCocycleSum_coboundary {α : Type} (R : α × α → α × α) (f : α → BitVec 2) :
    ∀ (k : Nat) (u : α × α),
      orbitCocycleSum R (coboundary2 R f) k u =
        (f (iterPair R k u).1 + f (iterPair R k u).2) - (f u.1 + f u.2)
  | 0, u => by
    simp [orbitCocycleSum, iterPair]
  | k + 1, u => by
    simp only [orbitCocycleSum, iterPair]
    rw [orbitCocycleSum_coboundary R f k (R u)]
    unfold coboundary2
    exact bv2_tel (f u.1 + f u.2) (f (R u).1 + f (R u).2)
      (f (iterPair R k (R u)).1 + f (iterPair R k (R u)).2)

theorem coboundary_vanishes_torus2k {α : Type} (R : α × α → α × α) (f : α → BitVec 2)
    (k : Nat) (u : α × α) (hFix : iterPair R k u = u) :
    orbitCocycleSum R (coboundary2 R f) k u = 0#2 := by
  rw [orbitCocycleSum_coboundary, hFix, bv2_sub_self]

private theorem orbitCocycleSum_eq_coboundary {α : Type} (R : α × α → α × α)
    (phi : α × α → BitVec 2) (f : α → BitVec 2) (hf : ∀ u, phi u = coboundary2 R f u) :
    ∀ (k : Nat) (u : α × α),
      orbitCocycleSum R phi k u = orbitCocycleSum R (coboundary2 R f) k u
  | 0, _ => rfl
  | k + 1, u => by
    simp only [orbitCocycleSum, hf u, orbitCocycleSum_eq_coboundary R phi f hf k (R u)]

theorem q6_non_coboundary_of_witness (R : Fin 6 × Fin 6 → Fin 6 × Fin 6)
    (phi : Fin 6 × Fin 6 → BitVec 2) (k : Nat) (u : Fin 6 × Fin 6)
    (hFix : iterPair R k u = u) (hNz : orbitCocycleSum R phi k u ≠ 0#2) :
    ¬ Is2CoboundaryMod4_6 R phi := by
  rintro ⟨f, hf⟩
  have hEq := orbitCocycleSum_eq_coboundary R phi f hf k u
  rw [coboundary_vanishes_torus2k R f k u hFix] at hEq
  exact hNz hEq

def iterPairPhiSum (R : Fin 6 × Fin 6 → Fin 6 × Fin 6) (phi : Fin 6 × Fin 6 → BitVec 2) :
    Nat → Fin 6 × Fin 6 → Fin 6 × Fin 6 × BitVec 2
  | 0, u => (u.1, u.2, 0#2)
  | k + 1, u =>
    let res := iterPairPhiSum R phi k (R u)
    (res.1, res.2.1, phi u + res.2.2)

def iterPairInvPhiSum (RInv : Fin 6 × Fin 6 → Fin 6 × Fin 6)
    (phi : Fin 6 × Fin 6 → BitVec 2) :
    Nat → Fin 6 × Fin 6 → Fin 6 × Fin 6 × BitVec 2
  | 0, u => (u.1, u.2, 0#2)
  | k + 1, u =>
    let y := RInv u
    let res := iterPairInvPhiSum RInv phi k y
    (res.1, res.2.1, (0#2 - phi y) + res.2.2)

def cocycleCountTorus2Pos (R : Fin 6 × Fin 6 → Fin 6 × Fin 6)
    (phi : Fin 6 × Fin 6 → BitVec 2) (k : Nat) (w : BitVec 2) : Nat :=
  sumFin6 fun a => sumFin6 fun b =>
    let res := iterPairPhiSum R phi k (a, b)
    if (res.1 == a) && (res.2.1 == b) && (res.2.2 == w) then 1 else 0

def cocycleCountTorus2Neg (RInv : Fin 6 × Fin 6 → Fin 6 × Fin 6)
    (phi : Fin 6 × Fin 6 → BitVec 2) (k : Nat) (w : BitVec 2) : Nat :=
  sumFin6 fun a => sumFin6 fun b =>
    let res := iterPairInvPhiSum RInv phi k (a, b)
    if (res.1 == a) && (res.2.1 == b) && (res.2.2 == w) then 1 else 0

def fig8CocycleStep (R RInv : Fin 6 × Fin 6 → Fin 6 × Fin 6)
    (phi : Fin 6 × Fin 6 → BitVec 2) (a b c : Fin 6) (s0 : BitVec 2) :
    Fin 6 × Fin 6 × Fin 6 × BitVec 2 :=
  let u1 := (a, b)
  let o1 := R u1
  let s1 := s0 + phi u1
  let o2 := RInv (o1.2, c)
  let s2 := s1 + (0#2 - phi o2)
  (o1.1, o2.1, o2.2, s2)

def cocycleCountFig8 (R RInv : Fin 6 × Fin 6 → Fin 6 × Fin 6)
    (phi : Fin 6 × Fin 6 → BitVec 2) (w : BitVec 2) : Nat :=
  sumFin6 fun a => sumFin6 fun b => sumFin6 fun c =>
    let st1 := fig8CocycleStep R RInv phi a b c 0#2
    let st2 := fig8CocycleStep R RInv phi st1.1 st1.2.1 st1.2.2.1 st1.2.2.2
    if (st2.1 == a) && (st2.2.1 == b) && (st2.2.2.1 == c) && (st2.2.2.2 == w) then 1 else 0

private theorem fin6_cases (x : Fin 6) :
    x = 0 ∨ x = 1 ∨ x = 2 ∨ x = 3 ∨ x = 4 ∨ x = 5 := by
  rcases x with ⟨v, hv⟩
  match v, hv with
  | 0, _ => exact Or.inl rfl
  | 1, _ => exact Or.inr (Or.inl rfl)
  | 2, _ => exact Or.inr (Or.inr (Or.inl rfl))
  | 3, _ => exact Or.inr (Or.inr (Or.inr (Or.inl rfl)))
  | 4, _ => exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inl rfl))))
  | 5, _ => exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inr rfl))))

theorem q6_inv_left (u : Fin 6 × Fin 6) : q6RInv (q6R u) = u := by
  rcases u with ⟨a, b⟩
  rcases fin6_cases a with rfl | rfl | rfl | rfl | rfl | rfl <;>
    rcases fin6_cases b with rfl | rfl | rfl | rfl | rfl | rfl <;> decide

theorem q6_inv_right (u : Fin 6 × Fin 6) : q6R (q6RInv u) = u := by
  rcases u with ⟨a, b⟩
  rcases fin6_cases a with rfl | rfl | rfl | rfl | rfl | rfl <;>
    rcases fin6_cases b with rfl | rfl | rfl | rfl | rfl | rfl <;> decide

theorem q6_bijective : Bijective q6R := by
  refine ⟨fun a b h => ?_, fun b => ⟨q6RInv b, q6_inv_right b⟩⟩
  have h1 := congrArg q6RInv h
  simpa [q6_inv_left] using h1

theorem q6_r1 : IsR1_6 q6R := by
  intro x
  rcases fin6_cases x with rfl | rfl | rfl | rfl | rfl | rfl <;> decide

theorem q6_satisfies_ybe : IsSetTheoreticYB6 q6R := by
  intro a b c
  rcases fin6_cases a with rfl | rfl | rfl | rfl | rfl | rfl <;>
    rcases fin6_cases b with rfl | rfl | rfl | rfl | rfl | rfl <;>
    rcases fin6_cases c with rfl | rfl | rfl | rfl | rfl | rfl <;> decide

theorem q6_noninvolutive : q6R (q6R (0, 2)) ≠ (0, 2) := by
  decide

theorem q6_phi_is_2cocycle : Is2CocycleMod4_6 q6R q6Phi := by
  refine ⟨fun x => ?_, fun a b c => ?_⟩
  · rcases fin6_cases x with rfl | rfl | rfl | rfl | rfl | rfl <;> decide
  · rcases fin6_cases a with rfl | rfl | rfl | rfl | rfl | rfl <;>
      rcases fin6_cases b with rfl | rfl | rfl | rfl | rfl | rfl <;>
      rcases fin6_cases c with rfl | rfl | rfl | rfl | rfl | rfl <;> decide

theorem q6_phi_not_coboundary : ¬ Is2CoboundaryMod4_6 q6R q6Phi :=
  q6_non_coboundary_of_witness q6R q6Phi 3 (0, 2) (by decide) (by decide)

theorem q6_trefoil_chiral_separation :
    cocycleCountTorus2Pos q6R q6Phi 1 0#2 = 6 ∧
    cocycleCountTorus2Pos q6R q6Phi 2 0#2 = 6 ∧
    cocycleCountTorus2Pos q6R q6Phi 2 2#2 = 6 ∧
    cocycleCountTorus2Pos q6R q6Phi 3 1#2 = 24 ∧
    cocycleCountTorus2Neg q6RInv q6Phi 3 1#2 = 0 ∧
    cocycleCountTorus2Pos q6R q6Phi 3 3#2 = 0 ∧
    cocycleCountTorus2Neg q6RInv q6Phi 3 3#2 = 24 ∧
    cocycleCountFig8 q6R q6RInv q6Phi 0#2 = 6 := by
  decide

/-! ## Section 17: Exact `Z[i]` Quantum Dual-Unitary, Quantum Scar & Period-2 DTC -/

structure GaussInt where
  re : Int
  im : Int
  deriving DecidableEq, Repr

def gZero : GaussInt := ⟨0, 0⟩
def gOne : GaussInt := ⟨1, 0⟩
def gI : GaussInt := ⟨0, 1⟩
def gAdd (x y : GaussInt) : GaussInt := ⟨x.re + y.re, x.im + y.im⟩
def gMul (x y : GaussInt) : GaussInt := ⟨x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re⟩
def gConj (x : GaussInt) : GaussInt := ⟨x.re, -x.im⟩

def gSumFin4 (f : Fin 4 → GaussInt) : GaussInt :=
  gAdd (gAdd (f 0) (f 1)) (gAdd (f 2) (f 3))

def gSumBV4 (f : BitVec 4 → GaussInt) : GaussInt :=
  gSumFin4 fun hi : Fin 4 =>
    gSumFin4 fun lo : Fin 4 =>
      f (BitVec.ofNat 2 hi.val ++ BitVec.ofNat 2 lo.val)

def wCxDU (i j : Fin 4) : GaussInt :=
  match i.val, j.val with
  | 0, 0 => ⟨1, 0⟩  | 0, 1 => ⟨0, 1⟩  | 0, 2 => ⟨1, 0⟩  | 0, 3 => ⟨0, -1⟩
  | 1, 0 => ⟨1, 0⟩  | 1, 1 => ⟨0, 1⟩  | 1, 2 => ⟨-1, 0⟩ | 1, 3 => ⟨0, 1⟩
  | 2, 0 => ⟨1, 0⟩  | 2, 1 => ⟨0, -1⟩ | 2, 2 => ⟨1, 0⟩  | 2, 3 => ⟨0, 1⟩
  | 3, 0 => ⟨1, 0⟩  | 3, 1 => ⟨0, -1⟩ | 3, 2 => ⟨-1, 0⟩ | 3, 3 => ⟨0, -1⟩
  | _, _ => gZero

def wCxDUReshuffle (ac bd : Fin 4) : GaussInt :=
  let a := ac.val / 2
  let c := ac.val % 2
  let b := bd.val / 2
  let d := bd.val % 2
  wCxDU ⟨c * 2 + d, by omega⟩ ⟨a * 2 + b, by omega⟩

theorem wCxDU_is_unitary (j k : Fin 4) :
    gSumFin4 (fun i => gMul (gConj (wCxDU i j)) (wCxDU i k)) =
      if j == k then ⟨4, 0⟩ else gZero := by
  revert j k
  decide

theorem wCxDU_is_dual_unitary (j k : Fin 4) :
    gSumFin4 (fun i => gMul (gConj (wCxDUReshuffle i j)) (wCxDUReshuffle i k)) =
      if j == k then ⟨4, 0⟩ else gZero := by
  revert j k
  decide

theorem wCxDU_full_superposition (i j : Fin 4) : wCxDU i j ≠ gZero := by
  revert i j
  decide

def scarGate4 (sOut sIn : BitVec 4) : GaussInt :=
  let b0 := BitVec.extractLsb' 0 1 sIn
  let b3 := BitVec.extractLsb' 3 1 sIn
  if b0 != b3 then
    if sOut == sIn then ⟨2, 0⟩ else gZero
  else
    let b0' := BitVec.extractLsb' 0 1 sOut
    let b3' := BitVec.extractLsb' 3 1 sOut
    if (b0' == b0) && (b3' == b3) then
      let mOut := (BitVec.extractLsb' 2 1 sOut).toNat * 2 + (BitVec.extractLsb' 1 1 sOut).toNat
      let mIn := (BitVec.extractLsb' 2 1 sIn).toNat * 2 + (BitVec.extractLsb' 1 1 sIn).toNat
      wCxDU ⟨mOut, by omega⟩ ⟨mIn, by omega⟩
    else gZero

theorem scarGate4_is_unitary (j k : BitVec 4) :
    gSumBV4 (fun i => gMul (gConj (scarGate4 i j)) (scarGate4 i k)) =
      if j == k then ⟨4, 0⟩ else gZero := by
  revert j k
  decide

def neel0 (L : Nat) (i : Fin L) : Bool := (i.val % 2 == 1)
def neel1 (L : Nat) (i : Fin L) : Bool := (i.val % 2 == 0)

def flipAll {L : Nat} (s : Fin L → Bool) (i : Fin L) : Bool := !(s i)

theorem dtc_neel_local_protected {L : Nat} (hGe : 4 ≤ L) (hEven : L % 2 = 0) (i : Fin L) :
    neel0 L i ≠ neel0 L ⟨(i.val + 3) % L, Nat.mod_lt _ (by omega)⟩ ∧
    neel1 L i ≠ neel1 L ⟨(i.val + 3) % L, Nat.mod_lt _ (by omega)⟩ := by
  have hMod : (i.val + 3) % L % 2 = (i.val + 3) % 2 :=
    Nat.mod_mod_of_dvd (i.val + 3) (Nat.dvd_of_mod_eq_zero hEven)
  simp only [neel0, neel1, hMod]
  rcases Nat.mod_two_eq_zero_or_one i.val with h | h <;> simp [Nat.add_mod, h]

theorem dtc_neel_exact_period2 {L : Nat} (hGe : 4 ≤ L) :
    flipAll (neel0 L) = neel1 L ∧
    flipAll (neel1 L) = neel0 L ∧
    neel0 L ≠ neel1 L ∧
    flipAll (flipAll (neel0 L)) = neel0 L := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · funext i; simp only [flipAll, neel0, neel1]
    rcases Nat.mod_two_eq_zero_or_one i.val with h | h <;> simp [h]
  · funext i; simp only [flipAll, neel0, neel1]
    rcases Nat.mod_two_eq_zero_or_one i.val with h | h <;> simp [h]
  · intro h
    have h0 := congrFun h ⟨0, by omega⟩
    simp [neel0, neel1] at h0
  · funext i; simp [flipAll]

/-! ## Section 18: Fracton Universal Involution & Dipole-Carrier Reversible Toffoli -/

theorem fracton3_universal_involution (T : BitVec 6 → BitVec 6)
    (hBij : Bijective T) (hLap : SatisfiesDipoleLaplacian T rho2) (u : BitVec 6) :
    T (T u) = u := by
  have h2 : ∀ p q r : BitVec 6,
      (let dp := rho2 (site0_3 q) - rho2 (site0_3 p)
       rho2 (site1_3 q) - rho2 (site1_3 p) = -(2 * dp) ∧
       rho2 (site2_3 q) - rho2 (site2_3 p) = dp) →
      (let dr := rho2 (site0_3 r) - rho2 (site0_3 p)
       rho2 (site1_3 r) - rho2 (site1_3 p) = -(2 * dr) ∧
       rho2 (site2_3 r) - rho2 (site2_3 p) = dr) →
      r = p ∨ r = q ∨ q = p := by
    intro p q r hq hr
    have hp : p = pack3_2 (site0_3 p) (site1_3 p) (site2_3 p) := (pack3_2_sites p).symm
    have hq' : q = pack3_2 (site0_3 q) (site1_3 q) (site2_3 q) := (pack3_2_sites q).symm
    have hr' : r = pack3_2 (site0_3 r) (site1_3 r) (site2_3 r) := (pack3_2_sites r).symm
    have hphi : ∀ a b c a' b' c' : BitVec 2,
        let phi := rho2 a' - rho2 a
        rho2 b' - rho2 b = -(2 * phi) → rho2 c' - rho2 c = phi →
        (a' = a ∧ b' = b ∧ c' = c) ∨
        (rho2 a' = rho2 a + 1 ∧ rho2 b' = rho2 b - 2 ∧ rho2 c' = rho2 c + 1) ∨
        (rho2 a' = rho2 a - 1 ∧ rho2 b' = rho2 b + 2 ∧ rho2 c' = rho2 c - 1) := by
      decide
    have h1 := hphi (site0_3 p) (site1_3 p) (site2_3 p) (site0_3 q) (site1_3 q) (site2_3 q) hq.1 hq.2
    have h2' := hphi (site0_3 p) (site1_3 p) (site2_3 p) (site0_3 r) (site1_3 r) (site2_3 r) hr.1 hr.2
    have huniq : ∀ x y : BitVec 2, rho2 x = rho2 y → x = y := by decide
    have hno_both : ∀ b b1 b2 : BitVec 2,
        rho2 b1 = rho2 b - 2 → rho2 b2 = rho2 b + 2 → False := by decide
    rcases h1 with ⟨ha1, hb1, hc1⟩ | ⟨ha1, hb1, hc1⟩ | ⟨ha1, hb1, hc1⟩
    · exact Or.inr (Or.inr (by rw [hq', hp, ha1, hb1, hc1]))
    · rcases h2' with ⟨ha2, hb2, hc2⟩ | ⟨ha2, hb2, hc2⟩ | ⟨ha2, hb2, hc2⟩
      · exact Or.inl (by rw [hr', hp, ha2, hb2, hc2])
      · exact Or.inr (Or.inl (by
          rw [hr', hq']
          rw [huniq _ _ (ha2.trans ha1.symm), huniq _ _ (hb2.trans hb1.symm),
              huniq _ _ (hc2.trans hc1.symm)]))
      · exact False.elim (hno_both _ _ _ hb1 hb2)
    · rcases h2' with ⟨ha2, hb2, hc2⟩ | ⟨ha2, hb2, hc2⟩ | ⟨ha2, hb2, hc2⟩
      · exact Or.inl (by rw [hr', hp, ha2, hb2, hc2])
      · exact False.elim (hno_both _ _ _ hb2 hb1)
      · exact Or.inr (Or.inl (by
          rw [hr', hq']
          rw [huniq _ _ (ha2.trans ha1.symm), huniq _ _ (hb2.trans hb1.symm),
              huniq _ _ (hc2.trans hc1.symm)]))
  by_cases hfix : T u = u
  · rw [hfix, hfix]
  · have hTu := hLap (site0_3 u) (site1_3 u) (site2_3 u)
    rw [pack3_2_sites] at hTu
    have hTTu := hLap (site0_3 (T u)) (site1_3 (T u)) (site2_3 (T u))
    rw [pack3_2_sites] at hTTu
    have hTTu_lap : (let dp := rho2 (site0_3 (T (T u))) - rho2 (site0_3 u)
        rho2 (site1_3 (T (T u))) - rho2 (site1_3 u) = -(2 * dp) ∧
        rho2 (site2_3 (T (T u))) - rho2 (site2_3 u) = dp) := by
      dsimp only at hTu hTTu ⊢
      constructor <;> omega
    rcases h2 u (T u) (T (T u)) hTu hTTu_lap with hEq | hEq | hEq
    · exact hEq
    · exact False.elim (hfix (hBij.1 (T u) u hEq))
    · exact False.elim (hfix hEq)

def dipoleToffoli (u : BitVec 6) : BitVec 6 :=
  if u == pack3_2 1#2 2#2 1#2 then pack3_2 2#2 0#2 2#2
  else if u == pack3_2 2#2 0#2 2#2 then pack3_2 1#2 2#2 1#2
  else u

def encodeDipoleLogical (ca cc tgt : Bool) : BitVec 6 :=
  if ca && cc && tgt then pack3_2 2#2 0#2 2#2
  else pack3_2 (if ca then 1#2 else 0#2) (if tgt then 0#2 else 2#2) (if cc then 1#2 else 0#2)

def decodeDipoleTarget (u : BitVec 6) : Bool :=
  site1_3 u == 0#2

theorem dipoleToffoli_involutive (u : BitVec 6) : dipoleToffoli (dipoleToffoli u) = u := by
  revert u
  decide

theorem dipoleToffoli_bijective : Bijective dipoleToffoli := by
  refine ⟨fun a b h => ?_, fun b => ⟨dipoleToffoli b, dipoleToffoli_involutive b⟩⟩
  have h1 := congrArg dipoleToffoli h
  simpa [dipoleToffoli_involutive] using h1

theorem dipoleToffoli_satisfies_laplacian : SatisfiesDipoleLaplacian dipoleToffoli rho2 := by
  unfold SatisfiesDipoleLaplacian rho2
  decide

theorem dipole_static_control_nogo (T : BitVec 6 → BitVec 6)
    (hL : SatisfiesDipoleLaplacian T rho2) (u : BitVec 6)
    (h0 : site0_3 (T u) = site0_3 u) (h2 : site2_3 (T u) = site2_3 u) :
    T u = u := by
  have h := hL (site0_3 u) (site1_3 u) (site2_3 u)
  rw [pack3_2_sites] at h
  dsimp only at h
  rw [h0] at h
  have h1 : site1_3 (T u) = site1_3 u := by
    have hs : rho2 (site1_3 (T u)) = rho2 (site1_3 u) := by omega
    have hinj : ∀ x y : BitVec 2, rho2 x = rho2 y → x = y := by decide
    exact hinj _ _ hs
  rw [← pack3_2_sites (T u), h0, h1, h2, pack3_2_sites]

theorem dipoleToffoli_open_chain_Q_D {L : Nat} (pos : Nat) (hpos : pos + 2 < L)
    (s : Fin L → BitVec 2) :
    openCharge L rho2 (applyGate3 dipoleToffoli pos hpos s) = openCharge L rho2 s ∧
    openDipole L rho2 (applyGate3 dipoleToffoli pos hpos s) = openDipole L rho2 s :=
  ⟨fracton_open_chain_Q_conservation dipoleToffoli rho2 dipoleToffoli_satisfies_laplacian pos hpos s,
   fracton_open_chain_D_conservation dipoleToffoli rho2 dipoleToffoli_satisfies_laplacian pos hpos s⟩

theorem dipoleToffoli_single_charge_frozen {L : Nat} (pos : Nat) (hpos : pos + 2 < L)
    (s : Fin L → BitVec 2) (hQ : openCharge L rho2 s ≤ 1) :
    applyGate3 dipoleToffoli pos hpos s = s :=
  fracton_single_charge_frozen dipoleToffoli dipoleToffoli_satisfies_laplacian pos hpos s hQ

theorem dipoleToffoli_implements_ccx (ca cc tgt : Bool) :
    decodeDipoleTarget (dipoleToffoli (encodeDipoleLogical ca cc tgt)) =
      xor tgt (ca && cc) := by
  revert ca cc tgt
  decide

theorem encodeDipoleLogical_injective (c1 c2 t1 c3 c4 t2 : Bool)
    (h : encodeDipoleLogical c1 c2 t1 = encodeDipoleLogical c3 c4 t2) :
    c1 = c3 ∧ c2 = c4 ∧ t1 = t2 := by
  revert c1 c2 t1 c3 c4 t2 h
  decide

end Phy.Noether
