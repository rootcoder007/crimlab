/-
  crimlab: inst/lean/Spatial.lean

  Identities behind the spatial functions.

  * The weighted mean centre minimises the weighted sum of squared
    distances (cgstat(), sdell()): coordinate by coordinate,
        ∑ w (x − c)² = ∑ w (x − m)² + W (m − c)²     (parallel_axis)
    with m the weighted mean and W the total weight, so with non-negative
    weights the sum is smallest at c = m (mean_minimises), and in the plane
    the mean centre minimises the sum of squared distances
    (mean_centre_minimises).
  * The local Moran statistics of lmorn() add up to S0 times global
    Moran's I (lisa_sum), whatever the weights.
  * The mean resultant length of dirmn(), the length of the weighted mean
    of unit vectors, lies in [0, 1] (resultant_nonneg, resultant_le_one).

  * Under randomisation (every assignment of the values to the units
    equally likely) Geary's C has expectation 1 (geary_mean_one) and the
    general G of Getis and Ord has expectation S0 / (n (n − 1))
    (getis_mean), the centres of gearc()'s and gorg()'s z statistics. Both
    rest on perm_pair_sum: over all permutations, a pair of distinct units
    sees every ordered pair of distinct values equally often.

  Not proved here: the variances of C and G, and the nearest-neighbour
  expectations, which are statements about a Poisson process.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section
open Finset

namespace Crimlab.Spatial

variable {ι : Type*} [Fintype ι]

/-- The weighted sum of squared deviations of `x` from `c`. -/
def ssq (w x : ι → ℝ) (c : ℝ) : ℝ := ∑ i, w i * (x i - c) ^ 2

/-- The weighted mean of `x`. -/
def wmean (w x : ι → ℝ) : ℝ := (∑ i, w i * x i) / ∑ i, w i

lemma ssq_expand (w x : ι → ℝ) (a : ℝ) :
    ssq w x a = (∑ i, w i * x i ^ 2) - 2 * a * (∑ i, w i * x i) + a ^ 2 * ∑ i, w i := by
  unfold ssq
  rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl (fun i _ => by ring)

/-- The parallel-axis identity: deviations about any point are those about
the weighted mean plus the total weight times the squared offset. -/
theorem parallel_axis (w x : ι → ℝ) (hW : ∑ i, w i ≠ 0) (c : ℝ) :
    ssq w x c = ssq w x (wmean w x) + (∑ i, w i) * (wmean w x - c) ^ 2 := by
  rw [ssq_expand, ssq_expand]
  unfold wmean
  field_simp
  ring

/-- With non-negative weights the weighted mean minimises the weighted sum
of squared deviations. -/
theorem mean_minimises (w x : ι → ℝ) (hw : ∀ i, 0 ≤ w i) (hW : ∑ i, w i ≠ 0) (c : ℝ) :
    ssq w x (wmean w x) ≤ ssq w x c := by
  rw [parallel_axis w x hW c]
  have : 0 ≤ (∑ i, w i) * (wmean w x - c) ^ 2 :=
    mul_nonneg (Finset.sum_nonneg (fun i _ => hw i)) (sq_nonneg _)
  linarith

/-- In the plane: the weighted mean centre minimises the weighted sum of
squared distances. -/
theorem mean_centre_minimises (w x y : ι → ℝ) (hw : ∀ i, 0 ≤ w i) (hW : ∑ i, w i ≠ 0)
    (cx cy : ℝ) :
    ∑ i, w i * ((x i - wmean w x) ^ 2 + (y i - wmean w y) ^ 2) ≤
      ∑ i, w i * ((x i - cx) ^ 2 + (y i - cy) ^ 2) := by
  have hx := mean_minimises w x hw hW cx
  have hy := mean_minimises w y hw hW cy
  unfold ssq at hx hy
  have e1 : ∑ i, w i * ((x i - wmean w x) ^ 2 + (y i - wmean w y) ^ 2) =
      (∑ i, w i * (x i - wmean w x) ^ 2) + ∑ i, w i * (y i - wmean w y) ^ 2 := by
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl (fun i _ => by ring)
  have e2 : ∑ i, w i * ((x i - cx) ^ 2 + (y i - cy) ^ 2) =
      (∑ i, w i * (x i - cx) ^ 2) + ∑ i, w i * (y i - cy) ^ 2 := by
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl (fun i _ => by ring)
  rw [e1, e2]
  linarith

/-- The local Moran statistics `I_i = z_i / m₂ ∑_j w_ij z_j`, with
`m₂ = ∑ z² / n`, add up to `S0` times global Moran's
`I = (n / S0) ∑∑ w_ij z_i z_j / ∑ z²`. -/
theorem lisa_sum (W : ι → ι → ℝ) (z : ι → ℝ) (n S0 : ℝ) (hn : n ≠ 0) (hS : S0 ≠ 0)
    (hz : ∑ i, z i ^ 2 ≠ 0) :
    ∑ i, z i / ((∑ j, z j ^ 2) / n) * ∑ j, W i j * z j =
      S0 * ((n / S0) * (∑ i, ∑ j, W i j * z i * z j) / ∑ i, z i ^ 2) := by
  have h : ∑ i, z i / ((∑ j, z j ^ 2) / n) * ∑ j, W i j * z j =
      n / (∑ j, z j ^ 2) * ∑ i, ∑ j, W i j * z i * z j := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [Finset.mul_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl (fun j _ => ?_)
    field_simp
  rw [h]
  field_simp

/-- The mean resultant length of weighted directions is at most 1. -/
theorem resultant_le_one (w θ : ι → ℝ) (hw : ∀ i, 0 ≤ w i) (hW : 0 < ∑ i, w i) :
    ‖∑ i, (w i : ℂ) * Complex.exp (θ i * Complex.I)‖ / ∑ i, w i ≤ 1 := by
  rw [div_le_one hW]
  calc ‖∑ i, (w i : ℂ) * Complex.exp (θ i * Complex.I)‖
      ≤ ∑ i, ‖(w i : ℂ) * Complex.exp (θ i * Complex.I)‖ := norm_sum_le _ _
    _ = ∑ i, w i := by
        refine Finset.sum_congr rfl (fun i _ => ?_)
        rw [norm_mul, Complex.norm_exp_ofReal_mul_I, mul_one, Complex.norm_real,
          Real.norm_of_nonneg (hw i)]

/-- The mean resultant length is not negative. -/
theorem resultant_nonneg (w θ : ι → ℝ) (hW : 0 < ∑ i, w i) :
    0 ≤ ‖∑ i, (w i : ℂ) * Complex.exp (θ i * Complex.I)‖ / ∑ i, w i :=
  div_nonneg (norm_nonneg _) hW.le

/-! ### Averages over all permutations of the values

Under randomisation every assignment of the n values to the n units is
equally likely: the expectation of a statistic is its average over the
permutations σ of the units. -/

section Perm

variable [DecidableEq ι]

/-- The sum over permutations of `f (σ i) (σ j)` is the same for every pair
of distinct units. -/
lemma perm_pair_eq (f : ι → ι → ℝ) {i j k l : ι} (hij : i ≠ j) (hkl : k ≠ l) :
    ∑ σ : Equiv.Perm ι, f (σ i) (σ j) = ∑ σ : Equiv.Perm ι, f (σ k) (σ l) := by
  set s1 := Equiv.swap k i
  have h1 : s1 l ≠ i := by
    intro h
    apply hkl
    have : s1 l = s1 k := by rw [h]; exact (Equiv.swap_apply_left k i).symm
    exact (s1.injective this).symm
  set τ := Equiv.swap (s1 l) j * s1
  have hk : τ k = i := by
    simp only [τ, Equiv.Perm.mul_apply, s1, Equiv.swap_apply_left]
    exact Equiv.swap_apply_of_ne_of_ne (Ne.symm h1) hij
  have hl : τ l = j := by
    simp only [τ, Equiv.Perm.mul_apply, Equiv.swap_apply_left]
  calc ∑ σ : Equiv.Perm ι, f (σ i) (σ j)
      = ∑ σ : Equiv.Perm ι, f ((σ * τ) k) ((σ * τ) l) := by
        simp only [Equiv.Perm.mul_apply, hk, hl]
    _ = ∑ σ : Equiv.Perm ι, f (σ k) (σ l) :=
        Equiv.sum_comp (Equiv.mulRight τ) (fun σ => f (σ k) (σ l))

/-- For distinct units `i`, `j`: `n (n − 1)` times the sum over
permutations of `f (σ i) (σ j)` is the number of permutations times the
sum of `f` over ordered pairs of distinct values. -/
theorem perm_pair_sum (f : ι → ι → ℝ) {i j : ι} (hij : i ≠ j) :
    (Fintype.card ι : ℝ) * ((Fintype.card ι : ℝ) - 1) * ∑ σ : Equiv.Perm ι, f (σ i) (σ j) =
      (Fintype.card (Equiv.Perm ι) : ℝ) * ((∑ a, ∑ b, f a b) - ∑ a, f a a) := by
  set S := ∑ σ : Equiv.Perm ι, f (σ i) (σ j)
  -- each unit k against the others l ≠ k gives S
  have row : ∀ k : ι, (∑ l, ∑ σ : Equiv.Perm ι, f (σ k) (σ l)) - ∑ σ : Equiv.Perm ι, f (σ k) (σ k)
      = ((Fintype.card ι : ℝ) - 1) * S := by
    intro k
    rw [← Finset.sum_erase_eq_sub (Finset.mem_univ k)]
    rw [Finset.sum_congr rfl (fun l hl => (perm_pair_eq f hij (Finset.ne_of_mem_erase hl).symm).symm)]
    rw [Finset.sum_const, Finset.card_erase_of_mem (Finset.mem_univ k), Finset.card_univ, nsmul_eq_mul]
    have : 1 ≤ Fintype.card ι := Fintype.card_pos_iff.mpr ⟨i⟩
    push_cast [Nat.cast_sub this]
    ring
  have lhs : (∑ k, ((∑ l, ∑ σ : Equiv.Perm ι, f (σ k) (σ l)) - ∑ σ : Equiv.Perm ι, f (σ k) (σ k)))
      = (Fintype.card ι : ℝ) * ((Fintype.card ι : ℝ) - 1) * S := by
    rw [Finset.sum_congr rfl (fun k _ => row k), Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
    ring
  have rhs : (∑ k, ((∑ l, ∑ σ : Equiv.Perm ι, f (σ k) (σ l)) - ∑ σ : Equiv.Perm ι, f (σ k) (σ k)))
      = ∑ σ : Equiv.Perm ι, ((∑ a, ∑ b, f a b) - ∑ a, f a a) := by
    calc (∑ k, ((∑ l, ∑ σ : Equiv.Perm ι, f (σ k) (σ l)) - ∑ σ : Equiv.Perm ι, f (σ k) (σ k)))
        = ∑ k, ∑ σ : Equiv.Perm ι, ((∑ l, f (σ k) (σ l)) - f (σ k) (σ k)) := by
          refine Finset.sum_congr rfl (fun k _ => ?_)
          rw [Finset.sum_sub_distrib, Finset.sum_comm]
      _ = ∑ σ : Equiv.Perm ι, ∑ k, ((∑ l, f (σ k) (σ l)) - f (σ k) (σ k)) := Finset.sum_comm
      _ = ∑ σ : Equiv.Perm ι, ((∑ a, ∑ b, f a b) - ∑ a, f a a) := by
          refine Finset.sum_congr rfl (fun σ _ => ?_)
          rw [Finset.sum_sub_distrib]
          congr 1
          · rw [Equiv.sum_comp σ (fun a => ∑ l, f a (σ l))]
            exact Finset.sum_congr rfl (fun a _ => Equiv.sum_comp σ (fun b => f a b))
          · exact Equiv.sum_comp σ (fun a => f a a)
  rw [← lhs, rhs, Finset.sum_const, Finset.card_univ, nsmul_eq_mul]

/-- Weighted sums over pairs of units, averaged over permutations: with no
unit its own neighbour (`W i i = 0`), the permutations' sum of
`∑ W_ij g(σ i, σ j)` is `S0` times the number of permutations times the
mean of `g` over ordered pairs of distinct values. -/
lemma perm_weighted (W : ι → ι → ℝ) (hW : ∀ i, W i i = 0) (g : ι → ι → ℝ)
    (h2 : 2 ≤ Fintype.card ι) :
    ∑ σ : Equiv.Perm ι, ∑ i, ∑ j, W i j * g (σ i) (σ j) =
      (∑ i, ∑ j, W i j) * ((Fintype.card (Equiv.Perm ι) : ℝ) * ((∑ a, ∑ b, g a b) - ∑ a, g a a) /
        ((Fintype.card ι : ℝ) * ((Fintype.card ι : ℝ) - 1))) := by
  have hN : (Fintype.card ι : ℝ) * ((Fintype.card ι : ℝ) - 1) ≠ 0 := by
    have : (2 : ℝ) ≤ Fintype.card ι := by exact_mod_cast h2
    apply mul_ne_zero <;> linarith
  set c := (Fintype.card (Equiv.Perm ι) : ℝ) * ((∑ a, ∑ b, g a b) - ∑ a, g a a) /
        ((Fintype.card ι : ℝ) * ((Fintype.card ι : ℝ) - 1))
  calc ∑ σ : Equiv.Perm ι, ∑ i, ∑ j, W i j * g (σ i) (σ j)
      = ∑ i, ∑ σ : Equiv.Perm ι, ∑ j, W i j * g (σ i) (σ j) := Finset.sum_comm
    _ = ∑ i, ∑ j, ∑ σ : Equiv.Perm ι, W i j * g (σ i) (σ j) :=
        Finset.sum_congr rfl (fun i _ => Finset.sum_comm)
    _ = ∑ i, ∑ j, W i j * c := by
        refine Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => ?_))
        rw [← Finset.mul_sum]
        by_cases h : i = j
        · subst h
          simp [hW i]
        · congr 1
          have := perm_pair_sum g h
          rw [show c = _ from rfl, eq_div_iff hN]
          linear_combination this
    _ = (∑ i, ∑ j, W i j) * c := by
        rw [Finset.sum_mul]
        exact Finset.sum_congr rfl (fun i _ => (Finset.sum_mul _ _ _).symm)

/-- Geary's C has expectation 1 under randomisation: its average over all
permutations of the values is 1 (n ≥ 2, no unit its own neighbour, a
positive total weight S0, values not all equal). -/
theorem geary_mean_one (W : ι → ι → ℝ) (hW : ∀ i, W i i = 0) (x : ι → ℝ)
    (h2 : 2 ≤ Fintype.card ι) (hS : ∑ i, ∑ j, W i j ≠ 0)
    (hv : ∑ a, (x a - (∑ b, x b) / Fintype.card ι) ^ 2 ≠ 0) :
    (∑ σ : Equiv.Perm ι,
        ((Fintype.card ι : ℝ) - 1) * (∑ i, ∑ j, W i j * (x (σ i) - x (σ j)) ^ 2) /
          (2 * (∑ i, ∑ j, W i j) * ∑ a, (x a - (∑ b, x b) / Fintype.card ι) ^ 2)) /
      (Fintype.card (Equiv.Perm ι) : ℝ) = 1 := by
  set V := ∑ a, (x a - (∑ b, x b) / (Fintype.card ι : ℝ)) ^ 2
  have hN2 : (2 : ℝ) ≤ (Fintype.card ι : ℝ) := by exact_mod_cast h2
  have hN : (Fintype.card ι : ℝ) ≠ 0 := by linarith
  have hN1 : (Fintype.card ι : ℝ) - 1 ≠ 0 := by linarith
  have hP : (Fintype.card (Equiv.Perm ι) : ℝ) ≠ 0 := by
    exact_mod_cast (Fintype.card_pos (α := Equiv.Perm ι)).ne'
  -- the pair sum of squared differences is 2 (Fintype.card ι : ℝ) V
  have hT : (∑ a, ∑ b, (x a - x b) ^ 2) - ∑ a, (x a - x a) ^ 2 = 2 * (Fintype.card ι : ℝ) * V := by
    simp only [sub_self, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow,
      Finset.sum_const_zero, sub_zero]
    have e1 : ∀ a, ∑ b, (x a - x b) ^ 2 =
        (Fintype.card ι : ℝ) * x a ^ 2 - 2 * x a * ∑ b, x b + ∑ b, x b ^ 2 := by
      intro a
      have h : ∑ b, (x a - x b) ^ 2 = ∑ b, (x a ^ 2 - 2 * x a * x b + x b ^ 2) :=
        Finset.sum_congr rfl (fun b _ => by ring)
      rw [h, Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum,
        Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
    have e2 : V = ∑ b, x b ^ 2 - (∑ b, x b) ^ 2 / (Fintype.card ι : ℝ) := by
      simp only [V]
      have : ∀ a, (x a - (∑ b, x b) / (Fintype.card ι : ℝ)) ^ 2 =
          x a ^ 2 - 2 * ((∑ b, x b) / (Fintype.card ι : ℝ)) * x a + ((∑ b, x b) / (Fintype.card ι : ℝ)) ^ 2 := fun a => by ring
      rw [Finset.sum_congr rfl (fun a _ => this a), Finset.sum_add_distrib,
        Finset.sum_sub_distrib, ← Finset.mul_sum, Finset.sum_const, Finset.card_univ,
        nsmul_eq_mul]
      field_simp
      ring
    rw [Finset.sum_congr rfl (fun a _ => e1 a), Finset.sum_add_distrib,
      Finset.sum_sub_distrib, ← Finset.mul_sum, ← Finset.sum_mul, Finset.sum_const,
      Finset.card_univ, nsmul_eq_mul, e2]
    field_simp
    rw [← Finset.mul_sum]
    ring
  have hsum := perm_weighted W hW (fun a b => (x a - x b) ^ 2) h2
  rw [← Finset.sum_div, ← Finset.mul_sum, hsum, hT]
  field_simp

/-- The general G of Getis and Ord has expectation `S0 / (n (n − 1))` under
randomisation (n ≥ 2, no unit its own neighbour, `∑_{a ≠ b} x_a x_b ≠ 0`). -/
theorem getis_mean (W : ι → ι → ℝ) (hW : ∀ i, W i i = 0) (x : ι → ℝ)
    (h2 : 2 ≤ Fintype.card ι)
    (hD : (∑ a, ∑ b, x a * x b) - ∑ a, x a * x a ≠ 0) :
    (∑ σ : Equiv.Perm ι, (∑ i, ∑ j, W i j * (x (σ i) * x (σ j))) /
        ((∑ a, ∑ b, x a * x b) - ∑ a, x a * x a)) /
      (Fintype.card (Equiv.Perm ι) : ℝ) =
      (∑ i, ∑ j, W i j) / ((Fintype.card ι : ℝ) * ((Fintype.card ι : ℝ) - 1)) := by
  have hN2 : (2 : ℝ) ≤ Fintype.card ι := by exact_mod_cast h2
  have hN : (Fintype.card ι : ℝ) * ((Fintype.card ι : ℝ) - 1) ≠ 0 := by
    apply mul_ne_zero <;> linarith
  have hP : (Fintype.card (Equiv.Perm ι) : ℝ) ≠ 0 := by
    exact_mod_cast (Fintype.card_pos (α := Equiv.Perm ι)).ne'
  have hD' : (∑ a, ∑ b, x a * x b) - ∑ a, x a ^ 2 ≠ 0 := by simpa [pow_two] using hD
  have hsum := perm_weighted W hW (fun a b => x a * x b) h2
  rw [← Finset.sum_div, hsum]
  field_simp

end Perm

end Crimlab.Spatial
