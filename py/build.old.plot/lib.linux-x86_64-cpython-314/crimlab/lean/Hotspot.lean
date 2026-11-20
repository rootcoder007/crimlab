/-
  crimlab: inst/lean/Hotspot.lean

  Identities and bounds behind the spatial crime methods.

  * Empirical-Bayes shrinkage (ebshr()): with noise mass Se ≥ 0 and signal
    mass St ≥ 0, the loss of shrinking by B is (1 − B)² Se + B² St. The
    factor B* = Se / (Se + St) lies in [0, 1] (bstar_mem), minimises the
    loss (bstar_min, through loss_gap: the loss exceeds its minimum by
    (Se + St) (B − B*)²), and never does worse than the raw counts, B = 0
    (bstar_le_raw).
  * The scan statistic (kscan()): the Poisson log-likelihood ratio
    c log(c / e) + (N − c) log((N − c) / (N − e)) of a cylinder holding
    more than it expects is never negative (llr_nonneg).
  * The dissimilarity index (dissix()): between 0 and 1 for non-negative
    counts with positive totals (dissim_nonneg, dissim_le_one).
  * The effects of a lag model (spill()): when every row of W sums to 1 and
    S inverts I − ρ W from the left, the total effect of a covariate,
    S (β I + θ W) applied to the vector of ones, is (β + θ) / (1 − ρ) at
    every unit (total_effect).
  * Bivariate Moran's I (bvmor()) of a variable with itself is Moran's I
    (bv_self).

  Not proved here: the moments of Moran's I, the null distributions of the
  permutation tests, and the likelihood's maximum.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section
open Finset

namespace Crimlab.Hotspot

/-- The loss of shrinking by `B`. -/
def loss (Se St B : ℝ) : ℝ := (1 - B) ^ 2 * Se + B ^ 2 * St

/-- The shrinkage factor that minimises it. -/
def bstar (Se St : ℝ) : ℝ := Se / (Se + St)

theorem bstar_mem (Se St : ℝ) (he : 0 ≤ Se) (ht : 0 ≤ St) (hp : 0 < Se + St) :
    0 ≤ bstar Se St ∧ bstar Se St ≤ 1 := by
  unfold bstar
  refine ⟨div_nonneg he hp.le, ?_⟩
  rw [div_le_one hp]
  linarith

/-- The loss exceeds its value at `B*` by `(Se + St) (B − B*)²`. -/
theorem loss_gap (Se St B : ℝ) (hp : 0 < Se + St) :
    loss Se St B = loss Se St (bstar Se St) + (Se + St) * (B - bstar Se St) ^ 2 := by
  unfold loss bstar
  field_simp
  ring

theorem bstar_min (Se St B : ℝ) (hp : 0 < Se + St) :
    loss Se St (bstar Se St) ≤ loss Se St B := by
  rw [loss_gap Se St B hp]
  have : 0 ≤ (Se + St) * (B - bstar Se St) ^ 2 := mul_nonneg hp.le (sq_nonneg _)
  linarith

theorem bstar_le_raw (Se St : ℝ) (hp : 0 < Se + St) :
    loss Se St (bstar Se St) ≤ Se := by
  have h := bstar_min Se St 0 hp
  unfold loss at h ⊢
  simpa using h

/-- The scan's log-likelihood ratio is never negative when the cylinder
holds more than it expects. -/
theorem llr_nonneg (c e N : ℝ) (he : 0 < e) (hce : e < c) (hcN : c ≤ N) :
    0 ≤ c * Real.log (c / e) + (N - c) * Real.log ((N - c) / (N - e)) := by
  have hc : 0 < c := lt_trans he hce
  -- log x ≥ 1 − 1 / x
  have h1 : c * (1 - (c / e)⁻¹) ≤ c * Real.log (c / e) :=
    mul_le_mul_of_nonneg_left (Real.one_sub_inv_le_log_of_pos (div_pos hc he)) hc.le
  have h1' : c * (1 - (c / e)⁻¹) = c - e := by
    field_simp
  rcases eq_or_lt_of_le hcN with h | h
  · subst h
    simp only [sub_self, zero_mul, add_zero]
    linarith
  · have hNc : 0 < N - c := by linarith
    have hNe : 0 < N - e := by linarith
    have h2 : (N - c) * (1 - ((N - c) / (N - e))⁻¹) ≤ (N - c) * Real.log ((N - c) / (N - e)) :=
      mul_le_mul_of_nonneg_left (Real.one_sub_inv_le_log_of_pos (div_pos hNc hNe)) hNc.le
    have h2' : (N - c) * (1 - ((N - c) / (N - e))⁻¹) = e - c := by
      field_simp
      ring
    linarith

variable {ι : Type*} [Fintype ι]

/-- Duncan's dissimilarity of two groups' counts. -/
def dissim (a b : ι → ℝ) : ℝ :=
  (∑ i, |a i / (∑ j, a j) - b i / (∑ j, b j)|) / 2

theorem dissim_nonneg (a b : ι → ℝ) : 0 ≤ dissim a b := by
  unfold dissim
  have : 0 ≤ ∑ i, |a i / (∑ j, a j) - b i / (∑ j, b j)| :=
    Finset.sum_nonneg (fun i _ => abs_nonneg _)
  linarith

theorem dissim_le_one (a b : ι → ℝ) (ha : ∀ i, 0 ≤ a i) (hb : ∀ i, 0 ≤ b i)
    (hA : 0 < ∑ j, a j) (hB : 0 < ∑ j, b j) : dissim a b ≤ 1 := by
  unfold dissim
  have step : ∀ i, |a i / (∑ j, a j) - b i / (∑ j, b j)| ≤
      a i / (∑ j, a j) + b i / (∑ j, b j) := by
    intro i
    have u : 0 ≤ a i / (∑ j, a j) := div_nonneg (ha i) hA.le
    have v : 0 ≤ b i / (∑ j, b j) := div_nonneg (hb i) hB.le
    rw [abs_le]
    constructor <;> linarith
  have hs := Finset.sum_le_sum (fun i (_ : i ∈ (Finset.univ : Finset ι)) => step i)
  rw [Finset.sum_add_distrib, ← Finset.sum_div, ← Finset.sum_div, div_self hA.ne',
    div_self hB.ne'] at hs
  linarith

/-- The total effect of a lag model with row-standardised weights. -/
theorem total_effect {n : Type*} [Fintype n] [DecidableEq n] (W S : Matrix n n ℝ)
    (ρ β θ : ℝ) (hW : W.mulVec (fun _ => (1 : ℝ)) = fun _ => 1)
    (hS : S * (1 - ρ • W) = 1) (hρ : ρ ≠ 1) :
    S.mulVec ((β • (1 : Matrix n n ℝ) + θ • W).mulVec (fun _ => (1 : ℝ))) =
      fun _ => (β + θ) / (1 - ρ) := by
  set o : n → ℝ := fun _ => 1 with ho
  have h1 : (1 - ρ) ≠ 0 := sub_ne_zero.mpr (Ne.symm hρ)
  have hA : (1 - ρ • W).mulVec o = (1 - ρ) • o := by
    rw [Matrix.sub_mulVec, Matrix.one_mulVec, Matrix.smul_mulVec, hW]
    funext i
    simp [ho]
  have hSo : S.mulVec o = (1 - ρ)⁻¹ • o := by
    have : S.mulVec ((1 - ρ • W).mulVec o) = o := by
      rw [Matrix.mulVec_mulVec, hS, Matrix.one_mulVec]
    rw [hA, Matrix.mulVec_smul] at this
    calc S.mulVec o = (1 - ρ)⁻¹ • ((1 - ρ) • S.mulVec o) := by
          rw [smul_smul, inv_mul_cancel₀ h1, one_smul]
      _ = (1 - ρ)⁻¹ • o := by rw [this]
  have hB : (β • (1 : Matrix n n ℝ) + θ • W).mulVec o = (β + θ) • o := by
    rw [Matrix.add_mulVec, Matrix.smul_mulVec, Matrix.smul_mulVec,
      Matrix.one_mulVec, hW]
    funext i
    simp [ho]
  rw [hB, Matrix.mulVec_smul, hSo, smul_smul]
  funext i
  simp [ho]
  field_simp

/-- Bivariate Moran's I of a variable with itself is Moran's I. -/
theorem bv_self (n S0 Q V : ℝ) (hV : 0 ≤ V) :
    n / S0 * Q / Real.sqrt (V * V) = n / S0 * Q / V := by
  rw [Real.sqrt_mul_self hV]

end Crimlab.Hotspot
