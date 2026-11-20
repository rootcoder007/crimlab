/-
  crimlab: inst/lean/Series.lean

  Two facts the series functions rest on. xcser()'s correlations are
  sums of products of centred values over the square root of the product
  of their sums of squares: by Cauchy-Schwarz they lie in [-1, 1].
  trndrr()'s Mann-Kendall S is a sum of signs of differences: reversing
  the sign of every value (a rise read as a fall) negates S, so the test
  treats rises and falls alike.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section

namespace Crimlab.Series

open Finset

/-- Cauchy-Schwarz for the sums behind a correlation. -/
theorem cross_sq_le {ι : Type*} (s : Finset ι) (a b : ι → ℝ) :
    (∑ i ∈ s, a i * b i) ^ 2 ≤ (∑ i ∈ s, a i ^ 2) * ∑ i ∈ s, b i ^ 2 :=
  Finset.sum_mul_sq_le_sq_mul_sq s a b

/-- A correlation lies between -1 and 1. -/
theorem corr_le_one {ι : Type*} (s : Finset ι) (a b : ι → ℝ)
    (ha : 0 < ∑ i ∈ s, a i ^ 2) (hb : 0 < ∑ i ∈ s, b i ^ 2) :
    |(∑ i ∈ s, a i * b i) / Real.sqrt ((∑ i ∈ s, a i ^ 2) * ∑ i ∈ s, b i ^ 2)| ≤ 1 := by
  have hp : 0 < Real.sqrt ((∑ i ∈ s, a i ^ 2) * ∑ i ∈ s, b i ^ 2) :=
    Real.sqrt_pos.mpr (mul_pos ha hb)
  rw [abs_div, abs_of_pos hp, div_le_one hp]
  exact Real.abs_le_sqrt (cross_sq_le s a b)

/-- Negating every difference negates Mann-Kendall's S. -/
theorem mk_neg {ι : Type*} (s : Finset ι) (d : ι → ℝ) :
    ∑ i ∈ s, Real.sign (-d i) = -∑ i ∈ s, Real.sign (d i) := by
  rw [← Finset.sum_neg_distrib]
  exact Finset.sum_congr rfl (fun i _ => Real.sign_neg)

end Crimlab.Series
