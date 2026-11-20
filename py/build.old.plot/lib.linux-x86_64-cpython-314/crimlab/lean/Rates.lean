/-
  crimlab: inst/lean/Rates.lean

  The algebra of rate comparisons. With y events over exposure e in a
  group and y0 over e0 in the reference group, the rate ratio is
  (y / e) / (y0 / e0): for the two-group log-linear model, cntreg()'s
  estimate. Swapping the groups inverts the ratio, and inverts its Wald
  interval exp(b ± z s), with the ends exchanged; the rate difference
  changes sign (grprr()). A directly standardised rate, a weighted mean of
  the stratum rates with weights from 0 summing to 1, lies between the
  smallest and the largest of them (dastd()). The empirical Bayes ratio
  (o + m a) / (e + a) is the convex combination of the observed ratio o / e
  and the prior mean m with weights e / (e + a) and a / (e + a) (sirex()).
  Division is Lean's (x / 0 = 0), so the inversions hold for every real
  count.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section

namespace Crimlab.Rates

/-- The rate ratio of a group (y events over exposure e) against a
reference group (y0 over e0). -/
def rateRatio (y e y0 e0 : ℝ) : ℝ := (y / e) / (y0 / e0)

/-- The groups swapped: the rate ratio inverts. -/
theorem rr_swap (y e y0 e0 : ℝ) : rateRatio y0 e0 y e = (rateRatio y e y0 e0)⁻¹ := by
  unfold rateRatio
  rw [inv_div]

/-- The groups swapped (the log ratio b becomes -b): the lower end of the
interval is the inverse of the upper end before the swap. -/
theorem ci_swap (b z s : ℝ) : Real.exp (-b - z * s) = (Real.exp (b + z * s))⁻¹ := by
  rw [← Real.exp_neg]
  ring_nf

/-- The groups swapped: the rate difference changes sign. -/
theorem rd_swap (y e y0 e0 : ℝ) : y0 / e0 - y / e = -(y / e - y0 / e0) := by
  ring

/-- A directly standardised rate is at most the largest stratum rate. -/
theorem dsr_le_max {k : ℕ} (w r : Fin k → ℝ) (M : ℝ) (hw : ∀ i, 0 ≤ w i)
    (hs : ∑ i, w i = 1) (hM : ∀ i, r i ≤ M) : ∑ i, w i * r i ≤ M := by
  calc ∑ i, w i * r i ≤ ∑ i, w i * M :=
        Finset.sum_le_sum fun i _ => mul_le_mul_of_nonneg_left (hM i) (hw i)
    _ = M := by rw [← Finset.sum_mul, hs, one_mul]

/-- A directly standardised rate is at least the smallest stratum rate. -/
theorem dsr_ge_min {k : ℕ} (w r : Fin k → ℝ) (m : ℝ) (hw : ∀ i, 0 ≤ w i)
    (hs : ∑ i, w i = 1) (hm : ∀ i, m ≤ r i) : m ≤ ∑ i, w i * r i := by
  calc m = ∑ i, w i * m := by rw [← Finset.sum_mul, hs, one_mul]
    _ ≤ ∑ i, w i * r i :=
        Finset.sum_le_sum fun i _ => mul_le_mul_of_nonneg_left (hm i) (hw i)

/-- The empirical Bayes ratio is a convex combination of the observed
ratio and the prior mean. -/
theorem eb_convex (o e m a : ℝ) (he : 0 < e) (ha : 0 < a) :
    (o + m * a) / (e + a) = e / (e + a) * (o / e) + a / (e + a) * m := by
  have h : e + a ≠ 0 := by linarith
  field_simp

end Crimlab.Rates
