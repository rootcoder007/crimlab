/-
  crimlab: inst/lean/Rates.lean

  The algebra of rate comparisons. With y events over exposure e in a
  group and y0 over e0 in the reference group, the rate ratio is
  (y / e) / (y0 / e0): for the two-group log-linear model, cntreg()'s
  estimate. Swapping the groups inverts the ratio, and inverts its Wald
  interval exp(b ± z s), with the ends exchanged. Division is Lean's
  (x / 0 = 0), so the inversions hold for every real count.
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

end Crimlab.Rates
