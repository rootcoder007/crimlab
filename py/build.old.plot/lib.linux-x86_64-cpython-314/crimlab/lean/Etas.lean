/-
  crimlab: inst/lean/Etas.lean

  etasf()'s modified Omori kernel K (d + c)^-p is, with no marks, the
  Lomax kernel of hwkfit() scaled by the branching ratio: with p = a + 1
  and K = eta a c^a, K (d + c)^-(a+1) = eta (a / c) (1 + d / c)^-(a+1),
  so the two fits describe one process.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section

namespace Crimlab.Etas

/-- The Omori kernel is the Lomax density times the branching ratio. -/
theorem omori_lomax (eta a c d : ℝ) (hc : 0 < c) (hd : 0 ≤ d) :
    eta * a * c ^ a * (d + c) ^ (-(a + 1)) = eta * (a / c * (1 + d / c) ^ (-(a + 1))) := by
  have hdc : 0 < d + c := by linarith
  have h1 : 1 + d / c = (d + c) / c := by field_simp; ring
  rw [h1, Real.div_rpow hdc.le hc.le, Real.rpow_neg hc.le, Real.rpow_add_one hc.ne']
  have hca : 0 < c ^ a := Real.rpow_pos_of_pos hc a
  field_simp

end Crimlab.Etas
