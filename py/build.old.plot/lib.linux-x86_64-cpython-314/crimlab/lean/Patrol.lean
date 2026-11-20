/-
  crimlab: inst/lean/Patrol.lean

  The square-root law of sqrtlw(): the mean distance to the nearest of u
  available units in an area A is c sqrt(A / u). Scaling the area and
  the units by the same factor leaves it unchanged, and four times the
  units halve it.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section

namespace Crimlab.Patrol

/-- The mean response distance of the square-root law. -/
def dist (c a u : ℝ) : ℝ := c * Real.sqrt (a / u)

/-- Area and units scaled together: the distance is unchanged. -/
theorem scale_invariant (c a u k : ℝ) (hk : k ≠ 0) : dist c (k * a) (k * u) = dist c a u := by
  unfold dist
  rw [mul_div_mul_left a u hk]

/-- Four times the units: half the distance. -/
theorem four_units_half (c a u : ℝ) (ha : 0 ≤ a) (hu : 0 ≤ u) : dist c a (4 * u) = dist c a u / 2 := by
  unfold dist
  have h4 : Real.sqrt 4 = 2 := by
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num]
    exact Real.sqrt_sq (by norm_num)
  rw [show a / (4 * u) = a / u / 4 by ring, Real.sqrt_div (div_nonneg ha hu), h4]
  ring

end Crimlab.Patrol
