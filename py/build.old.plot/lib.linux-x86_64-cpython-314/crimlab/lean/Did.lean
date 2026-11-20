/-
  crimlab: inst/lean/Did.lean

  The 2 x 2 difference in differences. With a group g (0 or 1) and a
  period t (0 or 1), the two-way fixed-effects regression of an outcome on
  a constant, the group, the period and the treatment g t has four
  coefficients for four cells: it fits the cell means exactly, so its
  treatment coefficient is the difference of the differences of the cell
  means (didfe()'s estimate on such a panel; the self-test of
  src/t_did.cpp checks the same number). Adding any unit effect and any
  period effect to the outcome leaves the difference of differences
  unchanged, and on outcomes made of such effects plus an effect tau on
  the treated cell it equals tau.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

namespace Crimlab.Did

/-- The difference of differences of the four cell means: y g t is the
mean of group g in period t. -/
def dd (y : Fin 2 → Fin 2 → ℝ) : ℝ := (y 1 1 - y 1 0) - (y 0 1 - y 0 0)

/-- A saturated fit a + b g + c t + d g t that reproduces the four cell
means has d equal to their difference of differences. -/
theorem twfe_2x2 (y : Fin 2 → Fin 2 → ℝ) (a b c d : ℝ)
    (h00 : a = y 0 0) (h10 : a + b = y 1 0) (h01 : a + c = y 0 1)
    (h11 : a + b + c + d = y 1 1) : d = dd y := by
  unfold dd
  linarith

/-- Unit and period effects cancel in the difference of differences. -/
theorem dd_effects (y : Fin 2 → Fin 2 → ℝ) (u m : Fin 2 → ℝ) :
    dd (fun g t => y g t + u g + m t) = dd y := by
  unfold dd
  ring

/-- Outcomes made of unit and period effects and an effect tau on the
treated cell give tau. -/
theorem dd_effect (u m : Fin 2 → ℝ) (tau : ℝ) :
    dd (fun g t => u g + m t + if g = 1 ∧ t = 1 then tau else 0) = tau := by
  unfold dd
  simp
  ring

end Crimlab.Did
