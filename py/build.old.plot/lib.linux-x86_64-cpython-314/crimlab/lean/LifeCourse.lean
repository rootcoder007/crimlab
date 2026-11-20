/-
  crimlab: inst/lean/LifeCourse.lean

  The age-crime curve of a mixture of types (agemix()).

  With types k of shares w_k (summing to 1) and curves x_k, y_k at two
  ages, the aggregate is sum w_k x_k at each age, and the covariance of
  an individual's offending at the two ages, given only the type, is the
  covariance across types sum_k w_k (x_k - xbar)(y_k - ybar).

    * one_type_cov_zero: when every type has the same curve (one type, in
      effect) the covariance is 0, so a covariance away from 0 in panel
      data refutes a single type, whatever the aggregate curve.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

open Finset

namespace Crimlab.LifeCourse

/-- One curve shared by every type: no covariance between ages. -/
theorem one_type_cov_zero {n : ℕ} (w : Fin n → ℝ) (x0 y0 : ℝ) (hw : ∑ k, w k = 1) :
    ∑ k, w k * ((x0 - ∑ j, w j * x0) * (y0 - ∑ j, w j * y0)) = 0 := by
  have hx : ∑ j, w j * x0 = x0 := by rw [← sum_mul, hw, one_mul]
  have hy : ∑ j, w j * y0 = y0 := by rw [← sum_mul, hw, one_mul]
  simp [hx, hy]

end Crimlab.LifeCourse
