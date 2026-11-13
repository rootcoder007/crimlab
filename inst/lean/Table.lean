/-
  crimlab: inst/lean/Table.lean

  The algebra of a 2 x 2 table behind ctabx(). With cells a, b (first row)
  and c, d (second row), positive margins and total n, Pearson's statistic
  sum (o - e)^2 / e over the four cells, e the product of the margins over
  n, equals n (a d - b c)^2 over the product of the four margins: phi^2
  times n, so Cramer's V of a 2 x 2 table is |phi|. Swapping the rows
  leaves the statistic unchanged.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section

namespace Crimlab.Table

/-- Pearson's chi-square statistic of the 2 x 2 table [[a, b], [c, d]]. -/
def chi2 (a b c d : ℝ) : ℝ :=
  let n := a + b + c + d
  let e := fun (r k : ℝ) => r * k / n
  (a - e (a + b) (a + c)) ^ 2 / e (a + b) (a + c) +
  (b - e (a + b) (b + d)) ^ 2 / e (a + b) (b + d) +
  (c - e (c + d) (a + c)) ^ 2 / e (c + d) (a + c) +
  (d - e (c + d) (b + d)) ^ 2 / e (c + d) (b + d)

/-- The closed form of the 2 x 2 statistic. -/
theorem chi2_closed (a b c d : ℝ) (h1 : 0 < a + b) (h2 : 0 < c + d) (h3 : 0 < a + c)
    (h4 : 0 < b + d) :
    chi2 a b c d = (a + b + c + d) * (a * d - b * c) ^ 2 /
      ((a + b) * (c + d) * (a + c) * (b + d)) := by
  have hn : 0 < a + b + c + d := by linarith
  unfold chi2
  field_simp
  ring

/-- Swapping the rows leaves the statistic unchanged. -/
theorem chi2_swap (a b c d : ℝ) : chi2 c d a b = chi2 a b c d := by
  unfold chi2
  ring_nf

end Crimlab.Table
