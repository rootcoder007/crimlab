/-
  crimlab: inst/lean/Published.lean

  What two published cells allow about the change between them (pbbnd(),
  yoypx()). A cell rounded or suppressed stands for a true count n in an
  interval [nlo, nhi]; the earlier one for p in [plo, phi].

    * chg_envelope: the true change n - p lies in [nlo - phi, nhi - plo]
      (its ends are the changes at n = nlo, p = phi and at n = nhi, p =
      plo).
    * pct_envelope: with 0 < plo and 0 <= nlo, the true ratio n / p lies in
      [nlo / phi, nhi / plo], so the percent change lies in the interval
      pbbnd() reports.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

namespace Crimlab.Published

/-- The change between two values known only within intervals. -/
theorem chg_envelope (n p nlo nhi plo phi : ℝ) (h1 : nlo ≤ n) (h2 : n ≤ nhi)
    (h3 : plo ≤ p) (h4 : p ≤ phi) : nlo - phi ≤ n - p ∧ n - p ≤ nhi - plo := by
  constructor <;> linarith

/-- The ratio between two counts known only within intervals. -/
theorem pct_envelope (n p nlo nhi plo phi : ℝ) (h0 : 0 ≤ nlo) (hp : 0 < plo) (h1 : nlo ≤ n)
    (h2 : n ≤ nhi) (h3 : plo ≤ p) (h4 : p ≤ phi) : nlo / phi ≤ n / p ∧ n / p ≤ nhi / plo := by
  have hp' : 0 < p := lt_of_lt_of_le hp h3
  constructor
  · calc nlo / phi ≤ nlo / p := div_le_div_of_nonneg_left h0 hp' h4
      _ ≤ n / p := div_le_div_of_nonneg_right h1 hp'.le
  · calc n / p ≤ nhi / p := div_le_div_of_nonneg_right h2 hp'.le
      _ ≤ nhi / plo := div_le_div_of_nonneg_left (le_trans h0 (le_trans h1 h2)) hp h3

end Crimlab.Published
