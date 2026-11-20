/-
  crimlab: inst/lean/RelChk.lean

  The algebra relchk() rests on: how the contrasts of a 2 x 2 table change
  when it is relabelled. A group with a events and b non-events is compared
  with a reference group with c events and d non-events:
      OR = a d / (b c),  RR = (a / (a + b)) / (c / (c + d)),
      RD = a / (a + b) − c / (c + d).
  Swapping the reference group, or which outcome is the event, inverts the
  odds ratio; swapping both leaves it unchanged, and so does transposing
  the table (groups read from the columns). The risk ratio inverts with the
  reference group; the risk difference changes sign with either swap.
  Division is Lean's (x / 0 = 0), so the odds-ratio and risk-ratio
  identities hold for every real count; relchk() reports a zero denominator
  as Inf or NaN instead, where these identities are not claimed.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section

namespace Crimlab.RelChk

/-- The odds ratio of a group (a events, b non-events) against a reference
group (c events, d non-events). -/
def oddsRatio (a b c d : ℝ) : ℝ := (a * d) / (b * c)

/-- The risk ratio of the same comparison. -/
def riskRatio (a b c d : ℝ) : ℝ := (a / (a + b)) / (c / (c + d))

/-- The risk difference of the same comparison. -/
def riskDiff (a b c d : ℝ) : ℝ := a / (a + b) - c / (c + d)

/-- Another reference group: the odds ratio inverts. -/
theorem or_swap_ref (a b c d : ℝ) : oddsRatio c d a b = (oddsRatio a b c d)⁻¹ := by
  unfold oddsRatio
  rw [inv_div]
  ring

/-- The other outcome as the event: the odds ratio inverts. -/
theorem or_swap_event (a b c d : ℝ) : oddsRatio b a d c = (oddsRatio a b c d)⁻¹ := by
  unfold oddsRatio
  rw [inv_div]

/-- Both swapped: the odds ratio is unchanged. -/
theorem or_swap_both (a b c d : ℝ) : oddsRatio d c b a = oddsRatio a b c d := by
  unfold oddsRatio
  ring

/-- The table transposed (groups in the columns): the odds ratio is
unchanged. -/
theorem or_transpose (a b c d : ℝ) : oddsRatio a c b d = oddsRatio a b c d := by
  unfold oddsRatio
  ring

/-- Another reference group: the risk ratio inverts. -/
theorem rr_swap_ref (a b c d : ℝ) : riskRatio c d a b = (riskRatio a b c d)⁻¹ := by
  unfold riskRatio
  rw [inv_div]

/-- Another reference group: the risk difference changes sign. -/
theorem rd_swap_ref (a b c d : ℝ) : riskDiff c d a b = -riskDiff a b c d := by
  unfold riskDiff
  ring

/-- The other outcome as the event: the risk difference changes sign, for
groups that are not empty. -/
theorem rd_swap_event (a b c d : ℝ) (hab : a + b ≠ 0) (hcd : c + d ≠ 0) :
    riskDiff b a d c = -riskDiff a b c d := by
  unfold riskDiff
  have h1 : b + a ≠ 0 := by rwa [add_comm]
  have h2 : d + c ≠ 0 := by rwa [add_comm]
  field_simp
  ring

end Crimlab.RelChk
