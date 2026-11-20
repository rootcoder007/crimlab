/-
  crimlab: inst/lean/Judge.lean

  The judge-leniency instrument (jdgiv(), jmono()).

  With shares pc of compliers and pd of defiers whose mean effects are tc
  and td, the intention-to-treat effect is pc tc - pd td and the first
  stage pc - pd, so the Wald ratio is their quotient.

    * wald_no_defiers: with pd = 0 and pc > 0 the ratio is tc, the
      compliers' effect.
    * complier_effect: whatever pd, a first stage f = pc - pd != 0 and pc
      != 0 give tc = (wald f + pd td) / pc, the sensitivity jdgiv()
      reports.
    * slope_bound: nested judges j, k with leniencies differing by m >= 0
      differ in mean outcome by m times the mean effect e of the marginal
      cases; with outcomes in [lo, hi], |e| <= hi - lo, so the outcomes
      differ by at most (hi - lo) m. A larger difference refutes the
      nesting.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

namespace Crimlab.Judge

/-- No defiers: the Wald ratio is the compliers' mean effect. -/
theorem wald_no_defiers (pc tc td : ℝ) (hc : 0 < pc) :
    (pc * tc - 0 * td) / (pc - 0) = tc := by
  rw [zero_mul, sub_zero, sub_zero]
  exact mul_div_cancel_left₀ tc hc.ne'

/-- The compliers' effect from the Wald ratio, the first stage and the
defiers. -/
theorem complier_effect (pc pd tc td : ℝ) (hf : pc - pd ≠ 0) (hc : pc ≠ 0) :
    tc = ((pc * tc - pd * td) / (pc - pd) * (pc - pd) + pd * td) / pc := by
  field_simp
  ring

/-- Nested judges: the outcome difference is at most (hi - lo) times the
difference in leniency. -/
theorem slope_bound (m e lo hi : ℝ) (hm : 0 ≤ m) (he : |e| ≤ hi - lo) :
    |m * e| ≤ (hi - lo) * m := by
  rw [abs_mul, abs_of_nonneg hm, mul_comm]
  exact mul_le_mul_of_nonneg_right he hm

end Crimlab.Judge
