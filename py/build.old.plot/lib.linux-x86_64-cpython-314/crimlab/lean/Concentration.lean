/-
  crimlab: inst/lean/Concentration.lean

  The Gini coefficient of counts split at zero (cconc()).

  With D the sum over ordered pairs of units of |x_i - x_j| and T > 0 the
  total, the Gini coefficient of n units is D / (2 n T). Write the units
  as z at zero and k > 0 above it. A pair of zeros adds nothing to D, and
  a zero against a positive unit adds that unit's value twice (once in
  each order), so D = D_pos + 2 z T, with D_pos the sum over the positive
  units alone.

    * gini_zero_split: the Gini coefficient is s + (1 - s) G_pos, with
      s = z / (z + k) the share of units at zero and G_pos = D_pos / (2 k
      T) the Gini coefficient of the positive units.
    * gini_excess: against a null share of zeros s0 (exp(-mean) for
      Poisson counts) at the same G_pos, the excess of the observed Gini
      is (s - s0)(1 - G_pos): all of it comes from the zeros.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

namespace Crimlab.Concentration

/-- The Gini coefficient of z zeros and k positive units, from the positive
units' pairwise sum D and the total T. -/
theorem gini_zero_split (z k D T : ℝ) (hz : 0 ≤ z) (hk : 0 < k) (hT : 0 < T) :
    (D + 2 * z * T) / (2 * (z + k) * T) =
      z / (z + k) + (1 - z / (z + k)) * (D / (2 * k * T)) := by
  have h1 : z + k ≠ 0 := by linarith
  have h2 : k ≠ 0 := hk.ne'
  have h3 : T ≠ 0 := hT.ne'
  field_simp
  ring

/-- The split at two shares of zeros, the Gini of the rest held fixed. -/
theorem gini_excess (s s0 g : ℝ) :
    (s + (1 - s) * g) - (s0 + (1 - s0) * g) = (s - s0) * (1 - g) := by
  ring

end Crimlab.Concentration
