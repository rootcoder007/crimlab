/-
  crimlab: inst/lean/Dark.lean

  The dark figure of crime as a bound problem (crcap(), dkfig()).

    * petersen_bounds: two lists of sizes n1, n2 with overlap m > 0
      estimate a population N = theta n1 n2 / m, theta the factor by which
      the lists' dependence moves the Lincoln-Petersen estimate. With theta
      in [1/kappa, kappa], N lies in [n1 n2 / (kappa m), kappa n1 n2 / m].
    * rate_lower, rate_upper: a survey measures a true rate v in [0, 1]
      with false positives alpha <= amax < 1 and false negatives beta <=
      bmax < 1, so that it reports vobs = (1 - beta) v + alpha (1 - v).
      Then (vobs - amax) / (1 - amax) <= v and v <= vobs / (1 - bmax):
      every admissible error leaves the true rate inside the interval
      dkfig() reports.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

namespace Crimlab.Dark

/-- The capture-recapture estimate over a box of dependence. -/
theorem petersen_bounds (n1 n2 m kappa theta : ℝ) (hp : 0 ≤ n1 * n2 / m)
    (h1 : 1 / kappa ≤ theta) (h2 : theta ≤ kappa) :
    n1 * n2 / (kappa * m) ≤ theta * (n1 * n2 / m) ∧ theta * (n1 * n2 / m) ≤ kappa * (n1 * n2 / m) := by
  constructor
  · have : n1 * n2 / (kappa * m) = 1 / kappa * (n1 * n2 / m) := by
      rw [div_mul_div_comm, one_mul]
    rw [this]
    exact mul_le_mul_of_nonneg_right h1 hp
  · exact mul_le_mul_of_nonneg_right h2 hp

/-- False positives at most amax bound the true rate from below. -/
theorem rate_lower (v vobs alpha beta amax : ℝ) (hv0 : 0 ≤ v) (hv1 : v ≤ 1)
    (ha : alpha ≤ amax) (hm : amax < 1) (hb : 0 ≤ beta)
    (hobs : vobs = (1 - beta) * v + alpha * (1 - v)) :
    (vobs - amax) / (1 - amax) ≤ v := by
  rw [div_le_iff₀ (by linarith)]
  nlinarith [mul_nonneg hb hv0, mul_le_mul_of_nonneg_right ha (by linarith : (0 : ℝ) ≤ 1 - v)]

/-- False negatives at most bmax bound the true rate from above. -/
theorem rate_upper (v vobs alpha beta bmax : ℝ) (hv0 : 0 ≤ v) (hv1 : v ≤ 1)
    (ha0 : 0 ≤ alpha) (hb : beta ≤ bmax) (hm : bmax < 1)
    (hobs : vobs = (1 - beta) * v + alpha * (1 - v)) :
    v ≤ vobs / (1 - bmax) := by
  rw [le_div_iff₀ (by linarith)]
  nlinarith [mul_le_mul_of_nonneg_right hb hv0, mul_nonneg ha0 (by linarith : (0 : ℝ) ≤ 1 - v)]

end Crimlab.Dark
