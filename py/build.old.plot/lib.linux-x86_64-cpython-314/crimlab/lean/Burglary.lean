/-
  crimlab: inst/lean/Burglary.lean

  The Short-D'Orsogna-Brantingham model of sdbsim() and sdbthr(). With
  no burglaries the lattice update b' = (1 - omega dt) ((1 - eta) b +
  eta / 4 (sum of four neighbours)) of non-negative values at most M
  stays at most (1 - omega dt) M. In the continuum model the uniform
  steady state b = eps D gam / omega, A = a0 + b, rho = gam / A is a
  rest point, and the trace of the linearisation at k = 0, -omega +
  eps D rho - A, is negative; with eps = 0 the determinant's middle
  coefficient Q is negative, so no mode grows.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section

namespace Crimlab.Burglary

/-- Without burglaries the largest attractiveness shrinks by 1 - omega dt. -/
theorem sup_contract (w eta m b n1 n2 n3 n4 : ℝ) (hw0 : 0 ≤ w) (he0 : 0 ≤ eta) (he1 : eta ≤ 1)
    (hb : b ≤ m) (h1 : n1 ≤ m) (h2 : n2 ≤ m) (h3 : n3 ≤ m) (h4 : n4 ≤ m) :
    w * ((1 - eta) * b + eta / 4 * (n1 + n2 + n3 + n4)) ≤ w * m := by
  apply mul_le_mul_of_nonneg_left _ hw0
  nlinarith

/-- The uniform steady state is a rest point of both equations. -/
theorem steady (a0 om eps d gam : ℝ) (ha : 0 < a0) (ho : 0 < om) (he : 0 ≤ eps) (hd : 0 < d)
    (hg : 0 < gam) :
    -om * (eps * d * gam / om) + eps * d * (gam / (a0 + eps * d * gam / om)) * (a0 + eps * d * gam / om) = 0 ∧
      -(gam / (a0 + eps * d * gam / om)) * (a0 + eps * d * gam / om) + gam = 0 := by
  have hb : 0 ≤ eps * d * gam / om := by positivity
  have hA : a0 + eps * d * gam / om ≠ 0 := by linarith
  constructor
  · field_simp
    ring
  · field_simp
    ring

/-- The trace at k = 0 is negative (and it falls as k grows). -/
theorem trace_neg (a0 om eps d gam : ℝ) (ha : 0 < a0) (ho : 0 < om) (he : 0 ≤ eps) (hd : 0 < d)
    (hg : 0 < gam) :
    -om + eps * d * (gam / (a0 + eps * d * gam / om)) - (a0 + eps * d * gam / om) < 0 := by
  set bb := eps * d * gam / om with hbb
  have hb : 0 ≤ bb := by positivity
  have hA : 0 < a0 + bb := by linarith
  have key : eps * d * gam = om * bb := by rw [hbb]; field_simp
  have : eps * d * (gam / (a0 + bb)) = om * bb / (a0 + bb) := by
    rw [← key]; ring
  rw [this]
  have hlt : om * bb / (a0 + bb) < om := by
    rw [div_lt_iff₀ hA]
    nlinarith
  linarith

/-- Without repeat reinforcement Q is negative: no hot spots. -/
theorem no_hotspots_without_reinforcement (a0 om d gam eta z : ℝ) (ha : 0 < a0) (ho : 0 < om)
    (hd : 0 < d) (he : 0 ≤ eta) (hz : 0 < z) :
    d / z * (3 * 0 * d * (gam / a0) - om - eta * a0) < 0 := by
  have : 3 * 0 * d * (gam / a0) - om - eta * a0 < 0 := by nlinarith [mul_nonneg he ha.le]
  exact mul_neg_of_pos_of_neg (div_pos hd hz) this

end Crimlab.Burglary
