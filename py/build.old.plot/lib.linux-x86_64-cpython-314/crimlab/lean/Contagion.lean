/-
  crimlab: inst/lean/Contagion.lean

  The branching arithmetic of a self-exciting (Hawkes) process, as
  hwkfit() reports it in `branch`. Each event has on average eta
  offspring; a cluster started by one background event has eta^k
  members in generation k, so 1 / (1 - eta) in all when eta < 1. A
  stationary process with background rate nu has the rate r that
  balances r = nu + eta r, so r = nu / (1 - eta), and the share of its
  events that were triggered, (r - nu) / r, is eta.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section

namespace Crimlab.Contagion

/-- The generations of a cluster sum to 1 / (1 - eta). -/
theorem cluster_size (eta : ℝ) (h0 : 0 ≤ eta) (h1 : eta < 1) :
    ∑' k : ℕ, eta ^ k = (1 - eta)⁻¹ :=
  tsum_geometric_of_lt_one h0 h1

/-- The stationary rate: the one solution of r = nu + eta r. -/
theorem stationary_rate (nu eta r : ℝ) (h1 : eta < 1) (hr : r = nu + eta * r) :
    r = nu / (1 - eta) := by
  have h : (1 - eta) ≠ 0 := by linarith
  field_simp
  linarith

/-- The share of events triggered is the branching ratio. -/
theorem triggered_share (nu eta : ℝ) (hn : 0 < nu) (h1 : eta < 1) :
    (nu / (1 - eta) - nu) / (nu / (1 - eta)) = eta := by
  have h : (1 - eta) ≠ 0 := by linarith
  have hn' : nu ≠ 0 := ne_of_gt hn
  field_simp
  ring

end Crimlab.Contagion
