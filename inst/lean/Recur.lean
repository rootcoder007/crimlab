/-
  crimlab: inst/lean/Recur.lean

  The Kaplan-Meier curve of kmrec() is a running product of factors
  1 - d / n, with d events among n at risk (0 <= d <= n, n > 0). Each
  factor lies in [0, 1], so a curve that starts in [0, 1] stays there
  and never rises.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section

namespace Crimlab.Recur

/-- A factor of the product lies in [0, 1]. -/
theorem km_factor (d n : ℝ) (hd : 0 ≤ d) (hdn : d ≤ n) (hn : 0 < n) :
    0 ≤ 1 - d / n ∧ 1 - d / n ≤ 1 := by
  constructor
  · rw [sub_nonneg, div_le_one hn]
    exact hdn
  · have : 0 ≤ d / n := div_nonneg hd hn.le
    linarith

/-- One step of the curve: it stays non-negative and does not rise. -/
theorem km_step (s d n : ℝ) (hs : 0 ≤ s) (hd : 0 ≤ d) (hdn : d ≤ n) (hn : 0 < n) :
    0 ≤ s * (1 - d / n) ∧ s * (1 - d / n) ≤ s := by
  obtain ⟨h0, h1⟩ := km_factor d n hd hdn hn
  exact ⟨mul_nonneg hs h0, mul_le_of_le_one_right hs h1⟩

end Crimlab.Recur
