/-
  crimlab: inst/lean/Hazard.lean

  The built-in selection of period-by-period hazard ratios (Hernán 2010;
  Hernán & Robins, Causal Inference: What If, Fine Point 17.2).

  Two types of people: a share s "high-risk" with per-period
  reconviction probability h, the rest "low-risk" with probability l < h.
  The period-2 hazard among survivors is a weighted average of h and l
  whose weight on h is the surviving high-risk share (survivorHazard); it
  is increasing in that share (survivor_hazard_mono). An arm that depletes
  the high-risk type less in period 1 therefore has the higher period-2
  hazard (hr2_gt_one_of_depletion), with no individual effect in period 2
  at all; hr2_witness gives numbers. Period-by-period hazard ratios are
  not causal contrasts; risk differences at a fixed horizon are.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section

namespace Crimlab.Hazard

/-- Period-2 hazard among survivors when a share `w` of survivors is high-risk. -/
def survivorHazard (h l w : ℝ) : ℝ := w * h + (1 - w) * l

theorem survivor_hazard_mono (h l w₁ w₂ : ℝ) (hl : l < h) (hw : w₁ < w₂) :
    survivorHazard h l w₁ < survivorHazard h l w₂ := by
  unfold survivorHazard; nlinarith

/-- Surviving high-risk share after period 1 when high-risk people survive with probability `a`
and low-risk with probability `b`, from an initial high-risk share `s`. -/
def survivingShare (s a b : ℝ) : ℝ := s * a / (s * a + (1 - s) * b)

/-- If an arm depletes the high-risk type less (`a₂ > a₁`, same `b`), its surviving high-risk
share is larger, hence its period-2 hazard is larger, hence the period-2 hazard ratio exceeds one
with no individual effect in period 2 at all. -/
theorem hr2_gt_one_of_depletion (h l s a₁ a₂ b : ℝ) (hl : l < h) (hs0 : 0 < s) (hs1 : s < 1)
    (ha1 : 0 < a₁) (ha2 : a₁ < a₂) (hb : 0 < b) :
    survivorHazard h l (survivingShare s a₁ b) < survivorHazard h l (survivingShare s a₂ b) := by
  apply survivor_hazard_mono h l _ _ hl
  unfold survivingShare
  have h1s : 0 < 1 - s := by linarith
  have d1 : 0 < s * a₁ + (1 - s) * b := by have := mul_pos hs0 ha1; have := mul_pos h1s hb; linarith
  have d2 : 0 < s * a₂ + (1 - s) * b := by have := mul_pos hs0 (by linarith : 0 < a₂); have := mul_pos h1s hb; linarith
  rw [div_lt_div_iff₀ d1 d2]
  nlinarith [mul_pos hs0 h1s, mul_pos (mul_pos hs0 h1s) hb, mul_pos hs0 (mul_pos h1s hb)]

/-- A concrete witness: s = 1/2, h = 1/2, l = 1/10, treated arm survives high-risk at 4/5 vs 1/2,
low-risk at 9/10 in both: period-2 hazards 17/70 ≈ 0.243 (control) vs 49/170 ≈ 0.288 (treated), a
"harmful" ratio 343/289 ≈ 1.19 with zero effect. -/
theorem hr2_witness :
    survivorHazard (1/2) (1/10) (survivingShare (1/2) (1/2) (9/10)) < survivorHazard (1/2) (1/10) (survivingShare (1/2) (4/5) (9/10)) := by
  unfold survivorHazard survivingShare; norm_num

end Crimlab.Hazard
