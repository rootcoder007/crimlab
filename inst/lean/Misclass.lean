/-
  crimlab: inst/lean/Misclass.lean

  Nondifferential misclassification of a binary group attenuates the risk
  difference by sensitivity + specificity − 1; an ordinal classification
  can flip its sign (Bross 1954; Hernán & Robins, Causal Inference: What
  If, Fine Point 9.1, Technical Point 9.1).

  With true binary group A, outcome risks r1 = P(Y | A = 1), r0 = P(Y | A =
  0), prevalence π, sensitivity se and specificity sp, the measured risks are
  mixtures and the measured risk difference is
      RD* = (se + sp − 1) · RD · π(1 − π) / (d1 d0)             (`binary_attenuation`),
  with d1, d0 the measured-group shares, so when se + sp ≥ 1 the measured
  difference has the sign of the true one or is zero (`binary_sign`).
  For any se, sp in [0, 1] the factor lies in [−1, 1] (`factor_abs_le_one`), so
  the measured difference is never larger in size than the true one
  (`measured_abs_le`), and below se + sp = 1 it has the opposite sign
  (`binary_flip`).
  For a three-level classification with risks (0, 1/2, 1) and a
  non-monotone misclassification that sends levels 0 and 2 to the measured
  "high" group and level 1 to the measured "low" group, the measured
  contrast is 1/2 − 1/2 = 0 while the true top-vs-bottom contrast is 1
  (`ordinal_witness`); with the groups arranged the other way the measured
  contrast is negative (`ordinal_reversal`). The binary theorem is a binary
  theorem.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section

namespace Crimlab.Misclass

/-- The measured risk difference under nondifferential misclassification of a binary
exposure, written out from the two mixtures: outcome mass over group share in each
measured group. -/
def measuredRD (r1 r0 π se sp : ℝ) : ℝ :=
  (π * se * r1 + (1 - π) * (1 - sp) * r0) / (π * se + (1 - π) * (1 - sp))
    - (π * (1 - se) * r1 + (1 - π) * sp * r0) / (π * (1 - se) + (1 - π) * sp)

/-- `RD* = (se + sp − 1) RD π(1 − π) / (d1 d0)` when both measured groups are non-empty. -/
theorem binary_attenuation (r1 r0 π se sp : ℝ) (hd1 : π * se + (1 - π) * (1 - sp) ≠ 0)
    (hd0 : π * (1 - se) + (1 - π) * sp ≠ 0) :
    measuredRD r1 r0 π se sp
      = (se + sp - 1) * (r1 - r0) * (π * (1 - π)) / ((π * se + (1 - π) * (1 - sp)) * (π * (1 - se) + (1 - π) * sp)) := by
  unfold measuredRD
  field_simp
  ring

/-- With `se + sp ≥ 1` the measured difference has the sign of the true one, or is zero. -/
theorem binary_sign (r1 r0 π se sp : ℝ) (hπ : 0 < π) (hπ1 : π < 1)
    (hvalid : 1 ≤ se + sp) (hd1 : 0 < π * se + (1 - π) * (1 - sp))
    (hd0 : 0 < π * (1 - se) + (1 - π) * sp) (hrd : 0 < r1 - r0) :
    0 ≤ measuredRD r1 r0 π se sp := by
  rw [binary_attenuation r1 r0 π se sp hd1.ne' hd0.ne']
  have h1 : 0 ≤ se + sp - 1 := by linarith
  have h2 : 0 ≤ π * (1 - π) := by nlinarith
  exact div_nonneg (mul_nonneg (mul_nonneg h1 hrd.le) h2) (mul_pos hd1 hd0).le

/-- Three true levels with risks `0, 1/2, 1`; the measured "high" group holds levels 0 and 2
(equal weight, mean risk `1/2`), the measured "low" group holds level 1 (risk `1/2`): the
measured contrast is `0` while the true top-vs-bottom contrast is `1`. -/
theorem ordinal_witness : ((0 : ℝ) + 1) / 2 - 1 / 2 = 0 ∧ (1 : ℝ) - 0 = 1 := by norm_num

/-- Level 2 alone in the measured "low" group (risk `1`), levels 0 and 1 in "high" (mean risk
`1/4`): the measured contrast `1/4 − 1` is negative while the true contrast is positive. -/
theorem ordinal_reversal : ((0 : ℝ) + 1 / 2) / 2 - 1 < 0 ∧ (0 : ℝ) < 1 - 0 := by norm_num

/-- The attenuation factor `(se + sp − 1) π(1 − π) / (d1 d0)` lies in `[−1, 1]` for any
sensitivity and specificity in `[0, 1]`. -/
theorem factor_abs_le_one (π se sp : ℝ) (hπ : 0 < π) (hπ1 : π < 1)
    (hse0 : 0 ≤ se) (hse1 : se ≤ 1) (hsp0 : 0 ≤ sp) (hsp1 : sp ≤ 1)
    (hd1 : 0 < π * se + (1 - π) * (1 - sp)) (hd0 : 0 < π * (1 - se) + (1 - π) * sp) :
    |(se + sp - 1) * (π * (1 - π)) / ((π * se + (1 - π) * (1 - sp)) * (π * (1 - se) + (1 - π) * sp))| ≤ 1 := by
  have hd := mul_pos hd1 hd0
  rw [abs_div, abs_of_pos hd, div_le_one hd, abs_le]
  have a1 : 0 ≤ π * π * (se * (1 - se)) := by positivity
  have a2 : 0 ≤ (1 - π) * (1 - π) * (sp * (1 - sp)) := by
    have : 0 ≤ 1 - π := by linarith
    have : 0 ≤ 1 - se := by linarith
    have : 0 ≤ 1 - sp := by linarith
    positivity
  have a3 : 0 ≤ π * (1 - π) * ((1 - se) * (1 - sp)) := by
    have : 0 ≤ 1 - π := by linarith
    have : 0 ≤ 1 - se := by linarith
    have : 0 ≤ 1 - sp := by linarith
    positivity
  have a4 : 0 ≤ π * (1 - π) * (se * sp) := by
    have : 0 ≤ 1 - π := by linarith
    positivity
  constructor <;> nlinarith [a1, a2, a3, a4]

/-- `|RD*| ≤ |RD|` for any sensitivity and specificity in `[0, 1]`. -/
theorem measured_abs_le (r1 r0 π se sp : ℝ) (hπ : 0 < π) (hπ1 : π < 1)
    (hse0 : 0 ≤ se) (hse1 : se ≤ 1) (hsp0 : 0 ≤ sp) (hsp1 : sp ≤ 1)
    (hd1 : 0 < π * se + (1 - π) * (1 - sp)) (hd0 : 0 < π * (1 - se) + (1 - π) * sp) :
    |measuredRD r1 r0 π se sp| ≤ |r1 - r0| := by
  have h := factor_abs_le_one π se sp hπ hπ1 hse0 hse1 hsp0 hsp1 hd1 hd0
  rw [binary_attenuation r1 r0 π se sp hd1.ne' hd0.ne']
  have e : (se + sp - 1) * (r1 - r0) * (π * (1 - π)) / ((π * se + (1 - π) * (1 - sp)) * (π * (1 - se) + (1 - π) * sp))
      = (r1 - r0) * ((se + sp - 1) * (π * (1 - π)) / ((π * se + (1 - π) * (1 - sp)) * (π * (1 - se) + (1 - π) * sp))) := by ring
  rw [e, abs_mul]
  calc |r1 - r0| * _ ≤ |r1 - r0| * 1 := mul_le_mul_of_nonneg_left h (abs_nonneg _)
    _ = |r1 - r0| := mul_one _

/-- Below `se + sp = 1` the sign flips: a positive true difference is measured negative. -/
theorem binary_flip (r1 r0 π se sp : ℝ) (hπ : 0 < π) (hπ1 : π < 1)
    (hinv : se + sp < 1) (hd1 : 0 < π * se + (1 - π) * (1 - sp))
    (hd0 : 0 < π * (1 - se) + (1 - π) * sp) (hrd : 0 < r1 - r0) :
    measuredRD r1 r0 π se sp < 0 := by
  rw [binary_attenuation r1 r0 π se sp hd1.ne' hd0.ne']
  have h2 : 0 < π * (1 - π) := by nlinarith
  apply div_neg_of_neg_of_pos _ (mul_pos hd1 hd0)
  have : (se + sp - 1) * (r1 - r0) < 0 := mul_neg_of_neg_of_pos (by linarith) hrd
  exact mul_neg_of_neg_of_pos this h2

end Crimlab.Misclass
