/-
  crimlab: inst/lean/Fair.lean

  The algebra behind impfpr(), nsbnd() and rnkres().

  A confusion table (tp, fp, fn, tn) of positive masses gives the base rate
  p, the positive predictive value ppv, the false positive rate fpr and the
  false negative rate fnr. Chouldechova (2017):
      fpr = p / (1 − p) · (1 − ppv) / ppv · (1 − fnr),
  so two groups with equal ppv and equal fnr but different base rates have
  different fpr, whatever the classifier.

  Label noise: with α = P(recorded 1 | true 0) and β = P(recorded 0 |
  true 1), a true base rate p is recorded as p (1 − β) + (1 − p) α. With
  α in [αlo, αhi] and β in [βlo, βhi], αhi + βhi < 1, the true rate lies
  between (pobs − αhi) / (1 − αhi − βlo) and (pobs − αlo) / (1 − αlo −
  βhi), and both ends are reached; nsbnd() clamps them to [0, 1], where p
  lies anyway. Disjoint intervals order two groups' true rates for every
  admissible noise; overlapping ones (boxes from 0) admit noise in either
  order.

  A ranking by estimates within δ of the true scores: below a true gap of
  2δ an admissible error reverses a pair, above it none can.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section

namespace Crimlab.Fair

/-- A confusion table for one group: masses of the four cells, all positive. -/
structure Table where
  tp : ℝ
  fp : ℝ
  fn : ℝ
  tn : ℝ
  tp_pos : 0 < tp
  fp_pos : 0 < fp
  fn_pos : 0 < fn
  tn_pos : 0 < tn

namespace Table

variable (T : Table)

/-- Total mass. -/
def total : ℝ := T.tp + T.fp + T.fn + T.tn
/-- Base rate. -/
def p : ℝ := (T.tp + T.fn) / T.total
/-- Positive predictive value. -/
def ppv : ℝ := T.tp / (T.tp + T.fp)
/-- False positive rate. -/
def fpr : ℝ := T.fp / (T.fp + T.tn)
/-- False negative rate. -/
def fnr : ℝ := T.fn / (T.tp + T.fn)

theorem total_pos : 0 < T.total := by
  unfold total; have := T.tp_pos; have := T.fp_pos; have := T.fn_pos; have := T.tn_pos; linarith

theorem p_lt_one : T.p < 1 := by
  unfold p; rw [div_lt_one T.total_pos]; unfold total; have := T.fp_pos; have := T.tn_pos; linarith

theorem ppv_pos : 0 < T.ppv := by unfold ppv; have := T.tp_pos; have := T.fp_pos; positivity

theorem one_sub_ppv : 1 - T.ppv = T.fp / (T.tp + T.fp) := by
  unfold ppv; have := T.tp_pos; have := T.fp_pos; field_simp; ring

theorem one_sub_fnr : 1 - T.fnr = T.tp / (T.tp + T.fn) := by
  unfold fnr; have := T.tp_pos; have := T.fn_pos; field_simp; ring

theorem one_sub_p : 1 - T.p = (T.fp + T.tn) / T.total := by
  unfold p total; have := T.tp_pos; have := T.fp_pos; have := T.fn_pos; have := T.tn_pos
  field_simp; ring

/-- Chouldechova's identity: `fpr = p/(1-p) · (1-ppv)/ppv · (1-fnr)`. -/
theorem chouldechova : T.fpr = T.p / (1 - T.p) * ((1 - T.ppv) / T.ppv) * (1 - T.fnr) := by
  rw [T.one_sub_ppv, T.one_sub_fnr, T.one_sub_p]
  unfold fpr p ppv total
  have := T.tp_pos; have := T.fp_pos; have := T.fn_pos; have := T.tn_pos
  field_simp

end Table

/-- Equal ppv and equal fnr with different base rates force different fpr. -/
theorem impossibility (A B : Table) (hppv : A.ppv = B.ppv) (hfnr : A.fnr = B.fnr)
    (hp : A.p ≠ B.p) : A.fpr ≠ B.fpr := by
  intro hfpr
  rw [A.chouldechova, B.chouldechova, hppv, hfnr] at hfpr
  have hmid : 0 < (1 - B.ppv) / B.ppv := by
    rw [B.one_sub_ppv]; have := B.tp_pos; have := B.fp_pos; have := B.ppv_pos; positivity
  have hright : 0 < 1 - B.fnr := by
    rw [B.one_sub_fnr]; have := B.tp_pos; have := B.fn_pos; positivity
  have hodds : A.p / (1 - A.p) = B.p / (1 - B.p) := by
    have hk : (1 - B.ppv) / B.ppv * (1 - B.fnr) ≠ 0 := by positivity
    have := hfpr
    rw [mul_assoc, mul_assoc] at this
    exact mul_right_cancel₀ hk this
  have hA1 : 0 < 1 - A.p := by have := A.p_lt_one; linarith
  have hB1 : 0 < 1 - B.p := by have := B.p_lt_one; linarith
  rw [div_eq_div_iff hA1.ne' hB1.ne'] at hodds
  apply hp
  linarith

/-! ### Label noise -/

/-- The recorded rate of a true base rate `p` under noise rates `α`, `β`. -/
def observedRate (p α β : ℝ) : ℝ := p * (1 - β) + (1 - p) * α

/-- The true rate a recorded one gives under known noise rates. -/
def trueRate (pobs α β : ℝ) : ℝ := (pobs - α) / (1 - α - β)

/-- Known noise rates invert the recording exactly. -/
theorem trueRate_observedRate (p α β : ℝ) (h : 1 - α - β ≠ 0) :
    trueRate (observedRate p α β) α β = p := by
  unfold trueRate observedRate; field_simp; ring

/-- Noise rates in boxes `[αlo, αhi]`, `[βlo, βhi]` with `αhi + βhi < 1`:
the true rate lies between the two corner inversions. -/
theorem bounds_box (p α β αlo αhi βlo βhi : ℝ) (ha : αlo ≤ α) (ha' : α ≤ αhi)
    (hb : βlo ≤ β) (hb' : β ≤ βhi) (hsum : αhi + βhi < 1) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    (observedRate p α β - αhi) / (1 - αhi - βlo) ≤ p ∧
      p ≤ (observedRate p α β - αlo) / (1 - αlo - βhi) := by
  have h1 : 0 < 1 - αhi - βlo := by linarith
  have h2 : 0 < 1 - αlo - βhi := by linarith
  constructor
  · rw [div_le_iff₀ h1]
    unfold observedRate
    nlinarith [mul_nonneg hp0 (sub_nonneg.mpr hb), mul_nonneg (sub_nonneg.mpr hp1) (sub_nonneg.mpr ha')]
  · rw [le_div_iff₀ h2]
    unfold observedRate
    nlinarith [mul_nonneg hp0 (sub_nonneg.mpr hb'), mul_nonneg (sub_nonneg.mpr hp1) (sub_nonneg.mpr ha)]

/-- The lower end is reached at `α = αhi`, `β = βlo`. -/
theorem lower_attained (p αhi βlo : ℝ) (h : 1 - αhi - βlo ≠ 0) :
    (observedRate p αhi βlo - αhi) / (1 - αhi - βlo) = p :=
  trueRate_observedRate p αhi βlo h

/-- The upper end is reached at `α = αlo`, `β = βhi`. -/
theorem upper_attained (p αlo βhi : ℝ) (h : 1 - αlo - βhi ≠ 0) :
    (observedRate p αlo βhi - αlo) / (1 - αlo - βhi) = p :=
  trueRate_observedRate p αlo βhi h

/-- Disjoint intervals order the true rates for every admissible noise. -/
theorem compare_decided (pA pB αA βA αB βB αlo αhi βlo βhi : ℝ)
    (hαA : αlo ≤ αA) (hαA' : αA ≤ αhi) (hβA : βlo ≤ βA) (hβA' : βA ≤ βhi)
    (hαB : αlo ≤ αB) (hαB' : αB ≤ αhi) (hβB : βlo ≤ βB) (hβB' : βB ≤ βhi)
    (hsum : αhi + βhi < 1) (hA0 : 0 ≤ pA) (hA1 : pA ≤ 1) (hB0 : 0 ≤ pB) (hB1 : pB ≤ 1)
    (hgap : (observedRate pA αA βA - αlo) / (1 - αlo - βhi) <
      (observedRate pB αB βB - αhi) / (1 - αhi - βlo)) :
    pA < pB := by
  obtain ⟨_, hAhi⟩ := bounds_box pA αA βA αlo αhi βlo βhi hαA hαA' hβA hβA' hsum hA0 hA1
  obtain ⟨hBlo, _⟩ := bounds_box pB αB βB αlo αhi βlo βhi hαB hαB' hβB hβB' hsum hB0 hB1
  linarith

/-- Overlapping intervals (boxes from 0) admit noise under which either true
rate is the larger, reproducing the same recorded rates `qA`, `qB`. -/
theorem compare_undecided (qA qB αb βb : ℝ) (hαb : 0 ≤ αb) (hβb : 0 ≤ βb)
    (hsum : αb + βb < 1)
    (hover1 : (qA - αb) / (1 - αb) < qB / (1 - βb))
    (hover2 : (qB - αb) / (1 - αb) < qA / (1 - βb)) :
    (∃ pA pB : ℝ, observedRate pA αb 0 = qA ∧ observedRate pB 0 βb = qB ∧ pA < pB) ∧
      (∃ pA pB : ℝ, observedRate pA 0 βb = qA ∧ observedRate pB αb 0 = qB ∧ pB < pA) := by
  have h1 : 1 - αb ≠ 0 := (by linarith : (0 : ℝ) < 1 - αb).ne'
  have h2 : 1 - βb ≠ 0 := (by linarith : (0 : ℝ) < 1 - βb).ne'
  constructor
  · refine ⟨(qA - αb) / (1 - αb), qB / (1 - βb), ?_, ?_, hover1⟩
    · unfold observedRate; field_simp; ring
    · unfold observedRate; field_simp; ring
  · refine ⟨qA / (1 - βb), (qB - αb) / (1 - αb), ?_, ?_, hover2⟩
    · unfold observedRate; field_simp; ring
    · unfold observedRate; field_simp; ring

/-! ### The resolution of a noisy ranking -/

/-- Below `2δ` an admissible error reverses the order of two true scores. -/
theorem rank_reversal_exists (t₁ t₂ δ : ℝ) (hδ : 0 < δ) (_h12 : t₁ < t₂) (hgap : t₂ - t₁ < 2 * δ) :
    ∃ e₁ e₂ : ℝ, |e₁| ≤ δ ∧ |e₂| ≤ δ ∧ t₂ + e₂ < t₁ + e₁ := by
  refine ⟨δ, -δ, by rw [abs_of_pos hδ], by rw [abs_neg, abs_of_pos hδ], ?_⟩
  linarith

/-- Above `2δ` the order survives every admissible error. -/
theorem rank_stable_of_gap (t₁ t₂ δ e₁ e₂ : ℝ) (hgap : 2 * δ < t₂ - t₁) (h1 : |e₁| ≤ δ)
    (h2 : |e₂| ≤ δ) : t₁ + e₁ < t₂ + e₂ := by
  have := abs_le.1 h1
  have := abs_le.1 h2
  linarith [this.1, this.2]

end Crimlab.Fair
