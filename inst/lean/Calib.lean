/-
  crimlab: inst/lean/Calib.lean

  The pool-adjacent-violators pass of `clab_lrcal` (lrcal, cllrmin, eceplt).
  Blocks are (same-source count, size); the pass pushes the blocks of tied
  ratios in increasing order and pools the top of its stack into the new
  block while the top's share is at least the new one's.

    * push_pos, push_total: sizes stay positive and pooling conserves counts.
    * pav_increasing, above_iff: the pooled shares strictly increase from
      block to block (so the calibrated values are non-decreasing in the
      ratio).
    * pool_weighted_mean, pool_between: a pooled block's share is the
      size-weighted mean of the two, and lies between them.
    * gibbs, block_mean_optimal: the cost the kernel adds up for one block
      is lowest at the block's own share s/c among all constant values.

  Not proved: that the whole pass minimises the cost over every monotone
  recalibration (the optimality of PAV), and hence not that cllrmin's
  `min` is at most `cllr`.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section

namespace Crimlab.Calib

/-- One block: (same-source count, size). -/
abbrev Blk := ℝ × ℝ

/-- The kernel's stack step (`clab_lrcal`): the new block is pooled into the top while the
top's share is at least the new one's, tested as the kernel does, by cross-multiplying. The
stack is a list with the top first. -/
def push : List Blk → Blk → List Blk
  | [], b => [b]
  | t :: r, b => if b.1 * t.2 ≤ t.1 * b.2 then push r (t.1 + b.1, t.2 + b.2) else b :: t :: r

/-- The blocks after the whole pass, given the blocks of tied ratios in increasing order. -/
def pav (bs : List Blk) : List Blk := bs.foldl push []

/-- Top block's share strictly above the one below, cross-multiplied. -/
def Above (x y : Blk) : Prop := y.1 * x.2 < x.1 * y.2

theorem push_pos (st : List Blk) (b : Blk) (hs : ∀ x ∈ st, 0 < x.2) (hb : 0 < b.2) :
    ∀ x ∈ push st b, 0 < x.2 := by
  induction st generalizing b with
  | nil => simpa [push] using hb
  | cons t r ih =>
    simp only [push]
    split_ifs
    · exact ih _ (fun x hx => hs x (by simp [hx])) (by have := hs t (by simp); simp; linarith)
    · intro x hx
      simp only [List.mem_cons] at hx
      rcases hx with rfl | hx
      · exact hb
      · exact hs x (by simpa using hx)

theorem push_chain (st : List Blk) (b : Blk) (hc : st.IsChain Above) :
    (push st b).IsChain Above := by
  induction st generalizing b with
  | nil => simp [push]
  | cons t r ih =>
    simp only [push]
    split_ifs with h
    · exact ih _ hc.tail
    · refine List.isChain_cons_cons.2 ⟨?_, hc⟩
      unfold Above; exact not_le.mp h

theorem pav_aux (acc bs : List Blk) (ha : ∀ x ∈ acc, 0 < x.2) (hc : acc.IsChain Above)
    (hb : ∀ x ∈ bs, 0 < x.2) :
    (∀ x ∈ bs.foldl push acc, 0 < x.2) ∧ (bs.foldl push acc).IsChain Above := by
  induction bs generalizing acc with
  | nil => exact ⟨ha, hc⟩
  | cons b bs ih =>
    simp only [List.foldl_cons]
    exact ih _ (push_pos acc b ha (hb b (by simp))) (push_chain acc b hc)
      (fun x hx => hb x (by simp [hx]))

/-- The pooled block shares strictly fall from the top of the stack down: read in the order of
the ratios, the calibrated values strictly increase from one block to the next. -/
theorem pav_increasing (bs : List Blk) (hb : ∀ x ∈ bs, 0 < x.2) :
    (pav bs).IsChain Above ∧ ∀ x ∈ pav bs, 0 < x.2 :=
  (pav_aux [] bs (by simp) List.isChain_nil hb).symm

/-- With positive sizes, `Above x y` is `y`'s share below `x`'s. -/
theorem above_iff (x y : Blk) (hx : 0 < x.2) (hy : 0 < y.2) : Above x y ↔ y.1 / y.2 < x.1 / x.2 := by
  unfold Above
  rw [div_lt_div_iff₀ hy hx]

/-- Pooling conserves the counts: the stack's totals grow by exactly the block pushed. -/
theorem push_total (st : List Blk) (b : Blk) : (push st b).sum = st.sum + b := by
  induction st generalizing b with
  | nil => simp [push]
  | cons t r ih =>
    simp only [push]
    split_ifs
    · rw [ih]; simp only [List.sum_cons]; ext <;> simp <;> ring
    · simp only [List.sum_cons]; abel

/-- Pooling two adjacent blocks gives their size-weighted mean share. -/
theorem pool_weighted_mean (s₁ c₁ s₂ c₂ : ℝ) (h₁ : 0 < c₁) (h₂ : 0 < c₂) :
    (s₁ + s₂) / (c₁ + c₂) = (c₁ * (s₁ / c₁) + c₂ * (s₂ / c₂)) / (c₁ + c₂) := by
  field_simp

/-- ... which lies between the two shares. -/
theorem pool_between (s₁ c₁ s₂ c₂ : ℝ) (h₁ : 0 < c₁) (h₂ : 0 < c₂) (h : s₂ / c₂ ≤ s₁ / c₁) :
    s₂ / c₂ ≤ (s₁ + s₂) / (c₁ + c₂) ∧ (s₁ + s₂) / (c₁ + c₂) ≤ s₁ / c₁ := by
  rw [div_le_div_iff₀ h₂ h₁] at h
  constructor
  · rw [div_le_div_iff₀ h₂ (by linarith)]; nlinarith
  · rw [div_le_div_iff₀ (by linarith) h₁]; nlinarith

/-- Gibbs: for weights `a, b > 0` the weighted log-loss `a log P + b log (1 − P)` over
`P ∈ (0, 1)` is largest at `P = a / (a + b)`. -/
theorem gibbs (a b P : ℝ) (ha : 0 < a) (hb : 0 < b) (hP : 0 < P) (hP1 : P < 1) :
    a * Real.log P + b * Real.log (1 - P) ≤
      a * Real.log (a / (a + b)) + b * Real.log (b / (a + b)) := by
  have hab : 0 < a + b := by linarith
  have h1 := Real.log_le_sub_one_of_pos (show 0 < P * (a + b) / a by positivity)
  have h2 := Real.log_le_sub_one_of_pos (show 0 < (1 - P) * (a + b) / b by
    have : 0 < 1 - P := by linarith
    positivity)
  rw [Real.log_div (by positivity) ha.ne', Real.log_mul hP.ne' hab.ne'] at h1
  have h1P : 0 < 1 - P := by linarith
  rw [Real.log_div (mul_pos h1P hab).ne' hb.ne',
    Real.log_mul (by linarith) hab.ne'] at h2
  rw [Real.log_div ha.ne' hab.ne', Real.log_div hb.ne' hab.ne']
  have e1 : a * (P * (a + b) / a - 1) = P * (a + b) - a := by field_simp
  have e2 : b * ((1 - P) * (a + b) / b - 1) = (1 - P) * (a + b) - b := by field_simp
  have k1 := mul_le_mul_of_nonneg_left h1 ha.le
  have k2 := mul_le_mul_of_nonneg_left h2 hb.le
  rw [e1] at k1
  rw [e2] at k2
  nlinarith

/-- The kernel's cost of one block of `c` ratios, `s` of them same-source, all calibrated to
`pr`, with `ns` same-source and `nd` different-source ratios in all (the two sums of
`clab_lrcal` restricted to the block, before the common factor `1/(2 ln 2)`). -/
def blockCost (s c ns nd pr : ℝ) : ℝ :=
  s / ns * Real.log (1 + (1 - pr) / pr * (ns / nd)) + (c - s) / nd * Real.log (1 + pr / (1 - pr) * (nd / ns))

/-- The block's own share `s / c` is the best constant value for it: no other `pr ∈ (0, 1)`
gives that block a lower cost. -/
theorem block_mean_optimal (s c ns nd pr : ℝ) (hs : 0 < s) (hsc : s < c) (hns : 0 < ns)
    (hnd : 0 < nd) (hp : 0 < pr) (hp1 : pr < 1) :
    blockCost s c ns nd (s / c) ≤ blockCost s c ns nd pr := by
  have hc : 0 < c := by linarith
  -- the cost as a log-loss in P = pr nd / (pr nd + (1 - pr) ns)
  have key : ∀ x : ℝ, 0 < x → x < 1 →
      blockCost s c ns nd x = -(s / ns * Real.log (x * nd / (x * nd + (1 - x) * ns))
        + (c - s) / nd * Real.log (1 - x * nd / (x * nd + (1 - x) * ns))) := by
    intro x hx hx1
    have h1x : 0 < 1 - x := by linarith
    have hD : 0 < x * nd + (1 - x) * ns := by positivity
    unfold blockCost
    have e1 : 1 + (1 - x) / x * (ns / nd) = (x * nd / (x * nd + (1 - x) * ns))⁻¹ := by
      field_simp
    have e2 : 1 + x / (1 - x) * (nd / ns) = (1 - x * nd / (x * nd + (1 - x) * ns))⁻¹ := by
      field_simp; ring
    rw [e1, e2, Real.log_inv, Real.log_inv]
    ring
  have hP : ∀ x : ℝ, 0 < x → x < 1 → 0 < x * nd / (x * nd + (1 - x) * ns) ∧
      x * nd / (x * nd + (1 - x) * ns) < 1 := by
    intro x hx hx1
    have h1x : 0 < 1 - x := by linarith
    have hD : 0 < x * nd + (1 - x) * ns := by positivity
    refine ⟨by positivity, (div_lt_one hD).2 (by nlinarith)⟩
  have hsc' : 0 < s / c := by positivity
  have hsc1 : s / c < 1 := (div_lt_one hc).2 hsc
  rw [key _ hsc' hsc1, key _ hp hp1]
  have ha : 0 < s / ns := by positivity
  have hb : 0 < (c - s) / nd := div_pos (by linarith) hnd
  have g := gibbs (s / ns) ((c - s) / nd) _ ha hb (hP pr hp hp1).1 (hP pr hp hp1).2
  have hstar : s / c * nd / (s / c * nd + (1 - s / c) * ns) = s / ns / (s / ns + (c - s) / nd) := by
    field_simp
  have hne : s / ns + (c - s) / nd ≠ 0 := (add_pos ha hb).ne'
  have hstar2 : 1 - s / c * nd / (s / c * nd + (1 - s / c) * ns) = (c - s) / nd / (s / ns + (c - s) / nd) := by
    rw [hstar, one_sub_div hne]; congr 1; ring
  rw [hstar2, hstar]
  linarith

end Crimlab.Calib
