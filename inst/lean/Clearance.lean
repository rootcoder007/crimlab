/-
  csls: inst/lean/Clearance.lean

  Simpson's paradox in clearance rates, and Kitagawa's decomposition of a
  change in an aggregate rate (Kitagawa 1955).

  A police service records d_k offences of type k and clears n_k of them.
  Its aggregate clearance rate is ∑ n / ∑ d; with workload shares
  w_k = d_k / ∑ d and type rates r_k = n_k / d_k it is the share-weighted
  mean ∑ w_k r_k (agg_eq_weighted). Between two comparisons A and B:

    * kitagawa: the change ∑ w_B r_B − ∑ w_A r_A is exactly the rate
      effect ∑ (w_A + w_B)/2 · (r_B − r_A) plus the mix effect
      ∑ ((r_A + r_B)/2 − c) · (w_B − w_A), for every centre c, because the
      shares sum to one in both comparisons.
    * change_ge_mix: if no type's rate falls, the change is at least the
      mix effect; so an aggregate fall while every type clears better is
      the mix effect alone (workload moved toward types that clear less).
    * fixed_mix_no_reversal: with the workload mix held fixed, the
      aggregate cannot fall when no type's rate falls.
    * simpson_witness: numbers where both types clear better in B and the
      aggregate is lower in B.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section
open Finset

namespace Csls.Clearance

variable {ι : Type*} [Fintype ι]

/-- The aggregate rate: everything cleared over everything recorded. -/
def agg (n d : ι → ℝ) : ℝ := (∑ k, n k) / ∑ k, d k

/-- The aggregate rate is the workload-weighted mean of the type rates. -/
theorem agg_eq_weighted (n d : ι → ℝ) (hd : ∀ k, 0 < d k) :
    agg n d = ∑ k, (d k / ∑ j, d j) * (n k / d k) := by
  unfold agg
  rw [sum_div]
  apply sum_congr rfl
  intro k _
  have hk : d k ≠ 0 := (hd k).ne'
  rw [div_mul_div_comm, mul_comm (d k) (n k), mul_div_mul_right _ _ hk]

/-- Kitagawa's identity: rate effect plus mix effect, for any centre `c`. -/
theorem kitagawa (wa wb ra rb : ι → ℝ) (c : ℝ) (ha : ∑ k, wa k = 1) (hb : ∑ k, wb k = 1) :
    ∑ k, wb k * rb k - ∑ k, wa k * ra k =
      ∑ k, (wa k + wb k) / 2 * (rb k - ra k) + ∑ k, ((ra k + rb k) / 2 - c) * (wb k - wa k) := by
  have hc : ∑ k, c * (wb k - wa k) = 0 := by
    rw [← mul_sum, sum_sub_distrib, ha, hb]; ring
  have hm : ∑ k, ((ra k + rb k) / 2 - c) * (wb k - wa k) =
      ∑ k, (ra k + rb k) / 2 * (wb k - wa k) - ∑ k, c * (wb k - wa k) := by
    rw [← sum_sub_distrib]; apply sum_congr rfl; intro k _; ring
  rw [hm, hc, sub_zero, ← sum_sub_distrib, ← sum_add_distrib]
  apply sum_congr rfl; intro k _; ring

/-- If no type's rate falls, the aggregate change is at least the mix effect. -/
theorem change_ge_mix (wa wb ra rb : ι → ℝ) (c : ℝ) (ha : ∑ k, wa k = 1) (hb : ∑ k, wb k = 1)
    (hwa : ∀ k, 0 ≤ wa k) (hwb : ∀ k, 0 ≤ wb k) (hr : ∀ k, ra k ≤ rb k) :
    ∑ k, ((ra k + rb k) / 2 - c) * (wb k - wa k) ≤ ∑ k, wb k * rb k - ∑ k, wa k * ra k := by
  rw [kitagawa wa wb ra rb c ha hb]
  have h : 0 ≤ ∑ k, (wa k + wb k) / 2 * (rb k - ra k) :=
    sum_nonneg fun k _ => mul_nonneg (by linarith [hwa k, hwb k]) (by linarith [hr k])
  linarith

/-- With the workload mix fixed, the aggregate cannot fall when no type's rate falls. -/
theorem fixed_mix_no_reversal (w ra rb : ι → ℝ) (hw : ∀ k, 0 ≤ w k) (hr : ∀ k, ra k ≤ rb k) :
    ∑ k, w k * ra k ≤ ∑ k, w k * rb k :=
  sum_le_sum fun k _ => mul_le_mul_of_nonneg_left (hr k) (hw k)

/-- Simpson's paradox: both types clear better in B (19/20 > 90/100, 60/140 > 20/50), yet
B's aggregate (79/160) is below A's (110/150). -/
theorem simpson_witness :
    (90 : ℝ) / 100 < 19 / 20 ∧ (20 : ℝ) / 50 < 60 / 140 ∧
      ((19 : ℝ) + 60) / (20 + 140) < (90 + 20) / (100 + 50) := by
  norm_num

end Csls.Clearance
