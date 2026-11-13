/-
  crimlab: inst/lean/Learn.lean

  The identities behind the learners and the cross-fitted estimators.

  * soft_threshold_prox: soft-thresholding is the proximal map of t |.|:
    for t >= 0, soft z t minimises (b - z)^2 / 2 + t |b|.
  * cd_update_min: the coordinate-descent update of the elastic net,
    soft u l1 / (v + l2), minimises the one-dimensional objective
    (v + l2) / 2 b^2 - u b + l1 |b| when v + l2 > 0 and l1 >= 0.
  * dml_root: theta = -sum psi_b / sum psi_a solves the empirical linear
    score sum (psi_a theta + psi_b) = 0 when sum psi_a is not 0.
  * plr_orthogonal: the partialling-out score of the partially linear model
    is Neyman-orthogonal: on a finite weighted population, when the
    treatment residual D - m0 and the structural error Y - l0 - theta (D -
    m0) are orthogonal to every function of the covariates, the mean score
    along a perturbation (l0 + r h, m0 + r g) has no term linear in r.
  * aipw_double_robust: the augmented inverse probability term has
    conditional bias (mu - g) (e - m) / m, so it is unbiased when either
    the outcome regression (g = mu) or the propensity (m = e) is right.

  Not proved: rates of convergence, or anything about the learners'
  statistical accuracy.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section
open Finset

namespace Crimlab.Learn

/-- Soft-thresholding. -/
def soft (z t : ℝ) : ℝ := if t < z then z - t else if z < -t then z + t else 0

/-- Soft-thresholding is the proximal map of `t |.|`. -/
theorem soft_threshold_prox (z t b : ℝ) (ht : 0 ≤ t) :
    (soft z t - z) ^ 2 / 2 + t * |soft z t| ≤ (b - z) ^ 2 / 2 + t * |b| := by
  unfold soft
  split_ifs with h1 h2
  · rw [abs_of_pos (by linarith)]
    rcases le_or_gt 0 b with hb | hb
    · rw [abs_of_nonneg hb]; nlinarith [sq_nonneg (b - (z - t))]
    · rw [abs_of_neg hb]; nlinarith [sq_nonneg (b - (z - t))]
  · rw [abs_of_neg (by linarith)]
    rcases le_or_gt 0 b with hb | hb
    · rw [abs_of_nonneg hb]; nlinarith [sq_nonneg (b - (z + t))]
    · rw [abs_of_neg hb]; nlinarith [sq_nonneg (b - (z + t))]
  · simp only [abs_zero, mul_zero, add_zero]
    rcases le_or_gt 0 b with hb | hb
    · rw [abs_of_nonneg hb]; nlinarith [sq_nonneg b]
    · rw [abs_of_neg hb]; nlinarith [sq_nonneg b]

/-- The one-dimensional elastic net objective of a coordinate. -/
def cdobj (v u l1 l2 b : ℝ) : ℝ := (v + l2) / 2 * b ^ 2 - u * b + l1 * |b|

/-- The coordinate update `soft u l1 / (v + l2)` minimises it. -/
theorem cd_update_min (v u l1 l2 b : ℝ) (hv : 0 < v + l2) (hl : 0 ≤ l1) :
    cdobj v u l1 l2 (soft u l1 / (v + l2)) ≤ cdobj v u l1 l2 b := by
  -- with c = v + l2, the objective is c times the proximal objective at
  -- u / c with threshold l1 / c, minus a constant
  set c := v + l2 with hc
  have key : ∀ x, cdobj v u l1 l2 x =
      c * ((x - u / c) ^ 2 / 2 + l1 / c * |x|) - u ^ 2 / (2 * c) := by
    intro x
    unfold cdobj
    field_simp
    ring
  have hs : soft u l1 / c = soft (u / c) (l1 / c) := by
    unfold soft
    have hc0 : (0 : ℝ) < c := hv
    by_cases h1 : l1 < u
    · have : l1 / c < u / c := div_lt_div_of_pos_right h1 hc0
      rw [if_pos h1, if_pos this]; field_simp
    · by_cases h2 : u < -l1
      · have h3 : ¬ (l1 / c < u / c) := by
          intro h; exact h1 ((div_lt_div_iff_of_pos_right hc0).mp h)
        have h4 : u / c < -(l1 / c) := by
          rw [← neg_div]; exact div_lt_div_of_pos_right h2 hc0
        rw [if_neg h1, if_pos h2, if_neg h3, if_pos h4]; field_simp
      · have h3 : ¬ (l1 / c < u / c) := by
          intro h; exact h1 ((div_lt_div_iff_of_pos_right hc0).mp h)
        have h4 : ¬ (u / c < -(l1 / c)) := by
          intro h; rw [← neg_div] at h
          exact h2 ((div_lt_div_iff_of_pos_right hc0).mp h)
        rw [if_neg h1, if_neg h2, if_neg h3, if_neg h4]; simp
  rw [key, key, hs]
  have := soft_threshold_prox (u / c) (l1 / c) b (div_nonneg hl hv.le)
  nlinarith [this]

/-- The root of the empirical linear score. -/
theorem dml_root {ι : Type*} (s : Finset ι) (pa pb : ι → ℝ)
    (h : ∑ i ∈ s, pa i ≠ 0) :
    ∑ i ∈ s, (pa i * (-(∑ j ∈ s, pb j) / ∑ j ∈ s, pa j) + pb i) = 0 := by
  rw [sum_add_distrib, ← sum_mul]
  field_simp
  ring

/-- Neyman orthogonality of the partialling-out score: on a finite
population with weights `w`, if the treatment residual and the structural
error are orthogonal to every function of the covariates (`x` maps a unit
to its covariates), the mean score along `(l0 + r h, m0 + r g)` is its
value at `r = 0` plus a term in `r ^ 2` only. -/
theorem plr_orthogonal {Ω X : Type*} [Fintype Ω] (w : Ω → ℝ) (x : Ω → X)
    (y d : Ω → ℝ) (l0 m0 h g : X → ℝ) (θ : ℝ)
    (hd : ∀ f : X → ℝ, ∑ i, w i * f (x i) * (d i - m0 (x i)) = 0)
    (he : ∀ f : X → ℝ,
      ∑ i, w i * f (x i) * (y i - l0 (x i) - θ * (d i - m0 (x i))) = 0) :
    ∃ c : ℝ, ∀ r : ℝ,
      ∑ i, w i * ((y i - (l0 (x i) + r * h (x i)) -
          θ * (d i - (m0 (x i) + r * g (x i)))) * (d i - (m0 (x i) + r * g (x i)))) =
        ∑ i, w i * ((y i - l0 (x i) - θ * (d i - m0 (x i))) * (d i - m0 (x i))) +
          c * r ^ 2 := by
  refine ⟨∑ i, w i * ((h (x i) - θ * g (x i)) * g (x i)), fun r => ?_⟩
  have e1 := hd (fun z => h z - θ * g z)
  have e2 := he g
  have hexp : ∀ i, w i * ((y i - (l0 (x i) + r * h (x i)) -
        θ * (d i - (m0 (x i) + r * g (x i)))) * (d i - (m0 (x i) + r * g (x i)))) =
      w i * ((y i - l0 (x i) - θ * (d i - m0 (x i))) * (d i - m0 (x i))) -
        r * (w i * (h (x i) - θ * g (x i)) * (d i - m0 (x i))) -
        r * (w i * g (x i) * (y i - l0 (x i) - θ * (d i - m0 (x i)))) +
        r ^ 2 * (w i * ((h (x i) - θ * g (x i)) * g (x i))) := by
    intro i; ring
  simp only [hexp, sum_add_distrib, sum_sub_distrib, ← mul_sum, e1, e2]
  ring

/-- The conditional bias of the augmented inverse probability term for the
treated, given the covariates: with true outcome mean `μ` and propensity
`e`, the term `g + D (Y - g) / m` has conditional mean
`g + e (μ - g) / m`, which exceeds `μ` by `(μ - g) (e - m) / m`. -/
theorem aipw_bias (μ g e m : ℝ) (hm : m ≠ 0) :
    g + e * (μ - g) / m - μ = (μ - g) * (e - m) / m := by
  field_simp
  ring

/-- Double robustness: the bias vanishes when the outcome regression is
right or when the propensity is. -/
theorem aipw_double_robust (μ g e m : ℝ) (hm : m ≠ 0) (h : g = μ ∨ m = e) :
    g + e * (μ - g) / m = μ := by
  have := aipw_bias μ g e m hm
  rcases h with h | h
  · subst h; simp
  · subst h; rw [sub_self, mul_zero, zero_div] at this; linarith

end Crimlab.Learn
