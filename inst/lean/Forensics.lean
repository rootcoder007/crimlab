/-
  crimlab: the identities and bounds the forensic functions rest on.

  Compiles against Mathlib (Lean 4). Each theorem is named by proof() in the
  R package next to the function that relies on it.

  gtmp    Balding-Nichols match probabilities: the product rule at theta = 0;
          a homozygote match is never rarer with coancestry than without; the
          same is false for a heterozygote (a counterexample is given).
  lrpost  Bayes' rule in odds form, and the prosecutor's fallacy: a ratio
          below the prior odds against leaves the posterior below one half.
  cllr    The cost of a system that always reports a ratio of 1 is one bit.
  dbsrch  The np rule is an upper bound on the chance of a chance match
          somewhere in a database of n profiles.
  mixlr   Inclusion-exclusion for one unknown contributor to two alleles.
  kinidx  Parent and child: the index of two identical heterozygotes, and
          the paternity probability at even prior odds.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section

namespace Crimlab.Forensics

/-- Balding-Nichols match probability of a homozygote. -/
def bnHom (θ p : ℝ) : ℝ :=
  (2 * θ + (1 - θ) * p) * (3 * θ + (1 - θ) * p) / ((1 + θ) * (1 + 2 * θ))

/-- Balding-Nichols match probability of a heterozygote. -/
def bnHet (θ p q : ℝ) : ℝ :=
  2 * (θ + (1 - θ) * p) * (θ + (1 - θ) * q) / ((1 + θ) * (1 + 2 * θ))

theorem bn_hom_zero (p : ℝ) : bnHom 0 p = p ^ 2 := by
  unfold bnHom; ring

theorem bn_het_zero (p q : ℝ) : bnHet 0 p q = 2 * p * q := by
  unfold bnHet; ring

/-- Coancestry never makes a homozygote match rarer than the product rule. -/
theorem bn_hom_ge_product (θ p : ℝ) (hθ : 0 ≤ θ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    p ^ 2 ≤ bnHom θ p := by
  unfold bnHom
  have hd : 0 < (1 + θ) * (1 + 2 * θ) := by positivity
  rw [le_div_iff₀ hd]
  have h1 : 0 ≤ θ * (5 * p * (1 - p)) := by
    have : 0 ≤ 5 * p * (1 - p) := by nlinarith
    exact mul_nonneg hθ this
  have h2 : 0 ≤ θ * (θ * ((1 - p) * (6 + p))) := by
    have : 0 ≤ (1 - p) * (6 + p) := by nlinarith
    exact mul_nonneg hθ (mul_nonneg hθ this)
  nlinarith [h1, h2]

/-- The same is false for a heterozygote: with two common alleles coancestry
    makes the match rarer. -/
theorem bn_het_not_ge_product :
    bnHet (1 / 10) (1 / 2) (1 / 2) < 2 * (1 / 2 : ℝ) * (1 / 2) := by
  unfold bnHet; norm_num

/-- Posterior probability from prior probability π and likelihood ratio L. -/
def post (π L : ℝ) : ℝ := 1 / (1 + (1 - π) / (π * L))

theorem post_eq (π L : ℝ) (hπ : 0 < π) (hL : 0 < L) (hπ1 : π < 1) :
    post π L = π * L / (π * L + (1 - π)) := by
  unfold post
  have h : 0 < π * L := mul_pos hπ hL
  have h2 : 0 < π * L + (1 - π) := by nlinarith
  field_simp

/-- Prosecutor's fallacy: however large L is, if it is below the prior odds
    against, the hypothesis is still more likely false than true. -/
theorem post_lt_half (π L : ℝ) (hπ : 0 < π) (hL : 0 < L) (hπ1 : π < 1)
    (h : π * L < 1 - π) : post π L < 1 / 2 := by
  rw [post_eq π L hπ hL hπ1]
  have h2 : 0 < π * L + (1 - π) := by nlinarith [mul_pos hπ hL]
  rw [div_lt_iff₀ h2]
  nlinarith

/-- Even prior odds: the posterior is L / (L + 1), the probability of
    paternity reported from a paternity index. -/
theorem post_even_prior (L : ℝ) (hL : 0 < L) : post (1 / 2) L = L / (L + 1) := by
  rw [post_eq (1 / 2) L (by norm_num) hL (by norm_num)]
  have : 0 < L + 1 := by linarith
  field_simp
  ring

/-- A system that always reports a ratio of 1 costs one bit. -/
theorem cllr_neutral :
    (Real.logb 2 (1 + 1 / 1) + Real.logb 2 (1 + 1)) / 2 = 1 := by
  have h : Real.logb 2 2 = 1 := Real.logb_self_eq_one (by norm_num)
  norm_num [h]

/-- The np rule bounds the chance of at least one chance match among n
    independent profiles, each matching with probability p. -/
theorem np_rule (p : ℝ) (n : ℕ) (hp1 : p ≤ 1) :
    1 - (1 - p) ^ n ≤ n * p := by
  have h := one_add_mul_le_pow (a := -p) (by linarith) n
  have e : (1 + -p) = 1 - p := by ring
  rw [e] at h
  linarith

/-- One unknown contributor who must carry both of two alleles: the
    inclusion-exclusion sum is the heterozygote probability. -/
theorem mix_two_alleles (p q : ℝ) :
    (p + q) ^ 2 - p ^ 2 - q ^ 2 = 2 * p * q := by ring

/-- Parent and child, both heterozygous for the same two alleles: the chance
    of the pair over the chance of two unrelated people is (p+q)/(4pq). -/
theorem po_identical_het (p q : ℝ) (hp : 0 < p) (hq : 0 < q) :
    (2 * p * q * ((p + q) / 2)) / ((2 * p * q) * (2 * p * q)) = (p + q) / (4 * p * q) := by
  field_simp
  ring

end Crimlab.Forensics
