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

  Added:
  gtmp    Both match probabilities lie in (0, 1] (bn_hom_mem, bn_het_mem);
          the homozygote one is 1 at θ = 1 and tends to 1 as θ → 1.
  cllr    Positive for finite positive ratios (cllr_pos); 0 only as a limit
          (cllr_perfect_tendsto).
  dbsrch  Balding-Donnelly: uniform prior over N, n − 1 others excluded:
          1/(1 + (N − n)p) (balding_donnelly), at least the no-search value
          (post_ge_no_search); the exact chance of a chance match is
          1 − (1 − p)^n, in [0, 1], and the kernel's expm1/log1p form equals it.
  kinidx  The kernel's R1, R2 for every genotype-pair case and the IBD mixture
          (kin_mix); patidx: the trio cases 1/q, 1/(2q), 1/(p + q).
  mixlr   The kernel's inclusion-exclusion sum in closed form for one and two
          unknowns (not the general identity).
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

/-! ### Added: bounds, costs, searches, kinship, mixtures -/

open Finset

/-- A homozygote match probability lies in `(0, 1]`. -/
theorem bn_hom_mem (θ p : ℝ) (hθ : 0 ≤ θ) (hθ1 : θ < 1) (hp : 0 < p) (hp1 : p ≤ 1) :
    0 < bnHom θ p ∧ bnHom θ p ≤ 1 := by
  unfold bnHom
  have hd : 0 < (1 + θ) * (1 + 2 * θ) := by positivity
  have ha : 0 < (1 - θ) * p := mul_pos (by linarith) hp
  have ha1 : (1 - θ) * p ≤ 1 - θ := by nlinarith
  refine ⟨div_pos (by positivity) hd, (div_le_one hd).2 ?_⟩
  have h1 : 0 < 2 * θ + (1 - θ) * p := by linarith
  have h2 : 2 * θ + (1 - θ) * p ≤ 1 + θ := by linarith
  have h3 : 3 * θ + (1 - θ) * p ≤ 1 + 2 * θ := by linarith
  exact mul_le_mul h2 h3 (by linarith) (by linarith)

/-- A heterozygote match probability lies in `(0, 1]` when `p + q ≤ 1`. -/
theorem bn_het_mem (θ p q : ℝ) (hθ : 0 ≤ θ) (hθ1 : θ < 1) (hp : 0 < p) (hq : 0 < q)
    (hpq : p + q ≤ 1) : 0 < bnHet θ p q ∧ bnHet θ p q ≤ 1 := by
  unfold bnHet
  have hd : 0 < (1 + θ) * (1 + 2 * θ) := by positivity
  have hx : 0 < θ + (1 - θ) * p := by nlinarith
  have hy : 0 < θ + (1 - θ) * q := by nlinarith
  refine ⟨div_pos (by positivity) hd, (div_le_one hd).2 ?_⟩
  have hs : (θ + (1 - θ) * p) + (θ + (1 - θ) * q) ≤ 1 + θ := by nlinarith
  nlinarith [sq_nonneg ((θ + (1 - θ) * p) - (θ + (1 - θ) * q))]

/-- At full coancestry a homozygote always matches. -/
theorem bn_hom_one (p : ℝ) : bnHom 1 p = 1 := by unfold bnHom; norm_num

/-- ... and it tends to 1 as `θ → 1` from below. -/
theorem bn_hom_tendsto_one (p : ℝ) :
    Filter.Tendsto (fun θ => bnHom θ p) (nhdsWithin 1 (Set.Iio 1)) (nhds 1) := by
  have hc : ContinuousAt (fun θ => bnHom θ p) 1 := by
    unfold bnHom
    apply ContinuousAt.div (by fun_prop) (by fun_prop)
    norm_num
  have := hc.tendsto
  rw [bn_hom_one] at this
  exact tendsto_nhdsWithin_of_tendsto_nhds this

/-! cllr -/

/-- The cost of `m` same-source ratios `s` and `k` different-source ratios `d`, in bits. -/
def cllrOf {m k : ℕ} (s : Fin m → ℝ) (d : Fin k → ℝ) : ℝ :=
  ((∑ i, Real.logb 2 (1 + 1 / s i)) / m + (∑ j, Real.logb 2 (1 + d j)) / k) / 2

/-- Each term `log₂(1 + x)` is positive for `x > 0`. -/
theorem logb_one_add_pos (x : ℝ) (hx : 0 < x) : 0 < Real.logb 2 (1 + x) :=
  Real.logb_pos (by norm_num) (by linarith)

/-- Cllr is positive for any finite positive ratios: 0 is never reached. -/
theorem cllr_pos {m k : ℕ} (hm : 0 < m) (hk : 0 < k) (s : Fin m → ℝ) (d : Fin k → ℝ)
    (hs : ∀ i, 0 < s i) (hd : ∀ j, 0 < d j) : 0 < cllrOf s d := by
  unfold cllrOf
  have h1 : 0 < ∑ i, Real.logb 2 (1 + 1 / s i) := by
    have : Nonempty (Fin m) := ⟨⟨0, hm⟩⟩
    exact Finset.sum_pos (fun i _ => logb_one_add_pos _ (by have := hs i; positivity)) univ_nonempty
  have h2 : 0 < ∑ j, Real.logb 2 (1 + d j) := by
    have : Nonempty (Fin k) := ⟨⟨0, hk⟩⟩
    exact Finset.sum_pos (fun j _ => logb_one_add_pos _ (hd j)) univ_nonempty
  have : (0 : ℝ) < m := by exact_mod_cast hm
  have : (0 : ℝ) < k := by exact_mod_cast hk
  positivity

/-- ... but it is approached: same-source ratios `t` and different-source ratios `1/t` cost
`log₂(1 + 1/t)`, which tends to 0 as `t → ∞`. -/
theorem cllr_perfect_tendsto {m k : ℕ} (hm : 0 < m) (hk : 0 < k) :
    Filter.Tendsto (fun t : ℝ => cllrOf (fun _ : Fin m => t) (fun _ : Fin k => 1 / t))
      Filter.atTop (nhds 0) := by
  have e : ∀ t : ℝ, cllrOf (fun _ : Fin m => t) (fun _ : Fin k => 1 / t) = Real.logb 2 (1 + 1 / t) := by
    intro t
    unfold cllrOf
    simp only [sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul]
    have : (m : ℝ) ≠ 0 := by exact_mod_cast hm.ne'
    have : (k : ℝ) ≠ 0 := by exact_mod_cast hk.ne'
    field_simp
    ring
  simp_rw [e]
  have h1 : Filter.Tendsto (fun t : ℝ => 1 + 1 / t) Filter.atTop (nhds 1) := by
    simpa using (tendsto_const_nhds (x := (1:ℝ))).add tendsto_inv_atTop_zero
  have h2 := ((Real.continuousAt_logb (b := 2) (x := 1) (by norm_num)).tendsto).comp h1
  simpa [Real.logb_one, Function.comp_def] using h2
/-- Likelihood of the search result (index 0 matched, indices 1..n-1 excluded, the rest
untested) when source `i` left the stain. -/
def dbLik (n : ℕ) (p : ℝ) (i : ℕ) : ℝ :=
  if i = 0 then (1 - p) ^ (n - 1) else if i < n then 0 else p * (1 - p) ^ (n - 1)

/-- Posterior that the matching person is the source, from a uniform prior over `N`. -/
def dbPost (N n : ℕ) (p : ℝ) : ℝ :=
  (1 / N * dbLik n p 0) / ∑ i ∈ range N, 1 / N * dbLik n p i

/-- The evidence summed over the possible sources. -/
theorem db_sum (N n : ℕ) (p : ℝ) (hn : 1 ≤ n) (hnN : n ≤ N) :
    ∑ i ∈ range N, dbLik n p i = (1 - p) ^ (n - 1) * (1 + (N - n : ℕ) * p) := by
  rw [← sum_range_add_sum_Ico _ hnN]
  have h1 : ∑ i ∈ range n, dbLik n p i = (1 - p) ^ (n - 1) := by
    obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
    rw [sum_range_succ']
    have : ∀ i ∈ range m, dbLik (m + 1) p (i + 1) = 0 := by
      intro i hi
      simp only [mem_range] at hi
      simp [dbLik, show i + 1 < m + 1 by omega]
    rw [sum_eq_zero this]
    simp [dbLik]
  have h2 : ∑ i ∈ Ico n N, dbLik n p i = (N - n : ℕ) * (p * (1 - p) ^ (n - 1)) := by
    have : ∀ i ∈ Ico n N, dbLik n p i = p * (1 - p) ^ (n - 1) := by
      intro i hi
      simp only [mem_Ico] at hi
      simp [dbLik, show i ≠ 0 by omega, show ¬ i < n by omega]
    rw [sum_congr rfl this, sum_const, Nat.card_Ico, nsmul_eq_mul]
  rw [h1, h2]
  ring

/-- Balding-Donnelly: `P(source | search) = 1 / (1 + (N − n) p)`. -/
theorem balding_donnelly (N n : ℕ) (p : ℝ) (hn : 1 ≤ n) (hnN : n ≤ N) (hp0 : 0 ≤ p) (hp1 : p < 1) :
    dbPost N n p = 1 / (1 + (N - n : ℕ) * p) := by
  unfold dbPost
  rw [← mul_sum, db_sum N n p hn hnN]
  have hN : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
  have hq : 0 < (1 - p) ^ (n - 1) := pow_pos (by linarith) _
  have hc : 0 < 1 + ((N - n : ℕ) : ℝ) * p := by positivity
  have h0 : dbLik n p 0 = (1 - p) ^ (n - 1) := by simp [dbLik]
  rw [h0]
  field_simp

/-- The chance of at least one chance match among `n` independent profiles. -/
theorem any_match_eq (p : ℝ) (n : ℕ) :
    1 - ∏ _i ∈ range n, (1 - p) = 1 - (1 - p) ^ n := by
  rw [prod_const, card_range]

/-- It lies in `[0, 1]`. -/
theorem any_match_mem (p : ℝ) (n : ℕ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    0 ≤ 1 - (1 - p) ^ n ∧ 1 - (1 - p) ^ n ≤ 1 := by
  have h0 : 0 ≤ 1 - p := by linarith
  exact ⟨by linarith [pow_le_one₀ h0 (by linarith : 1 - p ≤ 1) (n := n)], by linarith [pow_nonneg h0 n]⟩

/-- The kernel's form `−expm1(n log1p(−p))` is the same number for `p < 1`. -/
theorem any_match_kernel (p : ℝ) (n : ℕ) (hp1 : p < 1) :
    -(Real.exp (n * Real.log (1 + -p)) - 1) = 1 - (1 - p) ^ n := by
  have h : 0 < 1 - p := by linarith
  rw [show 1 + -p = 1 - p by ring, ← Real.log_pow, Real.exp_log (pow_pos h n)]
  ring

/-! Kinship: the kernel's per-locus ratios, alleles as labels. -/
/-- The kernel's `R1` (`clab_kinidx`): one allele shared by descent. -/
def kinR1 (a b c d : ℕ) (pa pb : ℝ) : ℝ :=
  (if a = c ∨ a = d then 1 / 2 * (if c = d then 1 else 1 / 2) / pa else 0) +
  (if b = c ∨ b = d then 1 / 2 * (if c = d then 1 else 1 / 2) / pb else 0)
/-- The kernel's `R2`: both shared by descent. -/
def kinR2 (a b c d : ℕ) (pa pb : ℝ) : ℝ :=
  if (a = c ∧ b = d) ∨ (a = d ∧ b = c) then (if a = b then 1 / (pa * pa) else 1 / 2 / (pa * pb)) else 0
/-- The kernel's index `k0 + k1 R1 + k2 R2`. -/
def kin (k0 k1 k2 : ℝ) (a b c d : ℕ) (pa pb : ℝ) : ℝ := k0 + k1 * kinR1 a b c d pa pb + k2 * kinR2 a b c d pa pb

/-- Two identical homozygotes: `R1 = 1/p`. -/
theorem r1_hom_hom (a : ℕ) (p : ℝ) (hp : 0 < p) : kinR1 a a a a p p = 1 / p := by
  simp [kinR1]; field_simp; ring
/-- ... and `R2 = 1/p²`. -/
theorem r2_hom_hom (a : ℕ) (p : ℝ) : kinR2 a a a a p p = 1 / (p * p) := by simp [kinR2]
/-- Homozygote `aa` and heterozygote `ad`: `R1 = 1/(2p)`. -/
theorem r1_hom_het (a d : ℕ) (p : ℝ) (h : a ≠ d) (hp : 0 < p) :
    kinR1 a a a d p p = 1 / (2 * p) := by
  simp [kinR1, h]; field_simp; ring
/-- Heterozygote `ab` and homozygote `aa`: `R1 = 1/(2p)`. -/
theorem r1_het_hom (a b : ℕ) (p q : ℝ) (h : a ≠ b) (hp : 0 < p) :
    kinR1 a b a a p q = 1 / (2 * p) := by
  simp [kinR1, Ne.symm h]; field_simp
/-- Heterozygotes `ab` and `ad` sharing one allele: `R1 = 1/(4p)`. -/
theorem r1_het_share_one (a b d : ℕ) (p q : ℝ) (hab : a ≠ b) (had : a ≠ d) (hbd : b ≠ d) (hp : 0 < p) :
    kinR1 a b a d p q = 1 / (4 * p) := by
  simp [kinR1, had, Ne.symm hab, hbd]; field_simp; ring
/-- Two identical heterozygotes: `R1 = (p + q)/(4pq)`. -/
theorem r1_het_identical (a b : ℕ) (p q : ℝ) (hab : a ≠ b) (hp : 0 < p) (hq : 0 < q) :
    kinR1 a b a b p q = (p + q) / (4 * p * q) := by
  simp [kinR1, hab]; field_simp; ring
/-- ... and `R2 = 1/(2pq)`. -/
theorem r2_het_identical (a b : ℕ) (p q : ℝ) (hab : a ≠ b) : kinR2 a b a b p q = 1 / (2 * p * q) := by
  simp [kinR2, hab]; field_simp
/-- Different genotypes cannot share both alleles by descent: `R2 = 0`. -/
theorem r2_not_identical (a b c d : ℕ) (p q : ℝ) (h : ¬ ((a = c ∧ b = d) ∨ (a = d ∧ b = c))) :
    kinR2 a b c d p q = 0 := by simp [kinR2, h]

/-- The index is the IBD mixture: with `P(g₂ | g₁, IBD = i) = Pᵢ` and `P₀ = P(g₂)`, the ratio of
`∑ kᵢ Pᵢ` to `P(g₂)` is `k₀ + k₁ (P₁/P) + k₂ (P₂/P)`. -/
theorem kin_mix (k0 k1 k2 P P1 P2 : ℝ) (hP : P ≠ 0) :
    (k0 * P + k1 * P1 + k2 * P2) / P = k0 + k1 * (P1 / P) + k2 * (P2 / P) := by
  field_simp

/-- The parent-child case of the derivation: g₁ = aa, g₂ = ab (one allele by descent, `b` from
the population) gives `P₁ / P = p_b / (2 p_a p_b) = 1 / (2 p_a)`, the kernel's `R1`. -/
theorem r1_hom_het_derived (a d : ℕ) (p q : ℝ) (h : a ≠ d) (hp : 0 < p) (hq : 0 < q) :
    q / (2 * p * q) = kinR1 a a a d p p := by
  rw [r1_hom_het a d p h hp]; field_simp

/-! Paternity with the mother: the kernel's index. Child c₁c₂ (frequencies p₁ p₂), mother m₁m₂,
alleged father f₁f₂. -/
/-- One mother/father allele pair in the kernel's numerator: 1/4 if it makes the child. -/
def pt (g h c1 c2 : ℕ) : ℝ := if (g = c1 ∧ h = c2) ∨ (g = c2 ∧ h = c1) then 1 / 4 else 0
/-- One maternal allele in the kernel's denominator: the chance a random man completes the child. -/
def pd (g c1 c2 : ℕ) (p1 p2 : ℝ) : ℝ :=
  if c1 = c2 then (if g = c1 then 1 / 2 * p1 else 0)
  else (if g = c1 then 1 / 2 * p2 else 0) + (if g = c2 then 1 / 2 * p1 else 0)
/-- The kernel's paternity index with the mother (`clab_patidx`). -/
def pat (c1 c2 m1 m2 f1 f2 : ℕ) (p1 p2 : ℝ) : ℝ :=
  (pt m1 f1 c1 c2 + pt m1 f2 c1 c2 + pt m2 f1 c1 c2 + pt m2 f2 c1 c2) /
    (pd m1 c1 c2 p1 p2 + pd m2 c1 c2 p1 p2)

/-- Child `qq`, mother `qx`, father `qq`: `1/q`. -/
theorem pat_hom_father_hom (q' x : ℕ) (q : ℝ) (hx : x ≠ q') (hq : 0 < q) :
    pat q' q' q' x q' q' q q = 1 / q := by
  simp [pat, pt, pd, hx]; field_simp; ring
/-- Child `qq`, mother `qx`, father `qy` (y ≠ q): `1/(2q)`. -/
theorem pat_hom_father_het (q' x y : ℕ) (q : ℝ) (hx : x ≠ q') (hy : y ≠ q') (hq : 0 < q) :
    pat q' q' q' x q' y q q = 1 / (2 * q) := by
  simp [pat, pt, pd, hx, hy]; field_simp; ring
/-- Child `ab`, mother `ab`, father `ab`: `1/(p + q)`. -/
theorem pat_het_trio (a b : ℕ) (p q : ℝ) (hab : a ≠ b) (hp : 0 < p) (hq : 0 < q) :
    pat a b a b a b p q = 1 / (p + q) := by
  simp [pat, pt, pd, hab, Ne.symm hab]; field_simp; ring

/-! Mixtures: the kernel's inclusion-exclusion sum, closed forms. -/
/-- One unknown, one allele to carry, all alleles within total `T`: `T² − (T − p)² = p(2T − p)`. -/
theorem mix_one_unknown_one_allele (T p : ℝ) : T ^ 2 - (T - p) ^ 2 = p * (2 * T - p) := by ring
/-- One unknown carrying both of two alleles: `2pq`. -/
theorem mix_one_unknown_two (T p q : ℝ) :
    T ^ 2 - (T - p) ^ 2 - (T - q) ^ 2 + (T - p - q) ^ 2 = 2 * p * q := by ring
/-- One unknown cannot carry three alleles: the kernel's sum is 0 (its shortcut for `m > 2x`). -/
theorem mix_one_unknown_three_zero (T p q r : ℝ) :
    T ^ 2 - (T - p) ^ 2 - (T - q) ^ 2 - (T - r) ^ 2 + (T - p - q) ^ 2 + (T - p - r) ^ 2
      + (T - q - r) ^ 2 - (T - p - q - r) ^ 2 = 0 := by ring
/-- Two unknowns carrying four distinct alleles: `24 p₁p₂p₃p₄`. -/
theorem mix_two_unknowns_four (T a b c d : ℝ) :
    T ^ 4 - (T - a) ^ 4 - (T - b) ^ 4 - (T - c) ^ 4 - (T - d) ^ 4
      + (T - a - b) ^ 4 + (T - a - c) ^ 4 + (T - a - d) ^ 4 + (T - b - c) ^ 4 + (T - b - d) ^ 4 + (T - c - d) ^ 4
      - (T - a - b - c) ^ 4 - (T - a - b - d) ^ 4 - (T - a - c - d) ^ 4 - (T - b - c - d) ^ 4
      + (T - a - b - c - d) ^ 4 = 24 * a * b * c * d := by ring

/-- The search raises the matching person's posterior above the no-search value
`1/(1 + (N − 1)p)`. -/
theorem post_ge_no_search (N n : ℕ) (p : ℝ) (hn : 1 ≤ n) (hnN : n ≤ N) (hp0 : 0 ≤ p) :
    1 / (1 + (N - 1 : ℕ) * p) ≤ 1 / (1 + (N - n : ℕ) * p) := by
  apply one_div_le_one_div_of_le (by positivity)
  have : ((N - n : ℕ) : ℝ) ≤ (N - 1 : ℕ) := by exact_mod_cast (by omega)
  nlinarith

end Crimlab.Forensics
