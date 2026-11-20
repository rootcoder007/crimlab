/-
  crimlab: inst/lean/Volat.lean

  Placement volatility and alert complexity (plcvol(), alrtcx()), the two
  person-year measures defined in Appendix H of Ruhela's study of Ontario's
  restrictive confinement records:

    * placement volatility of a person-year whose records r = 1, ..., R
      (taken in a given order) carry the region at placement A_r and the
      region of the most recent placement B_r:
        v = sum_r 1(A_r /= B_r) + sum_{r > 1} 1(A_r /= A_{r-1});
      `within` is the first sum, `chg` the second, `vol` both;
    * alert complexity: the number of distinct alert configurations among
      the person-year's records (`acx`), and the treatment D = 1(acx >= 2).

  What is proved: the bounds and their sharpness; the values with no
  change and with one record; invariance of both measures to an injective
  relabelling of the categories, and that a merging relabelling can lower
  them; invariance of `within` and `acx` to the order of the records, and
  that `chg` is not order invariant (a counterexample), while on records
  sorted by region it is one less than the number of distinct regions and
  so depends only on the set of regions; the effect of adding a record and
  of joining two runs of records (a pooling of years for one person would
  add at most one change); the pooled mean as the weighted mean of the
  yearly means, and a counterexample showing that a pooled crude ratio of
  means can lie above every year's ratio; the difference between the
  placement-weighted and the unweighted mean as a covariance; for `acx`,
  the bounds, the characterisation of D = 1 as two records in different
  configurations, and two counterexamples showing that D = 1 does not
  require two alerts on one record and that a record with all three alerts
  gives D = 0.

  These are facts about the definitions. Nothing here says that a measure
  is a good proxy for instability or for vulnerability; that is a
  substantive question about the records.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

namespace Crimlab.Volat

open Finset

section Vol

variable {α : Type*} [DecidableEq α]

/-- Records whose two region fields differ. -/
def within (l : List (α × α)) : ℕ := l.countP (fun r => decide (r.1 ≠ r.2))

/-- Changes between consecutive entries of a list. -/
def chg : List α → ℕ
  | a :: b :: t => (if a = b then 0 else 1) + chg (b :: t)
  | _ => 0

/-- Placement volatility of a person-year's records, in the order given. -/
def vol (l : List (α × α)) : ℕ := within l + chg (l.map Prod.fst)

theorem chg_nil : chg ([] : List α) = 0 := rfl

theorem chg_single (a : α) : chg [a] = 0 := rfl

theorem chg_cons_cons (a b : α) (t : List α) :
    chg (a :: b :: t) = (if a = b then 0 else 1) + chg (b :: t) := rfl

/-- No placement: no volatility. -/
theorem vol_nil : vol ([] : List (α × α)) = 0 := rfl

/-- One record: volatility 1 when its two regions differ, else 0. -/
theorem vol_single (a b : α) : vol [(a, b)] = if a = b then 0 else 1 := by
  by_cases h : a = b <;> simp [vol, within, chg, h]

/-- A list whose entries are all equal has no change. -/
theorem chg_const (c : α) : ∀ m : List α, (∀ x ∈ m, x = c) → chg m = 0
  | [], _ => rfl
  | [_], _ => rfl
  | a :: b :: t, h => by
      have hab : a = b := (h a (by simp)).trans (h b (by simp)).symm
      rw [chg_cons_cons, if_pos hab, zero_add]
      exact chg_const c (b :: t) (fun x hx => h x (List.mem_cons_of_mem a hx))

/-- No change of placement: every record in one region `c`, recorded as
such in both fields, gives volatility 0. -/
theorem vol_stable (c : α) (l : List (α × α)) (h : ∀ r ∈ l, r.1 = c ∧ r.2 = c) :
    vol l = 0 := by
  have h1 : within l = 0 := by
    rw [within, List.countP_eq_zero]
    intro r hr
    simp [(h r hr).1, (h r hr).2]
  have h2 : chg (l.map Prod.fst) = 0 :=
    chg_const c _ (by
      intro x hx
      obtain ⟨r, hr, rfl⟩ := List.mem_map.mp hx
      exact (h r hr).1)
  simp [vol, h1, h2]

/-- No change between consecutive entries: every entry is the first. -/
theorem chg_eq_zero : ∀ (a : α) (t : List α), chg (a :: t) = 0 → ∀ x ∈ a :: t, x = a
  | a, [], _ => by simp
  | a, b :: t, h => by
      rw [chg_cons_cons] at h
      have hab : a = b := by
        by_contra hne
        rw [if_neg hne] at h
        omega
      rw [if_pos hab, zero_add] at h
      intro x hx
      rcases List.mem_cons.mp hx with rfl | hx
      · rfl
      · rw [hab]; exact chg_eq_zero b t h x hx

/-- Volatility is 0 exactly when every record names one region in both
fields. So the regional cluster (the regions a person-year's records
name) has two regions or more exactly when volatility is positive: a
model that conditions on the cluster conditions on whether the count is
zero. -/
theorem vol_eq_zero_iff (r : α × α) (l : List (α × α)) :
    vol (r :: l) = 0 ↔ ∀ s ∈ r :: l, s.1 = r.1 ∧ s.2 = r.1 := by
  constructor
  · intro h
    have hw : within (r :: l) = 0 := by unfold vol at h; omega
    have hc : chg ((r :: l).map Prod.fst) = 0 := by unfold vol at h; omega
    rw [within, List.countP_eq_zero] at hw
    intro s hs
    have h1 : s.1 = r.1 :=
      chg_eq_zero r.1 (l.map Prod.fst) (by simpa using hc) s.1
        (by simpa using List.mem_map_of_mem (f := Prod.fst) hs)
    have h2 : s.1 = s.2 := by simpa using hw s hs
    exact ⟨h1, h2 ▸ h1⟩
  · intro h
    exact vol_stable r.1 _ h

theorem within_le (l : List (α × α)) : within l ≤ l.length :=
  List.countP_le_length

theorem chg_le : ∀ m : List α, chg m ≤ m.length - 1
  | [] => Nat.zero_le _
  | [_] => Nat.zero_le _
  | a :: b :: t => by
      have ih := chg_le (b :: t)
      rw [chg_cons_cons]
      simp only [List.length_cons] at ih ⊢
      split_ifs <;> omega

/-- The bound: at most `2 R - 1` for `R` records. -/
theorem vol_le (l : List (α × α)) : vol l ≤ 2 * l.length - 1 := by
  have h1 := within_le l
  have h2 := chg_le (l.map Prod.fst)
  rw [List.length_map] at h2
  unfold vol
  rcases Nat.eq_zero_or_pos l.length with h | h <;> omega

/-- The bound is attained: three records alternating between two regions
whose two fields always differ. -/
theorem vol_le_sharp : vol [((0 : ℕ), (1 : ℕ)), (1, 0), (0, 1)] = 2 * 3 - 1 := by
  decide

/-- Relabelling the regions by an injective map leaves the changes alone. -/
theorem chg_map {β : Type*} [DecidableEq β] (f : α → β) (hf : Function.Injective f) :
    ∀ m : List α, chg (m.map f) = chg m
  | [] => rfl
  | [_] => rfl
  | a :: b :: t => by
      have ih := chg_map f hf (b :: t)
      simp only [List.map_cons] at ih ⊢
      rw [chg_cons_cons, chg_cons_cons, ih]
      simp only [hf.eq_iff]

/-- Volatility is invariant to an injective relabelling of the regions. -/
theorem vol_map {β : Type*} [DecidableEq β] (f : α → β) (hf : Function.Injective f)
    (l : List (α × α)) : vol (l.map (Prod.map f f)) = vol l := by
  have h1 : within (l.map (Prod.map f f)) = within l := by
    unfold within
    rw [List.countP_map]
    apply List.countP_congr
    intro r _
    simp [hf.eq_iff]
  have h2 : (l.map (Prod.map f f)).map Prod.fst = (l.map Prod.fst).map f := by
    simp [List.map_map, Function.comp_def]
  rw [vol, vol, h1, h2, chg_map f hf]

/-- Counterexample: a relabelling that merges two regions can lower
volatility (here from 1 to 0). -/
theorem vol_merge_lowers :
    vol [((0 : ℕ), (1 : ℕ))] = 1 ∧
      vol ([((0 : ℕ), (1 : ℕ))].map (Prod.map (fun _ => (0 : ℕ)) (fun _ => (0 : ℕ)))) = 0 := by
  decide

/-- The records whose fields differ do not depend on the order of the
records. -/
theorem within_perm {l₁ l₂ : List (α × α)} (h : l₁.Perm l₂) : within l₁ = within l₂ := by
  unfold within; exact h.countP_eq _

/-- Counterexample: the consecutive changes, and so volatility, depend on
the order of the records: the same three records give 2 or 1. -/
theorem vol_order_matters :
    ([((0 : ℕ), (0 : ℕ)), (1, 1), (0, 0)]).Perm [(0, 0), (0, 0), (1, 1)] ∧
      vol [((0 : ℕ), (0 : ℕ)), (1, 1), (0, 0)] = 2 ∧
      vol [((0 : ℕ), (0 : ℕ)), (0, 0), (1, 1)] = 1 := by
  decide

/-- Whatever the order, the changes are at least the distinct entries
less one. -/
theorem dedup_le_chg : ∀ m : List α, m.dedup.length ≤ chg m + 1
  | [] => by simp
  | [a] => by simp
  | a :: b :: t => by
      have ih := dedup_le_chg (b :: t)
      rw [chg_cons_cons]
      by_cases hm : a ∈ b :: t
      · rw [List.dedup_cons_of_mem hm]
        split_ifs <;> omega
      · have hab : a ≠ b := fun h => hm (by simp [h])
        rw [List.dedup_cons_of_notMem hm, List.length_cons, if_neg hab]
        omega

/-- On records sorted by region the changes are exactly the distinct
regions less one (the order the published code uses). -/
theorem chg_sorted [LinearOrder α] :
    ∀ m : List α, m ≠ [] → m.Pairwise (· ≤ ·) → chg m + 1 = m.dedup.length
  | [], h, _ => absurd rfl h
  | [a], _, _ => by simp [chg]
  | a :: b :: t, _, hs => by
      have hs' : (b :: t).Pairwise (· ≤ ·) := List.Pairwise.of_cons hs
      have ih := chg_sorted (b :: t) (List.cons_ne_nil b t) hs'
      have hab : a ≤ b := List.rel_of_pairwise_cons hs (by simp)
      rw [chg_cons_cons]
      by_cases h : a = b
      · rw [if_pos h, List.dedup_cons_of_mem (by simp [h])]
        omega
      · have hlt : a < b := lt_of_le_of_ne hab h
        have hm : a ∉ b :: t := by
          intro hm
          rcases List.mem_cons.mp hm with h' | h'
          · exact h h'
          · have := List.rel_of_pairwise_cons hs' h'
            exact absurd (lt_of_lt_of_le hlt this) (lt_irrefl a)
        rw [List.dedup_cons_of_notMem hm, List.length_cons, if_neg h]
        omega

/-- So on sorted records the changes do not depend on the order in which
the records arrived: two sorted rearrangements of the same records agree. -/
theorem chg_sorted_perm [LinearOrder α] {m₁ m₂ : List α} (hp : m₁.Perm m₂)
    (h₁ : m₁.Pairwise (· ≤ ·)) (h₂ : m₂.Pairwise (· ≤ ·)) : chg m₁ = chg m₂ := by
  rcases m₁ with _ | ⟨a, t⟩
  · rw [List.Perm.nil_eq hp]
  · have hne₂ : m₂ ≠ [] := by
      intro h; rw [h] at hp; exact List.cons_ne_nil a t (List.Perm.eq_nil hp)
    have e₁ := chg_sorted _ (List.cons_ne_nil a t) h₁
    have e₂ := chg_sorted _ hne₂ h₂
    have : (a :: t).dedup.length = m₂.dedup.length :=
      (List.Perm.dedup hp).length_eq
    omega

/-- Adding a record in front never lowers the changes. -/
theorem chg_le_cons (a : α) : ∀ m : List α, chg m ≤ chg (a :: m)
  | [] => Nat.zero_le _
  | b :: t => by rw [chg_cons_cons]; omega

/-- Adding a record never lowers volatility. -/
theorem vol_le_cons (r : α × α) (l : List (α × α)) : vol l ≤ vol (r :: l) := by
  have h1 : within l ≤ within (r :: l) := by
    unfold within; rw [List.countP_cons]; split <;> omega
  have h2 := chg_le_cons r.1 (l.map Prod.fst)
  simp only [vol, List.map_cons]
  omega

/-- Joining two runs of records (one person's records of two years, were
the identifiers kept across years): the changes add, plus at most one at
the join. -/
theorem chg_append : ∀ m₁ m₂ : List α,
    chg m₁ + chg m₂ ≤ chg (m₁ ++ m₂) ∧ chg (m₁ ++ m₂) ≤ chg m₁ + chg m₂ + 1
  | [], m₂ => by simp [chg]
  | [a], [] => by simp [chg]
  | [a], b :: t => by
      simp only [List.singleton_append, chg_single, zero_add, chg_cons_cons]
      constructor <;> split_ifs <;> omega
  | a :: b :: t, m₂ => by
      obtain ⟨h1, h2⟩ := chg_append (b :: t) m₂
      simp only [List.cons_append] at h1 h2 ⊢
      rw [chg_cons_cons, chg_cons_cons]
      constructor <;> split_ifs <;> omega

end Vol

section Means

/-- Pooling two years: the pooled mean is the mean of the yearly means
weighted by the yearly counts. -/
theorem pool_mean (s₁ s₂ n₁ n₂ : ℝ) (h₁ : 0 < n₁) (h₂ : 0 < n₂) :
    (s₁ + s₂) / (n₁ + n₂) = n₁ / (n₁ + n₂) * (s₁ / n₁) + n₂ / (n₁ + n₂) * (s₂ / n₂) := by
  have h₃ : n₁ + n₂ ≠ 0 := by positivity
  have := h₁.ne'
  have := h₂.ne'
  field_simp

/-- Counterexample: pooled over two years in which the ratio of the mean
counts of two groups is 2 in each year, the pooled ratio is 38/11, above
both, because the groups' shares differ between the years. A crude ratio
on pooled years is not an average of the yearly ratios. -/
theorem pool_ratio_outside :
    (2 : ℝ) / 1 = 2 ∧ (4 : ℝ) / 2 = 2 ∧
      ((1 * 2 + 9 * 4) / 10) / ((9 * 1 + 1 * 2) / 10 : ℝ) = 38 / 11 ∧ (2 : ℝ) < 38 / 11 := by
  norm_num

/-- With equal weights the weighted mean is the unweighted mean. -/
theorem wmean_const {ι : Type*} (s : Finset ι) (v : ι → ℝ) (c : ℝ) (hc : c ≠ 0)
    (hs : s.Nonempty) :
    (∑ i ∈ s, c * v i) / (∑ _i ∈ s, c) = (∑ i ∈ s, v i) / s.card := by
  have hn : (s.card : ℝ) ≠ 0 := by exact_mod_cast hs.card_pos.ne'
  rw [← mul_sum, sum_const, nsmul_eq_mul]
  field_simp <;> ring

/-- The weighted mean (weights `w`, say placements) less the unweighted
mean is the covariance of weights and values over the mean weight: the
weighted mean is the larger exactly when the two covary positively. -/
theorem wmean_sub (n sw swv sv : ℝ) (hn : n ≠ 0) (hw : sw ≠ 0) :
    swv / sw - sv / n = (swv / n - (sw / n) * (sv / n)) / (sw / n) := by
  field_simp

end Means

section Alerts

variable {γ : Type*} [DecidableEq γ]

/-- Alert complexity: the distinct configurations among a person-year's
records. -/
def acx (l : List γ) : ℕ := l.toFinset.card

theorem acx_nil : acx ([] : List γ) = 0 := by simp [acx]

/-- One record, or several in one configuration: complexity 1. -/
theorem acx_replicate (n : ℕ) (hn : 0 < n) (c : γ) : acx (List.replicate n c) = 1 := by
  unfold acx
  rw [List.toFinset_replicate_of_ne_zero hn.ne', card_singleton]

theorem acx_pos (l : List γ) (h : l ≠ []) : 0 < acx l := by
  unfold acx
  obtain ⟨a, ha⟩ := List.exists_mem_of_ne_nil l h
  exact card_pos.mpr ⟨a, List.mem_toFinset.mpr ha⟩

theorem acx_le_length (l : List γ) : acx l ≤ l.length := List.toFinset_card_le l

/-- At most the number of configurations allowed: with every record in a
set `S` of configurations (the six the release shows), at most `|S|`. -/
theorem acx_le_of_subset (l : List γ) (S : Finset γ) (h : ∀ c ∈ l, c ∈ S) : acx l ≤ S.card :=
  card_le_card (fun c hc => h c (List.mem_toFinset.mp hc))

/-- Complexity does not depend on the order of the records. -/
theorem acx_perm {l₁ l₂ : List γ} (h : l₁.Perm l₂) : acx l₁ = acx l₂ := by
  unfold acx; rw [List.toFinset_eq_of_perm _ _ h]

/-- Nor on how the configurations are labelled, for an injective
relabelling. -/
theorem acx_map {δ : Type*} [DecidableEq δ] (f : γ → δ) (hf : Function.Injective f)
    (l : List γ) : acx (l.map f) = acx l := by
  have h : (l.map f).toFinset = l.toFinset.image f := by
    ext x
    simp
  unfold acx
  rw [h, card_image_of_injective _ hf]

/-- More records never lower complexity. -/
theorem acx_mono {l₁ l₂ : List γ} (h : l₁ ⊆ l₂) : acx l₁ ≤ acx l₂ :=
  card_le_card (fun _ hc => List.mem_toFinset.mpr (h (List.mem_toFinset.mp hc)))

/-- Pooling two runs of records: at least the larger complexity, at most
the sum. -/
theorem acx_append (l₁ l₂ : List γ) :
    max (acx l₁) (acx l₂) ≤ acx (l₁ ++ l₂) ∧ acx (l₁ ++ l₂) ≤ acx l₁ + acx l₂ := by
  refine ⟨max_le (acx_mono (List.subset_append_left _ _))
    (acx_mono (List.subset_append_right _ _)), ?_⟩
  unfold acx; rw [List.toFinset_append]; exact card_union_le _ _

/-- The treatment D = 1, complexity of 2 or more, holds exactly when two
records of the person-year are in different configurations. -/
theorem complex_iff (l : List γ) : 2 ≤ acx l ↔ ∃ a ∈ l, ∃ b ∈ l, a ≠ b := by
  unfold acx
  rw [show (2 ≤ l.toFinset.card) ↔ 1 < l.toFinset.card from Iff.rfl, one_lt_card]
  simp only [List.mem_toFinset]

/-- The number of alerts on a record (mental health, suicide risk, suicide
watch). -/
def nal (c : Bool × Bool × Bool) : ℕ :=
  (if c.1 then 1 else 0) + (if c.2.1 then 1 else 0) + (if c.2.2 then 1 else 0)

/-- Counterexample to reading D = 1 as "two or more alerts at once": a
person-year with one record without alerts and one with a mental health
alert only has complexity 2 (D = 1), though no record carries more than
one alert. -/
theorem complex_without_two_alerts :
    acx [(false, false, false), (true, false, false)] = 2 ∧
      ∀ c ∈ [(false, false, false), (true, false, false)], nal c ≤ 1 := by
  decide

/-- And a person-year whose one record carries all three alerts has
complexity 1 (D = 0). -/
theorem three_alerts_not_complex :
    acx [(true, true, true)] = 1 ∧ nal (true, true, true) = 3 := by
  decide

end Alerts

end Crimlab.Volat
