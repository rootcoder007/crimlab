/-
  crimlab: inst/lean/Network.lean

  Identities and bounds behind the network functions.

    * handshake (ntdeg): in a finite simple graph the degrees sum to twice
      the number of edges.
    * modularity_le_one, modularity_ge_neg_half (ntmod): modularity at
      resolution 1 lies in [-1/2, 1]. A partition's modularity is
      ∑ e − ∑ (e + x)², where e c is the share of the total weight that
      lies inside community c (counted from both ends) and x c the share
      of edge ends in c whose other end lies outside; every leaving end of
      c is matched by an end in another community, so x c is at most the
      sum of the others' x (Brandes et al. 2008,
      doi:10.1109/TKDE.2007.190689, for the lower bound).
    * core_dense, core_nested (ntcore): the k-core, the union of all node
      sets in which every member has at least k neighbours inside, is such
      a set itself, and the (k + 1)-core lies in the k-core.
    * cheeger_lower (ntchg): the easy half of Cheeger's inequality,
      λ₂ / 2 ≤ φ(S) for every set S of positive volume on both sides
      (Chung 1997, Spectral Graph Theory, ch. 2), where λ₂ is any number
      with λ₂ · ∑ d f² ≤ ½ ∑∑ w (f u − f v)² for every f with ∑ d f = 0
      (the second eigenvalue of the normalised Laplacian is the largest
      such number). The test vector is 1/vol S on S and −1/vol Sᶜ off it.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section
open Finset

namespace Crimlab.Network

/-- Handshake: the degrees of a finite simple graph sum to twice its edges. -/
theorem handshake {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    [DecidableRel G.Adj] : ∑ v, G.degree v = 2 * #G.edgeFinset :=
  G.sum_degrees_eq_twice_card_edges

section Modularity

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Modularity at resolution 1 from the inside shares `e` and leaving shares `x`. -/
def modq (e x : ι → ℝ) : ℝ := ∑ c, e c - ∑ c, (e c + x c) ^ 2

omit [DecidableEq ι] in
/-- Modularity is at most 1. -/
theorem modularity_le_one (e x : ι → ℝ) (hx : ∀ c, 0 ≤ x c)
    (h1 : ∑ c, (e c + x c) = 1) : modq e x ≤ 1 := by
  unfold modq
  have hs : ∑ c, e c ≤ ∑ c, (e c + x c) :=
    sum_le_sum fun c _ => by linarith [hx c]
  have hq : 0 ≤ ∑ c, (e c + x c) ^ 2 := sum_nonneg fun c _ => sq_nonneg _
  linarith

/-- Modularity is at least −1/2. -/
theorem modularity_ge_neg_half (e x : ι → ℝ) (he : ∀ c, 0 ≤ e c) (hx : ∀ c, 0 ≤ x c)
    (h1 : ∑ c, (e c + x c) = 1) (hcut : ∀ c, x c ≤ ∑ d ∈ univ.erase c, x d) :
    -1 / 2 ≤ modq e x := by
  unfold modq
  set E := ∑ c, e c with hE
  set X := ∑ c, x c with hX
  have hEX : E + X = 1 := by rw [hE, hX, ← sum_add_distrib]; exact h1
  have hE0 : 0 ≤ E := sum_nonneg fun c _ => he c
  -- every community's share of ends is at most (1 + E) / 2
  have hm : ∀ c, e c + x c ≤ (1 + E) / 2 := by
    intro c
    have h2 := hcut c
    rw [sum_erase_eq_sub (mem_univ c)] at h2
    have h3 : e c ≤ E := single_le_sum (fun d _ => he d) (mem_univ c)
    linarith
  have hsq : ∑ c, (e c + x c) ^ 2 ≤ ∑ c, (e c + x c) * ((1 + E) / 2) :=
    sum_le_sum fun c _ => by
      rw [sq]
      exact mul_le_mul_of_nonneg_left (hm c) (by linarith [he c, hx c])
  rw [← sum_mul, h1, one_mul] at hsq
  linarith

end Modularity

section Core

variable {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V) [DecidableRel G.Adj]

/-- Every member of `S` has at least `k` neighbours in `S`. -/
def Dense (k : ℕ) (S : Finset V) : Prop := ∀ v ∈ S, k ≤ #(S.filter (G.Adj v))

open Classical in
/-- The k-core: every node that lies in some k-dense set. -/
def core (k : ℕ) : Finset V := univ.filter fun v => ∃ S, Dense G k S ∧ v ∈ S

/-- The k-core is k-dense: the union of k-dense sets is k-dense. -/
theorem core_dense (k : ℕ) : Dense G k (core G k) := by
  classical
  intro v hv
  simp only [core, mem_filter, mem_univ, true_and] at hv
  obtain ⟨S, hS, hvS⟩ := hv
  have hsub : S ⊆ core G k := fun u hu => by
    simp only [core, mem_filter, mem_univ, true_and]
    exact ⟨S, hS, hu⟩
  calc k ≤ #(S.filter (G.Adj v)) := hS v hvS
    _ ≤ #((core G k).filter (G.Adj v)) := card_le_card (filter_subset_filter _ hsub)

/-- Every k-dense set lies in the k-core. -/
theorem dense_sub_core (k : ℕ) (S : Finset V) (hS : Dense G k S) : S ⊆ core G k := by
  classical
  intro u hu
  simp only [core, mem_filter, mem_univ, true_and]
  exact ⟨S, hS, hu⟩

/-- The cores are nested: the (k + 1)-core lies in the k-core. -/
theorem core_nested (k : ℕ) : core G (k + 1) ⊆ core G k :=
  dense_sub_core G k _ fun v hv => Nat.le_of_succ_le (core_dense G (k + 1) v hv)

end Core

section Cheeger

variable {V : Type*} [Fintype V] [DecidableEq V]

/-- A node's strength. -/
def deg (w : V → V → ℝ) (v : V) : ℝ := ∑ u, w v u

/-- The volume of a node set. -/
def vol (w : V → V → ℝ) (S : Finset V) : ℝ := ∑ v ∈ S, deg w v

/-- The weight leaving a node set. -/
def cut (w : V → V → ℝ) (S : Finset V) : ℝ := ∑ u ∈ S, ∑ v ∈ Sᶜ, w u v

/-- The Laplacian quadratic form, ½ ∑∑ w (f u − f v)². -/
def energy (w : V → V → ℝ) (f : V → ℝ) : ℝ := (∑ u, ∑ v, w u v * (f u - f v) ^ 2) / 2

/-- The degree-weighted squared norm. -/
def dnorm (w : V → V → ℝ) (f : V → ℝ) : ℝ := ∑ v, deg w v * f v ^ 2

/-- Split a sum over all nodes into the part in `S` and the part outside. -/
lemma split (S : Finset V) (g : V → ℝ) : ∑ v, g v = ∑ v ∈ S, g v + ∑ v ∈ Sᶜ, g v :=
  (sum_add_sum_compl S g).symm

/-- The easy half of Cheeger's inequality: λ₂ / 2 is at most the conductance of every set. -/
theorem cheeger_lower (w : V → V → ℝ) (hs : ∀ u v, w u v = w v u) (hn : ∀ u v, 0 ≤ w u v)
    (lam : ℝ) (hlam : ∀ f : V → ℝ, ∑ v, deg w v * f v = 0 → lam * dnorm w f ≤ energy w f)
    (S : Finset V) (ha : 0 < vol w S) (hb : 0 < vol w Sᶜ) :
    lam / 2 ≤ cut w S / min (vol w S) (vol w Sᶜ) := by
  set a := vol w S with ha'
  set b := vol w Sᶜ with hb'
  set f : V → ℝ := fun v => if v ∈ S then 1 / a else -(1 / b) with hf
  have fin : ∀ v ∈ S, f v = 1 / a := fun v hv => by simp [hf, hv]
  have fout : ∀ v ∈ Sᶜ, f v = -(1 / b) := fun v hv => by
    simp [hf, (mem_compl.mp hv)]
  -- the test vector is orthogonal to the degrees
  have va : ∑ v ∈ S, deg w v = a := rfl
  have vb : ∑ v ∈ Sᶜ, deg w v = b := rfl
  have h0 : ∑ v, deg w v * f v = 0 := by
    have s1 : ∑ v ∈ S, deg w v * f v = (∑ v ∈ S, deg w v) * (1 / a) := by
      rw [sum_mul]; exact sum_congr rfl fun v hv => by rw [fin v hv]
    have s2 : ∑ v ∈ Sᶜ, deg w v * f v = (∑ v ∈ Sᶜ, deg w v) * -(1 / b) := by
      rw [sum_mul]; exact sum_congr rfl fun v hv => by rw [fout v hv]
    rw [split S, s1, s2, va, vb]
    field_simp
    ring
  -- its degree norm
  have hd : dnorm w f = 1 / a + 1 / b := by
    unfold dnorm
    have s1 : ∑ v ∈ S, deg w v * f v ^ 2 = (∑ v ∈ S, deg w v) * (1 / a) ^ 2 := by
      rw [sum_mul]; exact sum_congr rfl fun v hv => by rw [fin v hv]
    have s2 : ∑ v ∈ Sᶜ, deg w v * f v ^ 2 = (∑ v ∈ Sᶜ, deg w v) * (-(1 / b)) ^ 2 := by
      rw [sum_mul]; exact sum_congr rfl fun v hv => by rw [fout v hv]
    rw [split S, s1, s2, va, vb]
    field_simp
  -- its energy: only edges across the cut count
  have hcut' : ∑ u ∈ Sᶜ, ∑ v ∈ S, w u v = cut w S := by
    unfold cut
    rw [sum_comm]
    exact sum_congr rfl fun u _ => sum_congr rfl fun v _ => hs v u
  have he : energy w f = cut w S * (1 / a + 1 / b) ^ 2 := by
    unfold energy
    have inner : ∀ u, ∑ v, w u v * (f u - f v) ^ 2 =
        ∑ v ∈ S, w u v * (f u - f v) ^ 2 + ∑ v ∈ Sᶜ, w u v * (f u - f v) ^ 2 :=
      fun u => split S _
    rw [sum_congr rfl fun u _ => inner u, split S]
    have b1 : ∑ u ∈ S, (∑ v ∈ S, w u v * (f u - f v) ^ 2 +
        ∑ v ∈ Sᶜ, w u v * (f u - f v) ^ 2) = cut w S * (1 / a + 1 / b) ^ 2 := by
      unfold cut
      rw [sum_mul]
      refine sum_congr rfl fun u hu => ?_
      rw [sum_mul]
      have z : ∑ v ∈ S, w u v * (f u - f v) ^ 2 = 0 :=
        sum_eq_zero fun v hv => by rw [fin u hu, fin v hv]; ring
      rw [z, zero_add]
      exact sum_congr rfl fun v hv => by rw [fin u hu, fout v hv]; ring
    have b2 : ∑ u ∈ Sᶜ, (∑ v ∈ S, w u v * (f u - f v) ^ 2 +
        ∑ v ∈ Sᶜ, w u v * (f u - f v) ^ 2) = cut w S * (1 / a + 1 / b) ^ 2 := by
      rw [← hcut', sum_mul]
      refine sum_congr rfl fun u hu => ?_
      rw [sum_mul]
      have z : ∑ v ∈ Sᶜ, w u v * (f u - f v) ^ 2 = 0 :=
        sum_eq_zero fun v hv => by rw [fout u hu, fout v hv]; ring
      rw [z, add_zero]
      exact sum_congr rfl fun v hv => by rw [fout u hu, fin v hv]; ring
    rw [b1, b2]
    ring
  have hk : 0 < 1 / a + 1 / b := by positivity
  have hc0 : 0 ≤ cut w S := sum_nonneg fun u _ => sum_nonneg fun v _ => hn u v
  -- λ₂ ≤ cut (1/a + 1/b)
  have h1 : lam ≤ cut w S * (1 / a + 1 / b) := by
    have := hlam f h0
    rw [hd, he] at this
    nlinarith
  -- and 1/a + 1/b ≤ 2 / min a b
  have hmin : 0 < min a b := lt_min ha hb
  have h2 : 1 / a + 1 / b ≤ 2 / min a b := by
    have x1 : 1 / a ≤ 1 / min a b := one_div_le_one_div_of_le hmin (min_le_left a b)
    have x2 : 1 / b ≤ 1 / min a b := one_div_le_one_div_of_le hmin (min_le_right a b)
    have : 2 / min a b = 1 / min a b + 1 / min a b := by ring
    linarith
  have h3 : cut w S * (1 / a + 1 / b) ≤ cut w S * (2 / min a b) :=
    mul_le_mul_of_nonneg_left h2 hc0
  have h4 : cut w S * (2 / min a b) = 2 * (cut w S / min a b) := by ring
  linarith

end Cheeger

end Crimlab.Network
