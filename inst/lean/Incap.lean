/-
  crimlab: inst/lean/Incap.lean

  Incapacitation: what a sentence length buys in crimes not committed
  (Avi-Itzhak & Shinnar 1973; Shinnar & Shinnar 1975; Blumstein, Cohen &
  Nagin 1978), and the two assumptions behind it, replacement and
  desistance.

  An offender commits crimes at rate λ while free; each crime leads to
  arrest, conviction and a sentence of length S with probability q; so the
  expected free time between sentences is 1/(λq). Over many cycles the
  long-run crime rate is

      λ · (1/(λq)) / (1/(λq) + S)  =  λ / (1 + λ q S)            (steady_state_rate)

    * cycle_rate: on N cycles with free times f_i and crimes λ f_i the
      realised rate is λ f̄ / (f̄ + S).
    * prevented_share_eq, prevented_share_lt_one: the share of the
      free-offending rate removed is λqS/(1 + λqS), strictly below one.
    * rate_antitone_in_S, rate_antitone_in_q: longer sentences and higher
      certainty lower the rate.
    * marginal_prevention_eq: one more year of sentence prevents
      λ²q/((1 + λqS)(1 + λq(S+1))) crimes per year.
    * high_rate_more_prevented: a sentence rule uniform across offenders
      removes the larger share from the higher-rate group.

  Replacement (Reuter 1983; Cook 1986): a share r of the crimes an
  incarcerated offender would have committed is committed by someone else;
  the prevented share becomes (1 − r) times the one above (replaced_le,
  replaced_antitone, replaced_full).

  Desistance (Blumstein, Cohen, Roth & Visher 1986; Laub & Sampson 2003):
  with a non-increasing rate path λ(t), a sentence of S periods from career
  age t₀ prevents ∑_{s<S} λ(t₀+s) crimes, at most S·λ(t₀)
  (prevented_le_const) and at least S·λ(t₀+S) (prevented_ge_const); the
  same sentence later in a career prevents fewer (later_sentence_prevents_less);
  with replacement as well the entry-rate estimate S·λ(t₀) is still an
  upper bound (prevented_net_le).

  What is not proved is that λ is constant, or that incarceration leaves
  it unchanged: those are the empirical assumptions.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

noncomputable section
open Finset

namespace Crimlab.Incap

/-- Offending rate `λ`, sentence probability per crime `q`, sentence length `S`. -/
structure Career where
  (lam q S : ℝ)
  (lam_pos : 0 < lam)
  (q_pos : 0 < q)
  (q_le : q ≤ 1)
  (S_nonneg : 0 ≤ S)

variable (C : Career)

/-- Long-run crime rate under the incapacitation model. -/
def rate : ℝ := C.lam / (1 + C.lam * C.q * C.S)
/-- Share of free-state offending removed. -/
def preventedShare : ℝ := C.lam * C.q * C.S / (1 + C.lam * C.q * C.S)

theorem denom_pos : 0 < 1 + C.lam * C.q * C.S := by
  have := mul_nonneg (mul_nonneg C.lam_pos.le C.q_pos.le) C.S_nonneg; linarith

/-- `λ f̄ / (f̄ + S) = λ/(1 + λqS)` at the mean free time `f̄ = 1/(λq)`. -/
theorem steady_state_rate : C.lam * (1 / (C.lam * C.q)) / (1 / (C.lam * C.q) + C.S) = rate C := by
  unfold rate
  have hl := C.lam_pos.ne'
  have hq := C.q_pos.ne'
  have h2 := (denom_pos C).ne'
  field_simp <;> ring

/-- The realised rate on `N` cycles with free times `f_i` and crimes `λ f_i`. -/
theorem cycle_rate {N : ℕ} (f : Fin N → ℝ) (hN : 0 < N) (hf : ∀ i, 0 < f i) :
    (∑ i, C.lam * f i) / (∑ i, (f i + C.S)) = C.lam * ((∑ i, f i) / N) / ((∑ i, f i) / N + C.S) := by
  have hsum : 0 < ∑ i, f i := sum_pos (fun i _ => hf i) ⟨⟨0, hN⟩, mem_univ _⟩
  have hN' : (N : ℝ) ≠ 0 := by exact_mod_cast hN.ne'
  rw [← mul_sum, sum_add_distrib, sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul]
  have hden : (∑ i, f i) + N * C.S ≠ 0 := by
    have := mul_nonneg (Nat.cast_nonneg N : (0 : ℝ) ≤ N) C.S_nonneg
    linarith
  field_simp <;> ring

theorem rate_pos : 0 < rate C := div_pos C.lam_pos (denom_pos C)

theorem rate_le_lam : rate C ≤ C.lam := by
  unfold rate
  rw [div_le_iff₀ (denom_pos C)]
  have := mul_nonneg (mul_nonneg C.lam_pos.le C.q_pos.le) C.S_nonneg
  nlinarith [C.lam_pos]

/-- The prevented share is `1 − rate/λ`, lies in `[0, 1)`, and never reaches one. -/
theorem prevented_share_eq : preventedShare C = 1 - rate C / C.lam := by
  unfold preventedShare rate
  have hl := C.lam_pos.ne'
  have h2 := (denom_pos C).ne'
  field_simp <;> ring

theorem prevented_share_lt_one : preventedShare C < 1 := by
  unfold preventedShare
  rw [div_lt_one (denom_pos C)]; linarith

theorem prevented_share_nonneg : 0 ≤ preventedShare C :=
  div_nonneg (mul_nonneg (mul_nonneg C.lam_pos.le C.q_pos.le) C.S_nonneg) (denom_pos C).le

/-- Longer sentences lower the rate. -/
theorem rate_antitone_in_S (C₁ C₂ : Career) (hl : C₁.lam = C₂.lam) (hq : C₁.q = C₂.q) (hS : C₁.S ≤ C₂.S) :
    rate C₂ ≤ rate C₁ := by
  unfold rate
  rw [hl, hq]
  apply div_le_div_of_nonneg_left C₂.lam_pos.le (by rw [← hl, ← hq]; exact denom_pos C₁)
  have := mul_pos C₂.lam_pos C₂.q_pos
  nlinarith

/-- Higher certainty lowers the rate. -/
theorem rate_antitone_in_q (C₁ C₂ : Career) (hl : C₁.lam = C₂.lam) (hS : C₁.S = C₂.S) (hq : C₁.q ≤ C₂.q) :
    rate C₂ ≤ rate C₁ := by
  unfold rate
  rw [hl, hS]
  apply div_le_div_of_nonneg_left C₂.lam_pos.le (by rw [← hl, ← hS]; exact denom_pos C₁)
  have := mul_nonneg C₂.lam_pos.le C₂.S_nonneg
  nlinarith

/-- Crimes prevented per year by one more year of sentence, and that it falls in `S`. -/
def marginalPrevention : ℝ := rate C - C.lam / (1 + C.lam * C.q * (C.S + 1))

theorem marginal_prevention_eq :
    marginalPrevention C = C.lam ^ 2 * C.q / ((1 + C.lam * C.q * C.S) * (1 + C.lam * C.q * (C.S + 1))) := by
  unfold marginalPrevention rate
  have h1 := (denom_pos C).ne'
  have h2 : (1 + C.lam * C.q * (C.S + 1)) ≠ 0 := by
    have := mul_nonneg (mul_nonneg C.lam_pos.le C.q_pos.le) (by linarith [C.S_nonneg] : (0 : ℝ) ≤ C.S + 1)
    linarith
  field_simp
  ring

theorem marginal_prevention_pos : 0 < marginalPrevention C := by
  rw [marginal_prevention_eq]
  have := mul_nonneg (mul_nonneg C.lam_pos.le C.q_pos.le) C.S_nonneg
  have h2 : 0 < 1 + C.lam * C.q * (C.S + 1) := by nlinarith [mul_pos C.lam_pos C.q_pos]
  exact div_pos (mul_pos (pow_pos C.lam_pos 2) C.q_pos) (mul_pos (denom_pos C) h2)

/-! ### Heterogeneous offenders -/

/-- Share-weighted aggregate rate of `K` groups. -/
def groupRate {K : ℕ} (π : Fin K → ℝ) (Cs : Fin K → Career) : ℝ := ∑ k, π k * rate (Cs k)

/-- Under a uniform sentence the high-rate group has the larger prevented share. -/
theorem high_rate_more_prevented (C₁ C₂ : Career) (hq : C₁.q = C₂.q) (hS : C₁.S = C₂.S) (hl : C₁.lam ≤ C₂.lam) :
    preventedShare C₁ ≤ preventedShare C₂ := by
  unfold preventedShare
  rw [hq, hS]
  have d1 := denom_pos C₁; have d2 := denom_pos C₂
  rw [hq, hS] at d1
  rw [div_le_div_iff₀ d1 d2]
  have := mul_nonneg C₂.q_pos.le C₂.S_nonneg
  nlinarith


/-! ### Replacement -/

/-- Prevented share net of a replacement fraction `r`. -/
def replacedShare (C : Career) (r : ℝ) : ℝ := (1 - r) * preventedShare C

theorem replaced_le (C : Career) (r : ℝ) (hr : 0 ≤ r) : replacedShare C r ≤ preventedShare C := by
  unfold replacedShare
  exact mul_le_of_le_one_left (prevented_share_nonneg C) (by linarith)

theorem replaced_antitone (C : Career) {r₁ r₂ : ℝ} (h : r₁ ≤ r₂) :
    replacedShare C r₂ ≤ replacedShare C r₁ :=
  mul_le_mul_of_nonneg_right (by linarith) (prevented_share_nonneg C)

theorem replaced_full (C : Career) : replacedShare C 1 = 0 := by simp [replacedShare]

theorem replaced_nonneg (C : Career) (r : ℝ) (hr : r ≤ 1) : 0 ≤ replacedShare C r :=
  mul_nonneg (by linarith) (prevented_share_nonneg C)

/-! ### Desistance -/

/-- Crimes prevented by `S` periods served from career age `t₀` under the rate path `lam`. -/
def prevented (lam : ℕ → ℝ) (t₀ S : ℕ) : ℝ := ∑ s ∈ range S, lam (t₀ + s)

theorem prevented_le_const (lam : ℕ → ℝ) (hanti : Antitone lam) (t₀ S : ℕ) :
    prevented lam t₀ S ≤ S * lam t₀ := by
  unfold prevented
  calc ∑ s ∈ range S, lam (t₀ + s) ≤ ∑ s ∈ range S, lam t₀ :=
        sum_le_sum fun s _ => hanti (Nat.le_add_right t₀ s)
    _ = S * lam t₀ := by simp [sum_const, card_range, nsmul_eq_mul]

theorem prevented_ge_const (lam : ℕ → ℝ) (hanti : Antitone lam) (t₀ S : ℕ) :
    S * lam (t₀ + S) ≤ prevented lam t₀ S := by
  unfold prevented
  calc (S : ℝ) * lam (t₀ + S) = ∑ s ∈ range S, lam (t₀ + S) := by
        simp [sum_const, card_range, nsmul_eq_mul]
    _ ≤ ∑ s ∈ range S, lam (t₀ + s) := by
        apply sum_le_sum; intro s hs
        exact hanti (Nat.add_le_add_left (mem_range.mp hs).le t₀)

theorem later_sentence_prevents_less (lam : ℕ → ℝ) (hanti : Antitone lam) {t₀ t₁ : ℕ}
    (h : t₀ ≤ t₁) (S : ℕ) : prevented lam t₁ S ≤ prevented lam t₀ S := by
  unfold prevented
  exact sum_le_sum fun s _ => hanti (Nat.add_le_add_right h s)

theorem constant_rate (c : ℝ) (t₀ S : ℕ) : prevented (fun _ => c) t₀ S = S * c := by
  simp [prevented, sum_const, card_range, nsmul_eq_mul]

/-- Both corrections together: desistance and replacement. -/
def preventedNet (lam : ℕ → ℝ) (t₀ S : ℕ) (r : ℝ) : ℝ := (1 - r) * prevented lam t₀ S

theorem prevented_net_le (lam : ℕ → ℝ) (hanti : Antitone lam) (hnn : ∀ t, 0 ≤ lam t)
    (t₀ S : ℕ) (r : ℝ) (hr0 : 0 ≤ r) (hr1 : r ≤ 1) :
    preventedNet lam t₀ S r ≤ S * lam t₀ := by
  unfold preventedNet
  have h1 := prevented_le_const lam hanti t₀ S
  have h0 : 0 ≤ prevented lam t₀ S := sum_nonneg fun s _ => hnn _
  calc (1 - r) * prevented lam t₀ S ≤ 1 * prevented lam t₀ S :=
        mul_le_mul_of_nonneg_right (by linarith) h0
    _ ≤ S * lam t₀ := by linarith

end Crimlab.Incap
