/-
  crimlab: inst/lean/Custody.lean

  Lakner's stock and flow measures for a custodial population (stkflw(),
  lostay()) and the arithmetic of year-over-year tables (yoysmry(),
  cxyoy()).

    * adp and alos as Lakner (1976) defines them: the average daily
      population is the person-days over the days of the period (a stock),
      the average length of stay the person-days over the people released
      (a flow); adp = N alos / t when the N people are those whose days are
      counted (Little's law, L = lambda W, on a finite period: Little 1961);
    * the expansion estimator of person-days from C periodic headcounts,
      t times their mean, whose adp is the mean headcount, and which is
      exact when the population does not move between counts;
    * the worked examples of Lakner's manual: 13,500 person-days over a
      year is an adp of about 37; 12,150 days for 2,700 people an alos of
      4.5; an adp of 25 with 1,750 admissions an alos of about 5.2 days;
      a mean headcount of 34,935 / 255 = 137 over 365 days is 50,005
      person-days;
    * the stock and the flow can move in opposite directions: at the OTIS
      b02 counts for fiscal 2023 and 2025 the people fell while the
      length of stay and the average daily population rose;
    * a percent change between two years, its chaining, the change of a
      rate as a rate ratio less one, and what a first-to-last change
      column means when a year is added: it is the old column chained with
      the change into the new year, and equals the old column only when
      the new year's count equals the old last year's.

  Lean checks the arithmetic. Whether the people counted are those released
  within the period (Lakner's N'), and whether a period is long enough
  that stays are not cut at its ends, are facts about the records.
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

namespace Crimlab.Custody

/-- Average daily population: person-days over the days of the period. -/
noncomputable def adp (days t : ℝ) : ℝ := days / t

/-- Average length of stay: person-days over the people released. -/
noncomputable def alos (days n : ℝ) : ℝ := days / n

/-- The stock is the flow scaled by throughput: adp = N alos / t. -/
theorem adp_eq (days n t : ℝ) (hn : n ≠ 0) (ht : t ≠ 0) :
    adp days t = n * alos days n / t := by
  rw [adp, alos]; field_simp

/-- Admissions implied by a population and a stay: N = adp t / alos. -/
theorem people_eq (days n t : ℝ) (hn : n ≠ 0) (ht : t ≠ 0) (hd : days ≠ 0) :
    n = adp days t * t / alos days n := by
  rw [adp, alos]; field_simp

/-- Lakner's expansion estimator: person-days from `C` periodic
headcounts, the mean headcount times the days of the period. -/
noncomputable def xdays {C : ℕ} (N : Fin C → ℝ) (t : ℝ) : ℝ := (∑ c, N c) / C * t

/-- Its average daily population is the mean headcount. -/
theorem adp_xdays {C : ℕ} (N : Fin C → ℝ) (t : ℝ) (ht : t ≠ 0) :
    adp (xdays N t) t = (∑ c, N c) / C := by
  rw [adp, xdays]; field_simp

/-- With the same headcount at every count it is exact: the headcount
times the days. -/
theorem xdays_const {C : ℕ} (hC : 0 < C) (k t : ℝ) :
    xdays (fun _ : Fin C => k) t = k * t := by
  have : (C : ℝ) ≠ 0 := by exact_mod_cast hC.ne'
  rw [xdays, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  field_simp

/-- Lakner's first example: 13,500 person-days in a year, an adp of
about 37. -/
theorem lakner_adp_example : 36.98 < adp 13500 365 ∧ adp 13500 365 < 37 := by
  rw [adp]; constructor <;> norm_num

/-- 12,150 person-days for 2,700 people: an alos of 4.5 days. -/
theorem lakner_alos_example : alos 12150 2700 = 4.5 := by
  rw [alos]; norm_num

/-- An adp of 25 with 1,750 admissions in a year: a stay of about 5.2
days (25 365 / 1750 = 5.214...). -/
theorem lakner_stay_example : 5.2 < 25 * 365 / (1750 : ℝ) ∧ 25 * 365 / (1750 : ℝ) < 5.22 := by
  constructor <;> norm_num

/-- 255 counts summing to 34,935 over 365 days: 50,005 person-days. -/
theorem lakner_expansion_example : (34935 : ℝ) / 255 * 365 = 50005 := by norm_num

/-- At the OTIS b02 counts for 2023 and 2025 (115,674 and 126,121
person-days; 12,647 and 9,608 people; 365 days each) the people fell and
both the length of stay and the average daily population rose. -/
theorem stock_flow_diverge :
    (9608 : ℝ) < 12647 ∧ alos 115674 12647 < alos 126121 9608 ∧
      adp 115674 365 < adp 126121 365 := by
  simp only [alos, adp]; refine ⟨by norm_num, by norm_num, by norm_num⟩

/-- The relative change from `x` to `y` (times 100 a percent change). -/
noncomputable def rchg (x y : ℝ) : ℝ := y / x - 1

/-- Changes chain: 1 + change(x, z) = (1 + change(x, y)) (1 + change(y, z)). -/
theorem rchg_chain (x y z : ℝ) (hx : x ≠ 0) (hy : y ≠ 0) :
    1 + rchg x z = (1 + rchg x y) * (1 + rchg y z) := by
  simp only [rchg]; field_simp; ring

/-- Days decompose into people and stays: the change in days is the
change in people chained with the change in stays. -/
theorem days_chg (n₁ n₂ l₁ l₂ : ℝ) (hn : n₁ ≠ 0) (hl : l₁ ≠ 0) :
    1 + rchg (n₁ * l₁) (n₂ * l₂) = (1 + rchg n₁ n₂) * (1 + rchg l₁ l₂) := by
  simp only [rchg]; field_simp; ring

/-- The change of a rate is the rate ratio less one. -/
theorem rate_chg (n₁ e₁ n₂ e₂ : ℝ) :
    rchg (n₁ / e₁) (n₂ / e₂) = (n₂ / e₂) / (n₁ / e₁) - 1 := rfl

/-- A first-to-last column after a year is added is the old column
chained with the change into the new year. -/
theorem first_last_add (x₀ xₒ xₙ : ℝ) (h₀ : x₀ ≠ 0) (hₒ : xₒ ≠ 0) :
    rchg x₀ xₙ = (1 + rchg x₀ xₒ) * (1 + rchg xₒ xₙ) - 1 := by
  have := rchg_chain x₀ xₒ xₙ h₀ hₒ; linarith

/-- So it keeps its old value only if the new year's count equals the old
last year's: a stored first-to-last column is a different quantity once a
year is added, even when no earlier count changed. -/
theorem first_last_same_iff (x₀ xₒ xₙ : ℝ) (h₀ : x₀ ≠ 0) :
    rchg x₀ xₙ = rchg x₀ xₒ ↔ xₙ = xₒ := by
  simp only [rchg]
  constructor
  · intro h
    have : xₙ / x₀ = xₒ / x₀ := by linarith
    exact (div_left_inj' h₀).mp this
  · rintro rfl; rfl

end Crimlab.Custody
