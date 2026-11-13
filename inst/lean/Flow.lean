/-
  crimlab: inst/lean/Flow.lean

  Stocks and flows (stkflw(), bklog()).

    * adp_alos: person-days d over a period of t days are the average
      daily population times t and the people times their mean stay, so
      adp = people alos / t: any two of the measures give the others.
    * days_ratio: with days = people alos in two periods, the ratio of
      days is the ratio of people times the ratio of stays, so the log
      change in days splits exactly into the two log changes.
    * little: n cases over a horizon h whose times in the system add up
      to s have a time-average backlog s / h equal to the filing rate n / h
      times the mean time s / n (Little's law on a finite docket with
      every case disposed).
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.style.header false

namespace Crimlab.Flow

/-- The average daily population from people and stays. -/
theorem adp_alos (d adp people alos t : ℝ) (ht : t ≠ 0) (h1 : d = adp * t)
    (h2 : d = people * alos) : adp = people * alos / t := by
  rw [← h2, h1]
  field_simp

/-- Days split into people and stays. -/
theorem days_ratio (p0 p1 l0 l1 : ℝ) (hp : p0 ≠ 0) (hl : l0 ≠ 0) :
    (p1 * l1) / (p0 * l0) = (p1 / p0) * (l1 / l0) := by
  field_simp

/-- Little's law over a horizon. -/
theorem little (s n h : ℝ) (hn : n ≠ 0) (hh : h ≠ 0) : s / h = (n / h) * (s / n) := by
  field_simp

end Crimlab.Flow
