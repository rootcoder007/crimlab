// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_relchk and clab_catchk (src/lab.cpp). Part of
// the self-test (selftest.cpp).
#include <cmath>
#include <limits>

#include "crimlab.h"
#include "tchk.h"

namespace {

constexpr int cap = 3 * 6 * 6 * 10;  // readings of a table up to 6 x 6

// 1 for true, 0 for false: counts of failures without bool arithmetic
int one(bool b) { return b ? 1 : 0; }

struct reads {
  int grp[cap], g[cap], h[cap], e[cap], kind[cap], match[cap], nm;
  double v[cap];
};

int relchk(const double *x, int r, int c, double value, double tol, reads &o) {
  return clab_relchk(x, r, c, value, tol, o.grp, o.g, o.h, o.e, o.kind, o.v,
                     o.match, &o.nm);
}

// the index of reading (t, q, h, g, e) of an r x c table, in the documented
// order: orientation, kind, reference, group, level
int at(int r, int c, int t, int q, int h, int g, int e) {
  const int ng = t == 0 ? r : c, nl = t == 0 ? c : r;
  const int base = t == 0 ? 0 : 3 * r * (r - 1) * c;
  const int gi = g < h ? g : g - 1;
  return base + ((q * ng + h) * (ng - 1) + gi) * nl + e;
}

// a contrast from the four cells, by its textbook formula
double formula(int q, double a, double b, double c, double d) {
  const double inf = std::numeric_limits<double>::infinity();
  if (q == 0) return b * c > 0 ? (a * d) / (b * c) : (a * d > 0 ? inf : NAN);
  if (q == 1) {
    const double pg = a / (a + b), ph = c / (c + d);
    return ph > 0 ? pg / ph : (pg > 0 ? inf : NAN);
  }
  return a / (a + b) - c / (c + d);
}

void rel(tk &t) {
  static reads o;
  // 1. A: 20 events, 80 not; B: 10 and 90. Every reading by hand.
  {
    const double x[4] = {20, 10, 80, 90};
    tk_true(t, "rel.2x2.rc", relchk(x, 2, 2, 2.25, 0, o) == CLAB_OK);
    // rows: A against B with the event first, then second; B against A
    tk_eq(t, "rel.or.ab1", o.v[at(2, 2, 0, 0, 1, 0, 0)], 2.25, 1e-15);
    tk_eq(t, "rel.or.ab2", o.v[at(2, 2, 0, 0, 1, 0, 1)], 4.0 / 9, 1e-15);
    tk_eq(t, "rel.or.ba1", o.v[at(2, 2, 0, 0, 0, 1, 0)], 4.0 / 9, 1e-15);
    tk_eq(t, "rel.or.ba2", o.v[at(2, 2, 0, 0, 0, 1, 1)], 2.25, 1e-15);
    tk_eq(t, "rel.rr.ab1", o.v[at(2, 2, 0, 1, 1, 0, 0)], 2, 1e-15);
    tk_eq(t, "rel.rr.ab2", o.v[at(2, 2, 0, 1, 1, 0, 1)], 8.0 / 9, 1e-15);
    tk_eq(t, "rel.rr.ba1", o.v[at(2, 2, 0, 1, 0, 1, 0)], 0.5, 1e-15);
    tk_eq(t, "rel.rr.ba2", o.v[at(2, 2, 0, 1, 0, 1, 1)], 1.125, 1e-15);
    tk_eq(t, "rel.rd.ab1", o.v[at(2, 2, 0, 2, 1, 0, 0)], 0.1, 1e-15);
    tk_eq(t, "rel.rd.ab2", o.v[at(2, 2, 0, 2, 1, 0, 1)], -0.1, 1e-15);
    // columns: the event (30) against the rest (170), A as the level
    tk_eq(t, "rel.col.rr", o.v[at(2, 2, 1, 1, 1, 0, 0)], 17.0 / 12, 1e-15);
    tk_eq(t, "rel.col.rd", o.v[at(2, 2, 1, 2, 1, 0, 0)], 2.0 / 3 - 8.0 / 17, 1e-15);
    tk_eq(t, "rel.col.or", o.v[at(2, 2, 1, 0, 1, 0, 0)], 2.25, 1e-15);
    // 2.25 is the odds ratio of A against B for the first level, of B
    // against A for the second, in both orientations: four readings
    tk_ieq(t, "rel.nm", o.nm, 4);
    tk_true(t, "rel.which", o.match[at(2, 2, 0, 0, 1, 0, 0)] == 1 &&
                                o.match[at(2, 2, 0, 0, 0, 1, 1)] == 1 &&
                                o.match[at(2, 2, 1, 0, 1, 0, 0)] == 1 &&
                                o.match[at(2, 2, 1, 0, 0, 1, 1)] == 1);
    // the layout: 24 readings, orientation, kind, reference, group, level
    tk_true(t, "rel.layout", o.grp[11] == 0 && o.grp[12] == 1 && o.kind[4] == 1 &&
                                 o.h[2] == 1 && o.g[2] == 0 && o.e[3] == 1);
  }
  // 2. the algebra of a 2 x 2 table on random tables (Lean: RelChk.lean):
  // the odds ratio inverts when the reference or the event changes and is
  // unchanged when both do or the table is transposed; the risk ratio
  // inverts with the reference; the risk difference changes sign with
  // either
  {
    rng g = {7};
    int bad = 0;
    for (int rep = 0; rep < 500; ++rep) {
      double x[4];
      for (double &u : x) u = 1 + rsam(g, 0, 60) + (rep % 2 ? runi(g) : 0.0);
      relchk(x, 2, 2, 1, 0, o);
      for (int tt = 0; tt < 2; ++tt) {
        const double or1 = o.v[at(2, 2, tt, 0, 1, 0, 0)], or2 = o.v[at(2, 2, tt, 0, 1, 0, 1)];
        const double or3 = o.v[at(2, 2, tt, 0, 0, 1, 0)], or4 = o.v[at(2, 2, tt, 0, 0, 1, 1)];
        bad += one(std::fabs(or1 * or2 - 1) > 1e-14);
        bad += one(std::fabs(or1 * or3 - 1) > 1e-14);
        bad += one(std::fabs(or1 / or4 - 1) > 1e-14);
        bad += one(std::fabs(or1 / o.v[at(2, 2, 1 - tt, 0, 1, 0, 0)] - 1) > 1e-14);
        const double rr1 = o.v[at(2, 2, tt, 1, 1, 0, 0)], rr3 = o.v[at(2, 2, tt, 1, 0, 1, 0)];
        bad += one(std::fabs(rr1 * rr3 - 1) > 1e-14);
        const double rd1 = o.v[at(2, 2, tt, 2, 1, 0, 0)];
        bad += one(std::fabs(rd1 + o.v[at(2, 2, tt, 2, 0, 1, 0)]) > 1e-14);
        bad += one(std::fabs(rd1 + o.v[at(2, 2, tt, 2, 1, 0, 1)]) > 1e-14);
      }
    }
    tk_ieq(t, "rel.algebra", bad, 0);
  }
  // 3. brute force: random tables up to 6 x 6, zero cells included; every
  // reading against the formula on the four cells, read off a transposed
  // copy for the columns; matches counted again
  {
    rng g = {11};
    int bad = 0, badm = 0;
    for (int rep = 0; rep < 300; ++rep) {
      const int r = rsam(g, 2, 6), c = rsam(g, 2, 6);
      double x[36], tx[36];
      for (int k = 0; k < r * c; ++k) x[k] = rsam(g, 0, 4) == 0 ? 0 : rsam(g, 1, 30);
      for (int i = 0; i < r; ++i) x[i + r * (i % c)] += 1;  // rows not empty
      for (int j = 0; j < c; ++j) x[j % r + r * j] += 1;    // columns not empty
      for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j) tx[j + c * i] = x[i + r * j];
      const double value = 1.5;
      if (relchk(x, r, c, value, 0.25, o) != CLAB_OK) ++bad;
      int nm = 0;
      for (int tt = 0; tt < 2; ++tt) {
        const int ng = tt == 0 ? r : c, nl = tt == 0 ? c : r;
        const double *y = tt == 0 ? x : tx;
        for (int q = 0; q < 3; ++q)
          for (int h = 0; h < ng; ++h)
            for (int gg = 0; gg < ng; ++gg) {
              if (gg == h) continue;
              for (int e = 0; e < nl; ++e) {
                double tg = 0, th = 0;
                for (int j = 0; j < nl; ++j) {
                  tg += y[gg + ng * j];
                  th += y[h + ng * j];
                }
                const double a = y[gg + ng * e], cc = y[h + ng * e];
                const double want = formula(q, a, tg - a, cc, th - cc);
                const int k = at(r, c, tt, q, h, gg, e);
                const double got = o.v[k];
                const bool same = (std::isnan(want) && std::isnan(got)) || got == want ||
                                  std::fabs(got - want) <= 1e-14 * std::fabs(want);
                bad += one(!same || o.grp[k] != tt || o.kind[k] != q || o.h[k] != h ||
                       o.g[k] != gg || o.e[k] != e);
                const int hit = one(std::isfinite(want) && std::fabs(want - value) <= 0.25 + 1e-12);
                nm += hit;
                badm += one(hit != o.match[k]);
              }
            }
      }
      badm += one(nm != o.nm);
    }
    tk_ieq(t, "rel.brute", bad, 0);
    tk_ieq(t, "rel.brute.match", badm, 0);
  }
  // 4. rounding: a reading matches a value reported to two decimals when it
  // rounds to it either way at the half (2.345 to 2.34 or 2.35), and not
  // when it lies beyond the half
  {
    const double x1[4] = {469, 1, 200, 1}, x2[4] = {23451, 1, 10000, 1},
                 x3[4] = {23449, 1, 10000, 1};
    const int k = at(2, 2, 0, 0, 1, 0, 0);
    relchk(x1, 2, 2, 2.35, 0.005, o);
    tk_true(t, "rel.half.up", o.match[k] == 1);
    relchk(x1, 2, 2, 2.34, 0.005, o);
    tk_true(t, "rel.half.down", o.match[k] == 1);
    relchk(x1, 2, 2, 2.33, 0.005, o);
    tk_true(t, "rel.half.out", o.match[k] == 0);
    relchk(x2, 2, 2, 2.35, 0.005, o);
    tk_true(t, "rel.above.in", o.match[k] == 1);
    relchk(x2, 2, 2, 2.34, 0.005, o);
    tk_true(t, "rel.above.out", o.match[k] == 0);
    relchk(x3, 2, 2, 2.34, 0.005, o);
    tk_true(t, "rel.below.in", o.match[k] == 1);
    relchk(x3, 2, 2, 2.35, 0.005, o);
    tk_true(t, "rel.below.out", o.match[k] == 0);
    // tol 0: equal up to a few units in the last place
    relchk(x1, 2, 2, 2.345, 0, o);
    tk_true(t, "rel.exact", o.match[k] == 1);
    relchk(x1, 2, 2, 2.345 + 1e-12, 0, o);
    tk_true(t, "rel.exact.out", o.match[k] == 0);
  }
  // 5. zero cells: a positive count over 0 is +Inf, 0 over 0 NaN; neither
  // matches any value
  {
    const double x[6] = {0, 0, 5, 5, 5, 0};  // A: 0, 5; B: 0, 5; C: 5, 0
    relchk(x, 3, 2, 0, 1e300, o);
    tk_nan(t, "rel.zero.or", o.v[at(3, 2, 0, 0, 1, 0, 0)]);
    tk_nan(t, "rel.zero.rr", o.v[at(3, 2, 0, 1, 1, 0, 0)]);
    tk_eq(t, "rel.zero.rd", o.v[at(3, 2, 0, 2, 1, 0, 0)], 0, 0);
    tk_true(t, "rel.zero.inf", std::isinf(o.v[at(3, 2, 0, 0, 0, 2, 0)]) &&
                                   o.v[at(3, 2, 0, 0, 0, 2, 0)] > 0);
    tk_true(t, "rel.zero.nomatch", o.match[at(3, 2, 0, 0, 1, 0, 0)] == 0 &&
                                       o.match[at(3, 2, 0, 0, 0, 2, 0)] == 0);
    // a difference of exactly 0 matches 0 with no tolerance at all
    relchk(x, 3, 2, 0, 0, o);
    tk_true(t, "rel.zero.exact", o.match[at(3, 2, 0, 2, 1, 0, 0)] == 1);
  }
  // 6. guards
  {
    const double x[4] = {1, 2, 3, 4}, neg[4] = {1, -2, 3, 4}, row0[4] = {0, 2, 0, 4},
                 col0[4] = {0, 0, 3, 4}, nan[4] = {1, NAN, 3, 4};
    tk_ieq(t, "rel.guard.r", relchk(x, 1, 4, 1, 0, o), CLAB_EARG);
    tk_ieq(t, "rel.guard.c", relchk(x, 4, 1, 1, 0, o), CLAB_EARG);
    tk_ieq(t, "rel.guard.big",
           clab_relchk(x, 50000, 50000, 1, 0, nullptr, nullptr, nullptr, nullptr,
                       nullptr, nullptr, nullptr, nullptr),
           CLAB_EARG);
    tk_ieq(t, "rel.guard.neg", relchk(neg, 2, 2, 1, 0, o), CLAB_EDOM);
    tk_ieq(t, "rel.guard.nan", relchk(nan, 2, 2, 1, 0, o), CLAB_EDOM);
    tk_ieq(t, "rel.guard.row", relchk(row0, 2, 2, 1, 0, o), CLAB_EDOM);
    tk_ieq(t, "rel.guard.col", relchk(col0, 2, 2, 1, 0, o), CLAB_EDOM);
    tk_ieq(t, "rel.guard.tol", relchk(x, 2, 2, 1, -1, o), CLAB_EDOM);
    tk_ieq(t, "rel.guard.value", relchk(x, 2, 2, NAN, 0, o), CLAB_EDOM);
  }
}

void cat(tk &t) {
  int cnt[8], twin[8], s[4], work[8];
  // 1. by hand: levels a, b, c, A (A folds onto a), d unused; 2 missing
  {
    const int code[9] = {0, 1, -1, 0, 3, 1, 1, -1, 2};
    const int fold[5] = {0, 1, 2, 0, 4};
    tk_ieq(t, "cat.rc", clab_catchk(code, 9, 5, fold, 2, cnt, twin, s, work), CLAB_OK);
    tk_true(t, "cat.cnt", cnt[0] == 2 && cnt[1] == 3 && cnt[2] == 1 && cnt[3] == 1 &&
                              cnt[4] == 0);
    tk_true(t, "cat.twin", twin[0] == 3 && twin[3] == 0 && twin[1] == -1 &&
                               twin[2] == -1 && twin[4] == -1);
    tk_true(t, "cat.sum", s[0] == 2 && s[1] == 1 && s[2] == 2 && s[3] == 2);
  }
  // 2. three labels with one key: each points at the first, the first at
  // the second
  {
    const int code[3] = {0, 1, 2}, fold[3] = {1, 1, 1};
    clab_catchk(code, 3, 3, fold, 1, cnt, twin, s, work);
    tk_true(t, "cat.three", twin[0] == 1 && twin[1] == 0 && twin[2] == 0 && s[3] == 3 &&
                                s[2] == 0);
  }
  // 3. brute force on random codes against direct counts
  {
    rng g = {5};
    int bad = 0;
    for (int rep = 0; rep < 200; ++rep) {
      const int k = rsam(g, 0, 8), n = rsam(g, 0, 40), nmin = rsam(g, 1, 6);
      int code[40], fold[8];
      for (int i = 0; i < n; ++i) code[i] = k == 0 ? -1 : rsam(g, -1, k - 1);
      for (int j = 0; j < k; ++j) fold[j] = rsam(g, 0, k - 1);
      clab_catchk(code, n, k, fold, nmin, cnt, twin, s, work);
      int miss = 0, un = 0, lo = 0, tw = 0;
      for (int i = 0; i < n; ++i) miss += one(code[i] < 0);
      for (int j = 0; j < k; ++j) {
        int c = 0, first = -1, other = -1;
        for (int i = 0; i < n; ++i) c += one(code[i] == j);
        for (int m = 0; m < k; ++m) {
          if (fold[m] != fold[j]) continue;
          if (first < 0) first = m;
          if (m != j && other < 0) other = m;
        }
        const int want = first == j ? other : first;
        bad += one(c != cnt[j] || want != twin[j]);
        un += one(c == 0);
        lo += one(c > 0 && c < nmin);
        tw += one(want >= 0);
      }
      bad += one(miss != s[0] || un != s[1] || lo != s[2] || tw != s[3]);
    }
    tk_ieq(t, "cat.brute", bad, 0);
  }
  // 4. guards and empty input
  {
    const int code[2] = {0, 2}, fold[2] = {0, 1}, badf[2] = {0, 2};
    tk_ieq(t, "cat.guard.code", clab_catchk(code, 2, 2, fold, 1, cnt, twin, s, work), CLAB_EARG);
    tk_ieq(t, "cat.guard.fold", clab_catchk(code, 1, 2, badf, 1, cnt, twin, s, work), CLAB_EARG);
    tk_ieq(t, "cat.guard.nmin", clab_catchk(code, 1, 2, fold, 0, cnt, twin, s, work), CLAB_EARG);
    tk_ieq(t, "cat.guard.n", clab_catchk(code, -1, 2, fold, 1, cnt, twin, s, work), CLAB_EARG);
    tk_true(t, "cat.empty", clab_catchk(code, 0, 0, fold, 1, cnt, twin, s, work) == CLAB_OK &&
                                s[0] == 0 && s[1] == 0 && s[2] == 0 && s[3] == 0);
  }
}

}  // namespace

void t_lab(tk &t, int) {
  rel(t);
  cat(t);
}
