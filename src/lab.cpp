// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, labels: every reading of a contrast reported from a table
// (clab_relchk) and the audit of a categorical variable (clab_catchk).
#include <cfloat>
#include <climits>
#include <cmath>
#include <limits>

#include "crimlab.h"

namespace {

// num / den with a zero denominator made explicit: +Inf over a positive
// numerator, NaN over 0 (the numerators here are never negative).
double ratio(double num, double den) {
  if (den > 0) return num / den;
  return num > 0 ? std::numeric_limits<double>::infinity()
                 : std::numeric_limits<double>::quiet_NaN();
}

}  // namespace

extern "C" {

int clab_relchk(const double *x, int r, int c, double value, double tol,
                int *grp, int *g, int *h, int *e, int *kind, double *v,
                int *match, int *nmatch) {
  if (r < 2 || c < 2) return CLAB_EARG;
  // 3 r c (r + c - 2) readings must fit an int
  const double m = 3.0 * r * c * (static_cast<double>(r) + c - 2);
  if (m > INT_MAX) return CLAB_EARG;
  if (!std::isfinite(value) || !(tol >= 0) || !std::isfinite(tol)) return CLAB_EDOM;
  for (long long k = 0; k < static_cast<long long>(r) * c; ++k) {
    if (!(x[k] >= 0) || !std::isfinite(x[k])) return CLAB_EDOM;
  }
  // every row and every column total positive: no empty group or level
  for (int i = 0; i < r; ++i) {
    double s = 0;
    for (int j = 0; j < c; ++j) s += x[i + static_cast<long long>(r) * j];
    if (!(s > 0)) return CLAB_EDOM;
  }
  for (int j = 0; j < c; ++j) {
    double s = 0;
    for (int i = 0; i < r; ++i) s += x[i + static_cast<long long>(r) * j];
    if (!(s > 0)) return CLAB_EDOM;
  }
  int k = 0, nm = 0;
  for (int t = 0; t < 2; ++t) {
    // t = 0: groups in rows, levels in columns; t = 1: the transpose
    const int ng = t == 0 ? r : c, nl = t == 0 ? c : r;
    // cell (group i, level j) of the table read this way
    auto cell = [&](int i, int j) {
      return t == 0 ? x[i + static_cast<long long>(r) * j]
                    : x[j + static_cast<long long>(r) * i];
    };
    for (int q = 0; q < 3; ++q) {
      for (int hh = 0; hh < ng; ++hh) {
        for (int gg = 0; gg < ng; ++gg) {
          if (gg == hh) continue;
          double tg = 0, th = 0;
          for (int j = 0; j < nl; ++j) {
            tg += cell(gg, j);
            th += cell(hh, j);
          }
          for (int ee = 0; ee < nl; ++ee) {
            const double a = cell(gg, ee), cc = cell(hh, ee);
            const double b = tg - a, d = th - cc;
            double y;
            if (q == 0) {
              y = ratio(a * d, b * cc);
            } else if (q == 1) {
              y = ratio(a * th, cc * tg);
            } else {
              y = a / tg - cc / th;
            }
            grp[k] = t;
            kind[k] = q;
            g[k] = gg;
            h[k] = hh;
            e[k] = ee;
            v[k] = y;
            // within tol, up to a few units in the last place of either
            const double big = std::fmax(std::fabs(y), std::fabs(value));
            const bool hit = std::isfinite(y) &&
                             std::fabs(y - value) <= tol + 4 * DBL_EPSILON * big;
            match[k] = hit ? 1 : 0;
            nm += hit ? 1 : 0;
            ++k;
          }
        }
      }
    }
  }
  *nmatch = nm;
  return CLAB_OK;
}

int clab_catchk(const int *code, int n, int k, const int *fold, int nmin,
                int *cnt, int *twin, int *s, int *work) {
  if (n < 0 || k < 0 || nmin < 1) return CLAB_EARG;
  for (int j = 0; j < k; ++j) {
    if (fold[j] < 0 || fold[j] >= k) return CLAB_EARG;
    cnt[j] = 0;
    work[j] = -1;
  }
  int miss = 0;
  for (int i = 0; i < n; ++i) {
    if (code[i] == -1) {
      ++miss;
    } else if (code[i] < 0 || code[i] >= k) {
      return CLAB_EARG;
    } else {
      ++cnt[code[i]];
    }
  }
  // work[f]: the first level with key f; a level's twin is the first other
  // level with its key (the first level of a key points at the second)
  for (int j = 0; j < k; ++j) {
    twin[j] = -1;
    if (work[fold[j]] < 0) {
      work[fold[j]] = j;
    } else {
      twin[j] = work[fold[j]];
      if (twin[twin[j]] < 0) twin[twin[j]] = j;
    }
  }
  int unused = 0, low = 0, tw = 0;
  for (int j = 0; j < k; ++j) {
    unused += cnt[j] == 0 ? 1 : 0;
    low += cnt[j] > 0 && cnt[j] < nmin ? 1 : 0;
    tw += twin[j] >= 0 ? 1 : 0;
  }
  s[0] = miss;
  s[1] = unused;
  s[2] = low;
  s[3] = tw;
  return CLAB_OK;
}

}  // extern "C"
