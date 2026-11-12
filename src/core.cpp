// SPDX-License-Identifier: AGPL-3.0-or-later
//
// libcsls: the kernels. No R here, no allocation, no exceptions: every
// function reads the caller's arrays and writes the caller's arrays.
#include "csls.h"

#include <cmath>
#include <limits>

namespace {

inline bool pos(double x) noexcept { return std::isfinite(x) && x > 0.0; }
inline bool nonneg(double x) noexcept { return std::isfinite(x) && x >= 0.0; }
inline bool unit_open(double x) noexcept { return x > 0.0 && x < 1.0; }

// The risks for an odds ratio orr >= 1: hi among the group with share s of
// the population, lo among the rest (share t = 1 - s). With c = OR - 1 and
// hi = OR lo / (1 + c lo), the rate equation s hi + t lo = q is the quadratic
//   t c lo^2 + (1 + s c - q c) lo - q = 0,
// negative at lo = 0 and positive at lo = 1, so it has one root in (0, 1).
inline void orrsolve(double orr, double q, double s, double t, double *hi,
                     double *lo) noexcept {
  const double c = orr - 1.0;
  double root;
  if (c <= 1.0) {
    const double qa = t * c;
    const double qb = 1.0 + s * c - q * c;
    if (qa == 0.0) {
      root = q / qb;  // OR = 1: both risks are the overall rate
    } else {
      const double sq = std::sqrt(qb * qb + 4.0 * qa * q);
      // the form without cancellation for the sign of qb
      root = qb >= 0.0 ? 2.0 * q / (qb + sq) : (sq - qb) / (2.0 * qa);
    }
  } else {
    // the same quadratic divided by c, so that a huge OR cannot overflow
    const double qb = 1.0 / c + s - q;
    const double qc = q / c;
    const double sq = std::sqrt(qb * qb + 4.0 * t * qc);
    root = qb >= 0.0 ? 2.0 * qc / (qb + sq) : (sq - qb) / (2.0 * t);
  }
  *lo = root;
  *hi = orr * root / (1.0 + c * root);
}

// log(1 + 1 / x) for x > 0 without forming 1 / x when it would overflow.
inline double l1pinv(double x) noexcept {
  return x < 1.0 ? std::log1p(x) - std::log(x) : std::log1p(1.0 / x);
}

}  // namespace

extern "C" {

const char *csls_strerror(int code) {
  switch (code) {
    case CSLS_OK:
      return "no error";
    case CSLS_EARG:
      return "a size or an index is out of range";
    case CSLS_EDOM:
      return "a value is outside the domain of the method";
    default:
      return "unknown error";
  }
}

int csls_vod(const double *day, const double *night, int n, int ref,
             const double *shift, double *share_day, double *share_night,
             double *ratio, int *sign, double *adj) {
  if (n < 1 || ref < 0 || ref >= n) return CSLS_EARG;
  double td = 0.0, tn = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!nonneg(day[i]) || !pos(night[i])) return CSLS_EDOM;
    if (shift != nullptr && !pos(shift[i])) return CSLS_EDOM;
    td += day[i];
    tn += night[i];
  }
  // the odds are taken against the reference group, so it must have been
  // seen in daylight
  if (!(day[ref] > 0.0)) return CSLS_EDOM;
  for (int i = 0; i < n; ++i) {
    share_day[i] = day[i] / td;
    share_night[i] = night[i] / tn;
    ratio[i] = (day[i] / day[ref]) / (night[i] / night[ref]);
    const double d = share_day[i] - share_night[i];
    sign[i] = std::fabs(d) < 1e-12 ? 0 : (d > 0.0 ? 1 : -1);
    if (shift != nullptr) adj[i] = ratio[i] / (shift[i] / shift[ref]);
  }
  return CSLS_OK;
}

int csls_orrisk(double orr, double q, double s, double *a, double *b) {
  if (!pos(orr) || !unit_open(q) || !unit_open(s)) return CSLS_EDOM;
  // An odds ratio below 1 is solved as its reciprocal with the two groups
  // exchanged: OR - 1 then never rounds to -1, which made 1 + (OR - 1) b
  // cancel as b neared 1.
  if (orr < 1.0) {
    const double inv = 1.0 / orr;
    if (!std::isfinite(inv)) return CSLS_EDOM;
    orrsolve(inv, q, 1.0 - s, s, b, a);
  } else {
    orrsolve(orr, q, s, 1.0 - s, a, b);
  }
  return CSLS_OK;
}

int csls_outtst(const double *q, const double *w, int n, double t,
                double *stop_rate, double *hit_rate) {
  if (n < 1) return CSLS_EARG;
  if (!std::isfinite(t)) return CSLS_EDOM;
  double tw = 0.0, sw = 0.0, sq = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!(q[i] >= 0.0 && q[i] <= 1.0) || !nonneg(w[i])) return CSLS_EDOM;
    tw += w[i];
    if (q[i] > t) {
      sw += w[i];
      sq += w[i] * q[i];
    }
  }
  if (!(tw > 0.0)) return CSLS_EDOM;
  *stop_rate = sw / tw;
  *hit_rate = sw > 0.0 ? sq / sw : std::numeric_limits<double>::quiet_NaN();
  return CSLS_OK;
}

int csls_dspbm(const double *pop, const double *contact, const double *force,
               int n, int ref, double *cd, double *fd, double *rd) {
  if (n < 1 || ref < 0 || ref >= n) return CSLS_EARG;
  for (int i = 0; i < n; ++i) {
    if (!pos(pop[i]) || !pos(contact[i]) || !pos(force[i])) return CSLS_EDOM;
  }
  const double cr = contact[ref] / pop[ref];
  const double fr = force[ref] / contact[ref];
  const double rr = force[ref] / pop[ref];
  for (int i = 0; i < n; ++i) {
    cd[i] = (contact[i] / pop[i]) / cr;
    fd[i] = (force[i] / contact[i]) / fr;
    rd[i] = (force[i] / pop[i]) / rr;
  }
  return CSLS_OK;
}

int csls_dyadrt(const double *count, const int *off, const int *vic, int m,
                const double *pop, int g, double *rate_off,
                double *rate_pair, double *null_off, double *ratio) {
  if (m < 1 || g < 1) return CSLS_EARG;
  double total = 0.0;
  for (int j = 0; j < g; ++j) {
    if (!pos(pop[j])) return CSLS_EDOM;
    total += pop[j];
  }
  double offences = 0.0, pairs = 0.0;
  for (int k = 0; k < m; ++k) {
    if (off[k] < 0 || off[k] >= g || vic[k] < 0 || vic[k] >= g) {
      return CSLS_EARG;
    }
    if (!nonneg(count[k])) return CSLS_EDOM;
    offences += count[k];
    // pairs exposed: N p_off p_vic, written without the division twice over
    pairs += pop[off[k]] * pop[vic[k]] / total;
  }
  const double pooled = offences / pairs;
  for (int k = 0; k < m; ++k) {
    const double pair = pop[off[k]] * pop[vic[k]] / total;
    rate_off[k] = count[k] / pop[off[k]];
    rate_pair[k] = count[k] / pair;
    null_off[k] = pooled * pop[vic[k]] / total;
    ratio[k] = rate_pair[k] / pooled;
  }
  return CSLS_OK;
}

int csls_gtmp(const double *p, const double *q, const int *hom, int n,
              double theta, double *locus, double *log10mp) {
  if (n < 1) return CSLS_EARG;
  if (!(theta >= 0.0 && theta < 1.0)) return CSLS_EDOM;
  const double den = (1.0 + theta) * (1.0 + 2.0 * theta);
  double sum = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!(p[i] > 0.0 && p[i] <= 1.0)) return CSLS_EDOM;
    const double a = theta + (1.0 - theta) * p[i];
    if (hom[i] != 0) {
      locus[i] = (theta + a) * (2.0 * theta + a) / den;
    } else {
      if (!(q[i] > 0.0 && p[i] + q[i] <= 1.0 + 1e-12)) return CSLS_EDOM;
      locus[i] = 2.0 * a * (theta + (1.0 - theta) * q[i]) / den;
    }
    sum += std::log10(locus[i]);
  }
  *log10mp = sum;
  return CSLS_OK;
}

int csls_cllr(const double *same, int ns, const double *diff, int nd,
              double *cllr) {
  if (ns < 1 || nd < 1) return CSLS_EARG;
  double cs = 0.0, cd = 0.0;
  for (int i = 0; i < ns; ++i) {
    if (!pos(same[i])) return CSLS_EDOM;
    cs += l1pinv(same[i]);
  }
  for (int i = 0; i < nd; ++i) {
    if (!pos(diff[i])) return CSLS_EDOM;
    cd += std::log1p(diff[i]);
  }
  *cllr = 0.5 * (cs / ns + cd / nd) / std::log(2.0);
  return CSLS_OK;
}

int csls_lrcal(const double *lr, const int *same, int n, double *post,
               double *work, double *cllrmin) {
  if (n < 2) return CSLS_EARG;
  double *sum = work;      // same-source count of each block
  double *cnt = work + n;  // size of each block
  int k = 0;
  double ns = 0.0;
  for (int i = 0; i < n;) {
    if (!pos(lr[i]) || (i > 0 && lr[i] < lr[i - 1])) return CSLS_EDOM;
    double s = 0.0, c = 0.0;
    int j = i;
    for (; j < n && lr[j] == lr[i]; ++j) {  // tied ratios share one block
      s += same[j] != 0 ? 1.0 : 0.0;
      c += 1.0;
    }
    ns += s;
    sum[k] = s;
    cnt[k] = c;
    ++k;
    // pool while the share of same-source falls from one block to the next
    while (k > 1 && sum[k - 2] * cnt[k - 1] >= sum[k - 1] * cnt[k - 2]) {
      sum[k - 2] += sum[k - 1];
      cnt[k - 2] += cnt[k - 1];
      --k;
    }
    i = j;
  }
  const double nd = n - ns;
  if (!(ns > 0.0) || !(nd > 0.0)) return CSLS_EDOM;
  int i = 0;
  for (int b = 0; b < k; ++b) {
    const double pr = sum[b] / cnt[b];
    for (int j = 0; j < static_cast<int>(cnt[b]); ++j) post[i++] = pr;
  }
  // a same-source item sits in a block with pr > 0 and a different-source
  // item in one with pr < 1, so neither division is by zero
  double cs = 0.0, cd = 0.0;
  for (i = 0; i < n; ++i) {
    const double pr = post[i];
    if (same[i] != 0) {
      cs += std::log1p((1.0 - pr) / pr * (ns / nd));
    } else {
      cd += std::log1p(pr / (1.0 - pr) * (nd / ns));
    }
  }
  *cllrmin = 0.5 * (cs / ns + cd / nd) / std::log(2.0);
  return CSLS_OK;
}

int csls_mixp(const double *p, const int *need, int n, int x, double *out) {
  if (n < 1 || x < 0) return CSLS_EARG;
  int idx[20];
  int m = 0;
  double total = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!(p[i] > 0.0 && p[i] <= 1.0)) return CSLS_EDOM;
    total += p[i];
    if (need[i] != 0) {
      if (m == 20) return CSLS_EARG;
      idx[m++] = i;
    }
  }
  if (total > 1.0 + 1e-12) return CSLS_EDOM;
  if (m > 2 * x) {  // more alleles to carry than the unknowns have
    *out = 0.0;
    return CSLS_OK;
  }
  // inclusion-exclusion over the needed alleles left out
  double acc = 0.0;
  for (unsigned long t = 0; t < (1UL << m); ++t) {
    double rest = total;
    int bits = 0;
    for (int j = 0; j < m; ++j) {
      if ((t >> j) & 1UL) {
        rest -= p[idx[j]];
        ++bits;
      }
    }
    if (rest < 0.0) rest = 0.0;
    const double term = std::pow(rest, 2.0 * x);
    acc += (bits % 2 == 0) ? term : -term;
  }
  *out = acc > 0.0 ? acc : 0.0;
  return CSLS_OK;
}

int csls_kinidx(const int *g1, const int *g2, const double *f, int n,
                const double *k, double *idx) {
  if (n < 1) return CSLS_EARG;
  if (!nonneg(k[0]) || !nonneg(k[1]) || !nonneg(k[2]) ||
      std::fabs(k[0] + k[1] + k[2] - 1.0) > 1e-9) {
    return CSLS_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    const int a = g1[2 * i], b = g1[2 * i + 1];
    const int c = g2[2 * i], d = g2[2 * i + 1];
    const double pa = f[2 * i], pb = f[2 * i + 1];
    if (!(pa > 0.0 && pa <= 1.0) || !(pb > 0.0 && pb <= 1.0)) return CSLS_EDOM;
    if (a == b && pa != pb) return CSLS_EDOM;
    // one allele shared by descent: half the weight on each allele of g1
    const double w = c == d ? 1.0 : 0.5;
    double r1 = 0.0;
    if (a == c || a == d) r1 += 0.5 * w / pa;
    if (b == c || b == d) r1 += 0.5 * w / pb;
    // both shared: the genotypes must be the same
    const bool eq = (a == c && b == d) || (a == d && b == c);
    const double r2 = !eq ? 0.0 : (a == b ? 1.0 / (pa * pa) : 0.5 / (pa * pb));
    idx[i] = k[0] + k[1] * r1 + k[2] * r2;
  }
  return CSLS_OK;
}

int csls_patidx(const int *c, const int *m, const int *af, const double *f,
                int n, double *idx) {
  if (n < 1) return CSLS_EARG;
  for (int i = 0; i < n; ++i) {
    const int c1 = c[2 * i], c2 = c[2 * i + 1];
    const double p1 = f[2 * i], p2 = f[2 * i + 1];
    if (!(p1 > 0.0 && p1 <= 1.0) || !(p2 > 0.0 && p2 <= 1.0)) return CSLS_EDOM;
    if (c1 == c2 && p1 != p2) return CSLS_EDOM;
    double num = 0.0, den = 0.0;
    for (int s = 0; s < 2; ++s) {
      const int g = m[2 * i + s];  // the allele the mother passes on
      if (c1 == c2) {
        if (g == c1) den += 0.5 * p1;
      } else {
        if (g == c1) den += 0.5 * p2;
        if (g == c2) den += 0.5 * p1;
      }
      for (int t = 0; t < 2; ++t) {
        const int h = af[2 * i + t];
        if ((g == c1 && h == c2) || (g == c2 && h == c1)) num += 0.25;
      }
    }
    idx[i] = den > 0.0 ? num / den : std::numeric_limits<double>::quiet_NaN();
  }
  return CSLS_OK;
}

}  // extern "C"
