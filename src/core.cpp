// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab: the kernels. No R here, no allocation, no exceptions: every
// function reads the caller's arrays and writes the caller's arrays.
#include "crimlab.h"

#include <cmath>
#include <limits>

namespace {

//' @srrstats {G3.0} Computed doubles are compared with tolerances (1e-12
//'   for the sign in clab_vod(), 1e-9 for probabilities summing to 1). The
//'   exact comparisons are deliberate: whether a number is whole, ties
//'   between likelihood ratios in clab_lrcal(), the two frequencies given
//'   for a homozygote's one allele, outcomes of 0 and 1, and the guard
//'   against an exact zero divisor in orrsolve().
//' @noRd
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

const char *clab_strerror(int code) {
  switch (code) {
    case CLAB_OK:
      return "no error";
    case CLAB_EARG:
      return "a size or an index is out of range";
    case CLAB_EDOM:
      return "a value is outside the domain of the method";
    default:
      return "unknown error";
  }
}

int clab_vod(const double *day, const double *night, int n, int ref,
             const double *shift, double *share_day, double *share_night,
             double *ratio, int *sign, double *adj) {
  if (n < 1 || ref < 0 || ref >= n) return CLAB_EARG;
  double td = 0.0, tn = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!nonneg(day[i]) || !pos(night[i])) return CLAB_EDOM;
    if (shift != nullptr && !pos(shift[i])) return CLAB_EDOM;
    td += day[i];
    tn += night[i];
  }
  // the odds are taken against the reference group, so it must have been
  // seen in daylight
  if (!(day[ref] > 0.0)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    share_day[i] = day[i] / td;
    share_night[i] = night[i] / tn;
    ratio[i] = (day[i] / day[ref]) / (night[i] / night[ref]);
    const double d = share_day[i] - share_night[i];
    sign[i] = std::fabs(d) < 1e-12 ? 0 : (d > 0.0 ? 1 : -1);
    if (shift != nullptr) adj[i] = ratio[i] / (shift[i] / shift[ref]);
  }
  return CLAB_OK;
}

int clab_orrisk(double orr, double q, double s, double *a, double *b) {
  if (!pos(orr) || !unit_open(q) || !unit_open(s)) return CLAB_EDOM;
  // An odds ratio below 1 is solved as its reciprocal with the two groups
  // exchanged: OR - 1 then never rounds to -1, which made 1 + (OR - 1) b
  // cancel as b neared 1.
  if (orr < 1.0) {
    const double inv = 1.0 / orr;
    if (!std::isfinite(inv)) return CLAB_EDOM;
    orrsolve(inv, q, 1.0 - s, s, b, a);
  } else {
    orrsolve(orr, q, s, 1.0 - s, a, b);
  }
  return CLAB_OK;
}

int clab_outtst(const double *q, const double *w, int n, double t,
                double *stop_rate, double *hit_rate) {
  if (n < 1) return CLAB_EARG;
  if (!std::isfinite(t)) return CLAB_EDOM;
  double tw = 0.0, sw = 0.0, sq = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!(q[i] >= 0.0 && q[i] <= 1.0) || !nonneg(w[i])) return CLAB_EDOM;
    tw += w[i];
    if (q[i] > t) {
      sw += w[i];
      sq += w[i] * q[i];
    }
  }
  if (!(tw > 0.0)) return CLAB_EDOM;
  *stop_rate = sw / tw;
  *hit_rate = sw > 0.0 ? sq / sw : std::numeric_limits<double>::quiet_NaN();
  return CLAB_OK;
}

int clab_dspbm(const double *pop, const double *contact, const double *force,
               int n, int ref, double *cd, double *fd, double *rd) {
  if (n < 1 || ref < 0 || ref >= n) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!pos(pop[i]) || !pos(contact[i]) || !pos(force[i])) return CLAB_EDOM;
  }
  const double cr = contact[ref] / pop[ref];
  const double fr = force[ref] / contact[ref];
  const double rr = force[ref] / pop[ref];
  for (int i = 0; i < n; ++i) {
    cd[i] = (contact[i] / pop[i]) / cr;
    fd[i] = (force[i] / contact[i]) / fr;
    rd[i] = (force[i] / pop[i]) / rr;
  }
  return CLAB_OK;
}

int clab_dyadrt(const double *count, const int *off, const int *vic, int m,
                const double *pop, int g, double *rate_off,
                double *rate_pair, double *null_off, double *ratio) {
  if (m < 1 || g < 1) return CLAB_EARG;
  double total = 0.0;
  for (int j = 0; j < g; ++j) {
    if (!pos(pop[j])) return CLAB_EDOM;
    total += pop[j];
  }
  double offences = 0.0, pairs = 0.0;
  for (int k = 0; k < m; ++k) {
    if (off[k] < 0 || off[k] >= g || vic[k] < 0 || vic[k] >= g) {
      return CLAB_EARG;
    }
    if (!nonneg(count[k])) return CLAB_EDOM;
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
  return CLAB_OK;
}

int clab_gtmp(const double *p, const double *q, const int *hom, int n,
              double theta, double *locus, double *log10mp) {
  if (n < 1) return CLAB_EARG;
  if (!(theta >= 0.0 && theta < 1.0)) return CLAB_EDOM;
  const double den = (1.0 + theta) * (1.0 + 2.0 * theta);
  double sum = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!(p[i] > 0.0 && p[i] <= 1.0)) return CLAB_EDOM;
    const double a = theta + (1.0 - theta) * p[i];
    if (hom[i] != 0) {
      locus[i] = (theta + a) * (2.0 * theta + a) / den;
    } else {
      if (!(q[i] > 0.0 && p[i] + q[i] <= 1.0 + 1e-12)) return CLAB_EDOM;
      locus[i] = 2.0 * a * (theta + (1.0 - theta) * q[i]) / den;
    }
    sum += std::log10(locus[i]);
  }
  *log10mp = sum;
  return CLAB_OK;
}

int clab_cllr(const double *same, int ns, const double *diff, int nd,
              double *cllr) {
  if (ns < 1 || nd < 1) return CLAB_EARG;
  double cs = 0.0, cd = 0.0;
  for (int i = 0; i < ns; ++i) {
    if (!pos(same[i])) return CLAB_EDOM;
    cs += l1pinv(same[i]);
  }
  for (int i = 0; i < nd; ++i) {
    if (!pos(diff[i])) return CLAB_EDOM;
    cd += std::log1p(diff[i]);
  }
  *cllr = 0.5 * (cs / ns + cd / nd) / std::log(2.0);
  return CLAB_OK;
}

int clab_lrcal(const double *lr, const int *same, int n, double *post,
               double *work, double *cllrmin) {
  if (n < 2) return CLAB_EARG;
  double *sum = work;      // same-source count of each block
  double *cnt = work + n;  // size of each block
  int k = 0;
  double ns = 0.0;
  for (int i = 0; i < n;) {
    if (!pos(lr[i]) || (i > 0 && lr[i] < lr[i - 1])) return CLAB_EDOM;
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
  if (!(ns > 0.0) || !(nd > 0.0)) return CLAB_EDOM;
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
  return CLAB_OK;
}

int clab_mixp(const double *p, const int *need, int n, int x, double *out) {
  if (n < 1 || x < 0) return CLAB_EARG;
  int idx[20];
  int m = 0;
  double total = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!(p[i] > 0.0 && p[i] <= 1.0)) return CLAB_EDOM;
    total += p[i];
    if (need[i] != 0) {
      if (m == 20) return CLAB_EARG;
      idx[m++] = i;
    }
  }
  if (total > 1.0 + 1e-12) return CLAB_EDOM;
  if (m > 2 * x) {  // more alleles to carry than the unknowns have
    *out = 0.0;
    return CLAB_OK;
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
  return CLAB_OK;
}

int clab_kinidx(const int *g1, const int *g2, const double *f, int n,
                const double *k, double *idx) {
  if (n < 1) return CLAB_EARG;
  if (!nonneg(k[0]) || !nonneg(k[1]) || !nonneg(k[2]) ||
      std::fabs(k[0] + k[1] + k[2] - 1.0) > 1e-9) {
    return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    const int a = g1[2 * i], b = g1[2 * i + 1];
    const int c = g2[2 * i], d = g2[2 * i + 1];
    const double pa = f[2 * i], pb = f[2 * i + 1];
    if (!(pa > 0.0 && pa <= 1.0) || !(pb > 0.0 && pb <= 1.0)) return CLAB_EDOM;
    if (a == b && pa != pb) return CLAB_EDOM;
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
  return CLAB_OK;
}

int clab_patidx(const int *c, const int *m, const int *af, const double *f,
                int n, double *idx) {
  if (n < 1) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    const int c1 = c[2 * i], c2 = c[2 * i + 1];
    const double p1 = f[2 * i], p2 = f[2 * i + 1];
    if (!(p1 > 0.0 && p1 <= 1.0) || !(p2 > 0.0 && p2 <= 1.0)) return CLAB_EDOM;
    if (c1 == c2 && p1 != p2) return CLAB_EDOM;
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
  return CLAB_OK;
}

}  // extern "C"

// Policing and corrections.

extern "C" {

int clab_clrsmp(const double *clra, const double *reca, const double *clrb,
                const double *recb, int n, double *ra, double *rb,
                double *wa, double *wb, double *reff, double *meff,
                double *tot, int *cnt) {
  if (n < 1) return CLAB_EARG;
  double sca = 0.0, sra = 0.0, scb = 0.0, srb = 0.0;
  for (int k = 0; k < n; ++k) {
    if (!nonneg(clra[k]) || !pos(reca[k]) || !nonneg(clrb[k]) ||
        !pos(recb[k])) {
      return CLAB_EDOM;
    }
    sca += clra[k];
    sra += reca[k];
    scb += clrb[k];
    srb += recb[k];
  }
  const double agga = sca / sra;
  const double aggb = scb / srb;
  // the centre of the mix terms; their sum does not depend on it
  // (Crimlab.Clearance.kitagawa), only how it is shared out among the parts
  const double centre = 0.5 * (agga + aggb);
  double sre = 0.0, sme = 0.0;
  int up = 0, down = 0;
  for (int k = 0; k < n; ++k) {
    ra[k] = clra[k] / reca[k];
    rb[k] = clrb[k] / recb[k];
    wa[k] = reca[k] / sra;
    wb[k] = recb[k] / srb;
    reff[k] = 0.5 * (wa[k] + wb[k]) * (rb[k] - ra[k]);
    meff[k] = (0.5 * (ra[k] + rb[k]) - centre) * (wb[k] - wa[k]);
    sre += reff[k];
    sme += meff[k];
    if (rb[k] > ra[k]) ++up;
    if (rb[k] < ra[k]) ++down;
  }
  const double change = aggb - agga;
  tot[0] = agga;
  tot[1] = aggb;
  tot[2] = change;
  tot[3] = sre;
  tot[4] = sme;
  int type = 0;
  if ((up == n && change < 0.0) || (down == n && change > 0.0)) {
    type = 2;  // every part one way, the aggregate the other
  } else if (change != 0.0 && sre != 0.0 && (change > 0.0) != (sre > 0.0)) {
    type = 1;  // the aggregate moves against the parts' own rates
  }
  cnt[0] = up;
  cnt[1] = down;
  cnt[2] = type;
  return CLAB_OK;
}

int clab_expbnd(const double *y, const double *m, int n, int ref, double gamma,
                double *rate, double *lo, double *hi, double *ratio,
                double *rlo, double *rhi, int *sign) {
  if (n < 1 || ref < 0 || ref >= n) return CLAB_EARG;
  if (!(std::isfinite(gamma) && gamma >= 1.0)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!pos(y[i]) || !pos(m[i])) return CLAB_EDOM;
  }
  const double g2 = gamma * gamma;
  const double base = y[ref] / m[ref];
  for (int i = 0; i < n; ++i) {
    rate[i] = y[i] / m[i];
    lo[i] = rate[i] / gamma;
    hi[i] = gamma * rate[i];
    ratio[i] = rate[i] / base;
    rlo[i] = ratio[i] / g2;
    rhi[i] = g2 * ratio[i];
    sign[i] = rlo[i] > 1.0 ? 1 : (rhi[i] < 1.0 ? -1 : 0);
  }
  return CLAB_OK;
}

int clab_misrd(double r1, double r0, double p, double se, double sp,
               double *out) {
  if (!(r1 >= 0.0 && r1 <= 1.0) || !(r0 >= 0.0 && r0 <= 1.0) ||
      !unit_open(p) || !(se >= 0.0 && se <= 1.0) ||
      !(sp >= 0.0 && sp <= 1.0)) {
    return CLAB_EDOM;
  }
  // shares of the population in the measured groups
  const double d1 = p * se + (1.0 - p) * (1.0 - sp);
  const double d0 = p * (1.0 - se) + (1.0 - p) * sp;
  if (!(d1 > 0.0) || !(d0 > 0.0)) return CLAB_EDOM;
  // each measured group is a mixture of the true groups
  const double m1 = (p * se * r1 + (1.0 - p) * (1.0 - sp) * r0) / d1;
  const double m0 = (p * (1.0 - se) * r1 + (1.0 - p) * sp * r0) / d0;
  out[0] = r1 - r0;
  out[1] = m1 - m0;
  out[2] = (se + sp - 1.0) * p * (1.0 - p) / (d1 * d0);
  out[3] = d1;
  out[4] = d0;
  return CLAB_OK;
}

int clab_rtmhs(const double *x1, const double *x2, const double *w, int n,
               double c, double *out) {
  if (n < 1) return CLAB_EARG;
  if (!std::isfinite(c)) return CLAB_EDOM;
  double s1 = 0.0, m1 = 0.0, s2 = 0.0, m2 = 0.0;  // high: by x1, by x2
  double l1 = 0.0, n1 = 0.0, l2 = 0.0, n2 = 0.0;  // low: by x1, by x2
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x1[i]) || !std::isfinite(x2[i]) || !nonneg(w[i])) {
      return CLAB_EDOM;
    }
    const double d = x2[i] - x1[i];
    if (x1[i] > c) {
      s1 += w[i] * d;
      m1 += w[i];
    }
    if (x2[i] > c) {
      s2 -= w[i] * d;
      m2 += w[i];
    }
    if (x1[i] < c) {
      l1 += w[i] * d;
      n1 += w[i];
    }
    if (x2[i] < c) {
      l2 -= w[i] * d;
      n2 += w[i];
    }
  }
  if (!(m1 > 0.0)) return CLAB_EDOM;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  out[0] = m1;
  out[1] = s1 / m1;
  out[2] = m2 > 0.0 ? s2 / m2 : nan;
  // selection on x1 in the data together with its period-swapped copy
  out[3] = (s1 + s2) / (m1 + m2);
  out[4] = out[1] - out[3];
  out[5] = n1 > 0.0 ? l1 / n1 : nan;
  out[6] = n1 + n2 > 0.0 ? (l1 + l2) / (n1 + n2) : nan;
  return CLAB_OK;
}

int clab_incap(const double *lam, const double *q, const double *s, int n,
               const double *share, double *rate, double *prev, double *marg,
               double *agg) {
  if (n < 1) return CLAB_EARG;
  double tot = 0.0, free = 0.0, inc = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!pos(lam[i]) || !(q[i] > 0.0 && q[i] <= 1.0) || !nonneg(s[i])) {
      return CLAB_EDOM;
    }
    if (share != nullptr && !nonneg(share[i])) return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    const double k = lam[i] * q[i];
    const double d = 1.0 + k * s[i];
    rate[i] = lam[i] / d;
    prev[i] = k * s[i] / d;
    marg[i] = lam[i] * k / (d * (d + k));
    if (share != nullptr) {
      tot += share[i];
      free += share[i] * lam[i];
      inc += share[i] * rate[i];
    }
  }
  if (share != nullptr) {
    if (std::fabs(tot - 1.0) > 1e-9) return CLAB_EDOM;
    agg[0] = free;
    agg[1] = inc;
    agg[2] = 1.0 - inc / free;
  }
  return CLAB_OK;
}

int clab_incapr(const double *lam, int len, int t0, int s, double r,
                double *out) {
  if (t0 < 0 || s < 1 || len < 1 || t0 > len - 1 - s) return CLAB_EARG;
  if (!(r >= 0.0 && r <= 1.0)) return CLAB_EDOM;
  for (int i = 0; i < len; ++i) {
    if (!nonneg(lam[i]) || (i > 0 && lam[i] > lam[i - 1])) return CLAB_EDOM;
  }
  double now = 0.0, later = 0.0;
  for (int j = 0; j < s; ++j) {
    now += lam[t0 + j];
    later += lam[t0 + 1 + j];  // t0 + s <= len - 1, so always in range
  }
  out[0] = now;
  out[1] = (1.0 - r) * now;
  out[2] = s * lam[t0];
  out[3] = s * lam[t0 + s];
  out[4] = later;
  return CLAB_OK;
}

int clab_sellab(const double *y, const int *rel, const int *rule,
                const double *w, int n, double *out) {
  if (n < 1) return CLAB_EARG;
  double mr = 0.0, fr = 0.0, mm = 0.0, fi = 0.0, mi = 0.0, mu = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!nonneg(w[i])) return CLAB_EDOM;
    if (rel[i] != 0) {
      if (!(y[i] == 0.0 || y[i] == 1.0)) return CLAB_EDOM;
      mr += w[i];
      fr += w[i] * y[i];
    }
    if (rule[i] != 0) {
      mm += w[i];
      if (rel[i] != 0) {
        mi += w[i];
        fi += w[i] * y[i];
      } else {
        mu += w[i];  // released by the rule, detained by the judge
      }
    }
  }
  if (!(mr > 0.0) || !(mm > 0.0)) return CLAB_EDOM;
  out[0] = fr / mr;
  out[1] = mu / mm;
  out[2] = fi / mm;
  out[3] = (fi + mu) / mm;
  out[4] = out[3] - out[2];
  out[5] = mi > 0.0 ? fi / mi : std::numeric_limits<double>::quiet_NaN();
  return CLAB_OK;
}

int clab_losbnd(const double *t, int n, const double *a, int m, double *out) {
  if (n < 1 || m < 1) return CLAB_EARG;
  double st = 0.0, sa = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!nonneg(t[i])) return CLAB_EDOM;
    st += t[i];
  }
  for (int j = 0; j < m; ++j) {
    if (!nonneg(a[j])) return CLAB_EDOM;
    sa += a[j];
  }
  const double tbar = st / n;
  const double abar = sa / m;
  out[0] = tbar;
  out[1] = abar;
  out[2] = (st + sa) / (static_cast<double>(n) + m);
  out[3] = static_cast<double>(m) / (static_cast<double>(n) + m) *
           (abar - tbar);
  return CLAB_OK;
}

int clab_sntbnd(const double *y, const int *z, const double *w, int n,
                double *out) {
  if (n < 1) return CLAB_EARG;
  double tw = 0.0, j[2] = {0.0, 0.0}, p[2] = {0.0, 0.0};
  for (int i = 0; i < n; ++i) {
    if (!(y[i] == 0.0 || y[i] == 1.0) || !nonneg(w[i]) ||
        (z[i] != 0 && z[i] != 1)) {
      return CLAB_EDOM;
    }
    tw += w[i];
    j[z[i]] += w[i] * y[i];
    p[z[i]] += w[i];
  }
  if (!(tw > 0.0)) return CLAB_EDOM;
  for (int t = 0; t < 2; ++t) {
    j[t] /= tw;
    p[t] /= tw;
  }
  // P(y(t) = 1) lies in [P(y = 1, z = t), P(y = 1, z = t) + P(z != t)]
  const double lo0 = j[0], hi0 = j[0] + p[1];
  const double lo1 = j[1], hi1 = j[1] + p[0];
  const double nan = std::numeric_limits<double>::quiet_NaN();
  out[0] = j[0];
  out[1] = j[1];
  out[2] = p[0];
  out[3] = p[1];
  out[4] = lo0;
  out[5] = hi0;
  out[6] = lo1;
  out[7] = hi1;
  out[8] = lo1 - hi0;
  out[9] = hi1 - lo0;
  out[10] = p[0] > 0.0 && p[1] > 0.0 ? j[1] / p[1] - j[0] / p[0] : nan;
  return CLAB_OK;
}

int clab_hzdsel(double s, double h, double l, const double *ah,
                const double *al, double *out) {
  if (!unit_open(s) || !(l >= 0.0 && l < h && h <= 1.0)) return CLAB_EDOM;
  for (int k = 0; k < 2; ++k) {
    if (!(ah[k] > 0.0 && ah[k] <= 1.0) || !(al[k] > 0.0 && al[k] <= 1.0)) {
      return CLAB_EDOM;
    }
    // the high-risk share among the survivors of period 1, and their
    // period-2 hazard
    out[k] = s * ah[k] / (s * ah[k] + (1.0 - s) * al[k]);
    out[2 + k] = out[k] * h + (1.0 - out[k]) * l;
  }
  out[4] = out[3] / out[2];
  return CLAB_OK;
}

}  // extern "C"
