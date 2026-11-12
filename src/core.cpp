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
  // With c = OR - 1 and a = OR b / (1 + c b), the rate equation
  // s a + (1 - s) b = q is the quadratic
  //   (1 - s) c b^2 + (1 + s c - q c) b - q = 0,
  // negative at b = 0 and positive at b = 1, so it has one root in (0, 1).
  const double c = orr - 1.0;
  const double qa = (1.0 - s) * c;
  const double qb = 1.0 + s * c - q * c;
  double root;
  if (qa == 0.0) {
    root = q / qb;  // OR = 1: both risks are the overall rate
  } else {
    const double sq = std::sqrt(qb * qb + 4.0 * qa * q);
    // the form without cancellation for the sign of qb
    root = qb >= 0.0 ? 2.0 * q / (qb + sq) : (sq - qb) / (2.0 * qa);
  }
  *b = root;
  *a = orr * root / (1.0 + c * root);
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

}  // extern "C"
