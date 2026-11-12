// SPDX-License-Identifier: AGPL-3.0-or-later
//
// libcrimlab, policing: aoristic weights (Ratcliffe 2002) on wall-clock
// bins, with daylight-saving transitions supplied by the caller.
#include "crimlab.h"
#include "kutil.h"

namespace {

struct bins {
  double origin, bw, per;
  int nb;
  // the bin of wall-clock time t
  int at(double t) const {
    double u = std::fmod(t - origin, per);
    if (u < 0.0) u += per;
    int j = static_cast<int>(u / bw);
    return j < nb ? j : nb - 1;  // u / bw can round up to nb
  }
  // spread mass dens per second uniformly over the wall-clock interval
  // [a, a + len) into out[j * stride], j = 0..nb-1; at most nb + 1 steps
  void spread(double a, double len, double dens, double *out, int stride) const {
    if (len >= per) {
      const double q = std::floor(len / per);
      for (int j = 0; j < nb; ++j) out[static_cast<long long>(j) * stride] += q * bw * dens;
      len -= q * per;
    }
    double u = std::fmod(a - origin, per);
    if (u < 0.0) u += per;
    int j = static_cast<int>(u / bw);
    if (j >= nb) j = nb - 1;
    double pos = u - j * bw;
    if (pos < 0.0) pos = 0.0;
    while (len > 0.0) {
      double take = bw - pos;
      if (take > len) take = len;
      out[static_cast<long long>(j) * stride] += take * dens;
      len -= take;
      pos = 0.0;
      j = j + 1 == nb ? 0 : j + 1;
    }
  }
};

// the UTC offset in force at real time t, and the index of the next
// transition after t (nt when none): binary search
int tzidx(const double *tt, int nt, double t) {
  int lo = 0, hi = nt;  // first index with tt[idx] > t
  while (lo < hi) {
    const int mid = lo + (hi - lo) / 2;
    if (tt[mid] <= t) lo = mid + 1; else hi = mid;
  }
  return lo;
}

}  // namespace

extern "C" {

int clab_aor(const double *from, const double *to, int n, const double *tt,
             const double *toff, int nt, double off0, double origin,
             double bw, int nb, double lo, double hi, double *w, double *tot,
             int *cnt) {
  if (n < 0 || nb < 1 || nt < 0) return CLAB_EARG;
  if (!kutil::pos(bw) || !std::isfinite(origin) || !std::isfinite(off0) ||
      !kutil::nonneg(lo) || !(hi > 0.0))
    return CLAB_EDOM;
  const bins b{origin, bw, bw * nb, nb};
  if (!std::isfinite(b.per)) return CLAB_EDOM;
  for (int k = 0; k < nt; ++k) {
    if (!std::isfinite(tt[k]) || !std::isfinite(toff[k])) return CLAB_EDOM;
    if (k > 0 && !(tt[k] > tt[k - 1])) return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    if (std::isinf(from[i]) || std::isinf(to[i])) return CLAB_EDOM;
  }
  for (int k = 0; k < 5; ++k) cnt[k] = 0;
  if (w != nullptr) {
    for (long long k = 0; k < static_cast<long long>(n) * nb; ++k) w[k] = 0.0;
  }
  if (tot != nullptr) {
    for (int j = 0; j < nb; ++j) tot[j] = 0.0;
  }
  for (int i = 0; i < n; ++i) {
    const double s = from[i], e = to[i];
    if (std::isnan(s)) {  // no start: the event is left out
      ++cnt[0];
      continue;
    }
    int k = tzidx(tt, nt, s);
    const double l = e - s;
    int rule = -1;  // -1: spread, 1..3: all in the start bin, 4: uniform
    if (std::isnan(e)) rule = 1;
    else if (l < 0.0) rule = 2;
    else if (l == 0.0 || l < lo) rule = 3;
    else if (l >= hi) rule = 4;
    if (rule >= 1) ++cnt[rule];
    if (rule >= 1 && rule <= 3) {
      const int j = b.at(s + (k == 0 ? off0 : toff[k - 1]));
      if (w != nullptr) w[i + static_cast<long long>(n) * j] += 1.0;
      if (tot != nullptr) tot[j] += 1.0;
      continue;
    }
    if (rule == 4) {
      for (int j = 0; j < nb; ++j) {
        if (w != nullptr) w[i + static_cast<long long>(n) * j] += 1.0 / nb;
        if (tot != nullptr) tot[j] += 1.0 / nb;
      }
      continue;
    }
    // walk the real interval [s, e) piece by piece between transitions;
    // each piece is a translate of wall-clock time
    const double dens = 1.0 / l;
    double a = s;
    while (a < e) {
      const double o = k == 0 ? off0 : toff[k - 1];
      const double pe = k < nt && tt[k] < e ? tt[k] : e;
      if (w != nullptr) b.spread(a + o, pe - a, dens, w + i, n);
      if (tot != nullptr) b.spread(a + o, pe - a, dens, tot, 1);
      a = pe;
      ++k;
    }
  }
  return CLAB_OK;
}

}  // extern "C"
