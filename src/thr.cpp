// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, tests of discrimination in searches: the hit-rate (outcome)
// test of Knowles, Persico and Todd (2001), the threshold test with the
// discriminant distributions of Pierson, Corbett-Davies and Goel (2018)
// fitted by maximum likelihood, and the rank-order test of Anwar and Fang
// (2006).
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "qn.h"

namespace {

constexpr double rt2 = 1.41421356237309504880;
constexpr double ir2pi = 0.39894228040143267794;

double qup(double z) { return 0.5 * std::erfc(z / rt2); }
double dn(double z) { return ir2pi * std::exp(-0.5 * z * z); }
double lgs(double x) { return 1.0 / (1.0 + std::exp(-x)); }

// One cell of the model: its counts and the positions of its six
// parameters in u (-1: fixed at 0, the first area's offsets).
struct cell {
  double n, k, m;
  int ia, ja, ib, jb, ic, jc;
};

struct mod {
  const std::vector<cell> *c;
  double lc;  // the binomial coefficients' part of the log-likelihood
};

double at(const double *u, int i) { return i < 0 ? 0.0 : u[i]; }

// The search and hit probabilities of a cell at its natural parameters
// (prior phi, separation dl, threshold t): the signal is N(dl, 1) for the
// guilty and N(0, 1) for the rest, the risk its posterior, and a search
// one whose risk exceeds t. A = phi Q(x - dl), B = (1 - phi) Q(x).
void ab(double al, double dl, double ga, double &x, double &A, double &B) {
  const double ph = lgs(al);
  x = (ga - al) / dl + 0.5 * dl;
  A = ph * qup(x - dl);
  B = (1.0 - ph) * qup(x);
}

// The negative log-likelihood and its gradient over u.
double nll(const double *u, double *g, void *ctx) {
  const mod &M = *static_cast<const mod *>(ctx);
  const std::vector<cell> &cs = *M.c;
  int np = 0;
  for (const cell &e : cs) {
    const int mx[6] = {e.ia, e.ja, e.ib, e.jb, e.ic, e.jc};
    for (int q : mx) np = q + 1 > np ? q + 1 : np;
  }
  if (g != nullptr) {
    for (int q = 0; q < np; ++q) g[q] = 0.0;
  }
  kutil::ksum s;
  for (const cell &e : cs) {
    const double al = at(u, e.ia) + at(u, e.ja);
    const double be = at(u, e.ib) + at(u, e.jb);
    const double ga = at(u, e.ic) + at(u, e.jc);
    const double dl = std::exp(be);
    double x = 0.0, A = 0.0, B = 0.0;
    ab(al, dl, ga, x, A, B);
    const double r = 1.0 - A - B;
    const double km = e.k - e.m, nk = e.n - e.k;
    if ((e.m > 0 && !(A > 0)) || (km > 0 && !(B > 0)) || (nk > 0 && !(r > 0)) || !std::isfinite(dl)) {
      return HUGE_VAL;
    }
    double l = 0.0;
    if (e.m > 0) l += e.m * std::log(A);
    if (km > 0) l += km * std::log(B);
    if (nk > 0) l += nk * std::log(r);
    s.add(l);
    if (g == nullptr) continue;
    const double ph = lgs(al), f1 = dn(x - dl), f0 = dn(x);
    const double q1 = qup(x - dl), q0 = qup(x);
    // derivatives of x, then of A and B, by al, ga and dl
    const double xa = -1.0 / dl, xg = 1.0 / dl, xd = 0.5 - (ga - al) / (dl * dl);
    const double Aa = ph * (1 - ph) * q1 - ph * f1 * xa, Ba = -ph * (1 - ph) * q0 - (1 - ph) * f0 * xa;
    const double Ag = -ph * f1 * xg, Bg = -(1 - ph) * f0 * xg;
    const double Ad = -ph * f1 * (xd - 1.0), Bd = -(1 - ph) * f0 * xd;
    auto d = [&](double dA, double dB) {
      double v = 0.0;
      if (e.m > 0) v += e.m * dA / A;
      if (km > 0) v += km * dB / B;
      if (nk > 0) v -= nk * (dA + dB) / r;
      return -v;
    };
    const double da = d(Aa, Ba), dg = d(Ag, Bg), db = d(Ad, Bd) * dl;
    if (e.ia >= 0) g[e.ia] += da;
    if (e.ja >= 0) g[e.ja] += da;
    if (e.ib >= 0) g[e.ib] += db;
    if (e.jb >= 0) g[e.jb] += db;
    if (e.ic >= 0) g[e.ic] += dg;
    if (e.jc >= 0) g[e.jc] += dg;
  }
  const double v = -(s.val() + M.lc);
  return std::isfinite(v) ? v : HUGE_VAL;
}

// lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
double lch(double n, double k) { return std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1); }

}  // namespace

int clab_thrtst(const double *stops, const double *srch, const double *hits, const int *g,
                const int *a, int n, int ng, int na, int ref, int maxit, double *gs, double *hz,
                double *chi, double *thr, double *tse, double *cfit, double *lrt, int *conv) {
  if (n < 1 || ng < 2 || na < 1 || ref < 0 || ref >= ng || maxit < 0) return CLAB_EARG;
  if (static_cast<long long>(ng) * na > 50000000LL) return CLAB_EARG;
  const int nc = ng * na;
  std::vector<double> S(static_cast<std::size_t>(3 * nc), 0.0);
  for (int i = 0; i < n; ++i) {
    if (g[i] < 0 || g[i] >= ng || a[i] < 0 || a[i] >= na) return CLAB_EARG;
    const double s0 = stops[i], s1 = srch[i], s2 = hits[i];
    if (!kutil::nonneg(s0) || !kutil::nonneg(s1) || !kutil::nonneg(s2) || s1 > s0 || s2 > s1) {
      return CLAB_EDOM;
    }
    const std::size_t c = static_cast<std::size_t>(a[i] * ng + g[i]);
    S[3 * c] += s0;
    S[3 * c + 1] += s1;
    S[3 * c + 2] += s2;
  }
  // hit and search rates by group (gs: ng x 6 column-major: stops,
  // searches, hits, search rate, hit rate, total), z against ref
  for (int j = 0; j < ng; ++j) {
    double t0 = 0, t1 = 0, t2 = 0;
    for (int q = 0; q < na; ++q) {
      const std::size_t c = static_cast<std::size_t>(q * ng + j);
      t0 += S[3 * c];
      t1 += S[3 * c + 1];
      t2 += S[3 * c + 2];
    }
    gs[j] = t0;
    gs[ng + j] = t1;
    gs[2 * ng + j] = t2;
    gs[3 * ng + j] = t0 > 0 ? t1 / t0 : kutil::nan();
    gs[4 * ng + j] = t1 > 0 ? t2 / t1 : kutil::nan();
  }
  double K = 0, H = 0;
  int used = 0;
  for (int j = 0; j < ng; ++j) {
    K += gs[ng + j];
    H += gs[2 * ng + j];
    used += gs[ng + j] > 0;
  }
  if (used < ng) return CLAB_EDOM;
  // hz: ng x 4: difference of hit rates from ref, its standard error, z, p
  const double h0 = gs[4 * ng + ref], k0 = gs[ng + ref];
  for (int j = 0; j < ng; ++j) {
    const double h = gs[4 * ng + j], k = gs[ng + j];
    const double d = h - h0;
    const double se = std::sqrt(h * (1 - h) / k + h0 * (1 - h0) / k0);
    hz[j] = d;
    hz[ng + j] = j == ref ? kutil::nan() : se;
    hz[2 * ng + j] = j == ref || !(se > 0) ? kutil::nan() : d / se;
    hz[3 * ng + j] = std::isnan(hz[2 * ng + j]) ? kutil::nan() : 2 * qup(std::fabs(hz[2 * ng + j]));
  }
  // Pearson's chi-square of hits and misses among the searched by group
  {
    const double hp = H / K;
    kutil::ksum x2;
    for (int j = 0; j < ng; ++j) {
      const double k = gs[ng + j];
      if (!(k > 0)) continue;
      const double e1 = k * hp, e0 = k * (1 - hp);
      const double o1 = gs[2 * ng + j], o0 = k - o1;
      if (e1 > 0) x2.add((o1 - e1) * (o1 - e1) / e1);
      if (e0 > 0) x2.add((o0 - e0) * (o0 - e0) / e0);
    }
    chi[0] = x2.val();
    chi[1] = used - 1;
    chi[2] = hp > 0 && hp < 1 ? clab_pchisq(chi[0], chi[1], 0) : kutil::nan();
  }
  if (maxit == 0) return CLAB_OK;
  // the threshold model on the cells with stops
  std::vector<cell> cs;
  double lc = 0.0;
  for (int q = 0; q < na; ++q) {
    for (int j = 0; j < ng; ++j) {
      const std::size_t c = static_cast<std::size_t>(q * ng + j);
      if (!(S[3 * c] > 0)) continue;
      cell e{S[3 * c], S[3 * c + 1], S[3 * c + 2], 0, 0, 0, 0, 0, 0};
      cs.push_back(e);
      lc += lch(S[3 * c], S[3 * c + 1]) + lch(S[3 * c + 1], S[3 * c + 2]);
    }
  }
  const int nfull = 3 * (ng + na - 1);
  if (2 * static_cast<long long>(cs.size()) <= nfull) return CLAB_EDOM;
  std::vector<double> wa(static_cast<std::size_t>(na), 0.0);
  double tot = 0.0;
  for (int q = 0; q < na; ++q) {
    for (int j = 0; j < ng; ++j) wa[static_cast<std::size_t>(q)] += S[3 * static_cast<std::size_t>(q * ng + j)];
    tot += wa[static_cast<std::size_t>(q)];
  }
  for (double &w : wa) w /= tot;
  std::vector<double> best(static_cast<std::size_t>(nfull), 0.0), cov(static_cast<std::size_t>(nfull) * nfull);
  std::vector<double> same(static_cast<std::size_t>(nfull - ng + 1), 0.0);
  double ll[2] = {0, 0};
  for (int r = 1; r >= 0; --r) {
    // r 1: one threshold for all groups; r 0: a threshold per group, also
    // started from the fit of r 1
    const int cg = r == 0 ? ng : 1, np = 3 * (ng + na - 1) - (ng - cg);
    std::size_t z = 0;
    for (int q = 0; q < na; ++q) {
      for (int j = 0; j < ng; ++j) {
        if (!(S[3 * static_cast<std::size_t>(q * ng + j)] > 0)) continue;
        cell &e = cs[z++];
        e.ia = j;
        e.ja = q > 0 ? ng + q - 1 : -1;
        e.ib = ng + na - 1 + j;
        e.jb = q > 0 ? 2 * ng + na - 1 + q - 1 : -1;
        e.ic = 2 * (ng + na - 1) + (r == 0 ? j : 0);
        e.jc = q > 0 ? 2 * (ng + na - 1) + cg + q - 1 : -1;
      }
    }
    mod M{&cs, lc};
    std::vector<double> u(static_cast<std::size_t>(np)), w(static_cast<std::size_t>(np) * np + 5 * np);
    double fbest = HUGE_VAL;
    int cbest = 0, ibest = 0;
    // three starts, separations 1, 2 and 3, and for r 0 the fit of r 1
    for (int st = 1; st <= (r == 0 ? 4 : 3); ++st) {
      for (double &v : u) v = 0.0;
      if (st == 4) {
        const int h = 2 * (ng + na - 1);
        for (int q = 0; q < h; ++q) u[static_cast<std::size_t>(q)] = same[static_cast<std::size_t>(q)];
        for (int j = 0; j < ng; ++j) u[static_cast<std::size_t>(h + j)] = same[static_cast<std::size_t>(h)];
        for (int q = 1; q < na; ++q) {
          u[static_cast<std::size_t>(h + ng + q - 1)] = same[static_cast<std::size_t>(h + q)];
        }
      } else {
        const double hp = H / K;
        for (int j = 0; j < ng; ++j) {
          const double ph = std::fmin(0.9, std::fmax(0.01, gs[2 * ng + j] / std::fmax(gs[j], 1) + 0.01));
          u[static_cast<std::size_t>(j)] = std::log(ph / (1 - ph));
          u[static_cast<std::size_t>(ng + na - 1 + j)] = std::log(static_cast<double>(st));
        }
        const double t0 = std::fmin(0.9, std::fmax(0.01, 0.5 * hp));
        for (int j = 0; j < cg; ++j) u[static_cast<std::size_t>(2 * (ng + na - 1) + j)] = std::log(t0 / (1 - t0));
      }
      double fv = 0.0;
      int it = 0;
      const int c = qn::min(np, u.data(), nll, &M, maxit, 1e-7, &fv, &it, w.data(),
                            static_cast<double>(cs.size()) * 40.0);
      if (kutil::stop(0.0)) return CLAB_EINTR;
      if (c >= 0 && fv < fbest) {
        fbest = fv;
        cbest = c;
        ibest = it;
        for (int q = 0; q < np; ++q) best[static_cast<std::size_t>(q)] = u[static_cast<std::size_t>(q)];
      }
    }
    if (!std::isfinite(fbest)) return CLAB_EDOM;
    ll[r] = -fbest;
    conv[r] = cbest;
    conv[2 + r] = ibest;
    if (r == 1) {
      for (int q = 0; q < np; ++q) same[static_cast<std::size_t>(q)] = best[static_cast<std::size_t>(q)];
      continue;
    }
    // the full model: thresholds, their averages and the fitted rates
    std::vector<double> hw(static_cast<std::size_t>(3 * np));
    qn::hcov(np, best.data(), nll, &M, cov.data(), hw.data());
    std::vector<double> G(static_cast<std::size_t>(ng) * np, 0.0);
    for (int j = 0; j < ng; ++j) {
      double tb = 0.0;
      for (int q = 0; q < na; ++q) {
        const double cc = best[static_cast<std::size_t>(2 * (ng + na - 1) + j)] +
                          (q > 0 ? best[static_cast<std::size_t>(2 * (ng + na - 1) + ng + q - 1)] : 0.0);
        const double t = lgs(cc), d = wa[static_cast<std::size_t>(q)] * t * (1 - t);
        tb += wa[static_cast<std::size_t>(q)] * t;
        G[static_cast<std::size_t>(j) * np + 2 * (ng + na - 1) + j] += d;
        if (q > 0) G[static_cast<std::size_t>(j) * np + 2 * (ng + na - 1) + ng + q - 1] += d;
        // cfit: nc x 5: phi, separation, threshold, search and hit rate
        const double al = best[static_cast<std::size_t>(j)] + (q > 0 ? best[static_cast<std::size_t>(ng + q - 1)] : 0.0);
        const double be = best[static_cast<std::size_t>(ng + na - 1 + j)] +
                          (q > 0 ? best[static_cast<std::size_t>(2 * ng + na - 1 + q - 1)] : 0.0);
        double x = 0.0, A = 0.0, B = 0.0;
        ab(al, std::exp(be), cc, x, A, B);
        const std::size_t c = static_cast<std::size_t>(q * ng + j);
        cfit[c] = lgs(al);
        cfit[nc + c] = std::exp(be);
        cfit[2 * nc + c] = t;
        cfit[3 * nc + c] = A + B;
        cfit[4 * nc + c] = A + B > 0 ? A / (A + B) : kutil::nan();
      }
      thr[j] = tb;
    }
    // tse: ng x 3: se of the average, difference from ref and its se
    auto quad = [&](const double *x1, const double *x2) {
      double v = 0.0;
      for (int p1 = 0; p1 < np; ++p1) {
        for (int p2 = 0; p2 < np; ++p2) {
          v += (x1[p1] - (x2 ? x2[p1] : 0.0)) * cov[static_cast<std::size_t>(p2) * np + p1] *
               (x1[p2] - (x2 ? x2[p2] : 0.0));
        }
      }
      return v >= 0 ? std::sqrt(v) : kutil::nan();
    };
    const double *gr = &G[static_cast<std::size_t>(ref) * np];
    for (int j = 0; j < ng; ++j) {
      const double *gj = &G[static_cast<std::size_t>(j) * np];
      tse[j] = quad(gj, nullptr);
      tse[ng + j] = thr[j] - thr[ref];
      tse[2 * ng + j] = j == ref ? kutil::nan() : quad(gj, gr);
    }
  }
  lrt[0] = std::fmax(2 * (ll[0] - ll[1]), 0.0);
  lrt[1] = ng - 1;
  lrt[2] = clab_pchisq(lrt[0], lrt[1], 0);
  lrt[3] = ll[0];
  lrt[4] = ll[1];
  lrt[5] = static_cast<double>(cs.size());
  return CLAB_OK;
}

int clab_thrrk(const double *stops, const double *srch, const double *hits, const int *o,
               const int *g, int n, int no, int ng, double zc, double *rate, double *res) {
  if (n < 1 || no < 2 || ng < 2 || !(zc > 0)) return CLAB_EARG;
  if (static_cast<long long>(no) * ng > 50000000LL) return CLAB_EARG;
  const std::size_t nc = static_cast<std::size_t>(no) * ng;
  // rate: nc x 5 (cell o + no g): stops, searches, hits, search and hit rate
  for (std::size_t c = 0; c < 3 * nc; ++c) rate[c] = 0.0;
  for (int i = 0; i < n; ++i) {
    if (o[i] < 0 || o[i] >= no || g[i] < 0 || g[i] >= ng) return CLAB_EARG;
    if (!kutil::nonneg(stops[i]) || !kutil::nonneg(srch[i]) || !kutil::nonneg(hits[i]) ||
        srch[i] > stops[i] || hits[i] > srch[i]) {
      return CLAB_EDOM;
    }
    const std::size_t c = static_cast<std::size_t>(o[i]) + static_cast<std::size_t>(no) * g[i];
    rate[c] += stops[i];
    rate[nc + c] += srch[i];
    rate[2 * nc + c] += hits[i];
  }
  for (std::size_t c = 0; c < nc; ++c) {
    rate[3 * nc + c] = rate[c] > 0 ? rate[nc + c] / rate[c] : kutil::nan();
    rate[4 * nc + c] = rate[nc + c] > 0 ? rate[2 * nc + c] / rate[nc + c] : kutil::nan();
  }
  // res: one row per pair of officer groups (a < b) and measure (search,
  // hit): zmin, zmax over the suspects' groups of the z of a's rate less
  // b's, and the verdict (1: significant differences of both signs)
  int row = 0;
  const int np = no * (no - 1);
  for (int x = 0; x < no; ++x) {
    for (int y = x + 1; y < no; ++y) {
      for (int ms = 0; ms < 2; ++ms) {
        double lo = HUGE_VAL, hi = -HUGE_VAL;
        for (int j = 0; j < ng; ++j) {
          const std::size_t cx = static_cast<std::size_t>(x) + static_cast<std::size_t>(no) * j;
          const std::size_t cy = static_cast<std::size_t>(y) + static_cast<std::size_t>(no) * j;
          const double nx = rate[(ms ? nc : 0) + cx], ny = rate[(ms ? nc : 0) + cy];
          const double px = rate[(3 + ms) * nc + cx], py = rate[(3 + ms) * nc + cy];
          if (!(nx > 0) || !(ny > 0)) continue;
          const double se = std::sqrt(px * (1 - px) / nx + py * (1 - py) / ny);
          if (!(se > 0)) continue;
          const double zz = (px - py) / se;
          lo = std::fmin(lo, zz);
          hi = std::fmax(hi, zz);
        }
        const bool ok = lo <= hi;
        res[row] = ok ? lo : kutil::nan();
        res[np + row] = ok ? hi : kutil::nan();
        res[2 * np + row] = ok && lo < -zc && hi > zc ? 1.0 : 0.0;
        ++row;
      }
    }
  }
  return CLAB_OK;
}
