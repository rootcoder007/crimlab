// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, tests of discrimination in searches: the hit-rate (outcome)
// test of Knowles, Persico and Todd (2001), the threshold test with the
// discriminant distributions of Pierson, Corbett-Davies and Goel (2018)
// fitted by maximum likelihood or by Hamiltonian Monte Carlo under the
// priors of Simoiu, Corbett-Davies and Goel (2017), and the rank-order
// test of Anwar and Fang (2006), pairwise and joint.
#include <algorithm>
#include <cmath>
#include <cstdint>
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
               const int *g, int n, int no, int ng, double zc, double zj, double *rate, double *res,
               double *mono, double *jt) {
  if (n < 1 || no < 2 || ng < 2 || !(zc > 0) || !(zj > 0)) return CLAB_EARG;
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
  // mono: ng x 2 rows (suspect group j, measure ms: row 2 j + ms) x 3:
  // Pearson's chi-square of equal rates across the officer groups (Anwar
  // and Fang 2006, Section 3.2.3), its df and p
  const int nm = 2 * ng;
  for (int j = 0; j < ng; ++j) {
    for (int ms = 0; ms < 2; ++ms) {
      double N = 0, K = 0;
      int used = 0;
      for (int x = 0; x < no; ++x) {
        const std::size_t c = static_cast<std::size_t>(x) + static_cast<std::size_t>(no) * j;
        const double d = rate[(ms ? nc : 0) + c];
        if (!(d > 0)) continue;
        N += d;
        K += rate[(ms ? 2 * nc : nc) + c];
        ++used;
      }
      const int row = 2 * j + ms;
      const double pp = N > 0 ? K / N : kutil::nan();
      if (used < 2 || !(pp > 0 && pp < 1)) {
        mono[row] = mono[nm + row] = mono[2 * nm + row] = kutil::nan();
        continue;
      }
      double x2 = 0;
      for (int x = 0; x < no; ++x) {
        const std::size_t c = static_cast<std::size_t>(x) + static_cast<std::size_t>(no) * j;
        const double d = rate[(ms ? nc : 0) + c];
        if (!(d > 0)) continue;
        const double o1 = rate[(ms ? 2 * nc : nc) + c], e1 = d * pp, o0 = d - o1, e0 = d - e1;
        x2 += (o1 - e1) * (o1 - e1) / e1 + (o0 - e0) * (o0 - e0) / e0;
      }
      mono[row] = x2;
      mono[nm + row] = used - 1;
      mono[2 * nm + row] = clab_pchisq(x2, used - 1, 0);
    }
  }
  // jt: 2 rows (measure) x (1 + no): the joint verdict and, when there is
  // one, the common rank of each officer group (1 the highest rate). Each
  // one-sided z of eq. 7 (sample variances of the 0/1 outcomes) beyond zj
  // in a suspect group orders that pair; the null of no prejudice holds
  // when one ranking of the officer groups agrees with every such order
  // in every suspect group (Section 2.1.2), so it is rejected when the
  // orders form a cycle.
  for (int ms = 0; ms < 2; ++ms) {
    std::vector<char> ab(static_cast<std::size_t>(no) * no, 0);  // ab[x no + y]: x above y
    for (int j = 0; j < ng; ++j) {
      for (int x = 0; x < no; ++x) {
        for (int y = 0; y < no; ++y) {
          if (x == y) continue;
          const std::size_t cx = static_cast<std::size_t>(x) + static_cast<std::size_t>(no) * j;
          const std::size_t cy = static_cast<std::size_t>(y) + static_cast<std::size_t>(no) * j;
          const double nx = rate[(ms ? nc : 0) + cx], ny = rate[(ms ? nc : 0) + cy];
          if (!(nx > 1) || !(ny > 1)) continue;
          const double px = rate[(3 + ms) * nc + cx], py = rate[(3 + ms) * nc + cy];
          const double se = std::sqrt(px * (1 - px) / (nx - 1) + py * (1 - py) / (ny - 1));
          if (se > 0 && (px - py) / se > zj) ab[static_cast<std::size_t>(x) * no + y] = 1;
        }
      }
    }
    // Kahn's order: take the groups with nobody required above them, the
    // highest pooled rate first; a cycle leaves groups untaken
    std::vector<double> pr(static_cast<std::size_t>(no), 0.0);
    for (int x = 0; x < no; ++x) {
      double d = 0, k = 0;
      for (int j = 0; j < ng; ++j) {
        const std::size_t c = static_cast<std::size_t>(x) + static_cast<std::size_t>(no) * j;
        d += rate[(ms ? nc : 0) + c];
        k += rate[(ms ? 2 * nc : nc) + c];
      }
      pr[static_cast<std::size_t>(x)] = d > 0 ? k / d : -1.0;
    }
    std::vector<char> done(static_cast<std::size_t>(no), 0);
    std::vector<int> rk(static_cast<std::size_t>(no), 0);
    int taken = 0;
    for (; taken < no; ++taken) {
      int best = -1;
      for (int x = 0; x < no; ++x) {
        if (done[static_cast<std::size_t>(x)]) continue;
        bool free = true;
        for (int y = 0; y < no && free; ++y) {
          if (!done[static_cast<std::size_t>(y)] && ab[static_cast<std::size_t>(y) * no + x]) free = false;
        }
        if (free && (best < 0 || pr[static_cast<std::size_t>(x)] > pr[static_cast<std::size_t>(best)])) best = x;
      }
      if (best < 0) break;
      done[static_cast<std::size_t>(best)] = 1;
      rk[static_cast<std::size_t>(best)] = taken + 1;
    }
    const bool rej = taken < no;
    jt[ms] = rej ? 1.0 : 0.0;
    for (int x = 0; x < no; ++x) jt[2 * (1 + x) + ms] = rej ? kutil::nan() : rk[static_cast<std::size_t>(x)];
  }
  return CLAB_OK;
}

// ------------------------------------------------------------------------
// The Bayesian fit of the threshold model (clab_thrbay): the hierarchical
// priors of Simoiu, Corbett-Davies and Goel (2017, Section 2) on the
// discriminant distributions of Pierson, Corbett-Davies and Goel (2017),
// sampled by Hamiltonian Monte Carlo with the No-U-Turn sampler.
// ------------------------------------------------------------------------

namespace {

// The log-likelihood of one cell at (al, be, ga) = (logit phi, log delta,
// logit t) and its gradient; false where it is not finite.
bool cll(const cell &e, double al, double be, double ga, double &l, double *gr) {
  const double dl = std::exp(be);
  double x = 0.0, A = 0.0, B = 0.0;
  ab(al, dl, ga, x, A, B);
  const double r = 1.0 - A - B;
  const double km = e.k - e.m, nk = e.n - e.k;
  if ((e.m > 0 && !(A > 0)) || (km > 0 && !(B > 0)) || (nk > 0 && !(r > 0)) || !std::isfinite(dl)) return false;
  l = 0.0;
  if (e.m > 0) l += e.m * std::log(A);
  if (km > 0) l += km * std::log(B);
  if (nk > 0) l += nk * std::log(r);
  const double ph = lgs(al), f1 = dn(x - dl), f0 = dn(x);
  const double q1 = qup(x - dl), q0 = qup(x);
  const double xa = -1.0 / dl, xg = 1.0 / dl, xd = 0.5 - (ga - al) / (dl * dl);
  const double Aa = ph * (1 - ph) * q1 - ph * f1 * xa, Ba = -ph * (1 - ph) * q0 - (1 - ph) * f0 * xa;
  const double Ag = -ph * f1 * xg, Bg = -(1 - ph) * f0 * xg;
  const double Ad = -ph * f1 * (xd - 1.0), Bd = -(1 - ph) * f0 * xd;
  auto d = [&](double dA, double dB) {
    double v = 0.0;
    if (e.m > 0) v += e.m * dA / A;
    if (km > 0) v += km * dB / B;
    if (nk > 0) v -= nk * (dA + dB) / r;
    return v;
  };
  gr[0] = d(Aa, Ba);
  gr[1] = d(Ad, Bd) * dl;
  gr[2] = d(Ag, Bg);
  return std::isfinite(l) && std::isfinite(gr[0]) && std::isfinite(gr[1]) && std::isfinite(gr[2]);
}

// The posterior. Unconstrained parameters u: phi_r (R), delta_r (R), the
// standardised area terms of phi and delta (D1 each, D1 = D - 1: the
// largest area is pinned at 0), mu_phi, log sigma_phi, mu_delta, log
// sigma_delta, mu_t (R), log sigma_t (R) and the standardised logit
// thresholds of the cells (C). phi_d = mu_phi + sigma_phi z_d (so phi_d ~
// N(mu_phi, sigma_phi)), likewise delta_d, and logit t_rd = mu_t,r +
// sigma_t,r z_rd; N(0, 2) priors on phi_r, delta_r, the mu's, and
// half-normal N+(0, 2) on the sigma's (with the Jacobian of the log).
struct bay {
  const std::vector<cell> *c;
  std::vector<int> cg, cd;  // each cell's group and area index (-1 pinned)
  int R, D1, C;
  int iphr() const { return 0; }
  int idr() const { return R; }
  int izp() const { return 2 * R; }
  int izd() const { return 2 * R + D1; }
  int ih() const { return 2 * R + 2 * D1; }
  int imt() const { return ih() + 4; }
  int ist() const { return imt() + R; }
  int izt() const { return ist() + R; }
  int np() const { return izt() + C; }
};

double lpost(const bay &M, const double *u, double *g) {
  const int np = M.np();
  for (int q = 0; q < np; ++q) g[q] = 0.0;
  double lp = 0.0;
  // N(0, 2): -x^2 / 8
  auto n2 = [&](int q) {
    lp -= u[q] * u[q] / 8;
    g[q] -= u[q] / 4;
  };
  auto n1 = [&](int q) {
    lp -= u[q] * u[q] / 2;
    g[q] -= u[q];
  };
  // half-normal N+(0, 2) on s = exp(v): -s^2 / 8 + v
  auto hn = [&](int q) {
    const double s2 = std::exp(2 * u[q]);
    lp += -s2 / 8 + u[q];
    g[q] += -s2 / 4 + 1;
  };
  for (int j = 0; j < M.R; ++j) {
    n2(M.iphr() + j);
    n2(M.idr() + j);
    n2(M.imt() + j);
    hn(M.ist() + j);
  }
  for (int d = 0; d < M.D1; ++d) {
    n1(M.izp() + d);
    n1(M.izd() + d);
  }
  n2(M.ih());
  hn(M.ih() + 1);
  n2(M.ih() + 2);
  hn(M.ih() + 3);
  for (int z = 0; z < M.C; ++z) n1(M.izt() + z);
  const double sp = std::exp(u[M.ih() + 1]), sd = std::exp(u[M.ih() + 3]);
  kutil::ksum s;
  for (int z = 0; z < M.C; ++z) {
    const int j = M.cg[static_cast<std::size_t>(z)], d = M.cd[static_cast<std::size_t>(z)];
    const double zp = d >= 0 ? u[M.izp() + d] : 0.0, zd = d >= 0 ? u[M.izd() + d] : 0.0;
    const double st = std::exp(u[M.ist() + j]), zt = u[M.izt() + z];
    const double al = u[M.iphr() + j] + (d >= 0 ? u[M.ih()] + sp * zp : 0.0);
    const double be = u[M.idr() + j] + (d >= 0 ? u[M.ih() + 2] + sd * zd : 0.0);
    const double ga = u[M.imt() + j] + st * zt;
    double l = 0.0, gr[3];
    if (!cll((*M.c)[static_cast<std::size_t>(z)], al, be, ga, l, gr)) return -HUGE_VAL;
    s.add(l);
    g[M.iphr() + j] += gr[0];
    g[M.idr() + j] += gr[1];
    g[M.imt() + j] += gr[2];
    g[M.izt() + z] += gr[2] * st;
    g[M.ist() + j] += gr[2] * st * zt;
    if (d >= 0) {
      g[M.ih()] += gr[0];
      g[M.izp() + d] += gr[0] * sp;
      g[M.ih() + 1] += gr[0] * sp * zp;
      g[M.ih() + 2] += gr[1];
      g[M.izd() + d] += gr[1] * sd;
      g[M.ih() + 3] += gr[1] * sd * zd;
    }
  }
  const double v = lp + s.val();
  return std::isfinite(v) ? v : -HUGE_VAL;
}

// A standard normal draw (Box-Muller).
double rnorm(kutil::rng &r) {
  const double a = 1.0 - r.unif(), b = r.unif();
  return std::sqrt(-2.0 * std::log(a)) * std::cos(6.283185307179586477 * b);
}

// The quantile of type 7 of sorted x.
double q7(const std::vector<double> &x, double p) {
  const double h = (static_cast<double>(x.size()) - 1) * p;
  const std::size_t lo = static_cast<std::size_t>(std::floor(h));
  const std::size_t hi = lo + 1 < x.size() ? lo + 1 : lo;
  return x[lo] + (h - static_cast<double>(lo)) * (x[hi] - x[lo]);
}

// Split R-hat and the effective sample size of Gelman et al. (2013,
// Bayesian Data Analysis, 3rd ed., Sections 11.4 and 11.5) of m chains of
// nd draws each, chain-major in x: the chains are split in halves, and
// the autocorrelations summed in pairs up to the first pair below 0
// (Geyer's initial positive sequence).
void diag(const double *x, int m, int nd, double &rhat, double &ess) {
  const int h = nd / 2, mm = 2 * m;
  std::vector<double> mu(static_cast<std::size_t>(mm)), va(static_cast<std::size_t>(mm));
  auto at = [&](int c, int i) { return x[static_cast<std::size_t>(c / 2) * nd + (c % 2) * (nd - h) + i]; };
  for (int c = 0; c < mm; ++c) {
    double s = 0;
    for (int i = 0; i < h; ++i) s += at(c, i);
    mu[static_cast<std::size_t>(c)] = s / h;
    double v = 0;
    for (int i = 0; i < h; ++i) v += (at(c, i) - s / h) * (at(c, i) - s / h);
    va[static_cast<std::size_t>(c)] = v / (h - 1);
  }
  double W = 0, mb = 0;
  for (int c = 0; c < mm; ++c) {
    W += va[static_cast<std::size_t>(c)] / mm;
    mb += mu[static_cast<std::size_t>(c)] / mm;
  }
  double B = 0;
  for (int c = 0; c < mm; ++c) B += (mu[static_cast<std::size_t>(c)] - mb) * (mu[static_cast<std::size_t>(c)] - mb);
  B *= static_cast<double>(h) / (mm - 1);
  const double vp = (h - 1.0) / h * W + B / h;
  if (!(W > 0) || !(vp > 0)) {
    rhat = ess = kutil::nan();
    return;
  }
  rhat = std::sqrt(vp / W);
  auto rho = [&](int t) {
    double v = 0;
    for (int c = 0; c < mm; ++c) {
      for (int i = t; i < h; ++i) {
        const double d = at(c, i) - at(c, i - t);
        v += d * d;
      }
    }
    v /= static_cast<double>(mm) * (h - t);
    return 1.0 - v / (2 * vp);
  };
  // Geyer: the sums of the pairs (rho 2k, rho 2k + 1), k = 0, 1, ..., up
  // to the first negative one; tau = -1 + 2 sum
  double sum = 0;
  for (int t = 0; t + 1 < h; t += 2) {
    const double pr = rho(t) + rho(t + 1);
    if (pr < 0) break;
    sum += pr;
  }
  ess = static_cast<double>(mm) * h / (-1 + 2 * sum);
}


// One state of a trajectory: position, momentum, gradient, log posterior.
struct st {
  std::vector<double> u, p, g;
  double lp;
};

double kin(const std::vector<double> &p, const std::vector<double> &minv) {
  double k = 0;
  for (std::size_t q = 0; q < p.size(); ++q) k += 0.5 * minv[q] * p[q] * p[q];
  return k;
}

// No U-turn between the ends: (u+ - u-) . M^-1 p >= 0 at both.
bool noturn(const st &lo, const st &hi, const std::vector<double> &minv) {
  double a = 0, b = 0;
  for (std::size_t q = 0; q < lo.u.size(); ++q) {
    const double d = hi.u[q] - lo.u[q];
    a += d * minv[q] * lo.p[q];
    b += d * minv[q] * hi.p[q];
  }
  return a >= 0 && b >= 0;
}

struct tree {
  st lo, hi, pr;  // the ends and the proposal
  double n = 0, a = 0, na = 0;
  bool ok = true, div = false;
};

// The recursive doubling of Hoffman and Gelman (2014, Algorithm 6): a
// subtree of 2^j leapfrog steps in direction v from s.
void build(const bay &M, const std::vector<double> &minv, const st &s, double logu, int v, int j, double eps,
           double h0, kutil::rng &rg, tree &t) {
  if (j == 0) {
    st x = s;
    const std::size_t np = x.u.size();
    for (std::size_t q = 0; q < np; ++q) x.p[q] += 0.5 * v * eps * x.g[q];
    for (std::size_t q = 0; q < np; ++q) x.u[q] += v * eps * minv[q] * x.p[q];
    x.lp = lpost(M, x.u.data(), x.g.data());
    for (std::size_t q = 0; q < np; ++q) x.p[q] += 0.5 * v * eps * x.g[q];
    const double h = std::isfinite(x.lp) ? x.lp - kin(x.p, minv) : -HUGE_VAL;
    t.lo = x;
    t.hi = x;
    t.pr = x;
    t.n = logu <= h ? 1 : 0;
    t.ok = logu < h + 1000;
    t.div = !t.ok;
    t.a = std::isfinite(h) ? std::fmin(1.0, std::exp(h - h0)) : 0.0;
    t.na = 1;
    return;
  }
  build(M, minv, s, logu, v, j - 1, eps, h0, rg, t);
  if (!t.ok) return;
  tree u;
  build(M, minv, v < 0 ? t.lo : t.hi, logu, v, j - 1, eps, h0, rg, u);
  if (v < 0) {
    t.lo = u.lo;
  } else {
    t.hi = u.hi;
  }
  if (u.n > 0 && rg.unif() < u.n / (t.n + u.n)) t.pr = u.pr;
  t.a += u.a;
  t.na += u.na;
  t.div = t.div || u.div;
  t.ok = u.ok && noturn(t.lo, t.hi, minv);
  t.n += u.n;
}

// One draw of the No-U-Turn sampler from (u, gr, lp) with momentum p, at
// most 2^10 leapfrog steps; ac the mean acceptance statistic of the
// trajectory (for the step size), div whether it diverged.
void nuts(const bay &M, const std::vector<double> &minv, double eps, kutil::rng &rg, std::vector<double> &u,
          std::vector<double> &gr, double &lp, const std::vector<double> &p, double &ac, bool &div) {
  st s{u, p, gr, lp};
  const double h0 = lp - kin(p, minv);
  const double logu = h0 + std::log(1.0 - rg.unif());
  st lo = s, hi = s;
  double n = 1, sa = 0, sna = 0;
  div = false;
  for (int j = 0; j < 10; ++j) {
    const int v = rg.unif() < 0.5 ? -1 : 1;
    tree t;
    build(M, minv, v < 0 ? lo : hi, logu, v, j, eps, h0, rg, t);
    if (v < 0) {
      lo = t.lo;
    } else {
      hi = t.hi;
    }
    sa += t.a;
    sna += t.na;
    div = div || t.div;
    if (t.ok && rg.unif() < t.n / n) {
      u = t.pr.u;
      gr = t.pr.g;
      lp = t.pr.lp;
    }
    n += t.n;
    if (!t.ok || !noturn(lo, hi, minv)) break;
  }
  ac = sna > 0 ? sa / sna : 0.0;
}
}  // namespace

int clab_thrbay(const double *stops, const double *srch, const double *hits, const int *g, const int *a, int n,
                int ng, int na, int ref, int chains, int iter, unsigned int seed, double level, const double *init,
                double *post, double *cfit, double *info) {
  if (n < 1 || ng < 2 || na < 1 || ref < 0 || ref >= ng || chains < 2 || iter < 20 || !(level > 0 && level < 1)) {
    return CLAB_EARG;
  }
  if (static_cast<long long>(ng) * na > 1000000LL || static_cast<double>(chains) * iter > 1e8) return CLAB_EARG;
  const int nc = ng * na;
  std::vector<double> S(static_cast<std::size_t>(3 * nc), 0.0);
  for (int i = 0; i < n; ++i) {
    if (g[i] < 0 || g[i] >= ng || a[i] < 0 || a[i] >= na) return CLAB_EARG;
    const double s0 = stops[i], s1 = srch[i], s2 = hits[i];
    if (!kutil::nonneg(s0) || !kutil::nonneg(s1) || !kutil::nonneg(s2) || s1 > s0 || s2 > s1) return CLAB_EDOM;
    const std::size_t c = static_cast<std::size_t>(a[i] * ng + g[i]);
    S[3 * c] += s0;
    S[3 * c + 1] += s1;
    S[3 * c + 2] += s2;
  }
  // the area pinned at 0: the one with the most stops (the first of ties)
  std::vector<double> wa(static_cast<std::size_t>(na), 0.0);
  for (int q = 0; q < na; ++q) {
    for (int j = 0; j < ng; ++j) wa[static_cast<std::size_t>(q)] += S[3 * static_cast<std::size_t>(q * ng + j)];
  }
  int big = 0;
  for (int q = 1; q < na; ++q) big = wa[static_cast<std::size_t>(q)] > wa[static_cast<std::size_t>(big)] ? q : big;
  std::vector<int> di(static_cast<std::size_t>(na), -1);
  for (int q = 0, k = 0; q < na; ++q) {
    if (q != big) di[static_cast<std::size_t>(q)] = k++;
  }
  std::vector<cell> cs;
  bay M;
  M.c = &cs;
  std::vector<int> cix;  // the grid cell of each model cell
  std::vector<double> tk(static_cast<std::size_t>(ng), 0.0);
  for (int q = 0; q < na; ++q) {
    for (int j = 0; j < ng; ++j) {
      const std::size_t c = static_cast<std::size_t>(q * ng + j);
      tk[static_cast<std::size_t>(j)] += S[3 * c + 1];
      if (!(S[3 * c] > 0)) continue;
      cs.push_back(cell{S[3 * c], S[3 * c + 1], S[3 * c + 2], 0, 0, 0, 0, 0, 0});
      M.cg.push_back(j);
      M.cd.push_back(di[static_cast<std::size_t>(q)]);
      cix.push_back(static_cast<int>(c));
    }
  }
  for (double t : tk) {
    if (!(t > 0)) return CLAB_EDOM;
  }
  M.R = ng;
  M.D1 = na - 1;
  M.C = static_cast<int>(cs.size());
  // identified by the data: more rates (2 per cell) than the signal and
  // threshold terms, 2 (R + D - 1) + C
  if (M.C <= 2 * (ng + na - 1)) return CLAB_EDOM;
  const int np = M.np(), R = ng, D1 = na - 1, C = M.C;
  // reported quantities: thr (R), diff (R), phi_r, delta_r, mu_phi,
  // sigma_phi, mu_delta, sigma_delta, phi_d (D1), delta_d (D1), mu_t (R),
  // sigma_t (R), the cells' thresholds (C)
  const int nq = 6 * R + 4 + 2 * D1 + C;
  const int warm = iter / 2, nd = iter - warm;
  std::vector<double> dr(static_cast<std::size_t>(nq) * chains * nd);
  std::vector<double> cm(static_cast<std::size_t>(5) * C, 0.0);
  std::vector<double> u(static_cast<std::size_t>(np)), gr(static_cast<std::size_t>(np)), p(static_cast<std::size_t>(np)),
      minv(static_cast<std::size_t>(np));
  double ndiv = 0, accs = 0, epss = 0;
  // the share of each area's stops, for the groups' average thresholds
  std::vector<double> wg(static_cast<std::size_t>(R), 0.0);
  for (int z = 0; z < C; ++z) {
    wg[static_cast<std::size_t>(M.cg[static_cast<std::size_t>(z)])] += wa[static_cast<std::size_t>(cix[static_cast<std::size_t>(z)] / ng)];
  }
  // the start of the chains from init: the groups' terms at the pinned
  // area, the areas' terms as their mean offsets from it, standardised by
  // their spread, and the cells' thresholds by their group's
  std::vector<double> u0(static_cast<std::size_t>(np), 0.0);
  if (init) {
    auto iv = [&](int k, int j, int q) { return init[static_cast<std::size_t>(k) * nc + static_cast<std::size_t>(q * ng + j)]; };
    for (int k = 0; k < 2; ++k) {
      std::vector<double> off(static_cast<std::size_t>(D1), 0.0);
      for (int j = 0; j < R; ++j) u0[static_cast<std::size_t>((k ? M.idr() : M.iphr()) + j)] = iv(k, j, big);
      for (int q = 0; q < na; ++q) {
        const int dq = di[static_cast<std::size_t>(q)];
        if (dq < 0) continue;
        for (int j = 0; j < R; ++j) off[static_cast<std::size_t>(dq)] += (iv(k, j, q) - iv(k, j, big)) / R;
      }
      double m = 0, v = 0;
      for (double o : off) m += o / std::max(D1, 1);
      for (double o : off) v += (o - m) * (o - m) / std::max(D1 - 1, 1);
      const double s = std::fmax(std::sqrt(v), 0.05);
      u0[static_cast<std::size_t>(M.ih() + 2 * k)] = m;
      u0[static_cast<std::size_t>(M.ih() + 2 * k + 1)] = std::log(s);
      for (int d = 0; d < D1; ++d) u0[static_cast<std::size_t>((k ? M.izd() : M.izp()) + d)] = (off[static_cast<std::size_t>(d)] - m) / s;
    }
    std::vector<double> mt(static_cast<std::size_t>(R), 0.0), vt(static_cast<std::size_t>(R), 0.0), nt(static_cast<std::size_t>(R), 0.0);
    for (int z = 0; z < C; ++z) {
      const int j = M.cg[static_cast<std::size_t>(z)];
      mt[static_cast<std::size_t>(j)] += init[2 * static_cast<std::size_t>(nc) + static_cast<std::size_t>(cix[static_cast<std::size_t>(z)])];
      nt[static_cast<std::size_t>(j)] += 1;
    }
    for (int j = 0; j < R; ++j) mt[static_cast<std::size_t>(j)] /= nt[static_cast<std::size_t>(j)];
    for (int z = 0; z < C; ++z) {
      const int j = M.cg[static_cast<std::size_t>(z)];
      const double e = init[2 * static_cast<std::size_t>(nc) + static_cast<std::size_t>(cix[static_cast<std::size_t>(z)])] - mt[static_cast<std::size_t>(j)];
      vt[static_cast<std::size_t>(j)] += e * e / std::max(nt[static_cast<std::size_t>(j)] - 1, 1.0);
    }
    for (int j = 0; j < R; ++j) {
      const double s = std::fmax(std::sqrt(vt[static_cast<std::size_t>(j)]), 0.05);
      u0[static_cast<std::size_t>(M.imt() + j)] = mt[static_cast<std::size_t>(j)];
      u0[static_cast<std::size_t>(M.ist() + j)] = std::log(s);
    }
    for (int z = 0; z < C; ++z) {
      const int j = M.cg[static_cast<std::size_t>(z)];
      u0[static_cast<std::size_t>(M.izt() + z)] =
          (init[2 * static_cast<std::size_t>(nc) + static_cast<std::size_t>(cix[static_cast<std::size_t>(z)])] - mt[static_cast<std::size_t>(j)]) /
          std::exp(u0[static_cast<std::size_t>(M.ist() + j)]);
    }
    for (double v : u0) {
      if (!std::isfinite(v)) return CLAB_EDOM;
    }
  }
  for (int ch = 0; ch < chains; ++ch) {
    kutil::rng rg(static_cast<std::uint64_t>(seed) * 1000003ULL + static_cast<std::uint64_t>(ch) + 1);
    // the start: the given fit moved by up to 0.25, or drawn uniformly
    // in (-2, 2) as Stan does, until finite
    double lp = -HUGE_VAL;
    for (int tr = 0; tr < 100 && !std::isfinite(lp); ++tr) {
      for (std::size_t q = 0; q < u.size(); ++q) u[q] = init ? u0[q] + 0.5 * rg.unif() - 0.25 : 4 * rg.unif() - 2;
      lp = lpost(M, u.data(), gr.data());
    }
    if (!std::isfinite(lp)) return CLAB_EDOM;
    for (double &v : minv) v = 1.0;
    // the step size: dual averaging to an acceptance of 0.8 (Hoffman and
    // Gelman 2014) in the first 15% and the last 10% of the warm-up, the
    // diagonal of the inverse mass matrix from the draws in between
    double eps = 0.1, mu = std::log(10 * eps), hb = 0, lb = 0;
    int ta = 0;
    const int w1 = std::max(1, warm * 15 / 100), w2 = warm - std::max(1, warm / 10);
    std::vector<double> ws(static_cast<std::size_t>(np), 0.0), ws2(static_cast<std::size_t>(np), 0.0);
    int nw = 0;
    for (int it = 0; it < iter; ++it) {
      if (kutil::stop(64.0 * C * 50)) return CLAB_EINTR;
      for (int q = 0; q < np; ++q) p[static_cast<std::size_t>(q)] = rnorm(rg) / std::sqrt(minv[static_cast<std::size_t>(q)]);
      double ac = 0;
      bool div = false;
      nuts(M, minv, eps, rg, u, gr, lp, p, ac, div);
      if (it < warm) {
        // dual averaging (gamma 0.05, t0 10, kappa 0.75)
        ++ta;
        hb = (1 - 1.0 / (ta + 10)) * hb + (0.8 - ac) / (ta + 10);
        const double le = mu - std::sqrt(static_cast<double>(ta)) / 0.05 * hb;
        const double w = std::pow(static_cast<double>(ta), -0.75);
        lb = w * le + (1 - w) * lb;
        eps = std::exp(le);
        if (it >= w1 && it < w2) {
          for (int q = 0; q < np; ++q) {
            ws[static_cast<std::size_t>(q)] += u[static_cast<std::size_t>(q)];
            ws2[static_cast<std::size_t>(q)] += u[static_cast<std::size_t>(q)] * u[static_cast<std::size_t>(q)];
          }
          ++nw;
        }
        if (it == w2 - 1 && nw > 2) {
          // the window's variances, shrunk towards 1e-3 as Stan does
          for (int q = 0; q < np; ++q) {
            const double m1 = ws[static_cast<std::size_t>(q)] / nw;
            const double v = std::fmax((ws2[static_cast<std::size_t>(q)] / nw - m1 * m1) * nw / (nw - 1.0), 0.0);
            minv[static_cast<std::size_t>(q)] = (nw / (nw + 5.0)) * v + 1e-3 * (5.0 / (nw + 5.0));
          }
          mu = std::log(10 * eps);
          hb = lb = 0;
          ta = 0;
        }
        if (it == warm - 1) eps = std::exp(lb);
        continue;
      }
      ndiv += div;
      accs += ac;
      // the draw's reported quantities
      const int k = it - warm;
      auto put = [&](int r, double v) {
        dr[(static_cast<std::size_t>(r) * chains + ch) * nd + static_cast<std::size_t>(k)] = v;
      };
      const double sp = std::exp(u[static_cast<std::size_t>(M.ih() + 1)]), sd = std::exp(u[static_cast<std::size_t>(M.ih() + 3)]);
      std::vector<double> tb(static_cast<std::size_t>(R), 0.0);
      for (int z = 0; z < C; ++z) {
        const int j = M.cg[static_cast<std::size_t>(z)], d = M.cd[static_cast<std::size_t>(z)];
        const double al = u[static_cast<std::size_t>(M.iphr() + j)] +
                          (d >= 0 ? u[static_cast<std::size_t>(M.ih())] + sp * u[static_cast<std::size_t>(M.izp() + d)] : 0.0);
        const double be = u[static_cast<std::size_t>(M.idr() + j)] +
                          (d >= 0 ? u[static_cast<std::size_t>(M.ih() + 2)] + sd * u[static_cast<std::size_t>(M.izd() + d)] : 0.0);
        const double ga = u[static_cast<std::size_t>(M.imt() + j)] +
                          std::exp(u[static_cast<std::size_t>(M.ist() + j)]) * u[static_cast<std::size_t>(M.izt() + z)];
        const double t = lgs(ga);
        tb[static_cast<std::size_t>(j)] += wa[static_cast<std::size_t>(cix[static_cast<std::size_t>(z)] / ng)] * t;
        put(6 * R + 4 + 2 * D1 + z, t);
        double x = 0.0, A = 0.0, B = 0.0;
        ab(al, std::exp(be), ga, x, A, B);
        cm[static_cast<std::size_t>(z)] += lgs(al);
        cm[static_cast<std::size_t>(C + z)] += std::exp(be);
        cm[static_cast<std::size_t>(2 * C + z)] += t;
        cm[static_cast<std::size_t>(3 * C + z)] += A + B;
        cm[static_cast<std::size_t>(4 * C + z)] += A + B > 0 ? A / (A + B) : 0.0;
      }
      for (int j = 0; j < R; ++j) tb[static_cast<std::size_t>(j)] /= wg[static_cast<std::size_t>(j)];
      for (int j = 0; j < R; ++j) {
        put(j, tb[static_cast<std::size_t>(j)]);
        put(R + j, tb[static_cast<std::size_t>(j)] - tb[static_cast<std::size_t>(ref)]);
        put(2 * R + j, u[static_cast<std::size_t>(M.iphr() + j)]);
        put(3 * R + j, u[static_cast<std::size_t>(M.idr() + j)]);
        put(4 * R + 4 + 2 * D1 + j, u[static_cast<std::size_t>(M.imt() + j)]);
        put(5 * R + 4 + 2 * D1 + j, std::exp(u[static_cast<std::size_t>(M.ist() + j)]));
      }
      put(4 * R, u[static_cast<std::size_t>(M.ih())]);
      put(4 * R + 1, sp);
      put(4 * R + 2, u[static_cast<std::size_t>(M.ih() + 2)]);
      put(4 * R + 3, sd);
      for (int d = 0; d < D1; ++d) {
        put(4 * R + 4 + d, u[static_cast<std::size_t>(M.ih())] + sp * u[static_cast<std::size_t>(M.izp() + d)]);
        put(4 * R + 4 + D1 + d, u[static_cast<std::size_t>(M.ih() + 2)] + sd * u[static_cast<std::size_t>(M.izd() + d)]);
      }
    }
    epss += eps;
  }
  // post: nq x 6 (mean, sd, lower and upper quantiles, ess, rhat)
  const std::size_t tot = static_cast<std::size_t>(chains) * nd;
  std::vector<double> sv(tot);
  for (int r = 0; r < nq; ++r) {
    const double *x = &dr[static_cast<std::size_t>(r) * tot];
    double m = 0;
    for (std::size_t i = 0; i < tot; ++i) m += x[i];
    m /= static_cast<double>(tot);
    double v = 0;
    for (std::size_t i = 0; i < tot; ++i) v += (x[i] - m) * (x[i] - m);
    sv.assign(x, x + tot);
    std::sort(sv.begin(), sv.end());
    double rh = 0, es = 0;
    diag(x, chains, nd, rh, es);
    const bool zero = r == R + ref;
    post[r] = m;
    post[nq + r] = std::sqrt(v / static_cast<double>(tot - 1));
    post[2 * nq + r] = q7(sv, 0.5 - 0.5 * level);
    post[3 * nq + r] = q7(sv, 0.5 + 0.5 * level);
    post[4 * nq + r] = zero ? kutil::nan() : es;
    post[5 * nq + r] = zero ? kutil::nan() : rh;
  }
  for (std::size_t c = 0; c < 5 * static_cast<std::size_t>(nc); ++c) cfit[c] = kutil::nan();
  for (int z = 0; z < C; ++z) {
    for (int k = 0; k < 5; ++k) {
      cfit[static_cast<std::size_t>(k) * nc + static_cast<std::size_t>(cix[static_cast<std::size_t>(z)])] =
          cm[static_cast<std::size_t>(k) * C + static_cast<std::size_t>(z)] / static_cast<double>(tot);
    }
  }
  info[0] = ndiv;
  info[1] = accs / static_cast<double>(tot);
  info[2] = epss / chains;
  info[3] = big;
  return CLAB_OK;
}
