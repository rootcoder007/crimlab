// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, forensics: the gamma model of DNA mixture peak heights
// (Cowell, Graversen, Lauritzen and Mortera 2015; Bleka, Storvik and Gill
// 2016) and the likelihood ratio of two accounts of a mixture with the
// unknown contributors summed over and the parameters fitted by maximum
// likelihood (gmixlr()). The height of the peak at a position is gamma
// with shape (sum of the contributors' allele doses times their mixture
// proportions, degraded and moved by stutter) / omega^2 and scale mu
// omega^2; a peak below the analytical threshold is not seen (drop-out);
// a seen peak that no contributor explains is a drop-in.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include "crimlab.h"
#include "frx.h"
#include "kutil.h"
#include "qn.h"

namespace {

// log(exp(a) + exp(b))
double lse(double a, double b) {
  if (a == -HUGE_VAL) return b;
  if (b == -HUGE_VAL) return a;
  const double m = std::max(a, b);
  return m + std::log1p(std::exp(std::min(a, b) - m));
}

// the log density and the log distribution function of the gamma
double ldg(double x, double a, double s) {
  return (a - 1) * std::log(x) - x / s - std::lgamma(a) - a * std::log(s);
}

}  // namespace

extern "C" {

int clab_gmxlik(int nl, const int *np, const double *p, const double *y, const double *sz,
                const int *st, int nk, const int *kg, int nu, const double *cn, const double *cl,
                double theta, double at, double pc, double lambda, int deg, int stut,
                const double *par, double *ll) {
  if (nl < 1 || nk < 0 || nu < 0 || nk + nu < 1 || nk + nu > 64) return CLAB_EARG;
  const double mu = par[0], om = par[1], be = par[2], xi = par[3];
  const double *m = par + 4;
  if (!(mu > 0) || !(om > 0) || !std::isfinite(mu) || !std::isfinite(om) || !(at > 0) ||
      !(pc >= 0 && pc < 1) || !(lambda > 0) || !(theta >= 0 && theta < 1) ||
      (deg && !(be > 0 && be <= 1)) || (stut && !(xi >= 0 && xi < 1))) {
    return CLAB_EDOM;
  }
  for (int k = 0; k < nk + nu; ++k) {
    if (!(m[k] >= 0 && m[k] <= 1)) return CLAB_EDOM;
  }
  const double sc = mu * om * om;
  const int nc = nk + nu;
  std::size_t off = 0;
  for (int l = 0; l < nl; ++l) {
    const int n = np[l];
    if (n < 1) return CLAB_EARG;
    const double *pl = p + off, *yl = y + off, *zl = sz ? sz + off : nullptr, *cnl = cn + off;
    const int *sl = st ? st + off : nullptr;
    // the genotypes a contributor may have: pairs of carriable positions
    std::vector<int> ga, gb;
    for (int i = 0; i < n; ++i) {
      if (!(pl[i] >= 0) || (yl[i] > 0 && yl[i] < at) || (sl && (sl[i] < -1 || sl[i] >= n))) return CLAB_EDOM;
      for (int j = i; j < n && pl[i] > 0; ++j) {
        if (pl[j] > 0) {
          ga.push_back(i);
          gb.push_back(j);
        }
      }
    }
    const std::size_t ng = ga.size();
    if (nu > 0 && ng == 0) return CLAB_EDOM;
    std::vector<double> dg(static_cast<std::size_t>(n), 1.0);
    for (int i = 0; deg && i < n; ++i) dg[static_cast<std::size_t>(i)] = std::pow(be, (zl[i] - 125.0) / 100.0);
    std::vector<int> pick(static_cast<std::size_t>(nu), 0);
    std::vector<double> c(static_cast<std::size_t>(n)), cnt(static_cast<std::size_t>(n));
    double tot = -HUGE_VAL;
    for (;;) {
      if (kutil::stop(static_cast<double>(n) * nc)) return CLAB_EINTR;
      // the doses at each position
      std::fill(c.begin(), c.end(), 0.0);
      bool bad = false;
      for (int k = 0; k < nc; ++k) {
        int a, b;
        if (k < nk) {
          a = kg[(static_cast<std::size_t>(l) * nk + k) * 2];
          b = kg[(static_cast<std::size_t>(l) * nk + k) * 2 + 1];
          if (a < 0 || a >= n || b < 0 || b >= n) bad = true;
          if (bad) break;
        } else {
          a = ga[static_cast<std::size_t>(pick[static_cast<std::size_t>(k - nk)])];
          b = gb[static_cast<std::size_t>(pick[static_cast<std::size_t>(k - nk)])];
        }
        c[static_cast<std::size_t>(a)] += m[k];
        c[static_cast<std::size_t>(b)] += m[k];
      }
      if (bad) return CLAB_EARG;
      double lg = 0.0;
      int nd = 0;
      for (int i = 0; i < n; ++i) {
        double s = c[static_cast<std::size_t>(i)] * dg[static_cast<std::size_t>(i)];
        if (stut) {
          s *= 1 - xi;
          for (int j = 0; j < n; ++j) {
            if (sl[j] == i) s += xi * c[static_cast<std::size_t>(j)] * dg[static_cast<std::size_t>(j)];
          }
        }
        const double a = s / (om * om);
        if (yl[i] > 0) {
          if (a > 0) {
            lg += ldg(yl[i], a, sc);
          } else {
            lg += std::log(pc * pl[i] * lambda) - lambda * (yl[i] - at);
            ++nd;
          }
        } else if (a > 0) {
          lg += std::log(clab_pgamma(at / sc, a, 1));
        }
      }
      if (nd == 0) lg += std::log1p(-pc);
      // the chance of the unknowns' genotypes, the alleles drawn in turn
      // given the conditioning profiles (Balding and Nichols)
      double lp = 0.0, nn = cl[l];
      for (int i = 0; i < n; ++i) cnt[static_cast<std::size_t>(i)] = cnl[i];
      for (int k = 0; k < nu; ++k) {
        const int g = pick[static_cast<std::size_t>(k)];
        for (int h = 0; h < 2; ++h) {
          const int a = h ? gb[static_cast<std::size_t>(g)] : ga[static_cast<std::size_t>(g)];
          lp += std::log((cnt[static_cast<std::size_t>(a)] * theta + (1 - theta) * pl[a]) / (1 + (nn - 1) * theta));
          cnt[static_cast<std::size_t>(a)] += 1;
          nn += 1;
        }
        if (ga[static_cast<std::size_t>(g)] != gb[static_cast<std::size_t>(g)]) lp += std::log(2.0);
      }
      tot = lse(tot, lg + lp);
      int k = nu - 1;
      while (k >= 0 && static_cast<std::size_t>(++pick[static_cast<std::size_t>(k)]) == ng) pick[static_cast<std::size_t>(k--)] = 0;
      if (k < 0) break;
    }
    ll[l] = tot;
    off += static_cast<std::size_t>(n);
  }
  return CLAB_OK;
}

}  // extern "C"

namespace frx {

namespace {

// One hypothesis: the loci laid out for clab_gmxlik(), and its fit.
struct Hyp {
  int nk = 0, nu = 0;
  std::vector<int> kg;
  std::vector<std::string> who;
};

struct Lay {
  int nl = 0;
  std::vector<int> np, st;
  std::vector<double> p, y, sz, cn, cl;
  bool deg = false, stut = false;
  double theta = 0, at = 0, pc = 0, lambda = 0;
};

struct Ctx {
  const Lay *ly;
  const Hyp *h;
  std::vector<double> *ll;
  int rc;
};

// parameters from the free values u: log mu, log omega, [logit beta],
// [logit xi], then the nc - 1 log-ratios of the proportions to the last
void unpack(const Lay &ly, int nc, const double *u, double *par) {
  int j = 2;
  par[0] = std::exp(u[0]);
  par[1] = std::exp(u[1]);
  par[2] = ly.deg ? 1 / (1 + std::exp(-u[j++])) : 1.0;
  par[3] = ly.stut ? 1 / (1 + std::exp(-u[j++])) : 0.0;
  double mx = 0.0;
  for (int k = 0; k < nc - 1; ++k) mx = std::max(mx, u[j + k]);
  double s = 0.0;
  for (int k = 0; k < nc; ++k) {
    par[4 + k] = std::exp((k < nc - 1 ? u[j + k] : 0.0) - mx);
    s += par[4 + k];
  }
  for (int k = 0; k < nc; ++k) par[4 + k] /= s;
}

double nll(const Ctx &c, const double *u) {
  const int nc = c.h->nk + c.h->nu;
  std::vector<double> par(4 + static_cast<std::size_t>(nc));
  unpack(*c.ly, nc, u, par.data());
  const Lay &ly = *c.ly;
  const int rc = clab_gmxlik(ly.nl, ly.np.data(), ly.p.data(), ly.y.data(), ly.deg ? ly.sz.data() : nullptr,
                             ly.stut ? ly.st.data() : nullptr, c.h->nk, c.h->kg.empty() ? nullptr : c.h->kg.data(),
                             c.h->nu, ly.cn.data(), ly.cl.data(), ly.theta, ly.at, ly.pc, ly.lambda, ly.deg,
                             ly.stut, par.data(), c.ll->data());
  if (rc != CLAB_OK) {
    const_cast<Ctx &>(c).rc = rc;
    return HUGE_VAL;
  }
  double s = 0.0;
  for (double v : *c.ll) s += v;
  return std::isfinite(s) ? -s : HUGE_VAL;
}

int npar(const Lay &ly, int nc) { return 2 + ly.deg + ly.stut + nc - 1; }

// the objective of qn::min with its gradient by central differences
double fobj(const double *u, double *g, void *vc) {
  Ctx &c = *static_cast<Ctx *>(vc);
  const int p = npar(*c.ly, c.h->nk + c.h->nu);
  const double f0 = nll(c, u);
  if (g != nullptr && std::isfinite(f0)) {
    std::vector<double> w(u, u + p);
    for (int k = 0; k < p; ++k) {
      const double h = 1e-5 * (1 + std::fabs(u[k]));
      w[static_cast<std::size_t>(k)] = u[k] + h;
      const double a = nll(c, w.data());
      w[static_cast<std::size_t>(k)] = u[k] - h;
      const double b = nll(c, w.data());
      w[static_cast<std::size_t>(k)] = u[k];
      g[k] = std::isfinite(a) && std::isfinite(b) ? (a - b) / (2 * h) : 0.0;
    }
  }
  return f0;
}

}  // namespace

int gmixlr(const GmxIn &in, GmxOut &out, std::string &err) {
  out = GmxOut();
  auto bad = [&err](const std::string &s) {
    err = s;
    return CLAB_EARG;
  };
  if (!(in.nc >= 1 && in.nc <= 5 && in.nc == std::floor(in.nc))) return bad("`nc` must be one whole number from 1 to 5.");
  if (!(in.at > 0) || !std::isfinite(in.at)) return bad("`at` must be one number above 0.");
  if (!(in.pc >= 0 && in.pc < 1)) return bad("`pc` must be one number from 0 and below 1.");
  if (!(in.lambda > 0) || !std::isfinite(in.lambda)) return bad("`lambda` must be one number above 0.");
  if (!(in.theta >= 0 && in.theta < 1)) return bad("`theta` must be one number from 0 and below 1.");
  const int nc = static_cast<int>(in.nc);
  Lay ly;
  ly.deg = in.deg;
  ly.stut = in.stut;
  ly.theta = in.theta;
  ly.at = in.at;
  ly.pc = in.pc;
  ly.lambda = in.lambda;
  if (in.deg && !in.hassz) return bad("`deg = TRUE` needs a column size in `freq`.");
  // frequencies, locus by locus in their order
  const std::size_t nf = in.floc.s.size();
  std::vector<std::string> loc;
  for (std::size_t i = 0; i < nf; ++i) {
    if (in.floc.na[i] || in.fal.na[i]) return bad("`freq` has a missing locus or allele.");
    if (!(in.fq[i] > 0 && in.fq[i] <= 1)) return bad("`freq` must hold frequencies above 0 and at most 1.");
    if (in.deg && !(in.fsz[i] > 0 && std::isfinite(in.fsz[i]))) return bad("`freq$size` must be sizes above 0.");
    for (std::size_t j = 0; j < i; ++j) {
      if (in.floc.s[j] == in.floc.s[i] && in.fal.s[j] == in.fal.s[i]) {
        return bad("`freq` lists the allele " + in.fal.s[i] + " twice at locus " + in.floc.s[i] + ".");
      }
    }
    if (std::find(loc.begin(), loc.end(), in.floc.s[i]) == loc.end()) loc.push_back(in.floc.s[i]);
  }
  auto inf = [&in, nf](const std::string &l, const std::string &a) {
    for (std::size_t i = 0; i < nf; ++i) {
      if (in.floc.s[i] == l && in.fal.s[i] == a) return static_cast<int>(i);
    }
    return -1;
  };
  // peaks
  const std::size_t npk = in.ploc.s.size();
  bool anyp = false;
  for (std::size_t i = 0; i < npk; ++i) {
    if (in.ploc.na[i] || in.pal.na[i] || !std::isfinite(in.ph[i])) return bad("`peaks` has a missing locus, allele or height.");
    if (in.ph[i] < in.at) return bad("`peaks` has a height below `at`: every peak must reach the threshold.");
    if (inf(in.ploc.s[i], in.pal.s[i]) < 0) {
      return bad("`peaks` has the allele " + in.pal.s[i] + " at locus " + in.ploc.s[i] + ", which `freq` does not list.");
    }
    for (std::size_t j = 0; j < i; ++j) {
      if (in.ploc.s[j] == in.ploc.s[i] && in.pal.s[j] == in.pal.s[i]) {
        return bad("`peaks` gives the allele " + in.pal.s[i] + " twice at locus " + in.ploc.s[i] + ".");
      }
    }
    anyp = true;
  }
  if (!anyp) return bad("`peaks` must have a peak.");
  // reference profiles
  const std::size_t nr = in.rid.s.size();
  std::vector<std::string> ids;
  for (std::size_t i = 0; i < nr; ++i) {
    if (in.rid.na[i] || in.rloc.na[i] || in.ra1.na[i] || in.ra2.na[i]) {
      return bad("`ref` has a missing id, locus or allele.");
    }
    if (inf(in.rloc.s[i], in.ra1.s[i]) < 0 || inf(in.rloc.s[i], in.ra2.s[i]) < 0) {
      return bad("`ref` has an allele of " + in.rid.s[i] + " at locus " + in.rloc.s[i] + " that `freq` does not list.");
    }
    for (std::size_t j = 0; j < i; ++j) {
      if (in.rid.s[j] == in.rid.s[i] && in.rloc.s[j] == in.rloc.s[i]) {
        return bad("`ref` gives " + in.rid.s[i] + " twice at locus " + in.rloc.s[i] + ".");
      }
    }
    if (std::find(ids.begin(), ids.end(), in.rid.s[i]) == ids.end()) ids.push_back(in.rid.s[i]);
  }
  for (int h = 0; h < 2; ++h) {
    const std::vector<std::string> &w = h ? in.hd : in.hp;
    const char *an = h ? "hd" : "hp";
    for (std::size_t i = 0; i < w.size(); ++i) {
      if (std::find(ids.begin(), ids.end(), w[i]) == ids.end()) {
        return bad(std::string("`") + an + "` names " + w[i] + ", who has no profile in `ref`.");
      }
      if (std::find(w.begin(), w.begin() + static_cast<std::ptrdiff_t>(i), w[i]) != w.begin() + static_cast<std::ptrdiff_t>(i)) {
        return bad(std::string("`") + an + "` names " + w[i] + " twice.");
      }
    }
    if (static_cast<int>(w.size()) > nc) return bad(std::string("`") + an + "` names more contributors than `nc`.");
  }
  // the parameters given, if any
  std::vector<double> fixed;
  if (in.haspar) {
    std::vector<std::string> want = {"mu", "omega"};
    if (in.deg) want.push_back("beta");
    if (in.stut) want.push_back("xi");
    for (int k = 0; k < nc; ++k) want.push_back("m" + std::to_string(k + 1));
    std::string ws;
    for (const std::string &s : want) ws += (ws.empty() ? "" : ", ") + s;
    if (in.parnm.size() != want.size()) return bad("`par` must name " + ws + ".");
    fixed.assign(4 + static_cast<std::size_t>(nc), 0.0);
    fixed[2] = 1.0;
    double ms = 0.0;
    for (std::size_t j = 0; j < want.size(); ++j) {
      const std::size_t i = static_cast<std::size_t>(std::find(in.parnm.begin(), in.parnm.end(), want[j]) - in.parnm.begin());
      if (i == in.parnm.size()) return bad("`par` must name " + ws + ".");
      const double v = in.par[i];
      const std::string &k = want[j];
      if (k == "mu" || k == "omega") {
        if (!(v > 0) || !std::isfinite(v)) return bad("`par`: mu and omega must be above 0.");
        fixed[k == "mu" ? 0 : 1] = v;
      } else if (k == "beta") {
        if (!(v > 0 && v <= 1)) return bad("`par`: beta must be above 0 and at most 1.");
        fixed[2] = v;
      } else if (k == "xi") {
        if (!(v >= 0 && v < 1)) return bad("`par`: xi must be from 0 and below 1.");
        fixed[3] = v;
      } else {
        if (!(v >= 0 && v <= 1)) return bad("`par`: the mixture proportions must be from 0 to 1.");
        fixed[4 + std::stoul(k.substr(1)) - 1] = v;
        ms += v;
      }
    }
    if (std::fabs(ms - 1) > 1e-9) return bad("`par`: the mixture proportions must sum to 1.");
  }
  // positions locus by locus
  Hyp hy[2];
  for (int h = 0; h < 2; ++h) {
    hy[h].nk = static_cast<int>((h ? in.hd : in.hp).size());
    hy[h].nu = nc - hy[h].nk;
    hy[h].who = h ? in.hd : in.hp;
    for (int k = 0; k < hy[h].nu; ++k) hy[h].who.push_back("U" + std::to_string(k + 1));
  }
  double hsum = 0.0;
  int hn = 0;
  for (const std::string &l : loc) {
    std::vector<int> fi;  // the rows of freq of this locus
    for (std::size_t i = 0; i < nf; ++i) {
      if (in.floc.s[i] == l) fi.push_back(static_cast<int>(i));
    }
    auto num = [&in](int i, double &v) {
      const std::string &a = in.fal.s[static_cast<std::size_t>(i)];
      char *e = nullptr;
      v = std::strtod(a.c_str(), &e);
      return e != a.c_str() && *e == '\0' && std::isfinite(v);
    };
    std::vector<char> keep(fi.size(), 0);
    double ps = 0.0, lh = 0.0;
    for (std::size_t j = 0; j < fi.size(); ++j) {
      const std::size_t i = static_cast<std::size_t>(fi[j]);
      ps += in.fq[i];
      for (std::size_t q = 0; q < npk; ++q) {
        if (in.ploc.s[q] == l && in.pal.s[q] == in.fal.s[i]) {
          keep[j] = 1;
          lh += in.ph[q];
        }
      }
      for (std::size_t q = 0; q < nr; ++q) {
        if (in.rloc.s[q] == l && (in.ra1.s[q] == in.fal.s[i] || in.ra2.s[q] == in.fal.s[i])) keep[j] = 1;
      }
    }
    if (ps > 1 + 1e-9) return bad("the frequencies at locus " + l + " sum to more than 1.");
    if (lh > 0) {
      hsum += lh;
      ++hn;
    }
    std::vector<double> lab(fi.size(), 0.0);
    if (in.stut) {
      for (std::size_t j = 0; j < fi.size(); ++j) {
        if (!num(fi[j], lab[j])) {
          return bad("the stutter model needs numeric alleles; locus " + l + " has " + in.fal.s[static_cast<std::size_t>(fi[j])] + ".");
        }
      }
      // the parents of the peaks seen, then the stutter targets of all
      // kept, as alleles of their own
      std::vector<char> k0 = keep;
      for (std::size_t j = 0; j < fi.size(); ++j) {
        for (std::size_t q = 0; q < fi.size(); ++q) {
          if (k0[q] && std::fabs(lab[j] - 1 - lab[q]) < 1e-9) keep[j] = 1;
        }
      }
      k0 = keep;
      for (std::size_t j = 0; j < fi.size(); ++j) {
        for (std::size_t q = 0; q < fi.size(); ++q) {
          if (k0[q] && std::fabs(lab[q] - 1 - lab[j]) < 1e-9) keep[j] = 1;
        }
      }
    }
    // positions: the alleles kept, the rest as one allele Q, then the
    // stutter targets that are no allele
    std::vector<std::string> pl;
    std::vector<double> pp, pz, pn;
    double rest = ps, zs = 0.0;
    int nz = 0;
    for (std::size_t j = 0; j < fi.size(); ++j) {
      const std::size_t i = static_cast<std::size_t>(fi[j]);
      if (in.deg) {
        zs += in.fsz[i];
        ++nz;
      }
      if (!keep[j]) continue;
      pl.push_back(in.fal.s[i]);
      pp.push_back(in.fq[i]);
      pz.push_back(in.deg ? in.fsz[i] : 0.0);
      pn.push_back(lab[j]);
      rest -= in.fq[i];
    }
    const std::size_t nkept = pl.size();
    const double qf = 1 - (ps - rest);
    const bool q = qf > 1e-12;
    if (q) {
      pl.push_back("");
      pp.push_back(qf);
      pz.push_back(nz ? zs / nz : 0.0);
      pn.push_back(NAN);
    }
    std::vector<int> st(pl.size(), -1);
    if (in.stut) {
      for (std::size_t j = 0; j < nkept; ++j) {
        std::size_t t = 0;
        while (t < nkept && std::fabs(pn[t] - (pn[j] - 1)) >= 1e-9) ++t;
        if (t < nkept) {
          st[j] = static_cast<int>(t);
          continue;
        }
        // a target that is no allele: a position no one carries (its size
        // is not read: the stutter it receives is degraded as its parent)
        pl.push_back("");
        pp.push_back(0.0);
        pz.push_back(pz[j]);
        pn.push_back(pn[j] - 1);
        st.push_back(-1);
        st[j] = static_cast<int>(pl.size()) - 1;
      }
    }
    // heights, conditioning counts and the known contributors' codes
    const std::size_t n = pl.size();
    std::vector<double> yy(n, 0.0), cnt(n, 0.0);
    for (std::size_t q2 = 0; q2 < npk; ++q2) {
      if (in.ploc.s[q2] != l) continue;
      for (std::size_t j = 0; j < nkept; ++j) {
        if (pl[j] == in.pal.s[q2]) yy[j] = in.ph[q2];
      }
    }
    double tc = 0.0;
    auto code = [&pl, nkept](const std::string &a) {
      for (std::size_t j = 0; j < nkept; ++j) {
        if (pl[j] == a) return static_cast<int>(j);
      }
      return -1;
    };
    for (std::size_t r = 0; r < nr; ++r) {
      if (in.rloc.s[r] != l) continue;
      cnt[static_cast<std::size_t>(code(in.ra1.s[r]))] += 1;
      cnt[static_cast<std::size_t>(code(in.ra2.s[r]))] += 1;
      tc += 2;
    }
    for (int h = 0; h < 2; ++h) {
      for (int k = 0; k < hy[h].nk; ++k) {
        const std::string &w = hy[h].who[static_cast<std::size_t>(k)];
        std::size_t r = 0;
        while (r < nr && !(in.rid.s[r] == w && in.rloc.s[r] == l)) ++r;
        if (r == nr) return bad("`ref` gives no genotype for " + w + " at locus " + l + ".");
        hy[h].kg.push_back(code(in.ra1.s[r]));
        hy[h].kg.push_back(code(in.ra2.s[r]));
      }
    }
    // the unknowns' genotype combinations at this locus
    double ng = 0.0;
    for (std::size_t j = 0; j < n; ++j) ng += pp[j] > 0;
    ng = ng * (ng + 1) / 2;
    for (int h = 0; h < 2; ++h) {
      if (std::pow(ng, hy[h].nu) > 2e6) {
        return bad("too many genotypes of the unknown contributors at locus " + l + " (more than 2e6 combinations).");
      }
    }
    ly.np.push_back(static_cast<int>(n));
    ly.p.insert(ly.p.end(), pp.begin(), pp.end());
    ly.y.insert(ly.y.end(), yy.begin(), yy.end());
    ly.sz.insert(ly.sz.end(), pz.begin(), pz.end());
    ly.st.insert(ly.st.end(), st.begin(), st.end());
    ly.cn.insert(ly.cn.end(), cnt.begin(), cnt.end());
    ly.cl.push_back(tc);
    out.loc.push_back(l);
  }
  ly.nl = static_cast<int>(out.loc.size());
  // each hypothesis: the likelihood at par, or its maximum
  std::vector<double> llh[2];
  for (int h = 0; h < 2; ++h) {
    llh[h].assign(static_cast<std::size_t>(ly.nl), 0.0);
    Ctx c{&ly, &hy[h], &llh[h], CLAB_OK};
    std::vector<double> par;
    if (in.haspar) {
      par = fixed;
      const int rc = clab_gmxlik(ly.nl, ly.np.data(), ly.p.data(), ly.y.data(), ly.deg ? ly.sz.data() : nullptr,
                                 ly.stut ? ly.st.data() : nullptr, hy[h].nk, hy[h].kg.empty() ? nullptr : hy[h].kg.data(),
                                 hy[h].nu, ly.cn.data(), ly.cl.data(), ly.theta, ly.at, ly.pc, ly.lambda, ly.deg, ly.stut,
                                 par.data(), llh[h].data());
      if (rc != CLAB_OK) {
        err = clab_strerror(rc);
        return rc;
      }
      out.conv[h] = -1;
      out.iter[h] = 0;
    } else {
      const int p = npar(ly, nc);
      std::vector<double> best, u(static_cast<std::size_t>(p)), work(static_cast<std::size_t>(p * p + 5 * p));
      double bf = HUGE_VAL;
      const int nst = nc == 1 ? 1 : 3;
      for (int s = 0; s < nst; ++s) {
        u[0] = std::log(hsum / hn / 2);
        u[1] = std::log(0.3);
        int j = 2;
        if (ly.deg) u[static_cast<std::size_t>(j++)] = std::log(0.8 / 0.2);
        if (ly.stut) u[static_cast<std::size_t>(j++)] = std::log(0.05 / 0.95);
        for (int k = 0; k < nc - 1; ++k) {
          // equal proportions, then falling and rising ones
          const double w = s == 0 ? 0.0 : std::log(static_cast<double>(s == 1 ? nc - k : k + 1));
          u[static_cast<std::size_t>(j + k)] = w;
        }
        double fv = 0.0;
        int it = 0;
        const int cv = qn::min(p, u.data(), fobj, &c, 500, 1e-6, &fv, &it, work.data(),
                               static_cast<double>(p) * 1000.0);
        if (c.rc == CLAB_EINTR || kutil::stop(0.0)) {
          err = clab_strerror(CLAB_EINTR);
          return CLAB_EINTR;
        }
        if (cv >= 0 && fv < bf) {
          bf = fv;
          best = u;
          out.conv[h] = cv;
          out.iter[h] = it;
        }
      }
      if (best.empty()) {
        err = "the likelihood under " + std::string(h ? "hd" : "hp") + " is 0 at every start: the account cannot explain the peaks.";
        return CLAB_EDOM;
      }
      par.assign(4 + static_cast<std::size_t>(nc), 0.0);
      unpack(ly, nc, best.data(), par.data());
      nll(c, best.data());
    }
    double s = 0.0;
    for (double v : llh[h]) s += v;
    out.loglik[h] = s;
    out.mu[h] = par[0];
    out.omega[h] = par[1];
    out.beta[h] = ly.deg ? par[2] : NAN;
    out.xi[h] = ly.stut ? par[3] : NAN;
    out.who[h] = hy[h].who;
    out.m[h].assign(par.begin() + 4, par.end());
  }
  const double ln10 = std::log(10.0);
  for (int l = 0; l < ly.nl; ++l) out.llr.push_back((llh[0][static_cast<std::size_t>(l)] - llh[1][static_cast<std::size_t>(l)]) / ln10);
  out.log10lr = (out.loglik[0] - out.loglik[1]) / ln10;
  return CLAB_OK;
}

}  // namespace frx
