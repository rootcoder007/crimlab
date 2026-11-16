// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, forensics: pedigree likelihoods at autosomal markers with
// mutation, and the likelihood ratio of two pedigrees (pedlr()). The
// likelihood of a pedigree is summed over the alleles every untyped
// person may carry by variable elimination on the alleles, two variables
// per person (the allele from the father and from the mother): founders'
// alleles are drawn from the population frequencies, a child's allele
// from one of the two alleles of the parent with probability 1/2 each,
// then mutated by the parent's mutation matrix (Elston and Stewart 1971
// give the model; Egeland, Kling and Mostad 2016 the mutation models).
// Elimination works on any pedigree, with loops and inbreeding; its cost
// is the largest table it builds, which is capped.
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include "crimlab.h"
#include "frx.h"
#include "kutil.h"

namespace {

// the largest table elimination may build, in cells
constexpr double maxcell = 67108864.0;

// A factor: a table over vars (ascending), the last varying fastest.
struct Fac {
  std::vector<int> v;
  std::vector<double> x;
};

struct Elim {
  std::vector<std::vector<int>> dom;  // the allele codes each variable may take
  std::vector<Fac> f;
  std::vector<char> live;
  double lsc = 0.0;  // the logarithm of the scale divided out
};

// the factor over vars (any order) whose value at an assignment is val(a),
// a[k] the allele code of vars[k]
template <class F>
Fac mkfac(const Elim &e, std::vector<int> vars, F val) {
  Fac r;
  std::vector<int> s = vars;
  std::sort(s.begin(), s.end());
  r.v = s;
  std::size_t tot = 1;
  for (int v : s) tot *= e.dom[static_cast<std::size_t>(v)].size();
  r.x.assign(tot, 0.0);
  std::vector<std::size_t> dg(s.size(), 0);
  std::vector<int> a(vars.size());
  for (std::size_t t = 0; t < tot; ++t) {
    for (std::size_t k = 0; k < vars.size(); ++k) {
      const std::size_t j = static_cast<std::size_t>(std::find(s.begin(), s.end(), vars[k]) - s.begin());
      a[k] = e.dom[static_cast<std::size_t>(s[j])][dg[j]];
    }
    r.x[t] = val(a);
    for (std::size_t k = s.size(); k-- > 0;) {
      if (++dg[k] < e.dom[static_cast<std::size_t>(s[k])].size()) break;
      dg[k] = 0;
    }
  }
  return r;
}

// the cells of the largest first table clab_pedlik() builds, and of all of
// them with the allele lists, counted before any is made: a child's
// transmission from a parent is a table over the child's allele and the
// parent's two; d(i) is the number of alleles person i may carry
void cells0(int n, const int *fa, const int *mo, const int *g, int na, double &big, double &all) {
  auto d = [&](int i) { return g[2 * i] < 0 ? static_cast<double>(na) : g[2 * i] == g[2 * i + 1] ? 1.0 : 2.0; };
  big = all = 0.0;
  for (int i = 0; i < n; ++i) {
    const double di = d(i);
    all += 4 * di + 4;
    for (int par : {fa[i], mo[i]}) {
      const double t = par < 0 ? di : di * d(par) * d(par);
      big = std::max(big, t);
      all += t;
    }
  }
}

// the union of the scopes of the live factors that hold var v
std::vector<int> scope(const Elim &e, int v, std::vector<std::size_t> &hold) {
  std::vector<int> u;
  hold.clear();
  for (std::size_t i = 0; i < e.f.size(); ++i) {
    if (!e.live[i] || !std::binary_search(e.f[i].v.begin(), e.f[i].v.end(), v)) continue;
    hold.push_back(i);
    for (int w : e.f[i].v) u.push_back(w);
  }
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  return u;
}

// multiplies the factors holding v and sums v out; CLAB_EINTR when asked
// to stop
int elim1(Elim &e, int v, const std::vector<int> &u, const std::vector<std::size_t> &hold) {
  const std::size_t m = u.size();
  std::vector<std::size_t> sz(m);
  for (std::size_t k = 0; k < m; ++k) sz[k] = e.dom[static_cast<std::size_t>(u[k])].size();
  // strides of each factor and of the result over the positions of u
  const std::size_t nf = hold.size();
  std::vector<std::vector<std::size_t>> st(nf + 1, std::vector<std::size_t>(m, 0));
  for (std::size_t h = 0; h <= nf; ++h) {
    std::vector<int> fv;
    if (h < nf) {
      fv = e.f[hold[h]].v;
    } else {
      for (int w : u) {
        if (w != v) fv.push_back(w);
      }
    }
    std::size_t s = 1;
    for (std::size_t k = fv.size(); k-- > 0;) {
      const std::size_t j = static_cast<std::size_t>(std::find(u.begin(), u.end(), fv[k]) - u.begin());
      st[h][j] = s;
      s *= sz[j];
    }
  }
  Fac r;
  std::size_t rn = 1;
  for (std::size_t k = 0; k < m; ++k) {
    if (u[k] != v) {
      r.v.push_back(u[k]);
      rn *= sz[k];
    }
  }
  r.x.assign(rn, 0.0);
  std::size_t tot = 1;
  for (std::size_t k = 0; k < m; ++k) tot *= sz[k];
  std::vector<std::size_t> dg(m, 0), ix(nf + 1, 0);
  for (std::size_t t = 0; t < tot; ++t) {
    if ((t & 4095u) == 4095u && kutil::stop(4096.0 * static_cast<double>(nf + 1))) return CLAB_EINTR;
    double p = 1.0;
    for (std::size_t h = 0; h < nf && p != 0.0; ++h) p *= e.f[hold[h]].x[ix[h]];
    r.x[ix[nf]] += p;
    for (std::size_t k = m; k-- > 0;) {
      if (++dg[k] < sz[k]) {
        for (std::size_t h = 0; h <= nf; ++h) ix[h] += st[h][k];
        break;
      }
      dg[k] = 0;
      for (std::size_t h = 0; h <= nf; ++h) ix[h] -= st[h][k] * (sz[k] - 1);
    }
  }
  double mx = 0.0;
  for (double x : r.x) mx = std::max(mx, x);
  if (mx > 0.0) {
    for (double &x : r.x) x /= mx;
    e.lsc += std::log(mx);
  }
  for (std::size_t h : hold) e.live[h] = 0;
  e.f.push_back(r);
  e.live.push_back(1);
  return CLAB_OK;
}

}  // namespace

extern "C" {

int clab_mutmat(int model, int na, const double *lab, double rate, double range, double *m) CLAB_TRY {
  if (na < 1 || na > 46340 || model < 0 || model > 2) return CLAB_EARG;
  if (!(rate >= 0.0 && rate <= 1.0) || (model == 2 && !(range > 0.0 && range < 1.0))) return CLAB_EDOM;
  const std::size_t n = static_cast<std::size_t>(na);
  for (std::size_t i = 0; i < n; ++i) {
    double tot = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
      if (j == i) continue;
      if (model == 2) {
        if (!std::isfinite(lab[i]) || !std::isfinite(lab[j])) return CLAB_EDOM;
        tot += std::pow(range, std::fabs(lab[i] - lab[j]));
      }
    }
    for (std::size_t j = 0; j < n; ++j) {
      double v = 0.0;
      if (model == 0 || na == 1) {
        v = i == j ? 1.0 : 0.0;
      } else if (i == j) {
        v = 1.0 - rate;
      } else if (model == 1) {
        v = rate / (na - 1);
      } else {
        v = rate * std::pow(range, std::fabs(lab[i] - lab[j])) / tot;
      }
      m[i * n + j] = v;
    }
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_pedlik(int n, const int *fa, const int *mo, const int *g, int na, const double *p,
                const double *mm, const double *mf, double *ll) CLAB_TRY {
  if (n < 1 || n > INT_MAX / 2 || na < 1 || na > 46340) return CLAB_EARG;
  double ps = 0.0;
  for (int a = 0; a < na; ++a) {
    if (!(p[a] > 0.0 && p[a] <= 1.0)) return CLAB_EDOM;
    ps += p[a];
  }
  if (std::fabs(ps - 1.0) > 1e-9) return CLAB_EDOM;
  for (const double *mt : {mm, mf}) {
    if (mt == nullptr) continue;
    for (int a = 0; a < na; ++a) {
      double s = 0.0;
      for (int b = 0; b < na; ++b) {
        const double v = mt[static_cast<std::size_t>(a) * static_cast<std::size_t>(na) + static_cast<std::size_t>(b)];
        if (!(v >= 0.0 && v <= 1.0)) return CLAB_EDOM;
        s += v;
      }
      if (std::fabs(s - 1.0) > 1e-9) return CLAB_EDOM;
    }
  }
  // parents: both or neither, in range, no one their own ancestor
  std::vector<char> done(static_cast<std::size_t>(n), 0);
  for (int i = 0; i < n; ++i) {
    if ((fa[i] < 0) != (mo[i] < 0) || fa[i] >= n || mo[i] >= n || fa[i] == i || mo[i] == i ||
        (fa[i] >= 0 && fa[i] == mo[i])) {
      return CLAB_EARG;
    }
    if (g[2 * i] < -1 || g[2 * i] >= na || g[2 * i + 1] < -1 || g[2 * i + 1] >= na ||
        (g[2 * i] < 0) != (g[2 * i + 1] < 0)) {
      return CLAB_EARG;
    }
  }
  for (int left = n; left > 0;) {
    int got = 0;
    for (int i = 0; i < n; ++i) {
      if (done[static_cast<std::size_t>(i)]) continue;
      if (fa[i] < 0 || (done[static_cast<std::size_t>(fa[i])] && done[static_cast<std::size_t>(mo[i])])) {
        done[static_cast<std::size_t>(i)] = 1;
        ++got;
      }
    }
    if (got == 0) return CLAB_EARG;
    left -= got;
  }
  // the first tables within the bound of elimination, which would refuse
  // any of them anyway, before they are built
  double big = 0.0, all = 0.0;
  cells0(n, fa, mo, g, na, big, all);
  if (big > maxcell || all > 4 * maxcell) return CLAB_EARG;
  Elim e;
  const std::size_t nv = 2 * static_cast<std::size_t>(n);
  e.dom.resize(nv);
  for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i) {
    const int a = g[2 * i], b = g[2 * i + 1];
    for (std::size_t k = 0; k < 2; ++k) {
      std::vector<int> &d = e.dom[2 * i + k];
      if (a < 0) {
        for (int c = 0; c < na; ++c) d.push_back(c);
      } else {
        d.push_back(std::min(a, b));
        if (a != b) d.push_back(std::max(a, b));
      }
    }
  }
  const std::size_t sn = static_cast<std::size_t>(na);
  for (int i = 0; i < n; ++i) {
    const int pv = 2 * i, mv = 2 * i + 1;
    if (g[2 * i] >= 0 && g[2 * i] != g[2 * i + 1]) {
      e.f.push_back(mkfac(e, {pv, mv}, [](const std::vector<int> &a) { return a[0] != a[1] ? 1.0 : 0.0; }));
    }
    if (fa[i] < 0) {
      for (int v : {pv, mv}) {
        e.f.push_back(mkfac(e, {v}, [p](const std::vector<int> &a) { return p[a[0]]; }));
      }
      continue;
    }
    for (int s = 0; s < 2; ++s) {
      const int par = s == 0 ? fa[i] : mo[i];
      const double *mt = s == 0 ? mm : mf;
      e.f.push_back(mkfac(e, {s == 0 ? pv : mv, 2 * par, 2 * par + 1}, [mt, sn](const std::vector<int> &a) {
        const std::size_t x = static_cast<std::size_t>(a[0]), u = static_cast<std::size_t>(a[1]),
                          w = static_cast<std::size_t>(a[2]);
        if (mt == nullptr) return 0.5 * ((u == x) + (w == x));
        return 0.5 * (mt[u * sn + x] + mt[w * sn + x]);
      }));
    }
  }
  e.live.assign(e.f.size(), 1);
  std::vector<char> gone(nv, 0);
  std::vector<std::size_t> hold;
  for (std::size_t step = 0; step < nv; ++step) {
    // the variable whose elimination builds the smallest table
    int best = -1;
    double bc = 0.0;
    for (std::size_t v = 0; v < nv; ++v) {
      if (gone[v]) continue;
      const std::vector<int> u = scope(e, static_cast<int>(v), hold);
      double c = 1.0;
      for (int w : u) c *= static_cast<double>(e.dom[static_cast<std::size_t>(w)].size());
      if (best < 0 || c < bc) {
        best = static_cast<int>(v);
        bc = c;
      }
    }
    if (bc > maxcell) return CLAB_EARG;
    const std::vector<int> u = scope(e, best, hold);
    const int rc = elim1(e, best, u, hold);
    if (rc != CLAB_OK) return rc;
    gone[static_cast<std::size_t>(best)] = 1;
  }
  double pr = 1.0;
  for (std::size_t i = 0; i < e.f.size(); ++i) {
    if (e.live[i]) pr *= e.f[i].x[0];
  }
  *ll = pr > 0.0 ? std::log(pr) + e.lsc : -HUGE_VAL;
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

}  // extern "C"

namespace frx {

namespace {

bool nolab(const Lab &x, std::size_t i) {
  return x.na[i] || x.s[i].empty() || x.s[i] == "0";
}

// A pedigree table: parents as row numbers (-1 for a founder)
bool ptab(const std::string &an, const Lab &id, const Lab &fid, const Lab &mid, std::vector<int> &fa,
          std::vector<int> &mo, std::string &err) {
  const std::size_t n = id.s.size();
  fa.assign(n, -1);
  mo.assign(n, -1);
  if (n == 0) {
    err = "`" + an + "` must have a row.";
    return false;
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (id.na[i] || id.s[i].empty()) {
      err = "`" + an + "` has an id that is missing.";
      return false;
    }
    for (std::size_t j = 0; j < i; ++j) {
      if (id.s[j] == id.s[i]) {
        err = "`" + an + "` has the id " + id.s[i] + " twice.";
        return false;
      }
    }
  }
  auto find = [&id](const std::string &s) {
    return static_cast<int>(std::find(id.s.begin(), id.s.end(), s) - id.s.begin());
  };
  for (std::size_t i = 0; i < n; ++i) {
    const bool nf = nolab(fid, i), nm = nolab(mid, i);
    if (nf != nm) {
      err = "`" + an + "`: " + id.s[i] + " has one parent; give both or neither.";
      return false;
    }
    if (nf) continue;
    fa[i] = find(fid.s[i]);
    mo[i] = find(mid.s[i]);
    for (int k = 0; k < 2; ++k) {
      if ((k ? mo[i] : fa[i]) == static_cast<int>(n)) {
        err = "`" + an + "`: the " + (k ? "mother " + mid.s[i] : "father " + fid.s[i]) + " of " + id.s[i] +
              " is not in `" + an + "`.";
        return false;
      }
    }
    if (fa[i] == mo[i]) {
      err = "`" + an + "`: " + id.s[i] + " has one person as father and mother.";
      return false;
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      if (fa[i] >= 0 && mo[j] == fa[i]) {
        err = "`" + an + "`: " + id.s[static_cast<std::size_t>(fa[i])] + " is both a father and a mother.";
        return false;
      }
    }
  }
  // no one is their own ancestor
  std::vector<char> done(n, 0);
  for (std::size_t left = n; left > 0;) {
    std::size_t got = 0;
    for (std::size_t i = 0; i < n; ++i) {
      if (!done[i] && (fa[i] < 0 || (done[static_cast<std::size_t>(fa[i])] && done[static_cast<std::size_t>(mo[i])]))) {
        done[i] = 1;
        ++got;
      }
    }
    if (got == 0) {
      for (std::size_t i = 0; i < n; ++i) {
        if (!done[i]) {
          err = "`" + an + "` has a cycle: " + id.s[i] + " is among their own ancestors.";
          return false;
        }
      }
    }
    left -= got;
  }
  return true;
}

}  // namespace

int pedlr(const PedIn &in, PedOut &out, std::string &err) {
  out = PedOut();
  int model = -1;
  static const char *const mods[] = {"none", "equal", "step"};
  for (int k = 0; k < 3; ++k) {
    if (in.mut == mods[k]) model = k;
  }
  if (model < 0) {
    err = "`mut` must be one of none, equal, step.";
    return CLAB_EARG;
  }
  if ((in.rate.size() != 1 && in.rate.size() != 2) ||
      std::any_of(in.rate.begin(), in.rate.end(), [](double r) { return !(r >= 0.0 && r <= 1.0); })) {
    err = "`rate` must be one or two numbers from 0 to 1.";
    return CLAB_EARG;
  }
  if (!(in.range > 0.0 && in.range < 1.0)) {
    err = "`range` must be one number strictly between 0 and 1.";
    return CLAB_EARG;
  }
  std::vector<int> f1, m1, f2, m2;
  if (!ptab("ped", in.id, in.fid, in.mid, f1, m1, err)) return CLAB_EARG;
  // the genotypes: one row per typed person and locus
  const std::size_t ng = in.gid.s.size();
  std::vector<std::string> typed;
  for (std::size_t i = 0; i < ng; ++i) {
    if (in.gid.na[i] || in.gloc.na[i]) {
      err = "`geno` has a missing id or locus.";
      return CLAB_EARG;
    }
    if (in.ga1.na[i] != in.ga2.na[i]) {
      err = "`geno` gives one allele of two for " + in.gid.s[i] + " at locus " + in.gloc.s[i] + ".";
      return CLAB_EARG;
    }
    if (std::find(in.id.s.begin(), in.id.s.end(), in.gid.s[i]) == in.id.s.end()) {
      err = "`geno`: the id " + in.gid.s[i] + " is not in `ped`.";
      return CLAB_EARG;
    }
    for (std::size_t j = 0; j < i; ++j) {
      if (in.gid.s[j] == in.gid.s[i] && in.gloc.s[j] == in.gloc.s[i]) {
        err = "`geno` gives " + in.gid.s[i] + " twice at locus " + in.gloc.s[i] + ".";
        return CLAB_EARG;
      }
    }
    if (!in.ga1.na[i] && std::find(typed.begin(), typed.end(), in.gid.s[i]) == typed.end()) {
      typed.push_back(in.gid.s[i]);
    }
  }
  Lab aid, afid, amid;
  if (in.hasalt) {
    aid = in.aid;
    afid = in.afid;
    amid = in.amid;
  } else {
    // everyone typed, unrelated
    for (const std::string &s : typed) {
      aid.s.push_back(s);
      aid.na.push_back(0);
    }
    afid.s.assign(typed.size(), std::string());
    afid.na.assign(typed.size(), 1);
    amid = afid;
  }
  if (!ptab("alt", aid, afid, amid, f2, m2, err)) {
    if (in.hasalt) return CLAB_EARG;
    err = "`geno` must type someone.";
    return CLAB_EARG;
  }
  for (const std::string &s : typed) {
    if (std::find(aid.s.begin(), aid.s.end(), s) == aid.s.end()) {
      err = "`geno`: the id " + s + " is not in `alt`.";
      return CLAB_EARG;
    }
  }
  // the frequencies, locus by locus in their order
  const std::size_t nf = in.floc.s.size();
  std::vector<std::string> loc;
  for (std::size_t i = 0; i < nf; ++i) {
    if (in.floc.na[i] || in.fal.na[i]) {
      err = "`freq` has a missing locus or allele.";
      return CLAB_EARG;
    }
    if (!(in.fq[i] > 0.0 && in.fq[i] <= 1.0)) {
      err = "`freq` must hold frequencies above 0 and at most 1.";
      return CLAB_EARG;
    }
    for (std::size_t j = 0; j < i; ++j) {
      if (in.floc.s[j] == in.floc.s[i] && in.fal.s[j] == in.fal.s[i]) {
        err = "`freq` lists the allele " + in.fal.s[i] + " twice at locus " + in.floc.s[i] + ".";
        return CLAB_EARG;
      }
    }
    if (std::find(loc.begin(), loc.end(), in.floc.s[i]) == loc.end()) loc.push_back(in.floc.s[i]);
  }
  for (std::size_t i = 0; i < ng; ++i) {
    if (std::find(loc.begin(), loc.end(), in.gloc.s[i]) == loc.end()) {
      err = "`geno` has the locus " + in.gloc.s[i] + ", which `freq` does not list.";
      return CLAB_EARG;
    }
  }
  const double rm = in.rate.back(), rf = in.rate[0];
  const double ln10 = std::log(10.0);
  double tot = 0.0;
  for (const std::string &l : loc) {
    std::vector<std::string> al;
    std::vector<double> p, num;
    double ps = 0.0;
    for (std::size_t i = 0; i < nf; ++i) {
      if (in.floc.s[i] != l) continue;
      al.push_back(in.fal.s[i]);
      p.push_back(in.fq[i]);
      ps += in.fq[i];
    }
    if (ps > 1.0 + 1e-9) {
      err = "the frequencies at locus " + l + " sum to more than 1.";
      return CLAB_EARG;
    }
    if (model == 2) {
      if (ps < 1.0 - 1e-6) {
        err = "the stepwise model needs frequencies that sum to 1 at locus " + l + ".";
        return CLAB_EARG;
      }
      for (const std::string &a : al) {
        char *e = nullptr;
        const double v = std::strtod(a.c_str(), &e);
        if (e == a.c_str() || *e != '\0' || v != std::floor(v) || !std::isfinite(v)) {
          err = "the stepwise model needs whole-number alleles; locus " + l + " has " + a + ".";
          return CLAB_EARG;
        }
        num.push_back(v);
      }
    } else if (ps < 1.0 - 1e-9) {
      // the alleles not listed, as one more
      p.push_back(1.0 - ps);
      ps = 1.0;
    }
    for (double &x : p) x /= ps;
    const int na = static_cast<int>(p.size());
    // genotype codes of each person of each pedigree
    std::vector<int> g1(2 * in.id.s.size(), -1), g2(2 * aid.s.size(), -1);
    bool any = false;
    for (std::size_t i = 0; i < ng; ++i) {
      if (in.gloc.s[i] != l || in.ga1.na[i]) continue;
      int c[2];
      for (int k = 0; k < 2; ++k) {
        const std::string &a = k ? in.ga2.s[i] : in.ga1.s[i];
        const std::size_t j = static_cast<std::size_t>(std::find(al.begin(), al.end(), a) - al.begin());
        if (j == al.size()) {
          err = "`geno` has the allele " + a + " at locus " + l + ", which `freq` does not list.";
          return CLAB_EARG;
        }
        c[k] = static_cast<int>(j);
      }
      const std::size_t r1 = static_cast<std::size_t>(std::find(in.id.s.begin(), in.id.s.end(), in.gid.s[i]) - in.id.s.begin());
      const std::size_t r2 = static_cast<std::size_t>(std::find(aid.s.begin(), aid.s.end(), in.gid.s[i]) - aid.s.begin());
      g1[2 * r1] = g2[2 * r2] = c[0];
      g1[2 * r1 + 1] = g2[2 * r2 + 1] = c[1];
      any = true;
    }
    if (!any) continue;
    // the mutation matrices, na x na each, only with mutation and within
    // the bound of the tables
    std::vector<double> mmat, fmat;
    int rc = model > 0 && static_cast<double>(na) * na > maxcell ? CLAB_EARG : CLAB_OK;
    if (model > 0 && rc == CLAB_OK) {
      mmat.resize(static_cast<std::size_t>(na) * static_cast<std::size_t>(na));
      fmat.resize(mmat.size());
      rc = clab_mutmat(model, na, num.empty() ? nullptr : num.data(), rm, in.range, mmat.data());
      if (rc == CLAB_OK) rc = clab_mutmat(model, na, num.empty() ? nullptr : num.data(), rf, in.range, fmat.data());
    }
    double l1 = 0.0, l2 = 0.0;
    const double *pm = model > 0 ? mmat.data() : nullptr, *pf = model > 0 ? fmat.data() : nullptr;
    if (rc == CLAB_OK) {
      rc = clab_pedlik(static_cast<int>(in.id.s.size()), f1.data(), m1.data(), g1.data(), na, p.data(), pm, pf, &l1);
    }
    if (rc == CLAB_OK) {
      rc = clab_pedlik(static_cast<int>(aid.s.size()), f2.data(), m2.data(), g2.data(), na, p.data(), pm, pf, &l2);
    }
    if (rc == CLAB_EARG) {
      err = "the pedigrees are too large for exact computation at locus " + l + ".";
      return rc;
    }
    if (rc != CLAB_OK) {
      err = clab_strerror(rc);
      return rc;
    }
    out.loc.push_back(l);
    out.l1.push_back(l1 / ln10);
    out.l2.push_back(l2 / ln10);
    out.lr.push_back(std::exp(l1 - l2));
    tot += (l1 - l2) / ln10;
  }
  if (out.loc.empty()) {
    err = "`geno` must type someone.";
    return CLAB_EARG;
  }
  out.log10lr = tot;
  return CLAB_OK;
}

}  // namespace frx
