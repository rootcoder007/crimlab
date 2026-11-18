// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The drift battery (R/drift.R, src/drf.cpp) for clab_jrun(). The Monte
// Carlo p-value of a sparse table is drawn by R's own generator, so here
// drfthomog() refuses that branch.

bool drfin(Ctx &c, const char *nm, std::vector<double> &v) {
  CK(pnum(c, nm, v));
  std::vector<double> k;
  for (double x : v) {
    if (std::isnan(x)) continue;
    if (std::isinf(x)) return c.fail("`%s` must be finite (no Inf).", nm);
    k.push_back(x);
  }
  if (k.empty()) return c.fail("`%s` needs at least one non-missing value, got 0.", nm);
  v = k;
  return true;
}

bool drftks(Ctx &c, Val &out) {
  std::vector<double> x, y;
  CK(drfin(c, "x", x));
  CK(drfin(c, "y", y));
  double r[3];
  clab_ks(x.data(), ilen(x), y.data(), ilen(y), r);
  static const char *const nm[] = {"statistic", "p_value", "neff"};
  out = nvec(3, nm, r);
  return true;
}

bool drftpsi(Ctx &c, Val &out) {
  std::vector<double> x, y;
  CK(drfin(c, "x", x));
  CK(drfin(c, "y", y));
  double b = c.has("bins") ? val1(c["bins"]) : 10, e = c.has("eps") ? val1(c["eps"]) : 1e-6;
  if (std::isnan(b) || b < 2 || b != std::floor(b)) return c.fail("`bins` must be one whole number of at least 2.");
  if (b > 1e6) return c.fail("`bins` must be at most 1,000,000.");
  if (std::isnan(e) || e <= 0 || e >= 1) return c.fail("`eps` must be a positive floor below one (0 < eps < 1).");
  std::vector<double> p(sat<std::size_t>(b) + 2), q(p.size());
  const int k = clab_psibin(x.data(), ilen(x), y.data(), ilen(y), sat<int>(b), p.data(), q.data());
  double r[2];
  clab_psi(p.data(), q.data(), k, e, r);
  static const char *const nm[] = {"psi", "jsdiv"};
  out = nvec(2, nm, r);
  return true;
}

// .drtab(): text tabulated (names sorted by their bytes), or named counts
bool drtab(Ctx &c, const char *nm, std::vector<std::string> &lv, std::vector<double> &n) {
  lv.clear();
  n.clear();
  const Val &o = c[nm];
  if (o.t == Val::OBJ && isnum(o)) {
    lv = o.k;
    n = dbl(o);
  } else {
    bool txt = o.t == Val::STR || o.t == Val::ARR;
    for (const Val &e : o.a) txt = txt && (e.t == Val::STR || e.t == Val::NUL);
    if (!txt) return c.fail("a numeric `%s` must be named.", nm);
    std::vector<std::string> s;
    for (const Val &e : o.t == Val::ARR ? o.a : std::vector<Val>{o}) {
      if (e.t == Val::STR) s.push_back(e.s);
    }
    std::sort(s.begin(), s.end(), [](const std::string &a, const std::string &b) {
      return std::strcmp(a.c_str(), b.c_str()) < 0;
    });
    for (const std::string &t : s) {
      if (lv.empty() || lv.back() != t) {
        lv.push_back(t);
        n.push_back(0);
      }
      n.back() += 1;
    }
  }
  for (double v : n) {
    if (!std::isfinite(v) || v < 0) return c.fail("`%s` must be finite, non-negative counts.", nm);
  }
  return true;
}

void dralign(std::vector<std::string> &all, const std::vector<std::string> &lv, const std::vector<double> &n,
             std::vector<double> &o) {
  o.assign(all.size(), 0);
  for (std::size_t j = 0; j < all.size(); ++j) {
    const auto it = std::find(lv.begin(), lv.end(), all[j]);
    if (it != lv.end()) o[j] = n[sat<std::size_t>(it - lv.begin())];
  }
}

std::vector<std::string> unionv(const std::vector<std::string> &a, const std::vector<std::string> &b) {
  std::vector<std::string> u;
  for (const std::vector<std::string> *v : {&a, &b}) {
    for (const std::string &s : *v) {
      if (std::find(u.begin(), u.end(), s) == u.end()) u.push_back(s);
    }
  }
  return u;
}

Val chout(double st, double df, double p, double lp) {
  const double r[4] = {st, df, p, lp};
  static const char *const nm[] = {"statistic", "df", "p_value", "logpvalue"};
  return nvec(4, nm, r);
}

bool drftchisq(Ctx &c, Val &out) {
  std::vector<std::string> lo, le;
  std::vector<double> no, ne, o, e;
  CK(drtab(c, "obs", lo, no));
  CK(drtab(c, "expd", le, ne));
  std::vector<std::string> lv = unionv(lo, le);
  dralign(lv, lo, no, o);
  dralign(lv, le, ne, e);
  double so = 0, se = 0;
  for (std::size_t j = 0; j < lv.size(); ++j) {
    so += o[j];
    se += e[j];
  }
  if (so <= 0) return c.fail("`obs` has no counts.");
  if (se <= 0) return c.fail("`expd` has no counts.");
  double r[3];
  clab_chig(o.data(), e.data(), sat<int>(lv.size()), r);
  if (r[2] > 0) {
    out = chout(HUGE_VAL, std::max(1.0, r[1]), 0, -HUGE_VAL);
  } else if (r[1] < 1) {
    out = chout(r[0], 0, NA, NA);
  } else {
    const double p = clab_pchisq(r[0], r[1], 0);
    out = chout(r[0], r[1], p, std::log(p));
  }
  return true;
}

bool drfthomog(Ctx &c, Val &out) {
  std::vector<std::string> lx, ly;
  std::vector<double> nx, ny, x, y;
  CK(drtab(c, "x", lx, nx));
  CK(drtab(c, "y", ly, ny));
  std::vector<std::string> lv = unionv(lx, ly);
  dralign(lv, lx, nx, x);
  dralign(lv, ly, ny, y);
  double sx = 0, sy = 0;
  for (std::size_t j = 0; j < lv.size(); ++j) {
    sx += x[j];
    sy += y[j];
  }
  if (sx <= 0 || sy <= 0) return c.fail("both samples must have at least one observation.");
  double r[4];
  clab_chih(x.data(), y.data(), sat<int>(lv.size()), r);
  if (r[1] < 1 || (sx + sy) / r[3] < 2) {
    out = chout(r[0], std::max(r[1], 0.0), NA, NA);
    return true;
  }
  if (r[2] < 5) {
    return c.fail("an expected count is below 5, and the Monte Carlo p-value of such a table is drawn by R's "
                  "random number generator: call drfthomog() in R.");
  }
  const double p = clab_pchisq(r[0], r[1], 0);
  out = chout(r[0], r[1], p, std::log(p));
  return true;
}
