// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// exprspplr() (src/plr.cpp) with the built-in estimator, with the checks of
// R/aexp.R in their words. R draws the resamples and folds with its own
// generator; this is that generator: the Mersenne-Twister as set.seed()
// seeds it and sample() with the "Rejection" kind, R's defaults, so the
// draws are R's for the same seed.

struct Rmt {
  std::uint32_t mt[624];
  int mti = 624;
  static void init(std::int32_t s, std::uint32_t *m) {
    std::uint32_t sd = static_cast<std::uint32_t>(s);
    for (int j = 0; j < 50; ++j) sd = 69069 * sd + 1;
    sd = 69069 * sd + 1;  // the word that holds the position
    for (int j = 0; j < 624; ++j) {
      sd = 69069 * sd + 1;
      m[j] = sd;
    }
  }
  // set.seed(), and past 2147483647 the state crimlab and the reference
  // implementation build from the seed's two halves
  explicit Rmt(double seed) {
    if (std::fabs(seed) <= 2147483647.0) {
      init(static_cast<std::int32_t>(seed), mt);
      return;
    }
    std::uint32_t lo[624], hi[624];
    init(static_cast<std::int32_t>(std::fmod(seed, 2147483648.0)), lo);
    init(-static_cast<std::int32_t>(std::floor(seed / 2147483648.0)), hi);
    const std::uint32_t na = 0x80000000u;
    for (int j = 0; j < 624; ++j) {
      const std::uint32_t a = lo[j], b = hi[(j + 1) % 624], w = a ^ b;
      mt[j] = a == na || b == na || w == na ? 1u : w;
    }
  }
  double genrand() {
    static const std::uint32_t mag[2] = {0x0u, 0x9908b0dfu};
    std::uint32_t y;
    if (mti >= 624) {
      int kk;
      for (kk = 0; kk < 624 - 397; ++kk) {
        y = (mt[kk] & 0x80000000u) | (mt[kk + 1] & 0x7fffffffu);
        mt[kk] = mt[kk + 397] ^ (y >> 1) ^ mag[y & 0x1u];
      }
      for (; kk < 623; ++kk) {
        y = (mt[kk] & 0x80000000u) | (mt[kk + 1] & 0x7fffffffu);
        mt[kk] = mt[kk + (397 - 624)] ^ (y >> 1) ^ mag[y & 0x1u];
      }
      y = (mt[623] & 0x80000000u) | (mt[0] & 0x7fffffffu);
      mt[623] = mt[396] ^ (y >> 1) ^ mag[y & 0x1u];
      mti = 0;
    }
    y = mt[mti++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680u;
    y ^= (y << 15) & 0xefc60000u;
    y ^= (y >> 18);
    return static_cast<double>(y) * 2.3283064365386963e-10;
  }
  double unif() {
    const double i2 = 2.328306437080797e-10, x = genrand();
    if (x <= 0) return 0.5 * i2;
    if (1 - x <= 0) return 1 - 0.5 * i2;
    return x;
  }
  // R_unif_index(): rejection sampling from the next power of two
  double index(double dn) {
    if (dn <= 0) return 0;
    const int bits = static_cast<int>(std::ceil(std::log2(dn)));
    double dv;
    do {
      std::int_least64_t v = 0;
      for (int n = 0; n <= bits; n += 16) {
        const int v1 = static_cast<int>(std::floor(unif() * 65536));
        v = 65536 * v + v1;
      }
      dv = static_cast<double>(v & ((static_cast<std::int_least64_t>(1) << bits) - 1));
    } while (dn <= dv);
    return dv;
  }
};

// The rows of the complete cases of the columns cs (null and NaN missing).
inline std::vector<int> ccase(const std::vector<const Val *> &cs, std::size_t n) {
  std::vector<int> r;
  for (std::size_t i = 0; i < n; ++i) {
    bool ok = true;
    for (const Val *c : cs) {
      const Val &e = c->t == Val::ARR ? (i < c->a.size() ? c->a[i] : nul) : *c;
      double x = 0;
      ok = ok && e.t != Val::NUL && !(e.t == Val::NUM && std::isnan(e.n)) && !(snum(e, x) && std::isnan(x));
    }
    if (ok) r.push_back(static_cast<int>(i));
  }
  return r;
}

bool exprspplr(Ctx &c, Val &out) {
  std::string y, ex;
  if (!c.has("y")) return ckstop(c, "`y` must name one column of `data`");
  CK(colnm(c, "y", "", y));
  if (!c.has("expo")) return ckstop(c, "`expo` must name one column of `data`");
  CK(colnm(c, "expo", "", ex));
  std::vector<std::string> cf;
  {
    const Val &o = c["confs"];
    bool ok = o.t == Val::STR || (o.t == Val::ARR && !o.a.empty());
    if (o.t == Val::STR) cf.push_back(o.s);
    else if (o.t == Val::ARR) {
      for (const Val &e : o.a) {
        ok = ok && e.t == Val::STR;
        if (e.t == Val::STR) cf.push_back(e.s);
      }
    }
    for (const std::string &s : cf) ok = ok && !s.empty();
    if (!ok) return ckstop(c, "`confs` must name one or more columns of `data`");
  }
  double rs = 42;
  if (c.has("rstate")) {
    const double v = val1(c["rstate"]);
    bool lg = false;
    if (jcls(c["rstate"], lg) != nullptr || !std::isfinite(v) || v != std::nearbyint(v) || v < -2147483647.0 ||
        v >= 9007199254740992.0) {
      return ckstop(c, "`rstate` must be a single whole number from -2147483647 to 2^53 - 1");
    }
    rs = v;
  }
  double nb = 100;
  if (c.has("nboot")) CK(polsc(c, c["nboot"], "nboot", nb));
  if (nb != std::nearbyint(nb) || nb > 1e6) {
    return ckstop(c, "`nboot` must be a whole number of resamples, at most 1000000");
  }
  const Val &d = c["data"];
  std::vector<std::string> cols = {y, ex};
  cols.insert(cols.end(), cf.begin(), cf.end());
  std::string miss;
  std::vector<std::string> seen;
  for (const std::string &k : cols) {
    if (std::find(seen.begin(), seen.end(), k) != seen.end()) continue;
    seen.push_back(k);
    if (d.t != Val::OBJ || col(d, k.c_str()).t == Val::NUL) miss += (miss.empty() ? "" : ", ") + k;
  }
  if (!miss.empty()) return ckstop(c, "data is missing columns: %s", miss.c_str());
  std::size_t n = 0;
  for (std::size_t j = 0; j < d.a.size(); ++j) n = std::max(n, vlen(d.a[j]));
  std::vector<const Val *> cv;
  for (const std::string &k : cols) cv.push_back(&col(d, k.c_str()));
  const std::vector<int> keep = ccase(cv, n);
  std::vector<std::vector<double>> num(cols.size());
  for (std::size_t j = 0; j < cols.size(); ++j) {
    Dv t;
    bool lg = false;
    const Val &o = *cv[j];
    if (jcls(o, lg) != nullptr || readd(o, t) != NUM) {
      // a column of missing values only is logical in R
      return ckstop(c, "column `%s` must be numeric (a factor or character column is not a code)",
                    cols[j].c_str());
    }
    num[j] = t.v;
    if (num[j].size() == 1 && n > 1) num[j].assign(n, t.v[0]);
  }
  if (keep.size() < 10) {
    return ckstop(c, "need at least 10 complete rows for 5 folds; have %d", static_cast<int>(keep.size()));
  }
  const int p = static_cast<int>(cf.size());
  std::vector<double> xm(n * static_cast<std::size_t>(p));
  for (int j = 0; j < p; ++j) {
    for (std::size_t i = 0; i < n; ++i) xm[static_cast<std::size_t>(j) * n + i] = num[static_cast<std::size_t>(2 + j)][i];
  }
  std::vector<char> ok(n, 0);
  for (int r : keep) ok[static_cast<std::size_t>(r)] = 1;
  // the folds: sample(rep_len(1:5, m)) after set.seed(1), one set per size
  auto folds = [](int m) {
    Rmt g(1);
    std::vector<int> x(static_cast<std::size_t>(m)), f(static_cast<std::size_t>(m));
    for (int i = 0; i < m; ++i) x[static_cast<std::size_t>(i)] = i;
    int left = m;
    for (int i = 0; i < m; ++i) {
      const int j = static_cast<int>(g.index(left));
      f[static_cast<std::size_t>(i)] = x[static_cast<std::size_t>(j)] % 5 + 1;
      x[static_cast<std::size_t>(j)] = x[static_cast<std::size_t>(--left)];
    }
    return f;
  };
  double o[2];
  const int nn = static_cast<int>(n);
  int rc = clab_plr(num[0].data(), num[1].data(), p ? xm.data() : nullptr, nn, p, keep.data(),
                    static_cast<int>(keep.size()), folds(static_cast<int>(keep.size())).data(), 5, o);
  if (rc == CLAB_EDOM) {
    return ckstop(c, "the treatment has no variation left after partialling out the covariates");
  }
  CK(kchk(c, rc));
  const double ate = o[0], se = o[1];
  Rmt g(rs);
  std::vector<double> est;
  for (int b = 0; b < static_cast<int>(nb); ++b) {
    std::vector<int> rows;
    for (std::size_t i = 0; i < n; ++i) {
      const int r = static_cast<int>(g.index(static_cast<double>(n)));
      if (ok[static_cast<std::size_t>(r)]) rows.push_back(r);
    }
    if (rows.size() < 10) continue;
    rc = clab_plr(num[0].data(), num[1].data(), p ? xm.data() : nullptr, nn, p, rows.data(),
                  static_cast<int>(rows.size()), folds(static_cast<int>(rows.size())).data(), 5, o);
    if (rc == CLAB_OK) est.push_back(o[0]);
  }
  if (est.size() < 10) {
    return ckstop(c,
                  "Only %d successful bootstrap fits; need >=10. Inspect the model's convergence on your resamples",
                  static_cast<int>(est.size()));
  }
  const double pr[2] = {0.025, 0.975};
  double s[3];
  clab_bootci(est.data(), static_cast<int>(est.size()), pr, 2, s);
  out = vobj();
  put(out, "ate", vnum(ate));
  put(out, "seana", vnum(se));
  put(out, "seboot", vnum(s[0]));
  put(out, "cilobs", vnum(s[1]));
  put(out, "ciupbs", vnum(s[2]));
  put(out, "nboot", vnum(static_cast<double>(est.size())));
  put(out, "method", vstr("partially linear regression, DML2 cross-fitting, OLS nuisances + percentile bootstrap"));
  return true;
}
