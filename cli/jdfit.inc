// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// dstfit() (src/dfit.cpp), with the checks of R/dfit.R in their words. The
// bootstrap draws are R's (jaexp.inc), and the parameters are written as
// R's format() writes a numeric vector.

// format() of a numeric vector, seven significant digits: fixed notation
// with the decimals of the most precise element, padded to one width, or
// scientific where that is narrower.
inline std::vector<std::string> rfmtv(const std::vector<double> &x) {
  int mxsl = 1, rgt = 0, mxns = 1, neg = 0, mxe = 0;
  std::vector<int> kp(x.size()), ns(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double v = x[i];
    int e = 0, s = 1;
    if (v != 0 && std::isfinite(v)) {
      char b[64];
      std::snprintf(b, sizeof b, "%.6e", v);
      const double r7 = std::strtod(b, nullptr);
      for (s = 1; s <= 7; ++s) {
        std::snprintf(b, sizeof b, "%.*e", s - 1, v);
        if (std::strtod(b, nullptr) == r7) break;
      }
      e = std::atoi(std::strchr(b, 'e') + 1);
    }
    kp[i] = e;
    ns[i] = s;
    const int left = (e >= 0 ? e + 1 : 1) + (v < 0);
    mxsl = std::max(mxsl, left);
    rgt = std::max(rgt, std::max(0, s - e - 1));
    mxns = std::max(mxns, s);
    neg = neg || v < 0;
    mxe = std::max(mxe, std::abs(e));
  }
  const int wf = mxsl + (rgt > 0 ? rgt + 1 : 0);
  const int ws = neg + (mxns > 1 ? mxns + 1 : mxns) + (mxe >= 100 ? 5 : 4);
  std::vector<std::string> r;
  for (double v : x) {
    char b[64];
    if (wf <= ws) {
      std::snprintf(b, sizeof b, "%*.*f", wf, rgt, v);
    } else {
      std::snprintf(b, sizeof b, "%.*e", mxns - 1, v);
      std::string s = b;
      const std::size_t k = s.find('e');
      // R writes two exponent digits at least
      std::string ex = s.substr(k + 2);
      while (ex.size() > 2 && ex[0] == '0') ex.erase(0, 1);
      s = s.substr(0, k + 2) + ex;
      std::snprintf(b, sizeof b, "%*s", ws, s.c_str());
    }
    r.push_back(b);
  }
  return r;
}

bool dstfit(Ctx &c, Val &out) {
  static const char *const fam[] = {"normal", "lognormal", "exponential", "gamma", "weibull", "t",
                                    "pareto", "beta",      "poisson",     "negbin", "zip",     "zinb"};
  static const char *const pars[12][3] = {{"mean", "sd", ""},       {"meanlog", "sdlog", ""}, {"rate", "", ""},
                                          {"shape", "rate", ""},    {"shape", "scale", ""},   {"location", "scale", "df"},
                                          {"xmin", "alpha", ""},    {"shape1", "shape2", ""}, {"lambda", "", ""},
                                          {"size", "mu", ""},       {"pi", "lambda", ""},     {"pi", "size", "mu"}};
  const char *pos = "every value above 0", *cn = "counts (whole numbers, 0 or more)", *zz = "counts with at least one zero";
  const char *need[12] = {"numbers", pos, "no value below 0", pos, pos, "numbers", pos,
                          "every value between 0 and 1, exclusive", cn, cn, zz, zz};
  const Val &xo = c["x"];
  bool lg = false;
  if (xo.t == Val::NUL || jcls(xo, lg) != nullptr) return ckstop(c, "`x` must be a numeric vector");
  Dv d;
  readd(xo, d);
  int nmiss = 0;
  std::vector<double> x;
  for (double v : d.v) {
    if (std::isinf(v)) return ckstop(c, "`x` has an infinite value");
    if (std::isnan(v)) ++nmiss;
    else x.push_back(v);
  }
  std::sort(x.begin(), x.end());
  if (x.size() < 5) return ckstop(c, "`x` needs at least five values that are not missing");
  if (x.front() == x.back()) return ckstop(c, "every value of `x` is the same; there is no distribution to fit");
  double nb = 200;
  if (c.has("nboot")) {
    nb = val1(c["nboot"]);
    if (std::isfinite(nb) && nb != std::floor(nb)) {
      return ckstop(c, "`nboot` must be a whole number, not %s", rfmt(nb).c_str());
    }
    if (!(nb >= 0 && nb <= 1e5)) return ckstop(c, "`nboot` must be one whole number from 0 to 1e+05");
  }
  double lvl = 0.95;
  if (c.has("conflvl")) {
    lvl = val1(c["conflvl"]);
    if (!(lvl > 0 && lvl < 1)) return ckstop(c, "`conflvl` must be a single number strictly between 0 and 1");
  }
  bool posv = true, cnt = true, z0 = false, nn = true, unit = true;
  for (double v : x) {
    posv = posv && v > 0;
    nn = nn && v >= 0;
    cnt = cnt && v >= 0 && v == std::nearbyint(v);
    z0 = z0 || v == 0;
    unit = unit && v > 0 && v < 1;
  }
  const bool okf[12] = {true, posv, nn, posv, posv, true, posv, unit, cnt, cnt, cnt && z0, cnt && z0};
  std::vector<int> fs;
  if (!c.has("fams")) {
    for (int i = 0; i < 12; ++i) {
      if (okf[i] && (i >= 8 || !cnt || i == 0 || i == 5)) fs.push_back(i);
    }
  } else {
    std::vector<std::string> nm;
    const Val &fo = c["fams"];
    bool good = labs(fo, nm) && !nm.empty();
    for (const std::string &s : nm) {
      int hit = -1, many = 0;
      for (int i = 0; i < 12; ++i) {
        const std::string f = fam[i];
        if (f == s) {
          hit = i;
          many = 1;
          break;
        }
        if (!s.empty() && f.compare(0, s.size(), s) == 0) {
          hit = i;
          ++many;
        }
      }
      if (hit < 0 || many != 1) good = false;
      else if (std::find(fs.begin(), fs.end(), hit) == fs.end()) fs.push_back(hit);
    }
    if (!good) {
      std::string ch;
      for (int i = 0; i < 12; ++i) ch += std::string(i ? ", " : "") + "\"" + fam[i] + "\"";
      return ckstop(c, "`fams` must be one or more of %s, not %s", ch.c_str(), jdep(fo).c_str());
    }
    for (int i : fs) {
      if (!okf[i]) return ckstop(c, "`%s` needs %s", fam[i], need[i]);
    }
  }
  const int n = sat<int>(x.size());
  struct Row {
    int f, k;
    bool ok;
    double ll, aic, bic, ks, ad;
    std::vector<double> p;
  };
  std::vector<Row> rows;
  for (int f : fs) {
    Row r{f, clab_dfnp(f), false, NA, NA, NA, NA, NA, std::vector<double>(3)};
    if (clab_dfit(f, x.data(), n, r.p.data()) == CLAB_OK) {
      r.ok = true;
      r.p.resize(sat<std::size_t>(r.k));
      r.ll = clab_dfll(f, x.data(), n, r.p.data());
      r.aic = 2 * r.k - 2 * r.ll;
      r.bic = std::log(static_cast<double>(n)) * r.k - 2 * r.ll;
      double g[2];
      clab_dfgof(f, x.data(), n, r.p.data(), g);
      r.ks = g[0];
      r.ad = f >= 8 ? NA : g[1];
    }
    rows.push_back(r);
  }
  std::stable_sort(rows.begin(), rows.end(), [](const Row &a, const Row &b) {
    if (std::isnan(a.aic) || std::isnan(b.aic)) return !std::isnan(a.aic) && std::isnan(b.aic);
    return a.aic < b.aic;
  });
  Val rk = varr(), fm = varr(), np = varr(), ll = varr(), aic = varr(), bic = varr(), ks = varr(), ad = varr(),
      cv = varr(), ps = varr();
  Val fits = vobj();
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const Row &r = rows[i];
    rk.a.push_back(vnum(static_cast<double>(i) + 1));
    fm.a.push_back(vstr(fam[r.f]));
    np.a.push_back(r.ok ? vnum(r.k) : Val());
    ll.a.push_back(vnum(r.ll));
    aic.a.push_back(vnum(r.aic));
    bic.a.push_back(vnum(r.bic));
    ks.a.push_back(vnum(r.ks));
    ad.a.push_back(vnum(r.ad));
    cv.a.push_back(vbool(r.ok));
    std::string t;
    if (r.ok) {
      std::vector<double> sg;
      for (double v : r.p) {
        char b[64];
        std::snprintf(b, sizeof b, "%.4g", v);
        sg.push_back(std::strtod(b, nullptr));
      }
      const std::vector<std::string> fv = rfmtv(sg);
      Val fo = vobj();
      for (int j = 0; j < r.k; ++j) {
        t += std::string(j ? ", " : "") + pars[r.f][j] + " = " + fv[sat<std::size_t>(j)];
        put(fo, pars[r.f][j], vnum(r.p[sat<std::size_t>(j)]));
      }
      put(fits, fam[r.f], fo);
    }
    ps.a.push_back(vstr(t));
  }
  Val tab = vobj();
  put(tab, "rank", rk);
  put(tab, "family", fm);
  put(tab, "npar", np);
  put(tab, "loglik", ll);
  put(tab, "aic", aic);
  put(tab, "bic", bic);
  put(tab, "ks", ks);
  put(tab, "ad", ad);
  put(tab, "converged", cv);
  put(tab, "params", ps);
  // the bootstrap: R's generator from seed, nboot resamples per fitted family
  const int B = sat<int>(nb);
  const double a = (1 - lvl) / 2;
  Val bf = varr(), bp = varr(), be = varr(), bl = varr(), bu = varr();
  double sd = 1;
  if (B > 0 && c.has("seed")) {
    sd = val1(c["seed"]);
    if (!std::isfinite(sd) || sd != std::nearbyint(sd) || sd < -2147483647.0 || sd >= 9007199254740992.0) {
      return ckstop(c, "`seed` must be a single whole number from -2147483647 to 2^53 - 1");
    }
  }
  Rmt g(sd);
  for (const Row &r : rows) {
    if (!r.ok) continue;
    std::vector<std::vector<double>> est(sat<std::size_t>(r.k));
    std::vector<double> xb(sat<std::size_t>(n)), p(3);
    for (int b = 0; b < B; ++b) {
      for (int i = 0; i < n; ++i) xb[sat<std::size_t>(i)] = x[sat<std::size_t>(g.index(n))];
      if (clab_dfit(r.f, xb.data(), n, p.data()) != CLAB_OK) continue;
      for (int j = 0; j < r.k; ++j) est[sat<std::size_t>(j)].push_back(p[sat<std::size_t>(j)]);
    }
    for (int j = 0; j < r.k; ++j) {
      bf.a.push_back(vstr(fam[r.f]));
      bp.a.push_back(vstr(pars[r.f][j]));
      be.a.push_back(vnum(r.p[sat<std::size_t>(j)]));
      const std::vector<double> &e = est[sat<std::size_t>(j)];
      double q[3] = {NA, NA, NA};
      const double pr[2] = {a, 1 - a};
      if (!e.empty()) clab_bootci(e.data(), sat<int>(e.size()), pr, 2, q);
      bl.a.push_back(vnum(e.empty() ? NA : q[1]));
      bu.a.push_back(vnum(e.empty() ? NA : q[2]));
    }
  }
  Val boot = vobj();
  put(boot, "family", bf);
  put(boot, "param", bp);
  put(boot, "estimate", be);
  put(boot, "lower", bl);
  put(boot, "upper", bu);
  Val hill;
  if (n >= 10 && posv) {
    int k = sat<int>(std::floor(std::sqrt(static_cast<double>(n))));
    if (c.has("k")) {
      const double v = val1(c["k"]);
      if (std::isfinite(v) && v != std::floor(v)) return ckstop(c, "`k` must be a whole number, not %s", rfmt(v).c_str());
      if (!(v >= 2 && v <= n - 1)) return ckstop(c, "`k` must be one whole number from 2 to %d", n - 1);
      k = sat<int>(v);
    }
    long double s = 0;
    for (int i = 0; i < k; ++i) {
      s += std::log(x[sat<std::size_t>(n - 1 - i)] / x[sat<std::size_t>(n - 1 - k)]);
    }
    const double al = k / static_cast<double>(s), z = clab_qnorm(1 - (1 - lvl) / 2);
    hill = vobj();
    put(hill, "k", vnum(k));
    put(hill, "alpha", vnum(al));
    put(hill, "lower", vnum(al * (1 - z / std::sqrt(static_cast<double>(k)))));
    put(hill, "upper", vnum(al * (1 + z / std::sqrt(static_cast<double>(k)))));
  }
  out = vobj();
  put(out, "table", tab);
  put(out, "fits", fits);
  put(out, "boot", boot);
  put(out, "hill", hill);
  put(out, "x", unbox(vnums(x)));
  put(out, "nmiss", vnum(nmiss));
  put(out, "discrete", vbool(cnt));
  put(out, "conflvl", vnum(lvl));
  return true;
}
