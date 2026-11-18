// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The plots of R/plotx.R and the palettes of R/plot.R for clab_jrun():
// outside R nothing is drawn, and each plot function returns what the R
// function returns invisibly, the rows it draws, with the checks of R in
// their order and words. The drawing arguments (main, col, pal and those
// of `...`) are not taken. A date is the number of days since 1970-01-01,
// as R stores it.

namespace plt {

using ax::Cols;

// The values of an argument, NA for null.
inline std::vector<double> vals(const Val &o) { return ax::nums(o); }
// Its elements, one per value.
inline std::vector<Val> elems(const Val &o) {
  if (o.t == Val::ARR || o.t == Val::OBJ) return o.a;
  if (o.t == Val::NUL) return {};
  return {o};
}
inline std::size_t len(const Val &o) { return elems(o).size(); }

// .pxg(): x and its groups, the rows with a finite x and a group; the
// groups in order of first appearance.
struct Xg {
  std::vector<double> x;
  ax::Cols g, lv;
};
inline bool pxg(Ctx &c, Xg &d) {
  const Val &xo = c["x"], &go = c["grp"];
  if (!ax::numv(xo)) return c.fail("`x` must be a numeric vector.");
  const std::vector<double> x = vals(xo);
  std::vector<Val> g = elems(go);
  if (go.t == Val::NUL) g.assign(x.size(), vstr("all"));
  if (g.size() != x.size()) return c.fail("`grp` must have one entry per value of `x`.");
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (!std::isfinite(x[i]) || ax::vna(g[i])) continue;
    d.x.push_back(x[i]);
    d.g.push_back(ax::vlab(g[i]));
    if (!ax::has(d.lv, d.g.back())) d.lv.push_back(d.g.back());
  }
  if (d.x.size() < 2) return c.fail("`x` needs at least two finite values.");
  return true;
}

// .pser(): the series, ordered by series (first appearance) and period.
struct Ser {
  ax::Cols g, lv;
  std::vector<Val> prd;
  std::vector<double> y, pos;
  std::vector<int> gi;
};
inline bool pser(Ctx &c, Ser &d) {
  const Val &yo = c["y"], &po = c["prd"], &go = c["grp"];
  if (!ax::numv(yo)) return c.fail("`y` must be a numeric vector.");
  const std::vector<double> y = vals(yo);
  const std::vector<Val> p = elems(po);
  std::vector<Val> g = elems(go);
  if (go.t == Val::NUL) g.assign(y.size(), vstr("series"));
  if (p.size() != y.size() || g.size() != y.size()) {
    return c.fail("`prd` and `grp` must have one entry per value of `y`.");
  }
  bool num = true;
  for (const Val &e : p) {
    double v = 0.0;
    num = num && (e.t == Val::NUL || e.t == Val::NUM || (e.t == Val::STR && snum(e, v)));
  }
  Ser s;
  ax::Cols pl;
  for (std::size_t i = 0; i < y.size(); ++i) {
    if (ax::vna(g[i]) || ax::vna(p[i])) continue;
    s.g.push_back(ax::vlab(g[i]));
    if (!ax::has(s.lv, s.g.back())) s.lv.push_back(s.g.back());
    s.gi.push_back(sat<int>(std::find(s.lv.begin(), s.lv.end(), s.g.back()) - s.lv.begin()));
    s.prd.push_back(p[i]);
    s.y.push_back(y[i]);
    if (num) {
      s.pos.push_back(ax::nv(p[i]));
    } else {
      const std::string l = ax::vlab(p[i]);
      if (!ax::has(pl, l)) pl.push_back(l);
      s.pos.push_back(static_cast<double>(std::find(pl.begin(), pl.end(), l) - pl.begin() + 1));
    }
  }
  if (s.y.empty()) return c.fail("no value has both a period and a series.");
  std::vector<std::size_t> o(s.y.size());
  for (std::size_t i = 0; i < o.size(); ++i) o[i] = i;
  std::stable_sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) {
    return s.gi[a] != s.gi[b] ? s.gi[a] < s.gi[b] : s.pos[a] < s.pos[b];
  });
  d = Ser();
  d.lv = s.lv;
  for (std::size_t i : o) {
    d.g.push_back(s.g[i]);
    d.gi.push_back(s.gi[i]);
    d.prd.push_back(s.prd[i]);
    d.y.push_back(s.y[i]);
    d.pos.push_back(s.pos[i]);
  }
  return true;
}
inline Val sframe(const Ser &d) {
  Frame f;
  f.str("grp", d.g);
  Val p = varr();
  p.a = d.prd;
  f.add("prd", p);
  f.num("y", d.y);
  return f.o;
}
// The median of the values that are not NA (NA for none).
inline double med(const std::vector<double> &x) {
  std::vector<double> v;
  for (double e : x) {
    if (!std::isnan(e)) v.push_back(e);
  }
  return ax::amed(v);
}
// .pcc(): cconc(x, nsim = 0) of a share of a total.
inline bool pcc(Ctx &c, const Val &x, Val &r) {
  const std::vector<double> v = vals(x);
  bool ok = ax::numv(x) && !v.empty();
  double s = 0.0;
  for (double e : v) {
    ok = ok && std::isfinite(e) && e >= 0;
    s += e;
  }
  if (!ok || !(s > 0)) return c.fail("`x` must be finite values from 0 with a positive total.");
  Val a = vobj();
  put(a, "x", x);
  put(a, "nsim", vnum(0));
  return ax::ccall(c, "cconc", cconc, a, r);
}
// The civil date of a day number: year, month, day.
inline void civil(double day, std::int64_t &y, std::int64_t &m, std::int64_t &d) {
  const std::int64_t z = sat<std::int64_t>(std::floor(day)) + 719468;
  const std::int64_t era = (z >= 0 ? z : z - 146096) / 146097;
  const std::int64_t doe = z - era * 146097;
  const std::int64_t yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
  const std::int64_t doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
  const std::int64_t mp = (5 * doy + 2) / 153;
  d = doy - (153 * mp + 2) / 5 + 1;
  m = mp < 10 ? mp + 3 : mp - 9;
  y = yoe + era * 400 + (m <= 2 ? 1 : 0);
}
// The day number of a civil date.
inline double daynum(std::int64_t y, std::int64_t m, std::int64_t d) {
  y -= m <= 2 ? 1 : 0;
  const std::int64_t era = (y >= 0 ? y : y - 399) / 400;
  const std::int64_t yoe = y - era * 400;
  const std::int64_t doy = (153 * (m > 2 ? m - 3 : m + 9) + 2) / 5 + d - 1;
  const std::int64_t doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
  return static_cast<double>(era * 146097 + doe - 719468);
}
// ISO weekday, 1 Monday to 7 Sunday.
inline std::int64_t isowd(double day) {
  const std::int64_t z = sat<std::int64_t>(std::floor(day));
  return ((z + 3) % 7 + 7) % 7 + 1;
}

// The palettes of R/plot.R.
struct Pal {
  const char *name;
  std::vector<const char *> col;
};
const Pal pals[] = {
    {"okabe", {"#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#F0E442", "#000000"}},
    {"diverging", {"#2166AC", "#8C8C8C", "#D6604D"}},
    {"sequential", {"#C6DBEF", "#6BAED6", "#2171B5", "#08306B", "#041E42"}},
    {"grey", {"#252525", "#636363", "#969696", "#BDBDBD", "#D9D9D9"}}};
// colorRampPalette(col)(n): linear in RGB on the 0-1 scale, the result
// times 255 cut to a whole number, as rgb(maxColorValue = 255) does.
inline Cols ramp(const std::vector<const char *> &col, int n) {
  const std::size_t k = col.size();
  std::vector<double> at = {0.0};
  if (k > 1) {
    for (std::size_t i = 1; i + 1 < k; ++i) at.push_back(static_cast<double>(i) * (1.0 / static_cast<double>(k - 1)));
    at.push_back(1.0);
  }
  std::vector<std::array<double, 3>> rgb;
  for (const char *h : col) {
    std::array<double, 3> v{};
    for (std::size_t j = 0; j < 3; ++j) {
      v[j] = static_cast<double>(std::strtol(std::string(h + 1 + 2 * j, 2).c_str(), nullptr, 16)) / 255.0;
    }
    rgb.push_back(v);
  }
  const std::vector<double> t = ax::seqn(0, 1, n);
  ax::Cols out;
  for (double v : t) {
    std::size_t i = 0, j = k - 1;
    while (i + 1 < j) {
      const std::size_t ij = (i + j) / 2;
      if (v < at[ij]) {
        j = ij;
      } else {
        i = ij;
      }
    }
    char b[8];
    unsigned u[3] = {0, 0, 0};
    for (std::size_t q = 0; q < 3; ++q) {
      double y = 0.0;
      if (k == 1 || v == at[j]) {
        y = rgb[j][q];
      } else if (v == at[i]) {
        y = rgb[i][q];
      } else {
        y = rgb[i][q] + (rgb[j][q] - rgb[i][q]) * ((v - at[i]) / (at[j] - at[i]));
      }
      u[q] = sat<unsigned>(y * 255.0);
    }
    std::snprintf(b, sizeof b, "#%02X%02X%02X", u[0] & 0xFFu, u[1] & 0xFFu, u[2] & 0xFFu);
    out.push_back(b);
  }
  return out;
}

}  // namespace plt

// blpal()
bool blpal(Ctx &c, Val &out) {
  const Val &nm = c["name"];
  const plt::Pal *p = nm.t == Val::NUL ? &plt::pals[0] : nullptr;
  if (nm.t == Val::STR) {
    for (const plt::Pal &q : plt::pals) {
      if (nm.s == q.name) p = &q;
    }
  }
  if (p == nullptr) return c.fail("the palette must be one of okabe, diverging, sequential, grey.");
  ax::Cols base(p->col.begin(), p->col.end());
  if (!c.has("n")) {
    out = vstrs(base);
    return true;
  }
  const Val &no = c["n"];
  double n = 0.0;
  if (!(no.t == Val::NUM || (no.t == Val::ARR && no.a.size() == 1 && no.a[0].t == Val::NUM)) ||
      !std::isfinite(n = ax::nv(no.t == Val::NUM ? no : no.a[0])) || n < 1 || n != std::nearbyint(n)) {
    return c.fail("`n` must be one positive whole number.");
  }
  if (n > kutil::maxcell) return c.fail("`n` must be at most 67108864.");
  const int k = sat<int>(n);
  if (std::strcmp(p->name, "okabe") == 0 || std::strcmp(p->name, "grey") == 0) {
    ax::Cols r;
    for (int i = 0; i < k; ++i) r.push_back(base[zz(i) % base.size()]);
    out = unbox(vstrs(r));
  } else if (k == 1) {
    out = vstr(base[0]);
  } else {
    out = vstrs(plt::ramp(p->col, k));
  }
  return true;
}

// yoypals()
bool yoypals(Ctx &, Val &out) {
  out = vstrs({"diverging", "safe", "mono"});
  return true;
}

// plotcdfs()
bool plotcdfs(Ctx &c, Val &out) {
  plt::Xg d;
  CK(plt::pxg(c, d));
  ax::Cols g;
  std::vector<double> x, y;
  for (const std::string &l : d.lv) {
    std::vector<double> v;
    for (std::size_t i = 0; i < d.x.size(); ++i) {
      if (d.g[i] == l) v.push_back(d.x[i]);
    }
    std::sort(v.begin(), v.end());
    for (std::size_t i = 0; i < v.size(); ++i) {
      g.push_back(l);
      x.push_back(v[i]);
      y.push_back(static_cast<double>(i + 1) / static_cast<double>(v.size()));
    }
  }
  Frame f;
  f.str("grp", g);
  f.num("x", x);
  f.num("y", y);
  out = f.o;
  return true;
}

// plotlrnz()
bool plotlrnz(Ctx &c, Val &out) {
  const Val &x = c["x"];
  bool frame = x.t == Val::OBJ;
  for (const Val &e : x.a) frame = frame && e.t == Val::ARR;
  if (!frame) {
    Val r;
    CK(plt::pcc(c, x, r));
    out = col(r, "lorenz");
    return true;
  }
  const Val *ps = ax::fc(x, "pop"), *sh = ax::fc(x, "share");
  if (ps == nullptr || sh == nullptr) return c.fail("a frame `x` needs the columns pop and share.");
  std::vector<double> p = ax::colv(*ps), s = ax::colv(*sh);
  std::vector<std::size_t> o(p.size());
  for (std::size_t i = 0; i < o.size(); ++i) o[i] = i;
  std::stable_sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) {
    if (p[a] != p[b]) return p[a] < p[b];
    return s[a] < s[b];
  });
  std::vector<double> pp, ss;
  for (std::size_t i : o) {
    pp.push_back(p[i]);
    ss.push_back(s[i]);
  }
  bool ok = !pp.empty();
  for (std::size_t i = 0; ok && i < pp.size(); ++i) {
    ok = !std::isnan(pp[i]) && !std::isnan(ss[i]) && pp[i] >= 0 && pp[i] <= 1 && ss[i] >= 0 && ss[i] <= 1 &&
         (i == 0 || ss[i] >= ss[i - 1]);
  }
  if (!ok) return c.fail("`x` must be a Lorenz curve: shares in [0, 1], non-decreasing.");
  if (pp[0] != 0 || ss[0] != 0) {
    pp.insert(pp.begin(), 0.0);
    ss.insert(ss.begin(), 0.0);
  }
  Frame f;
  f.num("pop", pp);
  f.num("share", ss);
  out = f.o;
  return true;
}

// plotgain()
bool plotgain(Ctx &c, Val &out) {
  plt::Xg d;
  CK(plt::pxg(c, d));
  for (double v : d.x) {
    if (v < 0) return c.fail("`x` must not be negative (a share of a total).");
  }
  ax::Cols g;
  std::vector<double> pop, share;
  for (const std::string &l : d.lv) {
    std::vector<double> v;
    for (std::size_t i = 0; i < d.x.size(); ++i) {
      if (d.g[i] == l) v.push_back(d.x[i]);
    }
    Val r;
    CK(plt::pcc(c, vnums(v), r));
    const Val &lz = col(r, "lorenz");
    const std::vector<double> a = ax::colv(col(lz, "pop")), b = ax::colv(col(lz, "share"));
    for (std::size_t i = 0; i < a.size(); ++i) {
      g.push_back(l);
      pop.push_back(a[i]);
      share.push_back(b[i]);
    }
  }
  Frame f;
  f.str("grp", g);
  f.num("pop", pop);
  f.num("share", share);
  out = f.o;
  return true;
}

// plotmult()
bool plotmult(Ctx &c, Val &out) {
  plt::Ser d;
  CK(plt::pser(c, d));
  out = plt::sframe(d);
  return true;
}

// plothrz()
bool plothrz(Ctx &c, Val &out) {
  plt::Ser d;
  CK(plt::pser(c, d));
  double bands = 3;
  if (c.has("bands")) {
    const Val &b = c["bands"];
    const Val &e = b.t == Val::ARR && b.a.size() == 1 ? b.a[0] : b;
    bands = e.t == Val::NUM ? e.n : NaN;
    if (!(std::isfinite(bands) && bands >= 1 && bands <= 6 && bands == std::nearbyint(bands))) {
      return c.fail("`bands` must be one whole number from 1 to 6.");
    }
  }
  std::vector<double> dev(d.y.size()), md;
  for (std::size_t l = 0; l < d.lv.size(); ++l) {
    std::vector<double> v;
    for (std::size_t i = 0; i < d.y.size(); ++i) {
      if (d.gi[i] == sat<int>(l)) v.push_back(d.y[i]);
    }
    md.push_back(plt::med(v));
  }
  double top = 0.0;
  for (std::size_t i = 0; i < d.y.size(); ++i) {
    dev[i] = d.y[i] - md[zz(d.gi[i])];
    if (!std::isnan(dev[i])) top = std::max(top, std::fabs(dev[i]));
  }
  if (top == 0) top = 1;
  const double bw = top / bands;
  std::vector<double> band;
  for (double v : dev) band.push_back(std::isnan(v) ? NA : std::min(bands, std::ceil(std::fabs(v) / bw)));
  Val f = plt::sframe(d);
  put(f, "band", vnums(band));
  out = f;
  return true;
}

// plotbump()
bool plotbump(Ctx &c, Val &out) {
  plt::Ser d;
  CK(plt::pser(c, d));
  bool desc = true;
  if (c.has("desc") && !ax::tf(c["desc"], desc)) return c.fail("`desc` must be TRUE or FALSE.");
  const std::size_t n = d.y.size();
  std::vector<double> rk(n, NA);
  std::vector<char> done(n, 0);
  for (std::size_t i = 0; i < n; ++i) {
    if (done[i]) continue;
    std::vector<std::size_t> s;
    for (std::size_t j = i; j < n; ++j) {
      if (d.pos[j] == d.pos[i] && !done[j]) s.push_back(j);
    }
    // rank(v, ties.method = "min"), NA last
    std::vector<std::size_t> ok, na;
    for (std::size_t j : s) {
      done[j] = 1;
      (std::isnan(d.y[j]) ? na : ok).push_back(j);
    }
    for (std::size_t j : ok) {
      const double vj = desc ? -d.y[j] : d.y[j];
      double r = 1;
      for (std::size_t q : ok) r += (desc ? -d.y[q] : d.y[q]) < vj ? 1 : 0;
      rk[j] = r;
    }
    for (std::size_t q = 0; q < na.size(); ++q) rk[na[q]] = static_cast<double>(ok.size() + q + 1);
  }
  Val f = plt::sframe(d);
  put(f, "rank", vnums(rk));
  out = f;
  return true;
}

// plotfnl()
bool plotfnl(Ctx &c, Val &out) {
  const Val &ob = c["obs"], &ex = c["expd"], &ar = c["area"];
  Val lvo = c.has("lvls") ? c["lvls"] : vnums({0.95, 0.998});
  const std::pair<const Val *, const char *> nn[] = {{&ob, "obs"}, {&ex, "expd"}, {&lvo, "lvls"}};
  for (const auto &a : nn) {
    if (!ax::numv(*a.first) && !(a.first->t == Val::ARR && a.first->a.empty())) {
      return c.fail("`%s` must be numbers.", a.second);
    }
  }
  const std::vector<double> o0 = plt::vals(ob), e0 = plt::vals(ex);
  const std::size_t n = o0.size();
  if (n != e0.size() || n == 0) return c.fail("`obs` and `expd` must be the same length, from 1.");
  const std::vector<Val> al = plt::elems(ar);
  if (ar.t != Val::NUL && al.size() != n) return c.fail("`area` must have one label per area.");
  std::vector<double> o, e;
  std::vector<Val> a;
  for (std::size_t i = 0; i < n; ++i) {
    if (std::isnan(o0[i]) || std::isnan(e0[i])) continue;
    o.push_back(o0[i]);
    e.push_back(e0[i]);
    a.push_back(ar.t == Val::NUL ? Val() : vstr(ax::vlab(al[i])));
  }
  bool ok = !o.empty();
  for (std::size_t i = 0; ok && i < o.size(); ++i) {
    ok = std::isfinite(o[i]) && std::isfinite(e[i]) && o[i] >= 0 && e[i] > 0;
  }
  if (!ok) return c.fail("the counts must be finite, `obs` from 0 and `expd` positive.");
  std::vector<double> lv = plt::vals(lvo);
  bool lok = !lv.empty();
  for (double v : lv) lok = lok && !std::isnan(v) && v > 0 && v < 1;
  if (!lok) return c.fail("`lvls` must be probabilities strictly between 0 and 1.");
  double tg = 1;
  if (c.has("target")) {
    const Val &t = c["target"];
    const Val &te = t.t == Val::ARR && t.a.size() == 1 ? t.a[0] : t;
    tg = te.t == Val::NUM ? te.n : NaN;
    if (!(std::isfinite(tg) && tg >= std::numeric_limits<double>::min())) {
      return c.fail("`target` must be one positive number.");
    }
  }
  const double pm = *std::max_element(lv.begin(), lv.end()), al2 = (1 - pm) / 2;
  std::vector<double> r;
  std::vector<int> out1;
  for (std::size_t i = 0; i < o.size(); ++i) {
    const double ri = o[i] / e[i];
    const double lo = clab_qpois(al2, tg * e[i], 1) / e[i], hi = clab_qpois(1 - al2, tg * e[i], 1) / e[i];
    r.push_back(ri);
    out1.push_back(ri < lo || ri > hi ? 1 : 0);
  }
  Frame f;
  Val av = varr();
  av.a = a;
  f.add("area", av);
  f.num("observed", o);
  f.num("expected", e);
  f.num("ratio", r);
  f.lgl("outside", out1.data(), out1.size());
  out = f.o;
  return true;
}

// plotcal(): x days since 1970-01-01 or "YYYY-MM-DD".
bool plotcal(Ctx &c, Val &out) {
  const Val &xo = c["x"], &yo = c["y"];
  const std::vector<Val> xs = plt::elems(xo);
  const std::vector<double> y = plt::vals(yo);
  bool ok = ax::numv(yo) && xs.size() == y.size();
  std::vector<double> dt;
  for (const Val &e : xs) {
    double v = NA;
    int yy = 0, mm = 0, dd = 0;
    char tail = 0;
    if (e.t == Val::NUM) {
      v = e.n;
    } else if (e.t == Val::STR && std::sscanf(e.s.c_str(), "%d-%d-%d%c", &yy, &mm, &dd, &tail) == 3 && mm >= 1 &&
               mm <= 12 && dd >= 1 && dd <= 31) {
      v = plt::daynum(yy, mm, dd);
    } else if (e.t != Val::NUL) {
      ok = false;
    }
    // a day number is a date up to the year 9999, or it is no date (and
    // no whole number the calendar arithmetic below can hold)
    if (!std::isnan(v) && !(std::fabs(v) <= 2932896)) ok = false;
    dt.push_back(v);
  }
  if (!ok) return c.fail("`x` must be dates, one per value of the numbers `y`.");
  std::vector<double> d, yr, wk, wd, yv;
  for (std::size_t i = 0; i < dt.size(); ++i) {
    if (std::isnan(dt[i])) continue;
    std::int64_t a = 0, m = 0, q = 0;
    plt::civil(dt[i], a, m, q);
    const double j1 = plt::daynum(a, 1, 1);
    d.push_back(std::floor(dt[i]));
    yr.push_back(static_cast<double>(a));
    wd.push_back(static_cast<double>(plt::isowd(dt[i])));
    wk.push_back(std::floor((std::floor(dt[i]) - j1 + static_cast<double>(plt::isowd(j1)) - 1) / 7));
    yv.push_back(y[i]);
  }
  Frame f;
  f.num("date", d);
  f.num("year", yr);
  f.num("week", wk);
  f.num("wday", wd);
  f.num("y", yv);
  out = f.o;
  return true;
}
