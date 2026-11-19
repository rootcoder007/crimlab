// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// plotsvg(), second part: the plots R draws of results with plot(), by
// the shared shapes of R/plot.R (.pdots, .phist, .pgrp, .pbars, .pmat,
// .psurf) and one short builder for each result. Included by jsvg.inc
// inside namespace sv, after the drawing primitives.

// A number as R's format() shows it: `sig` significant digits.
inline std::string gfmt(double v, int sig = 4) {
  if (!std::isfinite(v)) return "NA";
  char b[64];
  std::snprintf(b, sizeof b, "%.*g", sig, v);
  std::string s = b;
  for (char &ch : s) {
    if (ch == ',') ch = '.';
  }
  return s;
}
// A column of TRUE and FALSE (null is FALSE).
inline std::vector<bool> bv(const Val &o, const char *k) {
  std::vector<bool> r;
  for (const Val &e : plt::elems(get(o, k))) r.push_back(e.t == Val::BOOL && e.b);
  return r;
}
// A matrix written as rows.
inline std::vector<std::vector<double>> mat(const Val &m) {
  std::vector<std::vector<double>> r;
  for (const Val &row : plt::elems(m)) r.push_back(ax::nums(row));
  return r;
}
inline bool has(const Val &o, const char *k) { return get(o, k).t != Val::NUL; }
// The text of one value.
inline std::string str(const Val &o) {
  const Val *e = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  return e->t == Val::STR ? e->s : e->t == Val::NUM ? gfmt(e->n, 7) : "";
}
inline double num1(const Val &o) {
  const std::vector<double> v = ax::nums(o);
  return v.size() == 1 ? v[0] : NA;
}
// The point shapes of .ppch(): disc, triangle, square, ring, triangle, square.
inline int pch(std::size_t g) {
  static const int p[] = {0, 1, 3, 2, 1, 3};
  return p[g % 6];
}
// The index of each label among the labels in order of appearance.
inline std::vector<std::size_t> idx(const ax::Cols &g, const ax::Cols &lv) {
  std::vector<std::size_t> r;
  for (const std::string &s : g) r.push_back(sat<std::size_t>(std::find(lv.begin(), lv.end(), s) - lv.begin()));
  return r;
}
// R's pretty() breaks for a histogram of about n cells.
inline std::vector<double> pretty(double lo, double hi, int n) {
  if (hi - lo <= 0) {
    lo -= 0.5;
    hi += 0.5;
  }
  const double ns = (hi - lo) / n, base = std::pow(10.0, std::floor(std::log10(ns)));
  double unit = base;
  for (double m : {2.0, 5.0, 10.0}) {
    if (m * base - ns < 1.5 * (ns - unit)) unit = m * base;
  }
  std::vector<double> b;
  for (double v = std::floor(lo / unit + 1e-7) * unit; v < hi + unit * (1 - 1e-7); v += unit) b.push_back(v);
  if (b.size() < 2) b.push_back(b.back() + unit);
  return b;
}
// An axis in log10 units: ticks at 1, 2 and 5 times the powers of ten in
// the range, labelled with the values themselves.
inline void logax(Ax &a, double lo, double hi) {
  std::vector<double> at, all;
  for (int e = sat<int>(std::floor(lo)) - 1; e <= sat<int>(std::ceil(hi)) + 1; ++e) {
    for (double m : {1.0, 2.0, 5.0}) {
      const double pos = e + std::log10(m);
      if (pos < lo - 1e-9 || pos > hi + 1e-9) continue;
      all.push_back(pos);
      if (m == 1.0) at.push_back(pos);
    }
  }
  if (all.size() <= 8 || at.size() < 2) at = all;
  std::vector<double> val;
  for (double pos : at) val.push_back(std::pow(10.0, pos));
  // a range within one step of 1, 2, 5: round values of the scale itself
  if (at.size() < 3) val = pretty(std::pow(10.0, lo), std::pow(10.0, hi), 5);
  a.own = true;
  a.at.clear();
  a.lab.clear();
  for (double v : val) {
    if (v <= 0 || std::log10(v) < lo - 1e-9 || std::log10(v) > hi + 1e-9) continue;
    a.at.push_back(std::log10(v));
    a.lab.push_back(gfmt(v, 3));
  }
}

// .pdots(): one row per label, the estimate as a point and its interval
// as a capped bar, groups told apart by colour and point shape.
struct Dots {
  ax::Cols lab, grp;
  std::vector<double> est, lo, hi;
  std::string main, xlab;
  double ref = NA;
  bool logx = false, key = true;
};
inline Pan dots(const Fig &f, Dots d) {
  const std::size_t n = d.lab.size();
  d.lo.resize(n, NA);
  d.hi.resize(n, NA);
  const bool grouped = !d.grp.empty();
  if (!grouped) d.grp.assign(n, "");
  const ax::Cols u = lvls(d.lab), lv = lvls(d.grp), cl = cols(f, sat<int>(lv.size()));
  const std::vector<std::size_t> gi = idx(d.grp, lv), li = idx(d.lab, u);
  const double k = static_cast<double>(u.size());
  auto tr = [&](double v) { return d.logx ? (v > 0 ? std::log10(v) : NA) : v; };
  std::vector<double> all;
  for (std::size_t i = 0; i < n; ++i) {
    all.push_back(tr(d.est[i]));
    all.push_back(tr(d.lo[i]));
    all.push_back(tr(d.hi[i]));
  }
  if (std::isfinite(d.ref)) all.push_back(tr(d.ref));
  Pan p = frame(d.main.c_str(), d.xlab.c_str(), "");
  rng(all, p.x0, p.x1);
  const double top = grouped && d.key ? 0.75 * static_cast<double>(lv.size()) + 0.3 : 0;
  p.y0 = 0.5;
  p.y1 = k + 0.5 + top;
  p.pad = false;
  const double dx = 0.04 * (p.x1 - p.x0);
  p.x0 -= dx;
  p.x1 += dx;
  if (d.logx) logax(p.ax, p.x0, p.x1);
  std::vector<double> at;
  for (std::size_t i = 0; i < u.size(); ++i) at.push_back(k - static_cast<double>(i));
  labs(p.ay, u, at);
  if (std::isfinite(d.ref)) p.it.push_back(ref(tr(d.ref), 0, tr(d.ref), k + 1 + top, 2));
  It sg, pt;
  sg.kind = SEG;
  sg.cap = 4;
  sg.hz = true;
  sg.lw = 2;
  pt.kind = PTS;
  pt.r = 3.8;
  bool any = false;
  for (std::size_t i = 0; i < n; ++i) {
    const double y = k - static_cast<double>(li[i]);
    if (std::isfinite(tr(d.lo[i])) && std::isfinite(tr(d.hi[i]))) {
      sg.x.insert(sg.x.end(), {tr(d.lo[i]), tr(d.hi[i])});
      sg.y.insert(sg.y.end(), {y, y});
      sg.cols.insert(sg.cols.end(), {cl[gi[i]], cl[gi[i]]});
      any = true;
    }
    pt.x.push_back(tr(d.est[i]));
    pt.y.push_back(y);
    pt.pchs.push_back(pch(gi[i]));
    pt.cols.push_back(cl[gi[i]]);
  }
  p.it.push_back(sg);
  p.it.push_back(pt);
  if (grouped && d.key) {
    for (std::size_t g = 0; g < lv.size(); ++g) p.key.push_back({lv[g], cl[g], any ? 1 : 0, pch(g)});
    p.corner = 1;
  }
  return p;
}

// .phist(): the values as a histogram (or a density with its 5% and 95%
// points), a red line at the observed one, notes in the top left corner.
inline Pan hist(Ctx &c, std::vector<double> z, double obs, const std::string &main, const char *xlab, bool density,
                const std::string &fill, const ax::Cols &note = {}) {
  std::sort(z.begin(), z.end());
  Pan p = frame(main.c_str(), xlab, density ? "Density" : "Count");
  std::vector<double> xs = {obs};
  if (density && z.size() > 1 && z.front() != z.back()) {
    std::vector<double> kx, ky;
    Ctx q = c;
    if (plt::dens(q, z, 1, 0.0, kx, ky)) {
      xs.insert(xs.end(), kx.begin(), kx.end());
      rng(xs, p.x0, p.x1);
      p.y1 = hi(ky);
      p.it.push_back(line(kx, ky, fill == "#D5D5D5" ? "#000000" : fill, 1, 2));
      for (double pr : {0.05, 0.95}) p.it.push_back(ref(qtl(z, pr), -p.y1, qtl(z, pr), 3 * p.y1, 2));
    }
  } else {
    const int nb = std::min(30, std::max(5, sat<int>(z.size() / 3)));
    const std::vector<double> b = pretty(z.empty() ? 0 : z.front(), z.empty() ? 1 : z.back(), nb);
    std::vector<double> cnt(b.size() - 1, 0.0);
    for (double v : z) {
      std::size_t j = 0;
      while (j + 2 < b.size() && v > b[j + 1]) ++j;
      cnt[j] += 1;
    }
    xs.insert(xs.end(), b.begin(), b.end());
    rng(xs, p.x0, p.x1);
    p.y1 = hi(cnt) * (note.empty() ? 1.05 : 1.45);
    p.pad = false;
    const double dx = 0.04 * (p.x1 - p.x0);
    p.x0 -= dx;
    p.x1 += dx;
    for (std::size_t j = 0; j < cnt.size(); ++j) {
      p.it.push_back(poly({b[j], b[j + 1], b[j + 1], b[j]}, {0, 0, cnt[j], cnt[j]}, fill, 1, "#FFFFFF"));
    }
  }
  It o = line({obs, obs}, {-p.y1, 3 * p.y1}, "#C0392B", 1, 2);
  p.it.push_back(o);
  for (std::size_t i = 0; i < note.size(); ++i) {
    p.it.push_back(text(p.x0 + 0.03 * (p.x1 - p.x0), p.y1 * (0.95 - 0.075 * static_cast<double>(i)), note[i], 0));
  }
  return p;
}

// The places of the values of an x column: numbers as they are; text as
// the numbers it spells or else its order, shown on the axis as written.
inline std::vector<double> xpos(const Val &col, Ax &a) {
  const std::vector<Val> e = plt::elems(col);
  bool numeric = true;
  for (const Val &v : e) numeric = numeric && (v.t == Val::NUM || v.t == Val::NUL);
  if (numeric) return ax::nums(col);
  ax::Cols lab;
  for (const Val &v : e) lab.push_back(ax::vlab(v));
  const ax::Cols u = lvls(lab);
  std::vector<double> at;
  bool spelt = true;
  for (const std::string &s : u) {
    char *end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    spelt = spelt && !s.empty() && *end == 0;
    at.push_back(v);
  }
  if (!spelt) {
    for (std::size_t i = 0; i < u.size(); ++i) at[i] = static_cast<double>(i) + 1;
  }
  labs(a, u, at, true);
  std::vector<double> x;
  for (std::size_t i : idx(lab, u)) x.push_back(at[i]);
  return x;
}

// .pgrp(): one line per group over x, with its own colour, line type and
// point shape; a band between lo and hi; type 'o' (line and points), 'l'
// (line) or 's' (steps); a dashed line at `ref`, a dotted one at `vref`.
struct Grp {
  std::vector<double> x, y, lo, hi;
  ax::Cols g;
  Ax xa;
  std::string main, xlab, ylab;
  char type = 'o';
  double ref = NA, vref = NA;
  int corner = 1;
  ax::Cols col;
};
inline Pan grp(const Fig &f, const Grp &d) {
  const ax::Cols lv = lvls(d.g), cl = d.col.empty() ? cols(f, sat<int>(lv.size())) : d.col;
  const std::vector<std::size_t> gi = idx(d.g, lv);
  Pan p = frame(d.main.c_str(), d.xlab.c_str(), d.ylab.c_str());
  p.ax = d.xa;
  rng(d.x, p.x0, p.x1);
  std::vector<double> all = d.y;
  all.insert(all.end(), d.lo.begin(), d.lo.end());
  all.insert(all.end(), d.hi.begin(), d.hi.end());
  if (std::isfinite(d.ref)) all.push_back(d.ref);
  rng(all, p.y0, p.y1);
  const double wx = p.x1 - p.x0, wy = p.y1 - p.y0;
  if (std::isfinite(d.ref)) p.it.push_back(ref(p.x0 - wx, d.ref, p.x1 + wx, d.ref, 2));
  if (std::isfinite(d.vref)) p.it.push_back(ref(d.vref, p.y0 - wy, d.vref, p.y1 + wy, 3));
  const bool band = !d.lo.empty();
  for (std::size_t k = 0; k < lv.size(); ++k) {
    std::vector<double> x, y, l, h;
    for (std::size_t i = 0; i < d.x.size(); ++i) {
      if (gi[i] != k) continue;
      x.push_back(d.x[i]);
      y.push_back(d.y[i]);
      if (band) {
        l.push_back(d.lo[i]);
        h.push_back(d.hi[i]);
      }
    }
    const std::string &c0 = cl[k % cl.size()];
    const int lty = sat<int>(k) + 1;
    if (band && d.type == 's') {
      for (const std::vector<double> *e : {&l, &h}) {
        It b = line(x, *e, c0, 3, 1.2);
        b.kind = STEP;
        p.it.push_back(b);
      }
    } else if (band) {
      std::vector<double> qx, qy;
      for (std::size_t i = 0; i < x.size(); ++i) {
        if (std::isfinite(l[i]) && std::isfinite(h[i])) {
          qx.push_back(x[i]);
          qy.push_back(l[i]);
        }
      }
      for (std::size_t i = x.size(); i-- > 0;) {
        if (std::isfinite(l[i]) && std::isfinite(h[i])) {
          qx.push_back(x[i]);
          qy.push_back(h[i]);
        }
      }
      if (!qx.empty()) p.it.push_back(poly(qx, qy, c0, 0.2, "none"));
    }
    It ln = line(x, y, c0, lty, 2);
    if (d.type == 's') ln.kind = STEP;
    p.it.push_back(ln);
    if (d.type == 'o') {
      It pt = line(x, y, c0);
      pt.kind = PTS;
      pt.pch = pch(k);
      p.it.push_back(pt);
    }
    if (lv.size() > 1) p.key.push_back({lv[k], c0, lty, d.type == 'o' ? pch(k) : -1});
  }
  p.corner = d.corner;
  return p;
}
// .plines(): several named series over one x.
inline Grp series(const Val &xcol, const ax::Cols &names, const std::vector<std::vector<double>> &ys) {
  Grp d;
  const std::vector<double> x = xpos(xcol, d.xa);
  for (std::size_t k = 0; k < ys.size(); ++k) {
    for (std::size_t i = 0; i < x.size() && i < ys[k].size(); ++i) {
      d.x.push_back(x[i]);
      d.y.push_back(ys[k][i]);
      d.g.push_back(names[k]);
    }
  }
  return d;
}

// .pbars(): one horizontal bar per label, from 0.
inline Pan bars(const Fig &f, const std::vector<double> &v, const ax::Cols &lab, const std::string &main,
                const char *xlab, double r = NA) {
  const double k = static_cast<double>(v.size());
  Pan p = frame(main.c_str(), xlab, "");
  std::vector<double> all = v;
  all.push_back(0);
  if (std::isfinite(r)) all.push_back(r);
  rng(all, p.x0, p.x1);
  p.y0 = 0.4;
  p.y1 = k + 0.6;
  p.pad = false;
  const double dx = 0.04 * (p.x1 - p.x0);
  p.x0 -= dx;
  p.x1 += dx;
  std::vector<double> at;
  const std::string c0 = cols(f, 1)[0];
  for (std::size_t i = 0; i < v.size(); ++i) {
    const double y = k - static_cast<double>(i);
    at.push_back(y);
    if (std::isfinite(v[i])) {
      It b = poly({0, v[i], v[i], 0}, {y - 0.35, y - 0.35, y + 0.35, y + 0.35}, c0, 1, "none");
      b.lw = 0;
      p.it.push_back(b);
    } else {
      p.it.push_back(text(0, y, "no finite value", 0));
    }
  }
  labs(p.ay, lab, at);
  if (std::isfinite(r)) p.it.push_back(ref(r, 0, r, k + 1, 2));
  return p;
}

// .pmat(): a matrix as a grid of shaded cells with a label in each;
// sequential shades, or diverging ones around `mid`.
inline Pan grid(const Fig &f, const std::vector<std::vector<double>> &m, ax::Cols rl, ax::Cols cn,
                const std::string &main, const char *xlab, const char *ylab, bool div, double mid,
                const std::vector<ax::Cols> &lab = {}) {
  const std::size_t nr = m.size(), nc = nr ? m[0].size() : 0;
  if (rl.size() != nr) {
    rl.clear();
    for (std::size_t i = 0; i < nr; ++i) rl.push_back(std::to_string(i + 1));
  }
  if (cn.size() != nc) {
    cn.clear();
    for (std::size_t j = 0; j < nc; ++j) cn.push_back(std::to_string(j + 1));
  }
  const ax::Cols sh = mp::pal(f.palette == "grey" || f.palette == "okabe" ? (div ? "diverging" : "sequential")
                                                                           : f.palette.c_str(), 9);
  std::vector<double> all;
  for (const std::vector<double> &r : m) all.insert(all.end(), r.begin(), r.end());
  double r0, r1;
  rng(all, r0, r1);
  if (div) {
    const double w = std::max(std::fabs(r0 - mid), std::fabs(r1 - mid));
    r0 = mid - w;
    r1 = mid + w;
  }
  Pan p = frame(main.c_str(), xlab, ylab);
  p.x0 = 0.5;
  p.x1 = static_cast<double>(nc) + 0.5;
  p.y0 = 0.5;
  p.y1 = static_cast<double>(nr) + 0.5;
  p.pad = false;
  std::vector<double> xat, yat;
  for (std::size_t j = 0; j < nc; ++j) xat.push_back(static_cast<double>(j) + 1);
  for (std::size_t i = 0; i < nr; ++i) yat.push_back(static_cast<double>(nr - i));
  labs(p.ax, cn, xat);
  labs(p.ay, rl, yat);
  for (std::size_t i = 0; i < nr; ++i) {
    for (std::size_t j = 0; j < nc; ++j) {
      const double x = static_cast<double>(j) + 1, y = static_cast<double>(nr - i), v = m[i][j];
      const double z = (v - r0) / (r1 - r0);
      const std::string c0 = std::isfinite(z) ? sh[sat<std::size_t>(std::min(8.0, std::max(0.0, std::floor(9 * z))))]
                                              : std::string("#FFFFFF");
      It q = poly({x - 0.5, x + 0.5, x + 0.5, x - 0.5}, {y - 0.5, y - 0.5, y + 0.5, y + 0.5}, c0, 1, "none");
      q.lw = 0;
      p.it.push_back(q);
      const bool dark = std::isfinite(z) && (z > 0.7 || (div && z < 0.3));
      p.it.push_back(text(x, y, lab.empty() ? gfmt(v, 2) : lab[i][j], 1, dark ? "#FFFFFF" : "#000000"));
    }
  }
  return p;
}

// .psurf(): a surface on a grid, z[i][j] at (x[i], y[j]), as shaded cells;
// with `div`, diverging shades around 0.
inline Pan surf(const Fig &f, const std::vector<double> &x, const std::vector<double> &y,
                const std::vector<std::vector<double>> &z, const std::string &main, const char *xlab,
                const char *ylab, bool div = false) {
  const ax::Cols sh = mp::pal(f.palette == "grey" || f.palette == "okabe" ? (div ? "diverging" : "sequential")
                                                                           : f.palette.c_str(), 9);
  std::vector<double> all;
  for (const std::vector<double> &r : z) all.insert(all.end(), r.begin(), r.end());
  double r0, r1;
  rng(all, r0, r1);
  if (div) {
    r1 = std::max(std::fabs(r0), std::fabs(r1));
    r0 = -r1;
  }
  const double hx = x.size() > 1 ? (x[1] - x[0]) / 2 : 0.5, hy = y.size() > 1 ? (y[1] - y[0]) / 2 : 0.5;
  Pan p = frame(main.c_str(), xlab, ylab);
  p.x0 = lo(x) - hx;
  p.x1 = hi(x) + hx;
  p.y0 = lo(y) - hy;
  p.y1 = hi(y) + hy;
  p.pad = false;
  for (std::size_t i = 0; i < x.size() && i < z.size(); ++i) {
    for (std::size_t j = 0; j < y.size() && j < z[i].size(); ++j) {
      if (!std::isfinite(z[i][j])) continue;
      const std::size_t s = sat<std::size_t>(std::min(8.0, std::floor(9 * (z[i][j] - r0) / (r1 - r0))));
      It q = poly({x[i] - hx, x[i] + hx, x[i] + hx, x[i] - hx}, {y[j] - hy, y[j] - hy, y[j] + hy, y[j] + hy}, sh[s],
                  1, sh[s]);
      q.lw = 0.6;
      p.it.push_back(q);
    }
  }
  return p;
}

// A result that lacks what its plot needs is refused in the words of .pck().
inline bool need(In s, const char *cls, std::initializer_list<const char *> keys) {
  std::string all;
  bool ok = s.d.t == Val::OBJ;
  for (const char *k : keys) {
    ok = ok && has(s.d, k);
    all += (all.empty() ? "" : ", ") + std::string(k);
  }
  return ok || s.c.fail("this is not a complete %s result to plot (it needs %s).", cls, all.c_str());
}
inline std::string per(const Val &d) { return "Rate per " + gfmt(num1(get(d, "per")), 7); }
inline ax::Cols seqlab(const char *pre, std::size_t n) {
  ax::Cols r;
  for (std::size_t i = 0; i < n; ++i) r.push_back(std::string(pre) + " " + std::to_string(i + 1));
  return r;
}
inline void push(In s, const Pan &p) { s.f.p.push_back(p); }

// ---- the results drawn as dots with intervals ---------------------------

// A table of the result: label, estimate, lower and upper columns.
inline Dots tab(const Val &t, const char *lab, const char *est, const char *lo_, const char *hi_, const char *main,
                const std::string &xlab, double r = NA) {
  Dots d;
  d.lab = sv_(t, lab);
  d.est = nv(t, est);
  if (lo_ != nullptr) d.lo = nv(t, lo_);
  if (hi_ != nullptr) d.hi = nv(t, hi_);
  d.main = main;
  d.xlab = xlab;
  d.ref = r;
  return d;
}
inline bool rcntreg(In s) {
  CK(need(s, "cntreg", {"irr"}));
  const Val &t = get(s.d, "irr");
  Dots d = tab(t, "term", "irr", "lo", "hi", "Rate ratios with their intervals", "Rate ratio (log scale)", 1);
  if (lvls(sv_(t, "fit")).size() > 1) d.grp = sv_(t, "fit");
  d.logx = true;
  push(s, dots(s.f, d));
  return true;
}
inline bool rdastd(In s) {
  CK(need(s, "dastd", {"rates"}));
  push(s, dots(s.f, tab(get(s.d, "rates"), "group", "rate", "lo", "hi", "Directly standardised rates", per(s.d))));
  return true;
}
inline bool rgrprr(In s) {
  CK(need(s, "grprr", {"rates"}));
  push(s, dots(s.f, tab(get(s.d, "rates"), "group", "rate", "lo", "hi", "Rates by group with their intervals",
                        str(get(s.d, "type")) == "rate" ? per(s.d) : "Share")));
  return true;
}
inline bool rsirex(In s) {
  CK(need(s, "sirex", {"sir"}));
  const Val &t = get(s.d, "sir");
  Dots d = tab(t, "group", "sir", "lo", "hi", "Standardised ratios with their intervals", "Observed / expected", 1);
  for (bool o : bv(t, "out")) d.grp.push_back(o ? "outside the funnel limits" : "inside the funnel limits");
  push(s, dots(s.f, d));
  return true;
}
inline bool rsvymx(In s) {
  CK(need(s, "svymx", {"est"}));
  const std::string st = str(get(s.d, "stat"));
  push(s, dots(s.f, tab(get(s.d, "est"), "dom", "est", "lo", "hi", "Design-based estimates with their intervals",
                        st.empty() ? "Estimate" : st)));
  return true;
}
inline bool rmndla(In s) {
  CK(need(s, "mndla", {"tab"}));
  push(s, dots(s.f, tab(get(s.d, "tab"), "group", "share", "lo", "hi", "Share of placements past the limit",
                        "Share")));
  return true;
}
inline bool rthrtst(In s) {
  CK(need(s, "thrtst", {"thr"}));
  push(s, dots(s.f, tab(get(s.d, "thr"), "group", "thr", "lo", "hi", "Search thresholds by group", "Threshold")));
  return true;
}
// est - z se and est + z se at the 95% level.
inline void zci(Dots &d, const std::vector<double> &se, double z = 1.959963984540054) {
  d.lo.clear();
  d.hi.clear();
  for (std::size_t i = 0; i < d.est.size(); ++i) {
    d.lo.push_back(d.est[i] - z * se[i]);
    d.hi.push_back(d.est[i] + z * se[i]);
  }
}
inline bool rsarml(In s) {
  CK(need(s, "sarml", {"coef"}));
  const Val &t = get(s.d, "coef");
  Dots d = tab(t, "term", "est", nullptr, nullptr, "Coefficients with 95% intervals", "Estimate", 0);
  zci(d, nv(t, "se"));
  push(s, dots(s.f, d));
  return true;
}
inline bool rspcure(In s) {
  CK(need(s, "spcure", {"coef"}));
  const Val &t = get(s.d, "coef");
  Dots d = tab(t, "term", "est", "lo", "hi", "Cure model coefficients with their intervals", "Estimate", 0);
  d.grp = sv_(t, "part");
  for (std::size_t i = 0; i < d.lab.size(); ++i) d.lab[i] = d.grp[i] + ": " + d.lab[i];
  push(s, dots(s.f, d));
  return true;
}
inline bool rmtrbnd(In s) {
  CK(need(s, "mtrbnd", {"bnd"}));
  push(s, dots(s.f, tab(get(s.d, "bnd"), "z", "mean", "lower", "upper",
                        "Observed means (points) and bounds by treatment", "Mean outcome")));
  return true;
}
inline bool rpbbnd(In s) {
  CK(need(s, "pbbnd", {"bnd"}));
  const Val &t = get(s.d, "bnd");
  Dots d = tab(t, "status", "value", "lower", "upper", "Published values (points) and their bounds", "Value");
  d.grp = d.lab;
  d.lab = seqlab("cell", d.est.size());
  push(s, dots(s.f, d));
  return true;
}
inline bool rdkfig(In s) {
  CK(need(s, "dkfig", {"bnd"}));
  const Val &t = get(s.d, "bnd");
  Dots d = tab(t, "v", "v", "lo", "hi", "Recorded rates (points) and bounds on the true rate", "Rate");
  const std::vector<double> r = nv(t, "r");
  for (std::size_t i = 0; i < d.lab.size(); ++i) d.lab[i] = "v " + gfmt(d.est[i], 6) + ", r " + gfmt(r[i], 6);
  push(s, dots(s.f, d));
  return true;
}
// The rows of a table in the order `o`.
inline Dots order(const Dots &d, const std::vector<std::size_t> &o) {
  Dots r = d;
  for (std::size_t i = 0; i < o.size(); ++i) {
    r.lab[i] = d.lab[o[i]];
    r.est[i] = d.est[o[i]];
  }
  return r;
}
inline std::vector<std::size_t> desc(const std::vector<double> &v) {
  std::vector<std::size_t> o(v.size());
  for (std::size_t i = 0; i < o.size(); ++i) o[i] = i;
  std::stable_sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return v[a] > v[b]; });
  return o;
}
inline bool rsmdbal(In s) {
  CK(need(s, "smdbal", {"tab"}));
  Dots d = tab(get(s.d, "tab"), "var", "smd", nullptr, nullptr, "Covariate balance",
               "Absolute standardised difference", num1(get(s.d, "thr")));
  for (double &v : d.est) v = std::fabs(v);
  push(s, dots(s.f, order(d, desc(d.est))));
  return true;
}
inline bool rjmono(In s) {
  CK(need(s, "jmono", {"judges"}));
  const Dots d = tab(get(s.d, "judges"), "judge", "rate", nullptr, nullptr, "Treatment rate by judge", "Rate");
  push(s, dots(s.f, order(d, desc(d.est))));
  return true;
}
inline bool rrelchk(In s) {
  CK(need(s, "relchk", {"reads"}));
  const Val &t = get(s.d, "reads");
  const double v = num1(get(s.d, "value"));
  Dots d = tab(t, "grp", "value", nullptr, nullptr, "", "Value", v);
  d.main = "Readings of the table; the claim is " + gfmt(v, 7);
  const ax::Cols ty = sv_(t, "type"), g = sv_(t, "group"), r = sv_(t, "ref"), e = sv_(t, "event");
  const std::vector<bool> m = bv(t, "match");
  for (std::size_t i = 0; i < d.lab.size(); ++i) {
    d.lab[i] += " " + ty[i] + ": " + g[i] + " vs " + r[i] + ", " + e[i];
    d.grp.push_back(m[i] ? "reproduces the claim" : "does not");
  }
  push(s, dots(s.f, d));
  return true;
}
inline bool rdissix(In s) {
  CK(need(s, "dissix", {"grp"}));
  const Val &t = get(s.d, "grp");
  const ax::Cols g = sv_(t, "group");
  Dots d;
  d.main = "Segregation indices by group";
  d.xlab = "Index";
  const char *const key[] = {"is", "h", "atk", "gini"}, *const nm[] = {"isolation", "entropy", "Atkinson", "Gini"};
  for (int k = 0; k < 4; ++k) {
    const std::vector<double> v = nv(t, key[k]);
    for (std::size_t i = 0; i < g.size(); ++i) {
      d.lab.push_back(g[i]);
      d.est.push_back(v[i]);
      d.grp.push_back(nm[k]);
    }
  }
  push(s, dots(s.f, d));
  return true;
}
inline bool rgmmpnl(In s) {
  CK(need(s, "gmmpnl", {"coef"}));
  const Val &t = get(s.d, "coef");
  const Dots a = tab(t, "term", "est", "lo", "hi", "Dynamic panel GMM", "Coefficient", 0);
  const ax::Cols dm = sv_(s.d, "dummies");
  Dots d = a;
  d.lab.clear();
  d.est.clear();
  d.lo.clear();
  d.hi.clear();
  for (std::size_t i = 0; i < a.lab.size(); ++i) {
    if (a.lab[i] == "(const)" || ax::has(dm, a.lab[i])) continue;
    d.lab.push_back(a.lab[i]);
    d.est.push_back(a.est[i]);
    d.lo.push_back(a.lo[i]);
    d.hi.push_back(a.hi[i]);
  }
  push(s, dots(s.f, d));
  return true;
}
inline bool rccepnl(In s) {
  CK(need(s, "ccepnl", {"mg"}));
  Dots d = tab(get(s.d, "mg"), "term", "est", "lo", "hi", "Common correlated effects", "Coefficient", 0);
  d.grp.assign(d.lab.size(), "CCEMG");
  if (has(s.d, "pooled")) {
    const Dots b = tab(get(s.d, "pooled"), "term", "est", "lo", "hi", "", "");
    d.lab.insert(d.lab.end(), b.lab.begin(), b.lab.end());
    d.est.insert(d.est.end(), b.est.begin(), b.est.end());
    d.lo.insert(d.lo.end(), b.lo.begin(), b.lo.end());
    d.hi.insert(d.hi.end(), b.hi.begin(), b.hi.end());
    d.grp.resize(d.lab.size(), "CCEP");
  }
  push(s, dots(s.f, d));
  return true;
}
inline bool rfepnl(In s) {
  CK(need(s, "fepnl", {"coef"}));
  push(s, dots(s.f, tab(get(s.d, "coef"), "term", "est", "lo", "hi", "Two-way fixed-effects estimates",
                        "Coefficient", 0)));
  if (has(s.d, "slopes")) {
    const Val &t = get(s.d, "slopes");
    Dots d = tab(t, "at", "slope", "lo", "hi", "Simple slopes", "Slope", 0);
    const std::vector<double> at = nv(t, "at");
    for (std::size_t i = 0; i < d.lab.size(); ++i) d.lab[i] = "at " + gfmt(at[i], 3);
    push(s, dots(s.f, d));
  }
  return true;
}
inline bool rhetfx(In s) {
  CK(need(s, "hetfx", {"tab", "test"}));
  push(s, dots(s.f, tab(get(s.d, "tab"), "level", "est", "lo", "hi", "Effect by level of the moderator",
                        "Coefficient (equal effects: p " + gfmt(num1(get(get(s.d, "test"), "p")), 2) + ")", 0)));
  return true;
}
inline bool rtypefx(In s) {
  CK(need(s, "typefx", {"tab"}));
  const Val &t = get(s.d, "tab");
  Dots d = tab(t, "y", "est", "lo", "hi", "Effect by outcome", "Coefficient", 0);
  const std::vector<double> pa = nv(t, "padj");
  for (std::size_t i = 0; i < d.lab.size(); ++i) d.lab[i] += " (adjusted p " + gfmt(pa[i], 2) + ")";
  push(s, dots(s.f, d));
  return true;
}
inline bool rplcday(In s) {
  CK(need(s, "plcday", {"tab"}));
  const Val &t = get(s.d, "tab");
  const Dots a = tab(t, "shift", "est", "lo", "hi", "Placebo: the exposure of other days",
                     "Coefficient of " + str(get(s.d, "var")), 0);
  const std::vector<double> sh = nv(t, "shift");
  std::vector<std::size_t> o(sh.size());
  for (std::size_t i = 0; i < o.size(); ++i) o[i] = i;
  std::stable_sort(o.begin(), o.end(), [&](std::size_t x, std::size_t y) { return sh[x] < sh[y]; });
  Dots d = a;
  d.lo.resize(o.size(), NA);
  d.hi.resize(o.size(), NA);
  for (std::size_t i = 0; i < o.size(); ++i) {
    const double v = sh[o[i]];
    d.lab[i] = v == 0 ? "same day" : gfmt(std::fabs(v), 7) + (v < 0 ? " days before" : " days after");
    d.est[i] = a.est[o[i]];
    d.lo[i] = a.lo[o[i]];
    d.hi[i] = a.hi[o[i]];
    d.grp.push_back(v == 0 ? "same day" : "other days");
  }
  push(s, dots(s.f, d));
  return true;
}

// ---- bars ---------------------------------------------------------------

inline bool rarwt(In s) {
  CK(need(s, "arwt", {"out"}));
  const Val &t = get(s.d, "out");
  push(s, bars(s.f, nv(t, "value"), sv_(t, "to"), "Values moved to each target area", "Value"));
  return true;
}
inline bool rcrsev(In s) {
  CK(need(s, "crsev", {"tab"}));
  const Val &t = get(s.d, "tab");
  push(s, bars(s.f, nv(t, "idx"), sv_(t, "unit"), "Severity-weighted crime index", "Index"));
  return true;
}
inline bool rgmixlr(In s) {
  CK(need(s, "gmixlr", {"locus"}));
  const Val &t = get(s.d, "locus");
  push(s, bars(s.f, nv(t, "log10lr"), sv_(t, "locus"),
               "log10 likelihood ratio by locus (total " + num(num1(get(s.d, "log10lr")), 2) + ")",
               "log10 likelihood ratio", 0));
  return true;
}
inline bool rpedlr(In s) {
  CK(need(s, "pedlr", {"locus"}));
  const Val &t = get(s.d, "locus");
  std::vector<double> v = nv(t, "lr");
  for (double &e : v) e = std::log10(e);
  push(s, bars(s.f, v, sv_(t, "locus"), "log10 likelihood ratio by locus", "log10 likelihood ratio", 0));
  return true;
}
inline bool rhypq(In s) {
  CK(need(s, "hypq", {"wl"}));
  const std::vector<double> v = nv(s.d, "wl");
  push(s, bars(s.f, v, seqlab("unit", v.size()), "Workload of each unit", "Share of time busy"));
  return true;
}
inline bool rlnsize(In s) {
  CK(need(s, "lnsize", {"count"}));
  const std::vector<double> v = nv(s.d, "count");
  ax::Cols lab = seqlab("member", v.size());
  const double su = num1(get(s.d, "susp"));
  if (std::isfinite(su) && su >= 1 && su <= static_cast<double>(v.size())) lab[sat<std::size_t>(su) - 1] += " (suspect)";
  double m = 0;
  for (double e : v) m += e / static_cast<double>(v.size());
  push(s, bars(s.f, v, lab, "Choices of mock witnesses (line: an even spread)", "Choices", m));
  return true;
}
inline bool rlostay(In s) {
  CK(need(s, "lostay", {"bands"}));
  const Val &t = get(s.d, "bands");
  const std::vector<double> l = nv(t, "lo"), h = nv(t, "hi");
  ax::Cols lab;
  for (std::size_t i = 0; i < l.size(); ++i) {
    lab.push_back(std::isfinite(h[i]) ? gfmt(l[i], 6) + " to " + gfmt(h[i], 6) : gfmt(l[i], 6) + " or more");
  }
  push(s, bars(s.f, nv(t, "share"), lab, "Lengths of stay", "Share of stays"));
  return true;
}
inline bool rpinpl(In s) {
  CK(need(s, "pinpl", {"regs"}));
  const Val &t = get(s.d, "regs");
  push(s, bars(s.f, nv(t, "units"), sv_(t, "region"), "Units in each region", "Units"));
  return true;
}
inline bool rlmdix(In s) {
  CK(need(s, "lmdix", {"dec", "factors"}));
  const Val &t = get(s.d, "dec");
  const ax::Cols fa = sv_(s.d, "factors"), un = sv_(t, "unit"), keep = sv_(*s.c.in, "units");
  std::vector<double> v;
  for (const std::string &k : fa) {
    const std::vector<double> col = nv(t, k.c_str());
    double sum = 0;
    for (std::size_t i = 0; i < col.size(); ++i) {
      if (std::isfinite(col[i]) && (keep.empty() || (i < un.size() && ax::has(keep, un[i])))) sum += col[i];
    }
    v.push_back(sum);
  }
  push(s, bars(s.f, v, fa, "LMDI contributions", "Contribution to the change in the count", 0));
  return true;
}

// ---- lines over a period or a quantity ----------------------------------

inline void titled(Grp &d, const std::string &main, const char *xlab, const char *ylab) {
  d.main = main;
  d.xlab = xlab;
  d.ylab = ylab;
}
inline bool rtrndrr(In s) {
  CK(need(s, "trndrr", {"yoy", "fitted", "rr"}));
  const Val &y = get(s.d, "yoy");
  Grp d = series(get(y, "x"), {"observed", "fitted"}, {nv(y, "y"), nv(s.d, "fitted")});
  titled(d, "Trend: rate ratio " + num(nv(get(s.d, "rr"), "rr")[0], 3) + " per period", "Period", "Value");
  push(s, grp(s.f, d));
  return true;
}
inline bool rstkflw(In s) {
  CK(need(s, "stkflw", {"flow"}));
  const Val &t = get(s.d, "flow");
  const char *const key[] = {"people", "days", "adp", "alos"};
  const char *const nm[] = {"People held", "Person-days", "Average daily population", "Average length of stay"};
  const ax::Cols cl = cols(s.f, 4);
  for (int k = 0; k < 4; ++k) {
    Grp d = series(get(t, "period"), {key[k]}, {nv(t, key[k])});
    titled(d, nm[k], "Period", "");
    d.col = {cl[sat<std::size_t>(k)]};
    push(s, grp(s.f, d));
  }
  return true;
}
inline bool rdktflw(In s) {
  CK(need(s, "dktflw", {"flow"}));
  const Val &t = get(s.d, "flow");
  const ax::Cols cl = cols(s.f, 2);
  Grp a = series(get(t, "period"), {"clearance"}, {nv(t, "clear")});
  titled(a, "Clearance rate (disposed / received)", "Period", "Clearance rate");
  a.lo = nv(t, "lo");
  a.hi = nv(t, "hi");
  a.ref = 1;
  a.col = {cl[0]};
  push(s, grp(s.f, a));
  Grp b = series(get(t, "period"), {"pending"}, {nv(t, "pend")});
  titled(b, "Pending cases", "Period", "Cases");
  b.col = {cl[1]};
  push(s, grp(s.f, b));
  return true;
}
inline bool rjpipe(In s) {
  CK(need(s, "jpipe", {"index", "rate"}));
  Grp d;
  d.x = xpos(get(s.d, "index"), d.xa);
  d.y = nv(s.d, "rate");
  d.g = sv_(s.d, "group");
  titled(d, "Rate between stages over time", "Period", "Rate");
  push(s, grp(s.f, d));
  return true;
}
inline bool rhsper(In s) {
  CK(need(s, "hsper", {"periods"}));
  const Val &t = get(s.d, "periods");
  std::vector<double> j;
  for (const Val &e : plt::elems(get(t, "jacc"))) {
    j.push_back(e.t == Val::NUM ? e.n : e.t == Val::STR ? std::strtod(e.s.c_str(), nullptr) : NA);
  }
  Grp d = series(get(t, "period"), {"share of events in hot cells", "hot cells kept (Jaccard)"},
                 {nv(t, "hotshr"), j});
  titled(d, "Hot spots over the periods", "Period", "Share");
  push(s, grp(s.f, d));
  return true;
}
inline bool rcumrsk(In s) {
  CK(need(s, "cumrsk", {"age", "cumrisk"}));
  Grp d = series(get(s.d, "age"), {"cumrisk"}, {nv(s.d, "cumrisk")});
  titled(d, "Cumulative risk by age", "Age", "Cumulative risk");
  d.lo = nv(s.d, "lo");
  d.hi = nv(s.d, "hi");
  push(s, grp(s.f, d));
  return true;
}
inline bool ragemix(In s) {
  CK(need(s, "agemix", {"curve"}));
  const Val &t = get(s.d, "curve");
  ax::Cols nm;
  std::vector<std::vector<double>> ys;
  for (std::size_t i = 0; i < t.k.size(); ++i) {
    if (t.k[i] == "age") continue;
    nm.push_back(t.k[i]);
    ys.push_back(ax::nums(t.a[i]));
  }
  Grp d = series(get(t, "age"), nm, ys);
  titled(d, "Age-crime curve and its offender types", "Age", "Rate");
  push(s, grp(s.f, d));
  return true;
}
inline bool rapcie(In s) {
  CK(need(s, "apcie", {"coef"}));
  const Val &t = get(s.d, "coef");
  const ax::Cols tm = sv_(t, "term"), lv = sv_(t, "level");
  const std::vector<double> e = nv(t, "est"), l = nv(t, "lo"), h = nv(t, "hi");
  for (const std::string &k : lvls(tm)) {
    if (k == "(Intercept)") continue;
    Val xs = varr();
    Grp d;
    for (std::size_t i = 0; i < tm.size(); ++i) {
      if (tm[i] != k) continue;
      xs.a.push_back(vstr(lv[i]));
      d.y.push_back(e[i]);
      d.lo.push_back(l[i]);
      d.hi.push_back(h[i]);
      d.g.push_back("est");
    }
    d.x = xpos(xs, d.xa);
    titled(d, "Effect of " + k, k.c_str(), "Effect");
    d.xlab = k;
    d.ref = 0;
    push(s, grp(s.f, d));
  }
  return true;
}
inline bool rkmrec(In s) {
  CK(need(s, "kmrec", {"km"}));
  const Val &t = get(s.d, "km");
  const ax::Cols g = sv_(t, "grp");
  const std::vector<double> tm = nv(t, "time"), sv = nv(t, "surv"), l = nv(t, "lo"), h = nv(t, "hi");
  Grp d;
  for (const std::string &k : lvls(g)) {
    d.x.push_back(0);
    d.y.push_back(1);
    d.lo.push_back(1);
    d.hi.push_back(1);
    d.g.push_back(k);
    for (std::size_t i = 0; i < g.size(); ++i) {
      if (g[i] != k) continue;
      d.x.push_back(tm[i]);
      d.y.push_back(sv[i]);
      d.lo.push_back(l[i]);
      d.hi.push_back(h[i]);
      d.g.push_back(k);
    }
  }
  titled(d, "Kaplan-Meier survival", "Time", "Share without the event");
  d.type = 's';
  push(s, grp(s.f, d));
  return true;
}
inline bool rlftab(In s) {
  CK(need(s, "lftab", {"tab"}));
  const Val &t = get(s.d, "tab");
  const std::vector<double> lw = nv(t, "lower"), up = nv(t, "upper"), sv = nv(t, "surv"), l = nv(t, "lo"),
                            h = nv(t, "hi");
  Grp d;
  d.x = {lw.empty() ? 0.0 : lw[0]};
  d.y = d.lo = d.hi = {1};
  for (std::size_t i = 0; i < up.size(); ++i) {
    if (!std::isfinite(up[i])) continue;
    d.x.push_back(up[i]);
    d.y.push_back(sv[i]);
    d.lo.push_back(l[i]);
    d.hi.push_back(h[i]);
  }
  d.g.assign(d.x.size(), "survival");
  titled(d, "Life table survival", "Time", "Share surviving");
  d.type = 's';
  push(s, grp(s.f, d));
  return true;
}
inline bool rlnroc(In s) {
  CK(need(s, "lnroc", {"roc"}));
  const Val &t = get(s.d, "roc");
  const ax::Cols g = sv_(t, "proc");
  const std::vector<double> fa = nv(t, "far"), hr = nv(t, "hr");
  Grp d;
  for (const std::string &k : lvls(g)) {
    std::vector<std::pair<double, double>> pts = {{0.0, 0.0}};
    for (std::size_t i = 0; i < g.size(); ++i) {
      if (g[i] == k) pts.push_back({fa[i], hr[i]});
    }
    std::sort(pts.begin(), pts.end());
    for (const std::pair<double, double> &q : pts) {
      d.x.push_back(q.first);
      d.y.push_back(q.second);
      d.g.push_back(k);
    }
  }
  titled(d, "Lineup ROC curves", "False identification rate", "Correct identification rate");
  d.corner = 3;
  push(s, grp(s.f, d));
  return true;
}
inline bool rlnkpr(In s) {
  CK(need(s, "lnkpr", {"pairs"}));
  const Val &t = get(s.d, "pairs");
  const std::vector<double> pr = nv(t, "prob");
  const std::vector<Val> lk = plt::elems(get(t, "linked"));
  std::vector<std::pair<double, bool>> v;
  double np = 0, nn = 0;
  for (std::size_t i = 0; i < pr.size() && i < lk.size(); ++i) {
    if (lk[i].t != Val::BOOL || !std::isfinite(pr[i])) continue;
    v.push_back({pr[i], lk[i].b});
    (lk[i].b ? np : nn) += 1;
  }
  if (np == 0 || nn == 0) return s.c.fail("an ROC curve needs both linked and unlinked pairs.");
  std::sort(v.begin(), v.end(), [](const std::pair<double, bool> &a, const std::pair<double, bool> &b) {
    return a.first > b.first;
  });
  Grp d;
  d.x = d.y = {0};
  double tp = 0, fp = 0, auc = 0;
  for (std::size_t i = 0; i < v.size(); ++i) {
    (v[i].second ? tp : fp) += 1;
    if (i + 1 < v.size() && v[i + 1].first == v[i].first) continue;
    auc += (fp / nn - d.x.back()) * (tp / np + d.y.back()) / 2;
    d.x.push_back(fp / nn);
    d.y.push_back(tp / np);
  }
  d.g.assign(d.x.size(), "model");
  titled(d, "Linkage ROC curve (area " + num(auc, 3) + ")", "Share of unlinked pairs above the threshold",
         "Share of linked pairs above the threshold");
  d.type = 'l';
  Pan p = grp(s.f, d);
  p.it.push_back(ref(0, 0, 1, 1, 3));
  push(s, p);
  return true;
}
inline bool rcconc(In s) {
  CK(need(s, "cconc", {"lorenz", "stat"}));
  const Val &t = get(s.d, "lorenz");
  Pan p = frame(("Lorenz curve (Gini " + num(num1(get(get(s.d, "stat"), "gini")), 3) + ")").c_str(),
                "Cumulative share of units, smallest first", "Cumulative share of the total");
  p.it.push_back(ref(0, 0, 1, 1, 2));
  p.it.push_back(line(nv(t, "pop"), nv(t, "share"), cols(s.f, 1)[0], 1, 2.2));
  push(s, p);
  return true;
}
inline bool rcrossk(In s) {
  CK(need(s, "crossk", {"k"}));
  const Val &t = get(s.d, "k");
  Grp d = series(get(t, "r"), {"estimated", "under independence"}, {nv(t, "k"), nv(t, "theo")});
  titled(d, "Cross K function, " + str(get(s.d, "i")) + " to " + str(get(s.d, "j")), "Distance", "K(r)");
  d.corner = 0;
  push(s, grp(s.f, d));
  return true;
}
inline bool rjdgiv(In s) {
  CK(need(s, "jdgiv", {"sens"}));
  const Val &t = get(s.d, "sens");
  Grp d;
  d.x = nv(t, "dsh");
  d.y = nv(t, "cshare");
  for (double v : nv(t, "deff")) d.g.push_back("direct effect " + gfmt(v, 7));
  titled(d, "Complier share by the share of defiers", "Share of defiers", "Complier share");
  push(s, grp(s.f, d));
  return true;
}
inline bool rrptvic(In s) {
  CK(need(s, "rptvic", {"tab"}));
  const Val &t = get(s.d, "tab");
  Grp d = series(get(t, "k"), {"observed", "expected by chance"}, {nv(t, "victims"), nv(t, "expected")});
  titled(d, "Targets by times victimised", "Times victimised", "Targets");
  push(s, grp(s.f, d));
  return true;
}

// ---- grids --------------------------------------------------------------

inline bool rctabx(In s) {
  CK(need(s, "ctabx", {"tab", "res"}));
  std::vector<ax::Cols> lab;
  for (const std::vector<double> &r : mat(get(s.d, "tab"))) {
    lab.push_back({});
    for (double v : r) lab.back().push_back(gfmt(v, 7));
  }
  push(s, grid(s.f, mat(get(s.d, "res")), {}, {}, "Residuals of the table (cells show counts)", "", "", true, 0,
               lab));
  return true;
}
inline bool rknxnr(In s) {
  CK(need(s, "knxnr", {"tab"}));
  const Val &t = get(s.d, "tab");
  const std::vector<double> dl = nv(t, "dlo"), dh = nv(t, "dhi"), tl = nv(t, "tlo"), th = nv(t, "thi"),
                            ra = nv(t, "ratio"), pv = nv(t, "p");
  ax::Cols rs, cs;
  for (std::size_t i = 0; i < ra.size(); ++i) {
    rs.push_back(gfmt(dl[i], 7) + " to " + gfmt(dh[i], 7));
    cs.push_back(gfmt(tl[i], 7) + " to " + gfmt(th[i], 7));
  }
  const ax::Cols ru = lvls(rs), cu = lvls(cs);
  std::vector<std::vector<double>> m(ru.size(), std::vector<double>(cu.size(), NA));
  std::vector<ax::Cols> lab(ru.size(), ax::Cols(cu.size(), ""));
  const std::vector<std::size_t> ri = idx(rs, ru), ci = idx(cs, cu);
  for (std::size_t i = 0; i < ra.size(); ++i) {
    m[ri[i]][ci[i]] = ra[i];
    lab[ri[i]][ci[i]] = num(ra[i], 2) + (pv[i] < 0.05 ? "*" : "");
  }
  push(s, grid(s.f, m, ru, cu, "Knox ratios of close pairs (* p < 0.05)", "Time apart", "Distance apart", true, 1,
               lab));
  return true;
}
inline bool rqdrtst(In s) {
  CK(need(s, "qdrtst", {"counts"}));
  push(s, grid(s.f, mat(get(s.d, "counts")), {}, {}, "Counts by quadrat", "Column", "Row", false, 0));
  return true;
}
inline bool rtrnmat(In s) {
  CK(need(s, "trnmat", {"prob"}));
  const ax::Cols st = get(s.d, "pi").t == Val::OBJ ? get(s.d, "pi").k : ax::Cols();
  push(s, grid(s.f, mat(get(s.d, "prob")), st, st, "Transition probabilities, row to column", "To", "From", false,
               0));
  return true;
}
inline bool rzcomp(In s) {
  CK(need(s, "zcomp", {"z"}));
  push(s, grid(s.f, mat(get(s.d, "z")), sv_(get(s.d, "tab"), "unit"), sv_(get(s.d, "cat"), "cat"),
               "Standardised counts by unit and category", "Category", "Unit", true, 0));
  return true;
}
// A surface result: x, y and the matrix z.
inline bool surfof(In s, const char *cls, const char *main, const Val &z, const char *xlab = "x",
                   const char *ylab = "y") {
  const std::vector<std::vector<double>> m = mat(z);
  std::vector<double> x = nv(s.d, "x"), y = nv(s.d, "y");
  if (m.empty() || m[0].size() < 2) {
    return s.c.fail("this %s result holds no grid to draw (its z is not a matrix).", cls);
  }
  if (x.size() != m.size()) {
    x.clear();
    for (std::size_t i = 0; i < m.size(); ++i) x.push_back(static_cast<double>(i) + 1);
  }
  if (y.size() != m[0].size()) {
    y.clear();
    for (std::size_t j = 0; j < m[0].size(); ++j) y.push_back(static_cast<double>(j) + 1);
  }
  push(s, surf(s.f, x, y, m, main, xlab, ylab));
  return true;
}
inline bool rkdsurf(In s) {
  CK(need(s, "kdsurf", {"x", "y", "z"}));
  return surfof(s, "kdsurf", "Kernel intensity surface", get(s.d, "z"));
}
inline bool rtpsurf(In s) {
  CK(need(s, "tpsurf", {"x", "y", "z"}));
  return surfof(s, "tpsurf", "Thin-plate spline surface", get(s.d, "z"));
}
inline bool rsdbsim(In s) {
  CK(need(s, "sdbsim", {"b"}));
  return surfof(s, "sdbsim", "Attractiveness at the end of the run", get(s.d, "b"), "column", "row");
}
inline bool rprohs(In s) {
  CK(need(s, "prohs", {"x", "y", "z", "hot"}));
  CK(surfof(s, "prohs", "Prospective risk (hot cells outlined)", get(s.d, "z")));
  Pan &p = s.f.p.back();
  const std::vector<double> x = nv(s.d, "x"), y = nv(s.d, "y");
  const double hx = x.size() > 1 ? (x[1] - x[0]) / 2 : 0.5, hy = y.size() > 1 ? (y[1] - y[0]) / 2 : 0.5;
  const std::vector<Val> rows_ = plt::elems(get(s.d, "hot"));
  for (std::size_t i = 0; i < rows_.size() && i < x.size(); ++i) {
    const std::vector<Val> r = plt::elems(rows_[i]);
    for (std::size_t j = 0; j < r.size() && j < y.size(); ++j) {
      if (r[j].t != Val::BOOL || !r[j].b) continue;
      It q = poly({x[i] - hx, x[i] + hx, x[i] + hx, x[i] - hx}, {y[j] - hy, y[j] - hy, y[j] + hy, y[j] + hy},
                  "none", 1, "#252525");
      q.lw = 1.5;
      p.it.push_back(q);
    }
  }
  return true;
}

// ---- scatters and the rest ----------------------------------------------

inline It points(const std::vector<double> &x, const std::vector<double> &y, const std::string &col, int shape,
                 double r = 3.6) {
  It pt = line(x, y, col);
  pt.kind = PTS;
  pt.pch = shape;
  pt.r = r;
  return pt;
}
inline bool rxcser(In s) {
  CK(need(s, "xcser", {"ccf", "cor"}));
  const Val &t = get(s.d, "ccf");
  const std::vector<double> lag = nv(t, "lag"), r = nv(t, "r");
  // the normal quantile of the level, from the interval of the correlation
  const double lev = num1(get(s.d, "level")), n = nv(get(s.d, "cor"), "n")[0];
  const double b = clab_qnorm(1 - (1 - lev) / 2) / std::sqrt(n);
  const std::string c0 = cols(s.f, 1)[0];
  Pan p = frame("Cross-correlation at lags", "Lag", "Correlation");
  rng(lag, p.x0, p.x1);
  std::vector<double> all = r;
  all.insert(all.end(), {b, -b});
  rng(all, p.y0, p.y1);
  const double wx = p.x1 - p.x0;
  p.it.push_back(ref(p.x0 - wx, 0, p.x1 + wx, 0, 1));
  for (double v : {-b, b}) p.it.push_back(ref(p.x0 - wx, v, p.x1 + wx, v, 2));
  It sg;
  sg.kind = SEG;
  sg.col = c0;
  sg.lw = 3;
  for (std::size_t i = 0; i < lag.size(); ++i) {
    sg.x.insert(sg.x.end(), {lag[i], lag[i]});
    sg.y.insert(sg.y.end(), {0.0, r[i]});
  }
  p.it.push_back(sg);
  p.it.push_back(points(lag, r, c0, 0));
  push(s, p);
  return true;
}
inline bool rhwkres(In s) {
  CK(need(s, "hwkres", {"u", "n"}));
  std::vector<double> u = nv(s.d, "u"), q;
  std::sort(u.begin(), u.end());
  const double n = static_cast<double>(u.size()), b = 1.36 / std::sqrt(n);
  for (std::size_t i = 0; i < u.size(); ++i) q.push_back((static_cast<double>(i) + 0.5) / n);
  Pan p = frame("Time-rescaling residuals (95% Kolmogorov band)", "uniform quantile", "rescaled time quantile");
  p.it.push_back(line({-1, 2}, {-1, 2}, "#000000", 1, 1));
  p.it.push_back(line({-1, 2}, {-1 + b, 2 + b}, "#000000", 2, 1));
  p.it.push_back(line({-1, 2}, {-1 - b, 2 - b}, "#000000", 2, 1));
  p.it.push_back(points(q, u, cols(s.f, 1)[0], 2));
  push(s, p);
  return true;
}
inline bool rmmshr(In s) {
  CK(need(s, "mmshr", {"group", "below", "at", "above"}));
  const ax::Cols g = sv_(s.d, "group"), cl = cols(s.f, 3);
  const char *const key[] = {"below", "at", "above"};
  const double k = static_cast<double>(g.size());
  Pan p = frame("Sentences against the minimum", "Share of sentences", "");
  p.x1 = 1;
  p.y0 = 0.4;
  p.y1 = k + 0.6 + 0.9;
  p.pad = false;
  std::vector<double> at, left(g.size(), 0.0);
  for (std::size_t i = 0; i < g.size(); ++i) at.push_back(k - static_cast<double>(i));
  double top = 1;
  for (int j = 0; j < 3; ++j) {
    const std::vector<double> v = nv(s.d, key[j]);
    for (std::size_t i = 0; i < g.size(); ++i) {
      if (!std::isfinite(v[i])) continue;
      const double y = at[i];
      p.it.push_back(poly({left[i], left[i] + v[i], left[i] + v[i], left[i]},
                          {y - 0.35, y - 0.35, y + 0.35, y + 0.35}, cl[sat<std::size_t>(j)], j == 0 ? 1 : 0.55,
                          "#252525"));
      left[i] += v[i];
      top = std::max(top, left[i]);
    }
    p.key.push_back({key[j], cl[sat<std::size_t>(j)], 1, -1});
  }
  p.x1 = top;
  p.corner = 1;
  labs(p.ay, g, at);
  push(s, p);
  return true;
}
inline bool rbufas(In s) {
  CK(need(s, "bufas", {"pts"}));
  const Val &t = get(s.d, "pts");
  const ax::Cols g = sv_(t, "grp"), lv = lvls(g), cl = cols(s.f, sat<int>(lv.size()));
  const std::vector<double> x = nv(t, "dsite"), y = nv(t, "dcomp");
  Pan p = frame("Places by distance to sites and comparison sites", "Distance to the nearest site",
                "Distance to the nearest comparison site");
  rng(x, p.x0, p.x1);
  rng(y, p.y0, p.y1);
  const Val &pa = get(s.d, "par");
  const double wy = p.y1 - p.y0;
  int lty = 2;
  for (const char *k : {"r", "buf"}) {
    const double v = num1(get(pa, k));
    if (std::isfinite(v)) p.it.push_back(ref(v, p.y0 - wy, v, p.y1 + wy, lty));
    ++lty;
  }
  for (std::size_t k = 0; k < lv.size(); ++k) {
    p.it.push_back(points(rows(x, g, lv[k]), rows(y, g, lv[k]), cl[k], pch(k)));
    p.key.push_back({lv[k], cl[k], 0, pch(k)});
  }
  p.corner = 1;
  push(s, p);
  return true;
}
// Points in two classes: discs for the marked ones, rings for the others.
inline void twoclass(Pan &p, const std::vector<double> &x, const std::vector<double> &y, const std::vector<bool> &m,
                     const ax::Cols &cl, int shape, const char *yes, const char *no, int corner) {
  for (int k = 1; k >= 0; --k) {
    std::vector<double> qx, qy;
    for (std::size_t i = 0; i < x.size(); ++i) {
      if ((i < m.size() && m[i]) == (k == 1)) {
        qx.push_back(x[i]);
        qy.push_back(y[i]);
      }
    }
    p.it.push_back(points(qx, qy, cl[sat<std::size_t>(k)], k ? shape : 2, 4));
    p.key.push_back({k ? yes : no, cl[sat<std::size_t>(k)], 0, k ? shape : 2});
  }
  p.corner = corner;
}
inline bool rcmpsc(In s) {
  CK(need(s, "cmpsc", {"segs"}));
  const Val &t = get(s.d, "segs");
  const std::vector<double> x = nv(t, "seg"), y = nv(t, "ccp");
  const double pos = num1(get(s.d, "pos"));
  char ttl[96];
  std::snprintf(ttl, sizeof ttl, "Segments congruent at the best shift: %d of %d", sat<int>(num1(get(s.d, "score"))),
                sat<int>(num1(get(s.d, "nseg"))));
  Pan p = frame(ttl, "Segment", "Position of the best match");
  rng(x, p.x0, p.x1);
  std::vector<double> all = y;
  all.push_back(pos);
  rng(all, p.y0, p.y1);
  const double wx = p.x1 - p.x0;
  p.it.push_back(ref(p.x0 - wx, pos, p.x1 + wx, pos, 2));
  twoclass(p, x, y, bv(t, "cmps"), cols(s.f, 2), 0, "congruent", "not congruent", 1);
  push(s, p);
  return true;
}
inline bool rebshr(In s) {
  CK(need(s, "ebshr", {"est"}));
  const Val &t = get(s.d, "est");
  const std::vector<double> y = nv(t, "y"), sh = nv(t, "shrunk");
  const ax::Cols cl = cols(s.f, 2);
  Pan p = frame("Counts and their empirical-Bayes shrinkage", "", "Value");
  p.x0 = 0.8;
  p.x1 = 2.2;
  std::vector<double> all = y;
  all.insert(all.end(), sh.begin(), sh.end());
  rng(all, p.y0, p.y1);
  labs(p.ax, {"observed", "shrunk"}, {1, 2});
  It sg;
  sg.kind = SEG;
  sg.col = "#969696";
  sg.lw = 1;
  for (std::size_t i = 0; i < y.size(); ++i) {
    sg.x.insert(sg.x.end(), {1.0, 2.0});
    sg.y.insert(sg.y.end(), {y[i], sh[i]});
  }
  p.it.push_back(sg);
  p.it.push_back(points(std::vector<double>(y.size(), 1.0), y, cl[0], 2));
  p.it.push_back(points(std::vector<double>(y.size(), 2.0), sh, cl[1], 0));
  push(s, p);
  return true;
}
inline bool ribench(In s) {
  CK(need(s, "ibench", {"tab"}));
  const Val &t = get(s.d, "tab");
  const std::vector<double> x = nv(t, "bench"), y = nv(t, "rate");
  const std::vector<bool> fl = bv(t, "flag");
  const ax::Cols off = sv_(t, "off");
  Pan p = frame("Officer rate against benchmark", "Benchmark", "Rate");
  std::vector<double> all = x;
  all.insert(all.end(), y.begin(), y.end());
  rng(all, p.x0, p.x1);
  p.y0 = p.x0;
  p.y1 = p.x1;
  const double w = p.x1 - p.x0;
  p.it.push_back(ref(p.x0 - w, p.x0 - w, p.x1 + w, p.x1 + w, 2));
  twoclass(p, x, y, fl, cols(s.f, 2), 1, "flagged", "not flagged", 0);
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (i < fl.size() && fl[i]) p.it.push_back(text(x[i], y[i] + 0.045 * w, off[i], 1));
  }
  push(s, p);
  return true;
}
inline bool rkscan(In s) {
  CK(need(s, "kscan", {"cl"}));
  const Val &t = get(s.d, "cl");
  const std::vector<double> x = nv(t, "x"), y = nv(t, "y"), r = nv(t, "r");
  const std::string c0 = cols(s.f, 1)[0];
  Pan p = frame("Space-time clusters (numbered by rank)", "x", "y");
  std::vector<double> ax_, ay_;
  for (std::size_t i = 0; i < x.size(); ++i) {
    ax_.insert(ax_.end(), {x[i] - r[i], x[i] + r[i]});
    ay_.insert(ay_.end(), {y[i] - r[i], y[i] + r[i]});
  }
  rng(ax_, p.x0, p.x1);
  rng(ay_, p.y0, p.y1);
  // one scale on both axes, so that a circle is round
  const double w = std::max(p.x1 - p.x0, (p.y1 - p.y0) * 1.5);
  const double mx = (p.x0 + p.x1) / 2, my = (p.y0 + p.y1) / 2;
  p.x0 = mx - w / 2;
  p.x1 = mx + w / 2;
  p.y0 = my - w / 3;
  p.y1 = my + w / 3;
  const double pi = 3.14159265358979323846;
  for (std::size_t i = 0; i < x.size(); ++i) {
    std::vector<double> qx, qy;
    for (int k = 0; k < 72; ++k) {
      qx.push_back(x[i] + r[i] * std::cos(2 * pi * k / 72));
      qy.push_back(y[i] + r[i] * std::sin(2 * pi * k / 72));
    }
    It q = poly(qx, qy, "none", 1, c0);
    q.lw = 2;
    p.it.push_back(q);
    p.it.push_back(text(x[i], y[i], std::to_string(i + 1), 1));
  }
  push(s, p);
  return true;
}

// ---- the policy results with more than one view -------------------------

inline bool rdidfe(In s) {
  CK(need(s, "didfe", {"coef"}));
  const Val &t = get(s.d, "coef");
  const std::vector<double> est = nv(t, "est"), lo_ = nv(t, "lo"), hi_ = nv(t, "hi");
  std::string ty = "dots";
  if (!str1((*s.c.in).t == Val::OBJ ? s.c["type"] : nul, ty) || (ty != "dots" && ty != "hist")) {
    return s.c.fail("`type` must be one of \"dots\", \"hist\".");
  }
  const bool loo = has(s.d, "loo");
  if (!loo) push(s, dots(s.f, tab(t, "term", "est", "lo", "hi", "Two-way fixed-effects estimate", "Effect", 0)));
  if (has(s.d, "bacon")) {
    const Val &b = get(s.d, "bacon");
    const ax::Cols g = sv_(b, "type"), lv = lvls(g), cl = cols(s.f, sat<int>(lv.size()));
    const std::vector<double> w = nv(b, "weight"), e = nv(b, "est");
    Pan p = frame("Bacon decomposition", "Weight", "2 x 2 difference-in-differences estimate");
    std::vector<double> ax_ = w, ay_ = e;
    ax_.push_back(0);
    ay_.push_back(est[0]);
    rng(ax_, p.x0, p.x1);
    rng(ay_, p.y0, p.y1);
    const double wx = p.x1 - p.x0, wy = p.y1 - p.y0;
    p.x1 += 0.06 * wx;
    p.y0 -= 0.08 * wy;
    p.y1 += 0.2 * static_cast<double>(lv.size() + 1) * wy;
    p.it.push_back(line({p.x0 - wx, p.x1 + wx}, {est[0], est[0]}, "#000000", 1, 1.5));
    for (std::size_t k = 0; k < lv.size(); ++k) {
      p.it.push_back(points(rows(w, g, lv[k]), rows(e, g, lv[k]), cl[k], pch(k), 4.2));
      p.key.push_back({lv[k], cl[k], 0, pch(k)});
    }
    p.key.push_back({"Two-way estimate (line) = " + gfmt(est[0], 4), "#000000", 1, -1});
    p.corner = 0;
    push(s, p);
  }
  if (loo) {
    const Val &l = get(s.d, "loo");
    const std::vector<double> le = nv(l, "est"), ls = nv(l, "se");
    if (ty == "hist") {
      std::vector<double> e;
      double m = 0, sd = 0;
      for (double v : le) {
        if (std::isfinite(v)) e.push_back(v);
      }
      for (double v : e) m += v / static_cast<double>(e.size());
      for (double v : e) sd += (v - m) * (v - m);
      sd = e.size() > 1 ? std::sqrt(sd / static_cast<double>(e.size() - 1)) : NA;
      push(s, hist(s.c, e, est[0], "Leave-one-out sensitivity", "Estimate leaving one unit out", false, "#D5D5D5",
                   {"N = " + std::to_string(e.size()),
                    "Range: [" + gfmt(lo(e), 3) + ", " + gfmt(hi(e), 3) + "]",
                    "Mean " + gfmt(m, 3) + ", SD " + gfmt(sd, 2), "Full sample = " + gfmt(est[0], 3)}));
    } else {
      const double z = clab_qnorm(1 - (1 - num1(get(s.d, "level"))) / 2);
      Dots d;
      d.main = "Estimate, each unit left out";
      d.xlab = "Effect";
      d.ref = est[0];
      d.key = false;
      d.lab = {"all units"};
      d.est = {est[0]};
      d.lo = {lo_[0]};
      d.hi = {hi_[0]};
      d.grp = {"all units"};
      const ax::Cols un = sv_(l, "unit");
      for (std::size_t i = 0; i < le.size(); ++i) {
        d.lab.push_back("without " + un[i]);
        d.est.push_back(le[i]);
        d.lo.push_back(le[i] - z * ls[i]);
        d.hi.push_back(le[i] + z * ls[i]);
        d.grp.push_back("one left out");
      }
      push(s, dots(s.f, d));
    }
  }
  return true;
}
inline bool rwcbx(In s) {
  CK(need(s, "wcbx", {"res", "mthd"}));
  const Val &r = get(s.d, "res");
  const std::string mt = str(get(s.d, "mthd"));
  std::string ty = "hist";
  if (!str1(s.c["type"], ty) || (ty != "hist" && ty != "density")) {
    return s.c.fail("`type` must be one of \"hist\", \"density\".");
  }
  if (mt == "perm") {
    const std::vector<std::vector<double>> nu = mat(get(s.d, "null"));
    const ax::Cols out = sv_(r, "outcome");
    const std::vector<double> st = nv(r, "stat"), pv = nv(r, "p");
    if (nu.empty() || nu[0].size() != st.size()) return s.c.fail("this wcbx result holds no permuted statistics.");
    const std::string fill = s.f.palette == "grey" ? "#D5D5D5" : cols(s.f, 1)[0];
    for (std::size_t j = 0; j < st.size(); ++j) {
      std::vector<double> z;
      for (const std::vector<double> &row : nu) {
        if (std::isfinite(row[j])) z.push_back(row[j]);
      }
      push(s, hist(s.c, z, st[j], out[j] + " (p = " + gfmt(pv[j], 3) + ")", "Permuted statistic", ty == "density",
                   fill));
    }
    return true;
  }
  if (mt != "power" || !has(r, "eff") || !has(r, "power")) {
    return s.c.fail("a wcbx result is drawn for mthd = \"perm\" (the permuted statistics) and \"power\" (the "
                    "power at each effect); this one is \"%s\".", mt.c_str());
  }
  const std::vector<double> ef = nv(r, "eff"), pw = nv(r, "power");
  const std::string c0 = cols(s.f, 1)[0];
  Pan p = frame("Power of the design", "Effect size", "Power");
  rng(ef, p.x0, p.x1);
  p.it.push_back(line(ef, pw, c0, 1, 2));
  p.it.push_back(points(ef, pw, c0, 0));
  const double mde = num1(get(s.d, "mde"));
  if (std::isfinite(mde)) {
    const double wx = p.x1 - p.x0;
    for (std::size_t i = 0; i + 1 < ef.size(); ++i) {
      if (mde < ef[i] || mde > ef[i + 1] || ef[i + 1] == ef[i]) continue;
      const double v = pw[i] + (pw[i + 1] - pw[i]) * (mde - ef[i]) / (ef[i + 1] - ef[i]);
      p.it.push_back(ref(p.x0 - wx, v, p.x1 + wx, v, 2));
      break;
    }
    p.it.push_back(ref(mde, -1, mde, 2, 3));
  }
  push(s, p);
  return true;
}
inline bool ryoypx(In s) {
  std::string what = "pctchg";
  if (!str1(s.c["what"], what) || (what != "pctchg" && what != "value")) {
    return s.c.fail("`what` must be one of \"pctchg\", \"value\".");
  }
  CK(need(s, "yoypx", {"prd", what == "value" ? "value" : "pctchg"}));
  const std::vector<double> prd = nv(s.d, "prd");
  ax::Cols g = sv_(s.d, "group");
  if (g.size() != prd.size()) g.assign(prd.size(), "all");
  if (what == "value") {
    const std::vector<double> v = nv(s.d, "value");
    Grp d;
    for (std::size_t i = 0; i < prd.size(); ++i) {
      if (!std::isfinite(v[i])) continue;
      d.x.push_back(prd[i]);
      d.y.push_back(v[i]);
      d.g.push_back(g[i]);
    }
    titled(d, "Value by period", "Period", "Value");
    push(s, grp(s.f, d));
    return true;
  }
  const std::vector<double> pc = nv(s.d, "pctchg"), l = nv(s.d, "pctlo"), h = nv(s.d, "pctup");
  ax::Cols gg;
  std::vector<double> x, y, ql, qh;
  for (std::size_t i = 0; i < prd.size(); ++i) {
    if (!std::isfinite(pc[i])) continue;
    x.push_back(prd[i]);
    y.push_back(pc[i]);
    ql.push_back(l[i]);
    qh.push_back(h[i]);
    gg.push_back(g[i]);
  }
  const ax::Cols lv = lvls(gg), cl = cols(s.f, sat<int>(lv.size()));
  const std::vector<std::size_t> gi = idx(gg, lv);
  const double k = static_cast<double>(lv.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    x[i] += (static_cast<double>(gi[i]) + 1 - (k + 1) / 2) * 0.6 / std::max(k, 1.0);
  }
  Pan p = frame("Percent change from the period before", "Period", "Percent change");
  rng(x, p.x0, p.x1);
  std::vector<double> all = y;
  all.insert(all.end(), ql.begin(), ql.end());
  all.insert(all.end(), qh.begin(), qh.end());
  all.push_back(0);
  rng(all, p.y0, p.y1);
  const double wx = p.x1 - p.x0;
  p.it.push_back(ref(p.x0 - wx, 0, p.x1 + wx, 0, 2));
  It sg, pt;
  sg.kind = SEG;
  sg.cap = 4;
  sg.lw = 2;
  pt.kind = PTS;
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (std::isfinite(ql[i]) && std::isfinite(qh[i])) {
      sg.x.insert(sg.x.end(), {x[i], x[i]});
      sg.y.insert(sg.y.end(), {ql[i], qh[i]});
      sg.cols.insert(sg.cols.end(), {cl[gi[i]], cl[gi[i]]});
    }
    pt.x.push_back(x[i]);
    pt.y.push_back(y[i]);
    pt.pchs.push_back(pch(gi[i]));
    pt.cols.push_back(cl[gi[i]]);
  }
  p.it.push_back(sg);
  p.it.push_back(pt);
  if (lv.size() > 1) {
    for (std::size_t j = 0; j < lv.size(); ++j) p.key.push_back({lv[j], cl[j], 0, pch(j)});
  }
  p.corner = 1;
  push(s, p);
  return true;
}
