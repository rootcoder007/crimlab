// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// plotsvg() for clab_jrun(): the plots of R/plotx.R drawn outside R. It
// runs the plot function named by `plot` on `args` (that function's own
// arguments, checked by it) and draws the rows it returns the way the R
// function draws them (for evtst, stgdid, synthc and sdidx: the result,
// the way plot() of it draws it), as one standalone SVG document: a white page, a
// black box, a bold title, the palette of plotopt() (grey by default),
// no grid unless asked for. The result is the document (`svg`) and the
// rows drawn (`data`). Coordinates are written by integer arithmetic and
// labels with the decimal point ".", so the text does not depend on the
// locale.

namespace sv {

enum Kind { LINE, STEP, POLY, PTS, SEG, TXT };

// One thing drawn, in data coordinates. SEG: pairs of points, one segment
// each (cap > 0: a bar of that half-width in pixels across both ends; hz
// for a segment that runs across the page).
// PTS: pch 0 a disc, 1 a triangle, 2 a ring, 3 a square, 4 a diamond.
// TXT: anc 0 starts at the point, 1 is centred on it, 2 ends at it; with
// `edge` it starts right of the box, level with the point.
struct It {
  Kind kind = LINE;
  std::vector<double> x, y;
  std::string col = "#000000", fill = "none", txt;
  double lw = 1.6, op = 1, r = 3.4, cap = 0;
  int lty = 1, pch = 0, anc = 0;
  bool edge = false, hz = false;
  std::vector<int> pchs;
  ax::Cols cols;
};
struct Key {
  std::string lab, col;
  int lty, pch;
};
// An axis: `own` gives the places and their labels, `none` draws nothing,
// `ticks` false leaves the labels without tick marks.
struct Ax {
  std::vector<double> at;
  ax::Cols lab;
  bool own = false, none = false, ticks = true;
};
struct Pan {
  std::string main, xlab, ylab;
  double x0 = 0, x1 = 1, y0 = 0, y1 = 1;
  Ax ax, ay;
  std::vector<It> it;
  std::vector<Key> key;
  int corner = 0;
  double mr = 0;      // room right of the box, for labels at the line ends
  bool box = true, pad = true;
};
struct Fig {
  std::string main, palette = "grey";
  std::vector<Pan> p;
  bool grid = false;
  double w = 720, h = 480;
};

// A pixel coordinate with one decimal.
inline std::string px(double v) {
  const long long t = std::llround(v * 10);
  const long long a = t < 0 ? -t : t;
  std::string s = (t < 0 ? "-" : "") + std::to_string(a / 10);
  if (a % 10) s += "." + std::to_string(a % 10);
  return s;
}
// A number with `dec` decimals, "." the decimal point in every locale.
inline std::string num(double v, int dec) {
  char b[64];
  std::snprintf(b, sizeof b, "%.*f", dec, v);
  std::string s = b;
  for (char &ch : s) {
    if (ch == ',') ch = '.';
  }
  if (s == "-0") s = "0";
  return s;
}
inline std::string esc(const std::string &t) {
  std::string s;
  for (char ch : t) {
    if (ch == '&') s += "&amp;";
    else if (ch == '<') s += "&lt;";
    else if (ch == '>') s += "&gt;";
    else if (ch == '"') s += "&quot;";
    else s += ch;
  }
  return s;
}
// R's line types 1 to 6 as dash patterns.
inline const char *dash(int lty) {
  static const char *const d[] = {"", "7 5", "2 4", "2 4 7 4", "11 5", "8 3 3 3"};
  return d[sat<std::size_t>((lty < 1 ? 1 : lty) - 1) % 6];
}
// Ticks at 1, 2 or 5 times a power of 10, about five of them.
inline void ticks(double lo, double hi, Ax &a) {
  if (hi < lo) std::swap(lo, hi);
  const double raw = (hi - lo) / 5;
  double st = std::pow(10.0, std::floor(std::log10(raw)));
  const double e = raw / st;
  st *= e >= 7.5 ? 10 : e >= 3.5 ? 5 : e >= 1.5 ? 2 : 1;
  const int dec = std::max(0, sat<int>(-std::floor(std::log10(st) + 1e-9)));
  for (double v = std::ceil(lo / st - 1e-9) * st; v <= hi + st * 1e-9; v += st) {
    a.at.push_back(std::fabs(v) < st * 1e-9 ? 0.0 : v);
    a.lab.push_back(num(a.at.back(), dec));
  }
}

// The box of one panel in pixels and the map from data to it.
struct Box {
  double l, t, r, b, x0, x1, y0, y1;
  double X(double v) const { return l + (v - x0) / (x1 - x0) * (r - l); }
  double Y(double v) const { return b - (v - y0) / (y1 - y0) * (b - t); }
};

inline void mark(std::string &o, int pch, double x, double y, double r, const std::string &col) {
  const std::string c = " stroke=\"" + col + "\"", f = " fill=\"" + col + "\"";
  if (pch == 1) {
    o += "<path d=\"M" + px(x) + " " + px(y - r * 1.25) + "L" + px(x + r * 1.15) + " " + px(y + r * 0.85) + "L" +
         px(x - r * 1.15) + " " + px(y + r * 0.85) + "Z\"" + f + "/>";
  } else if (pch == 3) {
    o += "<rect x=\"" + px(x - r) + "\" y=\"" + px(y - r) + "\" width=\"" + px(2 * r) + "\" height=\"" + px(2 * r) +
         "\"" + f + "/>";
  } else if (pch == 4) {
    o += "<path d=\"M" + px(x) + " " + px(y - r * 1.3) + "L" + px(x + r * 1.3) + " " + px(y) + "L" + px(x) + " " +
         px(y + r * 1.3) + "L" + px(x - r * 1.3) + " " + px(y) + "Z\"" + f + "/>";
  } else {
    o += "<circle cx=\"" + px(x) + "\" cy=\"" + px(y) + "\" r=\"" + px(r) + "\"" +
         (pch == 2 ? " fill=\"#FFFFFF\" stroke-width=\"1.5\"" + c : f) + "/>";
  }
}

inline void path(std::string &o, const It &i, const Box &b) {
  std::string d;
  bool up = true;
  for (std::size_t k = 0; k < i.x.size(); ++k) {
    if (!std::isfinite(i.x[k]) || !std::isfinite(i.y[k])) {
      up = true;
      continue;
    }
    if (!up && i.kind == STEP) d += "H" + px(b.X(i.x[k])) + "V" + px(b.Y(i.y[k]));
    else d += (up ? "M" : "L") + px(b.X(i.x[k])) + " " + px(b.Y(i.y[k]));
    up = false;
  }
  if (d.empty()) return;
  if (i.kind == POLY) d += "Z";
  o += "<path d=\"" + d + "\" fill=\"" + i.fill + "\"";
  if (i.kind == POLY && i.op < 1) o += " fill-opacity=\"" + num(i.op, 2) + "\"";
  o += " stroke=\"" + i.col + "\" stroke-width=\"" + num(i.lw, 1) + "\"";
  if (*dash(i.lty)) o += std::string(" stroke-dasharray=\"") + dash(i.lty) + "\"";
  o += " stroke-linejoin=\"round\"/>";
}

inline void item(std::string &o, const It &i, const Box &b) {
  if (i.kind == LINE || i.kind == STEP || i.kind == POLY) return path(o, i, b);
  for (std::size_t k = 0; k < i.x.size(); ++k) {
    if (!std::isfinite(i.x[k]) || !std::isfinite(i.y[k])) continue;
    const std::string &col = i.cols.empty() ? i.col : i.cols[k % i.cols.size()];
    const double x = i.edge ? b.r + 8 : b.X(i.x[k]), y = b.Y(i.y[k]);
    if (i.kind == PTS) mark(o, i.pchs.empty() ? i.pch : i.pchs[k % i.pchs.size()], x, y, i.r, col);
    if (i.kind == TXT) {
      o += "<text x=\"" + px(x) + "\" y=\"" + px(y + 4) + "\" font-size=\"12\" fill=\"" + col + "\"" +
           (i.anc == 1 ? " text-anchor=\"middle\"" : i.anc == 2 ? " text-anchor=\"end\"" : "") + ">" + esc(i.txt) +
           "</text>";
    }
    if (i.kind == SEG && k % 2 == 0 && k + 1 < i.x.size()) {
      const double x2 = b.X(i.x[k + 1]), y2 = b.Y(i.y[k + 1]);
      std::string d = "M" + px(x) + " " + px(y) + "L" + px(x2) + " " + px(y2);
      if (i.cap > 0 && i.hz) {
        d += "M" + px(x) + " " + px(y - i.cap) + "V" + px(y + i.cap) + "M" + px(x2) + " " + px(y2 - i.cap) + "V" +
             px(y2 + i.cap);
      } else if (i.cap > 0) {
        d += "M" + px(x - i.cap) + " " + px(y) + "H" + px(x + i.cap) + "M" + px(x2 - i.cap) + " " + px(y2) + "H" +
             px(x2 + i.cap);
      }
      o += "<path d=\"" + d + "\" fill=\"none\" stroke=\"" + col + "\" stroke-width=\"" + num(i.lw, 1) + "\"";
      if (*dash(i.lty)) o += std::string(" stroke-dasharray=\"") + dash(i.lty) + "\"";
      o += "/>";
    }
  }
}

inline void axes(std::string &o, const Pan &p, const Box &b, bool grid) {
  const std::string tk = " stroke=\"#000000\" stroke-width=\"1\"/>";
  for (std::size_t k = 0; k < p.ax.at.size() && !p.ax.none; ++k) {
    const double x = b.X(p.ax.at[k]);
    if (x < b.l - 0.5 || x > b.r + 0.5) continue;
    if (grid) o += "<path d=\"M" + px(x) + " " + px(b.t) + "V" + px(b.b) + "\" stroke=\"#DDDDDD\" stroke-width=\"0.8\"/>";
    if (p.ax.ticks) o += "<path d=\"M" + px(x) + " " + px(b.b) + "v6\"" + tk;
    o += "<text x=\"" + px(x) + "\" y=\"" + px(b.b + 21) + "\" text-anchor=\"middle\">" + esc(p.ax.lab[k]) + "</text>";
  }
  for (std::size_t k = 0; k < p.ay.at.size() && !p.ay.none; ++k) {
    const double y = b.Y(p.ay.at[k]);
    if (y < b.t - 0.5 || y > b.b + 0.5) continue;
    if (grid) o += "<path d=\"M" + px(b.l) + " " + px(y) + "H" + px(b.r) + "\" stroke=\"#DDDDDD\" stroke-width=\"0.8\"/>";
    if (p.ay.ticks) o += "<path d=\"M" + px(b.l) + " " + px(y) + "h-6\"" + tk;
    o += "<text x=\"" + px(b.l - 10) + "\" y=\"" + px(y + 4) + "\" text-anchor=\"end\">" + esc(p.ay.lab[k]) + "</text>";
  }
}

inline void legend(std::string &o, const Pan &p, const Box &b) {
  if (p.key.empty()) return;
  std::size_t wd = 0;
  for (const Key &k : p.key) wd = std::max(wd, k.lab.size());
  const double w = 44 + 6.4 * static_cast<double>(wd), h = 17.0 * static_cast<double>(p.key.size());
  const double x = p.corner % 2 ? b.r - w - 6 : b.l + 10, y0 = p.corner >= 2 ? b.b - h - 6 : b.t + 8;
  for (std::size_t k = 0; k < p.key.size(); ++k) {
    const Key &e = p.key[k];
    const double y = y0 + 17.0 * static_cast<double>(k) + 9;
    if (e.lty > 0) {
      o += "<path d=\"M" + px(x) + " " + px(y) + "h26\" stroke=\"" + e.col + "\" stroke-width=\"2\"";
      if (*dash(e.lty)) o += std::string(" stroke-dasharray=\"") + dash(e.lty) + "\"";
      o += "/>";
    }
    if (e.pch >= 0) mark(o, e.pch, x + 13, y, 3.4, e.col);
    o += "<text x=\"" + px(x + 34) + "\" y=\"" + px(y + 4) + "\" font-size=\"12\">" + esc(e.lab) + "</text>";
  }
}

// One panel in the cell (cx, cy, cw, ch); `id` names its clip path.
inline void panel(std::string &o, Pan p, double cx, double cy, double cw, double ch, int id, bool grid) {
  if (p.pad) {
    const double dx = 0.04 * (p.x1 - p.x0), dy = 0.04 * (p.y1 - p.y0);
    p.x0 -= dx;
    p.x1 += dx;
    p.y0 -= dy;
    p.y1 += dy;
  }
  if (!p.ax.own && !p.ax.none) ticks(p.x0, p.x1, p.ax);
  if (!p.ay.own && !p.ay.none) ticks(p.y0, p.y1, p.ay);
  std::size_t wd = 0;
  for (const std::string &s : p.ay.lab) wd = std::max(wd, s.size());
  const double ml = (p.ay.none ? 14 : 20 + 6.6 * static_cast<double>(wd)) + (p.ylab.empty() ? 0 : 20);
  const Box b = {cx + ml, cy + 34, cx + cw - 14 - p.mr, cy + ch - (p.xlab.empty() ? 30 : 48), p.x0, p.x1, p.y0, p.y1};
  const std::string cid = "c" + std::to_string(id);
  o += "<clipPath id=\"" + cid + "\"><rect x=\"" + px(b.l) + "\" y=\"" + px(b.t) + "\" width=\"" + px(b.r - b.l) +
       "\" height=\"" + px(b.b - b.t) + "\"/></clipPath>";
  axes(o, p, b, grid);
  o += "<g clip-path=\"url(#" + cid + ")\">";
  for (const It &i : p.it) {
    if (i.kind != TXT) item(o, i, b);
  }
  o += "</g>";
  for (const It &i : p.it) {
    if (i.kind == TXT) item(o, i, b);
  }
  if (p.box) {
    o += "<rect x=\"" + px(b.l) + "\" y=\"" + px(b.t) + "\" width=\"" + px(b.r - b.l) + "\" height=\"" +
         px(b.b - b.t) + "\" fill=\"none\" stroke=\"#000000\" stroke-width=\"1\"/>";
  }
  legend(o, p, b);
  const double mx = (b.l + b.r) / 2, my = (b.t + b.b) / 2;
  o += "<text x=\"" + px(mx) + "\" y=\"" + px(cy + 22) + "\" text-anchor=\"middle\" font-size=\"14\" "
       "font-weight=\"bold\">" + esc(p.main) + "</text>";
  if (!p.xlab.empty()) {
    o += "<text x=\"" + px(mx) + "\" y=\"" + px(b.b + 40) + "\" text-anchor=\"middle\">" + esc(p.xlab) + "</text>";
  }
  if (!p.ylab.empty()) {
    o += "<text transform=\"translate(" + px(cx + 15) + " " + px(my) + ") rotate(-90)\" text-anchor=\"middle\">" +
         esc(p.ylab) + "</text>";
  }
}

// The document: the panels in a grid, the title of the figure above them.
inline std::string draw(const Fig &f) {
  const int n = sat<int>(f.p.size()), nc = sat<int>(std::ceil(std::sqrt(static_cast<double>(n))));
  const int nr = (n + nc - 1) / nc;
  const double top = f.main.empty() || n == 1 ? 0 : 26;
  std::string o = "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 " + px(f.w) + " " + px(f.h) + "\" width=\"" +
                  px(f.w) + "\" height=\"" + px(f.h) + "\" font-family=\"Helvetica, Arial, sans-serif\" "
                  "font-size=\"12.5\" fill=\"#000000\"><rect width=\"100%\" height=\"100%\" fill=\"#FFFFFF\"/>";
  for (int i = 0; i < n; ++i) {
    const double cw = f.w / nc, ch = (f.h - top) / nr;
    panel(o, f.p[sat<std::size_t>(i)], cw * (i % nc), top + ch * (i / nc), cw, ch, i, f.grid);
  }
  if (top > 0) {
    o += "<text x=\"" + px(f.w / 2) + "\" y=\"19\" text-anchor=\"middle\" font-size=\"15\" font-weight=\"bold\">" +
         esc(f.main) + "</text>";
  }
  return o + "</svg>";
}

// ---- reading arguments and rows ----------------------------------------

inline const Val &get(const Val &o, const char *k) {
  for (std::size_t i = 0; i < o.k.size(); ++i) {
    if (o.k[i] == k) return o.a[i];
  }
  return nul;
}
inline std::vector<double> nv(const Val &o, const char *k) { return ax::nums(get(o, k)); }
inline ax::Cols sv_(const Val &o, const char *k) {
  ax::Cols r;
  for (const Val &e : plt::elems(get(o, k))) r.push_back(ax::vlab(e));
  return r;
}
// The values of the rows whose label is `l`.
inline std::vector<double> rows(const std::vector<double> &v, const ax::Cols &g, const std::string &l) {
  std::vector<double> r;
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (g[i] == l) r.push_back(v[i]);
  }
  return r;
}
inline ax::Cols lvls(const ax::Cols &g) {
  ax::Cols l;
  for (const std::string &s : g) {
    if (!ax::has(l, s)) l.push_back(s);
  }
  return l;
}
inline double lo(const std::vector<double> &v) {
  double r = HUGE_VAL;
  for (double e : v) {
    if (std::isfinite(e)) r = std::min(r, e);
  }
  return r;
}
inline double hi(const std::vector<double> &v) {
  double r = -HUGE_VAL;
  for (double e : v) {
    if (std::isfinite(e)) r = std::max(r, e);
  }
  return r;
}
// .prange(): the range of the finite values, (0, 1) when there is none, a
// single value widened by 0.5 either way.
inline void rng(const std::vector<double> &v, double &a, double &b) {
  a = lo(v);
  b = hi(v);
  if (!(a <= b)) {
    a = 0;
    b = 1;
  }
  if (b - a == 0) {
    a -= 0.5;
    b += 0.5;
  }
}
inline It line(std::vector<double> x, std::vector<double> y, const std::string &col, int lty = 1, double lw = 2) {
  It i;
  i.x = std::move(x);
  i.y = std::move(y);
  i.col = col;
  i.lty = lty;
  i.lw = lw;
  return i;
}
inline It poly(std::vector<double> x, std::vector<double> y, const std::string &fill, double op,
               const std::string &border) {
  It i = line(std::move(x), std::move(y), border, 1, 1);
  i.kind = POLY;
  i.fill = fill;
  i.op = op;
  return i;
}
inline It text(double x, double y, const std::string &t, int anc, const std::string &col = "#000000") {
  It i;
  i.kind = TXT;
  i.x = {x};
  i.y = {y};
  i.txt = t;
  i.anc = anc;
  i.col = col;
  return i;
}
// A straight reference line across the panel.
inline It ref(double x0, double y0, double x1, double y1, int lty) { return line({x0, x1}, {y0, y1}, "#636363", lty, 1); }
inline ax::Cols cols(const Fig &f, int n) { return mp::pal(f.palette.c_str(), std::max(n, 1)); }
// An axis of labels at 1, 2, ... or at the places given.
inline void labs(Ax &a, const ax::Cols &l, const std::vector<double> &at, bool tk = false) {
  a.own = true;
  a.ticks = tk;
  a.at = at;
  a.lab = l;
}

// One text option, or "" when it is left out; false when it is not text.
inline bool str1(const Val &o, std::string &s) {
  if (o.t == Val::NUL) return true;
  const Val *e = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  if (e->t != Val::STR) return false;
  s = e->s;
  return true;
}

// What a builder is given: the arguments of the plot function, the rows it
// returned, and the figure to fill.
struct In {
  Ctx &c;
  const Val &a, &d;
  Fig &f;
};

inline Pan frame(const char *main, const char *xlab, const char *ylab) {
  Pan p;
  p.main = main;
  p.xlab = xlab;
  p.ylab = ylab;
  return p;
}

// cconc()'s Gini coefficient of the values x (NA when it cannot be had).
inline double gini(Ctx &c, const Val &x) {
  Ctx q = c;
  Val r;
  if (!plt::pcc(q, x, r)) return NA;
  const Val &g = get(get(r, "stat"), "gini");
  return g.t == Val::NUM ? g.n : NA;
}

inline bool blrnz(In s) {
  const double g = gini(s.c, get(s.a, "x"));
  char t[64];
  std::snprintf(t, sizeof t, "Lorenz curve (Gini %s)", num(g, 3).c_str());
  Pan p = frame(std::isfinite(g) ? t : "Lorenz curve", "Cumulative share of units, smallest first",
                "Cumulative share of the total");
  p.it.push_back(ref(0, 0, 1, 1, 2));
  p.it.push_back(line(nv(s.d, "pop"), nv(s.d, "share"), cols(s.f, 1)[0], 1, 2.2));
  s.f.p.push_back(p);
  return true;
}

inline bool bgain(In s) {
  const ax::Cols g = sv_(s.d, "grp"), lv = lvls(g), cl = cols(s.f, sat<int>(lv.size()));
  const std::vector<double> pop = nv(s.d, "pop"), sh = nv(s.d, "share");
  Pan p = frame("Concentration curves", "Cumulative share of units, smallest first",
                "Cumulative share of the total");
  p.it.push_back(ref(0, 0, 1, 1, 2));
  Ctx q = s.c;
  q.in = &s.a;
  plt::Xg d;
  CK(plt::pxg(q, d));
  for (std::size_t i = 0; i < lv.size(); ++i) {
    const int lty = sat<int>(i) + 1;
    p.it.push_back(line(rows(pop, g, lv[i]), rows(sh, g, lv[i]), cl[i], lty, 2.2));
    p.key.push_back({lv[i] + " (Gini " + num(gini(s.c, vnums(plt::ingrp(d, lv[i]))), 3) + ")", cl[i], lty, -1});
  }
  s.f.p.push_back(p);
  return true;
}

inline bool bcdfs(In s) {
  const ax::Cols g = sv_(s.d, "grp"), lv = lvls(g), cl = cols(s.f, sat<int>(lv.size()));
  const std::vector<double> x = nv(s.d, "x"), y = nv(s.d, "y");
  Pan p = frame("Empirical distribution functions", "Value", "Share at or below");
  p.x0 = lo(x);
  p.x1 = hi(x);
  for (std::size_t i = 0; i < lv.size(); ++i) {
    std::vector<double> v = rows(x, g, lv[i]), q = rows(y, g, lv[i]);
    v.insert(v.begin(), v[0]);
    q.insert(q.begin(), 0.0);
    It l = line(v, q, cl[i], sat<int>(i) + 1, 2);
    l.kind = STEP;
    p.it.push_back(l);
    if (lv.size() > 1) p.key.push_back({lv[i], cl[i], sat<int>(i) + 1, -1});
  }
  p.corner = 3;
  s.f.p.push_back(p);
  return true;
}

// The groups of `x`/`grp` as the plot functions read them.
inline bool groups(In s, plt::Xg &d) {
  Ctx q = s.c;
  q.in = &s.a;
  if (plt::pxg(q, d)) return true;
  s.c.err = q.err;
  return false;
}

inline bool brdg(In s) {
  plt::Xg d;
  CK(groups(s, d));
  const ax::Cols g = sv_(s.d, "grp"), cl = cols(s.f, sat<int>(d.lv.size()));
  const std::vector<double> x = nv(s.d, "x"), y = nv(s.d, "y");
  const double top = hi(y), step = 0.7 * top;
  Pan p = frame("Distribution by group", "Value", "");
  p.x0 = lo(x);
  p.x1 = hi(x);
  p.y1 = step * static_cast<double>(d.lv.size() - 1) + top;
  std::vector<double> at;
  for (std::size_t k = d.lv.size(); k-- > 0;) {
    std::vector<double> vx = rows(x, g, d.lv[k]), vy = rows(y, g, d.lv[k]);
    const double b = step * static_cast<double>(k);
    for (double &e : vy) e += b;
    if (vx.empty()) continue;
    std::vector<double> qx = vx, qy = vy;
    qx.push_back(vx.back());
    qy.push_back(b);
    qx.push_back(vx.front());
    qy.push_back(b);
    p.it.push_back(poly(qx, qy, cl[k], 0.5, "#FFFFFF"));
    p.it.push_back(line(vx, vy, cl[k], sat<int>(k) + 1, 1.8));
  }
  for (std::size_t k = 0; k < d.lv.size(); ++k) at.push_back(step * static_cast<double>(k));
  labs(p.ay, d.lv, at);
  s.f.p.push_back(p);
  return true;
}

// R's quantile() (type 7) of sorted values.
inline double qtl(const std::vector<double> &v, double pr) {
  const double h = static_cast<double>(v.size() - 1) * pr, fl = std::floor(h);
  const std::size_t i = sat<std::size_t>(fl);
  return i + 1 < v.size() ? v[i] + (h - fl) * (v[i + 1] - v[i]) : v[i];
}

inline bool bvln(In s) {
  plt::Xg d;
  CK(groups(s, d));
  const ax::Cols g = sv_(s.d, "grp"), cl = cols(s.f, sat<int>(d.lv.size()));
  const std::vector<double> x = nv(s.d, "x"), y = nv(s.d, "y");
  Pan p = frame("Violins: density, median and quartiles", "", "Value");
  p.x0 = 0.5;
  p.x1 = static_cast<double>(d.lv.size()) + 0.5;
  // the data and the tails of the densities
  p.y0 = std::min(lo(d.x), lo(x));
  p.y1 = std::max(hi(d.x), hi(x));
  std::vector<double> at;
  for (std::size_t k = 0; k < d.lv.size(); ++k) {
    const double c0 = static_cast<double>(k) + 1;
    at.push_back(c0);
    const std::vector<double> kx = rows(x, g, d.lv[k]), ky = rows(y, g, d.lv[k]);
    if (!kx.empty()) {
      const double mx = hi(ky);
      std::vector<double> qx, qy;
      for (std::size_t i = 0; i < kx.size(); ++i) {
        qx.push_back(c0 - 0.42 * ky[i] / mx);
        qy.push_back(kx[i]);
      }
      for (std::size_t i = kx.size(); i-- > 0;) {
        qx.push_back(c0 + 0.42 * ky[i] / mx);
        qy.push_back(kx[i]);
      }
      p.it.push_back(poly(qx, qy, cl[k], 0.5, cl[k]));
    }
    std::vector<double> v = plt::ingrp(d, d.lv[k]);
    std::sort(v.begin(), v.end());
    p.it.push_back(line({c0, c0}, {qtl(v, 0.25), qtl(v, 0.75)}, "#252525", 1, 5));
    It m;
    m.kind = PTS;
    m.x = {c0};
    m.y = {qtl(v, 0.5)};
    m.pch = 2;
    m.col = "#252525";
    m.r = 4;
    p.it.push_back(m);
  }
  labs(p.ax, d.lv, at);
  s.f.p.push_back(p);
  return true;
}

// The whole number `k` of the arguments, or its default.
inline double opt(const Val &a, const char *k, double def) {
  const std::vector<double> v = nv(a, k);
  return v.size() == 1 && std::isfinite(v[0]) ? v[0] : def;
}

inline bool bhex(In s) {
  const std::vector<double> ax_ = nv(s.a, "x"), ay_ = nv(s.a, "y"), x = nv(s.d, "x"), y = nv(s.d, "y"),
                            n = nv(s.d, "count");
  std::vector<double> fx, fy;
  for (std::size_t i = 0; i < ax_.size() && i < ay_.size(); ++i) {
    if (std::isfinite(ax_[i]) && std::isfinite(ay_[i])) {
      fx.push_back(ax_[i]);
      fy.push_back(ay_[i]);
    }
  }
  const double bins = opt(s.a, "bins", 30), mx = hi(n);
  char t[64];
  std::snprintf(t, sizeof t, "Counts in hexagons (most: %d)", sat<int>(mx));
  Pan p = frame(t, "x", "y");
  rng(fx, p.x0, p.x1);
  rng(fy, p.y0, p.y1);
  const ax::Cols sh = mp::pal(s.f.palette == "grey" ? "sequential" : s.f.palette.c_str(), 9);
  const double pi = 3.14159265358979323846;
  for (std::size_t i = 0; i < x.size(); ++i) {
    std::vector<double> qx, qy;
    for (int k = 0; k < 6; ++k) {
      const double an = pi / 6 + k * pi / 3;
      qx.push_back(x[i] + std::cos(an) / std::sqrt(3.0) / bins * (p.x1 - p.x0));
      qy.push_back(y[i] + std::sin(an) / std::sqrt(3.0) / bins * (p.y1 - p.y0));
    }
    const int lv = std::min(9, 1 + sat<int>(std::floor(8 * (n[i] - 1) / std::max(1.0, mx - 1))));
    It h = poly(qx, qy, sh[sat<std::size_t>(lv - 1)], 1, "none");
    h.lw = 0;
    p.it.push_back(h);
  }
  s.f.p.push_back(p);
  return true;
}

inline bool bfnl(In s) {
  const std::vector<double> ex = nv(s.d, "expected"), ra = nv(s.d, "ratio");
  const ax::Cols ar = sv_(s.d, "area");
  const std::vector<Val> ou = plt::elems(get(s.d, "outside"));
  std::vector<double> lv = nv(s.a, "lvls"), e = ex, all = ra;
  if (lv.empty()) lv = {0.95, 0.998};
  std::sort(lv.begin(), lv.end());
  const double tg = opt(s.a, "target", 1);
  std::sort(e.begin(), e.end());
  e.erase(std::unique(e.begin(), e.end()), e.end());
  Pan p = frame("Funnel plot", "Expected count", "Observed / expected");
  rng(e, p.x0, p.x1);
  for (std::size_t j = 0; j < lv.size(); ++j) {
    std::vector<double> bl, bh;
    for (double v : e) {
      const double al = (1 - lv[j]) / 2;
      bl.push_back(clab_qpois(al, tg * v, 1) / v);
      bh.push_back(clab_qpois(1 - al, tg * v, 1) / v);
    }
    all.insert(all.end(), bl.begin(), bl.end());
    all.insert(all.end(), bh.begin(), bh.end());
    const int lty = j + 1 == lv.size() ? 1 : 2;
    p.it.push_back(line(e, bl, "#636363", lty, 1));
    p.it.push_back(line(e, bh, "#636363", lty, 1));
  }
  all.push_back(tg);
  rng(all, p.y0, p.y1);
  p.it.push_back(ref(p.x0 - (p.x1 - p.x0), tg, p.x1 + (p.x1 - p.x0), tg, 1));
  const ax::Cols c2 = mp::pal(s.f.palette == "grey" ? "grey" : "diverging", s.f.palette == "grey" ? 2 : 3);
  It pt;
  pt.kind = PTS;
  pt.x = ex;
  pt.y = ra;
  pt.r = 4;
  for (std::size_t i = 0; i < ex.size(); ++i) {
    const bool out = i < ou.size() && ou[i].t == Val::BOOL && ou[i].b;
    pt.pchs.push_back(out ? 1 : 0);
    pt.cols.push_back(out ? c2.back() : c2[0]);
    if (out && i < ar.size() && get(s.a, "area").t != Val::NUL) {
      p.it.push_back(text(ex[i], ra[i] + 0.045 * (p.y1 - p.y0), ar[i], 1));
    }
  }
  p.it.push_back(pt);
  s.f.p.push_back(p);
  return true;
}

inline bool bcal(In s) {
  const std::vector<double> yr = nv(s.d, "year"), wk = nv(s.d, "week"), wd = nv(s.d, "wday"), y = nv(s.d, "y");
  std::vector<double> yrs = yr;
  std::sort(yrs.begin(), yrs.end());
  yrs.erase(std::unique(yrs.begin(), yrs.end()), yrs.end());
  double r0, r1;
  rng(y, r0, r1);
  const ax::Cols sh = mp::pal(s.f.palette == "grey" ? "sequential" : s.f.palette.c_str(), 9);
  Pan p = frame(("Daily values, " + num(r0, r0 == std::floor(r0) ? 0 : 2) + " to " +
                 num(r1, r1 == std::floor(r1) ? 0 : 2) + " (darker is higher)").c_str(), "", "");
  p.x0 = 0;
  p.x1 = 54;
  p.y1 = 8 * static_cast<double>(yrs.size());
  p.box = p.pad = false;
  std::vector<double> yat;
  ax::Cols ylab;
  for (std::size_t i = 0; i < yrs.size(); ++i) {
    const double top = 8 * static_cast<double>(yrs.size() - i);
    yat.push_back(top - 4);
    ylab.push_back(num(yrs[i], 0));
    for (std::size_t k = 0; k < y.size(); ++k) {
      if (yr[k] != yrs[i]) continue;
      const double a = wk[k], b = top - wd[k];
      const std::string c0 = std::isfinite(y[k]) ? sh[sat<std::size_t>(std::min(8.0, std::floor(8 * (y[k] - r0) / (r1 - r0))))]
                                                 : std::string("#E0E0E0");
      It q = poly({a, a + 0.92, a + 0.92, a}, {b, b, b + 0.92, b + 0.92}, c0, 1, "none");
      q.lw = 0;
      p.it.push_back(q);
    }
  }
  labs(p.ay, ylab, yat);
  labs(p.ax, {"Jan", "Apr", "Jul", "Oct", "Dec"}, {0, 13, 26, 39, 52});
  s.f.p.push_back(p);
  return true;
}

// The series of `y`/`prd`/`grp` as the plot functions read them, and the
// x axis: the periods as written when they are text.
inline bool series(In s, plt::Ser &d) {
  Ctx q = s.c;
  q.in = &s.a;
  if (plt::pser(q, d)) return true;
  s.c.err = q.err;
  return false;
}
inline void prdax(const plt::Ser &d, Pan &p) {
  p.x0 = lo(d.pos);
  p.x1 = hi(d.pos);
  if (p.x1 - p.x0 == 0) {
    p.x0 -= 0.5;
    p.x1 += 0.5;
  }
  bool txt = false;
  for (const Val &e : d.prd) txt = txt || e.t == Val::STR;
  if (!txt) return;
  std::vector<double> at;
  ax::Cols lab;
  for (std::size_t i = 0; i < d.pos.size(); ++i) {
    if (std::find(at.begin(), at.end(), d.pos[i]) != at.end()) continue;
    at.push_back(d.pos[i]);
    lab.push_back(ax::vlab(d.prd[i]));
  }
  labs(p.ax, lab, at, true);
}
// The rows of series k: positions and values.
inline void ofser(const plt::Ser &d, int k, const std::vector<double> &v, std::vector<double> &x,
                  std::vector<double> &y) {
  for (std::size_t i = 0; i < d.gi.size(); ++i) {
    if (d.gi[i] == k) {
      x.push_back(d.pos[i]);
      y.push_back(v[i]);
    }
  }
}

inline bool bmult(In s) {
  plt::Ser d;
  CK(series(s, d));
  const ax::Cols cl = cols(s.f, sat<int>(d.lv.size()));
  double y0, y1;
  rng(d.y, y0, y1);
  for (std::size_t k = 0; k < d.lv.size(); ++k) {
    Pan p = frame(d.lv[k].c_str(), "", "");
    prdax(d, p);
    p.y0 = y0;
    p.y1 = y1;
    std::vector<double> x, y;
    ofser(d, sat<int>(k), d.y, x, y);
    p.it.push_back(line(x, y, cl[k], 1, 2));
    It pt = line(x, y, cl[k]);
    pt.kind = PTS;
    p.it.push_back(pt);
    s.f.p.push_back(p);
  }
  return true;
}

inline bool bbump(In s) {
  plt::Ser d;
  CK(series(s, d));
  const std::vector<double> rk = nv(s.d, "rank");
  const ax::Cols cl = cols(s.f, sat<int>(d.lv.size()));
  const double nr = hi(rk);
  Pan p = frame("Ranks over time", "Period", "Rank");
  prdax(d, p);
  p.y0 = nr + 0.5;
  p.y1 = 0.5;
  p.pad = false;
  const double dx = 0.04 * (p.x1 - p.x0);
  p.x0 -= dx;
  p.x1 += dx;
  std::size_t wd = 0;
  for (const std::string &l : d.lv) wd = std::max(wd, l.size());
  p.mr = 14 + 6.8 * static_cast<double>(wd);
  std::vector<double> at;
  ax::Cols lab;
  for (int i = 1; i <= sat<int>(nr); ++i) {
    at.push_back(i);
    lab.push_back(std::to_string(i));
  }
  labs(p.ay, lab, at);
  for (std::size_t k = 0; k < d.lv.size(); ++k) {
    std::vector<double> x, y;
    ofser(d, sat<int>(k), rk, x, y);
    if (x.empty()) continue;
    p.it.push_back(line(x, y, cl[k], sat<int>(k) + 1, 3));
    It pt = line(x, y, cl[k]);
    pt.kind = PTS;
    pt.pch = sat<int>(k % 5);
    pt.r = 4.4;
    p.it.push_back(pt);
    It lb = text(x.back(), y.back(), d.lv[k], 0, cl[k]);
    lb.edge = true;
    p.it.push_back(lb);
  }
  s.f.p.push_back(p);
  return true;
}

inline bool bhrz(In s) {
  plt::Ser d;
  CK(series(s, d));
  const int nb = sat<int>(opt(s.a, "bands", 3)), ng = sat<int>(d.lv.size());
  std::vector<double> dev(d.y.size(), NA);
  double top = 0;
  for (int k = 0; k < ng; ++k) {
    std::vector<double> x, y;
    ofser(d, k, d.y, x, y);
    const double m = plt::med(y);
    for (std::size_t i = 0; i < d.gi.size(); ++i) {
      if (d.gi[i] != k || !std::isfinite(d.y[i])) continue;
      dev[i] = d.y[i] - m;
      top = std::max(top, std::fabs(dev[i]));
    }
  }
  if (top == 0) top = 1;
  const double bw = top / nb;
  ax::Cols up = mp::pal("sequential", nb + 1);
  up.erase(up.begin());
  const ax::Cols dn = plt::ramp({"#FBD5CC", "#B2182B"}, std::max(nb, 2));
  char t[96];
  std::snprintf(t, sizeof t, "Distance from each series' median, %d bands (above blue, below red)", nb);
  Pan p = frame(t, "Period", "");
  prdax(d, p);
  p.y1 = ng;
  p.pad = false;
  std::vector<double> at;
  for (int k = 0; k < ng; ++k) {
    const double b0 = ng - k - 1;
    at.push_back(b0 + 0.45);
    std::vector<double> x, v;
    ofser(d, k, dev, x, v);
    for (int b = 1; b <= nb; ++b) {
      for (int sg = 1; sg >= -1; sg -= 2) {
        std::vector<double> qx, qy;
        for (std::size_t i = 0; i < x.size(); ++i) {
          if (!std::isfinite(v[i])) continue;
          const double h = std::min(std::max(std::fabs(v[i]) - (b - 1) * bw, 0.0), bw) / bw * 0.9;
          qx.push_back(x[i]);
          qy.push_back(b0 + (sg * v[i] > 0 ? h : 0));
        }
        if (qx.empty()) continue;
        qx.push_back(qx.back());
        qy.push_back(b0);
        qx.push_back(qx.front());
        qy.push_back(b0);
        It q = poly(qx, qy, sg > 0 ? up[sat<std::size_t>(b - 1)] : dn[sat<std::size_t>(b - 1)], 1, "none");
        q.lw = 0;
        p.it.push_back(q);
      }
    }
  }
  labs(p.ay, d.lv, at);
  s.f.p.push_back(p);
  return true;
}

// .ppath(): estimates over event time with capped interval bars, rings
// before the adoption and discs after it.
inline Pan evpan(const Fig &f, const std::string &main, const std::vector<double> &t, const std::vector<double> &est,
                 const std::vector<double> &l, const std::vector<double> &h) {
  const ax::Cols cl = cols(f, 2);
  Pan p = frame(main.c_str(), "Periods since adoption", "Effect");
  rng(t, p.x0, p.x1);
  std::vector<double> all = est;
  all.insert(all.end(), l.begin(), l.end());
  all.insert(all.end(), h.begin(), h.end());
  all.push_back(0);
  rng(all, p.y0, p.y1);
  const double wx = p.x1 - p.x0, wy = p.y1 - p.y0;
  p.it.push_back(ref(p.x0 - wx, 0, p.x1 + wx, 0, 2));
  p.it.push_back(ref(-0.5, p.y0 - wy, -0.5, p.y1 + wy, 3));
  for (int post = 0; post < 2; ++post) {
    It sg, pt;
    sg.kind = SEG;
    sg.cap = 4;
    pt.kind = PTS;
    pt.pch = post ? 0 : 2;
    sg.col = pt.col = cl[sat<std::size_t>(post)];
    for (std::size_t i = 0; i < t.size(); ++i) {
      if ((t[i] > -0.5) != (post == 1)) continue;
      pt.x.push_back(t[i]);
      pt.y.push_back(est[i]);
      if (std::isfinite(l[i]) && std::isfinite(h[i]) && h[i] > l[i]) {
        sg.x.insert(sg.x.end(), {t[i], t[i]});
        sg.y.insert(sg.y.end(), {l[i], h[i]});
      }
    }
    p.it.push_back(sg);
    p.it.push_back(pt);
  }
  p.key = {{"before", cl[0], 0, 2}, {"after", cl[1], 0, 0}};
  p.corner = 2;
  return p;
}

inline bool bevg(In s) {
  const ax::Cols g = sv_(s.d, "panel");
  const std::vector<double> t = nv(s.d, "t"), e = nv(s.d, "est"), l = nv(s.d, "lo"), h = nv(s.d, "hi");
  for (const std::string &k : lvls(g)) {
    s.f.p.push_back(evpan(s.f, k, rows(t, g, k), rows(e, g, k), rows(l, g, k), rows(h, g, k)));
  }
  return true;
}

inline bool bcoh(In s) {
  const std::vector<double> co = nv(s.d, "cohort"), t = nv(s.d, "evtm"), v = nv(s.d, "value");
  std::vector<double> lv = co;
  std::sort(lv.begin(), lv.end());
  lv.erase(std::unique(lv.begin(), lv.end()), lv.end());
  const ax::Cols cl = cols(s.f, sat<int>(lv.size()));
  const Val &no = get(s.a, "norm");
  const bool nrm = !(no.t == Val::BOOL && !no.b) && !(no.t == Val::ARR && no.a.size() == 1 && !no.a[0].b);
  Pan p = frame("Outcome by adoption cohort", "Periods since adoption",
                nrm ? "Mean outcome (1 = period before adoption)" : "Mean outcome");
  rng(t, p.x0, p.x1);
  rng(v, p.y0, p.y1);
  const double wy = p.y1 - p.y0;
  p.it.push_back(poly({-1.5, -0.5, -0.5, -1.5}, {p.y0 - wy, p.y0 - wy, p.y1 + wy, p.y1 + wy}, "#9A9A9A", 0.25, "none"));
  for (std::size_t i = 0; i < lv.size(); ++i) {
    std::vector<double> x, y;
    for (std::size_t k = 0; k < co.size(); ++k) {
      if (co[k] == lv[i]) {
        x.push_back(t[k]);
        y.push_back(v[k]);
      }
    }
    p.it.push_back(line(x, y, cl[i], sat<int>(i) + 1, 2));
    p.key.push_back({"Adopted at " + num(lv[i], lv[i] == std::floor(lv[i]) ? 0 : 2), cl[i], sat<int>(i) + 1, -1});
  }
  p.corner = 2;
  s.f.p.push_back(p);
  return true;
}

// plot() of evtst() and stgdid(): the effects by time since adoption.
inline bool bevt(In s) {
  const Val &e = get(s.d, "est").t == Val::OBJ ? get(s.d, "est") : get(s.d, "dynamic");
  Pan p = evpan(s.f, get(s.d, "dynamic").t == Val::OBJ ? "Effect by time since adoption" : "Event study",
                nv(e, "evtm"), nv(e, "est"), nv(e, "lo"), nv(e, "hi"));
  p.ylab = "Effect, with its interval";
  s.f.p.push_back(p);
  return true;
}

// .ptwo(): the treated path and its comparison, both black in the grey
// palette and told apart by the line type.
inline ax::Cols two(const Fig &f) {
  return f.palette == "grey" ? ax::Cols{"#000000", "#000000"} : cols(f, 2);
}
// .padopt(): halfway between the last period before the adoption and the
// first after it (NA when there is none).
inline double adopt(const std::vector<double> &t, double last) {
  double nx = HUGE_VAL;
  for (double v : t) {
    if (v > last) nx = std::min(nx, v);
  }
  return std::isfinite(last) && std::isfinite(nx) ? (last + nx) / 2 : NA;
}
// Paths over the periods, one line type each, with a dotted line at the
// adoption `v` and, with `zero`, a dashed line at 0.
inline Pan paths(const char *main, const char *ylab, const std::vector<double> &t,
                 const std::vector<std::vector<double>> &ys, const ax::Cols &lab, const ax::Cols &cl, double v,
                 bool zero) {
  Pan p = frame(main, "Period", ylab);
  rng(t, p.x0, p.x1);
  std::vector<double> all;
  for (const std::vector<double> &y : ys) all.insert(all.end(), y.begin(), y.end());
  if (zero) all.push_back(0);
  rng(all, p.y0, p.y1);
  const double wx = p.x1 - p.x0, wy = p.y1 - p.y0;
  if (zero) p.it.push_back(ref(p.x0 - wx, 0, p.x1 + wx, 0, 2));
  if (std::isfinite(v)) p.it.push_back(ref(v, p.y0 - wy, v, p.y1 + wy, 3));
  for (std::size_t i = 0; i < ys.size(); ++i) {
    const int lty = sat<int>(i) + 1;
    p.it.push_back(line(t, ys[i], cl[i % cl.size()], lty, 2));
    if (ys.size() > 1) p.key.push_back({lab[i], cl[i % cl.size()], lty, -1});
  }
  return p;
}

// plot() of synthc(): the treated unit and its synthetic control, and the
// gap between them.
inline bool bsyn(In s) {
  const Val &pa = get(s.d, "path");
  const std::vector<double> t = nv(pa, "time");
  // the last period before the first adoption among the units
  const double f0 = lo(nv(s.a, "first"));
  double last = -HUGE_VAL;
  for (double e : t) {
    if (e < f0) last = std::max(last, e);
  }
  const double v = adopt(t, last);
  const ax::Cols cl = two(s.f);
  const Val &u = get(s.d, "unit");
  const std::string un = u.t == Val::STR ? u.s : u.t == Val::ARR && u.a.size() == 1 ? ax::vlab(u.a[0]) : "";
  s.f.p.push_back(paths(("Synthetic control: " + un).c_str(), "Outcome", t, {nv(pa, "treated"), nv(pa, "synth")},
                        {"Treated", "Synthetic"}, cl, v, false));
  s.f.p.push_back(paths("Treated minus synthetic", "Gap", t, {nv(pa, "gap")}, {"gap"}, cl, v, true));
  return true;
}

// plot() of sdidx(): the treated path and the weighted comparison.
inline bool bsdid(In s) {
  const Val &pa = get(s.d, "path");
  const std::vector<double> t = nv(pa, "time");
  s.f.p.push_back(paths("Synthetic difference in differences", "Outcome", t, {nv(pa, "treated"), nv(pa, "synth")},
                        {"Treated", "Weighted comparison"}, two(s.f), adopt(t, hi(nv(get(s.d, "timew"), "time"))),
                        false));
  return true;
}

#include "jsvg2.inc"
#include "jsvg3.inc"

// The plots: the function that makes the rows and the one that draws them.
struct Plot {
  const char *name;
  bool (*fn)(Ctx &, Val &);
  bool (*draw)(In);
};
const Plot plots[] = {
    {"plotlrnz", plotlrnz, blrnz}, {"plotgain", plotgain, bgain}, {"plotcdfs", plotcdfs, bcdfs},
    {"plotrdg", plotrdg, brdg},    {"plotvln", plotvln, bvln},    {"plothex", plothex, bhex},
    {"plotfnl", plotfnl, bfnl},    {"plotcal", plotcal, bcal},    {"plotmult", plotmult, bmult},
    {"plotbump", plotbump, bbump}, {"plothrz", plothrz, bhrz},    {"plotevgrid", plotevgrid, bevg},
    {"plotcoh", plotcoh, bcoh},  {"evtst", evtst, bevt},        {"stgdid", stgdid, bevt},
    {"synthc", synthc, bsyn},     {"sdidx", sdidx, bsdid},
    {"cntreg", cntreg, rcntreg},
    {"dastd", dastd, rdastd},
    {"grprr", grprr, rgrprr},
    {"sirex", sirex, rsirex},
    {"svymx", svymx, rsvymx},
    {"mndla", mndla, rmndla},
    {"thrtst", thrtst, rthrtst},
    {"sarml", sarml, rsarml},
    {"spcure", spcure, rspcure},
    {"mtrbnd", mtrbnd, rmtrbnd},
    {"pbbnd", pbbnd, rpbbnd},
    {"dkfig", dkfig, rdkfig},
    {"smdbal", smdbal, rsmdbal},
    {"jmono", jmono, rjmono},
    {"relchk", relchk, rrelchk},
    {"dissix", dissix, rdissix},
    {"gmmpnl", gmmpnl, rgmmpnl},
    {"ccepnl", ccepnl, rccepnl},
    {"fepnl", fepnl, rfepnl},
    {"hetfx", hetfx, rhetfx},
    {"typefx", typefx, rtypefx},
    {"plcday", plcday, rplcday},
    {"arwt", arwt, rarwt},
    {"crsev", crsev, rcrsev},
    {"gmixlr", gmixlr, rgmixlr},
    {"pedlr", pedlr, rpedlr},
    {"hypq", hypq, rhypq},
    {"lnsize", lnsize, rlnsize},
    {"lostay", lostay, rlostay},
    {"pinpl", pinpl, rpinpl},
    {"lmdix", lmdix, rlmdix},
    {"trndrr", trndrr, rtrndrr},
    {"stkflw", stkflw, rstkflw},
    {"dktflw", dktflw, rdktflw},
    {"jpipe", jpipe, rjpipe},
    {"hsper", hsper, rhsper},
    {"cumrsk", cumrsk, rcumrsk},
    {"agemix", agemix, ragemix},
    {"apcie", apcie, rapcie},
    {"kmrec", kmrec, rkmrec},
    {"lftab", lftab, rlftab},
    {"lnroc", lnroc, rlnroc},
    {"lnkpr", lnkpr, rlnkpr},
    {"cconc", cconc, rcconc},
    {"crossk", crossk, rcrossk},
    {"jdgiv", jdgiv, rjdgiv},
    {"rptvic", rptvic, rrptvic},
    {"ctabx", ctabx, rctabx},
    {"knxnr", knxnr, rknxnr},
    {"qdrtst", qdrtst, rqdrtst},
    {"trnmat", trnmat, rtrnmat},
    {"zcomp", zcomp, rzcomp},
    {"kdsurf", kdsurf, rkdsurf},
    {"tpsurf", tpsurf, rtpsurf},
    {"sdbsim", sdbsim, rsdbsim},
    {"prohs", prohs, rprohs},
    {"xcser", xcser, rxcser},
    {"hwkres", hwkres, rhwkres},
    {"mmshr", mmshr, rmmshr},
    {"bufas", bufas, rbufas},
    {"cmpsc", cmpsc, rcmpsc},
    {"ebshr", ebshr, rebshr},
    {"ibench", ibench, ribench},
    {"kscan", kscan, rkscan},
    {"didfe", didfe, rdidfe},
    {"wcbx", wcbx, rwcbx},
    {"yoypx", yoypx, ryoypx},
    {"indgap", indgap, rindgap},
    {"dualcont", dualcont, rdualcont},
    {"commidx", commidx, rcommidx},
    {"stdprev", stdprev, rstdprev},
    {"ordcmp", ordcmp, rordcmp},
    {"welprof", welprof, rwelprof},
    {"healchg", healchg, rhealchg},
    {"holbal", holbal, rholbal},
    {"trchlthgap", trchlthgap, rtrchlthgap},
    {"cnscwb", cnscwb, rcnscwb},
    {"cwbscr", cwbscr, rcwbscr},
    {"ppmlfe", ax::ppmlfe, rppmlfe},
    {"expgmm", ax::expgmm, rexpgmm},
    {"ivfe", ax::ivfe, rivfe},
    {"ppmlcf", ax::ppmlcf, rppmlcf},
    {"updwn", ax::updwn, rupdwn},
    {"runfx", runfx, rrunfx},
    {"dosebin", dosebin, rdosebin},
    {"jntx", jntx, rjntx},
    {"ptmap", ptmap, rptmap},
    {"dbscn", dbscn, rdbscn},
    {"tipplt", tipplt, rtipplt},
    {"eceplt", eceplt, receplt},
    {"rdloc", rdloc, rrdloc},
    {"rddens", rddens, rrddens}};

}  // namespace sv

// plotsvg()
bool plotsvg(Ctx &c, Val &out) {
  std::string nm, main, xlab, ylab, pal = "grey";
  if (!sv::str1(c["plot"], nm) || nm.empty()) return c.fail("`plot` must be the name of one function that has a plot.");
  const sv::Plot *pl = nullptr;
  std::string all;
  for (const sv::Plot &p : sv::plots) {
    if (nm == p.name) pl = &p;
    all += (all.empty() ? "" : ", ") + std::string(p.name);
  }
  if (pl == nullptr) return c.fail("`plot` must be one of %s.", all.c_str());
  const Val &a = c["args"];
  if (a.t != Val::OBJ) return c.fail("`args` must be the named arguments of %s().", pl->name);
  // a data frame among them has columns of one length, as jrun asks of
  // every argument
  for (std::size_t i = 0; i < a.a.size(); ++i) {
    const Val &v = a.a[i];
    if (v.t != Val::OBJ || v.a.empty() || v.a[0].t != Val::ARR) continue;
    for (const Val &x : v.a) {
      if (x.t == Val::ARR && x.a.size() != v.a[0].a.size()) {
        return c.fail("the columns of `%s` differ in length.", a.k[i].c_str());
      }
    }
  }
  if (!sv::str1(c["main"], main) || !sv::str1(c["xlab"], xlab) || !sv::str1(c["ylab"], ylab)) {
    return c.fail("`main`, `xlab` and `ylab` must each be one text.");
  }
  if (!sv::str1(c["palette"], pal) ||
      !(pal == "okabe" || pal == "diverging" || pal == "sequential" || pal == "grey")) {
    return c.fail("the palette must be one of okabe, diverging, sequential, grey.");
  }
  sv::Fig f;
  f.palette = pal;
  const Val &go = c["grid"], &wo = c["width"], &ho = c["height"];
  if (go.t != Val::NUL) {
    const Val *e = go.t == Val::ARR && go.a.size() == 1 ? &go.a[0] : &go;
    if (e->t != Val::BOOL) return c.fail("`grid` must be TRUE or FALSE.");
    f.grid = e->b;
  }
  if (wo.t != Val::NUL) {
    if (!plt::num1(wo, 200, 4000, false)) return c.fail("`width` must be one number from 200 to 4000.");
    f.w = ax::nums(wo)[0];
  }
  if (ho.t != Val::NUL) {
    if (!plt::num1(ho, 150, 4000, false)) return c.fail("`height` must be one number from 150 to 4000.");
    f.h = ax::nums(ho)[0];
  }
  // the plot function, on its own arguments and under its own name
  Ctx q = c;
  q.fn = pl->name;
  q.in = &a;
  Val d;
  if (!pl->fn(q, d)) {
    c.err = q.err;
    return false;
  }
  c.warn = q.warn;
  CK(pl->draw({c, a, d, f}));
  if (f.p.size() == 1) {
    if (!main.empty()) f.p[0].main = main;
    if (!xlab.empty()) f.p[0].xlab = xlab;
    if (!ylab.empty()) f.p[0].ylab = ylab;
  } else {
    f.main = main;
  }
  out = vobj();
  put(out, "svg", vstr(sv::draw(f)));
  put(out, "data", d);
  return true;
}
