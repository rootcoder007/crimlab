// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// More plots of R/plotx.R for clab_jrun(), as in jplt.inc: nothing is
// drawn, and each returns the rows the R function returns invisibly, with
// the checks of R in their order and words. The densities of plotrdg()
// and plotvln() are clab_dens() (src/dens.cpp), stats::density() with R's
// bandwidth rules, binning, FFT and interpolation: equal to R's to
// rounding, not bit for bit. The attribute `summary` of plotvln() is not
// part of the result.

namespace plt {

// .pbw(): the rule (0 for a number) and the bandwidth.
inline bool pbw(Ctx &c, int &rule, double &bw) {
  rule = 1;
  bw = 0.0;
  if (!c.has("bw")) return true;
  const Val &b = c["bw"];
  static const char *const rs[] = {"nrd0", "nrd", "ucv", "bcv", "sj", "sj-ste", "sj-dpi"};
  static const int rn[] = {1, 2, 3, 4, 5, 5, 6};
  if (b.t == Val::STR) {
    for (std::size_t i = 0; i < 7; ++i) {
      if (b.s == rs[i]) {
        rule = rn[i];
        return true;
      }
    }
  } else if (ax::numv(b)) {
    const std::vector<double> v = vals(b);
    if (v.size() == 1 && std::isfinite(v[0]) && v[0] > 0) {
      rule = 0;
      bw = v[0];
      return true;
    }
  }
  return c.fail("`bw` must be one positive number or one of nrd0, nrd, ucv, bcv, sj, sj-ste, sj-dpi.");
}

// stats::density(v, bw = bw): its 512 points and densities, or R's
// message (without the name of this function, as R gives it).
inline bool dens(Ctx &c, const std::vector<double> &v, int rule, double bw, std::vector<double> &x,
                 std::vector<double> &y) {
  x.assign(512, 0.0);
  y.assign(512, 0.0);
  double st[2] = {0.0, 0.0};
  int why = 0;
  const int rc = clab_dens(v.data(), nullptr, sat<std::int64_t>(v.size()), rule, bw, 1.0, 0, 512, 3.0, 0, 0.0,
                           0, 0.0, x.data(), y.data(), st, &why);
  if (rc == CLAB_EDOM) {
    static const char *const m[] = {"",
                                    "need at least 2 points to select a bandwidth automatically",
                                    "data are constant in bandwidth calculation",
                                    "sample is too sparse to find TD",
                                    "sample is too sparse to find alph2",
                                    "no solution in the specified range of bandwidths",
                                    "non-finite 'bw'",
                                    "'bw' is not positive.",
                                    "non-finite 'from'",
                                    "",
                                    "non-finite 'to'",
                                    "'xmin' not less than 'xmax'"};
    c.err = m[why >= 0 && why <= 11 ? why : 0];
    return false;
  }
  if (rc == CLAB_OK && why == 9) c.warn.push_back("minimum occurred at one end of the range");
  return kchk(c, rc);
}

// The values of x in group l.
inline std::vector<double> ingrp(const Xg &d, const std::string &l) {
  std::vector<double> v;
  for (std::size_t i = 0; i < d.x.size(); ++i) {
    if (d.g[i] == l) v.push_back(d.x[i]);
  }
  return v;
}

// The densities of the groups with two values or more, as rows.
inline bool grpdens(Ctx &c, const Xg &d, int rule, double bw, Frame &f, bool &any) {
  ax::Cols g;
  std::vector<double> xs, ys;
  any = false;
  for (const std::string &l : d.lv) {
    const std::vector<double> v = ingrp(d, l);
    if (v.size() < 2) continue;
    std::vector<double> x, y;
    CK(dens(c, v, rule, bw, x, y));
    any = true;
    for (std::size_t i = 0; i < x.size(); ++i) {
      g.push_back(l);
      xs.push_back(x[i]);
      ys.push_back(y[i]);
    }
  }
  f.str("grp", g);
  f.num("x", xs);
  f.num("y", ys);
  return true;
}

}  // namespace plt

// plotrdg()
bool plotrdg(Ctx &c, Val &out) {
  plt::Xg d;
  CK(plt::pxg(c, d));
  int rule = 1;
  double bw = 0.0;
  CK(plt::pbw(c, rule, bw));
  Frame f;
  bool any = false;
  CK(plt::grpdens(c, d, rule, bw, f, any));
  if (!any) return c.fail("no group has two values.");
  out = f.o;
  return true;
}

// plotvln()
bool plotvln(Ctx &c, Val &out) {
  plt::Xg d;
  CK(plt::pxg(c, d));
  int rule = 1;
  double bw = 0.0;
  CK(plt::pbw(c, rule, bw));
  Frame f;
  bool any = false;
  CK(plt::grpdens(c, d, rule, bw, f, any));
  if (!any) {
    // R sets the attribute summary on NULL
    c.err = "attempt to set an attribute on NULL";
    return false;
  }
  out = f.o;
  return true;
}

namespace plt {

// .prange() of finite values: their range, widened by 0.5 each side when
// it is one value.
inline void prange(const std::vector<double> &v, double &lo, double &hi) {
  lo = *std::min_element(v.begin(), v.end());
  hi = *std::max_element(v.begin(), v.end());
  if (hi - lo == 0) {
    lo -= 0.5;
    hi += 0.5;
  }
}
// R's x %% 2 of a whole number: 0 or 1.
inline double mod2(double r) { return r - 2 * std::floor(r / 2); }

}  // namespace plt

// plothex(): the hexagons in order of row, then column, as R orders them.
bool plothex(Ctx &c, Val &out) {
  const Val &xo = c["x"], &yo = c["y"];
  if (!ax::numv(xo)) return c.fail("`x` must be a numeric vector.");
  const std::vector<double> x0 = plt::vals(xo);
  std::size_t nf = 0;
  for (double v : x0) nf += std::isfinite(v) ? 1 : 0;
  if (nf < 2) return c.fail("`x` needs at least two finite values.");
  const std::vector<double> y0 = ax::numv(yo) ? plt::vals(yo) : std::vector<double>();
  if (!ax::numv(yo) || y0.size() != x0.size()) return c.fail("`y` must be numbers, one per value of `x`.");
  std::vector<double> x, y;
  for (std::size_t i = 0; i < x0.size(); ++i) {
    if (!std::isfinite(x0[i]) || !std::isfinite(y0[i])) continue;
    x.push_back(x0[i]);
    y.push_back(y0[i]);
  }
  if (x.size() < 2) return c.fail("needs at least two points with both values.");
  double bins = 30;
  if (c.has("bins")) {
    const std::vector<double> b = ax::numv(c["bins"]) ? plt::vals(c["bins"]) : std::vector<double>();
    const bool ok = b.size() == 1 && std::isfinite(b[0]) && b[0] >= 2 && b[0] <= 500 && b[0] == std::round(b[0]);
    if (!ok) return c.fail("`bins` must be one whole number from 2 to 500.");
    bins = b[0];
  }
  double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
  plt::prange(x, x1, x2);
  plt::prange(y, y1, y2);
  const double h = std::sqrt(3.0) / 2;
  std::map<std::pair<double, double>, int> cnt;
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double sx = (x[i] - x1) / (x2 - x1) * bins, sy = (y[i] - y1) / (y2 - y1) * bins;
    const double r1 = std::nearbyint(sy / h), c1 = std::nearbyint(sx - plt::mod2(r1) / 2);
    double r2 = std::floor(sy / h);
    if (r2 == r1) r2 += 1;
    const double c2 = std::nearbyint(sx - plt::mod2(r2) / 2);
    const double d1 = (sx - c1 - plt::mod2(r1) / 2) * (sx - c1 - plt::mod2(r1) / 2) + (sy - r1 * h) * (sy - r1 * h);
    const double d2 = (sx - c2 - plt::mod2(r2) / 2) * (sx - c2 - plt::mod2(r2) / 2) + (sy - r2 * h) * (sy - r2 * h);
    ++cnt[d1 <= d2 ? std::make_pair(r1, c1) : std::make_pair(r2, c2)];
  }
  std::vector<double> hx, hy;
  std::vector<int> n;
  for (const auto &e : cnt) {
    const double rr = e.first.first, cc = e.first.second;
    hx.push_back(x1 + (cc + plt::mod2(rr) / 2) / bins * (x2 - x1));
    hy.push_back(y1 + rr * h / bins * (y2 - y1));
    n.push_back(e.second);
  }
  Frame f;
  f.num("x", hx);
  f.num("y", hy);
  f.ints("count", n.data(), n.size());
  out = f.o;
  return true;
}

// plotevgrid(): x an object of panels, each a result of evtst() or
// stgdid() (its `est` or `dynamic` rows taken, as R takes them) or those
// rows themselves; the points by panel, in the order given.
bool plotevgrid(Ctx &c, Val &out) {
  const Val &x = c["x"];
  bool ok = x.t == Val::OBJ && !x.a.empty();
  std::vector<const Val *> fr;
  for (std::size_t i = 0; ok && i < x.a.size(); ++i) {
    const Val &p = x.a[i];
    const Val *d = p.t == Val::OBJ ? ax::fc(p, "dynamic") : nullptr;
    if (d == nullptr || d->t != Val::OBJ) d = p.t == Val::OBJ ? ax::fc(p, "est") : nullptr;
    if (d == nullptr || d->t != Val::OBJ) d = &p;
    for (const char *k : {"evtm", "est", "lo", "hi"}) {
      const Val *cl = d->t == Val::OBJ ? ax::fc(*d, k) : nullptr;
      ok = ok && cl != nullptr && isnum(*cl);
    }
    ok = ok && !x.k[i].empty();
    fr.push_back(d);
  }
  if (!ok) return c.fail("`x` must be a named list of evtst or stgdid results.");
  ax::Cols pn;
  std::vector<double> t, e, lo, hi;
  for (std::size_t i = 0; i < fr.size(); ++i) {
    const std::vector<double> a = ax::colv(*ax::fc(*fr[i], "evtm")), b = ax::colv(*ax::fc(*fr[i], "est")),
                              l = ax::colv(*ax::fc(*fr[i], "lo")), h = ax::colv(*ax::fc(*fr[i], "hi"));
    if (b.size() != a.size() || l.size() != a.size() || h.size() != a.size()) {
      return c.fail("`x` must be a named list of evtst or stgdid results.");
    }
    for (std::size_t j = 0; j < a.size(); ++j) {
      pn.push_back(x.k[i]);
      t.push_back(a[j]);
      e.push_back(b[j]);
      lo.push_back(l[j]);
      hi.push_back(h[j]);
    }
  }
  Frame f;
  f.str("panel", pn);
  f.num("t", t);
  f.num("est", e);
  f.num("lo", lo);
  f.num("hi", hi);
  out = f.o;
  return true;
}

namespace plt {

// The colour names of R's colors(), in one text (spaces between).
const char *const rcols =
    "white aliceblue antiquewhite antiquewhite1 antiquewhite2 antiquewhite3 antiquewhite4 aquamarine "
    "aquamarine1 aquamarine2 aquamarine3 aquamarine4 azure azure1 azure2 azure3 azure4 beige bisque "
    "bisque1 bisque2 bisque3 bisque4 black blanchedalmond blue blue1 blue2 blue3 blue4 blueviolet brown "
    "brown1 brown2 brown3 brown4 burlywood burlywood1 burlywood2 burlywood3 burlywood4 cadetblue "
    "cadetblue1 cadetblue2 cadetblue3 cadetblue4 chartreuse chartreuse1 chartreuse2 chartreuse3 "
    "chartreuse4 chocolate chocolate1 chocolate2 chocolate3 chocolate4 coral coral1 coral2 coral3 coral4 "
    "cornflowerblue cornsilk cornsilk1 cornsilk2 cornsilk3 cornsilk4 cyan cyan1 cyan2 cyan3 cyan4 "
    "darkblue darkcyan darkgoldenrod darkgoldenrod1 darkgoldenrod2 darkgoldenrod3 darkgoldenrod4 "
    "darkgray darkgreen darkgrey darkkhaki darkmagenta darkolivegreen darkolivegreen1 darkolivegreen2 "
    "darkolivegreen3 darkolivegreen4 darkorange darkorange1 darkorange2 darkorange3 darkorange4 "
    "darkorchid darkorchid1 darkorchid2 darkorchid3 darkorchid4 darkred darksalmon darkseagreen "
    "darkseagreen1 darkseagreen2 darkseagreen3 darkseagreen4 darkslateblue darkslategray darkslategray1 "
    "darkslategray2 darkslategray3 darkslategray4 darkslategrey darkturquoise darkviolet deeppink "
    "deeppink1 deeppink2 deeppink3 deeppink4 deepskyblue deepskyblue1 deepskyblue2 deepskyblue3 "
    "deepskyblue4 dimgray dimgrey dodgerblue dodgerblue1 dodgerblue2 dodgerblue3 dodgerblue4 firebrick "
    "firebrick1 firebrick2 firebrick3 firebrick4 floralwhite forestgreen gainsboro ghostwhite gold gold1 "
    "gold2 gold3 gold4 goldenrod goldenrod1 goldenrod2 goldenrod3 goldenrod4 gray gray0 gray1 gray2 "
    "gray3 gray4 gray5 gray6 gray7 gray8 gray9 gray10 gray11 gray12 gray13 gray14 gray15 gray16 gray17 "
    "gray18 gray19 gray20 gray21 gray22 gray23 gray24 gray25 gray26 gray27 gray28 gray29 gray30 gray31 "
    "gray32 gray33 gray34 gray35 gray36 gray37 gray38 gray39 gray40 gray41 gray42 gray43 gray44 gray45 "
    "gray46 gray47 gray48 gray49 gray50 gray51 gray52 gray53 gray54 gray55 gray56 gray57 gray58 gray59 "
    "gray60 gray61 gray62 gray63 gray64 gray65 gray66 gray67 gray68 gray69 gray70 gray71 gray72 gray73 "
    "gray74 gray75 gray76 gray77 gray78 gray79 gray80 gray81 gray82 gray83 gray84 gray85 gray86 gray87 "
    "gray88 gray89 gray90 gray91 gray92 gray93 gray94 gray95 gray96 gray97 gray98 gray99 gray100 green "
    "green1 green2 green3 green4 greenyellow grey grey0 grey1 grey2 grey3 grey4 grey5 grey6 grey7 grey8 "
    "grey9 grey10 grey11 grey12 grey13 grey14 grey15 grey16 grey17 grey18 grey19 grey20 grey21 grey22 "
    "grey23 grey24 grey25 grey26 grey27 grey28 grey29 grey30 grey31 grey32 grey33 grey34 grey35 grey36 "
    "grey37 grey38 grey39 grey40 grey41 grey42 grey43 grey44 grey45 grey46 grey47 grey48 grey49 grey50 "
    "grey51 grey52 grey53 grey54 grey55 grey56 grey57 grey58 grey59 grey60 grey61 grey62 grey63 grey64 "
    "grey65 grey66 grey67 grey68 grey69 grey70 grey71 grey72 grey73 grey74 grey75 grey76 grey77 grey78 "
    "grey79 grey80 grey81 grey82 grey83 grey84 grey85 grey86 grey87 grey88 grey89 grey90 grey91 grey92 "
    "grey93 grey94 grey95 grey96 grey97 grey98 grey99 grey100 honeydew honeydew1 honeydew2 honeydew3 "
    "honeydew4 hotpink hotpink1 hotpink2 hotpink3 hotpink4 indianred indianred1 indianred2 indianred3 "
    "indianred4 ivory ivory1 ivory2 ivory3 ivory4 khaki khaki1 khaki2 khaki3 khaki4 lavender "
    "lavenderblush lavenderblush1 lavenderblush2 lavenderblush3 lavenderblush4 lawngreen lemonchiffon "
    "lemonchiffon1 lemonchiffon2 lemonchiffon3 lemonchiffon4 lightblue lightblue1 lightblue2 lightblue3 "
    "lightblue4 lightcoral lightcyan lightcyan1 lightcyan2 lightcyan3 lightcyan4 lightgoldenrod "
    "lightgoldenrod1 lightgoldenrod2 lightgoldenrod3 lightgoldenrod4 lightgoldenrodyellow lightgray "
    "lightgreen lightgrey lightpink lightpink1 lightpink2 lightpink3 lightpink4 lightsalmon lightsalmon1 "
    "lightsalmon2 lightsalmon3 lightsalmon4 lightseagreen lightskyblue lightskyblue1 lightskyblue2 "
    "lightskyblue3 lightskyblue4 lightslateblue lightslategray lightslategrey lightsteelblue "
    "lightsteelblue1 lightsteelblue2 lightsteelblue3 lightsteelblue4 lightyellow lightyellow1 "
    "lightyellow2 lightyellow3 lightyellow4 limegreen linen magenta magenta1 magenta2 magenta3 magenta4 "
    "maroon maroon1 maroon2 maroon3 maroon4 mediumaquamarine mediumblue mediumorchid mediumorchid1 "
    "mediumorchid2 mediumorchid3 mediumorchid4 mediumpurple mediumpurple1 mediumpurple2 mediumpurple3 "
    "mediumpurple4 mediumseagreen mediumslateblue mediumspringgreen mediumturquoise mediumvioletred "
    "midnightblue mintcream mistyrose mistyrose1 mistyrose2 mistyrose3 mistyrose4 moccasin navajowhite "
    "navajowhite1 navajowhite2 navajowhite3 navajowhite4 navy navyblue oldlace olivedrab olivedrab1 "
    "olivedrab2 olivedrab3 olivedrab4 orange orange1 orange2 orange3 orange4 orangered orangered1 "
    "orangered2 orangered3 orangered4 orchid orchid1 orchid2 orchid3 orchid4 palegoldenrod palegreen "
    "palegreen1 palegreen2 palegreen3 palegreen4 paleturquoise paleturquoise1 paleturquoise2 "
    "paleturquoise3 paleturquoise4 palevioletred palevioletred1 palevioletred2 palevioletred3 "
    "palevioletred4 papayawhip peachpuff peachpuff1 peachpuff2 peachpuff3 peachpuff4 peru pink pink1 "
    "pink2 pink3 pink4 plum plum1 plum2 plum3 plum4 powderblue purple purple1 purple2 purple3 purple4 "
    "red red1 red2 red3 red4 rosybrown rosybrown1 rosybrown2 rosybrown3 rosybrown4 royalblue royalblue1 "
    "royalblue2 royalblue3 royalblue4 saddlebrown salmon salmon1 salmon2 salmon3 salmon4 sandybrown "
    "seagreen seagreen1 seagreen2 seagreen3 seagreen4 seashell seashell1 seashell2 seashell3 seashell4 "
    "sienna sienna1 sienna2 sienna3 sienna4 skyblue skyblue1 skyblue2 skyblue3 skyblue4 slateblue "
    "slateblue1 slateblue2 slateblue3 slateblue4 slategray slategray1 slategray2 slategray3 slategray4 "
    "slategrey snow snow1 snow2 snow3 snow4 springgreen springgreen1 springgreen2 springgreen3 "
    "springgreen4 steelblue steelblue1 steelblue2 steelblue3 steelblue4 tan tan1 tan2 tan3 tan4 thistle "
    "thistle1 thistle2 thistle3 thistle4 tomato tomato1 tomato2 tomato3 tomato4 turquoise turquoise1 "
    "turquoise2 turquoise3 turquoise4 violet violetred violetred1 violetred2 violetred3 violetred4 wheat "
    "wheat1 wheat2 wheat3 wheat4 whitesmoke yellow yellow1 yellow2 yellow3 yellow4 yellowgreen ";

// R's StrMatch(): equal ignoring case and blanks.
inline bool strmatch(const char *s, const char *t, std::size_t ns) {
  const char *const e = s + ns;
  for (;;) {
    if (s == e && *t == '\0') return true;
    if (s != e && *s == ' ') {
      ++s;
      continue;
    }
    if (*t == ' ') {
      ++t;
      continue;
    }
    if (s == e || *t == '\0') return false;
    if (std::tolower(static_cast<unsigned char>(*s++)) != std::tolower(static_cast<unsigned char>(*t++))) return false;
  }
}
// .iscol(): what grDevices::col2rgb() takes, as R's str2col() reads it: #
// and 3, 4, 6 or 8 hex digits; text that starts with a digit and is one
// number; "NA", "transparent" or a name of colors(); a number from 0 to
// 2147483647.
inline bool iscol(const Val &o) {
  const Val *v = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  if (v->t == Val::NUM) return std::isfinite(v->n) && v->n >= 0 && v->n <= 2147483647.0;
  if (v->t != Val::STR) return false;
  const std::string &s = v->s;
  if (!s.empty() && s[0] == '#') {
    bool hex = s.size() == 4 || s.size() == 5 || s.size() == 7 || s.size() == 9;
    for (std::size_t i = 1; hex && i < s.size(); ++i) hex = std::isxdigit(static_cast<unsigned char>(s[i])) != 0;
    return hex;
  }
  if (!s.empty() && std::isdigit(static_cast<unsigned char>(s[0]))) {
    char *end = nullptr;
    std::strtod(s.c_str(), &end);
    return end != nullptr && *end == '\0';
  }
  if (s == "NA" || s == "transparent") return true;
  for (const char *p = rcols; *p != '\0';) {
    const char *q = std::strchr(p, ' ');
    const std::size_t n = q == nullptr ? std::strlen(p) : sat<std::size_t>(q - p);
    if (strmatch(p, s.c_str(), n)) return true;
    p += n + (q == nullptr ? 0 : 1);
  }
  return false;
}
// .pnum1(v, lo, hi, int)
inline bool num1(const Val &o, double lo, double hi, bool whole) {
  if (!ax::numv(o)) return false;
  const std::vector<double> v = vals(o);
  return v.size() == 1 && std::isfinite(v[0]) && v[0] >= lo && v[0] <= hi && (!whole || v[0] == std::round(v[0]));
}

}  // namespace plt

// plotopt(): outside R there is no session to keep options in, so the
// options given are checked as R checks them and returned with the
// defaults for the others: the set a caller passes on. `panel` = NA
// cannot be given (null is an argument left out).
bool plotopt(Ctx &c, Val &out) {
  static const char *const nm[] = {"palette", "grid", "cex", "las", "bg", "panel"};
  static const char *const msg[] = {"", "TRUE or FALSE", "one positive number", "0, 1, 2 or 3", "one colour",
                                    "one colour, or NA"};
  out = vobj();
  put(out, "palette", vstr("okabe"));
  put(out, "grid", vbool(true));
  put(out, "cex", vnum(1));
  put(out, "las", vnum(1));
  put(out, "bg", vstr("white"));
  put(out, "panel", vstr("#EEF1F6"));
  for (std::size_t i = 0; i < c.in->k.size(); ++i) {
    const std::string &k = c.in->k[i];
    const Val &v = c.in->a[i];
    if (v.t == Val::NUL) continue;
    std::size_t j = 0;
    while (j < 6 && k != nm[j]) ++j;
    if (j == 6) return c.fail("the options are palette, grid, cex, las, bg, panel.");
    const Val *s = v.t == Val::ARR && v.a.size() == 1 ? &v.a[0] : &v;
    bool ok = true;
    if (j == 0) {
      ok = s->t == Val::STR && (s->s == "okabe" || s->s == "diverging" || s->s == "sequential" || s->s == "grey");
      if (!ok) return c.fail("the palette must be one of okabe, diverging, sequential, grey.");
    }
    if (j == 1) ok = s->t == Val::BOOL;
    if (j == 2) ok = plt::num1(v, std::numeric_limits<double>::min(), HUGE_VAL, false);
    if (j == 3) ok = plt::num1(v, 0, 3, true);
    if (j >= 4) ok = plt::iscol(v);
    if (!ok) return c.fail("`%s` must be %s.", nm[j], msg[j]);
    for (std::size_t m = 0; m < out.k.size(); ++m) {
      if (out.k[m] == k) out.a[m] = *s;
    }
  }
  return true;
}
