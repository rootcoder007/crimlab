// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// clab_jrun(): every function of the Python package by name, its named
// arguments a JSON object, its result JSON. The checks, their order and
// their words, and the results follow the Python binding
// (py/ext/clab.cpp) and through it the R binding: a data frame is an
// object of columns, a named vector an object of numbers. A JSON object
// given for a per-group argument names its groups, as a named vector does
// in R; other named arguments are matched to the first one by name.
//
// Adding a function: its entry point below and one line in the table at
// the end. The command-line program and every binding pick it up from
// there.
#include "jrun.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "crimlab.h"
#include "json.h"

extern "C" int clab_selftest(char *sink, int cap, int *nfail);

namespace {

using namespace clabj;


// Building results.

Val vnum(double x) {
  Val v;
  v.t = Val::NUM;
  v.n = x;
  return v;
}
Val vbool(bool b) {
  Val v;
  v.t = Val::BOOL;
  v.b = b;
  return v;
}
Val vstr(const std::string &s) {
  Val v;
  v.t = Val::STR;
  v.s = s;
  return v;
}
Val varr() {
  Val v;
  v.t = Val::ARR;
  return v;
}
Val vnums(const double *x, std::size_t n) {
  Val v = varr();
  for (std::size_t i = 0; i < n; ++i) v.a.push_back(vnum(x[i]));
  return v;
}
Val vnums(const std::vector<double> &x) { return vnums(x.data(), x.size()); }
Val vints(const std::vector<int> &x) {
  Val v = varr();
  for (int i : x) v.a.push_back(vnum(i));
  return v;
}
Val vstrs(const std::vector<std::string> &x) {
  Val v = varr();
  for (const std::string &s : x) v.a.push_back(vstr(s));
  return v;
}
void put(Val &o, const char *k, Val v) {
  o.t = Val::OBJ;
  o.k.push_back(k);
  o.a.push_back(std::move(v));
}
// A named vector.
Val nvec(int n, const char *const *nm, const double *x) {
  Val o;
  for (int i = 0; i < n; ++i) put(o, nm[i], vnum(x[i]));
  return o;
}
// True at the reference group i (from 1) of n.
Val refcol(int n, int i) {
  Val v = varr();
  for (int k = 0; k < n; ++k) v.a.push_back(vbool(k == i - 1));
  return v;
}

// One call: its function, its arguments and its message.

struct Ctx {
  const char *fn;
  const Val *in;
  std::string err;
  const Val &operator[](const char *k) const {
    for (std::size_t i = 0; i < in->k.size(); ++i) {
      if (in->k[i] == k) return in->a[i];
    }
    return nul;
  }
  // null is the same as left out
  bool has(const char *k) const { return (*this)[k].t != Val::NUL; }
  bool fail(const char *f, ...) {
    std::va_list ap, aq;
    va_start(ap, f);
    va_copy(aq, ap);
    std::vector<char> b(256);
    const int n = std::vsnprintf(b.data(), b.size(), f, ap);
    if (n >= static_cast<int>(b.size())) {
      b.resize(static_cast<std::size_t>(n) + 1);
      std::vsnprintf(b.data(), b.size(), f, aq);
    }
    va_end(aq);
    va_end(ap);
    err = std::string(fn) + ": " + b.data();
    return false;
  }
};

#define CK(x)               \
  do {                      \
    if (!(x)) return false; \
  } while (0)

bool kchk(Ctx &c, int code) { return code == CLAB_OK || c.fail("%s", clab_strerror(code)); }

// The arithmetic helpers of the R binding (src/rchk.cpp).
double rd(double x) {
  volatile double v = x;
  return v;
}

double nsum(const double *x, std::size_t n) {
  double s = 0.0, c = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double t = s + x[i];
    c += std::fabs(s) >= std::fabs(x[i]) ? (s - t) + x[i] : (x[i] - t) + s;
    s = t;
  }
  return std::isfinite(s) ? s + c : s;
}

double nmean(const double *x, std::size_t n) {
  const double dn = static_cast<double>(n);
  double s = 0.0;
  for (std::size_t i = 0; i < n; ++i) s += x[i];
  s /= dn;
  if (!std::isfinite(s)) return s;
  double t = 0.0;
  for (std::size_t i = 0; i < n; ++i) t += x[i] - s;
  return s + t / dn;
}

double rpow(double x, double y) { return y == 2.0 ? x * x : std::pow(x, y); }

// as.character() of one double, as R 4.3 and later writes it.
std::string rfmt(double x) {
  if (std::isnan(x)) return "NaN";
  if (std::isinf(x)) return x > 0 ? "Inf" : "-Inf";
  if (x == 0) return "0";
  char b[64];
  const bool whole = x == std::floor(x) && std::fabs(x) <= 9007199254740992.0;
  const int top = whole ? 17 : 15;
  std::snprintf(b, sizeof b, "%.*e", top - 1, x);
  const double v = std::strtod(b, nullptr);
  int d = 1;
  for (; d < top; ++d) {
    std::snprintf(b, sizeof b, "%.*e", d - 1, x);
    if (std::strtod(b, nullptr) == v) break;
  }
  std::snprintf(b, sizeof b, "%.*e", d - 1, x);
  const int e = std::atoi(std::strchr(b, 'e') + 1);
  const int neg = x < 0;
  const int sci = neg + d + (d > 1) + (std::abs(e) >= 100 ? 5 : 4);
  const int rgt = std::max(0, d - e - 1);
  const int fix = neg + (e >= 0 ? e + 1 : 1) + (rgt > 0 ? rgt + 1 : 0);
  if (fix <= sci) std::snprintf(b, sizeof b, "%.*f", rgt, x);
  return b;
}

// Reading arguments.

enum Got { NUM, NOTNUM, MATRIX };

// Numbers and, when the argument was an object, their names.
struct Dv {
  std::vector<double> v;
  std::vector<std::string> nm;
  std::size_t n() const { return v.size(); }
  double operator[](std::size_t i) const { return v[i]; }
  const double *p() const { return v.data(); }
};

// A number, or an array or object of numbers and nulls (NaN).
Got readd(const Val &o, Dv &d) {
  d = Dv();
  if (o.t == Val::NUM) {
    d.v.push_back(o.n);
    return NUM;
  }
  if (o.t != Val::ARR && o.t != Val::OBJ) return NOTNUM;
  for (const Val &e : o.a) {
    if (e.t == Val::NUM) d.v.push_back(e.n);
    else if (e.t == Val::NUL) d.v.push_back(NaN);
    else return e.t == Val::ARR && o.t == Val::ARR ? MATRIX : NOTNUM;
  }
  if (o.t == Val::OBJ) d.nm = o.k;
  return NUM;
}

// true, false and null (-1).
Got readb(const Val &o, std::vector<int> &b) {
  b.clear();
  if (o.t == Val::BOOL) {
    b.push_back(o.b);
    return NUM;
  }
  if (o.t != Val::ARR && o.t != Val::OBJ) return NOTNUM;
  for (const Val &e : o.a) {
    if (e.t == Val::BOOL) b.push_back(e.b);
    else if (e.t == Val::NUL) b.push_back(-1);
    else return e.t == Val::ARR && o.t == Val::ARR ? MATRIX : NOTNUM;
  }
  return NUM;
}

bool notvec(Ctx &c, const char *nm) {
  return c.fail("`%s` must be a vector, not a matrix with more than one row and column.", nm);
}

// .num(): finite numbers, len of them (any number when len < 0), all at
// least lo (above lo when open).
bool num(Ctx &c, const Val &o, const char *nm, Dv &d, long len = -1, double lo = -HUGE_VAL,
         bool open = false) {
  const Got g = readd(o, d);
  if (g == NOTNUM || d.n() < 1) return c.fail("`%s` must be finite numbers with no missing values.", nm);
  if (g == MATRIX) return notvec(c, nm);
  if (d.n() > INT_MAX) return c.fail("`%s` is too long for crimlab: more than %d values.", nm, INT_MAX);
  for (double x : d.v) {
    if (!std::isfinite(x)) return c.fail("`%s` must be finite numbers with no missing values.", nm);
  }
  if (len >= 0 && d.n() != static_cast<std::size_t>(len)) {
    return c.fail("`%s` must have length %ld, not %ld.", nm, len, static_cast<long>(d.n()));
  }
  for (double x : d.v) {
    if (open ? x <= lo : x < lo) {
      return c.fail("`%s` must be %s%s.", nm, open ? "greater than " : "at least ", rfmt(lo).c_str());
    }
  }
  return true;
}

// A per-group argument.
bool gnum(Ctx &c, const Val &o, const char *nm, Dv &d, long len = -1, double lo = -HUGE_VAL,
          bool open = false) {
  if (readd(o, d) == NOTNUM) {
    return c.fail("`%s` must be finite numbers, one per group: an array of numbers or an object "
                  "from group to number.", nm);
  }
  return num(c, o, nm, d, len, lo, open);
}

// x[1] when o is one number, NaN otherwise; type: o was not numbers at all.
double val1(const Val &o) {
  Dv d;
  return readd(o, d) == NUM && d.n() == 1 ? d[0] : NaN;
}

// .one(): one number strictly between lo and hi.
bool one(Ctx &c, const Val &o, const char *nm, double lo, double hi, double &v) {
  v = val1(o);
  if (std::isnan(v) || v <= lo || v >= hi) {
    return c.fail("`%s` must be one number strictly between %s and %s.", nm, rfmt(lo).c_str(),
                  rfmt(hi).c_str());
  }
  return true;
}

// .p01(): one number from 0 to 1.
bool p01(Ctx &c, const Val &o, const char *nm, double &v) {
  v = val1(o);
  if (std::isnan(v) || v < 0 || v > 1) return c.fail("`%s` must be one number from 0 to 1.", nm);
  return true;
}

// .whole(): one whole number from lo up to INT_MAX.
bool whole(Ctx &c, const Val &o, const char *nm, int lo, int &r) {
  const double v = val1(o);
  if (std::isnan(v) || v != std::floor(v) || v < lo || v > INT_MAX) {
    return c.fail("`%s` must be one whole number from %d.", nm, lo);
  }
  r = static_cast<int>(v);
  return true;
}

// .lgl(): true or false, n of them (any number when n < 0), none missing.
bool lgl(Ctx &c, const Val &o, const char *nm, long n, std::vector<int> &b) {
  const Got g = readb(o, b);
  bool na = false;
  for (int v : b) na = na || v < 0;
  if (g == NOTNUM || b.empty() || (n >= 0 && b.size() != static_cast<std::size_t>(n)) || na) {
    return c.fail("`%s` must be TRUE or FALSE, one per unit, with no missing values.", nm);
  }
  if (g == MATRIX) return notvec(c, nm);
  return true;
}

// The group names: those of the first argument given as an object. Every
// other argument given as an object is put in their order and must name
// the same groups. Null entries are arguments left out.
// ponytail: quadratic matching, fine for the group counts these methods see
bool align(Ctx &c, std::vector<std::pair<const char *, Dv *>> as, std::vector<std::string> &g) {
  g.clear();
  const char *gn = nullptr;
  for (auto &a : as) {
    if (a.second == nullptr || a.second->nm.empty()) continue;
    Dv &d = *a.second;
    if (gn == nullptr) {
      gn = a.first;
      g = d.nm;
      std::vector<std::string> s = g;
      std::sort(s.begin(), s.end());
      const auto t = std::adjacent_find(s.begin(), s.end());
      if (t != s.end()) return c.fail("`%s` names a group twice: %s.", gn, t->c_str());
      continue;
    }
    std::vector<double> v(g.size());
    bool ok = d.nm.size() == g.size();
    for (std::size_t i = 0; ok && i < g.size(); ++i) {
      const auto j = std::find(d.nm.begin(), d.nm.end(), g[i]);
      ok = j != d.nm.end();
      if (ok) v[i] = d.v[static_cast<std::size_t>(j - d.nm.begin())];
    }
    if (!ok) return c.fail("`%s` must name the same groups as `%s`.", a.first, gn);
    d.v = v;
    d.nm = g;
  }
  return true;
}

// .grp(): the group names, or "1", "2", ...
std::vector<std::string> grp(std::vector<std::string> g, std::size_t n) {
  if (g.empty()) {
    for (std::size_t i = 0; i < n; ++i) g.push_back(std::to_string(i + 1));
  }
  return g;
}

// .ref(): the group named or numbered (from 1) by r; 1 when left out.
bool ref(Ctx &c, const Val &r, const std::vector<std::string> &g, int &i) {
  const int n = static_cast<int>(g.size());
  if (r.t == Val::NUL) {
    i = 1;
    return true;
  }
  if (r.t == Val::STR) {
    for (int k = 0; k < n; ++k) {
      if (g[static_cast<std::size_t>(k)] == r.s) {
        i = k + 1;
        return true;
      }
    }
    return c.fail("`ref` names no group: %s.", r.s.c_str());
  }
  const double v = val1(r);
  if (std::isnan(v) || v != std::floor(v) || v < 1 || v > n) {
    return c.fail("`ref` must be one group name or a position from 1 to %d.", n);
  }
  i = static_cast<int>(v);
  return true;
}

// Some p[i] > 0 (where keep[i], when keep is given).
bool pos(const Dv &p, const std::vector<int> *keep = nullptr) {
  for (std::size_t i = 0; i < p.n(); ++i) {
    if ((keep == nullptr || (*keep)[i]) && p[i] > 0) return true;
  }
  return false;
}

void ones(Dv &d, std::size_t n) { d.v.assign(n, 1.0); }

// The weights w, or rep(1, n) when left out.
bool wts(Ctx &c, Dv &w, std::size_t n) {
  if (!c.has("w")) {
    ones(w, n);
    return true;
  }
  return num(c, c["w"], "w", w, static_cast<long>(n), 0.0);
}

// as.character() of one label: text, a number or a bool; false when
// missing or anything else.
bool lab(const Val &v, std::string &s) {
  switch (v.t) {
    case Val::STR: s = v.s; return true;
    case Val::NUM: s = rfmt(v.n); return true;
    case Val::BOOL: s = v.b ? "TRUE" : "FALSE"; return true;
    default: return false;
  }
}

// Labels: one, or an array of them.
bool labs(const Val &v, std::vector<std::string> &out) {
  out.clear();
  std::string s;
  if (lab(v, s)) {
    out.push_back(s);
    return true;
  }
  if (v.t != Val::ARR) return false;
  for (const Val &e : v.a) {
    if (!lab(e, s)) return false;
    out.push_back(s);
  }
  return true;
}

// Disparity.

bool vodtst(Ctx &c, Val &out) {
  Dv d, ni, sh;
  CK(gnum(c, c["day"], "day", d, -1, 0.0));
  const int n = static_cast<int>(d.n());
  CK(gnum(c, c["night"], "night", ni, n, 0.0, true));
  const bool has = c.has("shift");
  if (has) CK(gnum(c, c["shift"], "shift", sh, n, 0.0, true));
  std::vector<std::string> g;
  CK(align(c, {{"day", &d}, {"night", &ni}, {"shift", has ? &sh : nullptr}}, g));
  g = grp(g, d.n());
  int i = 1;
  CK(ref(c, c["ref"], g, i));
  if (d[static_cast<std::size_t>(i - 1)] <= 0) return c.fail("the reference group has no daylight stops.");
  const std::size_t un = d.n();
  std::vector<double> sd(un), sn(un), ra(un), adj(has ? un : 0);
  std::vector<int> sg(un);
  CK(kchk(c, clab_vod(d.p(), ni.p(), n, i - 1, has ? sh.p() : nullptr, sd.data(), sn.data(),
                      ra.data(), sg.data(), has ? adj.data() : nullptr)));
  put(out, "group", vstrs(g));
  put(out, "ref", refcol(n, i));
  put(out, "shday", vnums(sd));
  put(out, "shnight", vnums(sn));
  put(out, "ratio", vnums(ra));
  put(out, "sign", vints(sg));
  if (has) put(out, "adj", vnums(adj));
  return true;
}

bool orrisk(Ctx &c, Val &out) {
  Dv o;
  CK(num(c, c["or"], "or", o, 1, 0.0, true));
  const double v = o[0];
  if (!std::isfinite(1 / v)) return c.fail("`or` is too small: its reciprocal is not finite.");
  double q = 0, s = 0;
  CK(one(c, c["rate"], "rate", 0, 1, q));
  CK(one(c, c["share"], "share", 0, 1, s));
  double r[4];
  CK(kchk(c, clab_orrisk(v, q, s, r, r + 1)));
  r[2] = r[0] / r[1];
  r[3] = v;
  static const char *const nm[] = {"exposed", "unexposed", "rr", "or"};
  out = nvec(4, nm, r);
  return true;
}

bool outtst(Ctx &c, Val &out) {
  Dv q, t, w;
  CK(num(c, c["q"], "q", q, -1, 0.0));
  for (double x : q.v) {
    if (x > 1) return c.fail("`q` must be at most 1.");
  }
  CK(num(c, c["t"], "t", t, 1));
  CK(wts(c, w, q.n()));
  if (!pos(w)) return c.fail("`w` must have a positive total.");
  double r[2];
  CK(kchk(c, clab_outtst(q.p(), w.p(), static_cast<int>(q.n()), t[0], r, r + 1)));
  static const char *const nm[] = {"stop", "hit"};
  out = nvec(2, nm, r);
  return true;
}

bool dspbm(Ctx &c, Val &out) {
  Dv p, ct, f;
  CK(gnum(c, c["pop"], "pop", p, -1, 0.0, true));
  const int n = static_cast<int>(p.n());
  CK(gnum(c, c["contact"], "contact", ct, n, 0.0, true));
  CK(gnum(c, c["force"], "force", f, n, 0.0, true));
  std::vector<std::string> g;
  CK(align(c, {{"pop", &p}, {"contact", &ct}, {"force", &f}}, g));
  g = grp(g, p.n());
  int i = 1;
  CK(ref(c, c["ref"], g, i));
  const std::size_t un = p.n();
  std::vector<double> cd(un), fd(un), rdv(un);
  CK(kchk(c, clab_dspbm(p.p(), ct.p(), f.p(), n, i - 1, cd.data(), fd.data(), rdv.data())));
  put(out, "group", vstrs(g));
  put(out, "ref", refcol(n, i));
  put(out, "cd", vnums(cd));
  put(out, "fd", vnums(fd));
  put(out, "rd", vnums(rdv));
  return true;
}

// count: rows of offending groups, columns of victim groups, both in the
// order of pop.
bool dyadrt(Ctx &c, Val &out) {
  Dv p, x;
  CK(gnum(c, c["pop"], "pop", p, -1, 0.0, true));
  const std::size_t g = p.n();
  std::vector<std::string> nm;
  CK(align(c, {{"pop", &p}}, nm));
  nm = grp(nm, g);
  const Val &t = c["count"];
  bool ok = t.t == Val::ARR && t.a.size() == g;
  for (std::size_t i = 0; ok && i < g; ++i) ok = t.a[i].t == Val::ARR && t.a[i].a.size() == g;
  if (!ok) return c.fail("`count` must be a %d by %d matrix.", static_cast<int>(g), static_cast<int>(g));
  Val flat = varr();  // column by column
  for (std::size_t j = 0; j < g; ++j) {
    for (std::size_t i = 0; i < g; ++i) flat.a.push_back(t.a[i].a[j]);
  }
  CK(num(c, flat, "count", x, -1, 0.0));
  if (!pos(x)) return c.fail("`count` holds no offences.");
  const std::size_t m = x.n();
  std::vector<int> off(m), vic(m);
  std::vector<std::string> so(m), sv(m);
  for (std::size_t k = 0; k < m; ++k) {
    off[k] = static_cast<int>(k % g);
    vic[k] = static_cast<int>(k / g);
    so[k] = nm[k % g];
    sv[k] = nm[k / g];
  }
  std::vector<double> ro(m), rp(m), no(m), ra(m);
  CK(kchk(c, clab_dyadrt(x.p(), off.data(), vic.data(), static_cast<int>(m), p.p(),
                         static_cast<int>(g), ro.data(), rp.data(), no.data(), ra.data())));
  put(out, "off", vstrs(so));
  put(out, "vic", vstrs(sv));
  put(out, "count", vnums(x.v));
  put(out, "rtoff", vnums(ro));
  put(out, "rtpair", vnums(rp));
  put(out, "nulloff", vnums(no));
  put(out, "ratio", vnums(ra));
  return true;
}

// Forensics.

bool gtmp(Ctx &c, Val &out) {
  Dv p;
  CK(gnum(c, c["p"], "p", p, -1, 0.0, true));
  const std::size_t n = p.n();
  for (double x : p.v) {
    if (x > 1) return c.fail("`p` must be at most 1.");
  }
  std::vector<int> h0;
  if (c.has("hom")) {
    const Got gh = readb(c["hom"], h0);
    bool na = false;
    for (int v : h0) na = na || v < 0;
    if (gh == NOTNUM || na || (h0.size() != 1 && h0.size() != n)) {
      return c.fail("`hom` must be TRUE or FALSE, once or once per locus.");
    }
    if (gh == MATRIX) return notvec(c, "hom");
  } else {
    h0.assign(1, 0);
  }
  std::vector<int> h(n);
  for (std::size_t i = 0; i < n; ++i) h[i] = h0[h0.size() == 1 ? 0 : i];
  Dv qd;
  if (c.has("q")) {
    const Got gq = readd(c["q"], qd);
    if (gq == NOTNUM || qd.n() != n) return c.fail("`q` must be numbers, one per locus.");
    if (gq == MATRIX) return notvec(c, "q");
  } else {
    qd.v.assign(n, NaN);
  }
  std::vector<std::string> g;
  CK(align(c, {{"p", &p}, {"q", &qd}}, g));
  std::vector<double> q = qd.v;
  bool bad = false;
  for (std::size_t i = 0; i < n; ++i) {
    if (h[i]) q[i] = 0;
    bad = bad || std::isnan(q[i]) || (!h[i] && (q[i] <= 0 || p[i] + q[i] > 1 + 1e-12));
  }
  if (bad) return c.fail("`q` must be above 0 with `p + q` at most 1 at every heterozygous locus.");
  double th = 0;
  if (c.has("theta")) {
    th = val1(c["theta"]);
    if (std::isnan(th) || th < 0 || th >= 1) {
      return c.fail("`theta` must be one number from 0 up to, not including, 1.");
    }
  }
  std::vector<double> mp(n);
  double l10 = 0;
  CK(kchk(c, clab_gtmp(p.p(), q.data(), h.data(), static_cast<int>(n), th, mp.data(), &l10)));
  Val loc, hl = varr();
  for (int v : h) hl.a.push_back(vbool(v != 0));
  put(loc, "locus", vstrs(grp(g, n)));
  put(loc, "hom", hl);
  put(loc, "mp", vnums(mp));
  put(out, "locus", loc);
  put(out, "log10mp", vnum(l10));
  put(out, "log10lr", vnum(-l10));
  return true;
}

bool cllr(Ctx &c, Val &out) {
  Dv s, d;
  CK(num(c, c["same"], "same", s, -1, 0.0, true));
  CK(num(c, c["diff"], "diff", d, -1, 0.0, true));
  double r = 0;
  CK(kchk(c, clab_cllr(s.p(), static_cast<int>(s.n()), d.p(), static_cast<int>(d.n()), &r)));
  out = vnum(r);
  return true;
}

bool lrpost(Ctx &c, Val &out) {
  Dv l, p;
  CK(num(c, c["lr"], "lr", l, -1, 0.0, true));
  CK(num(c, c["prior"], "prior", p, -1, 0.0, true));
  bool bad = p.n() != 1 && p.n() != l.n();
  for (std::size_t i = 0; i < p.n() && !bad; ++i) bad = p[i] >= 1;
  if (bad) return c.fail("`prior` must be below 1, once or once per ratio.");
  std::vector<double> o(l.n());
  for (std::size_t i = 0; i < l.n(); ++i) {
    const double q = p[p.n() == 1 ? 0 : i];
    o[i] = 1 / (1 + (1 - q) / (q * l[i]));
  }
  out = vnums(o);
  return true;
}

// .mix(): one side of a mixture ratio.
bool mix(Ctx &c, const Dv &p, const char *an, const char *bn, double &r) {
  std::vector<std::string> need;
  bool bad = !labs(c[an], need);
  for (std::size_t i = 0; i < need.size() && !bad; ++i) {
    bad = std::find(p.nm.begin(), p.nm.end(), need[i]) == p.nm.end();
  }
  std::sort(need.begin(), need.end());
  need.erase(std::unique(need.begin(), need.end()), need.end());
  if (bad || need.size() > 20) return c.fail("`%s` must be names of `p`, at most 20.", an);
  const double v = val1(c[bn]);
  if (std::isnan(v) || v < 0 || v != std::floor(v) || v > 50) {
    return c.fail("`%s` must be a whole number from 0 to 50.", bn);
  }
  std::vector<int> in(p.n());
  for (std::size_t i = 0; i < p.n(); ++i) {
    in[i] = std::binary_search(need.begin(), need.end(), p.nm[i]);
  }
  return kchk(c, clab_mixp(p.p(), in.data(), static_cast<int>(p.n()), static_cast<int>(v), &r));
}

// p: an object from allele to frequency.
bool mixlr(Ctx &c, Val &out) {
  Dv p;
  CK(num(c, c["p"], "p", p, -1, 0.0, true));
  std::vector<std::string> s = p.nm;
  std::sort(s.begin(), s.end());
  if (p.nm.empty() || std::adjacent_find(s.begin(), s.end()) != s.end() ||
      nsum(p.p(), p.n()) > 1 + 1e-12) {
    return c.fail("`p` must be uniquely named with a sum of at most 1.");
  }
  double r[3];
  CK(mix(c, p, "np", "xp", r[0]));
  CK(mix(c, p, "nd", "xd", r[1]));
  r[2] = r[0] / r[1];
  static const char *const nm[] = {"num", "den", "lr"};
  out = nvec(3, nm, r);
  return true;
}

bool dbsrch(Ctx &c, Val &out) {
  Dv m, nn;
  CK(num(c, c["mp"], "mp", m, 1, 0.0, true));
  const double mv = m[0];
  if (mv > 1) return c.fail("`mp` must be at most 1.");
  CK(num(c, c["n"], "n", nn, 1, 1.0));
  const double nv = nn[0];
  double r[3];
  const double e = nv * mv;
  r[0] = e < 1 ? e : 1;
  r[1] = -std::expm1(nv * std::log1p(-mv));
  const bool has = c.has("pop");
  if (has) {
    Dv pp;
    CK(num(c, c["pop"], "pop", pp, 1, nv));
    r[2] = 1 / (1 + rd((pp[0] - nv) * mv));
  }
  static const char *const nm[] = {"np", "any", "post"};
  out = nvec(has ? 3 : 2, nm, r);
  return true;
}

// .cal(): the ratios checked, and calibrated in the order given.
struct Cal {
  Dv s, d;
  std::vector<double> cs, cd;
  double min = 0;
};

bool cal(Ctx &c, Cal &k) {
  CK(num(c, c["same"], "same", k.s, -1, 0.0, true));
  CK(num(c, c["diff"], "diff", k.d, -1, 0.0, true));
  const std::size_t ns = k.s.n(), nd = k.d.n(), n = ns + nd;
  if (n > INT_MAX) return c.fail("`same` is too long for crimlab: more than %d values.", INT_MAX);
  std::vector<double> lr(k.s.v);
  lr.insert(lr.end(), k.d.v.begin(), k.d.v.end());
  std::vector<std::size_t> o(n);
  for (std::size_t i = 0; i < n; ++i) o[i] = i;
  std::stable_sort(o.begin(), o.end(), [&lr](std::size_t i, std::size_t j) { return lr[i] < lr[j]; });
  std::vector<int> is(n);
  std::vector<double> ls(n), post(n), work(2 * n);
  for (std::size_t i = 0; i < n; ++i) {
    is[i] = o[i] < ns;
    ls[i] = lr[o[i]];
  }
  CK(kchk(c, clab_lrcal(ls.data(), is.data(), static_cast<int>(n), post.data(), work.data(), &k.min)));
  k.cs.resize(ns);
  k.cd.resize(nd);
  for (std::size_t i = 0; i < n; ++i) {
    const double x = post[i] / (1 - post[i]) * static_cast<double>(nd) / static_cast<double>(ns);
    if (o[i] < ns) k.cs[o[i]] = x;
    else k.cd[o[i] - ns] = x;
  }
  return true;
}

bool cllrmin(Ctx &c, Val &out) {
  Cal k;
  CK(cal(c, k));
  double r[3];
  CK(kchk(c, clab_cllr(k.s.p(), static_cast<int>(k.s.n()), k.d.p(), static_cast<int>(k.d.n()), r)));
  r[1] = k.min;
  r[2] = r[0] - k.min;
  static const char *const nm[] = {"cllr", "min", "cal"};
  out = nvec(3, nm, r);
  return true;
}

bool lrcal(Ctx &c, Val &out) {
  Cal k;
  CK(cal(c, k));
  put(out, "same", vnums(k.cs));
  put(out, "diff", vnums(k.cd));
  return true;
}

// one curve of tipplt(): sorted log10 ratios and the share at or above
Val curve(const Dv &x) {
  const std::size_t n = x.n();
  std::vector<double> l(n), s(n);
  for (std::size_t i = 0; i < n; ++i) l[i] = std::log10(x[i]);
  std::sort(l.begin(), l.end());
  for (std::size_t i = 0; i < n; ++i) s[i] = static_cast<double>(n - i) / static_cast<double>(n);
  Val out;
  put(out, "log10lr", vnums(l));
  put(out, "share", vnums(s));
  return out;
}

bool tipplt(Ctx &c, Val &out) {
  Dv s, d;
  CK(num(c, c["same"], "same", s, -1, 0.0, true));
  CK(num(c, c["diff"], "diff", d, -1, 0.0, true));
  put(out, "same", curve(s));
  put(out, "diff", curve(d));
  return true;
}

// log(1 + e^z), without overflow for large z
double softplus(double z) { return z > 0 ? z + std::log1p(std::exp(-z)) : std::log1p(std::exp(z)); }

// One empirical cross-entropy, as src/bind.cpp computes it.
double ece(const double *s, std::size_t ns, const double *d, std::size_t nd, double lo,
           double odds, double pr, double *t) {
  const double lg = lo * std::log(10.0);
  const bool fin = std::isfinite(odds) && odds > 0;
  for (std::size_t k = 0; k < ns; ++k) {
    const double x = fin && std::isfinite(s[k]) && s[k] > 0 ? s[k] * odds : 0;
    if (std::isfinite(x) && x > 0) {
      t[k] = x < 1 ? std::log1p(x) - std::log(x) : std::log1p(1 / x);
    } else {
      t[k] = softplus(-(std::log(s[k]) + lg));
    }
  }
  const double a = nmean(t, ns);
  for (std::size_t k = 0; k < nd; ++k) {
    const double x = fin && std::isfinite(d[k]) && d[k] > 0 ? d[k] * odds : 0;
    t[k] = std::isfinite(x) && x > 0 ? std::log1p(x) : softplus(std::log(d[k]) + lg);
  }
  const double b = nmean(t, nd);
  return (rd(pr * a) + rd((1 - pr) * b)) / std::log(2.0);
}

// lo: log10 prior odds; by default -2.5 to 2.5 by 0.05, as in Python.
bool eceplt(Ctx &c, Val &out) {
  Cal k;
  CK(cal(c, k));
  Dv l;
  if (c.has("lo")) {
    CK(num(c, c["lo"], "lo", l));
  } else {
    // the product rounded before the sum, as Python computes it: a fused
    // multiply-add would give other grid points (1.4e-16 for 0)
    for (int i = 0; i <= 100; ++i) l.v.push_back(std::min(-2.5 + rd(i * 0.05), 2.5));
  }
  const std::size_t nl = l.n(), ns = k.s.n(), nd = k.d.n();
  std::vector<double> e(nl), m(nl), r(nl), t(ns > nd ? ns : nd);
  const double unit = 1.0;
  for (std::size_t i = 0; i < nl; ++i) {
    const double odds = rpow(10.0, l[i]);
    const double pr = std::isinf(odds) ? 1.0 : odds / (1 + odds);
    e[i] = ece(k.s.p(), ns, k.d.p(), nd, l[i], odds, pr, t.data());
    m[i] = ece(k.cs.data(), ns, k.cd.data(), nd, l[i], odds, pr, t.data());
    r[i] = ece(&unit, 1, &unit, 1, l[i], odds, pr, t.data());
  }
  put(out, "lo", vnums(l.v));
  put(out, "ece", vnums(e));
  put(out, "min", vnums(m));
  put(out, "ref", vnums(r));
  return true;
}

// A genotype table: one pair of alleles, an array of pairs (one per
// locus) or an object from locus to pair. Labels (num false) or numbers
// (num true), two per row, and the loci when the table names them.
struct Tab {
  std::vector<std::string> lab;
  std::vector<double> x;
  std::vector<std::string> loci;
  std::size_t rows() const { return (lab.empty() ? x.size() : lab.size()) / 2; }
};

bool tab1(const Val &v, bool num, Tab &t) {
  t = Tab();
  auto pair = [&](const Val &r) {
    if (r.t != Val::ARR || r.a.size() != 2) return false;
    for (const Val &e : r.a) {
      std::string s;
      if (num ? e.t != Val::NUM : !lab(e, s)) return false;
      if (num) t.x.push_back(e.n);
      else t.lab.push_back(s);
    }
    return true;
  };
  if (v.t == Val::OBJ) {
    t.loci = v.k;
  } else if (v.t != Val::ARR || v.a.empty()) {
    return false;
  } else if (v.a[0].t != Val::ARR) {
    return pair(v);
  }
  for (const Val &r : v.a) {
    if (!pair(r)) return false;
  }
  return true;
}

// A failed read leaves t empty: no loci without their rows.
bool tab(const Val &v, bool num, Tab &t) {
  if (tab1(v, num, t)) return true;
  t = Tab();
  return false;
}

bool gtab(Ctx &c, const char *an, const char *nm, Tab &t) {
  if (!tab(c[an], false, t)) {
    return c.fail("`%s` must be a genotype: a pair of allele labels, an array of pairs (one per "
                  "locus) or an object from locus to pair.", nm);
  }
  return true;
}

bool nrows(Ctx &c, const char *nm, const Tab &t, std::size_t n) {
  if (t.rows() != n) return c.fail("`%s` must have %d rows, not %d.", nm, static_cast<int>(n), static_cast<int>(t.rows()));
  return true;
}

// The tables matched by locus: the loci of the first that names them;
// every other table that names its loci is put in their order.
bool loci(Ctx &c, std::vector<std::pair<const char *, Tab *>> ts, std::vector<std::string> &loc) {
  loc.clear();
  const char *ln = nullptr;
  for (auto &a : ts) {
    Tab &t = *a.second;
    if (t.loci.empty()) continue;
    if (ln == nullptr) {
      ln = a.first;
      loc = t.loci;
      std::vector<std::string> s = loc;
      std::sort(s.begin(), s.end());
      const auto d = std::adjacent_find(s.begin(), s.end());
      if (d != s.end()) return c.fail("`%s` names a locus twice: %s.", ln, d->c_str());
      continue;
    }
    bool ok = t.loci.size() == loc.size();
    Tab u = t;
    for (std::size_t i = 0; ok && i < loc.size(); ++i) {
      const auto j = std::find(t.loci.begin(), t.loci.end(), loc[i]);
      ok = j != t.loci.end();
      if (!ok) break;
      const std::size_t k = static_cast<std::size_t>(j - t.loci.begin());
      for (std::size_t h = 0; h < 2; ++h) {
        if (t.lab.empty()) u.x[2 * i + h] = t.x[2 * k + h];
        else u.lab[2 * i + h] = t.lab[2 * k + h];
      }
    }
    if (!ok) return c.fail("`%s` must name the same loci as `%s`.", a.first, ln);
    u.loci = loc;
    t = u;
  }
  return true;
}

// The alleles of the tables as codes, one per distinct label, from 1.
std::vector<std::vector<int>> codes(const std::vector<const Tab *> &ts) {
  std::vector<std::string> u;
  std::vector<std::vector<int>> out;
  for (const Tab *t : ts) {
    std::vector<int> v;
    for (const std::string &s : t->lab) {
      auto j = std::find(u.begin(), u.end(), s);
      if (j == u.end()) j = u.insert(u.end(), s);
      v.push_back(static_cast<int>(j - u.begin()) + 1);
    }
    out.push_back(v);
  }
  return out;
}

// .gf(): the frequencies of the alleles of the table g, read by tab()
// into f (ok false when freq is no table of pairs of numbers).
bool gf(Ctx &c, const Tab &f, bool ok, const std::vector<int> &g, const char *gn) {
  const std::size_t nr = g.size() / 2;
  bool bad = !ok || f.x.size() != 2 * nr;
  for (std::size_t i = 0; !bad && i < f.x.size(); ++i) bad = std::isnan(f.x[i]) || f.x[i] <= 0 || f.x[i] > 1;
  if (bad) {
    return c.fail("`freq` must hold the frequencies of the alleles of `%s`, same shape, above 0 "
                  "and at most 1.", gn);
  }
  for (std::size_t i = 0; i < nr; ++i) {
    if (g[2 * i] == g[2 * i + 1] && f.x[2 * i] != f.x[2 * i + 1]) {
      return c.fail("`freq` gives two frequencies for one allele of a homozygote.");
    }
  }
  return true;
}

// .kin(): the identity-by-descent probabilities of a relationship.
bool kin(Ctx &c, double *k) {
  static const char *const nm[] = {"parent", "sibling", "half", "cousin"};
  static const double ks[4][3] = {{0, 1, 0}, {0.25, 0.5, 0.25}, {0.5, 0.5, 0}, {0.75, 0.25, 0}};
  const Val &rel = c["rel"];
  if (rel.t == Val::NUL || rel.t == Val::STR) {
    const std::string s = rel.t == Val::NUL ? "parent" : rel.s;
    for (int i = 0; i < 4; ++i) {
      if (s == nm[i]) {
        std::memcpy(k, ks[i], sizeof ks[i]);
        return true;
      }
    }
  } else {
    Dv v;
    if (readd(rel, v) == NUM && v.n() == 3) {
      bool ok = true;
      for (int i = 0; i < 3; ++i) ok = ok && !std::isnan(v[static_cast<std::size_t>(i)]) && v[static_cast<std::size_t>(i)] >= 0;
      if (ok && std::fabs(nsum(v.p(), 3) - 1) < 1e-9) {
        std::memcpy(k, v.p(), 3 * sizeof(double));
        return true;
      }
    }
  }
  return c.fail("`rel` must be one of parent, sibling, half, cousin, or three probabilities that "
                "sum to 1.");
}

// .idx(): the per-locus indices with their loci, their log10 total and the
// posterior at even prior odds.
void idx(const std::vector<std::string> &loc, const std::vector<double> &ix, Val &out) {
  const std::size_t nr = ix.size();
  std::vector<double> lg(nr);
  for (std::size_t i = 0; i < nr; ++i) lg[i] = std::log10(ix[i]);
  const double l10 = nsum(lg.data(), nr);
  Val l;
  put(l, "locus", vstrs(grp(loc, nr)));
  put(l, "idx", vnums(ix));
  put(out, "locus", l);
  put(out, "log10idx", vnum(l10));
  put(out, "post", vnum(1 / (1 + rpow(10.0, -l10))));
}

bool kinidx(Ctx &c, Val &out) {
  Tab a, b, f;
  CK(gtab(c, "g1", "g1", a));
  CK(gtab(c, "g2", "g2", b));
  CK(nrows(c, "g2", b, a.rows()));
  const bool fok = tab(c["freq"], true, f);
  std::vector<std::string> loc;
  CK(loci(c, {{"g1", &a}, {"g2", &b}, {"freq", &f}}, loc));
  const auto k = codes({&a, &b});
  CK(gf(c, f, fok, k[0], "g1"));
  double r[3];
  CK(kin(c, r));
  std::vector<double> ix(a.rows());
  CK(kchk(c, clab_kinidx(k[0].data(), k[1].data(), f.x.data(), static_cast<int>(a.rows()), r, ix.data())));
  idx(loc, ix, out);
  return true;
}

bool patidx(Ctx &c, Val &out) {
  const bool hasm = c.has("mother");
  // without the mother: kinidx(child, father, freq, "parent"), its
  // messages naming g1 and g2
  const char *cn = hasm ? "child" : "g1", *fn = hasm ? "father" : "g2";
  Tab a, b, m, f;
  CK(gtab(c, "child", cn, a));
  CK(gtab(c, "father", fn, b));
  CK(nrows(c, fn, b, a.rows()));
  if (hasm) {
    CK(gtab(c, "mother", "mother", m));
    CK(nrows(c, "mother", m, a.rows()));
  }
  const bool fok = tab(c["freq"], true, f);
  std::vector<std::string> loc;
  CK(loci(c, {{cn, &a}, {fn, &b}, {"mother", &m}, {"freq", &f}}, loc));
  const auto k = hasm ? codes({&a, &b, &m}) : codes({&a, &b});
  CK(gf(c, f, fok, k[0], cn));
  const int nr = static_cast<int>(a.rows());
  std::vector<double> ix(a.rows());
  if (hasm) {
    CK(kchk(c, clab_patidx(k[0].data(), k[2].data(), k[1].data(), f.x.data(), nr, ix.data())));
  } else {
    const double r[3] = {0, 1, 0};
    CK(kchk(c, clab_kinidx(k[0].data(), k[1].data(), f.x.data(), nr, r, ix.data())));
  }
  idx(loc, ix, out);
  return true;
}

// Policing.

bool clrsmp(Ctx &c, Val &out) {
  Dv ca, ra, cb, rb;
  CK(gnum(c, c["clra"], "clra", ca, -1, 0.0));
  const int n = static_cast<int>(ca.n());
  CK(gnum(c, c["reca"], "reca", ra, n, 0.0, true));
  CK(gnum(c, c["clrb"], "clrb", cb, n, 0.0));
  CK(gnum(c, c["recb"], "recb", rb, n, 0.0, true));
  std::vector<std::string> g;
  CK(align(c, {{"clra", &ca}, {"reca", &ra}, {"clrb", &cb}, {"recb", &rb}}, g));
  std::vector<double> v[6];
  for (auto &x : v) x.resize(ca.n());
  double tot[8];
  int cnt[3];
  CK(kchk(c, clab_clrsmp(ca.p(), ra.p(), cb.p(), rb.p(), n, v[0].data(), v[1].data(), v[2].data(),
                         v[3].data(), v[4].data(), v[5].data(), tot, cnt)));
  static const char *const pnm[] = {"ratea", "rateb", "sharea", "shareb", "rateeff", "mixeff"};
  Val parts;
  put(parts, "part", vstrs(grp(g, ca.n())));
  for (int k = 0; k < 6; ++k) put(parts, pnm[k], vnums(v[k]));
  for (int k = 0; k < 3; ++k) tot[5 + k] = cnt[k];
  static const char *const tnm[] = {"ratea", "rateb", "change", "rateeff",
                                    "mixeff", "up",    "down",   "simpson"};
  put(out, "parts", parts);
  put(out, "total", nvec(8, tnm, tot));
  return true;
}

bool expbnd(Ctx &c, Val &out) {
  Dv y, m, gm;
  CK(gnum(c, c["y"], "y", y, -1, 0.0, true));
  const int n = static_cast<int>(y.n());
  CK(gnum(c, c["m"], "m", m, n, 0.0, true));
  CK(num(c, c["gamma"], "gamma", gm, 1, 1.0));
  std::vector<std::string> g;
  CK(align(c, {{"y", &y}, {"m", &m}}, g));
  g = grp(g, y.n());
  int i = 1;
  CK(ref(c, c["ref"], g, i));
  std::vector<double> v[6];
  for (auto &x : v) x.resize(y.n());
  std::vector<int> sg(y.n());
  CK(kchk(c, clab_expbnd(y.p(), m.p(), n, i - 1, gm[0], v[0].data(), v[1].data(), v[2].data(),
                         v[3].data(), v[4].data(), v[5].data(), sg.data())));
  static const char *const nm[] = {"rate", "lo", "hi", "ratio", "rlo", "rhi"};
  put(out, "group", vstrs(g));
  put(out, "ref", refcol(n, i));
  for (int k = 0; k < 6; ++k) put(out, nm[k], vnums(v[k]));
  put(out, "sign", vints(sg));
  return true;
}

bool misrd(Ctx &c, Val &out) {
  double r1 = 0, r0 = 0, p = 0, s = 0, sp = 0;
  CK(p01(c, c["r1"], "r1", r1));
  CK(p01(c, c["r0"], "r0", r0));
  CK(one(c, c["prev"], "prev", 0, 1, p));
  CK(p01(c, c["se"], "se", s));
  CK(p01(c, c["sp"], "sp", sp));
  const double d1 = rd(p * s) + rd((1 - p) * (1 - sp));
  if (d1 <= 0 || d1 >= 1) return c.fail("`se` and `sp` leave a recorded group empty.");
  double r[5];
  CK(kchk(c, clab_misrd(r1, r0, p, s, sp, r)));
  static const char *const nm[] = {"true", "measured", "factor", "d1", "d0"};
  out = nvec(5, nm, r);
  return true;
}

bool rtmhs(Ctx &c, Val &out) {
  Dv x1, x2, ct, w;
  CK(num(c, c["x1"], "x1", x1));
  const int n = static_cast<int>(x1.n());
  CK(num(c, c["x2"], "x2", x2, n));
  CK(num(c, c["cut"], "cut", ct, 1));
  CK(wts(c, w, x1.n()));
  bool any = false;
  for (std::size_t i = 0; i < x1.n() && !any; ++i) any = x1[i] > ct[0] && w[i] > 0;
  if (!any) return c.fail("no place with weight exceeds `cut`.");
  double r[7];
  CK(kchk(c, clab_rtmhs(x1.p(), x2.p(), w.p(), n, ct[0], r)));
  static const char *const nm[] = {"nsel", "change", "mirror", "sym", "excess", "low", "lowsym"};
  out = nvec(7, nm, r);
  return true;
}

// Corrections.

bool incap(Ctx &c, Val &out) {
  Dv l, q, s, sh;
  CK(gnum(c, c["lam"], "lam", l, -1, 0.0, true));
  const std::size_t n = l.n();
  CK(gnum(c, c["q"], "q", q, -1, 0.0, true));
  for (double x : q.v) {
    if (x > 1) return c.fail("`q` must be at most 1.");
  }
  CK(gnum(c, c["s"], "s", s, -1, 0.0));
  if (q.n() != 1 && q.n() != n) return c.fail("`q` must have length 1 or %d.", static_cast<int>(n));
  if (s.n() != 1 && s.n() != n) return c.fail("`s` must have length 1 or %d.", static_cast<int>(n));
  const bool has = c.has("share");
  if (has) {
    CK(gnum(c, c["share"], "share", sh, static_cast<long>(n), 0.0));
    if (std::fabs(nsum(sh.p(), n) - 1) > 1e-9) return c.fail("`share` must sum to 1.");
  }
  std::vector<std::string> g;
  CK(align(c, {{"lam", &l}, {"q", q.n() == n ? &q : nullptr}, {"s", s.n() == n ? &s : nullptr},
               {"share", has ? &sh : nullptr}}, g));
  g = grp(g, n);
  const std::size_t nr = has ? n + 1 : n;
  std::vector<double> v[6];
  for (auto &x : v) x.resize(nr);
  for (std::size_t i = 0; i < n; ++i) {
    v[0][i] = l[i];
    v[1][i] = q[q.n() == 1 ? 0 : i];
    v[2][i] = s[s.n() == 1 ? 0 : i];
  }
  double agg[3];
  CK(kchk(c, clab_incap(v[0].data(), v[1].data(), v[2].data(), static_cast<int>(n),
                        has ? sh.p() : nullptr, v[3].data(), v[4].data(), v[5].data(), agg)));
  if (has) {
    v[0][n] = agg[0];
    v[1][n] = NA;
    v[2][n] = NA;
    v[3][n] = agg[1];
    v[4][n] = agg[2];
    v[5][n] = NA;
    g.push_back("all");
  }
  static const char *const nm[] = {"lam", "q", "s", "rate", "prev", "marg"};
  put(out, "group", vstrs(g));
  for (int k = 0; k < 6; ++k) put(out, nm[k], vnums(v[k]));
  return true;
}

bool incapr(Ctx &c, Val &out) {
  Dv l;
  CK(num(c, c["lam"], "lam", l, -1, 0.0));
  for (std::size_t i = 1; i < l.n(); ++i) {
    if (l[i] - l[i - 1] > 0) return c.fail("`lam` must never increase.");
  }
  int sv = 0, tv = 0;
  CK(whole(c, c["s"], "s", 1, sv));
  CK(whole(c, c["t0"], "t0", 0, tv));
  if (static_cast<double>(tv) + sv > static_cast<double>(l.n()) - 1) {
    return c.fail("`lam` must have length at least t0 + s + 1.");
  }
  double rv = 0;
  if (c.has("r")) CK(p01(c, c["r"], "r", rv));
  double o[5];
  CK(kchk(c, clab_incapr(l.p(), static_cast<int>(l.n()), tv, sv, rv, o)));
  static const char *const nm[] = {"prev", "net", "upper", "lower", "later"};
  out = nvec(5, nm, o);
  return true;
}

bool sellab(Ctx &c, Val &out) {
  std::vector<int> rel, rule;
  CK(lgl(c, c["rel"], "rel", -1, rel));
  const std::size_t n = rel.size();
  CK(lgl(c, c["rule"], "rule", static_cast<long>(n), rule));
  Dv y;
  const Got gy = readd(c["y"], y);
  if (gy == NOTNUM || y.n() != n) return c.fail("`y` must be numbers, one per defendant.");
  if (gy == MATRIX) return notvec(c, "y");
  std::vector<double> yv = y.v;
  for (std::size_t i = 0; i < n; ++i) {
    if (rel[i] && !(yv[i] == 0 || yv[i] == 1)) {
      return c.fail("`y` must be 0 or 1 for every released defendant.");
    }
    if (!rel[i]) yv[i] = 0;
  }
  Dv w;
  CK(wts(c, w, n));
  if (!pos(w, &rel)) return c.fail("the judge released nobody.");
  if (!pos(w, &rule)) return c.fail("the rule releases nobody.");
  double r[6];
  CK(kchk(c, clab_sellab(yv.data(), rel.data(), rule.data(), w.p(), static_cast<int>(n), r)));
  static const char *const nm[] = {"obs", "unobs", "lower", "upper", "width", "naive"};
  out = nvec(6, nm, r);
  return true;
}

bool losbnd(Ctx &c, Val &out) {
  Dv d, o;
  CK(num(c, c["done"], "done", d, -1, 0.0));
  CK(num(c, c["open"], "open", o, -1, 0.0));
  double r[4];
  CK(kchk(c, clab_losbnd(d.p(), static_cast<int>(d.n()), o.p(), static_cast<int>(o.n()), r)));
  static const char *const nm[] = {"mean", "age", "lower", "bias"};
  out = nvec(4, nm, r);
  return true;
}

bool sntbnd(Ctx &c, Val &out) {
  Dv y;
  CK(num(c, c["y"], "y", y, -1, 0.0));
  const std::size_t n = y.n();
  for (double x : y.v) {
    if (!(x == 0 || x == 1)) return c.fail("`y` must be 0 or 1.");
  }
  // the labels of z, and its kind: 'l' (all bools), 'n' (numbers, bools
  // among them) or 's'
  const Val &zv = c["z"];
  std::vector<std::string> z;
  bool bad = zv.t != Val::ARR || zv.a.size() != n;
  bool alll = true, alln = true;
  for (std::size_t i = 0; !bad && i < n; ++i) {
    std::string s;
    bad = !lab(zv.a[i], s);
    alll = alll && zv.a[i].t == Val::BOOL;
    alln = alln && zv.a[i].t != Val::STR;
    z.push_back(s);
  }
  if (bad) return c.fail("`z` must have length %d and no missing values.", static_cast<int>(n));
  const char kind = alll ? 'l' : alln ? 'n' : 's';
  // lv <- sort(unique(z)): exactly two labels, in code-point order
  std::vector<std::string> lv = z;
  std::sort(lv.begin(), lv.end());
  lv.erase(std::unique(lv.begin(), lv.end()), lv.end());
  if (lv.size() != 2) return c.fail("`z` must take exactly two values.");
  std::size_t tr = 1;
  if (!c.has("treat")) {
    const bool b = kind == 'n' && lv[0] == "0" && lv[1] == "1";
    if (kind != 'l' && !b) {
      return c.fail("`treat` must name the treatment, %s or %s: it is chosen by default only for a "
                    "logical `z` (TRUE) or a `z` of 0 and 1 (1).",
                    lv[0].c_str(), lv[1].c_str());
    }
  } else {
    std::string t;
    if (!lab(c["treat"], t) || (t != lv[0] && t != lv[1])) {
      return c.fail("`treat` must be one of the values of `z`.");
    }
    tr = t == lv[0] ? 0 : 1;
  }
  Dv w;
  CK(wts(c, w, n));
  if (!pos(w)) return c.fail("`w` must have a positive total.");
  std::vector<int> zi(n);
  for (std::size_t i = 0; i < n; ++i) zi[i] = z[i] == lv[tr];
  double r[11];
  CK(kchk(c, clab_sntbnd(y.p(), zi.data(), w.p(), static_cast<int>(n), r)));
  const double col[4][3] = {{r[2], r[3], NA},
                            {r[4], r[6], r[8]},
                            {r[5], r[7], r[9]},
                            {r[0] / r[2], r[1] / r[3], r[10]}};
  static const char *const nm[] = {"share", "lower", "upper", "observed"};
  put(out, "term", vstrs({lv[1 - tr], lv[tr], "effect"}));
  for (int k = 0; k < 4; ++k) put(out, nm[k], vnums(col[k], 3));
  return true;
}

bool hzdsel(Ctx &c, Val &out) {
  double s = 0, h = 0, l = 0;
  CK(one(c, c["s"], "s", 0, 1, s));
  CK(p01(c, c["h"], "h", h));
  CK(p01(c, c["l"], "l", l));
  if (l >= h) return c.fail("`l` must be below `h`.");
  Dv sa, sb;
  CK(gnum(c, c["surh"], "surh", sa, 2, 0.0, true));
  CK(gnum(c, c["surl"], "surl", sb, 2, 0.0, true));
  for (double x : sa.v) {
    if (x > 1) return c.fail("`surh` must be at most 1.");
  }
  for (double x : sb.v) {
    if (x > 1) return c.fail("`surl` must be at most 1.");
  }
  double r[5];
  CK(kchk(c, clab_hzdsel(s, h, l, sa.p(), sb.p(), r)));
  static const char *const nm[] = {"share0", "share1", "haz0", "haz1", "ratio"};
  out = nvec(5, nm, r);
  return true;
}

// .pick() (src/spin.cpp): one of the n strings opt, into k.
bool pick(Ctx &c, const Val &o, const char *nm, const char *const *opt, int n, int &k) {
  if (o.t == Val::STR) {
    for (k = 0; k < n; ++k) {
      if (o.s == opt[k]) return true;
    }
  }
  std::string all;
  for (int j = 0; j < n; ++j) all += std::string(j ? ", \"" : "\"") + opt[j] + "\"";
  return c.fail("`%s` must be one of %s.", nm, all.c_str());
}

// byname() (src/spin.cpp): when every argument given is an object, all of
// one length with names that do not repeat, and two or more are given,
// each later one is put in the order of the first, which must name the
// same units.
bool byname(Ctx &c, std::vector<std::pair<const char *, Dv *>> as) {
  std::vector<std::size_t> g;
  for (std::size_t j = 0; j < as.size(); ++j) {
    if (c.has(as[j].first)) g.push_back(j);
  }
  if (g.size() < 2) return true;
  const Dv &f = *as[g[0]].second;
  for (std::size_t j : g) {
    std::vector<std::string> s = as[j].second->nm;
    if (s.empty() || s.size() != f.n() || as[j].second->n() != f.n()) return true;
    std::sort(s.begin(), s.end());
    if (std::adjacent_find(s.begin(), s.end()) != s.end()) return true;
  }
  for (std::size_t q = 1; q < g.size(); ++q) {
    Dv &d = *as[g[q]].second;
    std::vector<double> v(f.n());
    for (std::size_t i = 0; i < f.n(); ++i) {
      const auto it = std::find(d.nm.begin(), d.nm.end(), f.nm[i]);
      if (it == d.nm.end()) {
        return c.fail("`%s` must name the same units as `%s`; `%s` has no %s.", as[g[q]].first,
                      as[g[0]].first, as[g[q]].first, f.nm[i].c_str());
      }
      v[i] = d.v[static_cast<std::size_t>(it - d.nm.begin())];
    }
    d.v = v;
    d.nm = f.nm;
  }
  return true;
}

// Fairness: the false positive rates the other two rates imply.
bool impfpr(Ctx &c, Val &out) {
  static const char *const un[] = {"p", "ppv", "fnr"};
  Dv v[3];
  for (int t = 0; t < 3; ++t) readd(c[un[t]], v[t]);
  CK(byname(c, {{un[0], &v[0]}, {un[1], &v[1]}, {un[2], &v[2]}}));
  std::size_t n = 1;
  for (int t = 0; t < 3; ++t) {
    const std::vector<double> keep = v[t].v;
    const std::vector<std::string> nm = v[t].nm;
    CK(num(c, c[un[t]], un[t], v[t]));
    v[t].v = keep;  // in the order byname() put them
    v[t].nm = nm;
    n = std::max(n, v[t].n());
  }
  for (int t = 0; t < 3; ++t) {
    if (v[t].n() != 1 && v[t].n() != n) {
      return c.fail("`p`, `ppv` and `fnr` must have one length (or length 1).");
    }
  }
  std::vector<double> x[3];
  for (int t = 0; t < 3; ++t) {
    x[t].resize(n);
    for (std::size_t i = 0; i < n; ++i) x[t][i] = v[t][v[t].n() == 1 ? 0 : i];
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (!(x[0][i] > 0.0 && x[0][i] < 1.0)) return c.fail("`p` must be base rates in (0, 1).");
    if (!(x[1][i] > 0.0 && x[1][i] <= 1.0)) return c.fail("`ppv` must be in (0, 1].");
    if (!(x[2][i] >= 0.0 && x[2][i] < 1.0)) return c.fail("`fnr` must be in [0, 1).");
  }
  std::vector<double> o(n);
  CK(kchk(c, clab_impfpr(x[0].data(), x[1].data(), x[2].data(), static_cast<int>(n), o.data())));
  for (int t = 0; t < 3; ++t) {
    if (v[t].n() == n && !v[t].nm.empty()) {
      out = Val();
      out.t = Val::OBJ;
      for (std::size_t i = 0; i < n; ++i) put(out, v[t].nm[i].c_str(), vnum(o[i]));
      return true;
    }
  }
  out = vnums(o);
  return true;
}

// Policing: Kolesar and Blum's square-root law.
bool sqrtlw(Ctx &c, Val &out) {
  static const char *const ms[] = {"euclid", "grid"};
  Dv ar, un, bu;
  CK(num(c, c["area"], "area", ar, 1, 0.0, true));
  CK(num(c, c["units"], "units", un, 1, 0.0, true));
  if (c.has("busy")) {
    CK(num(c, c["busy"], "busy", bu, 1, 0.0));
  } else {
    bu.v.push_back(0.0);
  }
  if (!(bu[0] < 1)) return c.fail("`busy` must be below 1: some units must be free.");
  int mt = 0;
  if (c.has("metric")) CK(pick(c, c["metric"], "metric", ms, 2, mt));
  double d, k;
  CK(kchk(c, clab_sqrtlw(ar[0], un[0], bu[0], mt, &d, &k)));
  double tm = NA;
  if (c.has("speed")) {
    Dv sp;
    CK(num(c, c["speed"], "speed", sp, 1, 0.0, true));
    tm = d / sp[0];
  }
  out = Val();
  out.t = Val::OBJ;
  put(out, "dist", vnum(d));
  put(out, "const", vnum(k));
  put(out, "avail", vnum(un[0] * (1 - bu[0])));
  put(out, "time", vnum(tm));
  return true;
}

// The self-test (src/selftest.cpp): checks run, checks failed, report.
bool selftest(Ctx &, Val &out) {
  std::vector<char> sink(1 << 16);
  int nf = 0;
  const int n = clab_selftest(sink.data(), static_cast<int>(sink.size()), &nf);
  put(out, "n", vnum(n));
  put(out, "nfail", vnum(nf));
  put(out, "report", vstr(sink.data()));
  return true;
}

// The machine-checked theorems behind each function, from src/bind.cpp
// (cli/thm.inc, made by cli/thm.sh).
struct Thm {
  const char *fn, *mod, *thm;
};
const Thm thms[] = {
#include "thm.inc"
};

bool proof(Ctx &c, Val &out) {
  std::vector<std::string> fs;
  const Val &f = c["fn"];
  if (c.has("fn")) {
    bool ok = f.t == Val::STR || f.t == Val::ARR;
    if (f.t == Val::STR) fs.push_back(f.s);
    for (std::size_t i = 0; ok && i < f.a.size(); ++i) {
      ok = f.a[i].t == Val::STR;
      if (ok) fs.push_back(f.a[i].s);
    }
    if (!ok) return c.fail("`fn` must be function names.");
    std::vector<std::string> miss;
    for (const std::string &s : fs) {
      bool known = false;
      for (const Thm &t : thms) known = known || s == t.fn;
      if (!known && std::find(miss.begin(), miss.end(), s) == miss.end()) miss.push_back(s);
    }
    if (!miss.empty()) {
      std::string m = miss[0];
      for (std::size_t i = 1; i < miss.size(); ++i) m += ", " + miss[i];
      return c.fail("no proof recorded for: %s.", m.c_str());
    }
  }
  Val fn = varr(), th = varr(), fl = varr();
  for (const Thm &t : thms) {
    if (!fs.empty() && std::find(fs.begin(), fs.end(), t.fn) == fs.end()) continue;
    fn.a.push_back(vstr(t.fn));
    th.a.push_back(vstr(std::string("Crimlab.") + t.mod + "." + t.thm));
    fl.a.push_back(vstr(std::string("lean/") + t.mod + ".lean"));
  }
  put(out, "fn", fn);
  put(out, "theorem", th);
  put(out, "file", fl);
  return true;
}

// The families of the R package's other functions, each in its file.
#include "jx.inc"
#include "jfair.inc"
#include "jtrace.inc"
#include "jaor.inc"
#include "jsens.inc"
#include "jhot.inc"
#include "jtrj.inc"
#include "jnet.inc"

// The table: name, arguments (? for an optional one), entry point.
using Fn = bool (*)(Ctx &, Val &);
struct Ent {
  const char *name, *args;
  Fn f;
};
const Ent ents[] = {
    {"vodtst", "day night ref? shift?", vodtst},
    {"orrisk", "or rate share", orrisk},
    {"outtst", "q t w?", outtst},
    {"dspbm", "pop contact force ref?", dspbm},
    {"dyadrt", "count pop", dyadrt},
    {"gtmp", "p q? hom? theta?", gtmp},
    {"cllr", "same diff", cllr},
    {"lrpost", "lr prior", lrpost},
    {"mixlr", "p np xp nd xd", mixlr},
    {"dbsrch", "mp n pop?", dbsrch},
    {"cllrmin", "same diff", cllrmin},
    {"lrcal", "same diff", lrcal},
    {"tipplt", "same diff", tipplt},
    {"eceplt", "same diff lo?", eceplt},
    {"kinidx", "g1 g2 freq rel?", kinidx},
    {"patidx", "child father freq mother?", patidx},
    {"clrsmp", "clra reca clrb recb", clrsmp},
    {"expbnd", "y m gamma ref?", expbnd},
    {"misrd", "r1 r0 prev se sp", misrd},
    {"rtmhs", "x1 x2 cut w?", rtmhs},
    {"incap", "lam q s share?", incap},
    {"incapr", "lam t0 s r?", incapr},
    {"sellab", "y rel rule w?", sellab},
    {"losbnd", "done open", losbnd},
    {"sntbnd", "y z treat? w?", sntbnd},
    {"hzdsel", "s h l surh surl", hzdsel},
    {"impfpr", "p ppv fnr", impfpr},
    {"sqrtlw", "area units busy? metric? speed?", sqrtlw},
    {"fair", "score y grp cut? base? conf?", fair},
    {"fairdf", "pred grp? y? base? conf? nboot? seed?", fairdf},
    {"tmpaud", "pred grp? period? city? base?", tmpaud},
    {"lipaud", "h d? pairs? lip?", lipaud},
    {"nsbnd", "pobs a b", nsbnd},
    {"rnkres", "est hw", rnkres},
    {"nordet", "x off radius p? draw? seed?", nordet},
    {"fairsim", "n? groups? props? nareas? base? bias? seed?", fairsim},
    {"oaxbd", "y x grp ref? conf?", oaxbd},
    {"tlvc", "x item", tlvc},
    {"tlvlr", "xa xb na nb bg", tlvlr},
    {"skp", "x s2 alpha? nsim? seed?", skp},
    {"lwss", "x y f? iter? delta?", lwss},
    {"chmb", "a b wopt wval coarse?", chmb},
    {"tfrp", "par deffect? tdist? ymax? nsim? seed? hl? nq?", tfrp},
    {"aorwt", "from to bw? nb? origin? off? tt? toff? lo? hi? each?", aorwt},
    {"bdig", "x nd? dp? sgn?", bdig},
    {"bnft", "dig nd? nsim? seed?", bnft},
    {"rptc", "unit year rep y0 y1 minrep?", rptc},
    {"rswp", "d gam? alpha?", rswp},
    {"rshl", "d gam?", rshl},
    {"gstr", "cnt ex? rad? nsim? seed?", gstr},
    {"sptst", "a b nboot? seed? alpha? scale?", sptst},
    {"zgll", "y x z? par prob off?", zgll},
    {"zgem", "y x z? par prob off? maxit? tol?", zgem},
    {"mvtll", "y x beta sig prob", mvtll},
    {"mvtem", "y x beta sig prob maxit? tol?", mvtem},
    {"ntdeg", "from to w? nodes? dir?", ntdeg},
    {"ntbtw", "from to nodes? dir? norm?", ntbtw},
    {"ntcls", "from to nodes? dir?", ntcls},
    {"ntevc", "from to w? nodes? tol? maxit?", ntevc},
    {"ntcore", "from to nodes?", ntcore},
    {"ntcmp", "from to nodes? dir?", ntcmp},
    {"ntlvn", "from to w? nodes? res? seed?", ntlvn},
    {"ntmod", "from to memb w? nodes? res?", ntmod},
    {"ntbrk", "from to grp nodes? dir?", ntbrk},
    {"ntchg", "from to w? nodes? tol? maxit?", ntchg},
    {"coofg", "inc off", coofg},
    {"selftest", "", selftest},
    {"proof", "fn?", proof},
};

std::vector<std::string> words(const char *s) {
  std::vector<std::string> w;
  std::string t;
  for (const char *p = s;; ++p) {
    if (*p == ' ' || *p == '\0') {
      if (!t.empty()) w.push_back(t);
      t.clear();
      if (*p == '\0') return w;
    } else {
      t += *p;
    }
  }
}

int run(const char *fn, const char *in, std::string &res) {
  const Ent *e = nullptr;
  for (const Ent &x : ents) {
    if (fn != nullptr && std::strcmp(fn, x.name) == 0) e = &x;
  }
  if (e == nullptr) {
    res = std::string("crimlab: no function ") + (fn ? fn : "(null)") + ".";
    return CLAB_EARG;
  }
  Val v;
  v.t = Val::OBJ;
  if (in != nullptr && !Parser(in, std::strlen(in)).doc(v)) v.t = Val::NUL;
  Ctx c{e->name, &v, std::string()};
  if (v.t != Val::OBJ) {
    c.fail("the arguments must be one JSON object.");
    res = c.err;
    return CLAB_EARG;
  }
  const std::vector<std::string> as = words(e->args);
  for (std::size_t i = 0; i < v.k.size(); ++i) {
    const std::string &k = v.k[i];
    bool ok = false;
    for (const std::string &a : as) ok = ok || a == k || a == k + "?";
    if (!ok) {
      c.fail("unknown argument `%s`; the arguments are: %s.", k.c_str(), e->args[0] ? e->args : "none");
      res = c.err;
      return CLAB_EARG;
    }
    if (std::find(v.k.begin(), v.k.begin() + static_cast<std::ptrdiff_t>(i), k) != v.k.begin() + static_cast<std::ptrdiff_t>(i)) {
      c.fail("`%s` is given twice.", k.c_str());
      res = c.err;
      return CLAB_EARG;
    }
  }
  for (const std::string &a : as) {
    if (a.back() != '?' && !c.has(a.c_str())) {
      c.fail("`%s` is missing.", a.c_str());
      res = c.err;
      return CLAB_EARG;
    }
  }
  Val out;
  if (!e->f(c, out)) {
    res = c.err;
    return CLAB_EARG;
  }
  res.clear();
  wval(res, out);
  return CLAB_OK;
}

}  // namespace

extern "C" int clab_jrun(const char *fn, const char *in, char *out, int cap, int *need) {
  std::string res;
  int rc;
#if defined(__cpp_exceptions) || defined(_CPPUNWIND)
  // no exception crosses the C interface (std::bad_alloc is the only one)
  try {
    rc = run(fn, in, res);
  } catch (...) {
    res = "crimlab: out of memory.";
    rc = CLAB_EARG;
  }
#else
  rc = run(fn, in, res);
#endif
  const std::size_t n = res.size() + 1;
  if (need != nullptr) *need = n > INT_MAX ? INT_MAX : static_cast<int>(n);
  if (out == nullptr || cap < 1) return CLAB_ESIZE;
  if (n > static_cast<std::size_t>(cap)) {
    out[0] = '\0';
    return CLAB_ESIZE;
  }
  std::memcpy(out, res.c_str(), n);
  return rc;
}

extern "C" const char *clab_jfns(void) {
  static const std::string s = [] {
    Val o;
    o.t = Val::OBJ;
    for (const Ent &e : ents) put(o, e.name, vstrs(words(e.args)));
    std::string r;
    wval(r, o);
    return r;
  }();
  return s.c_str();
}
