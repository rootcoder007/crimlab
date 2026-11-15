// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The values and helpers of the Indigenous health, published-table and
// census functions (igh.h, igk.h), and the table of the functions.
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>

#include "crimlab.h"
#include "igk.h"

namespace igh {

const V nulv;

std::size_t V::len() const { return t == DBL ? d.size() : t == STR ? s.size() : t == LIST ? a.size() : i.size(); }

double na() {
  const std::uint64_t u = 0x7ff00000000007a2ULL;
  double x;
  std::memcpy(&x, &u, sizeof x);
  return x;
}

bool isna(double x) {
  std::uint64_t u;
  std::memcpy(&u, &x, sizeof u);
  return std::isnan(x) && (u & 0xffffffffULL) == 1954;
}

namespace {
std::string vf(const char *fmt, std::va_list ap) {
  std::va_list aq;
  va_copy(aq, ap);
  std::vector<char> b(512);
  const int n = std::vsnprintf(b.data(), b.size(), fmt, ap);
  if (n >= static_cast<int>(b.size())) {
    b.resize(static_cast<std::size_t>(n) + 1);
    std::vsnprintf(b.data(), b.size(), fmt, aq);
  }
  va_end(aq);
  return b.data();
}
}  // namespace

void stop(const char *fmt, ...) {
  std::va_list ap;
  va_start(ap, fmt);
  Stop e{vf(fmt, ap)};
  va_end(ap);
  throw e;
}

std::string sf(const char *fmt, ...) {
  std::va_list ap;
  va_start(ap, fmt);
  std::string s = vf(fmt, ap);
  va_end(ap);
  return s;
}

const V &Cx::operator[](const char *k) const {
  for (std::size_t j = 0; j < an->size(); ++j) {
    if ((*an)[j] == k) return (*av)[j];
  }
  return nulv;
}

std::string trim(const std::string &s) {
  const char *ws = " \t\r\n";
  const std::size_t b = s.find_first_not_of(ws);
  if (b == std::string::npos) return "";
  const std::size_t e = s.find_last_not_of(ws);
  return s.substr(b, e - b + 1);
}

std::string lower(const std::string &s) {
  std::string o = s;
  for (char &ch : o) {
    if (ch >= 'A' && ch <= 'Z') ch = static_cast<char>(ch - 'A' + 'a');
  }
  return o;
}

std::string rch(double x) {
  if (isna(x)) return "NA";
  if (std::isnan(x)) return "NaN";
  if (std::isinf(x)) return x > 0 ? "Inf" : "-Inf";
  if (x == 0) return "0";
  char b[64];
  const bool whole = x == std::floor(x) && std::fabs(x) <= 9007199254740992.0;
  const int top = whole ? 17 : 15;
  std::snprintf(b, sizeof b, "%.*e", top - 1, x);
  const double v = std::strtod(b, nullptr);
  int dg = 1;
  for (; dg < top; ++dg) {
    std::snprintf(b, sizeof b, "%.*e", dg - 1, x);
    if (std::strtod(b, nullptr) == v) break;
  }
  std::snprintf(b, sizeof b, "%.*e", dg - 1, x);
  const int e = std::atoi(std::strchr(b, 'e') + 1);
  const int neg = x < 0;
  const int sci = neg + dg + (dg > 1) + (std::abs(e) >= 100 ? 5 : 4);
  const int rgt = std::max(0, dg - e - 1);
  const int fix = neg + (e >= 0 ? e + 1 : 1) + (rgt > 0 ? rgt + 1 : 0);
  if (fix <= sci) std::snprintf(b, sizeof b, "%.*f", rgt, x);
  return b;
}

std::string quoted(const std::vector<std::string> &v) {
  std::string o;
  for (std::size_t j = 0; j < v.size(); ++j) o += (j ? ", \"" : "\"") + v[j] + "\"";
  return o;
}

std::string tostr(const std::vector<std::string> &v) {
  std::string o;
  for (std::size_t j = 0; j < v.size(); ++j) o += (j ? ", " : "") + v[j];
  return o;
}

Lab lab(const V &v) {
  Lab l;
  const std::size_t n = v.len();
  l.s.resize(n);
  l.na.assign(n, 0);
  for (std::size_t j = 0; j < n; ++j) {
    switch (v.t) {
      case V::STR:
        l.na[j] = v.na.empty() ? 0 : v.na[j];
        if (!l.na[j]) l.s[j] = v.s[j];
        break;
      case V::FAC:
        if (v.i[j] == INT_MIN) l.na[j] = 1;
        else l.s[j] = v.s[static_cast<std::size_t>(v.i[j] - 1)];
        break;
      case V::INT:
        if (v.i[j] == INT_MIN) l.na[j] = 1;
        else l.s[j] = std::to_string(v.i[j]);
        break;
      case V::LGL:
        if (v.i[j] == INT_MIN) l.na[j] = 1;
        else l.s[j] = v.i[j] ? "TRUE" : "FALSE";
        break;
      case V::DBL:
        if (isna(v.d[j])) l.na[j] = 1;
        else l.s[j] = rch(v.d[j]);
        break;
      default:
        break;
    }
  }
  return l;
}

Lab tlab(const V &v) {
  Lab l = lab(v);
  for (std::size_t j = 0; j < l.n(); ++j) {
    if (!l.na[j]) l.s[j] = trim(l.s[j]);
  }
  return l;
}

bool dflt(const V &v) { return v.t == V::NUL && !v.fn && !v.env && !v.odd && !v.given; }

void nofn(const V &v, const char *what) {
  if (v.fn) stop("`%s` must be data, not a function", what);
  if (v.env) stop("`%s` must be data, not an environment", what);
}

namespace {
// R's as.numeric() of one text: NA when it is not a number
double rnum(const std::string &s0) {
  const std::string s = trim(s0);
  if (s.empty() || s == "NA") return na();
  const char *p = s.c_str();
  char *e = nullptr;
  double x;
  if (s == "Inf" || s == "inf" || s == "+Inf") return HUGE_VAL;
  if (s == "-Inf" || s == "-inf") return -HUGE_VAL;
  if (s == "NaN") return std::nan("");
  x = std::strtod(p, &e);
  if (e == p || *e != '\0') return na();
  // strtod reads "infinity" and "nan(...)", as.numeric() does not
  for (char ch : s) {
    if (std::isalpha(static_cast<unsigned char>(ch)) && ch != 'e' && ch != 'E' && ch != 'x' && ch != 'X' &&
        !std::isxdigit(static_cast<unsigned char>(ch))) {
      return na();
    }
  }
  return x;
}
}  // namespace

std::vector<double> num(const V &v, const char *what, bool integer) {
  if (v.fn || v.env || v.odd || v.t == V::LIST || v.t == V::FAC) {
    const char *k = v.fn ? "a function" : v.env ? "an environment" : v.t == V::FAC ? "a factor"
                    : v.t == V::LIST ? (v.df ? "a data frame" : "a list") : "this kind of value";
    stop("`%s` must be a%s number, not %s", what, integer ? " whole" : "", k);
  }
  std::vector<double> o(v.len());
  if (v.t == V::STR) {
    for (std::size_t j = 0; j < o.size(); ++j) {
      const bool isn = !v.na.empty() && v.na[j];
      o[j] = isn ? na() : rnum(v.s[j]);
      if (!isn && isna(o[j]) && !trim(v.s[j]).empty() && v.s[j] != "NA") {
        stop("`%s` must be numeric, not text such as \"%s\"", what, v.s[j].c_str());
      }
    }
  } else if (v.t == V::LGL) {
    for (std::size_t j = 0; j < o.size(); ++j) {
      if (v.i[j] != INT_MIN) {
        stop("`%s` must be a%s number, not %s", what, integer ? " whole" : "", v.i[j] ? "TRUE" : "FALSE");
      }
      o[j] = na();
    }
  } else if (v.t == V::INT) {
    for (std::size_t j = 0; j < o.size(); ++j) o[j] = v.i[j] == INT_MIN ? na() : v.i[j];
  } else if (v.t == V::DBL) {
    o = v.d;
  }
  if (integer) {
    for (double x : o) {
      if (std::isfinite(x) && x != std::round(x)) stop("`%s` must be whole numbers", what);
    }
  }
  return o;
}

double scalar(const V &v, const char *what, double lo, double hi, bool integer, double def) {
  if (dflt(v)) return def;
  const std::vector<double> x = num(v, what, integer);
  if (x.size() != 1 || !std::isfinite(x[0]) || x[0] < lo || x[0] > hi) {
    std::string rg;
    if (std::isfinite(lo) && std::isfinite(hi)) rg = " from " + rch(lo) + " to " + rch(hi);
    else if (std::isfinite(lo)) rg = " of at least " + rch(lo);
    stop("`%s` must be one %s%s", what, integer ? "whole number" : "number", rg.c_str());
  }
  return x[0];
}

double conf(const V &v, const char *what, double def) {
  if (dflt(v)) return def;
  const std::vector<double> x = num(v, what);
  if (x.size() != 1 || std::isnan(x[0]) || x[0] <= 0 || x[0] >= 1) {
    stop("`%s` must be a single number strictly between 0 and 1", what);
  }
  return x[0];
}

std::string str1(const V &v, const char *what) {
  if (v.t != V::STR || v.s.size() != 1 || (!v.na.empty() && v.na[0]) || trim(v.s[0]).empty()) {
    stop("`%s` must be a single non-empty string", what);
  }
  return v.s[0];
}

std::string pick(const V &v, const char *what, const std::vector<std::string> &ch) {
  if (v.t == V::NUL) return ch[0];
  if (v.t == V::STR && v.s == ch) return ch[0];
  const bool one = v.t == V::STR && v.s.size() == 1 && (v.na.empty() || !v.na[0]);
  if (one) {
    for (const std::string &c : ch) {
      if (c == v.s[0]) return c;
    }
  }
  std::string sh;
  for (std::size_t j = 0; j < ch.size(); ++j) sh += (j ? ", \"" : "\"") + ch[j] + "\"";
  stop("`%s` must be one of %s, not %s", what, sh.c_str(),
       one ? ("\"" + v.s[0] + "\"").c_str() : "that value");
}

bool isstr(const V &v) {
  if (v.t != V::STR || v.s.empty()) return false;
  for (char x : v.na) {
    if (x) return false;
  }
  return true;
}

bool isdf(const V &v) {
  if (v.t != V::LIST) return false;
  if (v.df) return true;
  if (!v.js || v.nm.size() != v.a.size() || v.a.empty()) return false;
  std::size_t n = 0;
  for (std::size_t j = 0; j < v.a.size(); ++j) {
    const V &e = v.a[j];
    if (e.t == V::LIST || (e.t == V::NUL)) return false;
    if (j && e.len() != n) return false;
    n = e.len();
  }
  return true;
}

std::size_t nrow(const V &v) {
  if (v.t != V::LIST) return v.len();
  if (v.a.empty()) return v.nr > 0 ? static_cast<std::size_t>(v.nr) : 0;
  return v.a[0].len();
}

int coli(const V &df, const std::string &nm) {
  for (std::size_t j = 0; j < df.nm.size(); ++j) {
    if (df.nm[j] == nm) return static_cast<int>(j);
  }
  return -1;
}

const V *col(const V &df, const std::string &nm) {
  const int j = coli(df, nm);
  return j < 0 ? nullptr : &df.a[static_cast<std::size_t>(j)];
}

std::string colnm(const V &df, const V &x, const char *what, const char *fn) {
  const std::string s = str1(x, what);
  if (!col(df, s)) {
    (void)fn;
    stop("column `%s` not found in `data`", s.c_str());
  }
  return s;
}

V frame(const std::vector<std::string> &nm, std::vector<V> cols) {
  V f;
  f.t = V::LIST;
  f.df = true;
  f.nm = nm;
  f.a = std::move(cols);
  f.nr = f.a.empty() ? 0 : static_cast<int>(f.a[0].len());
  f.cls = {"data.frame"};
  return f;
}

V rows(const V &x, const std::vector<std::size_t> &r) {
  V o = x;
  if (x.t == V::LIST && x.df) {
    for (V &e : o.a) e = rows(e, r);
    o.nr = static_cast<int>(r.size());
    o.rn.clear();
    return o;
  }
  switch (x.t) {
    case V::DBL:
      o.d.clear();
      for (std::size_t k : r) o.d.push_back(x.d[k]);
      break;
    case V::INT:
    case V::LGL:
    case V::FAC:
      o.i.clear();
      for (std::size_t k : r) o.i.push_back(x.i[k]);
      break;
    case V::STR:
      o.s.clear();
      o.na.clear();
      for (std::size_t k : r) {
        o.s.push_back(x.s[k]);
        o.na.push_back(x.na.empty() ? 0 : x.na[k]);
      }
      break;
    default:
      break;
  }
  return o;
}

bool inh(const V &v, const char *c) {
  return std::find(v.cls.begin(), v.cls.end(), std::string(c)) != v.cls.end();
}

V dv(const std::vector<double> &x) {
  V v;
  v.t = V::DBL;
  v.d = x;
  return v;
}
V dv1(double x) { return dv(std::vector<double>(1, x)); }
V iv(const std::vector<int> &x) {
  V v;
  v.t = V::INT;
  v.i = x;
  return v;
}
V sv(const std::vector<std::string> &x) {
  V v;
  v.t = V::STR;
  v.s = x;
  v.na.assign(x.size(), 0);
  return v;
}
V sv1(const std::string &x) { return sv(std::vector<std::string>(1, x)); }
V svna(const std::vector<std::string> &x, const std::vector<char> &nas) {
  V v = sv(x);
  v.na = nas;
  return v;
}
V lv(const std::vector<int> &x) {
  V v;
  v.t = V::LGL;
  v.i = x;
  return v;
}
V fac(const std::vector<std::string> &x, const std::vector<std::string> &lev) {
  V v;
  v.t = V::FAC;
  v.s = lev;
  v.cls = {"factor"};
  for (const std::string &e : x) {
    const auto f = std::find(lev.begin(), lev.end(), e);
    v.i.push_back(f == lev.end() ? INT_MIN : static_cast<int>(f - lev.begin()) + 1);
  }
  return v;
}
V nlist(const std::vector<std::string> &nm, std::vector<V> el) {
  V v;
  v.t = V::LIST;
  v.nm = nm;
  v.a = std::move(el);
  return v;
}
V named(V v, const std::vector<std::string> &nm) {
  v.nm = nm;
  return v;
}

bool gone(const std::string &q) { return q == "F" || q == "x" || q == ".." || q == "..."; }

std::vector<std::string> qflags(const Cx &c, const Lab &q) {
  static const std::set<std::string> ok = {"ok", "E", "F", "x", "..", "..."};
  std::vector<std::string> o(q.n());
  std::vector<std::string> odd;
  for (std::size_t j = 0; j < q.n(); ++j) {
    o[j] = q.na[j] ? "" : trim(q.s[j]);
    if (o[j].empty()) o[j] = "ok";
    if (!ok.count(o[j]) && std::find(odd.begin(), odd.end(), o[j]) == odd.end()) odd.push_back(o[j]);
  }
  if (!odd.empty()) {
    c.w(sf("quality flag(s) %s not recognised: treated as ok; recode them with qltyrcd()",
           quoted(odd).c_str()));
    for (std::string &x : o) {
      if (!ok.count(x)) x = "ok";
    }
  }
  return o;
}

std::vector<std::string> quality(const Cx &c, const V &data, const std::string &qn, std::size_t n) {
  Lab q;
  if (qn.empty()) {
    q.s.assign(n, "ok");
    q.na.assign(n, 0);
  } else {
    q = lab(*col(data, qn));
  }
  return qflags(c, q);
}

void contains(const Cx &c, const char *fn, std::vector<double> &est, const std::vector<double> &lo,
              const std::vector<double> &hi, const std::vector<std::string> &labs, const std::vector<char> *use) {
  std::vector<std::string> bad;
  std::vector<std::size_t> at;
  for (std::size_t j = 0; j < est.size(); ++j) {
    if (use && !(*use)[j]) continue;
    if (std::isfinite(lo[j]) && std::isfinite(hi[j]) && std::isfinite(est[j]) && (lo[j] > est[j] || hi[j] < est[j])) {
      at.push_back(j);
      if (std::find(bad.begin(), bad.end(), labs[j]) == bad.end()) bad.push_back(labs[j]);
    }
  }
  if (at.empty()) return;
  c.w(sf("%s: the interval does not contain the estimate for: %s; not compared", fn, tostr(bad).c_str()));
  for (std::size_t j : at) est[j] = na();
}

void pooledtype(const V &pooled) {
  nofn(pooled, "pooled");
  if (pooled.t != V::STR) stop("`pooled` must be the labels of the rows that pool other groups");
  for (char x : pooled.na) {
    if (x) stop("`pooled` must be the labels of the rows that pool other groups");
  }
}

std::vector<int> pooledrows(const std::vector<std::string> &g, const V &pooled) {
  std::vector<int> o(g.size(), 0);
  if (pooled.t != V::NUL) {
    pooledtype(pooled);
    for (std::size_t j = 0; j < g.size(); ++j) {
      for (const std::string &p : pooled.s) o[j] = o[j] || trim(g[j]) == trim(p);
    }
    return o;
  }
  static const char *const pre[] = {"total", "indigenous identity", "aboriginal identity", "indigenous population"};
  for (std::size_t j = 0; j < g.size(); ++j) {
    const std::string l = trim(lower(g[j]));
    for (const char *p : pre) o[j] = o[j] || l.compare(0, std::strlen(p), p) == 0;
  }
  return o;
}

void pooledcheck(const V &pooled, const std::vector<std::string> &g) {
  if (pooled.t == V::NUL) return;
  pooledtype(pooled);
  std::vector<std::string> miss;
  for (const std::string &p : pooled.s) {
    const std::string tp = trim(p);
    bool hit = false;
    for (const std::string &x : g) hit = hit || trim(x) == tp;
    if (!hit && std::find(miss.begin(), miss.end(), tp) == miss.end()) miss.push_back(tp);
  }
  if (!miss.empty()) stop("`pooled` names %s, not a group of `data`", quoted(miss).c_str());
}

double qn(double conflvl) { return clab_qnorm(1 - (1 - conflvl) / 2); }

// R's mean(): a long double sum, divided, then corrected by the mean of
// the deviations
double rmean(const std::vector<double> &x) {
  const long double n = static_cast<long double>(x.size());
  long double s = 0;
  for (double v : x) s += v;
  s /= n;
  if (std::isfinite(static_cast<double>(s))) {
    long double t = 0;
    for (double v : x) t += v - s;
    s += t / n;
  }
  return static_cast<double>(s);
}

double rsum(const std::vector<double> &x) {
  long double s = 0;
  for (double v : x) s += v;
  return static_cast<double>(s);
}

double quant7(std::vector<double> x, double p) {
  std::vector<double> y;
  for (double v : x) {
    if (!std::isnan(v)) y.push_back(v);
  }
  if (y.empty()) return na();
  std::sort(y.begin(), y.end());
  const double n = static_cast<double>(y.size());
  const double idx = 1 + (n - 1) * p;
  const double lo = std::floor(idx), hi = std::ceil(idx);
  const double ql = y[static_cast<std::size_t>(lo) - 1], qh = y[static_cast<std::size_t>(hi) - 1];
  const double h = idx - lo;
  // as quantile(type = 7): (1 - h) * qs + h * x[hi], only where h is not 0
  double q = ql;
  if (h > 0) q = (1 - h) * ql + h * qh;
  return q;
}

double median(std::vector<double> x) {
  std::sort(x.begin(), x.end());
  const std::size_t n = x.size();
  if (n % 2) return x[n / 2];
  return rmean({x[n / 2 - 1], x[n / 2]});
}

std::vector<double> rank(const std::vector<double> &x) {
  const std::size_t n = x.size();
  std::vector<std::size_t> o(n);
  for (std::size_t j = 0; j < n; ++j) o[j] = j;
  std::stable_sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return x[a] < x[b]; });
  std::vector<double> r(n);
  for (std::size_t j = 0; j < n;) {
    std::size_t k = j;
    while (k + 1 < n && x[o[k + 1]] == x[o[j]]) ++k;
    const double av = (static_cast<double>(j) + static_cast<double>(k)) / 2 + 1;
    for (std::size_t m = j; m <= k; ++m) r[o[m]] = av;
    j = k + 1;
  }
  return r;
}

// cor() of two complete vectors as R's cov.c computes it: long double
// means with their correction pass, then the sums of products
double cor(const std::vector<double> &x, const std::vector<double> &y) {
  const std::size_t n = x.size();
  const long double ln = static_cast<long double>(n);
  long double mx = 0, my = 0;
  for (std::size_t j = 0; j < n; ++j) {
    mx += x[j];
    my += y[j];
  }
  mx /= ln;
  my /= ln;
  long double tx = 0, ty = 0;
  for (std::size_t j = 0; j < n; ++j) {
    tx += x[j] - mx;
    ty += y[j] - my;
  }
  mx += tx / ln;
  my += ty / ln;
  long double sxy = 0, sxx = 0, syy = 0;
  for (std::size_t j = 0; j < n; ++j) {
    const long double a = x[j] - mx, b = y[j] - my;
    sxy += a * b;
    sxx += a * a;
    syy += b * b;
  }
  const double cxx = static_cast<double>(sxx / (ln - 1)), cyy = static_cast<double>(syy / (ln - 1));
  const double cxy = static_cast<double>(sxy / (ln - 1));
  if (!(cxx > 0) || !(cyy > 0)) return na();
  double r = cxy / (std::sqrt(cxx) * std::sqrt(cyy));
  if (r > 1) r = 1;
  if (r < -1) r = -1;
  return r;
}

// round(x, dig) as R 4.0.0 and later: the nearer of the two candidates at
// dig digits, ties to the even last digit
double rround(double x, int dig) {
  if (std::isnan(x) || std::isinf(x) || x == 0) return x;
  double sgn = 1;
  if (x < 0) {
    sgn = -1;
    x = -x;
  }
  if (dig == 0) return sgn * std::nearbyint(x);
  const double l10x = 0.301029995663981195 * (0.5 + std::logb(x));
  if (l10x + dig > 15) return sgn * x;
  double p10 = 1;
  for (int k = 0; k < std::abs(dig); ++k) p10 *= 10;
  double xd, xu, i10;
  if (dig > 0) {
    const double x10 = p10 * x;
    i10 = std::floor(x10);
    xd = i10 / p10;
    xu = std::ceil(x10) / p10;
  } else {
    const double x10 = x / p10;
    i10 = std::floor(x10);
    xd = i10 * p10;
    xu = std::ceil(x10) * p10;
  }
  const double du = xu - x, dd = x - xd;
  return sgn * ((du < dd || (du == dd && std::fmod(i10, 2.) == 1)) ? xu : xd);
}

double plogis(double x) { return 1 / (1 + std::exp(-x)); }
double qlogis(double p) { return std::log(p / (1 - p)); }

std::vector<std::string> uniq(const std::vector<std::string> &x) {
  std::vector<std::string> o;
  std::set<std::string> seen;
  for (const std::string &s : x) {
    if (seen.insert(s).second) o.push_back(s);
  }
  return o;
}

// R's generator. set.seed(): the seed scrambled 50 times by the linear
// congruence, then each of the 625 words of the state; the first word is
// the position, 624.
Rng::Rng(std::uint32_t seed) {
  for (int j = 0; j < 50; ++j) seed = 69069 * seed + 1;
  for (int j = 0; j < 625; ++j) {
    seed = 69069 * seed + 1;
    mt_[j] = seed;
  }
  mti_ = 624;
}

double Rng::genrand() {
  static const std::uint32_t mag01[2] = {0x0u, 0x9908b0dfu};
  std::uint32_t *mt = mt_ + 1;
  const int N = 624, M = 397;
  std::uint32_t y;
  if (mti_ >= N) {
    int kk;
    for (kk = 0; kk < N - M; kk++) {
      y = (mt[kk] & 0x80000000u) | (mt[kk + 1] & 0x7fffffffu);
      mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    for (; kk < N - 1; kk++) {
      y = (mt[kk] & 0x80000000u) | (mt[kk + 1] & 0x7fffffffu);
      mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    y = (mt[N - 1] & 0x80000000u) | (mt[0] & 0x7fffffffu);
    mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];
    mti_ = 0;
  }
  y = mt[mti_++];
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680u;
  y ^= (y << 15) & 0xefc60000u;
  y ^= (y >> 18);
  return static_cast<double>(y) * 2.3283064365386963e-10;
}

double Rng::unif() {
  const double i2 = 2.328306437080797e-10;
  const double v = genrand();
  if (v <= 0.0) return 0.5 * i2;
  if (1.0 - v <= 0.0) return 1.0 - 0.5 * i2;
  return v;
}

double Rng::norm() {
  const double big = 134217728;
  double u = unif();
  u = static_cast<int>(big * u) + unif();
  return clab_qnorm(u / big);
}

double Rng::expo() {
  // q[k] = sum of ln(2)^i / i! for i = 1 .. k + 1
  double q[16];
  {
    double t = 1, s = 0;
    const double l2 = 0.693147180559945309417232121458;
    for (int k = 0; k < 16; ++k) {
      t *= l2 / (k + 1);
      s += t;
      q[k] = s;
    }
  }
  double a = 0.;
  double u = unif();
  while (u <= 0. || u >= 1.) u = unif();
  for (;;) {
    u += u;
    if (u > 1.) break;
    a += q[0];
  }
  u -= 1.;
  if (u <= q[0]) return a + u;
  int k = 0;
  double ustar = unif(), umin = ustar;
  do {
    ustar = unif();
    if (umin > ustar) umin = ustar;
    k++;
  } while (u > q[k]);
  return a + umin * q[0];
}

double Rng::rbits(int bits) {
  std::int64_t v = 0;
  for (int n = 0; n <= bits; n += 16) {
    const int v1 = static_cast<int>(std::floor(unif() * 65536));
    v = 65536 * v + v1;
  }
  const std::int64_t one64 = 1;
  return static_cast<double>(v & ((one64 << bits) - 1));
}

double Rng::index(double dn) {
  if (dn <= 0) return 0.0;
  const int bits = static_cast<int>(std::ceil(std::log2(dn)));
  double dv;
  do {
    dv = rbits(bits);
  } while (dn <= dv);
  return dv;
}

// rpois() of R's nmath (Ahrens and Dieter 1982, doi:10.1145/355993.355997):
// table lookup below 10, the normal approximation with its corrections
// from 10 up
double Rng::pois(double mu) {
  static const double a0 = -0.5, a1 = 0.3333333, a2 = -0.2500068, a3 = 0.2000118, a4 = -0.1661269,
                      a5 = 0.1421878, a6 = -0.1384794, a7 = 0.1250060;
  static const double one_7 = 0.1428571428571428571, one_12 = 0.0833333333333333333,
                      one_24 = 0.0416666666666666667;
  static const double fact[] = {1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.};
  const double sq2pi = 0.398942280401432677939946059934;
  double del, difmuk = 0., E = 0., fk = 0., fx, fy, g, omega, px, py, t, u = 0., v, x, xx, pois = -1.;
  int j, k, kflag;
  bool newbig = false;
  if (!std::isfinite(mu) || mu < 0) return std::nan("");
  if (mu <= 0.) return 0.;
  const bool big = mu >= 10.;
  if (!(big && mu == muprev_)) {
    if (big) {
      newbig = true;
      muprev_ = mu;
      s_ = std::sqrt(mu);
      d_ = 6. * mu * mu;
      bigl_ = std::floor(mu - 1.1484);
    } else {
      if (mu != muprev_) {
        muprev_ = mu;
        m_ = std::max(1, static_cast<int>(mu));
        l_ = 0;
        q_ = p0_ = p_ = std::exp(-mu);
      }
      for (;;) {
        u = unif();
        if (u <= p0_) return 0.;
        if (l_ > 0) {
          j = (u > 0.458) ? std::min(l_, m_) : 1;
          for (k = j; k <= l_; k++) {
            if (u <= pp_[k]) return static_cast<double>(k);
          }
          if (l_ == 35) continue;
        }
        l_++;
        for (k = l_; k <= 35; k++) {
          p_ *= mu / k;
          q_ += p_;
          pp_[k] = q_;
          if (u <= q_) {
            l_ = k;
            return static_cast<double>(k);
          }
        }
        l_ = 35;
      }
    }
  }
  g = mu + s_ * norm();
  if (g >= 0.) {
    pois = std::floor(g);
    if (pois >= bigl_) return pois;
    fk = pois;
    difmuk = mu - fk;
    u = unif();
    if (d_ * u >= difmuk * difmuk * difmuk) return pois;
  }
  if (newbig || mu != muprev2_) {
    muprev2_ = mu;
    omega_ = sq2pi / s_;
    b1_ = one_24 / mu;
    b2_ = 0.3 * b1_ * b1_;
    c3_ = one_7 * b1_ * b2_;
    c2_ = b2_ - 15. * c3_;
    c1_ = b1_ - 6. * b2_ + 45. * c3_;
    c0_ = 1. - b1_ + 3. * b2_ - 324. * c3_;
    c_ = 0.1069 / mu;
  }
  omega = omega_;
  bool fstep = false;
  if (g >= 0.) {
    kflag = 0;
    fstep = true;
  }
  for (;;) {
    if (!fstep) {
      E = expo();
      u = 2 * unif() - 1.;
      t = 1.8 + (u >= 0 ? std::fabs(E) : -std::fabs(E));
      if (!(t > -0.6744)) continue;
      pois = std::floor(mu + s_ * t);
      fk = pois;
      difmuk = mu - fk;
      kflag = 1;
    }
    fstep = false;
    if (pois < 10) {
      px = -mu;
      py = std::pow(mu, pois) / fact[static_cast<int>(pois)];
    } else {
      del = one_12 / fk;
      del = del * (1. - 4.8 * del * del);
      v = difmuk / fk;
      if (std::fabs(v) <= 0.25) {
        px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v + a0) - del;
      } else {
        px = fk * std::log(1. + v) - difmuk - del;
      }
      py = sq2pi / std::sqrt(fk);
    }
    x = (0.5 - difmuk) / s_;
    xx = x * x;
    fx = -0.5 * xx;
    fy = omega * (((c3_ * xx + c2_) * xx + c1_) * xx + c0_);
    if (kflag > 0) {
      if (c_ * std::fabs(u) <= py * std::exp(px + E) - fy * std::exp(fx + E)) break;
    } else if (fy - u * fy <= py * std::exp(px - fx)) {
      break;
    }
  }
  return pois;
}

bool seedarg(const V &v, std::uint32_t &s) {
  if (v.t == V::NUL) return false;
  const bool numt = v.t == V::DBL || v.t == V::INT;
  const double x = numt && v.len() == 1 ? (v.t == V::DBL ? v.d[0] : (v.i[0] == INT_MIN ? na() : v.i[0])) : na();
  if (!numt || !std::isfinite(x) || x != std::round(x) || x < -2147483647.0 || x > 2147483647.0) {
    stop("`seed` must be a single whole number from -2147483647 to 2147483647");
  }
  s = static_cast<std::uint32_t>(static_cast<std::int32_t>(x));
  return true;
}

}  // namespace igh
