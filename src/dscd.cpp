// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Data description between the arguments and the results (dscd.h): the
// sketches, the column profiles, the outlier, run and correlation tables,
// the cleaned names, the missingness map and patterns, and the groups of a
// distinctions test. The counting and the arithmetic are in dsc.cpp; the R
// binding (bdscd.cpp) and clab_jrun() (cli/jdsc.inc) both call these.
#include "dscd.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

using std::int64_t;
using std::size_t;

constexpr double qnan = std::numeric_limits<double>::quiet_NaN();

inline size_t z(int64_t i) { return static_cast<size_t>(i); }

// the finite values of x[0..n)
std::vector<double> fin(const double *x, int64_t n) {
  std::vector<double> f;
  for (int64_t i = 0; i < n; ++i) {
    if (std::isfinite(x[i])) f.push_back(x[i]);
  }
  return f;
}

// the sketch of finite values f in bins characters
std::string sketch(const std::vector<double> &f, int bins, bool uni) {
  const size_t b = static_cast<size_t>(bins);
  if (f.empty()) return std::string(b, ' ');
  const auto mm = std::minmax_element(f.begin(), f.end());
  if (*mm.first == *mm.second) {
    return std::string(b / 2, ' ') + (uni ? "\xe2\x96\x88" : "#") + std::string(b - b / 2 - 1, ' ');
  }
  std::vector<int> cnt(b, 0);
  if (clab_hbin(f.data(), static_cast<int>(f.size()), bins, cnt.data()) != CLAB_OK) return std::string(b, ' ');
  static const char *const up[] = {"\xe2\x96\x81", "\xe2\x96\x82", "\xe2\x96\x83", "\xe2\x96\x84",
                                   "\xe2\x96\x85", "\xe2\x96\x86", "\xe2\x96\x87", "\xe2\x96\x88"};
  static const char *const as[] = {".", ":", "-", "=", "#"};
  const int nl = uni ? 8 : 5;
  const int mx = *std::max_element(cnt.begin(), cnt.end());
  std::string r;
  for (int k : cnt) {
    if (k == 0) {
      r += ' ';
      continue;
    }
    const double h = std::ceil(static_cast<double>(k) / mx * nl);
    const int i = std::max(1, static_cast<int>(h));
    r += uni ? up[i - 1] : as[i - 1];
  }
  return r;
}

std::string low(std::string s) {
  for (char &c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

std::string upc(std::string s) {
  for (char &c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

bool alnum(char c) {
  return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9');
}

// one name of clncolnms() before the numbering of repeats
std::string cname(const std::string &s0, int cs, const std::string &sep) {
  // a control character (a tab, a line break) is a word break, which the
  // transliteration would delete
  std::string s(s0);
  for (char &c : s) {
    if (static_cast<unsigned char>(c) < 0x20 || c == '\x7f') c = ' ';
  }
  std::vector<char> b(s.size() * 2 + 1, 0);
  const int m = s.size() > INT_MAX / 2 ? 0 : clab_ascii(s.data(), static_cast<int>(s.size()), b.data());
  std::string u;
  for (int i = 0; i < m; ++i) {
    const char c = b[static_cast<size_t>(i)];
    if (c == '%') u += " pct ";
    else if (c == '&') u += " and ";
    else if (c == '#') u += " n ";
    else u += c;
  }
  std::vector<std::string> w(1);
  for (size_t i = 0; i < u.size(); ++i) {
    if (alnum(u[i])) w.back() += u[i];
    else if (!w.back().empty()) w.emplace_back();
    const bool lo = (u[i] >= 'a' && u[i] <= 'z') || (u[i] >= '0' && u[i] <= '9');
    if (lo && i + 1 < u.size() && u[i + 1] >= 'A' && u[i + 1] <= 'Z') w.emplace_back();
  }
  if (w.back().empty()) w.pop_back();
  if (w.empty()) w.emplace_back("x");
  std::string o;
  for (size_t i = 0; i < w.size(); ++i) {
    const std::string cap = upc(w[i].substr(0, 1)) + low(w[i].substr(1));
    const std::string jn = i ? sep : std::string();
    if (cs == 0) o += jn + low(w[i]);
    else if (cs == 3) o += jn + upc(w[i]);
    else if (cs == 4) o += jn + w[i];
    else if (cs == 1) o += i ? cap : low(w[i]);
    else o += cap;
  }
  if (!o.empty() && o[0] >= '0' && o[0] <= '9') o = "x" + o;
  return o;
}

// the bytes of the first k characters of UTF-8 text s, one per element
std::vector<std::string> chars(const std::string &s, int k) {
  std::vector<std::string> r;
  for (size_t i = 0; i < s.size() && static_cast<int>(r.size()) < k;) {
    const unsigned char c = static_cast<unsigned char>(s[i]);
    const size_t len = c >= 0xF0 ? 4 : c >= 0xE0 ? 3 : c >= 0xC0 ? 2 : 1;
    r.push_back(s.substr(i, len));
    i += len;
  }
  return r;
}

std::string rep(const std::string &s, int k) {
  std::string r;
  for (int i = 0; i < k; ++i) r += s;
  return r;
}

}  // namespace

namespace dscd {

std::vector<std::string> cnames(const std::vector<std::string> &nm, int cs, const std::string &sep) {
  std::vector<std::string> out;
  for (const std::string &s : nm) out.push_back(cname(s, cs, sep));
  // each repeated name, in the order of its first repeat, numbered from 2;
  // a number whose name is taken (by a name given, or one numbered
  // before) is skipped, so that the names are unique
  // ponytail: quadratic in the names, which are columns
  std::vector<std::string> rp;
  for (size_t i = 0; i < out.size(); ++i) {
    const bool dup = std::find(out.begin(), out.begin() + static_cast<std::ptrdiff_t>(i), out[i]) !=
                     out.begin() + static_cast<std::ptrdiff_t>(i);
    if (dup && std::find(rp.begin(), rp.end(), out[i]) == rp.end()) rp.push_back(out[i]);
  }
  std::vector<std::string> tk(out);
  for (const std::string &r : rp) {
    bool first = true;
    long k = 2;
    for (std::string &o : out) {
      if (o != r) continue;
      if (first) {
        first = false;
        continue;
      }
      std::string c = r + sep + std::to_string(k);
      while (std::find(tk.begin(), tk.end(), c) != tk.end()) c = r + sep + std::to_string(++k);
      ++k;
      o = c;
      tk.push_back(c);
    }
  }
  return out;
}

std::vector<std::string> mgrid(const int *na, size_t n, size_t p, const std::vector<std::string> &nm, int height,
                               int width, bool uni) {
  const std::string h = uni ? "\xe2\x94\x80" : "-";
  const std::string top = rep(h, 2) + " Missingness map " + rep(h, 47), bot = rep(h, 66);
  if (n == 0) return {top, "  (no rows)", bot};
  std::vector<std::string> lines = {top};
  std::vector<std::vector<std::string>> ch;
  size_t mx = 0;
  for (const std::string &s : nm) {
    ch.push_back(chars(s, width));
    mx = std::max(mx, ch.back().size());
  }
  for (size_t r = 0; r < mx; ++r) {
    std::string l = "         ";
    for (const auto &c : ch) l += r < c.size() ? c[r] : " ";
    lines.push_back(l);
  }
  static const char *const su[] = {" ", "\xe2\x96\x91", "\xe2\x96\x92", "\xe2\x96\x93", "\xe2\x96\x88"};
  static const char *const sa[] = {" ", ".", "-", "+", "#"};
  const char *const *sh = uni ? su : sa;
  const int b = static_cast<int>(std::min(static_cast<size_t>(height), n)), pp = static_cast<int>(p);
  std::vector<int> st(static_cast<size_t>(b) + 1, 0), g((static_cast<size_t>(b) + 1) * p, 0);
  const int k = clab_mgrid(na, static_cast<int>(n), pp, b, st.data(), g.data());
  for (int i = 0; i < k; ++i) {
    char lb[32] = {0};
    std::snprintf(lb, sizeof lb, "  %6d ", st[static_cast<size_t>(i)]);
    std::string l = lb;
    for (size_t j = 0; j < p; ++j) l += sh[g[j * static_cast<size_t>(k) + static_cast<size_t>(i)]];
    lines.push_back(l);
  }
  lines.emplace_back("");
  lines.push_back(std::string("  legend: '") + sh[0] + "' none  '" + sh[4] + "' all missing");
  lines.push_back(bot);
  return lines;
}

pat mpat(const int *na, size_t n, size_t p, const std::vector<std::string> &nm, int maxpat) {
  pat r;
  if (n == 0) return r;
  std::vector<std::string> code(n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < p; ++j) code[i] += na[j * n + i] ? 'X' : '.';
  }
  std::vector<std::string> lv = code;
  std::sort(lv.begin(), lv.end());
  lv.erase(std::unique(lv.begin(), lv.end()), lv.end());
  std::vector<int> cd(n, 0), lev(lv.size() + 2, 0), cnt(lv.size() + 2, 0);
  for (size_t i = 0; i < n; ++i) {
    cd[i] = static_cast<int>(std::lower_bound(lv.begin(), lv.end(), code[i]) - lv.begin()) + 1;
  }
  const int k = clab_freq(cd.data(), static_cast<int>(n), static_cast<int>(lv.size()), 1, maxpat, 0, lev.data(),
                          cnt.data());
  for (int i = 0; i < k; ++i) {
    const std::string &s = lv[static_cast<size_t>(lev[static_cast<size_t>(i)] - 1)];
    const int c = cnt[static_cast<size_t>(i)];
    int m = 0;
    std::string cols;
    for (size_t j = 0; j < s.size(); ++j) {
      if (s[j] != 'X') continue;
      ++m;
      cols += (cols.empty() ? "" : ", ") + nm[j];
    }
    r.code.push_back(s);
    r.rows.push_back(c);
    r.pct.push_back(100 * static_cast<double>(c) / static_cast<double>(n));
    r.nmiss.push_back(m);
    r.cols.push_back(cols);
  }
  return r;
}

sel dtsel(const std::vector<std::string> &g, const std::vector<char> &gna, const double *y, const double *lo,
          const double *hi, const std::vector<std::string> &q, const std::vector<char> &qna, bool hasp,
          const std::vector<std::string> &pooled, double z) {
  const size_t n = g.size();
  sel r;
  r.y.assign(y, y + n);
  r.s.assign(n, qnan);
  r.use.assign(n, 0);
  bool nabad = false;
  for (size_t i = 0; i < n; ++i) {
    const double l = lo[i], h = hi[i];
    if (std::isfinite(l) && std::isfinite(h) && std::isfinite(y[i]) && (l > y[i] || h < y[i])) {
      if (gna[i] && !nabad) {
        nabad = true;
        r.bad.emplace_back();
        r.bna.push_back(1);
      } else if (!gna[i] && std::find(r.bad.begin(), r.bad.end(), g[i]) == r.bad.end()) {
        r.bad.push_back(g[i]);
        r.bna.push_back(0);
      }
      r.y[i] = qnan;
    }
    r.s[i] = (h - l) / (2 * z);
  }
  static const char *const known[] = {"ok", "E", "F", "x", "..", "..."};
  for (size_t i = 0; i < n; ++i) {
    std::string f = q.empty() || qna[i] || q[i].empty() ? std::string("ok") : q[i];
    if (std::find(std::begin(known), std::end(known), f) == std::end(known)) {
      if (std::find(r.odd.begin(), r.odd.end(), f) == r.odd.end()) r.odd.push_back(f);
      f = "ok";
    }
    bool pan = false;
    if (gna[i]) {
      pan = false;
    } else if (hasp) {
      pan = std::find(pooled.begin(), pooled.end(), g[i]) != pooled.end();
    } else {
      const std::string t = low(g[i]);
      pan = t.rfind("total", 0) == 0 || t.rfind("indigenous identity", 0) == 0 ||
            t.rfind("aboriginal identity", 0) == 0 || t.rfind("indigenous population", 0) == 0;
    }
    const bool sup = f == "F" || f == "x" || f == ".." || f == "...";
    r.use[i] = !pan && !sup && std::isfinite(r.y[i]) && std::isfinite(r.s[i]) && r.s[i] > 0;
  }
  return r;
}

std::string mcnote(const double *o) {
  std::string note;
  if (o[1] <= 0) {
    char b[64] = {0};
    std::snprintf(b, sizeof b, "%.15g", o[1]);
    note = std::string("df is ") + b +
           ": with a single missingness pattern there is nothing to test. A complete table trivially satisfies MCAR.";
  }
  if (o[7] > 0) {
    note += (note.empty() ? "" : " ") + std::to_string(static_cast<int64_t>(o[7])) +
            " pattern(s) were skipped: their observed variables were collinear.";
  }
  return note;
}

}  // namespace dscd

extern "C" {

int clab_spark(const double *x, int64_t n, int bins, int uni, char *out, int64_t *len) CLAB_TRY {
  *len = 0;
  if (n < 0 || n > INT_MAX || bins < 1 || bins > kutil::maxcell) return CLAB_EARG;
  const std::string s = sketch(fin(x, n), bins, uni != 0);
  std::memcpy(out, s.data(), s.size());
  *len = static_cast<int64_t>(s.size());
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_prof(const double *x, int64_t n, const double *pr, int np, double *st, double *q) CLAB_TRY {
  if (n < 0 || n > INT_MAX || np < 0) return CLAB_EARG;
  st[0] = 0;
  for (int k = 1; k < 12; ++k) st[k] = qnan;
  for (int j = 0; j < np; ++j) q[j] = qnan;
  std::vector<double> ok;
  for (int64_t i = 0; i < n; ++i) {
    if (!std::isnan(x[i])) ok.push_back(x[i]);
  }
  if (ok.empty()) return CLAB_OK;
  double nz = 0, ng = 0, ni = 0;
  for (double v : ok) {
    nz += v == 0 ? 1 : 0;
    ng += v < 0 ? 1 : 0;
    ni += std::isinf(v) ? 1 : 0;
  }
  std::vector<double> u = ok;
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  st[0] = static_cast<double>(u.size());
  st[1] = nz;
  st[2] = ng;
  st[3] = ni;
  const std::vector<double> f = fin(ok.data(), static_cast<int64_t>(ok.size()));
  if (f.empty()) return CLAB_OK;
  double s[7] = {0, 0, 0, 0, 0, 0, 0};
  clab_pstat(f.data(), static_cast<int>(f.size()), s);
  const auto mm = std::minmax_element(f.begin(), f.end());
  st[4] = s[0];
  st[5] = s[1];
  st[6] = s[3];
  st[7] = s[4];
  st[8] = *mm.first;
  st[9] = *mm.second;
  st[10] = s[2];
  double no = 0;
  for (double v : f) no += v < s[5] || v > s[6] ? 1 : 0;
  st[11] = no;
  if (np > 0) clab_q7(f.data(), static_cast<int>(f.size()), pr, np, q);
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_mhtab(const double *x, int64_t n, int p, int robust, double alpha, int *row, double *d, double *pv,
               double *lp, int *ol, int *why) CLAB_TRY {
  *why = 0;
  if (n < 0 || n > INT_MAX || p < 1) return CLAB_EARG;
  const size_t nn = z(n), pp = static_cast<size_t>(p);
  std::vector<size_t> ok;
  for (size_t i = 0; i < nn; ++i) {
    bool all = true;
    for (size_t j = 0; j < pp && all; ++j) all = !std::isnan(x[j * nn + i]);
    if (all) ok.push_back(i);
  }
  const size_t m = ok.size();
  if (m < pp + 1) {
    *why = 1;
    return CLAB_EDOM;
  }
  std::vector<double> xc(m * pp, 0.0), d2(m, 0.0), dd(m, 0.0);
  for (size_t j = 0; j < pp; ++j) {
    for (size_t i = 0; i < m; ++i) xc[j * m + i] = x[j * nn + ok[i]];
  }
  const int rc = clab_maha(xc.data(), static_cast<int>(m), p, robust, d2.data());
  if (rc != CLAB_OK) {
    *why = rc == CLAB_EDOM ? 2 : 0;
    return rc;
  }
  for (size_t i = 0; i < m; ++i) dd[i] = std::sqrt(d2[i]);
  std::vector<size_t> o(m, 0);
  std::iota(o.begin(), o.end(), size_t{0});
  std::stable_sort(o.begin(), o.end(), [&](size_t a, size_t b) { return dd[a] > dd[b]; });
  size_t k = 0;
  for (size_t i : o) {
    row[k] = static_cast<int>(ok[i]) + 1;
    d[k] = dd[i];
    pv[k] = clab_pchisq(d2[i], p, 0);
    lp[k] = clab_lpchisq(d2[i], p);
    ol[k] = pv[k] < alpha ? 1 : 0;
    ++k;
  }
  size_t c = 0;
  for (size_t i = 0; i < nn; ++i) {
    if (c < m && ok[c] == i) {
      ++c;
      continue;
    }
    row[k] = static_cast<int>(i) + 1;
    d[k] = pv[k] = lp[k] = qnan;
    ol[k] = -1;
    ++k;
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_mruns(const int *na, int64_t n, int p, int minrun, const int *rank, int *col, int *st, int *en,
               int64_t *k) CLAB_TRY {
  *k = 0;
  if (n < 0 || n > INT_MAX || p < 0 || minrun < 1) return CLAB_EARG;
  struct Run {
    int c, s, e;
  };
  std::vector<Run> rs;
  std::vector<int> a(z(n) + 1, 0), b(z(n) + 1, 0);
  for (int j = 0; j < p; ++j) {
    const int m = clab_nrun(na + static_cast<size_t>(j) * z(n), static_cast<int>(n), minrun, a.data(), b.data());
    for (int i = 0; i < m; ++i) rs.push_back({j, a[static_cast<size_t>(i)], b[static_cast<size_t>(i)]});
  }
  std::stable_sort(rs.begin(), rs.end(), [&](const Run &u, const Run &v) {
    const int lu = u.e - u.s, lv = v.e - v.s;
    return lu != lv ? lu > lv : rank[u.c] < rank[v.c];
  });
  for (size_t i = 0; i < rs.size(); ++i) {
    col[i] = rs[i].c + 1;
    st[i] = rs[i].s;
    en[i] = rs[i].e;
  }
  *k = static_cast<int64_t>(rs.size());
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_mcdrop(const double *x, int64_t n, int p, double *y, int64_t *m, int *keep, int *k, int *why) CLAB_TRY {
  *m = 0;
  *k = 0;
  *why = 0;
  if (n < 0 || n > INT_MAX || p < 0) return CLAB_EARG;
  const size_t nn = z(n);
  std::vector<size_t> cols, rows;
  for (size_t j = 0; j < static_cast<size_t>(p); ++j) {
    bool all = true;
    for (size_t i = 0; i < nn && all; ++i) all = std::isnan(x[j * nn + i]);
    if (!all) cols.push_back(j);
  }
  for (size_t i = 0; i < nn; ++i) {
    bool all = true;
    for (size_t j : cols) all = all && std::isnan(x[j * nn + i]);
    if (!all) rows.push_back(i);
  }
  bool inf = false;
  for (size_t c = 0; c < cols.size(); ++c) {
    for (size_t r = 0; r < rows.size(); ++r) {
      const double v = x[cols[c] * nn + rows[r]];
      inf = inf || std::isinf(v);
      y[c * rows.size() + r] = v;
    }
  }
  for (size_t c = 0; c < cols.size(); ++c) keep[c] = static_cast<int>(cols[c]) + 1;
  *m = static_cast<int64_t>(rows.size());
  *k = static_cast<int>(cols.size());
  if (inf) *why = 1;
  else if (cols.size() < 2 || rows.size() < 2) *why = 2;
  return *why ? CLAB_EDOM : CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_topc(const double *r, int64_t k, double minabs, int64_t top, int *idx, int64_t *m) CLAB_TRY {
  *m = 0;
  if (k < 0 || k > INT_MAX || top < 0) return CLAB_EARG;
  std::vector<int> o;
  for (int64_t i = 0; i < k; ++i) {
    if (!std::isnan(r[i]) && std::fabs(r[i]) >= minabs) o.push_back(static_cast<int>(i));
  }
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) {
    return std::fabs(r[static_cast<size_t>(a)]) > std::fabs(r[static_cast<size_t>(b)]);
  });
  if (static_cast<int64_t>(o.size()) > top) o.resize(z(top));
  for (size_t i = 0; i < o.size(); ++i) idx[i] = o[i] + 1;
  *m = static_cast<int64_t>(o.size());
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

}  // extern "C"
