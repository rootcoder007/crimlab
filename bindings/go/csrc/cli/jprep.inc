// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// clprep(): the preset layer of the R package (R/clprep.R) for clab_jrun().
// The presets are the rows of inst/presets/presets.tsv, embedded as
// pre.inc (made by pre.sh; CI checks that it is up to date). A data frame
// is a JSON object of columns; each column is numbers, text or true and
// false, null missing. Every recode, filter and recipe of the R function
// is done here in the same order and with the same messages, so the
// result equals R's on every preset (cli/test.tsv holds a case for each).
// Text is ordered by its bytes, as R orders it under LC_COLLATE=C: where
// R runs in another locale, the levels of a dum recode and the order of
// the person-year recipes can differ for labels that sort differently
// there (mixed case, punctuation). Letters are lowered as ASCII.

namespace pp {

const char *const ptsv[][9] = {
#include "pre.inc"
};

// A column: numbers, text, true/false or a factor. Missing: NA (or NaN)
// in d, sna[i] for text.
struct Col {
  enum K { NUM, STR, LGL, FAC } k = LGL;
  std::vector<double> d;         // NUM; LGL 0 or 1; FAC the code from 1
  std::vector<std::string> s;    // STR
  std::vector<char> sna;         // STR: 1 where missing
  std::vector<std::string> lev;  // FAC
  std::size_t size() const { return k == STR ? s.size() : d.size(); }
  bool miss(std::size_t i) const { return k == STR ? sna[i] != 0 : std::isnan(d[i]); }
};

struct Tbl {
  std::vector<std::string> nm;
  std::vector<Col> c;
  std::size_t n = 0;
  const Col *col(const std::string &k) const {
    for (std::size_t i = 0; i < nm.size(); ++i) {
      if (nm[i] == k) return &c[i];
    }
    return nullptr;
  }
  // d$k <- v: replaced where it is, otherwise added last
  void set(const std::string &k, Col v) {
    for (std::size_t i = 0; i < nm.size(); ++i) {
      if (nm[i] == k) {
        c[i] = std::move(v);
        return;
      }
    }
    nm.push_back(k);
    c.push_back(std::move(v));
  }
};

Col ncol(std::size_t n, double v = NA) {
  Col r;
  r.k = Col::NUM;
  r.d.assign(n, v);
  return r;
}

// as.character(x[i]) of a value that is not missing
std::string lab(const Col &x, std::size_t i) {
  switch (x.k) {
    case Col::STR: return x.s[i];
    case Col::NUM: return rfmt(x.d[i]);
    case Col::LGL: return x.d[i] != 0 ? "TRUE" : "FALSE";
    case Col::FAC: return x.lev[static_cast<std::size_t>(x.d[i]) - 1];
  }
  return std::string();
}

// the rows idx of x, in that order
Col rows(const Col &x, const std::vector<std::size_t> &idx) {
  Col r;
  r.k = x.k;
  r.lev = x.lev;
  for (std::size_t i : idx) {
    if (x.k == Col::STR) {
      r.s.push_back(x.s[i]);
      r.sna.push_back(x.sna[i]);
    } else {
      r.d.push_back(x.d[i]);
    }
  }
  return r;
}

std::vector<std::size_t> which(const std::vector<char> &keep) {
  std::vector<std::size_t> w;
  for (std::size_t i = 0; i < keep.size(); ++i) {
    if (keep[i]) w.push_back(i);
  }
  return w;
}

Tbl frows(const Tbl &f, const std::vector<std::size_t> &idx) {
  Tbl g;
  g.nm = f.nm;
  for (const Col &x : f.c) g.c.push_back(rows(x, idx));
  g.n = idx.size();
  return g;
}

// trimws(): spaces, tabs, carriage returns and newlines at both ends
std::string trim(const std::string &s) {
  const char *ws = " \t\r\n";
  const std::size_t a = s.find_first_not_of(ws);
  if (a == std::string::npos) return std::string();
  return s.substr(a, s.find_last_not_of(ws) - a + 1);
}

std::string lower(std::string s) {
  for (char &ch : s) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  return s;
}

// as.numeric() of one text: NA for blank text and "NA" (no warning), the
// number when the whole text (spaces aside) is one as R reads it (decimal,
// hexadecimal, Inf, NaN), otherwise NA with bad set
double rnum(const std::string &s, bool &bad) {
  const std::string t = trim(s);
  if (t.empty() || t == "NA") return NA;
  char *e = nullptr;
  const double v = std::strtod(t.c_str(), &e);
  if (e == t.c_str() || *e != '\0') {
    bad = true;
    return NA;
  }
  return v;
}

// as.numeric(x); warn set when text failed to read as a number
Col asnum(const Col &x, bool &warn) {
  Col r = ncol(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (x.miss(i)) {
      r.d[i] = x.k == Col::NUM ? x.d[i] : NA;
    } else if (x.k == Col::STR) {
      r.d[i] = rnum(x.s[i], warn);
    } else {
      r.d[i] = x.d[i];
    }
  }
  return r;
}

// split(): strsplit(s, sep, fixed = TRUE)[[1]], whose last piece is
// dropped when it is empty
std::vector<std::string> split(const std::string &s, char sep) {
  std::vector<std::string> out;
  std::size_t a = 0;
  for (;;) {
    const std::size_t b = s.find(sep, a);
    if (b == std::string::npos) {
      if (a < s.size()) out.push_back(s.substr(a));
      return out;
    }
    out.push_back(s.substr(a, b - a));
    a = b + 1;
  }
}

std::string tostr(const std::vector<std::string> &v) {
  std::string s;
  for (std::size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + v[i];
  return s;
}

// factor(x): levels the sorted distinct values (numbers by value, text by
// its bytes, FALSE before TRUE; a factor keeps the order of its levels)
Col factor(const Col &x) {
  Col r;
  r.k = Col::FAC;
  const std::size_t n = x.size();
  std::vector<std::size_t> o;
  for (std::size_t i = 0; i < n; ++i) {
    if (!x.miss(i)) o.push_back(i);
  }
  std::stable_sort(o.begin(), o.end(), [&x](std::size_t a, std::size_t b) {
    return x.k == Col::STR ? x.s[a] < x.s[b] : x.d[a] < x.d[b];
  });
  std::vector<std::string> keys;
  for (std::size_t j = 0; j < o.size(); ++j) {
    const bool fresh = j == 0 || (x.k == Col::STR ? x.s[o[j]] != x.s[o[j - 1]] : x.d[o[j]] != x.d[o[j - 1]]);
    if (fresh) r.lev.push_back(lab(x, o[j]));
  }
  r.d.assign(n, NA);
  for (std::size_t i = 0; i < n; ++i) {
    if (x.miss(i)) continue;
    const std::string l = lab(x, i);
    r.d[i] = static_cast<double>(std::find(r.lev.begin(), r.lev.end(), l) - r.lev.begin()) + 1;
  }
  return r;
}

// .flag(): 1 for yes, y, true, t or 1 (any case, spaces aside), TRUE, a
// positive number; 0 otherwise
Col flag(const Col &x) {
  Col r = ncol(x.size(), 0);
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (x.k == Col::LGL) {
      r.d[i] = x.d[i];
    } else if (x.k == Col::NUM) {
      r.d[i] = !std::isnan(x.d[i]) && x.d[i] > 0;
    } else if (!x.miss(i)) {
      const std::string l = lower(trim(lab(x, i)));
      r.d[i] = l == "yes" || l == "y" || l == "true" || l == "t" || l == "1";
    }
  }
  return r;
}

// the number of UTF-8 characters of s, and the byte where character k
// (from 0) starts
std::size_t cpos(const std::string &s, std::size_t k) {
  std::size_t i = 0;
  for (std::size_t c = 0; i < s.size() && c < k; ++c) {
    ++i;
    while (i < s.size() && (static_cast<unsigned char>(s[i]) & 0xc0) == 0x80) ++i;
  }
  return i;
}

// the values of one argument: a vector, a dum matrix (columns m named mn),
// or NULL
struct V {
  bool nul = false, mat = false;
  Col c;
  std::vector<Col> m;
  std::vector<std::string> mn;
  std::size_t nr = 0;  // the rows of a matrix
  std::size_t len() const { return nul ? 0 : mat ? nr * m.size() : c.size(); }
};

struct Pr {
  Ctx &c;
  bool fail(const std::string &m) { return c.fail("%s", m.c_str()); }
  void coerce() { c.warn.push_back("NAs introduced by coercion"); }
};

// .remap(): the mapping rule s applied to the labels of x
bool remap(Pr &p, const Col &x, const std::string &s, Col &out) {
  const std::size_t n = x.size();
  std::vector<std::string> l(n), o(n);
  std::vector<char> ln(n), on(n, 1);
  for (std::size_t i = 0; i < n; ++i) {
    ln[i] = x.miss(i);
    if (!ln[i]) l[i] = lower(trim(lab(x, i)));
  }
  bool hasd = false;
  std::string dflt;
  for (const std::string &rule : split(s, ';')) {
    const std::vector<std::string> kv = split(rule, '=');
    const std::string key = kv.empty() ? std::string() : kv[0];
    const std::string to = kv.size() > 1 ? kv[1] : std::string();
    if (key == "*") {
      hasd = true;
      dflt = to;
      continue;
    }
    const bool has = key.compare(0, 4, "has:") == 0;
    const std::string pat = has ? lower(key.substr(4)) : std::string();
    const std::vector<std::string> keys = split(key, '|');
    for (std::size_t i = 0; i < n; ++i) {
      if (!on[i] || ln[i]) continue;
      const bool hit = has ? l[i].find(pat) != std::string::npos
                           : std::find(keys.begin(), keys.end(), l[i]) != keys.end();
      if (hit) {
        o[i] = to;
        on[i] = 0;
      }
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (!on[i] || ln[i]) continue;
    if (!hasd) return p.fail("a value of the column has no recode: " + l[i] + ".");
    if (dflt != "NA") {
      o[i] = dflt;
      on[i] = 0;
    }
  }
  // numbers when every value reads as one, text otherwise
  Col num = ncol(n);
  bool isn = true;
  for (std::size_t i = 0; i < n && isn; ++i) {
    if (on[i]) continue;
    bool bad = false;
    num.d[i] = rnum(o[i], bad);
    isn = !std::isnan(num.d[i]);
  }
  if (isn) {
    out = num;
  } else {
    out = Col();
    out.k = Col::STR;
    out.s = o;
    out.sna = on;
  }
  return true;
}

// .recode(): one argument's values from column x (null: no such column)
// under the rule s, for n rows
bool prc(Pr &p, const Col *x, const std::string &s, std::size_t n, V &v) {
  v = V();
  Col none;
  const Col &cx = x ? *x : none;
  if (s.compare(0, 6, "const:") == 0) {
    const std::string w = s.substr(6);
    if (w == "TRUE" || w == "FALSE") {
      v.c.k = Col::LGL;
      v.c.d.assign(1, w == "TRUE");
      return true;
    }
    bool bad = false;
    const double d = rnum(w, bad);
    if (std::isnan(d)) {
      v.c.k = Col::STR;
      v.c.s.assign(1, w);
      v.c.sna.assign(1, 0);
    } else {
      v.c = ncol(1, d);
    }
    return true;
  }
  if (s == "one") {
    v.c = ncol(n, 1);
    return true;
  }
  if (s.empty()) {
    if (x == nullptr) v.nul = true;
    else v.c = cx;
    return true;
  }
  bool warn = false;
  if (s == "fac") {
    v.c = factor(cx);
  } else if (s == "num" || s == "pct") {
    v.c = asnum(cx, warn);
    if (warn) p.coerce();
    if (s == "pct") {
      for (double &d : v.c.d) d /= 100;
    }
  } else if (s == "flag") {
    v.c = flag(cx);
  } else if (s == "dum") {
    const Col f = factor(cx);
    v.mat = true;
    v.nr = f.size();
    for (std::size_t k = 1; k < f.lev.size(); ++k) {
      Col m = ncol(f.size());
      for (std::size_t i = 0; i < f.size(); ++i) {
        if (!std::isnan(f.d[i])) m.d[i] = f.d[i] == static_cast<double>(k + 1);
      }
      v.m.push_back(m);
      v.mn.push_back(f.lev[k]);
    }
  } else if (s.compare(0, 4, "sub:") == 0) {
    const std::vector<std::string> ab = split(s.substr(4), ',');
    bool bad = false;
    const double a = ab.size() > 0 ? std::trunc(rnum(ab[0], bad)) : NA;
    const double b = ab.size() > 1 ? std::trunc(rnum(ab[1], bad)) : NA;
    v.c.k = Col::STR;
    for (std::size_t i = 0; i < cx.size(); ++i) {
      const bool m = cx.miss(i) || std::isnan(a) || std::isnan(b);
      std::string t;
      if (!m) {
        const std::string l = lab(cx, i);
        const double st = a < 1 ? 1 : a;
        if (st <= b) {
          const std::size_t i0 = cpos(l, static_cast<std::size_t>(st) - 1);
          const std::size_t i1 = cpos(l, static_cast<std::size_t>(b));
          if (i0 < i1) t = l.substr(i0, i1 - i0);
        }
      }
      v.c.s.push_back(t);
      v.c.sna.push_back(m);
    }
  } else {
    return remap(p, cx, s, v.c);
  }
  return true;
}

// .keep(): the rows with none of cols missing; an error names the columns
// data lacks
bool keep(Pr &p, Tbl &d, const std::vector<std::string> &cols, bool drop) {
  std::vector<std::string> miss;
  for (const std::string &k : cols) {
    if (!d.col(k) && std::find(miss.begin(), miss.end(), k) == miss.end()) miss.push_back(k);
  }
  if (!miss.empty()) return p.fail("`data` has no column " + tostr(miss) + ".");
  if (!drop) return true;
  std::vector<char> ok(d.n, 1);
  for (const std::string &k : cols) {
    const Col *x = d.col(k);
    for (std::size_t i = 0; i < d.n; ++i) ok[i] = ok[i] && !x->miss(i);
  }
  d = frows(d, which(ok));
  return true;
}

const char *const otc[] = {"UniqueIndividual_ID", "EndFiscalYear", "Region_AtTimeOfPlacement",
                           "Region_MostRecentPlacement", "MentalHealth_Alert", "SuicideRisk_Alert",
                           "SuicideWatch_Alert"};

// quantile(x, q, type = 7) of the values of x that are not missing (NA
// when there are none)
double q7(std::vector<double> x, double q) {
  x.erase(std::remove_if(x.begin(), x.end(), [](double e) { return std::isnan(e); }), x.end());
  const std::size_t n = x.size();
  if (n == 0) return NA;
  std::sort(x.begin(), x.end());
  const double ix = 1 + static_cast<double>(n - 1) * q;
  const double lo = std::floor(ix), hi = std::ceil(ix);
  const double a = x[static_cast<std::size_t>(lo) - 1], b = x[static_cast<std::size_t>(hi) - 1];
  if (!(ix > lo && b != a)) return a;
  const double h = ix - lo;
  return (1 - h) * a + h * b;
}

// .pyrows(): the rows of a placement table with none of need missing,
// ordered by person and fiscal year; their groups (person-years, in that
// order), each row's alert combination and region changes
struct Pyr {
  Tbl d;
  std::vector<std::vector<std::size_t>> g;
  std::vector<double> combo, vmrow;
};

bool cmpv(const Col &x, std::size_t a, std::size_t b, int &r) {
  // r < 0 when row a sorts first; missing values never reach here
  if (x.k == Col::STR) {
    const int c = x.s[a].compare(x.s[b]);
    r = c < 0 ? -1 : c > 0;
  } else {
    r = x.d[a] < x.d[b] ? -1 : x.d[a] > x.d[b];
  }
  return true;
}

bool pyrows(Pr &p, Tbl d, bool extra, Pyr &y) {
  std::vector<std::string> need(otc, otc + 7);
  if (extra) {
    need.push_back("Gender");
    need.push_back("Age_Category");
  }
  if (!keep(p, d, need, true)) return false;
  const Col &id = *d.col(otc[0]), &yr = *d.col(otc[1]);
  std::vector<std::size_t> o(d.n);
  for (std::size_t i = 0; i < d.n; ++i) o[i] = i;
  std::stable_sort(o.begin(), o.end(), [&id, &yr](std::size_t a, std::size_t b) {
    int r = 0;
    cmpv(id, a, b, r);
    if (r == 0) cmpv(yr, a, b, r);
    return r < 0;
  });
  y.d = frows(d, o);
  const Col &id2 = *y.d.col(otc[0]), &yr2 = *y.d.col(otc[1]);
  const Col mh = flag(*y.d.col(otc[4])), sr = flag(*y.d.col(otc[5])), sw = flag(*y.d.col(otc[6]));
  const Col &ra = *y.d.col(otc[2]), &rb = *y.d.col(otc[3]);
  std::vector<std::string> key;
  for (std::size_t i = 0; i < y.d.n; ++i) key.push_back(lab(id2, i) + "\r" + lab(yr2, i));
  std::vector<std::string> seen;
  for (std::size_t i = 0; i < y.d.n; ++i) {
    y.combo.push_back(4 * mh.d[i] + 2 * sr.d[i] + sw.d[i]);
    const std::string a = lab(ra, i), b = lab(rb, i);
    double v = a != b;
    if (i > 0 && key[i] == key[i - 1] && a != lab(ra, i - 1)) v += 1;
    y.vmrow.push_back(v);
    const std::size_t j = static_cast<std::size_t>(std::find(seen.begin(), seen.end(), key[i]) - seen.begin());
    if (j == seen.size()) {
      seen.push_back(key[i]);
      y.g.emplace_back();
    }
    y.g[j].push_back(i);
  }
  return true;
}

std::size_t ndist(const std::vector<double> &v, const std::vector<std::size_t> &ix) {
  std::vector<double> u;
  for (std::size_t i : ix) {
    if (std::find(u.begin(), u.end(), v[i]) == u.end()) u.push_back(v[i]);
  }
  return u.size();
}

// the first row of each person-year: the columns named
Tbl pfirst(const Pyr &y, std::initializer_list<const char *> cols) {
  std::vector<std::size_t> f;
  for (const auto &g : y.g) f.push_back(g[0]);
  Tbl o;
  for (const char *k : cols) {
    o.nm.push_back(k);
    o.c.push_back(rows(*y.d.col(k), f));
  }
  o.n = f.size();
  return o;
}

// .recipe(): the derived columns a preset needs
bool recipe(Pr &p, const std::string &name, const std::string &arg, Tbl &d) {
  static const char *const known[] = {"long", "acvm", "pyr", "keep", "acnv", "pairb", "pairc", "gent", "pop"};
  if (std::find_if(std::begin(known), std::end(known), [&name](const char *k) { return name == k; }) ==
      std::end(known)) {
    return p.fail("unknown recipe " + name + ".");
  }
  if (name == "long") {
    const std::vector<std::string> ab = split(arg, ',');
    if (!keep(p, d, ab, false)) return false;
    Tbl o;
    std::vector<std::size_t> twice;
    for (int r = 0; r < 2; ++r) {
      for (std::size_t i = 0; i < d.n; ++i) twice.push_back(i);
    }
    for (std::size_t j = 0; j < d.nm.size(); ++j) {
      if (std::find(ab.begin(), ab.end(), d.nm[j]) != ab.end()) continue;
      if (o.col(d.nm[j])) continue;
      o.nm.push_back(d.nm[j]);
      o.c.push_back(rows(d.c[j], twice));
    }
    // c(x, y): text when either is text, numbers when either is numbers
    const Col &x = *d.col(ab[0]), &z = *d.col(ab[1]);
    Col cnt;
    const bool str = x.k == Col::STR || x.k == Col::FAC || z.k == Col::STR || z.k == Col::FAC;
    cnt.k = str ? Col::STR : (x.k == Col::NUM || z.k == Col::NUM) ? Col::NUM : Col::LGL;
    for (const Col *e : {&x, &z}) {
      for (std::size_t i = 0; i < d.n; ++i) {
        if (str) {
          cnt.sna.push_back(e->miss(i));
          cnt.s.push_back(e->miss(i) ? std::string() : lab(*e, i));
        } else {
          cnt.d.push_back(e->miss(i) ? NA : e->d[i]);
        }
      }
    }
    o.n = 2 * d.n;
    o.set("count", cnt);
    Col w = ncol(2 * d.n, 0);
    for (std::size_t i = 0; i < d.n; ++i) w.d[i] = 1;
    o.set("with", w);
    d = o;
    return true;
  }
  if (name == "keep") return keep(p, d, split(arg, ','), true);
  if (name == "acvm" || name == "pyr" || name == "acnv") {
    Pyr y;
    if (!pyrows(p, d, name != "pyr", y)) return false;
    const std::size_t ng = y.g.size();
    Tbl o = name == "pyr" ? pfirst(y, {otc[0], otc[1]})
                            : pfirst(y, {otc[0], otc[1], "Gender", "Age_Category"});
    Col a = ncol(ng), b = ncol(ng), h = ncol(ng);
    for (std::size_t j = 0; j < ng; ++j) {
      const std::vector<std::size_t> &g = y.g[j];
      if (name == "acvm") {
        a.d[j] = static_cast<double>(ndist(y.combo, g));
        b.d[j] = 0;
        for (std::size_t i : g) b.d[j] += y.vmrow[i];
        h.d[j] = a.d[j] >= 2;
      } else if (name == "acnv") {
        const Col mh = flag(*y.d.col(otc[4])), sr = flag(*y.d.col(otc[5])), sw = flag(*y.d.col(otc[6]));
        const Col &ra = *y.d.col(otc[2]), &rb = *y.d.col(otc[3]);
        a.d[j] = 0;
        b.d[j] = 0;
        for (std::size_t i : g) {
          a.d[j] = std::max(a.d[j], mh.d[i] + sr.d[i] + sw.d[i]);
          b.d[j] = std::max(b.d[j], static_cast<double>(lab(ra, i) != lab(rb, i)));
        }
        h.d[j] = a.d[j] >= 2;
      } else {
        a.d[j] = static_cast<double>(g.size());
        std::vector<std::string> u;
        for (std::size_t i : g) {
          for (const char *k : {otc[2], otc[3]}) {
            const std::string l = lab(*y.d.col(k), i);
            if (std::find(u.begin(), u.end(), l) == u.end()) u.push_back(l);
          }
        }
        b.d[j] = static_cast<double>(u.size());
        h.d[j] = static_cast<double>(ndist(y.combo, g));
      }
    }
    const bool v = name == "acvm";
    o.set(v ? "ac" : name == "acnv" ? "acn" : "np", a);
    o.set(v ? "vm" : name == "acnv" ? "vmany" : "nreg", b);
    o.set(v ? "hiac" : name == "acnv" ? "hiacn" : "ac", h);
    d = o;
    return true;
  }
  if (name == "pairb") {
    if (!keep(p, d, {otc[0], otc[1], otc[4], otc[5], otc[6], "Number_Of_Placements"}, true)) return false;
    const Col mh = flag(*d.col(otc[4])), sr = flag(*d.col(otc[5])), sw = flag(*d.col(otc[6]));
    bool w = false;
    const Col np = asnum(*d.col("Number_Of_Placements"), w);
    Col tb = ncol(d.n), yb = ncol(d.n);
    for (std::size_t i = 0; i < d.n; ++i) {
      tb.d[i] = mh.d[i] + sr.d[i] + sw.d[i] >= 2;
      yb.d[i] = !std::isnan(np.d[i]) && np.d[i] >= 2;
    }
    d.set("tb", tb);
    d.set("yb", yb);
    return true;
  }
  if (name == "pairc") {
    if (!keep(p, d, {otc[2], otc[3], "NumberConsecutiveDays_Segregation"}, true)) return false;
    const Col &ra = *d.col(otc[2]), &rb = *d.col(otc[3]);
    bool w = false;
    Col y = asnum(*d.col("NumberConsecutiveDays_Segregation"), w);
    Col tc = ncol(d.n);
    for (std::size_t i = 0; i < d.n; ++i) {
      tc.d[i] = lab(ra, i) != lab(rb, i);
      if (std::isnan(y.d[i])) y.d[i] = 0;
    }
    const double cap = q7(y.d, 0.99);
    for (double &e : y.d) e = std::isnan(cap) ? cap : std::min(e, cap);
    d.set("tc", tc);
    d.set("yc", y);
    return true;
  }
  if (name == "gent") {
    const std::vector<std::string> a = split(arg, ',');
    std::vector<std::string> four(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(std::min<std::size_t>(4, a.size())));
    while (four.size() < 4) four.push_back("NA");
    if (!keep(p, d, four, false)) return false;
    bool bad = false;
    const double q1 = a.size() > 4 ? rnum(a[4], bad) : 0.5;
    const double q2 = a.size() > 4 ? (a.size() > 5 ? rnum(a[5], bad) : NA) : 0.667;
    std::vector<Col> v;
    double cut[4];
    for (int k = 0; k < 4; ++k) {
      bool w = false;
      v.push_back(asnum(*d.col(four[static_cast<std::size_t>(k)]), w));
      if (w) p.coerce();
      const double q = k < 2 ? q1 : q2;
      cut[k] = std::isnan(q) ? NA : q7(v.back().d, q);
    }
    // R's three-valued & and <=, >=; NA where an operand is NA
    auto cmp = [](double x, double c, bool le) { return std::isnan(x) || std::isnan(c) ? -1 : le ? x <= c : x >= c; };
    auto and3 = [](int x, int y) { return x == 0 || y == 0 ? 0 : x < 0 || y < 0 ? -1 : 1; };
    Col g;
    g.k = Col::FAC;
    g.lev = {"ineligible", "eligible", "gentrified"};
    g.d.assign(d.n, NA);
    for (std::size_t i = 0; i < d.n; ++i) {
      const int el = and3(cmp(v[0].d[i], cut[0], true), cmp(v[1].d[i], cut[1], true));
      const int gr = and3(cmp(v[2].d[i], cut[2], false), cmp(v[3].d[i], cut[3], false));
      const int eg = and3(el, gr);
      if (el >= 0 && eg >= 0) g.d[i] = 1 + el + eg;
    }
    d.set("gent", g);
    return true;
  }
  // pop: count / rate * per
  const std::vector<std::string> a = split(arg, ';');
  std::vector<std::string> two(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(std::min<std::size_t>(2, a.size())));
  while (two.size() < 2) two.push_back("NA");
  if (!keep(p, d, two, false)) return false;
  bool w1 = false, w2 = false, w3 = false;
  const Col cn = asnum(*d.col(two[0]), w1);
  if (w1) p.coerce();
  const Col rt = asnum(*d.col(two[1]), w2);
  if (w2) p.coerce();
  const double per = a.size() > 2 ? rnum(a[2], w3) : NA;
  if (w3) p.coerce();
  Col pop = ncol(d.n);
  for (std::size_t i = 0; i < d.n; ++i) pop.d[i] = cn.d[i] / rt.d[i] * per;
  d.set("pop", pop);
  return true;
}

// A JSON value as a column: numbers (with "NaN", "Inf", "-Inf"), true and
// false, or text when any element is text (numbers and true/false then as
// R writes them); null is missing
Col jcol(const Val &v) {
  const std::vector<Val> one1(1, v);
  const std::vector<Val> &e = v.t == Val::ARR ? v.a : one1;
  bool str = false, nums = false, lgl = false;
  for (const Val &x : e) {
    double t = 0;
    if (x.t == Val::NUM || snum(x, t)) nums = true;
    else if (x.t == Val::BOOL) lgl = true;
    else if (x.t != Val::NUL) str = true;
  }
  Col c;
  c.k = str ? Col::STR : nums ? Col::NUM : Col::LGL;
  for (const Val &x : e) {
    double t = 0;
    const bool special = snum(x, t);
    if (c.k == Col::STR) {
      c.sna.push_back(x.t == Val::NUL);
      c.s.push_back(x.t == Val::STR ? x.s : x.t == Val::NUM ? rfmt(x.n) : x.t == Val::BOOL ? (x.b ? "TRUE" : "FALSE") : "");
    } else {
      c.d.push_back(x.t == Val::NUL ? NA : x.t == Val::BOOL ? x.b : special ? t : x.n);
    }
  }
  (void)lgl;
  return c;
}

// A JSON object of columns of one length as a data frame
bool jframe(const Val &v, Tbl &f) {
  if (v.t != Val::OBJ) return false;
  f = Tbl();
  for (std::size_t j = 0; j < v.a.size(); ++j) {
    if (v.a[j].t == Val::OBJ) return false;
    const std::size_t len = v.a[j].t == Val::ARR ? v.a[j].a.size() : 1;
    if (j > 0 && len != f.n) return false;
    f.n = len;
    for (const Val &e : v.a[j].t == Val::ARR ? v.a[j].a : std::vector<Val>()) {
      if (e.t == Val::ARR || e.t == Val::OBJ) return false;
    }
    f.nm.push_back(v.k[j]);
    f.c.push_back(jcol(v.a[j]));
  }
  return true;
}

Val jv(const Col &c, std::size_t i) {
  if (c.miss(i)) {
    if (c.k == Col::NUM && !isna(c.d[i])) return vnum(NaN);
    return Val();
  }
  switch (c.k) {
    case Col::NUM: return vnum(c.d[i]);
    case Col::LGL: return vbool(c.d[i] != 0);
    default: return vstr(lab(c, i));
  }
}

Val jarr(const Col &c) {
  Val a = varr();
  for (std::size_t i = 0; i < c.size(); ++i) a.a.push_back(jv(c, i));
  return a;
}

// the rows of the preset table: one per preset and argument
typedef std::vector<std::array<std::string, 9>> Rows;
const char *const pcol[] = {"preset", "generic", "arg", "column", "recode", "filter", "units", "note", "source"};

}  // namespace pp

bool clprep(Ctx &c, Val &out) {
  using namespace pp;
  Pr p{c};
  Rows tab;
  for (const auto &r : ptsv) {
    std::array<std::string, 9> x;
    for (int j = 0; j < 9; ++j) x[static_cast<std::size_t>(j)] = r[j];
    tab.push_back(x);
  }
  const Val &pv = c["preset"];
  std::string pre;
  Tbl pt;
  if (jframe(pv, pt)) {
    // .ptab(): a preset of one's own
    if (!pt.col("preset") || !pt.col("generic") || !pt.col("arg") || pt.n == 0) {
      return c.fail("a preset table needs the columns preset, generic and arg, and a row.");
    }
    tab.assign(pt.n, std::array<std::string, 9>());
    for (std::size_t j = 0; j < 9; ++j) {
      const Col *x = pt.col(pcol[j]);
      for (std::size_t i = 0; x && i < pt.n; ++i) tab[i][j] = x->miss(i) ? "" : lab(*x, i);
    }
    for (const auto &r : tab) {
      if (r[0] != tab[0][0]) return c.fail("a preset table must hold one preset.");
    }
    pre = tab[0][0];
  } else if (pv.t == Val::NUL) {
    if (c.has("data")) return c.fail("`preset` must name a preset.");
    // .plist(): one row per preset
    std::vector<std::string> nm, gen, args, src;
    for (const auto &r : tab) {
      const std::size_t k = static_cast<std::size_t>(std::find(nm.begin(), nm.end(), r[0]) - nm.begin());
      if (k == nm.size()) {
        nm.push_back(r[0]);
        gen.push_back(r[1]);
        src.push_back(r[8]);
        args.push_back(r[2]);
      } else if (std::find_if(tab.begin(), tab.end(), [&r](const std::array<std::string, 9> &q) {
                   return &q < &r && q[0] == r[0] && q[2] == r[2];
                 }) == tab.end()) {
        args[k] += " " + r[2];
      }
    }
    put(out, "preset", vstrs(nm));
    put(out, "generic", vstrs(gen));
    put(out, "args", vstrs(args));
    put(out, "source", vstrs(src));
    return true;
  } else {
    if (pv.t != Val::STR) return c.fail("`preset` must be one name.");
    pre = pv.s;
  }
  // .prows()
  Rows r;
  for (const auto &x : tab) {
    if (x[0] == pre) r.push_back(x);
  }
  if (r.empty()) return c.fail("no preset named %s; clprep() lists them.", pre.c_str());
  Tbl d;
  if (!jframe(c["data"], d)) return c.fail("`data` must be a data frame.");
  Rows a;
  for (const auto &x : r) {
    if (x[2] == "recipe") {
      if (!recipe(p, x[3], x[4], d)) return false;
    } else {
      a.push_back(x);
    }
  }
  // the columns read, and the rows with none of them missing
  std::vector<std::string> used, miss;
  auto use = [&used](const std::string &k) {
    if (std::find(used.begin(), used.end(), k) == used.end()) used.push_back(k);
  };
  for (const auto &x : a) {
    if (!x[3].empty()) use(x[3]);
  }
  std::vector<std::pair<std::string, std::string>> flt;
  for (const auto &x : a) {
    if (x[5].empty()) continue;
    for (const std::string &f : split(x[5], ';')) {
      const std::vector<std::string> kv = split(f, '=');
      flt.emplace_back(kv.empty() ? std::string() : kv[0], kv.size() > 1 ? kv[1] : std::string());
      use(flt.back().first);
    }
  }
  for (const std::string &k : used) {
    if (!d.col(k)) miss.push_back(k);
  }
  if (!miss.empty()) {
    return c.fail("`data` has no column %s (preset %s).", tostr(miss).c_str(), pre.c_str());
  }
  const std::size_t n = d.n;
  std::vector<char> ok(n, 1);
  for (const std::string &k : used) {
    const Col *x = d.col(k);
    for (std::size_t i = 0; i < n; ++i) ok[i] = ok[i] && !x->miss(i);
  }
  for (const auto &f : flt) {
    const Col *x = d.col(f.first);
    for (std::size_t i = 0; i < n; ++i) ok[i] = ok[i] && trim(lab(*x, i)) == f.second;
  }
  std::vector<V> vals(a.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (!prc(p, a[i][3].empty() ? nullptr : d.col(a[i][3]), a[i][4], n, vals[i])) return false;
  }
  for (const V &v : vals) {
    if (v.len() != n || (v.mat && v.m.empty())) continue;
    const Col &x = v.mat ? v.m[0] : v.c;
    for (std::size_t i = 0; i < n; ++i) ok[i] = ok[i] && !x.miss(i);
  }
  const std::vector<std::size_t> kept = which(ok);
  // .pargs(): each argument's values on the rows kept
  Val args;
  args.t = Val::OBJ;
  std::vector<std::string> seen;
  for (std::size_t i = 0; i < a.size(); ++i) {
    const std::string &an = a[i][2];
    if (std::find(seen.begin(), seen.end(), an) != seen.end()) continue;
    seen.push_back(an);
    std::vector<std::size_t> j;
    for (std::size_t k = i; k < a.size(); ++k) {
      if (a[k][2] == an) j.push_back(k);
    }
    auto sub = [&](const Col &x) { return x.size() == n ? rows(x, kept) : x; };
    const V &v0 = vals[j[0]];
    if (an == "adj" || j.size() > 1 || v0.mat) {
      std::vector<std::string> nm;
      std::vector<Col> cs;
      for (std::size_t k : j) {
        const V &v = vals[k];
        if (v.mat) {
          for (std::size_t m = 0; m < v.m.size(); ++m) {
            nm.push_back(a[k][3] + v.mn[m]);
            cs.push_back(rows(v.m[m], kept));
          }
        } else if (!v.nul) {
          nm.push_back(a[k][3]);
          cs.push_back(sub(v.c));
        }
      }
      std::size_t nr = 0;
      for (const Col &x : cs) nr = std::max(nr, x.size());
      Val df;
      df.t = Val::OBJ;
      for (std::size_t m = 0; m < cs.size(); ++m) {
        Val col = jarr(cs[m]);
        if (cs[m].size() == 1 && nr > 1) col.a = std::vector<Val>(nr, col.a[0]);
        put(df, nm[m].c_str(), col);
      }
      put(args, an.c_str(), df);
    } else if (!v0.nul) {
      put(args, an.c_str(), unbox(jarr(sub(v0.c))));
    }
  }
  put(out, "args", args);
  put(out, "preset", vstr(pre));
  Val gen, src, units, note = varr();
  if (!a.empty()) {
    gen = vstr(a[0][1]);
    src = vstr(a[0][8]);
  }
  for (const auto &x : a) {
    if (!x[6].empty() && units.t == Val::NUL) units = vstr(x[6]);
    if (!x[7].empty()) note.a.push_back(vstr(x[7]));
  }
  put(out, "generic", gen);
  put(out, "source", src);
  put(out, "units", units);
  put(out, "note", unbox(note));
  put(out, "drop", vnum(static_cast<double>(n - kept.size())));
  return true;
}
