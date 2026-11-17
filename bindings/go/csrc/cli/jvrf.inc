// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// vrfpol() and polrpttxt() of R/vrfp.R for clab_jrun().
//
// polrpttxt(): the text of a vrfpol() report given as JSON (as R's report
// is written), line for line, with R's checks and words. Numbers are
// written as R's sprintf() and format() write them (rfmt7()). A line R
// builds from a part that is missing is left out, as sprintf() of NULL
// leaves it out.
//
// vrfpol(): the report R gives, member for member, from the stages of
// japol.inc (crfno2(), crfpm25(), attrfrac(), mortdisp(), polbrdn(),
// expconcidx()) run on the exposure R resolves. The demo data are drawn
// from R's generator after set.seed(20260417) (igh::Rng), as R draws them.
// The exposure file is read as .vpcsv() and read.csv() read it: plain or
// gzip text in UTF-8, with or without a byte-order mark, lines ending in
// LF, CRLF or CR, fields split at commas with double quotes; a column's
// type is the one type.convert() gives it. What R's reader decodes and
// this one does not is refused in words: a bzip2 or xz file, UTF-16 or
// UTF-32 text, and text that is not UTF-8 (R reads it as Latin-1). JSON
// null is an argument left out, so NA cannot be given for ref, expomean,
// expoprev, brate or pop as it can in R. R's notes on standard error
// (rows dropped, units converted) are not repeated.

namespace vrf {

// R's text of a number in sprintf(): NA, NaN, Inf, -Inf, else f.
inline std::string pf(const char *f, double x) {
  if (isna(x)) return "NA";
  if (std::isnan(x)) return "NaN";
  if (std::isinf(x)) return x > 0 ? "Inf" : "-Inf";
  char b[400];
  std::snprintf(b, sizeof b, f, x);
  return b;
}

// format(x) with R's 7 significant digits (dig of them with digits =
// dig): the fewest digits up to 7 that show x, in fixed notation unless
// scientific is narrower; with big, fixed always and commas between
// thousands (format(x, big.mark = ",", scientific = FALSE)).
inline std::string rfmt7(double x, bool big, int dig = 7) {
  if (isna(x)) return "NA";
  if (std::isnan(x)) return "NaN";
  if (std::isinf(x)) return x > 0 ? "Inf" : "-Inf";
  char b[400];
  std::snprintf(b, sizeof b, "%.*e", dig - 1, x);
  const std::string s = b;
  const std::size_t e = s.find('e');
  const int kp = std::atoi(s.c_str() + e + 1);
  std::string dg;
  for (std::size_t i = 0; i < e; ++i) {
    if (std::isdigit(static_cast<unsigned char>(s[i]))) dg += s[i];
  }
  int nsig = dig;
  while (nsig > 1 && dg[static_cast<std::size_t>(nsig) - 1] == '0') --nsig;
  const int neg = x < 0 ? 1 : 0;
  const int rgt = std::max(0, nsig - kp - 1);
  const int wf = neg + (kp >= 0 ? kp + 1 : 1) + (rgt > 0 ? rgt + 1 : 0);
  const int we = neg + (nsig > 1 ? nsig + 1 : 1) + (kp >= 100 || kp <= -100 ? 5 : 4);
  if (!big && wf > we) {
    std::snprintf(b, sizeof b, "%.*e", nsig - 1, x);
    return b;
  }
  std::snprintf(b, sizeof b, "%.*f", rgt, x);
  std::string o = b;
  if (!big) return o;
  const std::size_t dot = o.find('.') == std::string::npos ? o.size() : o.find('.');
  for (std::size_t k = dot; k > static_cast<std::size_t>(neg) + 3; k -= 3) o.insert(k - 3, ",");
  return o;
}

// A member of an object (null when absent).
inline const Val &mem(const Val &o, const char *k) {
  const Val *v = o.t == Val::OBJ ? ax::fc(o, k) : nullptr;
  return v == nullptr ? nul : *v;
}
// A number of a member, or false (the line is left out) when it is null.
inline bool num(const Val &o, const char *k, double &x) {
  const Val &v = mem(o, k);
  const Val &e = v.t == Val::ARR && v.a.size() == 1 ? v.a[0] : v;
  if (e.t != Val::NUM && e.t != Val::STR) return false;
  x = e.t == Val::NUM ? e.n : (e.s == "NaN" ? NaN : e.s == "Inf" ? HUGE_VAL : e.s == "-Inf" ? -HUGE_VAL : NA);
  return true;
}
// The texts of a member, one per element (none when null).
inline std::vector<std::string> txt(const Val &o, const char *k) {
  const Val &v = mem(o, k);
  std::vector<std::string> r;
  for (const Val &e : v.t == Val::ARR ? v.a : std::vector<Val>{v}) {
    if (e.t == Val::STR) r.push_back(e.s);
    else if (e.t == Val::NUM) r.push_back(rfmt(e.n));
    else if (e.t == Val::BOOL) r.push_back(e.b ? "TRUE" : "FALSE");
    else if (e.t == Val::NUL && v.t == Val::ARR) r.push_back("NA");
  }
  return r;
}
inline std::string up(std::string s) {
  for (char &ch : s) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
  return s;
}

// .vphead()
inline void head(const Val &r, std::vector<std::string> &l) {
  const std::string bar(66, '=');
  l.push_back(bar);
  const auto cm = txt(r, "command"), pl = txt(r, "pollutant"), oc = txt(r, "outcome");
  if (!cm.empty() && !pl.empty() && !oc.empty()) l.push_back("  " + cm[0] + " -- " + up(pl[0]) + " -> " + oc[0]);
  for (const std::string &s : txt(r, "region")) l.push_back("  region: " + s);
  for (const std::string &s : txt(r, "years")) l.push_back("  years:  " + s);
  for (const std::string &s : txt(r, "dsou")) l.push_back("  data:   " + s);
  l.push_back(bar);
  l.push_back("");
  l.push_back("Inputs");
  const Val &in = mem(r, "inputs");
  double x = 0;
  if (num(in, "expomean", x)) l.push_back("  exposure mean:      " + pf("%.3f", x));
  if (num(in, "expoprev", x)) l.push_back("  exposure prevalence:" + pf("%.3f", x));
  if (num(in, "baserateper100k", x)) l.push_back("  baseline rate/100k: " + pf("%.2f", x));
  if (num(in, "population", x)) l.push_back("  population:         " + rfmt7(x, true));
  if (num(in, "refconc", x)) l.push_back("  reference conc:     " + rfmt7(x, false));
  if (num(in, "drorows", x) && x > 0) l.push_back("  rows dropped:       " + pf("%.0f", x) + " (no exposure value)");
  l.push_back("");
  l.push_back("Assumption log");
}

// .vppipe()
inline void pipe(const Val &p, std::vector<std::string> &l) {
  const Val &ci = mem(p, "pafci"), &crf = mem(p, "crf"), &dp = mem(p, "displaced"), &bd = mem(p, "burden");
  std::string cis;
  double lo = 0, hi = 0, x = 0;
  if (ci.t != Val::NUL) {
    if (!num(ci, "lower", lo) || !num(ci, "upper", hi)) lo = hi = NA;
    cis = "  (95% CI " + pf("%.4f", lo) + " to " + pf("%.4f", hi) + ", from the risk's limits)";
  }
  l.push_back("");
  l.push_back("Concentration-response");
  if (num(crf, "rr", x)) l.push_back("  RR:       " + pf("%.4f", x));
  for (const std::string &s : txt(crf, "citation")) l.push_back("  source:   " + s);
  l.push_back("");
  if (num(p, "paf", x)) l.push_back("Attributable fraction (PAF): " + pf("%.4f", x) + cis);
  l.push_back("");
  const auto un = txt(dp, "unit");
  l.push_back(un.size() == 1 && un[0] == "deaths" ? "Mortality displaced" : "Cases displaced");
  if (num(dp, "deadisp", x)) l.push_back("  expected avoided " + (un.empty() ? "deaths" : un[0]) + ": " + pf("%.1f", x));
  l.push_back("");
  l.push_back("Burden of pollution");
  const auto bu = txt(mem(bd, "extra"), "unit");
  if (num(bd, "attrcases", x)) l.push_back("  attributable " + (bu.empty() ? "deaths" : bu[0]) + ":   " + pf("%.1f", x));
  const Val &eq = mem(p, "equity");
  if (eq.t != Val::NUL) {
    l.push_back("");
    l.push_back("Equity analysis");
    if (num(eq, "concidx", x)) l.push_back("  concentration index: " + pf("%.4f", x));
  }
}

}  // namespace vrf

// polrpttxt()
bool polrpttxt(Ctx &c, Val &out) {
  using namespace vrf;
  const Val &r = c["report"];
  const Val &st = mem(r, "status");
  if (r.t != Val::OBJ || st.t != Val::STR) return ckstop(c, "`report` must be a report from vrfpol()");
  if (st.s == "error") {
    const auto e = txt(r, "error");
    out = vstr("ERROR: " + (e.empty() ? std::string() : e[0]) + "\n");
    return true;
  }
  const Val &in = mem(r, "inputs"), &as = mem(r, "assumptions"), &pp = mem(r, "pipeline");
  bool ok = in.t == Val::OBJ && (as.t == Val::ARR || as.t == Val::OBJ) && (st.s != "ok" || pp.t == Val::OBJ);
  for (const Val &a : as.a) {
    ok = ok && a.t == Val::OBJ && (mem(a, "ok").t == Val::BOOL || mem(a, "ok").t == Val::NUL) && ax::fc(a, "ok");
  }
  if (!ok) {
    return ckstop(c, "`report` is not a complete vrfpol() report (its inputs, assumptions or pipeline are missing)");
  }
  std::vector<std::string> l;
  head(r, l);
  for (const Val &a : as.a) {
    const Val &k = mem(a, "ok");
    const auto an = txt(a, "assumption"), nt = txt(a, "note");
    if (an.empty() || nt.empty()) continue;
    l.push_back(std::string("  [") + (k.t == Val::BOOL ? (k.b ? "PASS" : "FAIL") : "NA") + "] " + an[0] + " -- " + nt[0]);
  }
  std::string s;
  if (st.s != "assumption_failure") pipe(pp, l);
  for (std::size_t i = 0; i < l.size(); ++i) s += (i ? "\n" : "") + l[i];
  s += st.s == "assumption_failure" ? "\n\nSTATUS: assumption_failure (pipeline skipped)\n" : "\n\nSTATUS: ok\n";
  out = vstr(s);
  return true;
}

// vrfpol()

bool u8ok(const std::string &s);

namespace vrf {

// One value of the type t: itself, or an array of one.
inline const Val *one(const Val &o, Val::T t) {
  const Val *s = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  return s->t == t ? s : nullptr;
}

// trimws() and tolower() (ASCII letters) of R.
inline std::string trim(const std::string &s) {
  const std::size_t a = s.find_first_not_of(" \t\r\n");
  return a == std::string::npos ? std::string() : s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}
inline std::string lower(std::string s) {
  for (char &ch : s) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  return s;
}

// Inflating deflate data (RFC 1951) as puff.c of zlib does: the bits of
// s from p, the bytes out appended to o (from st on for this member).
struct Zin {
  const unsigned char *s;
  std::size_t n, p;
  std::string &o;
  std::size_t lim, st;
  unsigned long bb = 0;
  int bc = 0;
  bool eof = false, big = false;
  Zin(const unsigned char *s0, std::size_t n0, std::size_t p0, std::string &o0, std::size_t l0)
      : s(s0), n(n0), p(p0), o(o0), lim(l0), st(o0.size()) {}
  int bits(int need) {
    unsigned long v = bb;
    while (bc < need) {
      if (p >= n) {
        eof = true;
        return 0;
      }
      v |= static_cast<unsigned long>(s[p++]) << bc;
      bc += 8;
    }
    bb = v >> need;
    bc -= need;
    return static_cast<int>(v & ((1UL << need) - 1));
  }
};

// A canonical Huffman code: the codes of each length, the symbols in order.
struct Huf {
  int cnt[16];
  int sym[320];
};

// The code of the lengths len (n of them): 0 when complete, above 0 when
// incomplete, below 0 when over-subscribed.
inline int mkhuf(Huf &h, const int *len, int n) {
  for (int &k : h.cnt) k = 0;
  for (int i = 0; i < n; ++i) h.cnt[len[i]]++;
  if (h.cnt[0] == n) return 0;
  int left = 1;
  for (int l = 1; l < 16; ++l) {
    left <<= 1;
    left -= h.cnt[l];
    if (left < 0) return left;
  }
  int offs[16];
  offs[1] = 0;
  for (int l = 1; l < 15; ++l) offs[l + 1] = offs[l] + h.cnt[l];
  for (int i = 0; i < n; ++i) {
    if (len[i] != 0) h.sym[offs[len[i]]++] = i;
  }
  return left;
}

inline int hdec(Zin &z, const Huf &h) {
  int code = 0, first = 0, index = 0;
  for (int l = 1; l < 16; ++l) {
    code |= z.bits(1);
    if (z.eof) return -1;
    const int k = h.cnt[l];
    if (code - k < first) return h.sym[index + (code - first)];
    index += k;
    first += k;
    first <<= 1;
    code <<= 1;
  }
  return -1;
}

inline bool hcodes(Zin &z, const Huf &lc, const Huf &dc) {
  static const int lb[29] = {3,  4,  5,  6,  7,  8,  9,  10, 11,  13,  15,  17,  19,  23, 27,
                             31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258};
  static const int le[29] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0};
  static const int db[30] = {1,   2,   3,   4,   5,   7,    9,    13,   17,   25,   33,   49,    65,    97,    129,
                             193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577};
  static const int de[30] = {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13};
  for (;;) {
    int s = hdec(z, lc);
    if (s < 0) return false;
    if (s == 256) return true;
    if (s < 256) {
      z.o.push_back(static_cast<char>(s));
    } else {
      s -= 257;
      if (s >= 29) return false;
      const int len = lb[s] + z.bits(le[s]);
      const int d = hdec(z, dc);
      if (d < 0 || d >= 30) return false;
      const std::size_t dist = static_cast<std::size_t>(db[d] + z.bits(de[d]));
      if (z.eof || dist > z.o.size() - z.st) return false;
      for (int k = 0; k < len; ++k) z.o.push_back(z.o[z.o.size() - dist]);
    }
    if (z.o.size() > z.lim) {
      z.big = true;
      return false;
    }
  }
}

inline bool stored(Zin &z) {
  z.bb = 0;
  z.bc = 0;
  if (z.n - z.p < 4) {
    z.eof = true;
    return false;
  }
  const unsigned len = static_cast<unsigned>(z.s[z.p]) | static_cast<unsigned>(z.s[z.p + 1]) << 8;
  const unsigned nl = static_cast<unsigned>(z.s[z.p + 2]) | static_cast<unsigned>(z.s[z.p + 3]) << 8;
  z.p += 4;
  if (len != (~nl & 0xffffu)) return false;
  if (z.n - z.p < len) {
    z.eof = true;
    return false;
  }
  z.o.append(reinterpret_cast<const char *>(z.s + z.p), len);
  z.p += len;
  z.big = z.o.size() > z.lim;
  return !z.big;
}

inline bool fixd(Zin &z) {
  int l[288];
  Huf lc, dc;
  for (int i = 0; i < 288; ++i) l[i] = i < 144 ? 8 : i < 256 ? 9 : i < 280 ? 7 : 8;
  mkhuf(lc, l, 288);
  for (int i = 0; i < 30; ++i) l[i] = 5;
  mkhuf(dc, l, 30);
  return hcodes(z, lc, dc);
}

inline bool dyn(Zin &z) {
  static const int ord[19] = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};
  int ln[320];
  const int nlen = z.bits(5) + 257, ndist = z.bits(5) + 1, ncode = z.bits(4) + 4;
  if (z.eof || nlen > 286 || ndist > 30) return false;
  int i = 0;
  for (; i < ncode; ++i) ln[ord[i]] = z.bits(3);
  for (; i < 19; ++i) ln[ord[i]] = 0;
  if (z.eof) return false;
  Huf lc, dc;
  if (mkhuf(lc, ln, 19) != 0) return false;
  i = 0;
  while (i < nlen + ndist) {
    const int s = hdec(z, lc);
    if (s < 0) return false;
    if (s < 16) {
      ln[i++] = s;
      continue;
    }
    int len = 0, rep = 0;
    if (s == 16) {
      if (i == 0) return false;
      len = ln[i - 1];
      rep = 3 + z.bits(2);
    } else if (s == 17) {
      rep = 3 + z.bits(3);
    } else {
      rep = 11 + z.bits(7);
    }
    if (z.eof || i + rep > nlen + ndist) return false;
    while (rep-- > 0) ln[i++] = len;
  }
  if (ln[256] == 0) return false;
  int e = mkhuf(lc, ln, nlen);
  if (e < 0 || (e > 0 && nlen - lc.cnt[0] != 1)) return false;
  e = mkhuf(dc, ln + nlen, ndist);
  if (e < 0 || (e > 0 && ndist - dc.cnt[0] != 1)) return false;
  return hcodes(z, lc, dc);
}

inline std::uint32_t le32(const unsigned char *t) {
  return static_cast<std::uint32_t>(t[0]) | static_cast<std::uint32_t>(t[1]) << 8 |
         static_cast<std::uint32_t>(t[2]) << 16 | static_cast<std::uint32_t>(t[3]) << 24;
}

// The data of every member of a gzip file (RFC 1952) into o, each checked
// against the CRC-32 and length that end it: "" or .vpunzip()'s refusal.
// A file cut short, or bytes after the last member, fail the check as
// they fail it in R (.vpgzend()).
inline std::string gunzip(const std::string &b, std::size_t lim, std::string &o) {
  const std::string bad = "is truncated or damaged: its data do not match the length and CRC-32 at its end";
  const unsigned char *s = reinterpret_cast<const unsigned char *>(b.data());
  const std::size_t n = b.size();
  std::size_t p = 0;
  while (p < n) {
    if (n - p < 10 || s[p] != 0x1f || s[p + 1] != 0x8b || s[p + 2] != 8) return bad;
    const int fl = s[p + 3];
    std::size_t q = p + 10;
    if (fl & 4) {
      if (n - q < 2) return bad;
      q += 2 + (static_cast<std::size_t>(s[q]) | static_cast<std::size_t>(s[q + 1]) << 8);
    }
    for (int f : {8, 16}) {
      if (!(fl & f)) continue;
      while (q < n && s[q] != 0) ++q;
      ++q;
    }
    if (fl & 2) q += 2;
    if (q > n) return bad;
    Zin z(s, n, q, o, lim);
    int last = 0;
    do {
      last = z.bits(1);
      const int t = z.bits(2);
      if (z.eof) return bad;
      const bool ok = t == 0 ? stored(z) : t == 1 ? fixd(z) : t == 2 ? dyn(z) : false;
      if (z.big) return "decompresses to more than the size limit";
      if (!ok) return bad;
    } while (!last);
    p = z.p;
    if (n - p < 8) return bad;
    const std::size_t m = o.size() - z.st;
    if (clab_crc32(reinterpret_cast<const unsigned char *>(o.data()) + z.st, m, 0) != le32(s + p) ||
        static_cast<std::uint32_t>(m) != le32(s + p + 4)) {
      return bad;
    }
    p += 8;
  }
  return "";
}

// -1: no such path, 1: a folder, 0: anything else (file.exists(),
// dir.exists()).
inline int ftype(const std::string &f) {
#ifdef _WIN32
  struct _stat64 st;
  if (_stat64(f.c_str(), &st) != 0) return -1;
  return (st.st_mode & _S_IFMT) == _S_IFDIR ? 1 : 0;
#else
  struct stat st;
  if (stat(f.c_str(), &st) != 0) return -1;
  return S_ISDIR(st.st_mode) ? 1 : 0;
#endif
}

inline bool slurp(const std::string &f, std::string &b) {
  std::FILE *h = std::fopen(f.c_str(), "rb");
  if (h == nullptr) return false;
  char t[65536];
  std::size_t k;
  while ((k = std::fread(t, 1, sizeof t, h)) > 0) b.append(t, k);
  const bool ok = std::ferror(h) == 0;
  std::fclose(h);
  return ok;
}

// The table read.csv() reads: the header, then each column's cells.
struct Tab {
  std::vector<std::string> nm;
  std::vector<std::vector<std::string>> v;
};

inline bool blank(const char *p) {
  for (; *p != '\0'; ++p) {
    if (!std::isspace(static_cast<unsigned char>(*p))) return false;
  }
  return true;
}
inline bool blank(const std::string &s) { return blank(s.c_str()); }

// R_strtod() with nothing but white space after it (type.convert(),
// as.numeric()): "NA" is NA.
inline bool rnum(const std::string &s, double &x) {
  const char *p = s.c_str();
  while (std::isspace(static_cast<unsigned char>(*p))) ++p;
  if (std::strncmp(p, "NA", 2) == 0) {
    x = NA;
    return blank(p + 2);
  }
  char *q = nullptr;
  x = std::strtod(p, &q);
  return q != p && blank(q);
}

inline bool islg(const std::string &s, bool &b) {
  static const char *const w[8] = {"T", "True", "TRUE", "true", "F", "False", "FALSE", "false"};
  for (int i = 0; i < 8; ++i) {
    if (s == w[i]) {
      b = i < 4;
      return true;
    }
  }
  return false;
}

// type.convert() of a column with na.strings "NA": 0 logical, 1 numeric,
// 2 character. Blank cells are missing in the first two.
inline int ctype(const std::vector<std::string> &v) {
  bool lg = true, nm = true;
  for (const std::string &s : v) {
    if (s == "NA" || blank(s)) continue;
    bool b = false;
    double x = 0;
    lg = lg && islg(s, b);
    nm = nm && rnum(s, x);
  }
  return lg ? 0 : nm ? 1 : 2;
}

// as.numeric() of a column of type ty.
inline std::vector<double> asnum(const std::vector<std::string> &v, int ty) {
  std::vector<double> r;
  for (const std::string &s : v) {
    bool b = false;
    double x = NA;
    if (s == "NA") x = NA;
    else if (ty == 0) x = islg(s, b) ? (b ? 1.0 : 0.0) : NA;
    else if (!rnum(s, x)) x = NA;
    r.push_back(x);
  }
  return r;
}

// .vpshape() and the read of .vpparse(): the records of the text, split
// as scan() splits them (sep ",", quote "\"", a doubled quote inside
// quotes is one, LF, CRLF and CR end a line, empty lines are skipped),
// each the same number of fields as the header.
inline std::string rdtab(const std::string &t, Tab &tb) {
  struct Rec {
    std::vector<std::string> f;
    int end;
  };
  std::vector<Rec> rs;
  bool eofq = false;
  int nl = 0;  // the lines count.fields() counts
  std::size_t i = 0;
  const std::size_t n = t.size();
  while (i < n && !eofq) {
    Rec r;
    std::string cur;
    bool q = false, any = false;
    int lines = 1;
    for (;;) {
      if (i >= n) {
        eofq = q;
        break;
      }
      char ch = t[i++];
      if (ch == '\r') {
        if (i < n && t[i] == '\n') ++i;
        ch = '\n';
      }
      if (q) {
        if (ch == '"' && i < n && t[i] == '"') {
          cur += '"';
          ++i;
        } else if (ch == '"') {
          q = false;
        } else {
          lines += ch == '\n' ? 1 : 0;
          cur += ch;
        }
        continue;
      }
      if (ch == '\n') break;
      any = true;
      if (ch == '"') {
        q = true;
      } else if (ch == ',') {
        r.f.push_back(cur);
        cur.clear();
      } else {
        cur += ch;
      }
    }
    // a record a quote leaves open ends at the end of the text, as
    // count.fields() ends it: its fields so far, on the last line
    if (!any) continue;
    r.f.push_back(cur);
    nl += lines;
    r.end = nl;
    rs.push_back(r);
  }
  const std::string open =
      "could not be read cleanly (EOF within quoted string): no rows are used from a file the parser stopped "
      "part-way through";
  // .vpshape() reads the header and the first row (readtablehead()) before
  // anything else: a quote left open there stops that read
  if (eofq && rs.size() <= 2) return open;
  if (rs.empty()) return "";
  const std::vector<std::string> &hd = rs[0].f;
  std::vector<std::string> key, dup;
  for (const std::string &h : hd) {
    std::string k;
    for (char ch : trim(h)) k += static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    if (std::find(key.begin(), key.end(), k) != key.end() && std::find(dup.begin(), dup.end(), h) == dup.end()) {
      dup.push_back(h);
    }
    key.push_back(k);
  }
  if (!dup.empty()) {
    std::string l;
    for (std::size_t k = 0; k < dup.size(); ++k) l += (k ? ", '" : "'") + dup[k] + "'";
    return "has the column name " + l + " more than once: say which column is meant";
  }
  for (const Rec &r : rs) {
    if (r.f.size() != hd.size()) {
      return "line " + std::to_string(r.end) + " has " + std::to_string(r.f.size()) + " fields where the header has " +
             std::to_string(hd.size()) + ": no rows are used from a file whose rows do not line up";
    }
  }
  if (eofq) return open;
  tb.nm = hd;
  tb.v.assign(hd.size(), std::vector<std::string>());
  for (std::size_t k = 1; k < rs.size(); ++k) {
    for (std::size_t j = 0; j < hd.size(); ++j) tb.v[j].push_back(rs[k].f[j]);
  }
  return "";
}

// .vpload() with .vpcsv(): the table of the file, or the refusal.
inline std::string vpload(const std::string &f, Tab &tb) {
  const int k = ftype(f);
  if (k < 0) return "exposure CSV not found: " + f;
  if (k > 0) return "exposure CSV is a directory, not a file: " + f;
  const std::string empty =
      f + " is empty (expected a column exposure in ug/m3, or a NAPS pull's value + unit)";
  std::string b;
  if (!slurp(f, b)) return f + ": cannot open the connection";
  if (b.empty()) return empty;
  static const unsigned char xz[6] = {0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00};
  const unsigned char *u = reinterpret_cast<const unsigned char *>(b.data());
  const std::size_t n = b.size();
  if ((n >= 3 && b.compare(0, 3, "BZh") == 0) || (n >= 6 && std::memcmp(u, xz, 6) == 0)) {
    return f + ": is " + (u[0] == 'B' ? "bzip2" : "xz") +
           "-compressed, which only R decompresses: decompress it first, or gzip it";
  }
  std::string d;
  if (n >= 2 && u[0] == 0x1f && u[1] == 0x8b) {
    const double lim = std::min(std::min(1073741824.0, std::max(67108864.0, 200.0 * static_cast<double>(n))),
                                2147483646.0);
    const std::string e = gunzip(b, static_cast<std::size_t>(lim), d);
    if (!e.empty()) return f + ": " + e;
  } else {
    d.swap(b);
  }
  if (d.empty()) return empty;
  int bom = 0;
  const int r = clab_txtenc(reinterpret_cast<const unsigned char *>(d.data()), static_cast<std::int64_t>(d.size()), &bom);
  if (r < 0) return f + ": contains NUL bytes and is not UTF-16 text: not a text CSV";
  if (r > 0) {
    static const char *const en[4] = {"UTF-32LE", "UTF-32BE", "UTF-16LE", "UTF-16BE"};
    return f + ": is " + en[r - 1] + " text, which only R decodes: save it as UTF-8";
  }
  if (d.compare(0, 3, "\xef\xbb\xbf") == 0) d.erase(0, 3);
  if (!u8ok(d)) return f + ": is not UTF-8 text, which only R decodes (as Latin-1): save it as UTF-8";
  if (blank(d)) return empty;
  const std::string e = rdtab(d, tb);
  if (!e.empty()) return f + ": " + e;
  if (tb.v.empty() || tb.v[0].empty()) return empty;
  return "";
}

// The exposure of a run: its mean and share exposed, the data behind
// them (exposure, and income for the equity stage) and where they came
// from.
struct Src {
  double em = NA, ep = NA, drop = 0;
  std::string ds;
  bool eq = false;
  std::vector<double> ex, inc;
};

inline int coln(const Tab &t, const char *k) {
  for (std::size_t j = 0; j < t.nm.size(); ++j) {
    if (t.nm[j] == k) return static_cast<int>(j);
  }
  return -1;
}

// .vpnaps(): a NAPS pull's values in ug/m3 by their units, or the refusal.
inline std::string naps(const Tab &t, int iv, int iu, const std::string &poll, std::vector<double> &ex) {
  const int ty = ctype(t.v[iv]), tu = ctype(t.v[iu]);
  ex = asnum(t.v[iv], ty);
  static const char *const ug[4] = {"ug/m3", "\xc2\xb5g/m3", "\xc2\xb5g/m\xc2\xb3", "ug/m\xc2\xb3"};
  std::vector<std::string> bad;
  bool any = false;
  for (std::size_t i = 0; i < ex.size(); ++i) {
    const std::string &c = t.v[iu][i];
    // as.character() of the unit; NA stays NA
    bool b = false, na0 = c == "NA" || (tu != 2 && blank(c));
    double x = 0;
    std::string w = c;
    if (!na0 && tu == 0) w = islg(c, b) && b ? "TRUE" : "FALSE";
    if (!na0 && tu == 1) {
      rnum(c, x);
      na0 = isna(x);
      w = rfmt(x);
    }
    w = lower(trim(w));
    bool isug = false;
    for (const char *g : ug) isug = isug || (!na0 && w == g);
    const bool ppb = !na0 && (w == "ppb" || w == "ppbv");
    if (!std::isnan(ex[i]) && !isug && !(ppb && poll == "no2")) {
      any = true;
      if (!na0 && std::find(bad.begin(), bad.end(), w) == bad.end()) bad.push_back(w);
    }
    // ug/m3 per ppb of NO2 at 25 C and 1 atm (WHO 2021)
    if (ppb) ex[i] *= 1.88;
  }
  if (!any) return "";
  // sort() of the units: byte order here, the locale's order in R
  std::sort(bad.begin(), bad.end());
  std::string l;
  for (std::size_t k = 0; k < bad.size(); ++k) l += (k ? ", " : "") + bad[k];
  return "exposure unit " + l + " is not ug/m3 for " + poll;
}

// .vpread() with .vptext() and .vprows(): the exposures (and incomes) of
// the file into s, or the refusal.
inline std::string vpread(const std::string &f, const std::string &poll, Src &s) {
  Tab t;
  std::string e = vpload(f, t);
  if (!e.empty()) return e;
  const int ie = coln(t, "exposure"), iv = coln(t, "value");
  std::vector<double> ex;
  const std::vector<std::string> *raw = nullptr;
  if (ie < 0 && iv >= 0) {
    const int iu = coln(t, "unit");
    if (iu < 0) return "a 'value' column needs a 'unit' column (NAPS pulls carry one; NO2 is in ppb)";
    e = naps(t, iv, iu, poll, ex);
    if (!e.empty()) return e;
    if (ctype(t.v[iv]) == 2) raw = &t.v[iv];
  } else if (ie < 0) {
    return "CSV missing 'exposure' column (ug/m3); a NAPS pull's 'value' column also works.";
  } else {
    const int ty = ctype(t.v[ie]);
    ex = asnum(t.v[ie], ty);
    if (ty == 2) raw = &t.v[ie];
  }
  if (raw != nullptr) {
    static const char *const tk[8] = {"", "na", "n/a", "nan", "null", "-", ".", "missing"};
    std::size_t k = 0;
    std::vector<std::string> ux;
    for (std::size_t i = 0; i < ex.size(); ++i) {
      const std::string &c = (*raw)[i];
      if (!std::isnan(ex[i]) || c == "NA") continue;
      const std::string w = trim(c), l = lower(w);
      bool tok = false;
      for (const char *x : tk) tok = tok || l == x;
      if (tok) continue;
      ++k;
      if (ux.size() < 3 && std::find(ux.begin(), ux.end(), w) == ux.end()) ux.push_back(w);
    }
    if (k > 0) {
      std::string l;
      for (std::size_t j = 0; j < ux.size(); ++j) l += (j ? ", \"" : "\"") + ux[j] + "\"";
      return std::to_string(k) + " exposure value(s) are text, not numbers (" + l +
             "): a censored value such as \"<2\" is a measurement, and dropping it would bias the mean and the "
             "prevalence; replace them (for example with half the detection limit) or remove those rows";
    }
  }
  std::size_t drop = 0, neg = 0;
  for (double x : ex) {
    if (std::isnan(x)) ++drop;
    else if (x < 0) ++neg;
  }
  if (neg > 0) return std::to_string(neg) + " exposure value(s) are negative: a concentration cannot be negative";
  const int ii = coln(t, "income");
  if (ii >= 0) {
    const int ty = ctype(t.v[ii]);
    const std::vector<double> iv2 = asnum(t.v[ii], ty == 2 ? 1 : ty);
    bool allna = true;
    for (std::size_t i = 0; i < iv2.size(); ++i) {
      allna = allna && (ty == 2 ? t.v[ii][i] == "NA" : std::isnan(iv2[i]));
    }
    if (!allna && ty != 1) {
      return "the income column must be numeric (a quintile or an amount); text cannot be ranked";
    }
    if (!allna) {
      s.eq = true;
      s.inc = iv2;
    }
  }
  s.ex = ex;
  s.drop = static_cast<double>(drop);
  s.ds = f;
  return "";
}

// .vpexposed(): the mean and the share of the rows above the reference.
inline void exposed(Src &s, double ref) {
  std::vector<double> x, ab;
  for (double v : s.ex) {
    if (!std::isnan(v)) x.push_back(v);
  }
  for (double v : x) {
    if (!std::isfinite(ref) || v > ref) ab.push_back(v);
  }
  s.em = !ab.empty() ? rmean(ab) : !x.empty() ? rmean(x) : NA;
  s.ep = x.empty() ? NA : static_cast<double>(ab.size()) / static_cast<double>(x.size());
}

// revsort() of R's sort.c: a[1..n] into decreasing order by heapsort, ib
// alongside (the order of ties is the one sample() uses).
inline void revsort(double *a, int *ib, int n) {
  if (n <= 1) return;
  int l = (n >> 1) + 1, ir = n, i = 0, j = 0, ii = 0;
  double ra = 0;
  for (;;) {
    if (l > 1) {
      l = l - 1;
      ra = a[l];
      ii = ib[l];
    } else {
      ra = a[ir];
      ii = ib[ir];
      a[ir] = a[1];
      ib[ir] = ib[1];
      if (--ir == 1) {
        a[1] = ra;
        ib[1] = ii;
        return;
      }
    }
    i = l;
    j = l << 1;
    while (j <= ir) {
      if (j < ir && a[j] > a[j + 1]) ++j;
      if (ra > a[j]) {
        a[i] = a[j];
        ib[i] = ib[j];
        i = j;
        j += j;
      } else {
        j = ir + 1;
      }
    }
    a[i] = ra;
    ib[i] = ii;
  }
}

// .vpdemo(): rlnorm() and sample(1:5, n, TRUE, prob) after
// set.seed(20260417), as R draws them (sample()'s ProbSampleReplace()).
inline void demo(const std::string &poll, Src &s) {
  igh::Rng g(20260417u);
  const double ml = std::log(poll == "no2" ? 22.0 : poll == "pm25" ? 9.0 : 15.0), sl = poll == "pm25" ? 0.30 : 0.35;
  const std::size_t n = 1000;
  s.ex.assign(n, 0);
  s.inc.assign(n, 0);
  for (double &x : s.ex) x = std::exp(ml + sl * g.norm());
  double p[6] = {0, 0.18, 0.22, 0.22, 0.20, 0.18};
  int pm[6] = {0, 1, 2, 3, 4, 5};
  double sm = 0;
  for (int i = 1; i <= 5; ++i) sm += p[i];
  for (int i = 1; i <= 5; ++i) p[i] /= sm;
  revsort(p, pm, 5);
  for (int i = 2; i <= 5; ++i) p[i] += p[i - 1];
  static const double mult[5] = {1.25, 1.12, 1.00, 0.93, 0.85};
  for (std::size_t i = 0; i < n; ++i) {
    const double u = g.unif();
    int j = 1;
    while (j < 5 && u > p[j]) ++j;
    s.inc[i] = pm[j];
    s.ex[i] *= mult[pm[j] - 1];
  }
  s.eq = true;
  s.ds = "demo (synthetic)";
}

inline Val row(const char *a, bool ok, const std::string &note) {
  Val r = vobj();
  put(r, "assumption", vstr(a));
  put(r, "ok", vbool(ok));
  put(r, "note", vstr(note));
  return r;
}

// .vpasmlist(): the assumption log, and whether every one holds.
inline bool asm7(const std::string &poll, double em, double ep, double br, double pop, double ref, bool ra, bool hi,
                 Val &as) {
  auto f7 = [](double x) { return rfmt7(x, false); };
  const bool fe = std::isfinite(em), fb = std::isfinite(br), fp = std::isfinite(pop);
  std::string n1;
  if (!fe || em < 0) n1 = "expomean=" + f7(em) + " is not a finite non-negative concentration";
  else if (em <= ref) n1 = "mean " + rfmt7(em, false, 10) + " is at or below the counterfactual floor " + f7(ref) + ": no excess risk, PAF and cases are 0.";
  else n1 = "mean " + rfmt7(em, false, 10) + " vs ref " + f7(ref);
  const double top = hi ? 1e5 : 1e4;
  std::string nr;
  if (fb && br < 1e-6 && br > 0) {
    nr = "brate=" + f7(br) + " per 100k is below one case per 100 billion person-years: not a usable rate (check basrateper)";
  } else if (fb && br < 1 && br > 0 && !ra) {
    nr = "brate=" + f7(br) +
         " per 100k is below one case per 100,000. If it is a rate per person-year (0.012 is 1,200 per 100k), pass "
         "basrateper = 1; if the outcome really is this rare, say so with rare = TRUE";
  } else if (fb && br > top && br <= 1e5) {
    nr = "brate=" + rfmt7(br, true) +
         " per 100k a year is above a whole population's all-cause death rate (an age-specific rate for the oldest "
         "ages can be this high): if the outcome really is this common, say so with hirate = TRUE";
  } else if (fb && br == 0) {
    nr = "brate=0: a measured zero rate, so no cases";
  } else {
    nr = "brate=" + f7(br) + " per 100k per year";
  }
  const std::string rn = std::string("brate finite: zero, or between ") + (ra ? "1e-6" : "1") + " and " +
                         rfmt7(top, true) + " per 100k" + (ra ? " (a rare outcome)" : "");
  const bool ok[7] = {fe && em >= 0,
                      fe && em <= 1e4,
                      std::isfinite(ep) && ep >= 0 && ep <= 1,
                      fb && br >= 0 && br <= top && (br == 0 || br >= 1 || (ra && br >= 1e-6)),
                      fp && pop > 0 && pop <= 1e10,
                      std::isfinite(ref) && ref >= 0,
                      poll == "no2" || poll == "pm25"};
  as = varr();
  as.a.push_back(row("exposure finite and non-negative", ok[0], n1));
  as.a.push_back(row("exposure within the concentration-response range (at most 10,000 ug/m3)", ok[1],
                     "exposure_mean=" + f7(em)));
  as.a.push_back(row("prevalence in [0,1]", ok[2], "exposure_prevalence=" + f7(ep)));
  as.a.push_back(row(rn.c_str(), ok[3], nr));
  as.a.push_back(row("population finite, positive and at most 10 billion", ok[4],
                     "population=" + rfmt7(pop, true) + (fp && pop > 1e10 ? " is more people than live on Earth" : "")));
  as.a.push_back(row("reference finite and non-negative", ok[5], "reference=" + f7(ref)));
  as.a.push_back(row("pollutant supported by the CRFs", ok[6],
                     "Current CRFs: NO2 (log-linear), PM2.5 (log-linear all-cause; Burnett IER for IHD and stroke). "
                     "Other pollutants reject."));
  bool all = true;
  for (bool b : ok) all = all && b;
  return all;
}

// .vppipeline(): the stages of japol.inc in R's order; false with the
// first stage's refusal in err.
inline bool stages(const std::string &poll, const std::string &y, const Src &s, double br, double pop, double ref,
                   bool ra, bool hi, Val &pp, std::string &err) {
  const bool pm = poll == "pm25";
  const double em = s.em, ep = s.ep;
  Val a1 = vobj();
  put(a1, "expo", vnum(em));
  put(a1, "y", vstr(y));
  put(a1, "refconc", vnum(ref));
  Ctx c1{pm ? "crfpm25" : "crfno2", &a1, std::string(), {}};
  Val crf;
  if (!(pm ? crfpm25(c1, crf) : crfno2(c1, crf))) {
    err = c1.err;
    return false;
  }
  const double rr = col(crf, "rr").n, lo = col(crf, "rrlo").n, up = col(crf, "rrup").n;
  Val a2 = vobj();
  put(a2, "rr", vnum(rr));
  put(a2, "expoprev", vnum(ep));
  if (std::isfinite(lo)) put(a2, "rrlo", vnum(lo));
  if (std::isfinite(up)) put(a2, "rrup", vnum(up));
  Ctx c2{"attrfrac", &a2, std::string(), {}};
  Val paf, ci;
  if (!attrfrac(c2, paf)) {
    err = c2.err;
    return false;
  }
  if (std::isfinite(lo) && std::isfinite(up)) {
    ci = vobj();
    put(ci, "lower", vnum(clab_levin(lo, ep)));
    put(ci, "upper", vnum(clab_levin(up, ep)));
  }
  const double bpp = br / 1e5, del = std::max(0.0, em - ref), bpu = std::log(rr) / std::max(del, 1e-9);
  Val a3 = vobj();
  put(a3, "expdel", vnum(del));
  put(a3, "pop", vnum(pop * ep));
  put(a3, "brate", vnum(bpp));
  put(a3, "betaperunit", vnum(bpu));
  put(a3, "rare", vbool(ra));
  put(a3, "hirate", vbool(hi));
  Ctx c3{"mortdisp", &a3, std::string(), {}};
  Val dd;
  if (!mortdisp(c3, dd)) {
    err = c3.err;
    return false;
  }
  Val a4 = vobj();
  put(a4, "expomean", vnum(em));
  put(a4, "expoprev", vnum(ep));
  put(a4, "brate", vnum(bpp));
  put(a4, "pop", vnum(pop));
  put(a4, "poll", vstr(pm ? "PM2.5" : "NO2"));
  put(a4, "y", vstr(y));
  put(a4, "refconc", vnum(ref));
  put(a4, "rare", vbool(ra));
  put(a4, "hirate", vbool(hi));
  Ctx c4{"polbrdn", &a4, std::string(), {}};
  Val bd;
  if (!polbrdn(c4, bd)) {
    err = c4.err;
    return false;
  }
  // .vpequity(): with incomes, two complete rows and an exposure above 0
  Val eq;
  if (s.eq) {
    std::size_t k = 0;
    bool pos = false;
    Val ec = varr(), ic = varr();
    for (std::size_t i = 0; i < s.ex.size(); ++i) {
      k += std::isfinite(s.ex[i]) && std::isfinite(s.inc[i]) ? 1 : 0;
      pos = pos || s.ex[i] > 0;
      // a missing value is null, which expconcidx() drops with its row
      ec.a.push_back(std::isnan(s.ex[i]) ? Val() : vnum(s.ex[i]));
      ic.a.push_back(std::isnan(s.inc[i]) ? Val() : vnum(s.inc[i]));
    }
    if (k >= 2 && pos) {
      Val d = vobj();
      put(d, "exposure", ec);
      put(d, "income", ic);
      Val a5 = vobj();
      put(a5, "data", d);
      put(a5, "expo", vstr("exposure"));
      put(a5, "income", vstr("income"));
      Ctx c5{"expconcidx", &a5, std::string(), {}};
      if (!expconcidx(c5, eq)) {
        err = c5.err;
        return false;
      }
    }
  }
  Val dp = vobj();
  put(dp, "deadisp", dd);
  put(dp, "expdel", vnum(del));
  put(dp, "betaperunit", vnum(bpu));
  put(dp, "unit", vstr(y == "chiast" ? "incident cases" : "deaths"));
  pp = vobj();
  put(pp, "crf", crf);
  put(pp, "paf", paf);
  put(pp, "pafci", ci);
  put(pp, "displaced", dp);
  put(pp, "burden", bd);
  put(pp, "equity", eq);
  return true;
}

// vrfpol(): R's checks in their order; an argument the pipeline cannot
// use is an error report, not a refusal, as in R.
inline bool run(Ctx &c, Val &out) {
  const Val *ps = one(c["poll"], Val::STR);
  if (ps == nullptr) return ckstop(c, "`poll` must be one pollutant name (\"PM2.5\" or \"NO2\")");
  const std::string pk = polkey(c["poll"]);
  const std::string poll = pk == "NO2" ? "no2" : pk == "PM2.5" ? "pm25" : lower(ps->s);
  auto fail = [&](const std::string &m) {
    out = vobj();
    put(out, "status", vstr("error"));
    put(out, "error", vstr(m));
    return true;
  };
  auto flag = [&](const char *k, bool &b) {
    b = false;
    if (!c.has(k)) return true;
    const Val *v = one(c[k], Val::BOOL);
    if (v != nullptr) b = v->b;
    return v != nullptr;
  };
  bool dm = false, ra = false, hi = false;
  if (!flag("demo", dm)) return fail("demo must be TRUE or FALSE");
  // .vpsrc()
  const bool csv = c.has("expcsv"), hem = c.has("expomean"), hep = c.has("expoprev");
  if (dm && csv) {
    return fail("expcsv and demo = TRUE are two sources of exposure: pass one (the demo was used silently before)");
  }
  if ((dm || csv) && (hem || hep)) {
    return fail("expomean / expoprev cannot be combined with expcsv or demo: pass one source of exposure");
  }
  if (hem != hep) return fail("expomean and expoprev go together: pass both (a zero is a measurement), or neither");
  // .vpnums()
  static const char *const nm[6] = {"expomean", "expoprev", "brate", "basrateper", "pop", "ref"};
  const double def[6] = {NA, NA, 500, 1e5, 1e6, poll == "no2" ? 10 : 5.8};
  std::vector<double> v[6];
  for (int i = 0; i < 6; ++i) {
    if (!c.has(nm[i])) {
      v[i] = {def[i]};
      continue;
    }
    bool lg = false;
    const char *cl = jcls(c[nm[i]], lg);
    if (cl != nullptr) return fail(std::string(nm[i]) + " must be a single number, not " + cl);
    Dv d;
    readd(c[nm[i]], d);
    for (double &x : d.v) x = na(x);
    v[i] = d.v;
  }
  if (!flag("rare", ra)) return fail("rare must be TRUE or FALSE");
  if (!flag("hirate", hi)) return fail("hirate must be TRUE or FALSE");
  if (v[3].size() != 1 || !std::isfinite(v[3][0]) || v[3][0] <= 0) {
    return fail("basrateper must be a single positive number (100000 for a rate per 100,000, 1 per person-year)");
  }
  for (int i : {5, 0, 1, 2, 4}) {
    if (v[i].size() != 1 && (i > 1 || c.has(nm[i]))) {
      return fail(std::string(nm[i]) + " must be a single number, not " + std::to_string(v[i].size()) + " values");
    }
  }
  std::string y = "allcause";
  if (c.has("y")) {
    const Val *yv = one(c["y"], Val::STR);
    if (yv == nullptr) return fail("outcome must be a single outcome name");
    y = yv->s;
  }
  const double ref = v[5][0];
  // .vpsource()
  Src s;
  if (dm) {
    demo(poll, s);
    exposed(s, ref);
  } else if (csv) {
    const Val *f = one(c["expcsv"], Val::STR);
    if (f == nullptr) return fail("expcsv must be a single file path");
    const std::string m = vpread(f->s, poll, s);
    if (!m.empty()) return fail(m);
    exposed(s, ref);
  } else if (!hem) {
    return fail("no exposure data: pass expcsv, or expomean and expoprev, or demo = TRUE");
  } else {
    s.em = v[0][0];
    s.ep = v[1][0];
    s.ds = "scalar arguments";
  }
  // .vprun(), the rate per 100,000 from here on
  const double br = v[2][0] * 1e5 / v[3][0], pop = v[4][0];
  Val as;
  const bool ok = asm7(poll, s.em, s.ep, br, pop, ref, ra, hi, as);
  Val pp;
  if (ok) {
    std::string err;
    if (!stages(poll, y, s, br, pop, ref, ra, hi, pp, err)) return fail(err);
  } else {
    pp = vobj();
    put(pp, "skipped", vbool(true));
    put(pp, "reason", vstr("assumptions"));
  }
  Val in = vobj();
  put(in, "expomean", vnum(s.em));
  put(in, "expoprev", vnum(s.ep));
  put(in, "baserateper100k", vnum(br));
  put(in, "population", vnum(pop));
  put(in, "refconc", vnum(ref));
  put(in, "drorows", vnum(s.drop));
  out = vobj();
  put(out, "command", vstr("verify-pollution"));
  put(out, "pollutant", vstr(poll));
  put(out, "outcome", vstr(y));
  put(out, "region", c["reg"]);
  put(out, "years", c["years"]);
  put(out, "dsou", vstr(s.ds));
  put(out, "inputs", in);
  put(out, "assumptions", as);
  put(out, "status", vstr(ok ? "ok" : "assumption_failure"));
  put(out, "pipeline", pp);
  return true;
}

}  // namespace vrf

bool vrfpol(Ctx &c, Val &out) { return vrf::run(c, out); }
