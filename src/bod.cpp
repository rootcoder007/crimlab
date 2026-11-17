// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The readers of odls() and odtab(): the JSON reader of the core (src/json.h,
// the one the command-line program reads its arguments with) for the
// manifest, and a CSV reader (RFC 4180) for the tables.
//
// JSON:
// An array is a list, an object a named list (an empty one has empty
// names, as jsonlite gives it), a string a character value in UTF-8, a number a double, true
// and false logicals, null NULL. Strings must be valid UTF-8 without NUL.
// The containers being read are kept in one protected list, two slots a
// level (items and keys), grown by doubling, so the work and the memory
// are linear in the text and its depth is limited only by maxd.
#include <climits>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "json.h"
#include "rchk.h"

namespace {

bool utf8ok(const std::string &s) { return clab::u8ok(s.data(), s.size()); }

class RB {
 public:
  RB(SEXP stk, int maxd) : stk_(stk), top_(2 * (maxd + 1)) {}
  SEXP root() const { return VECTOR_ELT(stk_, top_); }
  bool open(bool obj) {
    const std::size_t l = len_.size();
    SET_VECTOR_ELT(stk_, 2 * l, Rf_allocVector(VECSXP, 4));
    SET_VECTOR_ELT(stk_, 2 * l + 1, obj ? Rf_allocVector(STRSXP, 4) : R_NilValue);
    len_.push_back(0);
    obj_.push_back(obj);
    return true;
  }
  bool close() {
    const std::size_t l = len_.size() - 1;
    const R_xlen_t n = len_[l];
    SEXP it = VECTOR_ELT(stk_, 2 * l), ks = VECTOR_ELT(stk_, 2 * l + 1);
    SEXP v = PROTECT(Rf_allocVector(VECSXP, n));
    for (R_xlen_t i = 0; i < n; ++i) SET_VECTOR_ELT(v, i, VECTOR_ELT(it, i));
    if (obj_[l]) {
      SEXP nm = PROTECT(Rf_allocVector(STRSXP, n));
      for (R_xlen_t i = 0; i < n; ++i) SET_STRING_ELT(nm, i, STRING_ELT(ks, i));
      Rf_setAttrib(v, R_NamesSymbol, nm);
      UNPROTECT(1);
    }
    SET_VECTOR_ELT(stk_, 2 * l, R_NilValue);
    SET_VECTOR_ELT(stk_, 2 * l + 1, R_NilValue);
    len_.pop_back();
    obj_.pop_back();
    add(v);
    UNPROTECT(1);
    return true;
  }
  bool key(std::string k) {
    if (!utf8ok(k)) return false;
    grow();
    const std::size_t l = len_.size() - 1;
    SET_STRING_ELT(VECTOR_ELT(stk_, 2 * l + 1), len_[l],
                   Rf_mkCharLenCE(k.data(), static_cast<int>(k.size()), CE_UTF8));
    return true;
  }
  bool str(std::string s) {
    if (!utf8ok(s) || s.size() > static_cast<std::size_t>(INT_MAX)) return false;
    SEXP v = PROTECT(Rf_ScalarString(Rf_mkCharLenCE(s.data(), static_cast<int>(s.size()), CE_UTF8)));
    add(v);
    UNPROTECT(1);
    return true;
  }
  bool num(double d) { return put(Rf_ScalarReal(d)); }
  bool boolean(bool b) { return put(Rf_ScalarLogical(b ? 1 : 0)); }
  bool null() {
    add(R_NilValue);
    return true;
  }

 private:
  SEXP stk_;
  R_xlen_t top_;
  std::vector<R_xlen_t> len_;
  std::vector<bool> obj_;
  bool put(SEXP v) {
    PROTECT(v);
    add(v);
    UNPROTECT(1);
    return true;
  }
  // room for one more item (and its key) at the top level
  void grow() {
    const std::size_t l = len_.size() - 1;
    SEXP it = VECTOR_ELT(stk_, 2 * l);
    const R_xlen_t n = XLENGTH(it);
    if (len_[l] < n) return;
    SEXP ni = PROTECT(Rf_allocVector(VECSXP, 2 * n));
    for (R_xlen_t i = 0; i < n; ++i) SET_VECTOR_ELT(ni, i, VECTOR_ELT(it, i));
    SET_VECTOR_ELT(stk_, 2 * l, ni);
    UNPROTECT(1);
    if (obj_[l]) {
      SEXP ks = VECTOR_ELT(stk_, 2 * l + 1), nk = PROTECT(Rf_allocVector(STRSXP, 2 * n));
      for (R_xlen_t i = 0; i < n; ++i) SET_STRING_ELT(nk, i, STRING_ELT(ks, i));
      SET_VECTOR_ELT(stk_, 2 * l + 1, nk);
      UNPROTECT(1);
    }
  }
  void add(SEXP v) {
    if (len_.empty()) {
      SET_VECTOR_ELT(stk_, top_, v);
      return;
    }
    const std::size_t l = len_.size() - 1;
    // an object's key has made room for its value already
    if (!obj_[l]) grow();
    SET_VECTOR_ELT(VECTOR_ELT(stk_, 2 * l), len_[l]++, v);
  }
};


// CSV (RFC 4180) from the bytes of a file: records end at LF or CRLF (the
// last may lack it), fields are separated by commas, a field in double
// quotes may hold commas, line ends and doubled quotes; a UTF-8 byte order
// mark is skipped. A blank line is a row of a one-column table (its value
// missing) and no row of a wider one, as read.csv() skips it. A field that
// is empty or NA is
// missing. Bytes that are not valid UTF-8 are written <xx>, as the
// package's messages write them. The first record is the header, its
// names as they are.
struct Csv {
  int code = 0;  // 0, or the refusal: 1 ragged, 2 a quote inside a field
                 // not quoted, 3 a quote never closed, 4 text after a
                 // closing quote, 5 a NUL byte, 6 no header
  long line = 0, nf = 0, nh = 0;
  std::vector<std::vector<std::string>> col;  // the header's fields first
  std::vector<std::vector<char>> na;
};

// the field f, valid UTF-8 (invalid bytes as <xx>)
std::string u8fix(const std::string &f) {
  if (utf8ok(f)) return f;
  std::string o;
  for (std::size_t i = 0; i < f.size();) {
    // the longest valid prefix from i: one character, or the byte escaped
    std::size_t k = 1;
    const unsigned char c = static_cast<unsigned char>(f[i]);
    const std::size_t n = c >= 0xf0 ? 4 : c >= 0xe0 ? 3 : c >= 0xc0 ? 2 : 1;
    if (n > 1 && i + n <= f.size() && utf8ok(f.substr(i, n))) k = n;
    if (k == 1 && c >= 0x80) {
      char b[8];
      std::snprintf(b, sizeof b, "<%02x>", c);
      o += b;
    } else {
      o.append(f, i, k);
    }
    i += k;
  }
  return o;
}

void csv(const char *p, const char *e, Csv &r) {
  if (e - p >= 3 && std::memcmp(p, "\xef\xbb\xbf", 3) == 0) p += 3;
  long line = 1, start = 1;
  std::vector<std::string> rec;
  std::vector<char> rna;
  std::string f;
  bool blank = false;  // the record is a blank line
  auto fail = [&](int c, long l) {
    r.code = c;
    r.line = l;
  };
  // the record just read: the header, or a row of as many fields
  auto put = [&]() {
    if (blank && r.nh > 1) {
      rec.clear();
      rna.clear();
      return;
    }
    if (r.col.empty()) {
      r.nh = static_cast<long>(rec.size());
      r.col.resize(rec.size());
      r.na.resize(rec.size());
      for (std::size_t j = 0; j < rec.size(); ++j) r.col[j].push_back(rec[j]);
      for (std::size_t j = 0; j < rec.size(); ++j) r.na[j].push_back(0);
    } else if (static_cast<long>(rec.size()) != r.nh) {
      fail(1, start);
      r.nf = static_cast<long>(rec.size());
      return;
    } else {
      for (std::size_t j = 0; j < rec.size(); ++j) {
        r.col[j].push_back(rec[j]);
        r.na[j].push_back(rna[j]);
      }
    }
    rec.clear();
    rna.clear();
  };
  while (p < e && r.code == 0) {
    // one field
    f.clear();
    blank = rec.empty() && (*p == '\n' || (*p == '\r' && p + 1 < e && p[1] == '\n'));
    if (*p == '"') {
      const long ql = line;
      ++p;
      for (;;) {
        if (p == e) return fail(3, ql);
        if (*p == '"') {
          if (p + 1 < e && p[1] == '"') {
            f += '"';
            p += 2;
            continue;
          }
          ++p;
          break;
        }
        if (*p == 0) return fail(5, line);
        if (*p == '\n') ++line;
        f += *p++;
      }
      if (p < e && *p != ',' && *p != '\n' && !(*p == '\r' && p + 1 < e && p[1] == '\n')) {
        return fail(4, line);
      }
    } else {
      while (p < e && *p != ',' && *p != '\n' && !(*p == '\r' && p + 1 < e && p[1] == '\n')) {
        if (*p == '"') return fail(2, line);
        if (*p == 0) return fail(5, line);
        f += *p++;
      }
    }
    rec.push_back(u8fix(f));
    rna.push_back(f.empty() || f == "NA");
    if (p < e && *p == ',') {
      ++p;
      if (p == e) {
        // a comma at the very end: one more, empty, field
        rec.emplace_back();
        rna.push_back(1);
      } else {
        continue;
      }
    } else if (p < e) {
      p += *p == '\r' ? 2 : 1;
      ++line;
    }
    put();
    start = line;
  }
  if (r.code == 0 && r.col.empty()) r.code = 6;
}

}  // namespace

extern "C" {

//' @noRd
SEXP C_clab_r_json(SEXP txt) {
  if (TYPEOF(txt) != RAWSXP) clab::fail("the text must be raw bytes.");
  const int maxd = 10000;
  SEXP stk = PROTECT(Rf_allocVector(VECSXP, 2 * (maxd + 1) + 1));
  RB b(stk, maxd);
  clabj::Parser p(reinterpret_cast<const char *>(RAW(txt)), static_cast<std::size_t>(XLENGTH(txt)));
  const bool ok = p.run(b, maxd);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(r, 0, Rf_ScalarLogical(ok ? 1 : 0));
  if (ok) SET_VECTOR_ELT(r, 1, b.root());
  UNPROTECT(2);
  return r;
}

}  // extern "C"

extern "C" {

//' @noRd
SEXP C_clab_r_csv(SEXP txt) {
  if (TYPEOF(txt) != RAWSXP) clab::fail("the text must be raw bytes.");
  Csv r;
  const char *p = reinterpret_cast<const char *>(RAW(txt));
  csv(p, p + XLENGTH(txt), r);
  SEXP o = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP st = PROTECT(Rf_allocVector(REALSXP, 4));
  REAL(st)[0] = r.code;
  REAL(st)[1] = static_cast<double>(r.line);
  REAL(st)[2] = static_cast<double>(r.nf);
  REAL(st)[3] = static_cast<double>(r.nh);
  SET_VECTOR_ELT(o, 0, st);
  UNPROTECT(1);
  if (r.code == 0) {
    const R_xlen_t nc = static_cast<R_xlen_t>(r.col.size());
    SEXP nm = PROTECT(Rf_allocVector(STRSXP, nc)), cl = PROTECT(Rf_allocVector(VECSXP, nc));
    for (R_xlen_t j = 0; j < nc; ++j) {
      const std::vector<std::string> &c = r.col[static_cast<std::size_t>(j)];
      SET_STRING_ELT(nm, j, Rf_mkCharLenCE(c[0].data(), static_cast<int>(c[0].size()), CE_UTF8));
      SEXP v = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(c.size()) - 1));
      for (std::size_t i = 1; i < c.size(); ++i) {
        SET_STRING_ELT(v, static_cast<R_xlen_t>(i) - 1,
                       r.na[static_cast<std::size_t>(j)][i]
                           ? NA_STRING
                           : Rf_mkCharLenCE(c[i].data(), static_cast<int>(c[i].size()), CE_UTF8));
      }
      SET_VECTOR_ELT(cl, j, v);
      UNPROTECT(1);
    }
    SET_VECTOR_ELT(o, 1, nm);
    SET_VECTOR_ELT(o, 2, cl);
    UNPROTECT(2);
  }
  UNPROTECT(1);
  return o;
}

}  // extern "C"
