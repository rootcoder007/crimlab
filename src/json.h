// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// JSON for clab_jrun(), the command-line program and the R package's
// manifest reader (src/bjsn.cpp): a value, a reader (strict RFC 8259) and a
// writer. The reader is iterative, so its depth is a limit of its caller's
// (64 for arguments, 10000 for a manifest), never of the stack; it hands
// what it reads to a builder (ValB here, an R one in src/bjsn.cpp).
#ifndef CLAB_JSON_H
#define CLAB_JSON_H

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "dfmt.h"

namespace clabj {

inline const double NaN = std::numeric_limits<double>::quiet_NaN();

// R's NA_real_: a NaN whose low word is 1954. Written as null.
inline double mkna() {
  const std::uint64_t u = 0x7ff00000000007a2ULL;
  double d;
  std::memcpy(&d, &u, sizeof d);
  return d;
}
inline const double NA = mkna();
inline bool isna(double x) {
  std::uint64_t u;
  std::memcpy(&u, &x, sizeof u);
  return std::isnan(x) && (u & 0xffffffffULL) == 1954;
}

// JSON.

struct Val {
  enum T { NUL, BOOL, NUM, STR, ARR, OBJ };
  T t = NUL;
  bool b = false;
  double n = 0;
  std::string s;
  std::vector<std::string> k;  // an object's keys
  std::vector<Val> a;          // an array's items, an object's values
};

inline const Val nul;

class Parser {
 public:
  Parser(const char *s, std::size_t n) : p_(s), e_(s + n) {}
  // One JSON text into v, at most 64 levels deep (arguments).
  bool doc(Val &v);
  // One JSON text to the builder b: open(obj), key(k), str(s), num(d),
  // boolean(b), null() and close(), each false to stop. At most maxd
  // levels of arrays and objects inside the outermost.
  template <class B>
  bool run(B &b, int maxd) {
    std::vector<char> st;  // the closing bracket of each open container
    for (;;) {
      ws();
      if (p_ == e_ || static_cast<int>(st.size()) > maxd) return false;
      const char c = *p_;
      if (c == '{' || c == '[') {
        ++p_;
        const char end = c == '{' ? '}' : ']';
        if (!b.open(c == '{')) return false;
        ws();
        if (p_ < e_ && *p_ == end) {
          ++p_;
          if (!b.close()) return false;
        } else {
          st.push_back(end);
          if (end == '}' && !key(b)) return false;
          continue;
        }
      } else if (!scalar(b)) {
        return false;
      }
      // after a value: a comma (and the next key) or closing brackets
      for (;;) {
        ws();
        if (st.empty()) return p_ == e_;
        if (p_ == e_) return false;
        if (*p_ == ',') {
          ++p_;
          if (st.back() == '}' && !key(b)) return false;
          break;
        }
        if (*p_ != st.back()) return false;
        ++p_;
        st.pop_back();
        if (!b.close()) return false;
      }
    }
  }

 private:
  const char *p_, *e_;

  void ws() {
    while (p_ < e_ && (*p_ == ' ' || *p_ == '\t' || *p_ == '\n' || *p_ == '\r')) ++p_;
  }
  bool lit(const char *w) {
    const std::size_t n = std::strlen(w);
    if (static_cast<std::size_t>(e_ - p_) < n || std::memcmp(p_, w, n) != 0) return false;
    p_ += n;
    return true;
  }
  bool hex4(unsigned &u) {
    if (e_ - p_ < 4) return false;
    u = 0;
    for (int i = 0; i < 4; ++i) {
      const char c = *p_++;
      u <<= 4;
      if (c >= '0' && c <= '9') u |= static_cast<unsigned>(c - '0');
      else if (c >= 'a' && c <= 'f') u |= static_cast<unsigned>(c - 'a' + 10);
      else if (c >= 'A' && c <= 'F') u |= static_cast<unsigned>(c - 'A' + 10);
      else return false;
    }
    return true;
  }
  static void utf8(std::string &s, unsigned u) {
    if (u < 0x80) {
      s += static_cast<char>(u);
    } else if (u < 0x800) {
      s += static_cast<char>(0xc0 | (u >> 6));
      s += static_cast<char>(0x80 | (u & 0x3f));
    } else if (u < 0x10000) {
      s += static_cast<char>(0xe0 | (u >> 12));
      s += static_cast<char>(0x80 | ((u >> 6) & 0x3f));
      s += static_cast<char>(0x80 | (u & 0x3f));
    } else {
      s += static_cast<char>(0xf0 | (u >> 18));
      s += static_cast<char>(0x80 | ((u >> 12) & 0x3f));
      s += static_cast<char>(0x80 | ((u >> 6) & 0x3f));
      s += static_cast<char>(0x80 | (u & 0x3f));
    }
  }
  bool str(std::string &s) {
    if (p_ == e_ || *p_ != '"') return false;
    ++p_;
    while (p_ < e_ && *p_ != '"') {
      const unsigned char c = static_cast<unsigned char>(*p_++);
      if (c < 0x20) return false;
      if (c != '\\') {
        s += static_cast<char>(c);
        continue;
      }
      if (p_ == e_) return false;
      const char x = *p_++;
      unsigned u = 0, w = 0;
      switch (x) {
        case '"': case '\\': case '/': s += x; break;
        case 'b': s += '\b'; break;
        case 'f': s += '\f'; break;
        case 'n': s += '\n'; break;
        case 'r': s += '\r'; break;
        case 't': s += '\t'; break;
        case 'u':
          if (!hex4(u) || (u >= 0xdc00 && u < 0xe000)) return false;
          if (u >= 0xd800 && u < 0xdc00) {
            if (!lit("\\u") || !hex4(w) || w < 0xdc00 || w >= 0xe000) return false;
            u = 0x10000 + ((u - 0xd800) << 10) + (w - 0xdc00);
          }
          utf8(s, u);
          break;
        default: return false;
      }
    }
    if (p_ == e_) return false;
    ++p_;
    return true;
  }
  bool digits() {
    const char *b = p_;
    while (p_ < e_ && *p_ >= '0' && *p_ <= '9') ++p_;
    return p_ > b;
  }
  bool num(double &d) {
    const char *b = p_;
    if (p_ < e_ && *p_ == '-') ++p_;
    if (p_ < e_ && *p_ == '0') ++p_;
    else if (!digits()) return false;
    if (p_ < e_ && *p_ == '.') {
      ++p_;
      if (!digits()) return false;
    }
    if (p_ < e_ && (*p_ == 'e' || *p_ == 'E')) {
      ++p_;
      if (p_ < e_ && (*p_ == '+' || *p_ == '-')) ++p_;
      if (!digits()) return false;
    }
    d = std::strtod(std::string(b, p_).c_str(), nullptr);
    return true;
  }
  template <class B>
  bool key(B &b) {
    ws();
    std::string k;
    if (!str(k)) return false;
    ws();
    if (p_ == e_ || *p_ != ':') return false;
    ++p_;
    return b.key(std::move(k));
  }
  template <class B>
  bool scalar(B &b) {
    const char c = *p_;
    if (c == '"') {
      std::string s;
      return str(s) && b.str(std::move(s));
    }
    if (c == 't' || c == 'f') return lit(c == 't' ? "true" : "false") && b.boolean(c == 't');
    if (c == 'n') return lit("null") && b.null();
    double d;
    return num(d) && b.num(d);
  }
};

// The builder of a Val: the open containers on a stack, each the last
// item of the one below it (only the top one grows, so none moves).
class ValB {
 public:
  explicit ValB(Val &root) : root_(root) {}
  bool open(bool obj) {
    Val &v = slot();
    v.t = obj ? Val::OBJ : Val::ARR;
    st_.push_back(&v);
    return true;
  }
  bool close() {
    st_.pop_back();
    return true;
  }
  bool key(std::string k) {
    st_.back()->k.push_back(std::move(k));
    return true;
  }
  bool str(std::string s) {
    Val &v = slot();
    v.t = Val::STR;
    v.s = std::move(s);
    return true;
  }
  bool num(double d) {
    Val &v = slot();
    v.t = Val::NUM;
    v.n = d;
    return true;
  }
  bool boolean(bool b) {
    Val &v = slot();
    v.t = Val::BOOL;
    v.b = b;
    return true;
  }
  bool null() {
    slot();
    return true;
  }

 private:
  Val &root_;
  std::vector<Val *> st_;
  Val &slot() {
    if (st_.empty()) return root_;
    st_.back()->a.emplace_back();
    return st_.back()->a.back();
  }
};

inline bool Parser::doc(Val &v) {
  ValB b(v);
  return run(b, 64);
}

inline void wstr(std::string &o, const std::string &s) {
  o += '"';
  for (const char ch : s) {
    const unsigned char c = static_cast<unsigned char>(ch);
    if (c == '"' || c == '\\') {
      o += '\\';
      o += ch;
    } else if (c < 0x20) {
      char b[8];
      std::snprintf(b, sizeof b, "\\u%04x", c);
      o += b;
    } else {
      o += ch;
    }
  }
  o += '"';
}

// The shortest text that reads back the same (src/dfmt.h).
inline void wnum(std::string &o, double x) {
  if (isna(x)) {
    o += "null";
  } else if (std::isnan(x)) {
    o += "\"NaN\"";
  } else if (std::isinf(x)) {
    o += x > 0 ? "\"Inf\"" : "\"-Inf\"";
  } else {
    char b[32];
    o.append(b, static_cast<std::size_t>(dfmt::put(x, b)));
  }
}

inline void wval(std::string &o, const Val &v) {
  switch (v.t) {
    case Val::NUL: o += "null"; break;
    case Val::BOOL: o += v.b ? "true" : "false"; break;
    case Val::NUM: wnum(o, v.n); break;
    case Val::STR: wstr(o, v.s); break;
    case Val::ARR:
    case Val::OBJ:
      o += v.t == Val::ARR ? '[' : '{';
      for (std::size_t i = 0; i < v.a.size(); ++i) {
        if (i) o += ',';
        if (v.t == Val::OBJ) {
          wstr(o, v.k[i]);
          o += ':';
        }
        wval(o, v.a[i]);
      }
      o += v.t == Val::ARR ? ']' : '}';
      break;
  }
}

}  // namespace clabj

#endif
