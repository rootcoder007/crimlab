// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// JSON for clab_jrun() and the command-line program: a value, a reader
// (strict RFC 8259, at most 64 levels deep) and a writer.
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
  bool doc(Val &v) {
    if (!val(v, 0)) return false;
    ws();
    return p_ == e_;
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
  bool val(Val &v, int depth) {
    ws();
    if (p_ == e_ || depth > 64) return false;
    const char c = *p_;
    if (c == '{' || c == '[') {
      ++p_;
      v.t = c == '{' ? Val::OBJ : Val::ARR;
      const char end = c == '{' ? '}' : ']';
      ws();
      if (p_ < e_ && *p_ == end) {
        ++p_;
        return true;
      }
      for (;;) {
        if (v.t == Val::OBJ) {
          ws();
          std::string k;
          if (!str(k)) return false;
          ws();
          if (p_ == e_ || *p_ != ':') return false;
          ++p_;
          v.k.push_back(std::move(k));
        }
        Val x;
        if (!val(x, depth + 1)) return false;
        v.a.push_back(std::move(x));
        ws();
        if (p_ == e_) return false;
        if (*p_ == ',') {
          ++p_;
          continue;
        }
        if (*p_ != end) return false;
        ++p_;
        return true;
      }
    }
    if (c == '"') {
      v.t = Val::STR;
      return str(v.s);
    }
    if (c == 't' || c == 'f') {
      v.t = Val::BOOL;
      v.b = c == 't';
      return lit(v.b ? "true" : "false");
    }
    if (c == 'n') return lit("null");
    v.t = Val::NUM;
    return num(v.n);
  }
};

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
