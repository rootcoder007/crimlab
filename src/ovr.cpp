// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Coding police oversight reports (R/ovr.R): the text of a report page,
// written-out dates, the fields of a schema read by rules, the count of
// subject officials, keyword indicators and the checks of a coded table.
// Everything is a single pass over bytes: no regular expressions, so no
// input can make a pattern backtrack or recurse. Text is UTF-8; letters
// outside ASCII are compared as their bytes.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

using std::size_t;
using std::string;

inline bool isdig(unsigned char c) { return c >= '0' && c <= '9'; }
inline bool isup(unsigned char c) { return c >= 'A' && c <= 'Z'; }
inline bool islow(unsigned char c) { return c >= 'a' && c <= 'z'; }
// a byte of a word: an ASCII letter or digit, or a byte of a character
// outside ASCII (a letter, in the reports' languages)
inline bool iswrd(unsigned char c) { return isdig(c) || isup(c) || islow(c) || c >= 0x80; }
inline unsigned char low(unsigned char c) { return isup(c) ? static_cast<unsigned char>(c + 32) : c; }
inline bool isws(unsigned char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v'; }

// case-insensitive (ASCII) match of lit at s[i]
bool ieq(const string &s, size_t i, const char *lit) {
  const size_t n = std::strlen(lit);
  if (i + n > s.size()) return false;
  for (size_t k = 0; k < n; ++k) {
    if (low(static_cast<unsigned char>(s[i + k])) != low(static_cast<unsigned char>(lit[k]))) return false;
  }
  return true;
}

// The HTML elements whose start and end break a line.
const char *const kBlock[] = {"p",       "div",    "br",     "h1",      "h2",     "h3",         "h4",
                              "h5",      "h6",     "tr",     "li",      "ul",     "ol",         "dl",
                              "dt",      "dd",     "table",  "section", "article", "header",    "footer",
                              "main",    "nav",    "aside",  "blockquote", "pre", "title",      "hr",
                              "figure",  "figcaption", "address", "form", "fieldset", "legend", "caption"};

// The name of the tag whose text starts at s[i] (after '<' and an
// optional '/'), lower case.
string tagname(const string &s, size_t i) {
  string t;
  while (i < s.size() && (isup(static_cast<unsigned char>(s[i])) || islow(static_cast<unsigned char>(s[i])) ||
                          isdig(static_cast<unsigned char>(s[i])))) {
    t += static_cast<char>(low(static_cast<unsigned char>(s[i])));
    ++i;
  }
  return t;
}

// The named entities of report pages, and the characters they stand for.
const char *const kEnt[][2] = {
    {"nbsp", " "},          {"amp", "&"},           {"quot", "\""},         {"lt", "<"},
    {"gt", ">"},            {"apos", "'"},          {"rsquo", "'"},         {"lsquo", "'"},
    {"ldquo", "\""},        {"rdquo", "\""},        {"ndash", "-"},         {"mdash", "--"},
    {"hellip", "..."},      {"thinsp", " "},        {"ensp", " "},          {"emsp", " "},
    {"shy", ""},            {"eacute", "\xc3\xa9"}, {"egrave", "\xc3\xa8"}, {"ecirc", "\xc3\xaa"},
    {"euml", "\xc3\xab"},   {"agrave", "\xc3\xa0"}, {"acirc", "\xc3\xa2"},  {"ccedil", "\xc3\xa7"},
    {"icirc", "\xc3\xae"},  {"iuml", "\xc3\xaf"},   {"ocirc", "\xc3\xb4"},  {"ouml", "\xc3\xb6"},
    {"ucirc", "\xc3\xbb"},  {"ugrave", "\xc3\xb9"}, {"uuml", "\xc3\xbc"},   {"auml", "\xc3\xa4"},
    {"Eacute", "\xc3\x89"}, {"Egrave", "\xc3\x88"}, {"Ecirc", "\xc3\x8a"},  {"Agrave", "\xc3\x80"},
    {"Acirc", "\xc3\x82"},  {"Ccedil", "\xc3\x87"}, {"Icirc", "\xc3\x8e"},  {"Ocirc", "\xc3\x94"},
    {"laquo", "\xc2\xab"},  {"raquo", "\xc2\xbb"},  {"oelig", "\xc5\x93"},  {"OElig", "\xc5\x92"},
    {"copy", "\xc2\xa9"},   {"reg", "\xc2\xae"},    {"deg", "\xc2\xb0"},    {"ordm", "\xc2\xba"},
    {"sect", "\xc2\xa7"},   {"para", "\xc2\xb6"}};

// The code point cp as UTF-8 appended to o; nothing for NUL, a surrogate
// or a value past U+10FFFF.
void utf8(unsigned long cp, string &o) {
  if (cp == 0 || cp > 0x10FFFF || (cp >= 0xD800 && cp <= 0xDFFF)) return;
  if (cp < 0x80) {
    o += static_cast<char>(cp);
  } else if (cp < 0x800) {
    o += static_cast<char>(0xC0 | (cp >> 6));
    o += static_cast<char>(0x80 | (cp & 0x3F));
  } else if (cp < 0x10000) {
    o += static_cast<char>(0xE0 | (cp >> 12));
    o += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
    o += static_cast<char>(0x80 | (cp & 0x3F));
  } else {
    o += static_cast<char>(0xF0 | (cp >> 18));
    o += static_cast<char>(0x80 | ((cp >> 12) & 0x3F));
    o += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
    o += static_cast<char>(0x80 | (cp & 0x3F));
  }
}

// The entity at s[i] ('&') decoded onto o; its length, 0 when none starts
// there. A numeric one stands for its code point (curly quotes and dashes
// as their ASCII forms, as the named ones).
size_t entity(const string &s, size_t i, string &o) {
  size_t k = i + 1;
  if (k < s.size() && s[k] == '#') {
    ++k;
    const bool hex = k < s.size() && (s[k] == 'x' || s[k] == 'X');
    if (hex) ++k;
    unsigned long cp = 0;
    const size_t d0 = k;
    while (k < s.size() && k - d0 < 8) {
      const unsigned char c = static_cast<unsigned char>(s[k]);
      int v = -1;
      if (isdig(c)) v = c - '0';
      else if (hex && low(c) >= 'a' && low(c) <= 'f') v = low(c) - 'a' + 10;
      if (v < 0) break;
      cp = cp * (hex ? 16 : 10) + static_cast<unsigned long>(v);
      ++k;
    }
    if (k == d0 || k >= s.size() || s[k] != ';') return 0;
    if (cp == 0x2018 || cp == 0x2019) o += '\'';
    else if (cp == 0x201C || cp == 0x201D) o += '"';
    else if (cp == 0x2013) o += '-';
    else if (cp == 0x2014) o += "--";
    else if (cp == 0xA0) o += ' ';
    else utf8(cp, o);
    return k + 1 - i;
  }
  for (const auto &e : kEnt) {
    const size_t n = std::strlen(e[0]);
    if (s.compare(k, n, e[0]) == 0 && k + n < s.size() && s[k + n] == ';') {
      o += e[1];
      return n + 2;
    }
  }
  return 0;
}

// Spaces collapsed, lines trimmed, at most one empty line in a row, none
// at either end.
string tidy(const string &t) {
  string o;
  o.reserve(t.size());
  int nl = 0;
  bool sp = false;
  for (size_t i = 0; i < t.size(); ++i) {
    const unsigned char c = static_cast<unsigned char>(t[i]);
    if (c == '\n') {
      ++nl;
      sp = false;
      continue;
    }
    if (c == ' ' || c == '\t' || c == '\r' || c == '\f' || c == '\v') {
      sp = true;
      continue;
    }
    if (!o.empty()) {
      if (nl) o.append(nl > 1 ? 2 : 1, '\n');
      else if (sp) o += ' ';
    }
    nl = 0;
    sp = false;
    o += static_cast<char>(c);
  }
  return o;
}

// HTML to text: comments, scripts and styles dropped, every tag a space
// except the starts and ends of block elements, which end a line; the
// page's own line breaks are spaces, as a browser shows them; entities
// and no-break spaces decoded.
string totext(const string &s) {
  string o;
  o.reserve(s.size());
  size_t i = 0;
  bool nogt = false;  // a '<' with no '>' after it: no later one has one either
  while (i < s.size()) {
    const unsigned char c = static_cast<unsigned char>(s[i]);
    if (c == '<' && !nogt) {
      if (s.compare(i, 4, "<!--") == 0) {
        const size_t e = s.find("-->", i + 4);
        o += ' ';
        i = e == string::npos ? s.size() : e + 3;
        continue;
      }
      const size_t gt = s.find('>', i + 1);
      if (gt == string::npos) {
        nogt = true;
        o += static_cast<char>(c);
        ++i;
        continue;
      }
      const bool close = i + 1 < s.size() && s[i + 1] == '/';
      const string nm = tagname(s, i + (close ? 2 : 1));
      if (nm.empty() && !(i + 1 < s.size() && (s[i + 1] == '!' || s[i + 1] == '?' || close))) {
        o += static_cast<char>(c);  // a '<' that opens no tag
        ++i;
        continue;
      }
      if (!close && (nm == "script" || nm == "style")) {
        const string end = "</" + nm;
        size_t e = gt + 1;
        while (e < s.size() && !(s[e] == '<' && ieq(s, e, end.c_str()))) ++e;
        const size_t g2 = e < s.size() ? s.find('>', e) : string::npos;
        o += ' ';
        i = g2 == string::npos ? s.size() : g2 + 1;
        continue;
      }
      bool blk = false;
      for (const char *b : kBlock) blk = blk || nm == b;
      o += blk ? '\n' : ' ';
      i = gt + 1;
      continue;
    }
    if (c == '&') {
      const size_t k = entity(s, i, o);
      if (k) {
        i += k;
        continue;
      }
    }
    if (c == 0xC2 && i + 1 < s.size() && static_cast<unsigned char>(s[i + 1]) == 0xA0) {
      o += ' ';
      i += 2;
      continue;
    }
    o += isws(c) ? ' ' : static_cast<char>(c);
    ++i;
  }
  return tidy(o);
}

// Plain text as the rules read it: line ends of any system as "\n",
// no-break spaces as spaces, and tidied as above.
string toplain(const string &s) {
  string o;
  o.reserve(s.size());
  for (size_t i = 0; i < s.size(); ++i) {
    const unsigned char c = static_cast<unsigned char>(s[i]);
    if (c == '\r') {
      o += '\n';
      if (i + 1 < s.size() && s[i + 1] == '\n') ++i;
    } else if (c == 0xC2 && i + 1 < s.size() && static_cast<unsigned char>(s[i + 1]) == 0xA0) {
      o += ' ';
      ++i;
    } else {
      o += static_cast<char>(c);
    }
  }
  return tidy(o);
}

// True when s holds markup: a '<' followed by a letter, '/' or '!'.
bool ishtml(const string &s) {
  for (size_t i = 0; i + 1 < s.size(); ++i) {
    const unsigned char c = static_cast<unsigned char>(s[i + 1]);
    if (s[i] == '<' && (isup(c) || islow(c) || c == '/' || c == '!')) return true;
  }
  return false;
}

}  // namespace

extern "C" int clab_ovrtxt(const char *s, int64_t n, int html, char *out, int64_t *len) CLAB_TRY {
  if (n < 0) return CLAB_EARG;
  const string in(s, static_cast<size_t>(n));
  const bool h = html < 0 ? ishtml(in) : html != 0;
  const string t = h ? totext(in) : toplain(in);
  std::memcpy(out, t.data(), t.size());
  *len = static_cast<int64_t>(t.size());
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)
