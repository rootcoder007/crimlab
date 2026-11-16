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

namespace {

// A token of text: a run of letters (w), of digits (d) or one other
// character (p), at s[a, e).
struct Tok {
  char k;
  size_t a, e;
};

std::vector<Tok> toks(const string &s, size_t a, size_t e) {
  std::vector<Tok> t;
  size_t i = a;
  while (i < e) {
    const unsigned char c = static_cast<unsigned char>(s[i]);
    if (isws(c)) {
      ++i;
      continue;
    }
    size_t j = i + 1;
    char k = 'p';
    if (isdig(c)) {
      k = 'd';
      while (j < e && isdig(static_cast<unsigned char>(s[j]))) ++j;
    } else if (iswrd(c)) {
      k = 'w';
      while (j < e && iswrd(static_cast<unsigned char>(s[j])) && !isdig(static_cast<unsigned char>(s[j]))) ++j;
    }
    t.push_back({k, i, j});
    i = j;
  }
  return t;
}

// A word in lower case: ASCII letters and the capitals of Latin-1 (as
// UTF-8) folded.
string fold(const string &s, const Tok &t) {
  string w = s.substr(t.a, t.e - t.a);
  for (size_t i = 0; i < w.size(); ++i) {
    const unsigned char c = static_cast<unsigned char>(w[i]);
    if (isup(c)) w[i] = static_cast<char>(c + 32);
    if (c == 0xC3 && i + 1 < w.size()) {
      const unsigned char d = static_cast<unsigned char>(w[i + 1]);
      if (d >= 0x80 && d <= 0x9E && d != 0x97) w[i + 1] = static_cast<char>(d + 0x20);
    }
  }
  return w;
}

// The month of a word (English or French, whole or abbreviated), 0 for
// none.
int month(const string &w) {
  static const char *const nm[] = {
      "january", "february", "march", "april", "may", "june", "july", "august", "september", "october",
      "november", "december", "jan", "feb", "mar", "apr", "jun", "jul", "aug", "sep", "sept", "oct", "nov", "dec",
      "janvier", "f\xc3\xa9vrier", "fevrier", "mars", "avril", "mai", "juin", "juillet", "ao\xc3\xbbt", "aout",
      "septembre", "octobre", "novembre", "d\xc3\xa9" "cembre", "decembre", "janv", "f\xc3\xa9vr", "fevr", "avr",
      "juil", "d\xc3\xa9" "c"};
  static const int mo[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 6, 7, 8, 9, 9, 10, 11, 12,
                           1, 2, 2, 3, 4, 5, 6, 7, 8, 8, 9, 10, 11, 12, 12, 1, 2, 2, 4, 7, 12};
  for (size_t i = 0; i < sizeof mo / sizeof mo[0]; ++i) {
    if (w == nm[i]) return mo[i];
  }
  return 0;
}

int num(const string &s, const Tok &t) {
  int v = 0;
  for (size_t i = t.a; i < t.e && i < t.a + 9; ++i) v = v * 10 + (s[i] - '0');
  return v;
}

bool realday(int y, int m, int d) {
  static const int dm[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  if (y < 1000 || y > 9999 || m < 1 || m > 12 || d < 1) return false;
  const bool leap = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
  return d <= dm[m - 1] + (m == 2 && leap ? 1 : 0);
}

bool ordsfx(const string &w) {
  return w == "st" || w == "nd" || w == "rd" || w == "th" || w == "er" || w == "e";
}

// The first date written in s[a, e) as YYYY-MM-DD into iso: "January 5,
// 2023", "Jan. 5th 2023", "5 January 2023", "5th of January, 2023", "1er
// janvier 2023", or 2023-01-05 itself. A day that does not exist is
// passed over.
bool finddate(const string &s, size_t a, size_t e, string &iso) {
  const std::vector<Tok> t = toks(s, a, e);
  const size_t n = t.size();
  auto isp = [&](size_t i, char ch) { return i < n && t[i].k == 'p' && s[t[i].a] == ch; };
  auto day = [&](size_t i) { return i < n && t[i].k == 'd' && t[i].e - t[i].a <= 2; };
  auto year = [&](size_t i) { return i < n && t[i].k == 'd' && t[i].e - t[i].a == 4; };
  for (size_t i = 0; i < n; ++i) {
    int y = 0, m = 0, d = 0;
    size_t j = i;
    if (t[i].k == 'w' && (m = month(fold(s, t[i]))) != 0) {
      j = i + 1;
      if (isp(j, '.')) ++j;
      if (!day(j)) continue;
      d = num(s, t[j++]);
      if (j < n && t[j].k == 'w' && ordsfx(fold(s, t[j]))) ++j;
      if (isp(j, ',')) ++j;
      if (!year(j)) continue;
      y = num(s, t[j]);
    } else if (day(i)) {
      d = num(s, t[i]);
      j = i + 1;
      if (j < n && t[j].k == 'w' && ordsfx(fold(s, t[j]))) ++j;
      if (j < n && t[j].k == 'w' && fold(s, t[j]) == "of") ++j;
      if (!(j < n && t[j].k == 'w' && (m = month(fold(s, t[j]))) != 0)) continue;
      ++j;
      if (isp(j, '.')) ++j;
      if (isp(j, ',')) ++j;
      if (!year(j)) continue;
      y = num(s, t[j]);
    } else if (year(i) && i + 4 < n && isp(i + 1, '-') && isp(i + 3, '-') && t[i + 2].k == 'd' &&
               t[i + 4].k == 'd' && t[i + 2].e - t[i + 2].a == 2 && t[i + 4].e - t[i + 4].a == 2 &&
               t[i + 1].a == t[i].e && t[i + 2].a == t[i + 1].e && t[i + 3].a == t[i + 2].e &&
               t[i + 4].a == t[i + 3].e) {
      y = num(s, t[i]);
      m = num(s, t[i + 2]);
      d = num(s, t[i + 4]);
    } else {
      continue;
    }
    if (!realday(y, m, d)) continue;
    char b[16];
    std::snprintf(b, sizeof b, "%04d-%02d-%02d", y, m, d);
    iso = b;
    return true;
  }
  return false;
}

}  // namespace

extern "C" int clab_ovrdte(const char *s, int64_t n, char *out) CLAB_TRY {
  if (n < 0) return CLAB_EARG;
  const string in(s, static_cast<size_t>(n));
  string iso;
  out[0] = '\0';
  if (finddate(in, 0, in.size(), iso)) std::memcpy(out, iso.c_str(), 11);
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

namespace {

// A field of a preset schema: its name, the rule that reads it, the cues
// (alternatives separated by '|'), the cues of the opposite answer
// (rule flag), the heading of the section it is read in ("" for the whole
// report) and what it records.
struct Fld {
  const char *name, *rule, *cue, *alt, *sect, *desc;
};

// The director's reports of Ontario's Special Investigations Unit, as
// published in English since 2020.
const Fld kSiu[] = {
    {"case_number", "after", "Case #", "", "", "the case number of the report"},
    {"police_service", "name", "Police Service|Provincial Police|Police", "", "",
     "the police service whose officials were investigated"},
    {"date_of_incident_iso", "date", "", "", "Incident Narrative", "the date of the incident"},
    {"date_siu_notified_iso", "date", "notified the SIU|contacted the SIU|SIU was notified|SIU was contacted", "",
     "", "the date the unit was notified"},
    {"date_of_director_decision_iso", "date", "Date:", "", "", "the date of the director's decision"},
    {"siu_investigators", "num", "Number of SIU Investigators assigned", "", "",
     "the number of investigators assigned"},
    {"siu_forensics_investigators", "num", "Number of SIU Forensic Investigators assigned", "", "",
     "the number of forensic investigators assigned"},
    {"number_of_witness_officials", "tag", "WO", "", "", "the number of witness officials (WO #1, ...)"},
    {"number_of_civilian_witnesses", "tag", "CW", "", "", "the number of civilian witnesses (CW #1, ...)"},
    {"number_of_subject_officials", "tag", "SO", "", "", "the number of subject officials (SO #1, ...)"},
    {"age_affected", "before", "year-old|year old", "", "", "the age of the affected person"},
    {"sex_gender_affected", "word", "year-old|year old", "", "", "the sex or gender of the affected person"},
    {"charges_recommended", "flag", "charged with|charges have been laid|criminal charges were laid",
     "no reasonable grounds|No charges|no charges|not be laid|decline to lay",
     "Analysis and Director's Decision|Analysis and Director\xe2\x80\x99s Decision",
     "whether the director found grounds for charges"},
    {"directors_name", "prev", "Director|Acting Director|Interim Director", "", "",
     "the name of the director who decided"},
    {"location_of_call", "after",
     "in the City of|in the Town of|in the Township of|in the Municipality of|in the Region of", "", "",
     "the municipality where the incident happened"},
    {"specific_injuries", "sent",
     "fracture|fractures|fractured|broken|laceration|lacerations|gunshot|stab|stabbed|burn|burns", "",
     "Nature of Injuries / Treatment|Nature of Injuries/Treatment", "the injuries the affected person had"},
    {"relevant_legislation", "lines", "Section", "", "Relevant Legislation",
     "the provisions the director considered"},
};

}  // namespace

extern "C" int clab_ovrsch(const char *preset, int j, const char **name, const char **rule, const char **cue,
                           const char **alt, const char **sect, const char **desc) {
  if (std::strcmp(preset, "siu") != 0 || j < 0 || j >= static_cast<int>(sizeof kSiu / sizeof kSiu[0])) {
    return CLAB_EARG;
  }
  const Fld &f = kSiu[j];
  *name = f.name;
  *rule = f.rule;
  *cue = f.cue;
  *alt = f.alt;
  *sect = f.sect;
  *desc = f.desc;
  return CLAB_OK;
}

namespace {

// s without its leading and trailing spaces
string trim(const string &s) {
  size_t a = 0, b = s.size();
  while (a < b && isws(static_cast<unsigned char>(s[a]))) ++a;
  while (b > a && isws(static_cast<unsigned char>(s[b - 1]))) --b;
  return s.substr(a, b - a);
}

// The alternatives of a cue ("a|b"), empty ones left out.
std::vector<string> alts(const char *cue) {
  std::vector<string> v;
  string w;
  for (const char *p = cue;; ++p) {
    if (*p == '|' || *p == '\0') {
      if (!w.empty()) v.push_back(w);
      w.clear();
      if (*p == '\0') break;
    } else {
      w += *p;
    }
  }
  return v;
}

// The length of the alternative of al that occurs at s[i] as whole words
// at its ends and inside [i, e), 0 when none does.
size_t at(const string &s, size_t i, size_t e, const std::vector<string> &al) {
  for (const string &w : al) {
    if (i + w.size() > e || s.compare(i, w.size(), w) != 0) continue;
    if (iswrd(static_cast<unsigned char>(w[0])) && i > 0 && iswrd(static_cast<unsigned char>(s[i - 1]))) continue;
    const size_t j = i + w.size();
    if (iswrd(static_cast<unsigned char>(w.back())) && j < s.size() && iswrd(static_cast<unsigned char>(s[j]))) {
      continue;
    }
    return w.size();
  }
  return 0;
}

// The line of s that holds position i: [*a, *b).
void line(const string &s, size_t i, size_t *a, size_t *b) {
  const size_t p = i == 0 ? string::npos : s.rfind('\n', i - 1);
  *a = p == string::npos ? 0 : p + 1;
  const size_t q = s.find('\n', i);
  *b = q == string::npos ? s.size() : q;
}

// A line that looks like a heading: one to eight words, the first
// starting with a capital, no digit, not ending in a stop or a comma.
// ponytail: a heuristic; a section of a layout it misjudges is cut short
// or runs on, and a schema can avoid sections for such fields.
bool heading(const string &l) {
  if (l.empty()) return false;
  const unsigned char c0 = static_cast<unsigned char>(l[0]);
  if (!(isup(c0) || c0 >= 0xC3)) return false;
  int words = 1;
  for (size_t i = 0; i < l.size(); ++i) {
    const unsigned char c = static_cast<unsigned char>(l[i]);
    if (isdig(c)) return false;
    if (c == ' ') ++words;
  }
  const char z = l.back();
  return words <= 8 && z != '.' && z != ',' && z != ';' && z != ':' && z != '!' && z != '?';
}

// The section under the last line that is one of the headings of sect, up
// to the next line that looks like a heading: [*a, *e). The whole text
// when sect is empty; false when no line is one of its headings.
bool section(const string &s, const char *sect, size_t *a, size_t *e) {
  const std::vector<string> hs = alts(sect);
  *a = 0;
  *e = s.size();
  if (hs.empty()) return true;
  bool hit = false;
  for (size_t i = 0; i < s.size();) {
    size_t la, lb;
    line(s, i, &la, &lb);
    const string l = trim(s.substr(la, lb - la));
    for (const string &h : hs) {
      if (l == h) {
        hit = true;
        *a = lb;
      }
    }
    i = lb + 1;
  }
  if (!hit) return false;
  for (size_t i = *a + 1; i < s.size();) {
    size_t la, lb;
    line(s, i, &la, &lb);
    if (heading(trim(s.substr(la, lb - la)))) {
      *e = la;
      break;
    }
    i = lb + 1;
  }
  return true;
}

// The sentence of s[a, e) that holds position p: from the last line end
// or sentence end (".", "!" or "?", a space and a capital or a digit)
// before p to the next one.
string sentence(const string &s, size_t p, size_t a, size_t e) {
  auto end = [&](size_t k) {
    return (s[k] == '.' || s[k] == '!' || s[k] == '?') && k + 2 < e && s[k + 1] == ' ' &&
           (isup(static_cast<unsigned char>(s[k + 2])) || isdig(static_cast<unsigned char>(s[k + 2])));
  };
  size_t b = a;
  for (size_t k = p; k > a; --k) {
    if (s[k - 1] == '\n') {
      b = k;
      break;
    }
    if (k + 1 <= p && end(k - 1)) {
      b = k + 1;
      break;
    }
  }
  size_t z = e;
  for (size_t k = p; k < e; ++k) {
    if (s[k] == '\n') {
      z = k;
      break;
    }
    if (end(k)) {
      z = k + 1;
      break;
    }
  }
  return trim(s.substr(b, z - b));
}

// A number at s[i]: digits, or a number word in English or French up to
// twenty; its value and end.
bool readnum(const string &s, size_t i, size_t e, long *v, size_t *end) {
  if (i >= e) return false;
  if (isdig(static_cast<unsigned char>(s[i]))) {
    size_t j = i;
    long x = 0;
    while (j < e && isdig(static_cast<unsigned char>(s[j])) && j - i < 9) x = x * 10 + (s[j++] - '0');
    if (j < e && isdig(static_cast<unsigned char>(s[j]))) return false;
    *v = x;
    *end = j;
    return true;
  }
  if (i > 0 && iswrd(static_cast<unsigned char>(s[i - 1]))) return false;
  size_t j = i;
  while (j < e && iswrd(static_cast<unsigned char>(s[j])) && !isdig(static_cast<unsigned char>(s[j]))) ++j;
  if (j == i) return false;
  string w = s.substr(i, j - i);
  for (char &c : w) c = static_cast<char>(low(static_cast<unsigned char>(c)));
  static const char *const en[] = {"zero",    "one",     "two",       "three",    "four",    "five",    "six",
                                   "seven",   "eight",   "nine",      "ten",      "eleven",  "twelve",  "thirteen",
                                   "fourteen", "fifteen", "sixteen",  "seventeen", "eighteen", "nineteen", "twenty"};
  static const char *const fr[] = {"z\xc3\xa9ro", "un",    "deux",     "trois",  "quatre", "cinq",   "six",
                                   "sept",        "huit",  "neuf",     "dix",    "onze",   "douze",  "treize",
                                   "quatorze",    "quinze", "seize",   "dix-sept", "dix-huit", "dix-neuf", "vingt"};
  for (int k = 0; k <= 20; ++k) {
    if (w == en[k] || w == fr[k] || (k == 1 && w == "une")) {
      *v = k;
      *end = j;
      return true;
    }
  }
  return false;
}

// The first position from i in [i, e) where an alternative occurs, and
// its length; npos when none.
size_t next(const string &s, size_t i, size_t e, const std::vector<string> &al, size_t *len) {
  for (; i < e; ++i) {
    const size_t n = at(s, i, e, al);
    if (n) {
      *len = n;
      return i;
    }
  }
  return string::npos;
}

// A word of a name ending at s[k - 1]: its start, or npos when the word
// does not start with a capital.
size_t capword(const string &s, size_t lo, size_t k) {
  size_t j = k;
  while (j > lo) {
    const unsigned char c = static_cast<unsigned char>(s[j - 1]);
    if (!(iswrd(c) || c == '.' || c == '\'' || c == '-')) break;
    --j;
  }
  if (j == k) return string::npos;
  const unsigned char c = static_cast<unsigned char>(s[j]);
  const bool cap = isup(c) || (c == 0xC3 && j + 1 < k && static_cast<unsigned char>(s[j + 1]) <= 0x9E);
  return cap ? j : string::npos;
}

// Rules, in the order of their codes (R/ovr.R): after, num, date, tag,
// sent, prev, lines, name, flag, before, word.
bool field(const string &s, int rule, const char *cue, const char *alt, const char *sect, string &v) {
  size_t a, e, len = 0;
  if (!section(s, sect, &a, &e)) return false;
  const std::vector<string> al = alts(cue), op = alts(alt);
  switch (rule) {
    case 0:
      for (size_t p = next(s, a, e, al, &len); p != string::npos; p = next(s, p + 1, e, al, &len)) {
        size_t q = p + len;
        while (q < e && (s[q] == ' ' || s[q] == ':')) ++q;
        size_t z = q;
        while (z < e && s[z] != ',' && s[z] != ';' && s[z] != '\n' &&
               !(s[z] == '.' && (z + 1 == e || s[z + 1] == ' ' || s[z + 1] == '\n'))) {
          ++z;
        }
        v = trim(s.substr(q, z - q));
        if (!v.empty()) return true;
      }
      return false;
    case 1:
      for (size_t p = next(s, a, e, al, &len); p != string::npos; p = next(s, p + 1, e, al, &len)) {
        size_t la, lb;
        line(s, p, &la, &lb);
        const size_t z = std::min({lb, e, p + len + 40});
        for (size_t k = p + len; k < z; ++k) {
          long x;
          size_t end;
          if (readnum(s, k, std::min(lb, e), &x, &end)) {
            v = std::to_string(x);
            return true;
          }
        }
      }
      return false;
    case 2:
      if (al.empty()) return finddate(s, a, e, v);
      for (size_t p = next(s, a, e, al, &len); p != string::npos; p = next(s, p + 1, e, al, &len)) {
        const string t = sentence(s, p, a, e);
        if (finddate(t, 0, t.size(), v)) return true;
      }
      return false;
    case 3: {
      long mx = 0;
      for (size_t p = next(s, a, e, al, &len); p != string::npos; p = next(s, p + 1, e, al, &len)) {
        size_t k = p + len;
        while (k < e && isws(static_cast<unsigned char>(s[k]))) ++k;
        if (k < e && s[k] == '#') {
          ++k;
        } else if (k + 1 < e && (s[k] == 'n' || s[k] == 'N') && s[k + 1] == 'o' &&
                   (k + 2 >= e || !iswrd(static_cast<unsigned char>(s[k + 2])))) {
          k += 2;
          if (k < e && s[k] == '.') ++k;
        } else if (k + 2 < e && s[k] == 'n' && static_cast<unsigned char>(s[k + 1]) == 0xC2 &&
                   (static_cast<unsigned char>(s[k + 2]) == 0xB0 || static_cast<unsigned char>(s[k + 2]) == 0xBA)) {
          k += 3;
        }
        while (k < e && isws(static_cast<unsigned char>(s[k]))) ++k;
        long x;
        size_t end;
        if (k < e && isdig(static_cast<unsigned char>(s[k])) && readnum(s, k, e, &x, &end)) mx = std::max(mx, x);
      }
      if (mx < 1) return false;
      v = std::to_string(mx);
      return true;
    }
    case 4: {
      const size_t p = next(s, a, e, al, &len);
      if (p == string::npos) return false;
      v = sentence(s, p, a, e);
      return true;
    }
    case 5: {
      size_t hit = string::npos;
      for (size_t i = a; i < e;) {
        size_t la, lb;
        line(s, i, &la, &lb);
        const string l = trim(s.substr(la, lb - la));
        for (const string &w : al) {
          if (l == w) hit = la;
        }
        i = lb + 1;
      }
      if (hit == string::npos) return false;
      for (size_t i = hit; i > a;) {
        size_t la, lb;
        line(s, i - 1, &la, &lb);
        const string l = trim(s.substr(la, lb - la));
        if (!l.empty()) {
          v = l;
          return true;
        }
        i = la;
      }
      return false;
    }
    case 6: {
      std::vector<string> got;
      for (size_t i = a; i < e;) {
        size_t la, lb;
        line(s, i, &la, &lb);
        const string l = trim(s.substr(la, std::min(lb, e) - la));
        if (!l.empty() && at(l, 0, l.size(), al) && std::find(got.begin(), got.end(), l) == got.end()) {
          got.push_back(l);
        }
        i = lb + 1;
      }
      if (got.empty()) return false;
      v.clear();
      for (size_t k = 0; k < got.size(); ++k) v += (k ? "; " : "") + got[k];
      return true;
    }
    case 7:
      for (size_t p = next(s, a, e, al, &len); p != string::npos; p = next(s, p + 1, e, al, &len)) {
        // the capitalised words before the cue, nearest first
        std::vector<size_t> ws;
        size_t k = p;
        for (int w = 0; w < 4 && k > a + 1 && s[k - 1] == ' '; ++w) {
          const size_t j = capword(s, a, k - 1);
          if (j == string::npos) break;
          ws.push_back(j);
          k = j;
        }
        // without a leading article
        while (!ws.empty()) {
          const size_t st = ws.back(), we = (ws.size() > 1 ? ws[ws.size() - 2] : p) - 1;
          const string w = s.substr(st, we - st);
          if (w != "The" && w != "A" && w != "An") break;
          ws.pop_back();
        }
        if (ws.empty()) continue;
        v = s.substr(ws.back(), p + len - ws.back());
        return true;
      }
      return false;
    case 8: {
      size_t n2 = 0;
      if (!op.empty() && next(s, a, e, op, &n2) != string::npos) {
        v = "FALSE";
        return true;
      }
      if (next(s, a, e, al, &len) == string::npos) return false;
      v = "TRUE";
      return true;
    }
    case 9:
      for (size_t p = next(s, a, e, al, &len); p != string::npos; p = next(s, p + 1, e, al, &len)) {
        size_t k = p;
        while (k > a && (s[k - 1] == ' ' || s[k - 1] == '-')) --k;
        size_t j = k;
        while (j > a && iswrd(static_cast<unsigned char>(s[j - 1]))) --j;
        long x;
        size_t end;
        if (j < k && readnum(s, j, k, &x, &end) && end == k) {
          v = std::to_string(x);
          return true;
        }
      }
      return false;
    case 10:
      for (size_t p = next(s, a, e, al, &len); p != string::npos; p = next(s, p + 1, e, al, &len)) {
        size_t q = p + len;
        while (q < e && (s[q] == ' ' || s[q] == '-')) ++q;
        size_t z = q;
        while (z < e && iswrd(static_cast<unsigned char>(s[z])) && !isdig(static_cast<unsigned char>(s[z]))) ++z;
        if (z > q) {
          v = s.substr(q, z - q);
          return true;
        }
      }
      return false;
    default:
      return false;
  }
}

}  // namespace

extern "C" int clab_ovrfld(const char *t, int64_t n, int rule, const char *cue, const char *alt, const char *sect,
                           char *out, int64_t *len) CLAB_TRY {
  if (n < 0 || rule < 0 || rule > 10) return CLAB_EARG;
  const string s(t, static_cast<size_t>(n));
  string v;
  *len = -1;
  if (!field(s, rule, cue, alt, sect, v)) return CLAB_OK;
  std::memcpy(out, v.data(), v.size());
  *len = static_cast<int64_t>(v.size());
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)
