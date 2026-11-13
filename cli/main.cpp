// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// crimlab: the functions of the package from the command line. The
// arguments come as a JSON object or a CSV table (a file, or standard
// input) and as --name value options; the result is printed as a table,
// or with -j as JSON. Everything else is clab_jrun() (jrun.cpp).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

#include "jrun.h"
#include "json.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#ifndef CLAB_VER
#define CLAB_VER "unknown"
#endif

using clabj::Val;

// The WebAssembly build (bindings/npm) exports the program as clab_cmain.
#ifdef CLAB_CMAIN
extern "C" int clab_cmain(int argc, char **argv);
#define main clab_cmain
#endif

namespace {

const char *const usage =
    "usage: crimlab FN [-j] [-f FILE] [--ARG VALUE ...]\n"
    "       crimlab -l | -h | -V\n"
    "\n"
    "Runs the crimlab function FN (the name of the R function, such as\n"
    "orrisk) and prints its result as a table, or as JSON with -j.\n"
    "\n"
    "Its arguments come from FILE (- for standard input), a JSON object of\n"
    "named arguments or a CSV table whose columns are arguments, and from\n"
    "--ARG VALUE options, which win. Without -f, the arguments the options\n"
    "leave out are read from standard input. A VALUE\n"
    "is JSON when it starts with [ { or \", a list when it holds a comma\n"
    "(1,2,3), else one number, true, false, NA or text. A CSV column named\n"
    "group names the rows: the other columns become per-group values.\n"
    "\n"
    "  -l  the functions and their arguments ([optional])\n"
    "  -h  this help\n"
    "  -V  the version\n"
    "\n"
    "example: crimlab orrisk --or 2 --rate 0.1 --share 0.3\n"
    "         echo '{\"done\": [30, 45, 60], \"open\": [100, 150]}' | crimlab losbnd\n";

void put(Val &o, const std::string &k, Val v) {
  for (std::size_t i = 0; i < o.k.size(); ++i) {
    if (o.k[i] == k) {
      o.a[i] = std::move(v);
      return;
    }
  }
  o.k.push_back(k);
  o.a.push_back(std::move(v));
}

std::string fmt(const Val &v) {
  switch (v.t) {
    case Val::NUL: return "NA";
    case Val::BOOL: return v.b ? "TRUE" : "FALSE";
    case Val::STR: return v.s;
    case Val::NUM: {
      char b[32];
      std::snprintf(b, sizeof b, "%.7g", v.n);
      return b;
    }
    default: return "";
  }
}

// One CSV cell or option value: a number, a bool, NA (empty or "NA"),
// else text.
Val cell(const std::string &s) {
  Val v;
  if (s.empty() || s == "NA") return v;
  if (s == "true" || s == "TRUE" || s == "false" || s == "FALSE") {
    v.t = Val::BOOL;
    v.b = s[0] == 't' || s[0] == 'T';
    return v;
  }
  const char c = s[0];
  if ((c >= '0' && c <= '9') || c == '-' || c == '+' || c == '.') {
    char *e = nullptr;
    const double d = std::strtod(s.c_str(), &e);
    if (*e == '\0' && std::isfinite(d)) {
      v.t = Val::NUM;
      v.n = d;
      return v;
    }
  }
  v.t = Val::STR;
  v.s = s;
  return v;
}

// An option's value.
bool optval(const std::string &s, Val &v) {
  if (!s.empty() && (s[0] == '[' || s[0] == '{' || s[0] == '"')) {
    return clabj::Parser(s.data(), s.size()).doc(v);
  }
  if (s.find(',') == std::string::npos) {
    v = cell(s);
    return true;
  }
  v = Val();
  v.t = Val::ARR;
  std::size_t b = 0;
  for (;;) {
    const std::size_t e = s.find(',', b);
    v.a.push_back(cell(s.substr(b, e == std::string::npos ? std::string::npos : e - b)));
    if (e == std::string::npos) return true;
    b = e + 1;
  }
}

// RFC 4180 CSV: rows of cells; a quoted cell is text.
bool csv(const std::string &s, std::vector<std::vector<Val>> &rows) {
  std::vector<Val> row;
  std::string f;
  bool quoted = false, any = false;
  auto endcell = [&]() {
    Val v = quoted ? Val() : cell(f);
    if (quoted) {
      v.t = Val::STR;
      v.s = f;
    }
    row.push_back(v);
    f.clear();
    quoted = false;
  };
  for (std::size_t i = 0; i < s.size(); ++i) {
    const char c = s[i];
    if (c == '"' && f.empty() && !quoted) {
      quoted = true;
      std::size_t j = i + 1;
      for (;; ++j) {
        if (j >= s.size()) return false;
        if (s[j] == '"') {
          if (j + 1 < s.size() && s[j + 1] == '"') {
            f += '"';
            ++j;
          } else {
            break;
          }
        } else {
          f += s[j];
        }
      }
      i = j;
      any = true;
    } else if (c == ',') {
      endcell();
      any = true;
    } else if (c == '\n' || c == '\r') {
      if (c == '\r' && i + 1 < s.size() && s[i + 1] == '\n') ++i;
      if (any || !f.empty()) {
        endcell();
        rows.push_back(row);
      }
      row.clear();
      any = false;
    } else {
      if (quoted) return false;
      f += c;
      any = true;
    }
  }
  if (any || !f.empty()) {
    endcell();
    rows.push_back(row);
  }
  return true;
}

// A CSV table as named arguments: one column each, trailing NA cells
// dropped, a column of one value a single value. A column named group
// names the rows, and the other columns become objects from group to value.
bool table(const std::string &s, Val &o, std::string &err) {
  std::vector<std::vector<Val>> rows;
  if (!csv(s, rows) || rows.empty()) {
    err = "the CSV input cannot be read.";
    return false;
  }
  const std::vector<Val> &h = rows[0];
  std::size_t gi = h.size();
  for (std::size_t j = 0; j < h.size(); ++j) {
    if (h[j].t != Val::STR || h[j].s.empty()) {
      err = "the first CSV row must name every column.";
      return false;
    }
    if (h[j].s == "group") gi = j;
  }
  for (std::size_t r = 1; r < rows.size(); ++r) {
    if (rows[r].size() != h.size()) {
      err = "CSV row " + std::to_string(r + 1) + " has " + std::to_string(rows[r].size()) +
            " cells, not " + std::to_string(h.size()) + ".";
      return false;
    }
  }
  for (std::size_t j = 0; j < h.size(); ++j) {
    if (j == gi) continue;
    Val col;
    col.t = gi < h.size() ? Val::OBJ : Val::ARR;
    for (std::size_t r = 1; r < rows.size(); ++r) {
      if (gi < h.size()) {
        const Val &n = rows[r][gi];
        col.k.push_back(n.t == Val::NUM ? fmt(n) : n.s);
      }
      col.a.push_back(rows[r][j]);
    }
    if (gi == h.size()) {
      while (!col.a.empty() && col.a.back().t == Val::NUL) col.a.pop_back();
      if (col.a.size() == 1) col = Val(col.a[0]);
    }
    put(o, h[j].s, col);
  }
  return true;
}

bool slurp(std::FILE *f, std::string &s) {
  char b[1 << 16];
  std::size_t n;
  while ((n = std::fread(b, 1, sizeof b, f)) > 0) s.append(b, n);
  return !std::ferror(f);
}

// Printing a result.

bool scalar(const Val &v) { return v.t != Val::ARR && v.t != Val::OBJ; }

bool frame(const Val &v) {
  if (v.t != Val::OBJ || v.a.empty()) return false;
  for (const Val &c : v.a) {
    if (c.t != Val::ARR || c.a.size() != v.a[0].a.size()) return false;
    for (const Val &x : c.a) {
      if (!scalar(x)) return false;
    }
  }
  return true;
}

void row(const std::vector<std::string> &cells, const std::vector<std::size_t> &w,
         const std::vector<bool> &left) {
  std::string l;
  for (std::size_t j = 0; j < cells.size(); ++j) {
    if (j) l += "  ";
    const std::string pad(w[j] - cells[j].size(), ' ');
    l += left[j] ? cells[j] + pad : pad + cells[j];
  }
  while (!l.empty() && l.back() == ' ') l.pop_back();
  std::printf("%s\n", l.c_str());
}

// Columns side by side: text left, numbers right.
void cols(const std::vector<std::string> &head, const std::vector<std::vector<Val>> &c) {
  const std::size_t nc = head.size(), nr = c.empty() ? 0 : c[0].size();
  std::vector<std::size_t> w(nc);
  std::vector<bool> left(nc, true);
  for (std::size_t j = 0; j < nc; ++j) {
    w[j] = head[j].size();
    for (const Val &x : c[j]) {
      w[j] = std::max(w[j], fmt(x).size());
      if (x.t == Val::NUM) left[j] = false;
    }
  }
  if (!head.empty() && !head[0].empty()) row(head, w, left);
  for (std::size_t i = 0; i < nr; ++i) {
    std::vector<std::string> r(nc);
    for (std::size_t j = 0; j < nc; ++j) r[j] = fmt(c[j][i]);
    row(r, w, left);
  }
}

void show(const Val &v) {
  if (scalar(v)) {
    std::printf("%s\n", fmt(v).c_str());
  } else if (v.t == Val::ARR) {
    for (const Val &x : v.a) show(x);
  } else if (frame(v)) {
    std::vector<std::vector<Val>> c;
    for (const Val &x : v.a) c.push_back(x.a);
    cols(v.k, c);
  } else {
    // named values first, as name and value; then each nested part
    std::vector<Val> nm, val;
    for (std::size_t i = 0; i < v.a.size(); ++i) {
      if (!scalar(v.a[i])) continue;
      Val k;
      k.t = Val::STR;
      k.s = v.k[i];
      nm.push_back(k);
      val.push_back(v.a[i]);
    }
    if (!nm.empty()) cols({"", ""}, {nm, val});
    for (std::size_t i = 0; i < v.a.size(); ++i) {
      if (scalar(v.a[i])) continue;
      std::printf("\n%s:\n", v.k[i].c_str());
      show(v.a[i]);
    }
  }
}

// The functions and their arguments, from clab_jfns().
Val fns() {
  Val v;
  const char *s = clab_jfns();
  clabj::Parser(s, std::strlen(s)).doc(v);
  return v;
}

int fail(const std::string &m, int rc) {
  std::fprintf(stderr, "crimlab: %s\n", m.c_str());
  return rc;
}

}  // namespace

int main(int argc, char **argv) {
#ifdef _WIN32
  // the same bytes on every system: \n ends a line
  _setmode(_fileno(stdout), _O_BINARY);
#endif
  if (argc < 2) {
    std::fputs(usage, stderr);
    return 2;
  }
  const std::string fn = argv[1];
  if (fn == "-h" || fn == "--help") {
    std::fputs(usage, stdout);
    return 0;
  }
  if (fn == "-V" || fn == "--version") {
    std::printf("crimlab %s\n", CLAB_VER);
    return 0;
  }
  const Val all = fns();
  if (fn == "-l") {
    for (std::size_t i = 0; i < all.k.size(); ++i) {
      std::string l = all.k[i];
      for (const Val &a : all.a[i].a) {
        l += a.s.back() == '?' ? " [" + a.s.substr(0, a.s.size() - 1) + "]" : " " + a.s;
      }
      std::printf("%s\n", l.c_str());
    }
    return 0;
  }
  std::size_t fi = all.k.size();
  for (std::size_t i = 0; i < all.k.size(); ++i) {
    if (all.k[i] == fn) fi = i;
  }
  if (fi == all.k.size()) return fail("no function " + fn + " (crimlab -l lists them).", 2);
  bool json = false;
  const char *file = nullptr;
  Val opts;
  opts.t = Val::OBJ;
  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "-j") {
      json = true;
    } else if (a == "-f" && i + 1 < argc) {
      file = argv[++i];
    } else if (a.size() > 2 && a.compare(0, 2, "--") == 0) {
      std::string k = a.substr(2), s;
      const std::size_t eq = k.find('=');
      if (eq != std::string::npos) {
        s = k.substr(eq + 1);
        k.resize(eq);
      } else if (i + 1 < argc) {
        s = argv[++i];
      } else {
        return fail("--" + k + " needs a value.", 2);
      }
      Val v;
      if (!optval(s, v)) return fail("--" + k + ": the value is not JSON: " + s, 2);
      put(opts, k, v);
    } else {
      return fail("unknown option " + a + " (crimlab -h for help).", 2);
    }
  }
  // standard input, without -f, when the options leave a needed argument out
  bool needs = false;
  for (const Val &a : all.a[fi].a) {
    needs = needs || (a.s.back() != '?' && std::find(opts.k.begin(), opts.k.end(), a.s) == opts.k.end());
  }
  if (file == nullptr && needs) file = "-";
  Val in;
  in.t = Val::OBJ;
  if (file != nullptr) {
    std::string s;
    std::FILE *f = std::strcmp(file, "-") == 0 ? stdin : std::fopen(file, "rb");
    if (f == nullptr) return fail(std::string("cannot open ") + file + ".", 2);
    const bool ok = slurp(f, s);
    if (f != stdin) std::fclose(f);
    if (!ok) return fail(std::string("cannot read ") + file + ".", 2);
    std::size_t b = s.find_first_not_of(" \t\r\n");
    if (b != std::string::npos && s.compare(b, 3, "\xef\xbb\xbf") == 0) b = s.find_first_not_of(" \t\r\n", b + 3);
    if (b != std::string::npos && s[b] == '{') {
      if (!clabj::Parser(s.data(), s.size()).doc(in) || in.t != Val::OBJ) {
        return fail("the JSON input is not one object.", 2);
      }
    } else if (b != std::string::npos) {
      std::string err;
      if (!table(s.substr(b), in, err)) return fail(err, 2);
    }
  }
  for (std::size_t i = 0; i < opts.k.size(); ++i) put(in, opts.k[i], opts.a[i]);
  std::string q;
  clabj::wval(q, in);
  std::vector<char> out(1 << 16);
  int need = 0;
  int rc = clab_jrun(fn.c_str(), q.c_str(), out.data(), static_cast<int>(out.size()), &need);
  if (rc == CLAB_ESIZE) {
    out.resize(static_cast<std::size_t>(need));
    rc = clab_jrun(fn.c_str(), q.c_str(), out.data(), need, &need);
  }
  if (rc != CLAB_OK) {
    std::fprintf(stderr, "%s\n", out.data());
    return 1;
  }
  if (json) {
    std::printf("%s\n", out.data());
    return 0;
  }
  Val r;
  clabj::Parser(out.data(), std::strlen(out.data())).doc(r);
  show(r);
  return 0;
}
