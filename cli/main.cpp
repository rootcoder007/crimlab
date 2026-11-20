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

// The notice of a newer release (below): not in the WebAssembly build,
// which has no network.
#if !defined(CLAB_CMAIN) && !defined(__wasi__) && !defined(__EMSCRIPTEN__)
#define CLAB_UPD 1
#include <ctime>
#include <sys/stat.h>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <direct.h>
#include <windows.h>
#include <winhttp.h>
#ifdef _MSC_VER
#pragma comment(lib, "winhttp")
#endif
#define CLAB_TTY() (_isatty(_fileno(stderr)) != 0)
#else
#include <dlfcn.h>
#include <fcntl.h>
#include <unistd.h>
#define CLAB_TTY() (isatty(2) != 0)
#endif
#endif

#ifndef CLAB_VER
#include "ver.h"
#define CLAB_VER CLAB_VERSION
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
    "       crimlab -l | -h | -V | --update\n"
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
    "  --update  ask whether a newer release exists\n"
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

// RFC 4180 CSV: rows of cells; a quoted cell is text. A NUL byte is
// refused, naming its row and column: the reader of a number stops at it,
// so "4<NUL>0" would be read as 4.
bool csv(const std::string &s, std::vector<std::vector<Val>> &rows, std::string &err) {
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
  auto nul = [&]() {
    err = "CSV row " + std::to_string(rows.size() + 1) + ", column " + std::to_string(row.size() + 1) +
          " holds a NUL byte.";
    return false;
  };
  for (std::size_t i = 0; i < s.size(); ++i) {
    const char c = s[i];
    if (c == '\0') return nul();
    if (c == '"' && f.empty() && !quoted) {
      quoted = true;
      std::size_t j = i + 1;
      for (;; ++j) {
        if (j >= s.size()) return false;
        if (s[j] == '\0') return nul();
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
  if (!csv(s, rows, err) || rows.empty()) {
    if (err.empty()) err = "the CSV input cannot be read.";
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

#ifdef CLAB_UPD
// A newer release. The page of the latest release answers with the address
// of its tag. Nothing is sent but that request.
#define CLAB_REL "https://github.com/rootcoder007/crimlab/releases/"

// Is version a later than version b (numbers between dots).
bool later(const std::string &a, const std::string &b) {
  std::size_t i = 0, j = 0;
  while (i < a.size() || j < b.size()) {
    unsigned long x = 0, y = 0;
    while (i < a.size() && a[i] >= '0' && a[i] <= '9') x = x * 10 + static_cast<unsigned long>(a[i++] - '0');
    while (j < b.size() && b[j] >= '0' && b[j] <= '9') y = y * 10 + static_cast<unsigned long>(b[j++] - '0');
    if (x != y) return x > y;
    ++i;
    ++j;
  }
  return false;
}

/* The address the page of the latest release redirects to, in buf (cap
 * bytes, ended with 0); the length, or 0 when it cannot be had within secs
 * seconds. Nothing is linked for this: Windows asks through WinHTTP, which
 * is part of it; elsewhere the libcurl of the system is loaded when there
 * is one, and the curl program is the last resort. */
#ifdef _WIN32
static size_t asknet(int secs, char *buf, size_t cap) {
  HINTERNET s, c = NULL, r = NULL;
  wchar_t w[200];
  DWORD n = sizeof w, nored = WINHTTP_DISABLE_REDIRECTS;
  size_t k = 0, i;
  const int ms = secs * 1000;
  s = WinHttpOpen(L"crimlab", WINHTTP_ACCESS_TYPE_DEFAULT_PROXY,
                  WINHTTP_NO_PROXY_NAME, WINHTTP_NO_PROXY_BYPASS, 0);
  if (s == NULL)
    return 0;
  WinHttpSetTimeouts(s, ms, ms, ms, ms);
  c = WinHttpConnect(s, L"github.com", INTERNET_DEFAULT_HTTPS_PORT, 0);
  if (c != NULL)
    r = WinHttpOpenRequest(c, L"HEAD", L"/rootcoder007/crimlab/releases/latest",
                           NULL, WINHTTP_NO_REFERER,
                           WINHTTP_DEFAULT_ACCEPT_TYPES, WINHTTP_FLAG_SECURE);
  if (r != NULL &&
      WinHttpSetOption(r, WINHTTP_OPTION_DISABLE_FEATURE, &nored,
                       sizeof nored) &&
      WinHttpSendRequest(r, WINHTTP_NO_ADDITIONAL_HEADERS, 0,
                         WINHTTP_NO_REQUEST_DATA, 0, 0, 0) &&
      WinHttpReceiveResponse(r, NULL) &&
      WinHttpQueryHeaders(r, WINHTTP_QUERY_LOCATION,
                          WINHTTP_HEADER_NAME_BY_INDEX, w, &n,
                          WINHTTP_NO_HEADER_INDEX)) {
    k = (size_t)n / sizeof w[0];
    if (k >= cap)
      k = 0;
    for (i = 0; i < k; i++) {
      /* an address is ASCII; anything else is not the answer looked for */
      if (w[i] < 0x21 || w[i] > 0x7e) {
        k = 0;
        break;
      }
      buf[i] = (char)w[i];
    }
  }
  if (r != NULL)
    WinHttpCloseHandle(r);
  if (c != NULL)
    WinHttpCloseHandle(c);
  WinHttpCloseHandle(s);
  buf[k] = '\0';
  return k;
}
#else
/* libcurl, by the numbers of its stable interface (curl/curl.h):
 * CURLOPT_URL 10002, NOBODY 44, TIMEOUT 13, NOSIGNAL 99, USERAGENT 10018,
 * CURLINFO_REDIRECT_URL 0x100000 + 31. */
typedef void *(*askinit)(void);
typedef int (*askset)(void *, int, ...);
typedef int (*askrun)(void *);
typedef int (*askinfo)(void *, int, ...);
typedef void (*askend)(void *);

static size_t asklib(int secs, char *buf, size_t cap) {
  static const char *const names[] = {"libcurl.so.4", "libcurl-gnutls.so.4",
                                      "libcurl.4.dylib", "libcurl.so"};
  void *h = NULL, *p[5], *c;
  askinit finit;
  askset fset;
  askrun frun;
  askinfo finfo;
  askend fend;
  const char *url = NULL;
  size_t i, k = 0;
  for (i = 0; h == NULL && i < sizeof names / sizeof names[0]; i++)
    h = dlopen(names[i], RTLD_NOW | RTLD_LOCAL);
  if (h == NULL)
    return 0;
  p[0] = dlsym(h, "curl_easy_init");
  p[1] = dlsym(h, "curl_easy_setopt");
  p[2] = dlsym(h, "curl_easy_perform");
  p[3] = dlsym(h, "curl_easy_getinfo");
  p[4] = dlsym(h, "curl_easy_cleanup");
  for (i = 0; i < 5; i++) {
    if (p[i] == NULL)
      return 0;
  }
  /* an object pointer to a function pointer, as POSIX has dlsym give it */
  memcpy(&finit, &p[0], sizeof finit);
  memcpy(&fset, &p[1], sizeof fset);
  memcpy(&frun, &p[2], sizeof frun);
  memcpy(&finfo, &p[3], sizeof finfo);
  memcpy(&fend, &p[4], sizeof fend);
  c = finit();
  if (c == NULL)
    return 0;
  fset(c, 10002, CLAB_REL "latest");
  fset(c, 44, 1L);
  fset(c, 13, (long)secs);
  fset(c, 99, 1L);
  fset(c, 10018, "crimlab");
  if (frun(c) == 0 && finfo(c, 0x100000 + 31, &url) == 0 && url != NULL) {
    k = strlen(url);
    if (k >= cap)
      k = 0;
    memcpy(buf, url, k);
  }
  fend(c);
  /* the library stays loaded: the program is about to end, and unloading
   * a TLS library is not safe everywhere */
  buf[k] = '\0';
  return k;
}

static size_t asknet(int secs, char *buf, size_t cap) {
  char cmd[256];
  size_t n = asklib(secs, buf, cap);
  FILE *f;
  if (n > 0)
    return n;
  /* the curl program; the command is a constant but for the seconds */
  snprintf(cmd, sizeof cmd,
           "curl -fsS -m %d -o /dev/null -w \"%%{redirect_url}\" " CLAB_REL
           "latest 2>/dev/null",
           secs);
  f = popen(cmd, "r");
  if (f == NULL)
    return 0;
  n = fread(buf, 1, cap - 1, f);
  pclose(f);
  buf[n] = '\0';
  return n;
}
#endif

// The version of the latest release, or "" when it cannot be had (no
// network within secs seconds, no way to ask, an answer of another form).
std::string latest(int secs) {
  char buf[160];
  std::string s(buf, asknet(secs, buf, sizeof buf));
  const std::string pre = CLAB_REL "tag/v";
  if (s.compare(0, pre.size(), pre) != 0) return "";
  s.erase(0, pre.size());
  if (s.empty() || s.size() > 24 || s.find_first_not_of("0123456789.") != std::string::npos) return "";
  return s;
}

// The file that keeps the last answer and its time, or "" without a
// place for it.
std::string updfile() {
#ifdef _WIN32
  const char *d = std::getenv("LOCALAPPDATA");
  if (d == nullptr || *d == 0) return "";
  const std::string dir = std::string(d) + "\\crimlab";
  _mkdir(dir.c_str());
  return dir + "\\update";
#else
  const char *x = std::getenv("XDG_CACHE_HOME"), *h = std::getenv("HOME");
  std::string dir;
  if (x != nullptr && *x != 0) {
    dir = x;
  } else if (h != nullptr && *h != 0) {
    dir = std::string(h) + "/.cache";
  } else {
    return "";
  }
  mkdir(dir.c_str(), 0700);
  dir += "/crimlab";
  mkdir(dir.c_str(), 0700);
  return dir + "/update";
#endif
}

// One line on standard error when a newer release exists: only at a
// terminal, never under CI or with CRIMLAB_NO_UPDATE_CHECK set, and the
// question is asked at most once a day (two seconds at most), the answer
// kept in a file between.
void notice() {
  const char *off = std::getenv("CRIMLAB_NO_UPDATE_CHECK");
  if ((off != nullptr && *off != 0) || std::getenv("CI") != nullptr || !CLAB_TTY()) return;
  const std::string f = updfile();
  if (f.empty()) return;
  long long t = 0;
  char v[32] = "";
  if (std::FILE *c = std::fopen(f.c_str(), "r")) {
    if (std::fscanf(c, "%lld %24[0-9.]", &t, v) < 1) t = 0;
    std::fclose(c);
  }
  const long long now = static_cast<long long>(std::time(nullptr));
  std::string ver = v;
  if (now < t || now - t >= 86400) {
    ver = latest(2);
    // a file of the user alone
#ifdef _WIN32
    std::FILE *c = std::fopen(f.c_str(), "w");
#else
    const int fd = open(f.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);
    std::FILE *c = fd >= 0 ? fdopen(fd, "w") : nullptr;
    if (c == nullptr && fd >= 0) close(fd);
#endif
    if (c != nullptr) {
      std::fprintf(c, "%lld %s\n", now, ver.c_str());
      std::fclose(c);
    }
  }
  if (!ver.empty() && later(ver, CLAB_VER)) {
    std::fprintf(stderr,
                 "crimlab: version %s is available (this is %s): " CLAB_REL "latest\n"
                 "crimlab: CRIMLAB_NO_UPDATE_CHECK=1 stops this notice.\n",
                 ver.c_str(), CLAB_VER);
  }
}
#endif

}  // namespace

static int mainb(int argc, char **argv) {
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
  if (fn == "--update") {
#ifdef CLAB_UPD
    const std::string v = latest(10);
    if (v.empty()) return fail("the latest release could not be asked for (no network, or nothing to ask with).", 1);
    if (later(v, CLAB_VER)) {
      std::printf("crimlab %s is available (this is %s): " CLAB_REL "latest\n", v.c_str(), CLAB_VER);
    } else {
      std::printf("crimlab %s is up to date.\n", CLAB_VER);
    }
    return 0;
#else
    return fail("this build cannot ask for the latest release.", 1);
#endif
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

// The program: mainb(), with a failed allocation anywhere in it (reading
// the input, the call, the printing) reported in words and as status 1,
// never an abort. Without exceptions (WebAssembly) clab_jrun()'s
// new-handler answers instead.
int main(int argc, char **argv) {
#if defined(__cpp_exceptions) || defined(_CPPUNWIND)
  try {
    const int rc = mainb(argc, argv);
#ifdef CLAB_UPD
    // after the result of a function, never after -h, -l, -V or --update
    if (rc == 0 && argc > 1 && argv[1][0] != '-') notice();
#endif
    return rc;
  } catch (...) {
    std::fputs("crimlab: not enough memory for this size.\n", stderr);
    return 1;
  }
#else
  return mainb(argc, argv);
#endif
}
