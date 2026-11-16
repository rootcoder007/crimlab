// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// What the footprint of a computation reads from the machine: where it
// runs (the environment, the time zone, the locale), the RAPL counters of
// the CPU packages, the kernel's CPU accounting and the memory size; the
// row of a location or of its country; and the assumptions the flags of
// clab_fpe() record. Files are read with pwr::rd() and listed with
// pwr::ls() (pwr.cpp); nothing here blocks on a FIFO.
#ifdef _WIN32
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#endif

#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif
#ifdef __APPLE__
#include <sys/sysctl.h>
#include <sys/types.h>
#endif

#include "crimlab.h"
#include "kutil.h"
#include "pwr.h"

namespace {

const double qnan = std::numeric_limits<double>::quiet_NaN();

// The table's location code for key (static storage), nullptr when key is
// empty or not a location the table knows.
const char *known(const char *key) {
  if (key == nullptr || *key == 0) return nullptr;
  const int i = clab_cimatch(key);
  if (i < 0 || i >= clab_cin()) return nullptr;
  const char *loc = nullptr, *iso3 = nullptr, *name = nullptr, *src = nullptr;
  double g = 0.0;
  int yr = 0;
  clab_cirow(i, &loc, &iso3, &name, &g, &yr, &src);
  return loc;
}

// The first line of the file at path (LF, CRLF or CR ending it, as R's
// readLines() takes them), false when the file is unreadable or empty.
bool line1(const std::string &path, std::string &l) {
  std::string all, err;
  l.clear();
  if (!pwr::rd(path, all, err) || all.empty()) return false;
  l = all.substr(0, all.find_first_of("\r\n"));
  return true;
}

// The text s as one number, space around it allowed; NaN when it is not.
double num(const std::string &s) {
  const char *b = s.c_str();
  char *e = nullptr;
  const double v = std::strtod(b, &e);
  if (e == b) return qnan;
  while (*e != 0 && std::isspace(static_cast<unsigned char>(*e))) ++e;
  return *e == 0 ? v : qnan;
}

bool exists(const std::string &p) {
#ifdef _WIN32
  return _access(p.c_str(), 0) == 0;
#else
  return access(p.c_str(), F_OK) == 0;
#endif
}

// "^[A-Za-z]{2,3}_[A-Z]{2}": the country of a locale name, "" when none.
std::string lcterr(const std::string &lc) {
  auto alpha = [](char ch) { return (ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z'); };
  auto up = [](char ch) { return ch >= 'A' && ch <= 'Z'; };
  std::size_t k = 0;
  while (k < lc.size() && k < 3 && alpha(lc[k])) ++k;
  if (k < 2 || k + 3 > lc.size() || lc[k] != '_' || !up(lc[k + 1]) || !up(lc[k + 2])) return "";
  return lc.substr(k + 1, 2);
}

// The method text into how (cap bytes, terminated), its length into len.
void putstr(const std::string &s, char *how, std::size_t cap, std::size_t *len) {
  *len = s.size();
  if (cap == 0 || how == nullptr) return;
  const std::size_t n = s.size() < cap - 1 ? s.size() : cap - 1;
  std::memcpy(how, s.data(), n);
  how[n] = 0;
}

}  // namespace

extern "C" {

int clab_cizone(const char *key, int *zone) CLAB_TRY {
  *zone = 0;
  int i = clab_cimatch(key);
  if (i >= 0) return i;
  const char *d = std::strchr(key, '-');
  if (d == nullptr) return -1;
  i = clab_cimatch(std::string(key, static_cast<std::size_t>(d - key)).c_str());
  if (i >= 0) *zone = 1;
  return i;
} CLAB_CATCH(-1)

int clab_fploc(const char *env, const char *cod, const char *tz, const char *lc, const char **loc,
               char *how, std::size_t cap, std::size_t *len, int *ign) CLAB_TRY {
  std::string h;
  const char *l = known(env);
  int ig = 0;
  if (l != nullptr) {
    h = "CRIMLAB_LOCATION environment variable";
  } else {
    if (*env != 0) ig |= 1;
    l = known(cod);
    if (l != nullptr) h = "CODECARBON_COUNTRY_ISO_CODE environment variable";
    else if (*cod != 0) ig |= 2;
  }
  const char *terr = l == nullptr && *tz != 0 ? clab_tzterr(tz) : nullptr;
  if (terr != nullptr) {
    l = known(terr);
    h = std::string("system time zone ") + tz;
    if (l == nullptr) {
      // a territory the carbon table has no row for (Monaco, Jersey, the
      // Antarctic stations): its administering or neighbouring country's
      // grid, or the world average, never the locale's country
      l = known(clab_terrproxy(terr));
      if (l == nullptr) l = "WORLD";
      h += std::string(" (territory ") + terr + " has no row of its own; " +
           (std::strcmp(l, "WORLD") == 0 ? std::string("world average") : std::string("grid of ") + l) + ")";
    }
  }
  if (l == nullptr) {
    l = known(lcterr(lc).c_str());
    if (l != nullptr) h = std::string("locale ") + lc;
  }
  if (l == nullptr) {
    l = "WORLD";
    h = "fallback, no location found in environment, time zone or locale";
  }
  *loc = l;
  *ign = ig;
  putstr(h, how, cap, len);
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_rapldom(const char *dir, char *buf, std::size_t cap, std::size_t *len, int *n) CLAB_TRY {
  std::string all;
  int k = 0;
  const std::string d = dir;
  for (const std::string &e : pwr::ls(d)) {
    if (e.size() <= 11 || e.compare(0, 11, "intel-rapl:") != 0) continue;
    bool dig = true;
    for (std::size_t j = 11; j < e.size(); ++j) dig = dig && e[j] >= '0' && e[j] <= '9';
    if (!dig) continue;
    const std::string p = d + "/" + e;
    std::string nm;
    if (!line1(p + "/name", nm) || nm.compare(0, 7, "package") != 0 || !exists(p + "/energy_uj")) continue;
    all += p;
    all.push_back(0);
    ++k;
  }
  *n = k;
  *len = all.size();
  if (all.size() <= cap && buf != nullptr) std::memcpy(buf, all.data(), all.size());
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

void clab_raplrd(const char *dom, double *uj, double *mx) CLAB_TRY {
  *uj = qnan;
  *mx = qnan;
  const std::string p = dom;
  std::string l;
  if (line1(p + "/energy_uj", l)) *uj = num(l);
  if (line1(p + "/max_energy_range_uj", l)) *mx = num(l);
} CLAB_CATCHV

int clab_cpuacc(const char *path, double *out) CLAB_TRY {
  out[0] = qnan;
  out[1] = qnan;
  std::string l;
  if (!line1(path, l) || l.compare(0, 4, "cpu ") != 0) return CLAB_EDOM;
  std::vector<double> v;
  std::size_t i = 0;
  bool first = true;
  while (i < l.size()) {
    while (i < l.size() && std::isspace(static_cast<unsigned char>(l[i]))) ++i;
    std::size_t j = i;
    while (j < l.size() && !std::isspace(static_cast<unsigned char>(l[j]))) ++j;
    if (j > i) {
      if (!first) v.push_back(num(l.substr(i, j - i)));
      first = false;
    }
    i = j;
  }
  if (v.size() < 4) return CLAB_EDOM;
  long double tot = 0.0L;
  for (const double x : v) {
    if (std::isnan(x)) return CLAB_EDOM;
    tot += x;
  }
  out[0] = v[3] + (v.size() >= 5 ? v[4] : 0.0);
  out[1] = static_cast<double>(tot);
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

double clab_memgb(void) CLAB_TRY {
  std::string all, err;
  if (exists("/proc/meminfo") && pwr::rd("/proc/meminfo", all, err)) {
    int hits = 0;
    std::string dig;
    std::size_t i = 0;
    while (i < all.size()) {
      std::size_t j = all.find_first_of("\r\n", i);
      if (j == std::string::npos) j = all.size();
      if (all.compare(i, 9, "MemTotal:") == 0) {
        ++hits;
        dig.clear();
        for (std::size_t k = i; k < j; ++k) {
          if (all[k] >= '0' && all[k] <= '9') dig.push_back(all[k]);
        }
      }
      i = j + 1;
    }
    const double kb = hits == 1 && !dig.empty() ? std::strtod(dig.c_str(), nullptr) : qnan;
    if (std::isfinite(kb) && kb > 0) return kb / 1048576.0;
  }
#ifdef __APPLE__
  std::uint64_t b = 0;
  std::size_t sz = sizeof b;
  if (sysctlbyname("hw.memsize", &b, &sz, nullptr, 0) == 0 && b > 0) {
    return static_cast<double>(b) / 1073741824.0;
  }
#endif
  return qnan;
} CLAB_CATCH(qnan)

const char *clab_fpnote(int i, int flags, int machine) {
  static const int bit[] = {1, 4, 2, 8, 16};
  static const char *const txt[] = {
      "machine CPU accounting became unreadable; process utilisation used",
      "RAPL counters became unreadable during the run; modelled instead",
      "a RAPL counter wrapped and its range was unknown; modelled instead",
      "CPU power unknown: 85 W x 50% (codecarbon fallback)",
      "codecarbon process load curve applied to the CPU power"};
  if (i < 0 || i > 4 || (flags & bit[i]) == 0) return nullptr;
  if (i == 4 && machine != 0) return "codecarbon machine load curve applied to the CPU power";
  return txt[i];
}

}  // extern "C"
