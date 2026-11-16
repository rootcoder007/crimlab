// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The energy and carbon footprint of a computation: the Green Algorithms
// model of energy from running time, utilisation and memory (Lannelongue,
// Grealey and Inouye 2021), codecarbon's fallbacks and load curves, the
// RAPL counters' wrap-around, the measured energy of other devices, and
// the lookups of the bundled carbon-intensity and time-zone tables
// (fpdat.cpp). No file is read here: the callers read the counters.
#include <cctype>
#include <cmath>
#include <cstring>
#include <string>

#include "crimlab.h"
#include "fpdat.h"
#include "kutil.h"

namespace {

// toupper(trimws(s)) for ASCII text (every name in the table is ASCII).
std::string ukey(const char *s) {
  std::string k(s);
  const char *ws = " \t\n\r\f\v";
  const std::size_t a = k.find_first_not_of(ws);
  if (a == std::string::npos) return "";
  k = k.substr(a, k.find_last_not_of(ws) - a + 1);
  for (char &ch : k) {
    if (ch >= 'a' && ch <= 'z') ch = static_cast<char>(ch - 'a' + 'A');
  }
  return k;
}

std::string upper(const char *s) {
  std::string k(s);
  for (char &ch : k) {
    if (ch >= 'a' && ch <= 'z') ch = static_cast<char>(ch - 'a' + 'A');
  }
  return k;
}

// Territories the zone table names that the carbon table has no row for,
// with the country whose grid they draw on, or "" for the world average.
const char *const proxy[][2] = {
    {"AD", "ES"}, {"AI", "AG"}, {"BL", "GP"}, {"BQ", "NL"}, {"CC", "AU"}, {"CX", "AU"}, {"GG", "FR"},
    {"IM", "GB"}, {"IO", "GB"}, {"JE", "FR"}, {"LI", "CH"}, {"MC", "FR"}, {"MF", "GP"}, {"MP", "GU"},
    {"NF", "AU"}, {"PN", "GB"}, {"SJ", "NO"}, {"SM", "IT"}, {"SX", "NL"}, {"UM", "US"}, {"VA", "IT"},
    {"WF", "NC"}, {"AQ", ""},   {"MH", ""},   {"NU", ""},   {"TF", ""},   {"TK", ""},   {"TV", ""}};

double rpw(double x, double y) { return y == 2.0 ? x * x : std::pow(x, y); }

}  // namespace

extern "C" {

int clab_cin(void) CLAB_TRY { return cfp::nci; } CLAB_CATCH(CLAB_ENOMEM)

void clab_cirow(int i, const char **loc, const char **iso3, const char **name, double *g, int *year,
                const char **src) CLAB_TRY {
  const cfp::Ci &r = cfp::ci[i];
  *loc = r.loc;
  *iso3 = r.iso3;
  *name = r.name;
  *g = r.g;
  *year = r.year;
  *src = cfp::src[r.src];
} CLAB_CATCHV

int clab_cimatch(const char *key) CLAB_TRY {
  const std::string k = ukey(key);
  for (int i = 0; i < cfp::nci; ++i) {
    if (k == cfp::ci[i].loc) return i;
  }
  for (int i = 0; i < cfp::nci; ++i) {
    if (k == cfp::ci[i].iso3) return i;
  }
  for (int i = 0; i < cfp::nci; ++i) {
    if (k == upper(cfp::ci[i].name)) return i;
  }
  return -1;
} CLAB_CATCH(CLAB_ENOMEM)

const char *clab_tzterr(const char *tz) {
  const char *k = tz;
  if (std::strncmp(k, "posix/", 6) == 0) k += 6;
  else if (std::strncmp(k, "right/", 6) == 0) k += 6;
  for (int i = 0; i < cfp::ntz; ++i) {
    if (std::strcmp(k, cfp::tz[i].tz) == 0) return cfp::tz[i].terr;
  }
  return nullptr;
}

const char *clab_terrproxy(const char *terr) {
  for (const auto &p : proxy) {
    if (std::strcmp(terr, p[0]) == 0) return p[1];
  }
  return "";
}

int clab_fpe(const double *in, int nr, const double *r0, const double *r1, const double *rmx, int nd,
             const double *dw, const double *de, const int *dk, double *out, double *dout) CLAB_TRY {
  if (nr < 0 || nd < 0) return CLAB_EARG;
  const double cpu = in[1], cores = in[2], ufix = in[4], pw = in[9], memgb = in[11], mpw = in[12],
               pue = in[13], g = in[14];
  int mode = static_cast<int>(in[3]), flags = 0;
  // the running time, floored at the CPU time over the cores
  double wall = in[0] > 0 ? in[0] : 0;
  if (std::isfinite(cpu) && cpu > 0) {
    const double fl = cpu / (cores > 1 ? cores : 1);
    if (fl > wall) wall = fl;
  }
  const double hours = wall / 3600;
  double u = ufix;
  if (mode == 1) {
    if (std::isnan(in[7]) || std::isnan(in[8])) {
      mode = 0;
      flags |= 1;
    } else {
      const double dt = in[8] - in[6];
      u = 0;
      if (dt > 0) {
        u = 1 - (in[7] - in[5]) / dt;
        u = u < 0 ? 0 : u;
        u = u > 1 ? 1 : u;
      }
    }
  }
  if (mode == 0) {
    u = 0;
    if (wall > 0) {
      u = cpu / (wall * cores);
      u = u < 0 ? 0 : u;
      u = u > 1 ? 1 : u;
    }
  }
  // the CPU: the counters' difference, wrapped through their range
  int rapl = in[15] != 0 && nr > 0;
  double ecpu = 0;
  if (rapl) {
    bool lost = false;
    for (int i = 0; i < nr; ++i) lost = lost || std::isnan(r1[i]);
    if (lost) {
      rapl = 0;
      flags |= 4;
    } else {
      long double s = 0;
      bool unk = false;
      for (int i = 0; i < nr; ++i) unk = unk || (r1[i] - r0[i] < 0 && std::isnan(rmx[i]));
      if (unk) {
        rapl = 0;
        flags |= 2;
      } else {
        for (int i = 0; i < nr; ++i) {
          double d = r1[i] - r0[i];
          if (d < 0) d += rmx[i];
          s += d;
        }
        ecpu = static_cast<double>(s) * 1e-6 / 3.6e6;
      }
    }
  }
  if (!rapl) {
    if (std::isnan(pw)) {
      flags |= 8;
      ecpu = hours * 85 * 0.5 / 1000;
    } else {
      double p = pw * u;
      if (in[10] != 0) {
        flags |= 16;
        if (mode == 1) {
          p = pw * (0.1 + 0.9 * rpw(u, 3));
        } else {
          const double lo = 0.1 * pw;
          p = p > lo ? p : lo;
        }
      }
      ecpu = hours * p / 1000;
    }
  }
  const double emem = hours * memgb * mpw / 1000;
  // other devices: measured energy, else average power over the run
  double edev = 0;
  for (int i = 0; i < nd; ++i) {
    double e = NAN;
    if (dk[i] == 0) {
      if (std::isfinite(de[i]) && de[i] >= 0) e = de[i];
      else if (std::isfinite(dw[i]) && dw[i] >= 0) e = hours * dw[i] / 1000;
    }
    dout[i] = e;
    if (!std::isnan(e)) edev += e;
  }
  const double total = (nd > 0 ? ecpu + emem + edev : ecpu + emem) * pue;
  out[0] = wall;
  out[1] = u;
  out[2] = mode;
  out[3] = rapl ? 0 : 1;
  out[4] = flags;
  out[5] = ecpu * pue;
  out[6] = emem * pue;
  out[7] = edev * pue;
  out[8] = total;
  out[9] = total * g;
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

void clab_fpeqv(double g, double *out) CLAB_TRY {
  out[0] = g / 175;
  out[1] = g / 917;
} CLAB_CATCHV

}  // extern "C"
