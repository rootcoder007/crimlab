// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The figures of hwpwr() over an interval (rows(), src/pwr.h): energy from
// counters (pairwise differences, with the counter's range added back at a
// wrap), from power sampled at every snapshot (the trapezoid rule), or
// estimated as a rated power times a utilisation, with the range such an
// estimate carries. An unknown figure stays unknown (NaN); the total names
// the devices it could not count.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "pwr.h"

namespace pwr {

namespace {

const double nan = std::nan("");

// The range of an estimate from a rated power and a utilisation, as factors
// of the estimate: the power-from-load model of CodeCarbon, compared with
// RAPL on four CPUs, was off by -60% to +190% (CodeCarbon documentation,
// "Accuracy"). A linear model of a rated power is not better.
constexpr double elo = 0.4, ehi = 2.9;

int korder(const std::string &k) {
  static const char *o[] = {"cpu", "memory", "gpu", "accelerator", "meter", "machine"};
  for (int i = 0; i < 6; ++i) {
    if (k == o[i]) return i;
  }
  return 6;
}

std::string pct(double x) {
  char b[32];
  std::snprintf(b, sizeof b, "%.1f%%", 100 * x);
  return b;
}

std::string wt(double x) {
  char b[32];
  std::snprintf(b, sizeof b, "%.3g W", x);
  return b;
}

void note(std::string &w, const std::string &n) {
  if (n.empty()) return;
  w += w.empty() ? n : "; " + n;
}

double med(std::vector<double> v) {
  if (v.empty()) return nan;
  std::sort(v.begin(), v.end());
  const std::size_t m = v.size() / 2;
  return v.size() % 2 ? v[m] : (v[m - 1] + v[m]) / 2;
}

}  // namespace

std::vector<Row> rows(const std::vector<Snap> &s, const Rated &rt) {
  std::vector<Row> out;
  if (s.empty()) return out;
  const std::size_t n = s.size();
  const double T = s.back().t - s.front().t;
  // the machine's busy share and the process's share of the busy time
  double su = nan, ps = nan;
  if (n > 1 && s.front().scpu >= 0 && s.back().scpu >= 0) {
    const double db = s.back().scpu - s.front().scpu, dt = s.back().stot - s.front().stot;
    if (dt > 0 && db >= 0) su = std::min(1.0, db / dt);
    if (db > 0 && s.front().pcpu >= 0 && s.back().pcpu >= 0) {
      ps = std::max(0.0, std::min(1.0, (s.back().pcpu - s.front().pcpu) / db));
    }
  }
  // devices by id, in the order of the last snapshot
  std::vector<std::string> ids;
  std::map<std::string, std::vector<const Dev *>> at;
  for (const Dev &d : s.back().d) {
    if (!at.count(d.id)) ids.push_back(d.id);
    at[d.id].assign(n, nullptr);
  }
  for (std::size_t k = 0; k < n; ++k) {
    for (const Dev &d : s[k].d) {
      auto it = at.find(d.id);
      if (it != at.end()) it->second[k] = &d;
    }
  }
  struct Acc {
    Row r;
    const Dev *last;
    double c0 = nan, dc = nan;  // first counter, its total change
    bool est = false, mir = false;
    double pown = nan;
  };
  std::vector<Acc> acc;
  for (const std::string &id : ids) {
    const std::vector<const Dev *> &v = at[id];
    const Dev &L = *v[n - 1];
    Acc a;
    a.last = &L;
    Row &r = a.r;
    r.kind = L.kind;
    r.vendor = L.vendor;
    r.name = L.name;
    r.scope = L.kind == "machine" ? "machine" : "device";
    r.how = "unknown";
    r.watts = r.joules = r.util = r.flo = r.fhi = nan;
    std::string why = L.why;
    // utilisation: the mean of the readings that have one
    double us = 0, pw = 0;
    int un = 0, pn = 0;
    for (const Dev *d : v) {
      if (d != nullptr && d->util >= 0) {
        us += d->util;
        ++un;
      }
      if (d != nullptr && d->pown >= 0) {
        pw += d->pown;
        ++pn;
      }
    }
    if (un) r.util = us / un;
    if (pn) a.pown = pw / pn;
    if (r.kind == "cpu" && std::isnan(r.util)) r.util = su;
    if (L.bad || L.src == 0) {
      // nothing measured
    } else if (n == 1 || T <= 0) {
      if (L.w >= 0) {
        r.watts = L.w;
        r.how = L.src == 1 ? "measured counter" : "measured power integrated";
        note(why, "power at one instant; energy needs an interval (dt > 0 or since)");
      } else {
        note(why, "a counter needs an interval: dt > 0 or since");
      }
    } else if (L.src == 1) {
      // pairwise differences of the counter
      std::vector<double> dd, tt;
      bool gap = false, viaw = false, longw = false;
      for (std::size_t k = 1; k < n; ++k) {
        const Dev *p = v[k - 1], *q = v[k];
        const double ti = s[k].t - s[k - 1].t;
        if (p == nullptr || q == nullptr || p->bad || q->bad || p->src != 1 || q->src != 1) {
          gap = true;
          continue;
        }
        double d = q->ctr - p->ctr;
        // a fall is a wrap when the energy it implies is at most half the
        // range, or no more than 1 kW (the scale of the note on long
        // intervals below) could use in the interval; a larger fall is a
        // reset of the counter
        if (d < 0 && q->wrap > 0 && d + q->wrap <= std::fmax(q->wrap / 2, 1000 * ti)) d += q->wrap;
        if (d < 0) {
          // a reset: the power readings at the two ends, when there are any
          if (p->w >= 0 && q->w >= 0) {
            d = (p->w + q->w) / 2 * ti;
            viaw = true;
          } else {
            gap = true;
            continue;
          }
        }
        if (q->wrap > 0 && ti > q->wrap / 1000) longw = true;
        dd.push_back(d);
        tt.push_back(ti);
        if (std::isnan(a.c0)) a.c0 = p->ctr;
      }
      double J = 0;
      for (double d : dd) J += d;
      bool drop = false;
      if (L.aon && !dd.empty() && J == 0) {
        drop = true;
        note(why, "the counter did not change over the interval; not reported");
      } else if (L.aon && dd.size() >= 2) {
        // an always-on counter that stops changing has frozen (a known
        // fault of the macOS energy counters): at the end it under-counts,
        // and a first step far above the rest carries energy from before
        // the interval; either way the figure is not reported. A pause in
        // the middle followed by a jump is caught up within the interval
        // and is spread over it.
        bool moved = false, mid = false;
        for (std::size_t k = 0; k + 1 < dd.size(); ++k) {
          if (dd[k] > 0) moved = true;
          if (dd[k] == 0 && dd[k + 1] > 0) mid = true;
        }
        std::vector<double> rest;
        for (std::size_t k = 1; k < dd.size(); ++k) rest.push_back(dd[k] / tt[k]);
        const double m = med(rest);
        if (moved && dd.back() == 0) {
          drop = true;
          note(why, "the counter stopped changing before the end of the interval; not reported");
        } else if (dd.size() >= 3 && m > 0 && dd[0] / tt[0] > 10 * m) {
          drop = true;
          note(why, "the first step caught up energy from before the interval; not reported");
        } else if (mid) {
          note(why, "the counter paused and caught up; its energy is spread over the pause");
        }
      }
      if (gap || drop || dd.empty()) {
        if (gap) note(why, "the counter was not readable at every snapshot or went back with no power read");
      } else {
        r.joules = J;
        r.watts = J / T;
        a.dc = J;
        r.how = "measured counter";
        if (viaw) note(why, "the counter was reset once; that step comes from the power readings");
        if (longw) note(why, "readings far enough apart for the counter to wrap more than once at high power");
      }
    } else {
      // power only: the trapezoid rule over the snapshots
      double J = 0, cov = 0;
      for (std::size_t k = 1; k < n; ++k) {
        const Dev *p = v[k - 1], *q = v[k];
        if (p == nullptr || q == nullptr || p->w < 0 || q->w < 0) continue;
        const double ti = s[k].t - s[k - 1].t;
        J += (p->w + q->w) / 2 * ti;
        cov += ti;
      }
      if (cov >= T * (1 - 1e-9) && cov > 0) {
        r.joules = J;
        r.watts = J / T;
        r.how = "measured power integrated";
        if (n == 2) note(why, "from two readings, at the start and the end");
      } else {
        note(why, "power was not read at every snapshot");
      }
    }
    // an estimate where nothing was measured: a rated power times the
    // utilisation, never for a device inside another figure
    if (std::isnan(r.joules) && !L.inpkg && r.kind != "machine" && !std::isnan(r.util) && T > 0) {
      double rp = L.rated;
      std::string from = "the rated power the machine reports";
      if (r.kind == "cpu" && rt.cpu >= 0) {
        rp = rt.cpu;
        from = "the rated power given";
      }
      if (r.kind == "gpu" && rt.gpu >= 0) {
        rp = rt.gpu;
        from = "the rated power given";
      }
      if (rp >= 0) {
        r.watts = rp * r.util;
        r.joules = r.watts * T;
        r.how = "rated x utilisation";
        r.flo = elo;
        r.fhi = ehi;
        a.est = true;
        note(why, "estimated from " + from + " (" + wt(rp) + ") x utilisation " + pct(r.util) +
                      "; the true figure may lie from 0.4 to 2.9 times it");
      }
    }
    r.why = why;
    acc.push_back(a);
  }
  // Mirrored CPU package counters (one counter exposed once per die, as
  // multi-die AMD parts do): equal at the start and equal in change.
  for (std::size_t i = 0; i < acc.size(); ++i) {
    Acc &a = acc[i];
    if (a.r.kind != "cpu" || std::isnan(a.dc)) continue;
    for (std::size_t j = 0; j < i; ++j) {
      const Acc &b = acc[j];
      if (b.r.kind != "cpu" || std::isnan(b.dc) || b.mir) continue;
      if (std::fabs(a.c0 - b.c0) <= std::max(1.0, 1e-6 * std::fabs(b.c0)) &&
          std::fabs(a.dc - b.dc) <= std::max(1e-3, 1e-6 * std::fabs(b.dc))) {
        a.mir = true;
        a.r.joules = a.r.watts = nan;
        a.r.how = "unknown";
        note(a.r.why, "the same counter as " + b.r.name + " (a per-die copy); counted once");
        break;
      }
    }
  }
  std::stable_sort(acc.begin(), acc.end(), [](const Acc &x, const Acc &y) { return korder(x.r.kind) < korder(y.r.kind); });
  // the total over the devices, and the process's shares
  double tj = 0, tw = 0, lo = 0, hi = 0;
  bool anyj = false, anyw = false, anye = false;
  std::string miss;
  std::vector<Row> proc;
  for (const Acc &a : acc) {
    out.push_back(a.r);
    const Row &r = a.r;
    const bool dev = r.kind == "cpu" || r.kind == "memory" || r.kind == "gpu" || r.kind == "accelerator";
    if (!dev || a.last->inpkg || a.mir) continue;
    if (!std::isnan(r.watts)) {
      tw += r.watts;
      anyw = true;
    }
    if (!std::isnan(r.joules)) {
      tj += r.joules;
      lo += a.est ? r.joules * elo : r.joules;
      hi += a.est ? r.joules * ehi : r.joules;
      anyj = true;
      anye = anye || a.est;
    } else {
      note(miss, r.name);
    }
    double sh = nan;
    std::string how;
    if (r.kind == "cpu" && !std::isnan(ps)) {
      sh = ps;
      how = "the process's share of the busy CPU time, " + pct(ps);
    } else if (r.kind == "gpu" && !std::isnan(a.pown)) {
      sh = a.pown;
      how = "the process's share of the GPU's compute samples, " + pct(a.pown);
    }
    if (!std::isnan(sh) && !std::isnan(r.joules)) {
      Row p = r;
      p.scope = "process";
      p.joules = r.joules * sh;
      p.watts = r.watts * sh;
      p.util = sh;
      p.why = how;
      proc.push_back(p);
    }
  }
  for (const Row &p : proc) out.push_back(p);
  Row t;
  t.kind = "total";
  t.vendor = "";
  t.name = "all devices";
  t.scope = "device";
  t.how = anye ? "sum, with estimates" : "sum";
  t.watts = anyw ? tw : nan;
  t.joules = anyj ? tj : nan;
  t.util = nan;
  t.flo = anye && tj > 0 ? lo / tj : nan;
  t.fhi = anye && tj > 0 ? hi / tj : nan;
  t.why = miss.empty() ? "" : "not counted: " + miss;
  out.push_back(t);
  return out;
}

}  // namespace pwr
