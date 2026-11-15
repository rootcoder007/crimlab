// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/pwr.cpp and src/pwrd.cpp on constructed snapshots and
// recorded text: no file, program or device is read.
#include <cmath>
#include <string>
#include <vector>

#include "pwr.h"
#include "tchk.h"

namespace {

pwr::Dev ctr(const char *id, const char *kind, double c, double wrap = 0) {
  pwr::Dev d;
  d.id = id;
  d.kind = kind;
  d.name = id;
  d.src = 1;
  d.ctr = c;
  d.wrap = wrap;
  return d;
}

pwr::Dev pw(const char *id, const char *kind, double w) {
  pwr::Dev d;
  d.id = id;
  d.kind = kind;
  d.name = id;
  d.src = 2;
  d.w = w;
  return d;
}

const pwr::Row *find(const std::vector<pwr::Row> &r, const std::string &name, const char *scope = "device") {
  for (const pwr::Row &x : r) {
    if (x.name == name && x.scope == scope) return &x;
  }
  return nullptr;
}

double J(const std::vector<pwr::Row> &r, const std::string &name, const char *scope = "device") {
  const pwr::Row *x = find(r, name, scope);
  return x == nullptr ? -1 : x->joules;
}

}  // namespace

void t_pwr(tk &T, int big) {
  (void)big;
  pwr::Rated none;
  // a counter with a wrap, a power-only device, the process share
  {
    std::vector<pwr::Snap> s(3);
    for (int k = 0; k < 3; ++k) {
      s[k].t = 10.0 + k;
      s[k].scpu = 100.0 + 50 * k;
      s[k].stot = 1000.0 + 400 * k;
      s[k].pcpu = 5.0 + 10 * k;
    }
    s[0].d = {ctr("pkg", "cpu", 90), pw("gpu", "gpu", 10)};
    s[1].d = {ctr("pkg", "cpu", 99), pw("gpu", "gpu", 20)};
    s[2].d = {ctr("pkg", "cpu", 4, 100), pw("gpu", "gpu", 30)};
    s[1].d[0].wrap = s[0].d[0].wrap = 100;
    const std::vector<pwr::Row> r = pwr::rows(s, none);
    // 9 J, then 99 -> 4 through the range of 100: 5 J
    tk_eq(T, "pw.wrap", J(r, "pkg"), 14, 1e-12);
    tk_eq(T, "pw.w", find(r, "pkg")->watts, 7, 1e-12);
    // trapezoids: 15 + 25
    tk_eq(T, "pw.trap", J(r, "gpu"), 40, 1e-12);
    // the process: 20 of 100 busy ticks
    tk_eq(T, "pw.proc", J(r, "pkg", "process"), 14 * 0.2, 1e-12);
    tk_eq(T, "pw.tot", J(r, "all devices"), 54, 1e-12);
    tk_true(T, "pw.tot.why", find(r, "all devices")->why.empty());
    // the utilisation of a CPU is the machine's busy share, 100 / 800
    tk_eq(T, "pw.util", find(r, "pkg")->util, 0.125, 1e-12);
  }
  // a reset with no power reading leaves the energy unknown, and the total
  // names it
  {
    std::vector<pwr::Snap> s(2);
    s[0].t = 0;
    s[1].t = 1;
    s[0].d = {ctr("nv", "gpu", 50)};
    s[1].d = {ctr("nv", "gpu", 3)};
    const std::vector<pwr::Row> r = pwr::rows(s, none);
    tk_nan(T, "pw.reset", J(r, "nv"));
    tk_true(T, "pw.miss", find(r, "all devices")->why == "not counted: nv");
    s[0].d[0].w = 40;
    s[1].d[0].w = 60;
    tk_eq(T, "pw.resetw", J(pwr::rows(s, none), "nv"), 50, 1e-12);
  }
  // mirrored per-die package counters count once
  {
    std::vector<pwr::Snap> s(2);
    s[0].t = 0;
    s[1].t = 2;
    s[0].d = {ctr("p0", "cpu", 1000), ctr("p1", "cpu", 1000.5), ctr("p2", "cpu", 3000)};
    s[1].d = {ctr("p0", "cpu", 1040), ctr("p1", "cpu", 1040.5), ctr("p2", "cpu", 3010)};
    const std::vector<pwr::Row> r = pwr::rows(s, none);
    tk_nan(T, "pw.mir", J(r, "p1"));
    tk_eq(T, "pw.mirtot", J(r, "all devices"), 50, 1e-12);
  }
  // an always-on counter: steady, paused and caught up, frozen at the end,
  // a first step that carries energy from before
  {
    const double seq[4][5] = {{0, 1, 2, 3, 4}, {0, 1, 1, 3, 4}, {0, 1, 2, 3, 3}, {0, 30, 31, 32, 33}};
    const double want[4] = {4, 4, NAN, NAN};
    for (int c = 0; c < 4; ++c) {
      std::vector<pwr::Snap> s(5);
      for (int k = 0; k < 5; ++k) {
        s[k].t = 0.25 * k;
        s[k].d = {ctr("cpu", "cpu", seq[c][k])};
        s[k].d[0].aon = 1;
      }
      const double j = J(pwr::rows(s, none), "cpu");
      if (std::isnan(want[c])) {
        tk_nan(T, "pw.frz", j);
      } else {
        tk_eq(T, "pw.frz", j, want[c], 1e-12);
      }
    }
  }
  // an estimate from a rated power, with its range; a device inside the
  // package is never estimated or counted
  {
    std::vector<pwr::Snap> s(2);
    s[0].t = 0;
    s[1].t = 4;
    pwr::Dev g;
    g.id = "amd";
    g.kind = "gpu";
    g.name = "amd";
    g.bad = 1;
    g.util = 0.5;
    pwr::Dev i = g;
    i.id = i.name = "igpu";
    i.inpkg = 1;
    s[0].d = s[1].d = {g, i};
    pwr::Rated rt;
    rt.gpu = 100;
    const std::vector<pwr::Row> r = pwr::rows(s, rt);
    tk_eq(T, "pw.est", J(r, "amd"), 200, 1e-12);
    tk_eq(T, "pw.flo", find(r, "amd")->flo, 0.4, 0);
    tk_nan(T, "pw.inpkg", J(r, "igpu"));
    tk_eq(T, "pw.etot", find(r, "all devices")->fhi, 2.9, 1e-12);
    tk_nan(T, "pw.noest", J(pwr::rows(s, none), "amd"));
  }
  // snapshot text both ways
  {
    pwr::Snap a, b;
    a.t = 12.5;
    a.pcpu = 3;
    a.d = {ctr("x", "cpu", 1.25, 262143.328850), pw("y", "gpu", 7.5)};
    a.d[1].why = "a\tnote";
    tk_true(T, "pw.get", pwr::get(pwr::put(a), b));
    tk_eq(T, "pw.getc", b.d[0].wrap, 262143.32885, 0);
    tk_true(T, "pw.gets", b.d[1].why == "a note");
    tk_true(T, "pw.bad", !pwr::get("hwpwr1\t1\n", b));
  }
  // nvidia-smi, ioreg and the EMI metadata
  {
    pwr::Snap s;
    pwr::smi("41.20, 7, NVIDIA A40\n[N/A], [N/A], Tesla K80\n", s);
    tk_ieq(T, "pw.smin", static_cast<int>(s.d.size()), 2);
    tk_eq(T, "pw.smiw", s.d[0].w, 41.2, 1e-12);
    tk_eq(T, "pw.smiu", s.d[0].util, 0.07, 1e-12);
    tk_true(T, "pw.smina", s.d[1].bad == 1 && s.d[1].name == "Tesla K80");
    const std::vector<double> u = pwr::ioacc("\"PerformanceStatistics\" = {\"Device Utilization %\"=14,\"x\"=1}");
    tk_true(T, "pw.ioacc", u.size() == 1 && std::fabs(u[0] - 0.14) < 1e-12);
    std::string why;
    const double w = pwr::iobat("  \"ExternalConnected\" = No\n  \"InstantAmperage\" = 18446744073709550616\n"
                                "  \"Voltage\" = 12000\n",
                                why);
    tk_eq(T, "pw.iobat", w, 12, 1e-9);
    tk_true(T, "pw.iobatm", pwr::iobat("  \"ExternalConnected\" = Yes\n  \"Amperage\" = 0\n  \"Voltage\" = 1\n", why) < 0);
    std::vector<unsigned char> m(68, 0);
    m[66] = 2;
    const char *nm[2] = {"RAPL_Package0_PKG", "RAPL_Package0_DRAM"};
    for (const char *c : nm) {
      const std::size_t n = std::char_traits<char>::length(c);
      for (int k = 0; k < 4; ++k) m.push_back(0);
      m.push_back(static_cast<unsigned char>(2 * (n + 1)));
      m.push_back(0);
      for (std::size_t k = 0; k <= n; ++k) {
        m.push_back(static_cast<unsigned char>(c[k]));
        m.push_back(0);
      }
    }
    std::vector<std::string> names;
    std::vector<int> pwh;
    tk_true(T, "pw.emi", pwr::emimeta(m.data(), static_cast<unsigned long>(m.size()), 2, names, pwh) && names.size() == 2 &&
                             names[1] == "RAPL_Package0_DRAM" && pwh[0] == 1);
    tk_true(T, "pw.emis", !pwr::emimeta(m.data(), static_cast<unsigned long>(m.size() - 3), 2, names, pwh));
  }
}
