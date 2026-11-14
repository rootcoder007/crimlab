// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Self-exciting point processes and series (src/btm.cpp): the residuals
// and simulation of a Hawkes process, and change points in a series.

const char *const hkinds[] = {"exp", "weibull", "lomax"};
const char *const hbases[] = {"const", "cycle"};
int hnpar(int kind, int base) { return (base ? 4 : 1) + 1 + (kind ? 2 : 1); }

// pos1(): one positive number.
bool pos1(Ctx &c, const char *nm, double def, double &v) {
  Dv d;
  Val d1;
  CK(num(c, std::isnan(def) ? c[nm] : arg(c, nm, d1, def), nm, d, 1, 0.0, true));
  v = d[0];
  return true;
}
// evtimes(): event times from 0 to hor, sorted, at least nmin of them.
bool evtimes(Ctx &c, double hor, int nmin, std::vector<double> &t) {
  Dv tv;
  CK(num(c, c["times"], "times", tv));
  t = tv.v;
  const int n = static_cast<int>(t.size());
  if (n < nmin) return c.fail("`times` must hold at least %d events, not %d.", nmin, n);
  for (std::size_t i = 0; i < t.size(); ++i) {
    if (t[i] < 0 || t[i] > hor) return c.fail("`times` must lie between 0 and `hor`.");
    if (i > 0 && t[i] < t[i - 1]) return c.fail("`times` must be sorted from earliest to latest.");
  }
  return true;
}
// parv(): the parameters of a kernel and a baseline.
bool parv(Ctx &c, const char *nm, int kind, int base, std::vector<double> &th) {
  const int np = hnpar(kind, base);
  Dv v;
  CK(num(c, c[nm], nm, v, np, -HUGE_VAL));
  th = v.v;
  if (!(th[static_cast<std::size_t>(np - (kind ? 3 : 2))] >= 0)) return c.fail("`%s`: eta must be at least 0.", nm);
  for (int j = np - (kind ? 2 : 1); j < np; ++j) {
    if (!(th[static_cast<std::size_t>(j)] > 0)) return c.fail("`%s`: the kernel parameters must be positive.", nm);
  }
  return true;
}

bool hwkres(Ctx &c, Val &out) {
  double h = 0.0, pr = 0.0;
  CK(pos1(c, "hor", NaN, h));
  CK(pos1(c, "per", 365.25, pr));
  int kd = 0, bd = 0;
  CK(pickd(c, "kernel", hkinds, 3, 0, kd));
  CK(pickd(c, "base", hbases, 2, 0, bd));
  std::vector<double> t, th;
  CK(evtimes(c, h, 1, t));
  CK(parv(c, "par", kd, bd, th));
  const std::size_t n = t.size();
  std::vector<double> tau(n), work(n < 100 ? 3 * (2 * n + 1) * (2 * n + 1) + n : n), u(n);
  double ks[2];
  clab_hwkres(t.data(), static_cast<int>(n), h, kd, bd, pr, th.data(), tau.data(), ks, work.data());
  for (std::size_t i = 0; i < n; ++i) u[i] = -std::expm1(-tau[i]);
  for (double &e : tau) e = na(e);
  for (double &e : u) e = na(e);
  out = vobj();
  put(out, "tau", vauto(tau));
  put(out, "u", vauto(u));
  put(out, "d", vnum(ks[0]));
  put(out, "p", vnum(ks[1]));
  put(out, "n", vnum(static_cast<double>(n)));
  return true;
}

bool hwksim(Ctx &c, Val &out) {
  double h = 0.0, pr = 0.0;
  CK(pos1(c, "hor", NaN, h));
  CK(pos1(c, "per", 365.25, pr));
  int kd = 0, bd = 0;
  CK(pickd(c, "kernel", hkinds, 3, 0, kd));
  CK(pickd(c, "base", hbases, 2, 0, bd));
  std::vector<double> th;
  CK(parv(c, "par", kd, bd, th));
  int cap = 0;
  CK(wholed(c, "nmax", 1e5, 1, cap));
  double w[4], sg = 0;
  const bool sp = c.has("win");
  if (sp) {
    CK(win(c, c["win"], w));
    CK(pos1(c, "sigma", NaN, sg));
    th[0] += std::log((w[1] - w[0]) * (w[3] - w[2]));
  }
  unsigned int sd = 0;
  CK(seed(c, sd));
  const std::size_t sc = static_cast<std::size_t>(cap);
  std::vector<double> t(sc), x(sp ? sc : 0), y(sp ? sc : 0), work(3 * sc);
  std::vector<int> gen(sc), iwork(2 * sc);
  int n = 0;
  if (clab_hwksim(h, kd, bd, pr, th.data(), sp ? w : nullptr, sg, sd, cap, t.data(), sp ? x.data() : nullptr,
                  sp ? y.data() : nullptr, gen.data(), &n, work.data(), iwork.data()) != CLAB_OK) {
    return c.fail("more than `nmax` = %d events: raise it, or lower eta (at 1 or above the process explodes).", cap);
  }
  const std::size_t N = static_cast<std::size_t>(n);
  Frame f;
  f.num("t", t.data(), N);
  if (sp) {
    f.num("x", x.data(), N);
    f.num("y", y.data(), N);
  }
  f.ints("gen", gen.data(), N);
  out = f.o;
  return true;
}

bool stpcp(Ctx &c, Val &out) {
  Dv yv;
  CK(num(c, c["y"], "y", yv));
  const int n = static_cast<int>(yv.n());
  int m = 0, np = 0;
  CK(wholed(c, "minseg", 2, 1, m));
  CK(wholed(c, "nperm", 9999, 1, np));
  if (n < 2 * m) return c.fail("`y` must hold at least 2 `minseg` = %d values, not %d.", 2 * m, n);
  unsigned int sd = 0;
  CK(seed(c, sd));
  double st[12];
  std::vector<double> w(static_cast<std::size_t>(n));
  std::vector<int> iw(static_cast<std::size_t>(n));
  clab_stpcp(yv.p(), n, m, np, sd, st, w.data(), iw.data());
  const int k[2] = {static_cast<int>(st[0]), static_cast<int>(st[7])};
  const double be[2] = {st[1], st[10]}, af[2] = {st[2], st[11]}, sa[2] = {st[3], st[8]}, pv[2] = {st[4], st[9]};
  Frame f;
  f.str("test", {"step", "pettitt"});
  f.ints("k", k, 2);
  for (const auto &e : {std::make_pair("before", be), std::make_pair("after", af), std::make_pair("stat", sa),
                        std::make_pair("p", pv)}) {
    std::vector<double> v(e.second, e.second + 2);
    for (double &x : v) x = na(x);
    f.num(e.first, v);
  }
  if (c.has("x")) {
    // stpcp() in R/trn.R: the periods at the change points
    const Val &x = c["x"];
    if (vlen(x) != static_cast<std::size_t>(n)) return c.fail("`x` must have one period per value, %d.", n);
    Val at = varr();
    for (int q : k) at.a.push_back(x.t == Val::ARR || x.t == Val::OBJ ? x.a[static_cast<std::size_t>(q) - 1] : x);
    f.add("at", at);
  }
  out = vobj();
  put(out, "tab", f.o);
  put(out, "nperm", vnum(st[5]));
  put(out, "exact", vbool(st[6] != 0));
  put(out, "n", vnum(n));
  return true;
}
