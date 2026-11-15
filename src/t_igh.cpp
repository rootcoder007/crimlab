// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of the Indigenous health, published-table and census functions
// (igh.h): R's generator after set.seed(1), limits from a coefficient of
// variation, a gap against a reference, the generations of exposure and
// the Community Well-Being score, each by hand.
#include <climits>
#include <cmath>
#include <string>
#include <vector>

#include "crimlab.h"
#include "igk.h"
#include "tchk.h"

namespace {

igh::V nums(const std::vector<double> &x) { return igh::dv(x); }

// a data frame as clab_jrun() reads one: an object of columns, no class
igh::V jdf(const std::vector<std::string> &nm, std::vector<igh::V> cols) {
  igh::V f = igh::frame(nm, std::move(cols));
  f.df = false;
  f.cls.clear();
  f.js = true;
  for (igh::V &c : f.a) c.js = true;
  return f;
}

// an object of one number per name, as JSON gives a named vector
igh::V jnamed(const std::vector<std::string> &nm, const std::vector<double> &x) {
  igh::V o;
  o.t = igh::V::LIST;
  o.js = true;
  o.nm = nm;
  for (double v : x) o.a.push_back(igh::dv1(v));
  return o;
}

int run(const char *fn, const std::vector<std::string> &an, const std::vector<igh::V> &av, igh::V &out) {
  std::string err;
  std::vector<std::string> warn;
  return igh::call(fn, an, av, out, err, warn);
}

}  // namespace

void t_igh(tk &T, int big) {
  (void)big;
  // set.seed(1); runif(1); rnorm(1) in R
  igh::Rng g(1);
  tk_eq(T, "ig.unif", g.unif(), 0.26550866314209998, 1e-15);
  igh::Rng h(1);
  tk_eq(T, "ig.norm", h.norm(), -0.62645381074233242, 1e-14);
  // limfrom(c(24.1, 31.5), cv = c(8, 12)): s = est cv / 100
  std::vector<std::string> an{"est", "cv"};
  std::vector<igh::V> av{nums({24.1, 31.5}), nums({8, 12})};
  igh::V out;
  std::string err;
  std::vector<std::string> warn;
  tk_ieq(T, "ig.lim", igh::call("limfrom", an, av, out, err, warn), CLAB_OK);
  const double z = clab_qnorm(0.975);
  tk_eq(T, "ig.limlo", out.a[1].d[0], 24.1 - z * 24.1 * 8 / 100, 1e-13);
  tk_eq(T, "ig.limse", out.a[3].d[1], 31.5 * 12 / 100, 1e-15);
  // indgap(): the difference of First Nations and the reference
  igh::V d = igh::frame({"g", "e", "l", "u"}, {igh::sv({"First Nations", "Non-Indigenous"}), nums({45.8, 61.5}),
                                               nums({43.9, 61}), nums({47.7, 62})});
  an = {"data", "grp", "est", "lower", "upper"};
  av = {d, igh::sv1("g"), igh::sv1("e"), igh::sv1("l"), igh::sv1("u")};
  tk_ieq(T, "ig.gap", igh::call("indgap", an, av, out, err, warn), CLAB_OK);
  tk_eq(T, "ig.gapdif", out.a[6].d[0], 45.8 - 61.5, 1e-13);
  tk_eq(T, "ig.gapratio", out.a[3].d[0], 45.8 / 61.5, 1e-15);
  // a missing column is refused in words
  av[1] = igh::sv1("zz");
  tk_ieq(T, "ig.gapcol", igh::call("indgap", an, av, out, err, warn), CLAB_EARG);
  // genexpo(): TRUE and NA is unknown by default
  an = {"parent", "gparent"};
  av = {igh::lv({1, 0, 1, INT_MIN}), igh::lv({1, 0, 0, 1})};
  tk_ieq(T, "ig.gen", igh::call("genexpo", an, av, out, err, warn), CLAB_OK);
  tk_ieq(T, "ig.gen2", out.i[0], 2);
  tk_ieq(T, "ig.gen1", out.i[2], 1);
  tk_ieq(T, "ig.genna", out.i[3] == INT_MIN, 1);
  // cwbscr(): income at the upper bound scores 100
  an = {"incpercap", "highsch20p", "uni25p", "notcro", "nomajrep", "partic", "empl"};
  av = {nums({75000}), nums({60}), nums({30}), nums({80}), nums({90}), nums({70}), nums({90})};
  tk_ieq(T, "ig.cwb", igh::call("cwbscr", an, av, out, err, warn), CLAB_OK);
  tk_eq(T, "ig.cwbinc", out.a[0].d[0], 100, 1e-12);
  tk_eq(T, "ig.cwbedu", out.a[1].d[0], 2.0 / 3 * 60 + 1.0 / 3 * 30, 1e-13);
  tk_eq(T, "ig.cwb", out.a[4].d[0], (100 + 50 + 85 + 80) / 4.0, 1e-12);
  // R's round(): 0.15 to one digit is 0.1 (0.15 is below 0.15 in binary)
  tk_eq(T, "ig.round", igh::rround(0.15, 1), 0.1, 0);
  tk_eq(T, "ig.round2", igh::rround(2.675, 2), 2.67, 0);
  // the shapes clab_jrun() gives: frames and named vectors as objects
  const igh::V wp = jdf({"group", "indicator", "side", "estimate", "lower", "upper", "usewithcau"},
                        {igh::sv({"A", "A"}), igh::sv({"w", "d"}), igh::sv({"wellness", "distress"}),
                         nums({60, 20}), nums({55, 15}), nums({65, 25}), igh::lv({0, 1})});
  tk_ieq(T, "ig.jfrm", run("frmaud", {"x"}, {wp}, out), CLAB_OK);
  tk_eq(T, "ig.jshr", out.a[1].d[0], 0.5, 0);
  tk_ieq(T, "ig.jdc", run("dualcont", {"x", "well", "distress"}, {wp, igh::sv1("w"), igh::sv1("d")}, out), CLAB_OK);
  tk_eq(T, "ig.jdcr", out.a[7].d[0], 3, 1e-15);
  tk_ieq(T, "ig.jdcno", run("dualcont", {"x", "well", "distress"}, {jdf({"a"}, {nums({1})}), igh::sv1("w"),
                                                                    igh::sv1("d")}, out), CLAB_EARG);
  const igh::V cd = jdf({"g", "i", "e", "l", "u"}, {igh::sv({"A", "A"}), igh::sv({"w", "d"}), nums({60, 20}),
                                                    nums({55, 15}), nums({65, 25})});
  tk_ieq(T, "ig.jci", run("commidx", {"data", "grp", "ind", "est", "lower", "upper", "w"},
                          {cd, igh::sv1("g"), igh::sv1("i"), igh::sv1("e"), igh::sv1("l"), igh::sv1("u"),
                           jnamed({"w", "d"}, {2, -1})}, out), CLAB_OK);
  tk_eq(T, "ig.jcidx", out.a[1].d[0], (2 * 60.0 - 20) / 3, 1e-13);
  const igh::V hb = jdf({"g", "d", "e", "l", "u"}, {igh::sv({"A", "A"}), igh::sv({"a", "b"}), nums({60, 40}),
                                                    nums({55, 35}), nums({65, 45})});
  tk_ieq(T, "ig.jhb", run("holbal", {"data", "grp", "domain", "est", "lower", "upper", "w", "reps"},
                          {hb, igh::sv1("g"), igh::sv1("d"), igh::sv1("e"), igh::sv1("l"), igh::sv1("u"),
                           jnamed({"a", "b"}, {1, 3}), nums({100})}, out), CLAB_OK);
  tk_eq(T, "ig.jhblv", out.a[1].d[0], 100 * (0.25 * 0.6 + 0.75 * 0.4), 1e-12);
  igh::V dom;
  dom.t = igh::V::LIST;
  dom.js = true;
  dom.nm = {"w"};
  dom.a.push_back(igh::sv1("spirit"));
  tk_ieq(T, "ig.jwp", run("welprof", {"data", "grp", "ind", "est", "lower", "upper", "well", "distress", "domain"},
                          {cd, igh::sv1("g"), igh::sv1("i"), igh::sv1("e"), igh::sv1("l"), igh::sv1("u"),
                           igh::sv1("w"), igh::sv1("d"), dom}, out), CLAB_OK);
  tk_ieq(T, "ig.jwpdom", out.a[3].s[0] == "spirit", 1);
  const igh::V pf = jdf({"geo", "dguid", "residence", "gender", "age", "group", "chaid", "characteristic", "parid",
                         "rootid", "sample", "value", "quality"},
                        {igh::sv({"A", "A"}), igh::sv({"X", "X"}), igh::sv({"Total", "Total"}),
                         igh::sv({"Total", "Total"}), igh::sv({"Total", "Total"}), igh::sv({"G", "G"}), nums({1, 2}),
                         igh::sv({"T", "c"}), nums({igh::na(), 1}), nums({1, 1}), igh::sv({"100%", "100%"}),
                         nums({10, 4}), igh::sv({"ok", "ok"})});
  tk_ieq(T, "ig.jps", run("profshr", {"x", "num"}, {pf, nums({2})}, out), CLAB_OK);
  tk_eq(T, "ig.jpsest", out.a[8].d[0], 40, 1e-12);
  tk_ieq(T, "ig.jpsno", run("profshr", {"x", "num"}, {jdf({"a"}, {nums({1})}), nums({2})}, out), CLAB_EARG);
  const igh::V sc = jdf({"GEO", "VALUE", "STATUS"}, {igh::sv({"C", "C"}), nums({1, 2}), igh::sv({"", "x"})});
  tk_ieq(T, "ig.jsc", run("sctidy", {"x"}, {sc}, out), CLAB_OK);
  tk_nan(T, "ig.jscx", out.a[1].d[1]);
  // a function where data belong is refused in words
  igh::V f;
  f.fn = true;
  tk_ieq(T, "ig.fn", run("trcprog", {"done"}, {f}, out), CLAB_EARG);
  // a name no function has, and exptrnd() with its default levels
  tk_ieq(T, "ig.nofn", run("nosuch", {}, {}, out), CLAB_EARG);
  tk_ieq(T, "ig.dose", run("exptrnd", {"cases", "n"}, {nums({1, 4}), nums({10, 10})}, out), CLAB_OK);
  tk_eq(T, "ig.dose1", out.a[0].d[1], 1, 0);
  // round() of R below zero and to tens; rpois() of R below 10 (set.seed(2);
  // rpois(12, 3.5))
  tk_eq(T, "ig.rneg", igh::rround(-2.675, 2), -2.67, 0);
  tk_eq(T, "ig.rten", igh::rround(1234.5, -2), 1200, 0);
  igh::Rng p(2);
  const double pv[12] = {2, 4, 4, 2, 7, 7, 1, 5, 3, 4, 4, 2};
  int same = 0;
  for (double v : pv) same += p.pois(3.5) == v;
  tk_ieq(T, "ig.pois", same, 12);
}
