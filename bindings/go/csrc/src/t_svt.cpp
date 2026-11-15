// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/svla.cpp, src/svt.cpp, src/svb.cpp and src/tsd.cpp: R's
// generator after set.seed(42) (the first uniform and normal R draws),
// vote agreement and the Mardia fit by hand, classical scaling of points
// on a line, a Procrustes rotation onto itself, the correlation of two
// categories over two units, the Lotka-Volterra rates of a flat series,
// R's gamma draws of a shape below 1, and the refusals.
#include <cmath>
#include <string>

#include "svk.h"
#include "tchk.h"

void t_svt(tk &T, int big) {
  (void)big;
  using svk::Node;
  // set.seed(42); runif(1) and set.seed(42); rnorm(1)
  svk::rng g(42);
  tk_eq(T, "sv.unif", g.unif(), 0.914806043496355, 1e-15);
  svk::rng h(42);
  tk_eq(T, "sv.norm", h.norm(), 1.37095844714667, 1e-14);
  // two legislators agree on 1 of 2 votes, the third is absent from one
  Node a = svk::list();
  Node x = svk::vec({1, 1, 0, 0, 1, NAN});
  x.dim = {3, 2};
  a.add("x", x);
  Node out;
  std::string err;
  tk_ieq(T, "sv.agr", svk::run("votagr", a, out, err), 1);
  tk_eq(T, "sv.agr12", out.at("agree")->d[3], 0.5, 0);
  tk_eq(T, "sv.agr13", out.at("agree")->d[6], 0, 0);
  // Mardia: positive eigenvalues over the sum of absolute values
  Node b = svk::list();
  b.add("ev", svk::vec({4, 2, -2}));
  tk_ieq(T, "sv.fit", svk::run("mdsfit", b, out, err), 1);
  tk_eq(T, "sv.fit1", out.at("fit")->d[0], 0.5, 0);
  tk_eq(T, "sv.cum", out.at("cumfit")->d[2], 0.75, 0);
  // three points at 0, 1, 3 on a line: one dimension holds them exactly
  Node c = svk::list();
  Node d = svk::vec({0, 1, 3, 1, 0, 2, 3, 2, 0});
  d.dim = {3, 3};
  c.add("D", d);
  c.add("ndim", svk::num(1));
  tk_ieq(T, "sv.cmd", svk::run("cmdsf", c, out, err), 1);
  tk_eq(T, "sv.cmds", out.at("stress")->d[0], 0, 1e-12);
  tk_eq(T, "sv.cmdf", out.at("fit")->d[0], 1, 1e-12);
  // a configuration rotated onto itself
  Node p = svk::list();
  Node m = svk::vec({0, 1, 2, 0, 1, 0});
  m.dim = {3, 2};
  p.add("X", m);
  p.add("target", m);
  tk_ieq(T, "sv.pro", svk::run("procrs", p, out, err), 1);
  tk_eq(T, "sv.prom", out.at("mse")->d[0], 0, 1e-24);
  tk_eq(T, "sv.pror", out.at("rotation")->d[0], 1, 1e-14);
  // two categories over two units: counts (2, 1) and (1, 2), r = -1
  Node u = svk::list();
  Node th, fr;
  th.k = Node::STR;
  th.s = {"a", "a", "b"};
  fr.k = Node::STR;
  fr.s = {"a", "b", "b"};
  u.add("theft", th);
  u.add("fraud", fr);
  Node cu = svk::list();
  cu.add("units", u);
  tk_ieq(T, "sv.cor", svk::run("catcor", cu, out, err), 1);
  tk_eq(T, "sv.corr", out.at("r")->d[2], -1, 1e-15);
  tk_eq(T, "sv.corc", out.at("counts")->d[2], 2, 0);
  // a flat series: no change, so alpha = gamma = 0.05, the period 2 pi / 0.05
  Node lv = svk::list();
  lv.add("y", svk::vec({10, 10, 10, 10, 10}));
  tk_ieq(T, "sv.lv", svk::run("lvfit", lv, out, err), 1);
  tk_eq(T, "sv.lva", out.at("alpha")->d[0], 0.05, 1e-15);
  tk_eq(T, "sv.lvb", out.at("beta")->d[0], 0.005, 1e-15);
  tk_eq(T, "sv.lvp", out.at("period")->d[0], 40 * 3.141592653589793, 1e-10);
  // refusals
  Node e = svk::list();
  e.add("D", svk::vec({1, 2}));
  e.add("ndim", svk::num(1));
  tk_ieq(T, "sv.sq", svk::run("cmdsf", e, out, err), 0);
  tk_ieq(T, "sv.none", svk::run("nosuch", e, out, err), 0);
  Node l3 = svk::list();
  l3.add("y", svk::vec({1, 2, 3}));
  tk_ieq(T, "sv.lv3", svk::run("lvfit", l3, out, err), 0);
  Node c1 = svk::list();
  Node u1 = svk::list();
  u1.add("theft", th);
  c1.add("units", u1);
  tk_ieq(T, "sv.cor1", svk::run("catcor", c1, out, err), 0);
  // a category given as numbers, not the units' names (catcor() in R
  // makes every category text first)
  Node c2 = svk::list();
  Node u2 = svk::list();
  u2.add("theft", th);
  u2.add("fraud", svk::vec({1, 2}));
  c2.add("units", u2);
  tk_ieq(T, "sv.cornum", svk::run("catcor", c2, out, err), 0);
  tk_true(T, "sv.cornumm", err == "`units` must hold character vectors.");
  // set.seed(42); rgamma(3, shape = 0.5, scale = 2): a shape below 1 by
  // Ahrens and Dieter's GS; a zero shape or scale is a point mass at 0
  svk::rng gg(42);
  tk_eq(T, "sv.gam1", gg.gamma(0.5, 2), 1.52103366962146858, 1e-14);
  tk_eq(T, "sv.gam2", gg.gamma(0.5, 2), 0.58739504553944721, 1e-14);
  tk_eq(T, "sv.gam3", gg.gamma(0.5, 2), 3.73268704722661271, 1e-14);
  tk_true(T, "sv.gam0", gg.gamma(0, 2) == 0 && gg.gamma(2, 0) == 0 && std::isnan(gg.gamma(-1, 2)));
  // dnorm() far beyond its underflow, and a one-string node
  tk_eq(T, "sv.dnfar", svk::dnorm(3e154), 0, 0);
  const Node sn = svk::str("theft");
  tk_true(T, "sv.str", sn.k == Node::STR && sn.s.size() == 1 && sn.s[0] == "theft");
}
