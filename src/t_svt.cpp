// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/svla.cpp, src/svt.cpp, src/svb.cpp and src/tsd.cpp: R's
// generator after set.seed(42) (the first uniform and normal R draws),
// vote agreement and the Mardia fit by hand, classical scaling of points
// on a line, a Procrustes rotation onto itself and the refusals.
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
  // refusals
  Node e = svk::list();
  e.add("D", svk::vec({1, 2}));
  e.add("ndim", svk::num(1));
  tk_ieq(T, "sv.sq", svk::run("cmdsf", e, out, err), 0);
  tk_ieq(T, "sv.none", svk::run("nosuch", e, out, err), 0);
}
