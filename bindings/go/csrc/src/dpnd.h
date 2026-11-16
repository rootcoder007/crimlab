// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The work of gmmpnl(), wsctst() and ccepnl() between their arguments and
// their results, shared by the R binding (bdpn.cpp) and clab_jrun()
// (cli/jenv.inc): the design of difference and system GMM, the calls of
// the kernels (dpn.cpp), the refusals in the words both give, and the
// results as plain arrays. Nothing here knows about R or JSON; each
// function returns an empty message or the refusal, and CLAB_EINTR when
// the host asked to stop.
#ifndef CLAB_DPND_H
#define CLAB_DPND_H

#include <string>
#include <vector>

namespace dpnd {

// A panel: the outcome y (n), the unit codes u (0..nu-1, labels ulev),
// the period codes t (0..nt-1) and the periods' values tv, sorted.
struct pnl {
  int n = 0;
  const double *y = nullptr;
  std::vector<int> u, t;
  int nu = 0, nt = 0;
  std::vector<double> tv;
  std::vector<std::string> ulev;
};

// The columns of a numeric matrix of n rows, with their names.
struct xm {
  std::vector<double> v;  // n x p, column-major
  std::vector<std::string> nm;
  int p = 0;
};

// The options of gmmpnl(), read and checked one by one by the caller.
struct spec {
  int ylag = 1, sys = 1, fod = 0, coll = 0, cons = 1, tdum = 1, hm = 0;
  int sig = 0;                  // the sigma^2 H blocks of the AR tests are needed
  std::vector<int> gv, lo, hi;  // the GMM-style variables (-1: y) and lags
  std::vector<int> iv;          // the columns of x that are IV-style
  int steps = 1, rob = 1, wc = -1, small = 0, nar = 2;  // wc -1: not given
  double level = 0.95;
};

// A table of results: one row per name.
struct tab {
  std::vector<std::string> name;
  std::vector<double> est, se, stat, df, df2, p, lo, hi;
};

// A gmmpnl() result.
struct gout {
  tab coef, tests, wald, xwald;  // xwald: the joint test as xtabond2 prints it
  bool tstat = false;  // t (small) rather than z
  double n = 0, groups = 0, inst = 0, omin = 0, oavg = 0, omax = 0, steps = 0, k = 0, sigma2 = 0;
  std::vector<double> vcov, resid;  // resid: y - x b in levels, NaN off them
  std::vector<std::string> inst_names, dropped, dummies, note;
  std::string type, trans;
};

// The checks of spec that do not depend on the data (`wc` against
// `robust` and `steps`); then the design, the kernel and the result.
std::string gcheck(const spec &sp);
std::string gfit(const pnl &pn, const xm &x, const xm &z, spec sp, gout &o, int &code);

// A wsctst() result: stat (F, df1, df2, p, coef, se, pairs, clusters,
// ndiff), and the differenced regression (b, se) when x has columns.
struct wout {
  double stat[9];
  std::vector<double> b, se;
};
std::string wfit(const pnl &pn, const xm &x, int cons, int raw, wout &o, int &code);

// The dynamic and other options of ccepnl(): ylag lags of y among the
// regressors; lags of the averages (-1: the default of Chudik and Pesaran
// 2015); bias 0 none, 1 the half-panel jackknife, 2 recursive mean
// adjustment; pool -1 (the default), 0 or 1; pvar 0 Pesaran's pooled
// variance, 1 that of Brown, Schmidt and Wooldridge (2026).
struct copt {
  int ylag = 0, lags = -1, bias = 0, pool = -1, pvar = 0;
};

// A ccepnl() result.
struct cres {
  tab mg, pooled, mg0;  // mg0: the uncorrected mean group estimate (bias)
  std::vector<double> bi, vmg, vp, w, vmg0;  // bi: units x k, column-major
  std::vector<std::string> nm;               // the terms: lags of y, then x
  double units = 0, periods = 0, k = 0, h = 0, n = 0, lags = 0, ylag = 0;
  bool pool = true, bias = false;
  double hb[3] = {0, 0, 0};  // heterogeneity bias: stat, df, p (pool)
  double nf[3] = {0, 0, 0};  // the number of factors: ER, threshold, averages
};
// w: empty or one positive weight per row; avg: all (when useall) or the
// names among "y" and the columns of x.
std::string cfit(const pnl &pn, const xm &x, const xm &fac, int cons, const std::vector<double> &w,
                 bool useall, const std::vector<std::string> &avg, const copt &co, double level,
                 cres &o, int &code);

}  // namespace dpnd

#endif
