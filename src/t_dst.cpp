// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the distributions (src/dst.cpp) against values of
// R 4.6 (qnorm(), pnorm(), pchisq(), qt() and the others) to a relative
// 1e-12, over far tails, tiny and huge shapes and many degrees of freedom.
// Part of the self-test (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {
// got / want within 1e-12 of 1: relative also for the tiny values
void rel(tk &t, const char *id, double got, double want) { tk_eq(t, id, got / want, 1.0, 1e-12); }
}  // namespace

void t_dst(tk &t, int big) {
  (void)big;
  tk_eq(t, "qnorm.975", clab_qnorm(0.975), 1.959963984540054, 1e-15);
  tk_eq(t, "qnorm.3", clab_qnorm(0.3), -0.52440051270804067, 1e-15);
  tk_eq(t, "qnorm.995", clab_qnorm(0.995), 2.5758293035488999, 1e-15);
  tk_eq(t, "qnorm.tail", clab_qnorm(1e-10), -6.3613409024040557, 1e-15);
  tk_eq(t, "qnorm.far", clab_qnorm(1e-300), -37.04709629936120052, 1e-15);
  tk_eq(t, "qnorm.half", clab_qnorm(0.5), 0, 0);
  tk_true(t, "qnorm.ends", clab_qnorm(0) == -HUGE_VAL && clab_qnorm(1) == HUGE_VAL);
  tk_nan(t, "qnorm.out", clab_qnorm(1.5));
  tk_nan(t, "qnorm.nan", clab_qnorm(NAN));
  tk_eq(t, "pnorm.196", clab_pnorm(1.959963984540054), 0.975, 1e-15);
  tk_eq(t, "pnorm.tail", clab_pnorm(-10), 7.6198530241605269e-24, 1e-13);
  rel(t, "pchisq.384", clab_pchisq(3.84, 1, 0), 0.050043521248705224);
  rel(t, "pchisq.tiny", clab_pchisq(2, 1e-10, 0), 1.0969196720337289e-11);
  rel(t, "pchisq.big", clab_pchisq(1e6 + 3000, 1e6, 0), 0.017016772933266314);
  rel(t, "qchisq.95", clab_qchisq(0.95, 3, 1), 7.8147279032511792);
  rel(t, "qchisq.tail", clab_qchisq(1e-20, 10, 0), 118.54144259104658);
  rel(t, "qchisq.small", clab_qchisq(0.025, 0.5, 1), 5.2732025912599957e-07);
  rel(t, "pnchisq", clab_pnchisq(5, 3, 2, 1), 0.59340518008315557);
  rel(t, "pgamma.smalla", clab_pgamma(1, 1e-10, 0), 2.1938393441796981e-11);
  rel(t, "pgamma.mode", clab_pgamma(1e5 - 300, 1e5, 1), 0.17141731451450293);
  rel(t, "qgamma.tail", clab_qgamma(1e-10, 0.5, 1), 7.8539816339744674e-21);
  rel(t, "pt.2", clab_pt(-2, 3, 1), 0.069662984279421566);
  rel(t, "pt.far", clab_pt(10, 1e6, 0), 7.6393053840891504e-24);
  rel(t, "pt.heavy", clab_pt(50, 0.5, 0), 0.045352606402966393);
  rel(t, "qt.975", clab_qt(0.975, 10, 1), 2.2281388519862739);
  rel(t, "qt.tail", clab_qt(1e-10, 3, 1), -2225.7692846830946);
  rel(t, "qt.big", clab_qt(0.995, 1e7, 1), 2.5758297952037479);
  rel(t, "pbeta", clab_pbeta(0.3, 2, 5, 1), 0.57982500000000015);
  rel(t, "pbeta.tinya", clab_pbeta(1e-8, 1e-5, 10, 0), 0.00015590504833246161);
  rel(t, "pbeta.bulk", clab_pbeta(0.3008, 1e6, 7e6 / 3, 0), 0.00072172106242052582);
  rel(t, "qbeta", clab_qbeta(0.025, 3, 7, 1), 0.0748546314196918);
  rel(t, "qbeta.up", clab_qbeta(1e-6, 0.5, 40, 0), 0.25986995007384672);
  rel(t, "pf", clab_pf(3, 2, 10, 0), 0.095367431640624972);
  rel(t, "pf.big", clab_pf(5, 2, 1e6, 0), 0.0067381154487430577);
  rel(t, "qf", clab_qf(0.95, 3, 20, 1), 3.0983912121407795);
  rel(t, "ppois", clab_ppois(3, 5, 1), 0.26502591529736175);
  rel(t, "ppois.up", clab_ppois(1000, 800, 0), 4.3800283088961139e-12);
  rel(t, "qpois", clab_qpois(0.975, 5, 1), 10);
  rel(t, "qpois.big", clab_qpois(0.025, 1e5, 1), 99381);
  rel(t, "pbinom", clab_pbinom(3, 10, 0.3, 1), 0.64961071840000018);
  rel(t, "dbinom", clab_dbinom(3, 10, 0.3), 0.26682793199999993);
  rel(t, "dbinom.big", clab_dbinom(300, 1000, 0.3), 0.027521003821268392);
  tk_true(t, "dst.edges", clab_pchisq(0, 3, 1) == 0 && clab_pchisq(1, 0, 1) == 1 && clab_qt(0.5, 4, 1) == 0 &&
                              clab_qchisq(1, 2, 1) == HUGE_VAL && clab_pbeta(2, 1, 1, 1) == 1);
  tk_nan(t, "dst.out", clab_pt(1, -1, 1));
  tk_nan(t, "dst.outq", clab_qbeta(1.5, 1, 1, 1));
}
