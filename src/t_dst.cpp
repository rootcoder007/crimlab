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

void t_dstbig(tk &t);

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
  // no successes: dbinom(0, 10, 0.05), dbinom(0, 10, 0.3), dbinom(0, 0, 0.4)
  rel(t, "dbinom.0s", clab_dbinom(0, 10, 0.05), 0.598736939238378896);
  rel(t, "dbinom.0l", clab_dbinom(0, 10, 0.3), 0.028247524900000005);
  tk_eq(t, "dbinom.00", clab_dbinom(0, 0, 0.4), 1, 0);
  tk_nan(t, "dbinom.n", clab_dbinom(1, 2.5, 0.5));
  tk_nan(t, "pbinom.n", clab_pbinom(1, 2.5, 0.5, 1));
  // the limits of the beta at a zero or infinite shape: point masses at
  // 0, 1, or both (pbeta() and qbeta() in their lower tails)
  tk_true(t, "pbeta.lim", clab_pbeta(0.3, 0, 2, 1) == 1 && clab_pbeta(0.3, 2, 0, 1) == 0 &&
                              clab_pbeta(0.3, 0, 0, 1) == 0.5 && clab_pbeta(0.3, HUGE_VAL, HUGE_VAL, 1) == 0 &&
                              clab_pbeta(0.7, HUGE_VAL, HUGE_VAL, 1) == 1 && clab_pbeta(0.3, 0, 2, 0) == 0);
  tk_true(t, "qbeta.lim", clab_qbeta(0.3, 0, 0, 1) == 0 && clab_qbeta(0.7, 0, 0, 1) == 1 &&
                              clab_qbeta(0.3, 0, 2, 1) == 0 && clab_qbeta(0.3, 2, 0, 1) == 1 &&
                              clab_qbeta(0.3, HUGE_VAL, HUGE_VAL, 1) == 0.5);
  // the point masses as R gives them, both tails: qbeta(c(0.3, 0.5, 0.7),
  // 0, 0) is 0, 0.5, 1 with lower.tail TRUE and FALSE alike; a = 0 is 0
  // and b = 0 is 1 inside (0, 1), the ends of p as for any shapes
  tk_true(t, "qbeta.00", clab_qbeta(0.5, 0, 0, 1) == 0.5 && clab_qbeta(0.5, 0, 0, 0) == 0.5 &&
                             clab_qbeta(0.3, 0, 0, 0) == 0 && clab_qbeta(0.7, 0, 0, 0) == 1 &&
                             clab_qbeta(0, 0, 0, 0) == 1 && clab_qbeta(1, 0, 0, 0) == 0);
  tk_true(t, "qbeta.0x", clab_qbeta(0.7, 0, 2, 1) == 0 && clab_qbeta(0.7, 0, 2, 0) == 0 &&
                             clab_qbeta(1, 0, 2, 1) == 1 && clab_qbeta(0, 0, 2, 0) == 1 &&
                             clab_qbeta(0.3, 2, 0, 0) == 1 && clab_qbeta(0, 2, 0, 1) == 0 &&
                             clab_qbeta(1, 2, 0, 0) == 0);
  // qgamma(c(0, 0.3, 1), 0) is 0, 0, Inf; upper tail Inf, 0, 0
  tk_true(t, "qgamma.0", clab_qgamma(0.3, 0, 1) == 0 && clab_qgamma(1, 0, 1) == HUGE_VAL &&
                             clab_qgamma(0, 0, 0) == HUGE_VAL && clab_qgamma(1, 0, 0) == 0 &&
                             clab_qchisq(1, 0, 1) == HUGE_VAL);
  // t beyond 1e50: the tail from its logarithm, pt(-1e60, 3), pt(-1e60, 0.5)
  rel(t, "pt.far", clab_pt(-1e60, 3, 1), 1.1026577908435548e-180);
  rel(t, "pt.farup", clab_pt(1e60, 3, 0), 1.1026577908435548e-180);
  rel(t, "pt.farh", clab_pt(-1e60, 0.5, 1), 3.2070097541422669e-31);
  // F with infinite numerator degrees of freedom: pf(2, Inf, 5) and both
  // infinite, a point mass at 1
  rel(t, "pf.inf", clab_pf(2, HUGE_VAL, 5, 1), 0.77649507112332272);
  rel(t, "pf.infup", clab_pf(2, HUGE_VAL, 5, 0), 0.22350492887667725);
  tk_true(t, "pf.infs", clab_pf(2, HUGE_VAL, HUGE_VAL, 1) == 1 && clab_pf(1, HUGE_VAL, HUGE_VAL, 1) == 0.5 &&
                            clab_pf(0.5, HUGE_VAL, HUGE_VAL, 1) == 0);
  // F quantiles for a huge denominator: qf(0.9, 3, 1e6), and 1 for both
  tk_eq(t, "qf.infs", clab_qf(0.9, HUGE_VAL, HUGE_VAL, 1), 1, 0);
  rel(t, "qf.big", clab_qf(0.9, 3, 1e6, 1), 2.08379621039010798);
  rel(t, "qf.bigup", clab_qf(0.9, 3, 1e6, 0), 0.19479145805172773);
  // Poisson quantiles whose normal start falls short: qpois(0.05, 3),
  // qpois(0.2, 3), qpois(0.8, 3, lower.tail = FALSE)
  tk_eq(t, "qpois.up1", clab_qpois(0.05, 3, 1), 1, 0);
  tk_eq(t, "qpois.up2", clab_qpois(0.2, 3, 1), 2, 0);
  tk_eq(t, "qpois.up3", clab_qpois(0.8, 3, 0), 2, 0);
  t_dstbig(t);
}

// Large parameters, against references of 34 to 40 digits from mpmath
// (dev/ref/dstbig.py): the noncentral chi-square by an integral over the
// normal part (another representation from src/dst.cpp's), its
// noncentrality bounds by the secant method on it (df 1 in closed form),
// the incomplete gamma and beta ratios by their defining integrals, the
// noncentral F by its Poisson mixture. Each value takes at most a few
// milliseconds, whatever the parameter's size.
void t_dstbig(tk &t) {
  rel(t, "pnchisq.int", clab_pnchisq(1002, 2, 1000, 1), 0.5063033675594192372927274);
  rel(t, "pnchisq.intlo", clab_pnchisq(685.614159608872569151571951806545257568359375, 2, 1000, 1),
      2.436332487074393327752203e-8);
  rel(t, "pnchisq.e15", clab_pnchisq(999999873508894.625, 1, 1e15, 1), 0.02275012941425419019968595);
  rel(t, "pnchisq.e15up", clab_pnchisq(1000000316227767, 1, 1e15, 0), 2.86652136445840593867498e-7);
  rel(t, "pnchisq.e8", clab_pnchisq(100040001.00010000169277191162109375, 1, 1e8, 0), 0.02275822991822743518459779);
  tk_nan(t, "pnchisq.inf", clab_pnchisq(5, 3, HUGE_VAL, 1));
  // the bounds of Cramer's V for a 2 x 2 table of 1e8 and of 1e15 counts
  rel(t, "ncpq.e8lo", clab_ncpq(4e6, 1, 0.975), 3992163.985520660478702051);
  rel(t, "ncpq.e8hi", clab_ncpq(4e6, 1, 0.025), 4007843.697396980910972984);
  rel(t, "ncpq.e15lo", clab_ncpq(4e13, 1, 0.975), 39999975208202.54927636026);
  rel(t, "ncpq.e15hi", clab_ncpq(4e13, 1, 0.025), 40000024791805.13364128564);
  tk_true(t, "ncpq.edge", clab_ncpq(0.5, 4, 0.975) == 0 && clab_ncpq(0, 1, 0.5) == 0);
  tk_nan(t, "ncpq.out", clab_ncpq(5, 1, 1));
  tk_nan(t, "ncpq.inf", clab_ncpq(HUGE_VAL, 1, 0.5));
  // the bounds fall as the level of the lower tail rises (a wider interval)
  tk_true(t, "ncpq.mono", clab_ncpq(4e6, 4, 0.995) < clab_ncpq(4e6, 4, 0.975) &&
                              clab_ncpq(4e6, 4, 0.025) < clab_ncpq(4e6, 4, 0.005) &&
                              clab_ncpq(500, 99, 0.995) < clab_ncpq(500, 99, 0.975));
  rel(t, "pgamma.e8", clab_pgamma(99970000.0, 1e8, 1), 0.00134871644916155059179);
  rel(t, "pgamma.e16", clab_pgamma(9999999700000000.0, 1e16, 1), 0.001349897913447472645822);
  rel(t, "pgamma.e16up", clab_pgamma(9999999700000000.0, 1e16, 0), 0.9986501020865525273542);
  tk_true(t, "pgamma.e300", std::fabs(clab_pgamma(1e300, 1e300, 1) - 0.5) < 1e-15);
  rel(t, "pbeta.e8", clab_pbeta(0.49992928932188135, 1e8, 1e8, 1), 0.02275013154326005409058);
  rel(t, "pbeta.e12", clab_pbeta(0.4999992928932188, 1e12, 1e12, 1), 0.02275013194857566620789);
  rel(t, "pbeta.e8up", clab_pbeta(0.01, 1990001, 198010000, 0), 5.714803963440756959975e-13);
  tk_true(t, "pbeta.e300", std::fabs(clab_pbeta(0.5, 1e300, 1e300, 1) - 0.5) < 1e-15 &&
                               clab_pbeta(0.1, 5e307, 28.5, 1) == 0 && clab_pbeta(0.1, 28.5, 5e307, 0) == 0);
  tk_true(t, "pgamma.max", clab_pgamma(1.7e308, 28.5, 0) == 0);
  rel(t, "pnf", clab_pnf(3.1588427192606483, 2, 57, 6, 0), 0.560173535744839529465);
  rel(t, "pnf.60", clab_pnf(3.1588427192606483, 2, 57, 60, 0), 0.9999999008156065140529);
  rel(t, "pnf.60lo", clab_pnf(3.1588427192606483, 2, 57, 60, 1), 9.91843934859471e-8);
  tk_true(t, "pnf.huge", clab_pnf(3.1588427192606483, 2, 57, 6e7, 0) == 1 &&
                             clab_pnf(3.1588427192606483, 2, 57, 1e308, 0) == 1 &&
                             clab_pnf(3.1588427192606483, 2, 57, HUGE_VAL, 0) == 1);
  tk_nan(t, "pnf.out", clab_pnf(1, 0, 5, 1, 0));
}
