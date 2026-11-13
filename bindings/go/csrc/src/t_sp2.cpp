// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/sp2.cpp and of the Manhattan options of clab_ripkm and
// clab_mantm: brute force, closed forms, hand-computed cases and seeded
// simulations with known truth.
#include <cmath>
#include <vector>

#include "cgeo.h"
#include "crimlab.h"
#include "tchk.h"

namespace {

using vd = std::vector<double>;
using vi = std::vector<int>;
const double PI = cgeo::pi;

double l1(const vd &x, const vd &y, int i, int j) {
  return std::fabs(x[i] - x[j]) + std::fabs(y[i] - y[j]);
}

void famtest(tk &T, int big) {
  rng g = {77};
  // nnidm: means against brute force; the expectations in closed form
  {
    const int n = 150, km = 3;
    vd x(n), y(n), out(6 * km + 4), dist(n * km), work(n * km + 2 * n + 2 * km);
    vi iw(2 * n + 2 + n * km + 2 * km);
    for (int i = 0; i < n; ++i) {
      x[i] = runi(g, 0, 2);
      y[i] = runi(g, 0, 1);
    }
    const double win[4] = {0, 2, 0, 1};
    tk_true(T, "nnidm.ok", clab_nnidm(x.data(), y.data(), n, win, km, 19, 3, dist.data(),
                                      iw.data(), work.data(), out.data()) == CLAB_OK);
    for (int k = 0; k < km; ++k) {
      double s = 0;
      for (int i = 0; i < n; ++i) {
        vd d;
        for (int j = 0; j < n; ++j)
          if (j != i) d.push_back(l1(x, y, i, j));
        for (std::size_t a = 0; a <= static_cast<std::size_t>(k); ++a)
          for (std::size_t b = a + 1; b < d.size(); ++b)
            if (d[b] < d[a]) std::swap(d[a], d[b]);
        s += d[static_cast<std::size_t>(k)];
      }
      tk_eq(T, "nnidm.mean", out[6 * k], s / n, 1e-12);
    }
    const double lam = n / 2.0;
    tk_eq(T, "nnidm.e1", out[1], 0.5 * std::sqrt(PI / (2 * lam)), 1e-14);
    tk_eq(T, "nnidm.e2", out[7], 1.5 * out[1], 1e-14);
    tk_eq(T, "nnidm.sd1", out[5], std::sqrt(1 / (2 * lam) - out[1] * out[1]), 1e-14);
    tk_true(T, "nnidm.donnelly", std::isnan(out[6 * km]) && std::isnan(out[6 * km + 3]));
    tk_true(T, "nnidm.p", out[3] > 0 && out[3] <= 1 && out[4] > 0 && out[4] <= 1);
  }
  // the expectation by simulation, far from the edge: Poisson points in a
  // square of side 40 at intensity 1, distances from the central square of
  // side 20 only (the nearest neighbour of a central point lies inside)
  {
    const int sets = big ? 40 : 8;
    double s = 0, s2 = 0;
    int m = 0;
    for (int t = 0; t < sets; ++t) {
      const int n = 1600;
      vd x(n), y(n);
      for (int i = 0; i < n; ++i) {
        x[i] = runi(g, 0, 40);
        y[i] = runi(g, 0, 40);
      }
      for (int i = 0; i < n; ++i) {
        if (x[i] < 10 || x[i] > 30 || y[i] < 10 || y[i] > 30) continue;
        double b = 1e300;
        for (int j = 0; j < n; ++j)
          if (j != i) b = std::fmin(b, l1(x, y, i, j));
        s += b;
        s2 += b * b;
        ++m;
      }
    }
    // E R_1 = 0.5 sqrt(pi / 2) at lam = 1, E R_1^2 = 1 / 2
    const double se = std::sqrt((4 - PI) / 8 / m);
    tk_lt(T, "nnidm.csr", std::fabs(s / m - 0.5 * std::sqrt(PI / 2)), 4 * se);
    tk_lt(T, "nnidm.csr2", std::fabs(s2 / m - 0.5), 0.02);
  }
  // ripkm: brute force for none, border and modified border; E K = 2 r^2
  {
    const int n = 80, m = 6;
    vd x(n), y(n), r(m), k(m), wk(4 * m);
    vi iw(2 * n + 2);
    for (int i = 0; i < n; ++i) {
      x[i] = runi(g, 0, 1);
      y[i] = runi(g, 0, 1);
    }
    for (int j = 0; j < m; ++j) r[j] = 0.03 * (j + 1);
    const double win[4] = {0, 1, 0, 1};
    for (int corr = 0; corr <= 2; ++corr) {
      clab_ripkm(x.data(), y.data(), n, win, r.data(), m, corr, 1, iw.data(), wk.data(), k.data());
      bool ok = true;
      for (int j = 0; j < m; ++j) {
        double s = 0, cen = 0;
        for (int i = 0; i < n; ++i) {
          const double b = std::fmin(std::fmin(x[i], 1 - x[i]), std::fmin(y[i], 1 - y[i]));
          const bool c = corr == 0 || b >= r[j];
          if (corr != 0 && c) cen += 1;
          for (int q = 0; q < n; ++q)
            if (q != i && c && l1(x, y, i, q) <= r[j]) s += 1;
        }
        const double e = 1 - 2 * r[j];
        const double want = corr == 0 ? s / (n * (n - 1.0))
                            : corr == 1 ? s / cen / n
                                        : s / (n * (n - 1.0) * e * e);
        ok = ok && std::fabs(k[j] - want) <= 1e-12 * (1 + want);
      }
      tk_true(T, "ripkm.brute", ok);
    }
    tk_ieq(T, "ripkm.iso", clab_ripkm(x.data(), y.data(), n, win, r.data(), m, 3, 1,
                                     iw.data(), wk.data(), k.data()),
           CLAB_EDOM);
    // metric 0 is clab_ripk
    vd k0(m);
    clab_ripk(x.data(), y.data(), n, win, r.data(), m, 3, iw.data(), wk.data(), k0.data());
    clab_ripkm(x.data(), y.data(), n, win, r.data(), m, 3, 0, iw.data(), wk.data(), k.data());
    tk_veq(T, "ripkm.metric0", k.data(), k0.data(), m, 0);
    // under CSR the border estimator is unbiased for 2 r^2
    const int sets = big ? 300 : 60;
    double s = 0, s2 = 0;
    const double rr = 0.1;
    for (int t = 0; t < sets; ++t) {
      for (int i = 0; i < n; ++i) {
        x[i] = runi(g, 0, 1);
        y[i] = runi(g, 0, 1);
      }
      double kk;
      clab_ripkm(x.data(), y.data(), n, win, &rr, 1, 1, 1, iw.data(), wk.data(), &kk);
      s += kk;
      s2 += kk * kk;
    }
    const double mu = s / sets, sd = std::sqrt((s2 / sets - mu * mu) / sets);
    tk_lt(T, "ripkm.csr", std::fabs(mu - 2 * rr * rr), 4 * sd);
  }
  // mantm: the Manhattan dissimilarity written out; metric 0 is mantst
  {
    const int n = 25;
    vd x(n), y(n), t(n);
    vi iw(n);
    for (int i = 0; i < n; ++i) {
      x[i] = runi(g);
      y[i] = runi(g);
      t[i] = x[i] + runi(g);
    }
    double a[5], b[5];
    clab_mantm(x.data(), y.data(), t.data(), n, 1, 0, 0, 0, 1, iw.data(), a);
    double sx = 0, st = 0, sxx = 0, stt = 0, sxt = 0, np = 0;
    for (int i = 1; i < n; ++i)
      for (int j = 0; j < i; ++j) {
        const double u = l1(x, y, i, j), v = std::fabs(t[i] - t[j]);
        sx += u;
        st += v;
        sxx += u * u;
        stt += v * v;
        sxt += u * v;
        np += 1;
      }
    const double cv = sxt / np - sx / np * st / np;
    const double rr = cv / std::sqrt((sxx / np - sx * sx / np / np) * (stt / np - st * st / np / np));
    tk_eq(T, "mantm.r", a[0], rr, 1e-12);
    clab_mantm(x.data(), y.data(), t.data(), n, 0, 0.1, 0, 9, 4, iw.data(), a);
    clab_mantst(x.data(), y.data(), t.data(), n, 0.1, 0, 9, 4, iw.data(), b);
    tk_veq(T, "mantm.metric0", a, b, 5, 0);
    tk_ieq(T, "mantm.metric", clab_mantm(x.data(), y.data(), t.data(), n, 2, 0, 0, 0, 1,
                                         iw.data(), a),
           CLAB_EARG);
  }
  // head-bang on a 3 x 3 lattice with a spike in the middle: the middle
  // falls to the median screens of its opposite pairs; corners, with no
  // triple of angle 135 degrees or more, keep their values
  {
    const int n = 9, nn = 8, tt = nn * (nn - 1) / 2;
    vd x(n), y(n), v(n, 1.0), out(n), w(n, 1.0);
    for (int i = 0; i < n; ++i) {
      x[i] = i % 3;
      y[i] = i / 3;
    }
    v[4] = 100;
    v[0] = 7;
    vi iw(2 * n + 2 + n * nn + 3 * tt);
    vd wk(n * nn + n + 5 * tt);
    tk_true(T, "hdbang.ok", clab_hdbang(x.data(), y.data(), v.data(), nullptr, n, nn, 135, 0, 1,
                                        iw.data(), wk.data(), out.data()) == CLAB_OK);
    // middle: triples (0,8), (1,7), (2,6), (3,5) at 180 degrees and the
    // eight at 135 (e.g. (0,7)); lows all 1, highs 1 except the pairs with
    // the corner 0: (0,8), (0,5), (0,7) give 7; median of the 12 highs = 1
    tk_eq(T, "hdbang.mid", out[4], 1, 0);
    tk_eq(T, "hdbang.corner", out[0], 7, 0);
    // an edge midpoint (1, 0): its triples are (0, 2) at 180 degrees and
    // (0, 5), (2, 3) at 135: lows 1, highs 7, 7, 1, so its value 1 stays
    tk_eq(T, "hdbang.edge", out[1], 1, 0);
    // weighted: equal weights give the same screens here (odd counts)
    vd ow(n);
    clab_hdbang(x.data(), y.data(), v.data(), w.data(), n, nn, 135, 0, 1, iw.data(), wk.data(),
                ow.data());
    tk_true(T, "hdbang.wsame", ow[4] == 1 && ow[0] == 7);
    // a constant surface is a fixed point; zero passes change nothing
    vd c(n, 3.5);
    clab_hdbang(x.data(), y.data(), c.data(), nullptr, n, nn, 135, 0, 4, iw.data(), wk.data(),
                out.data());
    tk_true(T, "hdbang.const", out[4] == 3.5 && out[0] == 3.5);
    clab_hdbang(x.data(), y.data(), v.data(), nullptr, n, nn, 135, 2, 0, iw.data(), wk.data(),
                out.data());
    tk_true(T, "hdbang.iter0", out[4] == 100);
    tk_ieq(T, "hdbang.ang", clab_hdbang(x.data(), y.data(), v.data(), nullptr, n, nn, 60, 0, 1,
                                       iw.data(), wk.data(), out.data()),
           CLAB_EDOM);
  }
  // correlated random walk: a square walk, a straight line, the expected
  // squared displacement against simulation
  {
    const double sx[5] = {0, 1, 1, 0, 0}, sy[5] = {0, 0, 1, 1, 0};
    double st[4], tu[3], o[5];
    clab_crwfit(sx, sy, 5, st, tu, o);
    tk_true(T, "crw.square", o[0] == 1 && o[1] == 1 && std::fabs(o[2]) < 1e-15 &&
                                 std::fabs(o[3] - 1) < 1e-15);
    tk_eq(T, "crw.head", o[4], -PI / 2, 1e-15);
    double ob[4], ms[4];
    clab_crwmsd(sx, sy, 5, 1, 1, 1, 4, ob, ms);
    tk_true(T, "crw.line", ms[0] == 1 && ms[2] == 9);
    tk_true(T, "crw.obs", ob[0] == 1 && ob[1] == 2 && ob[3] == 0);
    const double one = 1, zero = 0;
    double px[3], py[3];
    clab_crwsim(&one, 1, &zero, 1, 2, 3, 0, 3, 1, 1, px, py);
    tk_true(T, "crw.straight", px[2] == 5 && py[2] == 3);
    // steps 1 or 2, turns +-0.6: c = cos 0.6, s = 0; simulated mean R_k^2
    // within 4 standard errors of the expectation
    const double steps[2] = {1, 2}, turns[2] = {0.6, -0.6};
    const int ns = big ? 40000 : 8000, h = 6;
    vd qx(static_cast<std::size_t>(ns) * h), qy(static_cast<std::size_t>(ns) * h);
    clab_crwsim(steps, 2, turns, 2, 0, 0, 0.3, h, ns, 11, qx.data(), qy.data());
    double e[6], dummy[6];
    clab_crwmsd(sx, sy, 0, 1.5, 2.5, std::cos(0.6), h, dummy, e);
    for (int k = 0; k < h; ++k) {
      double s = 0, s2 = 0;
      for (int i = 0; i < ns; ++i) {
        const double d2 = qx[k * ns + i] * qx[k * ns + i] + qy[k * ns + i] * qy[k * ns + i];
        s += d2;
        s2 += d2 * d2;
      }
      const double mu = s / ns, se = std::sqrt((s2 / ns - mu * mu) / ns);
      tk_lt(T, "crw.ks83", std::fabs(mu - e[k]), 4 * se);
    }
    const double rep[3] = {0, 0, 1};
    tk_ieq(T, "crw.repeat", clab_crwfit(rep, rep, 3, st, tu, o), CLAB_EDOM);
  }
  // risk-adjusted NNH, by hand: one baseline unit at the origin (size 1,
  // bandwidth 1), points (0, 0), (0.3, 0), (5, 0), z = 0
  {
    const double x[3] = {0, 0.3, 5}, y[3] = {0, 0, 0}, b0 = 0, b1 = 1;
    double thr[3], cs[24], wk[9];
    int cl[3], ncl = -1, iw[12];
    tk_true(T, "rnnh.ok", clab_rnnhcl(x, y, 3, &b0, &b0, &b1, 1, 1, 0, 2, 0, 1, iw, wk, thr, cl,
                                      &ncl, cs) == CLAB_OK);
    double lam[3];
    for (int i = 0; i < 3; ++i) lam[i] = 3 * std::exp(-x[i] * x[i] / 2) / (2 * PI);
    for (int i = 0; i < 3; ++i) tk_eq(T, "rnnh.thr", thr[i], 0.5 / std::sqrt(lam[i]), 1e-14);
    // pair (0, 1): threshold 0.5 / sqrt(mean lam) = 0.732 > 0.3, linked;
    // the far point is not
    tk_true(T, "rnnh.cl", ncl == 1 && cl[0] == 0 && cl[1] == 0 && cl[2] == -1);
    tk_eq(T, "rnnh.mean", cs[1], 0.15, 1e-15);
    // a uniform baseline reproduces the random-expectation threshold:
    // a large grid of equal sizes about the points
    vd bx, by, bw;
    for (int a = -40; a <= 40; ++a)
      for (int c = -40; c <= 40; ++c) {
        bx.push_back(a * 0.25);
        by.push_back(c * 0.25);
        bw.push_back(1);
      }
    clab_rnnhcl(x, y, 3, bx.data(), by.data(), bw.data(), static_cast<int>(bx.size()), 0.5, 0,
                2, 0, 1, iw, wk, thr, cl, &ncl, cs);
    // lam = 3 / A with A = 81 * 81 * 0.0625 square units: 0.5 sqrt(A / 3)
    tk_eq(T, "rnnh.uniform", thr[0], 0.5 * std::sqrt(81 * 81 * 0.0625 / 3), 1e-6);
    // a two-point cluster has no ellipse with divisor offset 2
    clab_rnnhcl(x, y, 3, &b0, &b0, &b1, 1, 1, 0, 2, 2, 1, iw, wk, thr, cl, &ncl, cs);
    tk_true(T, "rnnh.noell", ncl == 1 && std::isnan(cs[5]));
    tk_ieq(T, "rnnh.bwid", clab_rnnhcl(x, y, 3, &b0, &b0, &b1, 1, 0, 0, 2, 0, 1, iw, wk, thr,
                                       cl, &ncl, cs),
           CLAB_EDOM);
    tk_ieq(T, "rnnh.minpts", clab_rnnhcl(x, y, 3, &b0, &b0, &b1, 1, 1, 0, 1, 0, 1, iw, wk, thr,
                                         cl, &ncl, cs),
           CLAB_EARG);
  }
}

}  // namespace

void t_sp2(tk &T, int big) { famtest(T, big); }
