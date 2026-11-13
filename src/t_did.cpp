// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the difference-in-differences kernels (did.cpp,
// stg.cpp): closed forms of the 2 x 2 design and the balanced within
// transformation, the identities of the decomposition and of two-way
// clustering, and the failure paths. Part of the self-test (selftest.cpp).
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "pev.h"
#include "pevs.h"
#include "tchk.h"

namespace {

void run(tk &tt) {
  // 1. the balanced within transformation is y - ybar_i - ybar_t + ybar
  {
    const int nu = 3, nt = 4, n = 12;
    double v[n], y[n], wk[2 * (nu + nt)];
    int u[n], t[n];
    for (int i = 0; i < n; ++i) {
      u[i] = i % nu;
      t[i] = i / nu;
      y[i] = v[i] = std::sin(1.0 + i) + 0.1 * i * i;
    }
    tk_ieq(tt, "did.within.rc", pev::within(v, n, 1, u, nu, t, nt, 1, nullptr, 1e-12, 100, wk), 1);
    double mu[nu] = {0, 0, 0}, mt[nt] = {0, 0, 0, 0}, m = 0;
    for (int i = 0; i < n; ++i) {
      mu[u[i]] += y[i] / nt;
      mt[t[i]] += y[i] / nu;
      m += y[i] / n;
    }
    for (int i = 0; i < n; ++i) tk_near(tt, "did.within", v[i], y[i] - mu[u[i]] - mt[t[i]] + m, 1e-13);
    // no sweep allowed: no convergence
    tk_ieq(tt, "did.within.max", pev::within(v, n, 1, u, nu, t, nt, 1, nullptr, 1e-12, 0, wk), -1);
  }
  // 2. a singular matrix has no Cholesky factor; collinear columns no fit
  {
    double a[4] = {1, 1, 1, 1};
    tk_true(tt, "did.chol", !pev::chol(a, 2));
    const double x[6] = {1, 2, 3, 2, 4, 6}, y[3] = {1, 2, 4};
    double b[2], xi[4], r[3], wk[6];
    tk_ieq(tt, "did.ols", pev::ols(x, y, 3, 2, nullptr, b, xi, r, wk), CLAB_EDOM);
  }
  // 3. the 2 x 2 design: the coefficient is (y11 - y10) - (y01 - y00) of
  // the cell means; out-of-range sizes and no convergence are refused
  {
    const int n = 8;
    const double y[n] = {1, 2, 3, 5, 2, 3, 6, 9};  // units 0..3, periods 0, 1
    const int u[n] = {0, 1, 2, 3, 0, 1, 2, 3}, t[n] = {0, 0, 0, 0, 1, 1, 1, 1};
    const double x[n] = {0, 0, 0, 0, 0, 0, 1, 1};
    double b[1], vc[1], st[5], wk[n * 3 + 1 + 2 * 6 + (n + 2) * 1 + 1];
    int iw[n + 16];
    tk_ieq(tt, "did.fe.rc", clab_didfe(y, x, n, 1, u, 4, t, 2, 1, nullptr, u, 4, nullptr, 0, 1e-12,
                                       100, b, vc, st, wk, iw),
           CLAB_OK);
    tk_near(tt, "did.fe.b", b[0], (7.5 - 4) - (2.5 - 1.5), 1e-13);
    tk_ieq(tt, "did.fe.k", static_cast<int>(st[2]), 3);
    // two-way clusters with the same clusters twice: V + V - V = V
    double v2[1], st2[5];
    tk_ieq(tt, "did.fe.rc2", clab_didfe(y, x, n, 1, u, 4, t, 2, 1, nullptr, u, 4, u, 4, 1e-12, 100,
                                        b, v2, st2, wk, iw),
           CLAB_OK);
    tk_near(tt, "did.fe.cl2", v2[0], vc[0], 1e-13);
    tk_ieq(tt, "did.fe.earg", clab_didfe(y, x, 1, 1, u, 4, t, 2, 1, nullptr, u, 4, nullptr, 0,
                                         1e-12, 100, b, vc, st, wk, iw),
           CLAB_EARG);
    // unbalanced (one row less), no sweep allowed
    tk_ieq(tt, "did.fe.max", clab_didfe(y, x, n - 1, 1, u, 4, t, 2, 1, nullptr, u, 4, nullptr, 0,
                                        1e-12, 0, b, vc, st, wk, iw),
           CLAB_EDOM);
    // a treatment the effects absorb
    const double xa[n] = {1, 1, 1, 1, 1, 1, 1, 1};
    tk_ieq(tt, "did.fe.abs", clab_didfe(y, xa, n, 1, u, 4, t, 2, 1, nullptr, u, 4, nullptr, 0,
                                        1e-12, 100, b, vc, st, wk, iw),
           CLAB_EDOM);
  }
  // 4. the Goodman-Bacon decomposition sums to the two-way estimate, its
  // weights to 1; with never-treated, early, late and always-treated units
  {
    const int nu = 6, nt = 5, n = nu * nt;
    const int g[nu] = {1, 1, 3, 0, 5, 5};
    double yw[n], y[n], x[n];
    int u[n], t[n];
    for (int s = 0; s < nt; ++s) {
      for (int i = 0; i < nu; ++i) {
        const int q = s * nu + i;
        const double on = s >= g[i] ? 1.0 : 0.0;
        yw[q] = y[q] = i + 0.3 * s + on * (1.0 + 0.5 * (s - g[i])) + std::cos(q);
        u[q] = i;
        t[q] = s;
        x[q] = on;
      }
    }
    int ty[30], tr[30], co[30], iw[2 * nu];
    double es[30], wt[30], wk[64];
    const int r = pev::bacon(yw, nu, nt, g, ty, tr, co, es, wt, wk, iw);
    double sw = 0, sb = 0;
    for (int q = 0; q < r; ++q) {
      sw += wt[q];
      sb += wt[q] * es[q];
    }
    double b[1], vc[1], st[5], fw[n * 3 + 1 + 22 + (n + 2) + 1];
    int fi[n + nu];
    clab_didfe(y, x, n, 1, u, nu, t, nt, 1, nullptr, u, nu, nullptr, 0, 1e-13, 100, b, vc, st, fw,
               fi);
    tk_near(tt, "did.bacon.w", sw, 1.0, 1e-14);
    tk_near(tt, "did.bacon.sum", sb, b[0], 1e-12);
    // three timing groups (0, 1, 3) and never: 2 against never, 1 earlier
    // against later (the always treated never is), 3 later against earlier
    tk_ieq(tt, "did.bacon.rows", r, 6);
  }
  // 5. matching: the nearest comparison on one feature, the caliper, a
  // singular covariance
  {
    const double f[5] = {0.0, 1.0, 0.2, 2.0, 5.0};
    const int trt[5] = {1, 1, 0, 0, 0}, ex[5] = {0, 0, 0, 0, 0};
    int mt[5], mc[5];
    double md[5], wu[5], wk[16];
    tk_ieq(tt, "did.mtch.n", pev::mtch(f, 5, 1, trt, ex, 1, 0.0, mt, mc, md, wu, wk), 2);
    tk_ieq(tt, "did.mtch.c0", mc[0], 2);
    tk_ieq(tt, "did.mtch.c1", mc[1], 2);
    tk_near(tt, "did.mtch.w", wu[2], 2.0, 1e-15);
    tk_near(tt, "did.mtch.w3", wu[3], 0.0, 1e-15);
    const double f2[10] = {0, 1, 2, 3, 4, 0, 2, 4, 6, 8};
    tk_ieq(tt, "did.mtch.sing", pev::mtch(f2, 5, 2, trt, ex, 1, 0.0, mt, mc, md, wu, wk), -1);
  }
  // 6. Callaway-Sant'Anna on the 2 x 2 design is the difference in
  // differences; DID_M with one switch too; Wooldridge's coefficient too
  {
    const double yw[8] = {1, 2, 3, 5, 2, 3, 6, 9};
    const int g[4] = {2, 2, 1, 1}, dw[8] = {0, 0, 0, 0, 0, 0, 1, 1};
    pev::stg s{};
    s.yw = yw;
    s.dw = dw;
    s.g = g;
    s.nu = 4;
    s.nt = 2;
    pev::cells(g, 4, 2, s.cg, s.ct, 0);
    std::vector<int> idx = {0, 1, 2, 3}, ev = {0};
    double v[3];
    pev::stat(s, idx, ev, v);
    tk_near(tt, "did.csa", v[0], 2.5, 1e-14);
    tk_near(tt, "did.csa.all", v[1], 2.5, 1e-14);
    s.meth = 1;
    pev::stat(s, idx, ev, v);
    tk_near(tt, "did.wool", v[0], 2.5, 1e-13);
    s.meth = 2;
    double d4[4];
    pev::stat(s, idx, ev, d4);
    tk_near(tt, "did.dcdh", d4[0], 2.5, 1e-14);
    tk_nan(tt, "did.dcdh.pl", d4[1]);
    // no unit outside the cohort: no comparison, no cells for Wooldridge
    std::vector<int> two = {2, 3};
    s.meth = 0;
    pev::stat(s, two, ev, v);
    tk_nan(tt, "did.csa.none", v[0]);
    s.meth = 1;
    pev::stat(s, two, ev, v);
    tk_nan(tt, "did.wool.none", v[0]);
  }
}

}  // namespace

void t_did(tk &t, int) { run(t); }
