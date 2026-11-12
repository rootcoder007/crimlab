// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests of Ripley's K kernels (src/rpk.cpp).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "cgeo.h"
#include "crimlab.h"
#include "tchk.h"

namespace {

// a message of the extended run (standalone, "ext"); silent by default
template <class... A>
void info(int big, const char *f, A... a) {
  if (big) std::printf(f, a...);
}

// Reference values from dev/doss/dst.R; the first line names the versions.
// R 4.5.0 spatstat.explore 3.8.3 spatstat.geom 3.8.3 
// Kest(X, r = seq(0, 0.25, by = 0.005), correction = ...)
const double kun[] = {
    0, 0, 0.00035794183445190159, 0.00071588366890380317,
    0.0010738255033557046, 0.0019686800894854585, 0.0026845637583892616,
    0.0030425055928411633, 0.0035794183445190158, 0.0042953020134228184,
    0.0060850111856823267, 0.0078747203579418341, 0.0096644295302013416,
    0.010917225950782998, 0.013601789709172259, 0.016107382550335572,
    0.01807606263982103, 0.020223713646532439, 0.02308724832214765,
    0.025771812080536912, 0.028277404921700222, 0.031140939597315436,
    0.03275167785234899, 0.035436241610738252, 0.038478747203579418,
    0.04116331096196868, 0.043668903803131993, 0.04742729306487696,
    0.052259507829977632, 0.058523489932885905, 0.062639821029082776,
    0.065682326621923942, 0.069261744966442954, 0.073020134228187913,
    0.077315436241610733, 0.082684563758389257, 0.088053691275167781,
    0.092885906040268459, 0.095928411633109625, 0.10183445190156599,
    0.10684563758389262, 0.11149888143176734, 0.11686800894854586,
    0.12205816554809844, 0.12742729306487696, 0.13530201342281878,
    0.14156599552572707, 0.14675615212527965, 0.15248322147651006,
    0.15767337807606263, 0.16357941834451903};
const double kbord[] = {
    0, 0, 0.00036781609195402299, 0.00074592074592074592,
    0.0010628019323671498, 0.0019607843137254902, 0.002765432098765432,
    0.0032080200501253134, 0.0038383838383838384, 0.0045801526717557254,
    0.0064615384615384613, 0.0083333333333333332, 0.010498687664041995,
    0.012063492063492064, 0.01504, 0.017344173441734417, 0.01988888888888889,
    0.022372881355932205, 0.025507246376811593, 0.028771929824561403,
    0.031904761904761908, 0.034474474474474472, 0.036386292834890963,
    0.038993710691823898, 0.042285714285714288, 0.045436893203883492,
    0.047973856209150324, 0.051881188118811879, 0.057508417508417509,
    0.063194444444444442, 0.067777777777777784, 0.07244444444444445,
    0.07644444444444444, 0.081515151515151513, 0.086588235294117646,
    0.09285714285714286, 0.098536585365853663, 0.1055, 0.10987341772151898,
    0.11765765765765766, 0.12292237442922374, 0.12857142857142856,
    0.13681159420289854, 0.14294117647058824, 0.15044776119402986, 0.15875,
    0.16550264550264551, 0.17530054644808743, 0.18311111111111111,
    0.1910344827586207, 0.19857142857142857};
const double kbm[] = {
    0, 0, 0.00036893613116048405, 0.00074926335119975213,
    0.0010462798094629351, 0.0019322096326688347, 0.0027479631949586655,
    0.0031907456411111629, 0.0038501442790908795, 0.0046336179018283104,
    0.0065936653714823972, 0.0085117848986838267, 0.010817874590543447,
    0.012540835062854746, 0.015775756019541801, 0.01821007742839563,
    0.020727092509993005, 0.023330200546934132, 0.026382740411222974,
    0.030029906307607865, 0.033308476261496391, 0.036323211402197776,
    0.037640094951006647, 0.040708037522961878, 0.044555409267525124,
    0.047861936720997124, 0.051011437264619396, 0.055693697418480216,
    0.061709126013792789, 0.067071713481432843, 0.073393115635515954,
    0.075050834854206272, 0.080837704209200473, 0.086054476242346536,
    0.090171546466525046, 0.097620500305064084, 0.10333115076117204,
    0.11032095686437517, 0.11599966857237551, 0.11899807448176955,
    0.12546607009694258, 0.12877664335094, 0.13827085268002307,
    0.14579260102549793, 0.15487867836861124, 0.15997165856442752,
    0.1682963651643741, 0.17700662919621718, 0.18657899265458433,
    0.19312405170077276, 0.19901565995525727};
// Kest(..., correction = c("none", ..., "isotropic"))$bord.modif times |W| / |W eroded by r| equals the dedicated value to 2.8e-17 
const double kiso[] = {
    0, 0, 0.00035794183445190159, 0.00071588366890380317,
    0.001131553211281532, 0.0020264077974112859, 0.002742291466315089,
    0.0031002333007669903, 0.0036371460524448427, 0.0043530297213486462,
    0.0062110052400696763, 0.0080552676751958086, 0.0098449768474553161,
    0.011097773268036971, 0.01386786778855423, 0.016492365585108238,
    0.018512430546981058, 0.020750739798642514, 0.023891400922635163,
    0.026817026814436179, 0.029389483991111132, 0.032679517015106276,
    0.03432906250004112, 0.037402088912753366, 0.040745318507344558,
    0.043699489074348306, 0.046621747365321695, 0.051053040007360016,
    0.056629685872206584, 0.064008114466838226, 0.068647431406149539,
    0.071959467462786225, 0.075828661594382701, 0.079956372104566908,
    0.084831768285303089, 0.091092088745125518, 0.097096036516082992,
    0.10262632289626245, 0.10598611210859024, 0.11303297100302644,
    0.11972290351650097, 0.12536498755278611, 0.13155500412866478,
    0.13804278113728294, 0.14477043203745404, 0.15454578521091297,
    0.16270564763165668, 0.1686923531660014, 0.17600508112189417,
    0.18287913813264559, 0.19022630034222468};

const double PI = cgeo::pi;
void patp(double *x, double *y, int n) {
  for (long long i = 0; i < n; ++i) {
    x[i] = static_cast<double>((i * i * 7919 + i * 31) % 1009) / 1009 * 2;
    y[i] = static_cast<double>((i * i * 104729 + i * 7) % 1013) / 1013;
  }
}

// the inside share of a circle by brute angular sampling: an independent
// route to the isotropic weight
double share(double px, double py, const double *w, double d, int m) {
  int in = 0;
  for (int s = 0; s < m; ++s) {
    const double a = 2 * PI * (s + 0.5) / m;
    const double qx = px + d * std::cos(a), qy = py + d * std::sin(a);
    if (qx >= w[0] && qx <= w[1] && qy >= w[2] && qy <= w[3]) ++in;
  }
  return static_cast<double>(in) / m;
}
// the same share in closed form, written separately from the kernel:
// the arcs outside each edge, less the double-counted corner arcs
double sharecf(double px, double py, const double *w, double d) {
  const double e[4] = {px - w[0], w[1] - px, py - w[2], w[3] - py};
  double out = 0;
  for (int s = 0; s < 4; ++s)
    if (e[s] < d) out += 2 * std::acos(e[s] / d);
  const int pr[4][2] = {{0, 2}, {0, 3}, {1, 2}, {1, 3}};
  for (auto &p : pr) {
    const double u = e[p[0]], v = e[p[1]];
    if (u * u + v * v < d * d)
      out -= std::acos(u / d) + std::acos(v / d) - PI / 2;
  }
  return 1 - out / (2 * PI);
}
// brute-force K by all pairs, the four estimators
void kbrute(const double *x, const double *y, int n, const double *w,
            const double *r, int m, int corr, double *k) {
  const double A = (w[1] - w[0]) * (w[3] - w[2]);
  for (int j = 0; j < m; ++j) {
    double s = 0, den = 0;
    for (int i = 0; i < n; ++i) {
      const double b = std::fmin(std::fmin(x[i] - w[0], w[1] - x[i]),
                                 std::fmin(y[i] - w[2], w[3] - y[i]));
      if ((corr == 1 || corr == 2) && b < r[j]) continue;
      den += 1;
      for (int t = 0; t < n; ++t) {
        if (t == i) continue;
        const double d = std::hypot(x[t] - x[i], y[t] - y[i]);
        if (d > r[j]) continue;
        s += corr == 3 && d > 0 ? 1 / sharecf(x[i], y[i], w, d) : 1;
      }
    }
    if (corr == 1)
      k[j] = s / den * A / n;
    else if (corr == 2)
      k[j] = s * A * A /
             (double(n) * (n - 1) * (w[1] - w[0] - 2 * r[j]) *
              (w[3] - w[2] - 2 * r[j]));
    else
      k[j] = s * A / (double(n) * (n - 1));
  }
}

void famtest(tk &T, int big) {
  const double win[4] = {0, 2, 0, 1};
  static int iw[200000];
  static double wk[200000], k[64], kb[64], out[400];
  double r[51];
  for (int j = 0; j < 51; ++j) r[j] = 0.005 * j;
  double x[150], y[150];
  patp(x, y, 150);
  // spatstat Kest parity, all four estimators
  const double *ref[4] = {kun, kbord, kbm, kiso};
  const char *id[4] = {"k.none", "k.border", "k.bordmodif", "k.iso"};
  for (int c = 0; c < 4; ++c) {
    tk_true(T, "k.ok", clab_ripk(x, y, 150, win, r, 51, c, iw, wk, k) == CLAB_OK);
    for (int j = 0; j < 51; ++j) tk_eq(T, id[c], k[j], ref[c][j], 1e-12);
  }
  // the isotropic weight: closed form against angular sampling
  {
    cgeo::rng g(3);
    // a midpoint rule on an indicator with at most 8 jumps: error below
    // 8 / m
    const int m = big ? 200000 : 20000;
    for (int s = 0; s < (big ? 40 : 10); ++s) {
      const double px = 2 * g.unif(), py = g.unif(), d = 0.6 * g.unif() + 0.01;
      tk_eq(T, "k.share", sharecf(px, py, win, d), share(px, py, win, d, m),
            8.0 / m);
    }
    // a corner point sees a quarter circle, an edge point half
    tk_eq(T, "k.corner", sharecf(0, 0, win, 0.3), 0.25, 1e-15);
    tk_eq(T, "k.edge", sharecf(1, 0, win, 0.3), 0.5, 1e-15);
  }
  // brute force on random patterns, irregular radii
  {
    cgeo::rng g(8);
    for (int rep = 0; rep < 8; ++rep) {
      const int n = 30 + 25 * rep;
      double u[250], v[250], rr[12];
      for (int i = 0; i < n; ++i) {
        u[i] = 2 * g.unif();
        v[i] = g.unif();
      }
      if (rep == 3) {  // coincident points
        u[1] = u[0];
        v[1] = v[0];
      }
      for (int j = 0; j < 12; ++j) rr[j] = 0.02 * j + 0.013 * g.unif();
      for (int c = 0; c < 4; ++c) {
        clab_ripk(u, v, n, win, rr, 12, c, iw, wk, k);
        kbrute(u, v, n, win, rr, 12, c, kb);
        for (int j = 0; j < 12; ++j) tk_eq(T, "k.brute", k[j], kb[j], 1e-12);
      }
    }
  }
  // under CSR the isotropic estimator averages pi r^2 (Ripley 1977)
  {
    cgeo::rng g(13);
    const int reps = big ? 300 : 60, n = 100;
    double rr[4] = {0.05, 0.1, 0.15, 0.2}, s1[4] = {0}, s2[4] = {0};
    for (int rep = 0; rep < reps; ++rep) {
      double u[100], v[100];
      for (int i = 0; i < n; ++i) {
        u[i] = 2 * g.unif();
        v[i] = g.unif();
      }
      clab_ripk(u, v, n, win, rr, 4, 3, iw, wk, k);
      for (int j = 0; j < 4; ++j) {
        s1[j] += k[j];
        s2[j] += k[j] * k[j];
      }
    }
    for (int j = 0; j < 4; ++j) {
      const double m = s1[j] / reps, sd = std::sqrt((s2[j] / reps - m * m) / reps);
      tk_true(T, "k.csr", std::fabs(m - PI * rr[j] * rr[j]) < 4 * sd);
    }
  }
  // invariance: translation, transposition of the window, permutation
  {
    double u[150], v[150], k2[51];
    for (int i = 0; i < 150; ++i) {
      u[149 - i] = y[i] + 7e4;
      v[149 - i] = x[i] - 2e4;
    }
    const double w2[4] = {7e4, 1 + 7e4, -2e4, 2 - 2e4};
    for (int c = 0; c < 4; ++c) {
      clab_ripk(x, y, 150, win, r, 51, c, iw, wk, k);
      clab_ripk(u, v, 150, w2, r, 51, c, iw, wk, k2);
      for (int j = 1; j < 51; ++j) tk_eq(T, "k.inv", k2[j], k[j], 1e-9);
    }
  }
  // envelopes and the global test: nominal size under CSR (nsim = 19:
  // p <= 0.05 exactly when the data beat every simulation), power against
  // clustering
  {
    cgeo::rng g(21);
    double rr[10];
    for (int j = 0; j < 10; ++j) rr[j] = 0.02 * (j + 1);
    int rej = 0;
    const int reps = big ? 400 : 60;
    for (int rep = 0; rep < reps; ++rep) {
      double u[60], v[60];
      for (int i = 0; i < 60; ++i) {
        u[i] = 2 * g.unif();
        v[i] = g.unif();
      }
      clab_ripenv(u, v, 60, win, rr, 10, rep % 2 ? 3 : 1, 19, 77 + rep, iw,
                  wk, out);
      if (out[41] <= 0.05) ++rej;
    }
    const double sd = std::sqrt(0.05 * 0.95 / reps);
    info(big, "  envelope test size under CSR: %.4f\n", rej / double(reps));
    tk_true(T, "env.size", std::fabs(rej / double(reps) - 0.05) < 4 * sd);
    double u[60], v[60];
    for (int i = 0; i < 60; ++i) {
      const int c = i % 3;
      u[i] = 0.4 + 0.6 * c + 0.06 * g.unif();
      v[i] = 0.5 + 0.06 * g.unif();
    }
    clab_ripenv(u, v, 60, win, rr, 10, 3, 99, 1, iw, wk, out);
    tk_true(T, "env.power", out[41] == 0.01 && out[0] > out[20]);
    // the envelope brackets every simulated curve: lo <= mean <= hi
    bool br = true;
    for (int j = 0; j < 10; ++j) br = br && out[10 + j] <= out[30 + j] && out[30 + j] <= out[20 + j];
    tk_true(T, "env.bracket", br);
  }
  // beyond half the short side the modified border estimator has no
  // eroded window, and the border one no centre: both NaN
  {
    const double rr[3] = {0.2, 0.5, 0.6};
    clab_ripk(x, y, 150, win, rr, 3, 2, iw, wk, k);
    tk_true(T, "k.bmnan", !std::isnan(k[0]) && std::isnan(k[1]) && std::isnan(k[2]));
    clab_ripk(x, y, 150, win, rr, 3, 1, iw, wk, k);
    tk_true(T, "k.bnan", !std::isnan(k[0]) && std::isnan(k[2]));
  }
  // guards
  const double wrev[4] = {2, 0, 0, 1}, whuge[4] = {-1e151, 2, 0, 1},
               wsmall[4] = {0, 1, 0, 1};
  tk_true(T, "k.wrev", clab_ripk(x, y, 150, wrev, r, 51, 0, iw, wk, k) == CLAB_EDOM);
  tk_true(T, "k.whuge", clab_ripk(x, y, 150, whuge, r, 51, 0, iw, wk, k) == CLAB_EDOM);
  tk_true(T, "k.wout", clab_ripk(x, y, 150, wsmall, r, 51, 0, iw, wk, k) == CLAB_EDOM);
  double rbad[2] = {0.1, 0.1};
  tk_true(T, "k.radii", clab_ripk(x, y, 150, win, rbad, 2, 0, iw, wk, k) == CLAB_EDOM);
  tk_true(T, "k.corr", clab_ripk(x, y, 150, win, r, 51, 4, iw, wk, k) == CLAB_EARG);
  tk_true(T, "k.n", clab_ripk(x, y, 1, win, r, 51, 0, iw, wk, k) == CLAB_EARG);
  tk_true(T, "env.nsim", clab_ripenv(x, y, 150, win, r, 51, 0, 0, 1, iw, wk, out) ==
                     CLAB_EARG);
  double rneg[1] = {-0.1};
  tk_true(T, "k.neg", clab_ripk(x, y, 150, win, rneg, 1, 0, iw, wk, k) == CLAB_EDOM);
}

#ifdef CLAB_MAIN
double now() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

// a heap buffer for the timings (kept out of the program's static size)
template <class T>
struct hb {
  T *p;
  explicit hb(long n) : p(static_cast<T *>(std::malloc(sizeof(T) * n))) {
    if (p == nullptr) std::abort();
  }
  ~hb() { std::free(p); }
  hb(const hb &) = delete;
  hb &operator=(const hb &) = delete;
};

double med5(double *v) { return tk_q7(v, 5, 0.5); }

void fambench() {
  hb<double> x_(1000000), y_(1000000), wk_(3000000);
  double *x = x_.p, *y = y_.p, *wk = wk_.p;

  double k[20], eo[64];
  hb<int> iw_(2000100);
  int *iw = iw_.p;
  cgeo::rng g(4);
  for (int i = 0; i < 1000000; ++i) {
    x[i] = 1e4 * g.unif();
    y[i] = 1e4 * g.unif();
  }
  const double win[4] = {0, 1e4, 0, 1e4};
  const int ns[3] = {1000, 100000, 1000000};
  for (int n : ns) {
    // radii up to the distance holding about 20 neighbours on average
    const double rmax = std::sqrt(20 * 1e8 / (PI * n));
    double r[20];
    for (int j = 0; j < 20; ++j) r[j] = rmax * (j + 1) / 20;
    double t[2][5];
    for (int s = 0; s < 5; ++s) {
      double a = now();
      clab_ripk(x, y, n, win, r, 20, 1, iw, wk, k);
      t[0][s] = now() - a;
      a = now();
      clab_ripk(x, y, n, win, r, 20, 3, iw, wk, k);
      t[1][s] = now() - a;
    }
    std::printf("bench rpk n = %d (about 20 neighbours within rmax): border "
                "%.5f isotropic %.5f s (medians of 5)\n",
                n, med5(t[0]), med5(t[1]));
  }
  double r[10];
  for (int j = 0; j < 10; ++j) r[j] = 50.0 * (j + 1);
  const double a = now();
  clab_ripenv(x, y, 1000, win, r, 10, 3, 999, 1, iw, wk, eo);
  std::printf("bench rpk n = 1000, 999 simulations: %.4f s\n", now() - a);
}

#endif

}  // namespace

void t_rpk(tk &T, int big) { famtest(T, big); }
#ifdef CLAB_MAIN
void b_rpk() { fambench(); }
#endif
