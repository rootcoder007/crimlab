// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The self-test of libcrimlab: every kernel checked through the C interface
// against hand formulas, brute-force enumeration and seeded simulation. No R
// here. Built with -DCLAB_MAIN it is a program of its own:
//   g++ -std=c++17 -DCLAB_MAIN src/core.cpp src/t_*.cpp src/selftest.cpp
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>

#include "tchk.h"

bool tk_eq(tk &t, const char *id, double got, double want, double tol) {
  return tk_veq(t, id, &got, &want, 1, tol);
}
bool tk_near(tk &t, const char *id, double got, double want, double tol) {
  const double a = std::fabs(want);
  const bool ok = (std::isnan(got) && std::isnan(want)) || got == want ||
                  std::fabs(got - want) <= tol * (a > 1.0 ? a : 1.0);
  return tk_rec(t, id, ok, want, got);
}
bool tk_true(tk &t, const char *id, bool c) {
  return tk_rec(t, id, c, 1.0, c ? 1.0 : 0.0);
}
bool tk_lt(tk &t, const char *id, double a, double b) {
  return tk_rec(t, id, a < b, b, a);
}
bool tk_le(tk &t, const char *id, double a, double b) {
  return tk_rec(t, id, a <= b, b, a);
}
bool tk_gt(tk &t, const char *id, double a, double b) {
  return tk_rec(t, id, a > b, b, a);
}
bool tk_ge(tk &t, const char *id, double a, double b) {
  return tk_rec(t, id, a >= b, b, a);
}
bool tk_nan(tk &t, const char *id, double x) {
  return tk_rec(t, id, std::isnan(x), NAN, x);
}
bool tk_ieq(tk &t, const char *id, int got, int want) {
  return tk_rec(t, id, got == want, want, got);
}
std::uint64_t rnext(rng &r) {
  std::uint64_t z = (r.s += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}
double runi(rng &r) {
  return (static_cast<double>(rnext(r) >> 11) + 0.5) * 0x1.0p-53;
}
double runi(rng &r, double a, double b) { return a + (b - a) * runi(r); }
int rsam(rng &r, int a, int b) {
  const std::uint64_t k = static_cast<std::uint64_t>(b - a) + 1U;
  return a + static_cast<int>(rnext(r) % k);
}
double rexp(rng &r, double rate) { return -std::log(runi(r)) / rate; }
int rpois(rng &r, double mu) {
  const double lim = std::exp(-mu);
  int k = 0;
  for (double p = runi(r); p > lim; p *= runi(r)) ++k;
  return k;
}

bool tk_rec(tk &t, const char *id, bool ok, double want, double got) {
  ++t.n;
  if (ok) return true;
  // the one failure path: a line in the sink, cut short when it is full
  ++t.nf;
  const std::size_t room = t.cap - t.len;
  const int k = std::snprintf(t.buf + t.len, room, "%s: want %.17g got %.17g\n",
                              id, want, got);
  t.len += std::min(static_cast<std::size_t>(k < 0 ? 0 : k), room - 1);
  return false;
}

bool tk_veq(tk &t, const char *id, const double *got, const double *want,
            int n, double tol) {
  double sd = 0.0, sy = 0.0;
  int m = 0, bad = 0;
  bool nanbad = false;
  for (int i = 0; i < n; ++i) {
    const bool gn = std::isnan(got[i]), wn = std::isnan(want[i]);
    if (gn != wn && !nanbad) bad = i;
    nanbad = nanbad || gn != wn;
    if (gn || wn || got[i] == want[i]) continue;
    if (m == 0 && !nanbad) bad = i;
    sd += std::fabs(got[i] - want[i]);
    sy += std::fabs(want[i]);
    ++m;
  }
  // the mean differences over the m elements that differ
  double d = m > 0 ? sd / m : 0.0;
  const double y = m > 0 ? sy / m : 0.0;
  if (std::isfinite(y) && y > tol) d /= y;
  const bool ok = !nanbad && (m == 0 || d < tol);
  return tk_rec(t, id, ok, want[bad], got[bad]);
}

double tk_q7(double *x, int n, double p) {
  std::sort(x, x + n);
  const double h = (n - 1) * p;
  const int lo = static_cast<int>(std::floor(h));
  const int hi = lo + 1 < n ? lo + 1 : lo;
  const double f = h - lo;
  return (1.0 - f) * x[lo] + f * x[hi];
}

namespace {

// The framework checks itself: a failing comparison of each kind must be
// recorded, with its id, in a sink of its own, and a full sink must cut the
// line short without writing past its end.
void t_self(tk &t) {
  char buf[64];
  tk s = {0, 0, buf, sizeof buf, 0};
  const double one = 1.0, two = 2.0, nan = NAN;
  tk_eq(s, "s1", 1.0, 1.0 + 1e-6, 1e-7);
  tk_veq(s, "s2", &nan, &one, 1, 1e-7);
  tk_lt(s, "s3", two, one);
  tk_ieq(s, "s4", 3, 4);
  tk_ieq(t, "self.n", static_cast<int>(s.n), 4);
  tk_ieq(t, "self.nf", static_cast<int>(s.nf), 4);
  tk_true(t, "self.id", buf[0] == 's' && buf[1] == '1');
  tk_ieq(t, "self.len", static_cast<int>(s.len), 63);
  tk_true(t, "self.nul", buf[63] == '\0');
  // and passes what testthat passes
  const double a[2] = {1.0, 1e300}, b[2] = {1.0 + 1e-9, 1e300};
  tk_true(t, "self.ok", tk_veq(s, "s5", a, b, 2, 1e-8));
  tk_ieq(t, "self.nf2", static_cast<int>(s.nf), 4);
  double q[5] = {5.0, 1.0, 4.0, 2.0, 3.0};
  tk_eq(t, "self.q7", tk_q7(q, 5, 0.9), 4.6, 1e-15);
  // the comparisons and the draws behave
  tk_le(t, "self.le", 1.0, 1.0);
  tk_ge(t, "self.ge", 1.0, 1.0);
  tk_gt(t, "self.gt", 2.0, 1.0);
  tk_nan(t, "self.nan", nan);
  tk_true(t, "self.near", tk_near(s, "s6", 0.5 + 1e-9, 0.5, 2e-9) &&
                              tk_near(s, "s7", 2e3 + 1e-6, 2e3, 1e-9) &&
                              !tk_near(s, "s8", 0.5 + 1e-8, 0.5, 2e-9));
  rng g = {1};
  double lo = 1, hi = 0, se = 0, sp = 0;
  int kmin = 9, kmax = -9;
  for (int i = 0; i < 1000; ++i) {
    const double u = runi(g);
    lo = std::min(lo, u);
    hi = std::max(hi, u);
    const int k = rsam(g, -3, 6);
    kmin = std::min(kmin, k);
    kmax = std::max(kmax, k);
    se += rexp(g, 2) / 1000;
    sp += rpois(g, 3) / 1000.0;
  }
  tk_true(t, "self.runi", lo > 0 && hi < 1 && runi(g, 2, 3) >= 2);
  tk_true(t, "self.rsam", kmin == -3 && kmax == 6);
  // means within five standard errors: 0.5 +- 0.08, 3 +- 0.28
  tk_eq(t, "self.rexp", se, 0.5, 0.16);
  tk_eq(t, "self.rpois", sp, 3, 0.093);
}

}  // namespace

//' @srrstats {G5.4} Every kernel is checked against values computed
//'   independently in src/t_*.cpp: by hand, from closed forms or by brute
//'   force.
//' @srrstats {G5.4a} The references include exhaustive enumerations: every
//'   paternity trio from three alleles (t_frns.cpp), every filling of the
//'   unobserved outcomes (sellab(), sntbnd() in t_corr.cpp).
//' @srrstats {G5.5} Random cases come from splitmix64 streams with fixed
//'   seeds (rng g = {...}), so every run draws the same numbers.
//' @noRd
extern "C" int clab_selftest(char *sink, int cap, int *nfail) {
  sink[0] = '\0';
  tk t = {0, 0, sink, static_cast<std::size_t>(cap), 0};
  t_self(t);
  t_disp(t);
  t_frns(t);
  t_pol(t);
  t_corr(t);
  t_geo(t, 0);
  t_cgr(t, 0);
  t_nnd(t, 0);
  t_rpk(t, 0);
  t_jtc(t, 0);
  t_hsp(t, 0);
  t_stm(t, 0);
  t_fair(t, 0);
  t_calaud(t, 0);
  t_fsim(t, 0);
  t_oax(t, 0);
  t_tlv(t, 0);
  t_aor(t, 0);
  t_bnf(t, 0);
  t_tfr(t, 0);
  t_sk(t, 0);
  t_rpt(t, 0);
  t_sens(t, 0);
  t_trj(t, 0);
  t_mvt(t, 0);
  t_lws(t, 0);
  t_chm(t, 0);
  t_gst(t, 0);
  t_spt(t, 0);
  t_orc(t, 0);
  t_net(t, 0);
  t_geo2(t, 0);
  t_swt(t, 0);
  t_svg(t, 0);
  t_lab(t, 0);
  t_cnt(t, 0);
  t_rat(t, 0);
  t_tab(t, 0);
  t_learn(t, 0);
  t_rf(t, 0);
  t_dml(t, 0);
  t_sp2(t, 0);
  t_did(t, 0);
  t_knx(t, 0);
  t_hot(t, 0);
  t_kscn(t, 0);
  t_sreg(t, 0);
  t_ppx(t, 0);
  t_geou(t, 0);
  t_buf(t, 0);
  t_hwk(t, 0);
  t_trn(t, 0);
  t_kmr(t, 0);
  t_hyq(t, 0);
  t_sdb(t, 0);
  t_mhk(t, 0);
  t_sev(t, 0);
  t_cnc(t, 0);
  t_dkf(t, 0);
  t_pub(t, 0);
  t_jdg(t, 0);
  t_lcr(t, 0);
  *nfail = static_cast<int>(t.nf);
  return static_cast<int>(t.n);
}

#ifdef CLAB_MAIN
// The standalone program: the self-test, or with "bench" the time per
// element of each kernel. Nothing here is in the R package's build.
#include <chrono>
#include <cstdlib>
#include <cstring>

#include "crimlab.h"

#ifdef CLAB_WRAP
// Linked with -Wl,--wrap=malloc,--wrap=calloc,--wrap=realloc: every
// allocation made from the objects of this program is counted.
extern "C" void *__real_malloc(std::size_t n);
extern "C" void *__real_calloc(std::size_t n, std::size_t m);
extern "C" void *__real_realloc(void *p, std::size_t n);
static long nalloc = 0;
extern "C" void *__wrap_malloc(std::size_t n) {
  ++nalloc;
  return __real_malloc(n);
}
extern "C" void *__wrap_calloc(std::size_t n, std::size_t m) {
  ++nalloc;
  return __real_calloc(n, m);
}
extern "C" void *__wrap_realloc(void *p, std::size_t n) {
  ++nalloc;
  return __real_realloc(p, n);
}
#endif

namespace {

double now() {
  using c = std::chrono::steady_clock;
  return std::chrono::duration<double>(c::now().time_since_epoch()).count();
}

// Nanoseconds per element of f() over n elements, repeated for 50 ms.
template <class F> double per(F f, double n) {
  long reps = 0;
  const double t0 = now();
  double t;
  do {
    f();
    ++reps;
  } while ((t = now() - t0) < 0.05);
  return t / reps / n * 1e9;
}

int bench() {
  const int big = 1000000;
  double *x = static_cast<double *>(std::malloc(4 * 2 * sizeof(double) * big));
  int *g = static_cast<int *>(std::malloc(3 * 2 * sizeof(int) * big));
  if (x == nullptr || g == nullptr) return 2;
  double *a = x, *b = x + 2 * big, *c = x + 4 * big, *o = x + 6 * big;
  int *i0 = g, *i1 = g + 2 * big, *i2 = g + 4 * big;
  int nint[2];
  std::printf("%-8s %12s %12s %12s  ns per element\n", "kernel", "n=1e3",
              "n=1e5", "n=1e6");
  const char *names[] = {"vod",    "outtst", "dspbm",  "dyadrt", "gtmp",
                         "cllr",   "lrcal",  "kinidx", "patidx", "clrsmp",
                         "expbnd", "rtmhs",  "incap",  "incapr", "sellab",
                         "losbnd", "sntbnd"};
  for (int k = 0; k < 17; ++k) {
    std::printf("%-8s", names[k]);
    for (int n : {1000, 100000, 1000000}) {
      for (int i = 0; i < 2 * n; ++i) {
        a[i] = 1 + i % 7;
        b[i] = 2 + i % 5;
        c[i] = 0.1 + 0.1 * (i % 2);
        i0[i] = i / 2 + i % 2;
        i1[i] = i / 2 + 2 * (i % 2);
        i2[i] = i % 2;
      }
      double ns = 0;
      switch (k) {
        case 0: ns = per([&] { clab_vod(a, b, n, 0, nullptr, o, o, o, i2, o); }, n); break;
        case 1:
          for (int i = 0; i < n; ++i) c[i] = (i % 100) / 100.0;
          ns = per([&] { clab_outtst(c, a, n, 0.5, o, o + 1); }, n);
          break;
        case 2: ns = per([&] { clab_dspbm(a, b, a, n, 0, o, o, o); }, n); break;
        case 3:
          for (int i = 0; i < n; ++i) {
            i0[i] = i % 100;
            i1[i] = (i / 100) % 100;
          }
          ns = per([&] { clab_dyadrt(a, i0, i1, n, b, 100, o, o, o, o); }, n);
          break;
        case 4: ns = per([&] { clab_gtmp(c, c, i2, n, 0.01, o, o + n); }, n); break;
        case 5: ns = per([&] { clab_cllr(a, n / 2, b, n / 2, o); }, n); break;
        case 6:
          for (int i = 0; i < n; ++i) a[i] = 1 + i;
          ns = per([&] { clab_lrcal(a, i2, n, o, b, o + n); }, n);
          break;
        case 7: {
          const double kk[3] = {0.25, 0.5, 0.25};
          ns = per([&] { clab_kinidx(i0, i1, c, n, kk, o); }, n);
          break;
        }
        case 8: ns = per([&] { clab_patidx(i0, i1, i0, c, n, o); }, n); break;
        case 9: ns = per([&] { clab_clrsmp(a, b, a, b, n, o, o, o, o, o, o, o, nint); }, n); break;
        case 10: ns = per([&] { clab_expbnd(a, b, n, 0, 1.5, o, o, o, o, o, o, i2 + n); }, n); break;
        case 11: ns = per([&] { clab_rtmhs(a, b, c, n, 3, o); }, n); break;
        case 12: ns = per([&] { clab_incap(a, c, b, n, nullptr, o, o, o, o); }, n); break;
        case 13:
          for (int i = 0; i < n; ++i) a[i] = 1.0 / (1 + i);
          ns = per([&] { clab_incapr(a, n, 0, n / 2, 0.5, o); }, n);
          break;
        case 14:
          for (int i = 0; i < n; ++i) a[i] = i % 3 == 0;
          ns = per([&] { clab_sellab(a, i2, i0, b, n, o); }, n);
          break;
        case 15: ns = per([&] { clab_losbnd(a, n / 2, b, n / 2, o); }, n); break;
        default:
          for (int i = 0; i < n; ++i) a[i] = i % 3 == 0;
          ns = per([&] { clab_sntbnd(a, i2, b, n, o); }, n);
          break;
      }
      std::printf(" %12.2f", ns);
    }
    std::printf("\n");
  }
  // the scalar kernels, per call; the mixture by the number of alleles to
  // carry, m: inclusion-exclusion has 2^m terms
  std::printf("orrisk   %10.2f ns per call\n",
              per([&] { clab_orrisk(3, 0.4, 0.5, o, o + 1); }, 1));
  std::printf("misrd    %10.2f ns per call\n",
              per([&] { clab_misrd(0.3, 0.1, 0.4, 0.8, 0.9, o); }, 1));
  const double ah[2] = {0.5, 0.8}, al[2] = {0.9, 0.9};
  std::printf("hzdsel   %10.2f ns per call\n",
              per([&] { clab_hzdsel(0.5, 0.5, 0.1, ah, al, o); }, 1));
  for (int i = 0; i < 20; ++i) {
    a[i] = 0.04;
    i0[i] = 1;
  }
  for (int m : {4, 8, 12, 16, 20}) {
    std::printf("mixp m=%-2d %10.2f ns per call\n", m,
                per([&] { clab_mixp(a, i0, m, 10, o); }, 1));
  }
  std::free(x);
  std::free(g);
  return 0;
}

}  // namespace

double tk_med5(double *t) {
  std::sort(t, t + 5);
  return t[2];
}

int main(int argc, char **argv) {
  if (argc > 1 && std::strcmp(argv[1], "bench") == 0) {
    const int rc = bench();
    b_geo();
    b_cgr();
    b_nnd();
    b_rpk();
    b_jtc();
    b_hsp();
    b_stm();
    b_fair();
    b_tlv();
    b_aor();
    b_bnf();
    b_tfr();
    b_sk();
    b_rpt();
    b_sens();
    b_trj();
    b_mvt();
    b_chm();
    b_gst();
    b_spt();
    return rc;
  }
  if (argc > 1 && std::strcmp(argv[1], "ext") == 0) {
    // the spatial and dossier families with their full-size simulations
    static char xs[1 << 16];
    tk t = {0, 0, xs, sizeof xs, 0};
    xs[0] = '\0';
    t_geo(t, 1);
    t_cgr(t, 1);
    t_nnd(t, 1);
    t_rpk(t, 1);
    t_jtc(t, 1);
    t_hsp(t, 1);
    t_stm(t, 1);
    t_fair(t, 1);
    t_calaud(t, 1);
    t_fsim(t, 1);
    t_oax(t, 1);
    t_tlv(t, 1);
    t_aor(t, 1);
    t_bnf(t, 1);
    t_tfr(t, 1);
    t_sk(t, 1);
    t_rpt(t, 1);
    t_sens(t, 1);
    t_trj(t, 1);
    t_mvt(t, 1);
    t_lws(t, 1);
    t_chm(t, 1);
    t_gst(t, 1);
    t_spt(t, 1);
    t_orc(t, 1);
    t_net(t, 1);
    t_geo2(t, 1);
    t_swt(t, 1);
    t_svg(t, 1);
    t_lab(t, 1);
    t_cnt(t, 1);
    t_rat(t, 1);
    t_tab(t, 1);
    t_learn(t, 1);
    t_rf(t, 1);
    t_dml(t, 1);
    t_sp2(t, 1);
    t_knx(t, 1);
    t_hot(t, 1);
    t_kscn(t, 1);
    t_sreg(t, 1);
    t_ppx(t, 1);
    t_geou(t, 1);
    t_buf(t, 1);
    t_hwk(t, 1);
    t_trn(t, 1);
    t_kmr(t, 1);
    t_hyq(t, 1);
    t_sdb(t, 1);
    t_mhk(t, 1);
    t_sev(t, 1);
    t_cnc(t, 1);
    t_dkf(t, 1);
    t_pub(t, 1);
    t_jdg(t, 1);
    t_lcr(t, 1);
    std::fputs(xs, stdout);
    std::printf("%ld checks, %ld failed\n", t.n, t.nf);
    return t.nf == 0 ? 0 : 1;
  }
  static char sink[4096];
  int nf = 0;
#ifdef CLAB_WRAP
  nalloc = 0;
#endif
  const int n = clab_selftest(sink, static_cast<int>(sizeof sink), &nf);
  std::fputs(sink, stdout);
  // a report cut short at the end of the sink has no final newline
  const std::size_t sl = std::strlen(sink);
  if (sl > 0 && sink[sl - 1] != '\n') std::fputc('\n', stdout);
#ifdef CLAB_WRAP
  std::printf("allocations during the self-test: %ld\n", nalloc);
#endif
  std::printf("%d checks, %d failed\n", n, nf);
  return nf == 0 ? 0 : 1;
}
#endif
