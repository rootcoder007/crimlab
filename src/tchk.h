// SPDX-License-Identifier: AGPL-3.0-or-later
//
// The self-test's check framework. A check has a short id; a failed check
// appends one line "id: want <expected> got <observed>" to a sink the caller
// owns. No allocation, no exceptions, no I/O. Random cases come from a
// splitmix64 stream with a fixed seed, so every run sees the same draws.
#ifndef CLAB_TCHK_H
#define CLAB_TCHK_H

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

// What the kernels and this test assume of the platform, checked here
// because core.cpp states it nowhere: IEEE-754 doubles and an int of at
// least 32 bits.
static_assert(std::numeric_limits<double>::is_iec559,
              "crimlab needs IEEE-754 double");
static_assert(std::numeric_limits<int>::digits >= 31,
              "crimlab needs an int of at least 32 bits");

struct tk {
  long n;          // checks run
  long nf;         // checks failed
  char *buf;       // failure lines, NUL-terminated
  std::size_t cap; // size of buf, at least 1
  std::size_t len; // bytes used, not counting the NUL
};

// Count one check; on failure record id, the expected and the observed value.
bool tk_rec(tk &t, const char *id, bool ok, double want, double got);

// testthat's expect_equal(got, want, tolerance = tol), as waldo compares
// numbers: elements equal or NaN in both are set aside; over the rest the mean
// absolute difference, divided by the mean absolute expected value when that
// exceeds tol, must be below tol.
bool tk_veq(tk &t, const char *id, const double *got, const double *want,
            int n, double tol);
bool tk_eq(tk &t, const char *id, double got, double want, double tol);
bool tk_true(tk &t, const char *id, bool c);
bool tk_lt(tk &t, const char *id, double a, double b);
bool tk_le(tk &t, const char *id, double a, double b);
bool tk_gt(tk &t, const char *id, double a, double b);
bool tk_ge(tk &t, const char *id, double a, double b);
bool tk_nan(tk &t, const char *id, double x);
bool tk_ieq(tk &t, const char *id, int got, int want);

// testthat's default tolerance, sqrt(.Machine$double.eps).
constexpr double tk_tol = 1.4901161193847656e-08;

// splitmix64.
struct rng {
  std::uint64_t s;
};
std::uint64_t rnext(rng &r);
// Uniform on (0, 1), never 0 or 1.
double runi(rng &r);
double runi(rng &r, double a, double b);
// A whole number from a to b; the bias of the modulus is below 2^-40 here.
int rsam(rng &r, int a, int b);
double rexp(rng &r, double rate);
// Poisson by multiplying uniforms (Knuth); means here stay far below 700.
int rpois(rng &r, double mu);

// The sample quantile of type 7, R's default; x is sorted in place.
double tk_q7(double *x, int n, double p);

// The four families, each adding its checks to t.
void t_disp(tk &t);
void t_frns(tk &t);
void t_pol(tk &t);
void t_corr(tk &t);

// The spatial families: big = 1 runs the extended simulations (the
// standalone program's "ext" mode), 0 the default checks; b_<family>
// prints the timings (bench mode, standalone only).
void t_geo(tk &t, int big);
void b_geo();
void t_cgr(tk &t, int big);
void b_cgr();
void t_nnd(tk &t, int big);
void b_nnd();
void t_rpk(tk &t, int big);
void b_rpk();
void t_jtc(tk &t, int big);
void b_jtc();
void t_hsp(tk &t, int big);
void b_hsp();

#endif
