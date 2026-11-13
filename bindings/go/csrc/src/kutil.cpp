// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab: the shared numerical helpers of kutil.h that are not inline.
#include "kutil.h"

#include "crimlab.h"

namespace {

// The host's check and how much work runs between two calls of it; per
// thread, what has been done since the last call and whether the host
// asked to stop.
int (*hook)(void) = nullptr;
double every = 4194304.0;  // 2^22 operations: a few milliseconds
thread_local double done = 0.0;
thread_local bool hit = false;

}  // namespace

namespace kutil {

bool stop(double w) {
  if (hit) return true;
  done += w;
  if (done < every) return false;
  done = 0.0;
  hit = hook != nullptr && hook() != 0;
  return hit;
}


// A tail of Binomial(n, p) starting at g: P(X >= g) when up, P(X <= g)
// otherwise, for g on the far side of the mean (the terms then decrease
// monotonically away from g). From the log pmf at g by the ratio
// recursion, summed until a term is below 1e-17 of the sum: O(sd) terms.
// Needs p > 0 (both callers have it); p >= 1 is a point mass at n.
double btail(double g, double n, double p, bool up) {
  if (p >= 1.0) return up ? 1.0 : (g >= n ? 1.0 : 0.0);
  const double lt = std::lgamma(n + 1.0) - std::lgamma(g + 1.0) - std::lgamma(n - g + 1.0) +
                    g * std::log(p) + (n - g) * std::log1p(-p);
  double t = std::exp(lt), s = 0.0;
  const double r = p / (1.0 - p);
  for (double k = g; up ? k <= n : k >= 0.0; k += up ? 1.0 : -1.0) {
    s += t;
    if (t < 1e-17 * s) break;
    t *= up ? (n - k) / (k + 1.0) * r : k / (n - k + 1.0) / r;
  }
  return s;
}

}  // namespace kutil

extern "C" {

void clab_onintr(int (*f)(void), double ev) {
  hook = f;
  every = ev >= 0.0 ? ev : 4194304.0;
  done = 0.0;
  hit = false;
}

int clab_intrd(int clr) {
  const int r = hit ? 1 : 0;
  if (clr != 0) {
    hit = false;
    done = 0.0;
  }
  return r;
}

}  // extern "C"
