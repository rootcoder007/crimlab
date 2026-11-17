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
thread_local bool oom = false;

}  // namespace

namespace kutil {

void oomset() noexcept { oom = true; }

bool oomtake() noexcept {
  const bool r = oom;
  oom = false;
  return r;
}

bool stop(double w) {
  if (hit) return true;
  done += w;
  if (done < every) return false;
  done = 0.0;
  hit = hook != nullptr && hook() != 0;
  return hit;
}

}  // namespace kutil

extern "C" {

void clab_onintr(int (*f)(void), double ev) CLAB_TRY {
  hook = f;
  every = ev >= 0.0 ? ev : 4194304.0;
  done = 0.0;
  hit = false;
} CLAB_CATCHV

int clab_intrd(int clr) CLAB_TRY {
  const int r = hit ? 1 : 0;
  if (clr != 0) {
    hit = false;
    done = 0.0;
  }
  return r;
} CLAB_CATCH(CLAB_ENOMEM)

}  // extern "C"
