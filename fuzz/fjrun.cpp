// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// libFuzzer target of clab_jrun(), the entry every binding and the
// command-line program call: the first line of the input is the function
// name, the rest its JSON arguments. Whatever the bytes, the call must
// return (no crash, no read or write out of bounds, no undefined
// behaviour, no leak) and leave a NUL-terminated answer of the size it
// reports.
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>

#include "../cli/jrun.h"

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t *d, std::size_t n) {
  const char *p = reinterpret_cast<const char *>(d);
  const void *nl = std::memchr(p, 0x0a, n);
  const std::size_t k = nl ? static_cast<std::size_t>(static_cast<const char *>(nl) - p) : n;
  const std::string fn(p, k);
  const std::string in(nl ? p + k + 1 : p + n, nl ? n - k - 1 : 0);
  static char out[1 << 16];
  int need = 0;
  const int rc = clab_jrun(fn.c_str(), in.c_str(), out, static_cast<int>(sizeof out), &need);
  (void)rc;
  if (need > 0 && need <= static_cast<int>(sizeof out)) {
    // the answer is a string of the size reported
    if (std::strlen(out) + 1 > static_cast<std::size_t>(need)) std::abort();
  }
  return 0;
}
