// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Randomisation lists (crimlab.h, rndalc()): simple randomisation, or
// permuted blocks whose sizes are drawn from a list, stratum by stratum
// (Rosenberger and Lachin 2016, ch. 3). One xorshift64* stream (kutil.h)
// for the whole list: a draw of an arm or of a block size, then a
// Fisher-Yates shuffle of the block. Work: O(n).
#include <climits>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

extern "C" int clab_rndalc(int ns, const int *cnt, int k, const int *ratio, const int *blk, int nb,
                           unsigned int seed, int *arm, int *bno) CLAB_TRY {
  if (ns < 1 || k < 2 || nb < 0) return CLAB_EARG;
  long sr = 0;
  for (int a = 0; a < k; ++a) {
    if (ratio[a] < 1 || ratio[a] > 1000000) return CLAB_EARG;
    sr += ratio[a];
  }
  if (sr > 1000000) return CLAB_EARG;
  for (int b = 0; b < nb; ++b) {
    if (blk[b] < 1 || blk[b] % sr != 0) return CLAB_EARG;
  }
  double tot = 0;
  for (int s = 0; s < ns; ++s) {
    if (cnt[s] < 0) return CLAB_EARG;
    tot += cnt[s];
  }
  if (tot > INT_MAX) return CLAB_EARG;
  kutil::rng r(seed);
  std::vector<int> bk;
  int at = 0;
  for (int s = 0; s < ns; ++s) {
    int left = cnt[s], no = 0;
    while (left > 0) {
      if (kutil::stop(nb ? static_cast<double>(blk[0]) : 1.0)) return CLAB_EINTR;
      if (nb == 0) {
        // arm a with chance ratio[a] / sr
        int u = r.below(static_cast<int>(sr)), a = 0;
        while (u >= ratio[a]) u -= ratio[a++];
        arm[at] = a;
        bno[at] = 0;
        ++at;
        --left;
        continue;
      }
      const int sz = blk[r.below(nb)], m = static_cast<int>(sz / sr);
      bk.clear();
      for (int a = 0; a < k; ++a) bk.insert(bk.end(), static_cast<std::size_t>(ratio[a]) * m, a);
      for (int i = sz - 1; i > 0; --i) {
        const int j = r.below(i + 1);
        const int t = bk[static_cast<std::size_t>(i)];
        bk[static_cast<std::size_t>(i)] = bk[static_cast<std::size_t>(j)];
        bk[static_cast<std::size_t>(j)] = t;
      }
      ++no;
      for (int i = 0; i < sz && left > 0; ++i) {
        arm[at] = bk[static_cast<std::size_t>(i)];
        bno[at] = no;
        ++at;
        --left;
      }
    }
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
