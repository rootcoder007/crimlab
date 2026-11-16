// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The data behind the plots of R/plotm.R that take more than a column of
// the fitted object: the ROC curve of a record-linkage model over every
// threshold its probabilities take, with the area under it.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

int clab_rocpts(std::size_t n, const double *prob, const int *lk, double *fpr, double *tpr, std::size_t *k,
                double *auc, int *na) CLAB_TRY {
  if (k == nullptr || auc == nullptr || na == nullptr) return CLAB_EARG;
  if (n > 0 && (prob == nullptr || lk == nullptr || fpr == nullptr || tpr == nullptr)) return CLAB_EARG;
  *k = 0;
  *auc = 0.0;
  na[0] = 0;
  na[1] = 0;
  std::vector<double> th;
  th.reserve(n);
  std::size_t cnt[2] = {0, 0};
  for (std::size_t i = 0; i < n; ++i) {
    const bool miss = std::isnan(prob[i]);
    if (!miss) th.push_back(prob[i]);
    if (lk[i] != 0 && lk[i] != 1) {
      // an NA among the flags is a missing value in both groups
      na[0] = 1;
      na[1] = 1;
      continue;
    }
    const std::size_t g = static_cast<std::size_t>(lk[i]);
    ++cnt[g];
    if (miss) na[g] = 1;
  }
  // sort(unique(prob), decreasing = TRUE)
  std::sort(th.begin(), th.end(), [](double a, double b) { return a > b; });
  th.erase(std::unique(th.begin(), th.end(), [](double a, double b) { return a == b; }), th.end());
  const std::size_t m = th.size();
  // the count of each group at or above each threshold, in one pass down
  // the probabilities sorted the same way
  std::vector<std::size_t> o;
  o.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isnan(prob[i]) && (lk[i] == 0 || lk[i] == 1)) o.push_back(i);
  }
  std::sort(o.begin(), o.end(), [prob](std::size_t a, std::size_t b) { return prob[a] > prob[b]; });
  std::size_t up[2] = {0, 0};
  std::size_t j = 0;
  fpr[0] = 0.0;
  tpr[0] = 0.0;
  for (std::size_t t = 0; t < m; ++t) {
    while (j < o.size() && prob[o[j]] >= th[t]) {
      ++up[static_cast<std::size_t>(lk[o[j]])];
      ++j;
    }
    // mean() of a logical vector: the count over n in long double, as R
    const long double f0 = static_cast<long double>(up[0]) / static_cast<long double>(cnt[0]);
    const long double f1 = static_cast<long double>(up[1]) / static_cast<long double>(cnt[1]);
    fpr[t + 1] = static_cast<double>(f0);
    tpr[t + 1] = static_cast<double>(f1);
  }
  *k = m + 1;
  // sum(diff(fpr) * (tpr[-1] + tpr[-k]) / 2), accumulated as R's sum()
  long double s = 0.0L;
  for (std::size_t t = 0; t < m; ++t) {
    const double d = fpr[t + 1] - fpr[t];
    const double h = tpr[t + 1] + tpr[t];
    const double a = d * h;
    s += a / 2;
  }
  *auc = static_cast<double>(s);
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)
