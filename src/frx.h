// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The checks and the arithmetic of pedlr() as C++ shared by
// the R binding (src/bfrx.cpp) and clab_jrun() (cli/jfrx.inc): the
// arguments come in as columns of labels and numbers, the results go out
// as vectors, and a refusal is a message without the function's name,
// which each front end puts in front. Defined in ped.cpp.
#ifndef CLAB_FRX_H
#define CLAB_FRX_H

#include <string>
#include <vector>

namespace frx {

// A column as text (as.character()), with its missing values.
struct Lab {
  std::vector<std::string> s;
  std::vector<char> na;
};

// pedlr(): the pedigree (ped) and the alternative (alt, or every typed
// person unrelated when hasalt is false) as id, father and mother; the
// genotypes as id, locus and two alleles; the frequencies as locus,
// allele and frequency; the mutation model by name, its rates (maternal
// then paternal, or one for both) and the range of the stepwise model.
struct PedIn {
  Lab id, fid, mid;
  bool hasalt = false;
  Lab aid, afid, amid;
  Lab gid, gloc, ga1, ga2;
  Lab floc, fal;
  std::vector<double> fq;
  std::string mut;
  std::vector<double> rate;
  double range = 0.5;
};

// By locus: the base-10 log likelihoods under ped (l1) and alt (l2) and
// their ratio; log10lr, the base-10 log of the product of the ratios.
struct PedOut {
  std::vector<std::string> loc;
  std::vector<double> l1, l2, lr;
  double log10lr = 0.0;
};

// CLAB_OK, or an error code with err the message.
int pedlr(const PedIn &in, PedOut &out, std::string &err);

}  // namespace frx

#endif
