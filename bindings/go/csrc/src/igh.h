// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// First Nations, Metis and Inuit health and wellbeing, Statistics Canada's
// published tables, the 2021 Census Indigenous Population Profile, the
// Truth and Reconciliation Commission's Calls to Action, and the
// simulated policy panel polex(): the checks and the arithmetic as C++
// shared by the R binding (src/bigh.cpp) and clab_jrun() (cli/jigh.inc).
// The arguments come in as values of one tree type (V) that either front
// end builds from its own (an R object, a JSON value), and the result goes
// out the same way; a refusal is a message without the function's name,
// which each front end puts in front. Defined in igv.cpp (the values and
// helpers), igh.cpp, igf.cpp, igm.cpp and igc.cpp (the functions).
#ifndef CLAB_IGH_H
#define CLAB_IGH_H

#include <string>
#include <utility>
#include <vector>

namespace igh {

// An R value: NULL, a double, integer, text, logical or factor vector, or
// a list (a data frame when df). Doubles keep R's NA (the bits of
// NA_real_) apart from NaN; integers, logicals and factor codes are NA at
// INT_MIN; text is NA where na is set. A factor's codes count from 1 into
// its levels, held in s. A matrix has nr rows and nc columns, filled by
// column; rn and cn are its dimnames. at holds the other attributes of a
// result (the R binding sets them; JSON has none). js marks a value read
// from JSON, which carries no class; given, an argument an R caller gave
// (an explicit NULL is not the default, as in R).
struct V {
  enum T { NUL, DBL, INT, STR, LGL, FAC, LIST } t = NUL;
  std::vector<double> d;
  std::vector<int> i;
  std::vector<std::string> s;
  std::vector<char> na;
  std::vector<std::string> nm;
  std::vector<V> a;
  bool df = false;
  std::vector<std::string> cls;
  std::vector<std::pair<std::string, V>> at;
  int nr = -1, nc = -1;
  std::vector<std::string> rn, cn;
  bool fn = false;   // a function was given
  bool env = false;  // an environment was given
  bool odd = false;  // another kind of R object (complex, raw, ...)
  bool js = false;
  bool given = false;  // an argument the R caller gave, NULL included
  std::size_t len() const;
};

// R's NA_real_.
double na();
bool isna(double x);

// The functions by name, and their arguments in the order of the R
// function (a trailing ? marks an optional one); nf of them.
struct Fn {
  const char *name;
  const char *args;
};
extern const Fn fns[];
extern const int nfns;

// Runs fn on the arguments an (names) av (values): CLAB_OK with out the
// result, or CLAB_EARG with err the message; warn receives the warnings R
// gives, in order. fn must be one of fns.
int call(const char *fn, const std::vector<std::string> &an, const std::vector<V> &av, V &out,
         std::string &err, std::vector<std::string> &warn);

}  // namespace igh

#endif
