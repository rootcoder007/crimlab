// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The helpers of the Indigenous health, published-table and census
// functions (igh.h): reading arguments in the words of the R package,
// labels as R's as.character() writes them, the arithmetic R's stats
// package does for them (quantiles, ranks, correlation, rounding) and a
// generator that draws what R's Mersenne-Twister, inversion and rejection
// sampling draw after set.seed(), so a seeded result is R's. Internal to
// igv.cpp and the files of the functions.
#ifndef CLAB_IGK_H
#define CLAB_IGK_H

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "igh.h"

namespace igh {

// A refusal: thrown, caught by call().
struct Stop {
  std::string msg;
};
[[noreturn]] void stop(const char *fmt, ...)
#if defined(__GNUC__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;
std::string sf(const char *fmt, ...)
#if defined(__GNUC__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;

// One call: the arguments by name and the warnings given.
struct Cx {
  const std::vector<std::string> *an;
  const std::vector<V> *av;
  std::vector<std::string> *warn;
  const V &operator[](const char *k) const;
  bool has(const char *k) const { return (*this)[k].t != V::NUL; }
  void w(const std::string &m) const { warn->push_back(m); }
};

extern const V nulv;

// Text as R's trimws(), tolower() of ASCII letters.
std::string trim(const std::string &s);
std::string lower(const std::string &s);
// as.character() of one double, as R 4.3 and later writes it.
std::string rch(double x);
// "\"a\", \"b\"" and "a, b" (toString()).
std::string quoted(const std::vector<std::string> &v);
std::string tostr(const std::vector<std::string> &v);

// Labels: as.character() of an atomic vector, with its missing values.
struct Lab {
  std::vector<std::string> s;
  std::vector<char> na;
  std::size_t n() const { return s.size(); }
};
Lab lab(const V &v);
// trimws(as.character(x)), NA kept as NA.
Lab tlab(const V &v);

// The checks of the R package's arguments, each in its words.
void nofn(const V &v, const char *what);
// the argument was left out (JSON null, or not given in R): its default
bool dflt(const V &v);
std::vector<double> num(const V &v, const char *what, bool integer = false);
double scalar(const V &v, const char *what, double lo, double hi, bool integer, double def);
double conf(const V &v, const char *what, double def);
std::string str1(const V &v, const char *what);
// one of the choices (NULL: the first)
std::string pick(const V &v, const char *what, const std::vector<std::string> &ch);
// a character vector with no NA (NULL allowed when null)
bool isstr(const V &v);

// Data frames.
bool isdf(const V &v);
std::size_t nrow(const V &v);
const V *col(const V &df, const std::string &nm);
int coli(const V &df, const std::string &nm);
// .rmbl_col(): the name of a column of data, as one string
std::string colnm(const V &df, const V &x, const char *what, const char *fn = nullptr);
// a data frame from columns
V frame(const std::vector<std::string> &nm, std::vector<V> cols);
// rows r of a data frame (or of a vector)
V rows(const V &x, const std::vector<std::size_t> &r);
// the class of a value, its inherits()
bool inh(const V &v, const char *c);

// Builders.
V dv(const std::vector<double> &x);
V dv1(double x);
V iv(const std::vector<int> &x);
V sv(const std::vector<std::string> &x);
V sv1(const std::string &x);
V svna(const std::vector<std::string> &x, const std::vector<char> &na);
V lv(const std::vector<int> &x);
V fac(const std::vector<std::string> &x, const std::vector<std::string> &lev);
V nlist(const std::vector<std::string> &nm, std::vector<V> el);
V named(V v, const std::vector<std::string> &nm);

// Quality flags: .rmbl_qflags() with its warning.
std::vector<std::string> qflags(const Cx &c, const Lab &q);
bool gone(const std::string &q);
// the quality column, or "ok" for each row
std::vector<std::string> quality(const Cx &c, const V &data, const std::string &qn, std::size_t n);
// .rmbl_contains(): estimates whose interval does not hold them are NA,
// with a warning naming them
void contains(const Cx &c, const char *fn, std::vector<double> &est, const std::vector<double> &lo,
              const std::vector<double> &hi, const std::vector<std::string> &labs, const std::vector<char> *use);
// .rmbl_pan_indigenous(): rows that pool groups
std::vector<int> pooledrows(const std::vector<std::string> &g, const V &pooled);
void pooledtype(const V &pooled);
void pooledcheck(const V &pooled, const std::vector<std::string> &g);

// Arithmetic of R's.
double qn(double conflvl);  // qnorm(1 - (1 - conflvl) / 2)
double rmean(const std::vector<double> &x);
double rsum(const std::vector<double> &x);
double quant7(std::vector<double> x, double p);
double median(std::vector<double> x);
std::vector<double> rank(const std::vector<double> &x);
double cor(const std::vector<double> &x, const std::vector<double> &y);
double rround(double x, int dig);
double plogis(double x);
double qlogis(double p);
std::vector<std::string> uniq(const std::vector<std::string> &x);

// R's random numbers: set.seed(seed) (Mersenne-Twister, inversion,
// rejection sampling), then unif_rand(), norm_rand(), exp_rand(),
// R_unif_index() and rpois() as R draws them.
class Rng {
 public:
  explicit Rng(std::uint32_t seed);
  double unif();
  double norm();
  double expo();
  double index(double dn);
  double pois(double mu);

 private:
  std::uint32_t mt_[625];
  int mti_;
  double genrand();
  double rbits(int bits);
  // rpois's saved state
  double muprev_ = 0, muprev2_ = 0, s_ = 0, d_ = 0, bigl_ = 0, omega_ = 0, b1_ = 0, b2_ = 0, c_ = 0, c0_ = 0,
         c1_ = 0, c2_ = 0, c3_ = 0, pp_[36] = {}, p0_ = 0, p_ = 0, q_ = 0;
  int l_ = 0, m_ = 0;
};
// the seed argument: NULL draws none (has false), else a whole number
bool seedarg(const V &v, std::uint32_t &s);

}  // namespace igh

#endif
