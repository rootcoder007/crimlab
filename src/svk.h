// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The scaling, ideal-point and count-dynamics kernels (svt.cpp, svb.cpp,
// tsd.cpp) take their arguments and give their results as one small tree
// of named nodes, so that the R entry (bsvt.cpp) and clab_jrun()
// (cli/jsvt.inc) each convert once for all of them. The dense algebra
// (symmetric eigenproblem, singular values, Cholesky, LU, least squares)
// and the random numbers live in svla.cpp. Nothing here knows about R.
#ifndef CLAB_SVK_H
#define CLAB_SVK_H

#include <cstdint>
#include <string>
#include <vector>

namespace svk {

// One argument or one part of a result. NUM: doubles (NaN for a missing
// value) with dim empty for a vector, two or three extents (column-major)
// for a matrix or an array; INT and LGL: whole numbers and truth values
// held as doubles; STR: strings; LIST: the named (or unnamed) nodes l.
struct Node {
  enum K { NUL, NUM, INT, LGL, STR, LIST } k = NUL;
  std::vector<double> d;
  std::vector<int> dim;
  std::vector<std::string> s;
  std::vector<std::string> nm;
  std::vector<Node> l;
  const Node *at(const char *name) const;
  Node &add(const char *name, Node v);
};

Node num(double x);
Node ival(double x);
Node lgl(bool b);
Node str(const std::string &x);
Node vec(const std::vector<double> &x);
Node ivec(const std::vector<double> &x);
Node list();

// A dense column-major matrix.
struct Mat {
  int r = 0, c = 0;
  std::vector<double> v;
  Mat() = default;
  Mat(int nr, int nc, double x = 0.0) : r(nr), c(nc), v(static_cast<std::size_t>(nr) * nc, x) {}
  double &operator()(int i, int j) { return v[static_cast<std::size_t>(j) * r + i]; }
  double operator()(int i, int j) const { return v[static_cast<std::size_t>(j) * r + i]; }
};
Node mnode(const Mat &m);
Node anode(const std::vector<double> &x, int a, int b, int c);

Mat mul(const Mat &a, const Mat &b);
Mat tmul(const Mat &a, const Mat &b);  // t(a) %*% b
Mat mult(const Mat &a, const Mat &b);  // a %*% t(b)
Mat tr(const Mat &a);
Mat eye(int n);
// Symmetric eigenproblem by cyclic Jacobi: values decreasing, vectors in
// the columns, each with its largest component (the first of equals)
// positive: the sign of an eigenvector is arbitrary.
void eigsym(const Mat &a, std::vector<double> &val, Mat &vec);
// Thin singular value decomposition a = u diag(d) v' by one-sided Jacobi,
// d decreasing; the sign of each pair as eigsym() fixes it on v.
void svd(const Mat &a, Mat &u, std::vector<double> &d, Mat &v);
// Upper Cholesky factor (a = r'r); false unless a is positive definite.
bool chol(const Mat &a, Mat &r);
// a^-1 b by LU with partial pivoting; false when a is singular.
bool solve(const Mat &a, const Mat &b, Mat &x);
bool inv(const Mat &a, Mat &x);
// The least-squares solution of x b = y (Householder QR, full column
// rank); false when x is rank deficient.
bool lsq(const Mat &x, const std::vector<double> &y, std::vector<double> &b);
// The Moore-Penrose inverse with singular values above tol kept.
Mat pinv(const Mat &a, double tol);
double det(const Mat &a);

// The random numbers of R: the Mersenne twister seeded as set.seed()
// seeds it, uniforms as unif_rand(), normals by inversion as norm_rand(),
// exponentials as exp_rand(), gammas by Ahrens and Dieter (1974, 1982) as
// rgamma(), and indices as R_unif_index() with rejection sampling; so a
// kernel given a seed draws exactly what R draws after set.seed(seed).
struct rng {
  std::uint32_t mt[624];
  int mti = 624;
  double aa = 0, aaa = 0, s = 0, s2 = 0, dd = 0, q0 = 0, b = 0, si = 0, c = 0;
  explicit rng(int seed);
  double unif();
  double norm();
  double expo();
  double gamma(double a, double scale);
  double chisq(double df) { return gamma(df / 2.0, 2.0); }
  double index(double dn);
};

// Distribution functions as R's: pnorm (with log), qnorm, dnorm, pgamma
// and qgamma (scale 1 / rate).
double pnorm(double x, bool lg = false);
double qnorm(double p);
double dnorm(double x, bool lg = false);
double pgamma(double x, double shape, double rate);
double qgamma(double p, double shape, double rate);
// stats::quantile(type = 7) and the sample standard deviation.
double quant(std::vector<double> x, double p);
double sdv(const double *x, int n, int step = 1);

// The kernel table: R name, arguments in R's order (space separated:
// "name" required, "name?" NULL by default, "name=v" with the default v,
// a number, TRUE, FALSE, numbers separated by commas, or a string) and
// the kernel. A kernel reads a (a LIST) and fills out, or
// writes the reason into err and returns false ("" for an interrupt).
using Kfn = bool (*)(const Node &a, Node &out, std::string &err);
struct Fn {
  const char *name, *args;
  Kfn f;
};
const Fn *find(const char *name);
extern const Fn table[];
extern const int ntable;
// Run the kernel name on a; false with err set on a refusal.
bool run(const char *name, const Node &a, Node &out, std::string &err);

// Argument readers: false with err set ("`nm` must be ...").
bool getm(const Node &a, const char *nm, Mat &m, std::string &err, bool na = false, bool opt = false);
bool getv(const Node &a, const char *nm, std::vector<double> &v, std::string &err, bool na = false,
          bool opt = false);
bool getd(const Node &a, const char *nm, double &x, std::string &err, double lo, double hi);
bool geti(const Node &a, const char *nm, int &x, std::string &err, int lo);
bool getb(const Node &a, const char *nm, bool &x, std::string &err);
bool has(const Node &a, const char *nm);
std::string msg(const char *f, ...);

// The EM fit of the binary IRT model (svt.cpp), also the start of the
// Clinton-Jackman-Rivers sampler (svb.cpp): priors xp and bp (mean,
// variance), start NULL or a list (intercepts, slopes, ideal points), ab
// for the largest absolute change as the convergence measure.
struct Em {
  Mat x, b;
  std::vector<double> a;
  double ll = 0;
  int it = 0;
  bool conv = false;
};
bool emfit(const Mat &Y, int D, int mi, double tol, const double *xp, const double *bp, const Node *start, bool ab,
           Em &r, std::string &err);

// The kernels' families, defined in svt.cpp, svb.cpp and tsd.cpp.
#define SVK_DECL(n) bool k_##n(const Node &a, Node &out, std::string &err);
#define SVK_ALL(X)                                                                                      \
  X(votagr) X(amscal) X(bbscal) X(optcls) X(dblctr) X(cmdsf) X(smacf) X(nmdsf) X(mdsfit) X(unfstr)     \
  X(mlsmu) X(smunf) X(ipsrec) X(nomutl) X(nompr) X(nomll) X(procrs) X(irtlik) X(irtpst) X(ordocl)       \
  X(anchv) X(indscl) X(nrmvec) X(cutlns) X(dwnom) X(nombt) X(irtem) X(scboot) X(wrdfsh) X(bamscl)       \
  X(bmdsf) X(bunfld) X(cjrirt) X(anomin) X(ordirt) X(dynirt) X(arimfc) X(sarfc) X(catcor) X(ousim)      \
  X(fpdens) X(lvfit)
SVK_ALL(SVK_DECL)
#undef SVK_DECL

}  // namespace svk

#endif
