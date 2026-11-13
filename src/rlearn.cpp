// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the learning and causal families: each checks its
// arguments, calls one kernel with scratch from R's transient memory and
// builds the R result. Nothing here calls back into the R evaluator.
#include <climits>
#include <cmath>
#include <cstring>

#include "crimlab.h"
#include "learn.h"
#include "rchk.h"

using clab::fail;
using clab::kchk;
using clab::num;

namespace rl {

// x as a double matrix: n rows, p columns, every value finite

SEXP mat(SEXP x, const char *nm, int *n, int *p) {
  // is.matrix(x) && is.numeric(x) && all(is.finite(x)), at least one row
  SEXP d = Rf_getAttrib(x, R_DimSymbol);
  if (!clab::isnum(x) || TYPEOF(d) != INTSXP || XLENGTH(d) != 2 ||
      INTEGER_RO(d)[0] < 1) {
    fail("`%s` must be a numeric matrix with at least one row.", nm);
  }
  SEXP v = PROTECT(clab::dbl(x));
  const double *a = REAL_RO(v);
  for (R_xlen_t i = 0; i < XLENGTH(v); ++i) {
    if (!std::isfinite(a[i])) {
      fail("`%s` must be finite numbers with no missing values.", nm);
    }
  }
  *n = INTEGER_RO(d)[0];
  *p = INTEGER_RO(d)[1];
  UNPROTECT(1);
  return v;
}

int *folds(SEXP f, const char *nm, int n, int *k) {
  // whole numbers 1..K, length n, every one of them used, K >= 2
  SEXP v = PROTECT(num(f, nm, n, 1.0));
  const double *a = REAL_RO(v);
  double mx = 1.0;
  for (int i = 0; i < n; ++i) {
    if (a[i] != std::floor(a[i])) fail("`%s` must be whole numbers.", nm);
    mx = std::fmax(mx, a[i]);
  }
  if (mx < 2.0 || mx > n) fail("`%s` must number at least two folds.", nm);
  const int kk = static_cast<int>(mx);
  int *o = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(n), sizeof(int)));
  char *seen = R_alloc(static_cast<std::size_t>(kk), 1);
  std::memset(seen, 0, static_cast<std::size_t>(kk));
  for (int i = 0; i < n; ++i) {
    o[i] = static_cast<int>(a[i]) - 1;
    seen[o[i]] = 1;
  }
  for (int j = 0; j < kk; ++j) {
    if (!seen[j]) fail("`%s` leaves fold %d empty.", nm, j + 1);
  }
  *k = kk;
  UNPROTECT(1);
  return o;
}

int pick(SEXP s, const char *nm, const char *const *ch, int k) {
  if (TYPEOF(s) == STRSXP && XLENGTH(s) == 1 && STRING_ELT(s, 0) != NA_STRING) {
    const char *c = CHAR(STRING_ELT(s, 0));
    for (int i = 0; i < k; ++i) {
      if (std::strcmp(c, ch[i]) == 0) return i;
    }
  }
  // the choices, quoted and joined
  char buf[256];
  std::size_t len = 0;
  for (int i = 0; i < k && len + 40 < sizeof buf; ++i) {
    len += static_cast<std::size_t>(std::snprintf(
        buf + len, sizeof buf - len, "%s\"%s\"", i == 0 ? "" : ", ", ch[i]));
  }
  fail("`%s` must be one of %s.", nm, buf);
}

// scratch in R's transient memory; R_alloc() refuses a block too large
double *dalloc(double n) {
  return reinterpret_cast<double *>(R_alloc(static_cast<std::size_t>(n), sizeof(double)));
}

int *ialloc(double n) {
  return reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(n), sizeof(int)));
}

void ychk(const double *y, int n, bool bin, const char *nm) {
  bool one = false, zero = false, diff = false;
  for (int i = 0; i < n; ++i) {
    if (bin && y[i] != 0.0 && y[i] != 1.0) {
      fail("`%s` must be 0 or 1 for a binary outcome.", nm);
    }
    one = one || y[i] == 1.0;
    zero = zero || y[i] == 0.0;
    diff = diff || y[i] != y[0];
  }
  if (!diff || (bin && !(one && zero))) fail("`%s` is constant.", nm);
}

// a matrix of n rows and p columns, column names those of x (if any)
SEXP rmat(int n, int p, SEXP x) {
  SEXP m = PROTECT(Rf_allocMatrix(REALSXP, n, p));
  SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
  if (!Rf_isNull(dn) && !Rf_isNull(VECTOR_ELT(dn, 1))) {
    SEXP d = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(d, 0, VECTOR_ELT(dn, 1));
    Rf_setAttrib(m, R_DimNamesSymbol, d);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return m;
}

}  // namespace rl

extern "C" {

SEXP C_clab_r_enpath(SEXP x, SEXP y, SEXP fam, SEXP alpha, SEXP lambda,
                    SEXP nlam, SEXP lmr, SEXP fold, SEXP tol, SEXP maxit,
                    SEXP scr) {
  clab::who = "enpath";
  int n, p;
  SEXP xv = PROTECT(rl::mat(x, "x", &n, &p));
  if (n < 2 || p < 1) fail("`x` must have at least two rows and one column.");
  SEXP yv = PROTECT(num(y, "y", n));
  static const char *const fams[] = {"gaussian", "binomial"};
  const int f = rl::pick(fam, "fam", fams, 2);
  rl::ychk(REAL_RO(yv), n, f == 1, "y");
  const double a = clab::p01(alpha, "alpha");
  int nl;
  double r = 0.0;
  SEXP lam;
  if (Rf_isNull(lambda)) {
    nl = clab::whole(nlam, "nlam", 1);
    r = Rf_isNull(lmr) ? (n > p ? 1e-4 : 1e-2) : clab::one(lmr, "lmr", 0, 1);
    lam = PROTECT(Rf_allocVector(REALSXP, nl));
  } else {
    lam = PROTECT(clab::dcopy(num(lambda, "lambda", -1, 0.0)));
    nl = clab::ilen(static_cast<double>(XLENGTH(lam)), "lambda");
    for (int k = 1; k < nl; ++k) {
      if (REAL_RO(lam)[k] > REAL_RO(lam)[k - 1]) {
        fail("`lambda` must be decreasing.");
      }
    }
  }
  int nf = 0;
  const int *fd = Rf_isNull(fold) ? nullptr : rl::folds(fold, "fold", n, &nf);
  const double th = REAL_RO(num(tol, "tol", 1, 0.0, true))[0];
  const int mi = clab::whole(maxit, "maxit", 1);
  clab::lgl(scr, "scr", 1);
  const int sc = LOGICAL_RO(scr)[0];
  SEXP a0 = PROTECT(Rf_allocVector(REALSXP, nl));
  SEXP beta = PROTECT(rl::rmat(p, nl, x));
  SEXP cv = PROTECT(Rf_allocVector(REALSXP, 2 * static_cast<R_xlen_t>(nl)));
  int best[2], nk = nl;
  const double dp = p, dn = n;
  double *w = rl::dalloc(dp * dp + 6 * dp + 5 * dn + dn * dp + nl * (dp + 1) +
                         nf * (nl + 1.0));
  int *iw = rl::ialloc(2 * dn + dp);
  const int rc = clab_enpath(REAL_RO(xv), n, p, REAL_RO(yv), f, a, REAL(lam), nl,
                            r, th, mi, sc, fd, nf, REAL(a0), REAL(beta), REAL(cv),
                            best, &nk, w, iw);
  if (rc == CLAB_EDOM) {
    fail("no fit: no convergence within `maxit` passes, or `y` constant "
         "once a fold is held out.");
  }
  kchk(rc, nullptr);
  // the lambdas fitted: a path chosen here may stop early
  SEXP lk = PROTECT(Rf_allocVector(REALSXP, nk));
  SEXP ak = PROTECT(Rf_allocVector(REALSXP, nk));
  SEXP bk = PROTECT(rl::rmat(p, nk, x));
  std::memcpy(REAL(lk), REAL_RO(lam), sizeof(double) * nk);
  std::memcpy(REAL(ak), REAL_RO(a0), sizeof(double) * nk);
  std::memcpy(REAL(bk), REAL_RO(beta), sizeof(double) * nk * p);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out, 0, lk);
  SET_VECTOR_ELT(out, 1, ak);
  SET_VECTOR_ELT(out, 2, bk);
  if (fd != nullptr) {
    SEXP cm = PROTECT(Rf_allocVector(REALSXP, nk));
    SEXP cs = PROTECT(Rf_allocVector(REALSXP, nk));
    std::memcpy(REAL(cm), REAL_RO(cv), sizeof(double) * nk);
    std::memcpy(REAL(cs), REAL_RO(cv) + nk, sizeof(double) * nk);
    SET_VECTOR_ELT(out, 3, cm);
    SET_VECTOR_ELT(out, 4, cs);
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(REAL_RO(lam)[best[0]]));
    SET_VECTOR_ELT(out, 6, Rf_ScalarReal(REAL_RO(lam)[best[1]]));
    UNPROTECT(2);
  }
  static const char *const nm[] = {"lambda", "a0",   "beta", "cvm",
                                   "cvsd",   "lmin", "l1se"};
  UNPROTECT(10);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_lgtreg(SEXP x, SEXP y, SEXP tol, SEXP maxit) {
  clab::who = "lgtreg";
  int n, p;
  SEXP xv = PROTECT(rl::mat(x, "x", &n, &p));
  SEXP yv = PROTECT(num(y, "y", n));
  rl::ychk(REAL_RO(yv), n, true, "y");
  if (n <= p + 1) fail("`x` needs more rows than columns plus one.");
  const double th = REAL_RO(num(tol, "tol", 1, 0.0, true))[0];
  const int mi = clab::whole(maxit, "maxit", 1);
  const double q = p + 1.0;
  double *w = rl::dalloc(2 * q * q + 4 * q);
  int *iw = rl::ialloc(n);
  SEXP cf = PROTECT(Rf_allocVector(REALSXP, p + 1));
  SEXP se = PROTECT(Rf_allocVector(REALSXP, p + 1));
  double dev;
  const int rc = clab_lgtreg(REAL_RO(xv), n, p, REAL_RO(yv), th, mi, REAL(cf),
                             REAL(se), &dev, w, iw);
  if (rc == CLAB_EDOM) {
    fail("no fit: collinear columns, or separated outcomes (no convergence "
         "within `maxit` steps).");
  }
  kchk(rc, nullptr);
  // names: "(Intercept)" and the column names of x, or x1, x2, ...
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, p + 1));
  SET_STRING_ELT(nms, 0, Rf_mkChar("(Intercept)"));
  SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
  SEXP cn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
  char b[32];
  for (int j = 0; j < p; ++j) {
    if (!Rf_isNull(cn)) {
      SET_STRING_ELT(nms, j + 1, STRING_ELT(cn, j));
    } else {
      std::snprintf(b, sizeof b, "x%d", j + 1);
      SET_STRING_ELT(nms, j + 1, Rf_mkChar(b));
    }
  }
  Rf_setAttrib(cf, R_NamesSymbol, nms);
  Rf_setAttrib(se, R_NamesSymbol, nms);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, cf);
  SET_VECTOR_ELT(out, 1, se);
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(dev));
  static const char *const nm[] = {"coef", "se", "dev"};
  UNPROTECT(6);
  return clab::nlist(out, nm);
}

}  // extern "C"

extern "C" {

SEXP C_clab_r_rfor(SEXP x, SEXP y, SEXP newx, SEXP ntree, SEXP mtry,
                   SEXP minnode, SEXP frac, SEXP repl, SEXP hon, SEXP imb,
                   SEXP seed) {
  clab::who = "rfor";
  int n, p, nt, pt;
  SEXP xv = PROTECT(rl::mat(x, "x", &n, &p));
  if (n < 2) fail("`x` must have at least two rows.");
  SEXP yv = PROTECT(num(y, "y", n));
  SEXP tv = PROTECT(rl::mat(newx, "newx", &nt, &pt));
  if (pt != p) fail("`newx` must have the %d columns of `x`.", p);
  const int ntr = clab::whole(ntree, "ntree", 1);
  // grf's rule; see the help page for the measurements behind it
  const int mt = Rf_isNull(mtry) ? static_cast<int>(std::fmin(std::ceil(std::sqrt(p) + 20.0), p))
                                 : clab::whole(mtry, "mtry", 1);
  if (mt > p) fail("`mtry` must be at most the %d columns of `x`.", p);
  const int mn = clab::whole(minnode, "minnode", 1);
  const double fr = REAL_RO(num(frac, "frac", 1, 0.0, true))[0];
  if (fr > 1.0) fail("`frac` must be at most 1.");
  clab::lgl(repl, "repl", 1);
  clab::lgl(hon, "hon", 1);
  const double ib = REAL_RO(num(imb, "imb", 1, 0.0))[0];
  if (ib >= 0.5) fail("`imb` must be below 0.5.");
  const double sd = REAL_RO(num(seed, "seed", 1, 0.0))[0];
  if (std::lround(fr * n) < (LOGICAL_RO(hon)[0] ? 2 : 1)) {
    fail("`frac` leaves too few rows in a bag.");
  }
  SEXP pr = PROTECT(Rf_allocVector(REALSXP, nt));
  SEXP ob = PROTECT(Rf_allocVector(REALSXP, n));
  double *w = rl::dalloc(static_cast<double>(lrn::fwork(n, p, nt)));
  int *iw = rl::ialloc(static_cast<double>(lrn::fiwork(n, p, nt)) + n + nt);
  kchk(clab_rfor(REAL_RO(xv), n, p, REAL_RO(yv), REAL_RO(tv), nt, ntr, mt, mn,
                 fr, LOGICAL_RO(repl)[0], LOGICAL_RO(hon)[0], ib,
                 static_cast<unsigned int>(std::fmod(sd, 4294967296.0)), REAL(pr),
                 REAL(ob), w, iw),
       nullptr);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, pr);
  SET_VECTOR_ELT(out, 1, ob);
  static const char *const nm[] = {"pred", "oob"};
  UNPROTECT(6);
  return clab::nlist(out, nm);
}

}  // extern "C"

namespace {

const char *const mdl[] = {"plr", "pliv", "irm", "iivm", "aipw",
                           "tmle", "did", "didcs", "ssm"};
// the nuisance names, model by model (the IV-type score adds "g")
const char *const nn[][5] = {
    {"l", "m", "g"},       {"l", "m", "r", "g"}, {"g0", "g1", "m"},
    {"g0", "g1", "m", "r0", "r1"}, {"g0", "g1", "m"}, {"g0", "g1", "m"},
    {"g0", "m"}, {"g00", "g01", "g10", "g11", "m"}, {"pi", "m", "g0", "g1"}};

const char *const lch[] = {"glm", "lasso", "ridge", "enet", "forest", "stack"};

// The settings shared by every learner of a call.
struct Lset {
  double alpha;
  int nlam, nfin;
  double tol;
  int scr, hon;
  unsigned int seed;
  int ntree, mtry;
};

// Learner k of lch (not the stack).
clab_lrn base(int k, const Lset &s) {
  clab_lrn r = {k == 0 ? 0 : (k == 4 ? 2 : 1),
                k == 1 ? 1.0 : (k == 2 ? 0.0 : s.alpha),
                s.nlam, s.nfin, s.seed, s.ntree, s.mtry, 0, s.tol, s.scr, s.hon,
                nullptr, 0, 0, 0};
  return r;
}

// The learners named by stk, each once, into sub (5 at most); their
// number.
int subs(SEXP stk, const char *nm, const Lset &s, clab_lrn *sub) {
  const R_xlen_t m = XLENGTH(stk);
  if (TYPEOF(stk) != STRSXP || m < 1 || m > 5) {
    fail("`%s` must name from one to five learners.", nm);
  }
  int seen[5] = {0, 0, 0, 0, 0};
  for (R_xlen_t j = 0; j < m; ++j) {
    SEXP e = PROTECT(Rf_ScalarString(STRING_ELT(stk, j)));
    const int k = rl::pick(e, nm, lch, 5);
    UNPROTECT(1);
    if (seen[k]++ > 0) fail("`%s` names \"%s\" twice.", nm, lch[k]);
    sub[j] = base(k, s);
  }
  return static_cast<int>(m);
}

// The learner named by l; a stack of the learners in stk (into sub) on
// nstk folds, disc as clab_stack().
clab_lrn learner(SEXP l, const char *nm, const Lset &s, SEXP stk, int nstk,
                 int disc, clab_lrn *sub) {
  const int k = rl::pick(l, nm, lch, 6);
  if (k < 5) return base(k, s);
  const int ns = subs(stk, "stk", s, sub);
  clab_lrn r = {3, s.alpha, s.nlam, s.nfin, s.seed, s.ntree, s.mtry, 0, s.tol, s.scr,
                s.hon, sub, ns, nstk, disc};
  return r;
}

// a 0/1 vector of length n
SEXP bin(SEXP v, const char *nm, int n) {
  SEXP a = PROTECT(num(v, nm, n));
  rl::ychk(REAL_RO(a), n, true, nm);
  UNPROTECT(1);
  return a;
}

}  // namespace

extern "C" {

SEXP C_clab_r_dmlx(SEXP y, SEXP d, SEXP x, SEXP model, SEXP fold, SEXP z,
                   SEXP s, SEXP t, SEXP score, SEXP lrn, SEXP lrnd,
                   SEXP alpha, SEXP nlam, SEXP nfin, SEXP tol, SEXP scr,
                   SEXP hon, SEXP ntree, SEXP mtry,
                   SEXP trim, SEXP clu, SEXP pred, SEXP stk, SEXP nstk,
                   SEXP disc, SEXP seed) {
  clab::who = "dmlx";
  const int md = rl::pick(model, "model", mdl, 9);
  int n, p;
  SEXP xv = PROTECT(rl::mat(x, "x", &n, &p));
  if (n < 4 || p < 1) fail("`x` must have at least four rows and one column.");
  // s: the selection indicator (ssm) or the period (didcs)
  SEXP sv = PROTECT(md == 8 ? bin(s, "s", n) : (md == 7 ? bin(t, "t", n) : R_NilValue));
  SEXP yv;
  if (md == 8) {
    // the outcome may be missing where it is not observed; never read there
    if (!clab::isnum(y) || XLENGTH(y) != n) {
      fail("`y` must be numbers, one per row of `x`.");
    }
    yv = PROTECT(clab::dcopy(clab::dbl(y)));
    double *a = REAL(yv);
    for (int i = 0; i < n; ++i) {
      if (REAL_RO(sv)[i] == 0.0) {
        a[i] = 0.0;
      } else if (!std::isfinite(a[i])) {
        fail("`y` must be finite where `s` is 1.");
      }
    }
  } else {
    yv = PROTECT(num(y, "y", n));
  }
  SEXP dv = PROTECT(md >= 2 ? bin(d, "d", n) : num(d, "d", n));
  SEXP zv = PROTECT(md == 1 ? num(z, "z", n) : (md == 3 ? bin(z, "z", n) : R_NilValue));
  int sc = 0;
  if (!Rf_isNull(score)) {
    static const char *const s01[] = {"po", "iv"}, *const s23[] = {"ate", "atte"},
                             *const s67[] = {"obs", "exp"};
    if (md == 0 || md == 1) {
      sc = rl::pick(score, "score", s01, 2);
    } else if (md == 2) {
      sc = rl::pick(score, "score", s23, 2);
    } else if (md == 6 || md == 7) {
      sc = rl::pick(score, "score", s67, 2);
    } else {
      fail("`score` must be NULL for model \"%s\".", mdl[md]);
    }
  }
  // folds: a vector, or a matrix with one column per repetition
  SEXP fd = Rf_getAttrib(fold, R_DimSymbol);
  const int nrep = TYPEOF(fd) == INTSXP && XLENGTH(fd) == 2 ? INTEGER_RO(fd)[1] : 1;
  if (nrep < 1 || XLENGTH(fold) != static_cast<R_xlen_t>(n) * nrep) {
    fail("`fold` must have one row per row of `x`.");
  }
  if (!clab::isnum(fold)) fail("`fold` must be finite numbers with no missing values.");
  SEXP fdb = PROTECT(clab::dbl(fold));
  SEXP col = PROTECT(Rf_allocVector(REALSXP, n));
  int nf = 0;
  int *fv = rl::ialloc(static_cast<double>(n) * nrep);
  for (int r = 0; r < nrep; ++r) {
    std::memcpy(REAL(col), REAL_RO(fdb) + static_cast<R_xlen_t>(r) * n, sizeof(double) * n);
    int k;
    const int *f = rl::folds(col, "fold", n, &k);
    if (r > 0 && k != nf) fail("`fold` must number the same folds in every column.");
    nf = k;
    std::memcpy(fv + static_cast<R_xlen_t>(r) * n, f, sizeof(int) * n);
  }
  UNPROTECT(2);
  int nc = 0, nc2 = 0;
  const int *cv = nullptr;
  int nrk = nrep;  // the repetitions the kernel sees
  if (!Rf_isNull(clu)) {
    if (md > 1) fail("`clu` is for models \"plr\" and \"pliv\" only.");
    SEXP cd = Rf_getAttrib(clu, R_DimSymbol);
    const bool two = TYPEOF(cd) == INTSXP && XLENGTH(cd) == 2 && INTEGER_RO(cd)[1] == 2;
    int *cw;
    if (two) {
      // two-way: a column per kind of cluster, and fold's two columns are
      // the folds of each kind
      if (INTEGER_RO(cd)[0] != n || nrep != 2) {
        fail("two-way `clu` needs one row per row of `x` and `fold` with "
             "two columns, the folds of each kind of cluster.");
      }
      SEXP cdb = PROTECT(clab::dbl(clu));
      SEXP c1 = PROTECT(Rf_allocVector(REALSXP, n));
      std::memcpy(REAL(c1), REAL_RO(cdb), sizeof(double) * n);
      const int *a = rl::folds(c1, "clu", n, &nc);
      std::memcpy(REAL(c1), REAL_RO(cdb) + n, sizeof(double) * n);
      const int *b = rl::folds(c1, "clu", n, &nc2);
      cw = rl::ialloc(2.0 * n);
      std::memcpy(cw, a, sizeof(int) * n);
      std::memcpy(cw + n, b, sizeof(int) * n);
      UNPROTECT(2);
      nrk = 1;
    } else {
      cw = rl::folds(clu, "clu", n, &nc);
    }
    cv = cw;
    // a fold per cluster: in every repetition, or of each kind (two-way)
    int *cf = rl::ialloc(nc > nc2 ? nc : nc2);
    for (int r = 0; r < nrep; ++r) {
      const int *cr = cw + (two ? static_cast<R_xlen_t>(r) * n : 0);
      const int ncr = two && r == 1 ? nc2 : nc;
      for (int c = 0; c < ncr; ++c) cf[c] = -1;
      for (int i = 0; i < n; ++i) {
        const int f = fv[static_cast<R_xlen_t>(r) * n + i];
        if (cf[cr[i]] >= 0 && cf[cr[i]] != f) {
          fail("`fold` must be constant within each cluster.");
        }
        cf[cr[i]] = f;
      }
    }
  }
  const double al = clab::p01(alpha, "alpha");
  const int nl = clab::whole(nlam, "nlam", 1);
  const int ni = clab::whole(nfin, "nfin", 2);
  const double sd = REAL_RO(num(seed, "seed", 1, 0.0))[0];
  const unsigned int us = static_cast<unsigned int>(std::fmod(sd, 4294967296.0));
  const int nt = clab::whole(ntree, "ntree", 1);
  const int mt = Rf_isNull(mtry) ? 0 : clab::whole(mtry, "mtry", 1);
  const double tl = REAL_RO(num(tol, "tol", 1, 0.0, true))[0];
  clab::lgl(scr, "scr", 1);
  clab::lgl(hon, "hon", 1);
  const int sr = LOGICAL_RO(scr)[0], hn = LOGICAL_RO(hon)[0];
  const int ks = clab::whole(nstk, "nstk", 2);
  clab::lgl(disc, "disc", 1);
  const int dc = LOGICAL_RO(disc)[0];
  clab_lrn sy[5], sdd[5];
  const Lset ls = {al, nl, ni, tl, sr, hn, us, nt, mt};
  Lset lsd = ls;
  lsd.seed = us + 1U;
  const clab_lrn ly = learner(lrn, "lrn", ls, stk, ks, dc, sy);
  const clab_lrn ld = learner(lrnd, "lrnd", lsd, stk, ks, dc, sdd);
  double tr;
  if (Rf_isNull(trim)) {
    tr = md == 4 ? 0.025 : (md == 5 ? 5.0 / std::sqrt(n) / std::log(n)
                                    : (md == 6 || md == 7 ? 0.01 : 1e-12));
  } else {
    tr = REAL_RO(num(trim, "trim", 1, 0.0))[0];
    if (tr >= 0.5) fail("`trim` must be below 0.5.");
  }
  const int k = clab_dmlxk(md, sc);
  SEXP pr;
  int fit = 1;
  if (Rf_isNull(pred)) {
    pr = PROTECT(Rf_alloc3DArray(REALSXP, n, k, nrk));
  } else {
    SEXP pd = Rf_getAttrib(pred, R_DimSymbol);
    const bool ok = TYPEOF(pd) == INTSXP && INTEGER_RO(pd)[0] == n &&
                    INTEGER_RO(pd)[1] == k &&
                    (XLENGTH(pd) == 2 ? nrk == 1 : XLENGTH(pd) == 3 && INTEGER_RO(pd)[2] == nrk);
    if (!clab::isnum(pred) || !ok) {
      fail("`pred` must be an array of %d rows, %d columns and one slice "
           "per column of `fold`.", n, k);
    }
    SEXP pv = PROTECT(clab::dcopy(clab::dbl(pred)));
    for (R_xlen_t i = 0; i < XLENGTH(pv); ++i) {
      if (!std::isfinite(REAL_RO(pv)[i])) fail("`pred` must be finite.");
    }
    pr = PROTECT(Rf_alloc3DArray(REALSXP, n, k, nrk));
    std::memcpy(REAL(pr), REAL_RO(pv), sizeof(double) * XLENGTH(pv));
    UNPROTECT(2);
    PROTECT(pr);
    fit = 0;
  }
  long niw;
  const long nw = clab_dmlxsz(md, n, p, nf, nrk, &ly, &ld, &niw);
  double *w = rl::dalloc(static_cast<double>(nw));
  int *iw = rl::ialloc(static_cast<double>(niw));
  SEXP est = PROTECT(Rf_allocVector(REALSXP, 2 + 2 * nrk));
  SEXP psia = PROTECT(Rf_allocMatrix(REALSXP, n, nrk));
  SEXP psib = PROTECT(Rf_allocMatrix(REALSXP, n, nrk));
  double *ps = rl::dalloc(2.0 * n * nrk);
  const int rc = clab_dmlx(md, sc, REAL_RO(xv), n, p, REAL_RO(yv), REAL_RO(dv),
                          Rf_isNull(zv) ? nullptr : REAL_RO(zv),
                          Rf_isNull(sv) ? nullptr : REAL_RO(sv), fv, nf, nrk,
                          cv, nc, nc2, &ly, &ld, tr, fit, REAL(pr), REAL(est), ps, w, iw);
  if (rc == CLAB_EDOM) {
    fail("no estimate: a learner failed, a training set was empty, or the "
         "score has no root.");
  }
  kchk(rc, nullptr);
  std::memcpy(REAL(psia), ps, sizeof(double) * n * nrk);
  std::memcpy(REAL(psib), ps + static_cast<R_xlen_t>(n) * nrk, sizeof(double) * n * nrk);
  // dimnames of pred: none, the nuisance names, none
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP cn = PROTECT(Rf_allocVector(STRSXP, k));
  // the experimental panel score learns g1 in place of m
  for (int j = 0; j < k; ++j) {
    SET_STRING_ELT(cn, j, Rf_mkChar(md == 6 && sc == 1 && j == 1 ? "g1" : nn[md][j]));
  }
  SET_VECTOR_ELT(dn, 1, cn);
  Rf_setAttrib(pr, R_DimNamesSymbol, dn);
  SEXP reps = PROTECT(Rf_allocVector(VECSXP, 2));
  for (int c = 0; c < 2; ++c) {
    SEXP v = Rf_allocVector(REALSXP, nrk);
    SET_VECTOR_ELT(reps, c, v);
    std::memcpy(REAL(v), REAL_RO(est) + 2 + c * nrk, sizeof(double) * nrk);
  }
  static const char *const rnm[] = {"coef", "se"};
  clab::frame(reps, rnm);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(REAL_RO(est)[0]));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(REAL_RO(est)[1]));
  SET_VECTOR_ELT(out, 2, reps);
  SET_VECTOR_ELT(out, 3, psia);
  SET_VECTOR_ELT(out, 4, psib);
  SET_VECTOR_ELT(out, 5, pr);
  static const char *const onm[] = {"coef", "se", "reps", "psia", "psib", "pred"};
  UNPROTECT(13);
  return clab::nlist(out, onm);
}

SEXP C_clab_r_suplrn(SEXP x, SEXP y, SEXP fam, SEXP lrns, SEXP nfold, SEXP fold,
                     SEXP disc, SEXP newx, SEXP alpha, SEXP nlam, SEXP nfin,
                     SEXP ntree, SEXP mtry, SEXP seed) {
  clab::who = "suplrn";
  int n, p;
  SEXP xv = PROTECT(rl::mat(x, "x", &n, &p));
  if (n < 4 || p < 1) fail("`x` must have at least four rows and one column.");
  static const char *const fams[] = {"gaussian", "binomial"};
  const int fm = rl::pick(fam, "fam", fams, 2);
  SEXP yv = PROTECT(num(y, "y", n));
  rl::ychk(REAL_RO(yv), n, fm == 1, "y");
  int nt = n;
  SEXP tv = xv;
  if (!Rf_isNull(newx)) {
    int pt;
    tv = rl::mat(newx, "newx", &nt, &pt);
    if (pt != p) fail("`newx` must have the %d columns of `x`.", p);
  }
  PROTECT(tv);
  const double sd = REAL_RO(num(seed, "seed", 1, 0.0))[0];
  const unsigned int us = static_cast<unsigned int>(std::fmod(sd, 4294967296.0));
  const Lset ls = {clab::p01(alpha, "alpha"), clab::whole(nlam, "nlam", 1),
                   clab::whole(nfin, "nfin", 2), 1e-7, 1, 0, us,
                   clab::whole(ntree, "ntree", 1),
                   Rf_isNull(mtry) ? 0 : clab::whole(mtry, "mtry", 1)};
  clab::lgl(disc, "disc", 1);
  clab_lrn sub[5];
  const int ns = subs(lrns, "lrns", ls, sub);
  int k;
  int *fv;
  if (Rf_isNull(fold)) {
    k = clab::whole(nfold, "nfold", 2);
    if (k > n) fail("`nfold` must be at most the %d rows of `x`.", n);
    fv = rl::ialloc(n);
    kchk(clab_kfold(n, k, us, fv, rl::ialloc(n)), nullptr);
  } else {
    fv = rl::folds(fold, "fold", n, &k);
  }
  const clab_lrn l = {3, ls.alpha, ls.nlam, ls.nfin, us, ls.ntree, ls.mtry, 0, ls.tol,
                      ls.scr, ls.hon, sub, ns, k, LOGICAL_RO(disc)[0]};
  long niw;
  const long nw = clab_stksz(n, p, nt, &l, &niw);
  double *w = rl::dalloc(static_cast<double>(nw));
  int *iw = rl::ialloc(static_cast<double>(niw));
  SEXP oof = PROTECT(Rf_allocMatrix(REALSXP, n, ns));
  SEXP risk = PROTECT(Rf_allocVector(REALSXP, ns));
  SEXP wt = PROTECT(Rf_allocVector(REALSXP, ns));
  SEXP pr = PROTECT(Rf_allocVector(REALSXP, nt));
  const int rc = clab_stack(REAL_RO(xv), n, p, REAL_RO(yv), fm, &l, fv, REAL_RO(tv), nt,
                            REAL(oof), REAL(risk), REAL(wt), REAL(pr), w, iw);
  if (rc == CLAB_EDOM) {
    fail("a learner failed on a fold (a fold too small for its inner folds, "
         "or a design it cannot fit).");
  }
  kchk(rc, nullptr);
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, ns));
  for (int j = 0; j < ns; ++j) SET_STRING_ELT(nm, j, STRING_ELT(lrns, j));
  Rf_setAttrib(risk, R_NamesSymbol, nm);
  Rf_setAttrib(wt, R_NamesSymbol, nm);
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn, 1, nm);
  Rf_setAttrib(oof, R_DimNamesSymbol, dn);
  SEXP fo = PROTECT(Rf_allocVector(INTSXP, n));
  for (int i = 0; i < n; ++i) INTEGER(fo)[i] = fv[i] + 1;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, wt);
  SET_VECTOR_ELT(out, 1, risk);
  SET_VECTOR_ELT(out, 2, fo);
  SET_VECTOR_ELT(out, 3, oof);
  SET_VECTOR_ELT(out, 4, pr);
  static const char *const onm[] = {"wt", "risk", "fold", "oof", "pred"};
  UNPROTECT(11);
  return clab::nlist(out, onm);
}

}  // extern "C"
