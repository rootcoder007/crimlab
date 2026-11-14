// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of pedlr(), gmixlr() and cmpsc(): the shapes of the
// arguments are checked here, in the words clab_jrun() uses too; the rest
// of the checks and the arithmetic are frx:: (frx.h) and the kernels.
// Every C++ object is gone before an error is raised: a refusal from frx
// is copied to a buffer and raised once its scope has closed.
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <initializer_list>
#include <string>
#include <vector>

#include "crimlab.h"
#include "frx.h"
#include "rchk.h"

namespace {

char msg[1024];

void say(const std::string &s) { std::snprintf(msg, sizeof msg, "%s", s.c_str()); }

// the column nm of the data frame df, or R_NilValue
SEXP dcol(SEXP df, const char *nm) {
  SEXP n = Rf_getAttrib(df, R_NamesSymbol);
  for (R_xlen_t i = 0; i < Rf_xlength(n); ++i) {
    if (std::strcmp(Rf_translateCharUTF8(STRING_ELT(n, i)), nm) == 0) return VECTOR_ELT(df, i);
  }
  return R_NilValue;
}

bool isdf(SEXP x, std::initializer_list<const char *> cols) {
  if (!Rf_inherits(x, "data.frame")) return false;
  for (const char *c : cols) {
    if (Rf_isNull(dcol(x, c))) return false;
  }
  return true;
}

// as.character() of a column; false (and msg) for dates and times
bool lab(SEXP x, const char *an, frx::Lab &l) {
  SEXP s = PROTECT(clab::labs(x));
  if (Rf_isNull(s) || TYPEOF(s) != STRSXP) {
    UNPROTECT(1);
    say(std::string("`") + an + "` must hold labels: text or numbers.");
    return false;
  }
  const R_xlen_t n = XLENGTH(s);
  l.s.resize(static_cast<std::size_t>(n));
  l.na.resize(static_cast<std::size_t>(n));
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP e = STRING_ELT(s, i);
    l.na[static_cast<std::size_t>(i)] = e == NA_STRING;
    if (e != NA_STRING) l.s[static_cast<std::size_t>(i)] = Rf_translateCharUTF8(e);
  }
  UNPROTECT(1);
  return true;
}

SEXP strs(const std::vector<std::string> &v) {
  SEXP s = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(v.size())));
  for (std::size_t i = 0; i < v.size(); ++i) {
    SET_STRING_ELT(s, static_cast<R_xlen_t>(i), Rf_mkCharCE(v[i].c_str(), CE_UTF8));
  }
  UNPROTECT(1);
  return s;
}

SEXP reals(const std::vector<double> &v) {
  SEXP s = Rf_allocVector(REALSXP, static_cast<R_xlen_t>(v.size()));
  for (std::size_t i = 0; i < v.size(); ++i) REAL(s)[i] = v[i];
  return s;
}

SEXP classed(SEXP x, const char *cl) {
  SEXP c = PROTECT(Rf_mkString(cl));
  Rf_classgets(x, c);
  UNPROTECT(1);
  return x;
}

// pedlr() once its arguments have their shapes: the result, or R_NilValue
// with msg set; code receives the error code
SEXP pedrun(SEXP ped, SEXP geno, SEXP freq, SEXP alt, SEXP mut, SEXP rate, SEXP range, int *code) {
  frx::PedIn in;
  const bool ok = lab(dcol(ped, "id"), "ped", in.id) && lab(dcol(ped, "fid"), "ped", in.fid) &&
                  lab(dcol(ped, "mid"), "ped", in.mid) &&
                  (Rf_isNull(alt) || (lab(dcol(alt, "id"), "alt", in.aid) && lab(dcol(alt, "fid"), "alt", in.afid) &&
                                      lab(dcol(alt, "mid"), "alt", in.amid))) &&
                  lab(dcol(geno, "id"), "geno", in.gid) && lab(dcol(geno, "locus"), "geno", in.gloc) &&
                  lab(dcol(geno, "a1"), "geno", in.ga1) && lab(dcol(geno, "a2"), "geno", in.ga2) &&
                  lab(dcol(freq, "locus"), "freq", in.floc) && lab(dcol(freq, "allele"), "freq", in.fal);
  if (!ok) {
    *code = CLAB_EARG;
    return R_NilValue;
  }
  in.hasalt = !Rf_isNull(alt);
  SEXP fq = PROTECT(Rf_coerceVector(dcol(freq, "freq"), REALSXP));
  in.fq.assign(REAL(fq), REAL(fq) + XLENGTH(fq));
  SEXP rt = PROTECT(Rf_coerceVector(rate, REALSXP));
  in.rate.assign(REAL(rt), REAL(rt) + XLENGTH(rt));
  UNPROTECT(2);
  in.mut = Rf_translateCharUTF8(STRING_ELT(mut, 0));
  in.range = clab::val1(range);
  frx::PedOut o;
  std::string err;
  *code = frx::pedlr(in, o, err);
  if (*code != CLAB_OK) {
    say(err);
    return R_NilValue;
  }
  SEXP loc = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(loc, 0, strs(o.loc));
  SET_VECTOR_ELT(loc, 1, reals(o.l1));
  SET_VECTOR_ELT(loc, 2, reals(o.l2));
  SET_VECTOR_ELT(loc, 3, reals(o.lr));
  static const char *const ln[] = {"locus", "l1", "l2", "lr"};
  clab::frame(loc, ln);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, loc);
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(o.log10lr));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(std::pow(10.0, o.log10lr)));
  SET_VECTOR_ELT(out, 3, Rf_mkString(in.mut.c_str()));
  static const char *const on[] = {"locus", "log10lr", "lr", "mut"};
  clab::nlist(out, on);
  UNPROTECT(2);
  return classed(out, "pedlr");
}

// labels from an atomic vector (NULL: none); false for anything else
bool names(SEXP x, std::vector<std::string> &out) {
  out.clear();
  if (Rf_isNull(x)) return true;
  if (!Rf_isVectorAtomic(x) || !(TYPEOF(x) == STRSXP || clab::isnum(x))) return false;
  frx::Lab l;
  if (!lab(x, "hp", l)) return false;
  for (std::size_t i = 0; i < l.s.size(); ++i) {
    if (l.na[i]) return false;
    out.push_back(l.s[i]);
  }
  return true;
}

SEXP nas(const double *x, int n) {
  SEXP s = Rf_allocVector(REALSXP, n);
  for (int i = 0; i < n; ++i) REAL(s)[i] = std::isnan(x[i]) ? NA_REAL : x[i];
  return s;
}

// gmixlr() once its arguments have their shapes
SEXP gmxrun(SEXP peaks, SEXP freq, SEXP ref, SEXP hp, SEXP nc, SEXP hd, SEXP at, SEXP pc, SEXP lambda,
            SEXP theta, SEXP stut, SEXP deg, SEXP par, int *code) {
  frx::GmxIn in;
  *code = CLAB_EARG;
  if (!names(hp, in.hp)) {
    say("`hp` must be the ids of profiles in `ref`.");
    return R_NilValue;
  }
  if (!names(hd, in.hd)) {
    say("`hd` must be NULL or the ids of profiles in `ref`.");
    return R_NilValue;
  }
  const bool ok = lab(dcol(peaks, "locus"), "peaks", in.ploc) && lab(dcol(peaks, "allele"), "peaks", in.pal) &&
                  lab(dcol(freq, "locus"), "freq", in.floc) && lab(dcol(freq, "allele"), "freq", in.fal) &&
                  (Rf_isNull(ref) || (lab(dcol(ref, "id"), "ref", in.rid) && lab(dcol(ref, "locus"), "ref", in.rloc) &&
                                      lab(dcol(ref, "a1"), "ref", in.ra1) && lab(dcol(ref, "a2"), "ref", in.ra2)));
  if (!ok) return R_NilValue;
  SEXP h = PROTECT(Rf_coerceVector(dcol(peaks, "height"), REALSXP));
  in.ph.assign(REAL(h), REAL(h) + XLENGTH(h));
  SEXP f = PROTECT(Rf_coerceVector(dcol(freq, "freq"), REALSXP));
  in.fq.assign(REAL(f), REAL(f) + XLENGTH(f));
  in.hassz = !Rf_isNull(dcol(freq, "size"));
  if (in.hassz) {
    SEXP z = PROTECT(Rf_coerceVector(dcol(freq, "size"), REALSXP));
    in.fsz.assign(REAL(z), REAL(z) + XLENGTH(z));
    UNPROTECT(1);
  }
  if (!Rf_isNull(par)) {
    SEXP pv = PROTECT(Rf_coerceVector(par, REALSXP));
    SEXP pn = Rf_getAttrib(par, R_NamesSymbol);
    in.haspar = true;
    for (R_xlen_t i = 0; i < XLENGTH(pv); ++i) {
      in.par.push_back(REAL(pv)[i]);
      in.parnm.push_back(Rf_translateCharUTF8(STRING_ELT(pn, i)));
    }
    UNPROTECT(1);
  }
  UNPROTECT(2);
  in.nc = clab::val1(nc);
  in.at = clab::val1(at);
  in.pc = clab::val1(pc);
  in.lambda = clab::val1(lambda);
  in.theta = clab::val1(theta);
  in.stut = Rf_asLogical(stut) == TRUE;
  in.deg = Rf_asLogical(deg) == TRUE;
  frx::GmxOut o;
  std::string err;
  *code = frx::gmixlr(in, o, err);
  if (*code != CLAB_OK) {
    say(err);
    return R_NilValue;
  }
  SEXP fit = PROTECT(Rf_allocVector(VECSXP, 8));
  SET_VECTOR_ELT(fit, 0, strs({"hp", "hd"}));
  SET_VECTOR_ELT(fit, 1, nas(o.loglik, 2));
  SET_VECTOR_ELT(fit, 2, nas(o.mu, 2));
  SET_VECTOR_ELT(fit, 3, nas(o.omega, 2));
  SET_VECTOR_ELT(fit, 4, nas(o.beta, 2));
  SET_VECTOR_ELT(fit, 5, nas(o.xi, 2));
  SEXP cv = Rf_allocVector(LGLSXP, 2);
  SET_VECTOR_ELT(fit, 6, cv);
  SEXP it = Rf_allocVector(REALSXP, 2);
  SET_VECTOR_ELT(fit, 7, it);
  for (int k = 0; k < 2; ++k) {
    LOGICAL(cv)[k] = o.conv[k] < 0 ? NA_LOGICAL : o.conv[k] == 1;
    REAL(it)[k] = o.iter[k];
  }
  static const char *const fn[] = {"hyp", "loglik", "mu", "omega", "beta", "xi", "conv", "iter"};
  clab::frame(fit, fn);
  std::vector<std::string> hy, who;
  std::vector<double> m;
  for (int k = 0; k < 2; ++k) {
    for (std::size_t j = 0; j < o.who[k].size(); ++j) {
      hy.push_back(k ? "hd" : "hp");
      who.push_back(o.who[k][j]);
      m.push_back(o.m[k][j]);
    }
  }
  SEXP mx = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(mx, 0, strs(hy));
  SET_VECTOR_ELT(mx, 1, strs(who));
  SET_VECTOR_ELT(mx, 2, reals(m));
  static const char *const mn[] = {"hyp", "who", "m"};
  clab::frame(mx, mn);
  SEXP lc = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(lc, 0, strs(o.loc));
  SET_VECTOR_ELT(lc, 1, reals(o.llr));
  static const char *const lnm[] = {"locus", "log10lr"};
  clab::frame(lc, lnm);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(std::pow(10.0, o.log10lr)));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(o.log10lr));
  SET_VECTOR_ELT(out, 2, fit);
  SET_VECTOR_ELT(out, 3, mx);
  SET_VECTOR_ELT(out, 4, lc);
  static const char *const on[] = {"lr", "log10lr", "fit", "mix", "locus"};
  clab::nlist(out, on);
  UNPROTECT(4);
  return classed(out, "gmixlr");
}

}  // namespace

extern "C" {

//' @srrstats {G2.3} The character argument mut (pedlr()) is checked
//'   against none, equal and step, and refused otherwise.
//' @srrstats {G2.3b} mut is matched exactly, case included.
//' @noRd
SEXP C_clab_r_pedlr(SEXP ped, SEXP geno, SEXP freq, SEXP alt, SEXP mut, SEXP rate, SEXP range) {
  clab::who = "pedlr";
  if (!isdf(ped, {"id", "fid", "mid"})) clab::fail("`ped` must be a data frame with columns id, fid and mid.");
  if (!isdf(geno, {"id", "locus", "a1", "a2"})) {
    clab::fail("`geno` must be a data frame with columns id, locus, a1 and a2.");
  }
  if (!isdf(freq, {"locus", "allele", "freq"})) {
    clab::fail("`freq` must be a data frame with columns locus, allele and freq.");
  }
  if (!clab::isnum(dcol(freq, "freq"))) clab::fail("`freq$freq` must be numbers.");
  if (!Rf_isNull(alt) && !isdf(alt, {"id", "fid", "mid"})) {
    clab::fail("`alt` must be NULL or a data frame with columns id, fid and mid.");
  }
  if (TYPEOF(mut) != STRSXP || XLENGTH(mut) != 1 || STRING_ELT(mut, 0) == NA_STRING) {
    clab::fail("`mut` must be one of none, equal, step.");
  }
  if (!clab::isnum(rate)) clab::fail("`rate` must be one or two numbers from 0 to 1.");
  msg[0] = '\0';
  int code = CLAB_OK;
  SEXP out = pedrun(ped, geno, freq, alt, mut, rate, range, &code);
  clab::ichk(code);
  if (msg[0] != '\0') clab::fail("%s", msg);
  return out;
}

SEXP C_clab_r_gmixlr(SEXP peaks, SEXP freq, SEXP ref, SEXP hp, SEXP nc, SEXP hd, SEXP at, SEXP pc,
                     SEXP lambda, SEXP theta, SEXP stut, SEXP deg, SEXP par) {
  clab::who = "gmixlr";
  if (!isdf(peaks, {"locus", "allele", "height"})) {
    clab::fail("`peaks` must be a data frame with columns locus, allele and height.");
  }
  if (!clab::isnum(dcol(peaks, "height"))) clab::fail("`peaks$height` must be numbers.");
  if (!isdf(freq, {"locus", "allele", "freq"})) {
    clab::fail("`freq` must be a data frame with columns locus, allele and freq.");
  }
  if (!clab::isnum(dcol(freq, "freq"))) clab::fail("`freq$freq` must be numbers.");
  if (!Rf_isNull(dcol(freq, "size")) && !clab::isnum(dcol(freq, "size"))) clab::fail("`freq$size` must be numbers.");
  if (!Rf_isNull(ref) && !isdf(ref, {"id", "locus", "a1", "a2"})) {
    clab::fail("`ref` must be NULL or a data frame with columns id, locus, a1 and a2.");
  }
  for (SEXP b : {stut, deg}) {
    if (TYPEOF(b) != LGLSXP || XLENGTH(b) != 1 || LOGICAL(b)[0] == NA_LOGICAL) {
      clab::fail("`%s` must be TRUE or FALSE.", b == stut ? "stut" : "deg");
    }
  }
  if (!Rf_isNull(par) && (!clab::isnum(par) || Rf_isNull(Rf_getAttrib(par, R_NamesSymbol)))) {
    clab::fail("`par` must be NULL or named numbers.");
  }
  msg[0] = '\0';
  int code = CLAB_OK;
  SEXP out = gmxrun(peaks, freq, ref, hp, nc, hd, at, pc, lambda, theta, stut, deg, par, &code);
  clab::ichk(code);
  if (msg[0] != '\0') clab::fail("%s", msg);
  return out;
}

SEXP C_clab_r_cmpsc(SEXP a, SEXP b, SEXP len, SEXP npeaks, SEXP tx, SEXP full) {
  clab::who = "cmpsc";
  for (SEXP x : {a, b}) {
    const char *an = x == a ? "a" : "b";
    clab::unlab(x, an);
    bool lna = TYPEOF(x) == LGLSXP;
    for (R_xlen_t i = 0; lna && i < XLENGTH(x); ++i) lna = LOGICAL(x)[i] == NA_LOGICAL;
    if (!(clab::isnum(x) || lna) || XLENGTH(x) < 1) {
      clab::fail("`%s` must be numbers (NA for a missing point), not infinite.", an);
    }
    clab::vec(x, an);
  }
  SEXP av = PROTECT(Rf_coerceVector(a, REALSXP)), bv = PROTECT(Rf_coerceVector(b, REALSXP));
  for (SEXP x : {av, bv}) {
    for (R_xlen_t i = 0; i < XLENGTH(x); ++i) {
      if (std::isinf(REAL(x)[i])) clab::fail("`%s` must be numbers (NA for a missing point), not infinite.", x == av ? "a" : "b");
    }
  }
  const int na = clab::ilen(static_cast<double>(XLENGTH(av)), "a"), nb = clab::ilen(static_cast<double>(XLENGTH(bv)), "b");
  const int ln = clab::whole(len, "len", 2);
  int f0 = 0, f1 = na - 1;
  while (f0 < na && ISNAN(REAL(av)[f0])) ++f0;
  while (f1 >= f0 && ISNAN(REAL(av)[f1])) --f1;
  if (f1 - f0 + 1 < ln) clab::fail("`a` must hold at least `len` points from its first value to its last.");
  if (!clab::isnum(npeaks) || XLENGTH(npeaks) < 1 || XLENGTH(npeaks) > 8) {
    clab::fail("`npeaks` must be one to eight whole numbers from 1.");
  }
  SEXP pv = PROTECT(Rf_coerceVector(npeaks, REALSXP));
  int pk[8];
  for (R_xlen_t i = 0; i < XLENGTH(pv); ++i) {
    const double v = REAL(pv)[i];
    if (!(v >= 1 && v <= INT_MAX && v == std::floor(v))) clab::fail("`npeaks` must be one to eight whole numbers from 1.");
    pk[i] = static_cast<int>(v);
  }
  const int t = clab::whole(tx, "tx", 0);
  if (TYPEOF(full) != LGLSXP || XLENGTH(full) != 1 || LOGICAL(full)[0] == NA_LOGICAL) {
    clab::fail("`full` must be TRUE or FALSE.");
  }
  const int ns = (f1 - f0 + 1) / ln;
  SEXP cc = PROTECT(Rf_allocVector(REALSXP, ns)), cm = PROTECT(Rf_allocVector(LGLSXP, ns));
  int off = 0, nseg = 0, sc = 0;
  double pos = 0;
  const int rc = clab_cmpsc(REAL(av), na, REAL(bv), nb, ln, pk, static_cast<int>(XLENGTH(pv)), t, LOGICAL(full)[0],
                            &off, &nseg, &sc, &pos, REAL(cc), LOGICAL(cm));
  clab::ichk(rc);
  if (rc == CLAB_EARG) clab::fail("`b` is shorter than a widened segment; give a longer reference or `full = TRUE`.");
  clab::kchk(rc, "cmpsc(a, b, len, npeaks, tx, full)");
  SEXP sg = PROTECT(Rf_allocVector(INTSXP, ns)), st = PROTECT(Rf_allocVector(INTSXP, ns));
  for (int s = 0; s < ns; ++s) {
    INTEGER(sg)[s] = s + 1;
    INTEGER(st)[s] = off + s * ln + 1;
    if (ISNAN(REAL(cc)[s])) REAL(cc)[s] = NA_REAL;
  }
  SEXP segs = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(segs, 0, sg);
  SET_VECTOR_ELT(segs, 1, st);
  SET_VECTOR_ELT(segs, 2, cc);
  SET_VECTOR_ELT(segs, 3, cm);
  static const char *const sn[] = {"seg", "start", "ccp", "cmps"};
  clab::frame(segs, sn);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, Rf_ScalarInteger(sc));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(nseg));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(nseg > 0 ? static_cast<double>(sc) / nseg : NA_REAL));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(ISNAN(pos) ? NA_REAL : pos));
  SET_VECTOR_ELT(out, 4, segs);
  static const char *const on[] = {"score", "nseg", "scaled", "pos", "segs"};
  clab::nlist(out, on);
  UNPROTECT(9);
  return classed(out, "cmpsc");
}

}  // extern "C"
