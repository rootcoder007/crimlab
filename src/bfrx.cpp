// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of pedlr(): the shapes of the
// arguments are checked here, in the words clab_jrun() uses too; the rest
// of the checks and the arithmetic are frx:: (frx.h) and the kernels.
// Every C++ object is gone before an error is raised: a refusal from frx
// is copied to a buffer and raised once its scope has closed.
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

}  // extern "C"
