// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// pedlr() for clab_jrun(): the shapes of the
// arguments are checked here in the words of the R binding
// (src/bfrx.cpp); the rest is frx:: (src/frx.h) and the kernels. A data
// frame is an object of columns of one length.

namespace fj {

// the column nm of o, an object, as labels (as.character()); false when
// o has no such column of plain values
bool jlab(const Val &o, const char *nm, frx::Lab &l) {
  l = frx::Lab();
  if (o.t != Val::OBJ) return false;
  for (std::size_t i = 0; i < o.k.size(); ++i) {
    if (o.k[i] != nm) continue;
    const std::vector<Val> one1(1, o.a[i]);
    const std::vector<Val> &e = o.a[i].t == Val::ARR ? o.a[i].a : one1;
    for (const Val &x : e) {
      if (x.t == Val::ARR || x.t == Val::OBJ) return false;
      l.na.push_back(x.t == Val::NUL);
      l.s.push_back(x.t == Val::STR ? x.s : x.t == Val::NUM ? rfmt(x.n) : x.t == Val::BOOL ? (x.b ? "TRUE" : "FALSE") : "");
    }
    return true;
  }
  return false;
}

// an object with the columns cols, all of one length
bool isdf(const Val &o, std::initializer_list<const char *> cols) {
  if (o.t != Val::OBJ) return false;
  std::size_t n = 0;
  for (std::size_t i = 0; i < o.a.size(); ++i) {
    const std::size_t m = o.a[i].t == Val::ARR ? o.a[i].a.size() : 1;
    if (i > 0 && m != n) return false;
    n = m;
  }
  frx::Lab l;
  for (const char *c : cols) {
    if (!jlab(o, c, l)) return false;
  }
  return true;
}

const Val &jget(const Val &o, const char *nm) {
  for (std::size_t i = 0; i < o.k.size(); ++i) {
    if (o.k[i] == nm) return o.a[i];
  }
  return nul;
}

}  // namespace fj

bool pedlr(Ctx &c, Val &out) {
  const Val &ped = c["ped"], &geno = c["geno"], &freq = c["freq"], &alt = c["alt"];
  if (!fj::isdf(ped, {"id", "fid", "mid"})) return c.fail("`ped` must be a data frame with columns id, fid and mid.");
  if (!fj::isdf(geno, {"id", "locus", "a1", "a2"})) {
    return c.fail("`geno` must be a data frame with columns id, locus, a1 and a2.");
  }
  if (!fj::isdf(freq, {"locus", "allele", "freq"})) {
    return c.fail("`freq` must be a data frame with columns locus, allele and freq.");
  }
  Dv fq;
  if (readd(fj::jget(freq, "freq"), fq) != NUM) return c.fail("`freq$freq` must be numbers.");
  if (alt.t != Val::NUL && !fj::isdf(alt, {"id", "fid", "mid"})) {
    return c.fail("`alt` must be NULL or a data frame with columns id, fid and mid.");
  }
  frx::PedIn in;
  const Val &mut = c["mut"];
  if (mut.t != Val::NUL && mut.t != Val::STR) return c.fail("`mut` must be one of none, equal, step.");
  in.mut = mut.t == Val::STR ? mut.s : "none";
  Dv rt;
  if (c.has("rate")) {
    if (readd(c["rate"], rt) != NUM) return c.fail("`rate` must be one or two numbers from 0 to 1.");
    in.rate = rt.v;
  } else {
    in.rate.assign(1, 0.001);
  }
  in.range = c.has("range") ? val1(c["range"]) : 0.5;
  fj::jlab(ped, "id", in.id);
  fj::jlab(ped, "fid", in.fid);
  fj::jlab(ped, "mid", in.mid);
  in.hasalt = alt.t != Val::NUL;
  if (in.hasalt) {
    fj::jlab(alt, "id", in.aid);
    fj::jlab(alt, "fid", in.afid);
    fj::jlab(alt, "mid", in.amid);
  }
  fj::jlab(geno, "id", in.gid);
  fj::jlab(geno, "locus", in.gloc);
  fj::jlab(geno, "a1", in.ga1);
  fj::jlab(geno, "a2", in.ga2);
  fj::jlab(freq, "locus", in.floc);
  fj::jlab(freq, "allele", in.fal);
  in.fq = fq.v;
  frx::PedOut o;
  std::string err;
  const int rc = frx::pedlr(in, o, err);
  if (rc != CLAB_OK) return c.fail("%s", err.c_str());
  Val l;
  put(l, "locus", vstrs(o.loc));
  put(l, "l1", vnums(o.l1));
  put(l, "l2", vnums(o.l2));
  put(l, "lr", vnums(o.lr));
  put(out, "locus", l);
  put(out, "log10lr", vnum(o.log10lr));
  put(out, "lr", vnum(std::pow(10.0, o.log10lr)));
  put(out, "mut", vstr(in.mut));
  return true;
}
