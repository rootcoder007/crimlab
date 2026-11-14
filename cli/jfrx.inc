// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// pedlr() and gmixlr() for clab_jrun(): the shapes of the
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

namespace fj {

// labels: one text or number, or an array of them (null or absent: none)
bool jnames(const Val &v, std::vector<std::string> &out) {
  out.clear();
  const std::vector<Val> one1(1, v);
  const std::vector<Val> &e = v.t == Val::ARR ? v.a : v.t == Val::NUL ? std::vector<Val>() : one1;
  for (const Val &x : e) {
    if (x.t == Val::STR) out.push_back(x.s);
    else if (x.t == Val::NUM) out.push_back(rfmt(x.n));
    else return false;
  }
  return true;
}

bool jflag(Ctx &c, const char *nm, bool &b) {
  const Val &v = c[nm];
  if (v.t == Val::NUL) {
    b = false;
    return true;
  }
  if (v.t != Val::BOOL) return c.fail("`%s` must be TRUE or FALSE.", nm);
  b = v.b;
  return true;
}

double jdef(Ctx &c, const char *nm, double d) { return c.has(nm) ? val1(c[nm]) : d; }

Val nav(const std::vector<double> &x) {
  Val a = varr();
  for (double v : x) a.a.push_back(vnum(std::isnan(v) ? NA : v));
  return a;
}

}  // namespace fj

bool gmixlr(Ctx &c, Val &out) {
  const Val &pk = c["peaks"], &freq = c["freq"], &ref = c["ref"];
  if (!fj::isdf(pk, {"locus", "allele", "height"})) {
    return c.fail("`peaks` must be a data frame with columns locus, allele and height.");
  }
  Dv ph, fq, fs;
  if (readd(fj::jget(pk, "height"), ph) != NUM) return c.fail("`peaks$height` must be numbers.");
  if (!fj::isdf(freq, {"locus", "allele", "freq"})) {
    return c.fail("`freq` must be a data frame with columns locus, allele and freq.");
  }
  if (readd(fj::jget(freq, "freq"), fq) != NUM) return c.fail("`freq$freq` must be numbers.");
  frx::GmxIn in;
  in.hassz = fj::jget(freq, "size").t != Val::NUL;
  if (in.hassz && readd(fj::jget(freq, "size"), fs) != NUM) return c.fail("`freq$size` must be numbers.");
  if (ref.t != Val::NUL && !fj::isdf(ref, {"id", "locus", "a1", "a2"})) {
    return c.fail("`ref` must be NULL or a data frame with columns id, locus, a1 and a2.");
  }
  CK(fj::jflag(c, "stut", in.stut));
  CK(fj::jflag(c, "deg", in.deg));
  const Val &pv = c["par"];
  if (pv.t != Val::NUL) {
    Dv d;
    if (pv.t != Val::OBJ || readd(pv, d) != NUM) return c.fail("`par` must be NULL or named numbers.");
    in.haspar = true;
    in.parnm = d.nm;
    in.par = d.v;
  }
  if (!fj::jnames(c["hp"], in.hp)) return c.fail("`hp` must be the ids of profiles in `ref`.");
  if (!fj::jnames(c["hd"], in.hd)) return c.fail("`hd` must be NULL or the ids of profiles in `ref`.");
  in.nc = val1(c["nc"]);
  in.at = fj::jdef(c, "at", 50);
  in.pc = fj::jdef(c, "pc", 0.05);
  in.lambda = fj::jdef(c, "lambda", 0.01);
  in.theta = fj::jdef(c, "theta", 0);
  fj::jlab(pk, "locus", in.ploc);
  fj::jlab(pk, "allele", in.pal);
  in.ph = ph.v;
  fj::jlab(freq, "locus", in.floc);
  fj::jlab(freq, "allele", in.fal);
  in.fq = fq.v;
  in.fsz = fs.v;
  if (ref.t != Val::NUL) {
    fj::jlab(ref, "id", in.rid);
    fj::jlab(ref, "locus", in.rloc);
    fj::jlab(ref, "a1", in.ra1);
    fj::jlab(ref, "a2", in.ra2);
  }
  frx::GmxOut o;
  std::string err;
  if (frx::gmixlr(in, o, err) != CLAB_OK) return c.fail("%s", err.c_str());
  put(out, "lr", vnum(std::pow(10.0, o.log10lr)));
  put(out, "log10lr", vnum(o.log10lr));
  Val f, mx, lc;
  put(f, "hyp", vstrs({"hp", "hd"}));
  put(f, "loglik", vnums(o.loglik, 2));
  put(f, "mu", vnums(o.mu, 2));
  put(f, "omega", vnums(o.omega, 2));
  put(f, "beta", fj::nav({o.beta[0], o.beta[1]}));
  put(f, "xi", fj::nav({o.xi[0], o.xi[1]}));
  Val cv = varr();
  for (int h = 0; h < 2; ++h) cv.a.push_back(o.conv[h] < 0 ? Val() : vbool(o.conv[h] == 1));
  put(f, "conv", cv);
  put(f, "iter", vnums(std::vector<double>{static_cast<double>(o.iter[0]), static_cast<double>(o.iter[1])}));
  std::vector<std::string> hy, who;
  std::vector<double> m;
  for (int h = 0; h < 2; ++h) {
    for (std::size_t k = 0; k < o.who[h].size(); ++k) {
      hy.push_back(h ? "hd" : "hp");
      who.push_back(o.who[h][k]);
      m.push_back(o.m[h][k]);
    }
  }
  put(mx, "hyp", vstrs(hy));
  put(mx, "who", vstrs(who));
  put(mx, "m", vnums(m));
  put(lc, "locus", vstrs(o.loc));
  put(lc, "log10lr", vnums(o.llr));
  put(out, "fit", f);
  put(out, "mix", mx);
  put(out, "locus", lc);
  return true;
}
