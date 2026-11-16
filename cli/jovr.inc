// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Sentencing arithmetic, synthetic example data, published tables and
// the coding of police oversight reports (src/psc.cpp, src/sfx.cpp,
// src/ovr.cpp) for clab_jrun(): the R wrappers of R/psc.R, R/sfx.R,
// R/tabrpt.R and R/ovr.R in C++, with their checks in their order and
// words, and their results written as dev/clref.R writes R's.

namespace ovr {

// .pscone(): one value, or one per value of `days`.
inline bool pone(Ctx &c, const Nv &x, std::size_t n, const char *w) {
  if (x.n() != 1 && x.n() != n) return ckstop(c, "`%s` must have one value or one per value of `days`", w);
  return true;
}

bool pscred(Ctx &c, Val &out) {
  Nv d, r, s;
  CK(cknum(c, c["days"], "days", d, false, 0));
  if (c.has("ratio")) {
    CK(cknum(c, c["ratio"], "ratio", r, false, 0, 1.5));
  } else {
    r.v.push_back(1.5);
  }
  const std::size_t n = d.n();
  CK(pone(c, r, n, "ratio"));
  const bool hs = c.has("sentence");
  if (hs) {
    CK(cknum(c, c["sentence"], "sentence", s, false, 0));
    CK(pone(c, s, n, "sentence"));
  }
  std::vector<double> cr(n), rm(n), rr(n), ss(n);
  const int rc = clab_pscred(d.p(), static_cast<int64_t>(n), r.p(), static_cast<int64_t>(r.n()),
                             hs ? s.p() : nullptr, hs ? static_cast<int64_t>(s.n()) : 0, cr.data(), rm.data());
  if (rc != CLAB_OK) return kchk(c, rc);
  for (std::size_t i = 0; i < n; ++i) {
    rr[i] = r[r.n() == 1 ? 0 : i];
    if (hs) ss[i] = s[s.n() == 1 ? 0 : i];
  }
  out = vobj();
  put(out, "days", vnums(d.v));
  put(out, "ratio", vnums(rr));
  put(out, "credit", vnums(cr));
  if (hs) {
    put(out, "sentence", vnums(ss));
    put(out, "remain", vnums(rm));
  }
  return true;
}

}  // namespace ovr
