// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The age-period-cohort model (src/bapc.cpp).

// A matrix argument as num() reads its values: column-major into v.
bool apmat(Ctx &c, const Val &o, const char *nm, double lo, bool open, int *dim, std::vector<double> &v) {
  dim[0] = dim[1] = -1;
  if (!shape(o, 2, dim, 0)) return false;
  v.assign(zz(dim[0]) * zz(dim[1]), 0.0);
  bool bad = false;
  fill(o, 2, dim, 0, 0, 1, v, bad, true);
  for (double e : v) {
    if (!std::isfinite(e)) return c.fail("`%s` must be finite numbers with no missing values.", nm);
  }
  for (double e : v) {
    if (open ? e <= lo : e < lo) {
      return c.fail("`%s` must be %s%s.", nm, open ? "greater than " : "at least ", rfmt(lo).c_str());
    }
  }
  return true;
}

bool apcie(Ctx &c, Val &out) {
  int dy[2], de[2];
  std::vector<double> yv, ev;
  if (!isnum(c["y"].t == Val::ARR && !c["y"].a.empty() ? c["y"].a[0] : Val())) {
    return c.fail("`y` must be a numeric matrix of counts, ages by periods.");
  }
  dy[0] = dy[1] = -1;
  if (!shape(c["y"], 2, dy, 0)) return c.fail("`y` must be a numeric matrix of counts, ages by periods.");
  const int na = dy[0], np = dy[1];
  if (na < 3 || np < 3) return c.fail("`y` must have three ages and three periods at least.");
  if (static_cast<double>(na) * np > 1e6) return c.fail("`y` has more than 1e6 cells.");
  CK(apmat(c, c["y"], "y", 0.0, false, dy, yv));
  double tot = 0;
  for (double e : yv) tot += e;
  if (!(tot > 0)) return c.fail("`y` has no events: every count is 0.");
  de[0] = de[1] = -1;
  if (!shape(c["expo"], 2, de, 0) || de[0] != na || de[1] != np) {
    return c.fail("`expo` must be a matrix of the dimensions of `y`, %d by %d.", na, np);
  }
  CK(apmat(c, c["expo"], "expo", 0.0, true, de, ev));
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  int mi = 0;
  CK(wholed(c, "maxit", 100, 1, mi));
  const int nc = na + np - 1, ne = 1 + na + np + nc;
  std::vector<double> est(zz(ne)), se(zz(ne)), st(5);
  const int rc = clab_apcie(yv.data(), ev.data(), na, np, mi, est.data(), se.data(), st.data());
  if (rc == CLAB_EINTR) return kchk(c, rc);
  if (rc != CLAB_OK) return c.fail("the Poisson model could not be fitted to these counts.");
  const double z = qnorm(0.5 + 0.5 * lv);
  std::vector<std::string> tm, lb;
  std::vector<double> lo, hi;
  tm.push_back("(Intercept)");
  lb.push_back("");
  for (int k = 1; k < ne; ++k) {
    const bool ag = k <= na, pe = !ag && k <= na + np;
    tm.push_back(ag ? "age" : pe ? "period" : "cohort");
    lb.push_back(std::to_string(ag ? k : pe ? k - na : k - na - np));
  }
  for (int k = 0; k < ne; ++k) {
    lo.push_back(est[zz(k)] - z * se[zz(k)]);
    hi.push_back(est[zz(k)] + z * se[zz(k)]);
  }
  Frame f;
  f.str("term", tm);
  f.str("level", lb);
  f.num("est", est);
  f.num("se", se);
  f.num("lo", lo);
  f.num("hi", hi);
  static const char *const sn[] = {"deviance", "df", "loglik", "iter", "conv"};
  const double dm[3] = {static_cast<double>(na), static_cast<double>(np), static_cast<double>(nc)};
  static const char *const dn[] = {"ages", "periods", "cohorts"};
  out = vobj();
  put(out, "coef", f.o);
  put(out, "fit", nvec(5, sn, st.data()));
  put(out, "dims", nvec(3, dn, dm));
  return true;
}
