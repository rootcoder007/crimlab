// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// polrpttxt() of R/vrfp.R for clab_jrun(): the text of a vrfpol() report
// given as JSON (as R's report is written), line for line, with R's checks
// and words. Numbers are written as R's sprintf() and format() write them
// (rfmt7()). A line R builds from a part that is missing is left out, as
// sprintf() of NULL leaves it out. vrfpol() itself is R's.

namespace vrf {

// R's text of a number in sprintf(): NA, NaN, Inf, -Inf, else f.
inline std::string pf(const char *f, double x) {
  if (isna(x)) return "NA";
  if (std::isnan(x)) return "NaN";
  if (std::isinf(x)) return x > 0 ? "Inf" : "-Inf";
  char b[400];
  std::snprintf(b, sizeof b, f, x);
  return b;
}

// format(x) with R's 7 significant digits: the fewest digits up to 7 that
// show x, in fixed notation unless scientific is narrower; with big, fixed
// always and commas between thousands (format(x, big.mark = ",",
// scientific = FALSE)).
inline std::string rfmt7(double x, bool big) {
  if (isna(x)) return "NA";
  if (std::isnan(x)) return "NaN";
  if (std::isinf(x)) return x > 0 ? "Inf" : "-Inf";
  char b[400];
  std::snprintf(b, sizeof b, "%.6e", x);
  const std::string s = b;
  const std::size_t e = s.find('e');
  const int kp = std::atoi(s.c_str() + e + 1);
  std::string dg;
  for (std::size_t i = 0; i < e; ++i) {
    if (std::isdigit(static_cast<unsigned char>(s[i]))) dg += s[i];
  }
  int nsig = 7;
  while (nsig > 1 && dg[static_cast<std::size_t>(nsig) - 1] == '0') --nsig;
  const int neg = x < 0 ? 1 : 0;
  const int rgt = std::max(0, nsig - kp - 1);
  const int wf = neg + (kp >= 0 ? kp + 1 : 1) + (rgt > 0 ? rgt + 1 : 0);
  const int we = neg + (nsig > 1 ? nsig + 1 : 1) + (kp >= 100 || kp <= -100 ? 5 : 4);
  if (!big && wf > we) {
    std::snprintf(b, sizeof b, "%.*e", nsig - 1, x);
    return b;
  }
  std::snprintf(b, sizeof b, "%.*f", rgt, x);
  std::string o = b;
  if (!big) return o;
  const std::size_t dot = o.find('.') == std::string::npos ? o.size() : o.find('.');
  for (std::size_t k = dot; k > static_cast<std::size_t>(neg) + 3; k -= 3) o.insert(k - 3, ",");
  return o;
}

// A member of an object (null when absent).
inline const Val &mem(const Val &o, const char *k) {
  const Val *v = o.t == Val::OBJ ? ax::fc(o, k) : nullptr;
  return v == nullptr ? nul : *v;
}
// A number of a member, or false (the line is left out) when it is null.
inline bool num(const Val &o, const char *k, double &x) {
  const Val &v = mem(o, k);
  const Val &e = v.t == Val::ARR && v.a.size() == 1 ? v.a[0] : v;
  if (e.t != Val::NUM && e.t != Val::STR) return false;
  x = e.t == Val::NUM ? e.n : (e.s == "NaN" ? NaN : e.s == "Inf" ? HUGE_VAL : e.s == "-Inf" ? -HUGE_VAL : NA);
  return true;
}
// The texts of a member, one per element (none when null).
inline std::vector<std::string> txt(const Val &o, const char *k) {
  const Val &v = mem(o, k);
  std::vector<std::string> r;
  for (const Val &e : v.t == Val::ARR ? v.a : std::vector<Val>{v}) {
    if (e.t == Val::STR) r.push_back(e.s);
    else if (e.t == Val::NUM) r.push_back(rfmt(e.n));
    else if (e.t == Val::BOOL) r.push_back(e.b ? "TRUE" : "FALSE");
    else if (e.t == Val::NUL && v.t == Val::ARR) r.push_back("NA");
  }
  return r;
}
inline std::string up(std::string s) {
  for (char &ch : s) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
  return s;
}

// .vphead()
inline void head(const Val &r, std::vector<std::string> &l) {
  const std::string bar(66, '=');
  l.push_back(bar);
  const auto cm = txt(r, "command"), pl = txt(r, "pollutant"), oc = txt(r, "outcome");
  if (!cm.empty() && !pl.empty() && !oc.empty()) l.push_back("  " + cm[0] + " -- " + up(pl[0]) + " -> " + oc[0]);
  for (const std::string &s : txt(r, "region")) l.push_back("  region: " + s);
  for (const std::string &s : txt(r, "years")) l.push_back("  years:  " + s);
  for (const std::string &s : txt(r, "dsou")) l.push_back("  data:   " + s);
  l.push_back(bar);
  l.push_back("");
  l.push_back("Inputs");
  const Val &in = mem(r, "inputs");
  double x = 0;
  if (num(in, "expomean", x)) l.push_back("  exposure mean:      " + pf("%.3f", x));
  if (num(in, "expoprev", x)) l.push_back("  exposure prevalence:" + pf("%.3f", x));
  if (num(in, "baserateper100k", x)) l.push_back("  baseline rate/100k: " + pf("%.2f", x));
  if (num(in, "population", x)) l.push_back("  population:         " + rfmt7(x, true));
  if (num(in, "refconc", x)) l.push_back("  reference conc:     " + rfmt7(x, false));
  if (num(in, "drorows", x) && x > 0) l.push_back("  rows dropped:       " + pf("%.0f", x) + " (no exposure value)");
  l.push_back("");
  l.push_back("Assumption log");
}

// .vppipe()
inline void pipe(const Val &p, std::vector<std::string> &l) {
  const Val &ci = mem(p, "pafci"), &crf = mem(p, "crf"), &dp = mem(p, "displaced"), &bd = mem(p, "burden");
  std::string cis;
  double lo = 0, hi = 0, x = 0;
  if (ci.t != Val::NUL) {
    if (!num(ci, "lower", lo) || !num(ci, "upper", hi)) lo = hi = NA;
    cis = "  (95% CI " + pf("%.4f", lo) + " to " + pf("%.4f", hi) + ", from the risk's limits)";
  }
  l.push_back("");
  l.push_back("Concentration-response");
  if (num(crf, "rr", x)) l.push_back("  RR:       " + pf("%.4f", x));
  for (const std::string &s : txt(crf, "citation")) l.push_back("  source:   " + s);
  l.push_back("");
  if (num(p, "paf", x)) l.push_back("Attributable fraction (PAF): " + pf("%.4f", x) + cis);
  l.push_back("");
  const auto un = txt(dp, "unit");
  l.push_back(un.size() == 1 && un[0] == "deaths" ? "Mortality displaced" : "Cases displaced");
  if (num(dp, "deadisp", x)) l.push_back("  expected avoided " + (un.empty() ? "deaths" : un[0]) + ": " + pf("%.1f", x));
  l.push_back("");
  l.push_back("Burden of pollution");
  const auto bu = txt(mem(bd, "extra"), "unit");
  if (num(bd, "attrcases", x)) l.push_back("  attributable " + (bu.empty() ? "deaths" : bu[0]) + ":   " + pf("%.1f", x));
  const Val &eq = mem(p, "equity");
  if (eq.t != Val::NUL) {
    l.push_back("");
    l.push_back("Equity analysis");
    if (num(eq, "concidx", x)) l.push_back("  concentration index: " + pf("%.4f", x));
  }
}

}  // namespace vrf

// polrpttxt()
bool polrpttxt(Ctx &c, Val &out) {
  using namespace vrf;
  const Val &r = c["report"];
  const Val &st = mem(r, "status");
  if (r.t != Val::OBJ || st.t != Val::STR) return ckstop(c, "`report` must be a report from vrfpol()");
  if (st.s == "error") {
    const auto e = txt(r, "error");
    out = vstr("ERROR: " + (e.empty() ? std::string() : e[0]) + "\n");
    return true;
  }
  const Val &in = mem(r, "inputs"), &as = mem(r, "assumptions"), &pp = mem(r, "pipeline");
  bool ok = in.t == Val::OBJ && (as.t == Val::ARR || as.t == Val::OBJ) && (st.s != "ok" || pp.t == Val::OBJ);
  for (const Val &a : as.a) {
    ok = ok && a.t == Val::OBJ && (mem(a, "ok").t == Val::BOOL || mem(a, "ok").t == Val::NUL) && ax::fc(a, "ok");
  }
  if (!ok) {
    return ckstop(c, "`report` is not a complete vrfpol() report (its inputs, assumptions or pipeline are missing)");
  }
  std::vector<std::string> l;
  head(r, l);
  for (const Val &a : as.a) {
    const Val &k = mem(a, "ok");
    const auto an = txt(a, "assumption"), nt = txt(a, "note");
    if (an.empty() || nt.empty()) continue;
    l.push_back(std::string("  [") + (k.t == Val::BOOL ? (k.b ? "PASS" : "FAIL") : "NA") + "] " + an[0] + " -- " + nt[0]);
  }
  std::string s;
  if (st.s != "assumption_failure") pipe(pp, l);
  for (std::size_t i = 0; i < l.size(); ++i) s += (i ? "\n" : "") + l[i];
  s += st.s == "assumption_failure" ? "\n\nSTATUS: assumption_failure (pipeline skipped)\n" : "\n\nSTATUS: ok\n";
  out = vstr(s);
  return true;
}
