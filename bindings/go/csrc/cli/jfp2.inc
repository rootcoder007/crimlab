// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// co2fp() of R/fp.R for clab_jrun(): outside R there is no expression to
// time, so the caller gives what R measures around it, and gets the same
// result but `value`: `durs` (wall seconds) and `cputms` (CPU seconds of
// the computation); for usage = "machine" the machine's CPU accounting
// before and after (`stat`: {before: [idle, total], after: [idle, total]},
// as clab_cpuacc() reads /proc/stat), else the process utilisation is used
// with R's note; for RAPL the package counters before and after (`rapl`:
// {before: [uj...], after: [uj...], max: [uj...]}), else the model (meth
// "rapl" without them is refused as R refuses unreadable counters). The
// arithmetic is clab_fpe(), as in R; the checks are R's in their order and
// words.

namespace fp2 {

const ax::Cols modes = {"process", "machine", "fixed"};

// .fpnum(): NaN for null.
inline bool fpnum(Ctx &c, const char *what, bool positive, double &v) {
  v = NAN;
  const Val &o = c[what];
  if (o.t == Val::NUL) return true;
  const std::vector<double> x = ax::numv(o) ? ax::nums(o) : std::vector<double>();
  const bool ok = x.size() == 1 && std::isfinite(x[0]) && x[0] >= 0 && x[0] <= 1e9 && (!positive || x[0] > 0);
  if (!ok) return ckstop(c, "`%s` must be a single %s number", what, positive ? "positive" : "non-negative");
  v = x[0];
  return true;
}

// .ciget(): the carbon intensity, the detected location noted first.
inline bool ciget(Ctx &c, Val &ci) {
  Val a = vobj();
  std::string note;
  if (!c.has("loc")) {
    Val det;
    CK(ax::ccall(c, "dtctloc", dtctloc, vobj(), det));
    const std::string l = col(det, "location").s;
    note = "location " + l + " via " + col(det, "method").s;
    put(a, "loc", vstr(l));
  } else {
    put(a, "loc", c["loc"]);
  }
  CK(ax::ccall(c, "co2int", co2int, a, ci));
  if (note.empty()) return true;
  ax::Cols n = {note};
  for (std::size_t i = 0; i < ci.k.size(); ++i) {
    if (ci.k[i] != "note") continue;
    const Val &v = ci.a[i];
    if (v.t == Val::STR) n.push_back(v.s);
    for (const Val &e : v.a) n.push_back(e.s);
    ci.a[i] = vstrs(n);
  }
  return true;
}

// The arguments of .fpargs() but the paths: in[] filled.
struct Fpa {
  double cpupwrw = NAN, cores = NAN, memgb = NAN, mpw = 0.3725, pue = 1.0, fixed = NAN;
  int mode = 0;
};
inline bool fpargs(Ctx &c, Fpa &a) {
  if (c.has("pue")) {
    const std::vector<double> x = ax::numv(c["pue"]) ? ax::nums(c["pue"]) : std::vector<double>();
    if (x.size() != 1 || !std::isfinite(x[0]) || x[0] < 1) return ckstop(c, "`pue` must be a single number >= 1");
    a.pue = x[0];
  }
  CK(fpnum(c, "cpupwrw", true, a.cpupwrw));
  CK(fpnum(c, "cores", true, a.cores));
  if (!std::isnan(a.cores) && a.cores != std::round(a.cores)) return ckstop(c, "`cores` must be a whole number");
  CK(fpnum(c, "memgb", false, a.memgb));
  if (c.has("mempwrwpergb")) {
    CK(fpnum(c, "mempwrwpergb", false, a.mpw));
  }
  const Val &u = c["usage"];
  if (u.t == Val::NUL) return true;
  const Val *s = u.t == Val::ARR && u.a.size() == 1 ? &u.a[0] : &u;
  if (s->t == Val::ARR) {
    return ckstop(c, "`usage` takes one value: \"process\", \"machine\" or a single number in [0, 1]");
  }
  if (s->t == Val::STR) return ckarg(c, "usage", {"process", "machine"}, a.mode);
  if (s->t == Val::NUM && std::isfinite(s->n) && s->n >= 0 && s->n <= 1) {
    a.mode = 2;
    a.fixed = s->n;
    return true;
  }
  return ckstop(c, "`usage` must be \"process\", \"machine\" or a single number in [0, 1]");
}

// .fpdev(): kinds as codes (1 counted already), figures NaN unmeasured.
struct Dev {
  bool has = false;
  std::vector<int> kind;
  std::vector<double> w, e;
  ax::Cols k, v, nm, how;
};
inline bool devnum(Ctx &c, const Val &d, const char *cl, std::size_t n, std::vector<double> &v) {
  const Val *o = ax::fc(d, cl);
  v.assign(n, NA);
  if (o == nullptr) return true;
  if (!ax::numv(*o)) return ckstop(c, "`devices$%s` must be numeric", cl);
  v = ax::colv(*o);
  for (double x : v) {
    if (std::isinf(x) || (!std::isnan(x) && x < 0)) {
      return ckstop(c, "`devices$%s` must be non-negative and finite (NA: not measured)", cl);
    }
  }
  return true;
}
inline bool fpdev(Ctx &c, Dev &dv) {
  const Val &d = c["devices"];
  if (d.t == Val::NUL) return true;
  dv.has = true;
  bool frame = d.t == Val::OBJ;
  for (const Val &e : d.a) frame = frame && e.t == Val::ARR;
  if (!frame) return ckstop(c, "`devices` must be a data frame, one row per device");
  std::string miss;
  for (const char *k : {"kind", "vendor", "name", "how"}) {
    if (ax::fc(d, k) == nullptr) miss += (miss.empty() ? "" : ", ") + std::string(k);
  }
  if (!miss.empty()) return ckstop(c, "`devices` is missing columns: %s", miss.c_str());
  if (ax::fc(d, "watts") == nullptr && ax::fc(d, "energy") == nullptr) {
    return ckstop(c, "`devices` needs a `watts` or an `energy` column");
  }
  const std::size_t n = ax::fc(d, "kind")->a.size();
  CK(devnum(c, d, "watts", n, dv.w));
  CK(devnum(c, d, "energy", n, dv.e));
  for (const Val &e : ax::fc(d, "kind")->a) dv.k.push_back(ax::vlab(e));
  for (const Val &e : ax::fc(d, "vendor")->a) dv.v.push_back(ax::vlab(e));
  for (const Val &e : ax::fc(d, "name")->a) dv.nm.push_back(ax::vlab(e));
  for (const Val &e : ax::fc(d, "how")->a) dv.how.push_back(ax::vlab(e));
  for (std::string s : dv.k) {
    s.erase(0, s.find_first_not_of(" \t\r\n"));
    s.erase(s.find_last_not_of(" \t\r\n") + 1);
    for (char &ch : s) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    dv.kind.push_back(s == "cpu" || s == "memory" ? 1 : 0);
  }
  return true;
}

// Two numbers of an object's member, false when it is absent or not so.
inline bool pair(const Val &o, const char *k, std::size_t n, std::vector<double> &v) {
  const Val *m = ax::fc(o, k);
  if (m == nullptr || !ax::numv(*m)) return false;
  v = ax::nums(*m);
  return n == 0 || v.size() == n;
}

// .fpdevout(): the devices as R reports them.
inline Val devout(const Dev &dv, const std::vector<double> &kwh) {
  ax::Cols note;
  std::vector<int> cnt;
  std::vector<double> k;
  for (std::size_t i = 0; i < dv.kind.size(); ++i) {
    std::string s = "no power or energy reading";
    if (!std::isnan(dv.w[i])) s = "average power over the run";
    if (!std::isnan(dv.e[i])) s = "measured energy";
    if (dv.kind[i] == 1) s = "already counted as CPU or memory energy";
    note.push_back(s);
    cnt.push_back(std::isnan(kwh[i]) ? 0 : 1);
    k.push_back(std::isnan(kwh[i]) ? NA : kwh[i]);
  }
  Frame f;
  f.str("kind", dv.k);
  f.str("vendor", dv.v);
  f.str("name", dv.nm);
  f.str("how", dv.how);
  f.num("kwh", k);
  f.lgl("counted", cnt.data(), cnt.size());
  f.str("note", note);
  return f.o;
}

// .fpcores() and .fpmem(): the cores and the memory, with their notes.
inline void cores(double &n, ax::Cols &asm_) {
  if (!std::isnan(n)) return;
#if defined(_SC_NPROCESSORS_ONLN)
  const auto k = sysconf(_SC_NPROCESSORS_ONLN);
  if (k > 0) n = static_cast<double>(k);
#endif
  if (std::isnan(n)) {
    n = 1;
    asm_.push_back("core count unknown; 1 assumed");
  }
}
inline std::string mem(double &gb) {
  if (!std::isnan(gb)) return "";
  gb = clab_memgb();
  if (std::isnan(gb)) {
    gb = 8;
    return "memory size unknown: 8 GB assumed";
  }
  char b[96];
  std::snprintf(b, sizeof b, "memory charged: the machine's total %.1f GB", gb);
  return b;
}

// .fpcites()
inline Val cites(const Val &ci) {
  ax::Cols s = {"Lannelongue, Grealey & Inouye (2021) Adv Sci 8:2100707 (energy model, 0.3725 W/GB)",
                "codecarbon (RAPL method; process/machine tracking; 85 W x 50% unknown-CPU fallback)"};
  const Val &yr = col(ci, "year");
  if (yr.t == Val::NUM) {
    char b[64];
    std::snprintf(b, sizeof b, "%d", sat<int>(yr.n));
    s.push_back("carbon intensity: " + col(ci, "source").s + " (" + col(ci, "name").s + ", " + b + ")");
  }
  return vstrs(s);
}

}  // namespace fp2

// co2fp()
bool co2fp(Ctx &c, Val &out) {
  using namespace fp2;
  int meth = 0, lc = 0;
  CK(ckarg(c, "meth", {"auto", "rapl", "model"}, meth));
  CK(ckarg(c, "loadcur", {"linear", "codecarbon"}, lc));
  Val ci;
  CK(ciget(c, ci));
  Fpa a;
  CK(fpargs(c, a));
  Dev dv;
  CK(fpdev(c, dv));
  // R reaches `expr` only after every check above, so a missing `durs` is refused here
  if (!c.has("durs")) return ckstop(c, "`durs` is missing");
  double durs = 0.0, cpu = 0.0;
  if (!hw::one(c["durs"], durs) || durs < 0) return ckstop(c, "`durs` must be a single non-negative number of seconds");
  if (c.has("cputms") && (!hw::one(c["cputms"], cpu) || cpu < 0)) {
    return ckstop(c, "`cputms` must be a single non-negative number of seconds");
  }
  // the counters: before, after and their ranges
  std::vector<double> r0, r1, rmx;
  const bool rp = meth != 2 && c["rapl"].t == Val::OBJ && pair(c["rapl"], "before", 0, r0) &&
                  pair(c["rapl"], "after", r0.size(), r1) && pair(c["rapl"], "max", r0.size(), rmx);
  if (meth == 1 && !rp) {
    return ckstop(c, "RAPL energy counters are not readable here (see raplok()); use meth = \"model\"");
  }
  if (!rp) r0.clear();
  ax::Cols asm_;
  const Val &cn = col(ci, "note");
  if (cn.t == Val::STR) asm_.push_back(cn.s);
  for (const Val &e : cn.a) asm_.push_back(e.s);
  std::vector<double> s0 = {NAN, NAN}, s1 = {NAN, NAN};
  if (a.mode == 1 && !(c["stat"].t == Val::OBJ && pair(c["stat"], "before", 2, s0) && pair(c["stat"], "after", 2, s1))) {
    a.mode = 0;
    s0 = s1 = {NAN, NAN};
    asm_.push_back("machine CPU accounting unreadable (no /proc/stat); process utilisation used");
  }
  cores(a.cores, asm_);
  const std::string mn = mem(a.memgb);
  const double in[16] = {durs,  cpu,    a.cores,  static_cast<double>(a.mode), a.fixed, s0[0], s0[1], s1[0], s1[1],
                         a.cpupwrw, static_cast<double>(lc), a.memgb, a.mpw, a.pue, col(ci, "g_per_kwh").n,
                         rp ? 1.0 : 0.0};
  const std::size_t nd = dv.kind.size();
  std::vector<double> o(10, 0.0), dout(nd, 0.0);
  CK(kchk(c, clab_fpe(in, sat<int>(r0.size()), r0.data(), r1.empty() ? r0.data() : r1.data(), rmx.data(),
                      sat<int>(nd), dv.w.data(), dv.e.data(), dv.kind.data(), o.data(), dout.data())));
  const int um = sat<int>(o[2]);
  for (int i = 0; i < 5; ++i) {
    const char *f = clab_fpnote(i, sat<int>(o[4]), um == 1 ? 1 : 0);
    if (f != nullptr) asm_.push_back(f);
  }
  if (!mn.empty()) asm_.push_back(mn);
  if (a.pue > 1) {
    char b[48];
    std::snprintf(b, sizeof b, "PUE %.2f applied", a.pue);
    asm_.push_back(b);
  }
  Val en = vobj();
  put(en, "cpu", vnum(o[5]));
  put(en, "memory", vnum(o[6]));
  put(en, "total", vnum(o[8]));
  if (dv.has) put(en, "devices", vnum(o[7]));
  out = vobj();
  put(out, "durs", vnum(o[0]));
  put(out, "cputms", vnum(cpu));
  put(out, "usage", vnum(o[1]));
  put(out, "usagemode", vstr(modes[sat<std::size_t>(um)]));
  put(out, "cores", vnum(a.cores));
  put(out, "memgb", vnum(a.memgb));
  put(out, "enekwh", en);
  put(out, "co2eg", vnum(o[9]));
  put(out, "co2int", ci);
  put(out, "method", vstr(sat<int>(o[3]) == 0 ? "rapl" : "model"));
  put(out, "assumptions", vstrs(asm_));
  put(out, "citations", cites(ci));
  if (dv.has) put(out, "devices", devout(dv, dout));
  return true;
}
