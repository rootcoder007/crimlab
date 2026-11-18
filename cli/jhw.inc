// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// hwpwr() of R/hwpwr.R for clab_jrun(), on the readers of src/pwr.h, with
// the checks of R in their order and words. The result is an object: the
// rows R returns (`rows`), the interval in seconds (`dt`) and the last
// reading as text (`snap`), which a later call takes as `since` (R keeps
// the last two as attributes of its data frame). `pid` defaults to the
// calling process. WebAssembly cannot wait: there `dt` must be 0, and two
// readings are taken with `since`.

namespace hw {

// One finite number.
inline bool one(const Val &o, double &v) {
  if (!ax::numv(o)) return false;
  const std::vector<double> x = ax::nums(o);
  if (x.size() != 1 || !std::isfinite(x[0])) return false;
  v = x[0];
  return true;
}

// Strings: one, or an array of them (none for null); false for anything
// else or a null among them.
inline bool strs(const Val &o, std::vector<std::string> &v) {
  v.clear();
  if (o.t == Val::NUL) return true;
  if (o.t == Val::STR) {
    v.push_back(o.s);
    return true;
  }
  if (o.t != Val::ARR) return false;
  for (const Val &e : o.a) {
    if (e.t != Val::STR) return false;
    v.push_back(e.s);
  }
  return true;
}

// Unique keys of an object, all among the names nm.
inline bool keys(const Val &o, const ax::Cols &nm) {
  for (std::size_t i = 0; i < o.k.size(); ++i) {
    if (!ax::has(nm, o.k[i])) return false;
    for (std::size_t j = 0; j < i; ++j) {
      if (o.k[j] == o.k[i]) return false;
    }
  }
  return true;
}

// .hwck(): dt, step, pid and skip.
inline bool ck(Ctx &c, double &dt, double &step, double &pid, std::vector<std::string> &skip) {
  dt = 1.0;
  step = 0.25;
  if (c.has("dt") && (!one(c["dt"], dt) || dt < 0 || dt > 604800)) {
    return c.fail("`dt` must be one number from 0 to 604800 (a week, in seconds).");
  }
  if (c.has("step") && (!one(c["step"], step) || step <= 0)) return c.fail("`step` must be one number above 0.");
#if !defined(_WIN32) && !defined(__wasi__)
  pid = static_cast<double>(getpid());
#else
  pid = 0.0;  // no process identifier there
#endif
  if (c.has("pid") && (!one(c["pid"], pid) || pid != std::round(pid) || pid < 0 || pid > 2147483647.0)) {
    return c.fail("`pid` must be a process id: a whole number from 0 to 2147483647.");
  }
  static const ax::Cols src = {"rapl", "gpu", "battery", "proc", "nvml", "smi", "ioreport", "ioreg", "emi"};
  bool ok = strs(c["skip"], skip);
  for (const std::string &s : skip) ok = ok && ax::has(src, s);
  if (!ok) {
    return c.fail("`skip` must name sources among rapl, gpu, battery, proc, nvml, smi, ioreport, ioreg, emi.");
  }
  return true;
}

// .hwpth(): the programs and roots read, the defaults overridden by paths.
inline bool pth(Ctx &c, pwr::Cfg &g) {
#ifdef __linux__
  g.smi = "nvidia-smi";
#endif
  if (std::ifstream("/usr/sbin/ioreg").good()) g.ioreg = "/usr/sbin/ioreg";
  const Val &p = c["paths"];
  if (p.t == Val::NUL || (p.t == Val::ARR && p.a.empty())) return true;
  // a named list: gcc 16 on aarch64 reports the destructor of a temporary
  // one here as freeing a pointer inside an object (-Wfree-nonheap-object)
  static const ax::Cols nm = {"root", "nvml", "smi", "ioreg"};
  bool ok = p.t == Val::OBJ && keys(p, nm);
  for (std::size_t i = 0; ok && i < p.k.size(); ++i) {
    const Val &e = p.a[i].t == Val::ARR && p.a[i].a.size() == 1 ? p.a[i].a[0] : p.a[i];
    ok = e.t == Val::STR;
    if (!ok) break;
    std::string &d = p.k[i] == "root" ? g.root : p.k[i] == "nvml" ? g.nvml : p.k[i] == "smi" ? g.smi : g.ioreg;
    d = e.s;
  }
  if (!ok) return c.fail("`paths` must be named strings among root, nvml, smi and ioreg.");
  return true;
}

// .hwrt(): the rated watts of the CPU and the GPU.
inline bool rt(Ctx &c, pwr::Rated &r) {
  const Val &o = c["rated"];
  if (o.t == Val::NUL) return true;
  static const ax::Cols nm = {"cpu", "gpu"};
  bool ok = o.t == Val::OBJ && !o.a.empty() && keys(o, nm);
  for (std::size_t i = 0; ok && i < o.k.size(); ++i) {
    double v = 0.0;
    ok = one(o.a[i], v) && v >= 0;
    if (ok && o.k[i] == "cpu") r.cpu = v;
    if (ok && o.k[i] == "gpu") r.gpu = v;
  }
  if (!ok) return c.fail("`rated` must be watts from 0 named cpu or gpu.");
  return true;
}

// .hwsin(): the snapshot text of since (a string, or an earlier result).
inline bool snc(Ctx &c, std::string &s, bool &has) {
  const Val &o = c["since"];
  has = o.t != Val::NUL;
  if (!has) return true;
  const Val *t = o.t == Val::OBJ ? ax::fc(o, "snap") : &o;
  if (t != nullptr && t->t == Val::ARR && t->a.size() == 1) t = &t->a[0];
  if (t == nullptr || t->t != Val::STR) return c.fail("`since` must be an earlier result of hwpwr().");
  s = t->s;
  return true;
}

// The seconds of a steady clock.
inline double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// .hwsnap(): one reading then one every step seconds for dt seconds, or
// the reading of since then one now.
inline std::vector<std::string> snaps(const pwr::Cfg &g, bool has, const std::string &since, double dt, double step) {
  std::vector<std::string> s;
  if (has) {
    s.push_back(since);
    s.push_back(pwr::put(pwr::snap(g)));
    return s;
  }
  s.push_back(pwr::put(pwr::snap(g)));
  const double end = now() + dt;
  while (dt > 0) {
    const double left = end - now();
    if (left <= 0) break;
#ifndef __wasi__
    std::this_thread::sleep_for(std::chrono::duration<double>(std::min(step, left)));
#else
    (void)step;  // nothing can wait there: dt is 0
#endif
    s.push_back(pwr::put(pwr::snap(g)));
  }
  return s;
}

// C_clab_r_pwrd(): the rows over the readings, with its messages.
inline bool rows(Ctx &c, const std::vector<std::string> &txt, const pwr::Rated &r, Val &out) {
  std::vector<pwr::Snap> s(txt.size());
  for (std::size_t i = 0; i < txt.size(); ++i) {
    if (!pwr::get(txt[i], s[i])) return c.fail("`since` is not a reading of hwpwr().");
    if (i > 0 && s[i].t < s[i - 1].t) {
      return c.fail("`since` is not an earlier reading of this session: the clock went back.");
    }
  }
  const std::vector<pwr::Row> rw = pwr::rows(s, r);
  ax::Cols kind, ven, name, scope, how, why;
  std::vector<double> w, j, u, lo, hi;
  auto na = [](double v) { return std::isnan(v) ? NA : v; };
  for (const pwr::Row &x : rw) {
    kind.push_back(x.kind);
    ven.push_back(x.vendor);
    name.push_back(x.name);
    scope.push_back(x.scope);
    how.push_back(x.how);
    why.push_back(x.why);
    w.push_back(na(x.watts));
    j.push_back(na(x.joules));
    u.push_back(na(x.util));
    lo.push_back(na(x.flo));
    hi.push_back(na(x.fhi));
  }
  Frame f;
  f.str("kind", kind);
  f.str("vendor", ven);
  f.str("name", name);
  f.str("scope", scope);
  f.num("watts", w);
  f.num("joules", j);
  f.num("util", u);
  f.str("how", how);
  f.num("flo", lo);
  f.num("fhi", hi);
  f.str("why", why);
  out = vobj();
  put(out, "rows", f.o);
  put(out, "dt", vnum(s.back().t - s.front().t));
  put(out, "snap", vstr(txt.back()));
  return true;
}

}  // namespace hw

// hwpwr()
bool hwpwr(Ctx &c, Val &out) {
  double dt = 1.0, step = 0.25, pid = 0.0;
  pwr::Cfg g;
  CK(hw::ck(c, dt, step, pid, g.skip));
  CK(hw::pth(c, g));
  pwr::Rated r;
  CK(hw::rt(c, r));
  std::string since;
  bool has = false;
  CK(hw::snc(c, since, has));
#ifdef __wasi__
  if (!has && dt > 0) return c.fail("`dt` above 0 cannot be waited for in WebAssembly: take two readings with `since`.");
#endif
  g.pid = static_cast<decltype(g.pid)>(pid);
  return hw::rows(c, hw::snaps(g, has, since, dt, step), r, out);
}
