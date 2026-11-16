// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The footprint family (src/fp.cpp), with the checks of R/fp.R in their
// words. co2fp() evaluates an R expression and so has no entry here.

// One row of the carbon-intensity table as co2int() returns it.
inline Val cirow(int i, const std::vector<std::string> &note) {
  const char *loc, *iso3, *name, *src;
  double g;
  int yr;
  clab_cirow(i, &loc, &iso3, &name, &g, &yr, &src);
  Val o = vobj();
  put(o, "g_per_kwh", vnum(g));
  put(o, "location", vstr(loc));
  put(o, "name", vstr(name));
  put(o, "year", vnum(yr));
  put(o, "basis", vstr("lifecycle"));
  put(o, "source", vstr(src));
  put(o, "note", unbox(vstrs(note)));
  return o;
}

inline std::string upr(std::string s) {
  for (char &ch : s) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
  return s;
}

bool co2int(Ctx &c, Val &out) {
  const Val &l = c["loc"];
  if (l.t == Val::NUM || (l.t == Val::ARR && !l.a.empty() && l.a[0].t == Val::NUM)) {
    const double v = val1(l);
    if (!std::isfinite(v) || v < 0) {
      return ckstop(c, "a numeric `loc` must be a single non-negative gCO2e/kWh value");
    }
    out = vobj();
    put(out, "g_per_kwh", vnum(v));
    put(out, "location", vstr("user"));
    put(out, "name", vstr("user-supplied"));
    put(out, "year", Val());
    put(out, "basis", vstr("user"));
    put(out, "source", vstr("user-supplied carbon intensity"));
    put(out, "note", varr());
    return true;
  }
  std::string loc = "WORLD";
  if (c.has("loc")) {
    const Val *s = l.t == Val::ARR && l.a.size() == 1 ? &l.a[0] : &l;
    if (s->t != Val::STR) {
      return ckstop(c, "`loc` must be a single location code or name, a number in gCO2e/kWh, or NULL");
    }
    loc = s->s;
  }
  std::vector<std::string> note;
  int zone = 0;
  const int i = clab_cizone(loc.c_str(), &zone);
  if (zone != 0) {
    const char *cl, *iso3, *name, *src;
    double g;
    int yr;
    clab_cirow(i, &cl, &iso3, &name, &g, &yr, &src);
    note.push_back("no row for zone " + upr(loc) + "; its country " + cl + " used");
  }
  if (i < 0) {
    return ckstop(c,
                  "Unknown location '%s'. Use a code, alpha-3 code or name from co2inttab(), or pass a number in "
                  "gCO2e/kWh",
                  loc.c_str());
  }
  out = cirow(i, note);
  return true;
}

bool co2inttab(Ctx &, Val &out) {
  const int n = clab_cin();
  Val cols[7];
  for (Val &v : cols) v = varr();
  for (int i = 0; i < n; ++i) {
    const char *loc, *iso3, *name, *src;
    double g;
    int yr;
    clab_cirow(i, &loc, &iso3, &name, &g, &yr, &src);
    cols[0].a.push_back(vstr(loc));
    cols[1].a.push_back(vstr(iso3));
    cols[2].a.push_back(vstr(name));
    cols[3].a.push_back(vnum(g));
    cols[4].a.push_back(vnum(yr));
    cols[5].a.push_back(vstr("lifecycle"));
    cols[6].a.push_back(vstr(src));
  }
  static const char *const nm[] = {"location", "iso3", "name", "g_per_kwh", "year", "basis", "source"};
  out = vobj();
  for (int j = 0; j < 7; ++j) put(out, nm[j], cols[j]);
  return true;
}

inline std::string envs(const char *v) {
  const char *s = std::getenv(v);
  return s == nullptr ? "" : s;
}

// The system time zone as R's Sys.timezone() finds it on Linux and macOS:
// TZ, else the zone /etc/localtime links to.
inline std::string systz() {
  std::string tz = envs("TZ");
  if (!tz.empty()) return tz;
#ifndef _WIN32
  char b[4096];
  const ssize_t n = readlink("/etc/localtime", b, sizeof b - 1);
  if (n > 0) {
    b[n] = 0;
    const char *z = std::strstr(b, "zoneinfo/");
    if (z != nullptr) return z + 9;
  }
  std::ifstream f("/etc/timezone");
  if (f && std::getline(f, tz)) return tz;
#endif
  return "";
}

inline bool locarg(Ctx &c, const char *nm, std::string def, std::string &v) {
  v = def;
  if (!c.has(nm)) return true;
  const Val &o = c[nm];
  if (o.t == Val::STR) v = o.s;
  else if (o.t == Val::ARR && o.a.size() == 1 && o.a[0].t == Val::STR) v = o.a[0].s;
  else if (o.t == Val::ARR && o.a.empty()) v = "";
  else return ckstop(c, "`%s` must be one string (or NULL)", nm);
  return true;
}

inline void locign(Ctx &c, const std::string &v, const char *var) {
  c.warn.push_back(std::string(c.fn) + ": " + var + "='" + v +
                   "' is not a location co2inttab() knows (code, alpha-3 or name); it is ignored.");
}

bool dtctloc(Ctx &c, Val &out) {
  std::string env, cod, tz, lc;
  CK(locarg(c, "envloc", envs("CRIMLAB_LOCATION"), env));
  CK(locarg(c, "codloc", envs("CODECARBON_COUNTRY_ISO_CODE"), cod));
  CK(locarg(c, "tz", systz(), tz));
  std::string dl = envs("LC_ALL");
  if (dl.empty()) dl = envs("LC_TIME");
  if (dl.empty()) dl = envs("LANG");
  if (dl.empty()) dl = "C";
  CK(locarg(c, "locale", dl, lc));
  const char *loc = nullptr;
  std::size_t len = 0;
  int ig = 0;
  std::vector<char> how(256);
  if (clab_fploc(env.c_str(), cod.c_str(), tz.c_str(), lc.c_str(), &loc, how.data(), how.size(), &len, &ig) ==
          CLAB_OK &&
      len >= how.size()) {
    how.resize(len + 1);
    clab_fploc(env.c_str(), cod.c_str(), tz.c_str(), lc.c_str(), &loc, how.data(), how.size(), &len, &ig);
  }
  if ((ig & 1) != 0) locign(c, env, "CRIMLAB_LOCATION");
  if ((ig & 2) != 0) locign(c, cod, "CODECARBON_COUNTRY_ISO_CODE");
  out = vobj();
  put(out, "location", vstr(loc));
  put(out, "method", vstr(how.data()));
  return true;
}

bool raplok(Ctx &c, Val &out) {
  std::string dir = "/sys/class/powercap";
  if (c.has("rapldir")) {
    const Val &o = c["rapldir"];
    const Val *s = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
    if (s->t != Val::STR) return ckstop(c, "`rapldir` must be a single directory path");
    dir = s->s;
  }
  CK(ckpath(c, dir, "rapldir"));
  std::size_t len = 0;
  int n = 0;
  std::vector<char> b(4096);
  int rc = clab_rapldom(dir.c_str(), b.data(), b.size(), &len, &n);
  while (rc == CLAB_OK && len > b.size()) {
    b.resize(len);
    rc = clab_rapldom(dir.c_str(), b.data(), b.size(), &len, &n);
  }
  bool ok = rc == CLAB_OK && n > 0;
  std::size_t at = 0;
  for (int k = 0; ok && k < n; ++k) {
    double uj = 0, mx = 0;
    clab_raplrd(b.data() + at, &uj, &mx);
    ok = ok && !std::isnan(uj);
    at += std::strlen(b.data() + at) + 1;
  }
  out = vbool(ok);
  return true;
}

bool fpeqv(Ctx &c, Val &out) {
  const double g = val1(c["co2eg"]);
  if (!std::isfinite(g) || g < 0) return ckstop(c, "`co2eg` must be a single non-negative number");
  double e[2];
  clab_fpeqv(g, e);
  out = vobj();
  put(out, "carkm", vnum(e[0]));
  put(out, "treemon", vnum(e[1]));
  Val s = vobj();
  put(s, "car", vstr("175 gCO2/km, European passenger car (Green Algorithms context.csv)"));
  put(s, "tree", vstr("917 gCO2 per tree-month (Green Algorithms context.csv)"));
  put(out, "sources", s);
  return true;
}
