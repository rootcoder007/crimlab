// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Text to plain ASCII (R/txt.R, src/txt.cpp) for clab_jrun(). JSON text is
// UTF-8 already, so asciifb() changes it only with force.

// .txtin(): text, numbers or true and false, null missing (miss[i]).
bool txtin(Ctx &c, const char *nm, std::vector<std::string> &s, std::vector<char> &miss) {
  s.clear();
  miss.clear();
  const Val &o = c[nm];
  if (o.t == Val::NUL) return true;
  auto one1 = [&](const Val &e) {
    std::string t;
    if (e.t == Val::NUL) {
      s.emplace_back();
      miss.push_back(1);
      return true;
    }
    if (!lab(e, t)) return false;
    s.push_back(t);
    miss.push_back(0);
    return true;
  };
  bool ok = true;
  if (o.t == Val::ARR || o.t == Val::OBJ) {
    for (const Val &e : o.a) ok = ok && one1(e);
  } else {
    ok = one1(o);
  }
  return ok || c.fail("`%s` must be text (or numbers, TRUE and FALSE).", nm);
}

std::string ascii1(const std::string &s) {
  std::string o(s.size() + 1, '\0');
  o.resize(sat<std::size_t>(clab_ascii(s.data(), sat<int>(s.size()), &o[0])));
  return o;
}

Val strsna(const std::vector<std::string> &s, const std::vector<char> &miss) {
  Val v = varr();
  for (std::size_t i = 0; i < s.size(); ++i) v.a.push_back(miss[i] ? Val() : vstr(s[i]));
  return unbox(v);
}

// not toascii: a macro of the Windows C headers
bool jtoasc(Ctx &c, Val &out) {
  std::vector<std::string> s;
  std::vector<char> m;
  CK(txtin(c, "x", s, m));
  for (std::string &t : s) t = ascii1(t);
  out = unbox(vstrs(s));
  return true;
}

bool asciifb(Ctx &c, Val &out) {
  std::vector<std::string> s;
  std::vector<char> m;
  CK(txtin(c, "x", s, m));
  const bool force = c["force"].t == Val::BOOL && c["force"].b;
  if (force) {
    for (std::string &t : s) t = ascii1(t);
    out = unbox(vstrs(s));
    return true;
  }
  out = strsna(s, m);
  return true;
}

bool wrtxtfb(Ctx &c, Val &out) {
  std::vector<std::string> s;
  std::vector<char> m;
  CK(txtin(c, "text", s, m));
  const Val &p = c["path"];
  if (p.t != Val::STR || p.s.empty()) return c.fail("`path` must be one non-empty file path.");
  CK(ckpath(c, p.s, "path"));
  std::string b;
  for (std::size_t i = 0; i < s.size(); ++i) {
    if (i) b += '\n';
    b += m[i] ? std::string("NA") : s[i];
  }
  b += '\n';
  std::FILE *f = std::fopen(p.s.c_str(), "wb");
  if (f == nullptr) return c.fail("cannot open file '%s'.", p.s.c_str());
  const bool ok = std::fwrite(b.data(), 1, b.size(), f) == b.size();
  if (std::fclose(f) != 0 || !ok) return c.fail("cannot write file '%s'.", p.s.c_str());
  out = vstr(p.s);
  return true;
}
