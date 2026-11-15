// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Scaling, ideal points and counts over time (src/svk.h, src/bsvt.cpp):
// the JSON arguments into the kernels' node tree, R's defaults from the
// kernel table, and the result written as R's JSON writes the R value
// (a matrix as an array of rows, an array as an array of matrices, a
// vector of one as a scalar, NA as null).

std::vector<std::string> words(const char *s);

bool svscal(const Val &e) { return e.t == Val::NUM || e.t == Val::NUL || e.t == Val::BOOL || e.t == Val::STR; }

bool tonode(const Val &v, svk::Node &n) {
  double x = 0;
  n = svk::Node();
  switch (v.t) {
    case Val::NUL:
      return true;
    case Val::BOOL:
      n = svk::lgl(v.b);
      return true;
    case Val::NUM:
      n = svk::num(v.n);
      return true;
    case Val::STR:
      n = snum(v, x) ? svk::num(x) : svk::str(v.s);
      return true;
    case Val::OBJ:
      n = svk::list();
      for (std::size_t i = 0; i < v.a.size(); ++i) {
        svk::Node e;
        if (!tonode(v.a[i], e)) return false;
        n.add(v.k[i].c_str(), std::move(e));
      }
      return true;
    case Val::ARR:
      break;
  }
  bool scal = true, arr = true, str = false, num = false, lg = false;
  for (const Val &e : v.a) {
    scal = scal && svscal(e);
    arr = arr && e.t == Val::ARR;
    if (e.t == Val::STR && !snum(e, x)) str = true;
    if (e.t == Val::NUM || (e.t == Val::STR && snum(e, x))) num = true;
    if (e.t == Val::BOOL) lg = true;
  }
  if (v.a.empty() || (scal && !str)) {
    n.k = lg && !num ? svk::Node::LGL : svk::Node::NUM;
    for (const Val &e : v.a) {
      if (e.t == Val::NUM) {
        n.d.push_back(e.n);
      } else if (e.t == Val::BOOL) {
        n.d.push_back(e.b ? 1.0 : 0.0);
      } else if (e.t == Val::STR && snum(e, x)) {
        n.d.push_back(x);
      } else {
        n.d.push_back(NaN);
      }
    }
    if (lg && num) return false;
    return true;
  }
  if (scal && str) {
    n.k = svk::Node::STR;
    for (const Val &e : v.a) {
      if (e.t != Val::STR) return false;
      n.s.push_back(e.s);
    }
    return true;
  }
  std::vector<svk::Node> el;
  for (const Val &e : v.a) {
    svk::Node t;
    if (!tonode(e, t)) return false;
    el.push_back(std::move(t));
  }
  // rows of one length make a matrix, matrices of one shape an array
  bool same = arr;
  for (const svk::Node &t : el) {
    same = same && (t.k == svk::Node::NUM || t.k == svk::Node::LGL) && t.dim.size() == el[0].dim.size() &&
           t.d.size() == el[0].d.size() && t.dim == el[0].dim && t.dim.size() <= 2;
  }
  if (same && !el.empty()) {
    const int nr = static_cast<int>(el.size());
    n.k = svk::Node::NUM;
    if (el[0].dim.empty()) {
      const int nc = static_cast<int>(el[0].d.size());
      n.d.assign(zz(nr) * zz(nc), 0.0);
      for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nc; ++j) n.d[zz(j) * zz(nr) + zz(i)] = el[zz(i)].d[zz(j)];
      }
      n.dim = {nr, nc};
    } else {
      const int b = el[0].dim[0], c = el[0].dim[1];
      n.d.assign(zz(nr) * zz(b) * zz(c), 0.0);
      for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < b; ++j) {
          for (int k = 0; k < c; ++k) n.d[zz(i) + zz(nr) * (zz(j) + zz(b) * zz(k))] = el[zz(i)].d[zz(j) + zz(b) * zz(k)];
        }
      }
      n.dim = {nr, b, c};
    }
    return true;
  }
  n = svk::list();
  n.l = std::move(el);
  return true;
}

Val svone(const svk::Node &n, double x) {
  if (std::isnan(x)) return Val();
  return n.k == svk::Node::LGL ? vbool(x != 0) : vnum(x);
}

Val toval(const svk::Node &n) {
  switch (n.k) {
    case svk::Node::NUL:
      return Val();
    case svk::Node::STR: {
      if (n.s.size() == 1) return vstr(n.s[0]);
      Val a = varr();
      for (const std::string &s : n.s) a.a.push_back(vstr(s));
      return a;
    }
    case svk::Node::LIST: {
      Val o = n.nm.empty() ? varr() : vobj();
      for (std::size_t i = 0; i < n.l.size(); ++i) {
        if (n.nm.empty()) {
          o.a.push_back(toval(n.l[i]));
        } else {
          put(o, n.nm[i].c_str(), toval(n.l[i]));
        }
      }
      return o;
    }
    default:
      break;
  }
  if (n.dim.size() == 2 || n.dim.size() == 3) {
    const int a = n.dim[0], b = n.dim[1], c = n.dim.size() == 3 ? n.dim[2] : 1;
    Val o = varr();
    for (int i = 0; i < a; ++i) {
      if (n.dim.size() == 2) {
        Val r = varr();
        for (int j = 0; j < b; ++j) r.a.push_back(svone(n, n.d[zz(j) * zz(a) + zz(i)]));
        o.a.push_back(r);
      } else {
        Val m = varr();
        for (int j = 0; j < b; ++j) {
          Val r = varr();
          for (int k = 0; k < c; ++k) r.a.push_back(svone(n, n.d[zz(i) + zz(a) * (zz(j) + zz(b) * zz(k))]));
          m.a.push_back(r);
        }
        o.a.push_back(m);
      }
    }
    return o;
  }
  if (n.d.size() == 1) return svone(n, n.d[0]);
  Val o = varr();
  for (double x : n.d) o.a.push_back(svone(n, x));
  return o;
}

// A default of the kernel table: TRUE, FALSE, numbers (comma separated)
// or a string.
svk::Node svdef(const std::string &s) {
  if (s == "TRUE" || s == "FALSE") return svk::lgl(s == "TRUE");
  std::vector<double> v;
  std::size_t p = 0;
  bool ok = !s.empty();
  while (ok && p <= s.size()) {
    const std::size_t q = s.find(',', p);
    const std::string t = s.substr(p, q == std::string::npos ? std::string::npos : q - p);
    char *end = nullptr;
    const double x = std::strtod(t.c_str(), &end);
    ok = !t.empty() && end && *end == '\0';
    v.push_back(x);
    if (q == std::string::npos) break;
    p = q + 1;
  }
  if (!ok) return svk::str(s);
  return v.size() == 1 ? svk::num(v[0]) : svk::vec(v);
}

bool svrun(Ctx &c, Val &out) {
  const svk::Fn *f = svk::find(c.fn);
  svk::Node a = svk::list();
  for (const std::string &t : words(f->args)) {
    const std::size_t eq = t.find('=');
    const std::string nm = t.substr(0, std::min(eq, t.find('?')));
    svk::Node n;
    if (c.has(nm.c_str())) {
      if (!tonode(c[nm.c_str()], n)) return c.fail("`%s` must be numbers, strings or a list of them.", nm.c_str());
    } else if (eq != std::string::npos) {
      n = svdef(t.substr(eq + 1));
    }
    a.add(nm.c_str(), std::move(n));
  }
  svk::Node r;
  std::string err;
  if (!svk::run(c.fn, a, r, err)) return c.fail("%s", err.empty() ? "interrupted." : err.c_str());
  out = toval(r);
  return true;
}

#define SVJ(n) \
  bool n(Ctx &c, Val &out) { return svrun(c, out); }
SVK_ALL(SVJ)
#undef SVJ
