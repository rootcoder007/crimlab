// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// prechk() of R/prechk.R for clab_jrun(), on the preset table of jprep.inc
// and the publishers' codes of inst/presets/codes.tsv (codes.inc, made by
// codes.sh; CI checks that it is up to date): the columns a preset reads
// (.pkneed(), .rin(), .fcols()), the labels outside the codes and those a
// mapping leaves without a value (.pkmsg(), .unk(), .known()), with R's
// checks and messages. A table that passes is returned as given. Labels
// are put in lower case byte by byte (ASCII letters only).

namespace pk {

using pp::Col;
using pp::Tbl;
using Strs = std::vector<std::string>;

const char *const ctsv[][4] = {
#include "codes.inc"
};

// The columns each recipe adds (.rout).
const std::vector<std::pair<std::string, Strs>> rout = {
    {"long", {"count", "with"}}, {"acvm", {"ac", "vm", "hiac"}}, {"acnv", {"acn", "vmany", "hiacn"}},
    {"pyr", {"np", "nreg", "ac"}}, {"pairb", {"tb", "yb"}},    {"pairc", {"tc", "yc"}},
    {"gent", {"gent"}},           {"pop", {"pop"}}};

inline bool in(const Strs &v, const std::string &s) { return std::find(v.begin(), v.end(), s) != v.end(); }
inline void add(Strs &v, const std::string &s) {
  if (!in(v, s)) v.push_back(s);
}
// trimws() and tolower()
inline std::string trim(const std::string &s) {
  const std::size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}
inline std::string low(std::string s) {
  for (char &ch : s) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  return s;
}
// x[1:k] of split pieces: NA ("" dropped later) past their end
inline Strs head(Strs v, std::size_t k) {
  v.resize(k);
  return v;
}

// .rin(): the columns a recipe reads.
inline Strs rin(const std::string &name, const std::string &arg) {
  const Strs ot(pp::otc, pp::otc + 7);
  if (name == "long" || name == "keep") return pp::split(arg, ',');
  if (name == "gent") return head(pp::split(arg, ','), 4);
  if (name == "pop") return head(pp::split(arg, ';'), 2);
  if (name == "pyr") return ot;
  if (name == "acvm" || name == "acnv") {
    Strs v = ot;
    v.push_back("Gender");
    v.push_back("Age_Category");
    return v;
  }
  // built by push_back: a braced list of strings here trips gcc 16's
  // -Wfree-nonheap-object once the caller's loop is inlined
  Strs v;
  if (name == "pairb") {
    for (std::size_t i : {0u, 1u, 4u, 5u, 6u}) v.push_back(ot[i]);
    v.push_back("Number_Of_Placements");
  } else if (name == "pairc") {
    v.push_back("Region_AtTimeOfPlacement");
    v.push_back("Region_MostRecentPlacement");
    v.push_back("NumberConsecutiveDays_Segregation");
  }
  return v;
}

// .fcols(): the columns the filters read ("NA" for an empty rule).
inline Strs fcols(const pp::Rows &r) {
  Strs v;
  for (const auto &x : r) {
    if (x[5].empty()) continue;
    for (const std::string &f : pp::split(x[5], ';')) {
      const Strs kv = pp::split(f, '=');
      v.push_back(kv.empty() ? "NA" : kv[0]);
    }
  }
  return v;
}

// .pkneed(): the recipes' inputs, then the columns and filters read that
// no recipe adds.
inline Strs need(const pp::Rows &r) {
  Strs ins, outs, read, out;
  for (const auto &x : r) {
    if (x[2] != "recipe") continue;
    for (const std::string &s : rin(x[3], x[4])) {
      if (!s.empty()) ins.push_back(s);
    }
    for (const auto &o : rout) {
      if (o.first == x[3]) outs.insert(outs.end(), o.second.begin(), o.second.end());
    }
  }
  Strs rd;
  for (const auto &x : r) {
    if (x[2] != "recipe" && !x[3].empty()) rd.push_back(x[3]);
  }
  for (const std::string &s : fcols(r)) rd.push_back(s);
  for (const std::string &s : rd) {
    if (!in(outs, s)) add(read, s);
  }
  for (const std::string &s : ins) add(out, s);
  for (const std::string &s : read) add(out, s);
  return out;
}

// .unk(): the labels of x outside known(), as one line, or "".
template <class F>
std::string unk(const Col &x, const std::string &col, const char *why, F known) {
  Strs u;
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (!x.miss(i)) add(u, trim(pp::lab(x, i)));
  }
  Strs bad;
  for (const std::string &s : u) {
    if (!known(low(s))) bad.push_back(s);
  }
  if (bad.empty()) return "";
  return "column " + col + " has labels " + why + ": " + pp::tostr(bad) + ".";
}

// .known(): whether a mapping rule covers a label in lower case.
inline bool known(const std::string &rule, const std::string &v) {
  for (const std::string &piece : pp::split(rule, ';')) {
    const Strs kv = pp::split(piece, '=');
    if (kv.empty()) continue;
    const std::string &k = kv[0];
    if (k.compare(0, 4, "has:") == 0) {
      if (v.find(low(k.substr(4))) != std::string::npos) return true;
    } else if (in(pp::split(k, '|'), v)) {
      return true;
    }
  }
  return false;
}

// A mapping without a `*` rule: (^|;)\*= absent.
inline bool nostar(const std::string &s) { return s.compare(0, 2, "*=") != 0 && s.find(";*=") == std::string::npos; }

}  // namespace pk

// prechk()
bool prechk(Ctx &c, Val &out) {
  using namespace pp;
  Rows tab;
  for (const auto &r : ptsv) {
    std::array<std::string, 9> x;
    for (std::size_t j = 0; j < 9; ++j) x[j] = r[j];
    tab.push_back(x);
  }
  const Val &pv = c["preset"];
  Tbl pt;
  Rows r;
  if (jframe(pv, pt)) {
    if (!pt.col("preset") || !pt.col("generic") || !pt.col("arg") || pt.n == 0) {
      return c.fail("a preset table needs the columns preset, generic and arg, and a row.");
    }
    r.assign(pt.n, std::array<std::string, 9>());
    for (std::size_t j = 0; j < 9; ++j) {
      const Col *x = pt.col(pcol[j]);
      for (std::size_t i = 0; x && i < pt.n; ++i) r[i][j] = x->miss(i) ? "" : lab(*x, i);
    }
    for (const auto &q : r) {
      if (q[0] != r[0][0]) return c.fail("a preset table must hold one preset.");
    }
  } else {
    const Val *s = pv.t == Val::ARR && pv.a.size() == 1 ? &pv.a[0] : &pv;
    if (s->t != Val::STR) return c.fail("`preset` must be one name or a preset table.");
    for (const auto &x : tab) {
      if (x[0] == s->s) r.push_back(x);
    }
    if (r.empty()) return c.fail("no preset named %s; clprep() lists them.", s->s.c_str());
  }
  Tbl d;
  if (!jframe(c["data"], d)) return c.fail("`data` must be a data frame.");
  const std::string name = r[0][0];
  // .pkcodes(): (column, value) pairs
  std::vector<std::pair<std::string, std::string>> codes;
  if (!c.has("codes")) {
    const std::string ag = name.substr(0, name.find('-'));
    for (const auto &k : pk::ctsv) {
      if (ag == k[0]) codes.emplace_back(k[1], k[2]);
    }
  } else {
    Tbl k;
    if (!jframe(c["codes"], k) || !k.col("column") || !k.col("value")) {
      return c.fail("`codes` must be a data frame with columns column and value.");
    }
    const Col &kc = *k.col("column"), &kv = *k.col("value");
    for (std::size_t i = 0; i < k.n; ++i) {
      codes.emplace_back(kc.miss(i) ? "NA" : lab(kc, i), kv.miss(i) ? "NA" : lab(kv, i));
    }
  }
  pk::Strs bad;
  for (const auto &x : r) {
    if (x[2] != "recipe" || x[3] == "keep") continue;
    bool known = false;
    for (const auto &o : pk::rout) known = known || o.first == x[3];
    if (!known) pk::add(bad, x[3]);
  }
  if (!bad.empty()) return c.fail("unknown recipe %s (preset %s).", tostr(bad).c_str(), name.c_str());
  const pk::Strs need = pk::need(r);
  pk::Strs msg, miss;
  for (const std::string &s : need) {
    if (!d.col(s)) miss.push_back(s);
  }
  if (!miss.empty()) msg.push_back("`data` has no column " + tostr(miss) + ".");
  pk::Strs cc;
  for (const auto &k : codes) pk::add(cc, k.first);
  for (const std::string &col : need) {
    if (!d.col(col) || !pk::in(cc, col)) continue;
    pk::Strs ok;
    for (const auto &k : codes) {
      if (k.first == col) ok.push_back(pk::low(pk::trim(k.second)));
    }
    const std::string m = pk::unk(*d.col(col), col, "not among its codes", [&ok](const std::string &l) {
      return pk::in(ok, l);
    });
    if (!m.empty()) msg.push_back(m);
  }
  for (const auto &x : r) {
    const std::string &s = x[4];
    const bool map = x[2] != "recipe" && s.find('=') != std::string::npos && s.compare(0, 6, "const:") != 0 &&
                     pk::nostar(s) && d.col(x[3]);
    if (!map) continue;
    const std::string m = pk::unk(*d.col(x[3]), x[3], "with no recode", [&s](const std::string &l) {
      return pk::known(s, l);
    });
    if (!m.empty()) msg.push_back(m);
  }
  if (!msg.empty()) {
    std::string t;
    for (const std::string &m : msg) t += "\n* " + m;
    return c.fail("the table does not fit preset %s:%s", name.c_str(), t.c_str());
  }
  out = c["data"];
  return true;
}
