// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// yoyhtml() of R/yoyh.R for clab_jrun(): the self-contained HTML page of a
// change table, line for line as R writes it (the table of .yshow(), the
// cards of yoysmry(), the notes, the same styles), with R's checks in
// their order and words. The version in the page's last line is the one
// the build was given (CLAB_VER), "unknown" without it. R keeps its own
// writer. Without `file` (which R requires) the page is returned as text.

#ifndef CLAB_VER
#define CLAB_VER "unknown"
#endif

namespace yh {

using Strs = std::vector<std::string>;

// .YOY_PALETTES: bad, good, flat, badbg, goodbg
const char *const pals[][6] = {{"diverging", "#c2492f", "#1b7f6b", "#6b6b6b", "#fbeae6", "#e6f4f1"},
                               {"safe", "#b8621b", "#18567d", "#6b6b6b", "#fdf0e4", "#e7eff6"},
                               {"mono", "#1a1a1a", "#1a1a1a", "#6b6b6b", "#f0f0f0", "#f7f7f7"}};

inline void utf8(std::string &o, std::uint32_t k) {
  if (k < 0x80) {
    o += static_cast<char>(k);
  } else if (k < 0x800) {
    o += static_cast<char>(0xC0 | (k >> 6));
    o += static_cast<char>(0x80 | (k & 0x3F));
  } else if (k < 0x10000) {
    o += static_cast<char>(0xE0 | (k >> 12));
    o += static_cast<char>(0x80 | ((k >> 6) & 0x3F));
    o += static_cast<char>(0x80 | (k & 0x3F));
  } else {
    o += static_cast<char>(0xF0 | (k >> 18));
    o += static_cast<char>(0x80 | ((k >> 12) & 0x3F));
    o += static_cast<char>(0x80 | ((k >> 6) & 0x3F));
    o += static_cast<char>(0x80 | (k & 0x3F));
  }
}
// The next code point of s at i (i moved past it); false on bytes that
// are not UTF-8.
inline bool next(const std::string &s, std::size_t &i, std::uint32_t &k) {
  const unsigned char b = static_cast<unsigned char>(s[i]);
  std::size_t n = 0;
  if (b >= 0x80) {
    if ((b >> 5) == 6) n = 1;
    else if ((b >> 4) == 14) n = 2;
    else if ((b >> 3) == 30) n = 3;
    else return false;
  }
  if (i + n >= s.size() && n > 0) return false;
  k = n == 0 ? b : static_cast<std::uint32_t>(b & (0x3F >> n));
  for (std::size_t j = 1; j <= n; ++j) {
    const unsigned char c = static_cast<unsigned char>(s[i + j]);
    if ((c & 0xC0) != 0x80) return false;
    k = (k << 6) | (c & 0x3F);
  }
  i += n + 1;
  return true;
}
// .yctl(s, ws = TRUE): control characters shown escaped; text that is not
// UTF-8 is left as it is.
inline std::string yctl(const std::string &s) {
  std::string o;
  std::uint32_t k = 0;
  for (std::size_t i = 0; i < s.size();) {
    if (!next(s, i, k)) return s;
    char b[16];
    if (k == 9 || k == 10) {
      o += static_cast<char>(k);
    } else if (k == 7 || k == 8 || k == 12 || k == 13 || k == 27) {
      o += k == 7 ? "\\a" : k == 8 ? "\\b" : k == 12 ? "\\f" : k == 13 ? "\\r" : "\\033";
    } else if (k < 32 || k == 127) {
      std::snprintf(b, sizeof b, "\\%03o", static_cast<unsigned>(k));
      o += b;
    } else if ((k >= 0x80 && k <= 0x9f) || k == 0x200e || k == 0x200f || (k >= 0x202a && k <= 0x202e) ||
               (k >= 0x2066 && k <= 0x2069)) {
      std::snprintf(b, sizeof b, "<U+%04X>", static_cast<unsigned>(k));
      o += b;
    } else {
      utf8(o, k);
    }
  }
  return o;
}
// .yesc(): control characters escaped, then the five HTML characters.
inline std::string esc(const std::string &s0) {
  std::string o;
  for (char ch : yctl(s0)) {
    if (ch == '&') o += "&amp;";
    else if (ch == '<') o += "&lt;";
    else if (ch == '>') o += "&gt;";
    else if (ch == '"') o += "&quot;";
    else if (ch == '\'') o += "&#39;";
    else o += ch;
  }
  return o;
}

const char *const dash = "\xe2\x80\x94";

// .yvd(): the decimals of the values shown.
inline int yvd(const Yt &t) {
  if (t.get("pctlo")) return 0;
  const Yc &v = *t.get("value");
  int vd = 0;
  for (std::size_t i = 0; i < t.n; ++i) {
    if (v.na[i] || !std::isfinite(v.d[i])) continue;
    int k = 0;
    while (k < 6 && std::fabs(v.d[i] - rround(v.d[i], k)) >= 1e-9) ++k;
    vd = std::max(vd, k);
  }
  return vd;
}
inline std::string ynum(const Yc &v, std::size_t i, int d) {
  return v.na[i] || std::isnan(v.d[i]) ? dash : ynum::fix(v.d[i], d);
}
inline std::string ypct(double v, int d, bool pp) { return std::isnan(v) ? dash : ynum::pct(v, d, pp ? "pp" : "%"); }

// .yshow(md = FALSE): the heads and the columns shown, as text.
struct Show {
  Strs hd;
  std::vector<Strs> cols;
  std::vector<char> gna;
};
inline Show yshow(const Yt &t, int dg, bool intv, const char *na = "") {
  Show s;
  const int vd = yvd(t);
  const bool pp = t.get("ppchg") != nullptr;
  const Yc &ch = *t.get(pp ? "ppchg" : "pctchg");
  const Yc *g = t.get("group"), *pl = t.get("perlbl"), *lo = t.get("pctlo"), *hi = t.get("pctup");
  const Yc &p = *t.get("prd");
  if (g) s.hd.push_back("group");
  for (const char *h : {"prd", "value", "previous", "change"}) s.hd.push_back(h);
  s.hd.push_back(pp ? "points" : "change %");
  if (intv && lo && hi) s.hd.push_back("interval");
  s.cols.assign(s.hd.size(), Strs());
  for (std::size_t i = 0; i < t.n; ++i) {
    std::size_t q = 0;
    if (g) {
      s.cols[q++].push_back(g->na[i] ? na : g->s[i]);
      s.gna.push_back(g->na[i]);
    }
    s.cols[q++].push_back(pl ? (pl->na[i] ? na : pl->s[i]) : (p.na[i] ? na : rfmt(p.d[i])));
    s.cols[q++].push_back(ynum(*t.get("value"), i, vd));
    s.cols[q++].push_back(ynum(*t.get("prev"), i, vd));
    s.cols[q++].push_back(ynum(*t.get("change"), i, vd));
    s.cols[q++].push_back(ypct(ch.na[i] ? NaN : ch.d[i], dg, pp));
    if (intv && lo && hi) {
      s.cols[q++].push_back(lo->na[i] || std::isnan(lo->d[i]) ? std::string(dash)
                                                               : "[" + ypct(lo->d[i], 0, false) + ", " +
                                                                     ypct(hi->na[i] ? NaN : hi->d[i], 0, false) + "]");
    }
  }
  return s;
}

// .ycls()
inline Strs ycls(const Yc &ch, int better) {
  Strs o(ch.n(), "flat");
  for (std::size_t i = 0; better != 0 && i < ch.n(); ++i) {
    const double v = ch.na[i] ? NaN : ch.d[i];
    if (v > 0) o[i] = better == 1 ? "good" : "bad";
    if (v < 0) o[i] = better == 1 ? "bad" : "good";
  }
  return o;
}

// .yhpct(): the cell of a change.
inline std::string yhpct(const std::string &v, double ch, const std::string &cl, double sc, bool bars) {
  const char *a = std::isnan(ch) ? "" : ch > 0 ? "\xe2\x96\xb2" : ch < 0 ? "\xe2\x96\xbc" : "\xe2\x80\x93";
  std::string bar;
  if (bars && std::isfinite(ch)) {
    char b[128];
    std::snprintf(b, sizeof b, "<span class=\"bar %s\" style=\"width:%.1f%%\"></span>", cl.c_str(),
                  std::min(100.0, 100 * std::fabs(ch) / sc));
    bar = b;
  }
  return "<td class=\"pct " + cl + "\"><span class=\"arw\">" + a + "</span><span class=\"val\">" + v + "</span>" +
         bar + "</td>";
}

// .yhrow(): row r with its note row.
inline std::string yhrow(std::size_t r, const Show &s, const Yc &ch, const Strs &cls, double sc, bool bars,
                         const Yc *fl) {
  std::string o = "<tr class=\"" + cls[r] + "\">";
  for (std::size_t j = 0; j < s.hd.size(); ++j) {
    const std::string &h = s.hd[j], &v = s.cols[j][r];
    if (h == "change %" || h == "points") {
      o += yhpct(v, ch.na[r] ? NaN : ch.d[r], cls[r], sc, bars);
      continue;
    }
    const char *k = h == "group" ? "grp" : h == "prd" ? "per" : h == "interval" ? "ci" : "num";
    o += std::string("<td class=\"") + k + "\">" + (h == "group" || h == "prd" ? esc(v) : v) + "</td>";
  }
  o += "</tr>";
  if (fl && !fl->na[r] && fl->s[r] != "no earlier period") {
    o += "<tr class=\"note\"><td colspan=\"" + std::to_string(s.hd.size()) + "\">note: " + esc(fl->s[r]) +
         "</td></tr>";
  }
  return o;
}

// .yhcard(): the card of row r of the summary.
inline std::string yhcard(const Yt &sm, std::size_t r) {
  auto num = [&](const char *k) { return sm.get(k)->na[r] ? NaN : sm.get(k)->d[r]; };
  auto txt = [&](const char *k) { return sm.get(k)->na[r] ? std::string("NA") : rfmt(sm.get(k)->d[r]); };
  auto cnt = [&](const char *k) {
    return sm.get(k)->na[r] ? std::string("NA") : std::to_string(static_cast<std::int64_t>(num(k)));
  };
  const Yc *g = sm.get("group");
  return "<div class=\"card\"><div class=\"lab\">" + esc(g ? g->s[r] : "overall") + "</div><div class=\"big\">" +
         ypct(num("totalpct"), 1, false) + "</div><div class=\"cap\">" + esc(txt("from") + " to " + txt("to")) +
         "</div><div class=\"cap\">" + esc("compound " + ypct(num("cagrpct"), 1, false) + " per period") +
         "</div><div class=\"cap\">" + esc(cnt("up") + " up, " + cnt("down") + " down, " + cnt("flat") + " unchanged") +
         "</div></div>";
}

// .yfirst(): the note on the first periods, or "".
inline std::string yfirst(const Yt &t, double lag) {
  const Yc *fl = t.get("flag");
  bool any = false;
  for (std::size_t i = 0; fl && i < t.n; ++i) any = any || (!fl->na[i] && fl->s[i] == "no earlier period");
  if (!any) return "";
  char b[48];
  std::snprintf(b, sizeof b, "%.7g periods", lag);
  return std::string("Note: the first ") + (lag == 1 ? "period" : b) + (t.get("group") ? " of each group" : "") +
         (lag == 1 ? " has" : " have") + " nothing to compare with, so " + (lag == 1 ? "it shows" : "they show") +
         " no change.";
}

// .ycss()
inline std::string ycss(const char *const *p) {
  return std::string(":root{--bad:") + p[1] + ";--good:" + p[2] + ";--flat:" + p[3] + ";--badbg:" + p[4] +
         ";--goodbg:" + p[5] +
         ";--ink:#1a1a1a;--dim:#5c5c5c;--rule:#e0ddd8;--ground:#faf9f7;"
         "--panel:#ffffff}"
         "@media (prefers-color-scheme:dark){:root:not([data-theme=light]){"
         "--ink:#eceae6;--dim:#a6a29b;--rule:#332f2b;--ground:#171513;"
         "--panel:#201d1a;--badbg:#331f1a;--goodbg:#142a26}}"
         "*{box-sizing:border-box}"
         "body{margin:0;background:var(--ground);color:var(--ink);"
         "font:14px/1.5 ui-sans-serif,-apple-system,'Segoe UI',Roboto,"
         "'Helvetica Neue',Arial,sans-serif}"
         "main{max-width:1000px;margin:0 auto;padding:2.5rem 1.25rem 4rem}"
         "h1{font-size:1.6rem;margin:0 0 .25rem;letter-spacing:-.01em}"
         "h2{font-size:1rem;text-transform:uppercase;letter-spacing:.08em;"
         "color:var(--dim);margin:2.5rem 0 .75rem;font-weight:600}"
         "p.sub{color:var(--dim);margin:0 0 1.75rem}"
         ".scroll{overflow-x:auto;border:1px solid var(--rule);border-radius:6px;"
         "background:var(--panel)}"
         "table{border-collapse:collapse;width:100%;min-width:520px}"
         "th{text-align:left;font-size:.72rem;text-transform:uppercase;"
         "letter-spacing:.06em;color:var(--dim);font-weight:600;"
         "padding:.65rem .75rem;border-bottom:1px solid var(--rule);"
         "white-space:nowrap}"
         "td{padding:.55rem .75rem;border-bottom:1px solid var(--rule);"
         "vertical-align:middle}"
         "tbody tr:last-child td{border-bottom:0}"
         "td.num,td.ci{text-align:right;font-variant-numeric:tabular-nums}"
         "td.per{font-variant-numeric:tabular-nums;white-space:nowrap}"
         "td.grp{font-weight:600}"
         "td.ci{color:var(--dim);font-size:.85rem;white-space:nowrap}"
         "td.pct{position:relative;text-align:right;white-space:nowrap;"
         "font-variant-numeric:tabular-nums}"
         "td.pct .arw{margin-right:.3rem;font-size:.8rem}"
         "td.pct.bad{color:var(--bad)}td.pct.good{color:var(--good)}"
         "td.pct.flat{color:var(--flat)}"
         ".bar{position:absolute;left:0;bottom:0;height:2px;opacity:.55}"
         ".bar.bad{background:var(--bad)}.bar.good{background:var(--good)}"
         ".bar.flat{background:var(--flat)}"
         "tr.bad td.pct{background:var(--badbg)}"
         "tr.good td.pct{background:var(--goodbg)}"
         "tr.note td{border-bottom:1px solid var(--rule);color:var(--dim);"
         "font-size:.8rem;padding:.2rem .75rem .5rem}"
         ".cards{display:grid;gap:.75rem;"
         "grid-template-columns:repeat(auto-fit,minmax(180px,1fr))}"
         ".card{border:1px solid var(--rule);border-radius:6px;padding:.9rem;"
         "background:var(--panel)}"
         ".card .lab{font-size:.72rem;text-transform:uppercase;"
         "letter-spacing:.06em;color:var(--dim);font-weight:600}"
         ".card .big{font-size:1.5rem;margin:.2rem 0 .3rem;"
         "font-variant-numeric:tabular-nums}"
         ".card .cap{color:var(--dim);font-size:.8rem}"
         "ul.notes{color:var(--dim);font-size:.85rem;padding-left:1.1rem}"
         "p.gen{color:var(--dim);font-size:.75rem;margin-top:2.5rem}";
}

// .y1file(): one path.
// Left out (not in R): the page is returned as text.
inline bool y1file(Ctx &c, std::string &f) {
  f.clear();
  if (!c.has("file")) return true;
  const Val &o = c["file"];
  const Val *s = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  if (s->t != Val::STR) return c.fail("`file` must be one file path.");
  f = s->s;
  return ckpath(c, f, "file");
}
// .yargs(): title, sub and notes.
inline bool yargs(Ctx &c, std::string &title, std::string &sub, Strs &notes) {
  CK(ytxt1(c, "title", "Year-over-year change", title));
  CK(ytxt1(c, "sub", "", sub));
  const Val &n = c["notes"];
  bool ok = n.t == Val::NUL || n.t == Val::STR || n.t == Val::ARR;
  if (n.t == Val::STR) notes.push_back(n.s);
  for (std::size_t i = 0; ok && n.t == Val::ARR && i < n.a.size(); ++i) {
    ok = n.a[i].t == Val::STR;
    if (ok) notes.push_back(n.a[i].s);
  }
  return ok || c.fail("`notes` must be text (a character vector), or NULL.");
}
// .ypick()
inline bool ypick(Ctx &c, const char *nm, const Strs &opt, const char *def, int &k) {
  std::string v = def;
  if (c.has(nm)) {
    const Val &o = c[nm];
    const Val *s = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
    v = s->t == Val::STR ? s->s : std::string("\x01");
  }
  for (std::size_t i = 0; i < opt.size(); ++i) {
    if (opt[i] == v) {
      k = static_cast<int>(i);
      return true;
    }
  }
  std::string l;
  for (std::size_t i = 0; i < opt.size(); ++i) l += (i ? ", \"" : "\"") + opt[i] + "\"";
  return c.fail("`%s` must be one of %s.", nm, l.c_str());
}

}  // namespace yh

// yoyhtml()
bool yoyhtml(Ctx &c, Val &out) {
  using namespace yh;
  Yt t;
  CK(ytab(c, t, true));
  std::string f, title, sub;
  CK(y1file(c, f));
  Strs notes;
  CK(yargs(c, title, sub, notes));
  int pal = 0, better = 0, dg = 1;
  CK(ypick(c, "pal", {"diverging", "safe", "mono"}, "diverging", pal));
  CK(ypick(c, "better", {"neutral", "higher", "lower"}, "neutral", better));
  CK(ydig(c, "digits", 1, dg, false));
  bool bars = true, intv = true;
  CK(yflag(c, "bars", true, bars));
  CK(yflag(c, "intv", true, intv));
  const Show s = yshow(t, dg, intv);
  const Yc &ch = *t.get(t.get("ppchg") ? "ppchg" : "pctchg");
  const Strs cls = ycls(ch, better);
  double sc = -HUGE_VAL;
  for (std::size_t i = 0; i < ch.n(); ++i) {
    if (!ch.na[i] && !std::isnan(ch.d[i])) sc = std::max(sc, std::fabs(ch.d[i]));
  }
  if (!std::isfinite(sc) || sc == 0) sc = 1;
  double lag = 1;
  CK(ylag(c, lag));
  Yt sm;
  CK(ysmry(c, t, lag, sm));
  Strs nt;
  if (!yfirst(t, lag).empty()) nt.push_back(yfirst(t, lag));
  nt.insert(nt.end(), notes.begin(), notes.end());
  std::string th;
  for (const std::string &h : s.hd) th += "<th>" + esc(h) + "</th>";
  std::string cards;
  for (std::size_t r = 0; r < sm.n; ++r) cards += yhcard(sm, r);
  Strs l = {"<!DOCTYPE html>", "<html lang=\"en\"><head>", "<meta charset=\"utf-8\">",
            "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
            "<title>" + esc(title) + "</title>", "<style>" + ycss(pals[pal]) + "</style>", "</head><body><main>",
            "<h1>" + esc(title) + "</h1>", "<p class=\"sub\">" + esc(sub) + "</p>",
            "<div class=\"scroll\"><table><thead><tr>", th, "</tr></thead><tbody>"};
  for (std::size_t r = 0; r < t.n; ++r) l.push_back(yhrow(r, s, ch, cls, sc, bars, t.get("flag")));
  l.push_back("</tbody></table></div>");
  l.push_back("<h2>Across the span</h2><div class=\"cards\">" + cards + "</div>");
  if (!nt.empty()) {
    l.push_back("<ul class=\"notes\">");
    for (const std::string &n : nt) l.push_back("<li>" + esc(n) + "</li>");
    l.push_back("</ul>");
  }
  l.push_back("<p class=\"gen\">" + esc(std::string("crimlab ") + CLAB_VER) + "</p>");
  l.push_back("</main></body></html>");
  std::string txt;
  for (const std::string &x : l) txt += x + "\n";
  return ydone(c, txt, f, out);
}
