// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// yoypdf() of R/yoyh.R for clab_jrun(): the PDF R writes itself (no
// device), byte for byte: the layout of .yplay(), the content streams of
// .ypstr() and .yprows() in Helvetica with its widths, and the file of
// .ypdfw(), with R's checks in their order and words. Without `file`
// (which R requires) the PDF, ASCII throughout, is returned as text. R
// keeps its own writer.

namespace yp {

using yh::Strs;

// Helvetica widths (1/1000 em) of the printable ASCII characters (.yhw).
const int hw[95] = {278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278, 278, 556, 556, 556,
                    556, 556, 556, 556, 556, 556, 556, 278, 278, 584, 584, 584, 556, 1015, 667, 667, 722, 722, 667,
                    611, 778, 722, 278, 500, 667, 556, 833, 722, 778, 667, 778, 722, 667, 611, 722, 667, 944, 667,
                    667, 611, 278, 278, 278, 469, 556, 333, 556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500,
                    222, 833, 556, 556, 556, 556, 333, 500, 278, 556, 500, 722, 500, 500, 500, 334, 260, 334, 584};

// .yasc(): dashes and arrows as ASCII, then toascii().
inline std::string yasc(std::string s) {
  static const char *const rp[][2] = {{"\xe2\x80\x94", "-"}, {"\xe2\x80\x93", "-"}, {"\xe2\x96\xb2", "+"},
                                      {"\xe2\x96\xbc", "-"}, {"\xc2\xb7", "|"},     {"\xe2\x94\x80", "-"},
                                      {"\xe2\x86\xb3", ">"}};
  for (const auto &p : rp) {
    const std::string a = p[0];
    for (std::size_t k = s.find(a); k != std::string::npos; k = s.find(a, k + 1)) s.replace(k, a.size(), p[1]);
  }
  return ascii1(s);
}
// .ywid()
inline double ywid(const std::string &s, double size) {
  double w = 0;
  for (unsigned char ch : s) {
    if (ch >= 32 && ch <= 126) w += hw[ch - 32];
  }
  return w * size / 1000;
}
// .ypdfs()
inline std::string ypdfs(const std::string &s) {
  std::string o;
  for (char ch : s) {
    if (ch == '\\' || ch == '(' || ch == ')') o += '\\';
    o += ch;
  }
  return o;
}
// .yrgb()
inline std::string yrgb(const char *hex) {
  std::string o;
  for (int i = 0; i < 3; ++i) {
    const double v = static_cast<double>(std::strtol(std::string(hex + 1 + 2 * i, 2).c_str(), nullptr, 16)) / 255;
    char b[32];
    std::snprintf(b, sizeof b, "%s%.3f", i ? " " : "", v);
    o += b;
  }
  return o;
}
// .yptxt()
inline std::string yptxt(double x0, double y0, const std::string &s, double size, const char *font = "F1",
                         const std::string &rgb = "0 0 0", bool right = false) {
  if (right) x0 = x0 - ywid(s, size);
  char b[96];
  std::snprintf(b, sizeof b, "BT /%s %g Tf %s rg %.2f %.2f Td (", font, size, rgb.c_str(), x0, y0);
  return b + ypdfs(s) + ") Tj ET";
}
inline std::string fmt(const char *f, double a, double b, double c, double d) {
  char o[160];
  std::snprintf(o, sizeof o, f, a, b, c, d);
  return o;
}

// The layout of .yplay() and what .ypstr() draws with it.
struct Lay {
  Strs hdr, nt;
  std::vector<Strs> body;
  Strs cls;
  double W = 0, H = 0, m = 0, rowh = 0;
  std::size_t nc = 0, numfrom = 0;
  std::vector<double> xs;
  std::vector<std::vector<std::size_t>> pages;
  const char *const *p = nullptr;
  std::string dim, title, sub;
};
inline void yplay(Lay &l, const Yt &t, const yh::Show &s, double width, double height) {
  for (const std::string &h : s.hd) l.hdr.push_back(yasc(h));
  l.W = 72 * width;
  l.H = 72 * height;
  l.m = 0.6 * 72 / 6;
  l.nc = l.hdr.size();
  l.numfrom = (t.get("group") ? 1 : 0) + 2;
  std::vector<double> w(l.nc, 1.0);
  for (std::size_t j = 0; j + 1 < l.numfrom && j < l.nc; ++j) w[j] = 1.3;
  long double sw = 0.0L, cs = 0.0L;
  for (double v : w) sw += v;
  const double tot = static_cast<double>(sw);
  l.xs.push_back(l.m);
  for (double v : w) {
    cs += v / tot;
    l.xs.push_back(l.m + (l.W - 2 * l.m) * static_cast<double>(cs));
  }
  const double pf = std::floor((0.82 * height - 0.22 * 1.1) / 0.22);
  const std::size_t per = static_cast<std::size_t>(std::max(5.0, pf));
  for (std::size_t i = 0; i < t.n; ++i) {
    if (i % per == 0) l.pages.emplace_back();
    l.pages.back().push_back(i);
  }
  if (l.pages.empty()) l.pages.emplace_back();
  l.rowh = 0.22 * 72;
}

// .ypstr(): the content stream of page k (from 0).
inline std::string ypstr(const Lay &l, std::size_t k) {
  const double m = l.m, W = l.W, H = l.H, top = H - m - 0.12 * H;
  const std::size_t np = l.pages.size();
  auto rule = [&](double y) { return l.dim + fmt(" RG 0.7 w %.2f %.2f m %.2f %.2f l S", m, y, W - m, y); };
  Strs o = {yptxt(m, H - m - 18, l.title, 18, "F2"), yptxt(m, H - m - 32, l.sub, 9, "F1", l.dim), rule(top + 10)};
  for (std::size_t j = 0; j < l.nc; ++j) {
    const bool r = j + 1 >= l.numfrom && l.hdr[j] != "note";
    o.push_back(yptxt(r ? l.xs[j + 1] - 3 : l.xs[j] + 3, top, l.hdr[j], 8, "F2", l.dim, r));
  }
  o.push_back(rule(top - 5));
  double y = top - l.rowh * 1.1;
  for (std::size_t i : l.pages[k]) {
    for (std::size_t j = 0; j < l.nc; ++j) {
      const bool r = j + 1 >= l.numfrom;
      const bool pc = l.hdr[j] == "change %" || l.hdr[j] == "points";
      const std::string col = pc ? yrgb(l.p[l.cls[i] == "bad" ? 1 : l.cls[i] == "good" ? 2 : 3]) : "0 0 0";
      o.push_back(yptxt(r ? l.xs[j + 1] - 3 : l.xs[j] + 3, y, l.body[j][i], 9, "F1", col, r));
    }
    y = y - l.rowh;
  }
  if (np > 1) {
    o.push_back(yptxt(W - m, m, "page " + std::to_string(k + 1) + " of " + std::to_string(np), 7, "F1", l.dim, true));
  }
  if (k + 1 == np && !l.nt.empty()) {
    double yy = m + 10 + 10 * static_cast<double>(l.nt.size());
    for (const std::string &q : l.nt) {
      o.push_back(yptxt(m, yy, q, 7.5, "F1", l.dim));
      yy = yy - 10;
    }
  }
  std::string s;
  for (std::size_t i = 0; i < o.size(); ++i) s += (i ? "\n" : "") + o[i];
  return s;
}

// .ypdfw(): catalog, pages, two fonts, a page and a stream per page, the
// information and the cross-reference table.
inline std::string ypdfw(const Strs &st, double W, double H, const std::string &title) {
  const std::size_t np = st.size();
  std::string kids;
  for (std::size_t k = 0; k < np; ++k) kids += (k ? " " : "") + std::to_string(5 + 2 * k) + " 0 R";
  Strs obj = {"<< /Type /Catalog /Pages 2 0 R >>",
              "<< /Type /Pages /Kids [" + kids + "] /Count " + std::to_string(np) + " >>",
              "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>",
              "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>"};
  for (std::size_t k = 0; k < np; ++k) {
    char b[256];
    std::snprintf(b, sizeof b,
                  "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 %.2f %.2f] /Resources << /Font << /F1 3 0 R /F2 4 0 R "
                  ">> >> /Contents %d 0 R >>",
                  W, H, static_cast<int>(6 + 2 * k));
    obj.push_back(b);
    obj.push_back("<< /Length " + std::to_string(st[k].size()) + " >>\nstream\n" + st[k] + "\nendstream");
  }
  obj.push_back("<< /Title (" + ypdfs(title) + ") /Producer (crimlab) >>");
  const std::string head = "%PDF-1.4\n";
  std::string body, xref = "xref\n0 " + std::to_string(obj.size() + 1) + "\n0000000000 65535 f \n";
  for (std::size_t i = 0; i < obj.size(); ++i) {
    char b[32];
    std::snprintf(b, sizeof b, "%010d 00000 n \n", static_cast<int>(head.size() + body.size()));
    xref += b;
    body += std::to_string(i + 1) + " 0 obj\n" + obj[i] + "\nendobj\n";
  }
  return head + body + xref + "trailer\n<< /Size " + std::to_string(obj.size() + 1) + " /Root 1 0 R /Info " +
         std::to_string(obj.size()) + " 0 R >>\nstartxref\n" + std::to_string(head.size() + body.size()) +
         "\n%%EOF\n";
}

// .ypsize()
inline bool ypsize(Ctx &c, double &w, double &h) {
  w = 11;
  h = 8.5;
  for (const char *nm : {"width", "height"}) {
    double &v = nm[0] == 'w' ? w : h;
    if (!c.has(nm)) continue;
    const std::vector<double> x = ax::numv(c[nm]) ? ax::nums(c[nm]) : std::vector<double>();
    if (x.size() != 1 || !std::isfinite(x[0]) || x[0] < 2 || x[0] > 200) {
      return c.fail("`%s` must be one number from 2 to 200.", nm);
    }
    v = x[0];
  }
  return true;
}

}  // namespace yp

// yoypdf()
bool yoypdf(Ctx &c, Val &out) {
  using namespace yh;
  Yt t;
  CK(ytab(c, t, true));
  std::string f, title, sub;
  CK(y1file(c, f));
  Strs notes;
  CK(yargs(c, title, sub, notes));
  double lag = 1;
  CK(ylag(c, lag));
  int pal = 0, better = 0, dg = 1;
  CK(ypick(c, "pal", {"diverging", "safe", "mono"}, "diverging", pal));
  CK(ypick(c, "better", {"neutral", "higher", "lower"}, "neutral", better));
  CK(ydig(c, "digits", 1, dg, false));
  bool intv = true;
  CK(yflag(c, "intv", true, intv));
  double w = 11, h = 8.5;
  CK(yp::ypsize(c, w, h));
  const Show s = yshow(t, dg, intv, "NA");
  yp::Lay l;
  yp::yplay(l, t, s, w, h);
  for (const Strs &col : s.cols) {
    Strs b;
    for (const std::string &v : col) b.push_back(yp::yasc(v));
    l.body.push_back(b);
  }
  l.cls = ycls(*t.get(t.get("ppchg") ? "ppchg" : "pctchg"), better);
  if (!yfirst(t, lag).empty()) l.nt.push_back(yp::yasc(yfirst(t, lag)));
  for (const std::string &n : notes) l.nt.push_back(yp::yasc(n));
  l.p = pals[pal];
  l.dim = yp::yrgb(l.p[3]);
  l.title = yp::yasc(title);
  l.sub = yp::yasc(sub);
  Strs st;
  for (std::size_t k = 0; k < l.pages.size(); ++k) st.push_back(yp::ypstr(l, k));
  return ydone(c, yp::ypdfw(st, l.W, l.H, l.title), f, out);
}
