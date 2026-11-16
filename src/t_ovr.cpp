// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for psc.cpp, sfx.cpp and ovr.cpp: the credit at one
// and at one and a half days a day by hand, the sentence left; draws of
// synfix() against R's (set.seed(1); sample.int(3, 4, TRUE);
// sample.int(1, 4, TRUE); sample.int(7, 4, TRUE), less one: 0 2 0 1,
// 0 0 0 0, 1 2 2 0) and the presets; the text of a page, dates, the
// schema preset, fields, subject officials, keywords and checks of
// oversight reports; and the failure paths. Part of the self-test
// (selftest.cpp).
#include <cmath>
#include <cstdint>
#include <cstring>

#include "crimlab.h"
#include "tchk.h"

namespace {

void psc(tk &tt) {
  // 10 days at 1.5 is 15, 45 days at 1.5 is 67.5; a sentence of 60 days
  // less 15 leaves 45, less 67.5 leaves nothing
  const double d[2] = {10, 45}, r[1] = {1.5}, s[1] = {60};
  double c[2], m[2];
  tk_ieq(tt, "psc.rc", clab_pscred(d, 2, r, 1, s, 1, c, m), CLAB_OK);
  tk_near(tt, "psc.c0", c[0], 15, 0);
  tk_near(tt, "psc.c1", c[1], 67.5, 0);
  tk_near(tt, "psc.m0", m[0], 45, 0);
  tk_near(tt, "psc.m1", m[1], 0, 0);
  const double r2[2] = {1, 1.5};
  clab_pscred(d, 2, r2, 2, nullptr, 0, c, m);
  tk_near(tt, "psc.r1", c[0], 10, 0);
  const double big[1] = {1.6}, neg[1] = {-1};
  tk_ieq(tt, "psc.ratio", clab_pscred(d, 2, big, 1, nullptr, 0, c, m), CLAB_EDOM);
  tk_ieq(tt, "psc.days", clab_pscred(neg, 1, r, 1, nullptr, 0, c, m), CLAB_EDOM);
  tk_ieq(tt, "psc.len", clab_pscred(d, 2, r2, 2, s, 3, c, m), CLAB_EARG);
}

void sfx(tk &tt) {
  // three columns of four rows; the stream runs on from one column to the
  // next
  const double dn[3] = {3, 1, 7};
  int c[12];
  tk_ieq(tt, "sfx.rc", clab_synfix(dn, 3, 4, 1, c), CLAB_OK);
  tk_ieq(tt, "sfx.c0", c[0] * 1000 + c[1] * 100 + c[2] * 10 + c[3], 201);
  tk_ieq(tt, "sfx.c1", c[4] + c[5] + c[6] + c[7], 0);
  tk_ieq(tt, "sfx.c2", c[8] * 1000 + c[9] * 100 + c[10] * 10 + c[11], 1220);
  const char *k, *nm, *lv;
  int ty, lo, hi;
  tk_ieq(tt, "sfx.col", clab_sfxcol(0, &k, &nm, &ty, &lo, &hi, &lv), CLAB_OK);
  tk_ieq(tt, "sfx.rng", ty * 10000 + lo, 12018);
  tk_ieq(tt, "sfx.end", clab_sfxcol(-1, &k, &nm, &ty, &lo, &hi, &lv), CLAB_EARG);
  const double bad[1] = {0};
  tk_ieq(tt, "sfx.dn", clab_synfix(bad, 1, 4, 1, c), CLAB_EARG);
}

void txt(tk &tt) {
  // two paragraphs and an entity: "a", an empty line, "b&c"
  const char *h = "<p>a</p>\n<p>b&amp;c</p>";
  char o[64];
  int64_t len = 0;
  tk_ieq(tt, "ovr.txt.rc", clab_ovrtxt(h, static_cast<int64_t>(std::strlen(h)), -1, o, &len), CLAB_OK);
  tk_ieq(tt, "ovr.txt.len", static_cast<int>(len), 6);
  tk_true(tt, "ovr.txt", std::memcmp(o, "a\n\nb&c", 6) == 0);
  // as plain text the tags stay
  clab_ovrtxt(h, static_cast<int64_t>(std::strlen(h)), 0, o, &len);
  tk_ieq(tt, "ovr.txt.plain", static_cast<int>(len), static_cast<int>(std::strlen(h)));
}

void dte(tk &tt) {
  // the impossible day passed over, the next date read
  const char *a = "on February 30, 2023 or 5th of March, 2023";
  char o[16];
  tk_ieq(tt, "ovr.dte.rc", clab_ovrdte(a, static_cast<int64_t>(std::strlen(a)), o), CLAB_OK);
  tk_true(tt, "ovr.dte", std::strcmp(o, "2023-03-05") == 0);
  clab_ovrdte("May 2023", 8, o);
  tk_ieq(tt, "ovr.dte.none", o[0], 0);
}

void sch(tk &tt) {
  const char *v[6];
  tk_ieq(tt, "ovr.sch.rc", clab_ovrsch("siu", 0, &v[0], &v[1], &v[2], &v[3], &v[4], &v[5]), CLAB_OK);
  tk_true(tt, "ovr.sch.rule", std::strcmp(v[1], "after") == 0);
  tk_ieq(tt, "ovr.sch.end", clab_ovrsch("siu", 17, &v[0], &v[1], &v[2], &v[3], &v[4], &v[5]), CLAB_EARG);
  tk_ieq(tt, "ovr.sch.none", clab_ovrsch("x", 0, &v[0], &v[1], &v[2], &v[3], &v[4], &v[5]), CLAB_EARG);
}

void fld(tk &tt) {
  // a count after its label, the highest numbered tag, a flag
  const char *t = "Number of X assigned: 4\nWO #2 and WO\n#5\nno grounds";
  const int64_t n = static_cast<int64_t>(std::strlen(t));
  char o[256];
  int64_t len = 0;
  tk_ieq(tt, "ovr.fld.rc", clab_ovrfld(t, n, 1, "Number of X assigned", "", "", o, &len), CLAB_OK);
  tk_true(tt, "ovr.fld.num", len == 1 && o[0] == '4');
  clab_ovrfld(t, n, 3, "WO", "", "", o, &len);
  tk_true(tt, "ovr.fld.tag", len == 1 && o[0] == '5');
  clab_ovrfld(t, n, 8, "grounds to", "no grounds", "", o, &len);
  tk_true(tt, "ovr.fld.flag", len == 5 && std::memcmp(o, "FALSE", 5) == 0);
  clab_ovrfld(t, n, 3, "SO", "", "", o, &len);
  tk_true(tt, "ovr.fld.none", len == -1);
  tk_ieq(tt, "ovr.fld.rule", clab_ovrfld(t, n, 11, "", "", "", o, &len), CLAB_EARG);
}

void soc(tk &tt) {
  const char *t = "The unit designated two subject officials.";
  int c = 0, r = 0;
  char w[256];
  tk_ieq(tt, "ovr.soc.rc",
         clab_ovrsoc(t, static_cast<int64_t>(std::strlen(t)), "SO", "subject official", "WO", "Subject Officials", &c,
                     &r, w, sizeof w),
         CLAB_OK);
  tk_ieq(tt, "ovr.soc.n", c * 10 + r, 23);
  const char *u = "nothing here";
  clab_ovrsoc(u, 12, "SO", "subject official", "WO", "Subject Officials", &c, &r, w, sizeof w);
  tk_ieq(tt, "ovr.soc.na", c * 10 + r, -10 + 6);
}

void kwd(tk &tt) {
  // "Black" twice as a word, not in "Blackburn"; a stem with '*'
  const char *t = "Black, black and Blackburn; schizophrenia";
  const int64_t n = static_cast<int64_t>(std::strlen(t));
  int64_t k = 0;
  tk_ieq(tt, "ovr.kwd.rc", clab_ovrkwd(t, n, "Black", &k), CLAB_OK);
  tk_ieq(tt, "ovr.kwd.word", static_cast<int>(k), 2);
  clab_ovrkwd(t, n, "schizophren*", &k);
  tk_ieq(tt, "ovr.kwd.stem", static_cast<int>(k), 1);
  clab_ovrkwd(t, n, "schizophren", &k);
  tk_ieq(tt, "ovr.kwd.whole", static_cast<int>(k), 0);
  const char *g, *w;
  tk_ieq(tt, "ovr.kwd.voc", clab_ovrvoc("mhr", 0, &g, &w), CLAB_OK);
  tk_ieq(tt, "ovr.kwd.none", clab_ovrvoc("x", 0, &g, &w), CLAB_EARG);
}

void chk(tk &tt) {
  int ok = 0;
  tk_ieq(tt, "ovr.chk.rc", clab_ovrchk("2024-02-29", 10, 0, "", &ok), CLAB_OK);
  tk_ieq(tt, "ovr.chk.leap", ok, 1);
  clab_ovrchk("2023-02-29", 10, 0, "", &ok);
  tk_ieq(tt, "ovr.chk.date", ok, 0);
  clab_ovrchk("12", 2, 1, "", &ok);
  tk_ieq(tt, "ovr.chk.count", ok, 1);
  clab_ovrchk("Yes", 3, 2, "", &ok);
  tk_ieq(tt, "ovr.chk.flag", ok, 1);
  clab_ovrchk("23-OFD-001", 10, 3, "99-AA-999|99-AAA-999", &ok);
  tk_ieq(tt, "ovr.chk.shape", ok, 1);
  tk_ieq(tt, "ovr.chk.kind", clab_ovrchk("x", 1, 4, "", &ok), CLAB_EARG);
}

}  // namespace

void t_ovr(tk &t, int big) {
  (void)big;
  psc(t);
  sfx(t);
  txt(t);
  dte(t);
  sch(t);
  fld(t);
  soc(t);
  kwd(t);
  chk(t);
}
