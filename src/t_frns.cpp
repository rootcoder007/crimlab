// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Self-test of the forensic kernels: match probabilities, the cost of
// likelihood ratios and its calibrated minimum, mixtures, kinship and
// paternity indices.
#include <cmath>
#include <cstring>

#include "crimlab.h"
#include "tchk.h"

namespace {

void gtmp(tk &t) {
  const double p[2] = {0.1, 0.2}, q[2] = {0.3, 0};
  const int hom[2] = {0, 1};
  double loc[500], l10;
  tk_ieq(t, "frns3", clab_gtmp(p, q, hom, 2, 0, loc, &l10), CLAB_OK);
  const double e3[2] = {2 * 0.1 * 0.3, 0.2 * 0.2};
  tk_veq(t, "frns3", loc, e3, 2, 1e-10);
  tk_eq(t, "frns4", l10, std::log10(0.06 * 0.04), 1e-10);
  const double h = 0.5, z = 0;
  const int one = 1, nil = 0;
  clab_gtmp(&h, &z, &one, 1, 0, loc, &l10);
  tk_eq(t, "frns7", loc[0], 0.25, tk_tol);
  // Balding-Nichols by hand
  const double th = 0.03, pp = 0.1, qq = 0.3;
  const double d = (1 + th) * (1 + 2 * th);
  const double eh = (2 * th + (1 - th) * pp) * (3 * th + (1 - th) * pp) / d;
  const double et = 2 * (th + (1 - th) * pp) * (th + (1 - th) * qq) / d;
  const double p2[2] = {pp, pp}, q2[2] = {0, qq};
  const int h2[2] = {1, 0};
  clab_gtmp(p2, q2, h2, 2, th, loc, &l10);
  const double e18[2] = {eh, et};
  tk_veq(t, "frns18", loc, e18, 2, 1e-14);
  // shared ancestry makes a match less surprising
  tk_gt(t, "frns20", eh, pp * pp);
  tk_gt(t, "frns21", et, 2 * pp * qq);
  // 400 and 500 loci would underflow as a product; the logarithm does not
  double pv[500], qv[500];
  int hv[500];
  for (int i = 0; i < 500; ++i) {
    pv[i] = 1e-3;
    qv[i] = 0;
    hv[i] = 1;
  }
  clab_gtmp(pv, qv, hv, 400, 0, loc, &l10);
  tk_eq(t, "frns24", l10, -2400, 1e-10);
  for (int i = 0; i < 500; ++i) pv[i] = 1e-2;
  clab_gtmp(pv, qv, hv, 500, 0, loc, &l10);
  tk_eq(t, "edge16", l10, -2000, 1e-10);
  // coancestry can make a common heterozygote rarer
  clab_gtmp(&h, &h, &nil, 1, 0.1, loc, &l10);
  tk_lt(t, "frns2_2", loc[0], 0.5);
  clab_gtmp(&h, &z, &one, 1, 0.1, loc, &l10);
  tk_ge(t, "frns2_3", loc[0], 0.25);
  // a heterozygote of the only two alleles at a locus is allowed
  tk_ieq(t, "edge36", clab_gtmp(&h, &h, &nil, 1, 0, loc, &l10), CLAB_OK);
  tk_eq(t, "edge36", loc[0], 0.5, tk_tol);
}

double cllr(const double *s, int ns, const double *d, int nd) {
  double c = NAN;
  clab_cllr(s, ns, d, nd, &c);
  return c;
}

void cost(tk &t) {
  const double one[3] = {1, 1, 1};
  tk_eq(t, "frns28", cllr(one, 3, one, 2), 1, 1e-10);
  const double sg[2] = {1e12, 1e15}, dg[2] = {1e-12, 1e-15};
  tk_lt(t, "frns29", cllr(sg, 2, dg, 2), 1e-10);
  const double sb[2] = {0.01, 0.1}, db[2] = {100, 10};
  tk_gt(t, "frns30", cllr(sb, 2, db, 2), 1);
  const double s[3] = {50, 200, 8}, d[3] = {0.01, 0.2, 1.5};
  double ms = 0, md = 0;
  for (int i = 0; i < 3; ++i) {
    ms += std::log2(1 + 1 / s[i]) / 3;
    md += std::log2(1 + d[i]) / 3;
  }
  tk_eq(t, "frns33", cllr(s, 3, d, 3), 0.5 * (ms + md), 1e-14);
  // extreme ratios cost a finite number of bits
  const double e300 = 1e-300, b300 = 1e300, u = 1;
  tk_eq(t, "edge23", cllr(&e300, 1, &e300, 1), 0.5 * std::log2(1 + 1e300),
        1e-14);
  tk_eq(t, "edge24", cllr(&u, 1, &b300, 1), 0.5 * (1 + std::log2(1 + 1e300)),
        1e-14);
  // a subnormal same-source ratio: 1 / same overflows a double
  const double sub = 1e-310, want = 0.5 * (310 * std::log2(10.0) + 1);
  tk_eq(t, "edge63", cllr(&sub, 1, &u, 1), want, 1e-12);
  tk_eq(t, "edge64", cllr(&sub, 1, &u, 1), want, 1e-12);
  // both branches of the rewritten term agree where they meet
  const double lo = 1 - 1e-12, hi = 1 + 1e-12;
  tk_eq(t, "edge69", cllr(&lo, 1, &u, 1), cllr(&hi, 1, &u, 1), 1e-12);
}

// The calibrated minimum cost: same- and different-source ratios pooled and
// sorted, stably, as the R layer passes them to the kernel.
double cmin(const double *s, int ns, const double *d, int nd) {
  double lr[16], post[16], work[32], out = NAN;
  int is[16];
  const int n = ns + nd;
  for (int i = 0; i < n; ++i) {
    const double x = i < ns ? s[i] : d[i - ns];
    int j = i;
    for (; j > 0 && lr[j - 1] > x; --j) {
      lr[j] = lr[j - 1];
      is[j] = is[j - 1];
    }
    lr[j] = x;
    is[j] = i < ns;
  }
  clab_lrcal(lr, is, n, post, work, &out);
  return out;
}

void cal(tk &t) {
  double s[4] = {50, 200, 8}, d[4] = {0.01, 0.2, 1.5};
  tk_eq(t, "frns2_112", cllr(s, 3, d, 3), cllr(s, 3, d, 3), 1e-10);
  tk_eq(t, "frns2_113", cmin(s, 3, d, 3), 0, 1e-10);
  // scaling every ratio changes the cost, not the minimum
  const double s2[4] = {3, 0.7, 20, 5}, d2[4] = {0.2, 1, 4, 0.05};
  for (int i = 0; i < 4; ++i) {
    s[i] = s2[i] * 1000;
    d[i] = d2[i] * 1000;
  }
  const double am = cmin(s2, 4, d2, 4), ac = cllr(s2, 4, d2, 4);
  tk_eq(t, "frns2_120", am, cmin(s, 4, d, 4), 1e-10);
  tk_gt(t, "frns2_121", cllr(s, 4, d, 4), ac);
  tk_le(t, "frns2_122", am, ac);
  tk_gt(t, "frns2_123", am, 0);
  // by hand: order 0.05d 0.2d 0.7s 1d 3s 4d 5s 20s; blocks d d | s d | s d |
  // s s with same-source share 0, 0.5, 0.5, 1
  tk_eq(t, "frns2_126", am, 0.5 * (2 * std::log2(2.0) / 4 +
                                   2 * std::log2(2.0) / 4), 1e-10);
  // ties across the two sets are pooled
  const double two[2] = {2, 2}, rs[2] = {0.1, 0.2}, rd[2] = {5, 9};
  tk_eq(t, "frns2_136", cmin(two, 2, two, 2), 1, 1e-10);
  // an uninformative reversed system
  tk_eq(t, "frns2_138", cmin(rs, 2, rd, 2), 1, 1e-10);
  // by hand, two same-source against three different-source ratios: order
  // 1s 2d 3d 4s 5d pools to s d d (share 1/3) | s d (share 1/2)
  const double lr[5] = {1, 2, 3, 4, 5}, e[5] = {1.0 / 3, 1.0 / 3, 1.0 / 3,
                                                0.5, 0.5};
  const int is[5] = {1, 0, 0, 1, 0};
  double post[5], work[10], mn = NAN;
  clab_lrcal(lr, is, 5, post, work, &mn);
  tk_veq(t, "lrcal.post", post, e, 5, 1e-15);
  tk_eq(t, "lrcal.min", mn,
        0.5 * ((std::log2(7.0 / 3) + std::log2(5.0 / 3)) / 2 +
               (2 * std::log2(7.0 / 4) + std::log2(5.0 / 2)) / 3), 1e-14);
}

double mixp(const double *p, const int *need, int x) {
  double out = NAN;
  clab_mixp(p, need, 4, x, &out);
  return out;
}

void mix(tk &t) {
  const double p[4] = {0.1, 0.2, 0.05, 0.15};
  const int none[4] = {0, 0, 0, 0}, cd[4] = {0, 0, 1, 1}, all[4] = {1, 1, 1, 1},
            abc[4] = {1, 1, 1, 0}, c[4] = {0, 0, 1, 0}, a[4] = {1, 0, 0, 0};
  // mixlr(p, character(0), 0, c("c", "d"), 1): the numerator and
  // denominator; their ratio is formed in R and checked there
  tk_eq(t, "frns2_78", mixp(p, none, 0), 1, 1e-10);
  tk_eq(t, "frns2_78", mixp(p, cd, 1), 2 * 0.05 * 0.15, 1e-10);
  // nothing to explain, one unknown: random man not excluded
  tk_eq(t, "frns2_81", mixp(p, none, 1), 0.5 * 0.5, 1e-10);
  // two unknowns carry all four alleles: 4! ordered assignments
  tk_eq(t, "frns2_84", mixp(p, all, 2), 24 * 0.1 * 0.2 * 0.05 * 0.15, 1e-10);
  // three alleles cannot sit in one person
  tk_eq(t, "frns2_87", mixp(p, abc, 1), 0, 1e-10);
  const double s = 0.5;
  tk_eq(t, "frns2_91", mixp(p, c, 2),
        std::pow(s, 4) - std::pow(s - 0.05, 4), 1e-10);
  tk_eq(t, "frns2_93", mixp(p, a, 0), 0, 1e-10);
}

double kin(const int *g1, const int *g2, const double *f, double k0,
           double k1, double k2) {
  const double k[3] = {k0, k1, k2};
  double idx = NAN;
  clab_kinidx(g1, g2, f, 1, k, &idx);
  return idx;
}

void kinship(tk &t) {
  const double p = 0.1, q = 0.2, pq[2] = {p, q}, pp[2] = {p, p};
  const int g12[2] = {1, 2}, g11[2] = {1, 1}, g13[2] = {1, 3},
            g34[2] = {3, 4}, g21[2] = {2, 1}, g23[2] = {2, 3};
  tk_eq(t, "frns2_10", kin(g12, g12, pq, 0, 1, 0), (p + q) / (4 * p * q),
        1e-10);
  tk_eq(t, "frns2_12", kin(g11, g11, pp, 0, 1, 0), 1 / p, 1e-10);
  tk_eq(t, "frns2_13", kin(g11, g12, pp, 0, 1, 0), 1 / (2 * p), 1e-10);
  tk_eq(t, "frns2_15", kin(g12, g13, pq, 0, 1, 0), 1 / (4 * p), 1e-10);
  tk_eq(t, "frns2_17", kin(g12, g34, pq, 0, 1, 0), 0, 1e-10);
  tk_eq(t, "frns2_18", kin(g12, g34, pq, 0.25, 0.5, 0.25), 0.25, tk_tol);
  tk_eq(t, "frns2_19", kin(g12, g21, pq, 0.25, 0.5, 0.25),
        0.25 + 0.5 * (p + q) / (4 * p * q) + 0.25 / (2 * p * q), 1e-10);
  tk_eq(t, "frns2_22", kin(g11, g11, pp, 0, 0, 1), 1 / (p * p), 1e-10);
  // alleles a, b, c coded 1, 2, 3
  tk_eq(t, "frns2_24", kin(g12, g23, pq, 0.5, 0.5, 0), 0.5 + 0.5 / (4 * q),
        1e-10);
  tk_eq(t, "frns2_26", kin(g12, g12, pq, 1, 0, 0), 1, 1e-10);
  // no mother: the parent-child kinship index
  tk_eq(t, "frns2_51", kin(g12, g12, pq, 0, 1, 0), (p + q) / (4 * p * q),
        1e-10);
  // 500 loci each worth about 1e-3: every index finite
  int g[1000];
  double f[1000], idx[500];
  for (int i = 0; i < 1000; ++i) {
    g[i] = i / 2;
    f[i] = 1e-3;
  }
  const double k[3] = {0.9999, 0.0001, 0};
  clab_kinidx(g, g, f, 500, k, idx);
  bool fin = true;
  for (double x : idx) fin = fin && std::isfinite(x);
  tk_true(t, "edge13", fin);
}

double pat(const int *c, const int *m, const int *af, const double *f) {
  double idx = 0;
  clab_patidx(c, m, af, f, 1, &idx);
  return idx;
}

// The chance that a mother m and a father a have a child ch: the share of
// the four allele pairs that make the child's genotype.
double trans(const int *m, const int *a, const int *ch) {
  const int lo = ch[0] < ch[1] ? ch[0] : ch[1];
  const int hi = ch[0] < ch[1] ? ch[1] : ch[0];
  double k = 0;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      const int x = m[i], y = a[j];
      k += ((x < y ? x : y) == lo && (x < y ? y : x) == hi) ? 0.25 : 0;
    }
  }
  return k;
}

void paternity(tk &t) {
  const double p = 0.1, q = 0.2, pq[2] = {p, q}, pp[2] = {p, p};
  const int g12[2] = {1, 2}, g11[2] = {1, 1}, g22[2] = {2, 2},
            g23[2] = {2, 3}, g13[2] = {1, 3}, g33[2] = {3, 3};
  tk_eq(t, "frns2_39", pat(g12, g11, g22, pq), 1 / q, 1e-10);
  tk_eq(t, "frns2_41", pat(g12, g11, g23, pq), 1 / (2 * q), 1e-10);
  tk_eq(t, "frns2_43", pat(g12, g12, g12, pq), 1 / (p + q), 1e-10);
  tk_eq(t, "frns2_45", pat(g11, g12, g13, pp), 1 / (2 * p), 1e-10);
  tk_eq(t, "frns2_47", pat(g12, g11, g33, pq), 0, 1e-10);
  tk_nan(t, "frns2_49", pat(g12, g33, g12, pq));
  // brute force: every trio drawn from three alleles
  const double f[3] = {0.2, 0.3, 0.5};
  const int gs[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
  for (const auto &ch : gs) {
    for (const auto &m : gs) {
      for (const auto &a : gs) {
        double den = 0;
        for (const auto &g : gs) {
          const double pg = g[0] == g[1] ? f[g[0]] * f[g[0]]
                                         : 2 * f[g[0]] * f[g[1]];
          den += pg * trans(m, g, ch);
        }
        const double fc[2] = {f[ch[0]], f[ch[1]]};
        const double got = pat(ch, m, a, fc);
        if (den == 0) {
          tk_nan(t, "frns2_70", got);
        } else {
          tk_eq(t, "frns2_70", got, trans(m, a, ch) / den, 1e-10);
        }
      }
    }
  }
}

void guards(tk &t) {
  double o[4], w[4];
  const double h = 0.5, z = 0, six = 0.6, one = 1, m1 = -1;
  const int yes = 1, no = 0;
  tk_ieq(t, "crim36", clab_gtmp(&h, &z, &yes, 0, 0, o, o + 1), CLAB_EARG);
  tk_ieq(t, "crim38", clab_gtmp(&h, &z, &yes, 1, 1, o, o + 1), CLAB_EDOM);
  tk_ieq(t, "crim39", clab_gtmp(&z, &z, &yes, 1, 0, o, o + 1), CLAB_EDOM);
  tk_ieq(t, "crim40", clab_gtmp(&six, &six, &no, 1, 0, o, o + 1), CLAB_EDOM);
  tk_ieq(t, "crim41", clab_gtmp(&h, &z, &no, 1, 0, o, o + 1), CLAB_EDOM);
  tk_ieq(t, "crim42", clab_cllr(&one, 0, &one, 1, o), CLAB_EARG);
  tk_ieq(t, "crim43", clab_cllr(&z, 1, &one, 1, o), CLAB_EDOM);
  tk_ieq(t, "crim44", clab_cllr(&one, 1, &m1, 1, o), CLAB_EDOM);
  const double l21[2] = {2, 1}, l01[2] = {0, 1}, l12[2] = {1, 2};
  const int tf[2] = {1, 0}, tt[2] = {1, 1}, ff[2] = {0, 0};
  tk_ieq(t, "crim45", clab_lrcal(&one, &yes, 1, o, w, o + 1), CLAB_EARG);
  tk_ieq(t, "crim46", clab_lrcal(l21, tf, 2, o, w, o + 2), CLAB_EDOM);
  tk_ieq(t, "crim47", clab_lrcal(l01, tf, 2, o, w, o + 2), CLAB_EDOM);
  tk_ieq(t, "crim48", clab_lrcal(l12, tt, 2, o, w, o + 2), CLAB_EDOM);
  tk_ieq(t, "crim49", clab_lrcal(l12, ff, 2, o, w, o + 2), CLAB_EDOM);
  const double p01 = 0.1;
  double p21[21];
  int n21[21];
  for (int i = 0; i < 21; ++i) {
    p21[i] = 0.01;
    n21[i] = 1;
  }
  const double p66[2] = {0.6, 0.6};
  tk_ieq(t, "crim50", clab_mixp(&p01, &yes, 0, 1, o), CLAB_EARG);
  tk_ieq(t, "crim51", clab_mixp(&p01, &yes, 1, -1, o), CLAB_EARG);
  tk_ieq(t, "crim52", clab_mixp(p21, n21, 21, 11, o), CLAB_EARG);
  tk_ieq(t, "crim53", clab_mixp(&z, &yes, 1, 1, o), CLAB_EDOM);
  tk_ieq(t, "crim54", clab_mixp(p66, tf, 2, 1, o), CLAB_EDOM);
  const int g[2] = {0, 1}, g00[2] = {0, 0};
  const double f[2] = {0.1, 0.2}, f0[2] = {0, 0.2}, f15[2] = {0.1, 1.5};
  const double k[3] = {0, 1, 0}, kn[3] = {-1, 1, 1}, ks[3] = {0.5, 0.2, 0.2};
  tk_ieq(t, "crim58", clab_kinidx(g, g, f, 0, k, o), CLAB_EARG);
  tk_ieq(t, "crim60", clab_kinidx(g, g, f, 1, kn, o), CLAB_EDOM);
  tk_ieq(t, "crim61", clab_kinidx(g, g, f, 1, ks, o), CLAB_EDOM);
  tk_ieq(t, "crim62", clab_kinidx(g, g, f0, 1, k, o), CLAB_EDOM);
  tk_ieq(t, "crim63", clab_kinidx(g00, g, f, 1, k, o), CLAB_EDOM);
  tk_ieq(t, "crim64", clab_patidx(g, g, g, f, 0, o), CLAB_EARG);
  tk_ieq(t, "crim68", clab_patidx(g, g, g, f15, 1, o), CLAB_EDOM);
  tk_ieq(t, "crim69", clab_patidx(g00, g, g, f, 1, o), CLAB_EDOM);
  // the words for each code
  tk_true(t, "strerror", std::strcmp(clab_strerror(CLAB_OK), "no error") == 0);
  tk_true(t, "strerror",
          std::strcmp(clab_strerror(CLAB_EARG),
                      "a size or an index is out of range") == 0);
  tk_true(t, "strerror",
          std::strcmp(clab_strerror(CLAB_EDOM),
                      "a value is outside the domain of the method") == 0);
  tk_true(t, "strerror", std::strcmp(clab_strerror(99), "unknown error") == 0);
}

}  // namespace

void t_frns(tk &t) {
  gtmp(t);
  cost(t);
  cal(t);
  mix(t);
  kinship(t);
  paternity(t);
  guards(t);
}
