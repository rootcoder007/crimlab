// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_aorwt (src/aor.cpp). Part of the self-test (selftest.cpp).

#include <cmath>
#include <cstdint>

#include "crimlab.h"
#include "tchk.h"

namespace {

// aoristic 2.0.0, R 4.5.0, dev/ref/cpp/aor.R: 300 deterministic UTC
// events (i = 0..299; see section 3 below), aoristic.df(d, "x", "y",
// "from", "to"): ref_tot = colSums of the hour columns, ref_row0..2 = rows
// 1..3. hour1 is Sunday 00h (lubridate::wday 1).
const double ref_tot[168] = {1.7313913283399376, 2.9382950477747367, 1.9955372001258542, 2.713364801645711, 1.7287148471914218, 2.3025618182191181, 1.3304748829884907, 1.2674129101986078, 1.2741975983366987, 1.4909724538857558, 1.6193480076601678, 2.4420721733468951, 1.8018908339468667, 2.395195811976325, 2.401042680168751, 2.7004735405698974, 1.7809026631505021, 2.044278144800129, 1.6283224570296964, 1.3424636620778205, 1.2917136907318421, 1.0567885531379706, 1.054777398679446, 0.91912338888171152, 1.8310220691510057, 1.9164471525859155, 1.845188245222495, 1.169077319082451, 1.2632627687202105, 1.1071487244753495, 1.6726017638862587, 1.1122606490498441, 1.8464131292542534, 1.3925157390371283, 2.369779452819317, 1.7048773402052657, 3.7338338056659759, 1.8683744349993934, 2.5682632775790224, 1.4249924779942591, 2.2812217072869445, 1.1678272780823267, 2.4759830241079079, 1.1575473335702875, 1.8698446763478926, 1.8388525219281928, 1.7155215244408624, 1.6672912989224036, 2.0331004603437943, 2.2565491234395409, 2.2779269417926824, 2.0365969535428965, 1.9483870779450414, 1.8288139090543871, 1.6618733778439927, 1.6448700985909708, 1.4693720624630577, 1.8210759585093088, 1.9538281975565541, 1.9074045643367006, 0.94229684076549092, 2.0014169987243724, 0.7518382917395704, 1.6008933052321384, 1.0756086684722288, 1.2804062326392347, 1.5944901440330772, 1.3371987013678597, 1.7487040688820135, 1.941969470945637, 1.7969935589665345, 3.2359978166127603, 2.4100122821711745, 2.5439422636529452, 1.6265246575987242, 1.5338705450162342, 1.7591198909185932, 2.0093528085378631, 1.7237619176748871, 1.4952230375929467, 2.1632810621975378, 1.706274985803824, 1.5529855328972859, 1.9271823590449559, 2.3204052632206591, 1.8550464102838924, 2.9202635430316679, 1.4371879857647361, 2.6104320037522033, 1.163963627186223, 1.0710235091830882, 0.98127270004523015, 2.6101741094059752, 0.73073326665806337, 1.3540191954510983, 1.3390874336988914, 2.3993155493155491, 1.5654288944117059, 1.5087126528670975, 1.4816693533424043, 2.0462538926692155, 1.7253608520720749, 1.8533169651193326, 1.8866531954126684, 2.7065959171713549, 1.4949870075810712, 2.0178636425067622, 2.6500987742222422, 1.1680502794237644, 1.6544418175492679, 1.3858630499679132, 1.2944063206333083, 1.9507197439140218, 1.933167217901016, 1.6615192327408419, 1.7668398035695789, 1.928474789838716, 2.4009345263538324, 1.4200870925808622, 2.1633583393797045, 1.6736098865672622, 1.364264446330945, 1.6884744236238536, 2.1267602874518179, 1.0886272189314605, 2.5387859190119948, 1.1612155206323795, 0.8269037729694515, 1.2304474593494557, 1.0830404369424333, 1.3995532711554128, 2.4901051546608324, 2.3488942514421107, 2.1456109863283794, 2.9055871787237173, 2.0675923512530483, 2.0896410204924312, 1.647269254053964, 2.8866029128950386, 0.99386166913929752, 2.9909659162435447, 0.73360307167139771, 1.763075257237172, 1.4172645410162443, 2.4131242394668582, 1.7838534313788095, 2.300712404388976, 1.9201534954984525, 2.9799048894926812, 1.4354244996252545, 2.6572877332681588, 1.3909663509922161, 1.1155190381574944, 1.2010613383438105, 1.5791372201752969, 1.1256718672285471, 1.607434317851431, 1.3156403728377879, 1.4982675766932354, 1.7125617302332525, 2.1670308150794058, 0.93149341036477129, 1.6103007273722376, 1.531951591284529, 1.4375415748810132, 1.9051466114036482, 1.6127382688556358, 2.7969845525906827};
const double ref_row0[168] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
const double ref_row1[168] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01098901098901099, 0.65934065934065944, 0.32967032967032972, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
const double ref_row2[168] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.011111111111111112, 0.33333333333333337, 0.33333333333333337, 0.32222222222222224, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

const double SUN0 = 259200.0;  // 1970-01-04 00:00, a Sunday: bin 0
const double WEEK = 604800.0;
// Europe/London 2023: BST from 2023-03-26 01:00 UTC, GMT from 2023-10-29
// 01:00 UTC (R: as.POSIXct(..., tz = "UTC"))
const double LON_T[2] = {1679792400, 1698541200}, LON_O[2] = {3600, 0};

struct xrng {
  std::uint64_t s;
  double unif() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return static_cast<double>((s * 2685821657736338717ULL) >> 11) * 0x1p-53;
  }
};

void run(tk &tt) {
  double w[168], tot[168];
  int cnt[5];
  // 1. closed form (Ratcliffe's example shape): Monday 2023-01-02 10:30 to
  // 12:30 UTC puts .25, .5, .25 in Monday 10h, 11h, 12h (bins 34, 35, 36)
  {
    const double f = 1672617600 + 10.5 * 3600, t = f + 7200;
    tk_true(tt, "aor.hand.rc", clab_aorwt(&f, &t, 1, nullptr, nullptr, 0, 0, SUN0, 3600, 168,
                               120, WEEK, w, tot, cnt) == CLAB_OK);
    for (int j = 0; j < 168; ++j) {
      const double want = j == 34 || j == 36 ? .25 : (j == 35 ? .5 : 0);
      tk_near(tt, "aor.hand", w[j], want, 1e-15);
      tk_near(tt, "aor.hand.tot", tot[j], w[j], 0);
    }
  }
  // 2. across a clock change each real minute counts in the hour its wall
  // clock shows. Spring forward, 00:30 to 03:30 local (2 h elapsed): the
  // wall-clock hours are 00h (30 min), 02h (60) and 03h (30). Fall back,
  // 00:30 to 03:30 local (4 h elapsed): 00h 30 min, 01h 120 (twice), 02h
  // 60, 03h 30.
  {
    const double f = 1679790600, t = 1679797800;
    clab_aorwt(&f, &t, 1, LON_T, LON_O, 2, 0, SUN0, 3600, 168, 120, WEEK, w, nullptr, cnt);
    const double sp[4] = {.25, 0, .5, .25};
    for (int j = 0; j < 4; ++j) tk_near(tt, "aor.dst.spring", w[j], sp[j], 1e-15);
    const double f2 = 1698535800, t2 = 1698550200;
    clab_aorwt(&f2, &t2, 1, LON_T, LON_O, 2, 0, SUN0, 3600, 168, 120, WEEK, w, nullptr, cnt);
    const double fb[5] = {.125, .5, .25, .125, 0};
    for (int j = 0; j < 5; ++j) tk_near(tt, "aor.dst.fall", w[j], fb[j], 1e-15);
  }
  // 3. parity with aoristic 2.0.0 on 300 deterministic UTC events (whole
  // minutes, so its minute rounding is exact): totals and three rows
  {
    static double f[300], t[300], wr[300 * 168];
    for (int i = 0; i < 300; ++i) {
      f[i] = 1672617600 + 60.0 * ((i * 7919LL) % 200000);
      double sp = 2.0 + static_cast<double>((i * 104729LL) % 480);
      if (i % 50 == 7) sp = 8 * 1440;
      if (i % 50 == 21) sp = -30;
      if (i % 50 == 33) sp = 1;
      t[i] = i % 50 == 13 ? std::nan("") : f[i] + 60 * sp;
    }
    tk_true(tt, "aor.parity.rc", clab_aorwt(f, t, 300, nullptr, nullptr, 0, 0, SUN0, 3600, 168,
                                 120, WEEK, wr, tot, cnt) == CLAB_OK);
    for (int j = 0; j < 168; ++j) {
      tk_near(tt, "aor.parity.tot", tot[j], ref_tot[j], 1e-12);
      tk_near(tt, "aor.parity.row0", wr[0 + 300 * j], ref_row0[j], 1e-14);
      tk_near(tt, "aor.parity.row1", wr[1 + 300 * j], ref_row1[j], 1e-14);
      tk_near(tt, "aor.parity.row2", wr[2 + 300 * j], ref_row2[j], 1e-14);
    }
    tk_true(tt, "aor.parity.cnt", cnt[0] == 0 && cnt[1] == 6 && cnt[2] == 6 && cnt[3] == 6 && cnt[4] == 6);
    // every row sums to one
    for (int i = 0; i < 300; ++i) {
      double s = 0;
      for (int j = 0; j < 168; ++j) s += wr[i + 300 * j];
      tk_near(tt, "aor.rowsum", s, 1, 1e-13);
    }
  }
  // 4. brute force: random events on a whole-minute grid and random clock
  // changes (multiples of 30 min); each real minute is put in the bin of
  // its wall-clock time and divided by the length. 24 hourly bins and 7
  // daily bins.
  {
    xrng r{77};
    double ttr[6], to[6];
    double a = 1672617600;
    for (int k = 0; k < 6; ++k) {
      a += 60.0 * std::floor(r.unif() * 4000 + 30);
      ttr[k] = a;
      to[k] = 1800.0 * std::floor(r.unif() * 9 - 4);
    }
    const double off0 = 3600;
    for (int rep = 0; rep < 200; ++rep) {
      const double f = 1672617600 + 60.0 * std::floor(r.unif() * 26000);
      const double t = f + 60.0 * std::floor(r.unif() * 3000 + 1);
      for (int kind = 0; kind < 2; ++kind) {
        const double bw = kind == 0 ? 3600 : 86400;
        const int nb = kind == 0 ? 24 : 7;
        double got[24], want[24] = {0};
        clab_aorwt(&f, &t, 1, ttr, to, 6, off0, SUN0, bw, nb, 0, 1e300, got, nullptr, cnt);
        const double m = (t - f) / 60;
        for (double u = f; u < t; u += 60) {
          double o = off0;
          for (int k = 0; k < 6; ++k)
            if (ttr[k] <= u) o = to[k];
          double v = std::fmod(u + o - SUN0, bw * nb);
          if (v < 0) v += bw * nb;
          want[static_cast<int>(v / bw)] += 1 / m;
        }
        for (int j = 0; j < nb; ++j) tk_near(tt, "aor.brute", got[j], want[j], 1e-12);
      }
    }
  }
  // 5. the exact rule for long spans: 1.5 weeks with lo = 0, hi = Inf over
  // 7 daily bins starting Sunday 00:00 gives 1.5 days per day except the
  // first 3.5 days, which get 2: weights 2/10.5 and 1/10.5
  {
    const double f = SUN0 + 52 * WEEK, t = f + 1.5 * WEEK;
    double g[7];
    clab_aorwt(&f, &t, 1, nullptr, nullptr, 0, 0, SUN0, 86400, 7, 0, 1e300, g, nullptr, cnt);
    const double want[7] = {2 / 10.5, 2 / 10.5, 2 / 10.5, 1.5 / 10.5, 1 / 10.5, 1 / 10.5, 1 / 10.5};
    for (int j = 0; j < 7; ++j) tk_near(tt, "aor.long", g[j], want[j], 1e-15);
  }
  // 6. n = 0 and guards
  {
    tk_true(tt, "aor.n0", clab_aorwt(nullptr, nullptr, 0, nullptr, nullptr, 0, 0, SUN0, 3600, 168,
                          120, WEEK, nullptr, tot, cnt) == CLAB_OK && tot[5] == 0);
    const double f = 0, inf = HUGE_VAL, bad[2] = {5, 5}, o2[2] = {0, 0};
    tk_true(tt, "aor.guard.nb", clab_aorwt(&f, &f, 1, nullptr, nullptr, 0, 0, 0, 3600, 0, 0, 1, w, tot, cnt) == CLAB_EARG);
    tk_true(tt, "aor.guard.bw", clab_aorwt(&f, &f, 1, nullptr, nullptr, 0, 0, 0, 0, 24, 0, 1, w, tot, cnt) == CLAB_EDOM);
    tk_true(tt, "aor.guard.inf", clab_aorwt(&f, &inf, 1, nullptr, nullptr, 0, 0, 0, 3600, 24, 0, 1, w, tot, cnt) == CLAB_EDOM);
    tk_true(tt, "aor.guard.tt", clab_aorwt(&f, &f, 1, bad, o2, 2, 0, 0, 3600, 24, 0, 1, w, tot, cnt) == CLAB_EDOM);
    const double nanv = std::nan("");
    clab_aorwt(&nanv, &f, 1, nullptr, nullptr, 0, 0, 0, 3600, 24, 0, 1, w, tot, cnt);
    tk_true(tt, "aor.nostart", cnt[0] == 1 && tot[0] == 0);
  }
}


}  // namespace

void t_aor(tk &t, int) { run(t); }

#ifdef CLAB_MAIN
#include <chrono>
#include <cstdio>

// timing: n = 10,000 events, spans up to 8 h, row weights over 168 hours
// (aoristic::aoristic.df took 6.319 s, 1 run, arm64 Linux)
void b_aor() {
  const int n = 10000;
  static double f[n], t[n], w[n * 168];
  double tot[168], tm[5];
  int cnt[5];
  xrng r{3};
  for (int i = 0; i < n; ++i) {
    f[i] = 1672617600 + r.unif() * 3e7;
    t[i] = f[i] + r.unif() * 8 * 3600;
  }
  for (int k = 0; k < 5; ++k) {
    const auto a = std::chrono::steady_clock::now();
    clab_aorwt(f, t, n, LON_T, LON_O, 2, 0, SUN0, 3600, 168, 120, WEEK, w, tot, cnt);
    tm[k] = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
  }
  std::printf("BENCH clab_aorwt n=1e4 (row weights + totals): median %.6f s (k=5)\n", tk_med5(tm));
}
#endif
