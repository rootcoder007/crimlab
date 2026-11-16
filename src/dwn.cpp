// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Downwind exposure (crimlab.h, dwnexp()): whether places lie downwind of
// moving or fixed sources of pollution on each period, by three rules
// from the papers: a sector about the direction the wind blows to from
// each source, to a place's point (the centroid rule of Bao and Da 2026);
// the same sector, to a radius, against a place's polygon (their
// transport wedge); and the angle between the bearing from the sources'
// mean point to a place and the sources' mean wind direction, with the
// mean distance to the sources (Ayesh 2023); and the tailwind component
// of the wind along a flight course, by row or by group (Hener 2022).
// Bearings are degrees
// clockwise from north, planar or on the sphere. Nothing here allocates
// more than O(sources + places + periods).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

using std::size_t;
constexpr double pi = 3.14159265358979323846;
constexpr double rd = pi / 180.0;
constexpr double erad = 6371.0088;  // the mean radius of the earth, km

inline size_t z(int i) { return static_cast<size_t>(i); }

// The bearing from (x0, y0) to (x1, y1), degrees clockwise from north in
// [0, 360): planar, or the initial great-circle bearing with x the
// longitude and y the latitude in degrees.
double brg(double x0, double y0, double x1, double y1, int ll) {
  double b = 0.0;
  if (ll) {
    const double p0 = y0 * rd, p1 = y1 * rd, dl = (x1 - x0) * rd;
    b = std::atan2(std::sin(dl) * std::cos(p1),
                   std::cos(p0) * std::sin(p1) - std::sin(p0) * std::cos(p1) * std::cos(dl));
  } else {
    b = std::atan2(x1 - x0, y1 - y0);
  }
  b /= rd;
  return b < 0.0 ? b + 360.0 : b;
}

// The distance: planar, or the great-circle distance in km (haversine).
double dst(double x0, double y0, double x1, double y1, int ll) {
  if (!ll) return std::hypot(x1 - x0, y1 - y0);
  const double p0 = y0 * rd, p1 = y1 * rd;
  const double a = std::sin((p1 - p0) / 2), c = std::sin((x1 - x0) * rd / 2);
  const double h = a * a + std::cos(p0) * std::cos(p1) * c * c;
  return 2.0 * erad * std::asin(std::sqrt(std::min(1.0, h)));
}

// The angle between two directions, in [0, 180].
double adiff(double a, double b) {
  const double d = std::fmod(std::fabs(a - b), 360.0);
  return d > 180.0 ? 360.0 - d : d;
}

// A point within the sector of half-angle hf about th from (ax, ay), to
// the radius rad; the apex itself is within.
bool insec(double px, double py, double ax, double ay, double th, double hf, double rad, int ll) {
  const double d = dst(ax, ay, px, py, ll);
  if (d == 0.0) return true;
  if (d > rad) return false;
  return adiff(brg(ax, ay, px, py, ll), th) <= hf;
}

// The orientation of (c - a) against (b - a): +1, -1 or 0.
int ori(double ax, double ay, double bx, double by, double cx, double cy) {
  const double v = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
  return v > 0.0 ? 1 : (v < 0.0 ? -1 : 0);
}

bool onseg(double ax, double ay, double bx, double by, double px, double py) {
  return std::min(ax, bx) <= px && px <= std::max(ax, bx) && std::min(ay, by) <= py &&
         py <= std::max(ay, by);
}

// Whether the segments p1p2 and q1q2 meet (touching counts).
bool segx(const double *p, const double *q) {
  const int o1 = ori(p[0], p[1], p[2], p[3], q[0], q[1]);
  const int o2 = ori(p[0], p[1], p[2], p[3], q[2], q[3]);
  const int o3 = ori(q[0], q[1], q[2], q[3], p[0], p[1]);
  const int o4 = ori(q[0], q[1], q[2], q[3], p[2], p[3]);
  if (o1 != o2 && o3 != o4) return true;
  return (o1 == 0 && onseg(p[0], p[1], p[2], p[3], q[0], q[1])) ||
         (o2 == 0 && onseg(p[0], p[1], p[2], p[3], q[2], q[3])) ||
         (o3 == 0 && onseg(q[0], q[1], q[2], q[3], p[0], p[1])) ||
         (o4 == 0 && onseg(q[0], q[1], q[2], q[3], p[2], p[3]));
}

// Whether (x, y) is inside the ring of nv vertices (even-odd rule).
bool pip(const double *rx, const double *ry, int nv, double x, double y) {
  bool in = false;
  for (int i = 0, j = nv - 1; i < nv; j = i++) {
    if ((ry[i] > y) != (ry[j] > y) &&
        x < (rx[j] - rx[i]) * (y - ry[i]) / (ry[j] - ry[i]) + rx[i]) {
      in = !in;
    }
  }
  return in;
}

// Whether the segment e meets the arc of radius rad about (ax, ay)
// within hf of th (planar).
bool segarc(const double *e, double ax, double ay, double th, double hf, double rad) {
  const double dx = e[2] - e[0], dy = e[3] - e[1], fx = e[0] - ax, fy = e[1] - ay;
  const double a = dx * dx + dy * dy, b = 2 * (fx * dx + fy * dy), c = fx * fx + fy * fy - rad * rad;
  const double disc = b * b - 4 * a * c;
  if (a == 0.0 || disc < 0.0) return false;
  const double s = std::sqrt(disc);
  for (int k = 0; k < 2; ++k) {
    const double t = (-b + (k ? s : -s)) / (2 * a);
    if (t < 0.0 || t > 1.0) continue;
    if (adiff(brg(ax, ay, e[0] + t * dx, e[1] + t * dy, 0), th) <= hf) return true;
  }
  return false;
}

// Whether the ring meets the wedge of half-angle hf about th from (ax,
// ay), to the radius rad (planar): a vertex in the wedge, the apex in the
// ring, or an edge across a side or the arc.
bool wedge(const double *rx, const double *ry, int nv, double ax, double ay, double th, double hf,
           double rad) {
  double far = 0.0;
  for (int i = 0; i < nv; ++i) {
    if (insec(rx[i], ry[i], ax, ay, th, hf, rad, 0)) return true;
    far = std::max(far, std::hypot(rx[i] - ax, ry[i] - ay));
  }
  if (pip(rx, ry, nv, ax, ay)) return true;
  const double len = std::isfinite(rad) ? rad : 2.0 * far + 1.0;
  double sd[2][4];
  for (int k = 0; k < 2; ++k) {
    const double a = (th + (k ? hf : -hf)) * rd;
    sd[k][0] = ax;
    sd[k][1] = ay;
    sd[k][2] = ax + len * std::sin(a);
    sd[k][3] = ay + len * std::cos(a);
  }
  for (int i = 0, j = nv - 1; i < nv; j = i++) {
    const double e[4] = {rx[j], ry[j], rx[i], ry[i]};
    if (segx(e, sd[0]) || segx(e, sd[1])) return true;
    if (std::isfinite(rad) && segarc(e, ax, ay, th, hf, rad)) return true;
  }
  return false;
}

// Rules 0 and 1: per source, the places it reaches, counted by period.
int reach(int rule, int ns, const double *sx, const double *sy, const int *sp, const double *sto,
          int nr, const double *rx, const double *ry, const int *ro, double hf, double rad, int ll,
          double *out) {
  for (int s = 0; s < ns; ++s) {
    if (kutil::stop(static_cast<double>(nr) * (rule == 1 ? 8.0 : 1.0))) return CLAB_EINTR;
    double *o = out + z(sp[s]) * z(nr);
    for (int r = 0; r < nr; ++r) {
      bool hit = false;
      if (rule == 0) {
        hit = insec(rx[r], ry[r], sx[s], sy[s], sto[s], hf, rad, ll);
      } else {
        hit = wedge(rx + ro[r], ry + ro[r], ro[r + 1] - ro[r], sx[s], sy[s], sto[s], hf, rad);
      }
      if (hit) o[r] += 1.0;
    }
  }
  return CLAB_OK;
}

// Rule 2: per period, the angle between the bearing from the sources'
// mean point to each place and the sources' mean direction (the mean
// wind vector), plain or about the circle, and the mean distance from
// the place to the sources not at its own location.
int angle(int ns, const double *sx, const double *sy, const int *sp, const double *sto, int nr,
          const double *rx, const double *ry, int np, int ll, int wrap, double *out, double *out2) {
  std::vector<double> mx(z(np), 0.0), my(z(np), 0.0), sn(z(np), 0.0), cs(z(np), 0.0);
  std::vector<int> cnt(z(np), 0);
  for (int s = 0; s < ns; ++s) {
    const size_t p = z(sp[s]);
    mx[p] += sx[s];
    my[p] += sy[s];
    sn[p] += std::sin(sto[s] * rd);
    cs[p] += std::cos(sto[s] * rd);
    ++cnt[p];
  }
  for (int p = 0; p < np; ++p) {
    if (kutil::stop(static_cast<double>(nr) * (cnt[z(p)] + 1))) return CLAB_EINTR;
    double *o = out + z(p) * z(nr), *o2 = out2 + z(p) * z(nr);
    if (cnt[z(p)] == 0) {
      std::fill(o, o + nr, NAN);
      std::fill(o2, o2 + nr, NAN);
      continue;
    }
    const double ax = mx[z(p)] / cnt[z(p)], ay = my[z(p)] / cnt[z(p)];
    double dir = std::atan2(sn[z(p)], cs[z(p)]) / rd;
    if (dir < 0.0) dir += 360.0;
    for (int r = 0; r < nr; ++r) {
      const double b = brg(ax, ay, rx[r], ry[r], ll);
      o[r] = wrap ? adiff(b, dir) : std::fabs(b - dir);
      double sd = 0.0;
      int k = 0;
      for (int s = 0; s < ns; ++s) {
        if (sp[s] != p) continue;
        const double d = dst(sx[s], sy[s], rx[r], ry[r], ll);
        if (d > 0.0) {
          sd += d;
          ++k;
        }
      }
      o2[r] = k > 0 ? sd / k : NAN;
    }
  }
  return CLAB_OK;
}

bool allfin(const double *v, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    if (!std::isfinite(v[i])) return false;
  }
  return true;
}

}  // namespace

extern "C" int clab_dwnexp(int rule, int ns, const double *sx, const double *sy, const int *sp,
                           const double *sto, int nr, const double *rx, const double *ry,
                           const int *ro, int np, double half, double rad, int ll, int wrap,
                           double *out, double *out2) CLAB_TRY {
  if (rule < 0 || rule > 2 || ns < 1 || nr < 1 || np < 1 || !(half >= 0.0 && half <= 180.0) ||
      !(rad > 0.0) || (rule == 1 && (ll || !ro)) || (rule == 2 && !out2)) {
    return CLAB_EARG;
  }
  if (static_cast<double>(nr) * np > kutil::maxcell) return CLAB_EARG;
  for (int s = 0; s < ns; ++s) {
    if (sp[s] < 0 || sp[s] >= np) return CLAB_EARG;
  }
  const int nv = rule == 1 ? ro[nr] : nr;
  if (rule == 1) {
    if (ro[0] != 0) return CLAB_EARG;
    for (int r = 0; r < nr; ++r) {
      if (ro[r + 1] - ro[r] < 3) return CLAB_EARG;
    }
  }
  if (!allfin(sx, z(ns)) || !allfin(sy, z(ns)) || !allfin(sto, z(ns)) || !allfin(rx, z(nv)) ||
      !allfin(ry, z(nv))) {
    return CLAB_EDOM;
  }
  const size_t cells = z(nr) * z(np);
  if (rule == 2) return angle(ns, sx, sy, sp, sto, nr, rx, ry, np, ll, wrap, out, out2);
  std::fill(out, out + cells, 0.0);
  return reach(rule, ns, sx, sy, sp, sto, nr, rx, ry, ro, half, rad, ll, out);
}
CLAB_CATCH(CLAB_ENOMEM)

extern "C" int clab_tailwd(int n, const double *dir, const double *spd, double course,
                           const int *grp, int ng, int agg, double *out) CLAB_TRY {
  if (n < 1 || !std::isfinite(course) || (grp && (ng < 1 || agg < 0 || agg > 1))) return CLAB_EARG;
  if (!grp) {
    for (int i = 0; i < n; ++i) {
      out[i] = spd[i] * std::cos((dir[i] + 180.0 - course) * rd);
    }
    return CLAB_OK;
  }
  for (int i = 0; i < n; ++i) {
    if (grp[i] < 0 || grp[i] >= ng) return CLAB_EARG;
  }
  std::vector<double> sm(z(ng), 0.0);
  std::vector<int> cnt(z(ng), 0);
  std::fill(out, out + ng, NAN);
  for (int i = 0; i < n; ++i) {
    const double t = spd[i] * std::cos((dir[i] + 180.0 - course) * rd);
    if (std::isnan(t)) continue;
    const size_t g = z(grp[i]);
    if (agg == 0) {
      out[g] = cnt[g] == 0 ? t : std::max(out[g], t);
    } else {
      sm[g] += t;
    }
    ++cnt[g];
  }
  if (agg == 1) {
    for (int g = 0; g < ng; ++g) out[g] = cnt[z(g)] ? sm[z(g)] / cnt[z(g)] : NAN;
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
