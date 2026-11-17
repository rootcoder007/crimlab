# References for the distributions of src/dst.cpp at large parameters
# (the values in src/t_dst.cpp, t_dstbig()), with mpmath at 34 to 40
# digits. Lines on standard input, one value out per line:
#   pnchisq x df ncp    P(X <= x) and P(X > x), noncentral chi-square, by
#                       the integral over the normal part,
#                       int phi(z) P(chi2_{df-1} <= x - (z + sqrt(ncp))^2) dz
#                       (src/dst.cpp integrates over the chi part instead)
#   ncpq st df q        the ncp with P(X <= st) = q (secant method on the
#                       above; df 1 in closed form)
#   pbeta x a b         I_x(a, b) and 1 - I_x(a, b) by the defining integral
#   pgamma x a          P(a, x) and Q(a, x) by the defining integral
#   pnf x m n ncp       P(F > x), noncentral F, by its Poisson mixture of
#                       beta tails
# Numbers are read as the exact values of their doubles. Run with
# python3 -I dev/ref/dstbig.py < in > out (mpmath 1.4).
import sys
from decimal import Decimal

import mpmath as mp


def d(v):
    return mp.mpf(Decimal(float(v)))


def pnchisq(x, df, lam):
    mp.mp.dps = 34
    r, m, sx = mp.sqrt(lam), df - 1, mp.sqrt(x)
    lo, hi = max(-r - sx, mp.mpf(-45)), min(-r + sx, mp.mpf(45))
    if lo >= hi:
        return mp.mpf(0)
    if m == 0:
        return mp.ncdf(sx - r) - mp.ncdf(-sx - r)
    f = lambda z: mp.npdf(z) * mp.gammainc(m / 2, 0, max((x - (z + r) ** 2) / 2, 0), regularized=True)
    # a fine grid of break points, and more where (z + r)^2 = x - m
    pts = [lo + (hi - lo) * k / 400 for k in range(401)]
    if x > m:
        c = mp.sqrt(x - m) - r
        w = mp.sqrt(2 * m) / (2 * mp.sqrt(x - m)) + mp.mpf(1e-30)
        pts += [c + j * w for j in range(-60, 61) if lo < c + j * w < hi]
    pts = sorted(set(pts))
    return mp.fsum(mp.quad(f, [pts[i], pts[i + 1]]) for i in range(len(pts) - 1))


def ncpq(st, df, q):
    mp.mp.dps = 40
    if df == 1:
        g = lambda lam: mp.ncdf(mp.sqrt(st) - mp.sqrt(lam)) - mp.ncdf(-mp.sqrt(st) - mp.sqrt(lam)) - q
    else:
        g = lambda lam: pnchisq(st, df, lam) - q
    z = mp.sqrt(2) * mp.erfinv(2 * q - 1)
    l0 = max(st - df - z * mp.sqrt(2 * df + 4 * st), mp.mpf(1))
    return mp.findroot(g, (l0, l0 * (1 + mp.mpf(10) ** -6)), solver="secant", tol=mp.mpf(10) ** -60)


def tail(logf, x, sd, sign, lim):
    # the integral of exp(logf) from x toward lim, in pieces of sd / 4
    pts = [x + sign * sd * k / 4 for k in range(0, 4 * 60 + 1)]
    pts = [p for p in pts if (p > lim if sign < 0 else p < lim)] + [lim]
    pts = sorted(set(pts))
    return mp.fsum(mp.quad(lambda t: mp.exp(logf(t)), [pts[i], pts[i + 1]]) for i in range(len(pts) - 1))


def pbeta(x, a, b):
    mp.mp.dps = 40
    lb = mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b)
    logf = lambda t: (a - 1) * mp.log(t) + (b - 1) * mp.log1p(-t) - lb
    m, sd = a / (a + b), mp.sqrt(a * b / ((a + b) ** 2 * (a + b + 1)))
    if x <= m:
        return tail(logf, x, sd, -1, max(mp.mpf(0), x - 60 * sd))
    return 1 - tail(logf, x, sd, 1, min(mp.mpf(1), x + 60 * sd))


def pgamma(x, a):
    mp.mp.dps = 40
    lg = mp.loggamma(a)
    logf = lambda t: (a - 1) * mp.log(t) - t - lg
    sd = mp.sqrt(a)
    if x <= a - 1:
        return tail(logf, x, sd, -1, max(mp.mpf(0), x - 60 * sd))
    return 1 - tail(logf, x, sd, 1, x + 60 * sd)


def pnf(x, m, n, ncp):
    mp.mp.dps = 30
    y, mu = m * x / (m * x + n), ncp / 2
    j0, s = int(mu), mp.mpf(0)
    for step in (1, -1):
        j = j0 if step == 1 else j0 - 1
        while j >= 0:
            wj = mp.exp(-mu + j * mp.log(mu) - mp.loggamma(j + 1))
            if wj < mp.mpf(10) ** -45 and abs(j - j0) > 5:
                break
            s += wj * mp.betainc(m / 2 + j, n / 2, y, 1, regularized=True)
            j += step
    return s


for line in sys.stdin:
    w = line.split()
    v = [d(u) for u in w[1:]]
    if w[0] == "pnchisq":
        p = pnchisq(v[0], int(v[1]), v[2])
        out = [p, 1 - p]
    elif w[0] == "ncpq":
        out = [ncpq(v[0], int(v[1]), v[2])]
    elif w[0] == "pbeta":
        p = pbeta(*v)
        out = [p, 1 - p]
    elif w[0] == "pgamma":
        p = pgamma(*v)
        out = [p, 1 - p]
    else:
        out = [pnf(*v)]
    print(" ".join(w), " ".join(mp.nstr(o, 22) for o in out), flush=True)
