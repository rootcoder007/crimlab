# The conditional maximum likelihood odds ratio of a 2 x 2 table and its
# exact interval (fisher.test()'s definitions) at 40 digits, for the
# references of src/t_tab.cpp: the noncentral hypergeometric terms from
# log-gamma values of 40 digits over the whole support, the roots by
# bisection on the log odds ratio. Input lines "a b c d conf" (the table
# by column: a, b the first column); output the estimate and the ends.
# Run with python3 -I dev/ref/fshor.py < in > out (mpmath 1.4).
import sys

import mpmath as mp

mp.mp.dps = 40


def run(a, b, c, d, conf):
    m, nn, k = a + b, c + d, a + c
    lo, hi = max(0, k - nn), min(k, m)
    xs = range(lo, hi + 1)
    ld = [mp.loggamma(m + 1) - mp.loggamma(x + 1) - mp.loggamma(m - x + 1) + mp.loggamma(nn + 1)
          - mp.loggamma(k - x + 1) - mp.loggamma(nn - k + x + 1) for x in xs]

    def w(lt):
        v = [ld[i] + x * lt for i, x in enumerate(xs)]
        top = max(v)
        return [mp.exp(u - top) for u in v]

    def mean(lt):
        e = w(lt)
        return mp.fsum(e[i] * x for i, x in enumerate(xs)) / mp.fsum(e)

    def tail(lt, up):
        e = w(lt)
        return mp.fsum(e[i] for i, x in enumerate(xs) if (x >= a if up else x <= a)) / mp.fsum(e)

    def bis(g, lo, hi):
        # the root of the increasing g on [lo, hi], 130 halvings
        for _ in range(130):
            mid = (lo + hi) / 2
            if g(mid) < 0:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2

    al = (1 - mp.mpf(conf)) / 2
    est = bis(lambda lt: mean(lt) - a, mp.mpf(-20), mp.mpf(20))
    lo_ = bis(lambda lt: tail(lt, True) - al, est - 10, est)
    hi_ = bis(lambda lt: al - tail(lt, False), est, est + 10)
    return [mp.exp(est), mp.exp(lo_), mp.exp(hi_)]


for line in sys.stdin:
    v = line.split()
    out = run(int(v[0]), int(v[1]), int(v[2]), int(v[3]), mp.mpf(v[4]))
    print(" ".join(v), " ".join(mp.nstr(o, 22) for o in out), flush=True)
