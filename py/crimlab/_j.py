# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""The functions that reach the kernels through clab_jrun() (cli/jrun.cpp).

Each takes the arguments of the R function of its name, with the same
defaults (None stands for an argument left out, R's default), and checks
them in the order and words of R. An argument whose R name is a Python
keyword (``from``, ``lambda``) takes a final underscore. The result is that of R: a data frame a Table, a named vector a
Vector, a list a Record, an unnamed vector a list (one value: a number), a
matrix a list of its rows; numbers are floats."""

from __future__ import annotations

from typing import Any

from ._util import jcall


def fair(score: Any, y: Any, grp: Any, cut: Any = 0.5, base: Any = 1, conf: Any = 0.95) -> Any:
    """Fairness metrics of a binary risk tool by group.

    The arguments and the result are those of R's ``fair()``: see its help
    page.

    Examples
    --------
    >>> r = fair(score=[0.2, 0.7, 0.9, 0.4, 0.6, 0.1, 0.8, 0.3], y=[0, 1, 1, 0, 1, 0, 0,
    ... 1], grp=["a", "a", "a", "a", "b", "b", "b", "b"])
    >>> list(r)
    ['cnt', 'met', 'rat', 'lo', 'hi']
    """
    return jcall("fair", score=score, y=y, grp=grp, cut=cut, base=base, conf=conf)


def fairdf(pred: Any, grp: Any = None, y: Any = None, base: Any = None, conf: Any = 0.95,
           nboot: Any = 0, seed: Any = None) -> Any:
    """Differences between groups in a risk tool's decisions.

    The arguments and the result are those of R's ``fairdf()``: see its help
    page.

    Examples
    --------
    >>> r = fairdf(pred=[1, 0, 1, 1, 0, 0, 1, 0], grp=["a", "a", "a", "a", "b", "b", "b",
    ... "b"], y=[1, 0, 1, 0, 0, 1, 1, 0])
    >>> list(r)
    ['tab', 'sum', 'base']
    """
    return jcall("fairdf", pred=pred, grp=grp, y=y, base=base, conf=conf, nboot=nboot, seed=seed)


def tmpaud(pred: Any, grp: Any = None, period: Any = None, city: Any = None,
           base: Any = None) -> Any:
    """Disparity metrics over time and place.

    The arguments and the result are those of R's ``tmpaud()``: see its help
    page.

    Examples
    --------
    >>> r = tmpaud(pred=[1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1], grp=["a", "b", "a", "b", "a",
    ... "b", "a", "b", "a", "b", "a", "b"], period=[1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2])
    >>> list(r)
    ['cells', 'cities', 'sum', 'base']
    """
    return jcall("tmpaud", pred=pred, grp=grp, period=period, city=city, base=base)


def lipaud(h: Any, d: Any = None, pairs: Any = None, lip: Any = 1) -> Any:
    """Individual fairness: a Lipschitz audit.

    The arguments and the result are those of R's ``lipaud()``: see its help
    page.

    Examples
    --------
    >>> r = lipaud(h=[0.1, 0.5, 0.4], pairs=[[1, 2], [2, 3]], d=[1, 2])
    >>> list(r)
    ['lreq', 'nviol', 'rate', 'ok', 'maxgap', 'worst', 'pairs']
    """
    return jcall("lipaud", h=h, d=d, pairs=pairs, lip=lip)


def nsbnd(pobs: Any, a: Any, b: Any) -> Any:
    """A true base rate behind a noisily recorded one.

    The arguments and the result are those of R's ``nsbnd()``: see its help
    page.

    Examples
    --------
    >>> r = nsbnd(pobs=[0.2, 0.25], a=[0, 0.1], b=0.2)
    >>> list(r)
    ['bnd', 'pairs']
    """
    return jcall("nsbnd", pobs=pobs, a=a, b=b)


def rnkres(est: Any, hw: Any) -> Any:
    """Which pairs a noisy ranking orders.

    The arguments and the result are those of R's ``rnkres()``: see its help
    page.

    Examples
    --------
    >>> r = rnkres(est=[1, 2, 3.5], hw=0.8)
    >>> list(r)
    ['id', 'npair', 'nunid', 'share', 'res']
    """
    return jcall("rnkres", est=est, hw=hw)


def nordet(x: Any, off: Any, radius: Any, p: Any = 0.85, draw: Any = False,
           seed: Any = None) -> Any:
    """Noisy-OR detection by patrols.

    The arguments and the result are those of R's ``nordet()``: see its help
    page.

    Examples
    --------
    >>> r = nordet(x=[[0, 0], [1, 1], [3, 3]], off={"x": [0, 2.5], "y": [0.5, 3]}, radius=1)
    >>> list(r)
    ['k', 'prob']
    """
    return jcall("nordet", x=x, off=off, radius=radius, p=p, draw=draw, seed=seed)


def fairsim(n: Any = 2000, groups: Any = None, props: Any = None, nareas: Any = 20,
            base: Any = 0.3, bias: Any = 0.5, seed: Any = None) -> Any:
    """Simulated records with a known disparity.

    The arguments and the result are those of R's ``fairsim()``: see its help
    page.

    Examples
    --------
    >>> r = fairsim(n=20, nareas=4, seed=1)
    >>> list(r)
    ['area', 'group', 'truth', 'det', 'risk']
    """
    return jcall("fairsim", n=n, groups=groups, props=props, nareas=nareas, base=base, bias=bias,
                 seed=seed)


def oaxbd(y: Any, x: Any, grp: Any, ref: Any = 1, conf: Any = 0.95) -> Any:
    """Oaxaca-Blinder decomposition of a gap between two groups.

    The arguments and the result are those of R's ``oaxbd()``: see its help
    page.

    Examples
    --------
    >>> r = oaxbd(y=[8.31, 8.69, 7.96, 6.43, 10.83, 8.04, 7.15, 9.85, 5.73, 10.2, 6.38,
    ... 7.06, 8.45, 8.97, 9.35, 8, 7.96, 7.36, 7.49, 8.5, 7.47, 4.6, 5.73, 5.59, 6.11, 4.6,
    ... 6.35, 5.62, 5.58, 3.4, 5.07, 4.89, 5.13, 6.85, 5.53, 6.2, 6.78, 6.27, 5.35, 6.12],
    ... x=[4.04, 4.71, 5.26, 3.85, 5.2, 5.03, 5.09, 6.12, 3.78, 6.27, 4.26, 3.87, 4.28,
    ... 5.25, 5.15, 4.69, 4.05, 4.35, 6.22, 5.2, 3.42, 3.06, 3.8, 2.33, 3.52, 3.26, 5.16,
    ... 5.01, 3.93, 2.86, 4.9, 4.85, 4.73, 4.74, 3.65, 4.71, 5.3, 4.04, 3.02, 4.79],
    ... grp=["a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a",
    ... "a", "a", "a", "a", "a", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b",
    ... "b", "b", "b", "b", "b", "b", "b", "b", "b"], conf=0.9)
    >>> list(r)
    ['groups', 'means', 'coef', 'dec', 'by']
    """
    return jcall("oaxbd", y=y, x=x, grp=grp, ref=ref, conf=conf)


def tlvc(x: Any, item: Any) -> Any:
    """Two-level model for multivariate trace evidence.

    The arguments and the result are those of R's ``tlvc()``: see its help
    page.

    Examples
    --------
    >>> r = tlvc(x=[[8.931, 13.104], [8.903, 12.902], [8.762, 13.092], [9.969, 10.554],
    ... [10.491, 11.067], [10.356, 11.176], [8.298, 8.685], [8.035, 8.549], [8.233, 8.872],
    ... [13.275, 5.544], [13.463, 6.051], [13.17, 5.563], [10.737, 12.388], [10.648,
    ... 12.256], [10.384, 12.101], [8.276, 9.948], [8.28, 9.549], [8.347, 10.203], [11.195,
    ... 9.999], [11.128, 10.403], [10.942, 10.063], [11.426, 11.746], [11.616, 12.01],
    ... [11.588, 11.701], [11.014, 11.391], [11.011, 11.7], [11.225, 11.553], [9.543,
    ... 11.188], [9.367, 11.203], [9.565, 11.07]], item=[1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4,
    ... 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10])
    >>> list(r)
    ['mu', 'u', 'c']
    """
    return jcall("tlvc", x=x, item=item)


def tlvlr(xa: Any, xb: Any, na: Any, nb: Any, bg: Any) -> Any:
    """Two-level model for multivariate trace evidence.

    The arguments and the result are those of R's ``tlvlr()``: see its help
    page.

    Examples
    --------
    >>> r = tlvlr(xa=[8.747, 13.024], xb=[8.797, 13.074], na=3, nb=3, bg={"mu": [10.2726,
    ... 10.521033333333333], "u": [[0.02122506666666663, 0.01757643333333333],
    ... [0.01757643333333333, 0.04655626666666665]], "c": [[2.552874898765432,
    ... -1.2094229691358025], [-1.2094229691358025, 4.371759675308642]]})
    >>> list(r)
    ['ln', 'log10']
    """
    return jcall("tlvlr", xa=xa, xb=xb, na=na, nb=nb, bg=bg)


def aorwt(from_: Any, to: Any, bw: Any = 3600, nb: Any = 168, origin: Any = 259200, off: Any = 0,
          tt: Any = None, toff: Any = None, lo: Any = 120, hi: Any = 604800,
          each: Any = False) -> Any:
    """Aoristic weights.

    The arguments and the result are those of R's ``aorwt()``: see its help
    page.

    Examples
    --------
    >>> r = aorwt(from_=[79200, 120600], to=[93600, 124200], bw=3600, nb=24, origin=0)
    >>> list(r)
    ['tot', 'w', 'cnt']
    """
    return jcall("aorwt", from_=from_, to=to, bw=bw, nb=nb, origin=origin, off=off, tt=tt,
                 toff=toff, lo=lo, hi=hi, each=each)


def bdig(x: Any, nd: Any = 1, dp: Any = -1, sgn: Any = "pos") -> Any:
    """Benford digit screens.

    The arguments and the result are those of R's ``bdig()``: see its help
    page.

    Examples
    --------
    >>> r = bdig(x=[123, 45.6, 0.29, -7, 1234], nd=2, dp=2, sgn="both")
    >>> list(r)
    ['dig', 'cnt']
    """
    return jcall("bdig", x=x, nd=nd, dp=dp, sgn=sgn)


def bnft(dig: Any, nd: Any = 1, nsim: Any = 0, seed: Any = None) -> Any:
    """Benford digit screens.

    The arguments and the result are those of R's ``bnft()``: see its help
    page.

    Examples
    --------
    >>> r = bnft(dig=[1, 2, 1, 3, 1, 4, 1, 5, 9, 2, 1, 1, 6])
    >>> list(r)
    ['obs', 'stat']
    """
    return jcall("bnft", dig=dig, nd=nd, nsim=nsim, seed=seed)


def skp(x: Any, s2: Any, alpha: Any = 0.05, nsim: Any = 999, seed: Any = None) -> Any:
    """Scott-Knott grouping of measurements.

    The arguments and the result are those of R's ``skp()``: see its help
    page.

    Examples
    --------
    >>> r = skp(x=[1.517975, 1.518007, 1.517967, 1.518064, 1.518013, 1.519467, 1.519519,
    ... 1.51953, 1.519523, 1.519488], s2=1.6000000000000003e-09, nsim=99, seed=1)
    >>> len(r)
    10
    """
    return jcall("skp", x=x, s2=s2, alpha=alpha, nsim=nsim, seed=seed)


def rptc(unit: Any, year: Any, rep: Any, y0: Any, y1: Any, minrep: Any = 12) -> Any:
    """Constant reporters of a reporting panel.

    The arguments and the result are those of R's ``rptc()``: see its help
    page.

    Examples
    --------
    >>> r = rptc(unit=["a", "a", "a", "b", "b", "b"], year=[2001, 2002, 2003, 2001, 2002,
    ... 2003], rep=[12, 12, 12, 5, 6, 7], y0=2001, y1=2003, minrep=1)
    >>> list(r)
    ['unit', 'const']
    """
    return jcall("rptc", unit=unit, year=year, rep=rep, y0=y0, y1=y1, minrep=minrep)


def rswp(d: Any, gam: Any = 1, alpha: Any = 0.05) -> Any:
    """Rosenbaum sensitivity bounds for matched pairs.

    The arguments and the result are those of R's ``rswp()``: see its help
    page.

    Examples
    --------
    >>> r = rswp(d=[1.2, 0.5, -0.3, 2, 1.1, 0.7, 0.9, -0.1], gam=[1, 1.5, 2])
    >>> list(r)
    ['p', 'gcrit', 't']
    """
    return jcall("rswp", d=d, gam=gam, alpha=alpha)


def rshl(d: Any, gam: Any = 1) -> Any:
    """Rosenbaum sensitivity bounds for matched pairs.

    The arguments and the result are those of R's ``rshl()``: see its help
    page.

    Examples
    --------
    >>> r = rshl(d=[1.2, 0.5, -0.3, 2, 1.1, 0.7, 0.9, -0.1])
    >>> list(r)
    ['gam', 'lo', 'hi']
    """
    return jcall("rshl", d=d, gam=gam)


def lwss(x: Any, y: Any, f: Any = None, iter: Any = 3, delta: Any = None) -> Any:
    """Lowess smoothing.

    The arguments and the result are those of R's ``lwss()``: see its help
    page.

    Examples
    --------
    >>> r = lwss(x=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], y=[2, 4, 3, 5, 7, 6, 8, 9, 8, 10])
    >>> len(r)
    10
    """
    return jcall("lwss", x=x, y=y, f=f, iter=iter, delta=delta)


def chmb(a: Any, b: Any, wopt: Any, wval: Any, coarse: Any = 0.25) -> Any:
    """The Chumbley score for two tool-mark profiles.

    The arguments and the result are those of R's ``chmb()``: see its help
    page.

    Examples
    --------
    >>> import math
    >>> a = [math.sin(i / 7) + math.sin(i / 3) for i in range(300)]
    >>> b = a[30:] + a[:30]
    >>> r = chmb(a, b, wopt=60, wval=20)
    >>> list(r)
    ['nsame', 'ndiff', 'u', 'p']
    """
    return jcall("chmb", a=a, b=b, wopt=wopt, wval=wval, coarse=coarse)


def gstr(cnt: Any, ex: Any = None, rad: Any = 1.5, nsim: Any = 0, seed: Any = None) -> Any:
    """Local Getis-Ord Gi* on a grid of counts.

    The arguments and the result are those of R's ``gstr()``: see its help
    page.

    Examples
    --------
    >>> r = gstr(cnt=[[1, 4, 3, 1, 1, 2], [1, 2, 1, 3, 1, 2], [2, 2, 9, 9, 0, 2], [4, 0, 9,
    ... 9, 1, 1], [1, 1, 3, 2, 4, 3], [4, 1, 6, 0, 1, 2]])
    >>> list(r)
    ['row', 'col', 'gz', 'pz', 'qz', 'pp', 'ze', 'pe', 'qe']
    """
    return jcall("gstr", cnt=cnt, ex=ex, rad=rad, nsim=nsim, seed=seed)


def sptst(a: Any, b: Any, nboot: Any = 200, seed: Any = None, alpha: Any = 0.05,
          scale: Any = 100) -> Any:
    """Andresen's spatial point pattern test.

    The arguments and the result are those of R's ``sptst()``: see its help
    page.

    Examples
    --------
    >>> r = sptst(a=[5, 3, 8, 0, 2, 6], b=[4, 6, 7, 1, 0, 9], nboot=50, seed=1)
    >>> list(r)
    ['units', 'glob']
    """
    return jcall("sptst", a=a, b=b, nboot=nboot, seed=seed, alpha=alpha, scale=scale)


def tfrp(par: Any, deffect: Any = True, tdist: Any = "nb", ymax: Any = 20, nsim: Any = 1000,
         seed: Any = None, hl: Any = 0.02, nq: Any = 256) -> Any:
    """Glass fragments recovered under transfer and persistence.

    The arguments and the result are those of R's ``tfrp()``: see its help
    page.

    Examples
    --------
    >>> r = tfrp(par=[80, 50, 0.5, 0.05, 0.8, 0.9, 0.1, 0.2, 0.4, 0.6, 0.05, 0.1, 0.5, 0.7,
    ... 1, 1, 0.5], ymax=10, nsim=100, seed=1)
    >>> list(r)
    ['prob', 'cut']
    """
    return jcall("tfrp", par=par, deffect=deffect, tdist=tdist, ymax=ymax, nsim=nsim, seed=seed,
                 hl=hl, nq=nq)


def zgll(y: Any, x: Any, z: Any = None, par: Any = None, prob: Any = None, off: Any = None) -> Any:
    """Group-based trajectory model for counts.

    The arguments and the result are those of R's ``zgll()``: see its help
    page.

    Examples
    --------
    >>> r = zgll(y=[[3, 2, 2, 3, 0, 2], [8, 4, 2, 0, 0, 0], [5, 2, 3, 1, 1, 3], [4, 4, 4,
    ... 0, 1, 0], [3, 2, 0, 1, 2, 0], [4, 3, 2, 0, 2, 1], [6, 1, 3, 1, 2, 1], [6, 3, 2, 2,
    ... 0, 0], [6, 4, 2, 3, 1, 0], [2, 1, 1, 1, 1, 0], [7, 2, 2, 1, 1, 0], [4, 4, 0, 3, 0,
    ... 1], [4, 2, 2, 3, 2, 0], [6, 6, 2, 2, 0, 0], [6, 2, 3, 0, 0, 0], [1, 1, 5, 8, 8,
    ... 10], [2, 3, 6, 5, 7, 6], [2, 8, 5, 3, 3, 7], [5, 4, 8, 6, 5, 7], [1, 0, 5, 2, 6,
    ... 8], [7, 3, 4, 3, 8, 7], [3, 2, 2, 5, 6, 4], [0, 4, 7, 5, 6, 7], [7, 3, 5, 5, 4, 6],
    ... [4, 3, 2, 6, 3, 10], [3, 3, 3, 5, 6, 5], [3, 1, 3, 3, 5, 11], [2, 5, 6, 4, 3, 6],
    ... [4, 2, 5, 7, 9, 6], [2, 6, 5, 6, 7, 11]], x=[[1, -1], [1, -0.6], [1, -0.2], [1,
    ... 0.2], [1, 0.6], [1, 1]], par=[[0.5, 1.5], [-1, 0.5]], prob=[0.5, 0.5])
    >>> list(r)
    ['ll', 'post']
    """
    return jcall("zgll", y=y, x=x, z=z, par=par, prob=prob, off=off)


def zgem(y: Any, x: Any, z: Any = None, par: Any = None, prob: Any = None, off: Any = None,
         maxit: Any = 500, tol: Any = 1e-8) -> Any:
    """Group-based trajectory model for counts.

    The arguments and the result are those of R's ``zgem()``: see its help
    page.

    Examples
    --------
    >>> r = zgem(y=[[3, 2, 2, 3, 0, 2], [8, 4, 2, 0, 0, 0], [5, 2, 3, 1, 1, 3], [4, 4, 4,
    ... 0, 1, 0], [3, 2, 0, 1, 2, 0], [4, 3, 2, 0, 2, 1], [6, 1, 3, 1, 2, 1], [6, 3, 2, 2,
    ... 0, 0], [6, 4, 2, 3, 1, 0], [2, 1, 1, 1, 1, 0], [7, 2, 2, 1, 1, 0], [4, 4, 0, 3, 0,
    ... 1], [4, 2, 2, 3, 2, 0], [6, 6, 2, 2, 0, 0], [6, 2, 3, 0, 0, 0], [1, 1, 5, 8, 8,
    ... 10], [2, 3, 6, 5, 7, 6], [2, 8, 5, 3, 3, 7], [5, 4, 8, 6, 5, 7], [1, 0, 5, 2, 6,
    ... 8], [7, 3, 4, 3, 8, 7], [3, 2, 2, 5, 6, 4], [0, 4, 7, 5, 6, 7], [7, 3, 5, 5, 4, 6],
    ... [4, 3, 2, 6, 3, 10], [3, 3, 3, 5, 6, 5], [3, 1, 3, 3, 5, 11], [2, 5, 6, 4, 3, 6],
    ... [4, 2, 5, 7, 9, 6], [2, 6, 5, 6, 7, 11]], x=[[1, -1], [1, -0.6], [1, -0.2], [1,
    ... 0.2], [1, 0.6], [1, 1]], par=[[0.7, 1.7], [-0.8, 0.7]], prob=[0.5, 0.5])
    >>> list(r)
    ['par', 'prob', 'post', 'll', 'iter']
    """
    return jcall("zgem", y=y, x=x, z=z, par=par, prob=prob, off=off, maxit=maxit, tol=tol)


def mvtll(y: Any, x: Any, beta: Any, sig: Any, prob: Any) -> Any:
    """Group-based multivariate trajectories.

    The arguments and the result are those of R's ``mvtll()``: see its help
    page.

    Examples
    --------
    >>> r = mvtll(y=[[[2.374, 2.402], [2.379, 0.611], [3.62, 0.594], [2.606, -1.277],
    ... [2.888, 0.384]], [[3.184, -0.039], [0.7850000000000001, -0.934], [2.944, 0.333],
    ... [2.941, -0.573], [3.8810000000000002, 1.682]], [[2.164, 0.69], [4.125, -1.254],
    ... [2.844, 1.063], [4.1, -1.225], [3.398, -0.636]], [[4.595, 0.028], [2.955, 0.291],
    ... [1.529, -0.304], [3.763, -0.473], [2.388, -0.462]], [[3.33, -0.743], [2.984,
    ... -0.443], [2.5220000000000002, 0.37], [2.835, -0.62], [3.341, 1.432]], [[2.18,
    ... 0.189], [3.944, 0.001], [3.418, 0.267], [2.747, 0.042], [1.871, -0.651]], [[0.487,
    ... -1.805], [0.821, 0.074], [1.359, -0.543], [0.697, -0.911], [1.433, -0.207]],
    ... [[0.738, 1.466], [0.594, -0.59], [-0.103, 1.208], [0.557, 0.158], [1.98, -0.393]],
    ... [[0.576, 0.153], [0.919, -0.569], [0.388, 1.16], [-0.689, -0.655], [-0.367,
    ... -0.32]], [[-0.305, 2.173], [0.782, -0.135], [-0.054, 0.7], [-0.707, 1.767],
    ... [-1.044, -0.279]], [[1.512, 0.476], [0.075, 1.178], [-1.377, 1.587], [0.365,
    ... 0.717], [0.57, 0.494]], [[0.39, -0.71], [-1.989, -1.524], [-0.415, 0.558], [0.769,
    ... 0.91], [-0.135, -0.177]]], x=[[1, 1], [1, 2], [1, 3], [1, 4], [1, 5]], beta=[[[0,
    ... 3], [0, 0]], [[0, 0], [0, 0]]], sig=[[[1, 1], [0, 0]], [[0, 0], [1, 1]]],
    ... prob=[0.5, 0.5])
    >>> list(r)
    ['ll', 'post']
    """
    return jcall("mvtll", y=y, x=x, beta=beta, sig=sig, prob=prob)


def mvtem(y: Any, x: Any, beta: Any, sig: Any, prob: Any, maxit: Any = 500, tol: Any = 1e-8) -> Any:
    """Group-based multivariate trajectories.

    The arguments and the result are those of R's ``mvtem()``: see its help
    page.

    Examples
    --------
    >>> r = mvtem(y=[[[2.374, 2.402], [2.379, 0.611], [3.62, 0.594], [2.606, -1.277],
    ... [2.888, 0.384]], [[3.184, -0.039], [0.7850000000000001, -0.934], [2.944, 0.333],
    ... [2.941, -0.573], [3.8810000000000002, 1.682]], [[2.164, 0.69], [4.125, -1.254],
    ... [2.844, 1.063], [4.1, -1.225], [3.398, -0.636]], [[4.595, 0.028], [2.955, 0.291],
    ... [1.529, -0.304], [3.763, -0.473], [2.388, -0.462]], [[3.33, -0.743], [2.984,
    ... -0.443], [2.5220000000000002, 0.37], [2.835, -0.62], [3.341, 1.432]], [[2.18,
    ... 0.189], [3.944, 0.001], [3.418, 0.267], [2.747, 0.042], [1.871, -0.651]], [[0.487,
    ... -1.805], [0.821, 0.074], [1.359, -0.543], [0.697, -0.911], [1.433, -0.207]],
    ... [[0.738, 1.466], [0.594, -0.59], [-0.103, 1.208], [0.557, 0.158], [1.98, -0.393]],
    ... [[0.576, 0.153], [0.919, -0.569], [0.388, 1.16], [-0.689, -0.655], [-0.367,
    ... -0.32]], [[-0.305, 2.173], [0.782, -0.135], [-0.054, 0.7], [-0.707, 1.767],
    ... [-1.044, -0.279]], [[1.512, 0.476], [0.075, 1.178], [-1.377, 1.587], [0.365,
    ... 0.717], [0.57, 0.494]], [[0.39, -0.71], [-1.989, -1.524], [-0.415, 0.558], [0.769,
    ... 0.91], [-0.135, -0.177]]], x=[[1, 1], [1, 2], [1, 3], [1, 4], [1, 5]], beta=[[[0,
    ... 3], [0, 0]], [[0, 0], [0, 0]]], sig=[[[1, 1], [0, 0]], [[0, 0], [1, 1]]],
    ... prob=[0.5, 0.5])
    >>> list(r)
    ['beta', 'sig', 'prob', 'post', 'll', 'iter']
    """
    return jcall("mvtem", y=y, x=x, beta=beta, sig=sig, prob=prob, maxit=maxit, tol=tol)


def ntdeg(from_: Any, to: Any, w: Any = None, nodes: Any = None, dir: Any = False) -> Any:
    """Degree and strength of each node.

    The arguments and the result are those of R's ``ntdeg()``: see its help
    page.

    Examples
    --------
    >>> r = ntdeg(from_=["a", "a", "b", "c", "d", "e", "e"], to=["b", "c", "c", "d", "e",
    ... "f", "d"])
    >>> list(r)
    ['node', 'deg', 'wdeg']
    """
    return jcall("ntdeg", from_=from_, to=to, w=w, nodes=nodes, dir=dir)


def ntbtw(from_: Any, to: Any, nodes: Any = None, dir: Any = False, norm: Any = False) -> Any:
    """Betweenness of each node.

    The arguments and the result are those of R's ``ntbtw()``: see its help
    page.

    Examples
    --------
    >>> r = ntbtw(from_=["a", "a", "b", "c", "d", "e", "e"], to=["b", "c", "c", "d", "e",
    ... "f", "d"])
    >>> list(r)
    ['node', 'btw']
    """
    return jcall("ntbtw", from_=from_, to=to, nodes=nodes, dir=dir, norm=norm)


def ntcls(from_: Any, to: Any, nodes: Any = None, dir: Any = False) -> Any:
    """Closeness of each node.

    The arguments and the result are those of R's ``ntcls()``: see its help
    page.

    Examples
    --------
    >>> r = ntcls(from_=["a", "a", "b", "c", "d", "e", "e"], to=["b", "c", "c", "d", "e",
    ... "f", "d"])
    >>> list(r)
    ['node', 'reach', 'cls', 'ncls', 'hrm']
    """
    return jcall("ntcls", from_=from_, to=to, nodes=nodes, dir=dir)


def ntevc(from_: Any, to: Any, w: Any = None, nodes: Any = None, tol: Any = 1e-13,
          maxit: Any = 100000) -> Any:
    """Eigenvector centrality of each node.

    The arguments and the result are those of R's ``ntevc()``: see its help
    page.

    Examples
    --------
    >>> r = ntevc(from_=["a", "a", "b", "c", "d", "e", "e"], to=["b", "c", "c", "d", "e",
    ... "f", "d"])
    >>> list(r)
    ['cent', 'value', 'iter', 'conv']
    """
    return jcall("ntevc", from_=from_, to=to, w=w, nodes=nodes, tol=tol, maxit=maxit)


def ntcore(from_: Any, to: Any, nodes: Any = None) -> Any:
    """k-core of each node.

    The arguments and the result are those of R's ``ntcore()``: see its help
    page.

    Examples
    --------
    >>> r = ntcore(from_=["a", "a", "b", "c", "d", "e", "e"], to=["b", "c", "c", "d", "e",
    ... "f", "d"])
    >>> list(r)
    ['node', 'core']
    """
    return jcall("ntcore", from_=from_, to=to, nodes=nodes)


def ntcmp(from_: Any, to: Any, nodes: Any = None, dir: Any = False) -> Any:
    """Connected components.

    The arguments and the result are those of R's ``ntcmp()``: see its help
    page.

    Examples
    --------
    >>> r = ntcmp(from_=["a", "a", "b", "c", "d", "e", "e", "x"], to=["b", "c", "c", "d",
    ... "e", "f", "d", "y"], dir=True)
    >>> list(r)
    ['node', 'weak', 'strong']
    """
    return jcall("ntcmp", from_=from_, to=to, nodes=nodes, dir=dir)


def ntlvn(from_: Any, to: Any, w: Any = None, nodes: Any = None, res: Any = 1,
          seed: Any = None) -> Any:
    """Communities by the Louvain method.

    The arguments and the result are those of R's ``ntlvn()``: see its help
    page.

    Examples
    --------
    >>> r = ntlvn(from_=["a", "a", "b", "c", "d", "e", "e"], to=["b", "c", "c", "d", "e",
    ... "f", "d"], seed=1)
    >>> list(r)
    ['memb', 'q', 'levels']
    """
    return jcall("ntlvn", from_=from_, to=to, w=w, nodes=nodes, res=res, seed=seed)


def ntmod(from_: Any, to: Any, memb: Any, w: Any = None, nodes: Any = None, res: Any = 1) -> Any:
    """Modularity of a partition.

    The arguments and the result are those of R's ``ntmod()``: see its help
    page.

    Examples
    --------
    >>> r = ntmod(from_=["a", "a", "b", "c", "d", "e", "e"], to=["b", "c", "c", "d", "e",
    ... "f", "d"], memb=[1, 1, 1, 2, 2, 2], res=0.5)
    >>> round(r, 6)
    0.607143
    """
    return jcall("ntmod", from_=from_, to=to, memb=memb, w=w, nodes=nodes, res=res)


def ntbrk(from_: Any, to: Any, grp: Any, nodes: Any = None, dir: Any = False) -> Any:
    """Brokerage roles.

    The arguments and the result are those of R's ``ntbrk()``: see its help
    page.

    Examples
    --------
    >>> r = ntbrk(from_=["a", "a", "b", "c", "d", "e", "e"], to=["b", "c", "c", "d", "e",
    ... "f", "d"], grp=["x", "x", "y", "y", "z", "z"], dir=True)
    >>> list(r)
    ['node', 'grp', 'coord', 'itin', 'repr', 'gate', 'liais', 'total']
    """
    return jcall("ntbrk", from_=from_, to=to, grp=grp, nodes=nodes, dir=dir)


def ntchg(from_: Any, to: Any, w: Any = None, nodes: Any = None, tol: Any = 1e-10,
          maxit: Any = 100000) -> Any:
    """Cheeger bounds on a network's bottleneck.

    The arguments and the result are those of R's ``ntchg()``: see its help
    page.

    Examples
    --------
    >>> r = ntchg(from_=["a", "a", "b", "c", "d", "e", "e"], to=["b", "c", "c", "d", "e",
    ... "f", "d"])
    >>> list(r)
    ['bound', 'set', 'iter', 'conv']
    """
    return jcall("ntchg", from_=from_, to=to, w=w, nodes=nodes, tol=tol, maxit=maxit)


def coofg(inc: Any, off: Any) -> Any:
    """Co-offending graph from an incident table.

    The arguments and the result are those of R's ``coofg()``: see its help
    page.

    Examples
    --------
    >>> r = coofg(inc=[1, 1, 2, 2, 2, 3, 3], off=["p", "q", "p", "r", "s", "q", "r"])
    >>> list(r)
    ['edges', 'nodes']
    """
    return jcall("coofg", inc=inc, off=off)


def relchk(x: Any, value: Any, type: Any = "or", group: Any = None, ref: Any = None,
           event: Any = None, digits: Any = None, tol: Any = None) -> Any:
    """Which labelling reproduces a reported contrast?.

    The arguments and the result are those of R's ``relchk()``: see its help
    page.

    Examples
    --------
    >>> r = relchk(x={"yes": [30, 20], "no": [70, 80]}, value=1.71)
    >>> list(r)
    ['verdict', 'ok', 'value', 'tol', 'digits', 'type', 'claim', 'reads']
    """
    return jcall("relchk", x=x, value=value, type=type, group=group, ref=ref, event=event,
                 digits=digits, tol=tol)


def catchk(x: Any, ref: Any = None, nmin: Any = 5) -> Any:
    """Audit categorical inputs before an analysis.

    The arguments and the result are those of R's ``catchk()``: see its help
    page.

    Examples
    --------
    >>> r = catchk(x=["u", "v", "v", "w", "w", "w"], ref="v", nmin=2)
    >>> list(r)
    ['ok', 'nmin', 'vars', 'levels']
    """
    return jcall("catchk", x=x, ref=ref, nmin=nmin)


def hvdist(from_: Any, to: Any, r: Any = 6371008.8) -> Any:
    """Great-circle distance.

    The arguments and the result are those of R's ``hvdist()``: see its help
    page.

    Examples
    --------
    >>> r = hvdist(from_=[-79.38, 43.65], to=[[-73.57, 45.5], [-123.12, 49.28]])
    >>> len(r)
    2
    """
    return jcall("hvdist", from_=from_, to=to, r=r)


def eqproj(x: Any, y: Any = None, lon0: Any = None, lat0: Any = None, r: Any = 6371008.8) -> Any:
    """Equal-area projection of longitude and latitude, and its inverse.

    The arguments and the result are those of R's ``eqproj()``: see its help
    page.

    Examples
    --------
    >>> r = eqproj(x={"lon": [-79.4, -79.38, -79.36, -79.39], "lat": [43.64, 43.66, 43.65,
    ... 43.67]})
    >>> list(r)
    ['x', 'y']
    """
    return jcall("eqproj", x=x, y=y, lon0=lon0, lat0=lat0, r=r)


def eqinv(x: Any, y: Any = None, lon0: Any = None, lat0: Any = None, r: Any = 6371008.8) -> Any:
    """Equal-area projection of longitude and latitude, and its inverse.

    The arguments and the result are those of R's ``eqinv()``: see its help
    page.

    Examples
    --------
    >>> r = eqinv(x=[0, 1000, -500], y=[0, 2000, 300], lon0=-79.4, lat0=43.6)
    >>> list(r)
    ['lon', 'lat']
    """
    return jcall("eqinv", x=x, y=y, lon0=lon0, lat0=lat0, r=r)


def cgstat(x: Any, y: Any = None, w: Any = None, k: Any = 1, lonlat: Any = False) -> Any:
    """Mean centre and dispersion of points.

    The arguments and the result are those of R's ``cgstat()``: see its help
    page.

    Examples
    --------
    >>> r = cgstat(x={"x": [0, 1, 2, 3, 4, 2.5], "y": [0, 2, 1, 3, 0.5, 4]})
    >>> list(r)
    ['w', 'mx', 'my', 'sx', 'sy', 'sdd', 'gx', 'gy', 'hx', 'hy']
    """
    return jcall("cgstat", x=x, y=y, w=w, k=k, lonlat=lonlat)


def cgmedn(x: Any, y: Any = None, lonlat: Any = False) -> Any:
    """Median centre of points.

    The arguments and the result are those of R's ``cgmedn()``: see its help
    page.

    Examples
    --------
    >>> r = cgmedn(x={"x": [0, 1, 2, 3, 4, 2.5], "y": [0, 2, 1, 3, 0.5, 4]})
    >>> list(r)
    ['x', 'y']
    """
    return jcall("cgmedn", x=x, y=y, lonlat=lonlat)


def webpt(x: Any, y: Any = None, w: Any = None, start: Any = None, tol: Any = 1e-10,
          maxit: Any = 1000, lonlat: Any = False) -> Any:
    """Centre of minimum distance.

    The arguments and the result are those of R's ``webpt()``: see its help
    page.

    Examples
    --------
    >>> r = webpt(x={"x": [0, 1, 2, 3, 4, 2.5], "y": [0, 2, 1, 3, 0.5, 4]})
    >>> list(r)
    ['x', 'y', 'obj', 'gap', 'iter', 'code']
    """
    return jcall("webpt", x=x, y=y, w=w, start=start, tol=tol, maxit=maxit, lonlat=lonlat)


def sdell(x: Any, y: Any = None, w: Any = None, k: Any = 1, c: Any = 1, lonlat: Any = False) -> Any:
    """Standard deviational ellipse.

    The arguments and the result are those of R's ``sdell()``: see its help
    page.

    Examples
    --------
    >>> r = sdell(x={"x": [0, 1, 2, 3, 4, 2.5], "y": [0, 2, 1, 3, 0.5, 4]})
    >>> list(r)
    ['mx', 'my', 'a', 'b', 'angle', 'area']
    """
    return jcall("sdell", x=x, y=y, w=w, k=k, c=c, lonlat=lonlat)


def cvhull(x: Any, y: Any = None, lonlat: Any = False) -> Any:
    """Convex hull of points.

    The arguments and the result are those of R's ``cvhull()``: see its help
    page.

    Examples
    --------
    >>> r = cvhull(x={"x": [0, 1, 2, 3, 4, 2.5], "y": [0, 2, 1, 3, 0.5, 4]})
    >>> list(r)
    ['idx', 'area', 'perim']
    """
    return jcall("cvhull", x=x, y=y, lonlat=lonlat)


def nnidx(x: Any, y: Any = None, win: Any = None, kmax: Any = 1, nsim: Any = 0, seed: Any = None,
          dist: Any = False, metric: Any = "euclid") -> Any:
    """Nearest-neighbour index.

    The arguments and the result are those of R's ``nnidx()``: see its help
    page.

    Examples
    --------
    >>> r = nnidx(x={"x": [1.849, 7.024, 5.733, 1.681, 9.438, 9.435, 1.292, 8.334, 4.68,
    ... 5.5, 5.527, 2.389, 7.605, 1.808, 4.053, 8.535, 9.764, 2.258, 4.448, 0.75, 6.619,
    ... 3.875, 8.369, 1.505, 3.473, 4.888, 1.492, 3.571, 9.626, 1.324], "y": [0.104, 1.646,
    ... 8.102, 8.689, 5.143, 6.272, 8.444, 2.849, 6.672, 1.505, 9.817, 2.97, 1.151, 1.632,
    ... 9.44, 7.949, 9.747, 3.491, 5.02, 8.104, 0.071, 0.147, 6.834, 9.297, 2.754, 8.119,
    ... 7.859, 9.889, 6.14, 7.102]}, win=[0, 10, 0, 10])
    >>> list(r)
    ['order', 'ce', 'dist']
    """
    return jcall("nnidx", x=x, y=y, win=win, kmax=kmax, nsim=nsim, seed=seed, dist=dist,
                 metric=metric)


def ripk(x: Any, y: Any = None, win: Any = None, r: Any = None, corr: Any = "iso",
         metric: Any = "euclid") -> Any:
    """Ripley's K function.

    The arguments and the result are those of R's ``ripk()``: see its help
    page.

    Examples
    --------
    >>> r = ripk(x={"x": [1.849, 7.024, 5.733, 1.681, 9.438, 9.435, 1.292, 8.334, 4.68,
    ... 5.5, 5.527, 2.389, 7.605, 1.808, 4.053, 8.535, 9.764, 2.258, 4.448, 0.75, 6.619,
    ... 3.875, 8.369, 1.505, 3.473, 4.888, 1.492, 3.571, 9.626, 1.324], "y": [0.104, 1.646,
    ... 8.102, 8.689, 5.143, 6.272, 8.444, 2.849, 6.672, 1.505, 9.817, 2.97, 1.151, 1.632,
    ... 9.44, 7.949, 9.747, 3.491, 5.02, 8.104, 0.071, 0.147, 6.834, 9.297, 2.754, 8.119,
    ... 7.859, 9.889, 6.14, 7.102]}, win=[0, 10, 0, 10], r=[0.5, 1, 2])
    >>> list(r)
    ['r', 'k', 'l']
    """
    return jcall("ripk", x=x, y=y, win=win, r=r, corr=corr, metric=metric)


def ripenv(x: Any, y: Any = None, win: Any = None, r: Any = None, corr: Any = "iso",
           nsim: Any = 99, seed: Any = None) -> Any:
    """Envelopes of Ripley's L function and a test of complete spatial randomness.

    The arguments and the result are those of R's ``ripenv()``: see its help
    page.

    Examples
    --------
    >>> r = ripenv(x={"x": [1.849, 7.024, 5.733, 1.681, 9.438, 9.435, 1.292, 8.334, 4.68,
    ... 5.5, 5.527, 2.389, 7.605, 1.808, 4.053, 8.535, 9.764, 2.258, 4.448, 0.75, 6.619,
    ... 3.875, 8.369, 1.505, 3.473, 4.888, 1.492, 3.571, 9.626, 1.324], "y": [0.104, 1.646,
    ... 8.102, 8.689, 5.143, 6.272, 8.444, 2.849, 6.672, 1.505, 9.817, 2.97, 1.151, 1.632,
    ... 9.44, 7.949, 9.747, 3.491, 5.02, 8.104, 0.071, 0.147, 6.834, 9.297, 2.754, 8.119,
    ... 7.859, 9.889, 6.14, 7.102]}, win=[0, 10, 0, 10], r=[0.5, 1, 2], nsim=19, seed=1)
    >>> list(r)
    ['env', 'test']
    """
    return jcall("ripenv", x=x, y=y, win=win, r=r, corr=corr, nsim=nsim, seed=seed)


def fzmode(x: Any, y: Any = None, r: Any = None, count: Any = None, metric: Any = "euclid",
           lonlat: Any = False) -> Any:
    """Fuzzy mode.

    The arguments and the result are those of R's ``fzmode()``: see its help
    page.

    Examples
    --------
    >>> r = fzmode(x={"x": [1.849, 7.024, 5.733, 1.681, 9.438, 9.435, 1.292, 8.334, 4.68,
    ... 5.5, 5.527, 2.389, 7.605, 1.808, 4.053, 8.535, 9.764, 2.258, 4.448, 0.75, 6.619,
    ... 3.875, 8.369, 1.505, 3.473, 4.888, 1.492, 3.571, 9.626, 1.324], "y": [0.104, 1.646,
    ... 8.102, 8.689, 5.143, 6.272, 8.444, 2.849, 6.672, 1.505, 9.817, 2.97, 1.151, 1.632,
    ... 9.44, 7.949, 9.747, 3.491, 5.02, 8.104, 0.071, 0.147, 6.834, 9.297, 2.754, 8.119,
    ... 7.859, 9.889, 6.14, 7.102]}, r=2)
    >>> len(r)
    30
    """
    return jcall("fzmode", x=x, y=y, r=r, count=count, metric=metric, lonlat=lonlat)


def kmpart(x: Any, y: Any = None, k: Any = None, init: Any = "pp", cell: Any = None,
           upd: Any = "mean", maxit: Any = 100, seed: Any = None, count: Any = None,
           lonlat: Any = False) -> Any:
    """K-means partition of points.

    The arguments and the result are those of R's ``kmpart()``: see its help
    page.

    Examples
    --------
    >>> r = kmpart(x={"x": [1.849, 7.024, 5.733, 1.681, 9.438, 9.435, 1.292, 8.334, 4.68,
    ... 5.5, 5.527, 2.389, 7.605, 1.808, 4.053, 8.535, 9.764, 2.258, 4.448, 0.75, 6.619,
    ... 3.875, 8.369, 1.505, 3.473, 4.888, 1.492, 3.571, 9.626, 1.324], "y": [0.104, 1.646,
    ... 8.102, 8.689, 5.143, 6.272, 8.444, 2.849, 6.672, 1.505, 9.817, 2.97, 1.151, 1.632,
    ... 9.44, 7.949, 9.747, 3.491, 5.02, 8.104, 0.071, 0.147, 6.834, 9.297, 2.754, 8.119,
    ... 7.859, 9.889, 6.14, 7.102]}, k=3, seed=1)
    >>> list(r)
    ['cl', 'cen', 'fit']
    """
    return jcall("kmpart", x=x, y=y, k=k, init=init, cell=cell, upd=upd, maxit=maxit, seed=seed,
                 count=count, lonlat=lonlat)


def nnhclu(x: Any, y: Any = None, area: Any = None, fixd: Any = None, p: Any = 0.5,
           minpts: Any = 10, minhi: Any = None, maxord: Any = 3, k: Any = 1, c: Any = 1,
           count: Any = None, lonlat: Any = False) -> Any:
    """Nearest-neighbour hierarchical clustering.

    The arguments and the result are those of R's ``nnhclu()``: see its help
    page.

    Examples
    --------
    >>> r = nnhclu(x={"x": [2.065, 1.837, 2.267, 2.179, 2.491, 2.207, 1.616, 1.936, 2.569,
    ... 2.533, 2.17, 2.005, 2.115, 1.986, 2.01, 2.051, 2.35, 1.987, 1.97, 1.915, 2.462,
    ... 2.05, 2.392, 2.386, 2.178, 3.886, 7.459, 8.954, 8.091, 8.185, 4.212, 1.767, 1.738,
    ... 8.925, 7.433, 5.61, 0.712, 8.536, 9.137, 2.254, 6.286, 0.691, 5.149, 8.054, 9.718,
    ... 3.429, 6.314, 4.101, 3.457, 8.249], "y": [2.149, 1.956, 1.667, 2.165, 2.415, 1.53,
    ... 2.685, 2.441, 1.658, 1.818, 1.587, 1.958, 2.158, 2.099, 2.024, 2.056, 1.76, 2.263,
    ... 1.923, 1.995, 2.008, 2.243, 1.56, 1.919, 1.811, 7.243, 9.805, 2.006, 5.523, 3.455,
    ... 4.084, 1.88, 3.456, 1.837, 8.131, 8.816, 7.006, 2.764, 4.32, 0.128, 5.63, 6.707,
    ... 3.136, 9.475, 1.688, 4.404, 5.716, 9.747, 7.648, 4.431]}, area=100, minpts=5)
    >>> list(r)
    ['cl', 'clusters']
    """
    return jcall("nnhclu", x=x, y=y, area=area, fixd=fixd, p=p, minpts=minpts, minhi=minhi,
                 maxord=maxord, k=k, c=c, count=count, lonlat=lonlat)


def stacel(x: Any, y: Any = None, win: Any = None, nx: Any = 20, ny: Any = 20, r: Any = None,
           top: Any = 25, minpts: Any = 5, k: Any = 1, c: Any = 1, count: Any = None) -> Any:
    """Search-circle hot-spot areas.

    The arguments and the result are those of R's ``stacel()``: see its help
    page.

    Examples
    --------
    >>> r = stacel(x={"x": [2.065, 1.837, 2.267, 2.179, 2.491, 2.207, 1.616, 1.936, 2.569,
    ... 2.533, 2.17, 2.005, 2.115, 1.986, 2.01, 2.051, 2.35, 1.987, 1.97, 1.915, 2.462,
    ... 2.05, 2.392, 2.386, 2.178, 3.886, 7.459, 8.954, 8.091, 8.185, 4.212, 1.767, 1.738,
    ... 8.925, 7.433, 5.61, 0.712, 8.536, 9.137, 2.254, 6.286, 0.691, 5.149, 8.054, 9.718,
    ... 3.429, 6.314, 4.101, 3.457, 8.249], "y": [2.149, 1.956, 1.667, 2.165, 2.415, 1.53,
    ... 2.685, 2.441, 1.658, 1.818, 1.587, 1.958, 2.158, 2.099, 2.024, 2.056, 1.76, 2.263,
    ... 1.923, 1.995, 2.008, 2.243, 1.56, 1.919, 1.811, 7.243, 9.805, 2.006, 5.523, 3.455,
    ... 4.084, 1.88, 3.456, 1.837, 8.131, 8.816, 7.006, 2.764, 4.32, 0.128, 5.63, 6.707,
    ... 3.136, 9.475, 1.688, 4.404, 5.716, 9.747, 7.648, 4.431]}, win=[0, 10, 0, 10], r=1,
    ... top=5, minpts=4)
    >>> list(r)
    ['cl', 'clusters']
    """
    return jcall("stacel", x=x, y=y, win=win, nx=nx, ny=ny, r=r, top=top, minpts=minpts, k=k, c=c,
                 count=count)


def mantst(x: Any, y: Any = None, t: Any = None, cs: Any = 0, ct: Any = 0, nsim: Any = 999,
           seed: Any = None, metric: Any = "euclid", lonlat: Any = False) -> Any:
    """Mantel's space-time test.

    The arguments and the result are those of R's ``mantst()``: see its help
    page.

    Examples
    --------
    >>> r = mantst(x={"x": [1.001, 3.426, 4.584, 1.422, 0.523, 3.505, 2.64, 4.04, 4.783,
    ... 0.552, 1.366, 2.453, 1.592, 2.796, 1.313, 1.009, 1.938, 4.439, 2.775, 4.211],
    ... "y": [4.451, 3.604, 1.057, 1.129, 0.7, 2.4, 2.187, 4.83, 0.71, 4.775, 2.224, 0.297,
    ... 1.376, 0.156, 0.072, 2.436, 2.976, 2.989, 1.988, 1.984]}, t=[2.06, 4.09, 5.91,
    ... 7.06, 7.67, 7.72, 11.43, 11.54, 16, 16.49, 16.99, 18.98, 22.8, 23.81, 24.48, 24.81,
    ... 27.66, 27.87, 27.99, 29.28], nsim=49, seed=1)
    >>> list(r)
    ['r', 'pairs', 'pup', 'plo', 'p2']
    """
    return jcall("mantst", x=x, y=y, t=t, cs=cs, ct=ct, nsim=nsim, seed=seed, metric=metric,
                 lonlat=lonlat)


def stmavg(x: Any, y: Any = None, t: Any = None, span: Any = 5, lonlat: Any = False) -> Any:
    """Spatial-temporal moving average.

    The arguments and the result are those of R's ``stmavg()``: see its help
    page.

    Examples
    --------
    >>> r = stmavg(x={"x": [1.001, 3.426, 4.584, 1.422, 0.523, 3.505, 2.64, 4.04, 4.783,
    ... 0.552, 1.366, 2.453, 1.592, 2.796, 1.313, 1.009, 1.938, 4.439, 2.775, 4.211],
    ... "y": [4.451, 3.604, 1.057, 1.129, 0.7, 2.4, 2.187, 4.83, 0.71, 4.775, 2.224, 0.297,
    ... 1.376, 0.156, 0.072, 2.436, 2.976, 2.989, 1.988, 1.984]}, t=[2.06, 4.09, 5.91,
    ... 7.06, 7.67, 7.72, 11.43, 11.54, 16, 16.49, 16.99, 18.98, 22.8, 23.81, 24.48, 24.81,
    ... 27.66, 27.87, 27.99, 29.28])
    >>> list(r)
    ['mx', 'my', 'n']
    """
    return jcall("stmavg", x=x, y=y, t=t, span=span, lonlat=lonlat)


def jtdecay(d: Any, type: Any, par: Any, tab: Any = None) -> Any:
    """Distance-decay functions for journey to crime.

    The arguments and the result are those of R's ``jtdecay()``: see its help
    page.

    Examples
    --------
    >>> r = jtdecay(d=[0, 0.5, 1, 2, 5], type="negexp", par=[1, 0.5])
    >>> len(r)
    5
    """
    return jcall("jtdecay", d=d, type=type, par=par, tab=tab)


def jtsurf(x: Any, y: Any = None, win: Any = None, nx: Any = None, ny: Any = None,
           type: Any = None, par: Any = None, tab: Any = None, metric: Any = "euclid",
           mode: Any = "sum") -> Any:
    """Journey-to-crime surface.

    The arguments and the result are those of R's ``jtsurf()``: see its help
    page.

    Examples
    --------
    >>> r = jtsurf(x=[2, 4, 3], y=[3, 1, 4], win=[0, 6, 0, 6], nx=4, ny=3, type="negexp",
    ... par=[1, 0.4])
    >>> list(r)
    ['raw', 'prob', 'peak']
    """
    return jcall("jtsurf", x=x, y=y, win=win, nx=nx, ny=ny, type=type, par=par, tab=tab,
                 metric=metric, mode=mode)


def jtcal(d: Any, t: Any, h: Any, w: Any = None, kern: Any = "normal", refl: Any = True) -> Any:
    """Empirical distance decay from solved trips.

    The arguments and the result are those of R's ``jtcal()``: see its help
    page.

    Examples
    --------
    >>> r = jtcal(d=[0.5, 1, 1.5, 2, 3, 4], t=[0, 1, 2, 3], h=0.8)
    >>> len(r)
    4
    """
    return jcall("jtcal", d=d, t=t, h=h, w=w, kern=kern, refl=refl)


def jtacc(prob: Any, win: Any, base: Any) -> Any:
    """Accuracy of a journey-to-crime surface on a solved series.

    The arguments and the result are those of R's ``jtacc()``: see its help
    page.

    Examples
    --------
    >>> r = jtacc(prob=[[0.06148360186816512, 0.08318485708583567, 0.061836270729800615],
    ... [0.08569135054028787, 0.12900129640940466, 0.08527923167686345],
    ... [0.10702908625632654, 0.10766356069336598, 0.08012021156147313],
    ... [0.07273905502291408, 0.0719294039310867, 0.054042074224476176]], win=[0, 6, 0, 6],
    ... base=[3, 2])
    >>> list(r)
    ['p', 'hitbest', 'hitworst', 'err', 'peak', 'cell']
    """
    return jcall("jtacc", prob=prob, win=win, base=base)


def bjtc(pj: Any, od: Any, dest: Any) -> Any:
    """Bayesian journey to crime.

    The arguments and the result are those of R's ``bjtc()``: see its help
    page.

    Examples
    --------
    >>> r = bjtc(pj=[0.2, 0.5, 0.3], od=[[2, 1], [1, 3], [0, 2]], dest=[1, 2])
    >>> list(r)
    ['pj', 'po', 'pod', 'prod', 'risk']
    """
    return jcall("bjtc", pj=pj, od=od, dest=dest)


def spwt(x: Any = None, y: Any = None, type: Any = "knn", k: Any = 4, d: Any = None,
         nb: Any = None, idw: Any = 0, style: Any = "row", metric: Any = "euclid",
         lonlat: Any = False) -> Any:
    """Spatial weights.

    The arguments and the result are those of R's ``spwt()``: see its help
    page.

    Examples
    --------
    >>> r = spwt(type="contig", nb=[2, [1, 3], [2, 4, 5], 3, 3])
    >>> list(r)
    ['n', 'i', 'j', 'd', 'w']
    """
    return jcall("spwt", x=x, y=y, type=type, k=k, d=d, nb=nb, idw=idw, style=style, metric=metric,
                 lonlat=lonlat)


def gearc(x: Any, w: Any, nsim: Any = 0, seed: Any = None) -> Any:
    """Geary's C.

    The arguments and the result are those of R's ``gearc()``: see its help
    page.

    Examples
    --------
    >>> pts = {"x": [0, 1, 2, 3, 4, 5], "y": [0, 1, 0, 1, 0, 1]}
    >>> w = spwt(pts, k=2)
    >>> r = gearc([1, 3, 2, 5, 4, 6], w)
    >>> list(r)
    ['c', 'vn', 'vr', 'zn', 'zr', 'pn', 'pr', 'pperm', 's0']
    """
    return jcall("gearc", x=x, w=w, nsim=nsim, seed=seed)


def gorg(x: Any, w: Any, nsim: Any = 0, seed: Any = None) -> Any:
    """The general G of Getis and Ord.

    The arguments and the result are those of R's ``gorg()``: see its help
    page.

    Examples
    --------
    >>> pts = {"x": [0, 1, 2, 3, 4, 5], "y": [0, 1, 0, 1, 0, 1]}
    >>> w = spwt(pts, k=2)
    >>> r = gorg([1, 3, 2, 5, 4, 6], w)
    >>> list(r)
    ['g', 'e', 'v', 'z', 'p', 'pperm']
    """
    return jcall("gorg", x=x, w=w, nsim=nsim, seed=seed)


def lmorn(x: Any, w: Any, nsim: Any = 0, seed: Any = None) -> Any:
    """Local Moran (LISA).

    The arguments and the result are those of R's ``lmorn()``: see its help
    page.

    Examples
    --------
    >>> pts = {"x": [0, 1, 2, 3, 4, 5], "y": [0, 1, 0, 1, 0, 1]}
    >>> w = spwt(pts, k=2)
    >>> r = lmorn([1, 3, 2, 5, 4, 6], w)
    >>> list(r)
    ['ii', 'e', 'v', 'z', 'p', 'pperm', 'quad']
    """
    return jcall("lmorn", x=x, w=w, nsim=nsim, seed=seed)


def locgi(x: Any, w: Any, star: Any = True) -> Any:
    """Local G and G* of Getis and Ord.

    The arguments and the result are those of R's ``locgi()``: see its help
    page.

    Examples
    --------
    >>> pts = {"x": [0, 1, 2, 3, 4, 5], "y": [0, 1, 0, 1, 0, 1]}
    >>> w = spwt(pts, k=2)
    >>> r = locgi([1, 3, 2, 5, 4, 6], w)
    >>> len(r)
    6
    """
    return jcall("locgi", x=x, w=w, star=star)


def dirmn(x: Any = None, y: Any = None, origin: Any = None, w: Any = None, angle: Any = None,
          lonlat: Any = False) -> Any:
    """Directional mean and circular variance.

    The arguments and the result are those of R's ``dirmn()``: see its help
    page.

    Examples
    --------
    >>> r = dirmn(angle=[0.1, 0.4, -0.2, 1])
    >>> list(r)
    ['theta', 'bearing', 'r', 'var', 'sd', 'w', 'n']
    """
    return jcall("dirmn", x=x, y=y, origin=origin, w=w, angle=angle, lonlat=lonlat)


def hdbang(x: Any, y: Any = None, v: Any = None, w: Any = None, nn: Any = 8, ang: Any = 135,
           ntrip: Any = 0, iter: Any = 1, lonlat: Any = False) -> Any:
    """Head-bang smoothing.

    The arguments and the result are those of R's ``hdbang()``: see its help
    page.

    Examples
    --------
    >>> r = hdbang(x={"x": [1.849, 7.024, 5.733, 1.681, 9.438, 9.435, 1.292, 8.334, 4.68,
    ... 5.5, 5.527, 2.389, 7.605, 1.808, 4.053, 8.535, 9.764, 2.258, 4.448, 0.75, 6.619,
    ... 3.875, 8.369, 1.505, 3.473, 4.888, 1.492, 3.571, 9.626, 1.324], "y": [0.104, 1.646,
    ... 8.102, 8.689, 5.143, 6.272, 8.444, 2.849, 6.672, 1.505, 9.817, 2.97, 1.151, 1.632,
    ... 9.44, 7.949, 9.747, 3.491, 5.02, 8.104, 0.071, 0.147, 6.834, 9.297, 2.754, 8.119,
    ... 7.859, 9.889, 6.14, 7.102]}, v=[1.901, 7.847, 9.784, 6.026, 12.01, 12.571, 5.514,
    ... 9.758, 8.016, 6.252, 10.436, 3.874, 8.181, 2.624, 8.773, 12.509, 14.637, 4.004,
    ... 6.958, 4.802, 6.654, 3.949, 11.786, 6.154, 4.85, 8.948, 5.422, 8.515, 12.696,
    ... 4.875], nn=5)
    >>> len(r)
    30
    """
    return jcall("hdbang", x=x, y=y, v=v, w=w, nn=nn, ang=ang, ntrip=ntrip, iter=iter,
                 lonlat=lonlat)


def crwalk(x: Any, y: Any = None, h: Any = 1, nsim: Any = 999, seed: Any = None,
           lonlat: Any = False) -> Any:
    """Correlated random walk of a series of events.

    The arguments and the result are those of R's ``crwalk()``: see its help
    page.

    Examples
    --------
    >>> r = crwalk(x={"x": [0.27, -0.36, 0.508, 2.235, 2.26, 2.628, 1.318, 2.057, 2.102,
    ... 1.054, 2.781, 1.603], "y": [0.653, 0.285, -0.315, -0.26, 1.447, 0.353, 0.064,
    ... 2.271, 2.79, 1.385, 3.4, 2.212]}, h=3, nsim=9, seed=1)
    >>> list(r)
    ['fit', 'msd', 'sim']
    """
    return jcall("crwalk", x=x, y=y, h=h, nsim=nsim, seed=seed, lonlat=lonlat)


def rnnhcl(x: Any, y: Any = None, base: Any = None, bw: Any = None, h: Any = None, p: Any = 0.5,
           minpts: Any = 10, k: Any = 1, c: Any = 1, lonlat: Any = False) -> Any:
    """Risk-adjusted nearest-neighbour hierarchical clustering.

    The arguments and the result are those of R's ``rnnhcl()``: see its help
    page.

    Examples
    --------
    >>> r = rnnhcl(x={"x": [2.065, 1.837, 2.267, 2.179, 2.491, 2.207, 1.616, 1.936, 2.569,
    ... 2.533, 2.17, 2.005, 2.115, 1.986, 2.01, 2.051, 2.35, 1.987, 1.97, 1.915, 2.462,
    ... 2.05, 2.392, 2.386, 2.178, 3.886, 7.459, 8.954, 8.091, 8.185, 4.212, 1.767, 1.738,
    ... 8.925, 7.433, 5.61, 0.712, 8.536, 9.137, 2.254, 6.286, 0.691, 5.149, 8.054, 9.718,
    ... 3.429, 6.314, 4.101, 3.457, 8.249], "y": [2.149, 1.956, 1.667, 2.165, 2.415, 1.53,
    ... 2.685, 2.441, 1.658, 1.818, 1.587, 1.958, 2.158, 2.099, 2.024, 2.056, 1.76, 2.263,
    ... 1.923, 1.995, 2.008, 2.243, 1.56, 1.919, 1.811, 7.243, 9.805, 2.006, 5.523, 3.455,
    ... 4.084, 1.88, 3.456, 1.837, 8.131, 8.816, 7.006, 2.764, 4.32, 0.128, 5.63, 6.707,
    ... 3.136, 9.475, 1.688, 4.404, 5.716, 9.747, 7.648, 4.431]}, base={"x": [1.849, 7.024,
    ... 5.733, 1.681, 9.438, 9.435, 1.292, 8.334, 4.68, 5.5, 5.527, 2.389, 7.605, 1.808,
    ... 4.053, 8.535, 9.764, 2.258, 4.448, 0.75, 6.619, 3.875, 8.369, 1.505, 3.473, 4.888,
    ... 1.492, 3.571, 9.626, 1.324], "y": [0.104, 1.646, 8.102, 8.689, 5.143, 6.272, 8.444,
    ... 2.849, 6.672, 1.505, 9.817, 2.97, 1.151, 1.632, 9.44, 7.949, 9.747, 3.491, 5.02,
    ... 8.104, 0.071, 0.147, 6.834, 9.297, 2.754, 8.119, 7.859, 9.889, 6.14, 7.102]}, h=2,
    ... minpts=5)
    >>> list(r)
    ['cl', 'clusters', 'thr']
    """
    return jcall("rnnhcl", x=x, y=y, base=base, bw=bw, h=h, p=p, minpts=minpts, k=k, c=c,
                 lonlat=lonlat)


def txcal(conc: Any, resp: Any, w: Any = "1/x^2", unk: Any = None, k: Any = None) -> Any:
    """Calibration line of an assay.

    The arguments and the result are those of R's ``txcal()``: see its help
    page.

    Examples
    --------
    >>> r = txcal(conc=[1, 2, 5, 10, 20], resp=[0.11, 0.2, 0.52, 0.98, 2.05], w="none",
    ... k=[3, 9])
    >>> list(r)
    ['fit', 'unk', 'w']
    """
    return jcall("txcal", conc=conc, resp=resp, w=w, unk=unk, k=k)


def lodimp(x: Any, lod: Any, meth: Any = "half") -> Any:
    """Values below the limit of detection.

    The arguments and the result are those of R's ``lodimp()``: see its help
    page.

    Examples
    --------
    >>> r = lodimp(x=[0.4, 0.05, None, 1.2, 0.08], lod=0.1)
    >>> list(r)
    ['tab', 'frac', 'nmiss']
    """
    return jcall("lodimp", x=x, lod=lod, meth=meth)


def pmrrat(cen: Any, per: Any, cut: Any = None) -> Any:
    """Postmortem redistribution ratio.

    The arguments and the result are those of R's ``pmrrat()``: see its help
    page.

    Examples
    --------
    >>> r = pmrrat(cen=[0.5, 1, 2], per=[0.4, 0.6, 0.5], cut=[1.5, 3])
    >>> list(r)
    ['cen', 'per', 'ratio', 'redist']
    """
    return jcall("pmrrat", cen=cen, per=per, cut=cut)


def etohlr(eth: Any, prop: Any = 0, but: Any = 0, etg: Any = None, thr: Any = 0) -> Any:
    """Ethanol in a decomposed body: drunk before death or made after.

    The arguments and the result are those of R's ``etohlr()``: see its help
    page.

    Examples
    --------
    >>> r = etohlr(eth=[0.5, 0.1])
    >>> list(r)
    ['eth', 'prop', 'but', 'etg', 'higher', 'verdict']
    """
    return jcall("etohlr", eth=eth, prop=prop, but=but, etg=etg, thr=thr)


def mtxrel(base: Any, dil: Any, days: Any = 0, stage: Any = 0, dmax: Any = 30, dcap: Any = 0.8,
           smax: Any = 10, scap: Any = 0.6) -> Any:
    """Reliability of specimen matrices in a submerged or decomposed body.

    The arguments and the result are those of R's ``mtxrel()``: see its help
    page.

    Examples
    --------
    >>> r = mtxrel(base={"blood": 0.9, "vitreous": 0.8, "urine": 0.6}, dil=[True, False,
    ... False], days=3, stage=2)
    >>> list(r)
    ['matrix', 'base', 'rel', 'dil', 'rank']
    """
    return jcall("mtxrel", base=base, dil=dil, days=days, stage=stage, dmax=dmax, dcap=dcap,
                 smax=smax, scap=scap)


def detrsp(x: Any, ben: Any, sanc: Any, p: Any) -> Any:
    """Offending chosen under certainty and severity of sanction.

    The arguments and the result are those of R's ``detrsp()``: see its help
    page.

    Examples
    --------
    >>> r = detrsp(x=[0, 1, 2, 3], ben=[0, 1.5, 2.5, 3], sanc=[0, 1, 3, 6], p=[0.1, 0.5, 0.9])
    >>> list(r)
    ['p', 'xopt', 'value']
    """
    return jcall("detrsp", x=x, ben=ben, sanc=sanc, p=p)


def dtrid(p: Any, s: Any, c: Any, y: Any = None, log: Any = False, rho: Any = None) -> Any:
    """Which deterrence effects a design identifies.

    The arguments and the result are those of R's ``dtrid()``: see its help
    page.

    Examples
    --------
    >>> r = dtrid(p=[0.1, 0.2, 0.3, 0.4, 0.2, 0.3], s=[1, 2, 1, 3, 2, 2], c=[0.5, 1, 1.5,
    ... 0.5, 2, 1])
    >>> list(r)
    ['rank', 'ident', 'const', 'coef', 'disc', 'cost']
    """
    return jcall("dtrid", p=p, s=s, c=c, y=y, log=log, rho=rho)


def inspg(a: Any, x: Any, steps: Any = 1000, dt: Any = 0.01) -> Any:
    """The criminal, citizen and inspector game.

    The arguments and the result are those of R's ``inspg()``: see its help
    page.

    Examples
    --------
    >>> r = inspg(a=0.5, x=[0.5, 0.3, 0.2], steps=5, dt=0.1)
    >>> list(r)
    ['traj', 'inv', 'drift', 'dvdt', 'field', 'fixed']
    """
    return jcall("inspg", a=a, x=x, steps=steps, dt=dt)


def insphs(temp: Any = None, cost: Any = None, steps: Any = 600, rate: Any = 0.05,
           x0: Any = None) -> Any:
    """Phase diagram of the inspection game.

    The arguments and the result are those of R's ``insphs()``: see its help
    page.

    Examples
    --------
    >>> r = insphs(temp=[0.1, 0.5], cost=[0.2, 0.8, 1], steps=20)
    >>> list(r)
    ['crime', 'coop', 'insp', 'temp', 'cost']
    """
    return jcall("insphs", temp=temp, cost=cost, steps=steps, rate=rate, x0=x0)


def colarr(a: Any, e: Any, bg: Any) -> Any:
    """Collider bias of arrest data.

    The arguments and the result are those of R's ``colarr()``: see its help
    page.

    Examples
    --------
    >>> r = colarr(a=[0.2, 0.5], e=0.3, bg=[0, 0.5])
    >>> list(r)
    ['a', 'e', 'bg', 'or', 'lor', 'c11', 'c10', 'c01', 'c00']
    """
    return jcall("colarr", a=a, e=e, bg=bg)


def crsev(x: Any, off: Any, w: Any, unit: Any = None, time: Any = None, pop: Any = None,
          per: Any = 1e5, base: Any = None, sent: Any = None) -> Any:
    """Severity-weighted crime index.

    The arguments and the result are those of R's ``crsev()``: see its help
    page.

    Examples
    --------
    >>> r = crsev(x=[2, 1, 3, 1], off=["theft", "assault", "theft", "robbery"],
    ... w={"assault": 5, "robbery": 4, "theft": 1})
    >>> list(r)
    ['tab', 'w', 'per', 'base']
    """
    return jcall("crsev", x=x, off=off, w=w, unit=unit, time=time, pop=pop, per=per, base=base,
                 sent=sent)


def zcomp(x: Any, unit: Any, cat: Any, w: Any = None) -> Any:
    """Composite index of standardised counts.

    The arguments and the result are those of R's ``zcomp()``: see its help
    page.

    Examples
    --------
    >>> r = zcomp(x=[3, 1, 4, 1, 5, 9, 2, 6], unit=["p", "q", "r", "s", "p", "q", "r",
    ... "s"], cat=["x", "x", "x", "x", "y", "y", "y", "y"])
    >>> list(r)
    ['tab', 'cat', 'z']
    """
    return jcall("zcomp", x=x, unit=unit, cat=cat, w=w)


def cconc(x: Any, w: Any = None, by: Any = None, fracs: Any = None, pop: Any = None,
          nsim: Any = 2000, seed: Any = None) -> Any:
    """Concentration of crime across places, people or units.

    The arguments and the result are those of R's ``cconc()``: see its help
    page.

    Examples
    --------
    >>> r = cconc(x=[0, 0, 1, 0, 3, 7, 0, 1, 2, 0, 12, 1], nsim=0)
    >>> list(r)
    ['stat', 'top', 'lorenz', 'null', 'by', 'scale', 'grow', 'places']
    """
    return jcall("cconc", x=x, w=w, by=by, fracs=fracs, pop=pop, nsim=nsim, seed=seed)


def hilltl(x: Any, xmin: Any = None, discrete: Any = True, approx: Any = False,
           mintl: Any = 3) -> Any:
    """Tail index of a heavy-tailed count.

    The arguments and the result are those of R's ``hilltl()``: see its help
    page.

    Examples
    --------
    >>> r = hilltl(x=[1, 2, 3], mintl=5)
    >>> list(r)
    ['alpha', 'se', 'xmin', 'ntl', 'ks', 'reliable', 'method']
    """
    return jcall("hilltl", x=x, xmin=xmin, discrete=discrete, approx=approx, mintl=mintl)


def crcap(x: Any, kappa: Any = 1, int: Any = None) -> Any:
    """Capture-recapture estimates of a hidden population.

    The arguments and the result are those of R's ``crcap()``: see its help
    page.

    Examples
    --------
    >>> r = crcap(x={"n1": 120, "n2": 90, "m": 30})
    >>> list(r)
    ['type', 'est', 'pairs', 'sens']
    """
    return jcall("crcap", x=x, kappa=kappa, int=int)


def dkfig(v: Any = None, r: Any = None, amax: Any = 0, bmax: Any = 0, thr: Any = None,
          va: Any = None, ra: Any = None, rb: Any = None, lam: Any = 1, k: Any = None) -> Any:
    """The dark figure of crime: bounds instead of a point.

    The arguments and the result are those of R's ``dkfig()``: see its help
    page.

    Examples
    --------
    >>> r = dkfig(va=0.2, ra=0.05, rb=0.02, lam=0.5, k=[1, 1, 2, 3])
    >>> list(r)
    ['bnd', 'cmp', 'hier']
    """
    return jcall("dkfig", v=v, r=r, amax=amax, bmax=bmax, thr=thr, va=va, ra=ra, rb=rb, lam=lam,
                 k=k)


def bndxw(from_: Any, to: Any, x: Any, unit: Any, w: Any = None) -> Any:
    """Crosswalk between boundary versions.

    The arguments and the result are those of R's ``bndxw()``: see its help
    page.

    Examples
    --------
    >>> r = bndxw(from_=["a", "b"], to=["X", "X"], x=[10, 20], unit=["a", "b"])
    >>> list(r)
    ['tab', 'map', 'lost', 'unmatched']
    """
    return jcall("bndxw", from_=from_, to=to, x=x, unit=unit, w=w)


def synexp(y: Any, x: Any, ax: Any, pop: Any) -> Any:
    """Synthetic exposure of areas from a survey.

    The arguments and the result are those of R's ``synexp()``: see its help
    page.

    Examples
    --------
    >>> r = synexp(y=[1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    ... 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1], x={"age": [36, 24, 51, 44, 29,
    ... 47, 28, 57, 50, 44, 36, 19, 34, 40, 21, 57, 15, 27, 27, 38, 25, 33, 43, 25, 45, 59,
    ... 19, 43, 39, 52, 24, 43, 57, 23, 33, 30, 57, 58, 33, 48], "male": [1, 0, 0, 1, 0, 0,
    ... 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1,
    ... 0, 1, 0, 1, 1, 1]}, ax={"age": [20, 40], "male": [0.5, 0.5]}, pop=[1000, 2000])
    >>> list(r)
    ['coef', 'areas', 'total', 'dev']
    """
    return jcall("synexp", y=y, x=x, ax=ax, pop=pop)


def clrilr(x: Any, type: Any = "clr", ps: Any = 1e-6) -> Any:
    """Log-ratio transforms of compositions.

    The arguments and the result are those of R's ``clrilr()``: see its help
    page.

    Examples
    --------
    >>> r = clrilr(x={"a": [1, 2], "b": [3, 4], "c": [5, 6]}, type="ilr", ps=0)
    >>> len(r)
    2
    """
    return jcall("clrilr", x=x, type=type, ps=ps)


def preslr(x: Any, m1: Any, s1: Any, m2: Any, s2: Any) -> Any:
    """Likelihood ratio of measured evidence.

    The arguments and the result are those of R's ``preslr()``: see its help
    page.

    Examples
    --------
    >>> r = preslr(x=[10, 12, 15], m1=11, s1=2, m2=14, s2=[2, 2, 3])
    >>> list(r)
    ['ll1', 'll2', 'llr', 'lr', 'l10', 'verbal']
    """
    return jcall("preslr", x=x, m1=m1, s1=s1, m2=m2, s2=s2)


def tpdmc(pres: Any = 0, decay: Any = 0.5, mum: Any = 0.5, states: Any = None, start: Any = 1,
          steps: Any = 0, seed: Any = None) -> Any:
    """Decay chain of a body.

    The arguments and the result are those of R's ``tpdmc()``: see its help
    page.

    Examples
    --------
    >>> r = tpdmc(states=["fresh", "bloat"], start="bloat", steps=10, seed=1)
    >>> list(r)
    ['P', 'N', 'absorb', 'time', 'start', 'path']
    """
    return jcall("tpdmc", pres=pres, decay=decay, mum=mum, states=states, start=start, steps=steps,
                 seed=seed)


def tpbhm(y: Any, x: Any, grp: Any = None, pm: Any = 0, psd: Any = 10, level: Any = 0.95) -> Any:
    """Bayesian preservation model.

    The arguments and the result are those of R's ``tpbhm()``: see its help
    page.

    Examples
    --------
    >>> r = tpbhm(y=[9.233, 9.184, 9.858, 9.722, 10.436, 8.813, 11.192, 9.982, 9.752,
    ... 9.637, 11.278, 9.531, 10.071, 9.734, 11.845, 9.161, 9.923, 7.382, 10.888, 9.293,
    ... 11.757, 10.182, 9.733, 10.926, 9.307, 12.682, 10.223, 9.293, 10.417, 10.37],
    ... x={"temp": [9.51, 13, 14.39, 6.71, 26.32, 11.41, 19.17, 13.21, 17.81, 24.86, 8.9,
    ... 9.77, 21.23, 25.87, 10.15, 8.03, 21.03, 29.96, 18.43, 14.29, 10.24, 15.9, 5.54,
    ... 19.37, 10.75, 15.7, 21.23, 7.7, 28.96, 27.32]})
    >>> list(r)
    ['coef', 'sigma', 'tau2', 'groups', 'fitted', 'level']
    """
    return jcall("tpbhm", y=y, x=x, grp=grp, pm=pm, psd=psd, level=level)


def enpath(x: Any, y: Any, fam: Any = "gaussian", alpha: Any = 1, lambda_: Any = None,
           nlam: Any = 100, lmr: Any = None, fold: Any = None, tol: Any = 1e-7, maxit: Any = 1e5,
           scr: Any = True) -> Any:
    """Elastic net by coordinate descent, with cross-validation.

    The arguments and the result are those of R's ``enpath()``: see its help
    page.

    Examples
    --------
    >>> r = enpath(x=[[0.019, -0.596, 1.087], [-0.184, -2.185, -0.763], [-1.371, -0.675,
    ... -0.829], [-0.599, -2.119, 0.834], [0.295, -1.265, -0.968], [0.39, -0.374, -0.029],
    ... [-1.208, -0.688, 0.233], [-0.364, -0.872, -0.301], [-1.627, -0.102, -0.678],
    ... [-0.256, -0.254, 0.655], [1.102, -1.854, -0.401], [0.756, -0.078, -0.335], [-0.238,
    ... 0.969, 1.368], [0.987, 0.185, 2.138], [0.741, -1.38, 0.506], [0.089, -1.436,
    ... 0.786], [-0.955, 0.362, -0.902], [-0.195, -1.759, 0.533], [0.926, -0.325, -0.646],
    ... [0.483, -0.652, 0.291]], y=[0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0,
    ... 1, 1], fam="binomial", alpha=0.5, lambda_=[0.2, 0.1, 0.05], fold=[1, 2, 3, 4, 1, 2,
    ... 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4])
    >>> list(r)
    ['lambda', 'a0', 'beta', 'cvm', 'cvsd', 'lmin', 'l1se']
    """
    return jcall("enpath", x=x, y=y, fam=fam, alpha=alpha, lambda_=lambda_, nlam=nlam, lmr=lmr,
                 fold=fold, tol=tol, maxit=maxit, scr=scr)


def lgtreg(x: Any, y: Any, tol: Any = 1e-10, maxit: Any = 50) -> Any:
    """Logistic regression.

    The arguments and the result are those of R's ``lgtreg()``: see its help
    page.

    Examples
    --------
    >>> r = lgtreg(x=[[0.019, -0.596, 1.087], [-0.184, -2.185, -0.763], [-1.371, -0.675,
    ... -0.829], [-0.599, -2.119, 0.834], [0.295, -1.265, -0.968], [0.39, -0.374, -0.029],
    ... [-1.208, -0.688, 0.233], [-0.364, -0.872, -0.301], [-1.627, -0.102, -0.678],
    ... [-0.256, -0.254, 0.655], [1.102, -1.854, -0.401], [0.756, -0.078, -0.335], [-0.238,
    ... 0.969, 1.368], [0.987, 0.185, 2.138], [0.741, -1.38, 0.506], [0.089, -1.436,
    ... 0.786], [-0.955, 0.362, -0.902], [-0.195, -1.759, 0.533], [0.926, -0.325, -0.646],
    ... [0.483, -0.652, 0.291]], y=[0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1])
    >>> list(r)
    ['coef', 'se', 'dev']
    """
    return jcall("lgtreg", x=x, y=y, tol=tol, maxit=maxit)


def rfor(x: Any, y: Any, newx: Any = None, ntree: Any = 500, mtry: Any = None, minnode: Any = 5,
         frac: Any = 1, repl: Any = True, hon: Any = False, imb: Any = 0, seed: Any = 1) -> Any:
    """Random forest for regression and probabilities.

    The arguments and the result are those of R's ``rfor()``: see its help
    page.

    Examples
    --------
    >>> r = rfor(x=[[0.019, -0.596, 1.087], [-0.184, -2.185, -0.763], [-1.371, -0.675,
    ... -0.829], [-0.599, -2.119, 0.834], [0.295, -1.265, -0.968], [0.39, -0.374, -0.029],
    ... [-1.208, -0.688, 0.233], [-0.364, -0.872, -0.301], [-1.627, -0.102, -0.678],
    ... [-0.256, -0.254, 0.655], [1.102, -1.854, -0.401], [0.756, -0.078, -0.335], [-0.238,
    ... 0.969, 1.368], [0.987, 0.185, 2.138], [0.741, -1.38, 0.506], [0.089, -1.436,
    ... 0.786], [-0.955, 0.362, -0.902], [-0.195, -1.759, 0.533], [0.926, -0.325, -0.646],
    ... [0.483, -0.652, 0.291]], y=[-0.054, 0.772, -1.283, 0.563, 1.247, 0.942, -0.643,
    ... -0.072, -1.407, -0.503, 2.143, 0.366, -1.037, 0.829, 0.984, 1.159, -1.58, 0.555,
    ... 0.773, 1.266], ntree=20, seed=2)
    >>> list(r)
    ['pred', 'oob']
    """
    return jcall("rfor", x=x, y=y, newx=newx, ntree=ntree, mtry=mtry, minnode=minnode, frac=frac,
                 repl=repl, hon=hon, imb=imb, seed=seed)


def crchyp(x: Any, y: Any, home: Any = None) -> Any:
    """Circle hypothesis of a crime series.

    The arguments and the result are those of R's ``crchyp()``: see its help
    page.

    Examples
    --------
    >>> r = crchyp(x=[0, 4, 2, 1], y=[0, 0, 3, 1], home=[2, 1])
    >>> list(r)
    ['centre', 'r', 'pair', 'share', 'mar']
    """
    return jcall("crchyp", x=x, y=y, home=home)


def hitscr(s: Any, home: Any) -> Any:
    """Hit score of a geographic profile.

    The arguments and the result are those of R's ``hitscr()``: see its help
    page.

    Examples
    --------
    >>> r = hitscr(s=[0.1, 0.5, 0.5, 0.9, 0.2], home=[2, 4])
    >>> list(r)
    ['home', 'above', 'ties', 'rank', 'hit', 'hitt']
    """
    return jcall("hitscr", s=s, home=home)


def rtmfit(y: Any, lay: Any, cell: Any = None, rad: Any = None, kind: Any = "prox",
           level: Any = 0.95) -> Any:
    """Risk terrain model.

    The arguments and the result are those of R's ``rtmfit()``: see its help
    page.

    Examples
    --------
    >>> r = rtmfit(y=[1, 0, 2, 0, 0, 5, 0, 1, 4, 0, 1, 2, 4, 3, 3, 2, 2, 1, 1, 2, 1, 3, 1,
    ... 1, 0], lay={"bars": [[2, 2], [4, 4]]}, cell=[[1, 1], [2, 1], [3, 1], [4, 1], [5,
    ... 1], [1, 2], [2, 2], [3, 2], [4, 2], [5, 2], [1, 3], [2, 3], [3, 3], [4, 3], [5, 3],
    ... [1, 4], [2, 4], [3, 4], [4, 4], [5, 4], [1, 5], [2, 5], [3, 5], [4, 5], [5, 5]],
    ... rad=1.5, kind="dens", level=0.9)
    >>> list(r)
    ['coef', 'cells', 'lay', 'dev', 'aic', 'conv', 'level']
    """
    return jcall("rtmfit", y=y, lay=lay, cell=cell, rad=rad, kind=kind, level=level)


def agemix(curves: Any = None, shares: Any = None, ages: Any = None, panel: Any = None,
           nperm: Any = 2000, seed: Any = None) -> Any:
    """The age-crime curve of a mixture of offender types.

    The arguments and the result are those of R's ``agemix()``: see its help
    page.

    Examples
    --------
    >>> r = agemix(panel=[[1, 0], [0, 0], [2, 1], [3, 2], [0, 1], [1, 1]], nperm=99, seed=1)
    >>> list(r)
    ['curve', 'cov', 'panel']
    """
    return jcall("agemix", curves=curves, shares=shares, ages=ages, panel=panel, nperm=nperm,
                 seed=seed)


def cumrsk(age: Any, ev: Any, pop: Any, atrisk: Any = None, dth: Any = None,
           level: Any = 0.95) -> Any:
    """Cumulative risk of a first event by age.

    The arguments and the result are those of R's ``cumrsk()``: see its help
    page.

    Examples
    --------
    >>> r = cumrsk(age=[10, 15, 20, 25], ev=[5, 20, 15, 5], pop=[1000, 1000, 950, 900])
    >>> list(r)
    ['age', 'hazard', 'atrisk', 'cumrisk', 'lo', 'hi']
    """
    return jcall("cumrsk", age=age, ev=ev, pop=pop, atrisk=atrisk, dth=dth, level=level)


def jpipe(num: Any, den: Any, win: Any = 12, lag: Any = 0, grp: Any = None) -> Any:
    """Rates between stages of the justice pipeline.

    The arguments and the result are those of R's ``jpipe()``: see its help
    page.

    Examples
    --------
    >>> r = jpipe(num=[5, 6, 7, 8, 6, 5], den=[10, 12, 11, 15, 12, 10], win=3)
    >>> list(r)
    ['group', 'index', 'num', 'den', 'rate']
    """
    return jcall("jpipe", num=num, den=den, win=win, lag=lag, grp=grp)


def smdbal(x: Any, d: Any, w: Any = None, ps: Any = None, thr: Any = 0.1) -> Any:
    """Covariate balance and overlap.

    The arguments and the result are those of R's ``smdbal()``: see its help
    page.

    Examples
    --------
    >>> r = smdbal(x={"age": [25, 56, 30, 31, 20, 37, 37, 34, 39, 44, 32, 27, 32, 40, 38,
    ... 33, 52, 43, 45, 37, 42, 60, 50, 37, 30, 37, 38, 41, 41, 44], "prior": [3, 2, 5, 2,
    ... 1, 1, 1, 3, 1, 3, 1, 3, 2, 1, 1, 5, 3, 3, 1, 3, 2, 4, 0, 1, 3, 1, 0, 4, 1, 0]},
    ... d=[0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0,
    ... 0, 0, 0])
    >>> list(r)
    ['tab', 'overlap', 'bal', 'ovl', 'thr']
    """
    return jcall("smdbal", x=x, d=d, w=w, ps=ps, thr=thr)


def evlx(est: Any, lo: Any = None, hi: Any = None, type: Any = "rr", rare: Any = False,
         sd: Any = None) -> Any:
    """E-values: how strong an unmeasured confounder must be.

    The arguments and the result are those of R's ``evlx()``: see its help
    page.

    Examples
    --------
    >>> r = evlx(est=1.8, lo=1.2, hi=2.6)
    >>> list(r)
    ['rr', 'evalue', 'evlim']
    """
    return jcall("evlx", est=est, lo=lo, hi=hi, type=type, rare=rare, sd=sd)


def jdgiv(z: Any = None, d: Any = None, y: Any = None, w: Any = None, dsh: Any = None,
          deff: Any = 0) -> Any:
    """The judge-leniency instrument.

    The arguments and the result are those of R's ``jdgiv()``: see its help
    page.

    Examples
    --------
    >>> r = jdgiv(d=[[0, 1], [1, 1], [1, 1], [0, 0], [1, 0]], y={"a": [1, 2, 3, 4, 5],
    ... "b": [2, 2, 4, 4, 5]})
    >>> list(r)
    ['type', 'est', 'shares', 'effects']
    """
    return jcall("jdgiv", z=z, d=d, y=y, w=w, dsh=dsh, deff=deff)


def mmshr(sent: Any, min: Any, grp: Any = None) -> Any:
    """Sentences at, below and above a mandatory minimum.

    The arguments and the result are those of R's ``mmshr()``: see its help
    page.

    Examples
    --------
    >>> r = mmshr(sent=[12, 24, 36, 12, 60], min=24)
    >>> list(r)
    ['group', 'n', 'at', 'below', 'above', 'mexcess']
    """
    return jcall("mmshr", sent=sent, min=min, grp=grp)


def hwkres(times: Any, hor: Any = None, par: Any = None, kernel: Any = "exp", base: Any = "const",
           per: Any = 365.25) -> Any:
    """Hawkes process: time-rescaling residuals.

    The arguments and the result are those of R's ``hwkres()``: see its help
    page.

    Examples
    --------
    >>> r = hwkres(times=[0.5, 1.2, 1.3, 2.8, 3.1, 3.15, 5.6, 7.2, 7.25, 7.3, 9.1], hor=10,
    ... par=[0, 0.4, 2])
    >>> list(r)
    ['tau', 'u', 'd', 'p', 'n']
    """
    return jcall("hwkres", times=times, hor=hor, par=par, kernel=kernel, base=base, per=per)


def hwksim(par: Any, hor: Any, kernel: Any = "exp", base: Any = "const", per: Any = 365.25,
           win: Any = None, sigma: Any = None, seed: Any = None, nmax: Any = 1e5) -> Any:
    """Hawkes process: simulation.

    The arguments and the result are those of R's ``hwksim()``: see its help
    page.

    Examples
    --------
    >>> r = hwksim(par=[0, 0.4, 2], hor=20, seed=1)
    >>> list(r)
    ['t', 'gen']
    """
    return jcall("hwksim", par=par, hor=hor, kernel=kernel, base=base, per=per, win=win,
                 sigma=sigma, seed=seed, nmax=nmax)


def stpcp(y: Any, x: Any = None, minseg: Any = 2, nperm: Any = 9999, seed: Any = None) -> Any:
    """Change point: a step in a series and Pettitt's test.

    The arguments and the result are those of R's ``stpcp()``: see its help
    page.

    Examples
    --------
    >>> r = stpcp(y=[3, 4, 3, 5, 4, 9, 10, 8, 11, 9], nperm=99, seed=1)
    >>> list(r)
    ['tab', 'nperm', 'exact', 'n']
    """
    return jcall("stpcp", y=y, x=x, minseg=minseg, nperm=nperm, seed=seed)


def arwt(x: Any, from_: Any, to: Any, area: Any = None, id: Any = None, top: Any = [10, 20]) -> Any:
    """Crosswalk of counts between two sets of areal units.

    The arguments and the result are those of R's ``arwt()``: see its help
    page.

    Examples
    --------
    >>> r = arwt(x=[5, 6], id=["p", "q"], from_=["p", "q"], to=["Z", "Z"])
    >>> list(r)
    ['out', 'map', 'exact', 'lost', 'conc']
    """
    return jcall("arwt", x=x, from_=from_, to=to, area=area, id=id, top=top)


def bandx(x: Any, cnts: Any = None, rule: Any = "midpoint", cap: Any = None, fac: Any = 2,
          floor: Any = None, clsup: Any = True, intsca: Any = True, drop: Any = False) -> Any:
    """Banded categories: bounds, values and what the open band costs.

    The arguments and the result are those of R's ``bandx()``: see its help
    page.

    Examples
    --------
    >>> r = bandx(x=["0-4", "5-9", "lots", "10+"], cnts=[3, 2, 0, 1], clsup=False)
    >>> list(r)
    ['bands', 'vals', 'sens', 'span', 'relspan']
    """
    return jcall("bandx", x=x, cnts=cnts, rule=rule, cap=cap, fac=fac, floor=floor, clsup=clsup,
                  intsca=intsca, drop=drop)


def bklog(arrive: Any, dispose: Any, horizon: Any = None, target: Any = None) -> Any:
    """Little's law on a court docket.

    The arguments and the result are those of R's ``bklog()``: see its help
    page.

    Examples
    --------
    >>> r = bklog(arrive=[1, 2, 4], dispose=[3, None, 6])
    >>> list(r)
    ['est', 'horizon', 'need']
    """
    return jcall("bklog", arrive=arrive, dispose=dispose, horizon=horizon, target=target)


def bufas(x: Any, y: Any = None, site: Any = None, comp: Any = None, r: Any = None, buf: Any = None,
          lonlat: Any = False) -> Any:
    """Treated, comparison and excluded places around sites.

    The arguments and the result are those of R's ``bufas()``: see its help
    page.

    Examples
    --------
    >>> r = bufas(x=[[2425, 4903], [4554, 790], [288, 4646], [3527, 180], [1567, 247], [2715,
    ... 3351], [3396, 2901], [2534, 2009], [915, 1177], [4361, 1051], [1954, 2184], [794, 4713],
    ... [3252, 3843], [1316, 511], [2380, 1791], [3345, 533], [1231, 2480], [1053, 216], [142,
    ... 3600], [571, 653], [4822, 1923], [4372, 3553], [834, 408], [50, 4226], [1810, 2330],
    ... [783, 1879], [4329, 1228], [2858, 2872], [3572, 1094], [2058, 3240]], site=[[1000,
    ... 1200]], comp=[[2500, 500]], r=900, buf=1000)
    >>> list(r)
    ['pts', 'n', 'par']
    """
    return jcall("bufas", x=x, y=y, site=site, comp=comp, r=r, buf=buf, lonlat=lonlat)


def bvmor(x: Any, y: Any, w: Any, nsim: Any = 999, seed: Any = None) -> Any:
    """Bivariate Moran's I.

    The arguments and the result are those of R's ``bvmor()``: see its help
    page.

    Examples
    --------
    >>> r = bvmor(x=[0.725, 0.857, 0.807, 0.628, 0.566, 0.747, -0.039, 0.741, 0.479, 0.341,
    ... 0.438, 0.295, -0.021, 1.027, 0.082, 0.159, 0.881, 0.793, 0.543, 0.85], y=[0.494, 1.427,
    ... 0.244, 0.632, -0.016, 0.874, 0.97, 1.723, 0.483, 0.632, 0.956, 0.282, -0.013, 0.368,
    ... 0.677, 0.702, 0.982, 0.4, 0.425, 0.628], w={"n": 20, "i": [1, 1, 1, 2, 2, 2, 3, 3, 3, 4,
    ... 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12,
    ... 13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 19, 20,
    ... 20, 20], "j": [4, 6, 9, 6, 8, 17, 5, 14, 15, 3, 6, 15, 3, 15, 16, 1, 4, 8, 16, 18, 20,
    ... 2, 6, 17, 1, 4, 6, 7, 18, 20, 12, 14, 19, 11, 14, 19, 5, 16, 18, 3, 12, 15, 3, 5, 14, 5,
    ... 7, 18, 2, 8, 11, 7, 16, 20, 11, 12, 14, 7, 10, 18], "d": [0.366742416417845,
    ... 0.28800173610587837, 0.19442736432920127, 0.3944641935588071, 0.18802393464662948,
    ... 0.2251421773013666, 0.13849548729110273, 0.11827087553578014, 0.05575840743780257,
    ... 0.27614670014323905, 0.11901680553602503, 0.2229125389025929, 0.13849548729110273,
    ... 0.18106904760339354, 0.11417968295629481, 0.28800173610587837, 0.11901680553602503,
    ... 0.23700843866833102, 0.18056023925549064, 0.13892803892663277, 0.0458802789878178,
    ... 0.18802393464662948, 0.23700843866833102, 0.2545270123189285, 0.19442736432920127,
    ... 0.329802971484491, 0.31624199594614244, 0.24263758983306769, 0.261663906567184,
    ... 0.2241428116179505, 0.22227235545609356, 0.23221541723150083, 0.049648766349225645,
    ... 0.22227235545609356, 0.1661715980545412, 0.18404890654388575, 0.2651886121235224,
    ... 0.25132051249350895, 0.3528639964632266, 0.11827087553578014, 0.1661715980545412,
    ... 0.15305227865013962, 0.05575840743780257, 0.18106904760339354, 0.15305227865013962,
    ... 0.11417968295629481, 0.18056023925549064, 0.1494824404403407, 0.2251421773013666,
    ... 0.2545270123189285, 0.26475082624989105, 0.13892803892663277, 0.1494824404403407,
    ... 0.17768511473953016, 0.049648766349225645, 0.18404890654388575, 0.18309833423600555,
    ... 0.0458802789878178, 0.2241428116179505, 0.17768511473953016], "w": [0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333]}, nsim=19, seed=1)
    >>> list(r)
    ['i', 'mean', 'sd', 'z', 'p', 'ptwo']
    """
    return jcall("bvmor", x=x, y=y, w=w, nsim=nsim, seed=seed)


def calaud(area: Any, risk: Any = None, outcome: Any = None, group: Any = None,
           population: Any = None, nperm: Any = 999, seed: Any = None) -> Any:
    """Calibration audit of a predictive-policing ranking.

    The arguments and the result are those of R's ``calaud()``: see its help
    page.

    Examples
    --------
    >>> r = calaud(area={"area": ["a1", "a1", "a1", "a1", "a1", "a1", "a1", "a1", "a2", "a2",
    ... "a2", "a2", "a2", "a2", "a2", "a2", "a3", "a3", "a3", "a3", "a3", "a3", "a3", "a3",
    ... "a4", "a4", "a4", "a4", "a4", "a4", "a4", "a4", "a5", "a5", "a5", "a5", "a5", "a5",
    ... "a5", "a5", "a6", "a6", "a6", "a6", "a6", "a6", "a6", "a6"], "risk": [0.522, 0.95,
    ... 0.427, 0.424, 0.945, 0.452, 0.833, 0.173, 0.934, 0.243, 0.328, 0.83, 0.652, 0.014,
    ... 0.821, 0.115, 0.101, 0.474, 0.232, 0.108, 0.143, 0.567, 0.368, 0.479, 0.568, 0.049,
    ... 0.148, 0.815, 0.161, 0.07, 0.179, 0.667, 0.481, 0.682, 0.804, 0.783, 0.513, 0.217,
    ... 0.178, 0.251, 0.05, 0.738, 0.42, 0.684, 0.732, 0.503, 0.776, 0.017], "group": ["g", "h",
    ... "g", "h", "g", "h", "g", "h", "g", "h", "g", "h", "g", "h", "g", "h", "g", "h", "g",
    ... "h", "g", "h", "g", "h", "g", "h", "g", "h", "g", "h", "g", "h", "g", "h", "g", "h",
    ... "g", "h", "g", "h", "g", "h", "g", "h", "g", "h", "g", "h"], "outcome": [1, 1, 0, 0, 0,
    ... 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1,
    ... 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0]}, nperm=0)
    >>> list(r)
    ['areas', 'groups', 'spear', 'ties', 'between', 'nperm']
    """
    return jcall("calaud", area=area, risk=risk, outcome=outcome, group=group,
                  population=population, nperm=nperm, seed=seed)


def cbatt(y: Any, d: Any, x: Any, est: Any = "ipw", tgt: Any = "ate", ps: Any = None, k: Any = 1,
          calip: Any = None, nsub: Any = 5, reps: Any = 499, seed: Any = None,
          level: Any = 0.95) -> Any:
    """Design-based causal effects of a 0/1 treatment.

    The arguments and the result are those of R's ``cbatt()``: see its help
    page.

    Examples
    --------
    >>> r = cbatt(y=[-0.242, -0.834, -3.562, 0.228, 3.251, 0.281, -0.55, -2.047, -0.356, 2.096,
    ... 1.928, -1.491, 1.843, -2.177, -0.69, -0.975, -0.08, 2.173, 0.214, -2.681, -1.929, 1.241,
    ... 3.57, 0.486, 2.846, 4.007, 0.227, 0.641, -3.443, -1.683, -0.276, 2.735, 3.284, 1.853,
    ... -0.8, -0.307, -0.635, 2.525, -0.112, -0.153], d=[0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1,
    ... 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0],
    ... x={"a": [-0.186, -1.229, -0.427, -0.596, 0.467, 0.422, -1.022, -0.622, 0.837, 0.962,
    ... 0.104, -0.06, 0.704, -0.755, -0.546, -0.455, -0.117, 0.548, 0.08, -0.232, -1.472, 0.177,
    ... 1.21, -0.34, -0.249, 1.152, -0.299, -0.557, -1.728, -0.765, -0.557, 2.135, 2.281, 1.197,
    ... -1.374, -1.065, -1.613, -1.098, -1.715, -1.4], "b": [-0.695, -0.455, 1.049, -0.332,
    ... -0.734, 0.823, 0.039, 0.819, 0.912, 1.323, 0.784, 0.365, -0.909, 1.323, -0.573, 0.939,
    ... 0.237, 1.032, -0.099, 1.885, 0.771, 0.616, -0.7, 1.252, 0.3, -0.024, 0.924, -1.056,
    ... 1.078, -0.335, 1.201, 1.273, 1.174, 0.446, -1.306, -0.386, -0.224, -1.827, -1.17,
    ... -0.824]}, est="ebal", tgt="att", reps=0)
    >>> list(r)
    ['res', 'ps', 'w', 'match', 'est', 'tgt', 'level']
    """
    return jcall("cbatt", y=y, d=d, x=x, est=est, tgt=tgt, ps=ps, k=k, calip=calip, nsub=nsub,
                  reps=reps, seed=seed, level=level)


def cntreg(y: Any, grp: Any = None, adj: Any = None, expo: Any = None, clus: Any = None,
           ref: Any = None, gee: Any = False, cadj: Any = True, level: Any = 0.95) -> Any:
    """Count regression: rate ratios between groups.

    The arguments and the result are those of R's ``cntreg()``: see its help
    page.

    Examples
    --------
    >>> r = cntreg(y=[5, 12], grp=["a", "b"], expo=[1000, 1500])
    >>> list(r)
    ['irr', 'coef', 'fit', 'ref', 'n', 'ncl', 'level', 'warnings']
    """
    return jcall("cntreg", y=y, grp=grp, adj=adj, expo=expo, clus=clus, ref=ref, gee=gee, cadj=cadj,
                  level=level)


def crossk(x: Any, y: Any = None, marks: Any = None, i: Any = None, j: Any = None, win: Any = None,
           r: Any = None, sigma: Any = None, nsim: Any = 99, seed: Any = None) -> Any:
    """Cross-type K function and the segregation test.

    The arguments and the result are those of R's ``crossk()``: see its help
    page.

    Examples
    --------
    >>> r = crossk(x=[[0.061, 0.39], [0.182, 0.07], [0.484, 0.708], [0.338, 0.693], [0.813,
    ... 0.096], [0.87, 0.037], [0.019, 0.401], [0.825, 0.012], [0.401, 0.568], [0.672, 0.945],
    ... [0.119, 0.343], [0.476, 0.244], [0.925, 0.542], [0.575, 0.118], [0.913, 0.375], [0.755,
    ... 0.328], [0.631, 0.947], [0.794, 0.911], [0.924, 0.762], [0.164, 0.758], [0.193, 0.998],
    ... [0.098, 0.58], [0.829, 0.213], [0.806, 0.134], [0.279, 0.792], [0.449, 0.982], [0.626,
    ... 0.982], [0.371, 0.125], [0.382, 0.618], [0.914, 0.858]], marks=["burglary", "burglary",
    ... "burglary", "burglary", "robbery", "robbery", "burglary", "robbery", "burglary",
    ... "robbery", "burglary", "burglary", "robbery", "robbery", "robbery", "robbery",
    ... "robbery", "robbery", "robbery", "burglary", "burglary", "burglary", "robbery",
    ... "robbery", "burglary", "burglary", "robbery", "burglary", "burglary", "robbery"],
    ... i="burglary", j="robbery", win=[0, 1, 0, 1], r=[0.05, 0.1, 0.2])
    >>> list(r)
    ['k', 'seg', 'i', 'j']
    """
    return jcall("crossk", x=x, y=y, marks=marks, i=i, j=j, win=win, r=r, sigma=sigma, nsim=nsim,
                  seed=seed)


def ctabx(row: Any, col: Any = None, w: Any = None, agg: Any = None, correct: Any = False,
          exact: Any = False, level: Any = 0.95) -> Any:
    """Cross-tabulation with tests and measures of association.

    The arguments and the result are those of R's ``ctabx()``: see its help
    page.

    Examples
    --------
    >>> r = ctabx(row=["N", "N", "E", "W", "E"], col=["N", "E", "E", "W", "E"])
    >>> list(r)
    ['tab', 'expd', 'res', 'stat', 'drop', 'agg', 'correct', 'level']
    """
    return jcall("ctabx", row=row, col=col, w=w, agg=agg, correct=correct, exact=exact, level=level)


def cusrv(y: Any, expd: Any = None, start: Any = None, k: Any = 1.04, h: Any = 2.26,
          trans: Any = "std", reset: Any = False) -> Any:
    """CUSUM surveillance of counts.

    The arguments and the result are those of R's ``cusrv()``: see its help
    page.

    Examples
    --------
    >>> r = cusrv(y=[10, 12, 9, 11, 10, 13, 19, 22, 25, 12], start=6)
    >>> list(r)
    ['i', 'y', 'expd', 'z', 's', 'alarm']
    """
    return jcall("cusrv", y=y, expd=expd, start=start, k=k, h=h, trans=trans, reset=reset)


def dastd(x: Any, n: Any, strat: Any, grp: Any = None, std: Any = None, ref: Any = None,
          per: Any = 1, level: Any = 0.95) -> Any:
    """Directly standardised rates.

    The arguments and the result are those of R's ``dastd()``: see its help
    page.

    Examples
    --------
    >>> r = dastd(x=[30, 22, 5, 41, 30, 12], n=[5000, 9000, 6000, 3000, 8000, 7000],
    ... strat=["18-29", "30-49", "50+", "18-29", "30-49", "50+"], grp=["A", "A", "A", "B", "B",
    ... "B"], std=[30, 45, 25], per=10000)
    >>> list(r)
    ['rates', 'cmp', 'ref', 'per', 'level']
    """
    return jcall("dastd", x=x, n=n, strat=strat, grp=grp, std=std, ref=ref, per=per, level=level)


def dbscn(x: Any, y: Any = None, eps: Any = None, minpts: Any = 5, border: Any = True,
          lonlat: Any = False) -> Any:
    """DBSCAN clusters of incidents.

    The arguments and the result are those of R's ``dbscn()``: see its help
    page.

    Examples
    --------
    >>> r = dbscn(x=[[-0.299, 0.114], [-0.294, -0.099], [0.025, 0], [-0.199, 0.225], [0, 0.288],
    ... [-0.086, -0.219], [-0.123, -0.023], [-0.405, 0.24], [-0.245, -0.094], [0.036, -0.01],
    ... [2.983, 0.992], [2.822, 0.69], [2.911, 1.233], [2.994, 0.945], [2.917, 0.906], [3.223,
    ... 0.752], [2.904, 0.998], [2.913, 0.84], [3.139, 0.893], [2.789, 1.258], [-2, 3], [5,
    ... -2]], eps=0.5, minpts=3, border=False)
    >>> list(r)
    ['cl', 'core', 'sizes']
    """
    return jcall("dbscn", x=x, y=y, eps=eps, minpts=minpts, border=border, lonlat=lonlat)


def didfe(y: Any, unit: Any, time: Any, trt: Any, adj: Any = None, clus: Any = None, w: Any = None,
          fe: Any = "both", fam: Any = "ols", expo: Any = None, loo: Any = False,
          bacon: Any = False, mtch: Any = 0, exact: Any = None, calip: Any = None,
          level: Any = 0.95) -> Any:
    """Difference in differences with two-way fixed effects.

    The arguments and the result are those of R's ``didfe()``: see its help
    page.

    Examples
    --------
    >>> r = didfe(y=[32, 32, 32, 28, 28, 29, 33, 33, 33, 29, 29, 30, 33, 34, 35, 30, 30, 31, 34,
    ... 35, 36, 35, 31, 32, 35, 36, 37, 37, 32, 33, 36, 37, 38, 38, 38, 38, 36, 38, 39, 39, 39,
    ... 39, 37, 39, 40, 40, 40, 40], unit=["A", "A", "A", "A", "A", "A", "B", "B", "B", "B",
    ... "B", "B", "C", "C", "C", "C", "C", "C", "D", "D", "D", "D", "D", "D", "E", "E", "E",
    ... "E", "E", "E", "F", "F", "F", "F", "F", "F", "G", "G", "G", "G", "G", "G", "H", "H",
    ... "H", "H", "H", "H"], time=[2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018,
    ... 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020,
    ... 2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016,
    ... 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020], trt=[0, 0, 0, 1, 1, 1, 0,
    ... 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
    ... 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], fam="pois", clus=["p", "q", "r", "s", "p", "q",
    ... "r", "s", "p", "q", "r", "s", "p", "q", "r", "s", "p", "q", "r", "s", "p", "q", "r",
    ... "s", "p", "q", "r", "s", "p", "q", "r", "s", "p", "q", "r", "s", "p", "q", "r", "s",
    ... "p", "q", "r", "s", "p", "q", "r", "s"])
    >>> list(r)
    ['coef', 'fit', 'vcov', 'loo', 'bacon', 'match', 'fam', 'level', 'warnings']
    """
    return jcall("didfe", y=y, unit=unit, time=time, trt=trt, adj=adj, clus=clus, w=w, fe=fe,
                  fam=fam, expo=expo, loo=loo, bacon=bacon, mtch=mtch, exact=exact, calip=calip,
                  level=level)


def dissix(x: Any, delta: Any = 0.5, w: Any = None) -> Any:
    """Indices of residential segregation.

    The arguments and the result are those of R's ``dissix()``: see its help
    page.

    Examples
    --------
    >>> r = dissix(x=[[90, 10], [80, 20], [10, 85], [5, 95], [50, 50]])
    >>> list(r)
    ['d', 'grp', 'dm', 'delta']
    """
    return jcall("dissix", x=x, delta=delta, w=w)


def dktflw(recv: Any = None, disp: Any = None, pend: Any = None, prd: Any = None, win: Any = 12,
           level: Any = 0.95, span: Any = None, mtime: Any = None, shares: Any = None,
           floors: Any = None, mix: Any = None, scores: Any = None) -> Any:
    """Case flow through a court.

    The arguments and the result are those of R's ``dktflw()``: see its help
    page.

    Examples
    --------
    >>> r = dktflw(mix=[[50, 30, 20], [55, 28, 17], [60, 25, 15]])
    >>> list(r)
    ['flow', 'mix', 'level']
    """
    return jcall("dktflw", recv=recv, disp=disp, pend=pend, prd=prd, win=win, level=level,
                  span=span, mtime=mtime, shares=shares, floors=floors, mix=mix, scores=scores)


def dmlx(y: Any, d: Any, x: Any, model: Any = "plr", fold: Any = None, z: Any = None, s: Any = None,
         t: Any = None, score: Any = None, lrn: Any = "lasso", lrnd: Any = None, alpha: Any = 0.5,
         nlam: Any = 100, nfin: Any = 5, tol: Any = 1e-07, scr: Any = True, hon: Any = False,
         ntree: Any = 500, mtry: Any = None, trim: Any = None, clu: Any = None, pred: Any = None,
         stk: Any = ["glm", "lasso", "forest"], nstk: Any = 5, disc: Any = False,
         seed: Any = 1) -> Any:
    """Double machine learning by cross-fitting.

    The arguments and the result are those of R's ``dmlx()``: see its help
    page.

    Examples
    --------
    >>> r = dmlx(y=[-1.418, -0.367, -0.158, 1.992, 3.463, -0.068, 0.348, -4.519, -2.436, 2.461,
    ... 0.649, -1.032, 1.474, -2.791, 0.979, 2.049, 0.208, 1.547, 2.809, 3.293, -5.781, 3.076,
    ... -1.021, 0.83, -0.302, 1.775, 0.01, -0.769, -1.261, -2.589, -2.647, -2.836, 1.756,
    ... -0.766, -2.634, -0.94, -1.4, 4.201, -1.312, -2.64, 0.291, -0.592, 2.629, 1.583, 4.009,
    ... -1.028, -0.528, 0.363, -3.881, 0.847, -1.43, 1.396, 0.975, 0.232, 4.456, -1.902, -1.319,
    ... 3.525, -1.358, -2.401], d=[-1.09, 0.273, 0.831, 2.059, 2.354, -1.357, 1.367, -0.19,
    ... -0.576, 2.23, 2.117, -1.109, 1.728, -1.086, 0.028, 1.863, 0.458, 1.798, -0.598, 2.72,
    ... -1.65, 1.573, -0.466, -0.684, 0.143, 3.925, -0.059, -0.204, -1.4, -1.057, -1.642,
    ... -0.748, -0.432, -1.143, -1.958, -0.99, -0.752, 2.211, 0.391, -0.276, 1.574, -0.199,
    ... 1.196, 2.24, 2.68, -1.18, -0.663, -0.956, -1.668, 0.266, 0.277, 0.931, 0.391, -1.286,
    ... 2.459, -1.746, -0.424, 2.493, 0.455, -0.849], x=[[-0.794, 0.882, 0.105], [0.197, 0.899,
    ... -1.049], [1.002, 1.503, 0.889], [1.289, 0.505, -0.556], [0.906, -0.341, 0.345], [0.494,
    ... 1.373, 0.557], [0.599, 0.819, 0.139], [-1.58, 1.426, -0.097], [1.001, 2.196, 0.568],
    ... [2.188, -0.896, 0.052], [-1.209, -1.69, -1.457], [-0.587, -0.467, 0.873], [1.056, -0.36,
    ... -0.178], [-0.317, 0.746, -0.345], [-0.055, 0.169, 0.383], [0.33, 1.65, -1.602], [0.663,
    ... 0.033, 1.214], [0.878, -1.032, 0.588], [0.203, -0.632, -0.594], [2.274, 0.578, -0.943],
    ... [-0.899, 1.842, 1.262], [2.139, -0.058, 0.267], [-1.166, -0.257, 1.31], [-0.047, -2.167,
    ... -0.862], [0.416, 0.11, 0.714], [1.72, 1.315, 0.093], [-0.784, -0.893, -0.831], [-1.304,
    ... -1.044, 0.781], [-0.452, -0.542, -0.018], [-1.774, -0.331, -0.373], [-0.05, 1.672,
    ... 0.695], [-1.037, -1.053, -1.264], [0.622, 0.081, 0.378], [0.127, 0.703, -0.534],
    ... [-1.741, 1.511, -0.812], [-0.573, 0.001, -1.479], [0.912, 1.431, -0.251], [1.234,
    ... -1.112, 0.203], [0.661, 0.568, -0.176], [-0.861, 1.338, -0.812], [-0.533, 0.073,
    ... -0.117], [0.321, 1.196, 1.274], [-0.004, -2.112, -1.196], [1.598, 0.671, -0.871],
    ... [1.327, -0.186, 0.434], [-0.119, -0.192, -0.426], [0.082, 0.787, 0.322], [-0.02, 0.562,
    ... 0.227], [-0.825, -0.019, -1.156], [0.702, -0.26, -0.321], [-0.196, -0.162, -1.072],
    ... [-0.158, -2.09, -0.997], [0.224, 0.449, 0.864], [0.305, -0.697, 0.241], [1.369, 0.782,
    ... -0.591], [-1.389, -0.265, -1.827], [-1.644, 0.02, -0.339], [1.426, 0.401, 1.611],
    ... [0.097, 1.967, -0.664], [0.312, 0.954, -0.565]], fold=[2, 3, 2, 3, 1, 1, 2, 2, 3, 2, 3,
    ... 2, 2, 2, 1, 3, 3, 3, 3, 1, 1, 3, 1, 3, 1, 1, 1, 2, 2, 1, 1, 3, 3, 1, 2, 1, 2, 1, 1, 2,
    ... 3, 3, 2, 1, 2, 1, 2, 2, 3, 1, 3, 2, 1, 3, 3, 1, 2, 3, 3, 2], lrn="glm", lrnd="glm")
    >>> list(r)
    ['coef', 'se', 'reps', 'psia', 'psib', 'pred']
    """
    return jcall("dmlx", y=y, d=d, x=x, model=model, fold=fold, z=z, s=s, t=t, score=score, lrn=lrn,
                  lrnd=lrnd, alpha=alpha, nlam=nlam, nfin=nfin, tol=tol, scr=scr, hon=hon,
                  ntree=ntree, mtry=mtry, trim=trim, clu=clu, pred=pred, stk=stk, nstk=nstk,
                  disc=disc, seed=seed)


def ebshr(y: Any, v: Any = None, w: Any = None, b: Any = None) -> Any:
    """Empirical-Bayes shrinkage of hot-spot counts.

    The arguments and the result are those of R's ``ebshr()``: see its help
    page.

    Examples
    --------
    >>> r = ebshr(y=[40, 12, 9, 25, 7, 31, 5, 18])
    >>> list(r)
    ['fit', 'est']
    """
    return jcall("ebshr", y=y, v=v, w=w, b=b)


def etasf(times: Any, hor: Any, mag: Any = None, m0: Any = None, start: Any = None,
          maxit: Any = 500) -> Any:
    """ETAS: a self-exciting process with marks and the Omori law.

    The arguments and the result are those of R's ``etasf()``: see its help
    page.

    Examples
    --------
    >>> r = etasf(times=[0.222, 6.498, 7.884, 14.087, 14.727, 20.374, 35.21, 35.283, 36.968,
    ... 38.448, 41.998, 42.381, 49.685, 51.543, 54.964, 58.55, 62.19, 64.88, 71.879, 72.566,
    ... 80.022, 80.399, 81.498, 83.703, 83.947, 84.585, 85.471, 88.323, 91.356, 92.289],
    ... hor=100)
    >>> list(r)
    ['par', 'se', 'cov', 'loglik', 'aic', 'branch', 'conv', 'iter', 'afree', 'n', 'hor']
    """
    return jcall("etasf", times=times, hor=hor, mag=mag, m0=m0, start=start, maxit=maxit)


def evtst(y: Any, unit: Any, time: Any, first: Any, leads: Any = 4, lags: Any = 4, ref: Any = -1,
          bin: Any = True, adj: Any = None, clus: Any = None, w: Any = None, ma: Any = 1,
          norm: Any = True, mbar: Any = None, target: Any = None, draws: Any = 5000,
          seed: Any = None, level: Any = 0.95) -> Any:
    """An event study: the effect period by period around adoption.

    The arguments and the result are those of R's ``evtst()``: see its help
    page.

    Examples
    --------
    >>> r = evtst(y=[31.8415, 32.4093, 32.1411, 27.7432, 28.0411, 29.2206, 32.657, 33.4894,
    ... 33.4121, 28.956, 29, 29.9634, 33.4202, 34.4906, 34.6503, 30.2121, 30.0386, 30.749,
    ... 34.1499, 35.4129, 35.8367, 35.4911, 31.1538, 31.5944, 34.8676, 36.2626, 36.9564,
    ... 36.7709, 32.3364, 32.512, 35.596, 37.0514, 37.9999, 38.0291, 37.5718, 37.5082, 36.3565,
    ... 37.7964, 38.9638, 39.2451, 38.8414, 38.5835, 37.1682, 38.5177, 39.8509, 40.4018,
    ... 40.1236, 39.7317], unit=["A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B",
    ... "C", "C", "C", "C", "C", "C", "D", "D", "D", "D", "D", "D", "E", "E", "E", "E", "E",
    ... "E", "F", "F", "F", "F", "F", "F", "G", "G", "G", "G", "G", "G", "H", "H", "H", "H",
    ... "H", "H"], time=[2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020,
    ... 2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016,
    ... 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018,
    ... 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020], first=[2018, 2018, 2018, 2018, 2018,
    ... 2018, 2018, 2018, 2018, 2018, 2018, 2018, 2018, 2018, 2018, 2018, 2018, 2018, 2019,
    ... 2019, 2019, 2019, 2019, 2019, 2019, 2019, 2019, 2019, 2019, 2019, None, None, None,
    ... None, None, None, None, None, None, None, None, None, None, None, None, None, None,
    ... None], leads=2, lags=2)
    >>> list(r)
    ['est', 'pre', 'vcov', 'coh', 'sens', 'brk', 'n', 'ncl', 'level']
    """
    return jcall("evtst", y=y, unit=unit, time=time, first=first, leads=leads, lags=lags, ref=ref,
                  bin=bin, adj=adj, clus=clus, w=w, ma=ma, norm=norm, mbar=mbar, target=target,
                  draws=draws, seed=seed, level=level)


def gmoran(x: Any, w: Any, nsim: Any = 0, seed: Any = None) -> Any:
    """Global Moran's I with its inference.

    The arguments and the result are those of R's ``gmoran()``: see its help
    page.

    Examples
    --------
    >>> r = gmoran(x={"y1": [3, 3, 5, 2, 3, 4, 1, 2, 1, 1, 2, 2, 2, 2, 3, 5, 4, 6, 2, 5], "y2":
    ... [4, 4, 6, 6, 5, 3, 7, 10, 2, 6, 5, 8, 6, 4, 7, 4, 4, 5, 7, 6]}, w={"n": 20, "i": [1, 1,
    ... 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10,
    ... 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 16, 17, 17, 17, 18,
    ... 18, 18, 19, 19, 19, 20, 20, 20], "j": [4, 6, 9, 6, 8, 17, 5, 14, 15, 3, 6, 15, 3, 15,
    ... 16, 1, 4, 8, 16, 18, 20, 2, 6, 17, 1, 4, 6, 7, 18, 20, 12, 14, 19, 11, 14, 19, 5, 16,
    ... 18, 3, 12, 15, 3, 5, 14, 5, 7, 18, 2, 8, 11, 7, 16, 20, 11, 12, 14, 7, 10, 18], "d":
    ... [0.366742416417845, 0.28800173610587837, 0.19442736432920127, 0.3944641935588071,
    ... 0.18802393464662948, 0.2251421773013666, 0.13849548729110273, 0.11827087553578014,
    ... 0.05575840743780257, 0.27614670014323905, 0.11901680553602503, 0.2229125389025929,
    ... 0.13849548729110273, 0.18106904760339354, 0.11417968295629481, 0.28800173610587837,
    ... 0.11901680553602503, 0.23700843866833102, 0.18056023925549064, 0.13892803892663277,
    ... 0.0458802789878178, 0.18802393464662948, 0.23700843866833102, 0.2545270123189285,
    ... 0.19442736432920127, 0.329802971484491, 0.31624199594614244, 0.24263758983306769,
    ... 0.261663906567184, 0.2241428116179505, 0.22227235545609356, 0.23221541723150083,
    ... 0.049648766349225645, 0.22227235545609356, 0.1661715980545412, 0.18404890654388575,
    ... 0.2651886121235224, 0.25132051249350895, 0.3528639964632266, 0.11827087553578014,
    ... 0.1661715980545412, 0.15305227865013962, 0.05575840743780257, 0.18106904760339354,
    ... 0.15305227865013962, 0.11417968295629481, 0.18056023925549064, 0.1494824404403407,
    ... 0.2251421773013666, 0.2545270123189285, 0.26475082624989105, 0.13892803892663277,
    ... 0.1494824404403407, 0.17768511473953016, 0.049648766349225645, 0.18404890654388575,
    ... 0.18309833423600555, 0.0458802789878178, 0.2241428116179505, 0.17768511473953016], "w":
    ... [0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333]})
    >>> list(r)
    ['var', 'i', 'e', 'vn', 'vr', 'zn', 'zr', 'pn', 'pr', 'pperm', 's0']
    """
    return jcall("gmoran", x=x, w=w, nsim=nsim, seed=seed)


def grprr(x: Any, n: Any = None, grp: Any = None, ref: Any = None, type: Any = "rate", per: Any = 1,
          level: Any = 0.95) -> Any:
    """Rates and proportions by group: ratios and differences.

    The arguments and the result are those of R's ``grprr()``: see its help
    page.

    Examples
    --------
    >>> r = grprr(x=[1, 0, 0, 1, 1, 0, 1, 1], grp=[1, 1, 1, 1, 2, 2, 2, 2], type="prop", ref=2,
    ... level=0.9)
    >>> list(r)
    ['rates', 'cmp', 'ref', 'type', 'per', 'level']
    """
    return jcall("grprr", x=x, n=n, grp=grp, ref=ref, type=type, per=per, level=level)


def hsper(x: Any, y: Any = None, prd: Any = None, cell: Any = None, top: Any = 0.05,
          org: Any = None, lonlat: Any = False) -> Any:
    """Hot spots and how long they last.

    The arguments and the result are those of R's ``hsper()``: see its help
    page.

    Examples
    --------
    >>> r = hsper(x=[[1.947, 2.061], [1.678, 1.896], [2.049, 2.076], [1.891, 1.612], [2.177,
    ... 1.712], [2.43, 2.326], [1.702, 2.121], [2.136, 2.176], [2.025, 2.545], [2.269, 2.039],
    ... [1.931, 1.4], [2.251, 2.1], [1.476, 2.351], [2.507, 2.618], [2.259, 1.587], [1.955,
    ... 1.655], [1.565, 1.788], [2.193, 1.684], [2.145, 1.806], [1.998, 1.944], [5.602, 1.148],
    ... [6.525, 4.828], [2.796, 9.792], [9.799, 8.115], [6.439, 5.429], [5.826, 0.724], [6.159,
    ... 4.665], [9.251, 3.399], [3.9, 6.899], [2.879, 5.142], [0.907, 5.149], [3.22, 5.455],
    ... [7.583, 4.475], [1.044, 0.839], [7.103, 9.301], [9.665, 0.164], [2.015, 4.141], [1.085,
    ... 2.27], [0.554, 0.996], [8.297, 4.829]], prd=[1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
    ... 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2], cell=2,
    ... org=[-1, -1])
    >>> list(r)
    ['periods', 'cells', 'hot', 'cell', 'k']
    """
    return jcall("hsper", x=x, y=y, prd=prd, cell=cell, top=top, org=org, lonlat=lonlat)


def hwkfit(times: Any, hor: Any, kernel: Any = "exp", base: Any = "const", per: Any = 365.25,
           start: Any = None, res: Any = None, seed: Any = None, maxit: Any = 500) -> Any:
    """Hawkes process: self-exciting event times.

    The arguments and the result are those of R's ``hwkfit()``: see its help
    page.

    Examples
    --------
    >>> r = hwkfit(times=[1.99, 5.17, 7.57, 8.92, 26.03, 28.44, 36.05, 40.55, 40.85, 41.53,
    ... 41.55, 41.65, 42.07, 42.55, 42.89, 43.38, 43.67, 44.08, 52.37, 61.04, 62.24, 67.55,
    ... 75.98, 80.23, 95.54], hor=100)
    >>> list(r)
    ['par', 'se', 'cov', 'loglik', 'aic', 'bic', 'branch', 'conv', 'iter', 'kernel', 'base',
    'per', 'n', 'hor', 'jit', 'times', 'cmp']
    """
    return jcall("hwkfit", times=times, hor=hor, kernel=kernel, base=base, per=per, start=start,
                  res=res, seed=seed, maxit=maxit)


def hypq(lam: Any, pref: Any, mu: Any = 1, queue: Any = "loss", tt: Any = None) -> Any:
    """Hypercube queue of patrol units.

    The arguments and the result are those of R's ``hypq()``: see its help
    page.

    Examples
    --------
    >>> r = hypq(lam=[0.2, 0.3], pref=[[1, 2], [2, 1]], queue="wait")
    >>> list(r)
    ['state', 'wl', 'disp', 'queue', 'pfull', 'lq', 'wq', 'rho', 'travel']
    """
    return jcall("hypq", lam=lam, pref=pref, mu=mu, queue=queue, tt=tt)


def jmono(judge: Any, d: Any, y: Any = None, grp: Any = None, w: Any = None, lo: Any = None,
          hi: Any = None, minn: Any = 30, adj: Any = "holm", alpha: Any = 0.05) -> Any:
    """Do judges' leniency rankings hold? Monotonicity checks.

    The arguments and the result are those of R's ``jmono()``: see its help
    page.

    Examples
    --------
    >>> r = jmono(judge=["J2", "J4", "J2", "J4", "J1", "J2", "J4", "J4", "J1", "J1", "J1", "J4",
    ... "J1", "J3", "J2", "J4", "J2", "J4", "J2", "J2", "J4", "J1", "J3", "J2", "J4", "J2",
    ... "J1", "J3", "J3", "J1", "J2", "J4", "J1", "J1", "J3", "J2", "J3", "J2", "J3", "J4",
    ... "J1", "J1", "J2", "J3", "J3", "J3", "J2", "J3", "J1", "J3", "J2", "J4", "J1", "J2",
    ... "J2", "J1", "J1", "J3", "J2", "J3", "J3", "J1", "J2", "J1", "J4", "J3", "J1", "J4",
    ... "J2", "J3", "J1", "J1", "J3", "J2", "J1", "J2", "J4", "J4", "J1", "J3"], d=[0, 0, 0, 1,
    ... 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0,
    ... 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1,
    ... 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0])
    >>> list(r)
    ['judges', 'homog', 'slope', 'sub', 'alpha']
    """
    return jcall("jmono", judge=judge, d=d, y=y, grp=grp, w=w, lo=lo, hi=hi, minn=minn, adj=adj,
                  alpha=alpha)


def kdsurf(x: Any, y: Any = None, w: Any = None, h: Any = None, win: Any = None, nx: Any = 100,
           ny: Any = None, edge: Any = "none", at: Any = None, cover: Any = [0.5, 0.9]) -> Any:
    """Kernel intensity surface of incidents.

    The arguments and the result are those of R's ``kdsurf()``: see its help
    page.

    Examples
    --------
    >>> r = kdsurf(x=[[2.756, 3.6], [2.129, 3.522], [2.044, 2.498], [1.94, 3.924], [1.403,
    ... 2.667], [2.306, 3.053], [1.891, 2.789], [1.909, 2.939], [2.467, 3.094], [2.411, 3.06],
    ... [2.696, 2.987], [1.762, 3.054], [2.325, 2.757], [2.696, 2.748], [1.445, 2.169], [1.168,
    ... 2.107], [4.705, 0.954], [0.773, 1.825], [0.775, 0.105], [0.434, 5.979], [0.319, 4.826],
    ... [3.191, 0.519], [0.674, 5.22], [4.459, 3.328], [4.388, 2.528]], h=0.4, at="points")
    >>> list(r)
    ['x', 'y', 'z', 'h', 'regions']
    """
    return jcall("kdsurf", x=x, y=y, w=w, h=h, win=win, nx=nx, ny=ny, edge=edge, at=at, cover=cover)


def kmrec(time: Any, event: Any = None, grp: Any = None, unit: Any = None, end: Any = None,
          level: Any = 0.95) -> Any:
    """Kaplan-Meier survival and recurrence.

    The arguments and the result are those of R's ``kmrec()``: see its help
    page.

    Examples
    --------
    >>> r = kmrec(time=[3, 8, 8, 12, 20, 24], event=[1, 1, 0, 1, 0, 1], grp=["a", "a", "a", "b",
    ... "b", "b"])
    >>> list(r)
    ['km', 'med', 'lr', 'n', 'level']
    """
    return jcall("kmrec", time=time, event=event, grp=grp, unit=unit, end=end, level=level)


def knxnr(x: Any, y: Any = None, t: Any = None, ds: Any = None, dt: Any = None, nsim: Any = 999,
          seed: Any = None, lonlat: Any = False) -> Any:
    """Knox test of near repeats.

    The arguments and the result are those of R's ``knxnr()``: see its help
    page.

    Examples
    --------
    >>> r = knxnr(x=[[1830, 759], [1874, 872], [572, 75], [1661, 1947], [1283, 864], [1038,
    ... 1915], [1473, 1776], [269, 1280], [1314, 1942], [1410, 1238], [915, 667], [1438, 693],
    ... [1869, 797], [511, 1569], [925, 78], [1880, 1498], [1956, 1355], [235, 343], [950, 522],
    ... [1121, 1029], [1808, 1351], [277, 1966], [1978, 1519], [1893, 1133], [165, 1699], [1028,
    ... 379], [780, 543], [1811, 1656], [894, 1386], [1672, 481], [1475, 86], [1622, 281], [776,
    ... 433], [1370, 959], [8, 395], [1666, 1439], [15, 16], [415, 751], [1813, 1029], [1224,
    ... 3]], t=[130, 150, 209, 215, 263, 144, 336, 351, 85, 264, 330, 220, 231, 342, 310, 212,
    ... 300, 42, 279, 228, 54, 29, 169, 284, 268, 298, 62, 345, 107, 54, 263, 118, 284, 144,
    ... 248, 283, 69, 11, 50, 248], ds=300, dt=30, nsim=0)
    >>> list(r)
    ['tab', 'n', 'pairs', 'nsim']
    """
    return jcall("knxnr", x=x, y=y, t=t, ds=ds, dt=dt, nsim=nsim, seed=seed, lonlat=lonlat)


def kscan(x: Any, y: Any = None, t: Any = None, r: Any = None, tw: Any = None, minc: Any = 2,
          top: Any = 1, nsim: Any = 999, seed: Any = None, lonlat: Any = False) -> Any:
    """Space-time scan statistic of incidents.

    The arguments and the result are those of R's ``kscan()``: see its help
    page.

    Examples
    --------
    >>> r = kscan(x=[[17.883, 0.069], [4.971, 10.623], [1.674, 10.486], [0.973, 4.264], [19.596,
    ... 14.339], [9.683, 19.227], [16.908, 10.365], [8.326, 3.491], [9.787, 11.251], [3.666,
    ... 15.185], [15.183, 13.339], [6.103, 4.497], [3.314, 6.917], [0.656, 6.397], [2.73,
    ... 18.098], [3.543, 3.984], [10.391, 13.619], [16.222, 2.75], [2.307, 2.14], [17.868,
    ... 1.857], [11.507, 18.329], [2.931, 5.541], [18.056, 17.716], [5.06, 15.457], [3.012,
    ... 15.901], [15.371, 4.113], [4.602, 0.964], [6.108, 0.776], [10.371, 5.691], [6.692,
    ... 6.976], [4.695, 15.191], [4.885, 15.013], [5.262, 15.104], [5.291, 15.738], [5.115,
    ... 14.754], [4.445, 14.366], [4.984, 15.082], [5.319, 14.794], [5.244, 15.134], [4.943,
    ... 14.756]], t=[92, 83, 8, 53, 9, 70, 91, 81, 24, 73, 72, 18, 78, 45, 2, 21, 60, 19, 41,
    ... 97, 86, 16, 81, 82, 45, 16, 93, 73, 36, 44, 41, 42, 42, 41, 41, 41, 43, 40, 41, 44],
    ... r=2, tw=10, nsim=0, minc=3)
    >>> list(r)
    ['cl', 'memb', 'n', 'nsim']
    """
    return jcall("kscan", x=x, y=y, t=t, r=r, tw=tw, minc=minc, top=top, nsim=nsim, seed=seed,
                  lonlat=lonlat)


def lftab(breaks: Any, events: Any, cens: Any = 0, level: Any = 0.95, ci: Any = "log-log") -> Any:
    """Actuarial life table of grouped times.

    The arguments and the result are those of R's ``lftab()``: see its help
    page.

    Examples
    --------
    >>> r = lftab(breaks=[0, 3, 6, 12, 24, float("inf")], events=[40, 30, 20, 10, 5], cens=[2,
    ... 2, 3, 3, 4])
    >>> list(r)
    ['tab', 'median', 'ci']
    """
    return jcall("lftab", breaks=breaks, events=events, cens=cens, level=level, ci=ci)


def lostay(days: Any, w: Any = None, breaks: Any = None, thr: Any = None, flag: Any = None,
           level: Any = 0.95) -> Any:
    """Lengths of stay: their distribution and its fit.

    The arguments and the result are those of R's ``lostay()``: see its help
    page.

    Examples
    --------
    >>> r = lostay(days=[3, 8, 20, 25, 40, 2], w=[2, 1, 1, 3, 1, 1], flag=[True, False, True,
    ... False, True, False], breaks=[0, 10, 50], thr=[5, 30], level=0.9)
    >>> list(r)
    ['stay', 'bands', 'thr', 'fit', 'par']
    """
    return jcall("lostay", days=days, w=w, breaks=breaks, thr=thr, flag=flag, level=level)


def mndla(days: Any, id: Any = None, by: Any = None, w: Any = None, thr: Any = 15,
          hours: Any = None, miss: Any = None, alert: Any = None, amin: Any = 2, elig: Any = None,
          excl: Any = None, level: Any = 0.95) -> Any:
    """Segregation placements under the Mandela Rules.

    The arguments and the result are those of R's ``mndla()``: see its help
    page.

    Examples
    --------
    >>> r = mndla(days=[5, 20, 30], hours=[1, 1.5, 3], miss=[100, 100, 40])
    >>> list(r)
    ['tab', 'cls', 'alerts', 'rule']
    """
    return jcall("mndla", days=days, id=id, by=by, w=w, thr=thr, hours=hours, miss=miss,
                  alert=alert, amin=amin, elig=elig, excl=excl, level=level)


def morsta(counts: Any, crit: Any = 0.05) -> Any:
    """Morisita's index of dispersion.

    The arguments and the result are those of R's ``morsta()``: see its help
    page.

    Examples
    --------
    >>> r = morsta(counts=[3, 0, 7, 12, 1, 4, 0, 9])
    >>> list(r)
    ['imor', 'mclu', 'muni', 'imst', 'chi', 'p']
    """
    return jcall("morsta", counts=counts, crit=crit)


def mvhwk(times: Any, type: Any, hor: Any, beta: Any = None, start: Any = None,
          maxit: Any = 500) -> Any:
    """Multivariate Hawkes process: types that excite each other.

    The arguments and the result are those of R's ``mvhwk()``: see its help
    page.

    Examples
    --------
    >>> r = mvhwk(times=[2.554, 7.062, 11.256, 13.158, 13.615, 16.466, 16.955, 21.318, 24.296,
    ... 30.104, 31.279, 32.149, 37.06, 41.255, 41.478, 51.995, 52.68, 53.769, 55.16, 55.244,
    ... 61.538, 62.168, 63.948, 65.733, 73.003, 74.389, 81.565, 84.657, 86.932, 94.857],
    ... type=["a", "b", "a", "b", "a", "b", "a", "b", "a", "b", "a", "b", "a", "b", "a", "b",
    ... "a", "b", "a", "b", "a", "b", "a", "b", "a", "b", "a", "b", "a", "b"], hor=100)
    >>> list(r)
    ['par', 'se', 'cov', 'mu', 'a', 'beta', 'rad', 'loglik', 'aic', 'conv', 'iter', 'n',
    'hor']
    """
    return jcall("mvhwk", times=times, type=type, hor=hor, beta=beta, start=start, maxit=maxit)


def pbbnd(x: Any, rnd: Any = None, kind: Any = "nearest", suplim: Any = None, supmar: Any = ["x",
          "X", "s", "..", "F", "suppressed", "n/a"], prev: Any = None, pop: Any = None,
          prepop: Any = None, per: Any = 1000) -> Any:
    """Bounds implied by rounding and suppression in a published table.

    The arguments and the result are those of R's ``pbbnd()``: see its help
    page.

    Examples
    --------
    >>> r = pbbnd(x=["1,240", "35", "<5", "x"], rnd=5, suplim=10)
    >>> list(r)
    ['bnd', 'chg']
    """
    return jcall("pbbnd", x=x, rnd=rnd, kind=kind, suplim=suplim, supmar=supmar, prev=prev, pop=pop,
                  prepop=prepop, per=per)


def pinpl(x: Any, y: Any = None, poly: Any = None, pop: Any = None) -> Any:
    """Points in polygons, and the population the units cover.

    The arguments and the result are those of R's ``pinpl()``: see its help
    page.

    Examples
    --------
    >>> r = pinpl(x=[[0.5, 0.5], [1.5, 0.5], [2.5, 0.2], [4, 0.5]], poly={"x": [0, 1, 1, 0, 1,
    ... 2, 2, 1, 2, 3, 3, 2], "y": [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1], "id": ["A", "A", "A",
    ... "A", "B", "B", "B", "B", "C", "C", "C", "C"]})
    >>> list(r)
    ['reg', 'nreg', 'regs', 'cov']
    """
    return jcall("pinpl", x=x, y=y, poly=poly, pop=pop)


def qdrtst(x: Any, y: Any = None, win: Any = None, nx: Any = None, ny: Any = None) -> Any:
    """Quadrat count test.

    The arguments and the result are those of R's ``qdrtst()``: see its help
    page.

    Examples
    --------
    >>> r = qdrtst(x=[[0.279, 0.624], [0.733, 0.802], [0.314, 0.654], [0.896, 0.679], [0.972,
    ... 0.411], [1.9, 0.663], [1.81, 0.027], [1.88, 0.816], [1.547, 0.075], [0.184, 0.901],
    ... [0.032, 0.071], [0.667, 0.191], [0.484, 0.777], [1.584, 0.203], [0.993, 0.404], [1.028,
    ... 0.446], [0.145, 0.606], [0.69, 0.576], [1.513, 0.06], [0.77, 0.183], [0.014, 0.296],
    ... [1.437, 0.376], [1.114, 0.713], [0.14, 0.125], [0.696, 0.429], [0.826, 0.03], [0.623,
    ... 0.454], [1.475, 0.794], [0.777, 0.471], [1.239, 0.236]], win=[0, 2, 0, 1], nx=3)
    >>> list(r)
    ['counts', 'stat']
    """
    return jcall("qdrtst", x=x, y=y, win=win, nx=nx, ny=ny)


def sarml(y: Any, x: Any, w: Any, model: Any = "lag", durbin: Any = False, method: Any = "ml",
          int: Any = True) -> Any:
    """Spatial regression: lag, error and combined models.

    The arguments and the result are those of R's ``sarml()``: see its help
    page.

    Examples
    --------
    >>> r = sarml(y=[-0.987, -0.03, 3.668, 0.758, 3.462, 0.386, 5.335, 4.206, -2.171, 2.116,
    ... 0.91, 1.111, 0.658, 1.598, 3.914, 5.665, -2.305, 2.894, 0.331, -1.417], x=[-0.55,
    ... -0.279, 1.097, 0.442, 0.241, -0.256, 0.931, 1.335, -0.869, 0.055, 0.049, -0.578, -0.999,
    ... -0.002, 0.656, 1.477, -1.909, -0.702, -0.311, -1.663], w=[[0, 0, 0, 0.3333333333333333,
    ... 0, 0.3333333333333333, 0, 0, 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0,
    ... 0, 0, 0, 0, 0.3333333333333333, 0, 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0,
    ... 0.3333333333333333, 0, 0, 0], [0, 0, 0, 0, 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0,
    ... 0.3333333333333333, 0.3333333333333333, 0, 0, 0, 0, 0], [0, 0, 0.3333333333333333, 0, 0,
    ... 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0, 0.3333333333333333, 0, 0, 0, 0, 0], [0, 0,
    ... 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3333333333333333,
    ... 0.3333333333333333, 0, 0, 0, 0], [0.3333333333333333, 0, 0, 0.3333333333333333, 0, 0, 0,
    ... 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ... 0, 0, 0, 0, 0, 0.3333333333333333, 0, 0.3333333333333333, 0, 0.3333333333333333], [0,
    ... 0.3333333333333333, 0, 0, 0, 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ... 0.3333333333333333, 0, 0, 0], [0.3333333333333333, 0, 0, 0.3333333333333333, 0,
    ... 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0,
    ... 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3333333333333333, 0,
    ... 0.3333333333333333], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3333333333333333, 0,
    ... 0.3333333333333333, 0, 0, 0, 0, 0.3333333333333333, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ... 0.3333333333333333, 0, 0, 0.3333333333333333, 0, 0, 0, 0, 0.3333333333333333, 0], [0, 0,
    ... 0, 0, 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3333333333333333, 0,
    ... 0.3333333333333333, 0, 0], [0, 0, 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0,
    ... 0.3333333333333333, 0, 0, 0.3333333333333333, 0, 0, 0, 0, 0], [0, 0, 0.3333333333333333,
    ... 0, 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0, 0.3333333333333333, 0, 0, 0, 0, 0, 0],
    ... [0, 0, 0, 0, 0.3333333333333333, 0, 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ... 0.3333333333333333, 0, 0], [0, 0.3333333333333333, 0, 0, 0, 0, 0, 0.3333333333333333, 0,
    ... 0, 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0,
    ... 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0, 0.3333333333333333, 0, 0, 0,
    ... 0.3333333333333333], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3333333333333333,
    ... 0.3333333333333333, 0, 0.3333333333333333, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0,
    ... 0.3333333333333333, 0, 0, 0.3333333333333333, 0, 0, 0, 0, 0, 0, 0, 0.3333333333333333,
    ... 0, 0]], method="2sls")
    >>> list(r)
    ['coef', 'fit', 'vcov', 'model', 'method', 'durbin', 'nlag', 'w']
    """
    return jcall("sarml", y=y, x=x, w=w, model=model, durbin=durbin, method=method, int=int)


def scrdsp(score: Any, group: Any = None, ref: Any = None, conf: Any = 0.95) -> Any:
    """A risk score's spread across groups.

    The arguments and the result are those of R's ``scrdsp()``: see its help
    page.

    Examples
    --------
    >>> r = scrdsp(score=[0.015, 0.873, -1.028, 0.686, 0.449, 0.407, 0.285, -0.624, 1.34, 0.811,
    ... 0.975, 0.4, 0.704, 0.406, 0.6, 0.232, 2.346, 0.854, 1.05, 1.834, 0.709, -0.085],
    ... group=["x", "x", "x", "x", "x", "x", "x", "x", "y", "y", "y", "y", "y", "y", "y", "y",
    ... "z", "z", "z", "z", "z", "z"])
    >>> list(r)
    ['tab', 'anova', 'welch', 'spread', 'ref']
    """
    return jcall("scrdsp", score=score, group=group, ref=ref, conf=conf)


def sdbsim(n: Any, steps: Any, mode: Any = "agent", a0: Any = 1 / 30, omega: Any = 1 / 15,
           eta: Any = 0.03, gam: Any = 0.002, theta: Any = 1, dt: Any = None, b: Any = None,
           rho: Any = None, dx: Any = 1, dif: Any = 1, eps: Any = 1, z: Any = 4, sup: Any = None,
           seed: Any = None) -> Any:
    """Hot-spot dynamics: the Short-D'Orsogna-Brantingham model.

    The arguments and the result are those of R's ``sdbsim()``: see its help
    page.

    Examples
    --------
    >>> r = sdbsim(n=6, steps=5, a0=0.03333333333333333, omega=0.06666666666666667, eta=0.03,
    ... gam=0.02, theta=0.5, seed=1)
    >>> list(r)
    ['b', 'rho', 'tot', 'burg', 'meanb', 'mode']
    """
    return jcall("sdbsim", n=n, steps=steps, mode=mode, a0=a0, omega=omega, eta=eta, gam=gam,
                  theta=theta, dt=dt, b=b, rho=rho, dx=dx, dif=dif, eps=eps, z=z, sup=sup,
                  seed=seed)


def sdbthr(a0: Any, omega: Any, eps: Any, dif: Any, gam: Any, eta: Any, z: Any = 1,
           k: Any = None) -> Any:
    """Hot-spot threshold of the Short-D'Orsogna-Brantingham model.

    The arguments and the result are those of R's ``sdbthr()``: see its help
    page.

    Examples
    --------
    >>> r = sdbthr(a0=1, omega=1, eps=2, dif=1, gam=1, eta=0.01)
    >>> list(r)
    ['steady', 'p', 'q', 'r', 'thr', 'unstable', 'kcrit', 'detcrit', 'grcrit', 'disp']
    """
    return jcall("sdbthr", a0=a0, omega=omega, eps=eps, dif=dif, gam=gam, eta=eta, z=z, k=k)


def sdidx(y: Any, unit: Any, time: Any, first: Any, reps: Any = 200, seed: Any = None,
          level: Any = 0.95) -> Any:
    """Synthetic difference in differences.

    The arguments and the result are those of R's ``sdidx()``: see its help
    page.

    Examples
    --------
    >>> r = sdidx(y=[21.8415, 22.1593, 21.6411, 20.9932, 21.0411, 18.9706, 20.157, 20.7394,
    ... 22.4121, 21.756, 21.6, 22.3634, 23.6202, 24.4906, 24.4503, 23.8121, 22.0386, 22.599,
    ... 23.8499, 24.9629, 25.2367, 24.7411, 24.2538, 24.5444, 23.8676, 25.1626, 25.7564,
    ... 25.4709, 24.9364, 25.012, 25.996, 27.3514, 25.9999, 25.9791, 25.4718, 25.3582, 26.1565,
    ... 27.5464, 28.6638, 28.8951, 25.8414, 25.5835, 26.1682, 27.5177, 28.8509, 29.4018,
    ... 29.1236, 28.7317], unit=["A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B",
    ... "B", "B", "B", "B", "C", "C", "C", "C", "C", "C", "C", "C", "D", "D", "D", "D", "D",
    ... "D", "D", "D", "E", "E", "E", "E", "E", "E", "E", "E", "F", "F", "F", "F", "F", "F",
    ... "F", "F"], time=[2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2010, 2011, 2012, 2013,
    ... 2014, 2015, 2016, 2017, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2010, 2011,
    ... 2012, 2013, 2014, 2015, 2016, 2017, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017,
    ... 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017], first=[2015, 2015, 2015, 2015, 2015,
    ... 2015, 2015, 2015, None, None, None, None, None, None, None, None, None, None, None,
    ... None, None, None, None, None, None, None, None, None, None, None, None, None, None,
    ... None, None, None, None, None, None, None, None, None, None, None, None, None, None,
    ... None], reps=20, seed=1)
    >>> list(r)
    ['est', 'meth', 'unitw', 'timew', 'path', 'level']
    """
    return jcall("sdidx", y=y, unit=unit, time=time, first=first, reps=reps, seed=seed, level=level)


def sirex(x: Any, n: Any, strat: Any, grp: Any, rates: Any = None, level: Any = 0.95) -> Any:
    """Indirect standardisation: standardised ratios with shrinkage.

    The arguments and the result are those of R's ``sirex()``: see its help
    page.

    Examples
    --------
    >>> r = sirex(x=[12, 3, 30, 9, 4, 1], n=[4000, 2000, 8000, 3000, 1500, 900], strat=["young",
    ... "old", "young", "old", "young", "old"], grp=["N", "N", "S", "S", "E", "E"])
    >>> list(r)
    ['sir', 'prior', 'level']
    """
    return jcall("sirex", x=x, n=n, strat=strat, grp=grp, rates=rates, level=level)


def spill(fit: Any = None, w: Any = None, rho: Any = None, b: Any = None, th: Any = None) -> Any:
    """Direct, indirect and total effects of a spatial lag model.

    The arguments and the result are those of R's ``spill()``: see its help
    page.

    Examples
    --------
    >>> r = spill(w={"n": 20, "i": [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7,
    ... 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15,
    ... 15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20], "j": [4, 6, 9, 6, 8,
    ... 17, 5, 14, 15, 3, 6, 15, 3, 15, 16, 1, 4, 8, 16, 18, 20, 2, 6, 17, 1, 4, 6, 7, 18, 20,
    ... 12, 14, 19, 11, 14, 19, 5, 16, 18, 3, 12, 15, 3, 5, 14, 5, 7, 18, 2, 8, 11, 7, 16, 20,
    ... 11, 12, 14, 7, 10, 18], "d": [0.366742416417845, 0.28800173610587837,
    ... 0.19442736432920127, 0.3944641935588071, 0.18802393464662948, 0.2251421773013666,
    ... 0.13849548729110273, 0.11827087553578014, 0.05575840743780257, 0.27614670014323905,
    ... 0.11901680553602503, 0.2229125389025929, 0.13849548729110273, 0.18106904760339354,
    ... 0.11417968295629481, 0.28800173610587837, 0.11901680553602503, 0.23700843866833102,
    ... 0.18056023925549064, 0.13892803892663277, 0.0458802789878178, 0.18802393464662948,
    ... 0.23700843866833102, 0.2545270123189285, 0.19442736432920127, 0.329802971484491,
    ... 0.31624199594614244, 0.24263758983306769, 0.261663906567184, 0.2241428116179505,
    ... 0.22227235545609356, 0.23221541723150083, 0.049648766349225645, 0.22227235545609356,
    ... 0.1661715980545412, 0.18404890654388575, 0.2651886121235224, 0.25132051249350895,
    ... 0.3528639964632266, 0.11827087553578014, 0.1661715980545412, 0.15305227865013962,
    ... 0.05575840743780257, 0.18106904760339354, 0.15305227865013962, 0.11417968295629481,
    ... 0.18056023925549064, 0.1494824404403407, 0.2251421773013666, 0.2545270123189285,
    ... 0.26475082624989105, 0.13892803892663277, 0.1494824404403407, 0.17768511473953016,
    ... 0.049648766349225645, 0.18404890654388575, 0.18309833423600555, 0.0458802789878178,
    ... 0.2241428116179505, 0.17768511473953016], "w": [0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
    ... 0.3333333333333333, 0.3333333333333333]}, rho=0.4, b={"gent": 0.1}, th=0.3)
    >>> list(r)
    ['term', 'direct', 'indirect', 'total', 'sed', 'sei', 'set']
    """
    return jcall("spill", fit=fit, w=w, rho=rho, b=b, th=th)


def stgdid(y: Any, unit: Any, time: Any, first: Any = None, trt: Any = None, meth: Any = "csa",
           ctrl: Any = "notyet", adj: Any = None, est: Any = "dr", clus: Any = None,
           reps: Any = 499, seed: Any = None, level: Any = 0.95) -> Any:
    """Effects of a policy adopted at different times.

    The arguments and the result are those of R's ``stgdid()``: see its help
    page.

    Examples
    --------
    >>> r = stgdid(y=[31.8415, 32.4093, 32.1411, 27.7432, 28.0411, 29.2206, 32.657, 33.4894,
    ... 33.4121, 28.956, 29, 29.9634, 33.4202, 34.4906, 34.6503, 30.2121, 30.0386, 30.749,
    ... 34.1499, 35.4129, 35.8367, 35.4911, 31.1538, 31.5944, 34.8676, 36.2626, 36.9564,
    ... 36.7709, 32.3364, 32.512, 35.596, 37.0514, 37.9999, 38.0291, 37.5718, 37.5082, 36.3565,
    ... 37.7964, 38.9638, 39.2451, 38.8414, 38.5835, 37.1682, 38.5177, 39.8509, 40.4018,
    ... 40.1236, 39.7317], unit=["A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B",
    ... "C", "C", "C", "C", "C", "C", "D", "D", "D", "D", "D", "D", "E", "E", "E", "E", "E",
    ... "E", "F", "F", "F", "F", "F", "F", "G", "G", "G", "G", "G", "G", "H", "H", "H", "H",
    ... "H", "H"], time=[2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020,
    ... 2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016,
    ... 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018,
    ... 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020], trt=[0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1,
    ... 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ... 0, 0, 0, 0, 0, 0, 0, 0], meth="dcdh", reps=0)
    >>> list(r)
    ['att', 'overall', 'dynamic', 'meth', 'ctrl', 'reps', 'level']
    """
    return jcall("stgdid", y=y, unit=unit, time=time, first=first, trt=trt, meth=meth, ctrl=ctrl,
                  adj=adj, est=est, clus=clus, reps=reps, seed=seed, level=level)


def sthwk(times: Any, x: Any, y: Any, hor: Any, win: Any = None, par: Any = None, at: Any = None,
          start: Any = None, maxit: Any = 500) -> Any:
    """Space-time Hawkes process.

    The arguments and the result are those of R's ``sthwk()``: see its help
    page.

    Examples
    --------
    >>> r = sthwk(times=[0.291, 0.674, 3.968, 4.298, 4.791, 9.136, 15.968, 16.367, 17.374,
    ... 20.606, 21.86, 26.067, 37.036, 42.381, 48.259, 55.174, 57.689, 57.762, 65.634, 66.902,
    ... 72.783, 86.366, 87.618, 91.117, 91.568, 92.343, 93.022, 93.628, 94.745, 96.139],
    ... x=[5.693, 2.573, 3.61, 4.35, 0.458, 9.092, 9.699, 0.498, 0.744, 0.78, 8.689, 9.458,
    ... 4.793, 5.385, 8.48, 5.345, 6.539, 3.609, 2.831, 5.23, 1.811, 4.912, 1.408, 6.66, 8.714,
    ... 3.395, 0.216, 9.304, 4.522, 9.998], y=[7.844, 2.65, 5.634, 1.992, 2.379, 5.466, 8.601,
    ... 5.197, 9.54, 3.039, 1.794, 2.748, 4.755, 2.586, 1.406, 2.705, 0.539, 3.712, 7.846,
    ... 0.769, 1.196, 4.889, 2.554, 0.142, 6.397, 6.044, 4.677, 5.382, 0.224, 0.98], hor=100,
    ... win=[0, 10, 0, 10])
    >>> list(r)
    ['par', 'se', 'cov', 'loglik', 'aic', 'conv', 'iter', 'fit', 'n', 'hor', 'win', 'area',
    'lam']
    """
    return jcall("sthwk", times=times, x=x, y=y, hor=hor, win=win, par=par, at=at, start=start,
                  maxit=maxit)


def stkflw(people: Any = None, days: Any = None, adp: Any = None, alos: Any = None, t: Any = 365,
           prd: Any = None, from_: Any = None, to: Any = None, expo: Any = None, per: Any = 1e+05,
           base: Any = "first", rel: Any = None, early: Any = None) -> Any:
    """Stock and flow of a population held over time.

    The arguments and the result are those of R's ``stkflw()``: see its help
    page.

    Examples
    --------
    >>> r = stkflw(adp=900, alos=61)
    >>> list(r)
    ['flow', 'mean', 'base', 'per']
    """
    return jcall("stkflw", people=people, days=days, adp=adp, alos=alos, t=t, prd=prd, from_=from_,
                  to=to, expo=expo, per=per, base=base, rel=rel, early=early)


def suplrn(x: Any, y: Any, fam: Any = "gaussian", lrns: Any = ["glm", "lasso", "forest"],
           nfold: Any = 5, fold: Any = None, disc: Any = False, newx: Any = None, alpha: Any = 0.5,
           nlam: Any = 100, nfin: Any = 5, ntree: Any = 500, mtry: Any = None,
           seed: Any = 1) -> Any:
    """Stacking: the super learner.

    The arguments and the result are those of R's ``suplrn()``: see its help
    page.

    Examples
    --------
    >>> r = suplrn(x=[[-0.794, 0.882, 0.105], [0.197, 0.899, -1.049], [1.002, 1.503, 0.889],
    ... [1.289, 0.505, -0.556], [0.906, -0.341, 0.345], [0.494, 1.373, 0.557], [0.599, 0.819,
    ... 0.139], [-1.58, 1.426, -0.097], [1.001, 2.196, 0.568], [2.188, -0.896, 0.052], [-1.209,
    ... -1.69, -1.457], [-0.587, -0.467, 0.873], [1.056, -0.36, -0.178], [-0.317, 0.746,
    ... -0.345], [-0.055, 0.169, 0.383], [0.33, 1.65, -1.602], [0.663, 0.033, 1.214], [0.878,
    ... -1.032, 0.588], [0.203, -0.632, -0.594], [2.274, 0.578, -0.943], [-0.899, 1.842, 1.262],
    ... [2.139, -0.058, 0.267], [-1.166, -0.257, 1.31], [-0.047, -2.167, -0.862], [0.416, 0.11,
    ... 0.714], [1.72, 1.315, 0.093], [-0.784, -0.893, -0.831], [-1.304, -1.044, 0.781],
    ... [-0.452, -0.542, -0.018], [-1.774, -0.331, -0.373], [-0.05, 1.672, 0.695], [-1.037,
    ... -1.053, -1.264], [0.622, 0.081, 0.378], [0.127, 0.703, -0.534], [-1.741, 1.511, -0.812],
    ... [-0.573, 0.001, -1.479], [0.912, 1.431, -0.251], [1.234, -1.112, 0.203], [0.661, 0.568,
    ... -0.176], [-0.861, 1.338, -0.812], [-0.533, 0.073, -0.117], [0.321, 1.196, 1.274],
    ... [-0.004, -2.112, -1.196], [1.598, 0.671, -0.871], [1.327, -0.186, 0.434], [-0.119,
    ... -0.192, -0.426], [0.082, 0.787, 0.322], [-0.02, 0.562, 0.227], [-0.825, -0.019, -1.156],
    ... [0.702, -0.26, -0.321], [-0.196, -0.162, -1.072], [-0.158, -2.09, -0.997], [0.224,
    ... 0.449, 0.864], [0.305, -0.697, 0.241], [1.369, 0.782, -0.591], [-1.389, -0.265, -1.827],
    ... [-1.644, 0.02, -0.339], [1.426, 0.401, 1.611], [0.097, 1.967, -0.664], [0.312, 0.954,
    ... -0.565]], y=[0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0,
    ... 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1,
    ... 1, 0, 0, 1, 0, 0], fam="binomial", lrns=["glm", "forest"], fold=[2, 3, 2, 3, 1, 1, 2, 2,
    ... 3, 2, 3, 2, 2, 2, 1, 3, 3, 3, 3, 1, 1, 3, 1, 3, 1, 1, 1, 2, 2, 1, 1, 3, 3, 1, 2, 1, 2,
    ... 1, 1, 2, 3, 3, 2, 1, 2, 1, 2, 2, 3, 1, 3, 2, 1, 3, 3, 1, 2, 3, 3, 2], ntree=20,
    ... newx=[[-0.794, 0.882, 0.105], [0.197, 0.899, -1.049], [1.002, 1.503, 0.889]], seed=2)
    >>> list(r)
    ['wt', 'risk', 'fold', 'oof', 'pred']
    """
    return jcall("suplrn", x=x, y=y, fam=fam, lrns=lrns, nfold=nfold, fold=fold, disc=disc,
                  newx=newx, alpha=alpha, nlam=nlam, nfin=nfin, ntree=ntree, mtry=mtry, seed=seed)


def svymx(y: Any, w: Any = None, strata: Any = None, psu: Any = None, fpc: Any = None,
          stat: Any = "mean", den: Any = None, by: Any = None, ref: Any = None,
          level: Any = 0.95) -> Any:
    """Design-based estimates from a complex survey.

    The arguments and the result are those of R's ``svymx()``: see its help
    page.

    Examples
    --------
    >>> r = svymx(y=[141, 68, 78, 90, 145, 70, 63, 101, 57, 108, 127, 81, 55, 114, 113, 132,
    ... 128, 105, 53, 75, 129, 78, 60, 102, 72, 101, 142, 79, 91, 139, 123, 129, 57, 119, 110,
    ... 129, 88, 67, 130, 119], den=[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    ... 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2], stat="ratio", fpc=[1000,
    ... 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
    ... 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
    ... 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000], level=0.9)
    >>> list(r)
    ['est', 'cmp', 'des', 'stat', 'level']
    """
    return jcall("svymx", y=y, w=w, strata=strata, psu=psu, fpc=fpc, stat=stat, den=den, by=by,
                  ref=ref, level=level)


def synthc(y: Any, unit: Any, time: Any, first: Any, plac: Any = True) -> Any:
    """Synthetic control for one treated unit.

    The arguments and the result are those of R's ``synthc()``: see its help
    page.

    Examples
    --------
    >>> r = synthc(y=[21.8415, 22.1593, 21.6411, 20.9932, 21.0411, 18.9706, 20.157, 20.7394,
    ... 22.4121, 21.756, 21.6, 22.3634, 23.6202, 24.4906, 24.4503, 23.8121, 22.0386, 22.599,
    ... 23.8499, 24.9629, 25.2367, 24.7411, 24.2538, 24.5444, 23.8676, 25.1626, 25.7564,
    ... 25.4709, 24.9364, 25.012, 25.996, 27.3514, 25.9999, 25.9791, 25.4718, 25.3582, 26.1565,
    ... 27.5464, 28.6638, 28.8951, 25.8414, 25.5835, 26.1682, 27.5177, 28.8509, 29.4018,
    ... 29.1236, 28.7317], unit=["A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B",
    ... "B", "B", "B", "B", "C", "C", "C", "C", "C", "C", "C", "C", "D", "D", "D", "D", "D",
    ... "D", "D", "D", "E", "E", "E", "E", "E", "E", "E", "E", "F", "F", "F", "F", "F", "F",
    ... "F", "F"], time=[2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2010, 2011, 2012, 2013,
    ... 2014, 2015, 2016, 2017, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2010, 2011,
    ... 2012, 2013, 2014, 2015, 2016, 2017, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017,
    ... 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017], first=[2015, 2015, 2015, 2015, 2015,
    ... 2015, 2015, 2015, None, None, None, None, None, None, None, None, None, None, None,
    ... None, None, None, None, None, None, None, None, None, None, None, None, None, None,
    ... None, None, None, None, None, None, None, None, None, None, None, None, None, None,
    ... None])
    >>> list(r)
    ['weights', 'path', 'fit', 'plac', 'unit']
    """
    return jcall("synthc", y=y, unit=unit, time=time, first=first, plac=plac)


def trndrr(y: Any, x: Any = None, expo: Any = None, seas: Any = None, alt: Any = "both",
           level: Any = 0.95) -> Any:
    """Trend in a series: rate ratio, Mann-Kendall and Theil-Sen.

    The arguments and the result are those of R's ``trndrr()``: see its help
    page.

    Examples
    --------
    >>> r = trndrr(y=[0, 0, 5])
    >>> list(r)
    ['rr', 'mk', 'sen', 'yoy', 'seas', 'fitted', 'n', 'alt', 'level', 'warnings']
    """
    return jcall("trndrr", y=y, x=x, expo=expo, seas=seas, alt=alt, level=level)


def trnmat(state: Any, id: Any, time: Any = None) -> Any:
    """Transition matrix between states.

    The arguments and the result are those of R's ``trnmat()``: see its help
    page.

    Examples
    --------
    >>> r = trnmat(state=["a", "b", "a", "b"], id=[1, 1, 2, 2])
    >>> list(r)
    ['count', 'prob', 'pi', 'summ']
    """
    return jcall("trnmat", state=state, id=id, time=time)


def wcbx(y: Any, unit: Any, time: Any, trt: Any = None, adj: Any = None, clus: Any = None,
         meth: Any = "wild", reps: Any = None, wtype: Any = "rad", ties: Any = "count",
         stat: Any = "t", ntrt: Any = 1, adopt: Any = None, eff: Any = None, alpha: Any = 0.05,
         power: Any = 0.8, seed: Any = None, level: Any = 0.95) -> Any:
    """Inference for a difference in differences with few clusters.

    The arguments and the result are those of R's ``wcbx()``: see its help
    page.

    Examples
    --------
    >>> r = wcbx(y=[31.8415, 32.4093, 32.1411, 27.7432, 28.0411, 29.2206, 32.657, 33.4894,
    ... 33.4121, 28.956, 29, 29.9634, 33.4202, 34.4906, 34.6503, 30.2121, 30.0386, 30.749,
    ... 34.1499, 35.4129, 35.8367, 35.4911, 31.1538, 31.5944, 34.8676, 36.2626, 36.9564,
    ... 36.7709, 32.3364, 32.512, 35.596, 37.0514, 37.9999, 38.0291, 37.5718, 37.5082, 36.3565,
    ... 37.7964, 38.9638, 39.2451, 38.8414, 38.5835, 37.1682, 38.5177, 39.8509, 40.4018,
    ... 40.1236, 39.7317], unit=["A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B",
    ... "C", "C", "C", "C", "C", "C", "D", "D", "D", "D", "D", "D", "E", "E", "E", "E", "E",
    ... "E", "F", "F", "F", "F", "F", "F", "G", "G", "G", "G", "G", "G", "H", "H", "H", "H",
    ... "H", "H"], time=[2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020,
    ... 2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016,
    ... 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020, 2015, 2016, 2017, 2018,
    ... 2019, 2020, 2015, 2016, 2017, 2018, 2019, 2020], meth="power", reps=20, adopt=2018,
    ... ntrt=2, eff=[0, 2, 5], seed=1)
    >>> list(r)
    ['res', 'meth', 'ncl', 'draws', 'enum', 'mde', 'size', 'level']
    """
    return jcall("wcbx", y=y, unit=unit, time=time, trt=trt, adj=adj, clus=clus, meth=meth,
                  reps=reps, wtype=wtype, ties=ties, stat=stat, ntrt=ntrt, adopt=adopt, eff=eff,
                  alpha=alpha, power=power, seed=seed, level=level)


def wystep(obs: Any, null: Any) -> Any:
    """Westfall-Young stepdown p-values.

    The arguments and the result are those of R's ``wystep()``: see its help
    page.

    Examples
    --------
    >>> r = wystep(obs=[2.5, 0.4, -2.1], null=[[-0.254, 0.876, 2.001], [-1.056, 0.639, 0.286],
    ... [0.686, -0.138, -0.032], [0.025, -0.501, -0.174], [-1.672, 0.067, -0.345], [-1.484,
    ... 2.346, -0.896], [0.43, -2.015, -1.157], [0.008, 1.783, 1.257], [0.886, 1.026, -0.123],
    ... [-0.41, -0.973, 1.191], [0.061, 1.255, -0.078], [-1.166, 0.45, -0.858], [-0.013, -1.785,
    ... 0.971], [1.196, 0.841, -0.233], [-0.021, 0.379, -0.21], [0.042, 0.373, 0.032], [0.699,
    ... -0.975, -0.039], [-0.071, -1.768, 0.033], [0.98, -0.863, -0.069], [0.479, -1.289,
    ... 0.298], [-0.071, -1.063, -1.381], [-1.29, 0.813, -0.297], [-0.009, -0.397, -0.847],
    ... [-1.202, -0.096, 0.05], [0.895, 1.374, -2.064], [-0.286, 0.352, 1.6], [-0.014, -0.279,
    ... 0.561], [-0.692, 1.804, 1.853], [1.128, 0.355, 0.072], [0.217, -0.042, -0.668]])
    >>> list(r)
    ['stat', 'p', 'padj']
    """
    return jcall("wystep", obs=obs, null=null)


def xcser(x: Any, y: Any, lagmax: Any = None, dif: Any = False, level: Any = 0.95) -> Any:
    """Correlation of two series, at lags.

    The arguments and the result are those of R's ``xcser()``: see its help
    page.

    Examples
    --------
    >>> r = xcser(x=[1, 3, 2, 5, 4, 6, 8, 7], y=[2, 2, 4, 3, 6, 5, 7, 9], lagmax=2, level=0.9)
    >>> list(r)
    ['cor', 'ccf', 'common', 'level']
    """
    return jcall("xcser", x=x, y=y, lagmax=lagmax, dif=dif, level=level)


def yoypx(x: Any, prd: Any, by: Any = None, lag: Any = 1, units: Any = "count", minbase: Any = 0,
          level: Any = 0.95, rnd: Any = None, kind: Any = "nearest", adj: Any = "BH",
          alpha: Any = 0.05) -> Any:
    """Year-over-year change in a published series.

    The arguments and the result are those of R's ``yoypx()``: see its help
    page.

    Examples
    --------
    >>> r = yoypx(x=[10.5, 12, 9, 15], prd=[1, 2, 3, 4], units="percent", minbase=10)
    >>> list(r)
    ['prd', 'value', 'prev', 'change', 'ppchg', 'flag']
    """
    return jcall("yoypx", x=x, prd=prd, by=by, lag=lag, units=units, minbase=minbase, level=level,
                  rnd=rnd, kind=kind, adj=adj, alpha=alpha)


def clprep(data: Any = None, preset: Any = None) -> Any:
    """Presets: a published dataset's columns onto a generic function's arguments.

    The arguments are those of R's ``clprep()``: see its help page. ``data``
    is a dict of columns; ``preset`` a name, or a preset of one's own as a
    dict of columns laid out as the preset table. The result holds
    ``args``, the generic's arguments (a dict to unpack into the generic
    with ``**``), and the preset's ``preset``, ``generic``, ``source``,
    ``units``, ``note`` and ``drop`` (the rows dropped), the attributes of
    the R result. Without arguments, the list of presets.

    Examples
    --------
    >>> d = {"EndFiscalYear": [2020, 2020, 2021, 2021], "Gender": ["Female", "Male"] * 2,
    ...      "NumberIndividuals_RestrictiveConfinement": [40, 310, 52, 335]}
    >>> p = clprep(d, "otis-c01-cnt")
    >>> list(p["args"]), p["args"]["grp"], p["drop"]
    (['y', 'grp', 'adj'], [1.0, 0.0, 1.0, 0.0], 0.0)
    """
    return jcall("clprep", data=data, preset=preset)


def pedlr(ped: Any, geno: Any, freq: Any, alt: Any = None, mut: Any = "none", rate: Any = 0.001,
          range: Any = 0.5) -> Any:
    """Pedigree likelihood ratio.

    The arguments and the result are those of R's ``pedlr()``: see its help
    page. ``ped``, ``alt``, ``geno`` and ``freq`` are dicts of columns.

    Examples
    --------
    >>> trio = {"id": ["AF", "M", "C"], "fid": [None, None, "AF"], "mid": [None, None, "M"]}
    >>> alt = {"id": ["AF", "M", "C", "X"], "fid": [None, None, "X", None],
    ...        "mid": [None, None, "M", None]}
    >>> geno = {"id": ["AF", "M", "C"], "locus": ["D3"] * 3, "a1": [14, 15, 15], "a2": [16, 17, 16]}
    >>> freq = {"locus": ["D3"] * 4, "allele": [14, 15, 16, 17], "freq": [0.1, 0.25, 0.25, 0.2]}
    >>> pedlr(trio, geno, freq, alt)["lr"]
    2.0
    """
    return jcall("pedlr", ped=ped, geno=geno, freq=freq, alt=alt, mut=mut, rate=rate, range=range)


def gmixlr(peaks: Any, freq: Any, ref: Any = None, hp: Any = None, nc: Any = None, hd: Any = None,
           at: Any = 50, pc: Any = 0.05, lambda_: Any = 0.01, theta: Any = 0, stut: Any = False,
           deg: Any = False, par: Any = None) -> Any:
    """Continuous DNA mixture likelihood ratio (gamma model).

    The arguments and the result are those of R's ``gmixlr()``: see its help
    page (``lambda`` is ``lambda_`` here). ``peaks``, ``freq`` and ``ref``
    are dicts of columns; ``par`` a dict of named parameters.

    Examples
    --------
    >>> peaks = {"locus": ["L", "L"], "allele": [1, 3], "height": [800, 900]}
    >>> freq = {"locus": ["L"] * 4, "allele": [1, 2, 3, 4], "freq": [0.1, 0.2, 0.3, 0.25]}
    >>> ref = {"id": ["S"], "locus": ["L"], "a1": [1], "a2": [3]}
    >>> r = gmixlr(peaks, freq, ref, hp="S", nc=1, pc=0, par={"mu": 1000, "omega": 0.2, "m1": 1})
    >>> round(r["lr"], 6)
    16.666667
    """
    return jcall("gmixlr", peaks=peaks, freq=freq, ref=ref, hp=hp, nc=nc, hd=hd, at=at, pc=pc,
                 lambda_=lambda_, theta=theta, stut=stut, deg=deg, par=par)


def cmpsc(a: Any, b: Any, len: Any = 50, npeaks: Any = (5, 3, 1), tx: Any = 25,
          full: Any = False) -> Any:
    """Congruent matching profile segments (CMPS) of two bullet signatures.

    The arguments and the result are those of R's ``cmpsc()``: see its help
    page.

    Examples
    --------
    >>> import math
    >>> a = [math.sin(i / 9) + math.sin(i / 4.1) + 0.3 * math.cos(i * 1.7) for i in range(600)]
    >>> r = cmpsc(a, a[30:])
    >>> r["pos"], r["nseg"]
    (-30.0, 12.0)
    """
    return jcall("cmpsc", a=a, b=b, len=len, npeaks=npeaks, tx=tx, full=full)


def thrtst(stops: Any, srch: Any, hits: Any, grp: Any, area: Any = None, off: Any = None, ref: Any = None,
           level: Any = 0.95, maxit: Any = 2000, rank: Any = "pair", fit: Any = "ml", chains: Any = 4,
           iter: Any = 1000, seed: Any = None) -> Any:
    """Tests of discrimination in searches: hit rates, thresholds, rank order.

    The arguments and the result are those of R's ``thrtst()``: see its help
    page.

    Examples
    --------
    >>> r = thrtst(stops=[100, 100, 100, 100], srch=[10, 20, 30, 10], hits=[5, 5, 5, 5],
    ... grp=[1, 1, 2, 2], off=["x", "y", "x", "y"])
    >>> list(r)
    ['hit', 'kpt', 'thr', 'lrt', 'fit', 'cells', 'rank', 'ref', 'mono', 'post']
    """
    return jcall("thrtst", stops=stops, srch=srch, hits=hits, grp=grp, area=area, off=off, ref=ref,
                 level=level, maxit=maxit, rank=rank, fit=fit, chains=chains, iter=iter, seed=seed)


def ibench(y: Any, off: Any, x: Any = None, minn: Any = 20, fdr: Any = 0.05, dr: Any = False,
           local: Any = False) -> Any:
    """Internal benchmarking of officers.

    The arguments and the result are those of R's ``ibench()``: see its help
    page.

    Examples
    --------
    >>> r = ibench(y=[1, 0, 1, 1, 0, 0, 1, 0, 1], off=["a", "a", "a", "a", "b", "b", "b", "b",
    ... "c"], minn=2)
    >>> list(r)
    ['tab', 'fdr', 'minn', 'nx', 'dr', 'null']
    """
    return jcall("ibench", y=y, off=off, x=x, minn=minn, fdr=fdr, dr=dr, local=local)


def lnroc(tp: Any, pick: Any, conf: Any, proc: Any = None, size: Any = None, far: Any = None, nboot: Any = 0,
          seed: Any = None, level: Any = 0.95) -> Any:
    """Eyewitness lineups: ROC curves, partial areas and confidence-accuracy.

    The arguments and the result are those of R's ``lnroc()``: see its help
    page.

    Examples
    --------
    >>> r = lnroc(tp=[1, 1, 1, 1, 0, 0, 0, 0], pick=["suspect", "suspect", "filler", "none",
    ... "suspect", "filler", "none", "none"], conf=[3, 1, 2, None, 2, 3, None, None])
    >>> list(r)
    ['roc', 'cac', 'pauc', 'far', 'size', 'nboot']
    """
    return jcall("lnroc", tp=tp, pick=pick, conf=conf, proc=proc, size=size, far=far, nboot=nboot, seed=seed,
                 level=level)


def lnsize(choice: Any, k: Any = None, susp: Any = 1, level: Any = 0.95) -> Any:
    """Fairness of a lineup from mock witnesses.

    The arguments and the result are those of R's ``lnsize()``: see its help
    page.

    Examples
    --------
    >>> r = lnsize(choice=[1, 1, 1, 2, 3, 3], k=3)
    >>> list(r)
    ['count', 'stat', 'susp']
    """
    return jcall("lnsize", choice=choice, k=k, susp=susp, level=level)


def rptvic(id: Any, t: Any = None, pop: Any = None, brk: Any = None) -> Any:
    """Repeat victimisation: how often the same targets are victimised again.

    The arguments and the result are those of R's ``rptvic()``: see its help
    page.

    Examples
    --------
    >>> r = rptvic(id=[1, 2, 1, 3, 1, 2], t=[1, 5, 4, 9, 20, 6], pop=10, brk=[0, 7, 30])
    >>> list(r)
    ['tab', 'stat', 'gap', 'lag']
    """
    return jcall("rptvic", id=id, t=t, pop=pop, brk=brk)


def wdqx(a: Any, b: Any, c: Any) -> Any:
    """Weighted displacement quotient: did crime move next door?

    The arguments and the result are those of R's ``wdqx()``: see its help
    page.

    Examples
    --------
    >>> r = wdqx(a=[100, 60], b=[50, 40], c=[200, 180])
    >>> list(r)
    ['sm', 'bdm', 'wdq', 'ge', 'ne', 'ea', 'eb', 'tne']
    """
    return jcall("wdqx", a=a, b=b, c=c)


def prohs(x: Any, y: Any = None, t: Any = None, at: Any = None, win: Any = None, cell: Any = None, bw: Any = None,
          tw: Any = 60, tunit: Any = 7, top: Any = 0.1, ahead: Any = 7) -> Any:
    """Prospective hot-spotting: where the next burglaries are likely.

    The arguments and the result are those of R's ``prohs()``: see its help
    page.

    Examples
    --------
    >>> r = prohs(x=[[0.5, 0.5], [0.6, 0.4], [1.6, 0.5]], t=[0, 8, 9], at=7, win=[0, 2, 0, 1],
    ... cell=1,
    ... bw=2, top=0.5)
    >>> list(r)
    ['x', 'y', 'z', 'hot', 'eval', 'par']
    """
    return jcall("prohs", x=x, y=y, t=t, at=at, win=win, cell=cell, bw=bw, tw=tw, tunit=tunit, top=top,
                 ahead=ahead)


def lnkpr(x: Any, y: Any = None, t: Any = None, series: Any = None, mo: Any = None) -> Any:
    """Crime linkage: which crimes belong to one series.

    The arguments and the result are those of R's ``lnkpr()``: see its help
    page.

    Examples
    --------
    >>> r = lnkpr(x=[[20, 45], [20, 60], [49, 41], [7, 41], [-8, 21], [26, 49], [-1, 28], [4, 45],
    ... [25, 46],
    ... [55, 35], [14, 22], [29, 45]], t=[35, 59, 27, 70, 30, 26, 54, 30, 5, 18, 8, -6],
    ... series=[1, 1, 1, 1, None, 2, 2, 2, 3, 3, None, 3])
    >>> list(r)
    ['coef', 'auc', 'pairs', 'dev', 'nmo']
    """
    return jcall("lnkpr", x=x, y=y, t=t, series=series, mo=mo)


def spcure(time: Any, event: Any, x: Any = None, dist: Any = "weibull", level: Any = 0.95,
           maxit: Any = 500) -> Any:
    """Split-population (cure) model of recidivism.

    The arguments and the result are those of R's ``spcure()``: see its help
    page.

    Examples
    --------
    >>> r = spcure(time=[1, 2, 3, 4, 50, 60, 70, 80], event=[1, 1, 1, 1, 0, 0, 0, 0], dist="exp")
    >>> list(r)
    ['coef', 'fit', 'dist', 'n']
    """
    return jcall("spcure", time=time, event=event, x=x, dist=dist, level=level, maxit=maxit)


def apcie(y: Any, expo: Any, level: Any = 0.95, maxit: Any = 100) -> Any:
    """Age-period-cohort model by the intrinsic estimator.

    The arguments and the result are those of R's ``apcie()``: see its help
    page.

    Examples
    --------
    >>> r = apcie(y=[[5, 9, 4], [8, 12, 10], [3, 6, 7]], expo=[[100, 100, 100], [100, 100, 100],
    ... [100, 100, 100]])
    >>> list(r)
    ['coef', 'fit', 'dims']
    """
    return jcall("apcie", y=y, expo=expo, level=level, maxit=maxit)


def hsacc(pred: Any, obs: Any, area: Any = None, top: Any = 0.1, prev: Any = None, ties: Any = "in") -> Any:
    """Accuracy of a hot-spot map.

    The arguments and the result are those of R's ``hsacc()``: see its help
    page.

    Examples
    --------
    >>> r = hsacc(pred=[9, 8, 7, 6, 5, 4, 3, 2, 1, 0], obs=[3, 0, 5, 1, 0, 0, 2, 9, 0, 0], top=[0.2,
    ... 0.3])
    >>> list(r)
    ['top', 'cells', 'area', 'ashare', 'n', 'N', 'hit', 'pai', 'nmax', 'pei', 'rri']
    """
    return jcall("hsacc", pred=pred, obs=obs, area=area, top=top, prev=prev, ties=ties)


def mtrbnd(y: Any, z: Any, ylo: Any = 0, yhi: Any = 1, assume: Any = "both", w: Any = None) -> Any:
    """Bounds on sentence effects under monotone response and selection.

    The arguments and the result are those of R's ``mtrbnd()``: see its help
    page.

    Examples
    --------
    >>> r = mtrbnd(y=[1, 0, 0, 1, 1, 0, 1, 1], z=[1, 1, 1, 2, 2, 2, 3, 3])
    >>> list(r)
    ['bnd', 'assume', 'range']
    """
    return jcall("mtrbnd", y=y, z=z, ylo=ylo, yhi=yhi, assume=assume, w=w)


def ltrlrng(x: Any, W: Any, model: Any = None, m: Any = 1, navsd: Any = 0) -> Any:
    """Lateral range curve of search theory: the chance of detection at each distance.

    The arguments and the result are those of R's ``ltrlrng()``: see its help
    page.

    Examples
    --------
    >>> round(ltrlrng(x=0.5, W=1), 6)
    0.272623
    """
    return jcall("ltrlrng", x=x, W=W, model=model, m=m, navsd=navsd)


def swpwdth(x: Any, p: Any) -> Any:
    """Sweep width: the area under a sampled lateral range curve.

    The arguments and the result are those of R's ``swpwdth()``: see its help
    page.

    Examples
    --------
    >>> swpwdth(x=[-1, 0, 1], p=[0, 1, 0])
    1.0
    """
    return jcall("swpwdth", x=x, p=p)


def srchpod(coverage: Any = None, model: Any = None, m: Any = 1, W: Any = None, effort: Any = None,
            area: Any = None) -> Any:
    """Probability of detection from coverage.

    The arguments and the result are those of R's ``srchpod()``: see its help
    page.

    Examples
    --------
    >>> round(srchpod(coverage=1, model="invcube"), 4)
    0.7899
    """
    return jcall("srchpod", coverage=coverage, model=model, m=m, W=W, effort=effort, area=area)


def pswppod(W: Any, S: Any, model: Any = None, m: Any = 1, navsd: Any = 0, npts: Any = 2001,
            ntrk: Any = 60) -> Any:
    """Probability of detection of a parallel sweep at track spacing S.

    The arguments and the result are those of R's ``pswppod()``: see its help
    page.

    Examples
    --------
    >>> round(pswppod(W=1, S=2, model="definite"), 3)
    0.5
    """
    return jcall("pswppod", W=W, S=S, model=model, m=m, navsd=navsd, npts=npts, ntrk=ntrk)


def srchscs(poc: Any, pod: Any) -> Any:
    """Probability of success of a search and the containment after it failed.

    The arguments and the result are those of R's ``srchscs()``: see its help
    page.

    Examples
    --------
    >>> round(srchscs(poc=[0.5, 0.3, 0.2], pod=[0.8, 0.5, 0])["pos"], 2)
    0.55
    """
    return jcall("srchscs", poc=poc, pod=pod)


def datumpmap(xedg: Any, yedg: Any, datum: Any = None, center: Any = (0, 0), sigma: Any = 1, end: Any = None,
              nodes: Any = 64, wts: Any = None) -> Any:
    """Datum probability map: the chance that each grid cell holds the object.

    The arguments and the result are those of R's ``datumpmap()``: see its help
    page.

    Examples
    --------
    >>> round(datumpmap(xedg=[-1, 0, 1], yedg=[-1, 0, 1])["outside"], 6)
    0.533935
    """
    return jcall("datumpmap", xedg=xedg, yedg=yedg, datum=datum, center=center, sigma=sigma, end=end,
                 nodes=nodes, wts=wts)


def optcrcsrch(sigma: Any, W: Any, effort: Any) -> Any:
    """Koopman's optimal search of a circular normal datum.

    The arguments and the result are those of R's ``optcrcsrch()``: see its
    help page.

    Examples
    --------
    >>> round(optcrcsrch(sigma=1, W=1, effort=3.141592653589793)["radius"] ** 2, 10)
    2.0
    """
    return jcall("optcrcsrch", sigma=sigma, W=W, effort=effort)


def opteffort(poc: Any, area: Any = None, W: Any = None, effort: Any = None,
              tol: Any = 1e-14) -> Any:
    """Optimal allocation of a fixed search effort over cells.

    The arguments and the result are those of R's ``opteffort()``: see its
    help page.

    Examples
    --------
    >>> [round(z, 8) for z in opteffort(poc=[0.5, 0.5], area=[1, 1], W=1, effort=2)["effort"]]
    [1.0, 1.0]
    """
    return jcall("opteffort", poc=poc, area=area, W=W, effort=effort, tol=tol)


def optsqr(effort: Any, model: Any = None, m: Any = 1, E: Any = 1) -> Any:
    """Optimal square searches of a circular normal datum (Frost, chapter 6).

    The arguments and the result are those of R's ``optsqr()``: see its help
    page.

    Examples
    --------
    >>> round(optsqr(effort=4)["factor"][0], 2)
    1.02
    """
    return jcall("optsqr", effort=effort, model=model, m=m, E=E)


def srchfig(fig: Any = None, num: Any = None) -> Any:
    """The figures of Frost's The Theory of Search, with a line of R for each.

    The arguments and the result are those of R's ``srchfig()``: see its help
    page.

    Examples
    --------
    >>> len(srchfig()["fig"])
    36
    """
    return jcall("srchfig", fig=fig, num=num)


def pgsig(x: Any, stab: Any = None, set: Any = None) -> Any:
    """Pasquill-Gifford dispersion coefficients by stability class.

    The arguments and the result are those of R's ``pgsig()``: see its help
    page.

    Examples
    --------
    >>> round(pgsig(x=1000)['sigmay'], 4)
    76.277
    """
    return jcall("pgsig", x=x, stab=stab, set=set)


def briggs(x: Any, u: Any, diam: Any, exitvel: Any, stktmp: Any, ambtmp: Any, stab: Any = None,
           dthdz: Any = None, g: Any = 9.80616) -> Any:
    """Briggs plume rise of a buoyant stack plume.

    The arguments and the result are those of R's ``briggs()``: see its help
    page.

    Examples
    --------
    >>> round(briggs(x=100, u=4, diam=1, exitvel=12, stktmp=420, ambtmp=285)['flux'], 4)
    9.4559
    """
    return jcall("briggs", x=x, u=u, diam=diam, exitvel=exitvel, stktmp=stktmp, ambtmp=ambtmp,
                 stab=stab, dthdz=dthdz, g=g)


def gplume(q: Any, u: Any, h: Any, rcpt: Any, stab: Any = None, set: Any = None, mixh: Any = None,
           nimg: Any = 3) -> Any:
    """Gaussian plume concentration at receptor points.

    The arguments and the result are those of R's ``gplume()``: see its help
    page.

    Examples
    --------
    >>> gplume(q=100, u=4, h=50, rcpt=[-10, 0, 0])
    0.0
    """
    return jcall("gplume", q=q, u=u, h=h, rcpt=rcpt, stab=stab, set=set, mixh=mixh, nimg=nimg)


def gspuff(mass: Any, u: Any, h: Any, rcpt: Any, t: Any, stab: Any = None, set: Any = None,
           mixh: Any = None, nimg: Any = 3) -> Any:
    """Gaussian puff concentration some time after an instantaneous release.

    The arguments and the result are those of R's ``gspuff()``: see its help
    page.

    Examples
    --------
    >>> gspuff(mass=0, u=3, h=10, rcpt=[1800, 0, 0], t=600)
    0.0
    """
    return jcall("gspuff", mass=mass, u=u, h=h, rcpt=rcpt, t=t, stab=stab, set=set, mixh=mixh,
                 nimg=nimg)


def advdif2d(c0: Any, u: Any, v: Any, kx: Any, ky: Any, dx: Any, dy: Any, dt: Any, nsteps: Any,
             source: Any = None, maxwork: Any = 1e8) -> Any:
    """Explicit advection-diffusion of a concentration field on a grid.

    The arguments and the result are those of R's ``advdif2d()``: see its help
    page.

    Examples
    --------
    >>> advdif2d(c0=[[0, 0, 0], [0, 1, 0], [0, 0, 0]], u=0, v=0, kx=0, ky=0, dx=1, dy=1, dt=1,
    ... nsteps=0)['mass']
    1.0
    """
    return jcall("advdif2d", c0=c0, u=u, v=v, kx=kx, ky=ky, dx=dx, dy=dy, dt=dt, nsteps=nsteps,
                 source=source, maxwork=maxwork)


def lagrprt(n: Any, x0: Any, y0: Any, u: Any, v: Any, kx: Any, ky: Any, dt: Any, nsteps: Any,
            seed: Any = 0, grid: Any = None, maxwork: Any = 1e8) -> Any:
    """Lagrangian random walk of released particles.

    The arguments and the result are those of R's ``lagrprt()``: see its help
    page.

    Examples
    --------
    >>> lagrprt(n=2, x0=0, y0=0, u=1, v=0, kx=0, ky=0, dt=1, nsteps=3)['mean_x']
    3.0
    """
    return jcall("lagrprt", n=n, x0=x0, y0=y0, u=u, v=v, kx=kx, ky=ky, dt=dt, nsteps=nsteps,
                 seed=seed, grid=grid, maxwork=maxwork)


def plmfld(q: Any, u: Any, h: Any, x: Any, y: Any, z: Any = 0, stab: Any = None, set: Any = None,
           mixh: Any = None, nimg: Any = 3) -> Any:
    """A Gaussian-plume concentration field on a grid.

    The arguments and the result are those of R's ``plmfld()``: see its help
    page.

    Examples
    --------
    >>> plmfld(q=100, u=5, h=50, x=[100, 200], y=[-50, 0, 50])['units']
    'concentration'
    """
    return jcall("plmfld", q=q, u=u, h=h, x=x, y=y, z=z, stab=stab, set=set, mixh=mixh, nimg=nimg)


def crfpm25(expo: Any, y: Any = None, refconc: Any = 5.8) -> Any:
    """Concentration-response function for PM2.5.

    The arguments and the result are those of R's ``crfpm25()``: see its help
    page.

    Examples
    --------
    >>> crfpm25(expo=5)['rr']
    1.0
    """
    return jcall("crfpm25", expo=expo, y=y, refconc=refconc)


def crfno2(expo: Any, y: Any = None, refconc: Any = 10, betaper10: Any = None) -> Any:
    """Concentration-response function for NO2.

    The arguments and the result are those of R's ``crfno2()``: see its help
    page.

    Examples
    --------
    >>> crfno2(expo=8)['rr']
    1.0
    """
    return jcall("crfno2", expo=expo, y=y, refconc=refconc, betaper10=betaper10)


def attrfrac(rr: Any, expoprev: Any, rrlo: Any = None, rrup: Any = None) -> Any:
    """Population attributable fraction (Levin).

    The arguments and the result are those of R's ``attrfrac()``: see its help
    page.

    Examples
    --------
    >>> attrfrac(rr=2, expoprev=1)
    0.5
    """
    return jcall("attrfrac", rr=rr, expoprev=expoprev, rrlo=rrlo, rrup=rrup)


def mortdisp(expdel: Any, pop: Any, brate: Any, betaperunit: Any, rare: Any = False,
             hirate: Any = False) -> Any:
    """Deaths avoided under a counterfactual exposure reduction.

    The arguments and the result are those of R's ``mortdisp()``: see its help
    page.

    Examples
    --------
    >>> mortdisp(expdel=0, pop=1e6, brate=0.008, betaperunit=0.0039)
    0.0
    """
    return jcall("mortdisp", expdel=expdel, pop=pop, brate=brate, betaperunit=betaperunit,
                 rare=rare, hirate=hirate)


def polbrdn(expomean: Any, expoprev: Any, brate: Any, pop: Any, poll: Any = None, y: Any = None,
            refconc: Any = None, rare: Any = False, hirate: Any = False) -> Any:
    """End-to-end pollution burden of one population.

    The arguments and the result are those of R's ``polbrdn()``: see its help
    page.

    Examples
    --------
    >>> polbrdn(expomean=25, expoprev=1, brate=0.008, pop=1e6)['basecases']
    8000.0
    """
    return jcall("polbrdn", expomean=expomean, expoprev=expoprev, brate=brate, pop=pop, poll=poll,
                 y=y, refconc=refconc, rare=rare, hirate=hirate)


def polbrdarea(areatab: Any, areacol: Any = None, expcol: Any = None, popcol: Any = None,
               basratecol: Any = None, poll: Any = None, y: Any = None, rare: Any = False,
               hirate: Any = False) -> Any:
    """Pollution burden by area, worst first.

    The arguments and the result are those of R's ``polbrdarea()``: see its help
    page.

    Examples
    --------
    >>> polbrdarea(areatab={'fsa': ['A1', 'B2'], 'exposure': [28, 18], 'population': [40000, 60000],
    ... 'brate': [0.008, 0.008]})['fsa']
    ['A1', 'B2']
    """
    return jcall("polbrdarea", areatab=areatab, areacol=areacol, expcol=expcol, popcol=popcol,
                 basratecol=basratecol, poll=poll, y=y, rare=rare, hirate=hirate)


def expconcidx(data: Any, expo: Any, income: Any) -> Any:
    """Concentration index of exposure by income.

    The arguments and the result are those of R's ``expconcidx()``: see its help
    page.

    Examples
    --------
    >>> round(expconcidx(data={'e': [30, 25, 20, 15, 10], 'i': [1, 2, 3, 4, 5]}, expo='e',
    ... income='i')['concidx'], 4)
    -0.2
    """
    return jcall("expconcidx", data=data, expo=expo, income=income)


def onlsmry(x: Any = None) -> Any:
    """Start a moment accumulator, optionally with a first block.

    The arguments and the result are those of R's ``onlsmry()``: see its help
    page.

    Examples
    --------
    >>> onlsmry(x=[1, 2, 3])['state'][:2]
    [3.0, 2.0]
    """
    return jcall("onlsmry", x=x)


def smryupd(acc: Any, x: Any) -> Any:
    """Add a block to a moment accumulator.

    The arguments and the result are those of R's ``smryupd()``: see its help
    page.

    Examples
    --------
    >>> smryupd(acc=onlsmry(), x=[4, 6])['state'][1]
    5.0
    """
    return jcall("smryupd", acc=acc, x=x)


def smrymrg(a: Any, b: Any) -> Any:
    """Merge two moment accumulators exactly.

    The arguments and the result are those of R's ``smrymrg()``: see its help
    page.

    Examples
    --------
    >>> smrymrg(a=onlsmry(x=[1]), b=onlsmry(x=[3]))['state'][1]
    2.0
    """
    return jcall("smrymrg", a=a, b=b)


def smrystat(acc: Any) -> Any:
    """The statistics of a moment accumulator.

    The arguments and the result are those of R's ``smrystat()``: see its help
    page.

    Examples
    --------
    >>> smrystat(acc=onlsmry(x=[1, 2, 3]))['variance']
    1.0
    """
    return jcall("smrystat", acc=acc)


def rsvidx(n: Any, k: Any, seed: Any = 42) -> Any:
    """Positions of a uniform reservoir sample of a stream.

    The arguments and the result are those of R's ``rsvidx()``: see its help
    page.

    Examples
    --------
    >>> rsvidx(n=3, k=10)
    [1.0, 2.0, 3.0]
    """
    return jcall("rsvidx", n=n, k=k, seed=seed)


def rsvsmp(x: Any, k: Any, seed: Any = 42) -> Any:
    """A uniform reservoir sample of a vector.

    The arguments and the result are those of R's ``rsvsmp()``: see its help
    page.

    Examples
    --------
    >>> rsvsmp(x=['a', 'b'], k=5)
    ['a', 'b']
    """
    return jcall("rsvsmp", x=x, k=k, seed=seed)


def dstnskch(x: Any, p: Any = 14, regstr: Any = None) -> Any:
    """HyperLogLog registers of a vector.

    The arguments and the result are those of R's ``dstnskch()``: see its help
    page.

    Examples
    --------
    >>> len(dstnskch(x=['a', 'b'], p=4))
    16
    """
    return jcall("dstnskch", x=x, p=p, regstr=regstr)


def dstncnt(regstr: Any) -> Any:
    """Estimated distinct values from HyperLogLog registers.

    The arguments and the result are those of R's ``dstncnt()``: see its help
    page.

    Examples
    --------
    >>> dstncnt(regstr=[0] * 16)
    0.0
    """
    return jcall("dstncnt", regstr=regstr)


def skchmrg(a: Any, b: Any) -> Any:
    """Merge two HyperLogLog sketches.

    The arguments and the result are those of R's ``skchmrg()``: see its help
    page.

    Examples
    --------
    >>> skchmrg(a=[1, 0], b=[0, 2])
    [1.0, 2.0]
    """
    return jcall("skchmrg", a=a, b=b)


def sctidy(x: Any) -> Any:
    """Read a Statistics Canada table as published.

    The arguments and the result are those of R's ``sctidy()``: see its help
    page.

    Examples
    --------
    >>> r = sctidy(x={"GEO": ["Canada", "Canada"], "VALUE": ["45.8", "3"], "STATUS": ["", "x"]})
    >>> list(r)
    ['GEO', 'VALUE', 'STATUS', 'quality']
    """
    return jcall("sctidy", x=x)


def qltyrcd(x: Any, suppr: Any = None, caution: Any = None) -> Any:
    """Recode another agency's quality flags to the package's.

    The arguments and the result are those of R's ``qltyrcd()``: see its help
    page.

    Examples
    --------
    >>> qltyrcd(x=["", "np", "*"], suppr="np", caution="*")
    ['ok', 'x', 'E']
    """
    return jcall("qltyrcd", x=x, suppr=suppr, caution=caution)


def limfrom(est: Any, se: Any = None, cv: Any = None, moe: Any = None, moelev: Any = 0.9,
            conflvl: Any = 0.95) -> Any:
    """Confidence limits from a standard error, a coefficient of variation or a margin of error.

    The arguments and the result are those of R's ``limfrom()``: see its help
    page.

    Examples
    --------
    >>> r = limfrom(est=[24.1, 31.5], cv=[8, 12])
    >>> list(r)
    ['estimate', 'lower', 'upper', 'se']
    """
    return jcall("limfrom", est=est, se=se, cv=cv, moe=moe, moelev=moelev, conflvl=conflvl)


def indgap(data: Any, grp: Any, est: Any, lower: Any, upper: Any, ref: Any = "Non-Indigenous",
           qual: Any = None, conflvl: Any = 0.95, pooled: Any = None) -> Any:
    """Compare published group estimates with a reference group.

    The arguments and the result are those of R's ``indgap()``: see its help
    page.

    Examples
    --------
    >>> d = {"g": ["First Nations", "Non-Indigenous"], "e": [45.8, 61.5], "l": [43.9, 61], "u": [47.7, 62]}
    >>> r = indgap(data=d, grp="g", est="e", lower="l", upper="u")
    >>> list(r)[:4]
    ['group', 'estimate', 'refest', 'ratio']
    """
    return jcall("indgap", data=data, grp=grp, est=est, lower=lower, upper=upper, ref=ref, qual=qual,
                 conflvl=conflvl, pooled=pooled)


def genexpo(parent: Any, gparent: Any, partial: Any = "na") -> Any:
    """Generations of residential school attendance in a family.

    The arguments and the result are those of R's ``genexpo()``: see its help
    page.

    Examples
    --------
    >>> genexpo(parent=[True, False, True], gparent=[True, False, False])
    [2.0, 0.0, 1.0]
    """
    return jcall("genexpo", parent=parent, gparent=gparent, partial=partial)


def exptrnd(cases: Any, n: Any, dose: Any = None, conflvl: Any = 0.95) -> Any:
    """A dose-response across generations of exposure.

    The arguments and the result are those of R's ``exptrnd()``: see its help
    page.

    Examples
    --------
    >>> r = exptrnd(cases=[60, 95, 70], n=[500, 520, 300])
    >>> list(r)
    ['dose', 'cases', 'n', 'proportion', 'orratio', 'lower', 'upper']
    """
    return jcall("exptrnd", cases=cases, n=n, dose=dose, conflvl=conflvl)


def trccmpl() -> Any:
    """Calls to Action completed, year by year.

    The result is that of R's ``trccmpl()``: see its help page.

    Examples
    --------
    >>> list(trccmpl())
    ['period', 'yearend', 'completed', 'cumulative', 'source']
    """
    return jcall("trccmpl")


def trcprog(done: Any = None, yearend: Any = None, total: Any = 94, start: Any = 2015.96,
            conflvl: Any = 0.95) -> Any:
    """How fast the Calls to Action are being completed.

    The arguments and the result are those of R's ``trcprog()``: see its help
    page.

    Examples
    --------
    >>> list(trcprog())[:3]
    ['completed', 'remaining', 'years']
    """
    return jcall("trcprog", done=done, yearend=yearend, total=total, start=start, conflvl=conflvl)


def mhcont(data: Any, grp: Any, ind: Any, est: Any, lower: Any, upper: Any, hope: Any = None,
           belong: Any = None, meaning: Any = None, purp: Any = None, qual: Any = None,
           conflvl: Any = 0.95) -> Any:
    """The four directions of the First Nations Mental Wellness Continuum, from published indicators.

    The arguments and the result are those of R's ``mhcont()``: see its help
    page.

    Examples
    --------
    >>> d = {"c": ["A", "A"], "i": ["Hopeful", "Meaning"], "p": [64, 70], "l": [59, 65], "u": [69, 75]}
    >>> r = mhcont(data=d, grp="c", ind="i", est="p", lower="l", upper="u", hope="Hopeful")
    >>> r["indicator"]
    ['hope', 'belonging', 'meaning', 'purpose']
    """
    return jcall("mhcont", data=data, grp=grp, ind=ind, est=est, lower=lower, upper=upper, hope=hope,
                 belong=belong, meaning=meaning, purp=purp, qual=qual, conflvl=conflvl)


def cwbscr(incpercap: Any, highsch20p: Any, uni25p: Any, notcro: Any, nomajrep: Any, partic: Any,
           empl: Any, bounds: Any = None) -> Any:
    """The Community Well-Being index of a community.

    The arguments and the result are those of R's ``cwbscr()``: see its help
    page.

    Examples
    --------
    >>> r = cwbscr(incpercap=18000, highsch20p=58, uni25p=6, notcro=71, nomajrep=62, partic=55, empl=78)
    >>> list(r)
    ['income', 'education', 'housing', 'labour', 'cwb']
    """
    return jcall("cwbscr", incpercap=incpercap, highsch20p=highsch20p, uni25p=uni25p, notcro=notcro,
                 nomajrep=nomajrep, partic=partic, empl=empl, bounds=bounds)


def connassoc(data: Any, place: Any, x: Any, y: Any, xlo: Any = None, xup: Any = None, ylo: Any = None,
              yup: Any = None, reps: Any = 10000, seed: Any = 1, conflvl: Any = 0.95) -> Any:
    """Connection across places: does one indicator rise where another does?.

    The arguments and the result are those of R's ``connassoc()``: see its
    help page.

    Examples
    --------
    >>> d = {"p": ["N", "S", "E", "W"], "a": [41, 18, 27, 33], "b": [68, 55, 60, 66]}
    >>> r = connassoc(data=d, place="p", x="a", y="b")
    >>> round(r["rho"], 6)
    1.0
    """
    return jcall("connassoc", data=data, place=place, x=x, y=y, xlo=xlo, xup=xup, ylo=ylo, yup=yup, reps=reps,
                 seed=seed, conflvl=conflvl)


def healchg(data: Any, grp: Any, ind: Any, prd: Any, est: Any, lower: Any, upper: Any, from_: Any,
            to: Any, well: Any = None, qual: Any = None, conflvl: Any = 0.95, distress: Any = None) -> Any:
    """Healing as a process: change in published indicators between two periods.

    The arguments and the result are those of R's ``healchg()``: see its help
    page.

    Examples
    --------
    >>> d = {"c": ["A", "A"], "i": ["Good", "Good"], "s": ["2018", "2023"], "p": [55, 61], "l": [50, 56],
    ...      "u": [60, 66]}
    >>> r = healchg(data=d, grp="c", ind="i", prd="s", est="p", lower="l", upper="u", from_="2018", to="2023")
    >>> r["change"]
    [6.0]
    """
    return jcall("healchg", data=data, grp=grp, ind=ind, prd=prd, est=est, lower=lower, upper=upper,
                 from_=from_, to=to, well=well, qual=qual, conflvl=conflvl, distress=distress)


def commidx(data: Any, grp: Any, ind: Any, est: Any, lower: Any, upper: Any, w: Any = None,
            qual: Any = None, conflvl: Any = 0.95) -> Any:
    """A community-weighted index of protective and risk factors.

    The arguments and the result are those of R's ``commidx()``: see its help
    page.

    Examples
    --------
    >>> d = {"c": ["A", "A"], "i": ["Language", "Food insecure"], "p": [38, 21], "l": [33, 17],
    ...      "u": [43, 25]}
    >>> r = commidx(data=d, grp="c", ind="i", est="p", lower="l", upper="u",
    ...             w={"Language": 2, "Food insecure": -1})
    >>> list(r)
    ['group', 'index', 'lower', 'upper', 'usewithcau']
    """
    return jcall("commidx", data=data, grp=grp, ind=ind, est=est, lower=lower, upper=upper, w=w, qual=qual,
                 conflvl=conflvl)


def frmaud(x: Any) -> Any:
    """How deficit-framed is a set of results?.

    The arguments and the result are those of R's ``frmaud()``: see its help
    page.

    Examples
    --------
    >>> r = frmaud(x={"group": ["A", "A", "A"], "side": ["wellness", "distress", "distress"]})
    >>> list(r)
    ['bygrp', 'shrdef', 'defonly']
    """
    return jcall("frmaud", x=x)


def scests(x: Any) -> Any:
    """The estimates of a Statistics Canada table, with their confidence limits.

    The arguments and the result are those of R's ``scests()``: see its help
    page.

    Examples
    --------
    >>> r = scests(x={"GEO": ["Canada", "Canada"], "Characteristics": ["Rate", "Coefficient of variation for rate"],
    ...               "VALUE": [12.5, 8], "STATUS": ["", ""]})
    >>> r["cisou"]
    ['cv']
    """
    return jcall("scests", x=x)


def ordcmp(data: Any, grp: Any, cat: Any, est: Any, lower: Any, upper: Any, lvls: Any, ref: Any,
           qual: Any = None, conflvl: Any = 0.95, pooled: Any = None) -> Any:
    """Compare rated scales (perceived mental health, life satisfaction) between groups.

    The arguments and the result are those of R's ``ordcmp()``: see its help
    page.

    Examples
    --------
    >>> d = {"g": ["A", "A", "B", "B"], "r": ["Low", "High", "Low", "High"], "p": [40, 60, 30, 70],
    ...      "l": [35, 55, 27, 67], "u": [45, 65, 33, 73]}
    >>> r = ordcmp(data=d, grp="g", cat="r", est="p", lower="l", upper="u", lvls=["Low", "High"], ref="B")
    >>> r["group"]
    ['A']
    """
    return jcall("ordcmp", data=data, grp=grp, cat=cat, est=est, lower=lower, upper=upper, lvls=lvls, ref=ref,
                 qual=qual, conflvl=conflvl, pooled=pooled)


def stdprev(est: Any, lower: Any, upper: Any, std: Any, grp: Any = None, qual: Any = None,
            conflvl: Any = 0.95, scale: Any = 100) -> Any:
    """Age-standardise survey percentages.

    The arguments and the result are those of R's ``stdprev()``: see its help
    page.

    Examples
    --------
    >>> r = stdprev(est=[12, 9, 6], lower=[9, 6, 3], upper=[15, 12, 9], std=[500, 300, 200])
    >>> list(r)
    ['group', 'estimate', 'lower', 'upper', 'se', 'usewithcau']
    """
    return jcall("stdprev", est=est, lower=lower, upper=upper, std=std, grp=grp, qual=qual, conflvl=conflvl,
                 scale=scale)


def welprof(data: Any, grp: Any, ind: Any, est: Any, lower: Any, upper: Any, well: Any, distress: Any,
            domain: Any = None, qual: Any = None) -> Any:
    """Wellness beside distress: a strengths-first profile of published indicators.

    The arguments and the result are those of R's ``welprof()``: see its help
    page.

    Examples
    --------
    >>> d = {"c": ["A", "A"], "i": ["Low", "Hopeful"], "p": [12, 64], "l": [9, 59], "u": [15, 69]}
    >>> r = welprof(data=d, grp="c", ind="i", est="p", lower="l", upper="u", well="Hopeful", distress="Low")
    >>> r["side"]
    ['wellness', 'distress']
    """
    return jcall("welprof", data=data, grp=grp, ind=ind, est=est, lower=lower, upper=upper, well=well,
                 distress=distress, domain=domain, qual=qual)


def holbal(data: Any, grp: Any, domain: Any, est: Any, lower: Any, upper: Any, w: Any = None,
           qual: Any = None, reps: Any = 10000, seed: Any = 1, conflvl: Any = 0.95) -> Any:
    """Holistic balance: how evenly wellness is held across the domains of a framework.

    The arguments and the result are those of R's ``holbal()``: see its help
    page.

    Examples
    --------
    >>> d = {"c": ["A", "A"], "d": ["body", "mind"], "p": [62, 58], "l": [57, 53], "u": [67, 63]}
    >>> r = holbal(data=d, grp="c", domain="d", est="p", lower="l", upper="u", reps=200)
    >>> r["weadom"]
    ['mind']
    """
    return jcall("holbal", data=data, grp=grp, domain=domain, est=est, lower=lower, upper=upper, w=w, qual=qual,
                 reps=reps, seed=seed, conflvl=conflvl)


def dualcont(x: Any, well: Any, distress: Any, conflvl: Any = 0.95) -> Any:
    """Wellness and distress as two continua.

    The arguments and the result are those of R's ``dualcont()``: see its
    help page.

    Examples
    --------
    >>> w = {"group": ["A", "A"], "indicator": ["Hopeful", "Low"], "side": ["wellness", "distress"],
    ...      "estimate": [64, 12], "lower": [59, 9], "upper": [69, 15], "usewithcau": [False, False]}
    >>> r = dualcont(x=w, well="Hopeful", distress="Low")
    >>> r["strratio"]
    [5.333333333333333]
    """
    return jcall("dualcont", x=x, well=well, distress=distress, conflvl=conflvl)


def cnsprof(x: Any) -> Any:
    """Read the 2021 Census Indigenous Population Profile.

    The arguments and the result are those of R's ``cnsprof()``: see its help
    page.

    Examples
    --------
    >>> x = {"GEO_NAME": ["A"], "CHARACTERISTIC_ID": [1], "CHARACTERISTIC_NAME": ["Total - Age groups"],
    ...      "NON_IND": [5000], "SYMBOL_NON_IND": [""]}
    >>> cnsprof(x)["group"]
    ['Non-Indigenous']
    """
    return jcall("cnsprof", x=x)


def profshr(x: Any, num: Any, den: Any = None, smpfrac: Any = 0.25, deff: Any = 1, conflvl: Any = 0.95,
            roubase: Any = 5) -> Any:
    """The share of a population with a characteristic, from the census profile.

    The arguments and the result are those of R's ``profshr()``: see its help
    page.

    Examples
    --------
    >>> x = cnsprof({"GEO_NAME": ["A", "A"], "CHARACTERISTIC_ID": [1, 2],
    ...              "CHARACTERISTIC_NAME": ["Total - Age groups", "  0 to 14 years"],
    ...              "NON_IND": [5000, 900], "SYMBOL_NON_IND": ["", ""]})
    >>> profshr(x, num=2)["estimate"]
    [18.0]
    """
    return jcall("profshr", x=x, num=num, den=den, smpfrac=smpfrac, deff=deff, conflvl=conflvl, roubase=roubase)


def cnsstd(x: Any, grp: Any = "Indigenous identity", geo: Any = "Canada", ages: Any = None,
           resid: Any = "Total - Residence on or off reserve") -> Any:
    """An Indigenous standard population from the census.

    The arguments and the result are those of R's ``cnsstd()``: see its help
    page.

    Examples
    --------
    >>> x = cnsprof({"GEO_NAME": ["A", "A"], "AGE_GROUP": ["0 to 14 years", "15 years and over"],
    ...              "CHARACTERISTIC_ID": [1, 1], "CHARACTERISTIC_NAME": ["Total - Age groups"] * 2,
    ...              "NON_IND": [900, 4100], "SYMBOL_NON_IND": ["", ""]})
    >>> cnsstd(x, grp="Non-Indigenous", geo="A", ages=["0 to 14 years", "15 years and over"])["weight"]
    [0.18, 0.82]
    """
    return jcall("cnsstd", x=x, grp=grp, geo=geo, ages=ages, resid=resid)


def cnscwb(x: Any, dwell: Any = None, groups: Any = None, bounds: Any = None,
           resid: Any = "Total - Residence on or off reserve", minpop: Any = 250) -> Any:
    """A Community Well-Being score from the census profile.

    The arguments and the result are those of R's ``cnscwb()``: see its help
    page.

    Examples
    --------
    >>> x = cnsprof({"GEO_NAME": ["A"], "CHARACTERISTIC_ID": [1], "CHARACTERISTIC_NAME": ["Total - Age groups"],
    ...              "NON_IND": [5000], "SYMBOL_NON_IND": [""]})
    >>> list(cnscwb(x, groups="Non-Indigenous"))[:4]
    ['geo', 'dguid', 'group', 'income']
    """
    return jcall("cnscwb", x=x, dwell=dwell, groups=groups, bounds=bounds, resid=resid, minpop=minpop)


def trcind() -> Any:
    """The indicators of Calls to Action 19 and 55, and the tables that measure them.

    The result is that of R's ``trcind()``: see its help page.

    Examples
    --------
    >>> list(trcind())
    ['call', 'indicator', 'table', 'measure', 'better', 'note']
    """
    return jcall("trcind")


def trchlthgap(health: Any = None, lifeexp: Any = None, infmort: Any = None, suiide: Any = None,
               ref: Any = "Non-Indigenous", agetot: Any = None) -> Any:
    """The Call to Action 19 health gaps, in one place.

    The arguments and the result are those of R's ``trchlthgap()``: see its
    help page.

    Examples
    --------
    >>> le = {"period": ["2017/2019"] * 2, "group": ["Inuit", "Non-Indigenous"], "estimate": [70.9, 81.8],
    ...       "lower": [69.1, 81.6], "upper": [72.7, 82.0], "quality": ["ok", "ok"]}
    >>> trchlthgap(lifeexp=le)["group"]
    ['Inuit']
    """
    return jcall("trchlthgap", health=health, lifeexp=lifeexp, infmort=infmort, suiide=suiide, ref=ref,
                 agetot=agetot)


def polex(seed: Any = 20240101) -> Any:
    """A simulated policy panel for the examples.

    The arguments and the result are those of R's ``polex()``: see its help
    page. The panel is simulated.

    Examples
    --------
    >>> p = polex()
    >>> len(p["district"])
    720
    """
    return jcall("polex", seed=seed)


def vgemp(x: Any, y: Any = None, v: Any = None, nbins: Any = 15, cutoff: Any = None) -> Any:
    """Empirical semivariogram.

    The arguments and the result are those of R's ``vgemp()``: see its help
    page.

    Examples
    --------
    >>> r = vgemp(x=[[0, 0], [1, 0], [0, 1], [1, 1], [2, 2]], v=[1, 2, 3, 5, 4])
    >>> list(r)
    ['dist', 'gamma', 'np']
    """
    return jcall("vgemp", x=x, y=y, v=v, nbins=nbins, cutoff=cutoff)


def vgfit(x: Any, y: Any = None, v: Any = None, model: Any = "exponential") -> Any:
    """Fit a variogram model by Gaussian maximum likelihood.

    The arguments and the result are those of R's ``vgfit()``: see its help
    page.

    Examples
    --------
    >>> r = vgfit(x=[[0, 0], [1, 0], [0, 1], [1, 1], [2, 2]], v=[1, 2, 3, 5, 4])
    >>> list(r)
    ['model', 'nugget', 'psill', 'range', 'loglik', 'conv', 'nfn']
    """
    return jcall("vgfit", x=x, y=y, v=v, model=model)


def okrge(x: Any, y: Any = None, v: Any = None, at: Any = None, vgm: Any = None) -> Any:
    """Ordinary kriging.

    The arguments and the result are those of R's ``okrge()``: see its help
    page.

    Examples
    --------
    >>> r = okrge(x=[[0, 0], [1, 0], [0, 1], [1, 1]], v=[1, 2, 3, 5], at=[[0.5, 0.5]],
    ... vgm={"model": "exponential", "nugget": 0, "psill": 1, "range": 1})
    >>> list(r)
    ['pred', 'var']
    """
    return jcall("okrge", x=x, y=y, v=v, at=at, vgm=vgm)


def bkrig(x: Any, y: Any = None, v: Any = None, blocks: Any = None, nq: Any = 25, nugget: Any = 0,
          sill: Any = 1, range: Any = 1) -> Any:
    """Ordinary block kriging.

    The arguments and the result are those of R's ``bkrig()``: see its help
    page.

    Examples
    --------
    >>> r = bkrig(x=[[0, 0], [1, 0], [0, 1], [1, 1]], v=[1, 2, 3, 5], blocks=[[[0, 0], [1, 1]]])
    >>> list(r)
    ['est', 'se']
    """
    return jcall("bkrig", x=x, y=y, v=v, blocks=blocks, nq=nq, nugget=nugget, sill=sill, range=range)


def krigcv(x: Any, y: Any = None, v: Any = None, nugget: Any = 0, sill: Any = 1, range: Any = 1) -> Any:
    """Leave-one-out cross-validation of ordinary kriging.

    The arguments and the result are those of R's ``krigcv()``: see its help
    page.

    Examples
    --------
    >>> r = krigcv(x=[[0, 0], [1, 0], [0, 1], [1, 1]], v=[1, 2, 3, 5])
    >>> list(r)
    ['mspe', 'rmspe', 'mae', 'resid']
    """
    return jcall("krigcv", x=x, y=y, v=v, nugget=nugget, sill=sill, range=range)


def trsurf(x: Any, y: Any = None, v: Any = None, order: Any = 2) -> Any:
    """Polynomial trend surface.

    The arguments and the result are those of R's ``trsurf()``: see its help
    page.

    Examples
    --------
    >>> r = trsurf(x=[[0, 0], [1, 0], [0, 1], [1, 1]], v=[1, 3, 0, 2], order=1)
    >>> list(r)
    ['coef', 'se', 'r2', 'order', 'n']
    """
    return jcall("trsurf", x=x, y=y, v=v, order=order)


def sglmx(v: Any, z: Any = None, x: Any = None, y: Any = None) -> Any:
    """Spatial linear model by profile maximum likelihood.

    The arguments and the result are those of R's ``sglmx()``: see its help
    page.

    Examples
    --------
    >>> r = sglmx(v=[1, 2, 3, 5, 4], x=[[0, 0], [1, 0], [0, 1], [1, 1], [2, 2]])
    >>> list(r)
    ['coef', 'se', 'sigma2', 'phi', 'tau2', 'n']
    """
    return jcall("sglmx", v=v, z=z, x=x, y=y)


def smixm(v: Any, z: Any = None, x: Any = None, y: Any = None) -> Any:
    """Spatial linear mixed model by REML.

    The arguments and the result are those of R's ``smixm()``: see its help
    page.

    Examples
    --------
    >>> r = smixm(v=[1, 2, 3, 5, 4], x=[[0, 0], [1, 0], [0, 1], [1, 1], [2, 2]])
    >>> list(r)
    ['coef', 'se', 'sigma2', 'tau2', 'phi', 'conv', 'n']
    """
    return jcall("smixm", v=v, z=z, x=x, y=y)


def spobj(data: Any, coords: Any = None, crs: Any = 4326) -> Any:
    """A spatial data object.

    The arguments and the result are those of R's ``spobj()``: see its help
    page.

    Examples
    --------
    >>> r = spobj(data={"lon": [-79.4, -79.3], "lat": [43.6, 43.7], "v": [1, 2]})
    >>> list(r)
    ['coords', 'data', 'cnames', 'crs', 'units', 'n']
    """
    return jcall("spobj", data=data, coords=coords, crs=crs)


def spxy(x: Any) -> Any:
    """The coordinates of a spatial data object.

    The arguments and the result are those of R's ``spxy()``: see its help
    page.

    Examples
    --------
    >>> r = spxy(x=spobj(data={"lon": [1, 2], "lat": [3, 4]}))
    >>> isinstance(r, list)
    True
    """
    return jcall("spxy", x=x)


def crsof(x: Any) -> Any:
    """The coordinate reference system of a spatial data object.

    The arguments and the result are those of R's ``crsof()``: see its help
    page.

    Examples
    --------
    >>> r = crsof(x=spobj(data={"lon": [1, 2], "lat": [3, 4]}))
    >>> list(r)
    ['epsg', 'wkt']
    """
    return jcall("crsof", x=x)


def lfspec(x: Any) -> Any:
    """The specification of a web map.

    The arguments and the result are those of R's ``lfspec()``: see its help
    page.

    Examples
    --------
    >>> r = lfspec(x=spobj(data={"lon": [1, 2], "lat": [3, 4]}))
    >>> list(r)
    ['provider', 'interactive', 'center', 'zoom', 'points']
    """
    return jcall("lfspec", x=x)


def spsamp(x: Any, size: Any = None, dens: Any = True, seed: Any = None) -> Any:
    """Sample points, by local density or uniformly.

    The arguments and the result are those of R's ``spsamp()``: see its help
    page.

    Examples
    --------
    >>> r = spsamp(x=spobj(data={"lon": [1, 2, 3], "lat": [3, 4, 6]}), size=2, seed=1)
    >>> isinstance(r, list)
    True
    """
    return jcall("spsamp", x=x, size=size, dens=dens, seed=seed)


def bayout(x: Any, freq: Any = 52, b: Any = 0, w: Any = 6, acty: Any = True, alpha: Any = 0.05, at: Any = None) -> Any:
    """Bayes outbreak detection.

    The arguments and the result are those of R's ``bayout()``: see its help
    page.

    Examples
    --------
    >>> r = bayout(x=[3, 5, 2, 4, 6, 3, 4, 12], at=8)
    >>> list(r)
    ['t', 'ub', 'alarm']
    """
    return jcall("bayout", x=x, freq=freq, b=b, w=w, acty=acty, alpha=alpha, at=at)


def bufexp(x: Any, y: Any = None, src: Any = None, r: Any = None, w: Any = None, lonlat: Any = False) -> Any:
    """Exposure within a buffer.

    The arguments and the result are those of R's ``bufexp()``: see its help
    page.

    Examples
    --------
    >>> r = bufexp(x=[[0, 0], [5, 5]], src=[[0.5, 0.5], [4, 4]], r=1.5)
    >>> list(r)
    ['count', 'wsum', 'near']
    """
    return jcall("bufexp", x=x, y=y, src=src, r=r, w=w, lonlat=lonlat)


def kerexp(x: Any, y: Any = None, src: Any = None, h: Any = None, w: Any = None, kern: Any = "gaussian") -> Any:
    """Exposure by kernel.

    The arguments and the result are those of R's ``kerexp()``: see its help
    page.

    Examples
    --------
    >>> r = kerexp(x=[[0.2, 0.3], [1, 1]], src=[[0.25, 0.35]], h=0.5)
    >>> isinstance(r, list)
    True
    """
    return jcall("kerexp", x=x, y=y, src=src, h=h, w=w, kern=kern)


def lerprc(w: Any, rho: Any, tau: Any = 1) -> Any:
    """The Leroux precision matrix.

    The arguments and the result are those of R's ``lerprc()``: see its help
    page.

    Examples
    --------
    >>> r = lerprc(w=[[0, 1, 0], [1, 0, 1], [0, 1, 0]], rho=0.5)
    >>> isinstance(r, list)
    True
    """
    return jcall("lerprc", w=w, rho=rho, tau=tau)


def bym2sc(w: Any, tau: Any = 1, phi: Any = 0.5, tol: Any = 1e-10) -> Any:
    """The BYM2 area structure.

    The arguments and the result are those of R's ``bym2sc()``: see its help
    page.

    Examples
    --------
    >>> r = bym2sc(w=[[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    >>> list(r)
    ['scale', 'ginv', 'cov']
    """
    return jcall("bym2sc", w=w, tau=tau, phi=phi, tol=tol)


def ecolrg(y: Any, e: Any, x: Any = None, family: Any = "poisson", z: Any = None, tol: Any = 1e-12,
           maxit: Any = 500) -> Any:
    """Ecological regression with expected counts.

    The arguments and the result are those of R's ``ecolrg()``: see its help
    page.

    Examples
    --------
    >>> r = ecolrg(y=[2, 0, 5, 3, 0, 7], e=[2.5, 1, 4, 3, 1.5, 5])
    >>> list(r)
    ['coef', 'fitted', 'loglik']
    """
    return jcall("ecolrg", y=y, e=e, x=x, family=family, z=z, tol=tol, maxit=maxit)


def spbias(d: Any, b: Any, ub: Any = -10) -> Any:
    """The spatial encoding of graph attention.

    The arguments and the result are those of R's ``spbias()``: see its help
    page.

    Examples
    --------
    >>> r = spbias(d=[[0, 1], [1, -1]], b=[1, 0.5])
    >>> list(r)
    ['bias', 'ub']
    """
    return jcall("spbias", d=d, b=b, ub=ub)


def spsir(i0: Any, w: Any = None, pop: Any = None, beta: Any = 0.5, gamma: Any = 0.2, kappa: Any = 0.1,
          tmax: Any = 100, dt: Any = 0.1) -> Any:
    """The metapopulation SIR model.

    The arguments and the result are those of R's ``spsir()``: see its help
    page.

    Examples
    --------
    >>> r = spsir(i0=[10, 0], w=[[0, 1], [1, 0]], tmax=2, dt=1)
    >>> list(r)
    ['attack', 'r0', 't', 'S', 'I', 'R']
    """
    return jcall("spsir", i0=i0, w=w, pop=pop, beta=beta, gamma=gamma, kappa=kappa, tmax=tmax, dt=dt)


def sfrail(time: Any, event: Any, x: Any = None, reg: Any = None, w: Any = None, tau: Any = 1,
           tol: Any = 1e-10) -> Any:
    """Spatial frailty model.

    The arguments and the result are those of R's ``sfrail()``: see its help
    page.

    Examples
    --------
    >>> r = sfrail(time=[1, 2, 3, 4, 5, 6], event=[1, 1, 0, 1, 1, 0], reg=[1, 1, 1, 2, 2, 2],
    ... w=[[0, 1], [1, 0]])
    >>> list(r)
    ['pll', 'shape', 'coef', 'frail', 'se', 'iter', 'tau']
    """
    return jcall("sfrail", time=time, event=event, x=x, reg=reg, w=w, tau=tau, tol=tol)


def sfcure(time: Any, event: Any, x: Any = None, reg: Any = None, w: Any = None, tau: Any = 1,
           tol: Any = 1e-10) -> Any:
    """Spatial promotion-time cure model.

    The arguments and the result are those of R's ``sfcure()``: see its help
    page.

    Examples
    --------
    >>> r = sfcure(time=[1, 2, 3, 4, 5, 6, 8, 9], event=[1, 1, 0, 1, 1, 0, 0, 0],
    ... reg=[1, 1, 1, 1, 2, 2, 2, 2], w=[[0, 1], [1, 0]])
    >>> list(r)
    ['pll', 'shape', 'lscale', 'coef', 'effect', 'se', 'cure', 'iter', 'tau']
    """
    return jcall("sfcure", time=time, event=event, x=x, reg=reg, w=w, tau=tau, tol=tol)


def evfilt(v: Any, w: Any, x: Any = None, method: Any = "moran", cand: Any = "positive", tol: Any = 0.1,
           maxv: Any = None) -> Any:
    """Eigenvector spatial filtering.

    The arguments and the result are those of R's ``evfilt()``: see its help
    page.

    Examples
    --------
    >>> r = evfilt(v=[1, 2, 2.5, 4, 3.5, 3, 1.5, 1], w=[[0, 1, 0, 0, 0, 0, 0, 0],
    ... [1, 0, 1, 0, 0, 0, 0, 0], [0, 1, 0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 1, 0, 0, 0],
    ... [0, 0, 0, 1, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 1],
    ... [0, 0, 0, 0, 0, 0, 1, 0]], method="aic")
    >>> list(r)
    ['sel', 'coef', 'imoran', 'aic', 'r2', 'fr2', 'resid', 'mi']
    """
    return jcall("evfilt", v=v, w=w, x=x, method=method, cand=cand, tol=tol, maxv=maxv)


def cmn(x: Any) -> Any:
    """Mean (C backend).

    The arguments and the result are those of R's ``cmn()``: see its help
    page.

    Examples
    --------
    >>> cmn(x=[1, 2, 3])
    2.0
    """
    return jcall("cmn", x=x)


def cvari(x: Any) -> Any:
    """Sample variance (C backend).

    The arguments and the result are those of R's ``cvari()``: see its help
    page.

    Examples
    --------
    >>> cvari(x=[1, 2, 3])
    1.0
    """
    return jcall("cvari", x=x)


def ccorr(x: Any, y: Any) -> Any:
    """Pearson correlation (C backend).

    The arguments and the result are those of R's ``ccorr()``: see its help
    page.

    Examples
    --------
    >>> ccorr(x=[1, 2, 3], y=[2, 4, 6])
    1.0
    """
    return jcall("ccorr", x=x, y=y)


def csdev(x: Any, ddof: Any = 1) -> Any:
    """Standard deviation (C backend).

    The arguments and the result are those of R's ``csdev()``: see its help
    page.

    Examples
    --------
    >>> csdev(x=[1, 2, 3])
    1.0
    """
    return jcall("csdev", x=x, ddof=ddof)


def cdistm(a: Any, b: Any) -> Any:
    """Euclidean distance (C backend).

    The arguments and the result are those of R's ``cdistm()``: see its help
    page.

    Examples
    --------
    >>> cdistm(a=[0, 0], b=[3, 4])
    5.0
    """
    return jcall("cdistm", a=a, b=b)


def cnormpdf(x: Any, mean: Any = 0, sd: Any = 1) -> Any:
    """Normal density (C backend).

    The arguments and the result are those of R's ``cnormpdf()``: see its help
    page.

    Examples
    --------
    >>> r = cnormpdf(x=[0, 1])
    >>> len(r)
    2
    """
    return jcall("cnormpdf", x=x, mean=mean, sd=sd)


def cnormlpdf(x: Any, mean: Any = 0, sd: Any = 1) -> Any:
    """Normal log-density (C backend).

    The arguments and the result are those of R's ``cnormlpdf()``: see its help
    page.

    Examples
    --------
    >>> r = cnormlpdf(x=[0, 1])
    >>> len(r)
    2
    """
    return jcall("cnormlpdf", x=x, mean=mean, sd=sd)


def cmom(x: Any) -> Any:
    """Mean, variance, skewness and kurtosis (C backend).

    The arguments and the result are those of R's ``cmom()``: see its help
    page.

    Examples
    --------
    >>> r = cmom(x=[2, 4, 4, 5, 9])
    >>> list(r)
    ['mean', 'variance', 'skewness', 'kurtosis']
    """
    return jcall("cmom", x=x)


def cqnt(x: Any, probs: Any = None) -> Any:
    """Type-7 quantiles (C backend).

    The arguments and the result are those of R's ``cqnt()``: see its help
    page.

    Examples
    --------
    >>> r = cqnt(x=[1, 2, 3, 4], probs=[0.25, 0.5])
    >>> list(r)
    ['25%', '50%']
    """
    return jcall("cqnt", x=x, probs=probs)


def cmed(x: Any) -> Any:
    """Median (C backend).

    The arguments and the result are those of R's ``cmed()``: see its help
    page.

    Examples
    --------
    >>> cmed(x=[1, 2, 3, 4])
    2.5
    """
    return jcall("cmed", x=x)


def cmad(x: Any, const: Any = 1.4826) -> Any:
    """Median absolute deviation (C backend).

    The arguments and the result are those of R's ``cmad()``: see its help
    page.

    Examples
    --------
    >>> cmad(x=[1, 2, 3, 4], const=1)
    1.0
    """
    return jcall("cmad", x=x, const=const)


def ciqr(x: Any) -> Any:
    """Interquartile range (C backend).

    The arguments and the result are those of R's ``ciqr()``: see its help
    page.

    Examples
    --------
    >>> ciqr(x=[1, 2, 3, 4])
    1.5
    """
    return jcall("ciqr", x=x)


def ctkyfnc(x: Any, k: Any = 1.5) -> Any:
    """Tukey's outlier fences (C backend).

    The arguments and the result are those of R's ``ctkyfnc()``: see its help
    page.

    Examples
    --------
    >>> r = ctkyfnc(x=[1, 2, 3, 4])
    >>> list(r)
    ['lower', 'upper']
    """
    return jcall("ctkyfnc", x=x, k=k)


def ctrmmean(x: Any, trim: Any = 0.1) -> Any:
    """Trimmed mean (C backend).

    The arguments and the result are those of R's ``ctrmmean()``: see its help
    page.

    Examples
    --------
    >>> ctrmmean(x=[1, 2, 3, 100], trim=0.25)
    2.5
    """
    return jcall("ctrmmean", x=x, trim=trim)


def cwinsmean(x: Any, trim: Any = 0.1) -> Any:
    """Winsorized mean (C backend).

    The arguments and the result are those of R's ``cwinsmean()``: see its help
    page.

    Examples
    --------
    >>> cwinsmean(x=[1, 2, 3, 100], trim=0.25)
    2.5
    """
    return jcall("cwinsmean", x=x, trim=trim)


def cwtd(x: Any, w: Any) -> Any:
    """Weighted mean and variance (C backend).

    The arguments and the result are those of R's ``cwtd()``: see its help
    page.

    Examples
    --------
    >>> r = cwtd(x=[10, 20, 30], w=[1, 1, 2])
    >>> list(r)
    ['mean', 'variance']
    """
    return jcall("cwtd", x=x, w=w)


def cspr(x: Any, y: Any) -> Any:
    """Spearman's rank correlation (C backend).

    The arguments and the result are those of R's ``cspr()``: see its help
    page.

    Examples
    --------
    >>> cspr(x=[1, 2, 3], y=[1, 4, 9])
    1.0
    """
    return jcall("cspr", x=x, y=y)


def cmidrk(x: Any) -> Any:
    """Midranks (C backend).

    The arguments and the result are those of R's ``cmidrk()``: see its help
    page.

    Examples
    --------
    >>> cmidrk(x=[3, 1, 2])
    [3.0, 1.0, 2.0]
    """
    return jcall("cmidrk", x=x)


def ccovm(x: Any) -> Any:
    """Covariance matrix (C backend).

    The arguments and the result are those of R's ``ccovm()``: see its help
    page.

    Examples
    --------
    >>> ccovm(x=[[1, 2], [2, 4], [3, 6]])
    [[1.0, 2.0], [2.0, 4.0]]
    """
    return jcall("ccovm", x=x)


def cbootmean(x: Any, B: Any = 1000, seed: Any = 42) -> Any:
    """Bootstrap replicate means (C backend).

    The arguments and the result are those of R's ``cbootmean()``: see its help
    page.

    Examples
    --------
    >>> r = cbootmean(x=[1, 2, 3], B=5, seed=1)
    >>> len(r)
    5
    """
    return jcall("cbootmean", x=x, B=B, seed=seed)


def cipw(treat: Any, ps: Any, trimlo: Any = 0.01, trimhi: Any = 0.99) -> Any:
    """Trimmed inverse-probability weights (C backend).

    The arguments and the result are those of R's ``cipw()``: see its help
    page.

    Examples
    --------
    >>> cipw(treat=[1, 0], ps=[0.5, 0.5])
    [2.0, 2.0]
    """
    return jcall("cipw", treat=treat, ps=ps, trimlo=trimlo, trimhi=trimhi)


def cgamcdf(shape: Any, x: Any) -> Any:
    """Regularized incomplete gamma function (C backend).

    The arguments and the result are those of R's ``cgamcdf()``: see its help
    page.

    Examples
    --------
    >>> cgamcdf(shape=1, x=[0])
    0.0
    """
    return jcall("cgamcdf", shape=shape, x=x)


def cunif(n: Any, seed: Any) -> Any:
    """Reproducible uniforms (splitmix64).

    The arguments and the result are those of R's ``cunif()``: see its help
    page.

    Examples
    --------
    >>> r = cunif(n=3, seed=42)
    >>> len(r)
    3
    """
    return jcall("cunif", n=n, seed=seed)


def toascii(x: Any) -> Any:
    """Transliterate text to plain ASCII.

    The arguments and the result are those of R's ``toascii()``: see its help
    page.

    Examples
    --------
    >>> toascii(x="Stra\u00dfe")
    'Strasse'
    """
    return jcall("toascii", x=x)


def asciifb(x: Any, force: Any = False) -> Any:
    """Text as it is, or its ASCII transliteration.

    The arguments and the result are those of R's ``asciifb()``: see its help
    page.

    Examples
    --------
    >>> asciifb(x="\u00c1ngela", force=True)
    'Angela'
    """
    return jcall("asciifb", x=x, force=force)


def wrtxtfb(text: Any, path: Any) -> Any:
    """Write text as UTF-8.

    The arguments and the result are those of R's ``wrtxtfb()``: see its help
    page.

    Examples
    --------
    >>> import os, tempfile
    >>> p = os.path.join(tempfile.mkdtemp(), "a.txt")
    >>> wrtxtfb(text=["one", "two"], path=p) == p
    True
    >>> open(p, "rb").read()
    b'one\\ntwo\\n'
    """
    return jcall("wrtxtfb", text=text, path=path)


def shrs(x: Any, count: Any = None, by: Any = None, total: Any = None, conflvl: Any = 0.95) -> Any:
    """Share of a total, in percent.

    The arguments and the result are those of R's ``shrs()``: see its help
    page.

    Examples
    --------
    >>> r = shrs(x=[3, 7])
    >>> list(r)
    ['count', 'total', 'share', 'lower', 'upper']
    """
    return jcall("shrs", x=x, count=count, by=by, total=total, conflvl=conflvl)


def deccodes(x: Any, vlabs: Any, keepna: Any = True) -> Any:
    """Decode numeric category codes by a code-to-label dictionary.

    The arguments and the result are those of R's ``deccodes()``: see its help
    page.

    Examples
    --------
    >>> deccodes(x=[2, 1], vlabs={"1": "a", "2": "b"})
    ['b', 'a']
    """
    return jcall("deccodes", x=x, vlabs=vlabs, keepna=keepna)


def grdrcd(x: Any, map: Any, keep: Any = None) -> Any:
    """Recode a categorical column by name.

    The arguments and the result are those of R's ``grdrcd()``: see its help
    page.

    Examples
    --------
    >>> grdrcd(x=["M", "F"], map={"M": "Male", "F": "Female"})
    ['Male', 'Female']
    """
    return jcall("grdrcd", x=x, map=map, keep=keep)


def relbl(x: Any, map: Any, keep: Any = None) -> Any:
    """Relabel a categorical variable by name, never by position.

    The arguments and the result are those of R's ``relbl()``: see its help
    page.

    Examples
    --------
    >>> relbl(x=["a", "b"], map={"a": "A", "b": "B"})
    ['A', 'B']
    """
    return jcall("relbl", x=x, map=map, keep=keep)


def grdlvls(x: Any, lvls: Any, ref: Any = None) -> Any:
    """A factor with declared, verified levels.

    The arguments and the result are those of R's ``grdlvls()``: see its help
    page.

    Examples
    --------
    >>> grdlvls(x=["b", "a"], lvls=["b", "a"])
    ['b', 'a']
    """
    return jcall("grdlvls", x=x, lvls=lvls, ref=ref)


def declbld(x: Any, vlabs: Any = None) -> Any:
    """Decode codes by their value labels.

    The arguments and the result are those of R's ``declbld()``: see its help
    page.

    Examples
    --------
    >>> declbld(x=[2, 1], vlabs={"1": "No", "2": "Yes"})
    ['Yes', 'No']
    """
    return jcall("declbld", x=x, vlabs=vlabs)


def grdbin(x: Any, col: Any) -> Any:
    """Refuse a categorical column where a 0/1 indicator is required.

    The arguments and the result are those of R's ``grdbin()``: see its help
    page.

    Examples
    --------
    >>> grdbin(x=[0, 1], col="d")
    True
    """
    return jcall("grdbin", x=x, col=col)


def vrfrcd(orig: Any, recd: Any, decl: Any) -> Any:
    """Prove a recode with a before/after cross-tabulation.

    The arguments and the result are those of R's ``vrfrcd()``: see its help
    page.

    Examples
    --------
    >>> r = vrfrcd(orig=["M"], recd=["Male"], decl={"M": "Male"})
    >>> list(r)
    ['original', 'recoded', 'Freq']
    """
    return jcall("vrfrcd", orig=orig, recd=recd, decl=decl)


def vrfmrg(x: Any, pub: Any, tol: Any = 0, strict: Any = True) -> Any:
    """Verify recoded counts against published counts.

    The arguments and the result are those of R's ``vrfmrg()``: see its help
    page.

    Examples
    --------
    >>> r = vrfmrg(x=["a", "b"], pub={"a": 1, "b": 1})
    >>> r["ok"]
    True
    """
    return jcall("vrfmrg", x=x, pub=pub, tol=tol, strict=strict)


def xfrvrf(imported: Any, soucnts: Any,
           codebook: Any = None, tol: Any = 0, strict: Any = True) -> Any:
    """Verify a categorical variable that crossed between programs.

    The arguments and the result are those of R's ``xfrvrf()``: see its help
    page.

    Examples
    --------
    >>> r = xfrvrf(imported=[1, 2], soucnts={"No": 1, "Yes": 1}, codebook={"1": "No", "2": "Yes"})
    >>> r["ok"]
    True
    """
    return jcall("xfrvrf", imported=imported, soucnts=soucnts, codebook=codebook, tol=tol, strict=strict)


def rgnmapintg(map: Any, unit: Any, reg: Any, regs: Any = None) -> Any:
    """Internal soundness of a region map.

    The arguments and the result are those of R's ``rgnmapintg()``: see its help
    page.

    Examples
    --------
    >>> r = rgnmapintg(map={"u": ["a", "a"], "r": ["1", ""]}, unit="u", reg="r")
    >>> r["observed"]
    [1.0, 1.0]
    """
    return jcall("rgnmapintg", map=map, unit=unit, reg=reg, regs=regs)


def rgnmapcmp(pub: Any, obs: Any, unit: Any, cols: Any = None) -> Any:
    """Compare a recomputed region map against a published one.

    The arguments and the result are those of R's ``rgnmapcmp()``: see its help
    page.

    Examples
    --------
    >>> r = rgnmapcmp(pub={"u": ["a"], "r": ["1"]}, obs={"u": ["a"], "r": ["2"]}, unit="u")
    >>> r["mismatched"]
    [0.0, 1.0]
    """
    return jcall("rgnmapcmp", pub=pub, obs=obs, unit=unit, cols=cols)


def rgnmapalt(map: Any, unit: Any, reg: Any, route: Any, known: Any = None) -> Any:
    """Check a region map against an independently derived assignment.

    The arguments and the result are those of R's ``rgnmapalt()``: see its help
    page.

    Examples
    --------
    >>> r = rgnmapalt(map={"u": ["a"], "r": ["1"]}, unit="u", reg="r", route={"a": "2"})
    >>> r["second"]
    ['2']
    """
    return jcall("rgnmapalt", map=map, unit=unit, reg=reg, route=route, known=known)


def drftks(x: Any, y: Any) -> Any:
    """Two-sample Kolmogorov-Smirnov test.

    The arguments and the result are those of R's ``drftks()``: see its help
    page.

    Examples
    --------
    >>> r = drftks(x=[1, 2, 3], y=[2, 3, 4])
    >>> list(r)
    ['statistic', 'p_value', 'neff']
    """
    return jcall("drftks", x=x, y=y)


def drftpsi(x: Any, y: Any, bins: Any = 10, eps: Any = 1e-06) -> Any:
    """Population stability index and Jensen-Shannon divergence.

    The arguments and the result are those of R's ``drftpsi()``: see its help
    page.

    Examples
    --------
    >>> r = drftpsi(x=[1, 2, 3, 4], y=[1, 2, 3, 4])
    >>> r["psi"]
    0.0
    """
    return jcall("drftpsi", x=x, y=y, bins=bins, eps=eps)


def drftchisq(obs: Any, expd: Any) -> Any:
    """Chi-square test of a categorical column against a known distribution.

    The arguments and the result are those of R's ``drftchisq()``: see its help
    page.

    Examples
    --------
    >>> r = drftchisq(obs={"a": 5, "b": 5}, expd={"a": 1, "b": 1})
    >>> r["statistic"]
    0.0
    """
    return jcall("drftchisq", obs=obs, expd=expd)


def drfthomog(x: Any, y: Any, seed: Any = 1) -> Any:
    """Chi-square test of homogeneity of two categorical samples.

    The arguments and the result are those of R's ``drfthomog()``: see its help
    page.

    Examples
    --------
    >>> r = drfthomog(x={"a": 30, "b": 20}, y={"a": 30, "b": 20})
    >>> r["statistic"]
    0.0
    """
    return jcall("drfthomog", x=x, y=y, seed=seed)


def fylbl(endyear: Any, sep: Any = "/", short: Any = True) -> Any:
    """Name a fiscal year by the years it spans.

    The arguments and the result are those of R's ``fylbl()``: see its help
    page.

    Examples
    --------
    >>> fylbl(endyear=[2023, 2024])
    ['2022/23', '2023/24']
    """
    return jcall("fylbl", endyear=endyear, sep=sep, short=short)


def yoylbl(x: Any, labels: Any) -> Any:
    """Label the periods of a change table.

    The arguments and the result are those of R's ``yoylbl()``: see its help
    page.

    Examples
    --------
    >>> t = {"prd": [2021, 2022], "value": [10, 12], "prev": [None, 10], "change": [None, 2], "pctchg": [None, 20]}
    >>> yoylbl(x=t, labels=["a", "b"])["perlbl"]
    ['a', 'b']
    """
    return jcall("yoylbl", x=x, labels=labels)


def yoysmry(x: Any, lag: Any = 1) -> Any:
    """Summarise a change table.

    The arguments and the result are those of R's ``yoysmry()``: see its help
    page.

    Examples
    --------
    >>> t = {"prd": [2021, 2022], "value": [10, 12], "prev": [None, 10], "change": [None, 2], "pctchg": [None, 20]}
    >>> [round(v, 9) for v in yoysmry(x=t)["totalpct"]]
    [20.0]
    """
    return jcall("yoysmry", x=x, lag=lag)


def yoycsv(x: Any, file: Any = None,
           digits: Any = None, na: Any = "", meta: Any = None, comment: Any = "#") -> Any:
    """Write a change table as CSV.

    The arguments and the result are those of R's ``yoycsv()``: see its help
    page.

    Examples
    --------
    >>> yoycsv(x={"code": ["046"], "n": [0.5]})
    'code,n\\n046,0.5\\n'
    """
    return jcall("yoycsv", x=x, file=file, digits=digits, na=na, meta=meta, comment=comment)


def yoytsv(x: Any, file: Any = None,
           digits: Any = None, na: Any = "", meta: Any = None, comment: Any = "#") -> Any:
    """Write a change table as TSV.

    The arguments and the result are those of R's ``yoytsv()``: see its help
    page.

    Examples
    --------
    >>> yoytsv(x={"code": ["046"], "n": [0.5]})
    'code\\tn\\n046\\t0.5\\n'
    """
    return jcall("yoytsv", x=x, file=file, digits=digits, na=na, meta=meta, comment=comment)


def yoyjs(x: Any, file: Any = None,
          digits: Any = None, pretty: Any = True, meta: Any = None) -> Any:
    """Write a change table as JSON.

    The arguments and the result are those of R's ``yoyjs()``: see its help
    page.

    Examples
    --------
    >>> yoyjs(x={"n": [0.5]}, pretty=False)
    '{"rows":[{"n":0.5}]}\\n'
    """
    return jcall("yoyjs", x=x, file=file, digits=digits, pretty=pretty, meta=meta)


def yoymd(x: Any, file: Any = None, digits: Any = 1, align: Any = True, meta: Any = None) -> Any:
    """Write a change table as Markdown.

    The arguments and the result are those of R's ``yoymd()``: see its help
    page.

    Examples
    --------
    >>> t = {"prd": [2021, 2022], "value": [10, 12], "prev": [None, 10], "change": [None, 2], "pctchg": [None, 20]}
    >>> yoymd(x=t).splitlines()[0]
    '| prd  | value | previous | change | change % | note |'
    """
    return jcall("yoymd", x=x, file=file, digits=digits, align=align, meta=meta)


def yoywr(x: Any, file: Any,
          format: Any = "auto", digits: Any = None, na: Any = None, meta: Any = None, comment: Any = None, pretty: Any = None, align: Any = None) -> Any:
    """Write a change table in the format its file name says.

    The arguments and the result are those of R's ``yoywr()``: see its help
    page.

    Examples
    --------
    >>> import os, tempfile
    >>> p = os.path.join(tempfile.mkdtemp(), "t.csv")
    >>> yoywr(x={"n": [1]}, file=p) == p
    True
    """
    return jcall("yoywr", x=x, file=file, format=format, digits=digits, na=na, meta=meta, comment=comment, pretty=pretty, align=align)


def co2int(loc: Any = None) -> Any:
    """Carbon intensity of electricity at a location, offline.

    The arguments and the result are those of R's ``co2int()``: see its help
    page.

    Examples
    --------
    >>> co2int(loc='CA-ON')['g_per_kwh']
    90.97
    """
    return jcall("co2int", loc=loc)


def co2inttab() -> Any:
    """The carbon-intensity table, every location it holds.

    The arguments and the result are those of R's ``co2inttab()``: see its
    help page.

    Examples
    --------
    >>> len(co2inttab()['location'])
    378
    """
    return jcall("co2inttab")


def dtctloc(envloc: Any = None, codloc: Any = None, tz: Any = None, locale: Any = None) -> Any:
    """The location of this computation, without network access.

    The arguments and the result are those of R's ``dtctloc()``: see its help
    page. Left out, an argument is read from the environment.

    Examples
    --------
    >>> dtctloc(envloc='', codloc='', tz='Europe/Paris')['location']
    'FR'
    """
    return jcall("dtctloc", envloc=envloc, codloc=codloc, tz=tz, locale=locale)


def raplok(rapldir: Any = None) -> Any:
    """Whether the CPU's RAPL energy counters are readable here.

    The arguments and the result are those of R's ``raplok()``: see its help
    page.

    Examples
    --------
    >>> raplok(rapldir='/no/such/dir')
    False
    """
    return jcall("raplok", rapldir=rapldir)


def fpeqv(co2eg: Any) -> Any:
    """Grams of CO2-equivalent in everyday terms.

    The arguments and the result are those of R's ``fpeqv()``: see its help
    page.

    Examples
    --------
    >>> round(fpeqv(co2eg=1750)['carkm'], 6)
    10.0
    """
    return jcall("fpeqv", co2eg=co2eg)


def exprspplr(data: Any, y: Any, expo: Any, confs: Any, nboot: Any = 100,
              rstate: Any = 42) -> Any:
    """Exposure-response estimate with a bootstrap sensitivity interval.

    The arguments and the result are those of R's ``exprspplr()`` with its
    built-in estimator: see its help page.

    Examples
    --------
    >>> import math
    >>> d = {'x': [math.sin(i) for i in range(40)],
    ...      'e': [math.sin(i) + math.cos(3 * i) for i in range(40)]}
    >>> d['y'] = [3 * e + 2 * x for e, x in zip(d['e'], d['x'])]
    >>> round(exprspplr(data=d, y='y', expo='e', confs='x', nboot=10)['ate'], 6)
    3.0
    """
    return jcall("exprspplr", data=data, y=y, expo=expo, confs=confs, nboot=nboot, rstate=rstate)


def inlhist(x: Any, bins: Any = 10) -> Any:
    """Inline histogram of a numeric vector.

    The arguments and the result are those of R's ``inlhist()``: see its help
    page.

    Examples
    --------
    >>> len(inlhist(x=[1, 2, 2, 3], bins=4))
    4
    """
    return jcall("inlhist", x=x, bins=bins)

def freqtbl(data: Any, col: Any = None, maxlvl: Any = 25, sort: Any = True) -> Any:
    """Frequency table of one column.

    The arguments and the result are those of R's ``freqtbl()``: see its help
    page.

    Examples
    --------
    >>> freqtbl(data=['a', 'b', 'b'])['n']
    [2.0, 1.0]
    """
    return jcall("freqtbl", data=data, col=col, maxlvl=maxlvl, sort=sort)

def cortbl(data: Any, meth: Any = None, minprs: Any = 3) -> Any:
    """Correlation of every pair of numeric columns.

    The arguments and the result are those of R's ``cortbl()``: see its help
    page.

    Examples
    --------
    >>> cortbl(data={'a': [1, 2, 3], 'b': [2, 4, 7]})['correlation']
    [1.0]
    """
    return jcall("cortbl", data=data, meth=meth, minprs=minprs)

def topcors(data: Any, n: Any = 10, meth: Any = None, minabs: Any = 0) -> Any:
    """The strongest pairwise correlations.

    The arguments and the result are those of R's ``topcors()``: see its help
    page.

    Examples
    --------
    >>> topcors(data={'a': [1, 2, 3], 'b': [3, 2, 1]})['abscor']
    [1.0]
    """
    return jcall("topcors", data=data, n=n, meth=meth, minabs=minabs)

def drpempty(data: Any, which: Any = None, quiet: Any = False) -> Any:
    """Drop rows and columns that are entirely missing.

    The arguments and the result are those of R's ``drpempty()``: see its help
    page.

    Examples
    --------
    >>> drpempty(data={'a': [1, None], 'b': [None, None]}, quiet=True)['a']
    [1.0]
    """
    return jcall("drpempty", data=data, which=which, quiet=quiet)

def drpconst(data: Any, naasval: Any = False, quiet: Any = False) -> Any:
    """Drop columns that hold one value.

    The arguments and the result are those of R's ``drpconst()``: see its help
    page.

    Examples
    --------
    >>> drpconst(data={'a': [1, 2], 'b': [7, 7]}, quiet=True)['a']
    [1.0, 2.0]
    """
    return jcall("drpconst", data=data, naasval=naasval, quiet=quiet)

def clncolnms(data: Any, case: Any = None, sep: Any = "_", quiet: Any = False) -> Any:
    """Tidy column names.

    The arguments and the result are those of R's ``clncolnms()``: see its help
    page.

    Examples
    --------
    >>> clncolnms(data=['Total Pop', '% change'], quiet=True)
    ['total_pop', 'pct_change']
    """
    return jcall("clncolnms", data=data, case=case, sep=sep, quiet=quiet)

def mahaoutl(data: Any, alpha: Any = 0.001, robust: Any = True) -> Any:
    """Multivariate outliers by Mahalanobis distance.

    The arguments and the result are those of R's ``mahaoutl()``: see its help
    page.

    Examples
    --------
    >>> mahaoutl(data={'x': [1, 2, 3, 4, 10], 'y': [2, 1, 4, 3, 0]})['row'][0]
    5.0
    """
    return jcall("mahaoutl", data=data, alpha=alpha, robust=robust)

def missruns(data: Any, minrun: Any = 2) -> Any:
    """Runs of consecutive missing values.

    The arguments and the result are those of R's ``missruns()``: see its help
    page.

    Examples
    --------
    >>> missruns(data={'a': [1, None, None, 4]})['length']
    [2.0]
    """
    return jcall("missruns", data=data, minrun=minrun)

def missgrid(data: Any, height: Any = 20, width: Any = 12) -> Any:
    """A text map of where the missing values are.

    The arguments and the result are those of R's ``missgrid()``: see its help
    page.

    Examples
    --------
    >>> len(missgrid(data={'a': [1, None]}, height=2))
    7
    """
    return jcall("missgrid", data=data, height=height, width=width)

def mcartst(data: Any, maxit: Any = 500, tol: Any = 1e-07) -> Any:
    """Little's test of missing completely at random.

    The arguments and the result are those of R's ``mcartst()``: see its help
    page.

    Examples
    --------
    >>> mcartst(data={'a': [1, 2, 3, 4], 'b': [2, 1, 4, 3]})['df']
    0.0
    """
    return jcall("mcartst", data=data, maxit=maxit, tol=tol)

def misspat(data: Any, maxpat: Any = 20) -> Any:
    """Which columns are missing together.

    The arguments and the result are those of R's ``misspat()``: see its help
    page.

    Examples
    --------
    >>> misspat(data={'a': [1, None], 'b': [1, 2]})['pattern']
    ['..', 'X.']
    """
    return jcall("misspat", data=data, maxpat=maxpat)

def duprows(data: Any, cols: Any = None, quiet: Any = False) -> Any:
    """Duplicated rows, with their groups.

    The arguments and the result are those of R's ``duprows()``: see its help
    page.

    Examples
    --------
    >>> duprows(data={'a': [1, 1, 2]}, quiet=True)['dupe_count']
    [2.0, 2.0]
    """
    return jcall("duprows", data=data, cols=cols, quiet=quiet)

def misssmry(data: Any) -> Any:
    """Missingness in one line per question.

    The arguments and the result are those of R's ``misssmry()``: see its help
    page.

    Examples
    --------
    >>> misssmry(data={'a': [1, None]})['n_missing']
    1.0
    """
    return jcall("misssmry", data=data)

def profcols(data: Any, qnts: Any = None) -> Any:
    """A report of every column of a data frame.

    The arguments and the result are those of R's ``profcols()``: see its help
    page.

    Examples
    --------
    >>> profcols(data={'v': [1, 2, 3]})['median']
    [2.0]
    """
    return jcall("profcols", data=data, qnts=qnts)

def dstntst(data: Any, grp: Any, est: Any, lower: Any, upper: Any, qual: Any = None,
            conflvl: Any = 0.95, pooled: Any = None) -> Any:
    """Whether distinct groups differ: Cochran's Q before pooling.

    The arguments and the result are those of R's ``dstntst()``: see its help
    page.

    Examples
    --------
    >>> d = {'g': ['a', 'b'], 'e': [1, 3], 'lo': [0, 2], 'hi': [2, 4]}
    >>> round(dstntst(data=d, grp='g', est='e', lower='lo', upper='hi')['pooled'], 6)
    2.0
    """
    return jcall("dstntst", data=data, grp=grp, est=est, lower=lower, upper=upper,
                 qual=qual, conflvl=conflvl, pooled=pooled)


def dstfit(x: Any, fams: Any = None, nboot: Any = 200, conflvl: Any = 0.95, seed: Any = 1,
           k: Any = None) -> Any:
    """Fit and rank candidate distributions by maximum likelihood.

    The arguments and the result are those of R's ``dstfit()``: see its help
    page.

    Examples
    --------
    >>> dstfit(x=[0.5, 1.2, 1.9, 2.4, 3.1, 4.0, 5.2, 7.7], fams='normal',
    ...        nboot=0)['fits']['normal']['mean']
    3.25
    """
    return jcall("dstfit", x=x, fams=fams, nboot=nboot, conflvl=conflvl, seed=seed, k=k)


def votagr(x: Any) -> Any:
    """Vote agreement between legislators.

    The arguments and the result are those of R's ``votagr()``: see its help
    page.

    Examples
    --------
    >>> r = votagr(x=[[0, 0, 0, 1, 0, 0], [0, 1, None, 0, 1, 1], [0, 0, 0, 1, 0, 1], [0, 0, 1, 0, 0,
    ... 0], [0, 1, 0, 1, 0, 0], [0, 1, 0, 1, 0, 1], [1, 0, 1, 0, 1, 0], [0, 1, 1, 0, 0, 1]])
    >>> list(r)
    ['agree', 'mean', 'n', 'm']
    """
    return jcall("votagr", x=x)


def amscal(Z: Any) -> Any:
    """Aldrich-McKelvey scaling.

    The arguments and the result are those of R's ``amscal()``: see its help
    page.

    Examples
    --------
    >>> r = amscal(Z=[[1, 3, 5, 6], [2, 3, 4, 7], [1, 2, 6, 6], [3, 4, 5, 7], [2, 4, 4, 6], [1, 3,
    ... 6, 7]])
    >>> list(r)
    ['zhat', 'alpha', 'beta', 'weights']
    """
    return jcall("amscal", Z=Z)


def bbscal(X: Any, ndim: Any = 2, minscale: Any = 8) -> Any:
    """Basic Space (blackbox) scaling.

    The arguments and the result are those of R's ``bbscal()``: see its help
    page.

    Examples
    --------
    >>> r = bbscal(X=[[1, 3, 5, 6, 3.5, 1.5], [2, 3, 4, 7, 2.5, 1.5], [1, 2, 6, 6, 2.5, 1.5], [3, 4,
    ... 5, 7, 3.5, 2.5], [2, 4, 4, 6, 4.5, 2.5], [1, 3, 6, 7, 2.5, 0.5]], ndim=1, minscale=4)
    >>> list(r)
    ['ideal', 'weights', 'eigen', 'singular', 'fit', 'cmeans', 'ndim']
    """
    return jcall("bbscal", X=X, ndim=ndim, minscale=minscale)


def optcls(votes: Any, ndim: Any = 1, maxit: Any = 500, nstart: Any = 10, seed: Any = 42) -> Any:
    """Optimal Classification scaling.

    The arguments and the result are those of R's ``optcls()``: see its help
    page.

    Examples
    --------
    >>> r = optcls(votes=[[0, 0, 0, 1, 0, 0], [0, 1, None, 0, 1, 1], [0, 0, 0, 1, 0, 1], [0, 0, 1,
    ... 0, 0, 0], [0, 1, 0, 1, 0, 0], [0, 1, 0, 1, 0, 1], [1, 0, 1, 0, 1, 0], [0, 1, 1, 0, 0,
    ... 1]], maxit=5, nstart=2)
    >>> list(r)
    ['ideal', 'normals', 'pre', 'errors', 'null', 'ndim']
    """
    return jcall("optcls", votes=votes, ndim=ndim, maxit=maxit, nstart=nstart, seed=seed)


def dblctr(D: Any) -> Any:
    """Double-centre a distance matrix.

    The arguments and the result are those of R's ``dblctr()``: see its help
    page.

    Examples
    --------
    >>> r = dblctr(D=[[0, 2.23606797749979, 3.1622776601683795, 5, 4.47213595499958],
    ... [2.23606797749979, 0, 2.23606797749979, 3.1622776601683795, 2.23606797749979],
    ... [3.1622776601683795, 2.23606797749979, 0, 2.23606797749979, 3.1622776601683795], [5,
    ... 3.1622776601683795, 2.23606797749979, 0, 2.23606797749979], [4.47213595499958,
    ... 2.23606797749979, 3.1622776601683795, 2.23606797749979, 0]])
    >>> list(r)
    ['B']
    """
    return jcall("dblctr", D=D)


def cmdsf(D: Any, ndim: Any = 2) -> Any:
    """Classical (metric) multidimensional scaling.

    The arguments and the result are those of R's ``cmdsf()``: see its help
    page.

    Examples
    --------
    >>> r = cmdsf(D=[[0, 2.23606797749979, 3.1622776601683795, 5, 4.47213595499958],
    ... [2.23606797749979, 0, 2.23606797749979, 3.1622776601683795, 2.23606797749979],
    ... [3.1622776601683795, 2.23606797749979, 0, 2.23606797749979, 3.1622776601683795], [5,
    ... 3.1622776601683795, 2.23606797749979, 0, 2.23606797749979], [4.47213595499958,
    ... 2.23606797749979, 3.1622776601683795, 2.23606797749979, 0]])
    >>> list(r)
    ['coords', 'eigen', 'stress', 'fit', 'B']
    """
    return jcall("cmdsf", D=D, ndim=ndim)


def smacf(
    D: Any,
    ndim: Any = 2,
    maxit: Any = 300,
    tol: Any = 1e-06,
    w: Any = None,
    init: Any = None,
) -> Any:
    """SMACOF stress minimisation.

    The arguments and the result are those of R's ``smacf()``: see its help
    page.

    Examples
    --------
    >>> r = smacf(D=[[0, 2.23606797749979, 3.1622776601683795, 5, 4.47213595499958],
    ... [2.23606797749979, 0, 2.23606797749979, 3.1622776601683795, 2.23606797749979],
    ... [3.1622776601683795, 2.23606797749979, 0, 2.23606797749979, 3.1622776601683795], [5,
    ... 3.1622776601683795, 2.23606797749979, 0, 2.23606797749979], [4.47213595499958,
    ... 2.23606797749979, 3.1622776601683795, 2.23606797749979, 0]], maxit=20)
    >>> list(r)
    ['coords', 'stress', 'iter', 'conv']
    """
    return jcall("smacf", D=D, ndim=ndim, maxit=maxit, tol=tol, w=w, init=init)


def nmdsf(D: Any, ndim: Any = 2, maxit: Any = 300, tol: Any = 1e-06) -> Any:
    """Nonmetric multidimensional scaling.

    The arguments and the result are those of R's ``nmdsf()``: see its help
    page.

    Examples
    --------
    >>> r = nmdsf(D=[[0, 2.23606797749979, 3.1622776601683795, 5, 4.47213595499958],
    ... [2.23606797749979, 0, 2.23606797749979, 3.1622776601683795, 2.23606797749979],
    ... [3.1622776601683795, 2.23606797749979, 0, 2.23606797749979, 3.1622776601683795], [5,
    ... 3.1622776601683795, 2.23606797749979, 0, 2.23606797749979], [4.47213595499958,
    ... 2.23606797749979, 3.1622776601683795, 2.23606797749979, 0]], maxit=10)
    >>> list(r)
    ['coords', 'stress', 'iter', 'conv']
    """
    return jcall("nmdsf", D=D, ndim=ndim, maxit=maxit, tol=tol)


def mdsfit(ev: Any) -> Any:
    """Fit of a scaling by its eigenvalues.

    The arguments and the result are those of R's ``mdsfit()``: see its help
    page.

    Examples
    --------
    >>> r = mdsfit(ev=[4, 2, -1, 0.5])
    >>> list(r)
    ['fit', 'cumfit', 'ev']
    """
    return jcall("mdsfit", ev=ev)


def unfstr(xr: Any, xs: Any, D: Any, w: Any = None) -> Any:
    """Stress of an unfolding.

    The arguments and the result are those of R's ``unfstr()``: see its help
    page.

    Examples
    --------
    >>> r = unfstr(xr=[[0, 0], [1, 1]], xs=[[1, 1], [0, 1], [2, 0]], D=[[1, None, 1], [2, 1, 1]])
    >>> list(r)
    ['stress']
    """
    return jcall("unfstr", xr=xr, xs=xs, D=D, w=w)


def mlsmu(D: Any, ndim: Any = 2, maxit: Any = 200, tol: Any = 1e-06, nstart: Any = 5) -> Any:
    """Least-squares metric unfolding.

    The arguments and the result are those of R's ``mlsmu()``: see its help
    page.

    Examples
    --------
    >>> r = mlsmu(D=[[1, 3, 5, 6], [2, 3, 4, 7], [1, 2, 6, 6], [3, 4, 5, 7], [2, 4, 4, 6], [1, 3, 6,
    ... 7]], maxit=10, nstart=2)
    >>> list(r)
    ['xr', 'xs', 'stress', 'iter', 'conv']
    """
    return jcall("mlsmu", D=D, ndim=ndim, maxit=maxit, tol=tol, nstart=nstart)


def smunf(D: Any, ndim: Any = 2, maxit: Any = 300, tol: Any = 1e-06) -> Any:
    """SMACOF unfolding.

    The arguments and the result are those of R's ``smunf()``: see its help
    page.

    Examples
    --------
    >>> r = smunf(D=[[1, 3, 5, 6], [2, 3, 4, 7], [1, 2, 6, 6], [3, 4, 5, 7], [2, 4, 4, 6], [1, 3, 6,
    ... 7]], maxit=10)
    >>> list(r)
    ['xr', 'xs', 'stress', 'iter', 'conv']
    """
    return jcall("smunf", D=D, ndim=ndim, maxit=maxit, tol=tol)


def ipsrec(xr: Any, xs: Any = None) -> Any:
    """Ideal points from an unfolding.

    The arguments and the result are those of R's ``ipsrec()``: see its help
    page.

    Examples
    --------
    >>> r = ipsrec(xr=[[1, 3], [2, 4]])
    >>> list(r)
    ['ideal']
    """
    return jcall("ipsrec", xr=xr, xs=xs)


def nomutl(x: Any, zy: Any, zn: Any, beta: Any = 15, w: Any = None) -> Any:
    """NOMINATE utilities and vote probabilities.

    The arguments and the result are those of R's ``nomutl()``: see its help
    page.

    Examples
    --------
    >>> r = nomutl(x=[[0, -0.5], [0.5, 0.2]], zy=[[0.3, 0.1]], zn=[[-0.3, 0]])
    >>> list(r)
    ['uyea', 'unay', 'udiff', 'prob']
    """
    return jcall("nomutl", x=x, zy=zy, zn=zn, beta=beta, w=w)


def nompr(x: Any, zy: Any, zn: Any, beta: Any = 15, w: Any = None) -> Any:
    """One NOMINATE vote probability.

    The arguments and the result are those of R's ``nompr()``: see its help
    page.

    Examples
    --------
    >>> r = nompr(x=0.1, zy=0.3, zn=-0.3)
    >>> list(r)
    ['prob']
    """
    return jcall("nompr", x=x, zy=zy, zn=zn, beta=beta, w=w)


def nomll(votes: Any, x: Any, zy: Any, zn: Any, beta: Any = 15, w: Any = None) -> Any:
    """NOMINATE log-likelihood and classification.

    The arguments and the result are those of R's ``nomll()``: see its help
    page.

    Examples
    --------
    >>> r = nomll(votes=[[1, None], [0, 1]], x=[[0.2], [-0.4]], zy=[[0.3], [-0.1]], zn=[[-0.3],
    ... [0.4]])
    >>> list(r)
    ['loglik', 'gmp', 'correct', 'total']
    """
    return jcall("nomll", votes=votes, x=x, zy=zy, zn=zn, beta=beta, w=w)


def procrs(X: Any, target: Any) -> Any:
    """Procrustes rotation.

    The arguments and the result are those of R's ``procrs()``: see its help
    page.

    Examples
    --------
    >>> r = procrs(X=[[0, 2.23606797749979], [2.23606797749979, 0], [3.1622776601683795,
    ... 2.23606797749979], [5, 3.1622776601683795], [4.47213595499958, 2.23606797749979]],
    ... target=[[2.23606797749979, 0], [0, 2.23606797749979], [2.23606797749979,
    ... 3.1622776601683795], [3.1622776601683795, 5], [2.23606797749979, 4.47213595499958]])
    >>> list(r)
    ['rotated', 'rotation', 'scale', 'mse']
    """
    return jcall("procrs", X=X, target=target)


def irtlik(votes: Any, x: Any, alpha: Any, beta: Any) -> Any:
    """Binary IRT likelihood.

    The arguments and the result are those of R's ``irtlik()``: see its help
    page.

    Examples
    --------
    >>> r = irtlik(votes=[[0, 0, 0, 1], [0, 1, None, 0], [0, 0, 0, 1]], x=[[-1], [0], [1]],
    ... alpha=[0, 0.5, -0.5, 0], beta=[[1], [2], [-1], [0.5]])
    >>> list(r)
    ['loglik', 'prob', 'correct', 'total', 'accuracy']
    """
    return jcall("irtlik", votes=votes, x=x, alpha=alpha, beta=beta)


def irtpst(chain: Any, std: Any = True) -> Any:
    """Posterior summaries of ideal-point draws.

    The arguments and the result are those of R's ``irtpst()``: see its help
    page.

    Examples
    --------
    >>> r = irtpst(chain=[[[0.010000000000000002, 1.6900000000000006], [0.25, 2.8900000000000006],
    ... [0.81, 4.41]], [[0.04000000000000001, 1.9600000000000004], [0.36, 3.240000000000001],
    ... [1, 4.840000000000001]], [[0.09000000000000002, 2.250000000000001], [0.4900000000000001,
    ... 3.6100000000000003], [1.2100000000000002, 5.290000000000001]], [[0.16000000000000003,
    ... 2.5600000000000005], [0.6400000000000001, 4], [1.4400000000000004, 5.76]]])
    >>> list(r)
    ['mean', 'sd', 'lo', 'hi', 'nsamp', 'std']
    """
    return jcall("irtpst", chain=chain, std=std)


def ordocl(Y: Any, ndim: Any = 2, maxit: Any = 500) -> Any:
    """Ordered Optimal Classification.

    The arguments and the result are those of R's ``ordocl()``: see its help
    page.

    Examples
    --------
    >>> r = ordocl(Y=[[1, 3, 5, 6], [2, 3, 4, 7], [1, 2, 6, 6], [3, 4, 5, 7], [2, 4, 4, 6], [1, 3,
    ... 6, 7]], ndim=1, maxit=5)
    >>> list(r)
    ['ideal', 'cuts', 'normals', 'correct', 'iter']
    """
    return jcall("ordocl", Y=Y, ndim=ndim, maxit=maxit)


def anchv(y: Any, V: Any, ncat: Any = 5) -> Any:
    """Anchoring vignettes.

    The arguments and the result are those of R's ``anchv()``: see its help
    page.

    Examples
    --------
    >>> r = anchv(y=[1, 3, 5, None, 2, 4], V=[[1, 3, 5], [2, 3, 4], [1, 2, 6], [3, 4, 5], [2, 4, 4],
    ... [1, 3, 6]])
    >>> list(r)
    ['scores', 'thresh', 'dif', 'order', 'n', 'nvig']
    """
    return jcall("anchv", y=y, V=V, ncat=ncat)


def indscl(D: Any, ndim: Any = 2, maxit: Any = 300, tol: Any = 1e-06) -> Any:
    """INDSCAL: scaling with individual differences.

    The arguments and the result are those of R's ``indscl()``: see its help
    page.

    Examples
    --------
    >>> r = indscl(D=[[[0, 2.23606797749979, 3.1622776601683795, 5, 4.47213595499958],
    ... [2.23606797749979, 0, 2.23606797749979, 3.1622776601683795, 2.23606797749979],
    ... [3.1622776601683795, 2.23606797749979, 0, 2.23606797749979, 3.1622776601683795], [5,
    ... 3.1622776601683795, 2.23606797749979, 0, 2.23606797749979], [4.47213595499958,
    ... 2.23606797749979, 3.1622776601683795, 2.23606797749979, 0]], [[0, 2.23606797749979,
    ... 3.1622776601683795, 5, 4.47213595499958], [2.6832815729997477, 0, 2.6832815729997477,
    ... 3.794733192202055, 2.6832815729997477], [2.8460498941515415, 2.012461179749811, 0,
    ... 2.012461179749811, 2.8460498941515415], [5.5, 3.478505426185218, 2.459674775249769, 0,
    ... 2.459674775249769], [4.47213595499958, 2.23606797749979, 3.1622776601683795,
    ... 2.23606797749979, 0]]], maxit=10)
    >>> list(r)
    ['config', 'weights', 'stress', 'iter', 'nind', 'nstim']
    """
    return jcall("indscl", D=D, ndim=ndim, maxit=maxit, tol=tol)


def nrmvec(x: Any, y: Any) -> Any:
    """Normal vector of an outside measure.

    The arguments and the result are those of R's ``nrmvec()``: see its help
    page.

    Examples
    --------
    >>> r = nrmvec(x=[[0, 2.23606797749979], [2.23606797749979, 0], [3.1622776601683795,
    ... 2.23606797749979], [5, 3.1622776601683795], [4.47213595499958, 2.23606797749979]], y=[1,
    ... 2, None, 4, 3])
    >>> list(r)
    ['normal', 'degrees', 'radians', 'r2', 'coef']
    """
    return jcall("nrmvec", x=x, y=y)


def cutlns(normals: Any, cuts: Any, xlim: Any = (-1, 1)) -> Any:
    """Cutting lines of votes in two dimensions.

    The arguments and the result are those of R's ``cutlns()``: see its help
    page.

    Examples
    --------
    >>> r = cutlns(normals=[[0.6, 0.8], [1, 0], [0, 0]], cuts=[0.1, -0.2, 0])
    >>> list(r)
    ['ends', 'angles', 'mids', 'nlines']
    """
    return jcall("cutlns", normals=normals, cuts=cuts, xlim=xlim)


def dwnom(votes: Any, ndim: Any = 2, maxit: Any = 100) -> Any:
    """DW-NOMINATE-style ideal points.

    The arguments and the result are those of R's ``dwnom()``: see its help
    page.

    Examples
    --------
    >>> r = dwnom(votes=[[0, 0, 0, 1, 0, 0], [0, 1, None, 0, 1, 1], [0, 0, 0, 1, 0, 1], [0, 0, 1, 0,
    ... 0, 0], [0, 1, 0, 1, 0, 0], [0, 1, 0, 1, 0, 1], [1, 0, 1, 0, 1, 0], [0, 1, 1, 0, 0, 1]],
    ... maxit=5)
    >>> list(r)
    ['ideal', 'weights', 'normals', 'mids', 'cuts', 'loglik', 'gmp', 'ndim']
    """
    return jcall("dwnom", votes=votes, ndim=ndim, maxit=maxit)


def nombt(votes: Any, x: Any, normals: Any, cuts: Any, nboot: Any = 100, seed: Any = 42) -> Any:
    """Parametric bootstrap of ideal points.

    The arguments and the result are those of R's ``nombt()``: see its help
    page.

    Examples
    --------
    >>> r = nombt(votes=[[0, 0, 0], [0, 1, None], [0, 0, 0], [0, 0, 1]], x=[[-1], [-0.5], [0.5],
    ... [1]], normals=[[1], [1], [-1]], cuts=[0, 0.2, -0.1], nboot=3)
    >>> list(r)
    ['se', 'mean', 'nboot']
    """
    return jcall("nombt", votes=votes, x=x, normals=normals, cuts=cuts, nboot=nboot, seed=seed)


def irtem(
    votes: Any,
    ndim: Any = 1,
    maxit: Any = 500,
    tol: Any = 1e-06,
    xprior: Any = (0, 1),
    bprior: Any = (0, 25),
    start: Any = None,
    conv: Any = 'cor',
) -> Any:
    """EM estimation of binary IRT ideal points.

    The arguments and the result are those of R's ``irtem()``: see its help
    page.

    Examples
    --------
    >>> r = irtem(votes=[[0, 0, 0, 1, 0, 0], [0, 1, None, 0, 1, 1], [0, 0, 0, 1, 0, 1], [0, 0, 1, 0,
    ... 0, 0], [0, 1, 0, 1, 0, 0], [0, 1, 0, 1, 0, 1], [1, 0, 1, 0, 1, 0], [0, 1, 1, 0, 0, 1]],
    ... maxit=20)
    >>> list(r)
    ['ideal', 'disc', 'diff', 'loglik', 'iter', 'conv']
    """
    return jcall(
        "irtem",
        votes=votes,
        ndim=ndim,
        maxit=maxit,
        tol=tol,
        xprior=xprior,
        bprior=bprior,
        start=start,
        conv=conv,
    )


def scboot(Z: Any, method: Any = 'am', nboot: Any = 200, seed: Any = 42) -> Any:
    """Bootstrap of scaling positions.

    The arguments and the result are those of R's ``scboot()``: see its help
    page.

    Examples
    --------
    >>> r = scboot(Z=[[1, 3, 5, 6], [2, 3, 4, 7], [1, 2, 6, 6], [3, 4, 5, 7], [2, 4, 4, 6], [1, 3,
    ... 6, 7]], nboot=3)
    >>> list(r)
    ['se', 'mean', 'lo', 'hi', 'nboot']
    """
    return jcall("scboot", Z=Z, method=method, nboot=nboot, seed=seed)


def wrdfsh(dtm: Any, maxit: Any = 100, tol: Any = 1e-06) -> Any:
    """Wordfish: positions from word counts.

    The arguments and the result are those of R's ``wrdfsh()``: see its help
    page.

    Examples
    --------
    >>> r = wrdfsh(dtm=[[5, 4, 3], [3, 1, 1], [0, 6, 0], [2, 2, 7]], maxit=10)
    >>> list(r)
    ['pos', 'beta', 'alpha', 'psi', 'loglik', 'iter']
    """
    return jcall("wrdfsh", dtm=dtm, maxit=maxit, tol=tol)


def bamscl(Z: Any, nsamp: Any = 1000, burn: Any = 200, polarity: Any = 1, seed: Any = 42) -> Any:
    """Bayesian Aldrich-McKelvey scaling.

    The arguments and the result are those of R's ``bamscl()``: see its help
    page.

    Examples
    --------
    >>> r = bamscl(Z=[[1, 3, 5, 6], [2, 3, 4, 7], [1, 2, 6, 6], [3, 4, 5, 7], [2, 4, 4, 6], [1, 3,
    ... 6, 7]], nsamp=5, burn=3)
    >>> list(r)
    ['zeta', 'sd', 'interval', 'a', 'b', 'tau', 'draws', 'nsamp']
    """
    return jcall("bamscl", Z=Z, nsamp=nsamp, burn=burn, polarity=polarity, seed=seed)


def bmdsf(
    D: Any,
    ndim: Any = 2,
    nsamp: Any = 1000,
    burn: Any = 200,
    sigma: Any = 1,
    seed: Any = 42,
) -> Any:
    """Bayesian metric multidimensional scaling.

    The arguments and the result are those of R's ``bmdsf()``: see its help
    page.

    Examples
    --------
    >>> r = bmdsf(D=[[0, 2.23606797749979, 3.1622776601683795, 5, 4.47213595499958],
    ... [2.23606797749979, 0, 2.23606797749979, 3.1622776601683795, 2.23606797749979],
    ... [3.1622776601683795, 2.23606797749979, 0, 2.23606797749979, 3.1622776601683795], [5,
    ... 3.1622776601683795, 2.23606797749979, 0, 2.23606797749979], [4.47213595499958,
    ... 2.23606797749979, 3.1622776601683795, 2.23606797749979, 0]], nsamp=3, burn=2)
    >>> list(r)
    ['pos', 'sd', 'dist', 'sigma', 'tau', 'draws', 'nsamp']
    """
    return jcall("bmdsf", D=D, ndim=ndim, nsamp=nsamp, burn=burn, sigma=sigma, seed=seed)


def bunfld(D: Any, ndim: Any = 2, nsamp: Any = 1000, burn: Any = 200, seed: Any = 42) -> Any:
    """Bayesian unfolding.

    The arguments and the result are those of R's ``bunfld()``: see its help
    page.

    Examples
    --------
    >>> r = bunfld(D=[[1, 3, 5, 6], [2, 3, 4, 7], [1, 2, 6, 6], [3, 4, 5, 7], [2, 4, 4, 6], [1, 3,
    ... 6, 7]], nsamp=3, burn=2)
    >>> list(r)
    ['stim', 'sd', 'ideal', 'dist', 'sigma', 'tau', 'nsamp']
    """
    return jcall("bunfld", D=D, ndim=ndim, nsamp=nsamp, burn=burn, seed=seed)


def cjrirt(
    votes: Any,
    ndim: Any = 1,
    nsamp: Any = 1000,
    burn: Any = 200,
    bvar: Any = 25,
    start: Any = None,
    seed: Any = 0,
    chains: Any = True,
) -> Any:
    """Clinton-Jackman-Rivers Bayesian IRT.

    The arguments and the result are those of R's ``cjrirt()``: see its help
    page.

    Examples
    --------
    >>> r = cjrirt(votes=[[0, 0, 0, 1, 0, 0], [0, 1, None, 0, 1, 1], [0, 0, 0, 1, 0, 1], [0, 0, 1,
    ... 0, 0, 0], [0, 1, 0, 1, 0, 0], [0, 1, 0, 1, 0, 1], [1, 0, 1, 0, 1, 0], [0, 1, 1, 0, 0,
    ... 1]], nsamp=3, burn=2, chains=False)
    >>> list(r)
    ['ideal', 'sd', 'alpha', 'beta', 'nsamp']
    """
    return jcall(
        "cjrirt",
        votes=votes,
        ndim=ndim,
        nsamp=nsamp,
        burn=burn,
        bvar=bvar,
        start=start,
        seed=seed,
        chains=chains,
    )


def anomin(
    votes: Any,
    ndim: Any = 1,
    nsamp: Any = 500,
    burn: Any = 100,
    seed: Any = 42,
    thin: Any = 1,
    lop: Any = 0.025,
    minvotes: Any = 20,
    polarity: Any = 1,
    fixalpha: Any = False,
) -> Any:
    """Alpha-NOMINATE.

    The arguments and the result are those of R's ``anomin()``: see its help
    page.

    Examples
    --------
    >>> r = anomin(votes=[[0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1], [0, 1, None, 0, 1, 1, 1, 0, None, 1,
    ... 0, 0], [0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0], [0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1], [0,
    ... 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1], [0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0], [1, 0, 1, 0, 1,
    ... 0, 0, 1, 0, 1, 0, 1], [0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0]], nsamp=2, burn=1,
    ... minvotes=5)
    >>> list(r)
    ['ideal', 'sd', 'yea', 'nay', 'alpha', 'interval', 'beta', 'draws', 'legis', 'votes',
    'ndim']
    """
    return jcall(
        "anomin",
        votes=votes,
        ndim=ndim,
        nsamp=nsamp,
        burn=burn,
        seed=seed,
        thin=thin,
        lop=lop,
        minvotes=minvotes,
        polarity=polarity,
        fixalpha=fixalpha,
    )


def ordirt(
    Y: Any,
    ndim: Any = 1,
    nsamp: Any = 500,
    burn: Any = 100,
    seed: Any = 42,
    lprec: Any = 0,
) -> Any:
    """Ordinal IRT (ordinal factor model).

    The arguments and the result are those of R's ``ordirt()``: see its help
    page.

    Examples
    --------
    >>> r = ordirt(Y=[[1, 3, 5, 6], [2, 3, 4, 7], [1, 2, 6, 6], [3, 4, 5, 7], [2, 4, 4, 6], [1, 3,
    ... 6, 7]], nsamp=3, burn=2)
    >>> list(r)
    ['ideal', 'sd', 'disc', 'icept', 'cuts', 'accept', 'nsamp']
    """
    return jcall("ordirt", Y=Y, ndim=ndim, nsamp=nsamp, burn=burn, seed=seed, lprec=lprec)


def dynirt(
    votes: Any,
    period: Any,
    nsamp: Any = 500,
    burn: Any = 100,
    seed: Any = 42,
    thin: Any = 1,
    tau2: Any = 1,
    c0: Any = -1,
    d0: Any = -1,
    e0: Any = 0,
    E0: Any = 1,
    a0: Any = 0,
    A0: Any = 0.1,
    b0: Any = 0,
    B0: Any = 0.1,
    anchor: Any = None,
) -> Any:
    """Dynamic IRT (ideal points that drift).

    The arguments and the result are those of R's ``dynirt()``: see its help
    page.

    Examples
    --------
    >>> r = dynirt(votes=[[0, 0, 0, 1, 0, 0], [0, 1, None, 0, 1, 1], [0, 0, 0, 1, 0, 1], [0, 0, 1,
    ... 0, 0, 0], [0, 1, 0, 1, 0, 0], [0, 1, 0, 1, 0, 1], [1, 0, 1, 0, 1, 0], [0, 1, 1, 0, 0,
    ... 1]], period=[1, 1, 2, 2, 3, 3], nsamp=3, burn=2)
    >>> list(r)
    ['theta', 'sd', 'alpha', 'beta', 'tau2', 'periods', 'nsamp', 'anchor']
    """
    return jcall(
        "dynirt",
        votes=votes,
        period=period,
        nsamp=nsamp,
        burn=burn,
        seed=seed,
        thin=thin,
        tau2=tau2,
        c0=c0,
        d0=d0,
        e0=e0,
        E0=E0,
        a0=a0,
        A0=A0,
        b0=b0,
        B0=B0,
        anchor=anchor,
    )


def arimfc(y: Any, h: Any = 12) -> Any:
    """ARIMA(1,1,1) forecast of a count series.

    The arguments and the result are those of R's ``arimfc()``: see its help
    page.

    Examples
    --------
    >>> r = arimfc(y=[52, 48, 61, 55, 47, 50, 66, 70, 58, 49, 45, 53, 57, 51, 63, 60, 49, 52, 71,
    ... 74, 61, 50, 47, 55, 59, 54, 66, 62, 51, 55, 73, 77, 63, 54, 49, 58], h=3)
    >>> list(r)
    ['fc', 'se', 'coef', 'sigma2', 'aic', 'bic', 'n']
    """
    return jcall("arimfc", y=y, h=h)


def sarfc(y: Any, h: Any = 12, order: Any = (1, 1, 1), seasonal: Any = (0, 1, 1, 12)) -> Any:
    """Seasonal ARIMA forecast with a hold-out check.

    The arguments and the result are those of R's ``sarfc()``: see its help
    page.

    Examples
    --------
    >>> r = sarfc(y=[52, 48, 61, 55, 47, 50, 66, 70, 58, 49, 45, 53, 57, 51, 63, 60, 49, 52, 71, 74,
    ... 61, 50, 47, 55, 59, 54, 66, 62, 51, 55, 73, 77, 63, 54, 49, 58], h=3, order=[0, 1, 1],
    ... seasonal=[0, 0, 0, 12])
    >>> list(r)
    ['fc', 'actual', 'coef', 'aic', 'bic', 'mape', 'rmse', 'ntrain']
    """
    return jcall("sarfc", y=y, h=h, order=order, seasonal=seasonal)


def catcor(units: Any) -> Any:
    """Correlation of counts per unit between categories.

    The arguments and the result are those of R's ``catcor()``: see its help
    page.

    Examples
    --------
    >>> r = catcor(units={'theft': ['a', 'b', 'b', 'c', 'd'], 'fraud': ['b', 'c', 'c', 'd'],
    ... 'arson': ['a', 'a', 'd']})
    >>> list(r)
    ['r', 'counts', 'units', 'cats']
    """
    return jcall("catcor", units=units)


def ousim(y: Any, npath: Any = 100, days: Any = 365, dt: Any = 1, seed: Any = 42) -> Any:
    """Ornstein-Uhlenbeck (Langevin) fit and simulation of a count series.

    The arguments and the result are those of R's ``ousim()``: see its help
    page.

    Examples
    --------
    >>> r = ousim(y=[52, 48, 61, 55, 47, 50, 66, 70, 58, 49, 45, 53, 57, 51, 63, 60, 49, 52, 71, 74,
    ... 61, 50, 47, 55, 59, 54, 66, 62, 51, 55, 73, 77, 63, 54, 49, 58], npath=2, days=3)
    >>> list(r)
    ['theta', 'mu', 'sigma', 'halflife', 'q', 'paths', 'n']
    """
    return jcall("ousim", y=y, npath=npath, days=days, dt=dt, seed=seed)


def fpdens(y: Any, ngrid: Any = 64, nstep: Any = 200) -> Any:
    """Fokker-Planck evolution of a count density.

    The arguments and the result are those of R's ``fpdens()``: see its help
    page.

    Examples
    --------
    >>> r = fpdens(y=[52, 48, 61, 55, 47, 50, 66, 70, 58, 49, 45, 53, 57, 51, 63, 60, 49, 52, 71,
    ... 74, 61, 50, 47, 55, 59, 54, 66, 62, 51, 55, 73, 77, 63, 54, 49, 58], ngrid=8, nstep=4)
    >>> list(r)
    ['theta', 'mu', 'sigma', 'grid', 'density', 'peak', 'statvar']
    """
    return jcall("fpdens", y=y, ngrid=ngrid, nstep=nstep)


def lvfit(y: Any, year: Any = None) -> Any:
    """Lotka-Volterra rates of a yearly count series.

    The arguments and the result are those of R's ``lvfit()``: see its help
    page.

    Examples
    --------
    >>> r = lvfit(y=[410, 455, 430, 470, 520, 498, 515], year=[2017, 2018, 2019, 2020, 2021, 2022,
    ... 2023])
    >>> list(r)
    ['alpha', 'beta', 'gamma', 'delta', 'period', 'year', 'prey', 'pred']
    """
    return jcall("lvfit", y=y, year=year)


def ppmlfe(data: Any, y: Any, x: Any, fe: Any, expo: Any = None, clus: Any = None, se: Any = "cluster",
           ssc: Any = "G", time: Any = None, lag: Any = None, nboot: Any = 199, seed: Any = None,
           sep: Any = "fe", level: Any = 0.95, tol: Any = 1e-12, maxit: Any = 100, twc: Any = "min",
           psd: Any = "none") -> Any:
    """Poisson pseudo-maximum likelihood with absorbed fixed effects.

    The arguments and the result are those of R's ``ppmlfe()``: see its help
    page. ``se="boot"`` needs ``seed`` outside R.

    Examples
    --------
    >>> d = {"y": [1, 2, 3, 5, 2, 3, 4, 7], "x": [0, 1, 2, 3, 0, 1, 2, 3], "a": ["p"] * 4 + ["q"] * 4}
    >>> r = ppmlfe(data=d, y="y", x="x", fe="a", se="hc1")
    >>> list(r)[:5]
    ['coef', 'vcov', 'n', 'nsep', 'xdrop']
    """
    return jcall("ppmlfe", data=data, y=y, x=x, fe=fe, expo=expo, clus=clus, se=se, ssc=ssc, time=time, lag=lag,
                 nboot=nboot, seed=seed, sep=sep, level=level, tol=tol, maxit=maxit, twc=twc, psd=psd)


def ivfe(data: Any, y: Any, endo: Any = None, inst: Any = None, x: Any = None, fe: Any = None, clus: Any = None,
         se: Any = "cluster", ssc: Any = "G", time: Any = None, lag: Any = None, nboot: Any = 199,
         seed: Any = None, level: Any = 0.95, tol: Any = 1e-12, twc: Any = "min", psd: Any = "none") -> Any:
    """Two-stage least squares with absorbed fixed effects.

    The arguments and the result are those of R's ``ivfe()``: see its help
    page (``fe`` is needed). ``se="boot"`` needs ``seed`` outside R.

    Examples
    --------
    >>> d = {"y": [1, 2, 3, 5, 2, 3, 4, 7], "x": [0, 1, 2, 3, 0, 1, 2, 3], "a": ["p"] * 4 + ["q"] * 4}
    >>> r = ivfe(data=d, y="y", x="x", fe="a", se="hc1")
    >>> list(r)[:5]
    ['coef', 'first', 'fstab', 'overid', 'gmm']
    """
    return jcall("ivfe", data=data, y=y, endo=endo, inst=inst, x=x, fe=fe, clus=clus, se=se, ssc=ssc, time=time,
                 lag=lag, nboot=nboot, seed=seed, level=level, tol=tol, twc=twc, psd=psd)


def ppmlcf(data: Any, y: Any, endo: Any, inst: Any, x: Any = None, fe: Any = None, expo: Any = None,
           clus: Any = None, how: Any = "cf", se: Any = "boot", ssc: Any = "G", nboot: Any = 199,
           seed: Any = None, level: Any = 0.95, tol: Any = 1e-12, maxit: Any = 100) -> Any:
    """Instrumental variables for a count outcome with fixed effects.

    The arguments and the result are those of R's ``ppmlcf()``: see its help
    page (``fe`` is needed). The bootstrap, the default, needs ``seed``
    outside R.

    Examples
    --------
    >>> d = {"y": [1, 2, 3, 5, 2, 3, 4, 7], "x": [0, 1, 2, 3, 0, 1, 2, 3], "z": [0, 1, 1, 3, 1, 0, 2, 3],
    ...      "a": ["p"] * 4 + ["q"] * 4}
    >>> r = ppmlcf(data=d, y="y", endo="x", inst="z", fe="a", se="naive")
    >>> r["coef"]["term"]
    ['x', 'vhat:x']
    """
    return jcall("ppmlcf", data=data, y=y, endo=endo, inst=inst, x=x, fe=fe, expo=expo, clus=clus, how=how, se=se,
                 ssc=ssc, nboot=nboot, seed=seed, level=level, tol=tol, maxit=maxit)


def windiv(dir: Any, width: Any = 90, start: Any = 0, area: Any = None, base: Any = 1, sect: Any = None,
           edges: Any = "halfopen", by: Any = None, avg: Any = "vector", share: Any = False,
           spd: Any = None) -> Any:
    """Wind-direction instruments.

    The arguments and the result are those of R's ``windiv()``: see its help
    page. The result is a matrix, a list of its rows, without names.

    Examples
    --------
    >>> windiv(dir=[10, 100, 200, 300])[1]
    [1.0, 0.0, 0.0]
    """
    return jcall("windiv", dir=dir, width=width, start=start, area=area, base=base, sect=sect, edges=edges, by=by,
                 avg=avg, share=share, spd=spd)


def dwnexp(sx: Any, sy: Any, sper: Any, sdir: Any, rx: Any, ry: Any, rid: Any = None, rule: Any = "sector",
           half: Any = 30, rad: Any = None, from_: Any = True, lonlat: Any = False, wrap: Any = True,
           dwt: Any = "mult") -> Any:
    """Downwind exposure: places downwind of sources of pollution.

    The arguments and the result are those of R's ``dwnexp()``: see its help
    page (``rad=None`` is R's ``Inf``; ``from_`` is R's ``from``).

    Examples
    --------
    >>> r = dwnexp(sx=[0], sy=[0], sper=None, sdir=[180], rx=[0, 60], ry=[100, 100])
    >>> [int(v) for v in r["any"]]
    [1, 0]
    """
    return jcall("dwnexp", sx=sx, sy=sy, sper=sper, sdir=sdir, rx=rx, ry=ry, rid=rid, rule=rule, half=half, rad=rad,
                 from_=from_, lonlat=lonlat, wrap=wrap, dwt=dwt)


def tailwd(dir: Any, spd: Any, course: Any = 250, by: Any = None, agg: Any = "max") -> Any:
    """The tailwind component of the wind along a flight course.

    The arguments and the result are those of R's ``tailwd()``: see its help
    page.

    Examples
    --------
    >>> [round(v, 6) for v in tailwd(dir=[70, 250], spd=[4, 4])]
    [4.0, -4.0]
    """
    return jcall("tailwd", dir=dir, spd=spd, course=course, by=by, agg=agg)


def invsn(tlow: Any, tup: Any, sign: Any = "upper", brk: Any = None, by: Any = None, rule: Any = "all") -> Any:
    """Temperature inversions.

    The arguments and the result are those of R's ``invsn()``: see its help
    page.

    Examples
    --------
    >>> list(invsn(tlow=[1, 2], tup=[3, 1]))
    ['diff', 'inv']
    >>> [int(v) for v in invsn(tlow=[1, 2, 1, 2], tup=[3, 4, 3, 1], by=[1, 1, 2, 2])["inv"]]
    [1, 0]
    """
    return jcall("invsn", tlow=tlow, tup=tup, sign=sign, brk=brk, by=by, rule=rule)


def updwn(data: Any, y: Any, ax: Any, ay: Any, day: Any, wind: Any, src: Any, area: Any = None, tol: Any = 60,
          band: Any = None, dist: Any = None, twin: Any = "drop", agg: Any = True, norm: Any = "mean",
          wx: Any = None, se: Any = "hc1", ssc: Any = "G", clus: Any = None, lag: Any = None, nboot: Any = 199,
          seed: Any = None, offs: Any = None, level: Any = 0.95) -> Any:
    """The upwind and downwind design around a road or a point source.

    The arguments and the result are those of R's ``updwn()``: see its help
    page (``band`` left out is R's ``Inf``).

    Examples
    --------
    >>> d = {"y": [3, 5, 4, 6, 6, 3, 7, 4], "ax": [-1] * 4 + [1] * 4, "ay": [0] * 8, "day": [1, 2, 3, 4] * 2,
    ...      "wind": [270, 90, 270, 90] * 2}
    >>> r = updwn(data=d, y="y", ax="ax", ay="ay", day="day", wind="wind",
    ...           src={"x0": [0], "y0": [-10], "x1": [0], "y1": [10]})
    >>> list(r)
    ['coef', 'class', 'n', 'fit', 'placebo', 'norm', 'tol', 'band']
    """
    return jcall("updwn", data=data, y=y, ax=ax, ay=ay, day=day, wind=wind, src=src, area=area, tol=tol, band=band,
                 dist=dist, twin=twin, agg=agg, norm=norm, wx=wx, se=se, ssc=ssc, clus=clus, lag=lag, nboot=nboot,
                 seed=seed, offs=offs, level=level)


def expjn(expo: Any, outc: Any, by: Any = ("area", "day")) -> Any:
    """Join an exposure table to counts by place and day.

    The arguments and the result are those of R's ``expjn()``: see its help
    page.

    Examples
    --------
    >>> r = expjn(expo={"area": [1, 1], "day": [1, 2], "pm": [5, 7]},
    ...           outc={"area": [1, 1], "day": [1, 3], "n": [2, 4]})
    >>> r["data"]["pm"], r["nout"]
    ([5.0], 1.0)
    """
    return jcall("expjn", expo=expo, outc=outc, by=by)


def smkiv(data: Any, plume: Any, unit: Any, day: Any, pm: Any = None, adj: Any = "raw", nbr: Any = None,
          z: Any = 1.64, seas: Any = None, dens: Any = None, cf: Any = "mean", by: Any = None) -> Any:
    """Smoke instruments from plume data.

    The arguments and the result are those of R's ``smkiv()``: see its help
    page. Outside R a day column is numbers, not dates, so ``pm`` needs
    ``seas``.

    Examples
    --------
    >>> smkiv(data={"pl": [0, 1, 0], "u": ["a", "a", "a"], "d": [1, 2, 3]}, plume="pl", unit="u", day="d")["smoke"]
    [0.0, 1.0, 0.0]
    """
    return jcall("smkiv", data=data, plume=plume, unit=unit, day=day, pm=pm, adj=adj, nbr=nbr, z=z, seas=seas,
                 dens=dens, cf=cf, by=by)


def exprun(data: Any, var: Any, unit: Any, day: Any, thr: Any = 0, cap: Any = 3, step: Any = 1) -> Any:
    """Runs of exposed days.

    The arguments and the result are those of R's ``exprun()``: see its help
    page.

    Examples
    --------
    >>> list(exprun(data={"s": [0, 1, 1, 0], "a": [1, 1, 1, 1], "d": [1, 2, 3, 4]}, var="s", unit="a", day="d"))
    ['run', 'run1', 'run2', 'run3p']
    """
    return jcall("exprun", data=data, var=var, unit=unit, day=day, thr=thr, cap=cap, step=step)


def lmdix(y: Any, q: Any, unit: Any, time: Any, type: Any = "time", at: Any = None, ref: Any = "mean",
          zero: Any = 0, fnm: Any = None) -> Any:
    """LMDI decomposition of a count over the factors of a Kaya identity.

    The arguments and the result are those of R's ``lmdix()``: see its help
    page (``q`` is a matrix as a list of rows, or a dict of columns).

    Examples
    --------
    >>> r = lmdix(y=[6, 12, 10, 20], q={"P": [3, 4, 5, 5]}, unit=["a", "a", "b", "b"], time=[1, 2, 1, 2])
    >>> list(r)
    ['dec', 'fac', 'ref', 'factors', 'type', 'refkind']
    """
    return jcall("lmdix", y=y, q=q, unit=unit, time=time, type=type, at=at, ref=ref, zero=zero, fnm=fnm)


def fepnl(y: Any, x: Any, unit: Any, time: Any, int: Any = None, clus: Any = None, at: Any = None,
          ssc: Any = "cr1", level: Any = 0.95) -> Any:
    """Two-way fixed effects with an interaction, for counts of offences.

    The arguments and the result are those of R's ``fepnl()``: see its help
    page (``x`` is a matrix as a list of rows, or a dict of columns).

    Examples
    --------
    >>> r = fepnl(y=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 2, 3.3, 4.7, 6, 7.3, 8.7, 10, 11.3, 12.7, 14],
    ...           x=[0] * 10 + [1, 2, 0, 1, 2, 0, 1, 2, 0, 1], unit=list(range(1, 11)) * 2, time=[1] * 10 + [2] * 10)
    >>> list(r)
    ['coef', 'slopes', 'fit', 'jb', 'vcov', 'resid', 'ssc', 'level']
    """
    return jcall("fepnl", y=y, x=x, unit=unit, time=time, int=int, clus=clus, at=at, ssc=ssc, level=level)


def pcdtst(x: Any, unit: Any, time: Any) -> Any:
    """Pesaran's CD test of cross-sectional dependence in a panel.

    The arguments and the result are those of R's ``pcdtst()``: see its help
    page.

    Examples
    --------
    >>> r = pcdtst(x=[1, 2, 3, 1, 3, 2, 3, 2, 1], unit=[1, 1, 1, 2, 2, 2, 3, 3, 3], time=[1, 2, 3] * 3)
    >>> list(r)
    ['cd', 'p', 'rho', 'absrho', 'pairs', 'units', 'periods']
    """
    return jcall("pcdtst", x=x, unit=unit, time=time)


def gmmpnl(y: Any, x: Any = None, unit: Any = None, time: Any = None, ylag: Any = 1, gmm: Any = None,
           iv: Any = None, z: Any = None, lagmin: Any = 2, lagmax: Any = None, collapse: Any = False,
           type: Any = "system", trans: Any = "fd", steps: Any = 1, robust: Any = True, wc: Any = None,
           small: Any = False, hmat: Any = "ab", cons: Any = True, tdum: Any = True, ar: Any = 2,
           level: Any = 0.95) -> Any:
    """Difference and system GMM for dynamic panels.

    The arguments and the result are those of R's ``gmmpnl()``: see its help
    page (``unit`` and ``time`` are needed; ``lagmax=None`` is R's ``Inf``).

    Examples
    --------
    >>> r = gmmpnl(y=[1, 3, 2, 5, 2, 4, 1, 3, 3, 1, 4, 2], unit=[1] * 4 + [2] * 4 + [3] * 4, time=[1, 2, 3, 4] * 3,
    ...            type="diff", tdum=False, collapse=True, lagmax=2, ar=0)
    >>> list(r)[:5]
    ['coef', 'tests', 'wald', 'xwald', 'fit']
    """
    return jcall("gmmpnl", y=y, x=x, unit=unit, time=time, ylag=ylag, gmm=gmm, iv=iv, z=z, lagmin=lagmin,
                 lagmax=lagmax, collapse=collapse, type=type, trans=trans, steps=steps, robust=robust, wc=wc,
                 small=small, hmat=hmat, cons=cons, tdum=tdum, ar=ar, level=level)


def wsctst(y: Any, x: Any = None, unit: Any = None, time: Any = None, cons: Any = False,
           ssc: Any = "cr1") -> Any:
    """Wooldridge's test for serial correlation in a panel.

    The arguments and the result are those of R's ``wsctst()``: see its help
    page (``unit`` and ``time`` are needed).

    Examples
    --------
    >>> r = wsctst(y=[0, 1, 3, 2, 1, 1, 2, 4, 0, 2, 1, 1], unit=[1] * 4 + [2] * 4 + [3] * 4, time=[1, 2, 3, 4] * 3)
    >>> list(r)
    ['stat', 'coef', 'ssc']
    """
    return jcall("wsctst", y=y, x=x, unit=unit, time=time, cons=cons, ssc=ssc)


def ccepnl(y: Any, x: Any, unit: Any, time: Any, fac: Any = None, cons: Any = True, w: Any = None,
           avg: Any = None, ylag: Any = 0, lags: Any = None, bias: Any = "none", pool: Any = None,
           pvar: Any = "pesaran", level: Any = 0.95) -> Any:
    """Common correlated effects for panels with unobserved common factors.

    The arguments and the result are those of R's ``ccepnl()``: see its help
    page.

    Examples
    --------
    >>> r = ccepnl(y=[1, 3, 2, 5, 2, 4, 1, 3, 3, 1, 4, 2, 5, 2, 3], x=[2, 1, 4, 3, 1, 5, 2, 2, 4, 3, 1, 2, 5, 3, 4],
    ...            unit=[1] * 5 + [2] * 5 + [3] * 5, time=[1, 2, 3, 4, 5] * 3)
    >>> list(r)
    ['mg', 'pooled', 'mg0', 'bi', 'vmg', 'vp', 'vmg0', 'hb', 'nfac', 'fit', 'w', 'bias', 'pvar', 'level']
    >>> d = ccepnl(y=[1, 3, 2, 5, 2, 4, 1, 3, 3, 1, 4, 2, 5, 2, 3, 2, 3, 4],
    ...            x={"x": [2, 1, 4, 3, 1, 5, 2, 2, 4, 3, 1, 2, 5, 3, 4, 1, 2, 3]}, unit=[1] * 6 + [2] * 6 + [3] * 6,
    ...            time=[1, 2, 3, 4, 5, 6] * 3, ylag=1, lags=0, avg=["x"], cons=False)
    >>> d["pooled"] is None, d["fit"]["lags"][0] == 0, d["mg"]["term"]
    (True, True, ['L1.y', 'x'])
    """
    return jcall("ccepnl", y=y, x=x, unit=unit, time=time, fac=fac, cons=cons, w=w, avg=avg, ylag=ylag, lags=lags,
                 bias=bias, pool=pool, pvar=pvar, level=level)
