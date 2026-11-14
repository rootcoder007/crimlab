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
