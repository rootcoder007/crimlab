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
