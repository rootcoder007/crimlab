/* SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
 *
 * libcrimlab: criminology and sociolegal statistics. The C interface.
 *
 * Every function is plain C from the outside: no exceptions cross it, no
 * memory is allocated for the caller, and nothing here knows about R. A
 * function returns CLAB_OK or an error code and writes its results into
 * arrays the caller owns. Indices are zero-based.
 */
#ifndef CLAB_H
#define CLAB_H

#ifdef __cplusplus
extern "C" {
#endif

//' @srrstats {G2.8} The kernels take double, int and logical arrays only;
//'   the binding converts every argument to one of these first.
//' @noRd
#define CLAB_OK 0
#define CLAB_EARG 1 /* a size or an index is out of range */
#define CLAB_EDOM 2 /* a value is outside the domain of the method */
#define CLAB_EINTR 3 /* the host's interrupt check asked to stop */

/* The words for an error code. */
const char *clab_strerror(int code);

/* Interrupts. Every long loop (resampling, permutation, bootstrap,
 * simulation, iterative fitting, scans over pairs) polls a check the host
 * gives here: f returns non-zero when the user asked to stop, and is
 * called once every ev element operations or so (ev < 0: the default,
 * 2^22, a few milliseconds; 0: at every poll). The kernel then returns
 * CLAB_EINTR, its outputs unspecified, having freed what it took. f NULL
 * (the default) never stops a kernel. Set it before any kernel runs: the
 * check and its interval are shared by all threads; whether a stop was
 * asked is kept per thread. */
void clab_onintr(int (*f)(void), double ev);

/* Non-zero when the check has asked to stop on this thread since the flag
 * was last cleared; clr non-zero clears it. A host clears it before each
 * call into the kernels. */
int clab_intrd(int clr);

/* Veil of darkness (Lean: Crimlab.Veil.odds_ratio_identifies_bias).
 * day, night: stops by group in daylight and in darkness, n groups; day
 * non-negative with day[ref] positive, night positive. shift is NULL or n
 * positive exposure factors. Out, n each: shday, shnight, ratio (the
 * visible-race factor against group ref), sign (-1, 0, 1: the sign of
 * shday - shnight) and, when shift is given, adj. */
int clab_vod(const double *day, const double *night, int n, int ref,
             const double *shift, double *shday, double *shnight,
             double *ratio, int *sign, double *adj);

/* The risks behind an odds ratio (Lean: Crimlab.RelRisk.or_eq_rr_mul).
 * orr > 0; q in (0, 1) is the overall rate and s in (0, 1) the exposed share.
 * Out: a (risk among the exposed) and b (risk among the unexposed), the one
 * pair in (0, 1) with odds ratio orr and s a + (1 - s) b = q. A risk closer
 * to 1 than a double can hold comes back as 1. 1 / orr must be finite. */
int clab_orrisk(double orr, double q, double s, double *a, double *b);

/* Outcome test for one group
 * (Lean: Crimlab.OutcomeTest.infra_marginality).
 * q: guilt probabilities in [0, 1]; w: non-negative weights with a positive
 * sum; t: the search threshold. Out: stop_rate (weight share with q > t) and
 * hit_rate (weighted mean of q among them; NaN when nobody is searched). */
int clab_outtst(const double *q, const double *w, int n, double t,
                double *stop_rate, double *hit_rate);

/* Benchmark decomposition (Lean: Crimlab.Benchmark.benchmark_product).
 * pop, contact, force: positive counts by group, n groups. Out, n each, all
 * against group ref: cd (contact per resident), fd (force per contact) and
 * rd (force per resident); rd = cd * fd. */
int clab_dspbm(const double *pop, const double *contact, const double *force,
               int n, int ref, double *cd, double *fd, double *rd);

/* Dyadic offending rates (Lean: Crimlab.Dyad.pair_exposure_rate_constant).
 * count[k] offences by group off[k] against group vic[k], m dyads; pop: g
 * positive group sizes. Out, m each: rtoff (per member of the offending
 * group), rtpair (per offender-victim pair exposed), nulloff (the rate
 * per offender under random mixing) and ratio (rtpair over the pooled
 * pair rate). */
int clab_dyadrt(const double *count, const int *off, const int *vic, int m,
                const double *pop, int g, double *rtoff,
                double *rtpair, double *nulloff, double *ratio);

/* Genotype match probability with the Balding-Nichols correction for
 * shared ancestry (NRC II recommendation 4.2). p, q: allele frequencies in
 * (0, 1] at n loci; hom[i] non-zero for a homozygote (q[i] is then not
 * read), otherwise p[i] + q[i] <= 1; theta in [0, 1). Out: locus, n match
 * probabilities, and log10mp, the base-10 logarithm of their product. */
int clab_gtmp(const double *p, const double *q, const int *hom, int n,
              double theta, double *locus, double *log10mp);

/* Log-likelihood-ratio cost of a set of likelihood ratios (Brummer and du
 * Preez 2006). same: ns positive ratios from same-source comparisons; diff:
 * nd positive ratios from different-source comparisons. Out: cllr, 0 for a
 * perfect system and 1 for one that always reports a ratio of 1. */
int clab_cllr(const double *same, int ns, const double *diff, int nd,
              double *cllr);

/* Calibration of likelihood ratios by pooling adjacent violators. lr: n
 * positive ratios in ascending order; same[i] non-zero for a same-source
 * comparison (both kinds must occur). work: 2n doubles of scratch. Out: post,
 * n calibrated probabilities of same source at the sample's own prior, and
 * cllrmin, the cost of the calibrated ratios: what the system would cost if
 * its ratios kept their order and lost their miscalibration. */
int clab_lrcal(const double *lr, const int *same, int n, double *post,
               double *work, double *cllrmin);

/* Mixture: the chance that x unknown contributors carry only alleles of the
 * mixture and between them every allele flagged in need (Weir et al. 1997,
 * no coancestry). p: n allele frequencies of the mixture, sum at most 1;
 * at most 20 alleles may be flagged. Out: the probability. */
int clab_mixp(const double *p, const int *need, int n, int x, double *out);

/* Kinship index by locus: the chance of two genotypes under a relationship
 * with identity-by-descent probabilities k[0..2] over their chance if
 * unrelated. g1, g2: 2n allele codes, two per locus; f: 2n frequencies of
 * the alleles of g1. Out: idx, n indices. */
int clab_kinidx(const int *g1, const int *g2, const double *f, int n,
                const double *k, double *idx);

/* Paternity index by locus for a mother, child and alleged father, without
 * mutation. c, m, af: 2n allele codes each; f: 2n frequencies of the
 * child's alleles. Out: idx, n indices; NaN where the mother cannot be the
 * mother. */
int clab_patidx(const int *c, const int *m, const int *af, const double *f,
                int n, double *idx);

/* Clearance-rate decomposition and Simpson's paradox (Lean:
 * Crimlab.Clearance.kitagawa, change_ge_mix). clra, reca, clrb, recb: n parts'
 * cleared (non-negative) and recorded (positive) counts in comparisons A
 * and B. Out, n each: ra, rb (rates), wa, wb (workload shares), reff and
 * meff (each part's rate and mix terms, the mix centred on the mean of the
 * two aggregate rates); tot, 5: the aggregate rates of A and B, their
 * change, the summed rate effect and the summed mix effect, change = sum
 * of the two; cnt, 3: parts whose rate rose, parts whose rate fell, and the
 * kind of reversal (0 none, 1 the change and the rate effect differ in
 * sign, 2 every part moved one way and the aggregate the other). */
int clab_clrsmp(const double *clra, const double *reca, const double *clrb,
                const double *recb, int n, double *ra, double *rb,
                double *wa, double *wb, double *reff, double *meff,
                double *tot, int *cnt);

/* Rates and disparities when exposure is known only within a factor of a
 * proxy (Lean: Crimlab.Exposure.rate_bounds, disparity_bounds,
 * disparity_sign_identified). y, m: n positive counts and proxies; gamma
 * >= 1. Out, n each: rate (y / m), lo, hi (its bounds), ratio (against
 * group ref), rlo, rhi (its bounds) and sign (1 when the true ratio exceeds
 * 1 for every admissible exposure, -1 when it is below 1, 0 otherwise). */
int clab_expbnd(const double *y, const double *m, int n, int ref, double gamma,
                double *rate, double *lo, double *hi, double *ratio,
                double *rlo, double *rhi, int *sign);

/* Risk difference under nondifferential misclassification of a binary
 * group (Lean: Crimlab.Misclass.binary_attenuation, binary_sign). r1, r0: true
 * risks in [0, 1]; p in (0, 1) the true prevalence; se, sp in [0, 1], with
 * both measured groups non-empty. Out, 5: the true difference, the measured
 * difference (from the mixtures), the attenuation factor
 * (se + sp - 1) p (1 - p) / (d1 d0), d1 and d0 (measured group shares). */
int clab_misrd(double r1, double r0, double p, double se, double sp,
               double *out);

/* Regression to the mean at places selected on x1 > c (Lean:
 * Crimlab.RegMean.selected_change_nonpos, low_selected_change_nonneg).
 * x1, x2: n finite counts; w: n non-negative weights, positive on the
 * selected places. Out, 7: the selected weight, their mean change, the
 * mirror change (selection on x2 > c, read backwards), the symmetrised
 * change (never positive), the excess of the first over it, the low
 * places' change (x1 < c) and its symmetrised value (never negative); NaN
 * where no place qualifies. */
int clab_rtmhs(const double *x1, const double *x2, const double *w, int n,
               double c, double *out);

/* Incapacitation (Lean: Crimlab.Incap.steady_state_rate,
 * prevented_share_eq, marginal_prevention_eq). lam > 0, q in (0, 1], s >= 0,
 * n groups; share is NULL or n non-negative shares summing to 1. Out, n
 * each: rate (lam / (1 + lam q s)), prev (the share prevented) and marg (the
 * crimes one more year prevents); agg, 3 when share is given: the free
 * rate, the incapacitated rate and the share prevented, all weighted. */
int clab_incap(const double *lam, const double *q, const double *s, int n,
               const double *share, double *rate, double *prev, double *marg,
               double *agg);

/* Incapacitation under desistance and replacement (Lean:
 * Crimlab.Incap.prevented_le_const, prevented_ge_const,
 * later_sentence_prevents_less, prevented_net_le). lam: len non-negative,
 * non-increasing rates by career age; a sentence of s >= 1 periods from age
 * t0 >= 0, with t0 + s <= len - 1; r in [0, 1] the replacement share. Out,
 * 5: crimes prevented, net of replacement, the bounds s lam[t0] and
 * s lam[t0 + s], and the crimes the same sentence prevents one period
 * later. */
int clab_incapr(const double *lam, int len, int t0, int s, double r,
                double *out);

/* Selective labels (Lean: Crimlab.Select.unobserved_bounds,
 * unobserved_width, unobserved_ends_attained). y: n outcomes, 0 or 1 where
 * rel[i] is non-zero (not read elsewhere); rel, rule: the judge's and the
 * rule's releases; w: non-negative weights, both release sets of positive
 * weight. Out, 6: the judge's observed failure rate, the unobserved share
 * of the rule's releases, the rule's failure-rate bounds, their width and
 * the naive rate on the rule's observed releases (NaN when none). */
int clab_sellab(const double *y, const int *rel, const int *rule,
                const double *w, int n, double *out);

/* Mean spell length from completed and open spells (Lean:
 * Crimlab.Censor.true_mean_ge, lower_bound_sub, bias_lower). t: n >= 1
 * completed durations; a: m >= 1 elapsed ages of open spells; all
 * non-negative. Out, 4: the completed mean, the open spells' mean age, the
 * lower bound on the cohort mean and (m / (n + m)) (abar - tbar). */
int clab_losbnd(const double *t, int n, const double *a, int m, double *out);

/* No-assumption bounds on a sentence's effect on a binary outcome (Lean:
 * Crimlab.Sentence.outcome_bounds, ate_width_one, ate_contains_zero). y: n
 * outcomes, 0 or 1; z: 0 (comparison) or 1 (treatment); w: non-negative
 * weights with a positive total. Out, 11: P(y = 1, z = 0), P(y = 1, z = 1),
 * P(z = 0), P(z = 1), the bounds on P(y(0) = 1), on P(y(1) = 1) and on the
 * effect, and the naive difference (NaN when an arm is empty). */
int clab_sntbnd(const double *y, const int *z, const double *w, int n,
                double *out);

/* Period-2 hazards after selective depletion (Lean:
 * Crimlab.Hazard.survivor_hazard_mono, hr2_gt_one_of_depletion). s in (0, 1)
 * the high-risk share; 0 <= l < h <= 1 the two types' hazards; ah, al: 2
 * period-1 survival probabilities in (0, 1] of each type in the control and
 * treated arms. Out, 5: the surviving high-risk share in each arm, the
 * period-2 hazard in each arm and their ratio (treated over control). */
int clab_hzdsel(double s, double h, double l, const double *ah,
                const double *al, double *out);

/* ------------------------------------------------------------------------
 * Spatial kernels: units. Every planar spatial kernel below takes
 * projected coordinates in one length unit (metres, say) and returns
 * lengths in that unit and areas in its square. None of them can tell
 * longitude and latitude from metres, so lon/lat input must go through
 * clab_eqproj first (or clab_hvdist for distances alone); a distance or
 * bandwidth given in metres against degrees is the classic silent error.
 * ------------------------------------------------------------------------ */

/* Great-circle distances on a sphere. lon1, lat1, lon2, lat2: n pairs of
 * points in degrees (latitude in [-90, 90], longitude any finite value);
 * r > 0 the sphere's radius (6371008.8 m is the mean Earth radius). Out:
 * d, n distances in the unit of r. The central angle is atan2 of the cross
 * and dot products of the two unit vectors, written so that nothing
 * cancels: accurate to a few ulps from coincident points to antipodes
 * (Vincenty 1975, doi:10.1179/sre.1975.23.176.88, the sphere case). O(n). */
int clab_hvdist(const double *lon1, const double *lat1, const double *lon2,
                const double *lat2, int n, double r, double *d);

/* Lambert azimuthal equal-area projection of the sphere about (lon0,
 * lat0) (Snyder 1987, doi:10.3133/pp1395, ch. 24). lon, lat: n points in
 * degrees; r > 0 the radius. Out: x (east), y (north), n each, in the unit
 * of r. Areas are exact, so intensities, nearest-neighbour expectations and
 * K functions computed on x, y are right; local lengths are scaled by
 * cos(c / 2) towards the centre and 1 / cos(c / 2) across (c the central
 * angle), so within 100 km of the centre on the Earth every length is right
 * to 3.1e-5 relative (2R sin(c / 2) against R c for distances from the
 * centre).
 * CLAB_EDOM for a point within about 1e-6 rad of the antipode. O(n). */
int clab_eqproj(const double *lon, const double *lat, int n, double lon0,
                double lat0, double r, double *x, double *y);

/* Centrographic summary (Bachi 1959; Neft 1966). x, y: n >= 1 points; w:
 * NULL or n finite non-negative weights with a positive sum W (frequency
 * weights); k in {0, 1, 2}: the dispersions divide by W - k, which must be
 * positive. k = 1 makes sx^2, sy^2 and sdd^2 unbiased for the variances
 * and for E|X - mu|^2 (each coordinate loses one degree of freedom to its
 * mean); k = 2 is the convention of some crime-mapping texts, which take
 * the two means as costing two degrees of freedom from each sum. Out, 10:
 * W, the (weighted) mean centre mx, my (the point minimising the weighted
 * sum of squared distances), sx, sy, the standard distance deviation sdd =
 * sqrt(sx^2 + sy^2 with the same divisor), the geometric means gx, gy and
 * the harmonic means hx, hy of the coordinates (NaN for a coordinate with
 * a non-positive value under positive weight: they are not defined there,
 * and they change with the origin of the coordinates). Moments are taken
 * about the refined mean with compensated sums. O(n). */
int clab_cgstat(const double *x, const double *y, const double *w, int n,
                int k, double *out);

/* Median centre: the median of x and the median of y (the mean of the two
 * middle values for even n). Not rotation invariant: it depends on the
 * axes. work: n doubles. Out, 2. O(n log n). */
int clab_cgmedn(const double *x, const double *y, int n, double *work,
                double *out);

/* Centre of minimum distance (the Weber point, spatial median): the point
 * c minimising sum w_i |p_i - c|. Weiszfeld's iteration (1937, translated
 * doi:10.1007/s10479-008-0352-z) with the modification of Vardi and Zhang
 * (2000, doi:10.1073/pnas.97.4.1423), which keeps each step defined and
 * the objective non-increasing when an iterate falls on a data point
 * (including coincident points, whose weights add). x, y, w as for
 * clab_cgstat; start: NULL (the weighted mean centre) or a starting point;
 * tol > 0: stop when a step is at most tol times the root mean squared
 * distance to the mean centre; maxit >= 0 steps at most. Out, 6: cx, cy,
 * the objective sum w_i |p_i - c|, the optimality gap max(0, |R| - eta)
 * (R the sum of unit vectors from c to the points, eta the weight at c;
 * 0 if and only if c is a minimiser), the number of steps and a code: 2 at
 * an exact minimiser, 1 when the step fell below tol, 0 when maxit was
 * reached. The minimiser is unique unless all points are collinear with
 * an even split of weight, when every point of a segment minimises. O(n)
 * per step. */
int clab_webpt(const double *x, const double *y, const double *w, int n,
               const double *start, double tol, int maxit, double *out);

/* Standard deviational ellipse, from the weighted covariance matrix C =
 * (sum w_i d_i d_i^T) / (W - k) of the points about their mean centre: the
 * axes lie along the eigenvectors of C and the semi-axes are c times the
 * square roots of its eigenvalues (Wang, Shi and Miao 2015,
 * doi:10.1371/journal.pone.0118537). This is the same ellipse as the
 * rotated-axes formulas of Lefever (1926) and Ebdon (1988), whose angle
 * minimises the sum of squared distances to the rotated axis, once the same
 * divisor is used; texts differ only in that divisor (k) and in the
 * multiplier c (some multiply by sqrt 2). Gong (2002,
 * doi:10.1111/j.1538-4632.2002.tb01082.x) shows that Lefever's curve of the
 * standard deviation by direction is not this ellipse; the ellipse is
 * returned, not the curve. For bivariate normal points the ellipse with c =
 * sqrt(-2 log(1 - p)) covers a share p (c = 1: 39.3%, c = 2: 86.5%). x, y,
 * w, k as for clab_cgstat; c > 0. Out, 6: mx, my, a >= b (semi-axes), the
 * angle of the major axis in [0, pi) anticlockwise from the x axis (0 when
 * a = b), and the area pi a b. O(n). */
int clab_sdell(const double *x, const double *y, const double *w, int n, int k,
               double c, double *out);

/* Convex hull by Andrew's monotone chain (Andrew 1979,
 * doi:10.1016/0020-0190(79)90072-3) with an exact orientation predicate
 * (floating-point filter, then an exact expansion; Shewchuk 1997,
 * doi:10.1007/PL00009321), so near-collinear points never make a wrong
 * turn. x, y: n >= 1 points with each coordinate 0 or of magnitude in
 * [1e-140, 1e150] (the range in which the predicate is exact). idx: n ints
 * of workspace; hull: 2n ints, of which the first *nh receive the vertex
 * indices anticlockwise from the lowest-x (then lowest-y) point. Points on
 * an edge are not vertices; of coincident points the one with the smallest
 * index is reported. Out, 2: area (shoelace) and perimeter; one distinct
 * point gives nh = 1 and 0, 0; collinear points give their two ends, area 0
 * and twice the segment's length. O(n log n). */
int clab_cvhull(const double *x, const double *y, int n, int *idx, int *hull,
                int *nh, double *out);


/* Nearest-neighbour index (Clark and Evans 1954, doi:10.2307/1931034) and
 * k-order nearest-neighbour distances in a rectangle. x, y: n >= 2 points
 * in planar units inside win = {x0, x1, y0, y1} (x0 < x1, y0 < y1); the
 * intensity is lam = n / area. kmax in [1, n - 1] orders; nsim >= 0
 * simulations of n uniform points in the rectangle for a Monte Carlo test
 * (0: none), drawn from seed. dist: NULL or n * kmax doubles receiving
 * every point's k-th nearest-neighbour distance (row i holds orders 1 to
 * kmax; coincident points are neighbours at distance 0). iwork: 2n + 2 +
 * 2 kmax ints; work: 3 kmax doubles, plus 2n when nsim > 0. Out, 6 kmax +
 * 4: for each order k, six values: the mean k-th distance, its expectation
 * under complete spatial randomness without edge effects, c_k / sqrt(lam)
 * with c_k = k (2k)! / (2^k k!)^2 (Thompson 1956, doi:10.2307/1933159),
 * their ratio, the Monte Carlo p-values (1 + #{sim <= obs}) / (nsim + 1)
 * against clustering and (1 + #{sim >= obs}) / (nsim + 1) against
 * regularity (NaN without simulations), and the standard deviation of one
 * k-th distance under CSR, sqrt(k / (pi lam) - c_k^2 / lam). Then, for
 * order 1: Donnelly's (1978) edge-corrected expectation 0.5 / sqrt(lam) +
 * (0.0514 + 0.0412 / sqrt(n)) P / n (P the perimeter; the constants as
 * spatstat uses them), the ratio against it, and the Clark-Evans z
 * statistics with the standard error sqrt((4 - pi) / (4 pi lam n)) about
 * the uncorrected and the corrected expectation. The z statistics ignore
 * that nearest-neighbour distances are dependent and, uncorrected, the
 * edge; the Monte Carlo p-value has exact size in the rectangle. The
 * index tests only the k-th neighbour distances, not complete spatial
 * randomness. Neighbours are found on a uniform grid of about n / 2 cells:
 * O(n kmax) expected for spread-out points, O(n log n)-like for clustered
 * ones; nsim + 1 times that with simulations. */
int clab_nnidx(const double *x, const double *y, int n, const double *win,
               int kmax, int nsim, unsigned int seed, double *dist,
               int *iwork, double *work, double *out);

/* Ripley's K function (Ripley 1976, doi:10.2307/3212829; 1977,
 * doi:10.1111/j.2517-6161.1977.tb01615.x) of n >= 2 points in the
 * rectangle win (as for clab_nnidx), at m >= 1 radii r (non-negative,
 * strictly increasing), with the estimator corr: 0 none, A sum_{i != j}
 * 1{d_ij <= r} / (n (n - 1)); 1 border (reduced sample: only points at
 * least r from the edge are centres, intensity n / A); 2 modified border,
 * A^2 sum over those centres / (n (n - 1) |W eroded by r|) (NaN once r
 * reaches half a side); 3 isotropic, each pair weighted by the inverse of
 * the share of the circle about the centre through the other point that
 * lies in the rectangle (exact for a rectangle: the arcs cut by each edge
 * less the overlaps at corners inside the circle), normalised by n (n - 1)
 * as spatstat does. Besag's L is sqrt(K / pi) (Besag 1977, in the
 * discussion of Ripley 1977). Pairs are found on a grid of cells of side
 * at least max(r), so the cost is O(n + pairs within max(r)) plus O(n m)
 * for the border estimators; no n x n matrix. iwork: 2n + 2 ints; work:
 * 4m doubles. Out: k, m values. Units: planar; A is in the square of the
 * coordinates' unit. */
int clab_ripk(const double *x, const double *y, int n, const double *win,
              const double *r, int m, int corr, int *iwork, double *work,
              double *k);

/* Simulation envelopes of L(r) - r = sqrt(K(r) / pi) - r and a global
 * Monte Carlo test of complete spatial randomness in the rectangle. Inputs
 * as for clab_ripk, nsim >= 1 patterns of n uniform points in win drawn
 * from seed. iwork: 2n + 2 ints; work: 5m + 2n doubles. Out, 4m + 2: the
 * observed L(r) - r, the pointwise minimum and maximum over the
 * simulations (pointwise envelopes: each radius alone is outside them
 * with probability 2 / (nsim + 1) under CSR, which is not a test of the
 * whole curve), the simulations' mean, then the maximum absolute
 * deviation T = max_r |L(r) - r| of the data and its Monte Carlo p-value
 * (1 + #{T_sim >= T}) / (nsim + 1), which has exact size (Ripley 1977;
 * Diggle 2013, ch. 6). Cost: nsim + 1 times clab_ripk. */
int clab_ripenv(const double *x, const double *y, int n, const double *win,
                const double *r, int m, int corr, int nsim, unsigned int seed,
                int *iwork, double *work, double *out);


/* ------------------------------------------------------------------------
 * Journey to crime. What these kernels do and do not license: given the
 * locations of a series of incidents attributed to one offender, they rank
 * candidate locations of the offender's base under the assumptions that
 * (a) the series really is one offender's, (b) the offender travelled
 * from a single fixed base inside the grid, (c) travel lengths follow the
 * decay function supplied, which describes offenders in general, not this
 * one, and (d) for the product surface, the incidents are independent
 * given the base. The surface is a search priority, not evidence that any
 * person lives anywhere, and on its own it can never identify an
 * individual. Accuracy on known solved series (clab_jtacc) is the only
 * guide to how much a ranking is worth in a given place.
 * ------------------------------------------------------------------------ */

/* Distance-decay functions. d: n distances >= 0 (planar units); type and
 * par: 0 linear, f = max(0, A + B d), par {A, B}; 1 negative exponential,
 * f = A exp(-b d), par {A > 0, b >= 0}; 2 normal, f = A phi((d - mu) / s)
 * / s, par {A > 0, mu, s > 0}; 3 lognormal, f = A phi((log d - mu) / s) /
 * (d s), f(0) = 0, par {A > 0, mu, s > 0} (the lognormal density of d;
 * some crime-mapping texts use log d^2 and 1 / d^2 instead); 4 truncated
 * negative exponential, f = P d / dp for d <= dp, P exp(-C (d - dp))
 * beyond, par {P > 0, dp > 0, C >= 0}, continuous with its peak P at dp;
 * 5 a table, par {delta > 0} with tab[0..nt) (nt >= 2, values >= 0) at
 * distances k delta, linear in between and 0 beyond (nt - 1) delta, e.g.
 * from clab_jtcal. tab and nt are read only for type 5. Out: f, n values.
 * The forms follow Brantingham and Brantingham (1981), Rossmo (2000) and
 * Levine (2009, doi:10.1002/jip.107). O(n). */
int clab_jtdecay(const double *d, int n, int type, const double *par,
                 const double *tab, int nt, double *f);

/* The journey-to-crime surface over a grid of candidate bases. x, y: m >=
 * 1 incident locations (planar units); the grid has nx * ny cells of size
 * dx by dy from the corner (x0, y0); cell c = iy * nx + ix has its centre
 * at (x0 + (ix + 1/2) dx, y0 + (iy + 1/2) dy). metric 0 Euclidean, 1
 * Manhattan (|dx| + |dy|, a grid of streets aligned with the axes). type,
 * par, tab, nt as for clab_jtdecay. mode 0: raw = sum_i f(d_ic), the
 * classic aggregate surface (Canter et al. 2000,
 * doi:10.1023/A:1007551316253; Levine 2009); mode 1: raw = sum_i log
 * f(d_ic), the log-likelihood of the base when f is the planar density of
 * an incident about the base and incidents are independent (O'Leary 2009,
 * doi:10.1002/jip.111): with a uniform prior over the cells, prob is then
 * the posterior over cells, and simulation shows its credible regions are
 * calibrated (src/t_jtc.cpp); the sum is not a likelihood. Out: raw, prob
 * (nx ny each; prob sums to 1, computed from raw - max in mode 1 so that it
 * cannot underflow) and peak, 3: the first cell of maximum raw, its
 * centre's x and y. CLAB_EDOM when the surface is 0 (mode 0) or -infinity
 * (mode 1) everywhere. O(m nx ny). */
int clab_jtsurf(const double *x, const double *y, int m, int nx, int ny,
                double x0, double y0, double dx, double dy, int metric,
                int mode, int type, const double *par, const double *tab,
                int nt, double *raw, double *prob, double *peak);

/* Empirical decay: a kernel density estimate of trip lengths from known
 * offender trips (Silverman 1986, doi:10.1201/9781315140919; Levine 2009).
 * d: m >= 1 trip lengths >= 0; w: NULL or m non-negative weights with a
 * positive sum; kern 0 normal, 1 quartic (15/16 (1 - u^2)^2), 2 triangular,
 * 3 uniform, 4 double exponential (exp(-|u|) / 2); h > 0 the bandwidth;
 * refl 1 reflects the kernels at 0 so that no mass falls below 0 (t must
 * then be >= 0), 0 does not. t: k evaluation distances. Out: dens, k values
 * of sum_i w_i [K((t - d_i) / h) + refl K((t + d_i) / h)] / (h sum w): a
 * density in t (it integrates to 1 over the line, or over [0, inf) when
 * reflected). Tabulated at t = k delta it is a type 5 decay function.
 * O(m k). */
int clab_jtcal(const double *d, const double *w, int m, int kern, double h,
               int refl, const double *t, int k, double *dens);

/* Accuracy of a surface for a solved series. prob: nx ny non-negative
 * cell values on the grid of clab_jtsurf; (hx, hy) the known base, inside
 * the grid. Out, 6: prob in the base's cell; the hit score as the share of
 * the grid searched up to and including the base's cell in the best order
 * of ties, (#{prob > p_base} + 1) / (nx ny), and in the worst,
 * #{prob >= p_base} / (nx ny) (Rossmo 2000; Canter et al. 2000); the error
 * distance from the centre of the peak cell (first maximum) to the base;
 * the peak cell; the base's cell. O(nx ny). */
int clab_jtacc(const double *prob, int nx, int ny, double x0, double y0,
               double dx, double dy, double hx, double hy, double *out);

/* Bayesian journey to crime (Levine 2009, doi:10.1002/jip.107; Levine and
 * Block 2011, doi:10.1080/00330124.2010.547152). pj: no >= 1 non-negative
 * journey-to-crime values over origin zones (e.g. prob of clab_jtsurf with
 * the zones as cells); od: no x nd origin-destination counts of solved
 * offences (row o, column j: offences in zone j by offenders based in
 * zone o; od[o nd + j]); dest: the m >= 1 destination zones of the series'
 * incidents (repeats count). Out, 5 no, each block normalised to sum 1:
 * P(Jtc); P(O), the origins of all offences (row sums); P(O|Jtc), the
 * origins of offences in the series' zones (sum over k of od[o,
 * dest[k]]); the product P(Jtc) P(O|Jtc), the numerator of Bayes' rule
 * with the journey-to-crime surface as prior and the origin distribution
 * as likelihood; and the risk, product / P(O) renormalised (0 where P(O) is
 * 0; the product is 0 there too). CLAB_EDOM when a block sums to 0 (e.g. no
 * solved offence in the series' zones). The O-D matrix describes past
 * solved offences and inherits their selection. O(no (nd + m)). */
int clab_bjtc(const double *pj, int no, const double *od, int nd,
              const int *dest, int m, double *out);


/* Fuzzy mode: for each point, the number of points (itself and coincident
 * ones included) within distance r >= 0 of it; r = 0 counts coincident
 * points, the mode of exact locations. Locations are counted once per
 * point, so the counts of nearby points overlap by design. x, y: n >= 1
 * points (planar units). iwork: 2n + 2 ints. Out: cnt, n counts. Grid
 * search: O(n + pairs within r). */
int clab_fzmode(const double *x, const double *y, int n, double r, int *iwork,
                int *cnt);

/* K-means partition of n >= 1 points into k in [1, n] groups by Lloyd's
 * iterations (Lloyd 1982, doi:10.1109/TIT.1982.1056489; MacQueen 1967):
 * assign each point to the nearest centre (lowest index on ties), move
 * each centre to its points' mean (upd 0: K-means, minimising the
 * within-group sum of squares) or to their Weber point (upd 1: the
 * variant in some crime-mapping programs, a spatial K-medians minimising
 * the sum of distances), until no point changes group or after maxit
 * moves. An emptied group keeps its centre. Seeds (init): 0 the k centres
 * given in cen; 1 k-means++ (Arthur and Vassilvitskii 2007, ACM-SIAM SODA
 * 1027-1035), drawn from seed (CLAB_EDOM if fewer than k distinct points);
 * 2 separated cells: square cells of side par[0] > 0 from the lower-left
 * corner of the points' bounding box, taken by decreasing count (then
 * cell number), a cell's centre becoming a seed when it is at least par[1]
 * >= 0 from every earlier seed (a common choice is par[1] = t 0.5
 * sqrt(A / n) with t near 4); CLAB_EDOM if fewer than k cells qualify.
 * Deterministic for given inputs and seed. cen: 2k doubles (x then y), in
 * for init 0, out always. cl: n group numbers. iwork: 2n ints; work: 2n
 * doubles. Out, 4 + 3k: the within-group sum of squares, the sum of
 * distances to the centres, the number of moves, 1 if converged; then for
 * each group its size, sum of squares and sum of distances. O(n k) per
 * iteration (O(n) more per group for upd 1). */
int clab_kmpart(const double *x, const double *y, int n, int k, int init,
                const double *par, unsigned int seed, int upd, int maxit,
                double *cen, int *cl, int *iwork, double *work, double *out);

/* Nearest-neighbour hierarchical clustering (Levine, the reference crime-
 * mapping manual, ch. 7; after Johnson 1967, doi:10.1007/BF02289588). At
 * order 1 the units are the n >= 2 points; two units are linked when their
 * distance is below the threshold: fixd if fixd > 0, otherwise the mean
 * nearest-neighbour distance under complete spatial randomness plus z
 * standard errors, 0.5 sqrt(A / m) + z sqrt((4 - pi) A / (4 pi)) / m (m
 * units, area > 0; z the normal quantile of the one-tailed probability
 * that a pair is that close by chance, e.g. z = 0 for 0.5). Units with a
 * linked neighbour are taken in decreasing number of links (then index);
 * each not yet clustered takes every unclustered unit linked to it; groups
 * with at least minpts (order 1) or minhi (higher orders) members are kept
 * and numbered by size. The kept groups' Weber points are the units of
 * the next order, with the threshold recomputed for their number, while
 * at least 4 units remain, up to maxord orders. The iterative reassignment
 * some descriptions add after the first pass is not done: its rule is not
 * stated unambiguously. cl: n maxord ints, cl[o n + i] the cluster of
 * point i at order o + 1 or -1. Out: *ncl clusters (at most n), cs: 11 per
 * cluster: order, points, member units, mean x, y, Weber x, y, and the
 * ellipse of the member units (as clab_sdell with divisor k in {0, 1, 2}
 * and multiplier c > 0: a, b, angle; NaN if the members are too few for
 * the divisor), and the threshold used. iwork: 6n + 2 ints; work: 4n
 * doubles. Grid search: O(n + linked pairs) per order. */
int clab_nnhclu(const double *x, const double *y, int n, double area,
                double fixd, double z, int minpts, int minhi, int maxord,
                int k, double c, int *iwork, double *work, int *cl, int *ncl,
                double *cs);

/* Search-circle hot-spot areas (Block and Block 1995, in Eck and Weisburd
 * (eds.), Crime and Place, Crime Prevention Studies 4; as described in
 * the reference crime-mapping manual, ch. 8). Circles of radius r > 0 are
 * centred on the nx x ny nodes (nx,
 * ny >= 2) of a rectangular grid spanning win = {x0, x1, y0, y1}, corners
 * included; the points within each are counted (closed disc); the circles
 * with at least 2 points are ranked by count (then node number) and the
 * first top >= 1 are kept; kept circles that share a point are merged,
 * repeatedly, and a hot cluster is the set of points in a merged group,
 * kept if it has at least minpts >= 2 points. Descriptions differ on the
 * circle's radius (the search radius or sqrt 2 times it) and on the grid
 * (20 x 20, rectangular or triangular): both are the caller's choice
 * here. Out: cl, n cluster numbers (by size, then best circle) or -1;
 * *ncl; cs, 8 per cluster: points, mean x, y, Weber x, y, ellipse a, b,
 * angle (as clab_sdell with k, c). iwork: 3n + 2 + nx ny + top + max(n,
 * 2 top) ints; work: nx ny + 2n doubles. O(nx ny log(nx ny) + node
 * neighbourhoods). */
int clab_stacel(const double *x, const double *y, int n, const double *win,
                int nx, int ny, double r, int top, int minpts, int k, double c,
                int *iwork, double *work, int *cl, int *ncl, double *cs);


/* Mantel's space-time test (Mantel 1967, Cancer Research 27: 209-220):
 * the Pearson correlation r over all n (n - 1) / 2 pairs of events
 * between their spatial dissimilarity X (the distance, or 1 / (distance +
 * cs) when cs > 0, Mantel's reciprocal transform) and their temporal one
 * Y (|t_i - t_j|, or 1 / (|t_i - t_j| + ct) when ct > 0). x, y: planar
 * locations, t: times, n in [3, 46340]. The test permutes the times over
 * the fixed locations, nsim >= 0 times from seed (Fisher-Yates); drawing
 * distances and intervals independently instead does not preserve the
 * dependence among the pairs of one distance matrix and gives no test of
 * stated size. Out, 5: r; the number of pairs; the permutation p-values
 * (1 + #{r_perm >= r}) / (nsim + 1), (1 + #{r_perm <= r}) / (nsim + 1)
 * and (1 + #{|r_perm| >= |r|}) / (nsim + 1) (NaN without permutations).
 * CLAB_EDOM when either dissimilarity is constant. iwork: n ints. Memory
 * O(n), no distance matrix; time O(n^2 (nsim + 1)). */
int clab_mantst(const double *x, const double *y, const double *t, int n,
                double cs, double ct, int nsim, unsigned int seed, int *iwork,
                double *out);

/* Spatial-temporal moving average: the events in time order (ties in
 * input order), and for each the mean centre of the span (odd, >= 1)
 * events centred on it, the window shortened at the ends of the series.
 * x, y: n >= 1 planar locations, t: times. iwork: n ints. Out, 4 n, one
 * row per event in time order: its input index, the window's mean x and
 * y, and the number of events in the window. O(n log n + n span). */
int clab_stmavg(const double *x, const double *y, const double *t, int n,
                int span, int *iwork, double *out);
/* Fairness of a binary risk tool by group (confusion-matrix metrics).
 * score: n scores (any value but NaN); a case is predicted positive when
 * score > cut (cut not NaN; +Inf predicts nobody positive, -Inf everybody
 * with a finite score). y: n outcomes, 0 or 1. grp: n group codes in
 * [0, g). base in [0, g): the group the ratios are taken against. z > 0:
 * the normal quantile of the intervals (1.959963984540054 for 95%); read
 * only when lo and hi are both non-NULL.
 * Out (k-th metric of group j at [k * g + j]):
 *   cnt, 5g: tp, fp, tn, fn and the group size, exact counts;
 *   met, 10g: accuracy, sensitivity (TPR), specificity (TNR), FPR, FNR,
 *     PPV, NPV, share predicted positive, number predicted positive, and
 *     the Matthews correlation; NaN where a rate has an empty denominator
 *     or MCC has an empty margin (0 / 0);
 *   rat, 10g: met over the base group's met (IEEE: x / 0 = Inf, 0 / 0 =
 *     NaN);
 *   lo, hi, 8g each (may be NULL): Wilson score intervals of the first
 *     eight metrics (Wilson 1927, doi:10.1080/01621459.1927.10502953), the
 *     inversion of the score test, computed without cancellation; lo = 0
 *     exactly at 0 successes and hi = 1 at all successes; NaN on an empty
 *     denominator.
 * A group predicted all positive or all negative is reported as such (a
 * one-sided prediction is never relabelled). The MCC numerator is formed
 * in 64-bit integers, exact for n < 2^31. Time O(n + g), no workspace. */
int clab_fair(const double *score, const int *y, const int *grp, int n, int g,
              double cut, int base, double z, double *cnt, double *met,
              double *rat, double *lo, double *hi);

/* Differences and ratios between groups of a binary decision
 * (clab_fairdf) and the same over time and place (clab_tmpaud).
 *
 * clab_fairdf: pred, n decisions 0 or 1 (1: predicted positive); y, NULL
 * or n outcomes 0 or 1; grp, n codes in [0, g), g >= 2, each code used.
 * base in [0, g), or -1 for the group with the largest share predicted
 * positive (the first on ties; the reference of the four-fifths rule).
 * z > 0 the normal quantile of the intervals. nboot >= 0 bootstrap
 * replicates from seed. Out (k-th value of group j at [k * g + j]):
 *   cnt, 6g: cases, predicted positive, true positives, of them predicted
 *     positive, true negatives, of them predicted positive (exact counts;
 *     the last four 0 without y);
 *   met, 23g: n, number and share predicted positive with its Wilson
 *     interval (0-4); the parity gap (share minus the base group's) with
 *     Newcombe's hybrid score interval (5-7, doi:10.1002/(SICI)1097-0258(
 *     19980430)17:8<873::AID-SIM779>3.0.CO;2-I); the ratio of shares (the
 *     disparate impact ratio, IEEE division) with the log interval of Katz
 *     et al. 1978 (8-10, doi:10.2307/2530610; NaN when either count is 0);
 *     with y, the true and false positive rates (11, 12), their gaps to
 *     the base group with Newcombe intervals (13-15, 16-18), the average
 *     odds difference (dtpr + dfpr) / 2 with the MOVER interval of Zou and
 *     Donner 2008 (19-21, doi:10.1002/sim.3095) and the larger in size of
 *     the two gaps (22); NaN where undefined, and in the base group's own
 *     intervals; columns 11-22 NaN without y;
 *   sum, 11: the smallest finite ratio of another group (NaN if none), the
 *     parity gap of largest size (signed, the first on ties), the Gini
 *     coefficient of the g shares, the bias amplification score (gap times
 *     Gini, Barman and Barman, arXiv:2603.18987), percentile bootstrap
 *     ends of the Gini coefficient and of the score (NaN when nboot = 0),
 *     the average odds difference and equalized-odds gap of largest size
 *     (NaN without y), and the base group.
 * The bootstrap redraws each group's decisions with replacement (n_j draws,
 * Binomial(n_j, share_j) positives) and keeps the base group. work: 3g +
 * 2 nboot doubles. CLAB_EDOM for a code without cases. O(n + g + n nboot).
 *
 * clab_tmpaud: pred and grp as above; cell, n codes in [0, nc) (one cell
 * per city and period); city, nc codes in [0, ncity), the city of each
 * cell. base in [0, g) or -1 for the group with the largest pooled share.
 * Out: cout, 6 per cell (6 c + k): cases, audited (1) or skipped (0: fewer
 * than two groups present, or the base group absent), the smallest finite
 * ratio, the parity gap of largest size, the Gini coefficient of the
 * present groups' shares and the bias amplification score (NaN when
 * skipped); yout, 9 per city: audited cells, the means of the four values
 * over them (finite ones only), the smallest and largest finite ratio,
 * their range and the number of cells with a ratio above 1; sum, 3: the
 * largest range, the spread of the cities' mean ratios (NaN with fewer
 * than two), the base group. work: 2 nc g + 2g doubles. O(n + nc g). */
int clab_fairdf(const int *pred, const int *y, const int *grp, int n, int g,
                int base, double z, int nboot, unsigned int seed, double *cnt,
                double *met, double *sum, double *work);
int clab_tmpaud(const int *pred, const int *grp, const int *cell, int n, int g,
                int nc, const int *city, int ncity, int base, double *cout,
                double *yout, double *sum, double *work);

/* Individual fairness against a Lipschitz constraint (Dwork et al. 2012,
 * doi:10.1145/2090236.2090255): |h_i - h_j| <= lip d_ij for each audited
 * pair. h: n finite scores; pi, pj: m pairs of indices in [0, n); d: m
 * distances >= 0, or NULL for all 1; lip >= 0. Out: gap (m, |h_i - h_j|),
 * rat (m, gap / d; Inf for a gap at distance 0, 0 for none), viol (m, 1
 * when gap > lip d: at distance 0 any gap), sum, 6: the smallest constant
 * the scores satisfy on these pairs (the largest rat), the violations,
 * their share, the largest gap, the first pair with the largest gap and
 * the first with the largest rat. O(m). */
int clab_lipaud(const double *h, int n, const int *pi, const int *pj,
                const double *d, int m, double lip, double *gap, double *rat,
                int *viol, double *sum);

/* The true base rate behind a noisy recorded one
 * (Lean: Crimlab.Fair.bounds_box). With false-positive noise a = P(recorded
 * 1 | true 0) in [alo, ahi] and false-negative noise b = P(recorded 0 |
 * true 1) in [blo, bhi], 0 <= alo <= ahi, 0 <= blo <= bhi, ahi + bhi < 1,
 * the recorded rate is pobs = p (1 - b) + (1 - p) a, so p = (pobs - a) /
 * (1 - a - b), which falls in a and rises in b. pobs: n rates in [0, 1].
 * Out, n each: lo = (pobs - a*) / (1 - a* - blo) with a* = min(ahi, pobs),
 * hi = (pobs - alo) / (1 - alo - b*) with b* = min(bhi, 1 - pobs): the
 * sharp interval, both ends attained; NaN when no admissible noise gives
 * pobs (pobs < alo or pobs > 1 - blo). Known rates (alo = ahi, blo = bhi)
 * give lo = hi = p. For each pair i < j in the order (0, 1), (0, 2), ...,
 * (1, 2), ...: ord, 1 when hi_i < lo_j (p_i < p_j under every admissible
 * noise, Crimlab.Fair.compare_decided), -1 when hi_j < lo_i, 0 when the
 * intervals overlap (either order is possible, compare_undecided), 2 when
 * either is undefined; bd, the largest bhi (the a box held) at which the
 * lower recorded rate's interval stays below the other's, capped at
 * 1 - ahi, NaN if none or the rates are equal. n (n - 1) / 2 must fit an
 * int. O(n^2). */
int clab_nsbnd(const double *pobs, int n, double alo, double ahi, double blo,
               double bhi, double *lo, double *hi, int *ord, double *bd);

/* Chouldechova's identity (Lean: Crimlab.Fair.chouldechova): the false
 * positive rate fpr = p / (1 - p) (1 - ppv) / ppv (1 - fnr) a base rate p in
 * (0, 1), positive predictive value ppv in (0, 1] and false negative rate
 * fnr in [0, 1) imply, for n triples. */
int clab_impfpr(const double *p, const double *ppv, const double *fnr, int n,
                double *fpr);

/* The resolution of a ranking by noisy scores
 * (Lean: Crimlab.Fair.rank_stable_of_gap): estimates est, n >= 2, each off
 * its true score by at most hw >= 0 (n half-widths). Out: id, n x n
 * (column-major): 1 when |est_i - est_j| > hw_i + hw_j (the order of i and
 * j survives every admissible noise), 0 when not, -1 on the diagonal; out,
 * 4: the pairs, the pairs not identified, their share, and twice the
 * largest half-width. n^2 must fit an int. O(n^2). */
int clab_rnkres(const double *est, const double *hw, int n, int *id,
                double *out);

/* Calibration audit of an area ranking. area: n codes in [0, na), na >= 3,
 * each code used; risk, outc: n finite values (an algorithm's predicted
 * risk and the recorded outcome of each record); grp: n codes in [0, g);
 * pop: NULL or n positive populations, one value per area. Each area gets
 * its record count, mean risk, outcome rate (mean outcome, or summed
 * outcome over population), group (its most frequent, the first code on
 * ties), the average ranks of risk and rate with rank 1 the largest, and
 * the gap rate rank - risk rank (positive: ranked riskier than its outcome
 * warrants). Out: ar, 7 na (k na + a, in that order); gr, 3g: areas, mean
 * gap, two-sided permutation p-value of each group's mean gap (NaN for a
 * group without areas or nperm = 0); st, 7: areas, Spearman's rho (Pearson
 * on the ranks; NaN when either is constant), S = (na^3 - na)(1 - rho) / 6,
 * the two-sided p-value of algorithm AS 89 (Best and Roberts 1975,
 * doi:10.2307/2347111; exact enumeration for na <= 9, the Edgeworth series
 * up to 1289; NaN with ties, beyond, or rho NaN), whether either ranking
 * has ties, T = sum_g n_g mean_g^2 (the between-group sum of squares of the
 * gaps) and its permutation p-value. The permutations shuffle the areas'
 * groups nperm times from seed; p = (1 + #{at least as extreme}) / (1 +
 * nperm). CLAB_EDOM when a population varies within an area. iwork: n +
 * 2 na ints; work: na + 3g doubles. O(n log n + nperm (na + g)). */
int clab_calaud(const int *area, const double *risk, const double *outc,
                const int *grp, const double *pop, int n, int na, int g,
                int nperm, unsigned int seed, double *ar, double *gr,
                double *st, int *iwork, double *work);

/* A score's spread across groups. x: n finite scores; grp: n codes in
 * [0, g), g >= 2, each used; ref in [0, g), or -1 for the group with the
 * lowest mean. Out: gs, 6g (k g + j): cases, mean, median, standard
 * deviation (NaN for one case), the gap mean - ref's mean, and its standard
 * error from the pooled within-group variance (NaN for ref); st, 9: the
 * one-way analysis of variance F with its two degrees of freedom (F NaN
 * when n = g), Welch's F with its degrees of freedom (Welch 1951,
 * doi:10.1093/biomet/38.3-4.330; NaN unless every group has two or more
 * cases and positive variance), the spread of the means (largest minus
 * smallest), ref, and the pooled within-group variance. iwork: n ints.
 * O(n log n). */
int clab_scrdsp(const double *x, const int *grp, int n, int g, int ref,
                double *gs, double *st, int *iwork);

/* Noisy-OR detection (Barman and Barman, arXiv:2603.18987): an incident at
 * (cx, cy), n of them, is detected by each of the k officers within radius
 * > 0 (Euclidean, inclusive) of it independently with probability p in
 * (0, 1], so with probability 1 - (1 - p)^k. Officers (ox, oy), m >= 0.
 * Out, n each: k, prob, and, when draw != 0, det (1 with probability prob,
 * from seed). O(n m). */
int clab_nordet(const double *cx, const double *cy, int n, const double *ox,
                const double *oy, int m, double radius, double p, int draw,
                unsigned int seed, int *k, double *prob, int *det);

/* Records with a known disparity (after Barman and Barman,
 * arXiv:2603.18987). n records; group j in [0, g), g >= 2, drawn by props
 * (g weights >= 0 with a positive sum, or NULL for equal shares); area in
 * [0, nareas), nareas >= g, uniform among the areas a with a mod g = j
 * (segregated areas); truth, 1 with probability base in [0, 1] whatever
 * the group; det, the recorded outcome, 1 with probability base for group
 * 0 and base (1 - bias), clamped to [0, 1], for the others, bias in
 * [-1, 1]; risk, normal with mean 250 (group 0) or 250 + 100 bias and sd
 * 40, clamped to [0, 500]. One draw sequence from seed (xorshift64*). */
int clab_fairsim(int n, int g, const double *props, int nareas, double base,
                 double bias, unsigned int seed, int *area, int *grp,
                 int *truth, int *det, double *risk);

/* The Oaxaca-Blinder decomposition (Lean: Crimlab.Oaxaca.threefold) of the
 * gap in mean outcome y between group A (code 1 - ref) and group B (code
 * ref). y: n finite outcomes; x: n x p finite covariates (column-major, a
 * constant is added); grp: n codes 0 or 1, each group with more than p + 1
 * cases. Least squares in each group (Householder QR; CLAB_EDOM when a
 * covariate is collinear with the others there). With D = xbar_A - xbar_B
 * and Delta = b_A - b_B (the constant first): out xm, 2(p + 1): the means
 * of A then B; b, bse, 2(p + 1): the coefficients and their standard
 * errors; dec, 14: ybar_A, ybar_B, the gap and its standard error, then
 * with its standard error each: D'b_B (explained at B's coefficients),
 * xbar_A'Delta (unexplained), D'b_A, xbar_B'Delta, D'Delta (the
 * interaction); the threefold decomposition is D'b_B + xbar_B'Delta +
 * D'Delta. var, 5(p + 1): for each coefficient D_j b_Bj and its standard
 * error, xbar_Aj Delta_j and its, and Delta_j (the change in that line of
 * the unexplained part when covariate j moves by 1 for everyone). Standard
 * errors by the delta method of Jann (2008,
 * doi:10.1177/1536867X0800800401): the groups, and means and coefficients
 * within one, independent; the second-order term left out. work:
 * (n + 5p + 6)(p + 1) + n doubles. O(n p^2). */
int clab_oaxbd(const double *y, const double *x, const int *grp, int n, int p,
               int ref, double *xm, double *b, double *bse, double *dec,
               double *var, double *work);

/* Two-level normal model for multivariate trace evidence (Aitken and Lucy
 * 2004, doi:10.1046/j.0035-9254.2003.05271.x): x_ij = theta_i + e_ij,
 * theta_i ~ N(mu, C) between items, e_ij ~ N(0, U) within an item.
 *
 * clab_tlvc: the background estimates. x: n x p observations, column-major;
 * item: n codes in [0, m), every item observed at least once; m >= 2 and
 * n > m (some item replicated). Out: mu (p), the mean of all n
 * observations; u (p x p), the pooled within-item covariance
 * sum (x - xbar_i)(x - xbar_i)' / (n - m); c (p x p), the ANOVA between-item
 * covariance (S_B - (m - 1) U) / (n - sum n_i^2 / n) with
 * S_B = sum n_i (xbar_i - mu)(xbar_i - mu)' (Searle, Casella and McCulloch
 * 1992, doi:10.1002/9780470316856). Both are unbiased for any replication
 * pattern; with equal n_i they are the balanced estimates of Aitken and
 * Lucy. c need not be positive definite (clab_tlvlr then refuses it).
 * Means take one corrective pass; cross-product sums are compensated.
 * work: m p + 2 m doubles. Time O(n p^2 + m p^2). */
int clab_tlvc(const double *x, int n, int p, const int *item, int m,
              double *mu, double *u, double *c, double *work);

/* clab_tlvlr: the likelihood ratio of same source against different
 * sources for two items with mean vectors xa, xb of na, nb > 0 replicates,
 * given u, c (p x p, lower triangles read) and mu (p). With
 * d = xa - xb and ybar = (na xa + nb xb) / (na + nb),
 *   ln LR = ln N(d; 0, U (1/na + 1/nb)) + ln N(ybar; mu, U/(na+nb) + C)
 *         - ln N(xa; mu, U/na + C) - ln N(xb; mu, U/nb + C),
 * evaluated with four Cholesky factors in log space, so atypical or far
 * apart items give a finite value instead of 0 / 0. Out, 2: ln LR and
 * log10 LR. CLAB_EDOM when one of the four matrices is not positive
 * definite. work: p p + 2 p doubles. Time O(p^3). */
int clab_tlvlr(const double *u, const double *c, const double *mu, int p,
               const double *xa, double na, const double *xb, double nb,
               double *out, double *work);

/* Aoristic weights (Ratcliffe 2002, doi:10.1023/A:1013240828824): an event
 * known only to lie in the real-time interval [from, to) is spread with
 * uniform density over that interval, and each instant is counted in the
 * bin of its local wall-clock time. Times are seconds (UTC, e.g. POSIX
 * time). The local clock is UTC + off0 before tt[0] and UTC + toff[k] from
 * tt[k] (nt strictly increasing transitions; nt = 0 for a fixed offset), so
 * a daylight-saving change moves the remaining time to the right wall-clock
 * hours: hours skipped in spring get nothing and the repeated hour in
 * autumn gets both. Bins: nb periodic bins of bw > 0 seconds, bin 0
 * starting at wall-clock time origin (e.g. 259200, Sunday 1970-01-04
 * 00:00, with bw = 3600 and nb = 168 for the hour of the week).
 * Rules: from NaN: the event is left out; to NaN, to < from, or a length
 * of 0 or below lo: the whole weight in the start's bin; a length of at
 * least hi: 1 / nb in every bin (lo = 120, hi = 604800 reproduce the usual
 * convention for whole-minute data; lo = 0, hi = Inf spread every interval
 * exactly). Infinite times are refused.
 * Out: w (n x nb column-major, may be NULL): each event's weights, summing
 * to 1; tot (nb, may be NULL): their sums over events; cnt (5): events left
 * out, without end, reversed, short, long. n = 0 is allowed.
 * Time O(n log(nt + 1) + sum of bins touched) <= O(n (log(nt + 1) + nt +
 * nb)); no workspace. */
int clab_aorwt(const double *from, const double *to, int n, const double *tt,
             const double *toff, int nt, double off0, double origin,
             double bw, int nb, double lo, double hi, double *w, double *tot,
             int *cnt);

/* Benford digit screens.
 *
 * clab_bdig: the first nd significant digits (1 <= nd <= 9) of each value,
 * as an integer in [10^(nd-1), 10^nd), padded with zeros (0.07 -> 70 for
 * nd = 2), or 0 when the value is set aside. sgn: 1 keeps positive values
 * only, -1 negative only, 0 both (by absolute value). dp in 0..22: the
 * data are recorded to dp decimals and the digits are those of the
 * integer round(|x| 10^dp), exact for every double nearest to such a
 * decimal while |x| 10^dp < 2^49; dp = -1: the digits of the shortest
 * decimal that rounds to x, decided exactly by comparing x with the
 * correctly rounded digit boundaries D 10^q (|q| <= 22, so |x| from about
 * 10^(nd - 23) to 10^(nd + 21)), so 0.29 gives 29 where truncating
 * 10^(1 - floor(log10 x)) x gives 28.
 * cnt, 4: zeros (also values that round to 0 at dp decimals), non-finite
 * values, values of the excluded sign, values out of the exact range.
 * Time O(n) (at most 9 boundary pairs a value); no workspace.
 *
 * clab_bnft: goodness of fit of the digits dig (n values from clab_bdig
 * with the same nd, 1 <= nd <= 6; zeros are skipped) to Benford's law
 * P(d) = log10(1 + 1/d) (Benford 1938; Nigrini 2012,
 * doi:10.1002/9781119203094). Out: obs (9 10^(nd-1) counts) and out, 6:
 * the number of digits used m, the chi-square statistic, its degrees of
 * freedom, the upper tail by the closed forms for integer degrees of
 * freedom (finite sums of Poisson-type terms in logs, no series
 * truncation), a Monte Carlo p-value (1 + #{chi* >= chi}) / (nsim + 1)
 * from nsim samples of m Benford digits (North et al. 2002,
 * doi:10.1086/341527; NaN when nsim = 0; seed gives the stream, xorshift64*
 * seeded by splitmix64), and the mean absolute deviation of the
 * proportions (Nigrini's MAD). m = 0 gives NaN statistics. work: 9
 * 10^(nd-1) doubles. Time O(n + 10^nd + nsim (m + 10^nd)). */
int clab_bdig(const double *x, int n, int nd, int dp, int sgn, int *dig,
              int *cnt);
int clab_bnft(const int *dig, int n, int nd, int nsim, unsigned int seed,
              double *obs, double *work, double *out);

/* Scott-Knott grouping of measurements (Scott and Knott 1974,
 * doi:10.2307/2529204), with the measurement variance known as in the
 * forensic glass version (Triggs, Curran, Buckleton and Walsh 1997,
 * doi:10.1016/S0379-0738(96)02037-3). x: n finite values (n >= 1; NaN is
 * refused, never dropped, so labels stay aligned with the input); s2 > 0:
 * the variance of one measurement, on the scale of x; alpha in (0, 1).
 * A segment of the sorted values is split at the cut with the largest
 * between-groups sum of squares B0 when
 * lambda = pi / (2 (pi - 2)) B0 / s2 is significant: the p-value
 * (1 + #{lambda* >= lambda}) / (nsim + 1) from nsim >= 1 samples of the
 * same size of independent N(0, s2) values (exact size under that null;
 * stops early once p must exceed alpha); seed fixes the stream. Segments
 * are examined with an explicit stack, no recursion. Out: grp, n group
 * labels 1..ng in increasing order of the values; ng. Ties are ordered by
 * input position. iwork: 3n ints; work: 2n doubles. Time O(n log n) for
 * the sort plus O(nsim k log k) per segment of size k examined. */
int clab_skp(const double *x, int n, double s2, double alpha, int nsim,
             unsigned int seed, int *grp, int *ng, int *iwork, double *work);

/* Constant reporters of a reporting panel (for example agencies in UCR
 * data): a unit is constant when, in every year y0..y1, at least one of its
 * records has rep >= minrep months reported. A year without a record of the
 * unit, or with too few months, disqualifies it, whatever happens to other
 * units in that year. unit: n codes in [0, a); year, rep: n integers
 * (records outside y0..y1 are ignored; duplicates are allowed). Out: flag,
 * a values (1 constant, 0 not) and nc, their sum. work: a (y1 - y0 + 1)
 * ints. Time O(n + a (y1 - y0 + 1)). Monotone: a lower minrep never
 * removes a unit. */
int clab_rptc(const int *unit, const int *year, const int *rep, int n, int a,
              int y0, int y1, int minrep, int *flag, int *nc, int *work);

/* Rosenbaum sensitivity bounds for matched pairs (Rosenbaum 1987,
 * doi:10.1093/biomet/74.1.13). d: n finite treated-minus-control
 * differences; gam: ng values Gamma >= 1, the odds of differential
 * assignment.
 *
 * clab_rswp: the Wilcoxon signed-rank statistic T (zeros dropped, midranks
 * for tied |d|) and, for each Gamma, the normal-approximation bounds on its
 * one-sided p-value: with p = Gamma / (1 + Gamma) (upper) or 1 / (1 + Gamma)
 * (lower), 1 - Phi((T - p S1) / sqrt(p (1 - p) S2)), S1 and S2 the sums of
 * the ranks and of their squares; tails by erfc, so small p-values keep
 * their digits. gcrit: the Gamma at which the upper bound equals
 * 1 - Phi(z) (z > 0), in closed form from the quadratic in p (below 1 when
 * the test is not significant without bias). Out: plo, pup (ng each),
 * gcrit, tsum (T). NaN when every difference is 0. work: n doubles; iwork:
 * n ints. Time O(n log n + ng).
 *
 * clab_rshl: the bounds on the Hodges-Lehmann estimate under Gamma
 * (Rosenbaum 1993, doi:10.1080/01621459.1993.10476405), exactly. With
 * M = n (n + 1) / 2 and T(tau) = #{Walsh averages (d_i + d_j) / 2, i <= j,
 * above tau} (the signed-rank statistic of d - tau where no difference is
 * tied or zero), the estimate at mu is (sup{tau: T > mu} +
 * inf{tau: T < mu}) / 2 = (W_(ceil a) + W_(floor a + 1)) / 2, a = M - mu;
 * the upper bound takes mu = M / (1 + Gamma), the lower mu = Gamma M /
 * (1 + Gamma), and Gamma = 1 gives the usual estimate. The ceilings are
 * decided exactly (fma) and each order statistic of the Walsh averages is
 * found by bisection on the doubles with an O(n) count, so no grid and
 * no tolerance. n <= 2^26. Out: lo, hi (ng each). work: n doubles. Time
 * O(n log n + ng n 64). */
int clab_rswp(const double *d, int n, const double *gam, int ng, double z,
              double *plo, double *pup, double *gcrit, double *tsum,
              double *work, int *iwork);
int clab_rshl(const double *d, int n, const double *gam, int ng, double *lo,
              double *hi, double *work);

/* Group-based trajectory model for counts (Nagin 1999,
 * doi:10.1037/1082-989X.4.2.139; zero-inflated Poisson within groups as in
 * Nielsen et al. 2014, doi:10.1080/03610926.2012.719986, and Lambert 1992,
 * doi:10.2307/1269547). Unit i
 * belongs to group k with probability prob[k]; given k, its counts y_it
 * are independent, Poisson with log mean off_it + x_t' beta_k, and, when
 * q > 0, zero with extra probability logistic(z_t' gamma_k).
 * y: n x nt counts, column-major; NaN or a negative value is missing and
 * contributes nothing; other values must be whole numbers. off: n x nt log
 * exposures or NULL. x: nt x p (p >= 1), z: nt x q (q >= 0; 0 gives a
 * Poisson model) designs, column-major, shared by all units (for example
 * B-spline bases of time). par: ng blocks of (beta (p), gamma (q)); prob:
 * ng positive shares summing to 1.
 *
 * clab_zgll: the log-likelihood sum_i log sum_k prob_k prod_t f_k(y_it),
 * with each cell in closed form, log(e^a + e^-lam) and log(1 + e^a)
 * without overflow, and log-sum-exp over groups. Out: w (n x ng posterior
 * memberships), ll. work: 2 nt doubles. Time O(n nt ng + nt ng (p + q)).
 *
 * clab_zgem: EM from the given par and prob (updated in place): the
 * E-step above; the M-step sets prob to the mean memberships and raises
 * each group's expected complete-data log-likelihood by up to four damped
 * Newton steps on (beta, gamma) with the analytic Hessian, accepting a step
 * only if it increases (a generalised EM), so the log-likelihood never
 * decreases. The derivatives are summed over units period by period and
 * then expanded over the design, since the design rows depend on the
 * period only. Stops when an
 * iteration gains at most tol (1 + |ll|) or after maxit iterations. Out:
 * par, prob, w, ll, iter. One start per call; several starts are the
 * caller's loop. work: 2 (p + q)^2 + 3 (p + q) + 7 nt doubles. Time per
 * iteration O(ng (n nt + nt (p + q)^2 + (p + q)^3)) times the Newton
 * steps. */
int clab_zgll(const double *y, int n, int nt, const double *off,
               const double *x, int p, const double *z, int q, int ng,
               const double *par, const double *prob, double *w, double *ll,
               double *work);
int clab_zgem(const double *y, int n, int nt, const double *off,
               const double *x, int p, const double *z, int q, int ng,
               double *par, double *prob, double *w, double *ll, int maxit,
               double tol, int *iter, double *work);

/* Group-based multivariate trajectories (Nagin 1999,
 * doi:10.1037/1082-989X.4.2.139; several indicators per unit as in
 * Magrini 2022, doi:10.1007/s10182-022-00437-9): unit i belongs to group k with probability prob[k]; given
 * k, its r indicators at period t are multivariate normal with mean
 * B_k' x_t and covariance S_k, independently over periods.
 * y: n x nt x r values, y[i + n (t + nt j)]; a period with any NaN
 * indicator is left out for that unit; infinite values are refused. x:
 * nt x p design (for example 1, t, t^2), column-major. beta: ng blocks of
 * B_k (p x r, column-major); sig: ng blocks of S_k (r x r, positive
 * definite); prob: ng positive shares summing to 1.
 *
 * clab_mvtll: the log-likelihood (log densities summed per unit, then
 * log-sum-exp over groups: no clamping, so the result does not depend on
 * the measurement units beyond the Jacobian) and the posterior memberships
 * w (n x ng). CLAB_EDOM when a covariance is not positive definite. work:
 * r^2 + r doubles. Time O(ng (r^3 + n nt r (r + p))).
 *
 * clab_mvtem: EM from the given values (updated in place). The M-step is
 * closed form: weighted least squares with one p x p Cholesky for all r
 * indicators (the design is shared), the weighted residual covariance, and
 * the mean memberships; the log-likelihood never decreases. Stops when an
 * iteration gains at most tol (1 + |ll|) or after maxit. CLAB_EDOM when a
 * weighted design is singular or a residual variance is zero to rounding
 * (its sd below 64 ulps of the indicator's root mean square, e.g. a
 * constant indicator). work: p^2 + p r + r^2 + r
 * doubles. Time per iteration O(ng (n nt r (r + p) + nt p (p + r) + p^3 +
 * r^3)). */
int clab_mvtll(const double *y, int n, int nt, int r, const double *x, int p,
               int ng, const double *beta, const double *sig,
               const double *prob, double *w, double *ll, double *work);
int clab_mvtem(const double *y, int n, int nt, int r, const double *x, int p,
               int ng, double *beta, double *sig, double *prob, double *w,
               double *ll, int maxit, double tol, int *iter, double *work);

/* Glass fragments recovered under the transfer, persistence and recovery
 * model of Curran, Triggs, Buckleton, Walsh and Hicks (1998,
 * doi:10.1016/S1355-0306(98)72068-4), in the parameterisation of its
 * simulation software. par, 17 values: lam, lsd (transfer amount
 * X ~ N(lam, lsd^2), used as |X|), d (distance: D ~ Gamma(d, 1) and the
 * intensity is |X| e^{1 - D/d} when deffect = 1, |X| when 0), Q (share of
 * high-persistence fragments), l0, u0 (first-hour loss of the others,
 * uniform), ls0, us0 (first-hour loss of the high-persistence ones), lj, uj
 * and lsj, usj (each later hour), lR, uR (share recovered in the
 * laboratory), lt, ut (time since the offence: fixed when equal, otherwise
 * by tdist: 0 negative binomial with size (lt + ut) / 2 and probability r,
 * 1 that law restricted to floor(lt)..ceil(ut), 2 uniform on those
 * integers), r. A time t gives max(t, 1) hours of loss.
 * Given X, D, the hour count and the loss and recovery rates, the count is
 * Poisson with mean |X| e^{1-D/d} R ((1-Q) prod (1 - a_j) + Q prod
 * (1 - b_j)) (binomial thinning of a Poisson, one rate per case and hour
 * shared by all its fragments). clab_tfrp averages that Poisson law over
 * the intensity M = |X| e^{1-D/d} deterministically: L = log M has a
 * smooth density (the folded normal on the log scale, convolved with the
 * distance term, the latter by a composite 8-point Gauss-Legendre rule
 * with nq panels in w = D^d), tabulated on a grid of step hl and summed by
 * the trapezoid rule, which converges geometrically here; intensities
 * below 1e-17 only add to zero counts, through P(0) = 1 + E[e^-V - 1]
 * (lsd = 0: the Gauss-Legendre rule in w, or one node); over the time
 * exactly (the negative
 * binomial is cut where less than 1e-16 is left); over the rates by nsim
 * Monte Carlo draws (xorshift64* from seed, drawn in a fixed order:
 * recovery, then each hour's two rates), or exactly by one draw when every
 * rate interval is a point. Out: prob, se (ymax + 1 each: the
 * probabilities of 0..ymax fragments and their Monte Carlo standard
 * errors, 0 when exact) and cut, 2: the time mass left out (the rest is
 * renormalised) and the intensity mass on the grid (the rest lies below
 * 1e-17). 1e-4 <= hl <= 1, nq >= 1; hl = 0.02 and nq = 256 resolve the
 * model's defaults (accuracy measured in src/t_tfr.cpp). clab_tfrw gives
 * the workspace length (doubles) for the same inputs. Time
 * O(G 8 nq + nsim H G v), G the grid points (about 47 / hl), H the hours
 * considered, v the typical count range of a node. */
int clab_tfrw(const double *par, int deffect, int tdist, int ymax, double hl,
              int nq, int *len);
int clab_tfrp(const double *par, int deffect, int tdist, int ymax, int nsim,
              unsigned int seed, double hl, int nq, double *prob, double *se,
              double *cut, double *work);

/* Lowess: robust locally weighted regression (Cleveland 1979,
 * doi:10.1080/01621459.1979.10481038), with the conventions of the classic
 * implementation: each fit uses the ceil-free neighbourhood of
 * max(2, min(n, floor(f n + 1e-7))) nearest points (and points tied at its
 * edge), tricube weights, a weighted linear fit (a weighted mean when the
 * points' spread is below 0.001 of the range of x), points within delta of
 * the last fitted one interpolated linearly (equal x share the fit), and
 * iter robustness steps with bisquare weights on six times the median
 * absolute residual (stopping early when that is below 1e-7 of the mean
 * absolute residual). x: n values in non-decreasing order (the caller
 * sorts), y: n values, all finite; f > 0; delta >= 0. Out: ys, n fitted
 * values. work: 3n doubles. Time O(n k) for k neighbourhood points, less
 * with delta. */
int clab_lwss(const double *x, const double *y, int n, double f, int iter,
              double delta, double *ys, double *work);

/* The Chumbley score for two striated tool-mark profiles a (na values) and
 * b (nb), deterministic version (Chumbley et al. 2010,
 * doi:10.1111/j.1556-4029.2010.01424.x; Hadler and Morris 2018,
 * doi:10.1111/1556-4029.13640). Each profile is trimmed to positions
 * round(0.01 n)..round(0.99 n) and detrended by its lowess smooth
 * (clab_lwss, f = coarse, 3 robustness steps, delta 0.01 (m - 1)). The
 * optimisation step finds the pair of windows of wopt points with the
 * largest Pearson correlation (first in column-major order of (window of
 * b, window of a) among ties); the validation step takes windows of wval
 * points at the same shift as that pair, stepping wval along the
 * diagonal on both sides, and at different shifts along the anti-diagonal,
 * as in the published deterministic algorithm, and compares the two sets
 * of correlations by the standardised Wilcoxon rank sum (midranks,
 * tie-corrected variance); p is its upper normal tail (by erfc).
 * Out: nsame, ndiff (window pairs used), u, p (NaN when either set is
 * empty). 2 <= wval <= wopt; coarse > 0. CLAB_EARG when a trimmed profile
 * is shorter than wopt; CLAB_EDOM for a non-finite value or when every
 * window is flat (spread below 64 ulps of the profile's magnitude). The optimisation step is O(na nb) by sliding cross
 * products (no window matrix); the validation step O(n). work:
 * 3 na + 3 nb + 7 max(na, nb) doubles. */
int clab_chmb(const double *a, int na, const double *b, int nb, int wopt,
              int wval, double coarse, int *nsame, int *ndiff, double *u,
              double *p, double *work);

/* Local Getis-Ord G* on a regular grid of nr x nc cells (column-major,
 * cell (r, c) at r + nr c). cnt: counts (whole numbers >= 0) per cell; ex:
 * the exposure of each cell (for points under complete spatial
 * randomness: its area inside the study window, e.g. the share of the
 * cell inside), 0 for a cell outside the window (left out, NaN results).
 * Neighbourhood: the cells whose centres lie within rad cell widths, the
 * cell itself included (rad = 1.5: the 3 x 3 block). nsim >= 0
 * permutations, seed for them. Out, 7 nr nc values, column by column:
 *   gz: the classic G* (Getis and Ord 1992, doi:10.1111/
 *       j.1538-4632.1992.tb00261.x; Ord and Getis 1995, doi:10.1111/
 *       j.1538-4632.1995.tb00912.x) with binary weights: (G - W xbar) /
 *       (s sqrt((m W - W^2) / (m - 1))), xbar and s over the m cells in the
 *       window (exposure ignored);
 *   pz: its two-sided normal p-value; qz: Benjamini-Hochberg q-values of
 *       pz (Benjamini and Hochberg 1995, doi:10.1111/
 *       j.2517-6161.1995.tb02031.x);
 *   pp: the conditional permutation p-value of the classic G* (the cell
 *       keeps its count, the other counts are permuted; Anselin 1995,
 *       doi:10.1111/j.1538-4632.1995.tb00338.x), two-sided,
 *       (1 + #{|G* - W xbar| >= observed}) / (nsim + 1); NaN when nsim = 0;
 *   ze: the exposure-aware z: given the total n, the count G of the
 *       neighbourhood is Binomial(n, P) under complete spatial randomness,
 *       P = the neighbourhood's share of the exposure: (G - n P) /
 *       sqrt(n P (1 - P));
 *   pe: its exact two-sided p-value (twice the smaller binomial tail,
 *       capped at 1), so cells overhanging the window neither inflate the
 *       global mean nor get flagged; qe: Benjamini-Hochberg q-values of pe.
 * No threshold is applied. work: nr nc doubles; iwork:
 * 2 (2 floor(rad) + 1)^2 + 3 nr nc ints. Time O(nr nc (k + tail) +
 * nsim nr nc k), k the neighbourhood size; no distance matrix. */
int clab_gstr(const double *cnt, const double *ex, int nr, int nc, double rad,
              int nsim, unsigned int seed, double *out, double *work,
              int *iwork);

/* Andresen's spatial point pattern test (Andresen 2009,
 * doi:10.1016/j.apgeog.2008.12.004), streamed. a, b: n units' counts of a
 * base and a test pattern (whole numbers >= 0; a unit with a NaN count in
 * either is left out entirely and counted, never treated as zero). Per
 * unit: the difference of shares (b_i / N_b - a_i / N_a) times scale (100:
 * percentage points), and its bootstrap p-value: nboot multinomial
 * resamples of each pattern from its own shares, differences centred on
 * their replicate mean, p = (1 + #{|centred| >= |observed|}) / (nboot + 1);
 * units empty in both patterns get no p. The replicates are generated
 * twice from the seed (xorshift64* seeded by splitmix64, draws by Walker
 * alias tables, all base draws then all test draws per replicate), once
 * for the means and once for the counts, so memory is O(n), not O(n
 * nboot). Out, 3n: differences, p-values, signed flags (p <= alpha: the
 * sign of the difference; 0 otherwise; NaN where there is no p). glob, 8:
 * units left out (NaN), units empty in both, eligible units, flagged
 * units, the S index 1 - flagged / eligible, its z against 1 - alpha with
 * sd sqrt(alpha (1 - alpha) / eligible), the exact sameness p-value
 * P(Binomial(eligible, alpha) >= flagged) (by tail sums), and nboot; NaN
 * when either pattern is empty (no inference). work: 6n doubles; iwork:
 * 4n ints. Time O(nboot (N_a + N_b + n)). */
int clab_sptst(const double *a, const double *b, int n, int nboot,
               unsigned int seed, double alpha, double scale, double *out,
               double *glob, double *work, int *iwork);

/* Networks (src/net.cpp). A graph is held in compressed sparse rows: n >= 1
 * nodes; node i's neighbours are adj[off[i]] .. adj[off[i + 1] - 1] with
 * weights wt at the same positions (wt NULL: every weight 1; given: positive
 * and finite). An undirected graph holds each edge in both rows; a directed
 * one in the row of its tail (out-neighbours), and its in-neighbours in a
 * second CSR built from the reversed edges. Results are by node index. */

/* The CSR of an edge list a[k] - b[k], k < m, over nodes 0 .. n - 1, with
 * weights w (NULL: 1 each). Self-loops are dropped and counted (nloop);
 * repeated edges (the same unordered pair when dir is 0, the same ordered
 * pair otherwise) are merged into one whose weight is their sum, counted in
 * ndup. Each row lists neighbours in the order their edges first appear.
 * off: n + 1; adj, wt: m entries when dir, 2m otherwise; iwork: n ints.
 * O(n + m). */
int clab_ntcsr(const int *a, const int *b, const double *w, int m, int n,
               int dir, int *off, int *adj, double *wt, int *iwork,
               int *nloop, int *ndup);

/* Degree (the row length) and strength (the row's summed weight). */
int clab_ntdeg(int n, const int *off, const double *wt, int *deg,
               double *str);

/* Betweenness by Brandes' algorithm (Brandes 2001,
 * doi:10.1080/0022250X.2001.9990249) on hop counts: the sum over ordered
 * pairs s != t of others of the share of shortest s-t paths through the
 * node, halved when dir is 0. norm divides by (n - 1)(n - 2) (directed) or
 * (n - 1)(n - 2) / 2 (undirected); NaN when n < 3. work: 2n doubles;
 * iwork: 2n ints. O(n m). */
int clab_ntbtw(int n, const int *off, const int *adj, int dir, int norm,
               double *bc, double *work, int *iwork);

/* Closeness on hop counts to the nodes each node reaches along its rows:
 * reach (how many), cls (1 / summed distance), ncls (reach / summed
 * distance), both NaN when nothing is reached, and hrm (the harmonic
 * centrality, the sum of 1 / distance; 0 when nothing is reached).
 * iwork: 2n ints. O(n m). */
int clab_ntcls(int n, const int *off, const int *adj, double *cls,
               double *ncls, double *hrm, int *reach, int *iwork);

/* Eigenvector centrality of an undirected graph by power iteration on
 * A + c I, c the largest weight (so that a bipartite component converges),
 * each connected component scaled to a maximum of 1 at every step, from
 * x = 1; stops when no entry changes by tol or more, or after maxit steps
 * (conv 0). val: the largest of the components' eigenvalues (Rayleigh
 * quotients); a component whose own leading eigenvalue is below it by more
 * than 1e-10 of it gets 0. An edgeless graph gives 1 everywhere and val 0.
 * work: 4n doubles; iwork: 2n ints. O(maxit m). */
int clab_ntevc(int n, const int *off, const int *adj, const double *wt,
               double tol, int maxit, double *x, double *val, int *iter,
               int *conv, double *work, int *iwork);

/* Coreness of each node of an undirected graph (Batagelj and Zaversnik
 * 2003, arXiv:cs/0310049). iwork: 3n ints. O(n + m). */
int clab_ntcore(int n, const int *off, const int *adj, int *core,
                int *iwork);

/* Connected components following the rows of off and, when ioff is not
 * NULL, those of ioff too (the weak components of a directed graph).
 * Components are numbered 0, 1, ... in the order of their first node; nc
 * counts them. iwork: n ints. O(n + m). */
int clab_ntwcc(int n, const int *off, const int *adj, const int *ioff,
               const int *iadj, int *comp, int *nc, int *iwork);

/* Strong components of a directed graph (Tarjan 1972,
 * doi:10.1137/0201010), without recursion, numbered as clab_ntwcc
 * numbers them. iwork: 6n ints. O(n + m). */
int clab_ntscc(int n, const int *off, const int *adj, int *comp, int *nc,
               int *iwork);

/* Modularity of the partition memb (codes 0 .. n - 1) of an undirected
 * graph at resolution res >= 0 (Newman and Girvan 2004,
 * doi:10.1103/PhysRevE.69.026113): the sum over communities of
 * in / 2m - res (tot / 2m)^2, in the weight inside (both directions), tot
 * the summed strength, 2m the total. NaN for a graph without edges.
 * work: 2n doubles. */
int clab_ntmod(int n, const int *off, const int *adj, const double *wt,
               const int *memb, double res, double *q, double *work);

/* The Louvain method (Blondel et al. 2008,
 * doi:10.1088/1742-5468/2008/10/P10008) at resolution res > 0. Each level
 * visits its nodes in node order (useseed 0) or in a permutation drawn from
 * seed (xorshift64*, one stream for the whole run) and moves each to the
 * neighbouring community of largest modularity gain (the first met, in row
 * order, among equal gains) when the gain beats staying by more than 1e-12
 * of the node's strength; passes repeat until none moves a node. A level
 * that moves nothing ends the run; otherwise its communities, numbered by
 * their first node, become the nodes of the next level. Out: memb
 * (communities numbered by first node), q (clab_ntmod of memb) and lev
 * (levels that moved a node). clab_ntlvnw gives the lengths of work
 * (doubles) and iwork (ints) for nnz = off[n]. */
int clab_ntlvnw(int n, int nnz, double *len, double *ilen);
int clab_ntlvn(int n, const int *off, const int *adj, const double *wt,
               double res, unsigned int seed, int useseed, int *memb,
               double *q, int *lev, double *work, int *iwork);

/* Gould-Fernandez brokerage (Gould and Fernandez 1989, doi:10.2307/270949):
 * for each node v, the two-paths i -> v -> j with i != j and no tie
 * i -> j, by role: cnt[v] coordinator (i, v, j in one group), cnt[n + v]
 * itinerant (i and j in one group, v in another), cnt[2n + v]
 * representative (i with v, j apart), cnt[3n + v] gatekeeper (v with j, i
 * apart), cnt[4n + v] liaison (three groups). off/adj: out-neighbours,
 * ioff/iadj: in-neighbours (the same arrays for an undirected graph); grp:
 * any int per node. iwork: n ints. O(sum over v of indeg(v) (outdeg(i) +
 * outdeg(v))). */
int clab_ntbrk(int n, const int *off, const int *adj, const int *ioff,
               const int *iadj, const int *grp, double *cnt, int *iwork);

/* The Cheeger bounds of a connected undirected graph with n >= 2 nodes
 * (Lean: Crimlab.Network.cheeger_lower). lambda2, the second-smallest
 * eigenvalue of the normalised Laplacian I - D^-1/2 A D^-1/2, by power
 * iteration on (I + D^-1/2 A D^-1/2) / 2 kept orthogonal to D^1/2 1, from a
 * fixed pseudo-random start, until the residual norm is below tol or maxit
 * steps (conv 0); lambda2 is then the Rayleigh quotient of the result. The
 * conductance h = min over node sets S of w(S, not S) / min(vol S, vol not
 * S) satisfies lambda2 / 2 <= h <= sqrt(2 lambda2) (Chung 1997, Spectral
 * Graph Theory, ch. 2). Out, 4: lambda2, lambda2 / 2, sqrt(2 lambda2) and
 * the least conductance of a sweep set (a prefix of the nodes ordered by
 * the eigenvector over D^1/2, which h never exceeds); set: 1 for the nodes
 * of that set. CLAB_EDOM when the graph is not connected. work: 4n doubles;
 * iwork: n ints. O(maxit m). */
int clab_ntchg(int n, const int *off, const int *adj, const double *wt,
               double tol, int maxit, double *out, int *set, int *iter,
               int *conv, double *work, int *iwork);

/* The co-offending graph of an incident table: k rows, row r naming
 * offender off[r] (0 .. no - 1) in incident inc[r] (0 .. ni - 1). An
 * offender named twice in one incident counts once. clab_coofw gives
 * npair, the number of offender pairs met over all incidents (an upper
 * bound on the edges; CLAB_EARG above INT_MAX / 4), and ilen, the ints
 * of iwork clab_coofg needs.
 * clab_coofg writes the one-mode projection: np pairs (pa < pb by code)
 * with pw, the number of incidents they share, in the order each pair is
 * first met (incidents by code, offenders within one by row); cap >=
 * npair. O(k + npair). */
int clab_coofw(const int *inc, const int *off, int k, int ni, int no,
               double *npair, double *ilen, int *iwork);
int clab_coofg(const int *inc, const int *off, int k, int ni, int no,
               int cap, int *pa, int *pb, double *pw, int *np, int *iwork);
/* ------------------------------------------------------------------------
 * Learning: the learners of the cross-fitted estimators (dml.cpp).
 * Matrices are column-major, n rows by p columns, as R stores them.
 * ------------------------------------------------------------------------ */

/* Elastic net path by cyclical coordinate descent (Friedman, Hastie and
 * Tibshirani 2010, doi:10.18637/jss.v033.i01), with the objective, the
 * standardisation and the lambda scale of glmnet: fam 0 (gaussian)
 * minimises (1 / 2n) |y - a0 - x b|^2 + lam pen(b), fam 1 (binomial, y in
 * {0, 1}) the mean negative log-likelihood plus lam pen(b), where pen(b) =
 * sum_j (1 - alpha) / 2 (s_j b_j)^2 + alpha |s_j b_j| with s_j the
 * population standard deviation of column j (a constant column keeps a
 * zero coefficient). For the gaussian family the problem is solved with y
 * scaled to unit variance and lam divided by that sd, as glmnet does, so
 * results match it for every alpha. alpha in [0, 1]; lam: nl lambdas,
 * decreasing, read when lmr <= 0; when 0 < lmr < 1 the path is written to
 * lam instead: nl values from the smallest lambda that zeroes every
 * coefficient down to lmr times it, evenly spaced on the log scale, the
 * path ending early by glmnet's rule (from the fifth lambda, once the share
 * of deviance explained grows by less than 1e-5 of itself or passes
 * 0.999). thresh > 0: convergence when no pass moves a standardised
 * coefficient by more than sqrt(thresh) (weighted, binomial); maxit
 * passes at most in all, else CLAB_EDOM. scr non-zero: binomial paths
 * screen columns by the sequential strong rule (Tibshirani et al. 2012,
 * doi:10.1111/j.1467-9868.2011.01004.x) before each lambda; the
 * optimality conditions are checked over every column either way, so the
 * solution is the same. fold: NULL, or n fold ids in [0,
 * nf), nf >= 2, every fold non-empty: then each fold is held out in turn,
 * the path refitted on the rest at the same lambdas, and cv (2 nl) gets
 * the mean held-out loss (squared error; binomial deviance with the
 * probability kept in [1e-5, 1 - 1e-5]) and its standard error, both as
 * cv.glmnet computes them from fold means weighted by fold size; best (2)
 * the indices of the lambda at the minimum and of the largest lambda
 * within one standard error of it; nfit, the number of lambdas fitted
 * (the first nfit entries of every output). Out: a0 (nl), beta (p nl) on
 * the original scale. work: p^2 + 6p + 5n + n p + nl (p + 1) +
 * nf (nl + 1) doubles; iwork: 2n + p ints. */
int clab_enpath(const double *x, int n, int p, const double *y, int fam,
               double alpha, double *lam, int nl, double lmr, double thresh,
               int maxit, int scr, const int *fold, int nf, double *a0,
               double *beta, double *cv, int *best, int *nfit, double *work,
               int *iwork);

/* Logistic regression by Newton's method on the centred design. y: n
 * outcomes in {0, 1}, both present; x: n by p, p >= 0, full column rank
 * with the intercept. Newton steps until none moves a coefficient by more
 * than tol (1 + |coefficient|), at most maxit. Out: coef and se (p + 1
 * each, intercept first; se from the inverse Fisher information) and dev,
 * the deviance. CLAB_EDOM for a singular information matrix (collinear
 * columns, n <= p + 1) or no convergence (separated data). work:
 * 2 (p + 1)^2 + 4 (p + 1) doubles; iwork: n ints. */
int clab_lgtreg(const double *x, int n, int p, const double *y, double tol,
                int maxit, double *coef, double *se, double *dev,
                double *work, int *iwork);


/* ------------------------------------------------------------------------
 * Spatial kernels added with the R spatial functions (src/geo2.cpp and
 * the families after it).
 * ------------------------------------------------------------------------ */

/* The inverse of clab_eqproj: the longitude and latitude (degrees) of n
 * planar points x, y (east and north, the unit of r) of the Lambert
 * azimuthal equal-area projection about (lon0, lat0) on the sphere of
 * radius r > 0. The unit vector of each point is built from s = rho / (2 r)
 * (cos c = 1 - 2 s^2, sin c / rho = sqrt(1 - s^2) / r), so neither a
 * division by rho nor an inverse sine is needed and the poles and the
 * centre are well conditioned (Snyder 1987, doi:10.3133/pp1395, ch. 24).
 * Longitudes come back in [-180, 180]. A point with a NaN coordinate gives
 * NaN; CLAB_EDOM for another non-finite one or one farther than 2 r from
 * the centre. lon, lat may be x, y. O(n). */
int clab_eqinv(const double *x, const double *y, int n, double lon0,
               double lat0, double r, double *lon, double *lat);

/* Spatial weights (src/swt.cpp). A neighbour structure is a set of
 * compressed sparse rows: the neighbours of unit i are ci[rp[i] ..
 * rp[i + 1]), with rp[0] = 0.
 *
 * clab_swknb: the k nearest neighbours (k in [1, n - 1]) of each of n >= 2
 * planar points, ties at equal distance broken by the lower index; metric
 * 0 Euclidean, 1 Manhattan (|dx| + |dy|). Out: nb (n k ints, row i nearest
 * first) and d (their distances). Grid of cells holding about k + 1 points
 * each, searched by rings: O(n k log k) expected for spread-out points.
 * iwork: 2n + 2 ints.
 *
 * clab_swband: the neighbours within dmax > 0 (closed) of each of n >= 1
 * points, metric as above, the unit itself excluded, each row in
 * increasing order of neighbour. rp: n + 1 ints; ci, d: cap entries; *nnz
 * receives the number of pairs. When it exceeds cap only the first cap are
 * stored: call again with cap >= *nnz (cap = 0 counts). Grid of cells of
 * side about dmax: O(n + pairs within 2 dmax). iwork: 2n + 2 ints.
 *
 * clab_swwt: weights for the rows rp of n units: 1 (pw = 0), or d^-pw for
 * pw > 0 (every distance d positive: CLAB_EDOM otherwise; d is not read
 * when pw = 0); rowstd 1 divides each row by its sum (an empty row stays
 * empty), 0 leaves them. Out: w, rp[n] values. O(nnz). */
int clab_swknb(const double *x, const double *y, int n, int k, int metric,
               int *nb, double *d, int *iwork);
int clab_swband(const double *x, const double *y, int n, double dmax,
                int metric, int cap, int *rp, int *ci, double *d, int *nnz,
                int *iwork);
int clab_swwt(const int *rp, const double *d, int n, double pw, int rowstd,
              double *w);

/* Spatial autocorrelation on weights (src/lsa.cpp). x: n values; rp, ci,
 * w: a weights structure whose rows hold distinct neighbours in increasing
 * order, none the unit itself, weights finite and >= 0 (CLAB_EDOM
 * otherwise). S0 = sum w_ij, S1 = sum_ij (w_ij + w_ji)^2 / 2, S2 = sum_i
 * (w_i. + w_.i)^2 (Cliff and Ord 1981). Permutation tests shuffle x over
 * the units nsim >= 0 times from seed.
 *
 * clab_gearc: Geary's C (Geary 1954, doi:10.2307/2986645), n >= 4,
 * C = (n - 1) sum_ij w_ij (x_i - x_j)^2 / (2 S0 sum (x_i - xbar)^2), with
 * its variance under normality and under randomisation (Cliff and Ord
 * 1981, with b2 = n sum z^4 / (sum z^2)^2), the z statistics (1 - C) / sd
 * and their upper normal tails (positive autocorrelation makes C small),
 * and the permutation p-value (1 + #{C* <= C}) / (nsim + 1) (NaN without
 * permutations). Out, 9: C, Vn, Vr, zn, zr, pn, pr, pperm, S0. CLAB_EDOM
 * for constant x or S0 = 0. work: n doubles. O(nnz log(deg) + nsim nnz).
 *
 * clab_gorg: the general G of Getis and Ord (1992,
 * doi:10.1111/j.1538-4632.1992.tb00261.x), x >= 0, n >= 4: G = sum_ij
 * w_ij x_i x_j / sum_{i != j} x_i x_j, E[G] = S0 / (n (n - 1)), and its
 * variance under randomisation from the moments of x (their B0 to B4).
 * Out, 6: G, E, V, z, the upper normal tail of z, and the permutation
 * p-value (1 + #{G* >= G}) / (nsim + 1). CLAB_EDOM when fewer than two
 * values are positive or S0 = 0. work: n doubles.
 *
 * clab_lmorn: local Moran (Anselin 1995,
 * doi:10.1111/j.1538-4632.1995.tb00338.x), n >= 3: with z = x - xbar and
 * m2 = sum z^2 / n, I_i = z_i / m2 sum_j w_ij z_j, E = -w_i. / (n - 1)
 * and the variance under total randomisation (b2 = (sum z^4 / n) / m2^2);
 * the two-sided normal p-value; the conditional permutation p-value
 * (unit i keeps its value, its neighbours' values are drawn without
 * replacement from the other n - 1), (1 + min(#{I* >= I}, #{I* <= I})) /
 * (nsim + 1); and the quadrant: 1 high-high, 2 low-low, 3 low-high, 4
 * high-low, 0 when z_i or the lag is 0. Out, 7n, column by column. iwork:
 * n ints; work: n doubles. O(nnz + n (n + nsim deg)).
 *
 * clab_locgi: the local G (star = 0, the unit excluded) or G* (star = 1,
 * the unit in its own neighbourhood with weight 1) of Getis and Ord (1992)
 * as z values (Ord and Getis 1995, doi:10.1111/j.1538-4632.1995.tb00912.x),
 * the means and variances taken over the other n - 1 values for G and over
 * all n for G*; NaN for a unit without weight. Out: n values. O(n + nnz). */
int clab_gearc(const double *x, int n, const int *rp, const int *ci,
               const double *w, int nsim, unsigned int seed, double *work,
               double *out);
int clab_gorg(const double *x, int n, const int *rp, const int *ci,
              const double *w, int nsim, unsigned int seed, double *work,
              double *out);
int clab_lmorn(const double *x, int n, const int *rp, const int *ci,
               const double *w, int nsim, unsigned int seed, int *iwork,
               double *work, double *out);
int clab_locgi(const double *x, int n, const int *rp, const int *ci,
               const double *w, int star, double *out);

/* Directional mean of n >= 1 angles a (radians) with weights w (NULL or n
 * finite values >= 0, positive sum): the mean direction atan2(S, C) in
 * [0, 2 pi) (NaN when the resultant is 0), the mean resultant length R in
 * [0, 1], the circular variance 1 - R and the circular standard deviation
 * sqrt(-2 log R) (Mardia and Jupp 2000, doi:10.1002/9780470316979, ch. 2),
 * and the total weight. Out, 5. O(n). */
int clab_dirmn(const double *a, const double *w, int n, double *out);
/* A line chart as a standalone SVG document (UTF-8), the picture of the
 * plots of every binding. k (1 to 3) curves: curve j has np[j] >= 1 points
 * (x[i], y[i]) from i = off[j], drawn as steps (step[j] 1: across to the
 * next x, then up or down) or as straight lines (0), dash[j] 0 solid, 1
 * dashed, 2 dotted, with lab[j] in the legend. title, xlab, ylab: text,
 * escaped here. vline 1 adds a dotted vertical line at x = vx. corner puts
 * the legend top left (0), top right (1), bottom left (2) or bottom right
 * (3) inside the plot box. The canvas
 * is 640 x 400 px; the plot box spans x 64 to 624 and y 40 to 348. The x
 * axis spans the data (and vx), the y axis the data and 0, each padded 4%
 * either way (a single value is first widened by 1 either way); ticks fall
 * at 1, 2 or 5 times a power of 10. Text is real text (a system sans-serif
 * stack); colours suit light and dark backgrounds (prefers-color-scheme).
 * Two passes: with cap 0, *len receives the length of the document
 * (without its NUL); with cap > *len the document is written into out,
 * NUL-terminated. Numbers are written without printf: the bytes do not
 * depend on the locale or the platform. CLAB_EDOM for a non-finite point
 * or vx; CLAB_EARG for k or corner out of range, an empty curve or a cap
 * too small. */
int clab_svgln(int k, const int *off, const int *np, const double *x,
               const double *y, const int *step, const int *dash,
               const char *const *lab, const char *title, const char *xlab,
               const char *ylab, int vline, double vx, int corner, char *out,
               int cap, int *len);
/* Every reading of a contrast reported from an r x c table of counts
 * (Lean: Crimlab.RelChk.or_swap_ref and the theorems beside it).
 * x: column-major (x[i + r j]: group i, outcome level j), r, c >= 2,
 * finite and >= 0, with every row and column total positive. A reading
 * compares group g with reference group h (g != h) for outcome level e as
 * the event, the other levels pooled: kind 0, the odds ratio a d / (b c);
 * 1, the risk ratio (a / (a + b)) / (c / (c + d)); 2, the risk difference
 * a / (a + b) - c / (c + d), where a and b are the events and the others in
 * group g, c and d in group h. A zero denominator gives +Inf over a
 * positive numerator and NaN over 0. The readings are taken with the
 * groups in rows (grp 0) and, the table transposed, in columns (grp 1), in
 * the order orientation, kind, h, g, e: m = 3 r c (r + c - 2) of them,
 * CLAB_EARG above INT_MAX. A reading matches value when it is finite and
 * within tol >= 0 of it, up to 4 units in the last place of the larger of
 * the two. Out, m each: grp, g, h, e, kind, v (the value), match (1 or 0);
 * nmatch, the number that match. CLAB_EDOM for a count, value or tol out
 * of range or an empty row or column. O(m (r + c)). */
int clab_relchk(const double *x, int r, int c, double value, double tol,
                int *grp, int *g, int *h, int *e, int *kind, double *v,
                int *match, int *nmatch);

/* The audit of one categorical variable. code: n level codes in [0, k), or
 * -1 for a missing value; fold: k keys in [0, k), equal for labels that
 * differ only by case or white space; nmin >= 1. Out: cnt, k counts; twin,
 * k levels: the first other level with the same key (the first level of a
 * key gets the second), or -1; s, 4: the missing values, the levels with
 * no count, the levels counted fewer than nmin times (but at least once)
 * and the levels with a twin. work: k ints. O(n + k). */
int clab_catchk(const int *code, int n, int k, const int *fold, int nmin,
                int *cnt, int *twin, int *s, int *work);
/* Random forest (Breiman 2001, doi:10.1023/A:1010933404324) for
 * regression, and for the probability of a 0/1 outcome (the leaf mean; the
 * variance split of a 0/1 outcome is the Gini split). x: n >= 2 rows, p
 * columns to grow on, y their outcomes; xt: nt >= 0 rows to predict, the
 * same columns. Each column is cut at its order statistics over x into at
 * most 256 bins (midpoints between neighbouring values when it has at
 * most 256 distinct ones, else 255 quantiles), so the forest is unchanged
 * by an increasing transform of a column. ntree >= 1 trees, each on a bag
 * of round(frac n) rows (0 < frac <= 1) drawn with replacement when repl,
 * without otherwise, from a stream seeded by seed and the tree's number;
 * at each node the best split by the reduction in the sum of squares over
 * mtry (1..p) columns drawn without replacement; nodes of at most minnode
 * rows, or with one outcome value, are leaves; a split must leave at least
 * imb (0 <= imb < 0.5) of its node's rows, and one row, on each side. hon
 * non-zero grows honest trees (Wager and Athey 2018,
 * doi:10.1080/01621459.2017.1319839): a shuffled bag's first half (round
 * frac n / 2 rows, at least 1) places the splits and its second half sets
 * the leaf means, a leaf with no row of it taking its deepest such
 * ancestor's mean. Out: pred (nt), the mean
 * over trees; oob (NULL or n), the mean over the trees whose bag left the
 * row out (NaN if none). work: p 255 + 3n + 512 + nt + 5 (2n + 1)
 * doubles; iwork: ceil((n + nt) p / 4) + 4n + nt + p + 1 ints. */
int clab_rfor(const double *x, int n, int p, const double *y,
              const double *xt, int nt, int ntree, int mtry, int minnode,
              double frac, int repl, int hon, double imb, unsigned int seed,
              double *pred, double *oob, double *work, int *iwork);

/* ------------------------------------------------------------------------
 * Causal: double machine learning by cross-fitting (Chernozhukov et al.
 * 2018, doi:10.1111/ectj.12097), with the learners above.
 * ------------------------------------------------------------------------ */

/* A nuisance learner. kind 0: least squares, or logistic regression for a
 * probability; kind 1: the elastic net path of clab_enpath() with mixing
 * alpha, nlam lambdas from the largest down to 1e-4 of it (1e-2 when the
 * training rows do not outnumber the columns), the lambda of least
 * cross-validated loss over nfold >= 2 inner folds (a shuffle of the
 * training rows by seed, dealt round robin); kind 2: the random forest of
 * clab_rfor() with ntree >= 1 trees, bags of the training rows' size drawn
 * with replacement, mtry columns per split (0: the floor of sqrt(p); more
 * than p: p) and nodes of at most minnode rows left unsplit (0: 5 for a
 * regression, 10 for a probability, as ranger); kind 3: the stack of
 * clab_stack() over the nsub learners sub (kinds 0 to 2), on nstk >= 2
 * inner folds (a shuffle of the training rows by seed, dealt round robin,
 * as clab_kfold()), with disc as there. */
typedef struct clab_lrn {
  int kind;
  double alpha;
  int nlam, nfold;
  unsigned int seed;
  int ntree, mtry, minnode;
  double tol;  /* elastic net: the convergence threshold of clab_enpath() */
  int scr;     /* elastic net: screening of binomial paths */
  int hon;     /* forest: honest trees on half-samples without replacement */
  const struct clab_lrn *sub; /* stack: the learners stacked */
  int nsub, nstk, disc;       /* stack: their number, the folds, disc */
} clab_lrn;

/* The number of nuisance functions of a model (0 for an unknown one). */
int clab_dmlxk(int model, int sc);

/* Scratch for clab_dmlx(): returns the doubles, sets *niw to the ints. */
long clab_dmlxsz(int model, int n, int p, int nf, int nrep,
                const clab_lrn *ly, const clab_lrn *ld, long *niw);

/* Cross-fitted double machine learning. model: 0 partially linear
 * regression (PLR), 1 partially linear IV (PLIV, one instrument z), 2
 * interactive regression (IRM; d in {0, 1}), 3 interactive IV (IIVM; d, z
 * in {0, 1}), 4 augmented inverse probability weighting (the IRM average
 * effect with the AIPW package's variance), 5 targeted maximum likelihood
 * (the tmle package's fluctuation of the IRM nuisances), 6 difference in
 * differences on a panel (y the change in the outcome), 7 difference in
 * differences on repeated cross-sections (s the period, 0 or 1), 8 sample
 * selection with outcomes missing at random (s, 0 or 1, marks an observed
 * outcome; y is not read where s is 0). sc: 0 for the partialling-out
 * score (PLR, PLIV) or the average effect (IRM), 1 for the IV-type score
 * or the effect on the treated. x: n by p covariates; z, s: NULL unless
 * the model reads them. fold: nrep columns of n fold ids in [0, nf); each
 * column a repetition. clu: NULL or n cluster ids in [0, nc), one-way
 * clustering (PLR and PLIV only; folds constant within a cluster); with
 * nc2 > 0, two-way clustering (nrep 1): clu holds 2n ids, the first n in
 * [0, nc), the next n in [0, nc2), and fold the folds of each kind of
 * cluster (2n); the test set of fold pair (k, l) is the rows in both, its
 * training set the rows in neither, and the variance is the two-way form
 * of Chiang et al. as DoubleML computes it (nc2 = 0: one-way). ly
 * learns the outcome nuisances, ld those of the treatment, instrument and
 * selection. trim: propensities are kept in [trim, 1 - trim] (TMLE: only
 * bounded below, as tmle). fit 1: the nuisances are learned and written to
 * pred; fit 0: pred holds them. pred: nrep blocks of clab_dmlxk() columns
 * of n, in the order of the model (PLR: l, m[, g]; PLIV: l, m, r[, g];
 * IRM, AIPW, TMLE: g0, g1, m; IIVM: g0, g1, m, r0, r1; DiD: g0, m; DiD on
 * cross-sections: g00, g01, g10, g11, m; selection: pi, m, g0, g1). Out:
 * est (2 + 2 nrep): the estimate (the median over repetitions) and its
 * standard error (DoubleML's aggregation), then each repetition's
 * estimates and standard errors; psi (2 n nrep): the score elements psi_a
 * then psi_b of each repetition. work, iwork: clab_dmlxsz(). CLAB_EDOM
 * when a learner fails, a training set is empty or a score has no root. */
int clab_dmlx(int model, int sc, const double *x, int n, int p,
             const double *y, const double *d, const double *z,
             const double *s, const int *fold, int nf, int nrep,
             const int *clu, int nc, int nc2, const clab_lrn *ly,
             const clab_lrn *ld, double trim, int fit, double *pred,
             double *est, double *psi, double *work, int *iwork);

/* K-fold assignment: the n rows shuffled by seed (Fisher-Yates, the
 * xorshift64* stream of the other kernels) and dealt round robin, so that
 * the k folds differ in size by at most one. Out: fold, n ids in [0, k).
 * The inner folds of the learners of clab_dmlx() are these. iwork: n
 * ints. */
int clab_kfold(int n, int k, unsigned int seed, int *fold, int *iwork);

/* Stacking: the super learner (van der Laan, Polley and Hubbard 2007,
 * doi:10.2202/1544-6115.1309). l: a learner of kind 3 (see clab_lrn), its
 * nsub learners fitted on x (n by p, column-major) to y (fam 0: numbers;
 * fam 1: 0 or 1, the learners then predicting probabilities). fold: n ids
 * in [0, l->nstk), each used, with at least one row outside each fold.
 * Each learner is fitted on the rows outside each fold and predicts that
 * fold; oof (n by nsub) holds these cross-fitted predictions and risk
 * (nsub) their mean squared error (the Brier score for fam 1). wt (nsub)
 * are the meta-weights: the non-negative least squares fit of y on oof
 * scaled to sum to 1 (SuperLearner's method.NNLS), or with disc non-zero,
 * or when that fit is all zero, 1 for the learner of least risk (the first
 * on ties) and 0 for the rest. pred (nt): the stack's prediction at the nt
 * rows of xt, sum_j wt[j] f_j(xt), each f_j fitted on all n rows with
 * l->seed (a learner of weight 0 is not fitted). The fold fits use seed
 * l->seed + 7919 (k + 1). With one learner the stack is that learner,
 * bit for bit. work: clab_stksz() doubles; iwork: the ints it sets. */
long clab_stksz(int n, int p, int nt, const clab_lrn *l, long *niw);
int clab_stack(const double *x, int n, int p, const double *y, int fam,
               const clab_lrn *l, const int *fold, const double *xt, int nt,
               double *oof, double *risk, double *wt, double *pred,
               double *work, int *iwork);


/* Mantel's test (clab_mantst) with the spatial dissimilarity built from
 * Euclidean (metric 0) or Manhattan (metric 1, |dx| + |dy|) distances;
 * clab_mantst is metric 0. */
int clab_mantm(const double *x, const double *y, const double *t, int n,
               int metric, double cs, double ct, int nsim, unsigned int seed,
               int *iwork, double *out);

/* Ripley's K (clab_ripk) with pairs counted by Euclidean (metric 0) or
 * Manhattan (metric 1) distance. Under complete spatial randomness the
 * Manhattan K is 2 r^2, the area of the Manhattan ball. The border and
 * modified border estimators carry over unchanged: a Manhattan ball of
 * radius r lies in the rectangle exactly when its centre is at least r
 * from every edge, the Euclidean condition. The isotropic weight (the
 * share of a circle) has no Manhattan counterpart here: corr 3 with
 * metric 1 is CLAB_EDOM. */
int clab_ripkm(const double *x, const double *y, int n, const double *win,
               const double *r, int m, int corr, int metric, int *iwork,
               double *work, double *k);

/* The nearest-neighbour index under Manhattan distance (src/sp2.cpp),
 * laid out as clab_nnidx. A Manhattan ball of radius r has area 2 r^2, so
 * under complete spatial randomness without edges the k-th neighbour
 * distance has E R_k = Gamma(k + 1/2) / (Gamma(k) sqrt(2 lam)) (sqrt(pi / 2)
 * times the Euclidean value) and E R_k^2 = k / (2 lam); order 1's z uses
 * Var R_1 = (4 - pi) / (8 lam). No edge-corrected expectation is published
 * for Manhattan distance: Donnelly's expectation, the ratio against it and
 * the corrected z are NaN. The Monte Carlo p-values (simulated uniform
 * points in the rectangle) need no expectation and keep exact size.
 * Neighbours by clab_swknb. iwork: 2n + 2 + n kmax + 2 kmax ints; work: n
 * kmax + 2n + 2 kmax doubles. */
int clab_nnidm(const double *x, const double *y, int n, const double *win,
               int kmax, int nsim, unsigned int seed, double *dist,
               int *iwork, double *work, double *out);

/* Head-bang smoothing (Hansen 1991, IEEE Transactions on Geoscience and
 * Remote Sensing 29: 369-378; the weighted screens after Mungiole, Pickle
 * and Simonson 1999, Statistics in Medicine 18: 3201-3209) of the values v at n >= 3 points. For each point i: its nn nearest
 * neighbours (2 <= nn <= n - 1, Euclidean, ties by index); the triples
 * (j, i, k) of two of them whose angle at i is at least ang degrees (90
 * to 180; coincident neighbours are skipped), so that i lies roughly
 * between j and k; of these the ntrip with the smallest |p_j - p_i| +
 * |p_k - p_i| (all when ntrip is 0). Each triple gives a low value
 * min(v_j, v_k) and a high one; the low screen L is the median of the
 * lows and the high screen H of the highs (with weights w, NULL or n
 * positive values, the weighted medians, each value weighted by its
 * unit's w: the smallest value whose cumulative weight reaches half the
 * total); the new value is the median of L, v_i and H. A point with no
 * triple keeps its value (Hansen's extrapolated end triples are not
 * formed). iter >= 0 passes, each from the previous values. iwork: 2n + 2
 * + n nn + 3T ints, work: n nn + n + 5T doubles, T = nn (nn - 1) / 2.
 * Out: n values. */
int clab_hdbang(const double *x, const double *y, const double *v,
                const double *w, int n, int nn, double ang, int ntrip,
                int iter, int *iwork, double *work, double *out);

/* Correlated random walk of n >= 3 successive event locations (Kareiva and
 * Shigesada 1983, doi:10.1007/BF00379695). clab_crwfit: step lengths
 * step (n - 1), turning angles turn (n - 2, in (-pi, pi], the change of
 * heading between successive steps) and out, 5: the mean step, the mean
 * squared step, the mean cosine c and mean sine s of the turns, and the
 * last heading; CLAB_EDOM for a repeated location (no heading).
 * clab_crwmsd: for k = 1 .. h, obs, the observed mean squared displacement
 * over k moves (over every start; NaN when k >= n), and msd, Kareiva and
 * Shigesada's expected squared displacement for independent steps and
 * turns with symmetric turning (s = 0): k E l^2 + 2 (E l)^2 c / (1 - c)
 * (k - (1 - c^k) / (1 - c)), and k E l^2 + (E l)^2 k (k - 1) at c = 1.
 * clab_crwsim: nsim walks of h steps from (x0, y0) with heading head, each
 * turn and each step drawn independently and uniformly from the observed
 * ones (the walk the expectation describes), from seed. Out: px, py, h
 * nsim each, step k of walk s at k nsim + s. */
int clab_crwfit(const double *x, const double *y, int n, double *step,
                double *turn, double *out);
int clab_crwmsd(const double *x, const double *y, int n, double ml,
                double ml2, double c, int h, double *obs, double *msd);
int clab_crwsim(const double *step, int ns, const double *turn, int nt,
                double x0, double y0, double head, int h, int nsim,
                unsigned int seed, double *px, double *py);

/* Risk-adjusted nearest-neighbour clustering (first order). The
 * random-expectation threshold of clab_nnhclu, 0.5 sqrt(A / n) + z
 * sqrt((4 - pi) A / (4 pi)) / n, is f / sqrt(lam) with lam = n / A and f =
 * 0.5 + z sqrt((4 - pi) / (4 pi n)); here lam is not uniform but the
 * intensity the n incidents would have if they followed the baseline:
 * lam_i = n b(p_i) / B, b the normal-kernel density (bandwidth bwid) of the
 * nb baseline units (bx, by) with sizes bw (B their sum). Points i and j
 * are linked when |p_i - p_j| < f / sqrt((lam_i + lam_j) / 2) (always when
 * that mean is 0). Grouping as in clab_nnhclu's first order; groups of at
 * least minpts are kept, numbered by size. Out: thr (n, f / sqrt(lam_i)),
 * cl (n, group or -1), *ncl, cs (8 per group: points, mean x, y, Weber x,
 * y, ellipse a, b, angle as clab_sdell with k, c). iwork: 4n ints; work: 3n
 * doubles. O(n^2 + n nb). */
int clab_rnnhcl(const double *x, const double *y, int n, const double *bx,
                const double *by, const double *bw, int nb, double bwid,
                double z, int minpts, int k, double c, int *iwork,
                double *work, double *thr, int *cl, int *ncl, double *cs);
/* ------------------------------------------------------------------------
 * Counts: log-linear regression of a count on a group, with exposure and
 * clusters.
 * ------------------------------------------------------------------------ */

/* Poisson (nb 0) or negative binomial (nb 1: variance mu + mu^2 / theta,
 * theta by maximum likelihood) regression of y on x with a log link, by
 * the steps of R's glm.fit (iteratively reweighted least squares from mu =
 * y + 0.1, until the relative change of the deviance is below eps, at most
 * maxit steps) and, for nb, of MASS::glm.nb (Venables and Ripley 2002;
 * theta and the coefficients in turn). y: n counts, at least 0, not all 0
 * (whole numbers or not: the likelihood uses lgamma(y + 1)); x: n by p,
 * column-major, the first column usually all 1; off: NULL or n offsets
 * (log exposure). A column whose part not explained by the columns before
 * it is below 1e-11 of its norm is aliased (alias[j] 1, beta[j] NaN, its
 * rows and columns of vb and vc NaN), as in R. cl: NULL or n cluster ids
 * in [0, ncl). Out: beta (p); vb (p x p), the model-based covariance
 * (X'WX)^-1; vc (p x p, not read when cl is NULL), the cluster sandwich
 * (X'WX)^-1 (sum_g U_g U_g') (X'WX)^-1 with no small-sample factor; st
 * (8): log-likelihood, AIC, theta and its standard error (NaN for
 * Poisson), steps taken, converged (1 or 0), rank, deviance. work: n p +
 * 6 n + p + p^2 + 2 ncl p doubles; iwork: p ints. */
int clab_cntreg(const double *y, const double *x, int n, int p, const double *off,
                const int *cl, int ncl, int nb, int maxit, double eps, double *beta,
                double *vb, double *vc, double *st, int *alias, double *work,
                int *iwork);

/* The Poisson log-linear model by generalised estimating equations
 * (Liang and Zeger 1986, doi:10.1093/biomet/73.1.13) with an exchangeable
 * working correlation, estimated as geepack does: phi = sum r^2 / n and
 * alpha = sum over pairs within a cluster of r_i r_j / (phi times the
 * number of pairs), r the Pearson residuals. y, x, off, cl as in
 * clab_cntreg (cl required); alias: p flags, the columns left out;
 * beta: p, in the start (the Poisson fit) and out the estimate, NaN where
 * aliased. Fisher scoring until no coefficient moves more than eps, at
 * most maxit steps. Out: vc (p x p), the robust (sandwich) covariance; st
 * (4): alpha, phi, steps, converged. CLAB_EDOM when the working
 * correlation is not positive definite. work: 2 n + 2 p^2 + 2 p + 2 ncl p
 * doubles; iwork: p + ncl + 1 + n ints. */
int clab_geeex(const double *y, const double *x, int n, int p, const double *off,
               const int *cl, int ncl, const int *alias, int maxit, double eps,
               double *beta, double *vc, double *st, double *work, int *iwork);

/* ------------------------------------------------------------------------
 * Rates by group: ratios and differences, direct and indirect
 * standardisation.
 * ------------------------------------------------------------------------ */

/* Rates (prop 0: x events over exposure n) or proportions (prop 1: x of
 * n, x <= n) of k groups against group ref, with z the normal quantile of
 * the intervals. x at least 0, n positive. Out, k each: rate (x / n); lo,
 * hi (Wilson's score interval for a proportion, rate exp(-+z / sqrt(x))
 * for a rate, NaN when x is 0); rep (the share of the events over the
 * share of n); ratio (rate over the reference's; Inf or NaN over a zero
 * rate) with rlo, rhi (log-normal: 1/x - 1/n + 1/x0 - 1/n0 for risks,
 * Katz et al. 1978, doi:10.2307/2530610; 1/x + 1/x0 for rates; NaN with a
 * zero count and for ref); diff (rate minus the reference's) with dlo,
 * dhi (Newcombe's 1998 hybrid score interval for risks,
 * doi:10.1002/(SICI)1097-0258(19980430)17:8<873::AID-SIM779>3.0.CO;2-I;
 * Wald for rates; NaN for ref). */
int clab_grprr(const double *x, const double *n, int k, int ref, int prop, double z,
               double *rate, double *lo, double *hi, double *rep, double *ratio,
               double *rlo, double *rhi, double *diff, double *dlo, double *dhi);

/* Directly standardised rates. m rows, one per group and stratum: x
 * events, n positive exposure, grp in [0, ng), str in [0, ns); w: ns
 * standard weights, at least 0 with a positive sum; every group has every
 * stratum of positive weight. Out, ng each: rate (sum of w x / n over sum
 * w), var (its Poisson variance), wmax (the largest w / n, for the
 * interval of Fay and Feuer 1997, doi:10.1002/(SICI)1097-0258(19970415)
 * 16:7<791::AID-SIM500>3.0.CO;2-#), ev (events); against group ref,
 * ratio, rlo, rhi (log-normal) and diff, dlo, dhi (normal), with z the
 * normal quantile (NaN for ref). seen: ng ns ints. */
int clab_dastd(const double *x, const double *n, const int *grp, const int *str, int m,
               int ng, int ns, const double *w, int ref, double z, double *rate,
               double *var, double *wmax, double *ev, double *ratio, double *rlo,
               double *rhi, double *diff, double *dlo, double *dhi, int *seen);

/* Indirect standardisation. Rows as for clab_dastd (any number per group
 * and stratum); sr: NULL for the stratum rates pooled over all groups, or
 * ns rates at least 0. Out, ng each: obs, expd (the sum of n times the
 * stratum rate, positive), sir (obs / expd), eb (shrunk towards the
 * overall ratio by Marshall's 1991 moment estimates of a gamma prior,
 * doi:10.2307/2347593) and shr (the share of the prior in it); prior (2):
 * the prior mean and variance (0 when the ratios vary no more than Poisson
 * noise, and then eb is the mean). work: 2 ns doubles. */
int clab_sirex(const double *x, const double *n, const int *grp, const int *str, int m,
               int ng, int ns, const double *sr, double *obs, double *expd, double *sir,
               double *eb, double *shr, double *prior, double *work);

/* ------------------------------------------------------------------------
 * Contingency tables.
 * ------------------------------------------------------------------------ */

/* The r x c table x (column-major counts from 0, every row and column
 * total positive). yates: 1 for the continuity correction of a 2 x 2
 * table's chi-square, as chisq.test(correct = TRUE). dg: NULL, or for each
 * row the column holding the same label (-1 for none), for the diagonal
 * share. Out: e and res (r x c), the expected counts and the Pearson
 * residuals (o - e) / sqrt(e); st (9): n, the chi-square statistic, its
 * degrees of freedom, the likelihood-ratio statistic G^2, Cramer's V
 * (from the uncorrected statistic), V with Bergsma's bias correction
 * (doi:10.1016/j.jkss.2012.10.002), the contingency coefficient, phi
 * (signed; 2 x 2 only, NaN otherwise) and the diagonal share (NaN when dg
 * is NULL or matches nothing). */
int clab_ctabx(const double *x, int r, int c, int yates, const int *dg, double *e,
               double *res, double *st);

/* Fisher's exact test of the r x c table x (whole counts, r, c >= 2):
 * every table with its margins is visited, at most cap of them
 * (CLAB_EDOM beyond), and the p-value sums the probabilities at most the
 * observed one's (to a relative 1e-7, as fisher.test()). Out: st (5): the
 * p-value; for a 2 x 2 table the conditional maximum likelihood odds ratio
 * and its exact interval at level conf (NaN otherwise); the tables
 * visited. rem: r + c doubles. */
int clab_fshx(const double *x, int r, int c, double cap, double conf, double *st,
              double *rem);

/* ------------------------------------------------------------------------
 * Policy evaluation on unit-by-period panels.
 * ------------------------------------------------------------------------ */

/* Two-way fixed-effects regression (didfe(), evtst()): y (n) on the p
 * columns of x (n x p, column-major) with unit effects (fe 1; u, nu codes)
 * and period effects (t, nt codes), by least squares after the within
 * transformation (alternating projections until the means left are below
 * tol of each column's largest value, at most maxit sweeps); w: NULL or n
 * positive weights. The variance is CR1 by the clusters cl (ncl), or
 * two-way with cl2 (ncl2; NULL for one-way): V1 + V2 - V12, each with the
 * factor G / (G - 1) (n - 1) / (n - K) of G the smaller number of
 * clusters, K = p + 1 + the levels less one of each fixed effect not
 * nested in a cluster ('fixest's "nested" count). Out: beta (p), vc (p x
 * p), st (5): the sweeps, G, K, the factor, the weighted residual sum of
 * squares. work: n (p + 2) + p p + 2 (nu + nt) + (n + 2 p) p + p doubles;
 * iwork: n + max(nu, nt, ncl ncl2) ints. CLAB_EDOM when a column of x is
 * absorbed by the effects (or the columns are collinear), or the
 * transformation does not converge. */
int clab_didfe(const double *y, const double *x, int n, int p, const int *u, int nu,
               const int *t, int nt, int fe, const double *w, const int *cl, int ncl,
               const int *cl2, int ncl2, double tol, int maxit, double *beta, double *vc,
               double *st, double *work, int *iwork);
/* Near repeats (src/knx.cpp). The Knox table of n >= 3 incidents at
 * planar points (x, y) and times t (all finite): every pair within ds[ks -
 * 1] falls in the distance band (ds[k - 1], ds[k]] and the time band
 * (dt[l - 1], dt[l]] (the first bands from 0, 0 included; edges positive
 * and increasing, CLAB_EDOM otherwise). Out, 3 ks kt + 1: the observed
 * pairs per cell (distance band fastest), the expectation under
 * exchangeable times S_k T_l / P (Knox 1964, doi:10.2307/2985220; S_k and
 * T_l the pairs in each band, P = n (n - 1) / 2), the permutation p-value
 * (1 + #{obs* >= obs}) / (nsim + 1) with the times shuffled nsim >= 0 times
 * from seed (NaN for nsim = 0), and the number of pairs within ds[ks - 1].
 * O(n log n + pairs (kt + nsim log kt)). */
int clab_knxnr(const double *x, const double *y, const double *t, int n,
               const double *ds, int ks, const double *dt, int kt, int nsim,
               unsigned int seed, double *out);

/* Hot spots (src/hot.cpp).
 *
 * clab_kdsurf: the kernel intensity sum_j m_j k(u - p_j) of n >= 1 planar
 * points with weights w (NULL for 1; finite, >= 0) and the product normal
 * kernel of standard deviations hx, hy > 0, at u on the grid ax x ay (grid
 * = 1: na x nb values, ax fastest) or at the na points (ax[a], ay[a])
 * (grid = 0, nb = na); loo = 1 (points only, na = n, the evaluation points
 * the data) leaves each point out of its own sum. edge 0 none; with the
 * rectangle win = {xmin, xmax, ymin, ymax}, 1 divides the value at u by
 * e(u), the kernel's mass in win (Diggle 1985, doi:10.2307/2347366; NaN
 * where it is 0), 2 divides each point's kernel by e(p_j) (Jones 1993,
 * doi:10.1007/BF00147776; CLAB_EDOM where it is 0). The grid is separable:
 * O((na + nb) n + na nb n) without exponentials in the double loop.
 *
 * clab_dbscn: DBSCAN (Ester et al. 1996) of n points: core points have at
 * least minpts points within eps > 0 (closed, the point included); clusters
 * grow from core points in index order; border = 1 puts a border point in
 * the first cluster reaching it, 0 leaves it noise. Out: cl (cluster from
 * 0, -1 noise), core (0/1), *ncl. O(n + pairs within eps).
 *
 * clab_hsper: hot-spot persistence. n points with periods prd in [0, np)
 * in square cells of side cell > 0 from the corner (x0, y0); the m active
 * cells in order of column then row (out: cid per point, ix, iy their
 * column and row, *nc = m; CLAB_EDOM for an index beyond 1e15). In each
 * period the k = max(1, ceil(top m)) busiest cells (0 < top <= 1) and those
 * tied with the k-th are hot (hot: m x np, cells fastest), unless the tie
 * makes more than 2k, when the tied cells are not. per, 9 per period: n,
 * active cells, the fewest cells holding half (NaN for none), that over m,
 * hot cells, tied in, tied out, the share of incidents in hot cells, and
 * the Jaccard overlap with the period before (NaN when either set is
 * empty). cs, 3 per cell: incidents, periods hot, the longest run of hot
 * periods. O(n log n + m np log m).
 *
 * clab_ebshr: linear empirical-Bayes shrinkage of n >= 2 values y with
 * weights w (NULL for 1; positive sum) and noise variance v >= 0: with ybar
 * and S the weighted mean and scatter, the signal max(S - v, 0) and B* = v /
 * (v + signal) (0 when both are 0). b in [0, 1], or NaN for B*. Out: fit,
 * 6 (ybar, S, v, signal, B*, b), shr (1 - b) y + b ybar and fall b (y -
 * ybar). Lean: Crimlab.Hotspot.bstar_min. O(n). */
int clab_kdsurf(const double *x, const double *y, const double *w, int n,
                double hx, double hy, const double *ax, int na,
                const double *ay, int nb, int grid, const double *win, int edge,
                int loo, double *out);
int clab_dbscn(const double *x, const double *y, int n, double eps, int minpts,
               int border, int *cl, int *core, int *ncl);
int clab_hsper(const double *x, const double *y, const int *prd, int n, int np,
               double cell, double x0, double y0, double top, int *cid,
               double *ix, double *iy, int *nc, int *hot, double *per,
               double *cs);
int clab_ebshr(const double *y, const double *w, int n, double v, double b,
               double *fit, double *shr, double *fall);

/* The space-time permutation scan (src/kscn.cpp; Kulldorff et al. 2005,
 * doi:10.1371/journal.pmed.0020059). n >= 3 incidents at planar points (x,
 * y) and times t, all finite. Cylinders: a circle of radius r[k] (positive,
 * increasing) about each distinct location (closed) times a window [u, u +
 * tw[w]) (tw positive) starting at each distinct time u. A cylinder with c
 * >= minc incidents, n_C in its circle and n_W in its window expects e =
 * n_C n_W / n; its score is the Poisson log-likelihood ratio c log(c / e)
 * + (n - c) log((n - c) / (n - e)) for c > e, 0 otherwise (Lean:
 * Crimlab.Hotspot.llr_nonneg). Clusters: the best cylinder of each circle,
 * by score, keeping those whose circles overlap none kept before, at most
 * top, positive scores only. Each p-value: (1 + #{max* >= score}) / (nsim +
 * 1), the largest score with the times shuffled nsim >= 0 times from seed
 * (NaN for nsim = 0). Out: 10 per cluster (centre index, radius, start,
 * end, c, e, n_C, n_W, score, p), *ncl clusters. O(L R (U + W U) (1 +
 * nsim)) for L locations, R radii, U times, W lengths. */
int clab_kscan(const double *x, const double *y, const double *t, int n,
               const double *r, int nr, const double *tw, int nw, int minc,
               int nsim, unsigned int seed, int top, double *out, int *ncl);

/* Spatial econometrics (src/sreg.cpp). Weights rp, ci, w: rows of
 * distinct neighbours, none the unit itself, weights finite and >= 0
 * (CLAB_EDOM otherwise); S0, S1, S2 of Cliff and Ord (1981).
 *
 * clab_gmoran: Moran's I (Moran 1950, doi:10.2307/2332142) of each of k
 * columns of x (n >= 4 rows, column-major): I = n / S0 sum w_ij z_i z_j /
 * sum z^2, E = -1 / (n - 1), the variances under normality and
 * randomisation (Cliff and Ord 1981), z statistics, upper normal tails, and
 * the permutation p-value (1 + #{I* >= I}) / (nsim + 1) (NaN for nsim = 0;
 * column c draws from seed + c). Out, 10 per column: I, E, Vn, Vr, zn, zr,
 * pn, pr, pperm, S0; NaN but S0 for a constant column. CLAB_EDOM for S0 =
 * 0. O(nnz log nnz + k nsim nnz).
 *
 * clab_bvmor: bivariate Moran (Wartenberg 1985,
 * doi:10.1111/j.1538-4632.1985.tb00849.x), n >= 3: n / S0 sum w_ij zx_i
 * zy_j / sqrt(sum zx^2 sum zy^2), with y's values permuted nsim times from
 * seed. Out, 6: I, the mean and sd of the permuted values, z, (1 + #{I* >=
 * I}) / (nsim + 1), (1 + min(#{I* >= I}, #{I* <= I})) / (nsim + 1) (NaN for
 * nsim = 0). CLAB_EDOM for a constant x or y, or S0 = 0.
 *
 * clab_sarml: spatial regression of y on X (n x p, column-major, p < n),
 * W the weights. model 0: y = rho W y + X b + e; 1: y = X b + u, u = lam W
 * u + e; 2: both. method 0: maximum likelihood (Anselin 1988,
 * doi:10.1007/978-94-015-7799-1), the likelihood concentrated in rho or lam
 * and maximised by Brent's method on (-1/m, 1/m), m the largest row sum of
 * W, log |I - a W| by dense LU; model 2 by Nelder-Mead in (rho, lam); the
 * covariance from the inverse of the expected information matrix. method 1
 * (model 0 only): two-stage least squares with instruments [X, inst]
 * (inst: n x ni, column-major; Kelejian and Prucha 1998,
 * doi:10.1023/A:1007707430416), s2 = RSS / n, the covariance with RSS / (n -
 * p - 1). Out: coef, 4 p + 4: b, rho,
 * lam (NaN when not in the model), their standard errors (p + 2), then
 * OLS b and its standard errors (p each); vc, (p + 2) x (p + 2) row-major
 * over (b, rho, lam), NaN outside the model; st, 5: log-likelihood (NaN for
 * 2SLS), s2, AIC, the OLS log-likelihood, m. CLAB_EDOM for a rank-deficient
 * X (after filtering) or a singular I - a W. O(n^3) per likelihood.
 *
 * clab_spill: direct, indirect and total effects of a lag model (LeSage and
 * Pace 2009, doi:10.1201/9781420064254) for n >= 2 units: with S = (I - rho
 * W)^-1, direct_j = tr(S (b_j I + th_j W)) / n, total_j = 1' S (b_j I +
 * th_j W) 1 / n (Lean: Crimlab.Hotspot.total_effect), indirect the
 * difference; th NULL for none. vc: NULL, or the covariance (1 + k or 1 +
 * 2k square, row-major) of (rho, b, th), for delta-method standard errors
 * with dS / drho = S W S. Out, 6 per coefficient: direct, indirect, total,
 * their standard errors (NaN without vc). CLAB_EDOM when I - rho W is
 * singular. O(n^3). */
int clab_gmoran(const double *x, int n, int k, const int *rp, const int *ci,
                const double *w, int nsim, unsigned int seed, double *out);
int clab_bvmor(const double *x, const double *y, int n, const int *rp,
               const int *ci, const double *w, int nsim, unsigned int seed,
               double *out);
int clab_sarml(const double *y, const double *x, int n, int p, const int *rp,
               const int *ci, const double *w, int model, int method,
               const double *inst, int ni, double *coef, double *vc,
               double *st);
int clab_spill(int n, const int *rp, const int *ci, const double *w, double rho,
               const double *b, const double *th, int k, const double *vc,
               double *out);

/* Point patterns and segregation (src/ppx.cpp).
 *
 * clab_crossk: the cross-type K function of n >= 2 points with types mk
 * (codes from 0) inside the rectangle win: K_ij(r) = A / (n_i n_j) sum over
 * a of type ti and b != a of type tj with d_ab <= r of Ripley's isotropic
 * weight at a (Ripley 1977, doi:10.1111/j.2517-6161.1977.tb01615.x),
 * n_i (n_i - 1) in the divisor when ti = tj; r: nr distances from 0,
 * increasing; out k, nr values. With sigma > 0, the segregation test of
 * Diggle, Zheng and Durr (2005, doi:10.1111/j.1467-9876.2005.05373.x): T =
 * sum_a sum_m (p_m(a) - pbar_m)^2, p_m(a) the leave-one-out share of type m
 * under a normal kernel of sd sigma, and the p-value of nsim relabellings
 * from seed; seg, 2: T, p (NaN for nsim = 0). CLAB_EDOM for a point outside
 * win, a type without points, or a point whose kernel weights all vanish.
 * O(n^2 (1 + nsim types)).
 *
 * clab_qdrtst: the counts of n points in nx x ny equal tiles of win (cnt:
 * ny x nx column-major, rows from the top; a point on the far edge in the
 * last tile) and Pearson's X2 against the mean count. st, 3: X2, nx ny - 1,
 * the points per unit area. O(n + nx ny).
 *
 * clab_morsta: Morisita's index of n >= 2 counts (from 0, total S > 1), I
 * = n (sum x^2 - S) / (S^2 - S), the uniform and clumped indices (chiu -
 * n + S) / (S - 1) and (chil - n + S) / (S - 1) at the chi-square quantiles
 * given, Smith-Gill's standardised index, and chi = I (S - 1) + n - S.
 * Out, 5: I, clumped, uniform, standardised, chi. O(n).
 *
 * clab_dissix: segregation of g >= 2 groups over n units (x: n x g
 * counts, column-major, every group with a positive total; units with no
 * one left out). d: g x g Duncan's D; grp, 4 per group against the rest:
 * IS, Theil's H, Atkinson's index of shape delta in (0, 1), Gini; with
 * weights rp, ci, w (NULL for none), dm: g x g Morrill's adjusted D (D
 * less the weighted mean |z_i - z_j| over neighbouring units, z the pair's
 * share; CLAB_EDOM when no weight joins two occupied units). Lean:
 * Crimlab.Hotspot.dissim_le_one. O(g^2 n + g n^2). */
int clab_crossk(const double *x, const double *y, const int *mk, int n, int ti,
                int tj, const double *win, const double *r, int nr,
                double sigma, int nsim, unsigned int seed, double *k,
                double *seg);
int clab_qdrtst(const double *x, const double *y, int n, const double *win, int nx,
                int ny, double *cnt, double *st);
int clab_morsta(const double *x, int n, double chiu, double chil, double *out);
int clab_dissix(const double *x, int n, int g, double delta, const int *rp,
                const int *ci, const double *w, double *d, double *grp,
                double *dm);

/* Geographic utilities (src/geou.cpp).
 *
 * clab_pinpl: the region of each of n points (px, py) among np polygons
 * given by nv >= 3 vertices (vx, vy) with polygon codes pid (0 to np - 1)
 * and ring codes rid, sorted by polygon then ring (CLAB_EARG otherwise),
 * each ring closed from its last vertex to its first. Even-odd rule over
 * all rings of a polygon (so holes and parts), a point on an edge inside
 * (Shimrat 1962, doi:10.1145/368637.368653). Out: reg, the first polygon
 * holding the point or -1; nreg, how many hold it. O(n nv) at worst,
 * polygons whose bounding box misses the point skipped.
 *
 * clab_arwt: an areal crosswalk. ns source values val, np pairs (from[q]
 * in [0, ns), to[q] in [0, nt)) with overlap sizes area (from 0): each
 * source is split among its pairs in proportion to area (evenly when its
 * areas sum to 0). Out: out (nt target totals), share (per pair), main (per
 * source: the target of its largest share, -1 without pairs), *lost (the
 * total of sources without pairs). O(ns + np). */
int clab_pinpl(const double *px, const double *py, int n, const double *vx,
               const double *vy, const int *pid, const int *rid, int nv, int np,
               int *reg, int *nreg);
int clab_arwt(const double *val, int ns, const int *from, const int *to,
              const double *area, int np, int nt, double *out, double *share,
              int *main, double *lost);

/* The buffer design of place-based evaluations (src/buf.cpp). For n
 * places, the distances ds and dc to the nearest of ns sites (sx, sy) and
 * nc comparison points (cx, cy); grp 0 (treated) when ds <= r, 1
 * (comparison) when dc <= r and ds > buf, 2 (excluded) otherwise; r > 0,
 * buf >= 0 (Bowers and Johnson 2003, doi:10.1023/A:1024909009240). O(n (ns
 * + nc)). */
int clab_bufas(const double *x, const double *y, int n, const double *sx,
               const double *sy, int ns, const double *cx, const double *cy,
               int nc, double r, double buf, int *grp, double *ds, double *dc);

/* ------------------------------------------------------------------------
 * Hawkes processes in time (hwk.cpp).
 * ------------------------------------------------------------------------ */

/* The model of the next three: n event times t, sorted, in [0, hor]; kind
 * 0 exponential kernel (rate beta), 1 Weibull (shape, scale), 2 Lomax
 * (shape, scale: density shape / scale (1 + d / scale)^-(shape + 1)); base
 * 0 a constant baseline e^a0, 1 exp(a0 + a1 t / hor + a2 sin(2 pi t / per)
 * + a3 cos(2 pi t / per)) (per > 0), its integral by Simpson's rule on
 * steps of at most per / 16. th: the baseline terms, eta >= 0, then the
 * kernel parameters (> 0): 3, 4, 6 or 7 values. Events at one instant do
 * not excite each other. CLAB_EARG for a bad size, window, option or
 * order of times; CLAB_EDOM for parameters outside their range. */

/* The log-likelihood at th and, when grad is not NULL, its gradient. O(n)
 * for the exponential kernel (Ozaki 1979), O(n^2) for the others. */
int clab_hwkll(const double *t, int n, double hor, int kind, int base, double per,
               const double *th, double *ll, double *grad);

/* The maximum-likelihood fit by BFGS with the analytic gradient over the
 * baseline terms, logit(eta) and the logarithms of the kernel parameters,
 * at most maxit steps per start, from start (th's layout, 0 < eta < 1) or,
 * when NULL, from two starts for the exponential kernel and two more (one
 * at the exponential fit) for the others, the best kept. n >= 2. Out: th;
 * cov, the inverse of the numerical Hessian of the gradient (NaN when not
 * positive definite); st (3): the log-likelihood, the steps and 1 when the
 * gradient vanished. CLAB_EDOM when no start gives a finite likelihood. */
int clab_hwkfit(const double *t, int n, double hor, int kind, int base, double per,
                const double *start, int maxit, double *th, double *cov, double *st);

/* Time-rescaling residuals: tau (n), the compensator's increments from 0;
 * ks (2), the Kolmogorov-Smirnov statistic of 1 - e^-tau against the
 * uniform and its p-value (exact below 100 events without ties). work: n
 * doubles, 3 (2 n + 1)^2 when n < 100. */
int clab_hwkres(const double *t, int n, double hor, int kind, int base, double per,
                const double *th, double *tau, double *ks, double *work);

/* Simulation on [0, hor] by the branching construction, from seed: the
 * background by thinning, Poisson(eta) offspring per event at lags drawn
 * from the kernel. win: NULL, or a rectangle (xmin, xmax, ymin, ymax) in
 * which the background falls uniformly (th's baseline is then the
 * window's whole rate) and offspring land at normal offsets of standard
 * deviation sig > 0; those outside it are dropped. Out, in time order: t,
 * x and y (with win), gen (0 for the background), n of them. CLAB_EDOM
 * when more than cap events arise. work: 3 cap doubles; iwork: 2 cap
 * ints. */
int clab_hwksim(double hor, int kind, int base, double per, const double *th,
                const double *win, double sig, unsigned int seed, int cap, double *t,
                double *x, double *y, int *gen, int *n, double *work, int *iwork);

/* ------------------------------------------------------------------------
 * Series: trend, change points, surveillance, cross-correlation (trn.cpp).
 * ------------------------------------------------------------------------ */

/* Mann-Kendall's test of a monotone trend in y over strictly increasing x
 * (n >= 3, finite) and the Theil-Sen slope (Sen 1968). alt: 0 two-sided,
 * 1 increasing, 2 decreasing. z: the normal quantile of the slope
 * interval. Out, st (10): S; its variance with the tie correction; tau
 * (S over n (n - 1) / 2); the p-value, exact over every ordering of the
 * observed values (ties kept) when n <= 8, otherwise normal with the
 * continuity correction; the slope (the median of the pairwise slopes);
 * the intercept (the median of y - slope x); the interval's ends at the
 * ranks floor((N - z sd) / 2) and ceil((N + z sd) / 2) + 1 of the N
 * sorted slopes, clamped to 1..N; and 1 where a clamp was needed. work:
 * n (n + 1) / 2 doubles; iwork: n ints. CLAB_EDOM for a non-finite value
 * or x not increasing. */
int clab_trndrr(const double *x, const double *y, int n, int alt, double z, double *st,
                double *work, int *iwork);

/* A single change of level in y (n finite values): the split after k
 * (minseg <= k <= n - minseg) with the largest T_k = k (n - k) / n
 * (mean before - mean after)^2, its p-value over the reorderings of y
 * (all n! of them when that is at most nperm, else nperm random ones from
 * seed: (1 + #{T* >= T}) / (1 + nperm)), and Pettitt's (1979) rank test.
 * Out, st (12): k, the means before and after, T, p, the orderings used,
 * 1 when exact; Pettitt's t, K = max |U_t|, its p-value 2 exp(-6 K^2 /
 * (n^3 + n^2)) (at most 1) and the means before and after t. work: n
 * doubles; iwork: n ints. */
int clab_stpcp(const double *y, int n, int minseg, int nperm, unsigned int seed, double *st,
               double *work, int *iwork);

/* CUSUM of counts y (>= 0) against expected counts m (> 0): z, the
 * transformed values (trans 0 (y - m) / sqrt(m), 1 Rossi's, 2 Anscombe's,
 * 3 y itself); s_i = max(0, s_(i-1) + z_i - k), alarm where s >= h > 0,
 * s set to 0 after an alarm when reset. */
int clab_cusrv(const double *y, const double *m, int n, double k, double h, int trans,
               int reset, double *z, double *s, int *alarm);

/* The cross-correlation of two series of n >= 2 finite values, neither
 * constant, at lags -lag..lag (0 <= lag < n): r[lag + k] estimates
 * cor(x[t + k], y[t]) as stats::ccf() does (sums over the overlap, the
 * full-series means and variances). */
int clab_xcser(const double *x, const double *y, int n, int lag, double *r);

/* ------------------------------------------------------------------------
 * Durations: Kaplan-Meier by group and the log-rank test (kmr.cpp).
 * ------------------------------------------------------------------------ */

/* time: n durations >= 0; event: 1 for an event, 0 for a censored
 * duration; grp: groups in [0, g). z: the normal quantile of the
 * intervals. Out, one row per group and distinct time (m of them, at most
 * n), groups in order: gout, tt, the number at risk nr, events nd,
 * censored nc, the survival s, its Greenwood standard error se and the
 * log-scale interval lo, hi (survival::survfit's "log"), NaN once s is 0;
 * med (2 g; the first g read), each group's median (survfit's rule: where
 * the curve sits at exactly 1/2, halfway to its next fall), NaN when it
 * never reaches 1/2; lr (2), the log-rank chi-square and its degrees of
 * freedom g - 1 (NaN for one group or a singular variance). work: 3 g +
 * (g - 1)^2 doubles; iwork: n ints. */
int clab_kmrec(const double *time, const int *event, const int *grp, int n, int g, double z,
               int *m, int *gout, double *tt, double *nr, double *nd, double *nc, double *s,
               double *se, double *lo, double *hi, double *med, double *lr, double *work,
               int *iwork);

/* ------------------------------------------------------------------------
 * Police operations research (hyq.cpp).
 * ------------------------------------------------------------------------ */

/* Larson's (1974) hypercube queue: n units (1 to 10) serving m atoms with
 * calls at rates lam (> 0), each answered by the first free unit in its
 * atom's order pref (m x n, column-major, row j a permutation of the
 * units 0..n-1), service at rate mu > 0. wait 0: a call finding every
 * unit busy is lost; 1: it waits (the infinite line, lam's total below n
 * mu), the queue lumped into the all-busy state. Out: pst (2^n), the
 * stationary probability of each set of busy units (bit u for unit u); wl
 * (n), each unit's workload; disp (n x m), the share of all calls from
 * atom j answered at once by unit u; st (3): the loss or waiting
 * probability, the mean queue (0 for loss) and rho = lam's total / (n
 * mu). work: 4^n + 2^n doubles. Exact by Gaussian elimination. */
int clab_hypq(const double *lam, const int *pref, int m, int n, double mu, int wait,
              double *pst, double *wl, double *disp, double *st, double *work);

/* Kolesar and Blum's (1973) square-root law: the mean distance from a
 * call to the nearest available unit, c sqrt(area / (units (1 - busy))),
 * with c = 1/2 for straight-line travel (metric 0) and sqrt(2 pi) / 4 on a
 * rectangular street grid (metric 1). area, units > 0; busy in [0, 1). */
int clab_sqrtlw(double area, double units, double busy, int metric, double *d, double *c);

/* ------------------------------------------------------------------------
 * Hot-spot dynamics: the Short-D'Orsogna-Brantingham model (sdb.cpp).
 * ------------------------------------------------------------------------ */

/* The model of Short et al. (2008) on an n x n torus (n >= 3, sites
 * column-major), steps steps of length dt. par (10): a0 > 0 (static
 * attractiveness), omega >= 0 (decay; omega dt <= 1), eta in [0, 1]
 * (spread), gam >= 0 (offenders arriving per site per unit time), theta
 * >= 0 (attractiveness added per burglary), dt > 0, then for mode 1 dx >
 * 0, dif > 0, eps >= 0, z > 0. b (n^2, in and out): the dynamic
 * attractiveness. mode 0, agents: rho (n^2) holds whole numbers of
 * offenders; each burgles with probability 1 - e^-(a0 + b) dt and leaves,
 * or moves to a neighbour chosen in proportion to its attractiveness;
 * b' = ((1 - eta) b + eta / 4 (sum of the four neighbours)) (1 - omega
 * dt) + theta E; Poisson(gam dt) offenders arrive at each site; from
 * seed. mode 1, the continuum limit by explicit Euler steps: b' = b + dt
 * (eta dif / z lap b - omega b + eps dif rho A), rho' = rho + dt (dif / z
 * (lap rho - 2 div(rho grad A / A)) - rho A + gam), A = a0 + b, no crime
 * where sup (NULL or n^2 flags) is set; CLAB_EDOM when dt 4 dif max(1,
 * eta) / z > dx^2 or a field turns negative. Out: burg and meanb (steps),
 * the burglaries (the crime rate times dt in mode 1) and the mean of b
 * after each step; tot (n^2), the burglaries at each site. work: 3 n^2
 * doubles; iwork: n^2 ints. */
int clab_sdbsim(int n, int steps, int mode, const double *par, double *b, double *rho,
                const int *sup, unsigned int seed, double *burg, double *meanb, double *tot,
                double *work, int *iwork);

/* Linear stability of the uniform state of the continuum model above
 * (Short, Bertozzi and Brantingham 2010). par (7): a0, omega, eps >= 0,
 * dif, gam, eta >= 0, z (the others > 0). At wavenumber-squared k the
 * linearisation has trace tr (negative at every k) and determinant det =
 * P k^2 - Q k + R; a mode grows exactly when det < 0, so the uniform
 * state forms hot spots when Q > 0 and Q^2 > 4 P R (Q > 0 when P = 0).
 * Out, st (12): P, Q, R; the steady B, A = a0 + B, rho; the two sides of
 * Q > 0, 3 eps dif rho and omega + eta A; 1 when unstable; the critical
 * k = Q / (2 P), det and the growth rate there (NaN unless P, Q > 0);
 * for the nk >= 0 values k (>= 0): tr, det and gr, the larger real part
 * of the eigenvalues. */
int clab_sdbthr(const double *par, const double *k, int nk, double *st, double *tr, double *det,
                double *gr);

/* ------------------------------------------------------------------------
 * Hawkes processes with marks, types and places (mhk.cpp).
 * ------------------------------------------------------------------------ */

/* Each of the three: n event times t, sorted, in [0, hor]; maxit 0
 * evaluates at start, otherwise BFGS from start or, when NULL, from the
 * kernel's starts, the best kept. Out: th; cov, the inverse of the
 * numerical Hessian of the analytic gradient (NaN when not positive
 * definite); grad (NULL or as th), the gradient of the log-likelihood at
 * th; st: the log-likelihood, the steps, 1 when the gradient vanished,
 * then as stated. CLAB_EARG for bad sizes or times; CLAB_EDOM for a start
 * outside the parameter space or no finite likelihood. */

/* Temporal ETAS (Ogata 1988): lambda(t) = mu + K sum over t_j < t of
 * e^(alpha dm_j) (t - t_j + c)^-p, dm the marks less the reference
 * magnitude; th (5): mu > 0, K >= 0, alpha (0 and not estimated when
 * afree is 0), c > 0, p > 1. st (4): ..., the branching ratio K mean(e^(alpha
 * dm)) c^(1 - p) / (p - 1). O(n^2). */
int clab_etasf(const double *t, const double *dm, int n, double hor, int afree,
               const double *start, int maxit, double *th, double *cov, double *st,
               double *grad);

/* Multivariate Hawkes with exponential kernels (Hawkes 1971): k types
 * (1 to 10), ty in [0, k); lambda_c(t) = mu_c + beta sum over l of
 * A[c, l] sum over earlier events of type l of e^(-beta (t - t_j)). th:
 * mu (k), A (k x k, column-major, A[c + l k] the events of type c each
 * event of type l triggers), then beta unless beta > 0 is given. st (4):
 * ..., the spectral radius of A (below 1 for a stationary process). O(n
 * k). work: 2 k + 8 np + np^2 doubles (np the parameters); iwork: k
 * ints. */
int clab_mvhwk(const double *t, const int *ty, int n, int k, double hor, double beta,
               const double *start, int maxit, double *th, double *cov, double *st,
               double *grad, double *work, int *iwork);

/* Space-time Hawkes (Mohler et al. 2011): lambda(t, x, y) = mu + eta
 * beta sum over earlier events of e^(-beta (t - t_j)) phi_sigma(x - x_j,
 * y - y_j), phi the isotropic normal density; mu per unit area and time
 * over a region of area area (the spatial kernel's mass outside it is
 * not subtracted). th (4): mu, eta (< 1 when fitted), beta, sigma. st
 * (3). qlam: the intensity at the nq points (qt, qx, qy) at th. O(n^2). */
int clab_sthwk(const double *t, const double *x, const double *y, int n, double hor,
               double area, const double *start, int maxit, double *th, double *cov,
               double *st, double *grad, const double *qt, const double *qx, const double *qy,
               int nq, double *qlam);
/* ---- Custody, courts and published tables (crsev, zcomp, ...) ---- */

/* Severity-weighted crime index (crsev): x[i] counts of offence off[i]
 * (0..no-1) in cell key[i] (0..nk-1), n rows (several rows of one cell
 * and offence add up); w: no weights from 0. pop: NULL or nk positive
 * populations; per > 0. Out, nk each: cnt (the counts summed), sev (the
 * weighted sum), rate and srate (cnt and sev per per of pop; NaN without
 * pop) and idx (srate, or sev without pop, as 100 times its ratio to cell
 * base's; NaN when base is -1 or the base cell's is 0). */
int clab_crsev(const double *x, const int *off, const int *key, int n, int no,
               int nk, const double *w, const double *pop, double per, int base,
               double *cnt, double *sev, double *rate, double *srate, double *idx);

/* Composite index of standardised counts (zcomp): x, nu x nc finite
 * (column-major, nu >= 2 units, nc categories); w: nc finite weights.
 * Out: mu, sd (nc), each category's mean and standard deviation (n - 1
 * denominator; 1 for a category that does not vary); z (nu x nc), (x -
 * mu) / sd; comp (nu), the weighted sum of each unit's z. */
int clab_zcomp(const double *x, int nu, int nc, const double *w, double *z,
               double *comp, double *mu, double *sd);

/* Concentration of a total across units (cconc; Lean:
 * Crimlab.Concentration.gini_zero_decomposition). x: n values from 0 with a
 * positive total; w: NULL or n positive multiplicities (a unit of weight
 * w counts as w units); fr: nf shares in [0, 1]. o: n ints of scratch, out
 * the units' order by value. Out: st (16): W (the units), T (the total),
 * the mean m, the variance (W - 1 denominator; NaN for W <= 1), the
 * dispersion index, the implied variance of the intensity max(v - m, 0)
 * and its root, the share of units at 0, the Poisson null share exp(-m),
 * their difference, the Gini coefficient, the Gini of the units above 0,
 * the null Gini exp(-m) + (1 - exp(-m)) times that, the variance with
 * denominator W, the bound sqrt(that / 2) / m and the Poisson bound 1 /
 * sqrt(2 m); lp, lv (n + 1): the Lorenz curve from (0, 0); top, topn (nf):
 * the share of the total held by the largest ceil(f W) units and that
 * number of units. */
int clab_cconc(const double *x, const double *w, int n, const double *fr, int nf, double *st,
               double *lp, double *lv, double *top, double *topn, int *o);

/* The Gini coefficient of n independent Poisson(mu) counts, nsim times,
 * from the generator seeded with seed (a whole number from 0); 0 when a
 * draw has no events. Out: g (nsim). work: n ints. */
int clab_cnull(double mu, int n, int nsim, double seed, double *g, int *work);

/* Distinct places among the first i events (cconc): code[i] the place of
 * event i (0..ng-1), n >= 2 events. Out: dist (n), the distinct places
 * so far; mhat, the M with M sum_{j < n} 1/(M + j) equal to the final
 * count (Inf when every event is at a new place, 0 when all are at one);
 * ex, lo, hi (n): the expected count M sum_{j <= i} 1/(M + j) under
 * preferential allocation and its bounds M log(1 + (i+1)/M) and 1 + M
 * log(1 + i/M) (NaN unless 0 < M < Inf); slope, the least-squares slope of
 * log dist on log i over the second half (NaN for fewer than 3 points).
 * seen: ng ints of scratch. */
int clab_cgrow(const int *code, int n, int ng, double *dist, double *mhat, double *ex, double *lo,
               double *hi, double *slope, int *seen);

/* Urban scaling (cconc): the least-squares fit of log x on log pop over
 * the units with x above 0 (x from 0, pop positive). Out (5): the exponent
 * beta, the prefactor exp(intercept), beta's standard error, R^2 and the
 * units used; CLAB_EDOM, with the count set, for fewer than 3 or no spread
 * in pop. */
int clab_urbsc(const double *x, const double *pop, int n, double *out);

/* The Hurwitz zeta function sum_{k >= 0} (q + k)^-s, s > 1, q > 0, by
 * Euler-Maclaurin summation. */
int clab_hzta(double s, double q, double *z);

/* Power-law tail index of the values x >= xmin (hilltl): x n positive
 * values, xmin > 0, at least 2 in the tail and not all equal (else
 * CLAB_EDOM, with out[3] set). discrete && !approx: the maximum
 * likelihood exponent of the discrete power law, P(X = k) = k^-a /
 * zeta(a, xmin), on [1.0001, 25]; otherwise the continuous Hill estimator
 * 1 + n / sum log(x / d), d = xmin - 1/2 for counts (xmin > 1/2) and xmin
 * else. Out (5): alpha, its standard error (the observed information by a
 * central difference, or (alpha - 1) / sqrt(n)), the Kolmogorov-Smirnov
 * distance of the fit, the tail size and 1 when the maximum lies at an
 * end of the search interval. work: n doubles. */
int clab_hilltl(const double *x, int n, double xmin, int discrete, int approx, double *out,
                double *work);

/* Capture-recapture (crcap; Lean: Crimlab.Dark.petersen_bounds). k = 3:
 * x = (n1, n2, m), two lists and their overlap, 0 < m <= min(n1, n2);
 * kappa >= 1 bounds the dependence of the lists (the odds of being on one
 * list given the other, over those given not). Out: st (6), Petersen's n1
 * n2 / m, its bounds over the dependence box (divided and multiplied by
 * kappa), Chapman's estimate and Seber's standard error of it, and the
 * distinct floor n1 + n2 - m. k = 7: x the cells 100, 010, 001, 110, 101,
 * 011, 111 of three lists (counts from 0, the three pairwise cells
 * positive). Out: st (3), the observed total, the never-listed cell
 * under no three-way interaction (m100 m010 m001 m111 / (m110 m101 m011))
 * and the total; pw (18): for the pairs 1-2, 1-3, 2-3 each list's total,
 * the overlap, the floor, Petersen and Chapman (six blocks of three);
 * sens (2 ni): for each three-way log interaction it[j] the missing cell
 * times exp(-it[j]) and the total. */
int clab_crcap(const double *x, int k, double kappa, const double *it, int ni, double *st,
               double *pw, double *sens);

/* Bounds on a true victimisation rate (dkfig; Lean:
 * Crimlab.Dark.rate_bounds). v: n survey rates in [0, 1], r: n recorded
 * rates in [0, 1]; the survey's false positives at most amax and false
 * negatives at most bmax (amax + bmax < 1). Out, n each: lo = max(r, (v -
 * amax) / (1 - amax)), hi = max(lo, min(1, v / (1 - bmax))), and brk, the
 * largest bmax under which the true rate stays below thr, 1 - v / thr
 * (thr NaN: NaN; otherwise thr > v > 0). */
int clab_dkfig(const double *v, const double *r, int n, double amax, double bmax, double thr,
               double *lo, double *hi, double *brk);

/* The dark figure when a survey covers only some offences (dkfig): va
 * survey offences of the covered type A (positive), ra of them recorded
 * (0 < ra <= va), rb recorded offences of the uncovered type B, recorded
 * at no less than lam (0 < lam <= 1) times A's rate. Out (7): A's
 * recording rate, the overall recorded share's upper bound (ra + rb) /
 * (va + rb), the total's bounds va + rb and va + rb / (lam ra / va), the
 * dark figure's bounds and their width. */
int clab_dkcmp(double va, double ra, double rb, double lam, double *out);

/* The hierarchy rule (dkfig): k[i] offences in incident i (whole numbers
 * from 1), n incidents. Out (6): incidents, offences, offences per
 * incident, extra offences per incident, and the bounds n and max(k) n on
 * the offences behind n recorded incidents. */
int clab_dkhier(const double *k, int n, double *out);

/* Values for banded categories (bandx): n bands with lower and upper
 * ends lo, hi (NaN where a label could not be read), olo[i] / oup[i]
 * non-zero for a band open below / above. An open top band ends at cap,
 * or (cap NaN) at fac > 0 times its start; an open bottom band starts at
 * flr, or (flr NaN) at min(0, its end); CLAB_EDOM when that puts its ends
 * in the wrong order. rule: 0 the midpoint, 1 the lower end, 2 the upper
 * end, 3 the geometric mean (the midpoint where the lower end is not
 * above 0). Out, n each: val, and the ends used blo, bhi. */
int clab_bandx(const double *lo, const double *hi, const int *olo, const int *oup, int n, int rule,
               double cap, double fac, double flr, double *val, double *blo, double *bhi);

/* Bounds of published cells (pbbnd): val, n numbers from 0 (NaN for
 * none); st: 0 a number, 1 suppressed, 2 missing; lim: for a suppressed
 * cell its "<k" limit, NaN to use sup. rnd: NaN, or the rounding base:
 * kind 0 to the nearest multiple (+- rnd / 2), 1 random rounding (+- (rnd
 * - 1)); the lower end is floored at 0. A suppressed cell lies in [0, k -
 * 1] (CLAB_EDOM without a positive limit). Out, n each: lo, hi, and sto:
 * 0 exact, 1 rounded, 2 suppressed, 3 missing (lo, hi NaN). */
int clab_pbbnd(const double *val, const int *st, const double *lim, int n, double rnd, int kind,
               double sup, double *lo, double *hi, int *sto);

/* The change between two published values known only within [nlo, nhi]
 * (now) and [plo, phi] (before), n of them (pbbnd, yoypx). pop, ppop: NULL
 * or n positive populations; per > 0. Out (8 n, n each): the change's
 * bounds nlo - phi and nhi - plo; the percent change's, 100 (nlo / phi -
 * 1) and 100 (nhi / plo - 1) (NaN at a zero denominator); the rate's per
 * per of pop now and before (NaN without a population). */
int clab_pbchg(const double *nlo, const double *nhi, const double *plo, const double *phi, int n,
               const double *pop, const double *ppop, double per, double *out);

/* Period over period change (yoypx): v, n values; prd, their finite
 * periods; grp, their groups; rows sorted by group, then period, a period
 * at most once in a group (CLAB_EARG otherwise). Out, n each: prev, the
 * value at period prd - lag (lag > 0) in the same group (NaN for none);
 * chg = v - prev; pct = 100 chg / prev (NaN for a zero prev). */
int clab_yoypx(const double *v, const double *prd, const int *grp, int n, double lag, double *prev,
               double *chg, double *pct);

/* Adjusted p-values (yoypx), as p.adjust(): p, n values in [0, 1]; meth 0
 * none, 1 Bonferroni, 2 Holm, 3 Hochberg, 4 Benjamini-Hochberg, 5
 * Benjamini-Yekutieli. Out: q (n). o: n ints of scratch. */
int clab_pbadj(const double *p, int n, int meth, double *q, int *o);

/* The judge-leniency instrument on observed data (jdgiv; Lean:
 * Crimlab.Judge.wald_with_defiers): z, d n values 0 or 1 (the lenient or
 * strict judge, the decision), y n finite outcomes, w NULL or n weights
 * from 0 with a positive total in each arm. Out: est (6), the first stage
 * E[d | z = 1] - E[d | z = 0], the intention-to-treat difference in y,
 * their ratio (Wald; NaN at a zero first stage), and the shares of
 * compliers, always-takers and never-takers if no one defies; sens (nd x
 * ne), the compliers' effect implied by a defier share dsh[a] with effect
 * deff[b]: (wald (c) + s e) / (c + s), c the complier share. */
int clab_jdgiv(const double *z, const double *d, const double *y, const double *w, int n,
               const double *dsh, int nd, const double *deff, int ne, double *est,
               double *sens);

/* The same on a known population (jdgiv): each case's decisions d0, d1
 * under the strict and the lenient judge (0 or 1) and outcomes y0, y1
 * without and with the decision. Out (12): the shares of compliers,
 * defiers, always- and never-takers, the mean effect y1 - y0 in each (NaN
 * for an empty type), the intention-to-treat effect, the first stage,
 * their ratio and the average effect. */
int clab_jdgpop(const double *d0, const double *d1, const double *y0, const double *y1,
                const double *w, int n, double *out);

/* Decisions by judge (jmono): jc the judge (0..nj-1) and d the decision (0
 * or 1) of n cases, y NULL or their outcomes, w NULL or weights from 0.
 * Out, nj each: cnt, the cases; one, the decisions 1; rate and outc, the
 * weighted means of d and y (NaN without y). st (4 + nj; the last nj
 * scratch): Pearson's chi-square of decision by judge on the counts (NaN
 * when every decision is the same or there is one judge), its degrees of
 * freedom, and the lowest and highest rate. CLAB_EDOM for a judge whose
 * weights sum to 0. */
int clab_jrate(const int *jc, const double *d, const double *y, const double *w, int n, int nj,
               double *cnt, double *one, double *rate, double *outc, double *st);

/* The many-judge slope test (jmono; Lean: Crimlab.Judge.slope_bound): nj
 * judges in ascending order of rate, outc their mean outcomes, outcomes in
 * [lo, hi]. For each pair a < b (in order (0,1), (0,2), ..., nj (nj - 1)
 * / 2 of them): dp and dy, the differences in rate and outcome; bnd, (hi -
 * lo) dp, the most nested judges allow; viol, |dy| > bnd; late, dy / dp
 * (NaN for dp = 0). */
int clab_jslope(const double *rate, const double *outc, int nj, double lo, double hi, double *dp,
                double *dy, double *bnd, int *viol, double *late);

/* Monotonicity of judges within groups of cases (jmono): jc, gc the judge
 * (0..nj-1) and the group (0..ng-1) and d the decision of n cases; jt the
 * judges' overall rates. Out: cn, cr (nj x ng), the cases and the rate of
 * each judge in each group; for each group sj, sc and sl (ng each), the
 * judges with at least minn cases there, their cases and the slope of
 * their group rates on their overall rates weighted by cases (NaN for
 * fewer than two judges or no spread); for each pair of those judges with
 * different overall rates, the stricter (higher overall rate) first: pg,
 * ps, pl the group and the two judges, and pv the pooled two-proportion z
 * of the lenient judge's group rate over the stricter's (0 when the
 * standard error is 0). With pg NULL only the number of pairs, npair, is
 * found. */
int clab_jsub(const int *jc, const int *gc, const double *d, int n, int nj, int ng, double minn,
              const double *jt, double *cn, double *cr, int *pg, int *ps, int *pl, double *pv,
              double *sl, int *sj, double *sc, int *npair);

/* Sentences against a mandatory minimum (mmshr): sent and mn, n finite
 * values; gc NULL or n groups 0..ng-1. Out (5 (ng + 1)): for each group
 * and then overall, the cases, the shares at, below and above the
 * minimum, and the mean excess of those above (0 for none). */
int clab_mmshr(const double *sent, const double *mn, const int *gc, int n, int ng, double *out);

/* The age-crime curve of a mixture of latent types (agemix; Lean:
 * Crimlab.LifeCourse.one_type_cov_zero): cv, na x nt finite curves (ages
 * by types, column-major); w, nt shares from 0 summing to 1 (within
 * 1e-8). Out: agg (na), the aggregate curve cv w; cov (na x na), the
 * covariance across types of the curves at two ages, sum_k w_k (cv_ak -
 * agg_a)(cv_bk - agg_b). */
int clab_agemix(const double *cv, int na, int nt, const double *w, double *agg, double *cov);

/* The panel test of one type (agemix): ya, yb, n >= 3 individuals'
 * finite outcomes at two ages. Out (2): the covariance mean(ya yb) -
 * mean(ya) mean(yb), and its two-sided permutation p-value (1 + #{|perm|
 * >= |obs|}) / (nperm + 1) from nperm shuffles of yb by the generator
 * seeded with seed. work: n doubles. */
int clab_agepan(const double *ya, const double *yb, int n, int nperm, double seed, double *out,
                double *work);

/* The cumulative risk of a first event by age (cumrsk): age, n increasing
 * ages; e, first events (from 0); pop, the population (positive); atrisk
 * NULL (the population less those with an event before) or those at risk;
 * dth NULL or death rates in [0, 1). Out, n each: h, the hazard e / nr; nr,
 * those at risk; cr, the cumulative risk 1 - prod (1 - h), or with deaths
 * sum of the surviving cohort's hazards; v, Greenwood's variance of log
 * survival (NaN with deaths). CLAB_EDOM when h > 1 or events and deaths
 * exceed the cohort. */
int clab_cumrsk(const double *age, const double *e, const double *pop, const double *atrisk,
                const double *dth, int n, double *h, double *nr, double *cr, double *v);

/* Rates between pipeline stages with a lag (jpipe): num, den, n counts
 * from 0 in time order; win >= 1, lag >= 0. Out, n each: onum, the sum of
 * num over the win periods ending at i; oden, the sum of den over the win
 * periods ending at i - lag (NaN until a full window); rate, their ratio
 * (NaN for a denominator of 0 or none). */
int clab_jpipe(const double *num, const double *den, int n, int win, int lag, double *onum,
               double *oden, double *rate);

/* Segregation placements under the Mandela Rules (mndla): days, n stays
 * (from 0); thr from 0, the days beyond which a stay is prolonged (Rule
 * 43: 15). With hours and miss (average hours out of cell from 0, the
 * share of days missing the 4 hours in [0, 100]) a stay meets Rule 44
 * when hours <= 2 and miss = 100: cls 2 when it also lasts more than thr
 * days, 1 when not, 0 otherwise (Sprott and Doob). Without them, cls 1
 * for a stay of more than thr days, or (nal non-NULL) with at least amin
 * alerts. excl[i] non-zero (contact that meets the rules) or elig[i] zero
 * sets cls 0. */
int clab_mndla(const double *days, const double *hours, const double *miss, const int *nal,
               int amin, const int *elig, const int *excl, int n, double thr, int *cls);

/* The Mandela classes by group (mndla): cls (0, 1, 2) and the group gc
 * (0..ng-1) of n stays, w NULL or their weights from 0, idc NULL or the
 * person of each. Out (7 (ng + 1), one row of 7 per group, then all):
 * the stays (weighted), those with cls > 0, cls 1 and cls 2; with idc,
 * the people, those with any stay of cls > 0, and those whose days add up
 * to more than thr (and, with elig, have a stay with elig non-zero). o: 2 n
 * ints of scratch. */
int clab_mndtab(const int *cls, const int *idc, const int *gc, const double *days, const double *w,
                const int *elig, int n, int ng, double thr, double *out, int *o);

/* Stock and flow over n periods (stkflw; Lean: Crimlab.Flow.adp_alos,
 * days_ratio): exactly two of pe (people admitted or released, positive),
 * dy (person-days, from 0), ad (average daily population, from 0) and al
 * (average length of stay, positive), but not dy with ad (CLAB_EARG
 * otherwise); t, n period lengths (positive). days = adp t = people alos
 * gives the rest. base 0: changes against the first period, 1: against the
 * one before. expo NULL or n positive populations, per > 0; rel NULL or n
 * releases (from 0) and early NaN or the days by which each is brought
 * forward. Out (16 n + 4, n each): people, days, adp, alos; their percent
 * changes (people, alos, days, adp); the share of the change in days on
 * the log scale due to people; flow and stock rates per per of expo and
 * their changes (NaN without expo); the stock an earlier release removes,
 * rel early / t, and that as a percentage of adp; then the means over the
 * periods of people, days, adp and alos. */
int clab_stkflw(const double *pe, const double *dy, const double *ad, const double *al,
                const double *t, int n, int base, const double *expo, double per,
                const double *rel, double early, double *out);

/* Lengths of stay (lostay): x, n stays from 0; w NULL or whole-number
 * frequencies from 0 with a positive total; br, nb increasing band edges
 * (bands (br[k], br[k + 1]]); thr, nt thresholds; flag NULL or n flags.
 * Out: st (9), the stays N, the days, the mean, the standard deviation (N
 * - 1), the median and quartiles (type 7 over the expanded stays), the
 * longest and the standard error of the mean; bt (2 (nb - 1)), the stays
 * in each band and the flagged ones; tt (nt), the share at or below each
 * threshold; fit (6), the AIC of the lognormal, exponential and Pareto
 * (from the shortest) fits of the positive stays by maximum likelihood
 * (NaN for fewer than 2, or no spread), the lognormal's mean and standard
 * deviation of log stay, and the Pareto exponent. o: n ints of scratch. */
int clab_lostay(const double *x, const double *w, int n, const double *br, int nb,
                const double *thr, int nt, const int *flag, double *st, double *bt, double *tt,
                double *fit, int *o);

/* Little's law on a docket (bklog; Lean: Crimlab.Flow.little): a, d the
 * arrival and disposition times of n cases (d NaN or after the horizon h >
 * 0: pending); every arrival in [0, h]. Out (7): the cases, those pending,
 * the filing rate n / h, the mean time to disposition of the disposed,
 * the time-average backlog (the occupancy over h, a pending case counted to
 * h), the occupancy, and with nothing pending the identity's residual L -
 * lambda W (NaN otherwise). */
int clab_bklog(const double *a, const double *d, int n, double h, double *out);

/* Case flow (dktflw): rc, ds the cases received and disposed in n
 * periods, pd NULL or those pending; win >= 1. Out (8 n, n each): the
 * clearance ds / rc (NaN for none received), the net rc - ds, the sums
 * over the win periods ending at each (NaN until a full window) of rc
 * and ds and their ratio, the change in pending, the periods to clear the
 * pending at the mean disposals of the last win periods (NaN without pd),
 * and the cumulative net. */
int clab_dktflw(const double *rc, const double *ds, const double *pd, int n, int win, double *out);

/* The mix of outcomes over periods (dktflw): x, nt periods by k outcomes
 * (column-major counts from 0, every period with a case); sc, nt scores;
 * a > 0 the chi-square quantile for Goodman's simultaneous intervals.
 * Out, nt x k each: sh, lo, hi; st (3): the chi-square of homogeneity
 * over the outcomes that occur, its degrees of freedom and Cramer's V;
 * tr (2 k): each outcome's Cochran-Armitage trend statistic (NaN for fewer
 * than 3 periods or an outcome always or never seen) and the change of its
 * share from the first period to the last. */
int clab_dkmix(const double *x, int nt, int k, const double *sc, double a, double *sh,
               double *lo, double *hi, double *st, double *tr);

/* The actuarial life table (lftab): k contiguous bands [lo, hi) (only the
 * last may end at Inf), d events and w censored cases in each (from 0);
 * z the normal quantile; ci 0 plain, 1 log, 2 log-log intervals. Out (8 k
 * + 1, k each): entering, at risk (w counted as half), the conditional
 * probability q, survival, its standard error (Greenwood), the interval's
 * ends (NaN at survival 0 or 1), the hazard, then the median (linear in
 * its band; Inf when survival stays above 1/2 or crosses it in an open
 * band). */
int clab_lftab(const double *lo, const double *hi, const double *d, const double *w, int k,
               double z, int ci, double *out);

/* Transitions between k states (trnmat): from, to the states of m
 * transitions. Out: cnt and p (k x k), the counts and the row
 * probabilities (rows of 0 for a state never left); pi (k), the
 * stationary distribution of p (NaN when a state has no exit or it is not
 * unique); st (3), the transitions, the share that stay in their state,
 * and the Theil index of the moves between states (NaN for none). work:
 * k (k + 1) doubles. */
int clab_trnmat(const int *from, const int *to, int m, int k, double *cnt, double *p, double *pi,
                double *st, double *work);

/* ------------------------------------------------------------------------
 * Geographic profiling and risk terrain: Canter's circle hypothesis, the
 * hit score of a profile, risk terrain layers.
 * ------------------------------------------------------------------------ */

/* Canter and Larkin's (1993) circle hypothesis for a series of n >= 2
 * crime sites (x, y): the circle whose diameter joins the two sites
 * farthest apart (the first such pair, i < j, in the order of the sites).
 * Out: st (4), the centre's x and y, the radius and the share of the
 * sites within it; pair (2), the two sites (from 0); for each of nh homes
 * (hx, hy; nh may be 0) mar, 1 when it lies within the circle (a
 * "marauder"), else 0. Within means a distance at most r + 1e-12 max(1,
 * r). */
int clab_crchyp(const double *x, const double *y, int n, const double *hx, const double *hy,
                int nh, double *st, int *pair, int *mar);

/* The hit score of a geographic profile (Rossmo 2000): s, the scores of n
 * search cells (higher searched first); for each of nh cells home[k] (from
 * 0), the cells scored above it (above), with the same score (ties, itself
 * included), hit = 100 above / n, hitt = 100 (above + ties) / n. */
int clab_hitscr(const double *s, int n, const int *home, int nh, int *above, int *ties,
                double *hit, double *hitt);

/* A risk terrain layer (Caplan, Kennedy and Miller 2011) at nc cells (cx,
 * cy) from nf features (fx, fy, nf >= 1) within distance rad > 0: kind 0,
 * 1 when a feature lies within rad (proximity); kind 1, the features
 * within rad (density); kind 2, the Gaussian kernel sum sum_f exp(-d^2 /
 * (2 rad^2)), standardised over the cells to mean 0 and standard deviation
 * 1 (all 0 when it does not vary, nc >= 2). Out: out (nc). */
int clab_rtmlay(const double *cx, const double *cy, int nc, const double *fx, const double *fy,
                int nf, double rad, int kind, double *out);

/* ------------------------------------------------------------------------
 * Forensic toxicology: calibration, left-censored values, postmortem
 * redistribution, ethanol and the reliability of specimen matrices.
 * ------------------------------------------------------------------------ */

/* A calibration line by weighted least squares of the responses y on the
 * concentrations x (x > 0) of n >= 3 calibrators with weights w > 0, and
 * the inverse prediction of nu unknown responses u (NaN stays NaN). Out:
 * st (6), the intercept, the slope, R^2 (weighted, as summary.lm gives
 * it), the residual standard error sigma (n - 2
 * degrees of freedom), the limit of detection klod sigma / |slope| and of
 * quantitation kloq sigma / |slope|; xh (nu), (u - intercept) / slope;
 * flag (nu), 0 below the limit of detection, 1 below that of
 * quantitation, 2 quantifiable, -1 for NaN. CLAB_EDOM when x does not
 * vary or the slope is 0. */
int clab_txcal(const double *x, const double *y, const double *w, int n, const double *u, int nu,
               double klod, double kloq, double *st, double *xh, int *flag);

/* Substitution of left-censored values: x (n; NaN for not measured, kept
 * as it is) below the limit of detection lod (nl = 1 or n, positive) is
 * replaced by lod / 2 (meth 0), lod / sqrt(2) (1) or lod (2). Out: out
 * (n), cens (n, 1 for a replaced value), st (2), the share of measured
 * values replaced (NaN when none is measured) and the number not
 * measured. */
int clab_lodimp(const double *x, int n, const double *lod, int nl, int meth, double *out,
                int *cens, double *st);

/* Postmortem redistribution from the central and peripheral (femoral)
 * blood concentrations c and p (n each, positive): r = c / p and its
 * class cls, 0 when r <= c1, 1 when r <= c2, 2 above (0 < c1 <= c2). */
int clab_pmrrat(const double *c, const double *p, int n, double c1, double c2, double *r,
                int *cls);

/* Ethanol in a decomposed body, antemortem or postmortem: eth (n, at least
 * 0), the higher alcohols pr and bu (n-propanol, n-butanol; NaN not
 * measured) and ethyl glucuronide etg (NaN not measured). Out: hi (n),
 * the larger higher alcohol (NaN when neither is measured); v (n), 0
 * antemortem ingestion (etg > 0), else 1 postmortem production (hi >
 * thr), else 2 indeterminate. */
int clab_etohlr(const double *eth, const double *pr, const double *bu, const double *etg, int n,
                double thr, double *hi, int *v);

/* The reliability of n specimen matrices after submersion and
 * decomposition: base (n, at least 0) the reliability of each in a fresh
 * body, dil (n, 1 when water dilutes it, as blood and urine). rel = base
 * (1 - min(stage / smax, scap)) times, for a diluted matrix, (1 -
 * min(days / dmax, dcap)); days, stage at least 0, dmax, smax positive,
 * dcap, scap from 0 to 1. Out: rel (n) and ord (n), the matrices from the
 * most reliable down (ties in their order). */
int clab_mtxrel(const double *base, const int *dil, int n, double days, double stage,
                double dmax, double dcap, double smax, double scap, double *rel, int *ord);

/* ------------------------------------------------------------------------
 * Taphonomy: compositional transforms, likelihood ratios of measured
 * evidence, the decay chain, a conjugate Bayesian preservation model.
 * ------------------------------------------------------------------------ */

/* Log-ratio transforms of n compositions of d >= 2 parts (x, n by d,
 * column-major, from 0; each row is closed to sum 1 after adding ps >=
 * 0 to every part, and must then hold no 0). type 0, the centred
 * log-ratio (out n by d): log x_j minus the row's mean log; type 1, the
 * isometric log-ratio of the default pivot basis (out n by d - 1):
 * sqrt(i / (i + 1)) (mean of log x_1..x_i - log x_{i+1}). */
int clab_clrilr(const double *x, int n, int d, double ps, int type, double *out);

/* The likelihood ratio of n measurements x under two normal models, H1
 * (means m1, standard deviations s1, k1 = 1 or n of each) and H2 (m2, s2,
 * k2). Out: out (5), the log-likelihoods under H1 and H2, their
 * difference, the likelihood ratio and its log10; band, the verbal band
 * of |log10 LR| (0 at 0, 1 up to 1, 2 up to 2, 3 up to 3, 4 up to 4, 5 up
 * to 6, 6 above); side, 1 when the ratio favours H1, 2 H2, 0 neither. */
int clab_preslr(const double *x, int n, const double *m1, const double *s1, int k1,
                const double *m2, const double *s2, int k2, double *out, int *band, int *side);

/* A decay chain of k >= 1 transient stages and two absorbing states
 * (skeletal, mummified): from each stage the body advances (to the next
 * stage, from the last to skeletal) with probability dr (1 - pres),
 * mummifies with probability mr pres (the two sum to at most 1), and
 * otherwise stays. pres, mr from 0 to 1, dr above 0 to 1;
 * CLAB_EDOM when nothing leaves a stage. Out: P ((k + 2) by (k + 2), the
 * transition matrix, column-major), N (k by k, the fundamental matrix
 * (I - Q)^-1), B (k by 2, the absorption probabilities, skeletal then
 * mummified), t (k, the expected steps to absorption). */
int clab_tpdmc(double pres, double dr, double mr, int k, double *P, double *N, double *B,
               double *t);

/* One path through a chain with m states of which the first k are
 * transient (P m by m, column-major, rows summing to 1), from state start
 * (< k), for at most steps steps or until absorbed, the draws from
 * xorshift64* seeded with seed. Out: path (steps + 1), len, the states
 * visited. */
int clab_tpdsim(const double *P, int m, int k, int start, int steps, unsigned long long seed,
                int *path, int *len);

/* A conjugate Bayesian linear model y = b0 + x b + e (n cases, q >= 1
 * columns of x, column-major; p = q + 1) with independent normal priors
 * (means m0, standard deviations s0 > 0, p each) and the noise variance
 * s2 fixed at the least-squares residual variance, rss / max(1, n - p)
 * (the variance of y when the fit is exact: a residual sum of squares at
 * most 1e-24 of the sum of y squared). Out: mu, sd (p), the posterior
 * means and standard deviations; st (2), s2 and tau2 (NaN without
 * groups). With g (n group codes from 0, ng >= 2 groups, each used; NULL
 * for none), from the residuals about the posterior means: gm, the mean
 * residual of each group, nj its cases, and lam = tau2 / (tau2 + s2 /
 * nj), with tau2 = max(0, var(gm) - mean(s2 / nj)) (ng each). CLAB_EDOM
 * when x has dependent columns (or fewer cases than columns) or y does
 * not vary. work: 2 p^2 + p doubles. */
int clab_tpbhm(const double *y, const double *x, int n, int q, const double *m0,
               const double *s0, const int *g, int ng, double *mu, double *sd, double *st,
               double *gm, double *lam, int *nj, double *work);

/* ------------------------------------------------------------------------
 * Deterrence theory: the offender's response to certainty and severity,
 * the identification of a sanction design, the inspection game, collider
 * bias in arrest data.
 * ------------------------------------------------------------------------ */

/* The offending an offender chooses from a menu: n >= 1 levels x (strictly
 * increasing), the benefit b of each and the sanction s (strictly
 * increasing). For each of np certainties p (at least 0) the levels
 * maximising b - p s: xo the largest of them (within 1e-12 of the
 * maximum), v the maximum. */
int clab_detrsp(const double *x, const double *b, const double *s, int n, const double *p,
                int np, double *xo, double *v);

/* Which partial effects a deterrence design identifies: the n regimes'
 * certainty p, severity s and celerity c (lg 1: log p and log s, both
 * positive) in the linear response y = b0 + b1 p + b2 s + b3 c. Out: rank,
 * the rank of the design (1, p, s, c) (a column whose part not explained
 * by those before it is below 1e-7 of its norm adds nothing); ident (3),
 * 1 when the coefficient of p, s, c is identified (full rank, or leaving
 * the column out lowers the rank); cst (3), 1 when the dimension does not
 * vary; with y (n, else NULL) and full rank, beta (4), the least-squares
 * coefficients. CLAB_EDOM for y with a rank below 4. work: 8 n doubles. */
int clab_dtrid(const double *p, const double *s, const double *c, int n, int lg, const double *y,
               int *rank, int *ident, int *cst, double *beta, double *work);

/* The criminal-citizen-inspector game of rock-paper-scissors type
 * (stake a, not 0): replicator dynamics x_i' = a x_i (x_{i+1} - x_{i+2})
 * from x0 (three shares from 0 summing to 1 within 1e-10), steps
 * fourth-order Runge-Kutta steps of size dt > 0. Out: traj ((steps + 1)
 * by 3, column-major), v (steps + 1), the invariant x1 x2 x3 along it; f0
 * (3), the field at x0; dv, dV/dt at x0 (0). */
int clab_inspg(double a, const double *x0, int steps, double dt, double *traj, double *v,
               double *f0, double *dv);

/* The inspection game of cooperators, criminals and inspectors over a
 * grid of temptations T (nt) and inspection costs g (ng): payoffs C (1, 0,
 * 1), P (T, 0, T - 1), O (1 - g, 1 - g, 1 - g) against (C, P, O); from x0
 * (three shares from 0, positive sum, normalised), steps discrete
 * replicator steps x += rate x (f - mean f), negative shares set to 0 and
 * the shares normalised. Out: out (3 nt ng, the
 * final shares of C, P and O, each nt by ng, column-major). */
int clab_insphs(const double *T, int nt, const double *g, int ng, int steps, double rate,
                const double *x0, double *out);

/* Collider bias of arrest data: two independent causes of arrest of
 * prevalences a and e (n each, strictly between 0 and 1) and a background
 * cause bg (from 0 to 1), arrest when any is present. Out: cells (4 n,
 * each n: the shares of the arrestees with both causes, a only, e only,
 * neither) and orr (n), the odds ratio of a and e among the arrestees,
 * which equals bg. */
int clab_colarr(const double *a, const double *e, const double *bg, int n, double *cells,
                double *orr);

/* ------------------------------------------------------------------------
 * Areas: crosswalks between boundary versions, synthetic exposures.
 * ------------------------------------------------------------------------ */

/* Values of n records x, each in source unit u (from 0, below ns; -1 for
 * a unit the crosswalk lacks, skipped), carried to nt target units by m
 * links (from[k] to to[k] with weight w[k] >= 0, the share of the source
 * unit's value that falls in the target: 1 for a source wholly inside its
 * target, an area or population share to split one). Out: out (nt), the
 * values of the targets; wsum (ns), each source's total weight (1 when
 * nothing is lost or counted twice); best (ns), the target of each
 * source's largest weight (the first of ties, -1 for none) and bw (ns),
 * that weight. */
int clab_bndxw(const double *x, const int *u, int n, const int *from, const int *to,
               const double *w, int m, int ns, int nt, double *out, double *wsum, int *best,
               double *bw);

/* Synthetic exposure of na areas from a logistic model fitted elsewhere
 * (clab_lgtreg): coef (p + 1, intercept first) applied to the areas'
 * covariates ax (na by p, column-major) gives each area's rate
 * 1 / (1 + exp(-eta)), times its population pop (at least 0) its
 * exposure. */
int clab_synexp(const double *coef, const double *ax, int na, int p, const double *pop,
                double *rate, double *expo);

#ifdef __cplusplus
}
#endif

#endif
