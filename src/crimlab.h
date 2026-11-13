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

/* The words for an error code. */
const char *clab_strerror(int code);

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

/* Centrographic summary (Bachi 1957; Neft 1962). x, y: n >= 1 points; w:
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
 * passes at most in all, else CLAB_EDOM. fold: NULL, or n fold ids in [0,
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
               int maxit, const int *fold, int nf, double *a0, double *beta,
               double *cv, int *best, int *nfit, double *work, int *iwork);

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
 * rows, or with one outcome value, are leaves. Out: pred (nt), the mean
 * over trees; oob (NULL or n), the mean over the trees whose bag left the
 * row out (NaN if none). work: p 255 + 3n + 512 + nt + 3 (2n + 1)
 * doubles; iwork: ceil((n + nt) p / 4) + 4n + nt + p + 1 ints. */
int clab_rfor(const double *x, int n, int p, const double *y,
              const double *xt, int nt, int ntree, int mtry, int minnode,
              double frac, int repl, unsigned int seed, double *pred,
              double *oob, double *work, int *iwork);

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
 * regression, 10 for a probability, as ranger). */
typedef struct {
  int kind;
  double alpha;
  int nlam, nfold;
  unsigned int seed;
  int ntree, mtry, minnode;
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
 * clustering (PLR and PLIV only; folds constant within a cluster). ly
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
             const int *clu, int nc, const clab_lrn *ly, const clab_lrn *ld,
             double trim, int fit, double *pred, double *est, double *psi,
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

#ifdef __cplusplus
}
#endif

#endif
