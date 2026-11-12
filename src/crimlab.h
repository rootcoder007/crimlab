/* SPDX-License-Identifier: AGPL-3.0-or-later
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

#define CLAB_OK 0
#define CLAB_EARG 1 /* a size or an index is out of range */
#define CLAB_EDOM 2 /* a value is outside the domain of the method */

/* The words for an error code. */
const char *clab_strerror(int code);

/* Veil of darkness (Lean: Crimlab.Veil.odds_ratio_identifies_bias).
 * day, night: stops by group in daylight and in darkness, n groups; day
 * non-negative with day[ref] positive, night positive. shift is NULL or n
 * positive exposure factors. Out, n each: share_day, share_night, ratio (the
 * visible-race factor against group ref), sign (-1, 0, 1: the sign of
 * share_day - share_night) and, when shift is given, adj. */
int clab_vod(const double *day, const double *night, int n, int ref,
             const double *shift, double *share_day, double *share_night,
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
 * positive group sizes. Out, m each: rate_off (per member of the offending
 * group), rate_pair (per offender-victim pair exposed), null_off (the rate
 * per offender under random mixing) and ratio (rate_pair over the pooled
 * pair rate). */
int clab_dyadrt(const double *count, const int *off, const int *vic, int m,
                const double *pop, int g, double *rate_off,
                double *rate_pair, double *null_off, double *ratio);

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

#ifdef __cplusplus
}
#endif

#endif
