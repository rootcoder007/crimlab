/* SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * libcsls: criminology and sociolegal statistics. The C interface.
 *
 * Every function is plain C from the outside: no exceptions cross it, no
 * memory is allocated for the caller, and nothing here knows about R. A
 * function returns CSLS_OK or an error code and writes its results into
 * arrays the caller owns. Indices are zero-based.
 */
#ifndef CSLS_H
#define CSLS_H

#ifdef __cplusplus
extern "C" {
#endif

#define CSLS_OK 0
#define CSLS_EARG 1 /* a size or an index is out of range */
#define CSLS_EDOM 2 /* a value is outside the domain of the method */

/* The words for an error code. */
const char *csls_strerror(int code);

/* Veil of darkness (Lean: Research.P2.Veil.odds_ratio_identifies_bias).
 * day, night: stops by group in daylight and in darkness, n groups; day
 * non-negative with day[ref] positive, night positive. shift is NULL or n
 * positive exposure factors. Out, n each: share_day, share_night, ratio (the
 * visible-race factor against group ref), sign (-1, 0, 1: the sign of
 * share_day - share_night) and, when shift is given, adj. */
int csls_vod(const double *day, const double *night, int n, int ref,
             const double *shift, double *share_day, double *share_night,
             double *ratio, int *sign, double *adj);

/* The risks behind an odds ratio (Lean: Research.P2.or_eq_rr_mul).
 * orr > 0; q in (0, 1) is the overall rate and s in (0, 1) the exposed share.
 * Out: a (risk among the exposed) and b (risk among the unexposed), the one
 * pair in (0, 1) with odds ratio orr and s a + (1 - s) b = q. */
int csls_orrisk(double orr, double q, double s, double *a, double *b);

/* Outcome test for one group
 * (Lean: Research.P2.OutcomeTest.infra_marginality).
 * q: guilt probabilities in [0, 1]; w: non-negative weights with a positive
 * sum; t: the search threshold. Out: stop_rate (weight share with q > t) and
 * hit_rate (weighted mean of q among them; NaN when nobody is searched). */
int csls_outtst(const double *q, const double *w, int n, double t,
                double *stop_rate, double *hit_rate);

/* Benchmark decomposition (Lean: Research.P2.benchmark_product).
 * pop, contact, force: positive counts by group, n groups. Out, n each, all
 * against group ref: cd (contact per resident), fd (force per contact) and
 * rd (force per resident); rd = cd * fd. */
int csls_dspbm(const double *pop, const double *contact, const double *force,
               int n, int ref, double *cd, double *fd, double *rd);

/* Dyadic offending rates (Lean: Research.P2.pair_exposure_rate_constant).
 * count[k] offences by group off[k] against group vic[k], m dyads; pop: g
 * positive group sizes. Out, m each: rate_off (per member of the offending
 * group), rate_pair (per offender-victim pair exposed), null_off (the rate
 * per offender under random mixing) and ratio (rate_pair over the pooled
 * pair rate). */
int csls_dyadrt(const double *count, const int *off, const int *vic, int m,
                const double *pop, int g, double *rate_off,
                double *rate_pair, double *null_off, double *ratio);

/* Genotype match probability with the Balding-Nichols correction for
 * shared ancestry (NRC II recommendation 4.2). p, q: allele frequencies in
 * (0, 1] at n loci; hom[i] non-zero for a homozygote (q[i] is then not
 * read), otherwise p[i] + q[i] <= 1; theta in [0, 1). Out: locus, n match
 * probabilities, and log10mp, the base-10 logarithm of their product. */
int csls_gtmp(const double *p, const double *q, const int *hom, int n,
              double theta, double *locus, double *log10mp);

/* Log-likelihood-ratio cost of a set of likelihood ratios (Brummer and du
 * Preez 2006). same: ns positive ratios from same-source comparisons; diff:
 * nd positive ratios from different-source comparisons. Out: cllr, 0 for a
 * perfect system and 1 for one that always reports a ratio of 1. */
int csls_cllr(const double *same, int ns, const double *diff, int nd,
              double *cllr);

/* Calibration of likelihood ratios by pooling adjacent violators. lr: n
 * positive ratios in ascending order; same[i] non-zero for a same-source
 * comparison (both kinds must occur). work: 2n doubles of scratch. Out: post,
 * n calibrated probabilities of same source at the sample's own prior, and
 * cllrmin, the cost of the calibrated ratios: what the system would cost if
 * its ratios kept their order and lost their miscalibration. */
int csls_lrcal(const double *lr, const int *same, int n, double *post,
               double *work, double *cllrmin);

/* Mixture: the chance that x unknown contributors carry only alleles of the
 * mixture and between them every allele flagged in need (Weir et al. 1997,
 * no coancestry). p: n allele frequencies of the mixture, sum at most 1;
 * at most 20 alleles may be flagged. Out: the probability. */
int csls_mixp(const double *p, const int *need, int n, int x, double *out);

/* Kinship index by locus: the chance of two genotypes under a relationship
 * with identity-by-descent probabilities k[0..2] over their chance if
 * unrelated. g1, g2: 2n allele codes, two per locus; f: 2n frequencies of
 * the alleles of g1. Out: idx, n indices. */
int csls_kinidx(const int *g1, const int *g2, const double *f, int n,
                const double *k, double *idx);

/* Paternity index by locus for a mother, child and alleged father, without
 * mutation. c, m, af: 2n allele codes each; f: 2n frequencies of the
 * child's alleles. Out: idx, n indices; NaN where the mother cannot be the
 * mother. */
int csls_patidx(const int *c, const int *m, const int *af, const double *f,
                int n, double *idx);

#ifdef __cplusplus
}
#endif

#endif
