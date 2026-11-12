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

#ifdef __cplusplus
}
#endif

#endif
