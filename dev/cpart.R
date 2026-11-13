# The search for other implementations of each exported function, written
# to inst/bench/counterparts.csv. The table below is the one hand-kept part:
# for each candidate found, whether it computes the same quantity and, in
# one line, the argument and output mapping (same = "yes") or the difference
# (same = "no"). Everything else is derived from it.
#   Rscript dev/cpart.R            (from the package root)
#
# Where the candidates came from (searched on the date in the table):
# tools::CRAN_package_db() titles and descriptions, searched for the method
# names and their authors (veil of darkness, outcome test, hit rate,
# infra-marginality, benchmark, interracial and dyadic rates, odds ratio to
# risk ratio, match probability, Balding, theta correction, kinship and
# paternity index, DNA mixture, LRmix, Cllr, log-likelihood-ratio cost,
# empirical cross-entropy, Tippett, posterior from a likelihood ratio,
# database search, Kitagawa, Das Gupta, rate decomposition, exposure and
# sensitivity bounds, misclassification, regression to the mean,
# incapacitation, selective labels, length of stay and censoring bounds,
# Manski and treatment-effect bounds, dynamic selection and frailty,
# brokerage, Cheeger and conductance, Louvain, co-offending and bipartite
# projection, k-core); the function index of every CRAN package that
# pkgcheck uses (the same terms as function names); the CRAN Task Views
# OfficialStatistics, Spatial, CausalInference, Epidemiology and
# NetworkAnalysis (their package lists); Bioconductor's software list for forensic
# genetics (no package there computes any of these quantities); and the
# packages named in the project dossier and its reference generators
# (dev/doss). Packages not on CRAN (DNAprofiles, euroformix, FamLink,
# hmeasure, the ROC package on GitHub) are listed when they were candidates.
#
# own = TRUE marks rmorie, the maintainer's own package: it is reported
# apart from the third-party packages and never counts as "another package".
m <- read.csv(text = '
fn,pkg,pkgfn,same,own,note
vodtst,rmorie,morie_veil_of_darkness,yes,TRUE,"day and night stops per group; shday, shnight, ratio = share_day, share_night, bias_ratio"
vodtst,searchlight,sl_veil_of_darkness,no,FALSE,"a logistic regression of stop records on darkness with clock-time controls (an odds ratio), not the share contrast from counts"
orrisk,rmorie,morie_relative_risk_from_or,yes,TRUE,"or, rate, share; exposed, unexposed = risks (found by root search, agree to 1e-6)"
orrisk,effectsize,oddsratio_to_riskratio,no,FALSE,"needs the risk of the unexposed (p0); orrisk solves for it from the overall rate and the exposed share"
orrisk,MKclass,or2rr,no,FALSE,"Zhang and Yu with p0 given, as effectsize; not the risks from rate and share"
orrisk,MKmisc,or2rr,no,FALSE,"as MKclass::or2rr"
orrisk,metainc,or2rr,no,FALSE,"as MKclass::or2rr"
orrisk,orsk,orsk,no,FALSE,"rebuilds a 2 by 2 table from group sizes and a published odds ratio with its interval"
outtst,rmorie,morie_outcome_test,yes,TRUE,"q, w, t; stop, hit = stop_rate, hit_rate"
outtst,searchlight,sl_hit_rates,no,FALSE,"hit rates per group from search records with Wilson intervals; no risk threshold"
outtst,fairness,pred_rate_parity,no,FALSE,"parity ratios between groups of the predictions of a classifier, not the outcome test"
dspbm,rmorie,morie_disparity_benchmark,yes,TRUE,"pop, contact, force; cd, fd, rd = the three disparities"
dspbm,searchlight,sl_benchmark,no,FALSE,"compares submitted record counts with published totals (coverage), not a disparity decomposition"
dyadrt,rmorie,morie_interracial_rates,yes,TRUE,"count matrix (rmorie: a named vector a_on_b), pop; the four rates"
gtmp,Familias,FamiliasPosterior,yes,FALSE,"kinship = theta; log10 of the ratio of the likelihoods of two unrelated people and of one, per locus"
gtmp,fbnet,getQSetRMP,no,FALSE,"genotype probabilities in a pedigree Bayesian network, not the Balding-Nichols match probability"
gtmp,DNAtools,Pnm_all,no,FALSE,"the distribution of the number of matching loci between database profiles"
gtmp,forensim,LR,no,FALSE,"a mixture ratio; a single-source profile with theta is not its documented use"
kinidx,Familias,FamiliasPosterior,yes,FALSE,"full-sibling pedigree against unrelated; log10 ratio summed over loci"
kinidx,forrel,kinshipLR,yes,FALSE,"nuclear pedigree against two singletons; lnLRtotal / log(10)"
kinidx,paramlink,LR,yes,FALSE,"as forrel (its predecessor, kept for legacy use): likelihoodsPerSystem ratio"
kinidx,relSim,lrSib,yes,FALSE,"profiles as 2 x loci allele codes, Freqs = list(freqs = per-locus frequencies); the product over loci (log10)"
kinidx,relSim,lrPC,yes,FALSE,"rel = parent; as lrSib"
kinidx,pedprobr,likelihood,no,FALSE,"the likelihood of one pedigree; the index is a ratio of two (forrel::kinshipLR forms it)"
kinidx,ribd,kappaIBD,no,FALSE,"relatedness coefficients of a pedigree, no genotypes"
kinidx,fbnet,reportLR,no,FALSE,"missing-person likelihood ratios in a pedigree Bayesian network"
kinidx,DNAprofiles,ki,no,FALSE,"not on CRAN (archived); not installed"
kinidx,FamLink,,no,FALSE,"not on CRAN; not installed"
patidx,Familias,FamiliasPosterior,yes,FALSE,"trio against mother, child and an untested father"
patidx,forrel,kinshipLR,yes,FALSE,"trio against the same with the alleged father a singleton"
patidx,paramlink,LR,yes,FALSE,"as forrel"
patidx,relSim,lrPC,yes,FALSE,"duo (mother = NULL): the parent-child ratio"
mixlr,forensim,LR,yes,FALSE,"Repliste = mixture alleles, Tp, Td, xp, xd; num, deno, LR. Not installable on l14 (Tk; see bench.md)"
mixlr,relSim,lrMix,yes,FALSE,"victim plus suspect against victim plus one unknown, theta 0, per locus"
mixlr,relMix,relMix,no,FALSE,"relationship inference from mixtures with drop-in and drop-out (a GUI-driven model)"
mixlr,KinMixLite,KinMix,no,FALSE,"peak-height model on Bayesian networks"
mixlr,DNAmixturesLite,,no,FALSE,"peak-height model"
mixlr,numberofalleles,pr_total_number_of_distinct_alleles,no,FALSE,"the distribution of the number of distinct alleles"
mixlr,simDNAmixtures,,no,FALSE,"simulates mixtures"
mixlr,euroformix,,no,FALSE,"not on CRAN; a peak-height model"
cllr,comparison,calc.ece,yes,FALSE,"prior 0.5: ece"
cllr,idiolect,Cllr,yes,FALSE,"unexported (idiolect:::Cllr, bundled from the ROC package): data frame of score = log(ratio) and target"
cllr,probcal,ece,no,FALSE,"binned expected calibration error of probabilities"
cllr,CalibratR,getECE,no,FALSE,"as probcal::ece"
cllr,ROC,,no,FALSE,"on GitHub only (davidavdav/ROC), not CRAN; its code is the one idiolect bundles"
cllrmin,comparison,calc.ece,yes,FALSE,"ece.cal at prior 0.5"
cllrmin,idiolect,Cllr,yes,FALSE,"unexported; opt = TRUE (pool-adjacent-violators without the Laplace points)"
lrcal,comparison,calibrate.set,yes,FALSE,"LR.cal.ss, LR.cal.ds, each sorted"
lrcal,idiolect,opt.llr,yes,FALSE,"unexported; laplace = FALSE; exp(opt.llr)"
lrcal,isotone,gpava,no,FALSE,"the pool-adjacent-violators fit itself; the ratios are its posteriors over prior odds"
lrcal,fdrtool,monoreg,no,FALSE,"as isotone::gpava"
eceplt,comparison,calc.ece,yes,FALSE,"prior grid; ece, ece.cal, ece.null; then plot()"
tipplt,ROC,,no,FALSE,"on GitHub only, not CRAN"
lrpost,idiolect,posterior,no,FALSE,"a table of posteriors for a fixed set of priors from one log10 ratio"
dbsrch,DNAtools,dbExpect,no,FALSE,"expected numbers of matching pairs within a database"
clrsmp,DasGuptR,dgnpop,yes,FALSE,"two factors (rate, share) by part, ratefunction sum(rate * share): the rate and composition effects in total"
clrsmp,DemoDecomp,horiuchi,yes,FALSE,"f = sum(rate * share); exact for this bilinear rate: per-part rate effects and both totals"
clrsmp,DemoDecomp,stepwise_replacement,no,FALSE,"an order-dependent split, not the symmetric one"
clrsmp,ddecompose,,no,FALSE,"distributional (Oaxaca-Blinder type) decomposition of regressions"
clrsmp,oaxaca,oaxaca,no,FALSE,"Blinder-Oaxaca decomposition of a regression"
expbnd,rmorie,morie_disparity_exposure_bounds,yes,TRUE,"y, m, gamma; rate, lo, hi, rlo, rhi"
expbnd,rbounds,,no,FALSE,"Rosenbaum bounds for matched designs"
expbnd,sensemakr,,no,FALSE,"sensitivity of regression estimates to omitted variables"
expbnd,tipr,,no,FALSE,"tipping points for unmeasured confounding"
misrd,rmorie,morie_misclassified_rd,yes,TRUE,"true, measured, factor, d1, d0"
misrd,episensr,misclass_cov,no,FALSE,"corrects an observed table for misclassification of a covariate (the inverse problem)"
misrd,EValue,misclassification,no,FALSE,"bounds and E-values for a misclassified estimate"
rtmhs,rmorie,morie_regression_to_mean,yes,TRUE,"nsel, change, mirror, sym, excess, low, lowsym"
rtmhs,regtomean,meechua_reg,no,FALSE,"the Mee-Chua test of a pre-post change in one group"
incap,rmorie,morie_incapacitation,yes,TRUE,"lam, q, s, share; rate, prev, marg"
incapr,rmorie,morie_incapacitation_career,yes,TRUE,"prev, net, upper, lower, later"
sellab,rmorie,morie_selective_labels,yes,TRUE,"obs, unobs, lower, upper, naive"
losbnd,rmorie,morie_backlog_censoring,yes,TRUE,"mean, age, lower, bias"
sntbnd,rmorie,morie_sentence_effect_bounds,yes,TRUE,"outcome and effect bounds, naive difference"
sntbnd,ATbounds,atebounds,no,FALSE,"bounds with covariates by limited information pooling"
sntbnd,bpbounds,,no,FALSE,"Balke-Pearl bounds: needs an instrument"
sntbnd,causaloptim,,no,FALSE,"bounds by linear programming on a causal graph; not installed"
hzdsel,rmorie,morie_hazard_selection,yes,TRUE,"share0, share1, haz0, haz1, ratio"
hzdsel,frailtypack,,no,FALSE,"fits frailty models to survival data"
ntdeg,igraph,degree + strength,yes,FALSE,"simplified graph, weights summed"
ntdeg,sna,degree,yes,FALSE,"gmode graph on the 0/1 adjacency; degree only"
ntdeg,tnet,degree_w,yes,FALSE,"symmetric weighted edge list; degree, output"
ntdeg,netmem,gen_degree,yes,FALSE,"digraph = FALSE, loops = FALSE; weighted with alpha = 1 for the strength"
ntdeg,NetworkToolbox,degree,yes,FALSE,"0/1 adjacency"
ntdeg,NetworkToolbox,strength,no,FALSE,"rounds the strength to two decimals"
ntbtw,igraph,betweenness,yes,FALSE,"weights = NA"
ntbtw,sna,betweenness,yes,FALSE,"gmode graph"
ntbtw,tnet,betweenness_w,yes,FALSE,"alpha = 0 on the symmetric edge list"
ntbtw,netmem,betweenness_centrality,yes,FALSE,"digraph = FALSE"
ntbtw,NetworkToolbox,betweenness,no,FALSE,"counts each pair of an undirected graph twice"
ntcls,igraph,closeness + harmonic_centrality,yes,FALSE,"weights = NA; cls and hrm"
ntcls,tnet,closeness_w,yes,FALSE,"alpha = 0, gconly = FALSE: closeness is the harmonic centrality"
ntcls,netmem,closeness_centrality,yes,FALSE,"digraph = FALSE (cls, reached nodes only) and harmonic = TRUE (hrm)"
ntcls,sna,closeness,no,FALSE,"cmode suminvundir is the harmonic centrality divided by n - 1; the default refuses unreachable pairs"
ntcls,centiserve,closeness.freeman,no,FALSE,"refuses a graph that is not connected"
ntcls,NetworkToolbox,closeness,no,FALSE,"another definition: values differ"
ntevc,igraph,eigen_centrality,yes,FALSE,"weighted, scaled to a maximum of 1"
ntevc,netmem,eigenvector_centrality,yes,FALSE,"weighted, scale max; the same on a connected graph (timed on the largest component)"
ntevc,sna,evcent,no,FALSE,"scaled to unit length, not to a maximum of 1"
ntevc,NetworkToolbox,eigenvector,no,FALSE,"differs by 1e-3 on a connected graph"
ntcore,igraph,coreness,yes,FALSE,""
ntcore,sna,kcores,yes,FALSE,"mode graph"
ntcore,netmem,k_core,yes,FALSE,""
ntcmp,igraph,components,yes,FALSE,"weak and strong; the same partition"
ntcmp,sna,component.dist,yes,FALSE,"connected weak and strong; the same partition"
ntcmp,netmem,components_id,yes,FALSE,"weak and strong; the same partition"
ntlvn,igraph,cluster_louvain,yes,FALSE,"a heuristic: partitions differ with node order; checked by the modularity reached"
ntlvn,modMax,louvain,no,FALSE,"unweighted; not installed (last release 2015)"
ntlvn,NetworkToolbox,louvain,no,FALSE,"repeats igraph::cluster_louvain and keeps the most frequent partition"
ntlvn,leiden,leiden,no,FALSE,"the Leiden method, not Louvain"
ntmod,igraph,modularity,yes,FALSE,"weighted"
ntmod,netmem,modularity_score,yes,FALSE,"weighted = TRUE"
ntbrk,sna,brokerage,yes,FALSE,"raw.nli: w_I, w_O, b_IO, b_OI, b_O = coord, itin, repr, gate, liais"
ntbrk,netmem,brokerage_roles,yes,FALSE,"coordinator, consultant, representative, gatekeeper, liaison"
ntchg,clustAnalytics,conductance,no,FALSE,"the conductance of given clusters, not the spectral bounds"
coofg,igraph,bipartite_projection,yes,FALSE,"multiplicity = TRUE (with building the graph)"
coofg,tnet,projecting_tm,yes,FALSE,"method sum on the offender-incident list"
coofg,netmem,matrix_projection,yes,FALSE,"matrix2 of a dense offender by incident matrix"
', stringsAsFactors = FALSE, strip.white = TRUE)
# network interfaces whose measures call igraph's (or sna's) routine, read
# in their sources: listed, not timed apart from the routine they call
wr <- c(ntdeg = "igraph::degree", ntbtw = "igraph::betweenness",
        ntcls = "igraph::closeness", ntevc = "igraph::eigen_centrality",
        ntcore = "igraph::coreness", ntcmp = "igraph::components",
        ntlvn = "igraph::cluster_louvain", ntmod = "igraph::modularity",
        ntbrk = "sna::brokerage")
for (p in c("netrics", "cograph")) {
  m <- rbind(m, data.frame(fn = names(wr), pkg = p, pkgfn = "", same = "no",
                           own = FALSE, note = paste("calls", wr)))
}
m <- rbind(m, data.frame(fn = "coofg", pkg = "manynet", pkgfn = "", same = "no",
                         own = FALSE, note = "calls igraph::bipartite_projection"))
m$date <- "2026-11-12"
# the functions this table covers: those of the disparity, forensic,
# policing, corrections and network families (the spatial functions and
# the learners have speed tables of their own)
fns <- sort(c(
  "cllr", "cllrmin", "clrsmp", "coofg", "dbsrch", "dspbm", "dyadrt",
  "eceplt", "expbnd", "gtmp", "hzdsel", "incap", "incapr", "kinidx",
  "losbnd", "lrcal", "lrpost", "misrd", "mixlr", "ntbrk", "ntbtw",
  "ntchg", "ntcls", "ntcmp", "ntcore", "ntdeg", "ntevc", "ntlvn",
  "ntmod", "orrisk", "outtst", "patidx", "proof", "rtmhs", "sellab",
  "sntbnd", "tipplt", "vodtst"
))
stopifnot(length(fns) == 38L, all(m$fn %in% fns))
miss <- setdiff(fns, m$fn)
# the functions with no candidate at all
m <- rbind(m, data.frame(fn = miss, pkg = "", pkgfn = "", same = "no",
                         own = FALSE, note = "no package found", date = m$date[1]))
m <- m[order(match(m$fn, fns), m$same != "yes", m$own), ]
utils::write.csv(m, "inst/bench/counterparts.csv", row.names = FALSE)
third <- unique(m$fn[m$same == "yes" & !m$own])
cat(sprintf("%d functions: %d with a third-party counterpart; none for %s\n",
            length(fns), length(third), paste(setdiff(fns, third), collapse = ", ")))
