# Reference values from pedprobr's likelihood() of pedigrees built with
# pedtools, mutation models from pedmut: for tests/testthat/test-pedlr.R,
# written to tests/testthat/ref/ped.txt. Each case holds the pedigree
# (id, fid, mid; founders 0), the genotypes, the frequencies, the
# mutation model, rate and range, and the likelihood pedprobr gives. From
# the package root, with the packages in the library:
#   Rscript dev/ref/ped.R
# the packages it needs (dev/ref/all.R skips it without them)
need <- c("pedtools", "pedprobr", "pedmut")
source(file.path("dev", "ref", "put.R"))
suppressPackageStartupMessages({
  library(pedtools)
  library(pedprobr)
})

one <- function(x, typed, a1, a2, als, frq, mod, rate, range = 0.5) {
  mm <- if (mod == "none") {
    NULL
  } else {
    pedmut::mutationModel(
      if (mod == "step") "stepwise" else "equal",
      alleles = als,
      afreq = frq,
      rate = rate,
      rate2 = 0,
      range = range
    )
  }
  m <- marker(x, alleles = als, afreq = frq, mutmod = mm)
  for (i in seq_along(typed)) {
    genotype(m, typed[i]) <- c(a1[i], a2[i])
  }
  d <- as.data.frame(x)[c("id", "fid", "mid")]
  d[] <- lapply(d, as.character)
  list(
    ped = d,
    geno = data.frame(id = typed, a1 = a1, a2 = a2),
    freq = data.frame(allele = als, freq = frq),
    mut = mod,
    rate = rate,
    range = range,
    lik = likelihood(x, m)
  )
}

als <- as.character(9:13)
frq <- c(0.1, 0.2, 0.3, 0.25, 0.15)
cases <- list(
  # full siblings, stepwise mutation
  sib = one(
    nuclearPed(2),
    c("3", "4"),
    c("9", "11"),
    c("12", "13"),
    als,
    frq,
    "step",
    0.005
  ),
  # grandparent and grandchild through an untyped parent, equal mutation;
  # the grandchild's allele 9 needs a mutation
  gp = one(
    linearPed(2),
    c("1", "5"),
    c("10", "12"),
    c("11", "9"),
    als,
    frq,
    "equal",
    0.01
  ),
  # a child of first cousins (a loop), with one cousin typed
  loop = one(
    cousinPed(1, child = TRUE),
    c("7", "9"),
    c("10", "11"),
    c("11", "13"),
    als,
    frq,
    "step",
    0.002,
    0.3
  ),
  # a paternity trio inconsistent at the marker, stepwise mutation
  trio = one(
    nuclearPed(1),
    c("1", "2", "3"),
    c("10", "9", "12"),
    c("11", "13", "13"),
    als,
    frq,
    "step",
    0.003,
    0.5
  ),
  # half siblings without mutation
  half = one(
    halfSibPed(),
    c("4", "5"),
    c("10", "10"),
    c("12", "13"),
    als,
    frq,
    "none",
    0
  )
)
put(cases, "ped", need)
