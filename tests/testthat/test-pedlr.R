# pedlr(): against kinidx() and patidx() where they apply, against
# enumeration of every allele each person may carry, against pedprobr
# (ref/ped.txt), and its refusals.

frq <- data.frame(
  locus = rep(c("A", "B"), c(4, 3)),
  allele = c(10, 11, 12, 13, 7, 8, 9),
  freq = c(0.1, 0.2, 0.3, 0.4, 0.25, 0.5, 0.25)
)
gt <- function(id, a1, a2) {
  data.frame(
    id = rep(id, each = 2),
    locus = rep(c("A", "B"), length(id)),
    a1 = a1,
    a2 = a2
  )
}
mat <- function(g, who) {
  g <- g[g$id == who, ]
  cbind(g$a1, g$a2)
}
fq <- function(g) {
  # the frequencies of g's alleles, two per locus, for kinidx()
  p <- stats::setNames(frq$freq, paste(frq$locus, frq$allele))
  cbind(p[paste(g$locus, g$a1)], p[paste(g$locus, g$a2)])
}

test_that("pedlr() equals kinidx() for parents, siblings, halves, cousins", {
  g <- gt(c("x", "y"), c(10, 7, 11, 7), c(12, 8, 12, 9))
  f1 <- fq(g[g$id == "x", ])
  peds <- list(
    parent = data.frame(
      id = c("x", "y", "m"),
      fid = c(NA, "x", NA),
      mid = c(NA, "m", NA)
    ),
    sibling = data.frame(
      id = c("f", "m", "x", "y"),
      fid = c(NA, NA, "f", "f"),
      mid = c(NA, NA, "m", "m")
    ),
    half = data.frame(
      id = c("f", "m1", "m2", "x", "y"),
      fid = c(NA, NA, NA, "f", "f"),
      mid = c(NA, NA, NA, "m1", "m2")
    ),
    cousin = data.frame(
      id = c("gf", "gm", "a", "b", "sa", "sb", "x", "y"),
      fid = c(NA, NA, "gf", "gf", NA, NA, "a", "sb"),
      mid = c(NA, NA, "gm", "gm", NA, NA, "sa", "b")
    )
  )
  for (rel in names(peds)) {
    r <- pedlr(peds[[rel]], g, frq)
    k <- kinidx(mat(g, "x"), mat(g, "y"), f1, rel = rel)
    expect_equal(r$locus$lr, k$locus$idx, tolerance = 1e-13, info = rel)
    expect_equal(r$log10lr, k$log10idx, tolerance = 1e-13, info = rel)
  }
})

test_that("pedlr() equals patidx() for a trio", {
  g <- gt(c("c", "m", "af"), c(11, 7, 10, 7, 12, 8), c(12, 9, 11, 8, 13, 9))
  trio <- data.frame(
    id = c("af", "m", "c"),
    fid = c(NA, NA, "af"),
    mid = c(NA, NA, "m")
  )
  alt <- data.frame(
    id = c("af", "m", "c", "u"),
    fid = c(NA, NA, "u", NA),
    mid = c(NA, NA, "m", NA)
  )
  r <- pedlr(trio, g, frq, alt)
  k <- patidx(
    mat(g, "c"),
    mat(g, "af"),
    fq(g[g$id == "c", ]),
    mother = mat(g, "m")
  )
  expect_equal(r$locus$lr, k$locus$idx, tolerance = 1e-13)
  # the paternal allele is 12 at A (0.5 / 0.3) and 9 at B (0.5 / 0.25)
  expect_equal(r$locus$lr, c(0.5 / 0.3, 0.5 / 0.25), tolerance = 1e-13)
})

# every allele pair of every person, founders from p, children from a
# parent's allele (probability 1/2) passed through the mutation matrix m
brute <- function(ped, typed, p, m) {
  n <- nrow(ped)
  k <- length(p)
  fa <- match(ped$fid, ped$id)
  mo <- match(ped$mid, ped$id)
  grid <- as.matrix(expand.grid(rep(list(seq_len(k)), 2 * n)))
  tot <- 0
  for (r in seq_len(nrow(grid))) {
    a <- grid[r, ]
    pr <- 1
    for (i in seq_len(n)) {
      x <- a[2 * i - 1]
      y <- a[2 * i]
      g <- typed[[ped$id[i]]]
      if (
        !is.null(g) &&
          !setequal(c(x, y), g) ||
          (!is.null(g) && (x == y) != (g[1] == g[2]))
      ) {
        pr <- 0
        break
      }
      if (is.na(fa[i])) {
        pr <- pr * p[x] * p[y]
      } else {
        pr <- pr *
          0.5 *
          (m[a[2 * fa[i] - 1], x] + m[a[2 * fa[i]], x]) *
          0.5 *
          (m[a[2 * mo[i] - 1], y] + m[a[2 * mo[i]], y])
      }
    }
    tot <- tot + pr
  }
  tot
}

test_that("pedlr() equals the enumeration of every allele, with mutation", {
  p <- c(0.2, 0.3, 0.5)
  f3 <- data.frame(locus = "L", allele = c(14, 15, 16), freq = p)
  ped <- data.frame(
    id = c("g", "h", "p", "s", "c"),
    fid = c(NA, NA, "g", NA, "p"),
    mid = c(NA, NA, "h", NA, "s")
  )
  geno <- data.frame(
    id = c("g", "c"),
    locus = "L",
    a1 = c(14, 16),
    a2 = c(15, 16)
  )
  w <- 0.4^abs(outer(14:16, 14:16, "-"))
  diag(w) <- 0
  for (mod in c("equal", "step")) {
    m <- if (mod == "equal") matrix(0.01 / 2, 3, 3) else 0.01 * w / rowSums(w)
    diag(m) <- 0.99
    r <- pedlr(ped, geno, f3, mut = mod, rate = 0.01, range = 0.4)
    want <- brute(ped, list(g = c(1, 2), c = c(3, 3)), p, m)
    expect_equal(10^r$locus$l1, want, tolerance = 1e-12, info = mod)
    expect_equal(
      10^r$locus$l2,
      2 * 0.2 * 0.3 * 0.25,
      tolerance = 1e-13,
      info = mod
    )
  }
  # a father 14/15 cannot have a child 16/16 without mutation
  pc <- data.frame(
    id = c("g", "s", "c"),
    fid = c(NA, NA, "g"),
    mid = c(NA, NA, "s")
  )
  expect_identical(pedlr(pc, geno, f3)$lr, 0)
  expect_gt(pedlr(pc, geno, f3, mut = "equal")$lr, 0)
})

test_that("pedlr() agrees with pedprobr's likelihoods", {
  # pedprobr's likelihood() of each pedigree (pedtools, pedmut mutation
  # models), stored in ref/ped.txt by dev/ref/ped.R (versions named there)
  ref <- refget("ped")
  for (k in setdiff(names(ref), "pkgs")) {
    x <- ref[[k]]
    geno <- data.frame(
      id = x$geno$id,
      locus = "L",
      a1 = x$geno$a1,
      a2 = x$geno$a2
    )
    freq <- data.frame(locus = "L", allele = x$freq$allele, freq = x$freq$freq)
    r <- pedlr(x$ped, geno, freq, mut = x$mut, rate = x$rate, range = x$range)
    expect_equal(10^r$locus$l1, x$lik, tolerance = 1e-12, info = k)
  }
})

test_that("pedlr() takes sex-specific rates and lumps the alleles not listed", {
  ped <- data.frame(
    id = c("f", "m", "c"),
    fid = c(NA, NA, "f"),
    mid = c(NA, NA, "m")
  )
  geno <- data.frame(
    id = c("f", "m", "c"),
    locus = "L",
    a1 = c(1, 3, 2),
    a2 = c(1, 3, 3)
  )
  fr <- data.frame(locus = "L", allele = 1:3, freq = c(0.2, 0.3, 0.4))
  # the father's 1 mutates (paternal rate 0.02, split over 2, 3 and the
  # rest) to 2 while the mother passes 3, or to 3 while her 3 mutates to 2
  # (maternal rate 0.004)
  r <- pedlr(ped, geno, fr, mut = "equal", rate = c(0.004, 0.02))
  want <- (0.02 / 3) * ((1 - 0.004) + 0.004 / 3) / (2 * 0.3 * 0.4)
  expect_equal(r$locus$lr, want, tolerance = 1e-13)
  # without mutation the rest changes nothing
  expect_identical(pedlr(ped, geno[c(1, 3), ], fr)$lr, 0)
})

test_that("pedlr() prints and refuses bad input", {
  ped <- data.frame(id = c("f", "c"), fid = c(NA, "f"), mid = c(NA, "m"))
  ok <- data.frame(
    id = c("f", "m", "c"),
    fid = c(NA, NA, "f"),
    mid = c(NA, NA, "m")
  )
  geno <- data.frame(id = c("f", "c"), locus = "L", a1 = c(1, 1), a2 = c(2, 3))
  fr <- data.frame(locus = "L", allele = 1:3, freq = c(0.2, 0.3, 0.5))
  expect_output(print(pedlr(ok, geno, fr)), "pedlr: 1 markers")
  cyc <- data.frame(
    id = c("a", "b", "c"),
    fid = c(NA, "c", "b"),
    mid = c(NA, "a", "a")
  )
  errs(
    "pedlr",
    list(
      list(
        "data frame with columns id, fid and mid",
        ped = list(1),
        geno = geno,
        freq = fr
      ),
      list("columns id, locus, a1 and a2", ped = ok, geno = 1, freq = fr),
      list(
        "columns locus, allele and freq",
        ped = ok,
        geno = geno,
        freq = fr[1:2]
      ),
      list(
        "must be numbers",
        ped = ok,
        geno = geno,
        freq = transform(fr, freq = "a")
      ),
      list(
        "one of none, equal, step",
        ped = ok,
        geno = geno,
        freq = fr,
        mut = "x"
      ),
      list(
        "one or two numbers from 0 to 1",
        ped = ok,
        geno = geno,
        freq = fr,
        mut = "equal",
        rate = 2
      ),
      list(
        "strictly between 0 and 1",
        ped = ok,
        geno = geno,
        freq = fr,
        range = 1
      ),
      list("mother m of c is not in `ped`", ped = ped, geno = geno, freq = fr),
      list(
        "one parent",
        ped = transform(ok, mid = c(NA, NA, NA)),
        geno = geno,
        freq = fr
      ),
      list(
        "cycle",
        ped = cyc,
        geno = data.frame(id = "a", locus = "L", a1 = 1, a2 = 1),
        freq = fr
      ),
      list(
        "both a father and a mother",
        ped = data.frame(
          id = c("a", "b", "c", "d"),
          fid = c(NA, NA, "a", "c"),
          mid = c(NA, NA, "b", "a")
        ),
        geno = geno[0, ],
        freq = fr
      ),
      list("twice", ped = rbind(ok, ok[1, ]), geno = geno, freq = fr),
      list(
        "not in `ped`",
        ped = ok,
        geno = transform(geno, id = c("f", "z")),
        freq = fr
      ),
      list(
        "does not list",
        ped = ok,
        geno = transform(geno, a2 = c(2, 9)),
        freq = fr
      ),
      list(
        "sum to more than 1",
        ped = ok,
        geno = geno,
        freq = transform(fr, freq = c(0.5, 0.5, 0.5))
      ),
      list(
        "needs frequencies that sum to 1",
        ped = ok,
        geno = geno,
        freq = fr[1:2, ],
        mut = "step"
      ),
      list(
        "whole-number alleles",
        ped = ok,
        geno = geno,
        freq = transform(fr, allele = c("1", "2", "x")),
        mut = "step"
      ),
      list(
        "one allele of two",
        ped = ok,
        geno = transform(geno, a2 = c(NA, 3)),
        freq = fr
      ),
      list("must type someone", ped = ok, geno = geno[0, ], freq = fr)
    )
  )
})
