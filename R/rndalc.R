#' Randomisation lists for a trial
#'
#' Allocates participants to arms by simple, permuted-block or stratified
#' block randomisation from a seed, so that the list can be made again
#' exactly from the seed and the settings, in R or on the command line.
#'
#' With `block = 0` each participant is allocated alone (simple
#' randomisation), to arm `a` with chance `ratio[a] / sum(ratio)`.
#' Otherwise blocks are made one after another, each of a size drawn with
#' equal chances from `block` (every size a multiple of `sum(ratio)`),
#' holding the arms in the ratio in an order drawn uniformly (a
#' Fisher-Yates shuffle). Block sizes drawn at random keep the next
#' allocation hard to guess at a site (Schulz and Grimes 2002). The list
#' stops at the places asked for, so the last block of a stratum may be
#' incomplete. With `strata`, each stratum has its own sequence of blocks
#' (Rosenberger and Lachin 2016, ch. 3).
#'
#' The draws come from the package's own generator (xorshift64* seeded
#' through splitmix64; Vigna 2016), the same on every platform; the seed
#' is returned with the list. Sealing the list (a commitment signed and
#' timestamped before enrolment) is left to other tools.
#'
#' @param n The number of participants: a whole number from 1 to
#'   10000000. With `strata` it may be left out (it is their total).
#' @param arms The names of the arms: two or more, each once.
#' @param ratio The allocation ratio: one whole number from 1 per arm.
#' @param block The block sizes drawn from: whole numbers, each a multiple
#'   of `sum(ratio)`; 0 for simple randomisation.
#' @param strata `NULL`, or the strata and their places: a named vector of
#'   whole numbers from 0.
#' @param seed `NULL` (drawn from R's generator, so [set.seed()] fixes it)
#'   or a whole number from 0 to 2147483647.
#' @param prefix The start of the participant ids: one string.
#' @return An object of class `rndalc`, a list: `alloc`, a data frame in
#'   enrolment order with `id`, `stratum` (with `strata`), `block` (the
#'   block's number within its stratum; `NA` for simple randomisation)
#'   and `arm`; `tab`, the participants of each arm by stratum (a matrix,
#'   one row `all` without strata); and `seed`, the seed used.
#' @references
#'   Rosenberger, W. F. and Lachin, J. M. (2016). *Randomization in
#'   Clinical Trials: Theory and Practice*, 2nd ed. Wiley.
#'   \doi{10.1002/9781118742112}
#'
#'   Schulz, K. F. and Grimes, D. A. (2002). Generation of allocation
#'   sequences in randomised trials: chance, not choice. *The Lancet*,
#'   359(9305), 515-519. \doi{10.1016/S0140-6736(02)07683-3}
#'
#'   Vigna, S. (2016). An experimental exploration of Marsaglia's xorshift
#'   generators, scrambled. *ACM Transactions on Mathematical Software*,
#'   42(4), 30. \doi{10.1145/2845077}
#' @examples
#' a <- rndalc(40, seed = 12)
#' a
#' head(a$alloc)
#' s <- rndalc(arms = c("drug", "placebo"), ratio = c(2, 1),
#'             strata = c(Toronto = 30, Ottawa = 24), seed = 7)
#' s$tab
#' rndalc(20, block = 0, seed = 1)$tab
#' @export
rndalc <- function(
  n = NULL,
  arms = c("A", "B"),
  ratio = rep(1, length(arms)),
  block = c(2, 3) * sum(ratio),
  strata = NULL,
  seed = NULL,
  prefix = "P"
) {
  f <- "rndalc"
  if (!is.character(arms) || length(arms) < 2L || anyNA(arms) ||
        !all(nzchar(arms)) || anyDuplicated(arms)) {
    .ckstop(f, "`arms` must name two or more arms, each once")
  }
  ratio <- .cknum(ratio, "ratio", f, scalar = FALSE, min = 1, max = 1e6,
                  whole = TRUE)
  if (length(ratio) != length(arms)) {
    .ckstop(f, "`ratio` must give one number per arm")
  }
  sr <- sum(ratio)
  block <- .cknum(block, "block", f, scalar = FALSE, min = 0, max = 1e6,
                  whole = TRUE)
  simple <- length(block) == 1L && block == 0
  if (!simple && any(block == 0 | block %% sr != 0)) {
    .ckstop(f, paste("`block` must be 0, or sizes that are multiples of",
                     "sum(ratio) (%s)"), format(sr))
  }
  nm <- NULL
  if (is.null(strata)) {
    n <- .cknum(n, "n", f, min = 1, max = 1e7, whole = TRUE)
    cnt <- n
  } else {
    cnt <- .cknum(strata, "strata", f, scalar = FALSE, min = 0, whole = TRUE)
    nm <- names(strata)
    if (is.null(nm) || anyNA(nm) || !all(nzchar(nm)) || anyDuplicated(nm)) {
      .ckstop(f, "`strata` must be named, each stratum once")
    }
    tot <- sum(cnt)
    if (tot < 1 || tot > 1e7) {
      .ckstop(f, "`strata` must hold from 1 to 10000000 places in all")
    }
    if (!is.null(n)) {
      n <- .cknum(n, "n", f, min = 1, max = 1e7, whole = TRUE)
      if (n != tot) {
        .ckstop(f, "`n` must equal the places of `strata` (%s)", format(tot))
      }
    }
    n <- tot
  }
  seed <- .xpseed(seed, f)
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix)) {
    .ckstop(f, "`prefix` must be one string")
  }
  k <- .Call(C_clab_r_rndalc, as.integer(cnt), as.integer(ratio),
             if (simple) integer(0) else as.integer(block), seed)
  n <- as.integer(n)
  id <- paste0(prefix, formatC(seq_len(n), width = nchar(as.character(n)),
                               flag = "0"))
  alloc <- data.frame(id = id, stringsAsFactors = FALSE)
  if (!is.null(nm)) alloc$stratum <- rep(nm, cnt)
  b <- k$block
  b[b == 0L] <- NA_integer_
  alloc$block <- b
  alloc$arm <- arms[k$arm + 1L]
  g <- if (is.null(nm)) rep(1L, n) else rep(seq_along(nm), cnt)
  ka <- length(arms)
  tab <- matrix(tabulate((g - 1L) * ka + k$arm + 1L, length(cnt) * ka),
                ncol = ka, byrow = TRUE,
                dimnames = list(if (is.null(nm)) "all" else nm, arms))
  structure(list(alloc = alloc, tab = tab, seed = seed), class = "rndalc")
}

#' @rdname rndalc
#' @param x An object of class `rndalc`.
#' @param ... Not used by `print()`.
#' @export
print.rndalc <- function(x, ...) {
  cat("rndalc: ", nrow(x[["alloc"]]), " participants, seed ",
      format(x[["seed"]], scientific = FALSE), "\n", sep = "")
  print(x[["tab"]])
  invisible(x)
}

# The seed of a kernel's generator: one whole number from 0 to
# 2147483647, drawn from R's generator when NULL (so set.seed() fixes it)
# and returned with the result.
.xpseed <- function(seed, f) {
  if (is.null(seed)) {
    return(as.numeric(sample.int(.Machine$integer.max, 1L)))
  }
  .cknum(seed, "seed", f, min = 0, max = .Machine$integer.max, whole = TRUE)
}
