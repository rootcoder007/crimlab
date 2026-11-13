# The cost of input handling: for each exported function, the time of the
# whole R call and of the kernel alone, at a small and a large size, with
# named and with unnamed inputs. The kernel time comes from a build made by
# dev/kcost.sh, which times every kernel call of the binding; the rest of
# the call (the R closure, .Call, argument checks, name matching, building
# the result) is "outside".
#   Rscript dev/cost.R OUT LIB [LIB0]
# LIB holds the instrumented build of the working tree, LIB0 (optional)
# one of an earlier revision: each is timed in its own R process and OUT
# (a CSV) gets one row per case with full and kern seconds per call for
# each, and out = 1 - kern / full. Columns ending in 0 are LIB0's.
a <- commandArgs(TRUE)
if (identical(a[1], "--one")) {
  .libPaths(c(a[2], .libPaths()))
  suppressMessages(library(crimlab))
  grDevices::pdf(NULL)
  kt <- get("C_clab_kt", asNamespace("crimlab"))
  nm <- function(k, x) stats::setNames(x, paste0("g", seq_len(k)))
  un <- function(k, x) x
  # the 1-minute load average, as uptime reports it
  ld <- function() {
    x <- tryCatch(readLines("/proc/loadavg", n = 1L), error = function(e) "",
                  warning = function(w) "")
    if (!length(x) || !nzchar(x)) {
      x <- gsub("[{}]", "", paste(suppressWarnings(system2("sysctl",
        c("-n", "vm.loadavg"), stdout = TRUE, stderr = FALSE)), collapse = " "))
    }
    as.numeric(strsplit(trimws(x), " +")[[1]][1])
  }
  # a case is timed again (up to 3 times) when the load rose by more than
  # 1.5 while it was timed; the loads before and after are kept
  tl <- function(f) {
    for (try in 1:3) {
      l0 <- ld()
      r <- tm(f)
      l1 <- ld()
      if (!(l1 - l0 > 1.5)) break
    }
    c(r, l0, l1)
  }
  # a run repeats the call until it lasts 0.2 s; the result is the median
  # of five runs, per call, for the call and for the kernel alone
  tm <- function(f) {
    f()
    k <- 1
    repeat {
      t <- system.time(for (i in seq_len(k)) f())[["elapsed"]]
      if (t >= 0.2) break
      k <- k * 4
    }
    r <- replicate(5, {
      .Call(kt, TRUE)
      t <- system.time(for (i in seq_len(k)) f())[["elapsed"]]
      c(t, .Call(kt, TRUE)) / k
    })
    apply(r, 1, stats::median)
  }
  # COSTDRY=1: each call made once, untimed (a check of the cases)
  if (nzchar(Sys.getenv("COSTDRY"))) tm <- function(f) {
    f()
    c(NA_real_, NA_real_)
  }
  # six alleles a locus; the child drawn from the mother and the father
  dna <- function(n, named) {
    set.seed(1)
    fq <- t(replicate(n, {
      f <- stats::runif(6)
      f / sum(f)
    }))
    s <- function() t(replicate(n, sample(6, 2, TRUE)))
    mo <- s()
    fa <- s()
    ch <- cbind(mo[cbind(seq_len(n), sample(2, n, TRUE))],
                fa[cbind(seq_len(n), sample(2, n, TRUE))])
    fr <- function(g) matrix(fq[cbind(rep(seq_len(n), 2), c(g))], n)
    g1 <- s()
    g <- list(g1 = g1, g2 = s(), f = fr(g1), ch = ch, fa = fa, mo = mo,
              fc = fr(ch))
    if (named) g <- lapply(g, function(m) {
      rownames(m) <- paste0("L", seq_len(n))
      m
    })
    g
  }
  graph <- function(n) {
    set.seed(n)
    w <- (1 - stats::runif(n))^(-1 / 1.5)
    f <- sample.int(n, 3 * n, TRUE, w)
    t <- sample.int(n, 3 * n, TRUE, w)
    lab <- paste0("n", seq_len(n))
    list(f = lab[f], t = lab[t], w = stats::runif(3 * n, 0.5, 2), nodes = lab)
  }
  # fn, n, of, and a maker: function(n, h) giving the call as a closure,
  # h = TRUE for named inputs (or, for the networks, the node list given)
  cs <- list(
    list("vodtst", c(5, 1e5), "groups", function(n, h) {
      g <- if (h) nm else un
      d <- g(n, stats::runif(n, 10, 1e3)); e <- g(n, stats::runif(n, 10, 1e3))
      function() vodtst(d, e)
    }),
    list("orrisk", 1, "calls", function(n, h) function() orrisk(3, 0.4, 0.5)),
    list("outtst", c(1e3, 1e6), "values", function(n, h) {
      q <- stats::runif(n); w <- stats::runif(n)
      if (h) names(q) <- names(w) <- paste0("p", seq_len(n))
      function() outtst(q, 0.5, w)
    }),
    list("dspbm", c(5, 1e5), "groups", function(n, h) {
      g <- if (h) nm else un
      p <- g(n, stats::runif(n, 1e3, 1e5)); cn <- g(n, stats::runif(n, 10, 1e3))
      f <- g(n, stats::runif(n, 1, 50))
      function() dspbm(p, cn, f)
    }),
    list("dyadrt", c(5, 300), "groups", function(n, h) {
      m <- matrix(stats::runif(n * n, 0, 100), n)
      p <- stats::runif(n, 1e3, 1e5)
      if (h) {
        names(p) <- paste0("g", seq_len(n))
        dimnames(m) <- list(names(p), names(p))
      }
      function() dyadrt(m, p)
    }),
    list("gtmp", c(20, 1e5), "loci", function(n, h) {
      p <- stats::runif(n, 0.02, 0.3); q <- stats::runif(n, 0.02, 0.3)
      hm <- stats::runif(n) < 0.2
      if (h) names(p) <- names(q) <- names(hm) <- paste0("L", seq_len(n))
      function() gtmp(p, q, hm, 0.01)
    }),
    list("kinidx", c(20, 1e5), "loci", function(n, h) {
      g <- dna(n, h)
      function() kinidx(g$g1, g$g2, g$f, "sibling")
    }),
    list("patidx", c(20, 1e5), "loci", function(n, h) {
      g <- dna(n, h)
      function() patidx(g$ch, g$fa, g$fc, g$mo)
    }),
    list("mixlr", c(4, 8), "alleles", function(n, h) {
      p <- stats::setNames(stats::runif(n, 0.02, 0.08), seq_len(n))
      al <- names(p)
      function() mixlr(p, character(0), 0, al[-(1:2)], (n - 2) / 2)
    }),
    list("cllr", c(1e3, 1e6), "ratios", function(n, h) {
      s <- exp(stats::rnorm(n / 2, 1, 1.5)); d <- exp(stats::rnorm(n / 2, -1, 1.5))
      function() cllr(s, d)
    }),
    list("cllrmin", c(1e3, 1e6), "ratios", function(n, h) {
      s <- exp(stats::rnorm(n / 2, 1, 1.5)); d <- exp(stats::rnorm(n / 2, -1, 1.5))
      function() cllrmin(s, d)
    }),
    list("lrcal", c(1e3, 1e6), "ratios", function(n, h) {
      s <- exp(stats::rnorm(n / 2, 1, 1.5)); d <- exp(stats::rnorm(n / 2, -1, 1.5))
      function() lrcal(s, d)
    }),
    list("eceplt", c(1e3, 1e5), "ratios", function(n, h) {
      s <- exp(stats::rnorm(n / 2, 1, 1.5)); d <- exp(stats::rnorm(n / 2, -1, 1.5))
      function() eceplt(s, d)
    }),
    list("tipplt", c(1e3, 1e5), "ratios", function(n, h) {
      s <- exp(stats::rnorm(n / 2, 1, 1.5)); d <- exp(stats::rnorm(n / 2, -1, 1.5))
      function() tipplt(s, d)
    }),
    list("lrpost", c(1e3, 1e6), "ratios", function(n, h) {
      l <- exp(stats::rnorm(n, 0, 5))
      if (h) names(l) <- paste0("c", seq_len(n))
      function() lrpost(l, 1e-3)
    }),
    list("dbsrch", 1, "calls", function(n, h) function() dbsrch(1e-6, 1e5, 3e7)),
    list("clrsmp", c(5, 1e5), "offence types", function(n, h) {
      g <- if (h) nm else un
      ra <- g(n, stats::runif(n, 10, 1e3)); rb <- g(n, stats::runif(n, 10, 1e3))
      ca <- g(n, ra * stats::runif(n)); cb <- g(n, rb * stats::runif(n))
      function() clrsmp(ca, ra, cb, rb)
    }),
    list("expbnd", c(5, 1e5), "groups", function(n, h) {
      g <- if (h) nm else un
      y <- g(n, stats::runif(n, 10, 1e3)); m <- g(n, stats::runif(n, 1e3, 1e4))
      function() expbnd(y, m, 1.5)
    }),
    list("misrd", 1, "calls", function(n, h) function() misrd(0.3, 0.1, 0.4, 0.8, 0.9)),
    list("rtmhs", c(1e3, 1e6), "places", function(n, h) {
      mu <- stats::rgamma(n, 2, 0.5); x1 <- stats::rpois(n, mu); x2 <- stats::rpois(n, mu)
      if (h) names(x1) <- names(x2) <- paste0("p", seq_len(n))
      cut <- stats::quantile(x1, 0.9, names = FALSE)
      function() rtmhs(x1, x2, cut)
    }),
    list("incap", c(5, 1e5), "groups (share)", function(n, h) {
      g <- if (h) nm else un
      sh <- stats::runif(n)
      l <- g(n, stats::runif(n, 0.5, 20)); sh <- g(n, sh / sum(sh))
      function() incap(l, 0.1, 2, share = sh)
    }),
    list("incap", c(5, 1e5), "groups (no share)", function(n, h) {
      l <- (if (h) nm else un)(n, stats::runif(n, 0.5, 20))
      function() incap(l, 0.1, 2)
    }),
    list("incapr", c(10, 1e5), "career periods", function(n, h) {
      l <- sort(stats::runif(n, 0, 20), decreasing = TRUE)
      function() incapr(l, 2, floor(n / 2), 0.25)
    }),
    list("sellab", c(1e3, 1e6), "defendants", function(n, h) {
      k <- stats::runif(n); y <- stats::rbinom(n, 1, k); rel <- k < 0.7; rule <- k < 0.8
      if (h) names(y) <- names(rel) <- names(rule) <- paste0("d", seq_len(n))
      function() sellab(y, rel, rule)
    }),
    list("losbnd", c(1e3, 1e6), "stays", function(n, h) {
      d <- stats::rexp(n / 2, 1 / 60); o <- stats::rexp(n / 2, 1 / 90)
      function() losbnd(d, o)
    }),
    list("sntbnd", c(1e3, 1e6), "people", function(n, h) {
      z <- sample(c("community", "custody"), n, TRUE, c(0.6, 0.4))
      y <- stats::rbinom(n, 1, ifelse(z == "custody", 0.55, 0.35))
      if (h) names(y) <- names(z) <- paste0("p", seq_len(n))
      function() sntbnd(y, z, "custody")
    }),
    list("hzdsel", 1, "calls", function(n, h) {
      a <- c(0.5, 0.8); b <- c(0.9, 0.9)
      if (h) names(a) <- names(b) <- c("control", "treated")
      function() hzdsel(0.5, 0.5, 0.1, a, b)
    }),
    list("proof", 1, "calls", function(n, h) function() proof(c("cllr", "gtmp")))
  )
  net <- list(
    ntdeg = function(g, h) function() ntdeg(g$f, g$t, w = g$w, nodes = if (h) g$nodes),
    ntbtw = function(g, h) function() ntbtw(g$f, g$t, nodes = if (h) g$nodes),
    ntcls = function(g, h) function() ntcls(g$f, g$t, nodes = if (h) g$nodes),
    ntevc = function(g, h) function() ntevc(g$f, g$t, w = g$w, nodes = if (h) g$nodes),
    ntcore = function(g, h) function() ntcore(g$f, g$t, nodes = if (h) g$nodes),
    ntcmp = function(g, h) function() ntcmp(g$f, g$t, nodes = if (h) g$nodes, dir = TRUE),
    ntlvn = function(g, h) function() ntlvn(g$f, g$t, w = g$w, nodes = if (h) g$nodes, seed = 1),
    ntmod = function(g, h) {
      memb <- stats::setNames(rep_len(1:20, length(g$nodes)), g$nodes)
      function() ntmod(g$f, g$t, memb, w = g$w, nodes = if (h) g$nodes)
    },
    ntbrk = function(g, h) {
      grp <- stats::setNames(rep_len(c("a", "b", "c"), length(g$nodes)), g$nodes)
      function() ntbrk(g$f, g$t, grp, nodes = if (h) g$nodes, dir = TRUE)
    },
    ntchg = function(g, h) {
      # the largest component: the bounds need a connected graph
      cm <- ntcmp(g$f, g$t)
      big <- cm$node[cm$weak == names(which.max(table(cm$weak)))]
      e <- g$f %in% big & g$t %in% big
      f <- g$f[e]
      t <- g$t[e]
      w <- g$w[e]
      function() ntchg(f, t, w = w, nodes = if (h) big, maxit = 2e3)
    },
    coofg = function(g, h) {
      set.seed(1)
      k <- sample(1:4, length(g$nodes), TRUE)
      inc <- rep(seq_along(k), k)
      off <- sample(g$nodes, sum(k), TRUE)
      function() coofg(inc, off)
    }
  )
  for (x in names(net)) {
    big <- if (x %in% c("ntbtw", "ntcls", "ntchg")) 1e4 else 1e5
    cs[[length(cs) + 1L]] <- list(x, c(1e3, big), "nodes", local({
      f <- net[[x]]
      function(n, h) f(graph(n), h)
    }))
  }
  rows <- list()
  for (k in cs) {
    for (n in k[[2]]) {
      # functions whose inputs carry no names to match: one way only
      one <- c("orrisk", "dbsrch", "misrd", "proof", "mixlr", "cllr",
               "cllrmin", "lrcal", "eceplt", "tipplt", "incapr", "losbnd",
               "coofg")
      for (h in if (k[[1]] %in% one) FALSE else c(FALSE, TRUE)) {
        set.seed(20261112)
        r <- tl(k[[4]](n, h))
        rows[[length(rows) + 1L]] <- data.frame(
          fn = k[[1]], of = k[[3]], n = n,
          named = if (k[[3]] == "nodes") (if (h) "nodes given" else "no nodes") else h,
          full = r[1], kern = r[2], l0 = r[3], l1 = r[4])
        cat(sprintf("%-8s %-18s %8g %-11s %10.3g %10.3g\n", k[[1]], k[[3]], n,
                    rows[[length(rows)]]$named, r[1], r[2]))
      }
    }
  }
  saveRDS(do.call(rbind, rows), a[3])
  quit(save = "no")
}

out <- a[1]
libs <- a[-1]
run <- function(lib) {
  f <- tempfile(fileext = ".rds")
  st <- system2(file.path(R.home("bin"), "Rscript"),
                c("dev/cost.R", "--one", lib, f))
  if (st != 0) stop("cost run failed for ", lib)
  readRDS(f)
}
r <- run(libs[1])
r$out <- 1 - r$kern / r$full
r$mach <- Sys.info()[["nodename"]]
sh <- function(cmd, args) {
  tryCatch(suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = FALSE)),
           error = function(e) character(0))
}
cpu <- c(sub("^Model name: *", "", grep("^Model name:", sh("lscpu", character()), value = TRUE)),
         sh("sysctl", c("-n", "machdep.cpu.brand_string")))
r$cpu <- paste(c(cpu, "")[1], Sys.info()[["machine"]])
if (length(libs) > 1) {
  r0 <- run(libs[2])
  r$full0 <- r0$full
  r$kern0 <- r0$kern
  r$out0 <- 1 - r0$kern / r0$full
  r$l00 <- r0$l0
  r$l10 <- r0$l1
}
r$rver <- R.version.string
r$date <- format(Sys.Date())
utils::write.csv(r, out, row.names = FALSE)
