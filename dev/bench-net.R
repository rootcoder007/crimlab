# Speed of the network functions against every package found that computes
# the same measure (dev/cpart.R records the search): igraph, sna, tnet and
# netmem. Fixed seeded graphs of 1e3, 1e4 and 1e5 nodes with heavy-tailed
# degrees (Chung-Lu: 3n edges, ends drawn with probability proportional to
# Pareto weights of exponent 2.5). Each counterpart's values are checked
# against crimlab's on the graph timed before it is timed; a disagreement
# fails the run. Nothing in the output is typed by hand.
#   Rscript dev/bench-net.R OUTDIR             one run, into OUTDIR/net.csv
#   Rscript dev/bench-net.R --merge RUNDIR ... the runs merged into
#                                              inst/bench/net.csv and the
#                                              README lists
# Runs are started by dev/benchrun.sh on a quiet machine, as dev/bench.R's.
# Columns: seconds for crimlab (from labels: each call builds the graph),
# for the reference on its own graph object already built (ref_s) and for
# the reference including building that object from the same labels
# (ref_build_s), each the median of three timings (one when a timing takes
# over 20 s; a fast call is repeated until a timing lasts 0.2 s and
# divided); speed = ref_s / crimlab_s (above 1: crimlab is faster); faster
# is "yes" or "no"; agree is "yes" when the values were checked equal
# (tolerance in tol), "not run" otherwise; note gives the reason a measure
# is not faster or was not run. The merge keeps, per row, the run with the
# median speed, and lo and hi, the smallest and largest speed over the runs.
#
# sna, tnet and netmem: sna reads an edge list (sparse); netmem and sna's
# brokerage need a dense nodes x nodes matrix and tnet's betweenness and
# closeness run an R loop over nodes, so these are timed on the 1e3-node
# graph only. Eigenvector centrality is timed on the largest component of
# each graph, where every package agrees on it (on a disconnected graph
# they treat the smaller components differently).
a <- commandArgs(TRUE)
merge <- identical(a[1], "--merge")
if (!merge) {
out <- if (length(a)) a[1] else "net-run"
suppressMessages({
  library(crimlab)
  library(igraph)
})
has <- function(p) requireNamespace(p, quietly = TRUE)
for (p in c("crimlab", "igraph", "sna", "tnet", "netmem")) {
  cat(sprintf("%-8s %s\n", p, tryCatch(paste(find.package(p),
    utils::packageVersion(p)), error = function(e) "not installed")))
}

tm <- function(expr) {
  e <- substitute(expr)
  p <- parent.frame()
  f <- function() {
    k <- 1
    repeat {
      t <- system.time(for (i in seq_len(k)) eval(e, p))[["elapsed"]]
      if (t >= 0.2) return(t / k)
      k <- k * 4
    }
  }
  t1 <- f()
  if (t1 > 20) return(t1)
  stats::median(c(t1, f(), f()))
}

graph <- function(n) {
  set.seed(n)
  w <- (1 - stats::runif(n))^(-1 / 1.5)
  m <- 3 * n
  f <- sample.int(n, m, replace = TRUE, prob = w)
  t <- sample.int(n, m, replace = TRUE, prob = w)
  lab <- paste0("n", seq_len(n))
  list(f = lab[f], t = lab[t], w = stats::runif(m, 0.5, 2), nodes = lab)
}

# equal to a relative tolerance (or both 0); NaN where both are NaN
eq <- function(x, y, tol) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  length(x) == length(y) && identical(is.nan(x), is.nan(y)) &&
    all(x == y | abs(x - y) <= tol * pmax(abs(x), abs(y)), na.rm = TRUE)
}
# the same partition, whatever the labels
part <- function(x, y) {
  length(x) == length(y) &&
    length(unique(paste(x, y))) == length(unique(x)) &&
    length(unique(x)) == length(unique(y))
}

rows <- list()
add <- function(meas, ref, n, m, cs, rs, rbs, ok, tol, note = "", l0 = NA,
                l1 = NA) {
  sp <- if (is.na(rs) || is.na(cs)) NA else rs / cs
  fast <- if (is.na(sp)) NA else if (sp > 1) "yes" else "no"
  if (identical(fast, "no") && note == "") {
    note <- if (!is.na(rbs) && cs < rbs) {
      paste("slower than the measure alone; faster than building the",
            "graph from labels and measuring")
    } else {
      "not investigated"
    }
  }
  agree <- if (is.na(ok)) "not run" else if (ok) "yes" else "no"
  rows[[length(rows) + 1L]] <<- data.frame(
    measure = meas, reference = ref, nodes = n, edges = m, crimlab_s = cs,
    ref_s = rs, ref_build_s = rbs, speed = sp, faster = fast, agree = agree,
    tol = tol, note = note, l0 = l0, l1 = l1
  )
  cat(sprintf("%-7s %-40s %7g %9.3g %9.3g %7.3g %s %s\n", meas, ref, n, cs,
              rs, sp, agree, note))
}
# a counterpart: check it and time it with crimlab's call (cs, a function
# returning seconds), or record why it was not timed (skip). A pair whose load
# rose by more than 1.5 while it was timed is timed again (up to 3 times).
cp <- function(meas, ref, n, m, cs, ok, tol, run, build = NA, note = "",
               skip = "") {
  if (!isTRUE(ok)) {
    add(meas, ref, n, m, NA, NA, NA, if (is.na(ok)) NA else FALSE, tol,
        if (is.na(ok)) skip else "values differ")
    return(invisible())
  }
  for (try in 1:3) {
    l0 <- ld()
    c1 <- cs()
    r <- run()
    l1 <- ld()
    if (!(l1 - l0 > 1.5)) break
  }
  add(meas, ref, n, m, c1, r, if (is.na(build)) NA else build + r, TRUE, tol,
      note, l0, l1)
}
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
load0 <- ld()
start <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")

# NETSIZES (e.g. "1e3") limits the sizes, for a check of the script
sizes <- as.numeric(strsplit(Sys.getenv("NETSIZES", "1e3 1e4 1e5"), " +")[[1]])
for (n in sizes) {
  g <- graph(n)
  m <- length(g$f)
  small <- n <= 1e3
  bld <- function(dir = FALSE) {
    x <- graph_from_data_frame(data.frame(g$f, g$t, weight = g$w),
                               directed = dir, vertices = g$nodes)
    simplify(x, edge.attr.comb = list(weight = "sum"))
  }
  tb <- tm(bld())
  ig <- bld()
  igd <- bld(TRUE)
  tbd <- tm(bld(TRUE))
  # sna's edge list: both directions of each undirected edge, node ids in
  # node order; tnet's: the same with weights
  el <- as_edgelist(ig, names = FALSE)
  sel <- rbind(el, el[, 2:1])
  attr(sel, "n") <- n
  dl <- as_edgelist(igd, names = FALSE)
  attr(dl, "n") <- n
  tel <- cbind(sel, rep(E(ig)$weight, 2))
  A <- if (small) as.matrix(as_adjacency_matrix(ig, sparse = FALSE))
  W <- if (small) as.matrix(as_adjacency_matrix(ig, attr = "weight", sparse = FALSE))
  Ad <- if (small) as.matrix(as_adjacency_matrix(igd, sparse = FALSE))
  if (small) {
    A[A != 0] <- 1
    Ad[Ad != 0] <- 1
  }
  dense <- "not run: needs a dense nodes x nodes matrix"
  slow <- "not run: an R loop over nodes"
  cat("n =", n, "built\n")

  x <- ntdeg(g$f, g$t, w = g$w, nodes = g$nodes)
  cs <- function() tm(ntdeg(g$f, g$t, w = g$w, nodes = g$nodes))
  cp("ntdeg", "igraph::degree + strength", n, m, cs,
     eq(x$deg, degree(ig), 0) && eq(x$wdeg, strength(ig), 1e-12), 1e-12,
     function() tm({degree(ig); strength(ig)}), tb)
  cp("ntdeg", "sna::degree", n, m, cs,
     eq(x$deg, sna::degree(sel, gmode = "graph"), 0), 0,
     function() tm(sna::degree(sel, gmode = "graph")))
  td <- if (n <= 1e4) tnet::degree_w(tel, measure = c("degree", "output"))
  cp("ntdeg", "tnet::degree_w", n, m, cs,
     if (n <= 1e4) eq(x$deg[td[, 1]], td[, "degree"], 0) &&
       eq(x$wdeg[td[, 1]], td[, "output"], 1e-12) else NA, 1e-12,
     function() tm(tnet::degree_w(tel, measure = c("degree", "output"))),
     skip = "not run: tnet reshapes the edge list in R, minutes at 1e5")
  cp("ntdeg", "netmem::gen_degree", n, m, cs,
     if (small) eq(x$deg, netmem::gen_degree(A, digraph = FALSE, loops = FALSE), 0) &&
       eq(x$wdeg, netmem::gen_degree(W, weighted = TRUE, digraph = FALSE,
                                     loops = FALSE, alpha = 1), 1e-12) else NA,
     1e-12, function() tm({
       netmem::gen_degree(A, digraph = FALSE, loops = FALSE)
       netmem::gen_degree(W, weighted = TRUE, digraph = FALSE, loops = FALSE,
                          alpha = 1)
     }), skip = dense)

  if (n <= 1e4) {
    x <- ntbtw(g$f, g$t, nodes = g$nodes)$btw
    cs <- function() tm(ntbtw(g$f, g$t, nodes = g$nodes))
    cp("ntbtw", "igraph::betweenness", n, m, cs,
       eq(x, betweenness(ig, weights = NA), 1e-9), 1e-9,
       function() tm(betweenness(ig, weights = NA)), tb)
    cp("ntbtw", "sna::betweenness", n, m, cs,
       eq(x, sna::betweenness(sel, gmode = "graph"), 1e-9), 1e-9,
       function() tm(sna::betweenness(sel, gmode = "graph")))
    tbw <- if (small) tnet::betweenness_w(tel, alpha = 0)
    cp("ntbtw", "tnet::betweenness_w", n, m, cs,
       if (small) eq(x[tbw[, 1]], tbw[, 2], 1e-9) else NA, 1e-9,
       function() tm(tnet::betweenness_w(tel, alpha = 0)), skip = slow)
    cp("ntbtw", "netmem::betweenness_centrality", n, m, cs,
       if (small) eq(x, netmem::betweenness_centrality(A, digraph = FALSE), 1e-9) else NA,
       1e-9, function() tm(netmem::betweenness_centrality(A, digraph = FALSE)),
       skip = dense)

    x <- ntcls(g$f, g$t, nodes = g$nodes)
    cs <- function() tm(ntcls(g$f, g$t, nodes = g$nodes))
    cp("ntcls", "igraph::closeness + harmonic_centrality", n, m, cs,
       eq(x$cls, closeness(ig, weights = NA), 1e-12) &&
         eq(x$hrm, harmonic_centrality(ig, weights = NA), 1e-12), 1e-12,
       function() tm({
         closeness(ig, weights = NA)
         harmonic_centrality(ig, weights = NA)
       }), tb)
    tc <- if (small) tnet::closeness_w(tel, alpha = 0, gconly = FALSE)
    cp("ntcls", "tnet::closeness_w", n, m, cs,
       if (small) eq(x$hrm[tc[, 1]], tc[, "closeness"], 1e-12) else NA, 1e-12,
       function() tm(tnet::closeness_w(tel, alpha = 0, gconly = FALSE)),
       skip = slow)
    cp("ntcls", "netmem::closeness_centrality", n, m, cs,
       if (small) suppressWarnings(
         eq(x$cls, netmem::closeness_centrality(A, digraph = FALSE), 1e-12) &&
         eq(x$hrm, netmem::closeness_centrality(A, digraph = FALSE,
                                                harmonic = TRUE), 1e-12)) else NA,
       1e-12, function() tm(suppressWarnings({
         netmem::closeness_centrality(A, digraph = FALSE)
         netmem::closeness_centrality(A, digraph = FALSE, harmonic = TRUE)
       })), skip = dense)
  } else {
    for (x in c("ntbtw", "ntcls")) {
      add(x, "igraph", n, m, NA, NA, NA, NA, NA,
          "not run: O(nodes x edges), about 3e10 steps for any package")
    }
  }

  # eigenvector centrality on the largest component
  big <- components(ig)
  keep <- names(big$membership)[big$membership == which.max(big$csize)]
  e <- g$f %in% keep & g$t %in% keep
  gc <- induced_subgraph(ig, keep)
  x <- ntevc(g$f[e], g$t[e], w = g$w[e], nodes = keep)$cent$evc
  cs <- function() tm(ntevc(g$f[e], g$t[e], w = g$w[e], nodes = keep))
  cp("ntevc", "igraph::eigen_centrality", n, sum(e), cs,
     eq(x, eigen_centrality(gc)$vector, 1e-6), 1e-6,
     function() tm(eigen_centrality(gc)))
  Wc <- if (small) as.matrix(as_adjacency_matrix(gc, attr = "weight", sparse = FALSE))
  cp("ntevc", "netmem::eigenvector_centrality", n, sum(e), cs,
     if (small) eq(x, netmem::eigenvector_centrality(Wc, digraph = FALSE,
                                                     weighted = TRUE)$vector, 1e-6) else NA,
     1e-6, function() tm(netmem::eigenvector_centrality(Wc, digraph = FALSE,
                                                        weighted = TRUE)),
     skip = dense)

  x <- ntcore(g$f, g$t, nodes = g$nodes)$core
  cs <- function() tm(ntcore(g$f, g$t, nodes = g$nodes))
  cp("ntcore", "igraph::coreness", n, m, cs, eq(x, coreness(ig), 0), 0,
     function() tm(coreness(ig)), tb)
  cp("ntcore", "sna::kcores", n, m, cs,
     eq(x, sna::kcores(sel, mode = "graph", cmode = "freeman"), 0), 0,
     function() tm(sna::kcores(sel, mode = "graph", cmode = "freeman")))
  cp("ntcore", "netmem::k_core", n, m, cs,
     if (small) eq(x, netmem::k_core(A), 0) else NA, 0,
     function() tm(netmem::k_core(A)), skip = dense)

  x <- ntcmp(g$f, g$t, nodes = g$nodes, dir = TRUE)
  cs <- function() tm(ntcmp(g$f, g$t, nodes = g$nodes, dir = TRUE))
  cp("ntcmp", "igraph::components weak + strong", n, m, cs,
     part(x$weak, components(igd, "weak")$membership) &&
       part(x$strong, components(igd, "strong")$membership), 0,
     function() tm({components(igd, "weak"); components(igd, "strong")}), tbd)
  cp("ntcmp", "sna::component.dist weak + strong", n, m, cs,
     if (small) part(x$weak, sna::component.dist(dl, connected = "weak")$membership) &&
       part(x$strong, sna::component.dist(dl, connected = "strong")$membership) else NA,
     0, function() tm({
       sna::component.dist(dl, connected = "weak")
       sna::component.dist(dl, connected = "strong")
     }), skip = paste("not run: sna finds strong components from a list of",
                      "every reachable pair, nodes^2 rows"))
  cp("ntcmp", "netmem::components_id weak + strong", n, m, cs,
     if (small) part(x$weak, netmem::components_id(Ad, "weak")$components) &&
       part(x$strong, netmem::components_id(Ad, "strong")$components) else NA, 0,
     function() tm({
       netmem::components_id(Ad, "weak")
       netmem::components_id(Ad, "strong")
     }), skip = dense)

  # Louvain is a heuristic: the check is the modularity reached
  x <- ntlvn(g$f, g$t, w = g$w, nodes = g$nodes, seed = 1)
  set.seed(1)
  il <- cluster_louvain(ig)
  q <- ntmod(g$f, g$t, x$memb$comm, w = g$w, nodes = g$nodes)
  cs <- function() tm(ntlvn(g$f, g$t, w = g$w, nodes = g$nodes, seed = 1))
  cp("ntlvn", "igraph::cluster_louvain", n, m, cs,
     abs(q - modularity(il)) <= 0.02, 0.02,
     function() tm({set.seed(1); cluster_louvain(ig)}), tb,
     note = "values: the modularity reached, within 0.02")

  memb <- rep_len(1:20, n)
  x <- ntmod(g$f, g$t, memb, w = g$w, nodes = g$nodes)
  cs <- function() tm(ntmod(g$f, g$t, memb, w = g$w, nodes = g$nodes))
  cp("ntmod", "igraph::modularity", n, m, cs, eq(x, modularity(ig, memb), 1e-12),
     1e-12, function() tm(modularity(ig, memb)), tb)
  cp("ntmod", "netmem::modularity_score", n, m, cs,
     if (small) eq(x, netmem::modularity_score(W, memb, weighted = TRUE), 1e-12) else NA,
     1e-12, function() tm(netmem::modularity_score(W, memb, weighted = TRUE)),
     skip = dense)

  grp <- rep_len(c("a", "b", "c"), n)
  x <- ntbrk(g$f, g$t, grp, nodes = g$nodes, dir = TRUE)
  xb <- as.matrix(x[c("coord", "itin", "repr", "gate", "liais")])
  cs <- function() tm(ntbrk(g$f, g$t, grp, nodes = g$nodes, dir = TRUE))
  cp("ntbrk", "sna::brokerage", n, m, cs,
     if (small) eq(xb, sna::brokerage(Ad, grp)$raw.nli[, 1:5], 0) else NA, 0,
     function() tm(sna::brokerage(Ad, grp)), skip = dense)
  cp("ntbrk", "netmem::brokerage_roles", n, m, cs,
     if (small) eq(xb, as.matrix(netmem::brokerage_roles(Ad, grp)[,
       c("coordinator", "consultant", "representative", "gatekeeper", "liaison")]), 0)
     else NA, 0, function() tm(netmem::brokerage_roles(Ad, grp)), skip = dense)

  if (n <= 1e4) {
    cv <- ntchg(g$f[e], g$t[e], w = g$w[e], maxit = 2e4)$conv
    add("ntchg", "none", n, sum(e),
        tm(ntchg(g$f[e], g$t[e], w = g$w[e], maxit = 2e4)), NA, NA, NA, NA,
        paste("no other package found; converged within 2e4 steps:", cv))
  }

  # co-offending: incidents of 1 to 4 offenders drawn by weight
  set.seed(n + 1)
  k <- sample(1:4, n, replace = TRUE)
  inc <- rep(seq_len(n), k)
  off <- sample(g$nodes, sum(k), replace = TRUE)
  x <- coofg(inc, off)$edges
  xk <- paste(pmin(x$from, x$to), pmax(x$from, x$to))
  xw <- stats::setNames(x$w, xk)
  pw <- function(a, b, w) {
    kk <- paste(pmin(a, b), pmax(a, b))
    length(kk) == length(xk) && setequal(kk, xk) && eq(xw[kk], w, 0)
  }
  bip <- function() {
    # an offender named twice in one incident is one tie, as in coofg()
    y <- graph_from_data_frame(unique(data.frame(paste0("i", inc), off)),
                               directed = FALSE)
    V(y)$type <- startsWith(V(y)$name, "n")
    bipartite_projection(y, which = "true", multiplicity = TRUE)
  }
  pb <- bip()
  pe <- as_edgelist(pb)
  cs <- function() tm(coofg(inc, off))
  cp("coofg", "igraph::bipartite_projection (with building)", n, length(inc), cs,
     pw(pe[, 1], pe[, 2], E(pb)$weight), 0, function() tm(bip()), 0)
  tp <- function() {
    o <- match(off, g$nodes)
    tnet::projecting_tm(unique(cbind(o, inc)), method = "sum")
  }
  tw <- if (n <= 1e4) tp()
  cp("coofg", "tnet::projecting_tm (with building)", n, length(inc), cs,
     if (n <= 1e4) {
       h <- tw[, 1] < tw[, 2]
       pw(g$nodes[tw[h, 1]], g$nodes[tw[h, 2]], tw[h, 3])
     } else NA, 0, function() tm(tp()), 0,
     skip = "not run: tnet reshapes the edge list in R, minutes at 1e5")
  mp <- function() {
    b <- unclass(table(off, inc))
    b[b > 1] <- 1
    netmem::matrix_projection(b)$matrix2
  }
  cp("coofg", "netmem::matrix_projection (with building)", n, length(inc), cs,
     if (small) {
       p <- mp()
       diag(p) <- 0
       i <- which(upper.tri(p) & p > 0, arr.ind = TRUE)
       pw(rownames(p)[i[, 1]], colnames(p)[i[, 2]], p[i])
     } else NA, 0, function() tm(mp()), 0, skip = dense)
  cat("n =", n, "done\n")
}

load1 <- ld()
tab <- do.call(rbind, rows)
tab$platform <- paste(R.version$platform, "R", getRversion(), "igraph",
                      packageVersion("igraph"), "sna", packageVersion("sna"),
                      "tnet", packageVersion("tnet"), "netmem",
                      packageVersion("netmem"))
sh <- function(cmd, args) {
  tryCatch(suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = FALSE)),
           error = function(e) character(0))
}
cpu <- c(sub("^Model name: *", "", grep("^Model name:", sh("lscpu", character()), value = TRUE)),
         sh("sysctl", c("-n", "machdep.cpu.brand_string")))
# named as in the other speed tables (dev/bench-sp.R)
tab$mach <- Sys.info()[["nodename"]]
tab$cpu <- paste(c(cpu, "")[1], Sys.info()[["machine"]])
tab$rver <- R.version.string
tab$load0 <- load0
tab$load1 <- load1
tab$date <- start
dir.create(out, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, file.path(out, "net.csv"), row.names = FALSE)
if (any(tab$agree == "no")) {
  print(tab[tab$agree == "no", c("measure", "reference", "nodes")])
  quit(status = 1L)
}
quit(status = 0L)
}

## merge ------------------------------------------------------------------

rs <- lapply(a[-1], function(d) utils::read.csv(file.path(d, "net.csv"),
                                                 stringsAsFactors = FALSE))
key <- function(t) paste(t$measure, t$reference, t$nodes, sep = "\r")
for (t in rs[-1]) {
  if (!identical(key(t), key(rs[[1]]))) stop("the runs do not hold the same rows")
}
sp <- sapply(rs, function(t) t$speed)
tab <- rs[[1]]
for (i in seq_len(nrow(tab))) {
  r <- sp[i, ]
  if (all(is.na(r))) next
  tab[i, ] <- rs[[order(r)[ceiling(sum(!is.na(r)) / 2)]]][i, ]
}
tab$lo <- suppressWarnings(apply(sp, 1, min, na.rm = TRUE))
tab$hi <- suppressWarnings(apply(sp, 1, max, na.rm = TRUE))
tab$lo[!is.finite(tab$lo)] <- NA
tab$hi[!is.finite(tab$hi)] <- NA
utils::write.csv(tab, "inst/bench/net.csv", row.names = FALSE)

if (file.exists("README.md")) {
  lst <- function(x) if (length(x)) paste0("`", x, "()`", collapse = ", ") else "none"
  ok <- tab$agree == "yes"
  th <- unique(tab$measure[ok])
  sl <- unique(tab$measure[ok & tab$speed < 1.1])
  # with building where the package has a graph object to build (igraph);
  # without where its input is the edge list or matrix itself
  bs <- ifelse(is.na(tab$ref_build_s), tab$ref_s, tab$ref_build_s) / tab$crimlab_s
  slb <- unique(tab$measure[ok & bs < 1.1])
  no <- setdiff(unique(tab$measure), th)
  txt <- c(
    "<!-- netbench start: written by dev/bench-net.R -->",
    paste0("Networks, against igraph, sna, tnet and netmem: faster, by at least ",
           "a tenth, than every package that computes the same measure, at ",
           "every size timed, counting only the measure on a graph the other ",
           "package has already built: ", lst(setdiff(th, sl)), ". Not faster ",
           "at some size by that count: ", lst(intersect(th, sl)), "; of these, ",
           "counting the building of the other package's graph from the same ",
           "labels, which crimlab's time includes, faster at every size: ",
           lst(setdiff(intersect(th, sl), slb)), ". No other package found: ",
           lst(no), "."),
    "<!-- netbench end -->"
  )
  rd <- readLines("README.md")
  i <- grep("<!-- netbench start", rd, fixed = TRUE)
  j <- grep("<!-- netbench end -->", rd, fixed = TRUE)
  if (length(i) == 1L && length(j) == 1L && j > i) {
    rd <- c(rd[seq_len(i - 1L)], txt, rd[-seq_len(j)])
  } else {
    k <- grep("<!-- bench end -->", rd, fixed = TRUE)
    rd <- append(rd, c("", txt), after = k)
  }
  writeLines(rd, "README.md")
}
