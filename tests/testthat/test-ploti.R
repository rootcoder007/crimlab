own <- function() {
  data.frame(
    community = rep(c("Lakeside", "Riverbend"), each = 4),
    item = rep(c("Hopeful", "Satisfied", "Rarely hopeful",
                 "Low satisfaction"), 2),
    pct = c(64, 71, 12, 18, 58, 66, 16, 22),
    lo = c(59, 66, 9, 14, 52, 60, 12, 17),
    hi = c(69, 76, 15, 22, 64, 72, 20, 27)
  )
}

test_that("a gap is drawn as its ratio and interval, flagged rows starred", {
  d <- data.frame(group = c("First Nations", "Inuit", "Non-Indigenous"),
                  pct = c(45.8, 47.8, 61.5), lo = c(43.9, 43.1, 61),
                  hi = c(47.7, 52.5, 62), q = c("", "E", ""))
  g <- indgap(d, "group", "pct", "lo", "hi", qual = "q")
  p <- drw(plot(g))
  expect_identical(p$est, g$ratio)
  expect_identical(p$lo, g$ratiolo)
  expect_identical(p$label, c("First Nations", "Inuit *"))
})

test_that("wellness is drawn beside distress, one row per group", {
  w <- welprof(own(), "community", "item", "pct", "lo", "hi",
               well = "Hopeful", distress = "Rarely hopeful")
  x <- dualcont(w, "Hopeful", "Rarely hopeful")
  p <- drw(plot(x))
  expect_identical(p$est, c(x$wellness, x$distress))
  expect_identical(p$label, rep(x$group, 2))
  expect_identical(unique(p$group),
                   c("wellness: Hopeful", "distress: Rarely hopeful"))
})

test_that("the key of a grouped chart leaves the rows below it", {
  w <- welprof(own(), "community", "item", "pct", "lo", "hi",
               well = "Hopeful", distress = "Rarely hopeful")
  x <- dualcont(w, "Hopeful", "Rarely hopeful")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(x)
  # two rows are at 1 and 2: the plot region goes on above them
  expect_gt(graphics::par("usr")[4L], 2.5 + 0.04 * 2)
})

test_that("levels with intervals: index, prevalence, rating and profile", {
  w <- commidx(own(), "community", "item", "pct", "lo", "hi",
               w = c(Hopeful = 2, Satisfied = 1, "Rarely hopeful" = -1,
                     "Low satisfaction" = -1))
  expect_identical(drw(plot(w))$est, w$index)
  s <- stdprev(est = c(12, 9, 6, 8, 6, 4), lower = c(9, 6, 3, 7, 5, 3),
               upper = c(15, 12, 9, 9, 7, 5), std = c(500, 300, 200),
               grp = rep(c("Nation A", "Everyone else"), each = 3))
  expect_identical(drw(plot(s))$hi, s$upper)
  f <- welprof(own(), "community", "item", "pct", "lo", "hi",
               well = c("Hopeful", "Satisfied"),
               distress = c("Rarely hopeful", "Low satisfaction"))
  p <- drw(plot(f))
  expect_identical(p$est, f$estimate)
  expect_setequal(p$group, c("Lakeside", "Riverbend"))
})

test_that("balance draws its two panels on one page", {
  h <- data.frame(community = rep(c("Lakeside", "Riverbend"), each = 4),
                  domain = rep(c("body", "mind", "spirit", "family"), 2),
                  pct = c(62, 58, 70, 75, 55, 61, 66, 80),
                  lo = c(57, 53, 65, 70, 49, 55, 60, 75),
                  hi = c(67, 63, 75, 80, 61, 67, 72, 85))
  b <- holbal(h, "community", "domain", "pct", "lo", "hi", reps = 200)
  expect_identical(pages(plot(b)), 1L)
  p <- drw(plot(b))
  expect_named(p, c("level", "balance"))
  expect_identical(p$balance$est, b$balance)
})

test_that("Community Well-Being scores are drawn by component", {
  s <- cwbscr(incpercap = c(18000, 31000), highsch20p = c(58, 81),
              uni25p = c(6, 19), notcro = c(71, 95), nomajrep = c(62, 90),
              partic = c(55, 68), empl = c(78, 92))
  p <- drw(plot(s))
  expect_identical(nrow(p), 10L)
  expect_identical(p$est[p$group == "score"], s$cwb)
})

test_that("an object without its columns is refused", {
  expect_error(plot.indgap(list(group = "a")), "not a complete indgap")
  expect_error(plot.cwbscr(data.frame(income = 1)), "not a complete cwbscr")
})

test_that("periods given as text are drawn, each at its place", {
  s <- stkflw(people = c(5200, 4800, 4300),
              days = c(310000, 318000, 330000), prd = 2021:2023)
  expect_type(s$flow$period, "character")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p <- plot(s)
  # the last panel spans the three periods: before, the axis ran 0 to 1
  expect_gt(diff(graphics::par("usr")[1:2]), 1.9)
  expect_identical(p$period, s$flow$period)
  f <- stkflw(people = c(5200, 4800, 4300),
              days = c(310000, 318000, 330000),
              prd = c("2020/21", "2021/22", "2022/23"))
  plot(f)
  expect_gt(diff(graphics::par("usr")[1:2]), 1.9)
})

test_that("a chart of many rows keeps the label of every row", {
  x <- data.frame(group = sprintf("group %02d", 1:24), est = 1:24)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  d <- .pdots(x$group, x$est, x$est - 1, x$est + 1, "", "Value")
  expect_identical(d$label, x$group)
})

test_that("a pooled total is drawn in its own colour and shape", {
  d <- data.frame(group = c("First Nations", "Inuit", "Total Indigenous",
                            "Non-Indigenous"),
                  pct = c(45.8, 47.8, 49, 61.5), lo = c(43.9, 43.1, 47.5, 61),
                  hi = c(47.7, 52.5, 50.5, 62))
  g <- indgap(d, "group", "pct", "lo", "hi", pooled = "Total Indigenous")
  p <- drw(plot(g))
  expect_identical(p$group, c("distinct group", "distinct group",
                              "pooled total"))
})

test_that("a rating, a change and a health gap are drawn", {
  r <- data.frame(group = rep(c("Nation A", "Everyone else"), each = 4),
                  rating = rep(c("Poor", "Fair", "Good", "Excellent"), 2),
                  pct = c(10, 25, 40, 25, 6, 20, 44, 30),
                  lo = c(7, 21, 35, 21, 5, 18, 42, 28),
                  hi = c(13, 29, 45, 29, 7, 22, 46, 32))
  o <- ordcmp(r, "group", "rating", "pct", "lo", "hi",
              lvls = c("Poor", "Fair", "Good", "Excellent"),
              ref = "Everyone else")
  expect_identical(drw(plot(o))$est, o$probbet)
  h <- data.frame(community = "Lakeside",
                  item = rep(c("Good mental health", "Anxiety"), 2),
                  survey = rep(c("2018", "2023"), each = 2),
                  pct = c(55, 14, 61, 11), lo = c(50, 11, 56, 8),
                  hi = c(60, 17, 66, 14))
  ch <- healchg(h, "community", "item", "survey", "pct", "lo", "hi",
                from = "2018", to = "2023", well = "Good mental health",
                distress = "Anxiety")
  expect_identical(drw(plot(ch))$est, ch$change)
  le <- data.frame(period = "2017/2019",
                   group = c("First Nations", "Inuit", "Non-Indigenous"),
                   estimate = c(74.2, 70.9, 81.8),
                   lower = c(73.6, 69.1, 81.6),
                   upper = c(74.8, 72.7, 82.0), quality = "ok")
  tg <- trchlthgap(health = NULL, lifeexp = le, infmort = NULL,
                   suiide = NULL)
  p <- drw(plot(tg))
  expect_identical(p$est, tg$disratio)
  expect_identical(unique(p$group), "life expectancy 2017/2019")
})

test_that("Community Well-Being from a census profile is drawn", {
  nm <- c(
    "Total - Age groups of the population - 100% data",
    paste("Number of total income recipients aged 15 years and over in",
          "private households in 2020 - 25% sample data"),
    "Average total income in 2020 among recipients ($)",
    paste("Total - Highest certificate, diploma or degree for the",
          "population aged 25 years and over in private households - 25%",
          "sample data"),
    "  No certificate, diploma or degree",
    "  Bachelor's degree or higher",
    "Participation rate",
    "Unemployment rate"
  )
  x <- data.frame(
    GEO_NAME = "Lakeside", DGUID = "X0001", GENDER = "Total - Gender",
    AGE_GROUP = c(rep("Total - Age groups", 3), rep("25 years and over", 3),
                  rep("25 to 64 years", 2)),
    RESIDENCE = "Total - Residence on or off reserve",
    CHARACTERISTIC_ID = 1:8, CHARACTERISTIC_NAME = nm,
    `IND_IDENT (2)` = c(1200, 900, 38000, 700, 150, 120, 62, 14),
    SYMBOL_IND_IDENT = "",
    NON_IND = c(5000, 4000, 52000, 3500, 400, 900, 66, 6),
    SYMBOL_NON_IND = "",
    check.names = FALSE
  )
  w <- cnscwb(cnsprof(x), groups = c("Indigenous identity",
                                     "Non-Indigenous"), minpop = 0)
  p <- drw(plot(w))
  expect_identical(p$est[p$group == "income"], w$income)
  expect_identical(unique(p$label),
                   paste("Lakeside", w$group, sep = ": "))
})
