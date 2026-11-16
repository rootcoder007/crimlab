test_that("palettes and session options", {
  expect_length(blpal(), 8L)
  expect_length(blpal("okabe", 11), 11L)
  expect_identical(blpal("grey", 2), c("#252525", "#636363"))
  expect_length(blpal("sequential", 7), 7L)
  expect_identical(blpal("diverging", 1), "#2166AC")
  expect_error(blpal("okabe", 0), "^blpal: `n` must be one positive")
  expect_error(blpal("neon"), "^blpal: the palette must be one of")
  expect_identical(yoypals(), c("diverging", "safe", "mono"))
  for (k in yoypals()) {
    expect_named(.YOY_PALETTES[[k]], c("bad", "good", "flat", "badbg",
                                        "goodbg"))
    expect_true(all(vapply(.YOY_PALETTES[[k]], .iscol, NA)))
  }
  old <- options(crimlab.plot = NULL)
  on.exit(options(old), add = TRUE)
  expect_identical(plotopt()$palette, "okabe")
  expect_invisible(plotopt(palette = "grey", grid = FALSE))
  expect_identical(plotopt()$palette, "grey")
  expect_false(plotopt()$grid)
  expect_error(plotopt(size = 3), "^plotopt: the options are")
  expect_error(plotopt(3), "^plotopt: the options are")
  expect_error(plotopt(palette = "neon"), "^plotopt: the palette must be")
  expect_error(plotopt(grid = NA), "`grid` must be TRUE or FALSE")
  expect_error(plotopt(cex = -1), "`cex` must be one positive number")
  expect_error(plotopt(las = 7), "`las` must be 0, 1, 2 or 3")
  expect_error(plotopt(bg = "notacolour"), "`bg` must be one colour")
  expect_error(plotopt(panel = 1:2), "`panel` must be one colour, or NA")
  expect_invisible(plotopt(panel = NA, bg = "#FFFFFF", las = 0, cex = 1.2))
  expect_true(is.na(plotopt()$panel))
  # set by hand: an unusable entry keeps its default
  options(crimlab.plot = list(grid = "x", cex = 2, palette = "neon"))
  expect_true(plotopt()$grid)
  expect_identical(plotopt()$cex, 2)
  expect_identical(plotopt()$palette, "okabe")
  # the grey grid on a plain panel, no grid at all, and the frame helpers
  options(crimlab.plot = list(panel = NA))
  expect_s3_class(drw(plotcdfs(c(1, 2, 3))), "data.frame")
  options(crimlab.plot = list(grid = FALSE))
  expect_s3_class(drw(plotcdfs(c(1, 2, 3))), "data.frame")
  expect_identical(.prange(NA), c(0, 1))
  expect_identical(.prange(2), c(1.5, 2.5))
  expect_identical(.pcols(3, col = "red"), rep("red", 3))
  expect_identical(.pcols(2, pal = "grey"), blpal("grey", 2))
  expect_error(.pcols(2, col = "nocolour"), "^plot: `col` must be colours")
  expect_error(.pmain(1, "x"), "^plot: `main` must be one string")
  expect_identical(.pmain("t", "x"), "t")
  m <- drw(.pmat(matrix(NA_real_, 2, 2), "empty", pal = "grey"))
  expect_identical(nrow(m), 4L)
  expect_identical(drw(.pmat(matrix(1, 1, 1), "one", div = TRUE))$value, 1)
})

test_that("the par settings come back after every plot", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(1, 1, 1, 1), cex = 1.3)
  plotmult(1:8, rep(1:4, 2), rep(c("a", "b"), each = 4))
  plot(didfe(did()$rate, did()$town, did()$year, did()$bylaw, bacon = TRUE))
  expect_identical(graphics::par("mar"), c(1, 1, 1, 1))
  expect_identical(graphics::par("cex"), 1.3)
  expect_identical(graphics::par("mfrow"), c(1L, 1L))
})

test_that("a plot of several panels draws them on one page", {
  d <- did()
  one <- d[d$town %in% c("A", LETTERS[9:12]), ]
  expect_identical(pages(plot(didfe(d$rate, d$town, d$year, d$bylaw,
                                    bacon = TRUE, loo = TRUE))), 1L)
  expect_identical(pages(plot(synthc(one$rate, one$town, one$year,
                                     one$first))), 1L)
  expect_identical(pages(plot(stkflw(people = c(5, 4, 3), days = c(9, 8, 9),
                                     prd = 1:3))), 1L)
  e <- evtst(d$rate, d$town, d$year, d$first, leads = 2, lags = 2)
  expect_identical(pages(plotevgrid(list(a = e, b = e, c = e))), 1L)
  expect_identical(pages(plotmult(1:9, rep(1:3, 3), rep(1:3, each = 3))),
                   1L)
})

test_that("blpal() refuses more colours than its bound", {
  expect_error(blpal(n = 1e9), "blpal: `n` must be at most 67108864.",
               fixed = TRUE)
})
