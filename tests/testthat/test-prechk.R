# prechk(): the columns a preset reads (directly, in filters and in
# recipes), the publisher's codes, mappings without a default, refusals.
# The tables are simulated.

ns <- asNamespace("crimlab")

b03 <- data.frame(
  EndFiscalYear = rep(2021:2022, each = 4),
  Region_AtTimeOfPlacement = rep(c("Central", "Toronto"), 4),
  Institution_AtTimeOfPlacement = rep(c("Quinte Detention Centre",
                                        "Toronto South Detention Centre"), 4),
  Alert_Type = "Mental Health Alert",
  Alert_Presence = rep(c("Yes", "Yes", "No", "No"), 2),
  Number_SegregationPlacements = c(31, 52, 140, 260, 28, 61, 133, 248),
  stringsAsFactors = FALSE
)

test_that("a table that fits passes unchanged", {
  expect_identical(prechk(b03, "otis-b03-cnt"), b03)
  expect_invisible(prechk(b03, "otis-b03-cnt"))
  # case, surrounding space and missing values are not unknown codes
  d <- b03
  d$Region_AtTimeOfPlacement[1:2] <- c(" toronto", NA)
  expect_identical(prechk(d, "otis-b03-cnt"), d)
  # filters read their column: b04 keeps the medians
  b04 <- data.frame(EndFiscalYear = 2022, Region_AtTimeOfPlacement = "Eastern",
                    Gender = c("Female", "All"), Measure = "Median",
                    NumberConsecutiveDays_Segregation = c(3, 4),
                    stringsAsFactors = FALSE)
  expect_identical(prechk(b04, "otis-b04-cnt"), b04)
  expect_error(prechk(b04[-4], "otis-b04-cnt"), "has no column Measure")
})

test_that("missing columns and unknown codes are named together", {
  d <- b03[-6]
  d$Region_AtTimeOfPlacement[c(1, 3)] <- c("Southern", "Southern")
  d$Alert_Type[2] <- "Flag X"
  e <- expect_error(prechk(d, "otis-b03-cnt"), "does not fit preset otis-b03")
  m <- conditionMessage(e)
  expect_match(m, "has no column Number_SegregationPlacements", fixed = TRUE)
  expect_match(m, paste("column Region_AtTimeOfPlacement has labels not",
                        "among its codes: Southern."), fixed = TRUE)
  expect_match(m, "column Alert_Type has labels not among its codes: Flag X.",
               fixed = TRUE)
})

test_that("recipes: the columns read, not the columns made", {
  # acvm reads the placement columns and makes vm and hiac; opre codes
  # some suicide-watch flags y, which clprep() reads as yes but the
  # publisher does not list
  expect_error(prechk(opre, "otis-a01-cbatt"),
               "column SuicideWatch_Alert has labels not among its codes: y.",
               fixed = TRUE)
  o <- opre
  o$SuicideWatch_Alert[which(o$SuicideWatch_Alert == "y")] <- "Yes"
  expect_identical(prechk(o, "otis-a01-cbatt"), o)
  expect_identical(prechk(opre, "otis-pc-cbatt"), opre)
  expect_error(prechk(o[-9], "otis-a01-cbatt"),
               "has no column SuicideWatch_Alert")
  expect_error(prechk(o[-10], "otis-pb-cbatt"),
               "has no column Number_Of_Placements")
  # long stacks the two count columns it names
  b07 <- data.frame(EndFiscalYear = 2022, Alert_Type = "Transgender Alert",
                    Gender = "Male",
                    Number_Segregation_Placements_With_Alert = 4,
                    Number_Segregation_Placements_Without_Alert = 9,
                    stringsAsFactors = FALSE)
  expect_identical(prechk(b07, "otis-b07-cnt"), b07)
  expect_error(prechk(b07[-5], "otis-b07-cnt"),
               "no column Number_Segregation_Placements_Without_Alert")
  expect_identical(ns$.rin("gent", "a,b,c,d,0.5,0.6"), c("a", "b", "c", "d"))
  expect_identical(ns$.rin("pop", "n;r,x;1000"), c("n", "r,x"))
  expect_identical(ns$.rin("pairc", ""), c(
    "Region_AtTimeOfPlacement", "Region_MostRecentPlacement",
    "NumberConsecutiveDays_Segregation"
  ))
})

test_that("a preset of one's own, with codes of its own", {
  # a made-up agency's layout, simulated: counts by sex and facility
  mine <- data.frame(
    preset = "dmo-cnt", generic = "cntreg", arg = c("y", "grp", "adj"),
    column = c("n", "sex", "unit"),
    recode = c("", "f|female=1;m|male=0", "fac"), stringsAsFactors = FALSE
  )
  d <- data.frame(n = c(4, 30, 6, 41), sex = c("F", "m", "f", "Male"),
                  unit = c("North", "North", "South", "South"),
                  stringsAsFactors = FALSE)
  expect_identical(prechk(d, mine), d)
  d$sex[2] <- "x"
  expect_error(prechk(d, mine), "column sex has labels with no recode: x.",
               fixed = TRUE)
  # a has: rule covers labels containing its text
  mine$recode[2] <- "has:fem=1;m|male=0"
  d$sex <- c("Female", "m", "female", "male")
  expect_identical(prechk(d, mine), d)
  # the agency's codes, given as a data frame
  cd <- data.frame(column = "unit", value = c("North", "South"),
                   stringsAsFactors = FALSE)
  expect_identical(prechk(d, mine, cd), d)
  d$unit[1] <- "East"
  expect_error(prechk(d, mine, cd), "not among its codes: East.", fixed = TRUE)
  # a constant and a mapping with a default are not checked for labels
  mine2 <- data.frame(preset = "dmo-c", generic = "cntreg",
                      arg = c("y", "grp", "gee"), column = c("n", "sex", ""),
                      recode = c("", "f=1;*=0", "const:TRUE"),
                      stringsAsFactors = FALSE)
  expect_identical(prechk(d, mine2), d)
})

test_that("the shipped codes", {
  k <- ns$.presets("codes.tsv")
  expect_named(k, c("agency", "column", "value", "source"))
  expect_false(anyDuplicated(k[c("agency", "column", "value")]) > 0)
  expect_true(all(nzchar(k$source)))
  o <- k$column[k$agency == "otis"]
  expect_true(all(c("Gender", "Region_AtTimeOfPlacement", "Alert_Type",
                    "Race", "Religion") %in% o))
  expect_setequal(k$value[k$column == "Gender"], c("Female", "Male", "All"))
})

test_that("refusals", {
  errs(
    "prechk",
    list(
      list("`preset` must be one name", b03, c("a", "b")),
      list("`preset` must be one name", b03, NA_character_),
      list("no preset named zz", b03, "zz"),
      list("`data` must be a data frame", 1:3, "otis-b03-cnt"),
      list("`codes` must be a data frame", b03, "otis-b03-cnt", list(1)),
      list("`codes` must be a data frame", b03, "otis-b03-cnt",
           data.frame(a = 1)),
      list("unknown recipe zz", b03,
           data.frame(preset = "x", generic = "ctabx", arg = "recipe",
                      column = "zz", stringsAsFactors = FALSE))
    )
  )
})

test_that("prechk() refuses a bad preset table in its own name", {
  expect_error(prechk(data.frame(a = 1), data.frame()), "^prechk: ")
})

test_that("text that is not UTF-8 is refused in prechk's words", {
  skip_if_not(isTRUE(l10n_info()[["UTF-8"]]), "not a UTF-8 session")
  expect_error(prechk(data.frame(a = "\xff"), "otis-b03-cnt"),
               "^prechk: `data` holds text that is not valid")
})
