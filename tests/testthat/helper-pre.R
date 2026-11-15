# A simulated table in the layout of the published OTIS placement tables
# (a01 and b01: one row per placement), made by arithmetic on the row
# number, no random draws. Twelve people, two fiscal years each, two
# placements a year; one missing gender and one missing suicide-watch
# flag. dev/ref/preset.R sources this file and runs the rmorie pair
# builders on it.
opre <- local({
  i <- 1:48
  p <- (i - 1L) %/% 4L + 1L
  reg <- c("Central", "Eastern", "Northern", "Toronto", "Western")
  ra <- reg[1L + i %% 5L]
  d <- data.frame(
    UniqueIndividual_ID = sprintf("S%02d", p),
    EndFiscalYear = 2021L + ((i - 1L) %/% 2L) %% 2L,
    Gender = c("Female", "Male")[1L + p %% 2L],
    Age_Category = c("18 to 24", "25 to 49", "50+")[1L + p %% 3L],
    Region_AtTimeOfPlacement = ra,
    Region_MostRecentPlacement = ifelse(i %% 3L == 0L, reg[1L + (i + 1L) %% 5L],
                                        ra),
    MentalHealth_Alert = c("Yes", "No")[1L + i %% 2L],
    SuicideRisk_Alert = c("Yes", "No", "No")[1L + i %% 3L],
    SuicideWatch_Alert = c("No", "No", "No", "Yes", "y")[1L + i %% 5L],
    Number_Of_Placements = i %% 4L,
    NumberConsecutiveDays_Segregation = (i * 37) %% 23 + i / 8,
    stringsAsFactors = FALSE
  )
  d$Gender[5] <- NA
  d$SuicideWatch_Alert[9] <- NA
  d
})
