test_that("ovrtxt() keeps blocks on lines of their own", {
  h <- paste0("<html><head><title>R</title><style>p {color: red}</style>",
              "</head><body><h2>The Team</h2>\n<p>Number of SIU ",
              "Investigators\nassigned: 3</p><!-- a <p> note -->",
              "<p>WO\n#1 <b>Interviewed</b></p><script>var a = '<p>';",
              "</script><p>WO&nbsp;#2</p></body></html>")
  # written out: the title, the heading and each paragraph on a line,
  # one empty line between blocks (a block ends a line where it starts
  # and where it ends), the page's own breaks spaces
  expect_identical(
    ovrtxt(h),
    paste("R", "", "The Team", "", "Number of SIU Investigators assigned: 3",
          "", "WO #1 Interviewed", "", "WO #2", sep = "\n")
  )
  expect_identical(ovrtxt("a<br/>b<BR>c"), "a\nb\nc")
})

test_that("ovrtxt() decodes entities and leaves plain text its lines", {
  expect_identical(
    ovrtxt("<p>&eacute;t&#233; &amp; &#x2019;s &bogus; &lt;b&gt;</p>"),
    "\u00e9t\u00e9 & 's &bogus; <b>"
  )
  expect_identical(ovrtxt("a  b\r\nc\n\n\n\nd ", html = FALSE),
                   "a b\nc\n\nd")
  # no markup: a '<' that opens no tag
  expect_identical(ovrtxt("1 < 2\nand 3 > 2"), "1 < 2\nand 3 > 2")
  expect_identical(ovrtxt("<p>a</p>", html = FALSE), "<p>a</p>")
  expect_identical(ovrtxt(c(NA, "x")), c(NA, "x"))
  expect_identical(ovrtxt(NA), NA_character_)
  expect_identical(ovrtxt(character(0)), character(0))
})

test_that("ovrtxt() refuses what is not text", {
  expect_error(ovrtxt(1),
               "ovrtxt: `x` must be text (a character vector), not numeric.",
               fixed = TRUE)
  expect_error(ovrtxt(list("a")),
               "ovrtxt: `x` must be text (a character vector), not list.",
               fixed = TRUE)
  expect_error(ovrtxt("a", html = "yes"),
               "ovrtxt: `html` must be NA, TRUE or FALSE.", fixed = TRUE)
  expect_error(ovrtxt("a", html = c(TRUE, FALSE)),
               "ovrtxt: `html` must be NA, TRUE or FALSE.", fixed = TRUE)
})

test_that("ovrdte() reads dates written out in English and French", {
  x <- c("January 5, 2023", "Jan. 5th 2023", "5 January 2023",
         "5th of January, 2023", "5 janvier 2023", "2023-01-05",
         "le 1er février 2021", "3 AOÛT 2017", "SEPT. 9, 2019")
  expect_identical(ovrdte(x),
                   c(rep("2023-01-05", 6), "2021-02-01", "2017-08-03",
                     "2019-09-09"))
  # a day that does not exist is passed over; the next date is read
  expect_identical(ovrdte("February 30, 2023, then March 2, 2023"),
                   "2023-03-02")
  expect_identical(ovrdte("February 29, 2024"), "2024-02-29")
  expect_identical(ovrdte(c("May 2023", "no date", NA, "2023-13-01")),
                   rep(NA_character_, 4))
  expect_identical(ovrdte(character(0)), character(0))
})

test_that("ovrdte() refuses what is not text", {
  expect_error(ovrdte(20230105),
               "ovrdte: `x` must be text (a character vector), not numeric.",
               fixed = TRUE)
})

test_that("ovrsch() gives the fields of the preset and their rules", {
  s <- ovrsch()
  expect_named(s, c("name", "rule", "cue", "alt", "sect", "desc"))
  expect_identical(nrow(s), 17L)
  expect_true(all(s$rule %in% c("after", "num", "date", "tag", "sent",
                                "prev", "lines", "name", "flag", "before",
                                "word")))
  expect_identical(s$cue[s$name == "number_of_witness_officials"], "WO")
  expect_identical(anyDuplicated(s$name), 0L)
  expect_identical(ovrsch("siu"), s)
})

test_that("ovrsch() refuses a preset it does not have", {
  expect_error(ovrsch("iio"),
               "ovrsch: `preset` must be one of \"siu\", not \"iio\".",
               fixed = TRUE)
})

# A director's report laid out as the unit publishes them (no real case):
# headings on lines of their own, tags split over lines, a signature.
ovr_rep <- paste0(
  "<html><body><h1>SIU Director's Report - Case # 23-OFD-001</h1>",
  "<h2>The Investigation</h2><h3>Notification of the SIU</h3>",
  "<p>On January 6, 2023, at 9 a.m. the Barrie Police Service contacted ",
  "the SIU with the\nfollowing information. On January 5, 2023, in the ",
  "City of Barrie, a 34-year-old man was injured during an arrest.</p>",
  "<h3>The Team</h3><p>Number of SIU Investigators assigned: 3</p>",
  "<p>Number of SIU Forensic Investigators assigned: one</p>",
  "<h3>Civilian Witnesses</h3><p>CW\n#1 Interviewed</p>",
  "<p>CW #2 Interviewed</p><h3>Subject Officials</h3>",
  "<p>SO #1 Declined interview</p><h3>Witness Officials</h3>",
  "<p>WO #1 Interviewed</p><p>WO #2 Interviewed</p><p>WO #3 Interviewed",
  "</p><h2>Incident Narrative</h2><p>On January 5, 2023, officers of the ",
  "Barrie Police Service responded to a disturbance.</p>",
  "<h2>Nature of Injuries / Treatment</h2><p>The man was taken to ",
  "hospital. He had a fractured left arm. He was released.</p>",
  "<h2>Relevant Legislation</h2>",
  "<p>Section 25, Criminal Code - Protection of persons</p>",
  "<p>Section 34, Criminal Code - Defence of person</p>",
  "<p>Section 25, Criminal Code - Protection of persons</p>",
  "<h2>Analysis and Director's Decision</h2><p>There are no reasonable ",
  "grounds to believe that the subject official committed an offence.</p>",
  "<p>Date: April 28, 2023</p><p>Joseph Martino</p><p>Director</p>",
  "<p>Special Investigations Unit</p></body></html>"
)

test_that("ovrpar() reads every field of the preset from a report", {
  p <- ovrpar(ovr_rep)
  expect_named(p, ovrsch()$name)
  # written out from the report above
  want <- list(
    case_number = "23-OFD-001",
    police_service = "Barrie Police Service",
    date_of_incident_iso = "2023-01-05",
    date_siu_notified_iso = "2023-01-06",
    date_of_director_decision_iso = "2023-04-28",
    siu_investigators = 3L,
    siu_forensics_investigators = 1L,
    number_of_witness_officials = 3L,
    number_of_civilian_witnesses = 2L,
    number_of_subject_officials = 1L,
    age_affected = 34L,
    sex_gender_affected = "man",
    charges_recommended = FALSE,
    directors_name = "Joseph Martino",
    location_of_call = "Barrie",
    specific_injuries = "He had a fractured left arm.",
    relevant_legislation = paste(
      "Section 25, Criminal Code - Protection of persons;",
      "Section 34, Criminal Code - Defence of person"
    )
  )
  expect_identical(as.list(p), want)
})

test_that("ovrpar() leaves unstated fields missing, one row per report", {
  p <- ovrpar(c(ovr_rep, "<p>not a report</p>", NA))
  expect_identical(nrow(p), 3L)
  expect_identical(p$siu_investigators, c(3L, NA, NA))
  expect_identical(p$charges_recommended, c(FALSE, NA, NA))
  expect_identical(p$directors_name, c("Joseph Martino", NA, NA))
  # a schema of one's own: another agency's wording
  s <- data.frame(name = c("officers", "file", "opened", "found"),
                  rule = c("tag", "after", "date", "flag"),
                  cue = c("Officer", "File:", "opened", "substantiated"),
                  alt = c("", "", "", "not substantiated"),
                  sect = c("", "", "", "Findings"))
  r <- paste("File: 2024-117", "Opened on 3 March 2024 after a complaint.",
             "Officer #1 and Officer #2 were interviewed.", "Findings",
             "The complaint was not substantiated.", sep = "\n")
  expect_identical(
    as.list(ovrpar(r, s)),
    list(officers = 2L, file = "2024-117", opened = NA_character_,
         found = FALSE)
  )
  s$cue[3] <- "Opened"
  expect_identical(ovrpar(r, s)$opened, "2024-03-03")
})

test_that("ovrpar() refuses a schema it cannot read", {
  s <- ovrsch()
  expect_error(ovrpar("a", s[, 1:3]),
               paste("ovrpar: `schema` must be a data frame with the columns",
                     "name, rule, cue, alt, sect."), fixed = TRUE)
  s2 <- s
  s2$rule[1] <- "regex"
  expect_error(ovrpar("a", s2),
               paste("ovrpar: `schema` has a rule that is not one of after,",
                     "num, date, tag, sent, prev, lines, name, flag, before,",
                     "word: \"regex\"."), fixed = TRUE)
  s3 <- s
  s3$name[2] <- s3$name[1]
  expect_error(ovrpar("a", s3),
               "ovrpar: the names of `schema` must be given once each.",
               fixed = TRUE)
  s4 <- s
  s4$cue[1] <- NA
  expect_error(ovrpar("a", s4),
               paste("ovrpar: the columns of `schema` must be text with no",
                     "missing values."), fixed = TRUE)
  expect_error(ovrpar(3),
               "ovrpar: `x` must be text (a character vector), not numeric.",
               fixed = TRUE)
})

test_that("ovrsoc() counts subject officials by the first rule that applies", {
  x <- c("Subject Officials\nSO #1 Interviewed\nSO #2 Declined interview",
         "SO #1 and SO #3 were designated.",
         "SO #3 was designated.",
         "The SIU designated three subject officers.",
         "The subject official declined an interview.",
         "There was no subject official in this case.",
         "Witness Officials\nWO #1 Interviewed\nWO #2 Interviewed",
         "No officers are named here.", NA)
  r <- ovrsoc(x)
  expect_named(r, c("count", "rule", "reason"))
  # written out: the rule each text meets first, and its count
  expect_identical(r$count, c(2L, 3L, NA, 3L, 1L, 0L, 0L, NA, NA))
  expect_identical(r$rule, c(1L, 2L, 6L, 3L, 4L, 5L, 5L, 6L, NA))
  expect_identical(r$reason[c(1, 4, 6, 7, 8)], c(
    "the roster under the heading: highest SO #2",
    "a plural: 'three subject officers'",
    "none: 'no subject official'",
    "witness officials only: WO #2",
    "unresolved"
  ))
  # a body's own words
  expect_identical(
    ovrsoc("Agents impliqués\nAI no 1\nAI no 2", tag = "AI",
           term = "agent impliqué", wtag = "AT",
           head = "Agents impliqués")$count,
    2L
  )
})

test_that("ovrsoc() refuses words it cannot look for", {
  expect_error(ovrsoc("a", tag = c("SO", "AI")),
               "ovrsoc: `tag` must be one text, not empty, with no \"|\".",
               fixed = TRUE)
  expect_error(ovrsoc("a", term = "a|b"),
               "ovrsoc: `term` must be text, not empty, with no \"|\".",
               fixed = TRUE)
  expect_error(ovrsoc("a", head = NA_character_),
               "ovrsoc: `head` must be text, not empty, with no \"|\".",
               fixed = TRUE)
  expect_error(ovrsoc(1),
               "ovrsoc: `x` must be text (a character vector), not numeric.",
               fixed = TRUE)
})

test_that("ovrkwd() marks the groups whose words a report holds", {
  x <- c("The man had a history of schizophrenia; a crisis team came.",
         "<p>An <b>Indigenous</b> woman, in a mental health crisis.</p>",
         "Blackburn Road; the EDPs unit; nothing of either kind.",
         NA)
  k <- ovrkwd(x)
  expect_named(k, c("mental_health", "race", "terms"))
  # written out: schizophren* starts a word of the first; the second has
  # both groups; "Blackburn" and "EDPs" are not the words "Black", "EDP"
  expect_identical(k$mental_health, c(TRUE, TRUE, FALSE, NA))
  expect_identical(k$race, c(FALSE, TRUE, FALSE, NA))
  expect_identical(k$terms, c("schizophren*", "mental health; Indigenous",
                              "", NA))
  w <- ovrkwd(c("A Taser and a baton.", "none"),
              list(force = c("taser", "baton"), dog = "police dog"))
  expect_identical(w$force, c(TRUE, FALSE))
  expect_identical(w$dog, c(FALSE, FALSE))
  expect_identical(w$terms, c("taser; baton", ""))
  v <- ovrkwd("baton", c("taser", "baton"))
  expect_named(v, c("taser", "baton", "terms"))
})

test_that("ovrkwd() refuses keyword groups it cannot read", {
  m <- paste("ovrkwd: `vocab` must be \"mhr\", a list of words named by",
             "group, or words.")
  expect_error(ovrkwd("a", 3), m, fixed = TRUE)
  expect_error(ovrkwd("a", list("taser")), m, fixed = TRUE)
  expect_error(ovrkwd("a", list(a = NA_character_)), m, fixed = TRUE)
  expect_error(ovrkwd("a", ""), m, fixed = TRUE)
})

test_that("ovrchk() flags cells not of the form of their field", {
  d <- data.frame(
    case_number = c("23-OFD-001", "23-OFD-1", NA, "24-TCI-140"),
    date_of_incident_iso = c("2023-01-05", "2023-02-30", "", "5 Jan 2023"),
    number_of_subject_officials = c(1, 2.5, NA, 100000),
    charges_recommended = c("FALSE", "maybe", "yes", NA),
    other = c("x", "y", "z", "w")
  )
  r <- ovrchk(d, shape = c(case_number = "99-AAA-999|99-AA-999"))
  expect_named(r, c("row", "issues", "what"))
  # written out: row 2 fails the shape, the date (February 30), the count
  # (2.5) and the flag; row 4 the date only; the missing cells pass
  expect_identical(r$issues, c(0L, 4L, 0L, 1L))
  expect_identical(r$what[2], paste(
    "date_of_incident_iso: date; number_of_subject_officials: count;",
    "charges_recommended: flag; case_number: shape"
  ))
  expect_identical(r$what[c(1, 4)], c("", "date_of_incident_iso: date"))
  # a table coded by ovrpar() passes its own check
  p <- ovrpar(c("<p>Date: May 5, 2023</p><p>WO #2</p>", "<p>none</p>"))
  expect_identical(ovrchk(p)$issues, c(0L, 0L))
})

test_that("ovrchk() refuses what it cannot check", {
  expect_error(ovrchk(list(a = 1)),
               "ovrchk: `data` must be a data frame, not list.", fixed = TRUE)
  expect_error(ovrchk(data.frame(a = 1), data.frame(x = 1)),
               paste("ovrchk: `schema` must be a data frame with the columns",
                     "name and rule."), fixed = TRUE)
  expect_error(ovrchk(data.frame(a = 1), shape = c(b = "999")),
               "ovrchk: `shape` must be patterns named by columns of `data`.",
               fixed = TRUE)
  expect_error(ovrchk(data.frame(a = 1), shape = "999"),
               "ovrchk: `shape` must be patterns named by columns of `data`.",
               fixed = TRUE)
})

test_that("keyword groups and words are refused in words", {
  errs("ovrkwd", list(
    list("a list of words named by group, or words", "a text",
         strrep("a", 1e5))
  ))
  errs("ovrsoc", list(list("`tag` holds text that is not valid", "a text",
                           tag = "\xff")))
})
