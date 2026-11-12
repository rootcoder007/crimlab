test_that("proof() selects rows by function name", {
  p <- proof()
  expect_named(p, c("fn", "theorem", "file"))
  expect_false(anyDuplicated(p[["theorem"]]) > 0)
  two <- proof(c("cllr", "dbsrch"))
  expect_identical(unique(two[["fn"]]), c("cllr", "dbsrch"))
  expect_identical(two, p[p[["fn"]] %in% c("cllr", "dbsrch"), ],
                   ignore_attr = "row.names")
  expect_identical(rownames(two), as.character(seq_len(nrow(two))))
  expect_true(all(startsWith(p[["theorem"]], "Crimlab.")))
  expect_error(proof(NA_character_), "function names")
  expect_error(proof(c("cllr", "nope", "zzz")), "nope, zzz")
})
