test_that("main-factor SIMPER matches vegan::simper() with the same seed", {
  dd <- make_comm()
  sim <- nested_simper(dd$comm, dd$main, dd$nested, permutations = 99,
                       cutoff = 1, seed = 5, verbose = FALSE)
  set.seed(5)
  ref <- vegan::simper(dd$comm, dd$main, permutations = 99)[[1]]
  o <- ref$ord
  expect_identical(sim$main$species, ref$species[o])
  expect_equal(sim$main$contribution, unname(ref$average[o]))
  expect_equal(sim$main$p, unname(ref$p[o]))
  expect_equal(sim$main$mean_a, unname(ref$ava[o]))
  expect_equal(sum(sim$main$pct), 100)
  expect_equal(sim$main$cum_pct, cumsum(sim$main$pct))
  expect_s3_class(sim$simper_objects$main, "simper")
})

test_that("nested SIMPER matches vegan::simper() on the subset", {
  dd <- make_comm()
  sim <- nested_simper(dd$comm, dd$main, dd$nested, permutations = 0,
                       cutoff = 1, verbose = FALSE)
  idx <- dd$main == "B"
  ref <- vegan::simper(dd$comm[idx, ], dd$nested[idx], permutations = 0)[[1]]
  b <- sim$nested[sim$nested$main_level == "B", ]
  expect_identical(b$species, ref$species[ref$ord])
  expect_equal(b$contribution, unname(ref$average[ref$ord]))
  expect_true(all(is.na(sim$nested$p)))
  expect_identical(names(sim$simper_objects$nested), c("A", "B"))
})

test_that("cutoff keeps taxa up to and including the one that crosses it", {
  dd <- make_comm()
  sim <- nested_simper(dd$comm, dd$main, dd$nested, permutations = 0,
                       cutoff = 0.5, verbose = FALSE)
  cp <- sim$main$cum_pct
  expect_gte(cp[length(cp)], 50)
  if (length(cp) > 1) expect_lt(cp[length(cp) - 1], 50)
  all_sp <- nested_simper(dd$comm, dd$main, dd$nested, permutations = 0,
                          cutoff = 1, verbose = FALSE)
  expect_identical(nrow(all_sp$main), ncol(dd$comm))
  expect_error(nested_simper(dd$comm, dd$main, dd$nested, cutoff = 0,
                             verbose = FALSE), "cutoff")
  expect_error(nested_simper(dd$comm, dd$main, dd$nested, cutoff = 1.2,
                             verbose = FALSE), "cutoff")
})

test_that("group names containing underscores are recovered correctly", {
  dd <- make_comm()
  nested <- ifelse(dd$nested == "x", "T_aphyllus", "T_verticillatus")
  main <- ifelse(dd$main == "A", "year_1", "year_2")
  sim <- nested_simper(dd$comm, main, nested, permutations = 0,
                       verbose = FALSE)
  expect_identical(unique(sim$main$group_a), "year_1")
  expect_identical(unique(sim$main$group_b), "year_2")
  expect_identical(unique(sim$nested$comparison),
                   "T_aphyllus vs T_verticillatus")
})

test_that("more than two nested groups give all pairwise comparisons", {
  dd <- make_comm()
  nested <- rep(c("x", "y", "z", "z", "x"), 4)
  sim <- nested_simper(dd$comm, dd$main, nested, permutations = 0,
                       cutoff = 1, verbose = FALSE)
  a <- sim$nested[sim$nested$main_level == "A", ]
  expect_setequal(unique(a$comparison), c("x vs y", "x vs z", "y vs z"))
  # group means must belong to the right group
  ref <- colMeans(dd$comm[dd$main == "A" & nested == "z", ])
  yz <- a[a$comparison == "y vs z", ]
  expect_equal(yz$mean_b, unname(ref[yz$species]))
})

test_that("skipped levels and empty rows are handled like nested_anosim()", {
  dd <- make_comm()
  main <- c(rep("A", 10), rep("B", 6), rep("C", 4))
  nested <- c(rep(c("x", "y"), each = 5), rep("x", 6), "x", "y", "y", "y")
  dd$comm[12, ] <- 0
  expect_warning(
    sim <- nested_simper(dd$comm, main, nested, permutations = 0,
                         verbose = FALSE),
    "1 empty sample\\(s\\) removed \\(rows: 12\\)"
  )
  expect_identical(sim$levels$note,
                   c("", "skipped: fewer than 2 nested groups",
                     "skipped: a nested group has < 2 samples"))
  expect_identical(unique(sim$nested$main_level), "A")
  expect_identical(sim$removed, 12L)
})

test_that("no testable nested level gives an empty nested table", {
  dd <- make_comm()
  sim <- nested_simper(dd$comm, dd$main, rep("only", 20), permutations = 0,
                       verbose = FALSE)
  expect_identical(nrow(sim$nested), 0L)
  expect_true("main_level" %in% names(sim$nested))
  expect_output(print(sim), "fewer than 2 nested groups")
})

test_that("print method works and returns the object invisibly", {
  dd <- make_comm()
  sim <- nested_simper(dd$comm, dd$main, dd$nested, permutations = 9,
                       verbose = FALSE)
  expect_output(out <- print(sim), "NESTED \\(TWO-WAY\\) SIMPER")
  expect_identical(out, sim)
  expect_output(print(sim), "A vs B")
})

test_that("strata restricts SIMPER permutations for the main factor only", {
  dd <- make_comm(seed = 9)
  cam <- rep(paste0("c", 1:10), 2)
  sim <- suppressMessages(nested_simper(dd$comm, dd$main, dd$nested, permutations = 99,
                       cutoff = 1, seed = 4, strata = cam, verbose = FALSE))
  set.seed(4)
  ref <- suppressMessages(vegan::simper(dd$comm, dd$main,
                       permutations = permute::how(nperm = 99,
                                                   blocks = factor(cam)))[[1]])
  expect_equal(sim$main$p, unname(ref$p[ref$ord]))
  free <- nested_simper(dd$comm, dd$main, dd$nested, permutations = 99,
                        cutoff = 1, seed = 4, verbose = FALSE)
  expect_identical(sim$nested, free$nested)
  expect_output(print(sim), "permutations within strata")
  # without permutations, strata is harmless
  expect_silent(nested_simper(dd$comm, dd$main, dd$nested, permutations = 0,
                              strata = cam, verbose = FALSE))
})

test_that("the report flags p-values computed with more than two groups", {
  dd <- make_comm()
  # main has 2 groups, nested has 2 per level: no note
  sim2 <- nested_simper(dd$comm, dd$main, dd$nested, permutations = 9,
                        verbose = FALSE)
  expect_false(any(grepl("more than two groups", capture.output(print(sim2)))))
  # 3 nested groups within each main level: note
  nested3 <- rep(c("x", "y", "z", "z", "x"), 4)
  sim3 <- nested_simper(dd$comm, dd$main, nested3, permutations = 9,
                        verbose = FALSE)
  expect_output(print(sim3), "more than two groups")
  # 3 main groups: note
  main3 <- rep(c("A", "B", "C", "A", "B"), each = 4)
  sim4 <- nested_simper(dd$comm, main3, rep("x", 20), permutations = 9,
                        verbose = FALSE)
  expect_output(print(sim4), "more than two groups")
  # no permutations, no p-values: no note
  sim5 <- nested_simper(dd$comm, dd$main, nested3, permutations = 0,
                        verbose = FALSE)
  expect_false(any(grepl("more than two groups", capture.output(print(sim5)))))
})
