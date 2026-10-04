test_that("main-factor test matches vegan::anosim() exactly with the same seed", {
  dd <- make_comm()
  res <- nested_anosim(dd$comm, dd$main, dd$nested, permutations = 199,
                       seed = 7, verbose = FALSE)
  set.seed(7)
  ref <- vegan::anosim(vegan::vegdist(dd$comm, "bray"), dd$main,
                       permutations = 199)
  expect_identical(res$main$statistic, ref$statistic)
  expect_identical(res$main$signif, ref$signif)
  expect_identical(res$main$perm, ref$perm)
})

test_that("with a single nested level the overall test equals the main test", {
  dd <- make_comm()
  res <- nested_anosim(dd$comm, dd$main, rep("only", 20), permutations = 199,
                       seed = 7, verbose = FALSE)
  set.seed(7)
  ref <- vegan::anosim(vegan::vegdist(dd$comm), dd$main, permutations = 199)
  expect_identical(res$overall$statistic, ref$statistic)
  expect_identical(res$overall$perm, ref$perm)
  expect_true(all(is.na(res$nested$R)))
  expect_true(all(res$nested$note == "skipped: fewer than 2 nested groups"))
})

test_that("nested tests match vegan::anosim() on the subset", {
  dd <- make_comm()
  res <- nested_anosim(dd$comm, dd$main, dd$nested, permutations = 199,
                       seed = 3, verbose = FALSE)
  idx <- dd$main == "B"
  set.seed(3)
  ref <- vegan::anosim(vegan::vegdist(dd$comm[idx, ]), dd$nested[idx],
                       permutations = 199)
  expect_equal(res$nested$R[res$nested$main_level == "B"], unname(ref$statistic))
  expect_equal(res$nested$p[res$nested$main_level == "B"], ref$signif)
  expect_s3_class(res$nested_objects[["B"]], "anosim")
})

test_that("empty rows are removed with a warning", {
  dd <- make_comm()
  dd$comm[c(3, 15), ] <- 0
  expect_warning(
    res <- nested_anosim(dd$comm, dd$main, dd$nested, permutations = 99,
                         verbose = FALSE),
    "2 empty sample\\(s\\) removed \\(rows: 3, 15\\)"
  )
  expect_identical(res$removed, c(3L, 15L))
  expect_identical(res$kept, setdiff(1:20, c(3L, 15L)))
  expect_identical(attr(res$dist, "Size"), 18L)
  expect_equal(sum(res$group_summary$N_samples), 18)
  expect_error(
    nested_anosim(dd$comm, dd$main, dd$nested, remove_empty = FALSE,
                  verbose = FALSE),
    "empty sample"
  )
})

test_that("nested levels with too few groups or samples are skipped", {
  dd <- make_comm()
  main <- c(rep("A", 10), rep("B", 6), rep("C", 4))
  nested <- c(rep(c("x", "y"), each = 5),   # A: testable
              rep("x", 6),                    # B: one group only
              "x", "y", "y", "y")             # C: group x has 1 sample
  res <- nested_anosim(dd$comm, main, nested, permutations = 99,
                       verbose = FALSE)
  tab <- res$nested
  expect_identical(tab$main_level, c("A", "B", "C"))
  expect_false(is.na(tab$R[1]))
  expect_identical(tab$note[2], "skipped: fewer than 2 nested groups")
  expect_identical(tab$note[3], "skipped: a nested group has < 2 samples")
  expect_true(all(is.na(tab[2:3, c("R", "p", "p_adj")])))
  expect_identical(names(res$nested_objects), "A")
})

test_that("p_adjust is applied to the nested p-values only", {
  dd <- make_comm()
  res <- nested_anosim(dd$comm, dd$main, dd$nested, permutations = 99,
                       p_adjust = "holm", seed = 1, verbose = FALSE)
  expect_equal(res$nested$p_adj, stats::p.adjust(res$nested$p, "holm"))
  res0 <- nested_anosim(dd$comm, dd$main, dd$nested, permutations = 99,
                        seed = 1, verbose = FALSE)
  expect_equal(res0$nested$p_adj, res0$nested$p)
  expect_error(
    nested_anosim(dd$comm, dd$main, dd$nested, p_adjust = "nonsense",
                  verbose = FALSE)
  )
})

test_that("input checks work", {
  dd <- make_comm()
  bad <- dd$comm; bad$site <- "a"
  expect_error(nested_anosim(bad, dd$main, dd$nested, verbose = FALSE),
               "only numeric")
  bad <- dd$comm; bad[1, 1] <- NA
  expect_error(nested_anosim(bad, dd$main, dd$nested, verbose = FALSE), "NA")
  expect_error(nested_anosim(dd$comm, dd$main[-1], dd$nested, verbose = FALSE),
               "same length")
  expect_error(nested_anosim(dd$comm, rep("A", 20), dd$nested, verbose = FALSE),
               "at least two levels")
})

test_that("print method works and returns the object invisibly", {
  dd <- make_comm()
  res <- nested_anosim(dd$comm, dd$main, dd$nested, permutations = 99,
                       p_adjust = "BH", verbose = FALSE)
  expect_output(out <- print(res), "NESTED \\(TWO-WAY\\) ANOSIM")
  expect_identical(out, res)
  expect_output(print(res), "p_adj")
})

test_that("example dataset has the documented structure", {
  expect_identical(dim(mistletoe_visitors), c(87L, 43L))
  expect_identical(names(mistletoe_visitors)[1:3],
                   c("Year_study", "Mistletoe", "Camera"))
  expect_false(any(c("Mimus_thenca", "Sephanoides_sephaniodes",
                     "not_identified") %in% names(mistletoe_visitors)))
  expect_equal(sum(rowSums(mistletoe_visitors[, -(1:3)]) == 0), 7)
})
