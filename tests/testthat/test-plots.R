test_that("plot_anosim_box() returns ggplot objects", {
  dd <- make_comm()
  res <- nested_anosim(dd$comm, dd$main, dd$nested, permutations = 99,
                       verbose = FALSE)
  for (w in c("overall", "main", "A", "B")) {
    expect_s3_class(plot_anosim_box(res, which = w), "ggplot")
  }
  p <- plot_anosim_box(res$main, palette = c("grey", "red", "blue"))
  expect_s3_class(p, "ggplot")
  expect_identical(levels(p$data$class)[1], "Between")
  expect_error(plot_anosim_box(res, which = "Z"), "'which' must be")
  expect_error(plot_anosim_box(list()), "'anosim' or a 'nested_anosim'")
})

test_that("plot_nmds_hulls() works with all input types", {
  dd <- make_comm()
  dd$comm[2, ] <- 0
  res <- suppressWarnings(
    nested_anosim(dd$comm, dd$main, dd$nested, permutations = 99,
                  verbose = FALSE)
  )
  grp <- paste(dd$main, dd$nested)

  out <- plot_nmds_hulls(res, groups = grp, trymax = 5, seed = 1)
  expect_named(out, c("plot", "mds", "scores", "stress"))
  expect_s3_class(out$plot, "ggplot")
  expect_s3_class(out$mds, "metaMDS")
  expect_identical(nrow(out$scores), 19L)
  expect_true(is.finite(out$stress))

  keep <- rowSums(dd$comm) > 0
  out2 <- plot_nmds_hulls(dd$comm[keep, ], groups = grp[keep], trymax = 5,
                          seed = 1)
  expect_s3_class(out2$plot, "ggplot")
  out3 <- plot_nmds_hulls(vegan::vegdist(dd$comm[keep, ]), groups = grp[keep],
                          trymax = 5, seed = 1, hull = FALSE)
  expect_s3_class(out3$plot, "ggplot")
  out4 <- plot_nmds_hulls(out$mds, groups = grp[keep])
  expect_identical(out4$stress, out$stress)

  expect_error(plot_nmds_hulls(out$mds, groups = grp[1:5]), "does not match")
})

test_that("ggsci palettes are resolved by name", {
  skip_if_not_installed("ggsci")
  cols <- .resolve_palette("startrek", 2)
  expect_length(cols, 2)
  expect_null(.resolve_palette(NULL, 3))
  expect_identical(.resolve_palette(c("red", "blue"), 3), c("red", "blue", "red"))
})
