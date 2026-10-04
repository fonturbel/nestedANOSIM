#' Two-way (hierarchical) SIMPER with tidy output
#'
#' Companion to [nested_anosim()]. Once an ANOSIM shows that groups differ,
#' SIMPER (similarity percentages; Clarke 1993) asks which taxa contribute most
#' to the average Bray-Curtis dissimilarity between them. This function wraps
#' [vegan::simper()] with the same design as [nested_anosim()]:
#'   1. SIMPER between the levels of the main factor.
#'   2. SIMPER between the levels of the nested factor, separately within each
#'      level of the main factor.
#'
#' Results are returned as tidy data frames (one row per taxon and pairwise
#' comparison), sorted by contribution and trimmed to the taxa that together
#' account for `cutoff` of the dissimilarity.
#'
#' SIMPER always uses Bray-Curtis dissimilarity. Its contributions are
#' descriptive: highly abundant and highly variable taxa tend to dominate them
#' (Warton et al. 2012), so read them alongside the group means (`mean_a`,
#' `mean_b`) and, if permutations are used, the p-values.
#'
#' @inheritParams nested_anosim
#' @param permutations Number of permutations for the taxon-level tests in
#'   [vegan::simper()] (default 999). Use 0 to skip them (no p-values).
#' @param cutoff Cumulative proportion of the dissimilarity to report (default
#'   0.7, i.e. taxa up to 70%). Taxa are added in decreasing order of
#'   contribution until this proportion is reached, including the taxon that
#'   crosses it. Use 1 to keep all taxa.
#' @param seed Optional integer. If given, the RNG seed is set to this value
#'   before each SIMPER, so every comparison is reproducible on its own.
#' @param verbose Logical; print the formatted report (default TRUE).
#'
#' @return An object of class `"nested_simper"`, a list with:
#'   \item{main}{tidy data.frame for the main-factor comparisons}
#'   \item{nested}{tidy data.frame for the nested-factor comparisons within
#'     each main-factor level (column `main_level`)}
#'   \item{levels}{data.frame with the group sizes of each main-factor level and
#'     a note when its nested comparison was skipped}
#'   \item{simper_objects}{list with the full [vegan::simper()] objects:
#'     `main`, and `nested` (one per main-factor level that was tested)}
#'   \item{removed}{indices of the samples removed (empty rows)}
#'   \item{kept}{indices of the samples retained (relative to the input)}
#'   \item{settings}{list with the main settings used}
#'
#'   The tidy tables have the columns `comparison`, `group_a`, `group_b`,
#'   `species`, `contribution` (average contribution to the Bray-Curtis
#'   dissimilarity), `sd`, `ratio` (`contribution / sd`), `mean_a`, `mean_b`
#'   (mean abundance in each group), `pct` (percentage of the average
#'   dissimilarity), `cum_pct`, `p` (`NA` without permutations) and
#'   `dissimilarity` (average Bray-Curtis dissimilarity between the two groups).
#'
#' @references
#' Clarke KR (1993) Non-parametric multivariate analyses of changes in community
#' structure. *Australian Journal of Ecology* 18: 117-143.
#'
#' Warton DI, Wright ST, Wang Y (2012) Distance-based multivariate analyses
#' confound location and dispersion effects. *Methods in Ecology and Evolution*
#' 3: 89-101.
#'
#' @examples
#' data(mistletoe_visitors)
#' comm <- mistletoe_visitors[, -(1:3)]
#' sim <- nested_simper(comm   = comm,
#'                      main   = mistletoe_visitors$Year_study,
#'                      nested = mistletoe_visitors$Mistletoe,
#'                      permutations = 99, seed = 123)
#' sim$main
#' subset(sim$nested, main_level == "First")
#'
#' @seealso [nested_anosim()], [vegan::simper()]
#' @importFrom vegan simper
#' @importFrom utils combn
#' @export
nested_simper <- function(comm, main, nested,
                          permutations = 999,
                          cutoff = 0.7,
                          remove_empty = TRUE,
                          seed = NULL,
                          verbose = TRUE) {

  if (!is.numeric(cutoff) || length(cutoff) != 1 || cutoff <= 0 || cutoff > 1) {
    stop("'cutoff' must be a single number in (0, 1].")
  }

  # ---- Input checks and removal of empty samples ---------------------------
  dat <- .prepare_input(comm, main, nested, remove_empty)
  comm <- dat$comm; main <- dat$main; nested <- dat$nested

  # Reset the seed before each SIMPER so each one is reproducible on its own
  run_simper <- function(x, grouping) {
    if (!is.null(seed)) set.seed(seed)
    vegan::simper(x, grouping, permutations = permutations)
  }

  # ---- 1. Main factor -------------------------------------------------------
  main_sim <- run_simper(comm, main)
  main_tab <- .tidy_simper(main_sim, main, cutoff)

  # ---- 2. Nested factor within each level of main --------------------------
  nested_sims <- list()
  nested_rows <- list()
  level_rows  <- list()
  for (lv in sort(unique(main))) {
    idx <- which(main == lv)
    sub_nested <- nested[idx]
    sizes <- table(sub_nested)
    note <- .skip_note(sizes)
    level_rows[[lv]] <- data.frame(
      main_level = lv, n_groups = length(sizes), n_samples = length(idx),
      group_sizes = paste(names(sizes), sizes, sep = "=", collapse = ", "),
      note = note, stringsAsFactors = FALSE)
    if (note == "") {
      s <- run_simper(comm[idx, , drop = FALSE], sub_nested)
      nested_sims[[lv]] <- s
      tab <- .tidy_simper(s, sub_nested, cutoff)
      nested_rows[[lv]] <- cbind(main_level = lv, tab, stringsAsFactors = FALSE)
    }
  }
  levels_tab <- do.call(rbind, level_rows)
  rownames(levels_tab) <- NULL
  nested_tab <- if (length(nested_rows)) do.call(rbind, nested_rows) else
    cbind(main_level = character(0), main_tab[0, ], stringsAsFactors = FALSE)
  rownames(nested_tab) <- NULL

  res <- structure(
    list(main = main_tab, nested = nested_tab, levels = levels_tab,
         simper_objects = list(main = main_sim, nested = nested_sims),
         removed = dat$removed, kept = dat$kept,
         settings = list(permutations = permutations, cutoff = cutoff,
                         seed = seed)),
    class = "nested_simper"
  )

  if (verbose) print(res)
  invisible(res)
}


# Internal: turn a vegan::simper() object into one tidy data.frame.
# vegan builds its comparisons as combn(unique(group), 2), in order of first
# appearance; the same call is used here to recover the two group names, since
# splitting the list names on "_" fails when group names contain "_".
.tidy_simper <- function(sim, group, cutoff) {
  pairs <- t(utils::combn(as.character(unique(group)), 2))
  out <- lapply(seq_along(sim), function(i) {
    s <- sim[[i]]
    o <- s$ord
    contrib <- unname(s$average[o])
    cum <- unname(s$cusum)
    df <- data.frame(
      comparison = paste(pairs[i, 1], "vs", pairs[i, 2]),
      group_a = pairs[i, 1],
      group_b = pairs[i, 2],
      species = s$species[o],
      contribution = contrib,
      sd = unname(s$sd[o]),
      ratio = unname(s$ratio[o]),
      mean_a = unname(s$ava[o]),
      mean_b = unname(s$avb[o]),
      pct = 100 * contrib / s$overall,
      cum_pct = 100 * cum,
      p = if (is.null(s$p)) NA_real_ else unname(s$p[o]),
      dissimilarity = s$overall,
      stringsAsFactors = FALSE
    )
    # Keep taxa up to and including the one that reaches the cutoff
    n_keep <- which(cum >= cutoff - sqrt(.Machine$double.eps))[1]
    if (is.na(n_keep)) n_keep <- nrow(df)
    df[seq_len(n_keep), , drop = FALSE]
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}


#' Print method for "nested_simper" objects
#'
#' @param x An object of class "nested_simper".
#' @param digits Number of digits to print.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
#' @method print nested_simper
print.nested_simper <- function(x, digits = 3, ...) {
  s <- x$settings
  cat("=== NESTED (TWO-WAY) SIMPER ===\n")
  cat("Bray-Curtis | permutations:", s$permutations,
      "| taxa shown up to", paste0(100 * s$cutoff, "%"), "of the dissimilarity\n")
  cat("Samples analysed:", length(x$kept))
  if (length(x$removed) > 0) cat(" (", length(x$removed), " empty removed)", sep = "")
  cat("\n")

  show_cols <- c("species", "pct", "cum_pct", "mean_a", "mean_b")
  if (s$permutations > 0) show_cols <- c(show_cols, "p")

  show <- function(tab) {
    for (cmp in unique(tab$comparison)) {
      t <- tab[tab$comparison == cmp, , drop = FALSE]
      cat("\n   ", cmp, " (average dissimilarity = ",
          round(t$dissimilarity[1], digits), ")\n", sep = "")
      t <- t[, show_cols, drop = FALSE]
      num <- vapply(t, is.numeric, logical(1))
      t[num] <- lapply(t[num], round, digits = digits)
      names(t)[names(t) == "mean_a"] <- "mean_A"
      names(t)[names(t) == "mean_b"] <- "mean_B"
      print(t, row.names = FALSE)
    }
  }

  cat("\n1. MAIN FACTOR\n")
  show(x$main)

  cat("\n2. NESTED FACTOR WITHIN EACH MAIN-FACTOR LEVEL\n")
  for (i in seq_len(nrow(x$levels))) {
    r <- x$levels[i, ]
    cat("\n ", r$main_level, " (", r$group_sizes, ")\n", sep = "")
    if (r$note != "") {
      cat("   ", r$note, "\n", sep = "")
    } else {
      show(x$nested[x$nested$main_level == r$main_level, , drop = FALSE])
    }
  }

  cat("\npct: % of the average dissimilarity | mean_A, mean_B: mean abundance ",
      "in the first and second group of each comparison.\n", sep = "")
  invisible(x)
}
