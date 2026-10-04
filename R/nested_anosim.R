#' Two-way (hierarchical) ANOSIM
#'
#' vegan::anosim() only handles one grouping factor. This function wraps it to
#' run a "two-way" ANOSIM in which a nested factor (e.g., mistletoe species)
#' is evaluated within each level of a main factor (e.g., study year):
#'   1. Overall test on the combined main.nested factor.
#'   2. Test of the main factor alone.
#'   3. Test of the nested factor separately within each level of the main factor
#'      (optionally with multiple-testing correction).
#'
#' @param comm        Community matrix or data.frame (rows = samples, columns = taxa).
#'                    Must be numeric and contain only the community data.
#' @param main        Vector (factor/character) with the main grouping factor; length = nrow(comm).
#' @param nested      Vector (factor/character) with the factor nested within `main`.
#' @param method      Dissimilarity index passed to vegan::vegdist() (default "bray").
#' @param binary      Logical; passed to vegdist() (presence/absence transformation).
#' @param permutations Number of permutations for anosim() (default 999).
#' @param remove_empty Logical; drop samples whose total abundance is zero (default TRUE).
#'                    Needed because Bray-Curtis is undefined for empty rows.
#' @param p_adjust    Method for p.adjust() applied to the nested tests (default "none";
#'                    e.g., "holm", "bonferroni", "BH").
#' @param seed        Optional integer. If given, the RNG seed is set to this value
#'                    before \emph{each} ANOSIM, so every test is reproducible on its
#'                    own (e.g., the main-factor test equals
#'                    \code{set.seed(seed); vegan::anosim(dist, main)}).
#' @param verbose     Logical; print the formatted report when the function runs (default TRUE).
#' @param ...         Further arguments passed to vegan::anosim() (e.g., parallel = 4).
#'
#' @return An object of class "nested_anosim", a list with:
#'   \item{overall}{anosim object for the combined main.nested factor}
#'   \item{main}{anosim object for the main factor}
#'   \item{nested}{data.frame with the nested-factor test for each main-factor level}
#'   \item{nested_objects}{list with the full anosim objects of the nested tests}
#'   \item{group_summary}{data.frame: N, mean richness and mean abundance per main.nested group}
#'   \item{dist}{the dissimilarity matrix used}
#'   \item{removed}{indices of the samples removed (empty rows)}
#'   \item{kept}{indices of the samples retained (relative to the input)}
#'   \item{settings}{list with the main settings used}
#'
#' @examples
#' data(mistletoe_visitors)
#' comm <- mistletoe_visitors[, -(1:3)]
#' res <- nested_anosim(comm   = comm,
#'                      main   = mistletoe_visitors$Year_study,
#'                      nested = mistletoe_visitors$Mistletoe,
#'                      permutations = 199, p_adjust = "holm", seed = 123)
#' res$nested          # table of nested tests
#' res$group_summary
#'
#' @importFrom vegan vegdist anosim
#' @importFrom stats as.dist p.adjust
#' @export
nested_anosim <- function(comm, main, nested,
                          method = "bray",
                          binary = FALSE,
                          permutations = 999,
                          remove_empty = TRUE,
                          p_adjust = "none",
                          seed = NULL,
                          verbose = TRUE,
                          ...) {

  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required. Install it with install.packages('vegan').")
  }

  # ---- Input checks ---------------------------------------------------------
  comm <- as.data.frame(comm)
  if (!all(vapply(comm, is.numeric, logical(1)))) {
    stop("'comm' must contain only numeric columns (community data).")
  }
  if (anyNA(comm)) stop("'comm' contains NA values.")
  n <- nrow(comm)
  if (length(main) != n || length(nested) != n) {
    stop("'main' and 'nested' must have the same length as nrow(comm).")
  }
  if (anyNA(main) || anyNA(nested)) stop("'main' and 'nested' must not contain NA.")

  main   <- as.character(main)
  nested <- as.character(nested)
  kept   <- seq_len(n)
  removed <- integer(0)

  # ---- Remove empty samples -------------------------------------------------
  empty <- which(rowSums(comm) == 0)
  if (length(empty) > 0) {
    if (remove_empty) {
      removed <- empty
      kept    <- setdiff(kept, empty)
      comm    <- comm[kept, , drop = FALSE]
      main    <- main[kept]
      nested  <- nested[kept]
      warning(length(empty), " empty sample(s) removed (rows: ",
              paste(empty, collapse = ", "), ").", call. = FALSE)
    } else {
      stop("Found ", length(empty), " empty sample(s); set remove_empty = TRUE ",
           "or remove them beforehand.")
    }
  }

  if (length(unique(main)) < 2) stop("'main' needs at least two levels.")

  # ---- Distance matrix ------------------------------------------------------
  d <- vegan::vegdist(comm, method = method, binary = binary)
  d_mat <- as.matrix(d)

  # Reset the seed before each test so each one is reproducible on its own
  run_anosim <- function(dis, grouping) {
    if (!is.null(seed)) set.seed(seed)
    vegan::anosim(dis, grouping, permutations = permutations, ...)
  }

  # ---- 1. Overall test (combined factor) -----------------------------------
  combined <- paste(main, nested, sep = ".")
  overall  <- run_anosim(d, combined)

  # ---- 2. Main factor -------------------------------------------------------
  main_test <- run_anosim(d, main)

  # ---- 3. Nested factor within each level of main --------------------------
  nested_objects <- list()
  rows <- list()
  for (lv in sort(unique(main))) {
    idx <- which(main == lv)
    sub_nested <- nested[idx]
    sizes <- table(sub_nested)
    out <- data.frame(main_level = lv,
                      n_groups = length(sizes),
                      n_samples = length(idx),
                      group_sizes = paste(names(sizes), sizes, sep = "=", collapse = ", "),
                      R = NA_real_, p = NA_real_,
                      note = "", stringsAsFactors = FALSE)

    if (length(sizes) < 2) {
      out$note <- "skipped: fewer than 2 nested groups"
    } else if (any(sizes < 2)) {
      out$note <- "skipped: a nested group has < 2 samples"
    } else {
      sub_d <- stats::as.dist(d_mat[idx, idx, drop = FALSE])
      a <- run_anosim(sub_d, sub_nested)
      nested_objects[[lv]] <- a
      out$R <- unname(a$statistic)
      out$p <- a$signif
    }
    rows[[lv]] <- out
  }
  nested_tab <- do.call(rbind, rows)
  rownames(nested_tab) <- NULL
  nested_tab$p_adj <- NA_real_
  ok <- !is.na(nested_tab$p)
  nested_tab$p_adj[ok] <- stats::p.adjust(nested_tab$p[ok], method = p_adjust)

  # ---- Group summary --------------------------------------------------------
  grp_levels <- sort(unique(combined))
  group_summary <- data.frame(
    Group = grp_levels,
    N_samples = as.vector(table(factor(combined, levels = grp_levels))),
    Mean_richness = as.vector(tapply(rowSums(comm > 0), factor(combined, levels = grp_levels), mean)),
    Mean_abundance = as.vector(tapply(rowSums(comm), factor(combined, levels = grp_levels), mean)),
    stringsAsFactors = FALSE
  )

  res <- structure(
    list(overall = overall, main = main_test,
         nested = nested_tab, nested_objects = nested_objects,
         group_summary = group_summary, dist = d,
         removed = removed, kept = kept,
         settings = list(method = method, binary = binary,
                         permutations = permutations, p_adjust = p_adjust,
                         seed = seed)),
    class = "nested_anosim"
  )

  if (verbose) print(res)
  invisible(res)
}


# Significance stars helper (internal)
.anosim_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p <= 0.001, "***",
                ifelse(p <= 0.01, "**",
                       ifelse(p <= 0.05, "*", "ns"))))
}


#' Print method for "nested_anosim" objects
#'
#' @param x An object of class "nested_anosim".
#' @param digits Number of digits to print.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
#' @method print nested_anosim
print.nested_anosim <- function(x, digits = 4, ...) {
  s <- x$settings
  cat("=== NESTED (TWO-WAY) ANOSIM ===\n")
  cat("Dissimilarity:", s$method, "| permutations:", s$permutations,
      "| nested p-adjust:", s$p_adjust, "\n")
  cat("Samples analysed:", length(x$kept))
  if (length(x$removed) > 0) cat(" (", length(x$removed), " empty removed)", sep = "")
  cat("\n\n")

  cat("1. OVERALL TEST (main.nested combined)\n")
  cat("   R =", round(x$overall$statistic, digits),
      "| p =", round(x$overall$signif, digits),
      .anosim_stars(x$overall$signif), "\n\n")

  cat("2. MAIN FACTOR\n")
  cat("   R =", round(x$main$statistic, digits),
      "| p =", round(x$main$signif, digits),
      .anosim_stars(x$main$signif), "\n\n")

  cat("3. NESTED FACTOR WITHIN EACH MAIN-FACTOR LEVEL\n")
  for (i in seq_len(nrow(x$nested))) {
    r <- x$nested[i, ]
    cat("   ", r$main_level, " (", r$group_sizes, ")\n", sep = "")
    if (is.na(r$R)) {
      cat("      ", r$note, "\n", sep = "")
    } else {
      cat("      R =", round(r$R, digits), "| p =", round(r$p, digits))
      if (s$p_adjust != "none") cat(" | p_adj =", round(r$p_adj, digits))
      cat(" ", .anosim_stars(if (s$p_adjust != "none") r$p_adj else r$p), "\n", sep = "")
    }
  }

  cat("\n4. GROUP SUMMARY\n")
  print(x$group_summary, row.names = FALSE)

  cat("\nR near 1: groups well separated | R near 0: no separation | R < 0: ",
      "within-group dissimilarities exceed between-group ones.\n", sep = "")
  cat("Significance: *** p<=0.001, ** p<=0.01, * p<=0.05, ns = not significant\n")
  invisible(x)
}
