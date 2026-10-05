# Internal helpers shared by nested_anosim() and nested_simper()

# Check the inputs and drop empty samples (rows with zero total abundance).
# `strata` (optional) is checked and trimmed in step with the samples.
# Returns the cleaned data plus the indices kept and removed.
.prepare_input <- function(comm, main, nested, remove_empty, strata = NULL) {
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
  if (!is.null(strata)) {
    if (length(strata) != n) stop("'strata' must have the same length as nrow(comm).")
    if (anyNA(strata)) stop("'strata' must not contain NA.")
    strata <- as.character(strata)
  }

  main    <- as.character(main)
  nested  <- as.character(nested)
  kept    <- seq_len(n)
  removed <- integer(0)

  empty <- which(rowSums(comm) == 0)
  if (length(empty) > 0) {
    if (remove_empty) {
      removed <- empty
      kept    <- setdiff(kept, empty)
      comm    <- comm[kept, , drop = FALSE]
      main    <- main[kept]
      nested  <- nested[kept]
      if (!is.null(strata)) strata <- strata[kept]
      warning(length(empty), " empty sample(s) removed (rows: ",
              paste(empty, collapse = ", "), ").", call. = FALSE)
    } else {
      stop("Found ", length(empty), " empty sample(s); set remove_empty = TRUE ",
           "or remove them beforehand.")
    }
  }

  if (length(unique(main)) < 2) stop("'main' needs at least two levels.")

  list(comm = comm, main = main, nested = nested, strata = strata,
       kept = kept, removed = removed)
}

# Reason for skipping a second-factor test within one main-factor level ("" = run it).
# `sizes` is the table of nested-group sizes within that level.
.skip_note <- function(sizes) {
  if (length(sizes) < 2) return("skipped: fewer than 2 nested groups")
  if (any(sizes < 2)) return("skipped: a nested group has < 2 samples")
  ""
}
