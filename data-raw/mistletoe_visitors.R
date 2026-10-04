# Build the `mistletoe_visitors` example dataset from the raw camera-trap matrix.
#
# Source: sqmat_all.csv, Las Chinchillas National Reserve project
# (doi:10.1016/j.jaridenv.2025.105518). Semicolon-separated; columns 1-3 are
# metadata (Year_study, Mistletoe, Camera), the rest are visitor taxa.
#
# The focal species (Mimus_thenca, Sephanoides_sephaniodes) and unidentified
# records (not_identified) are excluded, as in the original analysis.
# Empty samples are kept on purpose: nested_anosim() removes them with a warning.

raw <- utils::read.csv("data-raw/sqmat_all.csv", sep = ";",
                       check.names = FALSE, stringsAsFactors = FALSE)

drop <- c("Mimus_thenca", "Sephanoides_sephaniodes", "not_identified")
stopifnot(all(drop %in% names(raw)))
mistletoe_visitors <- raw[, setdiff(names(raw), drop)]
mistletoe_visitors <- mistletoe_visitors[, c("Year_study", "Mistletoe", "Camera",
                                             setdiff(names(mistletoe_visitors),
                                                     c("Year_study", "Mistletoe", "Camera")))]

# Checks
stopifnot(
  nrow(mistletoe_visitors) == 87,
  ncol(mistletoe_visitors) == 43,
  !anyNA(mistletoe_visitors),
  all(vapply(mistletoe_visitors[, -(1:3)], is.numeric, logical(1)))
)

usethis::use_data(mistletoe_visitors, overwrite = TRUE, compress = "xz")

# Plain-text copy for users who want to try reading from a file
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
utils::write.table(mistletoe_visitors, "inst/extdata/mistletoe_visitors.csv",
                   sep = ";", row.names = FALSE, quote = FALSE)
