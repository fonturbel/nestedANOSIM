#' Boxplot of ANOSIM rank dissimilarities (between vs. within groups)
#'
#' A ggplot2 version of \code{plot(vegan::anosim(...))}. It shows the ranks of
#' the dissimilarities between groups and within each group, with the R
#' statistic and p-value in the subtitle.
#'
#' @param x Either an object returned by \code{vegan::anosim()} or an object of
#'   class \code{"nested_anosim"} (see \code{\link{nested_anosim}}).
#' @param which Only used when \code{x} is a \code{"nested_anosim"} object.
#'   One of \code{"overall"}, \code{"main"}, or the name of a main-factor level
#'   (e.g., \code{"First"}) to plot the second-factor test run within that level.
#' @param title Plot title. If \code{NULL}, a title is built automatically.
#' @param xlab,ylab Axis labels.
#' @param palette \code{NULL} (default ggplot2 colors), a character vector of
#'   colors, or the name of a \pkg{ggsci} palette (e.g., \code{"startrek"}).
#' @param notch Logical; draw notched boxes (default FALSE).
#' @param base_size Base font size for \code{theme_bw()}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' data(mistletoe_visitors)
#' comm <- mistletoe_visitors[, -(1:3)]
#' res <- nested_anosim(comm, mistletoe_visitors$Year_study,
#'                      mistletoe_visitors$Mistletoe,
#'                      permutations = 199, seed = 1, verbose = FALSE)
#' plot_anosim_box(res, which = "main")
#' plot_anosim_box(res, which = "First", palette = c("grey60", "tomato", "steelblue"))
#'
#' # A plain vegan::anosim() object also works
#' plot_anosim_box(res$overall)
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot labs theme_bw theme element_blank
#'   element_text scale_fill_manual
#' @export
plot_anosim_box <- function(x, which = "overall", title = NULL,
                            xlab = "Category", ylab = "Rank distance",
                            palette = NULL, notch = FALSE, base_size = 14) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.")
  }

  # ---- Resolve which anosim object to plot ---------------------------------
  if (inherits(x, "nested_anosim")) {
    if (identical(which, "overall")) {
      a <- x$overall
      auto_title <- "Overall (main x second factor)"
    } else if (identical(which, "main")) {
      a <- x$main
      auto_title <- "Main factor"
    } else if (which %in% names(x$nested_objects)) {
      a <- x$nested_objects[[which]]
      auto_title <- paste("Within", which)
    } else {
      stop("'which' must be \"overall\", \"main\" or one of: ",
           paste(names(x$nested_objects), collapse = ", "), ".")
    }
  } else if (inherits(x, "anosim")) {
    a <- x
    auto_title <- "ANOSIM"
  } else {
    stop("'x' must be an 'anosim' or a 'nested_anosim' object.")
  }

  if (is.null(title)) title <- auto_title

  # ---- Data: rank dissimilarities by class (Between + each group) ----------
  df <- data.frame(class = a$class.vec, rank = a$dis.rank)
  lv <- c("Between", setdiff(unique(as.character(df$class)), "Between"))
  df$class <- factor(df$class, levels = lv)

  sub <- paste0("R = ", round(unname(a$statistic), 3),
                ", p = ", format.pval(a$signif, digits = 3, eps = 1e-3))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = class, y = rank, fill = class)) +
    ggplot2::geom_boxplot(notch = notch, alpha = 0.7, show.legend = FALSE) +
    ggplot2::labs(title = title, subtitle = sub, x = xlab, y = ylab) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"))

  cols <- .resolve_palette(palette, length(lv))
  if (!is.null(cols)) {
    p <- p + ggplot2::scale_fill_manual(values = stats::setNames(cols, lv))
  }
  p
}


# Internal: turn a palette spec into a vector of n colors (or NULL = default).
# A single string naming a ggsci palette (e.g. "startrek") is resolved via ggsci;
# anything else is treated as a vector of colors and recycled to length n.
.ggsci_palettes <- c("startrek", "npg", "lancet", "jco", "nejm", "aaas", "d3",
                     "futurama", "simpsons", "uchicago", "locuszoom", "igv",
                     "rickandmorty", "tron", "observable")

.resolve_palette <- function(palette, n) {
  if (is.null(palette)) return(NULL)
  if (length(palette) == 1 && palette %in% .ggsci_palettes) {
    if (!requireNamespace("ggsci", quietly = TRUE)) {
      stop("Palette '", palette, "' requires the 'ggsci' package.")
    }
    fn <- getExportedValue("ggsci", paste0("pal_", palette))
    cols <- fn()(max(n, 3))
    return(cols[seq_len(n)])
  }
  rep_len(palette, n)
}
