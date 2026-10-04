#' nMDS ordination with convex hulls for any number of groups
#'
#' Runs \code{vegan::metaMDS()} on a dissimilarity matrix (or on a community
#' matrix) and draws the ordination with ggplot2, one convex hull per group.
#' It generalizes the two-group code used in the one-way ANOSIM workflow to an
#' arbitrary number of groups.
#'
#' @param x A \code{dist} object, a community matrix/data.frame (rows =
#'   samples), a \code{"nested_anosim"} object (its stored dissimilarity matrix
#'   is used), or an already fitted \code{metaMDS} object.
#' @param groups Vector (factor/character) with the group of each sample, in the
#'   same order as the rows of \code{x}. For a \code{"nested_anosim"} object,
#'   \code{groups} must refer to the \emph{original} samples; samples removed
#'   as empty are dropped automatically.
#' @param method Dissimilarity index used when \code{x} is a community matrix
#'   (default \code{"bray"}).
#' @param k Number of nMDS dimensions (default 2).
#' @param trymax Maximum number of random starts for \code{metaMDS()}.
#' @param seed Optional integer for reproducible ordinations.
#' @param palette \code{NULL} (default ggplot2 colors), a character vector of
#'   colors, or the name of a \pkg{ggsci} palette (e.g., \code{"startrek"}).
#' @param hull Logical; draw convex hulls (default TRUE). Groups with fewer than
#'   three samples get no hull.
#' @param hull_alpha Transparency of the hulls.
#' @param point_size Size of the points.
#' @param show_stress Logical; print the stress value on the plot (default TRUE).
#' @param legend_title Legend title (\code{NULL} = none).
#' @param base_size Base font size for \code{theme_bw()}.
#'
#' @return A list with \code{plot} (ggplot object), \code{mds} (the
#'   \code{metaMDS} fit), \code{scores} (data.frame of site scores plus group)
#'   and \code{stress}.
#'
#' @examples
#' data(mistletoe_visitors)
#' comm <- mistletoe_visitors[, -(1:3)]
#' res <- nested_anosim(comm, mistletoe_visitors$Year_study,
#'                      mistletoe_visitors$Mistletoe,
#'                      permutations = 99, seed = 1, verbose = FALSE)
#' # `groups` refers to the original rows; empty samples are dropped for you
#' out <- plot_nmds_hulls(res, groups = mistletoe_visitors$Mistletoe,
#'                        trymax = 20, seed = 1)
#' out$plot
#' out$stress
#'
#' @importFrom ggplot2 ggplot aes geom_polygon geom_point coord_equal labs
#'   theme_bw theme element_blank element_text scale_fill_manual
#'   scale_colour_manual scale_shape_manual annotate
#' @export
plot_nmds_hulls <- function(x, groups, method = "bray", k = 2, trymax = 100,
                            seed = NULL, palette = NULL, hull = TRUE,
                            hull_alpha = 0.30, point_size = 4,
                            show_stress = TRUE, legend_title = NULL,
                            base_size = 14) {

  for (pkg in c("vegan", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for this function.")
    }
  }

  # ---- Resolve input -------------------------------------------------------
  groups <- as.character(groups)

  if (inherits(x, "nested_anosim")) {
    if (length(groups) == length(x$kept) + length(x$removed)) {
      groups <- groups[x$kept]          # drop samples removed as empty
    }
    x <- x$dist
  }

  if (inherits(x, "metaMDS")) {
    mds <- x
  } else {
    if (!inherits(x, "dist")) {
      x <- vegan::vegdist(as.data.frame(x), method = method)
    }
    if (!is.null(seed)) set.seed(seed)
    mds <- vegan::metaMDS(x, k = k, trymax = trymax, autotransform = FALSE,
                          trace = 0)
  }

  sc <- as.data.frame(vegan::scores(mds, display = "sites"))
  if (nrow(sc) != length(groups)) {
    stop("Length of 'groups' (", length(groups), ") does not match the number ",
         "of ordinated samples (", nrow(sc), ").")
  }
  if (ncol(sc) < 2) stop("At least two nMDS axes are needed to plot.")
  names(sc)[1:2] <- c("NMDS1", "NMDS2")
  sc$group <- factor(groups, levels = unique(groups))
  lv <- levels(sc$group)

  # ---- Convex hulls, one per group ----------------------------------------
  hulls <- NULL
  if (hull) {
    hl <- lapply(lv, function(g) {
      d <- sc[sc$group == g, , drop = FALSE]
      if (nrow(d) < 3) return(NULL)
      d[grDevices::chull(d$NMDS1, d$NMDS2), , drop = FALSE]
    })
    hulls <- do.call(rbind, hl)
  }

  # ---- Plot ----------------------------------------------------------------
  p <- ggplot2::ggplot()
  if (!is.null(hulls) && nrow(hulls) > 0) {
    p <- p + ggplot2::geom_polygon(
      data = hulls,
      ggplot2::aes(x = NMDS1, y = NMDS2, fill = group, group = group),
      alpha = hull_alpha)
  }
  p <- p +
    ggplot2::geom_point(
      data = sc,
      ggplot2::aes(x = NMDS1, y = NMDS2, shape = group, colour = group),
      size = point_size) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "nMDS 1", y = "nMDS 2",
                  fill = legend_title, colour = legend_title,
                  shape = legend_title) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank())
  if (is.null(legend_title)) {
    p <- p + ggplot2::theme(legend.title = ggplot2::element_blank())
  }

  # Shapes recycle through a fixed set so any number of groups works
  shapes <- rep_len(c(16, 17, 15, 18, 3, 4, 8, 7, 9, 10, 11, 12), length(lv))
  p <- p + ggplot2::scale_shape_manual(values = stats::setNames(shapes, lv))

  cols <- .resolve_palette(palette, length(lv))
  if (!is.null(cols)) {
    cols <- stats::setNames(cols, lv)
    p <- p + ggplot2::scale_fill_manual(values = cols) +
      ggplot2::scale_colour_manual(values = cols)
  }

  stress <- mds$stress
  if (show_stress && is.finite(stress)) {
    p <- p + ggplot2::annotate("text", x = -Inf, y = -Inf, hjust = -0.1,
                               vjust = -0.8,
                               label = paste0("Stress = ", round(stress, 3)),
                               size = base_size / 3)
  }

  list(plot = p, mds = mds, scores = sc, stress = stress)
}
