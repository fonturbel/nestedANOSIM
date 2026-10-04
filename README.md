# nestedANOSIM

[![R-CMD-check](https://github.com/fonturbel/nestedANOSIM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fonturbel/nestedANOSIM/actions/workflows/R-CMD-check.yaml)

### nestedANOSIM package version 0.0-1

A small toolkit for **hierarchical (two-way, nested) ANOSIM** on community data, plus ggplot2 graphics for ANOSIM and nMDS results. The `vegan` package only implements one-way ANOSIM (`vegan::anosim()`); this package wraps it to evaluate a nested factor (e.g., mistletoe species) within each level of a main factor (e.g., study year) in a single call.

This package was developed by Francisco E. Fontúrbel from code originally written for camera-trap community analyses at Las Chinchillas National Reserve, Chile. Software development was made using Claude Code (model Opus 5.5).

## What this package computes

- **Overall test**: ANOSIM on the combined `main.nested` factor.
- **Main factor test**: ANOSIM on the main factor alone.
- **Nested tests**: ANOSIM of the nested factor within each level of the main factor, with optional multiple-testing correction (`p_adjust`).
- **Group summary**: sample size, mean richness and mean abundance per `main.nested` group.
- **Plots**: rank-dissimilarity boxplots (`plot_anosim_box()`) and nMDS ordinations with convex hulls for any number of groups (`plot_nmds_hulls()`).

> **Note**: the nested tests are separate ANOSIMs run within each level of the main factor, not a variance partition as in a nested PERMANOVA (`vegan::adonis2()`). Report them as a hierarchical ANOSIM.

## Installation

```r
# install.packages("remotes")
remotes::install_github("fonturbel/nestedANOSIM")

# with the vignette
remotes::install_github("fonturbel/nestedANOSIM", build_vignettes = TRUE)
vignette("mistletoe-visitors", package = "nestedANOSIM")
```

## Dependencies

```r
install.packages(c("vegan", "ggplot2"))
install.packages("ggsci")   # optional, only for named palettes (e.g., "startrek")
```

| Package   | Used for                                                  |
| --------- | --------------------------------------------------------- |
| `vegan`   | dissimilarities, `anosim()`, `metaMDS()`                  |
| `ggplot2` | `plot_anosim_box()`, `plot_nmds_hulls()`                  |
| `ggsci`   | optional named color palettes                             |
| `testthat`, `rmarkdown`, `knitr` | tests and vignettes (dev only)     |

## Input format

- A **community matrix**: rows = samples (e.g., cameras), columns = taxa, numeric (abundances or counts). Keep only community columns; remove any taxa you want to exclude beforehand.
- A **main factor** and a **nested factor**: vectors of the same length as `nrow(comm)`.
- Samples with zero total abundance are removed automatically (with a warning), since Bray-Curtis is undefined for empty rows.
- With `seed`, the seed is reset before each ANOSIM, so every test is reproducible on its own. The main-factor test is identical to `set.seed(seed); vegan::anosim(dist, main)`.

## Quick start

```r
library(nestedANOSIM)

data(mistletoe_visitors)
comm <- mistletoe_visitors[, -(1:3)]   # community data only
meta <- mistletoe_visitors[, 1:3]      # Year_study, Mistletoe, Camera

res <- nested_anosim(comm   = comm,
                     main   = meta$Year_study,
                     nested = meta$Mistletoe,
                     method = "bray", permutations = 9999,
                     p_adjust = "holm", seed = 123)

res$overall            # anosim object, combined factor
res$main               # anosim object, main factor
res$nested             # table: nested test within each main-factor level
res$group_summary

plot_anosim_box(res, which = "main")
plot_anosim_box(res, which = "First")     # nested test within a main-factor level

# `groups` refers to the original rows; samples removed as empty are dropped
out <- plot_nmds_hulls(res, groups = meta$Mistletoe, palette = "startrek", seed = 1)
out$plot
out$stress
```

## Example data

The package includes a camera-trap dataset of animal visitors to the mistletoes *Tristerix aphyllus* and *T. verticillatus* in Las Chinchillas National Reserve, northern Chile (two sampling years), published in [doi:10.1016/j.jaridenv.2025.105518](https://doi.org/10.1016/j.jaridenv.2025.105518). It is available as `mistletoe_visitors` (87 camera-years x 40 taxa; the focal species and unidentified records are excluded) and as a semicolon-separated file at `system.file("extdata", "mistletoe_visitors.csv", package = "nestedANOSIM")`. The vignette reproduces the full analysis.

## Citation

Fontúrbel FE (2026) nestedANOSIM: hierarchical ANOSIM and graphics for community data. <https://github.com/fonturbel/nestedANOSIM>

If you use the example data, please also cite the paper above.

## License

GPL (>= 3)
