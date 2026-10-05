# nestedANOSIM 0.0-3

* New `strata` argument in `nested_anosim()` and `nested_simper()`: restricts
  the permutations of the main-factor test within blocks, for sampling units
  measured in every main-factor level (e.g., the same cameras every year). It
  is trimmed together with empty samples. Before, passing `strata` through
  `...` failed, and would have been applied to the wrong samples.
* `nested_anosim()` now stops if combining `main` and `nested` with "." gives
  ambiguous group labels.
* Documentation, vignette and plot titles now make clear that the second
  factor can be nested in or crossed with the main factor; in the example data
  mistletoe species is crossed with year. The vignette has a new study-design
  section and uses `strata = Camera`.
* `permute` added to Imports (already a dependency of `vegan`).

# nestedANOSIM 0.0-2

* New `nested_simper()`: SIMPER for the main factor and for the nested factor
  within each main-factor level, returned as tidy data frames and trimmed to a
  cumulative `cutoff` (default 70%), with a `print()` method.
* Input checks and empty-sample removal are now shared by `nested_anosim()`
  and `nested_simper()`.
* The vignette has a new SIMPER section.

# nestedANOSIM 0.0-1

* First version: `nested_anosim()`, `plot_anosim_box()`, `plot_nmds_hulls()`
  and the `mistletoe_visitors` example dataset.
