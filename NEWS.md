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
