# Approximate Modes of a Continuous Family at Many Parameter Settings

The grid search of
[`find_lp_anchor()`](https://statmodels7.github.io/distributions7/reference/find_lp_anchor.md),
run on every parameter setting at once: each refinement evaluates the
log-density on a grid of 129 points per setting in one call, and keeps
the three grid points around the largest value of each.

## Usage

``` r
quantile_anchors(distrib, theta)
```

## Arguments

- distrib:

  A continuous family.

- theta:

  A named list of parameter vectors of a common length.

## Value

A numeric vector of approximate modes, one per setting.

## See also

[`find_pdf_anchor()`](https://statmodels7.github.io/distributions7/reference/find_pdf_anchor.md),
the search for one setting;
[`distrib_quantile.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.continuous_distrib.md),
its consumer.
