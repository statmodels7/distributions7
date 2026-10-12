# Draw the Panel Matrix of a Multivariate Density

The shared engine of
[`plot.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md)
and the multivariate branch of
[`plot.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md).
It lays out one panel per ordered pair of the chosen coordinates: the
marginal density on the diagonal, contours of the bivariate marginal
below it, and the implied correlation printed above it. The density
panels come from
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md),
so nothing is integrated; the printed correlation comes from
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md).

## Usage

``` r
mv_pairs_panels(
  d,
  theta,
  which,
  n_grid,
  col_fit,
  data,
  col_data = "#4682B4",
  which_arg = "which",
  ...
)
```

## Arguments

- d:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list of parameters, already aligned by the caller.

- which:

  An integer vector of distinct coordinates, or `NULL` for all of them.
  An empty, non-integer, out-of-range or duplicated `which`, or more
  than three coordinates, is an error.

- n_grid:

  Points per axis, a single positive whole number.

- col_fit:

  Color of the fitted density, the contours and the printed
  correlations.

- data:

  An \\n \times p\\ matrix of observations with one column per
  coordinate of the distribution, whichever coordinates `which` selects,
  or `NULL` for a distribution with no data behind it.

- col_data:

  Color of the observed summary, used only when `data` is given.
  Defaults to `"#4682B4"`.

- which_arg:

  The name of the caller's argument that carries `which`, used in the
  error message for more than three coordinates. Defaults to `"which"`;
  [`plot.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md)
  passes `"mv_which"`.

- ...:

  Unused.

## Value

`NULL`, invisibly. Called for the plot it draws. The function sets
`mfrow`, `mar` and `mgp` for the panel matrix and restores the previous
values on exit.

## Details

When `data` is supplied the diagonal also carries a kernel density
estimate of the observed coordinate and the panels below the diagonal
carry the observations, so the fitted shape and the sample are read
against each other in one frame. The panels above the diagonal show only
the correlation. The kernel estimate
([`stats::density()`](https://rdrr.io/r/stats/density.html)) does not
depend on the fitted family and serves as a model-free comparison.

Each coordinate is drawn over three standard deviations either side of
the location, widened when `data` is supplied to include the
observations. The spread is taken from
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md),
so a Student \\t\\ with \\\nu \le 2\\ still has a finite range, where
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
returns an infinite matrix.

## See also

[`plot.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md)
and
[`plot.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md),
its two callers, and
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
for the panels.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0,
              sigma_log_L2 = 0, sigma_L2.1 = 0.6)

# Without data: the density alone.
op <- graphics::par(no.readonly = TRUE)
distributions7:::mv_pairs_panels(d, theta, NULL, 40, "#B22222", data = NULL)

graphics::par(op)

# With data: the same frame carries a kernel estimate and the observations.
set.seed(1)
op <- graphics::par(no.readonly = TRUE)
distributions7:::mv_pairs_panels(d, theta, NULL, 40, "#B22222",
                                 data = distrib_rng(d, 200, theta))

graphics::par(op)
```
