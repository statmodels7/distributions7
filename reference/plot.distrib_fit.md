# Plot a Fitted Distribution against the Data

Compares the fitted distribution with the sample it was estimated from.
For a continuous distribution the observations are summarized by a
kernel density estimate, with the fitted density drawn on top and a rug
of the data underneath. For a discrete one the observed relative
frequencies are drawn as bars with the fitted probability mass overlaid,
since a kernel density would misrepresent a discrete sample.

## Arguments

- x:

  A
  [`distrib_fit`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
  object.

- n_grid:

  Number of points at which the fitted density is evaluated (continuous
  distributions only). Defaults to 512.

- rug:

  Logical; draw a rug of the observations. Defaults to `TRUE` when there
  are at most 2000 of them.

- legend:

  Logical; add a legend. Defaults to `TRUE`.

- col_fit, col_data:

  Colors of the fitted curve and of the empirical summary.

- mv_which:

  For a multivariate fit, which coordinates to show. Defaults to all of
  them, and at most three are drawn: above that the panel matrix stops
  being readable.

- ...:

  Further arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html), for instance
  `main`, `xlab` or `xlim`.

## Value

`x`, invisibly.

## Details

A univariate fit is drawn as the fitted density over a histogram or, for
a discrete family, over the observed proportions. A multivariate one is
drawn as a panel matrix: the fitted marginal density and a kernel
estimate of the data on the diagonal, the fitted contours over the
observations below it, and the fitted correlation above.

## See also

[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
[`simulate.distrib_fit`](https://statmodels7.github.io/distributions7/reference/simulate.distrib_fit.md),
[`plot.multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md)

## Examples

``` r
set.seed(1)
y <- rgamma(300, shape = 4, rate = 2)
fit <- fit_distrib(gamma2_distrib(), y)
plot(fit)

```
