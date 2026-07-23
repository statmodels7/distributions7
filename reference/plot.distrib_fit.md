# Plot a Fitted Distribution against the Data

Compares the fitted distribution with the sample it was estimated from.
For a continuous distribution the observations are summarised by a
kernel density estimate, with the fitted density drawn on top and a rug
of the data underneath. For a discrete one the observed relative
frequencies are drawn as bars with the fitted probability mass overlaid,
since a kernel density would misrepresent a lattice-valued sample.

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

  Colours of the fitted curve and of the empirical summary.

- ...:

  Further arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html), for instance
  `main`, `xlab` or `xlim`.

## Value

`x`, invisibly.

## See also

[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
[`simulate.distrib_fit`](https://statmodels7.github.io/distributions7/reference/simulate.distrib_fit.md)

## Examples

``` r
set.seed(1)
y <- rgamma(300, shape = 4, rate = 2)
fit <- fit_distrib(gamma_distrib(), y)
plot(fit)

```
