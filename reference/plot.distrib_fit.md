# Plot a Fitted Distribution Against the Data

Compares the fitted distribution with the sample it was estimated from,
choosing the comparison from the family's support. A continuous family
is drawn as the fitted density over a kernel estimate of the data, with
a rug of the observations underneath. A discrete one is drawn as the
observed relative frequencies in bars with the fitted mass overlaid, a
kernel estimate being a misrepresentation of a lattice sample.

A multivariate fit is drawn as a panel matrix: the fitted marginal
density and a kernel estimate of the data on the diagonal, the fitted
contours over the observations below it, the fitted correlation above.

## Arguments

- x:

  A
  [`distrib_fit()`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
  object.

- n_grid:

  Number of points at which the fitted density is evaluated, a single
  positive integer. Defaults to 512. Read for a continuous family only;
  a discrete one is evaluated on its support.

- rug:

  Logical of length 1: draw a rug of the observations. Defaults to
  `TRUE` when there are at most 2000 of them, above which the rug
  becomes a solid band and says nothing.

- legend:

  Logical of length 1: add a legend. Defaults to `TRUE`.

- col_fit, col_data:

  Colors of the fitted curve and of the empirical summary, in any form
  [`grDevices::col2rgb()`](https://rdrr.io/r/grDevices/col2rgb.html)
  accepts.

- mv_which:

  For a multivariate fit, which coordinates to show, by index. Defaults
  to all of them, of which at most three are drawn: above that the panel
  matrix stops being readable. Ignored for a univariate fit.

- ...:

  Further arguments passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `xlab` or `xlim`.

## Value

`x`, invisibly. Called for the plot it draws.

## See also

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
for the fit;
[`simulate.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/simulate.distrib_fit.md)
for draws from it;
[`plot.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md),
the same panel matrix drawn from a distribution with no data beside it.

## Examples

``` r
set.seed(1)
y <- rgamma(300, shape = 4, rate = 2)
fit <- fit_distrib(gamma2_distrib(), y)
plot(fit)


# A count family is compared on its support, not through a kernel estimate.
fp <- fit_distrib(poisson_distrib(), rpois(300, 4))
plot(fp, main = "Poisson fit")

```
