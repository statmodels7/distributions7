# Panels of a Multivariate Density

Draws a multivariate distribution as a matrix of panels: the marginal
density of each coordinate on the diagonal, and contours of each
bivariate marginal below it.

## Arguments

- x:

  An object inheriting from class
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md).

- theta:

  A named list or vector of parameters. Generated at random when
  missing, as for a univariate distribution.

- which:

  The coordinates to show. Defaults to all of them.

- n_grid:

  Points per axis for the marginal densities and the contours.

- col_fit:

  Colour of the density and of the contours.

- ...:

  Passed to the underlying plotting calls.

## Value

`x`, invisibly.

## Details

The picture is built from marginals, so it exists exactly when
[`mv_marginal`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
does. Above about three coordinates the panel matrix stops being
readable, so a larger distribution is refused with the suggestion of
choosing coordinates rather than being drawn illegibly.

The contours are drawn at levels of the density itself rather than at
probability levels: the equal-density contour is what the density's
shape means, and computing a probability level would need the integral
this package refuses to approximate in several dimensions.

## See also

[`mv_marginal`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md),
[`plot.distrib_fit`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md)

## Examples

``` r
d <- mvgaussian_distrib(3)
theta <- as.list(stats::setNames(
  c(0, 1, -1, 0, 0, 0, 0.6, -0.3, 0.2), d@params
))
plot(d, theta)

```
