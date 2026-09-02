# Panels of a Multivariate Density

Draws a multivariate distribution as a matrix of panels: the marginal
density of each coordinate on the diagonal, contours of each bivariate
marginal below it, and the implied correlation printed above it. The
picture is built entirely from marginals, so it exists exactly where
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
does, which for the two families that ship is everywhere.

## Arguments

- x:

  An object inheriting from
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md).

- theta:

  A named list or vector of parameters, each component a single number.
  When missing, parameters are drawn by
  [`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md)
  and reported in a message. A component of length greater than one is
  an error.

- which:

  An integer vector of the coordinates to show, at most three. Defaults
  to `NULL`, which means all of them and is itself an error above three
  coordinates.

- n_grid:

  Points per axis for the marginal densities and the contours. Defaults
  to `80`. The bivariate panels evaluate the density on an
  `n_grid * n_grid` grid, so the cost is quadratic in it.

- col_fit:

  Color of the density curves, the contours and the printed
  correlations. Defaults to `"#B22222"`.

- ...:

  Passed to the underlying plotting calls.

## Value

`x`, invisibly. Called for the plot it draws.

## Why at most three coordinates

A panel matrix of \\k\\ coordinates has \\k^2\\ panels, and past three
the panels are too small to read. A larger distribution is rejected with
an error naming the count and suggesting `which`; drawing it illegibly
would be the worse answer.

## What the contours are levels of

The contours are levels of the DENSITY, not of the probability. An
equal-density contour is what the shape of the density means, and a
probability level would need the orthant integral this package does not
approximate in several dimensions. The levels are chosen from the
density's own range on the grid.

## One setting only

A univariate [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
method reads a `theta` component of length \\k\\ as \\k\\ settings drawn
over one another. Here the picture is already a matrix of panels, with
no axis left to overlay them on, so a component longer than one is
rejected by name.

## The range each panel is drawn over

Two and a half standard deviations either side of the location, with the
spread taken from the matrix the parametrization carries. A multivariate
Student t at \\\nu \le 2\\ has no variance at all, and that is precisely
the shape worth drawing, so
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
cannot be the source. For an elliptical family the two agree up to a
factor, which is all a plotting range needs.

## See also

[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md),
which supplies every panel,
[`plot.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md),
which draws the same matrix with the data on it, and
[`plot.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.continuous_distrib.md)
for the one-dimensional case, which does overlay several settings.

## Examples

``` r
d <- mvgaussian1_distrib(3)
theta <- as.list(stats::setNames(
  c(0, 1, -1, 0, 0, 0, 0.6, -0.3, 0.2), d@params
))

op <- graphics::par(no.readonly = TRUE)
plot(d, theta)

graphics::par(op)

# Two coordinates of a heavy-tailed family, whose contours are wider than a
# gaussian's at the same matrix.
t2 <- mvstudent_t1_distrib(3)
th2 <- as.list(stats::setNames(c(unlist(theta), 3), t2@params))
op <- graphics::par(no.readonly = TRUE)
plot(t2, th2, which = c(1, 2))

graphics::par(op)

# Four coordinates would be sixteen panels, so it is refused by name.
d4 <- mvgaussian1_distrib(4)
th4 <- as.list(stats::setNames(rep(0, d4@n_params), d4@params))
try(plot(d4, th4))
#> Error : A panel matrix of 4 coordinates has 16 panels and is not readable.
#>   Choose at most three with 'which', for instance which = c(1, 2, 3).
```
