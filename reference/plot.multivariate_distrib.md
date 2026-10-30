# Panels of a Multivariate Density

[`autoplot()`](https://statmodels7.github.io/distributions7/reference/autoplot.md)
returns a graph of a multivariate distribution as a matrix of panels,
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) prints that
graph: the marginal density of each coordinate on the diagonal, contours
of each bivariate marginal below it, and the implied correlation printed
above it. The density panels are built from marginals, so the graph
exists where
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
returns a marginal of several coordinates, which holds for the
multivariate gaussian and Student \\t\\. For the Dirichlet and the
multinomial, whose marginals are returned one coordinate at a time, the
graph signals the error of
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md).

## Usage

``` r
# S3 method for class 'multivariate_distrib'
plot(x, theta, which = NULL, n_grid = 80, legend = TRUE, ...)

autoplot_multivariate_distrib(
  object,
  theta,
  which = NULL,
  n_grid = 80,
  legend = TRUE,
  ...
)
```

## Arguments

- x, object:

  An object inheriting from
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md).

- theta:

  A named list or vector of parameters in which each component is a
  single number. When missing, parameters are drawn by
  [`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md)
  and reported in a message. A component of length greater than one is
  an error.

- which:

  An integer vector of the coordinates to show, at most three. Defaults
  to `NULL`, which means all of them and is itself an error above three
  coordinates.

- n_grid:

  Points per axis for the contours. The marginal density curves on the
  diagonal use `4 * n_grid` points. Defaults to `80`. The bivariate
  panels evaluate the density on an `n_grid * n_grid` grid, so the cost
  is quadratic in it.

- legend:

  Logical of length 1: draw the legend below the panel matrix. Defaults
  to `TRUE`.

- ...:

  Unused; an argument passed here signals an error.

## Value

[`autoplot()`](https://statmodels7.github.io/distributions7/reference/autoplot.md)
returns a patchwork object;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) prints the same
object and returns it invisibly.

## At most three coordinates

A panel matrix of \\k\\ coordinates has \\k^2\\ panels, and past three
the panels are too small to read. A larger distribution signals an error
that names the number of coordinates and suggests `which`.

## The contour levels

The contours are levels of the density and not of the probability
content of the enclosed region, which would require the mass of each
level set of the density. The levels are the values that
[`pretty()`](https://rdrr.io/r/base/pretty.html) selects, about six,
over the range of the density on the grid.

## One setting only

A univariate [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
method reads a `theta` component of length \\k\\ as \\k\\ settings drawn
over one another. Here the picture is already a matrix of panels and has
no axis on which to overlay them, so a component longer than one signals
an error that names it.

## The plotting range

Each coordinate is drawn over three standard deviations either side of
the location, and the spread is taken from the matrix that the
parametrization carries (see
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)).
A multivariate Student \\t\\ with \\\nu \le 2\\ has an infinite
variance, so
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
cannot be the source. For an elliptical family the two matrices agree up
to a factor, which is enough for a plotting range.

The panels are drawn with ggplot2 and joined with patchwork; both
packages are suggested, and the methods signal an error that names them
when they are not installed.

## See also

[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md),
which supplies the density panels,
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
plot(d, theta)


# Two coordinates of a heavy-tailed family are drawn, and its contours are
# wider than a gaussian's at the same matrix.
t2 <- mvstudent_t1_distrib(3)
th2 <- as.list(stats::setNames(c(unlist(theta), 3), t2@params))
plot(t2, th2, which = c(1, 2))


# Four coordinates would give sixteen panels, so the call signals an error.
d4 <- mvgaussian1_distrib(4)
th4 <- as.list(stats::setNames(rep(0, d4@n_params), d4@params))
try(plot(d4, th4))
#> Error : A panel matrix of 4 coordinates has 16 panels and is not readable.
#>   Choose at most three with 'which', for instance which = c(1, 2, 3).
```
