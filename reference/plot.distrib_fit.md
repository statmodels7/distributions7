# Plot a Fitted Distribution Against the Data

[`autoplot()`](https://statmodels7.github.io/distributions7/reference/autoplot.md)
returns a graph that compares the fitted distribution with the sample
from which it was estimated, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) prints that
graph. The comparison is chosen from the family's support. A continuous
family is drawn as the fitted density over a kernel estimate of the
data, with a rug of the observations underneath. A discrete one is drawn
as the observed relative frequencies in bars with the fitted
probabilities overlaid, because a kernel estimate misrepresents a sample
from a discrete distribution. A legend below the panel names each
element.

A multivariate fit is drawn as a panel matrix: the fitted marginal
density and a kernel estimate of the data on the diagonal, the fitted
contours over the observations below it, the fitted correlation above.

## Usage

``` r
# S3 method for class 'distrib_fit'
plot(x, n_grid = 512, rug = NULL, legend = TRUE, mv_which = NULL, ...)

autoplot_distrib_fit(
  object,
  n_grid = 512,
  rug = NULL,
  legend = TRUE,
  mv_which = NULL,
  ...
)
```

## Arguments

- x, object:

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
  becomes a solid band that conveys no information.

- legend:

  Logical of length 1: draw the legend below the panel. Defaults to
  `TRUE`.

- mv_which:

  For a multivariate fit, which coordinates to show, by index. Defaults
  to all of them. At most three coordinates can be drawn, because a
  larger panel matrix is not readable, so a fit with more than three
  coordinates signals an error unless `mv_which` selects at most three.
  The argument is ignored for a univariate fit.

- ...:

  Unused; an argument passed here signals an error.

## Value

[`autoplot()`](https://statmodels7.github.io/distributions7/reference/autoplot.md)
returns a ggplot2 object (a patchwork object for a multivariate fit);
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) prints the same
object and returns it invisibly.

## Details

The graph is drawn with ggplot2, which is suggested, and the panel
matrix of a multivariate fit is joined with patchwork, also suggested;
the methods signal an error that names them when they are not installed.
The title, the labels of the axes and the colours are changed on the
returned object.

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
p <- plot(fp)

p + ggplot2::labs(title = "Poisson fit")
```
