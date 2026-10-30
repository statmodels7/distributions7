# Plot Method for Discrete Distributions

[`autoplot()`](https://statmodels7.github.io/distributions7/reference/autoplot.md)
returns a graph of the probability mass of a discrete family as stems
over the integers in range, one set per parameter setting, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) prints that
graph. A component of `theta` given as a vector gives one setting per
element, exactly as in
[`plot.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.continuous_distrib.md).

## Usage

``` r
# S3 method for class 'discrete_distrib'
plot(x, theta, xlim = NULL, legend = TRUE, ...)

autoplot_discrete_distrib(object, theta, xlim = NULL, legend = TRUE, ...)
```

## Arguments

- x, object:

  An object inheriting from `discrete_distrib`.

- theta:

  A named list or numeric vector holding one entry per parameter of `x`.
  A component of length \\k \> 1\\ draws \\k\\ sets of stems, and every
  component must have length 1 or that same \\k\\. When `theta` is
  missing, random values for all parameters are drawn by
  [`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md).

- xlim:

  A numeric vector of length 2, or `NULL` (the default) to compute the
  range from the quantiles over every setting.

- legend:

  Logical of length 1: draw the legend below the panel when several
  settings are plotted. Defaults to `TRUE`.

- ...:

  Unused; an argument passed here signals an error.

## Value

[`autoplot()`](https://statmodels7.github.io/distributions7/reference/autoplot.md)
returns a ggplot2 object;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) prints the same
object and returns it invisibly.

## The two keys, and why they are not the continuous ones

Settings are separated by color and symbol together, and the symbol sits
at the top of each stem, where the point already is. The line type is
not used here, because it would dash the stems and a dashed stem reads
as a broken one, so at a support of any size the panel fills with
fragments that cross each other. The parameters that vary are named in a
legend below the panel and those held fixed are stated once in the
title.

Several settings are also shifted sideways by a fraction of the spacing,
so that equal masses at one support point stay countable instead of
standing one in front of another.

The graph is drawn with ggplot2, which is suggested; the methods signal
an error that names it when it is not installed.

## See also

[`plot.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.continuous_distrib.md)
for a density;
[`plot.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md)
to compare a fit with data;
[`plot_settings()`](https://statmodels7.github.io/distributions7/reference/plot_settings.md),
[`plot_keys()`](https://statmodels7.github.io/distributions7/reference/plot_keys.md)
and
[`plot_labels()`](https://statmodels7.github.io/distributions7/reference/plot_labels.md)
for the pieces.

## Examples

``` r
# The three rates have masses that are countable at every support point
# because the settings are offset sideways.
plot(poisson_distrib(), list(mu = c(1, 4, 10)))


# Overdispersion at a fixed mean: the mass spreads as theta falls.
plot(negbin2_distrib(), list(mu = 5, theta = c(0.5, 2, 50)))
```
