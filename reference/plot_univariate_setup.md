# The Settings of a Univariate Graph

Draws random parameters when `theta` is missing, reporting them in a
message, and returns the settings, the legend entries and the title of
the graph. It is the first step of the
[`autoplot()`](https://statmodels7.github.io/distributions7/reference/autoplot.md)
methods of
[`plot.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.continuous_distrib.md)
and
[`plot.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.discrete_distrib.md).

## Usage

``` r
plot_univariate_setup(x, theta)
```

## Arguments

- x:

  A univariate distribution object.

- theta:

  The caller's `theta`, possibly missing.

## Value

A list with the components `ps` (the value of
[`plot_settings()`](https://statmodels7.github.io/distributions7/reference/plot_settings.md))
and `labs` (the value of
[`plot_labels()`](https://statmodels7.github.io/distributions7/reference/plot_labels.md)).
