# Plot Method for Continuous Distributions

Draws the density of a continuous family over a grid of 1000 points, one
curve per parameter setting. A component of `theta` given as a vector
gives one curve per element, so
`plot(d, list(mu = 0, sigma = c(1, 2, 4)))` draws three densities that
share a mean and differ in scale.

## Arguments

- x:

  An object inheriting from `continuous_distrib`.

- theta:

  A named list or numeric vector holding one entry per parameter of `x`.
  A component of length \\k \> 1\\ draws \\k\\ curves, and every
  component must have length 1 or that same \\k\\. When `theta` is
  missing, random values for all parameters are drawn by
  [`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md).

- xlim:

  A numeric vector of length 2, or `NULL` (the default) to compute the
  range from the quantiles as above.

- legend:

  Logical of length 1: draw the key when several settings are plotted.
  Defaults to `TRUE`. A single setting is drawn without a key.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html),
  for instance `main`, `xlab` or `ylim`. A `col` or `lty` given here
  wins over the palette and is recycled over the curves.

## Value

`x`, invisibly. Called for the plot it draws.

## The two keys

Curves are separated by **color and line type together**, so that a
printed copy with no color is still readable. The parameters that vary
are named in a legend and those held fixed are stated once in the title,
so a reader sees what the panel is about without counting curves. The
legend goes in whichever top corner the mass leaves emptier. See
[`plot_settings()`](https://statmodels7.github.io/distributions7/reference/plot_settings.md)
for the rule on lengths and
[`plot_keys()`](https://statmodels7.github.io/distributions7/reference/plot_keys.md)
for the palette.

## The horizontal range

It covers every setting, not the first: from the smallest 0.5% quantile
to the largest 99.5% quantile over them, widened by a tenth and clamped
to the support. Give `xlim` to override it, which is useful for a
heavy-tailed family, whose 99.5% quantile can sit orders of magnitude
past anything worth drawing.

## See also

[`plot.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.discrete_distrib.md)
for a discrete family;
[`plot.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md)
to compare a fit with data;
[`plot_settings()`](https://statmodels7.github.io/distributions7/reference/plot_settings.md),
[`plot_keys()`](https://statmodels7.github.io/distributions7/reference/plot_keys.md)
and
[`plot_labels()`](https://statmodels7.github.io/distributions7/reference/plot_labels.md)
for the pieces.

## Examples

``` r
# One curve per value of the scale, at a shared mean.
plot(gaussian1_distrib(), list(mu = 0, sigma = c(1, 2, 4)))


# A heavy-tailed family needs its own window: the 99.5% quantile of a
# standard Cauchy is about ten times as far out as the window drawn here.
plot(cauchy_distrib(), list(mu = 0, sigma = 1), xlim = c(-6, 6))

```
