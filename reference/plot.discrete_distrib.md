# Plot Method for Discrete Distributions

Draws the probability mass of a discrete family as stems over the
integers in range, one set per parameter setting. A component of `theta`
given as a vector asks for one setting per element, exactly as in
[`plot.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.continuous_distrib.md).

## Arguments

- x:

  An object inheriting from `discrete_distrib`.

- theta:

  A named list or numeric vector holding one entry per parameter of `x`.
  A component of length \\k \> 1\\ draws \\k\\ sets of stems, and every
  component must have length 1 or that same \\k\\. Missing, a random
  parameter is drawn by
  [`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md).

- xlim:

  A numeric vector of length 2, or `NULL` (the default) to compute the
  range from the quantiles over every setting.

- legend:

  Logical of length 1: draw the key when several settings are plotted.
  Defaults to `TRUE`.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html),
  for instance `main`, `xlab` or `ylim`. A `col` or `pch` given here
  wins over the palette and is recycled over the settings; an `lty` is
  accepted and has no effect on the stems.

## Value

`x`, invisibly. Called for the plot it draws.

## The two keys, and why they are not the continuous ones

Settings are separated by **color and symbol together**, the symbol
sitting at the top of each stem where the point already is. The line
type is not used here: dashing a stem is what it would do, and a dashed
stem reads as a broken one, so at a support of any size the panel fills
with fragments that cross each other.

Several settings are also shifted sideways by a fraction of the spacing,
so that equal masses at one support point stay countable instead of
standing one in front of another.

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
# Three rates, whose masses are countable at every support point because
# the settings are offset sideways.
plot(poisson_distrib(), list(mu = c(1, 4, 10)))


# Overdispersion at a fixed mean: the mass spreads as theta falls.
plot(negbin2_distrib(), list(mu = 5, theta = c(0.5, 2, 50)))

```
