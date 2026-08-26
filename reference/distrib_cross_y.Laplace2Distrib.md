# Laplace Mixed Derivatives in Location and Rate

Closed form, away from the kink at \\y = \mu\\: with \\\ell^{(y)} =
-\lambda\\\mathrm{sign}(y-\mu)\\, the location component is 0 almost
everywhere and the rate component is \\-\mathrm{sign}(y-\mu)\\.

## Arguments

- distrib:

  A `Laplace2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `lambda`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with components `mu` and `lambda`, each a numeric vector of
length `length(y)`. The `mu` component is exactly zero almost everywhere
and the `lambda` component is \\\pm 1\\.

## See also

[`distrib_cross_y.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.LaplaceDistrib.md)
for the scale chart;
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
for the family;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
