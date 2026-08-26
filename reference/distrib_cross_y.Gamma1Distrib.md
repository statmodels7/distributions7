# Gamma Mixed Derivatives in Mean and Dispersion

Closed form. The shape is \\1/\phi\\ and the rate \\1/(\phi\mu)\\, so
the mean component is \\1/(\phi\mu^2)\\ and the dispersion one
\\1/(\phi^2\mu) - 1/(\phi^2 y)\\.

## Arguments

- distrib:

  A `Gamma1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `phi`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with components `mu` and `sigma2`, each a numeric vector of
length `length(y)`.

## See also

[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
for the family;
[`distrib_cross_y.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Gamma2Distrib.md)
for the dispersion chart;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
