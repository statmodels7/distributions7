# Lognormal Mixed Derivatives

Closed form. On the log scale the family is location-scale, and
\\\ell^{(y)} = -1/y - (\log y - \mu)/(\sigma^2 y)\\ gives \\1/(\sigma^2
y)\\ and \\(\log y - \mu)/(\sigma^4 y)\\.

## Arguments

- distrib:

  A `Lognormal1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `sigma2`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with components `mu` and `sigma`, each a numeric vector of
length `length(y)`.

## See also

[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
for the family;
[`distrib_cross_y.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Gaussian1Distrib.md),
the law on the log scale;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
