# Gaussian Mixed Derivatives in Mean and Precision

Closed form: \\\ell^{(y)} = -\tau(y-\mu)\\ gives \\\tau\\ and
\\-(y-\mu)\\.

## Arguments

- distrib:

  A `Gaussian3Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `tau`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with components `mu` and `tau`, each a numeric vector of
length `length(y)`.

## See also

[`distrib_cross_y.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Gaussian1Distrib.md)
for the scale chart;
[`distrib_cross_y.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Gaussian2Distrib.md)
for the variance chart;
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md)
for the family.
