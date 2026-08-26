# Gaussian Mixed Derivatives in Mean and Variance

Closed form: \\\ell^{(y)} = -(y-\mu)/\sigma^2\\ gives \\1/\sigma^2\\ and
\\(y-\mu)/\sigma^4\\.

## Arguments

- distrib:

  A `Gaussian2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `sigma2`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with components `mu` and `sigma2`, each a numeric vector of
length `length(y)`.

## See also

[`distrib_cross_y.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Gaussian1Distrib.md)
for the scale chart;
[`distrib_cross_y.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Gaussian3Distrib.md)
for the precision chart;
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
for the family.
