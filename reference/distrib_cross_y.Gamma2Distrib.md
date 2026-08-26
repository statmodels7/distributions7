# Gamma Mixed Derivatives in Mean and Variance

Closed form: \\\ell^{(y)} = (\mu^2/\sigma^2 - 1)/y - \mu/\sigma^2\\
gives \\2\mu/(\sigma^2 y) - 1/\sigma^2\\ and \\-\mu^2/(\sigma^4 y) +
\mu/\sigma^4\\.

## Arguments

- distrib:

  A `Gamma2Distrib` object.

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

[`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md)
for the family;
[`distrib_cross_y.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Gamma1Distrib.md)
for the variance chart;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
