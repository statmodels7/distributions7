# Exponential Mixed Derivatives

Closed form: the mean is a pure scale, so the identity \\-z\ell^{(yy)} -
\ell^{(y)}/\mu\\ applies with \\z = y/\mu\\ and returns \\1/\mu^2\\.

## Arguments

- distrib:

  An `ExponentialDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with the single component `mu`, a numeric vector of length
`length(y)` holding \\1/\mu^2\\ at every observation.

## See also

[`loc_scale_cross_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_block.md)
for the identity, whose scale half applies here with no location;
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
for the family;
[`distrib_cross_y.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Weibull1Distrib.md),
which adds a shape to it.
