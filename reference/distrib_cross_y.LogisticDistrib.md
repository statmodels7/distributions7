# Logistic Mixed Derivatives

Closed form from the location-scale identity: the location component is
\\-\ell^{(yy)}\\ and the scale one \\-z\ell^{(yy)} -
\ell^{(y)}/\sigma\\.

## Arguments

- distrib:

  A `LogisticDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `sigma`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with components `mu` and `sigma`, each a numeric vector of
length `length(y)`.

## See also

[`loc_scale_cross_y()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_y.md),
the shared body, and
[`loc_scale_cross_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_block.md)
for the identity;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
