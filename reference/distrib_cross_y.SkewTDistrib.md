# Skew t Mixed Derivatives

Closed form in the location and the scale, from the location-scale
identity. The two shape components are minus the mixed entries of
[`distrib_hessian.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewTDistrib.md),
since the response and the location enter the density only through their
difference, so that \\\partial_y\partial\_\theta\ell =
-\partial\_\mu\partial\_\theta\ell\\.

## Arguments

- distrib:

  A `SkewTDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `sigma`, `alpha` and `nu`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with components `mu`, `sigma`, `alpha` and `nu`, each a
numeric vector of length `length(y)`.

## See also

[`loc_scale_cross_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_block.md)
for the location-scale identity and
[`distrib_hessian.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewTDistrib.md)
for the mixed entries.
