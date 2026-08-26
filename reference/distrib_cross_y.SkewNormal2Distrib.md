# Skew Normal Mixed Derivatives in the Centered Parametrization

The direct parametrization's mixed block carried by the first-order
chain rule on the centered-to-direct map.

## Arguments

- distrib:

  A `SkewNormal2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `sigma` and `gamma1`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter of the centered
parametrization, keyed by `distrib@params`, each of length `length(y)`.

## See also

[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
for the family;
[`distrib_cross_y.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.SkewNormal1Distrib.md)
for the direct chart;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
