# Beta Mixed Derivatives in the Shapes

Closed form: \\\ell^{(y)} = (\alpha-1)/y - (\beta-1)/(1-y)\\ gives
\\1/y\\ and \\-1/(1-y)\\.

## Arguments

- distrib:

  A `Beta2Distrib` object.

- y:

  A numeric vector in the unit interval.

- theta:

  A list containing `alpha` and `beta`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with components `alpha` and `beta`, each a numeric vector
of length `length(y)`.

## See also

[`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md)
for the family;
[`distrib_cross_y.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Beta1Distrib.md)
for the mean chart;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
