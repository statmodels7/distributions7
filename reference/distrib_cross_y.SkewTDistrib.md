# Skew t Mixed Derivatives

Closed form in the location and the scale, from the location-scale
identity. The two shape components come from one central difference of
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md):
the density carries \\T\_{\nu+1}\\, whose derivative in the degrees of
freedom has no elementary form, the same obstruction that the family's
parameter derivatives meet.

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
numeric vector of length `length(y)`. The first two are exact and the
last two carry one stencil's error, measured at \\9.2\times10^{-10}\\
relative against Richardson.

## See also

[`partial_loc_scale_cross_y()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_cross_y.md),
the shared body;
[`distrib_cross_y.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.SkewNormal1Distrib.md),
where every component is exact;
[`numerical_cross_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross_y.md)
for the differenced half.
