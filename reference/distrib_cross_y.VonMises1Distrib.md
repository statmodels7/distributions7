# von Mises Mixed Derivatives

Closed form. With \\\ell^{(y)} = -\kappa\sin(y-\mu)\\, the location
component is \\\kappa\cos(y-\mu)\\ and the concentration one
\\-\sin(y-\mu)\\.

## Arguments

- distrib:

  A `VonMises1Distrib` object.

- y:

  A numeric vector of angles.

- theta:

  A list containing `mu` and `kappa`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with components `mu` and `kappa`, each a numeric vector of
length `length(y)`.

## See also

[`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md)
for the family;
[`distrib_cross_y.VonMises2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.VonMises2Distrib.md)
for the resultant-length chart;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
