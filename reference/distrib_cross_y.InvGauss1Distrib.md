# Inverse Gaussian Mixed Derivatives

Closed form. With \\\ell^{(y)} = -3/(2y) - 1/(2\phi\mu^2) + 1/(2\phi
y^2)\\, the components are \\1/(\phi\mu^3)\\ and \\1/(2\phi^2\mu^2) -
1/(2\phi^2 y^2)\\.

## Arguments

- distrib:

  An `InvGauss1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `phi`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with components `mu` and `sigma2`, each a numeric vector of
length `length(y)`.

## See also

[`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md)
for the family;
[`distrib_cross_y.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.InvGauss2Distrib.md)
for the shape chart;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
