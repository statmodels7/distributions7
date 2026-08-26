# Weibull Mixed Derivatives

Closed form at both parameters. The scale follows the pure-scale
identity and reduces to \\\sigma^2 w/(\mu y)\\ with \\w =
(y/\mu)^\sigma\\; the shape comes from differentiating \\\ell^{(y)} =
((\sigma-1) - \sigma w)/y\\, giving \\(1 - w - \sigma w L)/y\\ with \\L
= \log(y/\mu)\\.

## Arguments

- distrib:

  A `Weibull1Distrib` object.

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

[`distrib_cross_y.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.ExponentialDistrib.md),
the case \\\sigma = 1\\;
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
for the family;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
