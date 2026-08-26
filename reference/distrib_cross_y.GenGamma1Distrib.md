# Generalized Gamma Mixed Derivatives

Closed form. With \\w = (y/a)^p\\ and \\L = \log(y/a)\\, \\\ell^{(y)} =
((d-1) - pw)/y\\, so the components are \\p^2w/(ay)\\, \\1/y\\ and
\\-w(1 + pL)/y\\.

## Arguments

- distrib:

  A `GenGamma1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `a`, `d` and `p`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with components `a`, `d` and `p`, each a numeric vector of
length `length(y)`. Measured against Richardson on the analytic response
gradient the worst is \\2.8\times10^{-11}\\ relative.

## See also

[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
for the family and its parametrization;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
