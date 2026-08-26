# Generalized Pareto Mixed Derivatives

Closed form at both parameters, from \\\ell^{(y)} = -(\xi+1)/(\sigma
t)\\ with \\t = 1 + \xi y/\sigma\\: the scale gives \\(\xi+1)(t - \xi
z)/(\sigma^2 t^2)\\ and the shape \\-1/(\sigma t) + (\xi+1)y/(\sigma^2
t^2)\\. Neither carries a \\1/\xi\\, so the shape direction needs no
series at the exponential limit.

## Arguments

- distrib:

  A `GPDDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `sigma` and `xi`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with components `sigma` and `xi`, each a numeric vector of
length `length(y)`. Measured against Richardson on the analytic response
gradient the worst is \\4.9\times10^{-11}\\ relative.

## See also

[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
for the family and for the removable singularity at \\\xi = 0\\;
[`distrib_grad_y.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GPDDistrib.md)
for the quantity differentiated;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
