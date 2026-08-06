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

A named list with one numeric vector per parameter.
