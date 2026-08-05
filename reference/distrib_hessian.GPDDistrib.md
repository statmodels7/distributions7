# Generalized Pareto Analytical Observed Hessian

The second derivatives of the same expressions, kept short by \\t - \xi
z = 1\\, which makes \\\partial u/\partial\sigma\\ equal to \\-z/(\sigma
t^2)\\. The pure-\\\xi\\ component goes through a series near zero,
where its two singular terms cancel.

## Arguments

- distrib:

  A `GPDDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `sigma` and `xi`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of second-derivative components.

## See also

[`gpd_distrib`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
