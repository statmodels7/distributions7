# Beta Analytical Expected Hessian in Its Shapes

Equal to the observed Hessian: it is free of the data, so there is
nothing to average. Fisher scoring and Newton's method therefore take
the same step on the parameter scale here.

## Arguments

- distrib:

  A `Beta2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `alpha` and `beta`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is exact.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list of expected second derivatives.

## See also

[`beta2_distrib`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md)
