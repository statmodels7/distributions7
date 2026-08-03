# No Distribution Function in Several Dimensions

Refused. The distribution function of a multivariate law is an integral
over an orthant, which has no closed form for the gaussian and no
one-dimensional fallback to stand in for it, and the quadrature
registered on
[`continuous_distrib`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)
integrates over an interval.

## Arguments

- distrib:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- q:

  A numeric matrix of quantiles.

- theta:

  A named list of parameters.

- ...:

  Unused.

## Value

Never returns; raises an error.
