# Log-CDF Hessian for Discrete Distributions

Exact, by the same finite sum as the gradient.

## Arguments

- distrib:

  A `discrete_distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list keyed as
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md).
