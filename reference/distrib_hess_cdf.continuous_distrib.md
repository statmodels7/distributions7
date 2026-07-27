# Default Log-CDF Hessian for Continuous Distributions

Fallback: finite differences of the distribution function.

## Arguments

- distrib:

  A `continuous_distrib` object.

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
