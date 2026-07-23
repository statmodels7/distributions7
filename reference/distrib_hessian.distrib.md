# Default Numerical Hessian for `distrib` Objects

Fallback method: distributions that do not implement an analytical
Hessian get one computed by finite differences of the log-density (see
[`numerical_hessian`](https://statmodels7.github.io/distributions7/reference/numerical_hessian.md)).

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

## Value

A named list of Hessian component vectors in
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
order.
