# Default Fourth-Order Derivatives for `distrib` Objects

Fallback method: observed fourth derivatives via
[`numerical_deriv4`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md)
(second differences of the Hessian); expected fourth derivatives via the
[`expectation`](https://statmodels7.github.io/distributions7/reference/expectation.md)
of the observed ones.

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- expected:

  Logical; if `TRUE`, returns expected derivatives.

## Value

A named list of fourth-derivative component vectors.
