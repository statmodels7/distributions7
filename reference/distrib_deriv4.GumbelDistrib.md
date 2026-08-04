# Gumbel Analytical Fourth-Order Derivatives

Closed-form fourth-order derivatives of the Gumbel log-density
(observed, or expected when `expected = TRUE`), in the notation of
[`distrib_deriv3.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GumbelDistrib.md).

## Arguments

- distrib:

  A `GumbelDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- expected:

  Logical; if `TRUE`, returns the expected fourth derivatives.

## Value

A named list of fourth-derivative component vectors.

## See also

[`gumbel_distrib`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
