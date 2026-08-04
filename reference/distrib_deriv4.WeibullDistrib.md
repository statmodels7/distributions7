# Weibull Analytical Fourth-Order Derivatives

Closed-form fourth-order derivatives of the Weibull log-density
(observed, or expected when `expected = TRUE`), in the notation of
[`distrib_deriv3.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.WeibullDistrib.md).

## Arguments

- distrib:

  A `WeibullDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- expected:

  Logical; if `TRUE`, returns the expected fourth derivatives.

## Value

A named list of fourth-derivative component vectors.

## See also

[`weibull_distrib`](https://statmodels7.github.io/distributions7/reference/weibull_distrib.md)
