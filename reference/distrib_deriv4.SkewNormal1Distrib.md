# Skew Normal Analytical Fourth-Order Derivatives

Closed-form fourth-order derivatives of the skew normal log-density, in
the notation of
[`distrib_deriv3.SkewNormal1Distrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewNormal1Distrib.md).
The expected derivatives are approximated numerically, as at third
order.

## Arguments

- distrib:

  A `SkewNormal1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `sigma` and `alpha`.

- expected:

  Logical; if `TRUE`, the expectation is approximated numerically.

- approx:

  Strategy for the expectation; see
  [`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md).

- nsim:

  Monte Carlo sample size when `approx = "mc"`.

## Value

A named list of fourth-derivative component vectors.

## See also

[`skewnormal1_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
