# Skew Normal Analytical Third-Order Derivatives

Closed-form third-order derivatives of the skew normal log-density. With
\\t = \alpha z\\ and \\R\\ the inverse Mills ratio, the derivatives of
\\\log \Phi(t)\\ follow from \\R' = -R(t+R)\\ and stay polynomials in
\\t\\ and \\R\\, so every component is elementary. The expected
derivatives have no closed form (the same integrals as the expected
information) and come from
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md).

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
  [`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md).

- nsim:

  Monte Carlo sample size when `approx = "mc"`.

## Value

A named list of third-derivative component vectors.

## See also

[`skewnormal1_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
