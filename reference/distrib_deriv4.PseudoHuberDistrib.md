# Pseudo-Huber Analytical Fourth-Order Derivatives

Closed-form observed fourth-order derivatives of the Pseudo-Huber
log-density (see
[`distrib_deriv3.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.PseudoHuberDistrib.md)
for the Bessel handling). The expected fourth derivatives have no closed
form and are handled by the strategies in
[`expected_derivative_methods`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

## Arguments

- distrib:

  A `PseudoHuberDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu`, `sigma` and `nu`.

- expected:

  Logical; if `TRUE`, returns the (approximated) expected derivatives.

- approx, nsim:

  Passed to
  [`expected_derivative_methods`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
  when `expected = TRUE`.

## Value

A named list of fourth-derivative component vectors.

## See also

[`pseudohuber_distrib`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
