# Pseudo-Huber Analytical Third-Order Derivatives

Closed-form observed third-order derivatives of the Pseudo-Huber
log-density. Bessel functions enter only through the pure-\\\nu\\
component; the exponentially scaled forms are used so that large \\\nu\\
does not overflow. The expected third derivatives have no closed form,
so `expected = TRUE` is handled by the strategies in
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

A named list of third-derivative component vectors.

## See also

[`pseudohuber_distrib`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
