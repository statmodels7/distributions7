# Laplace Analytical Fourth-Order Derivatives

Closed-form fourth-order derivatives of the Laplace log-density, almost
everywhere (observed, or expected when `expected = TRUE`), in the
notation of
[`distrib_deriv3.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.LaplaceDistrib.md):
the non-zero components are \\\ell^{(\mu bbb)} = -6s/b^4\\ and
\\\ell^{(bbbb)} = 6/b^4 - 24a/b^5\\.

## Arguments

- distrib:

  A `LaplaceDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `b`.

- expected:

  Logical; if `TRUE`, returns the expected fourth derivatives.

## Value

A named list of fourth-derivative component vectors.

## See also

[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
