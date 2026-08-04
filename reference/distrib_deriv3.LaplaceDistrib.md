# Laplace Analytical Third-Order Derivatives

Closed-form third-order derivatives of the Laplace log-density, almost
everywhere (observed, or expected when `expected = TRUE`). With \\s =
\mathrm{sign}(y-\mu)\\ and \\a = \lvert y-\mu \rvert\\, the only
non-zero components are \\\ell^{(\mu b b)} = 2s/b^3\\ and \\\ell^{(bbb)}
= -2/b^3 + 6a/b^4\\; the kink at \\y = \mu\\ is the same one the
gradient carries, and `params_smooth` records it.

## Arguments

- distrib:

  A `LaplaceDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `b`.

- expected:

  Logical; if `TRUE`, returns the expected third derivatives.

## Value

A named list of third-derivative component vectors.

## See also

[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
