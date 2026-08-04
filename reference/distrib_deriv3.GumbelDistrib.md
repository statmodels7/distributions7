# Gumbel Analytical Third-Order Derivatives

Closed-form third-order derivatives of the Gumbel log-density (observed,
or expected when `expected = TRUE`). With \\z = (y-\mu)/\sigma\\ and \\w
= e^{-z}\\, every derivative is a polynomial in \\z\\ and \\z^j w\\; the
expected values use \\E\[z^k w\] = (-1)^k \Gamma^{(k)}(2)\\ and \\E\[z\]
= \gamma\\, since \\w\\ is standard exponential under the model.

## Arguments

- distrib:

  A `GumbelDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- expected:

  Logical; if `TRUE`, returns the expected third derivatives.

## Value

A named list of third-derivative component vectors.

## See also

[`gumbel_distrib`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
