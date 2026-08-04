# Weibull Analytical Third-Order Derivatives

Closed-form third-order derivatives of the Weibull log-density
(observed, or expected when `expected = TRUE`). With \\u =
(y/\mu)^{\sigma}\\ and \\L = \log(y/\mu)\\, every derivative is a
polynomial in \\u\\ and \\L u\\; the expected values use \\E\[u L^k\] =
\Gamma^{(k)}(2)/\sigma^k\\.

## Arguments

- distrib:

  A `WeibullDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- expected:

  Logical; if `TRUE`, returns the expected third derivatives.

## Value

A named list of third-derivative component vectors.

## See also

[`weibull_distrib`](https://statmodels7.github.io/distributions7/reference/weibull_distrib.md)
