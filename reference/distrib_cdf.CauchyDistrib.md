# Cauchy Cumulative Distribution Function

Computes the cumulative distribution function for the Cauchy
distribution: \$\$F(q; \mu, \sigma) = \dfrac{1}{2} +
\dfrac{1}{\pi}\arctan\left(\dfrac{q-\mu}{\sigma}\right)\$\$

## Arguments

- distrib:

  A `CauchyDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameters `mu` and `sigma`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of cumulative probabilities.

## See also

[`cauchy_distrib`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
