# Poisson Cumulative Distribution Function

Computes the cumulative distribution function for the Poisson
distribution: \$\$F(q; \mu) = \sum\_{k=0}^{\lfloor q \rfloor}
\dfrac{\mu^k e^{-\mu}}{k!}\$\$

## Arguments

- distrib:

  A `PoissonDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameter `mu`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of cumulative probabilities.

## See also

[`poisson_distrib`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
