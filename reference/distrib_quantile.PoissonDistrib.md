# Poisson Quantile Function

Computes the quantile function for the Poisson distribution, defined as
the generalized inverse of the CDF: \$\$Q(p; \mu) = \min\left\\y \in
\mathbb{N}\_0 : F(y; \mu) \ge p\right\\\$\$

## Arguments

- distrib:

  A `PoissonDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing the parameter `mu`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\,
  otherwise \\P(Y \> p)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of quantiles.

## See also

[`poisson_distrib`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
