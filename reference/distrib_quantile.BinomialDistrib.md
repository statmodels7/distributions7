# Binomial Quantile Function

Computes the quantile function for the Binomial distribution, the
generalized inverse of the CDF: \$\$Q(p; \mu, n) = \min\left\\y \in \\0,
1, \dots, n\\ : F(y; \mu, n) \ge p\right\\\$\$

## Arguments

- distrib:

  A `BinomialDistrib` object.

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

[`binomial_distrib`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
