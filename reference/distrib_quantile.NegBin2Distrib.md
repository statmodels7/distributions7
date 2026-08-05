# Negative Binomial Quantile Function

Computes the quantile function for the Negative Binomial distribution,
the generalized inverse of the CDF: \$\$Q(p; \mu, \theta) = \min\left\\y
\in \mathbb{N}\_0 : F(y; \mu, \theta) \ge p\right\\\$\$

## Arguments

- distrib:

  A `NegBin2Distrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing the parameters `mu` and `theta`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\,
  otherwise \\P(Y \> p)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of quantiles.

## See also

[`negbin2_distrib`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
