# Laplace Quantile Function

Computes the quantile function (inverse CDF) for the Laplace
distribution: \$\$Q(p; \mu, b) = \mu - b\\\mathrm{sign}(p -
\tfrac{1}{2})\\\log\left(1 - 2\left\|p - \tfrac{1}{2}\right\|\right)\$\$

## Arguments

- distrib:

  A `LaplaceDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing the parameters `mu` and `b`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\,
  otherwise \\P(Y \> p)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of quantiles.

## See also

[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
