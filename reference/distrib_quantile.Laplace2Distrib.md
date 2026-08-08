# Laplace Quantile Function in Location and Rate

\$\$Q(p; \mu, \lambda) = \mu - \dfrac{1}{\lambda}\\\mathrm{sign}(p -
\tfrac{1}{2})\\\log\left(1 - 2\left\|p - \tfrac{1}{2}\right\|\right)\$\$

## Arguments

- distrib:

  A `Laplace2Distrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing the parameters `mu` and `lambda`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\,
  otherwise \\P(Y \> p)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of quantiles.

## See also

[`laplace2_distrib`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
