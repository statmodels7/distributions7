# Cauchy Quantile Function

Computes the quantile function (inverse CDF) for the Cauchy
distribution: \$\$Q(p; \mu, \sigma) = \mu + \sigma
\tan\left(\pi\left(p - \dfrac{1}{2}\right)\right)\$\$

## Arguments

- distrib:

  A `CauchyDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing the parameters `mu` and `sigma`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\,
  otherwise \\P(Y \> p)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of quantiles.

## See also

[`cauchy_distrib`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
