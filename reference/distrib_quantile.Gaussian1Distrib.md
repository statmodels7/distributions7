# Gaussian Quantile Function

Computes the quantile function (inverse CDF) for the Gaussian
distribution: \$\$Q(p; \mu, \sigma) = \mu + \sigma \Phi^{-1}(p)\$\$

## Arguments

- distrib:

  A `Gaussian1Distrib` object.

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

[`gaussian1_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
