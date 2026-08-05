# Lognormal Quantile Function

Computes the quantile function (inverse CDF) for the Lognormal
distribution: \$\$Q(p; \mu, \sigma^2) = \exp\left(\mu +
\sqrt{\sigma^2}\\\Phi^{-1}(p)\right)\$\$ where \\\Phi^{-1}\\ is the
standard normal quantile function.

## Arguments

- distrib:

  A `Lognormal1Distrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing the parameters `mu` and `sigma2`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\,
  otherwise \\P(Y \> p)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of quantiles.

## See also

[`lognormal1_distrib`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
