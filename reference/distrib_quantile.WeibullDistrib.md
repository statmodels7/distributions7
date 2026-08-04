# Weibull Quantile Function

Computes the quantile function for the Weibull distribution: \$\$Q(p;
\mu, \sigma) = \mu \left\\-\log(1 - p)\right\\^{1/\sigma}\$\$

## Arguments

- distrib:

  A `WeibullDistrib` object.

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

[`weibull_distrib`](https://statmodels7.github.io/distributions7/reference/weibull_distrib.md)
