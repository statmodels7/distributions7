# Weibull Cumulative Distribution Function

Computes the cumulative distribution function for the Weibull
distribution: \$\$F(q; \mu, \sigma) = 1 -
\exp\left\\-(q/\mu)^{\sigma}\right\\\$\$ The survival function is exact
on the log scale, which is what a censored observation in the far tail
needs.

## Arguments

- distrib:

  A `WeibullDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameters `mu` and `sigma`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of cumulative probabilities.

## See also

[`weibull_distrib`](https://statmodels7.github.io/distributions7/reference/weibull_distrib.md)
