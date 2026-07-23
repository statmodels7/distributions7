# Lognormal Cumulative Distribution Function

Computes the cumulative distribution function for the Lognormal
distribution (with \\\sigma = \sqrt{\sigma^2}\\ on the log scale):
\$\$F(q; \mu, \sigma^2) = \Phi\left(\dfrac{\log q -
\mu}{\sqrt{\sigma^2}}\right), \quad q \> 0\$\$ where \\\Phi\\ is the
standard normal CDF.

## Arguments

- distrib:

  A `LognormalDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameters `mu` and `sigma2`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## See also

[`lognormal_distrib`](https://statmodels7.github.io/distributions7/reference/lognormal_distrib.md)
