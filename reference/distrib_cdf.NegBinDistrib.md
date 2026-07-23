# Negative Binomial Cumulative Distribution Function

Computes the cumulative distribution function for the Negative Binomial
distribution: \$\$F(q; \mu, \theta) = \sum\_{k=0}^{\lfloor q \rfloor}
P(Y=k; \mu, \theta)\$\$

## Arguments

- distrib:

  A `NegBinDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameters `mu` and `theta`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## See also

[`negbin_distrib`](https://statmodels7.github.io/distributions7/reference/negbin_distrib.md)
