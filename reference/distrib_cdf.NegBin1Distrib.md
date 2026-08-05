# NB1 Cumulative Distribution Function

[`pnbinom`](https://rdrr.io/r/stats/NegBinomial.html) at size
\\\mu/\theta\\ and probability \\1/(1+\theta)\\.

## Arguments

- distrib:

  A `NegBin1Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `theta`.

- lower.tail:

  Logical; if `TRUE` (default), \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned as logarithms.

## Value

A numeric vector of cumulative probabilities.

## See also

[`negbin1_distrib`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
