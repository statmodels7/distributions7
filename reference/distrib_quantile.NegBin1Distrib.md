# NB1 Quantile Function

[`qnbinom`](https://rdrr.io/r/stats/NegBinomial.html) at the same size
and probability.

## Arguments

- distrib:

  A `NegBin1Distrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing `mu` and `theta`.

- lower.tail:

  Logical; if `TRUE` (default), `p` is \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, `p` is given as its logarithm.

## Value

A numeric vector of quantiles.

## See also

[`negbin1_distrib`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
