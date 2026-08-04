# Chi-Squared Quantile Function

The inverse of [`pchisq`](https://rdrr.io/r/stats/Chisquare.html).

## Arguments

- distrib:

  A `ChisqDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing the parameter `mu`.

- lower.tail:

  Logical; if `TRUE` (default), `p` is \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, `p` is given as its logarithm.

## Value

A numeric vector of quantiles.

## See also

[`chisq_distrib`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md)
