# Pseudo-Huber Cumulative Distribution Function

Computes the cumulative distribution function for the Pseudo-Huber
distribution by numerical integration of the density. The distribution
is symmetric around \\\mu\\, so the upper tail is computed by
reflection.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameters `mu`, `sigma` and `nu`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of cumulative probabilities.

## See also

[`pseudohuber_distrib`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
