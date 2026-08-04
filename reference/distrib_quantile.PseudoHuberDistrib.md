# Pseudo-Huber Quantile Function

Computes the quantile function for the Pseudo-Huber distribution by
root-finding on the numerical CDF. Symmetry around \\\mu\\ is exploited:
\\Q(1/2) = \mu\\ and \\Q(p) = 2\mu - Q(1-p)\\ for \\p \> 1/2\\.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing the parameters `mu`, `sigma` and `nu`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\,
  otherwise \\P(Y \> p)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of quantiles.

## See also

[`pseudohuber_distrib`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
