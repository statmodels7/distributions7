# Gumbel Cumulative Distribution Function

Computes the cumulative distribution function for the Gumbel
distribution: \$\$F(q; \mu, \sigma) = \exp\left\\-e^{-z}\right\\, \qquad
z = (q - \mu)/\sigma\$\$

## Arguments

- distrib:

  A `GumbelDistrib` object.

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

## Details

The lower tail is exact on the log scale, since \\\log F = -e^{-z}\\,
and the upper tail uses `expm1`, so neither loses precision in its own
tail.

## See also

[`gumbel_distrib`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
