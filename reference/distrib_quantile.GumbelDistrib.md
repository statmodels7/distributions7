# Gumbel Quantile Function

Computes the quantile function for the Gumbel distribution: \$\$Q(p;
\mu, \sigma) = \mu - \sigma \log(-\log p)\$\$

## Arguments

- distrib:

  A `GumbelDistrib` object.

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

[`gumbel_distrib`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
