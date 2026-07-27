# Binomial Log-CDF Gradient

Closed form: \\\partial F(k)/\partial\mu = -n\\\mathrm{dbinom}(k; n-1,
\mu)\\.

## Arguments

- distrib:

  A `BinomialDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list with one element.

## See also

[`binomial_distrib`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
