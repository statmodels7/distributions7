# Bernoulli Log-CDF Gradient

Closed form, the binomial identity at \\n = 1\\: the derivative is
\\-1\\ at \\k = 0\\ and zero at \\k = 1\\.

## Arguments

- distrib:

  A `BernoulliDistrib` object.

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

[`bernoulli_distrib`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md)
