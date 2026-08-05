# Negative Binomial Log-CDF Gradient

Closed form in \\\mu\\, \\\partial F(k)/\partial\mu =
-f(k)(k+\theta)/(\theta+\mu)\\, which reduces to the Poisson identity as
\\\theta\to\infty\\. The \\\theta\\ direction is a derivative of the
incomplete beta in its parameter, has no elementary form, and keeps the
exact summation.

## Arguments

- distrib:

  A `NegBin2Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `theta`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.

## See also

[`negbin2_distrib`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
