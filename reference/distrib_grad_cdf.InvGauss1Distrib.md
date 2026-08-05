# Inverse Gaussian Log-CDF Gradient

Closed form, obtained by differentiating the elementary representation
\\F(y) = \Phi(a) + e^{2/(\phi\mu)}\Phi(b)\\. The exponential factor is
combined with \\\Phi(b)\\ on the log scale, since it overflows for small
\\\phi\mu\\ exactly where \\\Phi(b)\\ underflows.

## Arguments

- distrib:

  An `InvGauss1Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `phi`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.

## See also

[`invgauss1_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md)
