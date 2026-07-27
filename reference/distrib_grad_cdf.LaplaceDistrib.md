# Laplace Log-CDF Derivatives

Closed form, from the location-scale structure. Note that the second
derivatives inherit the kink at \\y = \mu\\, where
\\\partial\ell/\partial y\\ does not exist;
[`param_smoothness`](https://statmodels7.github.io/distributions7/reference/param_smoothness.md)
records this.

## Arguments

- distrib:

  A `LaplaceDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `b`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.

## See also

[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
