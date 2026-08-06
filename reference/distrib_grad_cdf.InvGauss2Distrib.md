# Inverse Gaussian Log-CDF Gradient in Mean and Shape

Closed form, by the chain rule on the dispersion parametrization's
gradient through \\\phi = 1/\lambda\\. The second order stays on the
fallback, the parent having no closed form there.

## Arguments

- distrib:

  An `InvGauss2Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `lambda`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.

## See also

[`invgauss2_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md)
