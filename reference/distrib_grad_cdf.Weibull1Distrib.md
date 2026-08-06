# Weibull Log-CDF Derivatives

Closed form; see
[`weibull_cdf_deriv`](https://statmodels7.github.io/distributions7/reference/weibull_cdf_deriv.md).

## Arguments

- distrib:

  A `Weibull1Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `sigma`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.

## See also

[`weibull1_distrib`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
