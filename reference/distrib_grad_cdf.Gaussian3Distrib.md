# Gaussian Log-CDF Derivatives in Mean and Precision

Closed form, by the chain rule on the scale parametrization's
derivatives through \\\sigma = \tau^{-1/2}\\.

## Arguments

- distrib:

  A `Gaussian3Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `tau`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.

## See also

[`gaussian3_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md)
