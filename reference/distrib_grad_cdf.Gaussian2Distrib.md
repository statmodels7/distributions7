# Gaussian Log-CDF Derivatives in Mean and Variance

Closed form, by the chain rule on the scale parametrization's
derivatives through \\\sigma = \sqrt{\sigma^2}\\.

## Arguments

- distrib:

  A `Gaussian2Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `sigma2`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.

## See also

[`gaussian2_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
