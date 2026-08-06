# Pseudo-Huber Log-CDF Hessian

Closed form in the location and scale block; the shape is differenced.
Differencing the cdf itself would be poor here, it being a quadrature.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu`, `sigma` and `nu`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list keyed as
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## See also

[`pseudohuber_distrib`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
