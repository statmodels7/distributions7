# Skew t Log-CDF Derivatives

Closed form in the location and scale; the shape and the degrees of
freedom are differenced.

## Arguments

- distrib:

  A `SkewTDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu`, `sigma`, `alpha` and `nu`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.

## See also

[`skewt_distrib`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
