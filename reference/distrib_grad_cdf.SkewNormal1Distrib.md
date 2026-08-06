# Skew Normal Log-CDF Derivatives

Closed form in the location and scale, the family being location-scale
in them; the shape direction is a derivative of Owen's T in its second
argument and is differenced.

## Arguments

- distrib:

  A `SkewNormal1Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu`, `sigma` and `alpha`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.

## See also

[`skewnormal1_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
