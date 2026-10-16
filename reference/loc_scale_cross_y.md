# Mixed Derivatives of a Location-Scale Family

The
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
body shared by every family that is location-scale in both its
parameters: it takes the two components of
[`loc_scale_cross_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_block.md)
and names them after the family's own parameters. Registered on the
logistic, the Cauchy, the Gumbel and the Laplace.

## Usage

``` r
loc_scale_cross_y(distrib, y, theta, scale = c("parameter", "link"), ...)
```

## Arguments

- distrib:

  An object inheriting from `distrib` with exactly two parameters, a
  location and a scale, in that order.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned by the generic.

- scale:

  Handled by the generic after dispatch; this body always returns the
  parameter scale.

- ...:

  Unused.

## Value

A named list of two numeric vectors, keyed by `distrib@params`, each of
length `length(y)`.

## See also

[`loc_scale_cross_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_block.md)
for the identity;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
