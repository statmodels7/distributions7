# Mixed Derivatives When Only Two Parameters Are Location-Scale

The
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
body for a family whose first two parameters are a location and a scale
and whose remaining ones are shapes: the first two components come from
[`loc_scale_cross_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_block.md)
in closed form, and the rest from one central difference of
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
through
[`numerical_cross_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross_y.md),
which is asked for those parameters alone.

Registered on the skew \\t\\, whose `alpha` and `nu` components are
differenced because its density carries \\T\_{\nu+1}\\, a Student \\t\\
distribution function whose derivative in the degrees of freedom has no
elementary form. Measured against Richardson on the analytic response
gradient, the worst of its four components is \\9.2\times10^{-10}\\
relative, which is one stencil's accuracy.

## Usage

``` r
partial_loc_scale_cross_y(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  ...
)
```

## Arguments

- distrib:

  An object inheriting from `distrib` with at least three parameters, a
  location and a scale first.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned by the generic.

- scale:

  Handled by the generic after dispatch; this body always returns the
  parameter scale.

- ...:

  Unused. Not forwarded to
  [`numerical_cross_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross_y.md),
  so a step given here does not reach the difference.

## Value

A named list with one numeric vector per parameter, keyed by
`distrib@params`, each of length `length(y)`. The first two components
are exact and the rest carry one stencil's error.

## See also

[`loc_scale_cross_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_block.md)
for the exact half;
[`numerical_cross_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross_y.md)
for the differenced half;
[`loc_scale_cross_y()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_y.md)
where every parameter is covered exactly.
