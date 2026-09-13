# Default Mixed Third-Response Derivatives for Continuous Distributions

Falls back to one central difference of
[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
in each parameter, through
[`numerical_cross3_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross3_y.md).

## Arguments

- distrib:

  A `continuous_distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  One of `"parameter"` or `"link"`, applied by the generic after
  dispatch.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector per parameter, keyed by
`distrib@params`.

## See also

[`distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross3_y.md)
for the generic,
[`numerical_cross3_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross3_y.md),
which does the work, and
[`loc_cross3_y()`](https://statmodels7.github.io/distributions7/reference/loc_cross3_y.md)
for the location families.
