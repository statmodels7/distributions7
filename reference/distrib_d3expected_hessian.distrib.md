# Default Third Derivative of the Expected Information

Takes one central difference of
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
per parameter and signals an error where the family has no analytic
second derivative.

## Arguments

- distrib:

  A distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  Either `"parameter"` or `"link"`.

- method:

  `NULL`: the generic handles an object before dispatching on the
  family. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused.

## Value

A named list keyed as
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md).
