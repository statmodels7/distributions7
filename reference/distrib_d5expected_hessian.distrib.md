# Default Fifth Derivative of the Expected Information

Applies one tensor stencil to the highest derivative of the expected
information that the family computes itself, through
[`numerical_d5expected()`](https://statmodels7.github.io/distributions7/reference/numerical_d5expected.md),
and signals an error where the family has no analytic second derivative.

## Arguments

- distrib:

  A distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  Either `"parameter"` or `"link"`, passed to
  [`numerical_d5expected()`](https://statmodels7.github.io/distributions7/reference/numerical_d5expected.md),
  which perturbs on that scale.

- method:

  `NULL`: the generic handles an object before dispatching on the
  family. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Passed to
  [`numerical_d5expected()`](https://statmodels7.github.io/distributions7/reference/numerical_d5expected.md),
  which accepts `accuracy`.

## Value

A named list keyed as
[`d5expected_names()`](https://statmodels7.github.io/distributions7/reference/d5expected_names.md).
