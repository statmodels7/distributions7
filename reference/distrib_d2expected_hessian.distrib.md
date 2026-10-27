# Default Second Derivative of the Expected Information

Signals an error: a family without its own method has no second
derivative of the expected information, because differencing the first
derivative would be a difference of a difference.

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

Does not return.
