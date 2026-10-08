# Default Fourth Derivative of the Expected Information

One second-order stencil on
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md);
an error where the family has no analytic second derivative.

## Arguments

- distrib:

  A distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  Either `"parameter"` or `"link"`.

- approx, nsim:

  Unused.

- ...:

  Unused.

## Value

A named list keyed as
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md).
