# Default Derivative of the Expected Information

One central difference of
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
per parameter, refused where that quantity is itself approximated.

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

  Passed through.

- ...:

  Unused.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md).
