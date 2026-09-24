# Derivatives of the Expected Information, Reparametrized Families

The parent's derivatives carried through the map by
[`reparam_dexpected()`](https://statmodels7.github.io/distributions7/reference/reparam_dexpected.md).
A parent without an analytic second derivative leaves the reparametrized
family without one, the error naming the parent.

## Arguments

- distrib:

  A reparametrized distribution.

- y:

  A numeric vector of observations.

- theta:

  A named list of the new parameters.

- scale:

  `"parameter"` or `"link"`.

- approx, nsim:

  Unused.

- ...:

  Unused.

- threads:

  Passed to the parent.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).
