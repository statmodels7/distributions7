# Default Fifth-Order Derivatives for `distrib` Objects

Defines the route that a family takes at the fifth order when it
registers no method of its own: one tensor-product difference of its
highest analytic order, through
[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md).

## Arguments

- distrib:

  An object inheriting from `distrib`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned by the generic.

- scale:

  Passed through to
  [`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md),
  which reads the fourth order on that scale and perturbs on it.

- ...:

  Passed to
  [`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md),
  which accepts `accuracy`.

## Value

A named list of fifth-derivative component vectors, each of length
`length(y)`, keyed lexicographically as
[`deriv_names(distrib@params, 5)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
gives them.

## Details

This is registered on the base class. The families that write the fifth
order out themselves (see
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md))
and the wrappers override it. Among the families that reach it, the
quantity that is differenced varies:
[`has_exact_deriv4()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv4.md)
reports whether that quantity is analytic, and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
emits no order-5 row where it is not.

## See also

[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md),
which does the differencing;
[`distrib_deriv4.distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.distrib.md)
for the order below.
