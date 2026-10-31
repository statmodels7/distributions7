# Default Sixth-Order Derivatives for `distrib` Objects

Defines the route of every family at the sixth order: one tensor-product
difference of its highest analytic order, through
[`numerical_deriv6()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv6.md).

## Arguments

- distrib:

  An object inheriting from `distrib`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned by the generic.

- scale:

  Passed through to
  [`numerical_deriv6()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv6.md),
  which reads the differenced order on that scale and perturbs on it.

- ...:

  Passed to
  [`numerical_deriv6()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv6.md),
  which accepts `accuracy`.

## Value

A named list of sixth-derivative component vectors, each of length
`length(y)`, keyed lexicographically as
[`deriv_names(distrib@params, 6)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
gives them.

## Details

This method is registered on the base class, and no family or wrapper
overrides it. The quantity that is differenced is the fifth order where
[`has_exact_deriv5()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv5.md)
is `TRUE`, and the order that
[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md)
differences otherwise.

## See also

[`numerical_deriv6()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv6.md),
which does the differencing;
[`distrib_deriv5.distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.distrib.md)
for the order below.
