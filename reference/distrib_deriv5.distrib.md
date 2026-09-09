# Default Fifth-Order Derivatives for `distrib` Objects

The route every family takes at the fifth order: one central difference
of the fourth, through
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

  Unused.

## Value

A named list of fifth-derivative component vectors, each of length
`length(y)`, keyed lexicographically as
[`deriv_names(distrib@params, 5)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
gives them.

## Details

This is registered on the base class and nothing overrides it, so it is
not a fallback in the usual sense – no family writes the fifth order out
yet. What differs between families is the quantity being differenced.
[`has_exact_deriv4()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv4.md)
answers whether that quantity is analytic, and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
reports the order-5 row as unchecked where it is not.

## See also

[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md),
which does the differencing;
[`distrib_deriv4.distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.distrib.md)
for the order below.
