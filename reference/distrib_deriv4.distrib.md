# Default Fourth-Order Derivatives for `distrib` Objects

Provides the fallback for a family that registers no fourth-order
method. Observed derivatives come from
[`numerical_deriv4()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md):
a second difference of the family's own
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
rather than a difference of the third order, or, for a family without
its own Hessian, one tensor-product stencil on its score or its
log-density. On both routes each component is a single stencil, so the
package's rule against nesting one difference inside another holds.
Expected ones come from
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
by the method that `method` names.

As at the order below, the families shipped in this package write the
fourth order out themselves. A second difference of the Hessian is less
accurate than the first difference of the third order, and its step is
chosen accordingly.

## Arguments

- distrib:

  An object inheriting from `distrib` that registers no method of its
  own.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned by the generic.

- expected:

  Logical of length 1. `FALSE`, the default, differences the Hessian
  twice; `TRUE` takes the expectation of the observed derivatives.

- scale:

  Handled by the generic after dispatch; this method always returns the
  parameter scale.

- method:

  `NULL` (the default) or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).
  It is read only when `expected` is `TRUE`: `NULL` selects
  [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md),
  and an object selects its own route. The generic applies an object
  before dispatching on the family, so that the object is used also by a
  family that has a closed form. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused.

## Value

A named list of fourth-derivative component vectors, each of length
`length(y)`, keyed lexicographically as
[`deriv_names(distrib@params, 4)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
gives them.

## See also

[`numerical_deriv4()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md),
which does the differencing;
[`distrib_deriv3.distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.distrib.md)
for the order below;
[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for what `method` selects.
