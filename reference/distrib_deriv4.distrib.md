# Default Fourth-Order Derivatives for `distrib` Objects

The fallback for a family that registers no fourth-order method.
Observed derivatives come from
[`numerical_deriv4()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md),
a **second** difference of
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
rather than a difference of the third order, so the package's rule
against nesting one difference inside another holds here; expected ones
from
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
at the strategy `approx` names.

As at the order below, no family shipped in this package reaches it.
Being a second difference it is the least accurate route the package
offers, and its step is chosen accordingly.

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

- approx:

  Which strategy takes the expectation, read only when `expected` is
  `TRUE`: `"integrate"` (the default at this order), `"bartlett"`,
  `"opg"` or `"mc"`. See
  [`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- nsim:

  Monte Carlo sample size, a single positive number, read only under
  `approx = "mc"`. Defaults to 10000.

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
for what `approx` selects.
