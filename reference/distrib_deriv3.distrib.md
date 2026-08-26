# Default Third-Order Derivatives for `distrib` Objects

The fallback for a family that registers no third-order method. Observed
derivatives come from
[`numerical_deriv3()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv3.md),
one central difference of
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
along each parameter; expected ones from
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
at the strategy `approx` names.

**No family shipped in this package reaches this method for its observed
derivatives.** All 46 write the third order out, 24 of them in compiled
kernels. It exists for a family defined outside the package, which gets
four orders from a density alone, and it is the reference the analytical
kernels are validated against.

## Arguments

- distrib:

  An object inheriting from `distrib` that registers no method of its
  own.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned by the generic.

- expected:

  Logical of length 1. `FALSE`, the default, differences the Hessian;
  `TRUE` takes the expectation of the observed derivatives through
  [`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md).

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

A named list of third-derivative component vectors, each of length
`length(y)`, keyed lexicographically as
[`deriv_names(distrib@params, 3)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
gives them.

## See also

[`numerical_deriv3()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv3.md),
which does the differencing;
[`distrib_deriv4.distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.distrib.md)
for the order above;
[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for what `approx` selects.
