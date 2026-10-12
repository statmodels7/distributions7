# Higher Response Derivatives of a Reparametrized Distribution

The parent's third and fourth response derivatives, read at the parent's
parameters. A reparametrization acts on \\\theta\\ and these derivatives
are taken in \\y\\, so the two do not interact: the chain rule does not
enter, and the result of the parent is returned exactly.

The same holds at the first and second orders. It does not hold for the
mixed block
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md),
which takes one derivative in each, so that the first-order Jacobian of
the map enters.

## Arguments

- distrib:

  A `ReparamContinuousDistrib` object, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters of the **new** parametrization, carried to
  the parent's by
  [`reparam_theta()`](https://statmodels7.github.io/distributions7/reference/reparam_theta.md)
  before the call.

- ...:

  Passed to the parent's method.

## Value

A numeric vector of length `length(y)`: the parent's derivative of that
order at the mapped parameters.

## See also

[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
for the wrapper;
[`distrib_cross_y.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.ReparamContinuousDistrib.md),
where the map does enter;
[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
for the generic.
