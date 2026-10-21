# Default Mixed Derivatives for Continuous Distributions

The fallback applies one tensor stencil to the highest quantity that the
family implements itself, through
[`mixed_tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/mixed_tensor_derivatives.md).
A family with an analytic response gradient or score differences it
once, in the parameter or in the response; a family with only a density
takes one four-point mixed stencil of the log-density.

Every continuous family shipped in this package registers a closed form,
so this method serves a family defined outside the package, which
obtains the mixed block only from its density.

## Arguments

- distrib:

  An object inheriting from `continuous_distrib` that registers no
  method of its own.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned by the generic.

- scale:

  Handled by the generic after dispatch; this method always returns the
  parameter scale.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter, keyed by
`distrib@params`. Each vector has length `length(y)`.

## See also

[`mixed_tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/mixed_tensor_derivatives.md),
which does the work;
[`numerical_cross_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross_y.md)
for the differenced part of a family's own method;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic;
[`distrib_cross_y.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Gaussian1Distrib.md)
for a closed form to compare against.
