# Default Mixed Derivatives for Continuous Distributions

The fallback: one central difference of
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
in each parameter, through
[`numerical_cross_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross_y.md)
at its default step of \\\varepsilon^{1/3} \approx 6.1\times10^{-6}\\.
The quantity differenced is the response gradient, so a family with an
analytic one pays for a single difference layer and reaches about
\\10^{-10}\\ relative accuracy; where the response gradient is itself a
difference, the two act on different variables and commute into one
four-point mixed stencil.

**No family shipped in this package reaches this method.** All 32
continuous families register a closed form, so this exists for a family
defined outside the package, which gets the mixed block for free from
its density alone.

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
`distrib@params`, each of length `length(y)`.

## See also

[`numerical_cross_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross_y.md),
which does the work;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic;
[`distrib_cross_y.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Gaussian1Distrib.md)
for a closed form to compare against.
