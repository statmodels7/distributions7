# Gaussian Third-Order Derivatives in Mean and Variance

Closed form, observed or expected.

## Arguments

- distrib:

  A `Gaussian2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `sigma2`.

- expected:

  Logical; if `TRUE`, returns the expected derivatives.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is closed form.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list of third-derivative components.

## See also

[`gaussian2_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
