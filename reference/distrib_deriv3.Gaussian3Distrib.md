# Gaussian Third-Order Derivatives in Mean and Precision

Closed form. Every third derivative is free of the response, so the
observed and the expected ones coincide and `expected` changes nothing.

## Arguments

- distrib:

  A `Gaussian3Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `tau`.

- expected:

  Logical; makes no difference here.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list of third-derivative components.

## See also

[`gaussian3_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md)
