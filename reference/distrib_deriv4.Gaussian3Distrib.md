# Gaussian Fourth-Order Derivatives in Mean and Precision

Closed form, and free of the response: the only non-zero component is
\\\ell^{(\tau\tau\tau\tau)} = -3/\tau^4\\.

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

A named list of fourth-derivative components.

## See also

[`gaussian3_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md)
