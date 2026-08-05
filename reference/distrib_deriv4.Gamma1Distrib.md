# Gamma Fourth-Order Derivatives in Mean and Dispersion

Closed form, observed or expected.

## Arguments

- distrib:

  A `Gamma1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `phi`.

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

A named list of fourth-derivative components.

## See also

[`gamma1_distrib`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
