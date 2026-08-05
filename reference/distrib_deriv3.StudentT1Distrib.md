# Student's t Analytical Third-Order Derivatives

Closed-form observed third-order derivatives of the Student's t
log-density. The expected third derivatives have no closed form, so
`expected = TRUE` falls back to the numerical
[`expectation`](https://statmodels7.github.io/distributions7/reference/expectation.md)
of the observed derivatives.

## Arguments

- distrib:

  A `StudentT1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu`, `sigma` and `nu`.

- expected:

  Logical; if `TRUE`, returns the (numerically integrated) expected
  third derivatives.

## Value

A named list of third-derivative component vectors.

## See also

[`student_t1_distrib`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
