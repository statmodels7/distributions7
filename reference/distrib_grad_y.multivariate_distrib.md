# Response Derivatives of a Multivariate Distribution

Refused on the base class rather than served numerically: the univariate
fallbacks difference along a line, and the derivative of a multivariate
log-density in its response is a vector (a matrix at second order) whose
shape the base class cannot guess. A family that has the closed form
registers it, as the gaussian and the Student t do.

## Arguments

- distrib:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- y:

  An \\n \times p\\ matrix of observations.

- theta:

  A named list of parameters.

- ...:

  Unused.

## Value

No return value; called for its error.
