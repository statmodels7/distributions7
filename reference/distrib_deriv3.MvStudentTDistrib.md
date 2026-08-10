# Multivariate Student t Third and Fourth Derivatives

Closed form. Every term of the log-density is elementary in \\\nu\\, so
the obstruction the univariate skew t meets – the derivative of a
Student t distribution function in its degrees of freedom – does not
arise here, and all four orders close.

## Arguments

- distrib:

  A `MvStudentTDistrib` object.

- y:

  An \\n \times p\\ matrix of observations.

- theta:

  A named list of parameters.

- expected:

  Logical; the expectation is approximated by sampling.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Strategy label; sampling is the only multivariate route.

- nsim:

  Monte Carlo sample size.

- ...:

  Unused.

## Value

A named list of third-derivative component vectors.

## See also

[`mvstudent_t_distrib`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md)
