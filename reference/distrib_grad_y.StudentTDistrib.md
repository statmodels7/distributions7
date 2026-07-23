# Student's t Response Derivatives

Closed-form derivatives of the Student's t log-density with respect to
the response. Let \\r = y - \mu\\ and \\d = \nu\sigma^2 + r^2\\:
\\\partial \ell / \partial y = -(\nu+1)r/d\\ and \\\partial^2 \ell /
\partial y^2 = (\nu+1)(r^2 - \nu\sigma^2)/d^2\\.

## Arguments

- distrib:

  A `StudentTDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu`, `sigma` and `nu`.

## Value

A numeric vector.

## See also

[`student_t_distrib`](https://statmodels7.github.io/distributions7/reference/student_t_distrib.md)
