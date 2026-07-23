# Student's t Response Second Derivative

Closed-form \\\partial^2 \ell / \partial y^2 = (\nu+1)(r^2 -
\nu\sigma^2)/(\nu\sigma^2 + r^2)^2\\, \\r = y - \mu\\.

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
