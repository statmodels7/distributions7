# Student's t Probability Density Function

Computes the probability density function for the Student's t
distribution: \$\$f(y; \mu, \sigma, \nu) =
\dfrac{\Gamma\left(\dfrac{\nu+1}{2}\right)}{\sigma\sqrt{\nu\pi}\\\Gamma\left(\dfrac{\nu}{2}\right)}
\left(1 + \dfrac{(y-\mu)^2}{\nu\sigma^2}\right)^{-\dfrac{\nu+1}{2}}\$\$

## Arguments

- distrib:

  A `StudentTDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu`, `sigma`, and `nu`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`student_t_distrib`](https://statmodels7.github.io/distributions7/reference/student_t_distrib.md)
