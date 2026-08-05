# Student's t Analytical Gradient

Computes the analytical gradient (first derivatives) of the Student's t
log-density with respect to the parameters \\\mu\\, \\\sigma\\, and
\\\nu\\.

\$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{(\nu+1)(y-\mu)}{\nu\sigma^2 + (y-\mu)^2}\$\$ \$\$\dfrac{\partial
\ell}{\partial \sigma} = \dfrac{\nu\left\[(y-\mu)^2 -
\sigma^2\right\]}{\sigma\left\[\nu\sigma^2 + (y-\mu)^2\right\]}\$\$
\$\$\dfrac{\partial \ell}{\partial \nu} = \dfrac{1}{2}\left\[
-\dfrac{1}{\nu} - \psi\left(\dfrac{\nu}{2}\right) +
\psi\left(\dfrac{\nu+1}{2}\right) +
\dfrac{(\nu+1)(y-\mu)^2}{\nu\left\[\nu\sigma^2 + (y-\mu)^2\right\]} -
\log\left(1 + \dfrac{(y-\mu)^2}{\nu\sigma^2}\right) \right\]\$\$

## Arguments

- distrib:

  A `StudentT1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu`, `sigma`, and `nu`.

## Value

A list containing the vectors of first derivatives.

## See also

[`student_t1_distrib`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
