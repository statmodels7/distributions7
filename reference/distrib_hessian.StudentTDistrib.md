# Student's t Analytical Observed Hessian

Computes the analytical observed Hessian (second derivatives) of the
Student's t log-density with respect to the parameters \\\mu\\,
\\\sigma\\, and \\\nu\\.

\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} =
\dfrac{(\nu+1)\left\[(y-\mu)^2 -
\nu\sigma^2\right\]}{\left\[\nu\sigma^2 + (y-\mu)^2\right\]^2}\$\$
\$\$\dfrac{\partial^2 \ell}{\partial \sigma^2} =
\dfrac{\nu\left\[\nu\sigma^4 - (3\nu+1)\sigma^2(y-\mu)^2 -
(y-\mu)^4\right\]}{\sigma^2\left\[\nu\sigma^2 + (y-\mu)^2\right\]^2}\$\$
\$\$\dfrac{\partial^2 \ell}{\partial \nu^2} = \dfrac{1}{4}\left\[
-\psi_1\left(\dfrac{\nu}{2}\right) +
\psi_1\left(\dfrac{\nu+1}{2}\right) + \dfrac{2\left(\nu\sigma^4 +
(y-\mu)^4\right)}{\nu\left\[\nu\sigma^2 + (y-\mu)^2\right\]^2}
\right\]\$\$ \$\$\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} =
-\dfrac{2\nu(\nu+1)\sigma(y-\mu)}{\left\[\nu\sigma^2 +
(y-\mu)^2\right\]^2}\$\$ \$\$\dfrac{\partial^2 \ell}{\partial \mu
\partial \nu} = \dfrac{(y-\mu)\left\[(y-\mu)^2 -
\sigma^2\right\]}{\left\[\nu\sigma^2 + (y-\mu)^2\right\]^2}\$\$
\$\$\dfrac{\partial^2 \ell}{\partial \sigma \partial \nu} =
\dfrac{(y-\mu)^2\left\[(y-\mu)^2 -
\sigma^2\right\]}{\sigma\left\[\nu\sigma^2 + (y-\mu)^2\right\]^2}\$\$

## Arguments

- distrib:

  A `StudentTDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu`, `sigma`, and `nu`.

## Value

A list containing the vectors of second derivatives.

## See also

[`student_t_distrib`](https://statmodels7.github.io/distributions7/reference/student_t_distrib.md)
