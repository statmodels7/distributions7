# Student's t Analytical Expected Hessian

Computes the analytical expected Hessian of the Student's t log-density
with respect to the parameters \\\mu\\, \\\sigma\\, and \\\nu\\.

\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{\nu+1}{\sigma^2(\nu+3)}\$\$
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right\]
= -\dfrac{2\nu}{\sigma^2(\nu+3)}\$\$
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \nu^2}\right\] =
\dfrac{1}{4}\left\[\psi_1\left(\dfrac{\nu+1}{2}\right) -
\psi_1\left(\dfrac{\nu}{2}\right)\right\] +
\dfrac{\nu+5}{2\nu(\nu+1)(\nu+3)}\$\$
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \sigma \partial
\nu}\right\] = \dfrac{2}{\sigma(\nu+1)(\nu+3)}\$\$

The parameter \\\mu\\ is orthogonal to \\\sigma\\ and \\\nu\\ (mixed
expected derivatives are 0).

## Arguments

- distrib:

  A `StudentTDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu`, `sigma`, and `nu`.

## Value

A list containing the vectors of expected second derivatives.

## See also

[`student_t_distrib`](https://statmodels7.github.io/distributions7/reference/student_t_distrib.md)
