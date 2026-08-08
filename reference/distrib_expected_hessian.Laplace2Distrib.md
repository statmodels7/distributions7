# Laplace Analytical Expected Hessian in Location and Rate

As for
[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md),
the second Bartlett identity fails in \\\mu\\ and the expected Hessian
is defined from the variance of the score:

\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\lambda^2, \qquad \mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu
\partial \lambda}\right\] = 0, \qquad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \lambda^2}\right\] = -\dfrac{1}{\lambda^2}\$\$

Because the closed form exists, the `approx` argument is ignored.

## Arguments

- distrib:

  A `Laplace2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `lambda`.

## Value

A list containing the vectors of expected second derivatives.

## See also

[`laplace2_distrib`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
