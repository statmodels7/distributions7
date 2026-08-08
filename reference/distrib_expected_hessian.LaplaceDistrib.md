# Laplace Analytical Expected Hessian (Fisher Information)

Computes the expected Hessian (negative Fisher information) of the
Laplace log-density. Because the log-likelihood is not differentiable in
\\\mu\\, the second Bartlett identity fails: \\\mathbb{E}\[\partial^2
\ell / \partial \mu^2\] = 0\\, yet the Fisher information for \\\mu\\ is
\\1/\sigma^2\\. The expected Hessian is therefore defined here from the
variance of the score, which is what the Fisher information *is*
whenever the score exists, whether or not the identity relating it to
\\-\mathbb{E}\[H\]\\ holds:

\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\sigma^2}, \qquad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu \partial \sigma}\right\] = 0, \qquad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right\] =
-\dfrac{1}{\sigma^2}\$\$

## Arguments

- distrib:

  A `LaplaceDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A list containing the vectors of expected second derivatives.

## See also

[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
