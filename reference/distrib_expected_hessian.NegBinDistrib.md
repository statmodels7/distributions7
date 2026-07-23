# Negative Binomial Analytical Expected Hessian

Computes the analytical expected Hessian of the Negative Binomial
log-probability with respect to the parameters \\\mu\\ and \\\theta\\.

\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{\theta}{\mu(\theta+\mu)}\$\$
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \theta^2}\right\]
= \mathbb{E}\[\psi_1(Y+\theta)\] - \psi_1(\theta) +
\dfrac{\mu}{\theta(\theta+\mu)}\$\$
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu \partial
\theta}\right\] = 0\$\$

The term \\\mathbb{E}\[\psi_1(Y+\theta)\]\\ has no closed form and is
evaluated by summing over the support up to a far-tail quantile (with a
tail-mass correction).

## Arguments

- distrib:

  A `NegBinDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `theta`.

## Value

A list containing the vectors of expected second derivatives.

## See also

[`negbin_distrib`](https://statmodels7.github.io/distributions7/reference/negbin_distrib.md)
