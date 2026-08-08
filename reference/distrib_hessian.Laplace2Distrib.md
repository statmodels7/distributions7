# Laplace Analytical Observed Hessian in Location and Rate

Almost everywhere, \$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = 0,
\qquad \dfrac{\partial^2 \ell}{\partial \mu \partial \lambda} =
\mathrm{sign}(y-\mu), \qquad \dfrac{\partial^2 \ell}{\partial \lambda^2}
= -\dfrac{1}{\lambda^2}\$\$

The log-likelihood is linear in \\\lambda\\ apart from the
\\\log\lambda\\ term, so every derivative in \\\lambda\\ beyond the
first is free of the data. As for
[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md),
the zero second derivative in \\\mu\\ means Newton-Raphson cannot update
\\\mu\\;
[`distrib_expected_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
supplies the Fisher information \\\lambda^2\\.

## Arguments

- distrib:

  A `Laplace2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `lambda`.

## Value

A list containing the vectors of second derivatives.

## See also

[`laplace2_distrib`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
