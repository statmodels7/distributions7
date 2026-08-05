# Inverse Gaussian Analytical Gradient in Mean and Shape

\$\$\dfrac{\partial\ell}{\partial\mu} = \dfrac{\lambda(y-\mu)}{\mu^3},
\qquad \dfrac{\partial\ell}{\partial\lambda} = \dfrac{1}{2\lambda} -
\dfrac{(y-\mu)^2}{2\mu^2 y}\$\$

## Arguments

- distrib:

  An `InvGauss2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `lambda`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of first derivatives.

## See also

[`invgauss2_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md)
