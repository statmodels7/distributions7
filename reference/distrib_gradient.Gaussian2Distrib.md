# Gaussian Analytical Gradient in Mean and Variance

With \\r = y - \mu\\ and \\v = \sigma^2\\,
\$\$\dfrac{\partial\ell}{\partial\mu} = \dfrac{r}{v}, \qquad
\dfrac{\partial\ell}{\partial v} = \dfrac{r^2 - v}{2v^2}\$\$

## Arguments

- distrib:

  A `Gaussian2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `sigma2`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of first derivatives.

## See also

[`gaussian2_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
