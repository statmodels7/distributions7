# Gaussian Analytical Gradient in Mean and Precision

With \\r = y - \mu\\, \$\$\dfrac{\partial\ell}{\partial\mu} = \tau r,
\qquad \dfrac{\partial\ell}{\partial\tau} = \dfrac{1}{2\tau} -
\dfrac{r^2}{2}\$\$

## Arguments

- distrib:

  A `Gaussian3Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `tau`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of first derivatives.

## See also

[`gaussian3_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md)
