# Gaussian Analytical Observed Hessian in Mean and Variance

\$\$\ell^{(\mu\mu)} = -\dfrac{1}{v}, \qquad \ell^{(\mu v)} =
-\dfrac{r}{v^2}, \qquad \ell^{(vv)} = \dfrac{1}{2v^2} -
\dfrac{r^2}{v^3}\$\$

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

A named list of second derivatives.

## See also

[`gaussian2_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
