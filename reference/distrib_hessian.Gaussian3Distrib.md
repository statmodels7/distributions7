# Gaussian Analytical Observed Hessian in Mean and Precision

\$\$\ell^{(\mu\mu)} = -\tau, \qquad \ell^{(\mu\tau)} = r, \qquad
\ell^{(\tau\tau)} = -\dfrac{1}{2\tau^2}\$\$ The mean block is free of
the data here, as it is in every parametrization of this family, and so
is the precision block.

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

A named list of second derivatives.

## See also

[`gaussian3_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md)
