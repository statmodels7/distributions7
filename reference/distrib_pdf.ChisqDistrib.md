# Chi-Squared Density

\$\$f(y; \mu) = \dfrac{y^{\mu/2 - 1} e^{-y/2}}{2^{\mu/2}\Gamma(\mu/2)},
\qquad y \> 0\$\$

## Arguments

- distrib:

  A `ChisqDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`chisq_distrib`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md)
