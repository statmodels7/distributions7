# Beta Probability Density Function

Computes the probability density function for the Beta distribution:
\$\$f(y; \mu, \phi) =
\dfrac{\Gamma(\phi)}{\Gamma(\mu\phi)\Gamma((1-\mu)\phi)} y^{\mu\phi-1}
(1-y)^{(1-\mu)\phi-1}\$\$

## Arguments

- distrib:

  A `BetaDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `phi`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`beta_distrib`](https://statmodels7.github.io/distributions7/reference/beta_distrib.md)
