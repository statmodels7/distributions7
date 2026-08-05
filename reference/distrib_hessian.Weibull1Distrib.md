# Weibull Analytical Observed Hessian

Closed-form second derivatives of the Weibull log-density:
\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} =
\dfrac{\sigma}{\mu^2}\left\\1 - (1 + \sigma) u\right\\, \qquad
\dfrac{\partial^2 \ell}{\partial \sigma^2} = -\dfrac{1}{\sigma^2} - u
(\log z)^2,\$\$ \$\$\dfrac{\partial^2 \ell}{\partial \mu \\ \partial
\sigma} = \dfrac{1}{\mu}\left(u - 1 + \sigma u \log z\right).\$\$

## Arguments

- distrib:

  A `Weibull1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- scale:

  Either `"parameter"` or `"link"`.

- ...:

  Unused.

## Value

A named list of second derivatives.

## See also

[`weibull1_distrib`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
