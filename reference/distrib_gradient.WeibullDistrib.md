# Weibull Analytical Gradient

Closed-form first derivatives of the Weibull log-density, written in \\u
= (y/\mu)^{\sigma}\\ and \\\log z = \log(y/\mu)\\: \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{\sigma}{\mu}(u - 1), \qquad \dfrac{\partial
\ell}{\partial \sigma} = \dfrac{1}{\sigma} + (1 - u)\log z\$\$

## Arguments

- distrib:

  A `WeibullDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- scale:

  Either `"parameter"` or `"link"`.

- ...:

  Unused.

## Value

A named list of first derivatives.

## See also

[`weibull_distrib`](https://statmodels7.github.io/distributions7/reference/weibull_distrib.md)
