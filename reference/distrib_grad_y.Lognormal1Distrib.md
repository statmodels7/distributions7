# Lognormal Response Derivatives

Closed-form derivatives of the Lognormal log-density with respect to the
response. With \\r = \log y - \mu\\: \\\partial \ell / \partial y =
-\dfrac{1}{y}\left(1 + \dfrac{r}{\sigma^2}\right)\\ and \\\partial^2
\ell / \partial y^2 = \dfrac{1}{y^2}\left(1 +
\dfrac{r-1}{\sigma^2}\right)\\.

## Arguments

- distrib:

  A `Lognormal1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma2`.

## Value

A numeric vector.

## See also

[`lognormal1_distrib`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
