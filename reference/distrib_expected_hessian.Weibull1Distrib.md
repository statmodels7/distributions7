# Weibull Analytical Expected Hessian

Closed form. Under the model \\u = (Y/\mu)^{\sigma}\\ is standard
exponential whatever the parameters are, so every expectation needed is
a derivative of \\\Gamma\\ at 2: \\E\[u\] = 1\\, \\E\[u\log u\] = 1 -
\gamma\\ and \\E\[u(\log u)^2\] = (1-\gamma)^2 + \pi^2/6 - 1\\. Hence
\$\$E\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{\sigma^2}{\mu^2}, \qquad E\left\[\dfrac{\partial^2
\ell}{\partial \sigma^2}\right\] =
-\dfrac{1}{\sigma^2}\left\\(1-\gamma)^2 + \dfrac{\pi^2}{6}\right\\,\$\$
\$\$E\left\[\dfrac{\partial^2 \ell}{\partial \mu \\ \partial
\sigma}\right\] = \dfrac{1 - \gamma}{\mu},\$\$ with \\\gamma\\ the
Euler-Mascheroni constant. Because the closed form exists, the `approx`
argument is ignored.

## Arguments

- distrib:

  A `Weibull1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- scale:

  Either `"parameter"` or `"link"`.

- approx:

  Ignored; the expectation is exact.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list of expected second derivatives.

## See also

[`weibull1_distrib`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
