# Gumbel Analytical Expected Hessian

Closed form. Under the model \\w = e^{-Z}\\ is standard exponential, so
every expectation needed is a derivative of \\\Gamma\\ at 2: \\E\[w\] =
1\\, \\E\[w \log w\] = 1 - \gamma\\ and \\E\[w(\log w)^2\] =
(1-\gamma)^2 + \pi^2/6 - 1\\, together with \\E\[Z\] = \gamma\\. Hence
\$\$E\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\sigma^2}, \qquad E\left\[\dfrac{\partial^2 \ell}{\partial
\mu \\ \partial \sigma}\right\] = \dfrac{1 - \gamma}{\sigma^2},\$\$
\$\$E\left\[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right\] =
-\dfrac{(1-\gamma)^2 + \pi^2/6}{\sigma^2},\$\$ with \\\gamma\\ the
Euler-Mascheroni constant. Because the closed form exists, the `approx`
argument is ignored.

## Arguments

- distrib:

  A `GumbelDistrib` object.

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

[`gumbel_distrib`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
