# Gamma Distribution Object (Mean-Variance Parameterization)

Creates a distribution object for the Gamma distribution parameterized
by mean (\\\mu\\) and variance (\\\sigma^2\\).

## Usage

``` r
gamma2_distrib(link_mu = log_link(), link_sigma2 = log_link())
```

## Arguments

- link_mu:

  A link function object for the mean parameter \\\mu\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

- link_sigma2:

  A link function object for the variance parameter \\\sigma^2\\.
  Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `Gamma2Distrib` (inheriting from
`continuous_distrib`) representing the Gamma distribution.

## Details

The Gamma distribution is given a mean/variance parameterization:
\\\mu\\ is the mean and \\\sigma^2\\ the variance. The standard shape
\\\alpha\\ and rate \\\lambda\\ are recovered as \$\$\alpha =
\dfrac{\mu^2}{\sigma^2}, \qquad \lambda = \dfrac{\mu}{\sigma^2}\$\$

**Probability density function:** \$\$f(y; \mu, \sigma^2) =
\dfrac{\lambda^\alpha}{\Gamma(\alpha)}\\ y^{\alpha-1} e^{-\lambda y},
\quad y \> 0\$\$

**Cumulative distribution function** (\\\gamma\\ the lower incomplete
gamma function): \$\$F(q; \mu, \sigma^2) = \dfrac{\gamma(\alpha, \lambda
q)}{\Gamma(\alpha)}\$\$

**Quantile function:** no closed form; the numerical inverse of the CDF.

**Score** (\\\psi\\ the digamma function): \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{-2\mu\psi(\alpha) + 2\mu\log\lambda + \mu +
2\mu\log y - y}{\sigma^2}\$\$ \$\$\dfrac{\partial \ell}{\partial
\sigma^2} = -\dfrac{\mu\left\[-\mu\psi(\alpha) + \mu + \mu(\log\lambda +
\log y) - y\right\]}{(\sigma^2)^2}\$\$

**Expected Hessian** (\\\psi_1\\ the trigamma function):
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
\dfrac{3\sigma^2 - 4\mu^2\psi_1(\alpha)}{(\sigma^2)^2}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu\\\partial
\sigma^2}\right\] = \dfrac{2\mu(\mu^2\psi_1(\alpha) -
\sigma^2)}{(\sigma^2)^3}\$\$ \$\$\mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial (\sigma^2)^2}\right\] =
-\dfrac{\mu^2(\mu^2\psi_1(\alpha) - \sigma^2)}{(\sigma^2)^4}\$\$ The
observed Hessian is available via
[`distrib_hessian.Gamma2Distrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gamma2Distrib.md).

**Moments:** mean \\\mu\\, variance \\\sigma^2\\, skewness
\\2\sqrt{\sigma^2}/\mu\\, excess kurtosis \\6\sigma^2/\mu^2\\.

**Parameter domains:**

- \\\mu \in (0, +\infty)\\

- \\\sigma^2 \in (0, +\infty)\\

Analytical third- and fourth-order derivatives
([`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md))
and response derivatives
([`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md))
are also available.

## See also

- [`distrib_pdf.Gamma2Distrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gamma2Distrib.md)
  for the probability density function.

- [`distrib_cdf.Gamma2Distrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gamma2Distrib.md)
  for the cumulative distribution function.

- [`distrib_quantile.Gamma2Distrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Gamma2Distrib.md)
  for the quantile function.

- [`distrib_rng.Gamma2Distrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Gamma2Distrib.md)
  for random number generation.

- [`distrib_gradient.Gamma2Distrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gamma2Distrib.md)
  for the analytical gradient.

- [`distrib_hessian.Gamma2Distrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gamma2Distrib.md)
  for the analytical observed Hessian.

- [`distrib_expected_hessian.Gamma2Distrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gamma2Distrib.md)
  for the analytical expected Hessian.

## Examples

``` r
d <- gamma2_distrib()
d@params
#> [1] "mu"     "sigma2"

theta <- list(mu = 2, sigma2 = 1)
distrib_pdf(d, c(0.5, 1, 2), theta)
#> [1] 0.1226265 0.3608941 0.3907336
distrib_gradient(d, c(0.5, 1, 2), theta)
#> $mu
#> [1] -3.5244707 -1.2518820  0.5207068
#> 
#> $sigma2
#> [1]  2.0244707  0.2518820 -0.5207068
#> 
```
