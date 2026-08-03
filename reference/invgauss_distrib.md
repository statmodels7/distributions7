# Inverse-Gaussian Distribution Object (Mean-Dispersion Parameterization)

Creates a distribution object for the Inverse-Gaussian distribution
parameterized by mean (\\\mu\\) and dispersion (\\\phi\\).

## Usage

``` r
invgauss_distrib(link_mu = log_link(), link_phi = log_link())
```

## Arguments

- link_mu:

  A link function object for the mean parameter \\\mu\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

- link_phi:

  A link function object for the dispersion parameter \\\phi\\. Defaults
  to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `InvGaussDistrib` (inheriting from
`continuous_distrib`) representing the Inverse-Gaussian distribution.

## Details

The Inverse-Gaussian distribution is given a mean/dispersion
parameterization, with mean \\\mu\\ and dispersion \\\phi\\.

**Probability density function:** \$\$f(y; \mu, \phi) =
\sqrt{\dfrac{1}{2\pi\phi y^3}} \exp\left\\-\dfrac{(y-\mu)^2}{2\phi\mu^2
y}\right\\, \quad y \> 0\$\$

**Cumulative distribution function** (\\\Phi\\ the standard normal CDF):
\$\$F(q; \mu, \phi) = \Phi\\\left(\sqrt{\tfrac{1}{\phi
q}}\left(\tfrac{q}{\mu}-1\right)\right) + e^{2/(\phi\mu)}\\
\Phi\\\left(-\sqrt{\tfrac{1}{\phi
q}}\left(\tfrac{q}{\mu}+1\right)\right)\$\$

**Quantile function:** no closed form; the numerical inverse of the CDF.

**Score:** \$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y -
\mu}{\phi\mu^3}, \qquad \dfrac{\partial \ell}{\partial \phi} =
\dfrac{(y - \mu)^2 - y\mu^2\phi}{2y\phi^2\mu^2}\$\$

**Observed Hessian:** \$\$\dfrac{\partial^2 \ell}{\partial \mu^2} =
-\dfrac{3y - 2\mu}{\phi\mu^4}, \quad \dfrac{\partial^2 \ell}{\partial
\phi^2} = \dfrac{\phi - 2(y-\mu)^2/(\mu^2 y)}{2\phi^3}, \quad
\dfrac{\partial^2 \ell}{\partial \mu\\\partial \phi} = -\dfrac{y -
\mu}{\phi^2\mu^3}\$\$

**Expected Hessian:** \$\$\mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu^2}\right\] = -\dfrac{1}{\phi\mu^3}, \quad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \phi^2}\right\] =
-\dfrac{1}{2\phi^2}, \quad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu\\\partial \phi}\right\] = 0\$\$

**Moments:** mean \\\mu\\, variance \\\phi\mu^3\\, skewness
\\3\sqrt{\phi\mu}\\, excess kurtosis \\15\phi\mu\\.

**Parameter domains:**

- \\\mu \in (0, +\infty)\\

- \\\phi \in (0, +\infty)\\

Analytical third- and fourth-order derivatives
([`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md))
and response derivatives
([`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md))
are also available.

## See also

- [`distrib_pdf.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.InvGaussDistrib.md)
  for the probability density function.

- [`distrib_cdf.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.InvGaussDistrib.md)
  for the cumulative distribution function.

- [`distrib_quantile.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.InvGaussDistrib.md)
  for the quantile function.

- [`distrib_rng.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.InvGaussDistrib.md)
  for random number generation.

- [`distrib_gradient.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.InvGaussDistrib.md)
  for the analytical gradient.

- [`distrib_hessian.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.InvGaussDistrib.md)
  for the analytical observed Hessian.

- [`distrib_expected_hessian.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.InvGaussDistrib.md)
  for the analytical expected Hessian.
