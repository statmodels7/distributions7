# Gaussian Distribution Object (Standard Deviation Parameterization)

Creates a distribution object for the Gaussian distribution
parameterized by mean (\\\mu\\) and standard deviation (\\\sigma\\).

## Usage

``` r
gaussian_distrib(link_mu = identity_link(), link_sigma = log_link())
```

## Arguments

- link_mu:

  A link function object for the location parameter \\\mu\\. Defaults to
  [`identity_link`](https://rdrr.io/pkg/linkfunctions7/man/identity_link.html).

- link_sigma:

  A link function object for the scale parameter \\\sigma\\. Defaults to
  [`log_link`](https://rdrr.io/pkg/linkfunctions7/man/log_link.html) to
  ensure positivity.

## Value

An S7 object of class `GaussianDistrib` (inheriting from
`continuous_distrib`) representing the Gaussian distribution.

## Details

The Gaussian (Normal) distribution is parameterized by its mean \\\mu\\
and standard deviation \\\sigma\\.

**Probability density function:** \$\$f(y; \mu, \sigma) =
\dfrac{1}{\sqrt{2\pi}\\\sigma}
\exp\left\\-\dfrac{1}{2}\left(\dfrac{y-\mu}{\sigma}\right)^2\right\\\$\$

**Cumulative distribution function** (\\\Phi\\ the standard normal CDF):
\$\$F(q; \mu, \sigma) = \Phi\left(\dfrac{q-\mu}{\sigma}\right)\$\$

**Quantile function:** \$\$Q(p; \mu, \sigma) = \mu +
\sigma\\\Phi^{-1}(p)\$\$

**Score** (gradient of the log-density): \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{y - \mu}{\sigma^2}, \qquad \dfrac{\partial
\ell}{\partial \sigma} = \dfrac{(y - \mu)^2 - \sigma^2}{\sigma^3}\$\$

**Observed Hessian:** \$\$\dfrac{\partial^2 \ell}{\partial \mu^2} =
-\dfrac{1}{\sigma^2}, \quad \dfrac{\partial^2 \ell}{\partial \sigma^2} =
\dfrac{\sigma^2 - 3(y-\mu)^2}{\sigma^4}, \quad \dfrac{\partial^2
\ell}{\partial \mu\\\partial \sigma} = -\dfrac{2(y-\mu)}{\sigma^3}\$\$

**Expected Hessian** (negative Fisher information):
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\sigma^2}, \quad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \sigma^2}\right\] = -\dfrac{2}{\sigma^2}, \quad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu\\\partial
\sigma}\right\] = 0\$\$

**Moments:** mean \\\mu\\, variance \\\sigma^2\\, skewness 0, excess
kurtosis 0.

**Parameter domains:**

- \\\mu \in (-\infty, +\infty)\\

- \\\sigma \in (0, +\infty)\\

Analytical third- and fourth-order derivatives
([`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md))
and derivatives with respect to the response
([`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md))
are also available.

## See also

- [`distrib_pdf.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GaussianDistrib.md)
  for the density function.

- [`distrib_cdf.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GaussianDistrib.md)
  for the cumulative distribution function.

- [`distrib_quantile.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GaussianDistrib.md)
  for the quantile function.

- [`distrib_rng.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.GaussianDistrib.md)
  for random number generation.

- [`distrib_gradient.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GaussianDistrib.md)
  for the analytical gradient.

- [`distrib_hessian.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GaussianDistrib.md)
  for the analytical observed Hessian.

- [`distrib_expected_hessian.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GaussianDistrib.md)
  for the analytical expected Hessian.
