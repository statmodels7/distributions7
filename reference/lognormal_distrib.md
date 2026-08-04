# Lognormal Distribution Object (Log-Scale Parameterization)

Creates a distribution object for the Lognormal distribution
parameterized by the mean (\\\mu\\) and the variance (\\\sigma^2\\) of
the log-transformed variable.

## Usage

``` r
lognormal_distrib(link_mu = identity_link(), link_sigma2 = log_link())
```

## Arguments

- link_mu:

  A link function object for the location parameter \\\mu\\. Defaults to
  [`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_sigma2:

  A link function object for the variance parameter \\\sigma^2\\.
  Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `LognormalDistrib` (inheriting from
`continuous_distrib`) representing the Lognormal distribution.

## Details

The Lognormal distribution describes a variable whose logarithm is
Gaussian with mean \\\mu\\ and variance \\\sigma^2\\ (both on the log
scale). Write \\r = \log y - \mu\\.

**Probability density function:** \$\$f(y; \mu, \sigma^2) =
\dfrac{1}{y\sqrt{2\pi\sigma^2}} \exp\left\\-\dfrac{(\log y -
\mu)^2}{2\sigma^2}\right\\, \quad y \> 0\$\$

**Cumulative distribution function** (\\\Phi\\ the standard normal CDF):
\$\$F(q; \mu, \sigma^2) = \Phi\left(\dfrac{\log q -
\mu}{\sqrt{\sigma^2}}\right)\$\$

**Quantile function:** \$\$Q(p; \mu, \sigma^2) = \exp\left(\mu +
\sqrt{\sigma^2}\\\Phi^{-1}(p)\right)\$\$

**Score:** \$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{r}{\sigma^2}, \qquad \dfrac{\partial \ell}{\partial \sigma^2} =
\dfrac{r^2 - \sigma^2}{2\sigma^4}\$\$

**Observed Hessian:** \$\$\dfrac{\partial^2 \ell}{\partial \mu^2} =
-\dfrac{1}{\sigma^2}, \quad \dfrac{\partial^2 \ell}{\partial
(\sigma^2)^2} = \dfrac{1}{2\sigma^4} - \dfrac{r^2}{\sigma^6}, \quad
\dfrac{\partial^2 \ell}{\partial \mu\\\partial \sigma^2} =
-\dfrac{r}{\sigma^4}\$\$

**Expected Hessian:** \$\$\mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu^2}\right\] = -\dfrac{1}{\sigma^2}, \quad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial (\sigma^2)^2}\right\]
= -\dfrac{1}{2\sigma^4}, \quad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu\\\partial \sigma^2}\right\] = 0\$\$

**Moments:** mean \\e^{\mu + \sigma^2/2}\\, variance
\\(e^{\sigma^2}-1)e^{2\mu+\sigma^2}\\, skewness
\\(e^{\sigma^2}+2)\sqrt{e^{\sigma^2}-1}\\, excess kurtosis
\\e^{4\sigma^2}+2e^{3\sigma^2}+3e^{2\sigma^2}-6\\.

**Parameter domains:**

- \\\mu \in (-\infty, +\infty)\\

- \\\sigma^2 \in (0, +\infty)\\

Analytical third- and fourth-order derivatives
([`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md))
and response derivatives
([`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md))
are also available.

## Examples

``` r
d <- lognormal_distrib()
d@params
#> [1] "mu"     "sigma2"

theta <- list(mu = 0, sigma2 = 1)
distrib_pdf(d, c(0.5, 1, 2), theta)
#> [1] 0.6274961 0.3989423 0.1568740
distrib_gradient(d, c(0.5, 1, 2), theta)
#> $mu
#> [1] -0.6931472  0.0000000  0.6931472
#> 
#> $sigma2
#> [1] -0.2597735 -0.5000000 -0.2597735
#> 
```
