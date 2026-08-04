# Cauchy Distribution Object

Creates a distribution object for the Cauchy distribution, parameterized
by location (\\\mu\\) and scale (\\\sigma\\).

## Usage

``` r
cauchy_distrib(link_mu = identity_link(), link_sigma = log_link())
```

## Arguments

- link_mu:

  A link function object for the location parameter \\\mu\\. Defaults to
  [`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_sigma:

  A link function object for the scale parameter \\\sigma\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `CauchyDistrib` (inheriting from
`continuous_distrib`) representing the Cauchy distribution.

## Details

The Cauchy distribution is a heavy-tailed location-scale distribution
with location \\\mu\\ and scale \\\sigma\\. Its moments (mean, variance,
and higher) are **undefined**.

**Probability density function:** \$\$f(y; \mu, \sigma) = \dfrac{1}{\pi
\sigma \left\[1 + \left(\dfrac{y-\mu}{\sigma}\right)^2\right\]}\$\$

**Cumulative distribution function:** \$\$F(q; \mu, \sigma) =
\dfrac{1}{2} +
\dfrac{1}{\pi}\arctan\left(\dfrac{q-\mu}{\sigma}\right)\$\$

**Quantile function:** \$\$Q(p; \mu, \sigma) = \mu + \sigma
\tan\left(\pi\left(p - \tfrac{1}{2}\right)\right)\$\$

**Score** (with \\d = \sigma^2 + (y-\mu)^2\\): \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{2(y-\mu)}{d}, \qquad \dfrac{\partial
\ell}{\partial \sigma} = \dfrac{(y-\mu)^2 - \sigma^2}{\sigma d}\$\$

**Observed Hessian:** \$\$\dfrac{\partial^2 \ell}{\partial \mu^2} =
\dfrac{2(y-\mu)^2 - 2\sigma^2}{d^2}, \quad \dfrac{\partial^2
\ell}{\partial \sigma^2} = \dfrac{\sigma^4 - 4\sigma^2 (y-\mu)^2 -
(y-\mu)^4}{\sigma^2 d^2}, \quad \dfrac{\partial^2 \ell}{\partial
\mu\\\partial \sigma} = -\dfrac{4\sigma(y-\mu)}{d^2}\$\$

**Expected Hessian:** \$\$\mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu^2}\right\] = \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \sigma^2}\right\] = -\dfrac{1}{2\sigma^2}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu\\\partial
\sigma}\right\] = 0\$\$

**Parameter domains:**

- \\\mu \in (-\infty, +\infty)\\

- \\\sigma \in (0, +\infty)\\

Analytical third- and fourth-order derivatives
([`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md))
and response derivatives
([`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md))
are also available.

## See also

- [`distrib_pdf.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.CauchyDistrib.md)
  for the density function.

- [`distrib_cdf.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.CauchyDistrib.md)
  for the cumulative distribution function.

- [`distrib_quantile.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.CauchyDistrib.md)
  for the quantile function.

- [`distrib_rng.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.CauchyDistrib.md)
  for random number generation.

- [`distrib_gradient.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.CauchyDistrib.md)
  for the analytical gradient.

- [`distrib_hessian.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.CauchyDistrib.md)
  for the analytical observed Hessian.

- [`distrib_expected_hessian.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.CauchyDistrib.md)
  for the analytical expected Hessian.

## Examples

``` r
d <- cauchy_distrib()
d@params
#> [1] "mu"    "sigma"

theta <- list(mu = 0, sigma = 1)
distrib_pdf(d, c(-1, 0, 1), theta)
#> [1] 0.1591549 0.3183099 0.1591549
distrib_gradient(d, c(-1, 0, 1), theta)
#> $mu
#> [1] -1  0  1
#> 
#> $sigma
#> [1]  0 -1  0
#> 
```
