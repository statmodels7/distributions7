# Logistic Distribution Object

Creates a distribution object for the Logistic distribution
parameterized by location (\\\mu\\) and scale (\\\sigma\\).

## Usage

``` r
logistic_distrib(link_mu = identity_link(), link_sigma = log_link())
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

An S7 object of class `LogisticDistrib` (inheriting from
`continuous_distrib`) representing the Logistic distribution.

## Details

The Logistic distribution is a location-scale distribution with location
\\\mu\\ and scale \\\sigma\\.

**Probability density function:** \$\$f(y; \mu, \sigma) =
\dfrac{\exp\left(-\dfrac{y-\mu}{\sigma}\right)}{\sigma\left\[1 +
\exp\left(-\dfrac{y-\mu}{\sigma}\right)\right\]^2}\$\$

**Cumulative distribution function:** \$\$F(q; \mu, \sigma) =
\dfrac{1}{1 + \exp\left(-\dfrac{q-\mu}{\sigma}\right)}\$\$

**Quantile function:** \$\$Q(p; \mu, \sigma) = \mu + \sigma
\log\left(\dfrac{p}{1-p}\right)\$\$

**Score:** \$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{1}{\sigma}
\tanh\left(\dfrac{y-\mu}{2\sigma}\right), \qquad \dfrac{\partial
\ell}{\partial \sigma} = -\dfrac{1}{\sigma}\left\[1 -
\dfrac{y-\mu}{\sigma}\tanh\left(\dfrac{y-\mu}{2\sigma}\right)\right\]\$\$

**Expected Hessian:** \$\$\mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu^2}\right\] = -\dfrac{1}{3\sigma^2}, \quad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right\] =
-\dfrac{3+\pi^2}{9\sigma^2}, \quad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu\\\partial \sigma}\right\] = 0\$\$ The observed
Hessian is available via
[`distrib_hessian.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LogisticDistrib.md).

**Moments:** mean \\\mu\\, variance \\\pi^2\sigma^2/3\\, skewness 0,
excess kurtosis \\6/5\\.

**Parameter domains:**

- \\\mu \in (-\infty, +\infty)\\

- \\\sigma \in (0, +\infty)\\

Response derivatives
([`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md))
and observed third- and fourth-order parameter derivatives
([`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md))
are available in closed form. The corresponding expected derivatives are
not: they are obtained through
[`expected_derivative_methods`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).
Seven of the nine components are nonetheless known exactly, and the
location-scale structure of the family forbids any of them from
depending on \\\mu\\: \$\$\mathbb{E}\left\[\dfrac{\partial^3
\ell}{\partial \mu^2 \partial \sigma}\right\] = \dfrac{1}{2\sigma^3},
\qquad \mathbb{E}\left\[\dfrac{\partial^3 \ell}{\partial
\sigma^3}\right\] = \dfrac{\pi^2+2}{2\sigma^3}, \qquad
\mathbb{E}\left\[\dfrac{\partial^4 \ell}{\partial \mu^4}\right\] =
\dfrac{1}{15\sigma^4},\$\$ the remaining ones in that list being zero by
symmetry.

## See also

- [`distrib_pdf.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.LogisticDistrib.md)
  for the probability density function.

- [`distrib_cdf.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.LogisticDistrib.md)
  for the cumulative distribution function.

- [`distrib_quantile.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.LogisticDistrib.md)
  for the quantile function.

- [`distrib_rng.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.LogisticDistrib.md)
  for random number generation.

- [`distrib_gradient.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.LogisticDistrib.md)
  for the analytical gradient.

- [`distrib_hessian.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LogisticDistrib.md)
  for the analytical observed Hessian.

- [`distrib_expected_hessian.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LogisticDistrib.md)
  for the analytical expected Hessian.

- [`distrib_deriv3.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.LogisticDistrib.md)
  and
  [`distrib_deriv4.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.LogisticDistrib.md)
  for the observed higher-order derivatives.
