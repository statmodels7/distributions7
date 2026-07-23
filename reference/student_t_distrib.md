# Student's t Distribution Object (Location-Scale Parameterization)

Creates a distribution object for the Student's t distribution
parameterized by location (\\\mu\\), scale (\\\sigma\\), and degrees of
freedom (\\\nu\\).

## Usage

``` r
student_t_distrib(
  link_mu = identity_link(),
  link_sigma = log_link(),
  link_nu = log_link()
)
```

## Arguments

- link_mu:

  A link function object for the location parameter \\\mu\\. Defaults to
  [`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_sigma:

  A link function object for the scale parameter \\\sigma\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

- link_nu:

  A link function object for the degrees of freedom parameter \\\nu\\.
  Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `StudentTDistrib` (inheriting from
`continuous_distrib`) representing the Student's t distribution.

## Details

The (location-scale) Student's t distribution has location \\\mu\\,
scale \\\sigma\\ and degrees of freedom \\\nu\\. Write \\d =
\nu\sigma^2 + (y-\mu)^2\\.

**Probability density function:** \$\$f(y; \mu, \sigma, \nu) =
\dfrac{\Gamma\left(\dfrac{\nu+1}{2}\right)}{\sigma\sqrt{\nu\pi}\\\Gamma\left(\dfrac{\nu}{2}\right)}
\left(1 + \dfrac{(y-\mu)^2}{\nu\sigma^2}\right)^{-\dfrac{\nu+1}{2}}\$\$

**Cumulative distribution function** (\\T\_\nu\\ the standard t CDF):
\$\$F(q; \mu, \sigma, \nu) =
T\_\nu\\\left(\dfrac{q-\mu}{\sigma}\right)\$\$

**Quantile function:** \$\$Q(p; \mu, \sigma, \nu) = \mu +
\sigma\\T\_\nu^{-1}(p)\$\$

**Score** (\\\psi\\ the digamma function): \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{(\nu+1)(y-\mu)}{d}, \qquad \dfrac{\partial
\ell}{\partial \sigma} = \dfrac{\nu\left\[(y-\mu)^2 -
\sigma^2\right\]}{\sigma d}\$\$ \$\$\dfrac{\partial \ell}{\partial \nu}
= \dfrac{1}{2}\left\[ -\dfrac{1}{\nu} -
\psi\left(\dfrac{\nu}{2}\right) + \psi\left(\dfrac{\nu+1}{2}\right) +
\dfrac{(\nu+1)(y-\mu)^2}{\nu d} - \log\left(1 +
\dfrac{(y-\mu)^2}{\nu\sigma^2}\right) \right\]\$\$

**Expected Hessian** (\\\psi_1\\ the trigamma function; \\\mu\\ is
orthogonal to \\\sigma, \nu\\): \$\$\mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu^2}\right\] = -\dfrac{\nu+1}{\sigma^2(\nu+3)}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right\] =
-\dfrac{2\nu}{\sigma^2(\nu+3)}\$\$
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \nu^2}\right\] =
\dfrac{1}{4}\left\[\psi_1\left(\dfrac{\nu+1}{2}\right) -
\psi_1\left(\dfrac{\nu}{2}\right)\right\] +
\dfrac{\nu+5}{2\nu(\nu+1)(\nu+3)}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \sigma\\\partial
\nu}\right\] = \dfrac{2}{\sigma(\nu+1)(\nu+3)}\$\$ The observed Hessian
is available via
[`distrib_hessian.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.StudentTDistrib.md).

**Moments** (defined for \\\nu\\ large enough): mean \\\mu\\
(\\\nu\>1\\), variance \\\sigma^2\nu/(\nu-2)\\ (\\\nu\>2\\), skewness 0
(\\\nu\>3\\), excess kurtosis \\6/(\nu-4)\\ (\\\nu\>4\\).

**Parameter domains:**

- \\\mu \in (-\infty, +\infty)\\

- \\\sigma \in (0, +\infty)\\

- \\\nu \in (0, +\infty)\\

Analytical third- and fourth-order observed derivatives
([`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md);
the expected ones use the numerical fallback) and response derivatives
([`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md))
are also available.

## See also

- [`distrib_pdf.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.StudentTDistrib.md)
  for the probability density function.

- [`distrib_cdf.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.StudentTDistrib.md)
  for the cumulative distribution function.

- [`distrib_quantile.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.StudentTDistrib.md)
  for the quantile function.

- [`distrib_rng.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.StudentTDistrib.md)
  for random number generation.

- [`distrib_gradient.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.StudentTDistrib.md)
  for the analytical gradient.

- [`distrib_hessian.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.StudentTDistrib.md)
  for the analytical observed Hessian.

- [`distrib_expected_hessian.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.StudentTDistrib.md)
  for the analytical expected Hessian.
