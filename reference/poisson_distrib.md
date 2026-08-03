# Poisson Distribution Object

Creates a distribution object for the Poisson distribution parameterized
by the mean parameter \\\mu\\.

## Usage

``` r
poisson_distrib(link_mu = log_link())
```

## Arguments

- link_mu:

  A link function object for the mean parameter \\\mu\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `PoissonDistrib` (inheriting from
`discrete_distrib`) representing the Poisson distribution.

## Details

The Poisson distribution models counts \\y \in \\0, 1, 2, \dots\\\\ with
mean (and variance) \\\mu\\.

**Probability mass function:** \$\$P(Y=y; \mu) = \dfrac{\mu^y
e^{-\mu}}{y!}\$\$

**Cumulative distribution function:** \$\$F(q; \mu) =
\sum\_{k=0}^{\lfloor q \rfloor} \dfrac{\mu^k e^{-\mu}}{k!}\$\$

**Quantile function:** the generalized inverse \\Q(p) = \min\\y \in
\mathbb{N}\_0 : F(y) \ge p\\\\.

**Score, observed and expected Hessian:** \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{y - \mu}{\mu}, \qquad \dfrac{\partial^2
\ell}{\partial \mu^2} = -\dfrac{y}{\mu^2}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\mu}\$\$

**Moments:** mean \\\mu\\, variance \\\mu\\, skewness \\1/\sqrt{\mu}\\,
excess kurtosis \\1/\mu\\.

**Parameter domains:**

- \\\mu \in (0, +\infty)\\

Analytical third- and fourth-order derivatives
([`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md))
are also available.

## See also

- [`distrib_pdf.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.PoissonDistrib.md)
  for the probability mass function.

- [`distrib_cdf.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.PoissonDistrib.md)
  for the cumulative distribution function.

- [`distrib_quantile.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.PoissonDistrib.md)
  for the quantile function.

- [`distrib_rng.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.PoissonDistrib.md)
  for random number generation.

- [`distrib_gradient.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PoissonDistrib.md)
  for the analytical gradient.

- [`distrib_hessian.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PoissonDistrib.md)
  for the analytical observed Hessian.

- [`distrib_expected_hessian.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.PoissonDistrib.md)
  for the analytical expected Hessian.
