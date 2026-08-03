# Negative Binomial Distribution Object (NB2)

Creates a distribution object for the Negative Binomial distribution
(NB2) parameterized by mean (\\\mu\\) and dispersion (\\\theta\\).

## Usage

``` r
negbin_distrib(link_mu = log_link(), link_theta = log_link())
```

## Arguments

- link_mu:

  A link function object for the mean parameter \\\mu\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- link_theta:

  A link function object for the dispersion parameter \\\theta\\.
  Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

## Value

An S7 object of class `NegBinDistrib` (inheriting from
`discrete_distrib`) representing the Negative Binomial distribution.

## Details

The Negative Binomial distribution (NB2) is given a mean/dispersion
parameterization, with mean \\\mu\\ and dispersion \\\theta\\ (smaller
\\\theta\\ means more overdispersion). Write \\s = \theta + \mu\\.

**Probability mass function:** \$\$P(Y=y; \mu, \theta) =
\dfrac{\Gamma(y+\theta)}{y!\\\Gamma(\theta)}
\left(\dfrac{\theta}{s}\right)^\theta \left(\dfrac{\mu}{s}\right)^y\$\$

**Cumulative distribution function:** \$\$F(q; \mu, \theta) =
\sum\_{k=0}^{\lfloor q \rfloor} P(Y=k; \mu, \theta)\$\$

**Quantile function:** the generalized inverse \\Q(p) = \min\\y \in
\mathbb{N}\_0 : F(y) \ge p\\\\.

**Score** (\\\psi\\ the digamma function): \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{\theta}{s}\left(\dfrac{y}{\mu} - 1\right),
\qquad \dfrac{\partial \ell}{\partial \theta} = \psi(y+\theta) -
\psi(\theta) + \log\left(\dfrac{\theta}{s}\right) + \dfrac{\mu -
y}{s}\$\$

**Observed Hessian** (\\\psi_1\\ the trigamma function):
\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{y+\theta}{s^2} -
\dfrac{y}{\mu^2}, \qquad \dfrac{\partial^2 \ell}{\partial \mu\\\partial
\theta} = \dfrac{y-\mu}{s^2}\$\$ \$\$\dfrac{\partial^2 \ell}{\partial
\theta^2} = \psi_1(y+\theta) - \psi_1(\theta) + \dfrac{\mu}{\theta s} +
\dfrac{y-\mu}{s^2}\$\$

**Expected Hessian:** \$\$\mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu^2}\right\] = -\dfrac{\theta}{\mu s}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \theta^2}\right\] =
\mathbb{E}\[\psi_1(Y+\theta)\] - \psi_1(\theta) + \dfrac{\mu}{\theta s},
\qquad \mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu\\\partial
\theta}\right\] = 0\$\$

**Moments:** mean \\\mu\\, variance \\\mu + \mu^2/\theta\\, skewness
\\(\theta + 2\mu)/\sqrt{\mu\theta(\theta+\mu)}\\, excess kurtosis
\\6/\theta + \theta/(\mu(\theta+\mu))\\.

**Parameter domains:**

- \\\mu \in (0, +\infty)\\

- \\\theta \in (0, +\infty)\\

Analytical third- and fourth-order derivatives
([`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md))
are also available.

## See also

- [`distrib_pdf.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.NegBinDistrib.md)
  for the probability mass function.

- [`distrib_cdf.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.NegBinDistrib.md)
  for the cumulative distribution function.

- [`distrib_quantile.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.NegBinDistrib.md)
  for the quantile function.

- [`distrib_rng.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.NegBinDistrib.md)
  for random number generation.

- [`distrib_gradient.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.NegBinDistrib.md)
  for the analytical gradient.

- [`distrib_hessian.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.NegBinDistrib.md)
  for the analytical observed Hessian.

- [`distrib_expected_hessian.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.NegBinDistrib.md)
  for the analytical expected Hessian.
