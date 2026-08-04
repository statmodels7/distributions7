# Bernoulli Distribution Object

Creates a distribution object for the Bernoulli distribution
parameterized by the probability of success \\\mu\\.

## Usage

``` r
bernoulli_distrib(link_mu = logit_link())
```

## Arguments

- link_mu:

  A link function object for the mean parameter \\\mu\\ (probability).
  Defaults to
  [`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html)
  to ensure the parameter stays within (0, 1).

## Value

An S7 object of class `BernoulliDistrib` (inheriting from
`discrete_distrib`) representing the Bernoulli distribution.

## Details

The Bernoulli distribution models a binary outcome \\y \in \\0, 1\\\\
with success probability \\\mu\\.

**Probability mass function:** \$\$P(Y=y; \mu) = \mu^y (1-\mu)^{1-y},
\quad y \in \\0, 1\\\$\$

**Cumulative distribution function:** \$\$F(q; \mu) = \begin{cases} 0 &
q \< 0 \\ 1-\mu & 0 \le q \< 1 \\ 1 & q \ge 1 \end{cases}\$\$

**Quantile function** (generalized inverse): \$\$Q(p; \mu) =
\begin{cases} 0 & p \le 1-\mu \\ 1 & p \> 1-\mu \end{cases}\$\$

**Score, observed and expected Hessian:** \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{y - \mu}{\mu(1-\mu)}, \qquad
\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} -
\dfrac{1-y}{(1-\mu)^2}, \qquad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu^2}\right\] = -\dfrac{1}{\mu(1-\mu)}\$\$

**Moments:** mean \\\mu\\, variance \\\mu(1-\mu)\\, skewness
\\(1-2\mu)/\sqrt{\mu(1-\mu)}\\, excess kurtosis \\(1 -
6\mu(1-\mu))/(\mu(1-\mu))\\.

**Parameter domains:**

- \\\mu \in (0, 1)\\

Analytical third- and fourth-order derivatives
([`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md))
are also available.

## See also

- [`distrib_pdf.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BernoulliDistrib.md)
  for the probability mass function.

- [`distrib_cdf.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.BernoulliDistrib.md)
  for the cumulative distribution function.

- [`distrib_quantile.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.BernoulliDistrib.md)
  for the quantile function.

- [`distrib_rng.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.BernoulliDistrib.md)
  for random number generation.

- [`distrib_gradient.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BernoulliDistrib.md)
  for the analytical gradient.

- [`distrib_hessian.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BernoulliDistrib.md)
  for the analytical observed Hessian.

- [`distrib_expected_hessian.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BernoulliDistrib.md)
  for the analytical expected Hessian.

## Examples

``` r
d <- bernoulli_distrib()
d@params
#> [1] "mu"

theta <- list(mu = 0.3)
distrib_pdf(d, c(0, 1), theta)
#> [1] 0.7 0.3
distrib_gradient(d, c(0, 1), theta)
#> $mu
#> [1] -1.428571  3.333333
#> 
```
