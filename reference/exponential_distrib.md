# Exponential Distribution Object

Creates a distribution object for the exponential distribution
parametrised by its mean \\\mu\\.

## Usage

``` r
exponential_distrib(link_mu = log_link())
```

## Arguments

- link_mu:

  A link function object for \\\mu\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `ExponentialDistrib`.

## Details

The exponential distribution models waiting times on \\y \> 0\\, with
mean \\\mu\\ and variance \\\mu^2\\: the coefficient of variation is
one, and fixing it is what distinguishes the family from the Gamma.

**Density:** \$\$f(y; \mu) = \dfrac{1}{\mu} e^{-y/\mu}\$\$

**Distribution function:** \$\$F(q; \mu) = 1 - e^{-q/\mu}\$\$

**Score, observed and expected Hessian:** \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{y - \mu}{\mu^2}, \qquad \dfrac{\partial^2
\ell}{\partial \mu^2} = \dfrac{1}{\mu^2} - \dfrac{2y}{\mu^3}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\mu^2}\$\$

Every order follows the same pattern, the log-density being a logarithm
plus a reciprocal: \$\$\ell^{(k)} = \dfrac{(-1)^k (k-1)!}{\mu^k} +
\dfrac{(-1)^{k+1} k!\\ y}{\mu^{k+1}}, \qquad \mathbb{E}\[\ell^{(k)}\] =
\dfrac{(-1)^k (k-1)! (1-k)}{\mu^k}\$\$ so the expected orders are closed
form as well, and vanish at \\k = 1\\.

**Moments:** mean \\\mu\\, variance \\\mu^2\\, skewness 2, excess
kurtosis 6.

**Parameter domains:**

- \\\mu \in (0, +\infty)\\

The family is the Weibull with unit shape, so
`fixed(weibull_distrib(), sigma = 1)` describes the same law and is used
in the tests as an independent implementation. It is **not** a Gamma
with a fixed parameter: this package writes the Gamma in \\(\mu,
\sigma^2)\\, whose shape is \\\mu^2/\sigma^2\\, so unit shape is the
relation \\\sigma^2 = \mu^2\\ between two parameters rather than a value
one of them can be held at.

## See also

[`gamma_distrib`](https://statmodels7.github.io/distributions7/reference/gamma_distrib.md),
[`weibull_distrib`](https://statmodels7.github.io/distributions7/reference/weibull_distrib.md),
[`geometric_distrib`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)

## Examples

``` r
d <- exponential_distrib()
d@params
#> [1] "mu"

theta <- list(mu = 2)
distrib_pdf(d, c(0.5, 1, 3), theta)
#> [1] 0.3894004 0.3032653 0.1115651
distrib_gradient(d, c(0.5, 1, 3), theta)
#> $mu
#> [1] -0.375 -0.250  0.250
#> 
c(mean = mean(d, theta), variance = variance(d, theta))
#>     mean variance 
#>        2        4 
```
