# Geometric Distribution Object

Creates a distribution object for the geometric distribution on \\\\0,
1, 2, \dots\\\\, parametrised by its mean \\\mu\\.

## Usage

``` r
geometric_distrib(link_mu = log_link())
```

## Arguments

- link_mu:

  A link function object for \\\mu\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `GeometricDistrib`.

## Details

The distribution counts the failures before the first success, so the
support includes zero and the success probability is \\p = 1/(1+\mu)\\.

**Probability mass function:** \$\$P(Y = y; \mu) =
\dfrac{1}{1+\mu}\left(\dfrac{\mu}{1+\mu}\right)^{y}\$\$

**Score, observed and expected Hessian:** \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{y-\mu}{\mu(1+\mu)}, \qquad
\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} +
\dfrac{y+1}{(1+\mu)^2}, \qquad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu^2}\right\] = -\dfrac{1}{\mu(1+\mu)}\$\$

Every order is the same expression evaluated at \\\mu\\ and at
\\1+\mu\\, \$\$\ell^{(j)} = (-1)^{j-1}(j-1)!\left(\dfrac{y}{\mu^{j}} -
\dfrac{y+1}{(1+\mu)^{j}}\right), \qquad \mathbb{E}\[\ell^{(j)}\] =
(-1)^{j-1}(j-1)!\left(\mu^{1-j} - (1+\mu)^{1-j}\right)\$\$ so the
expected orders are closed form and vanish at \\j = 1\\.

**Moments:** mean \\\mu\\, variance \\\mu(1+\mu)\\, so the family is
overdispersed relative to the Poisson at every mean.

**Parameter domains:**

- \\\mu \in (0, +\infty)\\

The family is the negative binomial at \\\theta = 1\\, so
`fixed(negbin2_distrib(), theta = 1)` describes the same law and is used
in the tests as an independent implementation.

## See also

[`negbin2_distrib`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md),
[`poisson_distrib`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md),
[`exponential_distrib`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)

## Examples

``` r
d <- geometric_distrib()
d@params
#> [1] "mu"

theta <- list(mu = 2)
distrib_pdf(d, 0:4, theta)
#> [1] 0.33333333 0.22222222 0.14814815 0.09876543 0.06584362
distrib_gradient(d, 0:4, theta)
#> $mu
#> [1] -0.3333333 -0.1666667  0.0000000  0.1666667  0.3333333
#> 
c(mean = mean(d, theta), variance = variance(d, theta))
#>     mean variance 
#>        2        6 
```
