# Poisson-Inverse Gaussian Distribution Object

Creates a Poisson-inverse Gaussian distribution in its mean-dispersion
parametrization: mean \\\mu\\ and variance \\\mu + \sigma\mu^2\\,
gamlss's `PIG` (Rigby and Stasinopoulos, 2005). The family is the
Poisson mixed over an inverse Gaussian rate, an overdispersed count
model with a heavier tail than the negative binomial at the same
variance.

## Usage

``` r
pig1_distrib(link_mu = log_link(), link_sigma = log_link())
```

## Arguments

- link_mu:

  The link for \\\mu\\; defaults to `log_link()`.

- link_sigma:

  The link for \\\sigma\\; defaults to `log_link()`.

## Value

A `Pig1Distrib` object.

## Details

The mass function carries the modified Bessel function \\K\_{y-1/2}\\,
which at half-integer order is a finite sum; the log-likelihood and its
derivatives to fourth order are exact, computed by a compiled kernel
that carries a bivariate jet through the closed expression. The expected
information has no closed form and goes through the summation strategies
of
[`expected_derivative_methods`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).
For the parametrization in which \\\mu\\ and the second parameter are
orthogonal, see
[`pig2_distrib`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md).

## The distribution

\$\$P(Y=y) =
\sqrt{\frac{2\alpha}{\pi}}\\\frac{\mu^{y}e^{1/\sigma}}{(\alpha\sigma)^{y}\\y!}\\K\_{y-1/2}(\alpha),
\qquad \alpha = \sqrt{\frac{1}{\sigma^{2}} + \frac{2\mu}{\sigma}}\$\$ on
\\y \in \\0, 1, \dots\\\\.

\$\$\mathbb{E}\[Y\] = \mu, \qquad \operatorname{Var}(Y) = \mu +
\sigma\mu^{2}\$\$

## References

Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive
models for location, scale and shape. *Applied Statistics* 54(3),
507–554.

Dean, C., Lawless, J. F., and Willmot, G. E. (1989). A mixed
Poisson-inverse-Gaussian regression model. *Canadian Journal of
Statistics* 17(2), 171–181.

## See also

[`pig2_distrib`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md),
[`negbin2_distrib`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)

## Examples

``` r
d <- pig1_distrib()
theta <- list(mu = 3, sigma = 0.8)
distrib_pdf(d, 0:5, theta)
#> [1] 0.17197629 0.21422781 0.17775288 0.12895666 0.08968701 0.06196187
mean(d, theta)
#> [1] 3
variance(d, theta)
#> [1] 10.2
```
