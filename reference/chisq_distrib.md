# Chi-Squared Distribution Object

Creates a distribution object for the chi-squared distribution,
parametrized by its mean \\\mu\\, which is the degrees of freedom.

## Usage

``` r
chisq_distrib(link_mu = log_link())
```

## Arguments

- link_mu:

  A link function object for \\\mu\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `ChisqDistrib`.

## Details

The degrees of freedom are treated as a continuous positive parameter,
which is what makes the family estimable; the mean is \\\mu\\ and the
variance \\2\mu\\.

**Density:** \$\$f(y; \mu) = \dfrac{y^{\mu/2 - 1}
e^{-y/2}}{2^{\mu/2}\Gamma(\mu/2)}\$\$

**Score and information:** \$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{\log y - \log 2 - \psi(\mu/2)}{2}, \qquad \dfrac{\partial^2
\ell}{\partial \mu^2} = -\dfrac{\psi'(\mu/2)}{4}\$\$

The family is a one-parameter exponential family in \\\log y\\, so from
the second order on the derivatives do not involve the response at all,
\$\$\ell^{(k)} = -\dfrac{\psi^{(k-2)}(\mu/2)}{2^{k}}, \qquad k \ge
2.\$\$ On the parameter scale the observed information is therefore
exactly the expected information, and the same holds at third and fourth
order: there is nothing to average. \\\mathbb{E}\[\log y\] =
\psi(\mu/2) + \log 2\\ is what makes the score have mean zero.

That coincidence does not carry to the scale a fit optimizes on. The
second-order chain rule of
[`distrib_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
adds a term \\h''(\eta)\\\partial\ell/\partial\mu\\ to the link-scale
Hessian, and the expected version drops it because the score has mean
zero, while a sample does not. Fisher scoring and Newton's method
therefore take different steps here, and agree at the optimum, where the
summed score vanishes.

**Moments:** mean \\\mu\\, variance \\2\mu\\, skewness
\\2\sqrt{2/\mu}\\, excess kurtosis \\12/\mu\\.

**Parameter domains:**

- \\\mu \in (0, +\infty)\\

The family is a Gamma with shape \\\mu/2\\ and scale 2, but it is
**not** a Gamma with a fixed parameter: this package writes the Gamma in
\\(\mu, \sigma^2)\\, and a scale of 2 is the relation \\\sigma^2 =
2\mu\\ between two parameters rather than a value one of them can be
held at.

## See also

[`gamma2_distrib`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md),
[`exponential_distrib`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)

## Examples

``` r
d <- chisq_distrib()
d@params
#> [1] "mu"

theta <- list(mu = 4)
distrib_pdf(d, c(1, 4, 9), theta)
#> [1] 0.15163266 0.13533528 0.02499524
distrib_gradient(d, c(1, 4, 9), theta)
#> $mu
#> [1] -0.5579658  0.1351814  0.5406465
#> 

# the observed and expected information coincide exactly
distrib_hessian(d, c(1, 4, 9), theta)
#> $mu_mu
#> [1] -0.1612335 -0.1612335 -0.1612335
#> 
distrib_expected_hessian(d, c(1, 4, 9), theta)
#> $mu_mu
#> [1] -0.1612335 -0.1612335 -0.1612335
#> 
```
