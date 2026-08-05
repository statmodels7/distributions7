# NB1 Negative Binomial Distribution Object

Creates a distribution object for the negative binomial whose variance
is **linear** in the mean, \\\operatorname{Var}(Y) = \mu(1+\theta)\\.

## Usage

``` r
negbin1_distrib(link_mu = log_link(), link_theta = log_link())
```

## Arguments

- link_mu:

  A link function object for \\\mu\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- link_theta:

  A link function object for \\\theta\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

## Value

An S7 object of class `NegBin1Distrib`.

## Details

Two negative binomials are in common use and they are **different
families**, not two parametrisations of one. Here the variance is
\\\mu(1+\theta)\\, growing in proportion to the mean, so the dispersion
relative to a Poisson is the same at every mean;
[`negbin2_distrib`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
has \\\mu + \mu^2/\theta\\, growing quadratically. Fitting one is not
fitting the other, and a likelihood ratio between them is not a test of
nested models.

The difference is visible in where the mean sits. The size is \\r =
\mu/\theta\\ and the success probability \\1/(1+\theta)\\, so \\\mu\\
appears **inside** the gamma functions, while in the quadratic form it
stays outside them and \\\theta\\ is the size.

**Probability mass function:** the negative binomial mass at size
\\\mu/\theta\\ and probability \\1/(1+\theta)\\.

**Score.** Writing \\P = \psi(y+r) - \psi(r) - \log(1+\theta)\\,
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{P}{\theta}, \qquad
\dfrac{\partial \ell}{\partial \theta} = -\dfrac{\mu}{\theta^2}P -
\dfrac{r}{1+\theta} + \dfrac{y}{\theta} - \dfrac{y}{1+\theta}\$\$ and
the Hessian is the same chain rule at second order.

**Expected information.** Every term carrying \\P\\ drops out, its
expectation vanishing by the first Bartlett identity, and only
\\\mathbb{E}\[\psi'(Y+r)\]\\ remains. That has no closed form and is
summed against the exact mass to a far-tail quantile.

**Parameter domains:**

- \\\mu \in (0, +\infty)\\

- \\\theta \in (0, +\infty)\\

As \\\theta \to 0\\ the family approaches the Poisson, as the quadratic
form does when \\\theta \to \infty\\.

## See also

[`negbin2_distrib`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
for the quadratic variance,
[`poisson_distrib`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md),
[`geometric_distrib`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)

## Examples

``` r
d <- negbin1_distrib()
d@params
#> [1] "mu"    "theta"

theta <- list(mu = 4, theta = 4)
distrib_pdf(d, 0:6, theta)
#> [1] 0.2000000 0.1600000 0.1280000 0.1024000 0.0819200 0.0655360 0.0524288
c(mean = mean(d, theta), variance = variance(d, theta))
#>     mean variance 
#>        4       20 

# the two negative binomials are different families: at the same (mu, theta)
# this one has variance mu(1+theta) = 20 and the other mu + mu^2/theta = 8
variance(negbin2_distrib(), list(mu = 4, theta = 4))
#> [1] 8
```
