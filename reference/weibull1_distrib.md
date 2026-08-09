# Weibull Distribution Object

Creates a distribution object for the Weibull distribution, parametrized
by a scale \\\mu\\ and a shape \\\sigma\\, both positive.

## Usage

``` r
weibull1_distrib(link_mu = log_link(), link_sigma = log_link())
```

## Arguments

- link_mu:

  A link function object for the scale \\\mu\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- link_sigma:

  A link function object for the shape \\\sigma\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

## Value

An S7 object of class
[`Weibull1Distrib`](https://statmodels7.github.io/distributions7/reference/Weibull1Distrib.md)
(inheriting from `continuous_distrib`).

## Details

**Parametrization.** \\\mu\\ is the **scale** and not the mean. The mean
is \\\mu\\\Gamma(1 + 1/\sigma)\\, which involves the shape, so a mean
parametrization would make every derivative a derivative of the gamma
function and its inverse. The scale-shape form keeps the whole family
elementary, and
[`mean.Weibull1Distrib`](https://statmodels7.github.io/distributions7/reference/mean.Weibull1Distrib.md)
reports the mean. This is the parametrization of `WEI` in gamlss.

**Probability density function:** \$\$f(y; \mu, \sigma) =
\dfrac{\sigma}{\mu} \left(\dfrac{y}{\mu}\right)^{\sigma - 1}
\exp\left\\-\left(\dfrac{y}{\mu}\right)^{\sigma}\right\\, \qquad y \>
0\$\$

**Cumulative distribution function:** \$\$F(q; \mu, \sigma) = 1 -
\exp\left\\-(q/\mu)^{\sigma}\right\\\$\$

**Quantile function:** \$\$Q(p; \mu, \sigma) =
\mu\left\\-\log(1-p)\right\\^{1/\sigma}\$\$

**Score**, with \\u = (y/\mu)^{\sigma}\\ and \\z = y/\mu\\:
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{\sigma}{\mu}(u - 1),
\qquad \dfrac{\partial \ell}{\partial \sigma} = \dfrac{1}{\sigma} + (1 -
u)\log z\$\$

**Observed Hessian:** \$\$\dfrac{\partial^2 \ell}{\partial \mu^2} =
\dfrac{\sigma}{\mu^2}\left\\1 - (1 + \sigma) u\right\\, \quad
\dfrac{\partial^2 \ell}{\partial \sigma^2} = -\dfrac{1}{\sigma^2} - u
(\log z)^2, \quad \dfrac{\partial^2 \ell}{\partial \mu \\ \partial
\sigma} = \dfrac{u - 1 + \sigma u \log z}{\mu}\$\$

**Expected Hessian:** see
[`distrib_expected_hessian.Weibull1Distrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Weibull1Distrib.md).
The substitution \\u \sim \mathrm{Exp}(1)\\ turns every expectation into
a derivative of \\\Gamma\\ at 2.

**Moments.** With \\g_k = \Gamma(1 + k/\sigma)\\, the mean is \\\mu
g_1\\ and the variance \\\mu^2 (g_2 - g_1^2)\\; the skewness and the
excess kurtosis follow from \\g_3\\ and \\g_4\\ and do not depend on
\\\mu\\.

**Special cases.** \\\sigma = 1\\ is the exponential distribution with
mean \\\mu\\, and \\\sigma = 2\\ the Rayleigh distribution. The hazard
is increasing for \\\sigma \> 1\\ and decreasing for \\\sigma \< 1\\,
which is what the family is used for.

**Higher orders.** Third and fourth derivatives are closed form,
observed and expected: with \\u = (y/\mu)^{\sigma}\\ and \\L =
\log(y/\mu)\\, every derivative is a polynomial in \\u\\ and \\Lu\\, and
every expectation is a derivative of \\\Gamma\\ at 2.

**Parameter Domains:**

- \\\mu \in (0, +\infty)\\

- \\\sigma \in (0, +\infty)\\

## References

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994). *Continuous
Univariate Distributions, Volume 1*, 2nd edition, chapter 21. Wiley.

## See also

[`weibull3_distrib`](https://statmodels7.github.io/distributions7/reference/weibull3_distrib.md),
[`gumbel_distrib`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)

## Examples

``` r
d <- weibull1_distrib()
d@params
#> [1] "mu"    "sigma"

theta <- list(mu = 2, sigma = 1.5)
distrib_pdf(d, c(0.5, 1, 2), theta)
#> [1] 0.3309363 0.3723917 0.2759096
distrib_gradient(d, c(0.5, 1, 2), theta)
#> $mu
#> [1] -0.656250 -0.484835  0.000000
#> 
#> $sigma
#> [1] -0.5463409  0.2185840  0.6666667
#> 

# the scale is not the mean
c(scale = theta$mu, mean = mean(d, theta))
#>    scale     mean 
#> 2.000000 1.805491 

# shape 1 is the exponential distribution
max(abs(distrib_pdf(d, c(0.5, 1, 2), list(mu = 2, sigma = 1)) -
        stats::dexp(c(0.5, 1, 2), rate = 1 / 2)))
#> [1] 0
```
