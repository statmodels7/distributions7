# Gumbel Distribution Object

Creates a distribution object for the Gumbel (type I extreme value)
distribution, in the form for **maxima**, with location \\\mu\\ and
scale \\\sigma\\.

## Usage

``` r
gumbel_distrib(link_mu = identity_link(), link_sigma = log_link())
```

## Arguments

- link_mu:

  A link function object for the location \\\mu\\. Defaults to
  [`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_sigma:

  A link function object for the scale \\\sigma\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class
[`GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/GumbelDistrib.md)
(inheriting from `continuous_distrib`).

## Details

The Gumbel distribution is the limit law of the maximum of a sample from
a light-tailed distribution, which is what it is used for. It is a
location-scale family on the whole real line, skewed to the right with a
fixed shape: unlike the gaussian, its skewness and kurtosis are
constants and cannot be fitted.

**Probability density function**, with \\z = (y-\mu)/\sigma\\: \$\$f(y;
\mu, \sigma) = \dfrac{1}{\sigma}\exp\left\\-z - e^{-z}\right\\\$\$

**Cumulative distribution function:** \$\$F(q; \mu, \sigma) =
\exp\left\\-e^{-z}\right\\\$\$

**Quantile function:** \$\$Q(p; \mu, \sigma) = \mu - \sigma\log(-\log
p)\$\$

**Score**, with \\w = e^{-z}\\: \$\$\dfrac{\partial \ell}{\partial \mu}
= \dfrac{1 - w}{\sigma}, \qquad \dfrac{\partial \ell}{\partial \sigma} =
\dfrac{z(1 - w) - 1}{\sigma}\$\$

**Observed Hessian:** \$\$\dfrac{\partial^2 \ell}{\partial \mu^2} =
-\dfrac{w}{\sigma^2}, \quad \dfrac{\partial^2 \ell}{\partial \mu \\
\partial \sigma} = -\dfrac{1 - w + zw}{\sigma^2}, \quad
\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{1 - 2z + 2zw - z^2
w}{\sigma^2}\$\$

**Expected Hessian:** see
[`distrib_expected_hessian.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GumbelDistrib.md).
Note that \\E\[\partial^2 \ell / \partial \mu \\ \partial \sigma\\\]
does not vanish: the location and the scale are **not** orthogonal here,
unlike in a symmetric location-scale family, because the density is
skewed.

**Moments:** mean \\\mu + \gamma\sigma\\, variance \\\pi^2\sigma^2/6\\,
skewness \\12\sqrt{6}\\\zeta(3)/\pi^3 \approx 1.1395\\ and excess
kurtosis \\12/5\\, the last two free of both parameters.

**Relation to the Weibull.** If \\Y\\ is Gumbel then \\e^{-Y}\\ is
Weibull, so
[`weibull1_distrib`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
is this family on the log scale and reversed; the two share the
expectations that produce their information matrices.

**Minima.** The distribution of minima is the reflection: fit this
family to \\-Y\\, or use
[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)
with
[`affine_transform`](https://statmodels7.github.io/distributions7/reference/affine_transform.md)`(scale = -1)`.

**Higher orders.** Third and fourth derivatives are closed form,
observed and expected: with \\z = (y-\mu)/\sigma\\ and \\w = e^{-z}\\,
every derivative is a polynomial in \\z\\ and \\z^j w\\, and every
expectation is a derivative of \\\Gamma\\ at 2.

**Parameter Domains:**

- \\\mu \in (-\infty, +\infty)\\

- \\\sigma \in (0, +\infty)\\

## References

Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme
Values*, chapter 3. Springer.

## Examples

``` r
d <- gumbel_distrib()
d@params
#> [1] "mu"    "sigma"

theta <- list(mu = 0, sigma = 1)
distrib_pdf(d, c(-1, 0, 1), theta)
#> [1] 0.1793741 0.3678794 0.2546464
distrib_gradient(d, c(-1, 0, 1), theta)
#> $mu
#> [1] -1.7182818  0.0000000  0.6321206
#> 
#> $sigma
#> [1]  0.7182818 -1.0000000 -0.3678794
#> 

# the mean is shifted by Euler's constant, and the shape is fixed
c(mean = mean(d, theta), skewness = skewness(d, theta))
#>      mean  skewness 
#> 0.5772157 1.1395471 

# exp(-Y) is Weibull: the two families are one on the log scale
stats::sd(exp(-distrib_rng(d, 1000, theta))) > 0
#> [1] TRUE
```
