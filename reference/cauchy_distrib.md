# Cauchy Distribution

Builds the distribution object for the Cauchy family with location
\\\mu\\ and scale \\\sigma \> 0\\, the Student's t on one degree of
freedom. The returned object carries closed-form derivatives of the
log-density to fourth order, in the parameters and in the response, so
every generic of the toolkit answers without a numerical fallback.

No moment of this family exists, so the two parameters are the median
and the half-interquartile range. Its tails decay like \\y^{-2}\\, and
the resulting score is bounded, which makes a Cauchy likelihood a
standard choice when the data carry gross outliers.

## Usage

``` r
cauchy_distrib(link_mu = identity_link(), link_sigma = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the location \\\mu\\.
  Defaults to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html),
  the location ranging over the whole line already.

- link_sigma:

  A `link` object from `linkfunctions7` for the scale \\\sigma\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

## Value

An S7 object of class `CauchyDistrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"cauchy"`, `dimension`
`"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "sigma")`,
`n_params` `2`, `params_bounds` the list of \\(-\infty, \infty)\\ and
\\(0, \infty)\\, and `link_params` the two links given here.

## The parametrization

The density on \\y \in (-\infty, \infty)\\ is \$\$f(y; \mu, \sigma) =
\dfrac{1}{\pi \sigma \left\[1 +
\left(\dfrac{y-\mu}{\sigma}\right)^2\right\]},\$\$ with \\\mu \in
(-\infty, \infty)\\ and \\\sigma \in (0, \infty)\\. The distribution
function is \\F(q) = 1/2 + \pi^{-1}\arctan((q-\mu)/\sigma)\\ and the
quantile function \\Q(p) = \mu + \sigma\tan(\pi(p - 1/2))\\.

## No moments

The density decays like \\y^{-2}\\, so \\\int \|y\|^p f(y)\\dy\\
diverges for every \\p \ge 1\\ and the family has no mean, no variance
and no higher moment. [`mean()`](https://rdrr.io/r/base/mean.html),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
return `NaN` directly, without attempting a quadrature over a divergent
integral;
[`skewness.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.CauchyDistrib.md)
gives the argument in full. What the family does have is a median, equal
to \\\mu\\, and quartiles at \\\mu \pm \sigma\\, so \\\sigma\\ is the
half-interquartile range.

## Derivatives

With \\r = y - \mu\\ and \\d = \sigma^2 + r^2\\ the score is
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{2r}{d}, \qquad
\dfrac{\partial \ell}{\partial \sigma} = \dfrac{r^2 - \sigma^2}{\sigma
d},\$\$ the observed Hessian \$\$\dfrac{\partial^2 \ell}{\partial \mu^2}
= \dfrac{2(r^2 - \sigma^2)}{d^2}, \quad \dfrac{\partial^2 \ell}{\partial
\sigma^2} = \dfrac{\sigma^4 - 4\sigma^2 r^2 - r^4}{\sigma^2 d^2}, \quad
\dfrac{\partial^2 \ell}{\partial \mu\\\partial \sigma} = -\dfrac{4\sigma
r}{d^2},\$\$ and its expectation \$\$\mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu^2}\right\] = \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \sigma^2}\right\] = -\dfrac{1}{2\sigma^2}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu\\\partial
\sigma}\right\] = 0.\$\$ The information in the location is half a
Gaussian's at the same scale. Every derivative is bounded in \\y\\, so
an expectation of any of them exists even though no moment of \\Y\\
does.

Third and fourth orders are closed form as well, observed and expected,
in
[`distrib_deriv3.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.CauchyDistrib.md)
and
[`distrib_deriv4.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.CauchyDistrib.md),
as are the derivatives in the response.

## Estimation

There is no closed-form estimate: the likelihood equations are
polynomial in \\\mu\\ of degree \\2n - 1\\ and are solved numerically.
Because the score redescends, the log-likelihood can carry several local
maxima, and
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
starts from a data-based value and takes Fisher scoring steps, whose
expected information is positive definite where the observed one need
not be.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the location and
\\\sigma \> 0\\ the scale. Neither is a moment. \\\eta\\ is a parameter
on the unconstrained scale of its link, with \\\theta = g^{-1}(\eta)\\.

## References

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994). *Continuous
Univariate Distributions*, Volume 1, 2nd edition, Chapter 16. Wiley, New
York.

## See also

[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md),
of which this is the case \\\nu = 1\\, and which estimates the degrees
of freedom instead of fixing them;
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
and
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
for the other robust location families;
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the light-tailed comparison;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[CauchyDistrib](https://statmodels7.github.io/distributions7/reference/CauchyDistrib.md)
for the class.

## Examples

``` r
d <- cauchy_distrib()
d
#> Distribution: Cauchy
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu    (location)           | Link: identity   | Domain: (-Inf, Inf)
#>   sigma (scale)              | Link: log        | Domain: (0, Inf)

# The density and the distribution function are R's own.
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
all.equal(distrib_pdf(d, y, th), dcauchy(y, 0.4, 1.5))
#> [1] TRUE

# No moment exists; the median and the half-interquartile range do.
c(mean = mean(d, th), variance = variance(d, th))
#>     mean variance 
#>      NaN      NaN 
distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#> [1] -1.1  0.4  1.9

# Fitting a contaminated sample: the Cauchy location tracks the bulk where
# the sample mean is dragged away by the outliers.
set.seed(5)
z <- c(rnorm(200, mean = 3, sd = 1), rnorm(10, mean = 60, sd = 1))
fit <- fit_distrib(d, z)
c(cauchy_mu = unname(coef(fit)["mu"]),
  median = median(z), mean = mean(z))
#> cauchy_mu    median      mean 
#>  2.957763  2.973282  5.748700 
```
