# Weibull Distribution, Scale and Shape

Builds the distribution object for the Weibull family parametrized by a
scale \\\mu \> 0\\ and a shape \\\sigma \> 0\\, with density \\f(y) =
(\sigma/\mu)(y/\mu)^{\sigma-1}\exp\\-(y/\mu)^{\sigma}\\\\ on \\y \> 0\\.
The returned object carries closed-form derivatives of the log-density
to fourth order, observed and expected, in the parameters and in the
response, and closed-form moments, so every generic of the toolkit
answers without a numerical fallback.

The two arguments choose the links that carry each parameter to the
unconstrained scale an optimizer works on. Both default to the
logarithm, both parameters being positive.

## Usage

``` r
weibull1_distrib(link_mu = log_link(), link_sigma = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the scale \\\mu\\. Defaults
  to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

- link_sigma:

  A `link` object from `linkfunctions7` for the shape \\\sigma\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  for the same reason.

## Value

An S7 object of class `Weibull1Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"weibull1"`, `dimension`
`"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "sigma")`,
`n_params` `2`, `params_bounds` the domain \\(0, \infty)\\ for both, and
`link_params` the two links given here.

## The parametrization

The density on \\y \in (0, \infty)\\ is \$\$f(y; \mu, \sigma) =
\dfrac{\sigma}{\mu} \left(\dfrac{y}{\mu}\right)^{\sigma - 1}
\exp\left\\-\left(\dfrac{y}{\mu}\right)^{\sigma}\right\\,\$\$ the
distribution function \\F(q) = 1 - \exp\\-(q/\mu)^{\sigma}\\\\ and the
quantile function \\Q(p) = \mu\\-\log(1-p)\\^{1/\sigma}\\. This is `WEI`
in gamlss.

**\\\mu\\ is the scale and not the mean.** The mean is \\\mu\\\Gamma(1 +
1/\sigma)\\, which involves the shape, so a mean parametrization would
make every derivative a derivative of the gamma function and of its
inverse. The scale-shape form keeps the whole family elementary;
[`mean.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Weibull1Distrib.md)
reports the mean, and
[`weibull3_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull3_distrib.md)
is the parametrization by a quantile for a reader who wants a location
that is one.

## Derivatives

Write \\z = y/\mu\\ and \\u = z^{\sigma}\\. The score is
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{\sigma}{\mu}(u - 1),
\qquad \dfrac{\partial \ell}{\partial \sigma} = \dfrac{1}{\sigma} + (1 -
u)\log z,\$\$ the observed Hessian \$\$\dfrac{\partial^2 \ell}{\partial
\mu^2} = \dfrac{\sigma}{\mu^2}\left\\1 - (1 + \sigma) u\right\\, \quad
\dfrac{\partial^2 \ell}{\partial \sigma^2} = -\dfrac{1}{\sigma^2} - u
(\log z)^2, \quad \dfrac{\partial^2 \ell}{\partial \mu \\ \partial
\sigma} = \dfrac{u - 1 + \sigma u \log z}{\mu},\$\$ and its expectation
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{\sigma^2}{\mu^2}, \quad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \sigma^2}\right\] = -\dfrac{(1-\gamma)^2 +
\pi^2/6}{\sigma^2}, \quad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu\\\partial \sigma}\right\] = \dfrac{1 -
\gamma}{\mu},\$\$ with \\\gamma\\ the Euler-Mascheroni constant. The
mixed entry is not zero, so the scale and the shape are asymptotically
correlated, unlike the mean and the standard deviation of a Gaussian.

Third and fourth orders are closed form as well, observed and expected,
in
[`distrib_deriv3.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Weibull1Distrib.md)
and
[`distrib_deriv4.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Weibull1Distrib.md),
as are the derivatives in the response,
[`distrib_grad_y.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Weibull1Distrib.md)
and
[`distrib_hess_y.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Weibull1Distrib.md),
and the mixed derivative
[`distrib_cross_y.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Weibull1Distrib.md).

## Why every expectation is elementary

Under the model \\u = (Y/\mu)^{\sigma}\\ is standard exponential
whatever the parameters are. Every expectation any order needs is
therefore a moment of \\u\\ against a power of \\\log u\\, and each of
those is a derivative of \\\Gamma\\ at 2: \\\mathbb{E}\[u\] = 1\\,
\\\mathbb{E}\[u\log u\] = 1 - \gamma\\, \\\mathbb{E}\[u(\log u)^2\] =
(1-\gamma)^2 + \pi^2/6 - 1\\. The Gumbel family shares the substitution,
\\\exp(-\text{Gumbel})\\ being Weibull, so
[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
rests on the same three moments.

## Moments and shape

With \\g_k = \Gamma(1 + k/\sigma)\\, the mean is \\\mu g_1\\ and the
variance \\\mu^2(g_2 - g_1^2)\\; the skewness and the excess kurtosis
follow from \\g_3\\ and \\g_4\\ and do not depend on \\\mu\\, the scale
being a multiplier. The hazard \\f/(1-F)\\ is
\\(\sigma/\mu)(y/\mu)^{\sigma-1}\\, increasing for \\\sigma \> 1\\ and
decreasing for \\\sigma \< 1\\. That monotone hazard is why the family
is used in survival work.

Two shapes are named families: \\\sigma = 1\\ is the exponential with
mean \\\mu\\, and \\\sigma = 2\\ the Rayleigh. Both are exact, not
limiting, and the example checks the first.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale. Neither estimate is
closed form: the shape solves \\1/\hat\sigma + \overline{\log y} = \sum
y_i^{\hat\sigma}\log y_i / \sum y_i^{\hat\sigma}\\ by iteration, and the
scale follows from it as \\\hat\mu = (n^{-1}\sum
y_i^{\hat\sigma})^{1/\hat\sigma}\\. The example checks that the fit
satisfies the second of these.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu \> 0\\ the scale,
\\\sigma \> 0\\ the shape and \\\gamma\\ the Euler-Mascheroni constant.
\\\eta\\ is a parameter on the unconstrained scale of its link, with
\\\theta = g^{-1}(\eta)\\.

## References

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994). *Continuous
Univariate Distributions*, Volume 1, 2nd edition, Chapter 21. Wiley, New
York.

Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive
models for location, scale and shape. *Journal of the Royal Statistical
Society, Series C*, **54**(3), 507-554.

## See also

[`weibull3_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull3_distrib.md)
for the same law parametrized by a quantile;
[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md),
which is this family after a negative logarithm;
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
and
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
for the special case and the generalization;
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
for another extreme-value family;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[Weibull1Distrib](https://statmodels7.github.io/distributions7/reference/Weibull1Distrib.md)
for the class.

## Examples

``` r
d <- weibull1_distrib()
d
#> Distribution: Weibull1
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu    (scale)              | Link: log        | Domain: (0, Inf)
#>   sigma (shape)              | Link: log        | Domain: (0, Inf)

# The density is stats::dweibull with the arguments in this order.
y <- c(0.5, 1.2, 3.0)
th <- list(mu = 2, sigma = 1.5)
all.equal(distrib_pdf(d, y, th), dweibull(y, shape = 1.5, scale = 2))
#> [1] TRUE

# The scale is not the mean: the mean carries a gamma function of the shape.
c(scale = th$mu, mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#>    scale     mean      var     skew     kurt 
#> 2.000000 1.805491 1.502761 1.071987 1.390404 

# Shape one is exactly the exponential with mean mu.
all.equal(distrib_pdf(d, y, list(mu = 2, sigma = 1)), dexp(y, rate = 1 / 2))
#> [1] TRUE

# The hazard increases with y above shape one and decreases below it.
haz <- function(s) {
  th <- list(mu = 2, sigma = s)
  distrib_pdf(d, c(0.5, 1, 2, 4), th) /
    distrib_cdf(d, c(0.5, 1, 2, 4), th, lower.tail = FALSE)
}
rbind(shape_0.5 = haz(0.5), shape_2.5 = haz(2.5))
#>              [,1]      [,2] [,3]      [,4]
#> shape_0.5 0.50000 0.3535534 0.25 0.1767767
#> shape_2.5 0.15625 0.4419417 1.25 3.5355339

# Fitting recovers the parameters, and the scale satisfies the profile
# identity given the fitted shape.
set.seed(9)
z <- distrib_rng(d, 2000, list(mu = 2, sigma = 1.5))
fit <- fit_distrib(d, z)
cf <- coef(fit)
c(cf, profile_mu = mean(z^cf[["sigma"]])^(1 / cf[["sigma"]]))
#>         mu      sigma profile_mu 
#>   2.024867   1.512341   2.024867 
```
