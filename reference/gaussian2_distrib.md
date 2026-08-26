# Gaussian Distribution, Mean and Variance

Builds the distribution object for the Gaussian (normal) family
parametrized by its mean \\\mu\\ and its variance \\\sigma^2 \> 0\\. The
returned object carries closed-form derivatives of the log-density to
fourth order, in the parameters and in the response, and closed-form
moments, so every generic of the toolkit answers without a numerical
fallback.

The two arguments choose the links that carry each parameter to the
unconstrained scale an optimizer works on. The defaults are the identity
for the mean, which is already free, and the logarithm for the variance,
which keeps it positive at every predictor.

## Usage

``` r
gaussian2_distrib(link_mu = identity_link(), link_sigma2 = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean \\\mu\\. Defaults
  to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html),
  the mean ranging over the whole line already.

- link_sigma2:

  A `link` object from `linkfunctions7` for the variance \\\sigma^2\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

## Value

An S7 object of class `Gaussian2Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"gaussian2"`, `dimension`
`"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "sigma2")`,
`n_params` `2`, `params_bounds` the list of \\(-\infty, \infty)\\ and
\\(0, \infty)\\, and `link_params` the two links given here.

## The parametrization

The density on \\y \in (-\infty, \infty)\\ is \$\$f(y; \mu, \sigma^2) =
\dfrac{1}{\sqrt{2\pi\sigma^{2}}}
\exp\left\\-\dfrac{(y-\mu)^{2}}{2\sigma^{2}}\right\\,\$\$ with \\\mu \in
(-\infty, \infty)\\ and \\\sigma^2 \in (0, \infty)\\. The mean is
\\\mu\\, the variance \\\sigma^2\\, and both the skewness and the excess
kurtosis are 0. The numbering follows the literature where it has one:
this parametrization is `NO2` in the gamlss family catalog.

This is the same law as
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
in different coordinates, \\\sigma^2\\ here being the square of the
\\\sigma\\ there. The two are separate families and not one family under
a link, because a link changes the scale a parameter is *modeled* on and
leaves the parameter what it was, while here the parameter, its
interpretation, its standard error and its confidence interval are all
about the variance.

## Derivatives

Writing \\r = y - \mu\\ and \\v = \sigma^2\\, the score is
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{r}{v}, \qquad
\dfrac{\partial \ell}{\partial v} = \dfrac{r^2 - v}{2v^2},\$\$ the
observed Hessian \$\$\ell^{(\mu\mu)} = -\dfrac{1}{v}, \quad \ell^{(\mu
v)} = -\dfrac{r}{v^2}, \quad \ell^{(vv)} = \dfrac{1}{2v^2} -
\dfrac{r^2}{v^3},\$\$ and its expectation
\$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] = -\dfrac{1}{v}, \quad
\mathbb{E}\left\[\ell^{(\mu v)}\right\] = 0, \quad
\mathbb{E}\left\[\ell^{(vv)}\right\] = -\dfrac{1}{2v^2}.\$\$ The zero
off-diagonal makes the mean and the variance orthogonal, so their
maximum likelihood estimates are asymptotically independent.
Orthogonality holds in all three parametrizations of this family, the
mean being orthogonal to any smooth function of the spread alone.

Third and fourth orders are closed form as well, observed and expected,
in
[`distrib_deriv3.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gaussian2Distrib.md)
and
[`distrib_deriv4.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gaussian2Distrib.md),
as are the derivatives in the response,
[`distrib_grad_y.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Gaussian2Distrib.md)
and
[`distrib_hess_y.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Gaussian2Distrib.md).

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale and reaches the
closed-form estimates \\\hat\mu = \bar y\\ and \\\hat\sigma^2 =
n^{-1}\sum (y_i - \bar y)^2\\, the maximum likelihood estimate with
divisor \\n\\. The example below checks both against the sample.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean and
\\\sigma^2 \> 0\\ the variance. \\\ell^{(ij)}\\ is a second derivative
of \\\ell\\ in parameters \\i\\ and \\j\\. \\\eta\\ is a parameter on
the unconstrained scale of its link, with \\\theta = g^{-1}(\eta)\\.

## References

Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive
models for location, scale and shape. *Journal of the Royal Statistical
Society: Series C* **54**, 507-554.

## See also

[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
and
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md)
for the same law in the standard deviation and in the precision;
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
for the Gaussian on a log scale;
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
and
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
for heavier tails;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[Gaussian2Distrib](https://statmodels7.github.io/distributions7/reference/Gaussian2Distrib.md)
for the class.

## Examples

``` r
d <- gaussian2_distrib()
d
#> Distribution: Gaussian2
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu     (mean)               | Link: identity   | Domain: (-Inf, Inf)
#>   sigma2 (variance)           | Link: log        | Domain: (0, Inf)

# The same law as gaussian1 at sigma = sqrt(sigma2).
y <- c(-1.2, 0.3, 2.5)
all.equal(distrib_pdf(d, y, list(mu = 1, sigma2 = 4)),
          distrib_pdf(gaussian1_distrib(), y, list(mu = 1, sigma = 2)))
#> [1] TRUE

# Moments in closed form: the variance is the parameter itself.
th <- list(mu = 1, sigma2 = 4)
c(mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#> mean  var skew kurt 
#>    1    4    0    0 

# Fitting recovers the closed-form maximum likelihood estimates, with the
# variance divided by n.
set.seed(7)
z <- distrib_rng(d, 400, list(mu = 3, sigma2 = 4))
fit <- fit_distrib(d, z)
rbind(fitted = coef(fit),
      closed = c(mu = mean(z), sigma2 = mean((z - mean(z))^2)))
#>              mu   sigma2
#> fitted 3.091666 4.047184
#> closed 3.091666 4.047184

# The estimate is a variance, so the interval it reports is an interval for
# the variance and stays positive.
confint(fit)
#>            2.5%    97.5%
#> mu     2.894517 3.288815
#> sigma2 3.523416 4.648812
```
