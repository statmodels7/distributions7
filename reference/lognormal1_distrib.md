# Lognormal Distribution, Log-Scale Parametrization

Builds the distribution object for the lognormal family on \\(0,
\infty)\\: the law of a variable whose logarithm is Gaussian with mean
\\\mu\\ and variance \\\sigma^2 \> 0\\. Both parameters are moments **of
the logarithm**, so neither is the mean or the variance of \\Y\\. The
returned object carries closed-form derivatives of the log-density to
fourth order, in the parameters and in the response, and closed-form
moments, so every generic of the toolkit answers without a numerical
fallback.

The two arguments choose the links that carry each parameter to the
unconstrained scale an optimizer works on. The defaults are the identity
for \\\mu\\, which ranges over the whole line, and the logarithm for
\\\sigma^2\\, which keeps it positive at every predictor.

## Usage

``` r
lognormal1_distrib(link_mu = identity_link(), link_sigma2 = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the location \\\mu\\.
  Defaults to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html),
  \\\mu\\ being a mean on the log scale and so already free.

- link_sigma2:

  A `link` object from `linkfunctions7` for the variance \\\sigma^2\\ of
  the logarithm. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line.

## Value

An S7 object of class `Lognormal1Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"lognormal1"`, `dimension`
`"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "sigma2")`,
`n_params` `2`, `params_bounds` the list of \\(-\infty, \infty)\\ and
\\(0, \infty)\\, and `link_params` the two links given here.

## The parametrization

The density on \\y \in (0, \infty)\\ is \$\$f(y; \mu, \sigma^2) =
\dfrac{1}{y\sqrt{2\pi\sigma^2}} \exp\left\\-\dfrac{(\log y -
\mu)^2}{2\sigma^2}\right\\,\$\$ the distribution function \\F(q) =
\Phi\\(\log q - \mu)/\sigma\\\\ and the quantile function \\Q(p) =
\exp\\\mu + \sigma\Phi^{-1}(p)\\\\, both closed form because the log
transformation is monotone.

On the original scale the mean is \\e^{\mu + \sigma^2/2}\\, the variance
\\(e^{\sigma^2}-1)e^{2\mu+\sigma^2}\\, the skewness
\\(e^{\sigma^2}+2)\sqrt{e^{\sigma^2}-1}\\ and the excess kurtosis
\\e^{4\sigma^2}+2e^{3\sigma^2}+3e^{2\sigma^2}-6\\. Three quantities are
worth keeping apart: the mode is \\e^{\mu-\sigma^2}\\, the median
\\e^{\mu}\\ and the mean \\e^{\mu+\sigma^2/2}\\, in that order.

[`lognormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal2_distrib.md)
carries the same law parametrized by the mean and the variance of \\Y\\
itself.

## Every derivative in the parameters is the Gaussian's

The log-density is the Gaussian's evaluated at \\\log y\\, minus \\\log
y\\. That last term is the Jacobian of the transformation and carries no
parameter, so it disappears from every derivative in \\\mu\\ and
\\\sigma^2\\. Writing \\r = \log y - \mu\\, the score is
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{r}{\sigma^2}, \qquad
\dfrac{\partial \ell}{\partial \sigma^2} = \dfrac{r^2 -
\sigma^2}{2\sigma^4},\$\$ and the expected Hessian is
\$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] = -\dfrac{1}{\sigma^2},
\quad \mathbb{E}\left\[\ell^{(\sigma^2\sigma^2)}\right\] =
-\dfrac{1}{2\sigma^4}, \quad
\mathbb{E}\left\[\ell^{(\mu\sigma^2)}\right\] = 0.\$\$ Every one of
these agrees with
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
at \\\log y\\, component for component and at all four orders, observed
and expected. The two parameters are orthogonal and the information does
not move with the mean of \\Y\\.

The derivatives **in the response** are where the two families differ,
the Jacobian being a function of \\y\\. They are still closed form; see
[`distrib_grad_y.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Lognormal1Distrib.md)
and
[`distrib_hess_y.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Lognormal1Distrib.md),
and note that the log-density is convex in the response above
\\e^{\mu+1-\sigma^2}\\.

## Estimation

Both maximum likelihood estimates are closed form and are the sample
moments of the logarithm: \$\$\hat\mu = \dfrac{1}{n}\sum_i \log y_i,
\qquad \hat\sigma^2 = \dfrac{1}{n}\sum_i (\log y_i - \hat\mu)^2,\$\$ the
second with divisor \\n\\.
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
reaches them on the link scale, and the example below checks both
against the sample.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean and
\\\sigma^2 \> 0\\ the variance **of \\\log Y\\**. \\\Phi\\ is the
standard normal distribution function and \\r = \log y - \mu\\ the
residual on the log scale. \\\eta\\ is a parameter on the unconstrained
scale of its link, with \\\theta = g^{-1}(\eta)\\.

## References

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994). *Continuous
Univariate Distributions*, Volume 1, 2nd edition, Chapter 14. Wiley, New
York.

## See also

[`lognormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal2_distrib.md)
for the same law in the mean and the variance of \\Y\\;
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
for the law of its logarithm;
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
and
[`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md)
for other positive families with a multiplicative variance function;
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
for the general wrapper this family is a named instance of;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[Lognormal1Distrib](https://statmodels7.github.io/distributions7/reference/Lognormal1Distrib.md)
for the class.

## Examples

``` r
d <- lognormal1_distrib()
d
#> Distribution: Lognormal1
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu     (mean (log scale))   | Link: identity   | Domain: (-Inf, Inf)
#>   sigma2 (variance (log scale)) | Link: log        | Domain: (0, Inf)

# The density is stats::dlnorm at meanlog = mu, sdlog = sqrt(sigma2).
y <- c(0.5, 1.6, 4)
th <- list(mu = 0.5, sigma2 = 0.36)
all.equal(distrib_pdf(d, y, th), dlnorm(y, 0.5, sqrt(0.36)))
#> [1] TRUE

# Mode, median and mean, in that order.
c(mode = exp(0.5 - 0.36), median = exp(0.5), mean = mean(d, th))
#>     mode   median     mean 
#> 1.150274 1.648721 1.973878 

# Moments on the original scale, all exponential in sigma2.
c(mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#>      mean       var      skew      kurt 
#>  1.973878  1.688335  2.260084 10.273355 

# Every derivative in the parameters is the Gaussian's at log y.
all.equal(distrib_gradient(d, y, th),
          distrib_gradient(gaussian2_distrib(), log(y), th))
#> [1] TRUE

# Fitting recovers the closed-form estimates, the sample moments of the
# logarithm with divisor n.
set.seed(6)
z <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, z)
rbind(fitted = coef(fit),
      closed = c(mu = mean(log(z)),
                 sigma2 = mean((log(z) - mean(log(z)))^2)))
#>               mu    sigma2
#> fitted 0.4928186 0.3600616
#> closed 0.4928186 0.3600616
```
