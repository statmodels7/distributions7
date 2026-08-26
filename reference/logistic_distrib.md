# Logistic Distribution

Builds the distribution object for the logistic family with mean \\\mu\\
and scale \\\sigma \> 0\\. The returned object carries closed-form
derivatives of the log-density to fourth order in the parameters and in
the response, and closed-form moments.

The family is symmetric about \\\mu\\ with variance \\\pi^2\sigma^2/3\\,
so \\\sigma\\ is a scale, not a standard deviation. Its distribution
function is the logistic sigmoid, the same curve `linkfunctions7` uses
as the inverse logit link.

## Usage

``` r
logistic_distrib(link_mu = identity_link(), link_sigma = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean \\\mu\\. Defaults
  to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html),
  the mean ranging over the whole line already.

- link_sigma:

  A `link` object from `linkfunctions7` for the scale \\\sigma\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

## Value

An S7 object of class `LogisticDistrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"logistic"`, `dimension`
`"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "sigma")`,
`n_params` `2`, `params_bounds` the list of \\(-\infty, \infty)\\ and
\\(0, \infty)\\, and `link_params` the two links given here.

## The parametrization

The density on \\y \in (-\infty, \infty)\\ is \$\$f(y; \mu, \sigma) =
\dfrac{\exp\left(-\dfrac{y-\mu}{\sigma}\right)}{\sigma \left\[1 +
\exp\left(-\dfrac{y-\mu}{\sigma}\right)\right\]^2},\$\$ with \\\mu \in
(-\infty, \infty)\\ and \\\sigma \in (0, \infty)\\. The distribution
function is \\F(q) = \[1 + e^{-(q-\mu)/\sigma}\]^{-1}\\ and the quantile
function \\Q(p) = \mu + \sigma \log(p/(1-p))\\.

The mean and the median are \\\mu\\, the variance is \\\pi^2
\sigma^2/3\\, the skewness is 0 and the excess kurtosis is \\6/5\\. A
logistic and a Gaussian matched on their variance are close in the body
and differ in the tails, the logistic's decaying exponentially and the
Gaussian's as a square exponential.

## Derivatives

With \\z = (y-\mu)/\sigma\\ the score is \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{1}{\sigma}\tanh\left(\dfrac{z}{2}\right),
\qquad \dfrac{\partial \ell}{\partial \sigma} =
-\dfrac{1}{\sigma}\left\[1 -
z\tanh\left(\dfrac{z}{2}\right)\right\],\$\$ and the expected Hessian is
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{3\sigma^2}, \quad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \sigma^2}\right\] = -\dfrac{3+\pi^2}{9\sigma^2}, \quad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu\\\partial
\sigma}\right\] = 0.\$\$

Third and fourth orders are closed form **observed** and numerical
**expected**: two of the nine expectations at those orders require
\\\int w^k \mathrm{sech}^4 w \tanh^2 w \\ dw\\, which has no elementary
form, so
[`distrib_deriv3.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.LogisticDistrib.md)
and
[`distrib_deriv4.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.LogisticDistrib.md)
route their `expected = TRUE` branch to
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md).
This is the only place in this family where a numerical route is taken.

## Estimation

There is no closed-form estimate;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale. The log-likelihood is
concave in \\\mu\\ for a fixed \\\sigma\\, so the location is well
determined.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean and
\\\sigma \> 0\\ the scale, with standard deviation \\\pi\sigma/\sqrt
3\\. \\z = (y-\mu)/\sigma\\ is the standardized residual. \\\eta\\ is a
parameter on the unconstrained scale of its link.

## References

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995). *Continuous
Univariate Distributions*, Volume 2, 2nd edition, Chapter 23. Wiley, New
York.

## See also

[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the light-tailed comparison;
[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md),
the asymmetric extreme-value relative;
[`linkfunctions7::logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html),
whose inverse is this distribution function;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[LogisticDistrib](https://statmodels7.github.io/distributions7/reference/LogisticDistrib.md)
for the class.

## Examples

``` r
d <- logistic_distrib()
d
#> Distribution: Logistic
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu    (mean)               | Link: identity   | Domain: (-Inf, Inf)
#>   sigma (scale)              | Link: log        | Domain: (0, Inf)

# The density and the distribution function are R's own.
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
all.equal(distrib_pdf(d, y, th), dlogis(y, 0.4, 1.5))
#> [1] TRUE

# sigma is a scale: the standard deviation is pi sigma / sqrt(3).
c(sd_from_variance = sqrt(variance(d, th)), pi * 1.5 / sqrt(3))
#> sd_from_variance                  
#>         2.720699         2.720699 

# Fitting recovers the parameters.
set.seed(9)
z <- distrib_rng(d, 3000, list(mu = 2, sigma = 1))
coef(fit_distrib(d, z))
#>        mu     sigma 
#> 2.0128244 0.9978246 

# Matched on variance, a logistic sits close to a Gaussian in the body and
# puts more mass in the tails.
g <- gaussian1_distrib()
s <- pi * 1.5 / sqrt(3)
rbind(logistic = distrib_cdf(d, 0.4 + c(1, 2, 3, 4) * s, th,
                             lower.tail = FALSE),
      gaussian = distrib_cdf(g, 0.4 + c(1, 2, 3, 4) * s,
                             list(mu = 0.4, sigma = s), lower.tail = FALSE))
#>               [,1]       [,2]        [,3]         [,4]
#> logistic 0.1401796 0.02589173 0.004314723 7.059941e-04
#> gaussian 0.1586553 0.02275013 0.001349898 3.167124e-05
```
