# Gaussian Distribution, Mean and Standard Deviation

Builds the distribution object for the Gaussian (normal) family
parametrized by its mean \\\mu\\ and its standard deviation \\\sigma \>
0\\. The returned object carries closed-form derivatives of the
log-density to fourth order, in the parameters and in the response, and
closed-form moments, so every generic of the toolkit answers without a
numerical fallback.

The two arguments choose the links that carry each parameter to the
unconstrained scale an optimizer works on. The defaults are the identity
for the mean, which is already free, and the logarithm for the standard
deviation, which keeps it positive at every predictor.

## Usage

``` r
gaussian1_distrib(link_mu = identity_link(), link_sigma = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean \\\mu\\. Defaults
  to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html),
  the mean ranging over the whole line already.

- link_sigma:

  A `link` object from `linkfunctions7` for the standard deviation
  \\\sigma\\. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

## Value

An S7 object of class `Gaussian1Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"gaussian1"`, `dimension`
`"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "sigma")`,
`n_params` `2`, `params_bounds` the list of \\(-\infty, \infty)\\ and
\\(0, \infty)\\, and `link_params` the two links given here.

## The parametrization

The density on \\y \in (-\infty, \infty)\\ is \$\$f(y; \mu, \sigma) =
\dfrac{1}{\sqrt{2\pi}\\\sigma}
\exp\left\\-\dfrac{1}{2}\left(\dfrac{y-\mu}{\sigma}\right)^2\right\\,\$\$
with \\\mu \in (-\infty, \infty)\\ and \\\sigma \in (0, \infty)\\. The
distribution function is \\F(q) = \Phi((q-\mu)/\sigma)\\ and the
quantile function \\Q(p) = \mu + \sigma \Phi^{-1}(p)\\, with \\\Phi\\
the standard normal distribution function.

The mean is \\\mu\\, the variance \\\sigma^2\\, and both the skewness
and the excess kurtosis are 0. Two sibling parametrizations of the same
law are
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md),
in the mean and the variance, and
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md),
in the mean and the precision.

## Derivatives

The score is \$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y -
\mu}{\sigma^2}, \qquad \dfrac{\partial \ell}{\partial \sigma} =
\dfrac{(y - \mu)^2 - \sigma^2}{\sigma^3},\$\$ the observed Hessian
\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{\sigma^2},
\quad \dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^2 -
3(y-\mu)^2}{\sigma^4}, \quad \dfrac{\partial^2 \ell}{\partial
\mu\\\partial \sigma} = -\dfrac{2(y-\mu)}{\sigma^3},\$\$ and its
expectation \$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial
\mu^2}\right\] = -\dfrac{1}{\sigma^2}, \quad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right\] =
-\dfrac{2}{\sigma^2}, \quad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu\\\partial \sigma}\right\] = 0.\$\$ The zero
off-diagonal makes the mean and the standard deviation orthogonal, so
their maximum likelihood estimates are asymptotically independent.

Third and fourth orders are closed form as well, observed and expected,
in
[`distrib_deriv3.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gaussian1Distrib.md)
and
[`distrib_deriv4.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gaussian1Distrib.md),
as are the derivatives in the response,
[`distrib_grad_y.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Gaussian1Distrib.md)
and
[`distrib_hess_y.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Gaussian1Distrib.md).

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale and reaches the
closed-form estimates \\\hat\mu = \bar y\\ and \\\hat\sigma^2 =
n^{-1}\sum (y_i - \bar y)^2\\, the maximum likelihood estimate of the
variance with divisor \\n\\. The example below checks both against the
sample.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean and
\\\sigma \> 0\\ the standard deviation. \\\eta\\ is a parameter on the
unconstrained scale of its link, with \\\theta = g^{-1}(\eta)\\.

## References

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994). *Continuous
Univariate Distributions*, Volume 1, 2nd edition, Chapter 13. Wiley, New
York.

## See also

[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
and
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md)
for the same law in the variance and in the precision;
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
[Gaussian1Distrib](https://statmodels7.github.io/distributions7/reference/Gaussian1Distrib.md)
for the class.

## Examples

``` r
d <- gaussian1_distrib()
d
#> Distribution: Gaussian1
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu    (mean)               | Link: identity   | Domain: (-Inf, Inf)
#>   sigma (standard deviation) | Link: log        | Domain: (0, Inf)

# The density and the distribution function are R's own.
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
all.equal(distrib_pdf(d, y, th), dnorm(y, 0.4, 1.5))
#> [1] TRUE

# Moments in closed form: skewness and excess kurtosis are 0.
c(mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#> mean  var skew kurt 
#> 0.40 2.25 0.00 0.00 

# Fitting recovers the closed-form maximum likelihood estimates, with the
# variance divided by n.
set.seed(7)
z <- distrib_rng(d, 400, list(mu = 3, sigma = 2))
fit <- fit_distrib(d, z)
rbind(fitted = coef(fit),
      closed = c(mu = mean(z), sigma = sqrt(mean((z - mean(z))^2))))
#>              mu    sigma
#> fitted 3.091666 2.011761
#> closed 3.091666 2.011761

# A different link changes the scale the optimizer moves on, not the answer.
d2 <- gaussian1_distrib(link_sigma = linkfunctions7::sqrt_link())
all.equal(coef(fit_distrib(d2, z)), coef(fit), tolerance = 1e-6)
#> [1] TRUE
```
