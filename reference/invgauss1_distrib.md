# Inverse Gaussian Distribution, Mean and Dispersion

Builds the distribution object for the inverse Gaussian family on \\(0,
\infty)\\ parametrized by its mean \\\mu \> 0\\ and a dispersion \\\phi
\> 0\\, so that \\\operatorname{Var}(Y) = \phi\mu^3\\. This is the
generalized linear model parametrization, with variance function
\\V(\mu) = \mu^3\\. The returned object carries closed-form derivatives
of the log-density to fourth order, in the parameters and in the
response, and closed-form moments, so every generic of the toolkit
answers without a numerical fallback.

The two arguments choose the links that carry each parameter to the
unconstrained scale an optimizer works on. Both default to the
logarithm, both parameters being positive.

## Usage

``` r
invgauss1_distrib(link_mu = log_link(), link_phi = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean \\\mu\\. Defaults
  to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  mean positive.

- link_phi:

  A `link` object from `linkfunctions7` for the dispersion \\\phi\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  for the same reason.

## Value

An S7 object of class `InvGauss1Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"invgauss1"`, `dimension`
`"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "phi")`,
`n_params` `2`, `params_bounds` the domain \\(0, \infty)\\ for both, and
`link_params` the two links given here.

## The parametrization

The density on \\y \in (0, \infty)\\ is \$\$f(y; \mu, \phi) =
\sqrt{\dfrac{1}{2\pi\phi y^3}} \exp\left\\-\dfrac{(y-\mu)^2}{2\phi\mu^2
y}\right\\,\$\$ the mean is \\\mu\\, the variance \\\phi\mu^3\\, the
skewness \\3\sqrt{\phi\mu}\\ and the excess kurtosis \\15\phi\mu\\. The
distribution function is elementary in the standard normal one, \$\$F(q)
= \Phi\left\\\sqrt{\dfrac{1}{\phi q}} \left(\dfrac{q}{\mu} -
1\right)\right\\ + e^{2/(\phi\mu)}\\\Phi\left\\-\sqrt{\dfrac{1}{\phi q}}
\left(\dfrac{q}{\mu} + 1\right)\right\\,\$\$ and the quantile function
is its numerical inverse.

The law is the first-passage time of a Brownian motion with drift, which
is why it is skewed however small the dispersion and why its variance
grows with the cube of the mean. Beside the gamma, whose variance
function is \\\mu^2\\, it puts more mass both near zero and far out.

This is the same law as
[`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md),
which carries the mean and the *variance*, the two being related by
\\\sigma^2 = \phi\mu^3\\. They are separate families because the second
parameter is a different quantity in each, with its own interpretation,
standard error and interval.

## Derivatives

The score is \$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y -
\mu}{\phi\mu^3}, \qquad \dfrac{\partial \ell}{\partial \phi} =
\dfrac{(y - \mu)^2 - y\mu^2\phi}{2y\phi^2\mu^2},\$\$ and the expected
Hessian is \$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] =
-\dfrac{1}{\phi\mu^3}, \qquad \mathbb{E}\left\[\ell^{(\phi\phi)}\right\]
= -\dfrac{1}{2\phi^2}, \qquad \mathbb{E}\left\[\ell^{(\mu\phi)}\right\]
= 0.\$\$ The zero off-diagonal makes the mean and the dispersion
orthogonal, so the mean equation can be fitted with the dispersion held
at any value without biasing it. Note that the expected curvature in
\\\phi\\ is \\-1/(2\phi^2)\\ whatever the mean, so the dispersion is
estimated with the same precision at every scale of the response.

The observed Hessian is a different matter, and the method page says so
at length: neither diagonal entry is negative at every observation. That
is the practical reason to fit this family by
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md).

Third and fourth orders are closed form as well, observed and expected,
in
[`distrib_deriv3.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.InvGauss1Distrib.md)
and
[`distrib_deriv4.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.InvGauss1Distrib.md),
as are the derivatives in the response,
[`distrib_grad_y.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.InvGauss1Distrib.md)
and
[`distrib_hess_y.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.InvGauss1Distrib.md).
The derivatives of the *distribution* function in the parameters are
closed form too, the expression above being elementary.

## Estimation

Both maximum likelihood estimates are available in closed form, which
few two-parameter families offer: \$\$\hat\mu = \bar y, \qquad \hat\phi
= \dfrac{1}{n}\sum_i \dfrac{1}{y_i} - \dfrac{1}{\bar y}.\$\$
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
reaches them numerically on the link scale, and the example below checks
both against the sample.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu \> 0\\ the mean
and \\\phi \> 0\\ the dispersion, with \\\operatorname{Var}(Y) =
\phi\mu^3\\. \\\Phi\\ is the standard normal distribution function.
\\\eta\\ is a parameter on the unconstrained scale of its link, with
\\\theta = g^{-1}(\eta)\\.

## References

Chhikara, R. S. and Folks, J. L. (1989). *The Inverse Gaussian
Distribution: Theory, Methodology, and Applications*. Marcel Dekker, New
York.

## See also

[`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md)
for the same law in the mean and the variance;
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
for the other positive generalized linear model family, with variance
function \\\mu^2\\;
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
for a third;
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
for the counts this law mixes a Poisson into;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[InvGauss1Distrib](https://statmodels7.github.io/distributions7/reference/InvGauss1Distrib.md)
for the class.

## Examples

``` r
d <- invgauss1_distrib()
d
#> Distribution: Invgauss1
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu  (mean)               | Link: log        | Domain: (0, Inf)
#>   phi (dispersion)         | Link: log        | Domain: (0, Inf)

# The density is statmod::dinvgauss at this parametrization.
y <- c(0.5, 1, 2)
th <- list(mu = 1, phi = 2)
all.equal(distrib_pdf(d, y, th),
          statmod::dinvgauss(y, mean = 1, dispersion = 2))
#> [1] TRUE

# Moments: variance phi mu^3, skewness 3 sqrt(phi mu), kurtosis 15 phi mu.
c(mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#>      mean       var      skew      kurt 
#>  1.000000  2.000000  4.242641 30.000000 
c(2 * 1^3, 3 * sqrt(2 * 1), 15 * 2 * 1)
#> [1]  2.000000  4.242641 30.000000

# Fitting recovers the closed-form maximum likelihood estimates.
set.seed(3)
z <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, z)
rbind(fitted = coef(fit),
      closed = c(mu = mean(z), phi = mean(1 / z) - 1 / mean(z)))
#>             mu      phi
#> fitted 1.00815 1.988966
#> closed 1.00815 1.988964

# The mean and the dispersion are orthogonal: the mixed entry is 0.
distrib_expected_hessian(d, 1, th)$mu_phi
#> [1] 0
```
