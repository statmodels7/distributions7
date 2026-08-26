# Beta Distribution, Mean and Precision

Builds the distribution object for the beta family on \\(0, 1)\\
parametrized by its mean \\\mu \in (0, 1)\\ and a precision \\\phi \>
0\\, so that the shapes are \\\alpha = \mu\phi\\ and \\\beta =
(1-\mu)\phi\\. This is the parametrization beta regression uses. The
returned object carries closed-form derivatives of the log-density to
fourth order, in the parameters and in the response, and closed-form
moments, so every generic of the toolkit answers without a numerical
fallback.

The two arguments choose the links that carry each parameter to the
unconstrained scale an optimizer works on. The defaults are the logit
for the mean, which lives in \\(0, 1)\\, and the logarithm for the
precision, which is positive.

## Usage

``` r
beta1_distrib(link_mu = logit_link(), link_phi = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean \\\mu\\. Defaults
  to
  [`linkfunctions7::logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html),
  which maps \\(0, 1)\\ onto the line and so keeps every fitted mean
  inside the unit interval.

- link_phi:

  A `link` object from `linkfunctions7` for the precision \\\phi\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line.

## Value

An S7 object of class `Beta1Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"beta1"`, `dimension`
`"univariate"`, `bounds` `c(0, 1)`, `params` `c("mu", "phi")`,
`n_params` `2`, `params_bounds` the list of \\(0, 1)\\ and \\(0,
\infty)\\, and `link_params` the two links given here.

## The parametrization

The density on \\y \in (0, 1)\\ is \$\$f(y; \mu, \phi) =
\dfrac{\Gamma(\phi)}{\Gamma(\alpha)\Gamma(\beta)}\\
y^{\alpha-1}(1-y)^{\beta-1}, \qquad \alpha = \mu\phi, \quad \beta =
(1-\mu)\phi,\$\$ the distribution function \\F(q) = I_q(\alpha, \beta)\\
with \\I\\ the regularized incomplete beta function, and the quantile
function its numerical inverse. The mean is \\\mu\\, the variance
\\\mu(1-\mu)/(\phi+1)\\ and the skewness
\\2(1-2\mu)\sqrt{\phi+1}/\\(\phi+2)\sqrt{\mu(1-\mu)}\\\\.

The name *precision* is earned: at a fixed mean the variance falls as
\\1/(\phi+1)\\, so \\\phi\\ says how tightly the mass concentrates. Two
settings are worth recognizing. At \\\mu = 1/2\\ and \\\phi = 2\\ both
shapes are 1 and the density is the uniform. Where a shape falls below
one, which happens when \\\phi\\ is small, the density is unbounded at
the corresponding endpoint.

This is the same law as
[`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md),
which carries the mean and the *variance*, the two being related by
\\\sigma^2 = \mu(1-\mu)/(\phi+1)\\. They are separate families because
the second parameter is a different quantity in each, with its own
interpretation, standard error and interval.

## Derivatives

With \\\psi\\ the digamma function and \\\psi_1\\ the trigamma, the
score is \$\$\dfrac{\partial \ell}{\partial \mu} = \phi\left\\
\log\dfrac{y}{1-y} - \psi(\alpha) + \psi(\beta)\right\\, \qquad
\dfrac{\partial \ell}{\partial \phi} = \psi(\phi) - \mu\psi(\alpha) -
(1-\mu)\psi(\beta) + \mu\log y + (1-\mu)\log(1-y),\$\$ and the expected
Hessian is \$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] =
-\phi^2\left\\\psi_1(\alpha) + \psi_1(\beta)\right\\, \qquad
\mathbb{E}\left\[\ell^{(\phi\phi)}\right\] = \psi_1(\phi) -
\mu^2\psi_1(\alpha) - (1-\mu)^2\psi_1(\beta),\$\$
\$\$\mathbb{E}\left\[\ell^{(\mu\phi)}\right\] =
-\phi\left\\\mu\psi_1(\alpha) - (1-\mu)\psi_1(\beta)\right\\.\$\$ The
mixed entry is zero only at \\\mu = 1/2\\, so the mean and the precision
are in general not orthogonal.

## Where the response stops entering

The data reach the log-density only through \\(\alpha-1)\log y +
(\beta-1)\log(1-y)\\, which is linear in the shapes and so at most
quadratic in \\(\mu, \phi)\\ through the bilinear map \\\alpha =
\mu\phi\\. Two consequences follow, and both are visible in the methods:

- the observed Hessian equals its expectation in `mu_mu` and `phi_phi`
  and differs in `mu_phi` by the log-odds residual \\\log\\y/(1-y)\\ -
  \psi(\alpha) + \psi(\beta)\\, whose expectation is zero;

- **every** third and fourth derivative is free of the response, so
  [`distrib_deriv3.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Beta1Distrib.md)
  and
  [`distrib_deriv4.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Beta1Distrib.md)
  return the same numbers whether or not `expected = TRUE` is asked for.

The derivatives of the *distribution* function in the parameters have no
elementary form, the derivative of an incomplete beta in its shapes
being hypergeometric, and are taken by finite difference on the analytic
cdf.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale. Neither estimate is
closed form: both estimating equations involve digamma functions of the
shapes. The method of moments supplies the starting values \\\hat\mu =
\bar y\\ and \\\hat\phi = \bar y(1-\bar y)/s^2 - 1\\, and the example
below shows them landing beside the estimates.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu \in (0,1)\\ the
mean and \\\phi \> 0\\ the precision. \\\alpha\\ and \\\beta\\ are the
implied shapes. \\\psi\\ is the digamma function and \\\psi_m\\ its
\\m\\th derivative. \\\eta\\ is a parameter on the unconstrained scale
of its link, with \\\theta = g^{-1}(\eta)\\.

## References

Ferrari, S. and Cribari-Neto, F. (2004). Beta regression for modelling
rates and proportions. *Journal of Applied Statistics* **31**, 799-815.

## See also

[`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md)
for the same law in the mean and the variance;
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
for the beta as a mixing law over a binomial probability;
[`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md)
for a response at the endpoints rather than between them;
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
and
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
for the wrappers a response with mass at 0 or 1 needs;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[Beta1Distrib](https://statmodels7.github.io/distributions7/reference/Beta1Distrib.md)
for the class.

## Examples

``` r
d <- beta1_distrib()
d
#> Distribution: Beta1
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu  (mean)               | Link: logit      | Domain: (0, 1)
#>   phi (precision)          | Link: log        | Domain: (0, Inf)

# The density is stats::dbeta at shape1 = mu phi, shape2 = (1 - mu) phi.
y <- c(0.2, 0.5, 0.8)
th <- list(mu = 0.4, phi = 5)
all.equal(distrib_pdf(d, y, th), dbeta(y, 2, 3))
#> [1] TRUE

# Moments: variance mu(1 - mu)/(phi + 1), so phi is a precision.
c(mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#>       mean        var       skew       kurt 
#>  0.4000000  0.0400000  0.2857143 -0.6428571 
0.4 * 0.6 / (5 + 1)
#> [1] 0.04

# Both shapes are 1 at mu = 1/2, phi = 2: the beta is the uniform.
distrib_pdf(d, y, list(mu = 0.5, phi = 2))
#> [1] 1 1 1

# Fitting recovers the parameters; the moment estimates start it off.
set.seed(6)
z <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, z)
rbind(fitted  = coef(fit),
      moments = c(mu = mean(z),
                  phi = mean(z) * (1 - mean(z)) / var(z) - 1))
#>                mu      phi
#> fitted  0.3980411 5.155223
#> moments 0.3976899 5.132059

# Every third derivative is free of the response, so it is a constant.
lapply(distrib_deriv3(d, y, th), unique)
#> $mu_mu_mu
#> [1] 31.25
#> 
#> $mu_mu_phi
#> [1] -4.045836
#> 
#> $mu_phi_phi
#> [1] 0.00385982
#> 
#> $phi_phi_phi
#> [1] 0.01036213
#> 
```
