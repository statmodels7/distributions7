
<!-- README.md is generated from README.Rmd. Please edit that file, then
     regenerate with devtools::build_readme(). Do not use knitr::knit(): it
     processes the code but leaves this YAML header in the output as literal
     text, which GitHub and pkgdown both render verbatim. -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/statmodels7/distributions7/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/statmodels7/distributions7/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/statmodels7/distributions7/graph/badge.svg)](https://app.codecov.io/gh/statmodels7/distributions7)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

# distributions7 <img src="man/figures/logo.png" align="right" height="139" alt="" />

Almost every R package that fits statistical models writes its own
distributions, privately: a `switch` on a character string, a list of
closures nothing outside can reach. The same Gamma is therefore
re-implemented again and again, each time only as far as that package
happened to need it: the density and the score, sometimes the Hessian,
rarely anything beyond.

`{distributions7}` writes them once, as objects. Each carries the usual
density, distribution and quantile functions, together with the **exact
derivatives of the log-likelihood up to fourth order**. The derivatives
are available with respect to the parameters, with respect to the
response, and with respect to the unconstrained parameters behind a link
function, which is the scale an optimiser works on.

It is the distribution layer of
[statmodels7](https://statmodels7.github.io), an S7 toolkit for
statistical modelling, and works alongside
[linkfunctions7](https://statmodels7.github.io/linkfunctions7/). The
mathematics behind every formula, from the likelihood theory and the
change of scale to the derivations for the zero-inflated, zero-adjusted,
truncated and transformed wrappers, is worked out in [the statmodels7
book](https://statmodels7.github.io/book/).

## Installation

``` r
# install.packages("pak")
pak::pak("statmodels7/distributions7")
```

## The usual functions

Fourteen distributions ship with the package. Each constructor takes the
link functions used for its parameters, and parameters travel as a named
list.

``` r
d <- gaussian_distrib()
theta <- list(mu = 2, sigma = 3)

distrib_pdf(d, c(0, 2, 4), theta)
#> [1] 0.1064827 0.1329808 0.1064827
distrib_cdf(d, c(0, 2, 4), theta)
#> [1] 0.2524925 0.5000000 0.7475075
distrib_quantile(d, c(0.025, 0.5, 0.975), theta)
#> [1] -3.879892  2.000000  7.879892
```

## Derivatives

The score, the observed and expected information, and the third and
fourth order derivatives are all closed form.

``` r
y <- distrib_rng(d, 5, theta)

distrib_gradient(d, y, theta)
#> $mu
#> [1] -0.20881794  0.06121444 -0.27854287  0.53176027  0.10983592
#> 
#> $sigma
#> [1] -0.2025185 -0.3220917 -0.1005749  0.5149736 -0.2971415
distrib_expected_hessian(d, 0, theta)
#> $mu_mu
#> [1] -0.1111111
#> 
#> $sigma_sigma
#> [1] -0.2222222
#> 
#> $mu_sigma
#> [1] 0
```

On the **link scale** the same generics return derivatives with respect
to the unconstrained parameters instead, computed through the chain rule
rather than by differentiating numerically:

``` r
distrib_gradient(d, y, theta, scale = "link")
#> $mu
#> [1] -0.20881794  0.06121444 -0.27854287  0.53176027  0.10983592
#> 
#> $sigma
#> [1] -0.6075556 -0.9662751 -0.3017248  1.5449208 -0.8914246
```

## Fitting

`fit_distrib()` maximises the likelihood on the link scale, where the
parameters are unconstrained, then reports estimates back on their
natural scale. Confidence intervals are built on the link scale and
mapped through $g^{-1}$, so they can never leave a parameter’s domain.

``` r
y <- distrib_rng(gamma_distrib(), 500, list(mu = 3, sigma2 = 2))
fit <- fit_distrib(gamma_distrib(), y)
fit
#> Maximum-likelihood fit: gamma
#> Observations: 500   Log-likelihood: -841.2   AIC: 1686   BIC: 1695
#> Method: Fisher scoring (converged in 17 iterations)
#> 
#> Parameter scale:
#>        Estimate Std. Error   2.5%  97.5%
#> mu       2.9570     0.0631 2.8359 3.0832
#> sigma2   1.9888     0.1480 1.7188 2.3011
#> 
#> Link scale (log, log):
#>        Estimate Std. Error   2.5%  97.5%
#> mu       1.0842     0.0213 1.0424 1.1260
#> sigma2   0.6875     0.0744 0.5416 0.8334
```

The fit knows what it was estimated from, so it can be checked against
the data and simulated from:

``` r
plot(fit)
```

<img src="man/figures/README-fit-plot-1.png" alt="" width="100%" />

``` r
sims <- simulate(fit, 100, seed = 1)
quantile(vapply(sims, median, numeric(1)), c(0.025, 0.975))
#>     2.5%    97.5% 
#> 2.570081 2.881259
```

The optimisation itself is delegated to
[optimizers7](https://statmodels7.github.io/optimizers7/), so `method`
takes any optimiser of that package as well as the three named
strategies, and the object brings its own stopping rule:

``` r
fit2 <- fit_distrib(gamma_distrib(), y,
                    method = optimizers7::lbfgs(
                      criterion = optimizers7::crit_grad(1e-12)))
c(fisher = as.numeric(logLik(fit)), lbfgs = as.numeric(logLik(fit2)))
#>    fisher     lbfgs 
#> -841.2326 -841.2326
```

## A user-defined distribution needs only its density

Everything above has a numerical fallback registered on the base
classes: the distribution function by quadrature, the quantile function
by root finding, the generator by Generalized Ratio-of-Uniforms, and
every derivative by finite differences. So a new distribution is a
subclass and one method.

``` r
MyLaplace <- S7::new_class("MyLaplace", parent = continuous_distrib)

S7::method(distrib_pdf, MyLaplace) <- function(distrib, y, theta, log = FALSE) {
  ld <- -log(2 * theta$b) - abs(y - theta$mu) / theta$b
  if (log) ld else exp(ld)
}

d2 <- MyLaplace(
  distrib_name = "my laplace", dimension = "univariate", bounds = c(-Inf, Inf),
  params = c("mu", "b"), params_interpretation = c(mu = "location", b = "scale"),
  n_params = 2, params_bounds = list(mu = c(-Inf, Inf), b = c(0, Inf)),
  link_params = list(
    mu = linkfunctions7::identity_link(),
    b  = linkfunctions7::log_link()
  )
)

# never implemented, yet available
distrib_cdf(d2, 1, list(mu = 0, b = 2))
#> [1] 0.6967347
distrib_gradient(d2, 1, list(mu = 0, b = 2))
#> $mu
#> [1] 0.5
#> 
#> $b
#> [1] -0.25
```

See `vignette("defining-a-distribution")` for the full treatment,
including what to do when a parameter is not differentiable.

## Validating a distribution

`check_distrib()` puts a distribution through thirteen numerical checks.
It verifies that the density integrates to one, that the distribution
function agrees with it, that the quantile function inverts it, and that
the generator follows it. It also checks every derivative order against
finite differences, on both scales.

``` r
invisible(check_distrib(d2, list(mu = 0, b = 2), nsim = 2e4))
#> Distribution: my laplace
#> Parameters:   mu = 0, b = 2
#> Observations: 100   Monte Carlo: 20000
#> 
#>   [OK  ] density integrates to 1                     1.81e-10
#>   [OK  ] density is non-negative                     5.00e-03
#>   [OK  ] cdf in [0,1] and non-decreasing             2.46e-02
#>   [OK  ] cdf agrees with the density                 8.42e-06
#>   [OK  ] quantile/cdf round-trip                     3.31e-10
#>   [OK  ] rng matches the cdf                         1.78e+00
#>   [OK  ] gradient vs finite differences              0.00e+00
#>   [OK  ] hessian vs finite differences               0.00e+00
#>   [OK  ] deriv3 vs finite differences                0.00e+00
#>   [OK  ] deriv4 vs finite differences                0.00e+00
#>   [OK  ] expected information vs Monte Carlo         1.97e+00
#>   [OK  ] response derivatives vs finite differences  0.00e+00
#>   [OK  ] link-scale gradient vs finite differences   3.29e-09
#> 
#> All 13 checks passed.
```

## What is in the box

|  |  |
|----|----|
| continuous | gaussian, cauchy, logistic, Student’s t, Laplace, pseudo-Huber, gamma, inverse gaussian, lognormal, beta |
| discrete | bernoulli, binomial, poisson, negative binomial |
| wrappers | `zero_inflated()`, `zero_adjusted()`, `truncated()`, `transformation()` with twelve transformers |
| tools | `fit_distrib()`, `check_distrib()`, `expectation()`, moments, `rng_grou()` |
