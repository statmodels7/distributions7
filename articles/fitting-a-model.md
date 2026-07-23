# Fitting a distribution to data

``` r

library(distributions7)
```

This vignette walks through fitting a distribution to a sample by
maximum likelihood, reading the result, checking it against the data,
and simulating from it. It assumes you already know what a distribution
object is; if not,
[`vignette("defining-a-distribution")`](https://statmodels7.github.io/distributions7/articles/defining-a-distribution.md)
is the place to start.

## A first fit

Every distribution has a constructor. It takes the link functions for
its parameters, but the defaults are sensible, so most of the time you
call it with no arguments:

``` r

d <- gamma_distrib()
d
#> Distribution: Gamma
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu     (mean)               | Link: log        | Domain: (0, Inf)
#>   sigma2 (variance)           | Link: log        | Domain: (0, Inf)
```

We simulate some data from a known parameter value and fit the
distribution back to it. Parameters travel as a named list.

``` r

y <- distrib_rng(d, 500, list(mu = 3, sigma2 = 2))
fit <- fit_distrib(d, y)
fit
#> Maximum-likelihood fit: gamma
#> Observations: 500   Log-likelihood: -844   AIC: 1692   BIC: 1700
#> Method: Fisher scoring (converged in 4 iterations)
#> 
#> Parameter scale (95% CI mapped from the link scale):
#>        Estimate Std. Error   2.5%  97.5%
#> mu       2.9706     0.0634 2.8488 3.0976
#> sigma2   2.0120     0.1498 1.7388 2.3281
#> 
#> Link scale (log, log):
#>        Estimate Std. Error
#> mu       1.0888     0.0214
#> sigma2   0.6991     0.0745
```

The estimates recover the values we generated from, and the print method
shows them twice: on the **parameter scale** you asked about, and on the
**link scale** where the fitting actually happened. More on that split
below.

## Reading the fit

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
returns an object you interrogate with the familiar extractors:

``` r

coef(fit)
#>       mu   sigma2 
#> 2.970589 2.012009
vcov(fit)
#>                 mu      sigma2
#> mu     0.004024019 0.005451015
#> sigma2 0.005451015 0.022443705
logLik(fit)
#> 'log Lik.' -844.0323 (df=2)
```

[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) also take a `scale`
argument, in case you want the quantities on the unconstrained scale:

``` r

coef(fit, scale = "link")
#>        mu    sigma2 
#> 1.0887603 0.6991339
```

The object carries confidence intervals, the information criteria, the
number of iterations, and the method that was actually used — all
visible in the print above.

## Why the link scale

A gamma’s mean and variance are positive; a beta’s probability lives in
$`(0, 1)`$; a Gaussian’s standard deviation is positive. Optimising
those constrained parameters directly means fighting the boundary.
Instead
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
optimises the **unconstrained** parameters behind each link —
$`\log\mu`$ and $`\log\sigma^2`$ for the gamma — where there is no
boundary to hit, using the exact score and information the package
provides on that scale.

The estimates are then mapped back, and so are the intervals. Because a
confidence interval is built symmetrically on the link scale and pushed
through $`g^{-1}`$, its endpoints can never leave the parameter’s
domain:

``` r

b <- bernoulli_distrib()
fit_distrib(b, rbinom(40, 1, 0.9))
#> Maximum-likelihood fit: bernoulli
#> Observations: 40   Log-likelihood: -10.66   AIC: 23.31   BIC: 25
#> Method: Fisher scoring (converged in 4 iterations)
#> 
#> Parameter scale (95% CI mapped from the link scale):
#>    Estimate Std. Error   2.5%  97.5%
#> mu    0.925     0.0416 0.7918 0.9756
#> 
#> Link scale (logit):
#>    Estimate Std. Error
#> mu   2.5123     0.6003
```

Even with a probability estimate close to one and only forty
observations, the interval stays inside $`(0, 1)`$ — a symmetric
interval on the probability scale would have spilled past it.

## Choosing how it is optimised

The default is Fisher scoring, which uses the *expected* information.
Two other methods are available, and Fisher scoring and Newton fall back
to BFGS if they fail to converge:

``` r

fit_newton <- fit_distrib(d, y, method = "newton")
fit_bfgs   <- fit_distrib(d, y, method = "bfgs")

c(scoring = as.numeric(logLik(fit)),
  newton  = as.numeric(logLik(fit_newton)),
  bfgs    = as.numeric(logLik(fit_bfgs)))
#>   scoring    newton      bfgs 
#> -844.0323 -844.0323 -844.0323
```

They agree on the maximum. Fisher scoring is the default because it uses
the expected information, which stays well-behaved even where the
observed Hessian does not — a point that matters for distributions with
a non-smooth parameter, such as the Laplace.

## Checking the fit against the data

A fit remembers the data it came from, so it can be drawn over it. For a
continuous distribution
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows a kernel
density of the sample with the fitted density on top:

``` r

plot(fit)
```

![](fitting-a-model_files/figure-html/unnamed-chunk-8-1.png)

For a discrete distribution it shows the observed frequencies as bars
with the fitted probability mass overlaid, because a kernel density
would misrepresent counts:

``` r

p <- poisson_distrib()
yp <- distrib_rng(p, 300, list(mu = 4))
plot(fit_distrib(p, yp))
```

![](fitting-a-model_files/figure-html/unnamed-chunk-9-1.png)

## Simulating from the fit

[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws new samples
from the fitted distribution. Each replicate has the same length as the
original data, which makes it a one-liner to build a parametric
bootstrap of any statistic:

``` r

sims <- simulate(fit, 200, seed = 1)
dim(sims)
#> [1] 500 200

boot_median <- vapply(sims, median, numeric(1))
quantile(boot_median, c(0.025, 0.975))
#>     2.5%    97.5% 
#> 2.597115 2.869105
```

It follows the
[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html) contract: a
supplied `seed` makes the draw reproducible and is reported back on the
result, and the caller’s random stream is left untouched.

## Validating a distribution numerically

Before trusting a fit — and especially before trusting a distribution
you wrote yourself —
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
puts it through thirteen numerical checks. It confirms that the density
integrates to one, that the distribution function agrees with the
density, that the quantile inverts it, that the generator follows it,
and that every analytical derivative matches finite differences, on both
the parameter and the link scale:

``` r

invisible(check_distrib(d, list(mu = 3, sigma2 = 2)))
#> Distribution: gamma
#> Parameters:   mu = 3, sigma2 = 2
#> Observations: 100   Monte Carlo: 200000
#> 
#>   [OK  ] density integrates to 1                     5.58e-10
#>   [OK  ] density is non-negative                     1.07e-02
#>   [OK  ] cdf in [0,1] and non-decreasing             2.46e-02
#>   [OK  ] cdf agrees with the density                 2.32e-11
#>   [OK  ] quantile/cdf round-trip                     5.55e-17
#>   [OK  ] rng matches the cdf                         1.16e+00
#>   [OK  ] gradient vs finite differences              6.49e-11
#>   [OK  ] hessian vs finite differences               5.78e-08
#>   [OK  ] deriv3 vs finite differences                1.75e-10
#>   [OK  ] deriv4 vs finite differences                7.01e-08
#>   [OK  ] expected information vs Monte Carlo         1.76e+00
#>   [OK  ] response derivatives vs finite differences  2.70e-08
#>   [OK  ] link-scale gradient vs finite differences   5.59e-08
#> 
#> All 13 checks passed.
```

If you had made a mistake in an analytical derivative, the corresponding
line would turn red, with the size of the discrepancy next to it. It is
the fastest way to catch a wrong formula.

## Where to go next

- [`vignette("defining-a-distribution")`](https://statmodels7.github.io/distributions7/articles/defining-a-distribution.md)
  — add your own distribution, which then works with everything shown
  here.
- [`?fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
  [`?check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
  — the full argument lists.
- [`?distrib_gradient`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
  — the derivatives, including what `scale = "link"` does.
