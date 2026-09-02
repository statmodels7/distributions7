# Interpretable Estimates of a Multivariate Fit

Reports the standard deviations, correlations and whatever else
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
declares for the fitted distribution, each with its standard error and
its confidence interval. A multivariate fit estimates the free values of
a parameters7 parametrization, and those are coordinates: the estimate
and standard error of `sigma_log_L2` answer no question anybody asked.
This function carries the fit's variance matrix onto the quantities that
do.

## Usage

``` r
mv_summary(object, level = object@level)
```

## Arguments

- object:

  A
  [`distrib_fit()`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
  of a MULTIVARIATE distribution. A univariate fit is rejected: there
  the parameters are already the interpretable quantities and
  [`confint.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md)
  is the function to call.

- level:

  The confidence level, a single number strictly inside \\(0, 1)\\.
  Defaults to the fit's own `level`. Anything else is an error.

## Value

A data frame with one row per quantity and four columns: `Estimate`,
`Std. Error`, and the two confidence limits, named after the percentiles
they are (`2.5%` and `97.5%` at the default level). Its row names are
the quantity names, and it carries the attribute `"block"`, a character
vector naming the group each row belongs to, by which
[`print()`](https://rdrr.io/r/base/print.html) lays the summary out.

## The delta method

\$\$\widehat{\mathrm{Var}}\\g(\hat\theta)\\ = J \\
\widehat{\mathrm{Var}}(\hat\theta) \\ J^\top, \qquad J = \partial
g/\partial\theta,\$\$ with \\J\\ taken from
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md),
in closed form for the two families that ship and by one central
difference for any other.

## Where each interval is built

On the scale the quantity declares, then mapped back, which is the
discipline
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
applies to a univariate parameter. A standard deviation is intervalled
on the log scale and a correlation on Fisher's \\z\\, so neither
interval can leave the set its quantity lives in. An interval on the raw
scale routinely puts a correlation above one.

The standard error is carried onto that scale by the transform's own
derivative, so the reported `Std. Error` column stays on the quantity's
natural scale while the limits are computed on the transformed one.

## Notation

\\\theta\\ is the parameter vector, \\g\\ the map to the reported
quantities, \\J\\ its Jacobian and \\z\\ Fisher's transform of a
correlation.

## See also

[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
for the quantities and the Jacobian,
[`confint.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md)
for the coordinates themselves, and
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
for the fit.

## Examples

``` r
set.seed(1)
d <- mvgaussian1_distrib(2)
truth <- list(mu1 = 0, mu2 = 1, sigma_log_L1 = 0,
              sigma_log_L2 = 0, sigma_L2.1 = 0.7)
y <- distrib_rng(d, 500, truth)
fit <- fit_distrib(d, y)

# The standard deviations and the correlation, which is what one reads.
mv_summary(fit)
#>            Estimate Std. Error      2.5%     97.5%
#> sd_v1     1.0109159 0.03196797 0.9501620 1.0755544
#> sd_v2     1.2469079 0.03943069 1.1719714 1.3266358
#> cor_v1_v2 0.5325916 0.03203598 0.4669039 0.5924342

# The coordinates the fit was estimated on, beside them.
confint(fit)
#>                      2.5%      97.5%
#> mu1          -0.065964967 0.11125314
#> mu2           0.860616204 1.07920478
#> sigma_log_L1 -0.051122766 0.07283624
#> sigma_log_L2 -0.008108583 0.11585042
#> sigma_L2.1    0.562845028 0.76534030

# Every interval stays inside the set its quantity lives in: the standard
# deviations are positive and the correlation is inside (-1, 1).
s <- mv_summary(fit)
c(sd_lower_positive = all(s[1:2, 3] > 0),
  cor_inside = s[3, 3] > -1 && s[3, 4] < 1)
#> sd_lower_positive        cor_inside 
#>              TRUE              TRUE 

# A wider level widens the interval and leaves the estimate alone.
mv_summary(fit, level = 0.99)
#>            Estimate Std. Error      0.5%     99.5%
#> sd_v1     1.0109159 0.03196797 0.9318363 1.0967065
#> sd_v2     1.2469079 0.03943069 1.1493676 1.3527258
#> cor_v1_v2 0.5325916 0.03203598 0.4450908 0.6100184

# A univariate fit is rejected by name.
try(mv_summary(fit_distrib(gaussian1_distrib(), rnorm(50))))
#> Error : mv_summary() is for a multivariate fit. For a univariate one the
#>   parameters are already the interpretable quantities; use confint().
```
