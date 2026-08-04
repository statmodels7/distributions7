# Interpretable Estimates of a Multivariate Fit

Reports the standard deviations, correlations and whatever else
[`mv_derived`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
declares for the fitted distribution, each with its standard error and
confidence interval.

## Usage

``` r
mv_summary(object, level = object@level)
```

## Arguments

- object:

  A
  [`distrib_fit`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
  of a multivariate distribution.

- level:

  The confidence level. Defaults to the fit's own.

## Value

A data frame with one row per quantity and the columns `Estimate`,
`Std. Error` and the two confidence limits, carrying the attribute
`"block"` that names the group each row belongs to.

## Details

A multivariate fit estimates the free values of a covstructs7 structure,
and those are coordinates rather than quantities: the estimate and
standard error of `sigma_log_L2` answer no question anybody asked. This
function carries the fit's variance matrix onto the quantities that do,
by the delta method, \$\$\widehat{\mathrm{Var}}\\g(\hat\theta)\\ = J \\
\widehat{\mathrm{Var}}(\hat\theta) \\ J^\top, \qquad J = \partial
g/\partial\theta,\$\$ with \\J\\ taken from
[`mv_derived`](https://statmodels7.github.io/distributions7/reference/mv_derived.md),
in closed form for the families that ship with the package.

Each interval is built on the scale the quantity declares and mapped
back, which is the same discipline
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
applies to a univariate parameter: a standard deviation is intervalled
on the log scale and a correlation on Fisher's \\z\\, so neither
interval can leave the set its quantity lives in. An interval on the raw
scale would routinely put a correlation above one.

## See also

[`mv_derived`](https://statmodels7.github.io/distributions7/reference/mv_derived.md),
[`confint.distrib_fit`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md)

## Examples

``` r
set.seed(1)
d <- mvgaussian_distrib(2)
y <- distrib_rng(d, 500, list(mu1 = 0, mu2 = 1, sigma_log_L1 = 0,
                              sigma_log_L2 = 0, sigma_L2.1 = 0.7))
fit <- fit_distrib(d, y)

# the standard deviations and the correlation, which is what one reads
mv_summary(fit)
#>            Estimate Std. Error      2.5%     97.5%
#> sd_v1     1.0109159 0.02109291 0.9704085 1.0531142
#> sd_v2     1.2469079 0.02698994 1.1951150 1.3009453
#> cor_v1_v2 0.5325916 0.02135920 0.4894342 0.5731465
```
