# Log-Likelihood of a Maximum-Likelihood Fit

Returns the maximized log-likelihood as a `logLik` object, so that
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html),
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) and anything else in
stats that reads one can be applied to a fit. The value is summed over
observations and is **not** divided by \\n\\, whatever scaling the
optimizer worked with.

The `df` attribute is the number of estimated parameters, which for
these fits is every parameter of the family: a fit estimates all of
them, so nothing is held. The `nobs` attribute is `object@n`, the row
count for a multivariate response and the length otherwise.

## Arguments

- object:

  A
  [`distrib_fit()`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
  object.

- ...:

  Unused, accepted for compatibility with
  [`stats::logLik()`](https://rdrr.io/r/stats/logLik.html).

## Value

An object of class `logLik`: a single number carrying the attributes
`df` (the parameter count) and `nobs` (the observation count).

## See also

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
for the fit; [`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html), which read this.

## Examples

``` r
set.seed(1)
d <- gaussian1_distrib()
y <- distrib_rng(d, 400, list(mu = 1, sigma = 2))
fit <- fit_distrib(d, y)

logLik(fit)
#> 'log Lik.' -831.9624 (df=2)
c(df = attr(logLik(fit), "df"), nobs = attr(logLik(fit), "nobs"))
#>   df nobs 
#>    2  400 

# The criteria the fit reports are the ones stats builds from this.
all.equal(AIC(logLik(fit)), fit@aic)
#> [1] TRUE
all.equal(BIC(logLik(fit)), fit@bic)
#> [1] TRUE

# Comparing two families on the same data. The Student t spends one more
# parameter, so AIC decides whether the tails are worth it.
ft <- fit_distrib(student_t1_distrib(), y)
c(gaussian = AIC(logLik(fit)), student_t = AIC(logLik(ft)))
#>  gaussian student_t 
#>  1667.925  1669.925 
```
