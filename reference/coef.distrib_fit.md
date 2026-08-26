# Estimates From a Maximum-Likelihood Fit

Returns the maximum likelihood estimates, on either of the two scales
the fit carries. The default is the parameter scale, \\\hat\theta\\,
which is what the family is interpreted in; `scale = "link"` gives
\\\hat\eta = g(\hat\theta)\\, the point the optimizer actually reached.

Neither is recomputed. Both were stored at the optimum, and each is the
image of the other under the family's links, so `coef(fit, "link")` is
`g(coef(fit))` component by component.

## Arguments

- object:

  A
  [`distrib_fit()`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
  object.

- scale:

  Either `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html); anything
  else signals an error.

- ...:

  Unused, accepted for compatibility with
  [`stats::coef()`](https://rdrr.io/r/stats/coef.html).

## Value

A named numeric vector of length `length(object@distrib@params)`, named
and ordered as `object@distrib@params`.

## See also

[`vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md)
and
[`confint.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md),
which take the same `scale`;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
for the fit itself.

## Examples

``` r
set.seed(1)
d <- gaussian1_distrib()
fit <- fit_distrib(d, distrib_rng(d, 400, list(mu = 1, sigma = 2)))

coef(fit)
#>       mu    sigma 
#> 1.076177 1.936665 
coef(fit, scale = "link")
#>        mu     sigma 
#> 1.0761773 0.6609674 

# The scale carries a log link and the location the identity, so the link
# scale is the logarithm of one estimate and the other unchanged.
all.equal(coef(fit, "link")[["sigma"]], log(coef(fit)[["sigma"]]))
#> [1] TRUE
all.equal(coef(fit, "link")[["mu"]], coef(fit)[["mu"]])
#> [1] TRUE
```
