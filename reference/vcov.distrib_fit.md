# Variance-Covariance Matrix of a Maximum-Likelihood Fit

Returns the estimated variance matrix of the estimates, on either scale.
The one the fit computes is on the link scale: the inverse of the
information at \\\hat\eta\\, the expected information where the fit used
it or the family writes it out, and the observed Hessian otherwise. The
parameter-scale matrix is its image under the delta method,
\$\$\widehat{\mathrm{Var}}(\hat\theta) =
J\\\widehat{\mathrm{Var}}(\hat\eta)\\J, \qquad J =
\mathrm{diag}\\\left(\frac{dg^{-1}}{d\eta}\Big\|\_{\hat\eta}\right),\$\$
the Jacobian being diagonal because each parameter carries its own
scalar link.

Every entry is `NA` when the information could not be evaluated or
inverted at the optimum. The estimates stand in that case and only the
uncertainty is missing.

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
  [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html).

## Value

A symmetric numeric matrix of dimension `length(object@distrib@params)`,
with both dimnames set to the parameter names. Its diagonal is the
square of `object@se` on the parameter scale and of `object@se_eta` on
the link scale.

## See also

[`coef.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/coef.distrib_fit.md)
and
[`confint.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md),
which take the same `scale`;
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the information itself.

## Examples

``` r
set.seed(1)
d <- gaussian1_distrib()
fit <- fit_distrib(d, distrib_rng(d, 400, list(mu = 1, sigma = 2)))

vcov(fit)
#>                mu       sigma
#> mu    0.009376678 0.000000000
#> sigma 0.000000000 0.004688339
sqrt(diag(vcov(fit)))          # the standard errors the fit reports
#>         mu      sigma 
#> 0.09683325 0.06847145 
all.equal(sqrt(diag(vcov(fit))), fit@se)
#> [1] TRUE

# The parameter-scale matrix is the link-scale one under the delta method.
J <- diag(c(1, coef(fit)[["sigma"]]))   # d theta / d eta: identity, then exp
all.equal(unname(vcov(fit)), J %*% vcov(fit, "link") %*% J)
#> [1] TRUE

# The location and the scale of a Gaussian are orthogonal, so the
# off-diagonal entry is zero rather than merely small.
vcov(fit)[["mu", "sigma"]]
#> [1] 0
```
