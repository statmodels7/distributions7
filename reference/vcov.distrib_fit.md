# Variance-Covariance Matrix of a Maximum-Likelihood Fit

Returns the estimated variance matrix of the estimates, on either scale.
The matrix that the fit computes is on the link scale and is the inverse
of the information at \\\hat\eta\\. The parameter-scale matrix is its
image under the delta method, \$\$\widehat{\mathrm{Var}}(\hat\theta) =
J\\\widehat{\mathrm{Var}}(\hat\eta)\\J, \qquad J =
\mathrm{diag}\\\left(\frac{dg^{-1}}{d\eta}\Big\|\_{\hat\eta}\right),\$\$
and the Jacobian is diagonal because each parameter carries its own
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

- information:

  Which information to invert: `"fit"` (the default) is the matrix that
  the fit stored, `"observed"` the Hessian at the estimates, and
  `"expected"` the expected information there. The last two are
  recomputed.

- approx:

  One of `"opg"`, `"bartlett"`, `"integrate"` or `"mc"`, read only when
  `information = "expected"` and the family has no closed form. Defaults
  to `"opg"`.

- nsim:

  Number of draws, read only by `approx = "mc"`.

- ...:

  Unused, accepted for compatibility with
  [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html).

## Value

A symmetric numeric matrix of dimension `length(object@distrib@params)`,
with both dimnames set to the parameter names. With
`information = "fit"` its diagonal is the square of `object@se` on the
parameter scale and of `object@se_eta` on the link scale.

## Which information

`information` selects the matrix that is inverted, and the one that the
fit stored is the default. That one is the expected information where
the family writes it out, where it costs one evaluation and has the
smaller variance, and the observed Hessian otherwise.

The alternative is available because the two are different estimators
and the choice is the reader's. They agree asymptotically but not in a
sample. The expected information is averaged over the model, and the
observed information is the curvature of the likelihood at the data in
hand.

Where a family does not write its expected information out,
`information = "expected"` reaches a fallback, and `approx` selects it.
The default `"opg"` is the outer product of the observed scores, which
costs one gradient; `"bartlett"` evaluates the expectation itself, a sum
over the support for a discrete family and a quadrature for a continuous
one, and is orders of magnitude dearer. A fit does not report standard
errors from the expected information for such a family, because the
outer product of the scores is a noisier estimate of the information
than the observed Hessian. The expensive route is available through
`approx`, and it is not the default.

A matrix requested through `information` is recomputed at the optimum,
so it costs one evaluation of that information and is not read from the
fit.

## See also

[`coef.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/coef.distrib_fit.md)
and
[`confint.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md),
which take the same `scale`;
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the information itself, and
[`expected_by_opg()`](https://statmodels7.github.io/distributions7/reference/expected_by_opg.md)
for what `approx = "opg"` computes.

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

# The two informations are different estimators of one quantity: close on
# a well-specified fit of this size, and not identical.
sqrt(diag(vcov(fit, information = "observed")))
#>         mu      sigma 
#> 0.09683325 0.06847145 
sqrt(diag(vcov(fit, information = "expected")))
#>         mu      sigma 
#> 0.09683325 0.06847145 
```
