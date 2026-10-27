# Variance-Covariance Matrix of a Maximum-Likelihood Fit

Returns the estimated variance matrix of the estimates, on either scale.
The matrix is computed on the link scale as the inverse of an
information matrix at \\\hat\eta\\. The parameter-scale matrix is its
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

  `NULL` (the default),
  [`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md)
  or
  [`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md),
  the information matrix to invert; `NULL` selects the default of the
  family described above.
  [`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md)
  carries the `method` that obtains the expectation. Any other value
  signals an error. See
  [information_objects](https://statmodels7.github.io/distributions7/reference/information_objects.md).

- ...:

  Unused, accepted for compatibility with
  [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html).

## Value

A symmetric numeric matrix of dimension `length(object@distrib@params)`,
with both dimnames set to the parameter names. With the default
`information` its diagonal is the square of `object@se` on the parameter
scale and of `object@se_eta` on the link scale.

## Which information

The `information` argument selects the matrix that is inverted, after
the fit and independently of the `method` that computed the estimates.
It is `NULL`,
[`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md)
or
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md)
(see
[information_objects](https://statmodels7.github.io/distributions7/reference/information_objects.md)):

- `NULL`, the default:
  [`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md)
  when every parameter of the family has
  [`params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.md)
  of at least 2, and
  [`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md)
  with `method = NULL` otherwise, for the whole matrix. A family with a
  kink in a parameter (the Laplace, the elastic net) has a singular or
  undefined observed information in that parameter, and the closed-form
  expected information gives finite standard errors.

- [`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md):
  the negative Hessian of the log-likelihood at the data. It exists for
  every family, and it describes the curvature of the likelihood in the
  sample at hand.

- [`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md):
  the Fisher information. With `method = NULL` the closed form of the
  family is used where there is one, and
  [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
  otherwise. An object in `method` is used for every family, one of
  [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md),
  [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
  and
  [`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
  (see
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)).

The two matrices agree asymptotically and differ in a sample. The
expected information is averaged over the model, and the observed
information is the curvature of the likelihood in the sample at hand.
For a family with no closed-form expected information,
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
is the outer product of the observed scores, which costs one gradient
and is a noisier estimate of the information than the observed Hessian;
`expect_integrate(identity = "score")` evaluates the expectation itself,
a sum over the support for a discrete family and a quadrature for a
continuous one, and is dearer.

The fit stores the matrix of the default information and returns it as
stored. Any other information is recomputed at the optimum, so it costs
one evaluation of that matrix and no refit.

## See also

[`coef.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/coef.distrib_fit.md)
and
[`confint.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md),
which take the same `scale`;
[information_objects](https://statmodels7.github.io/distributions7/reference/information_objects.md)
for the two information objects;
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the information itself, and
[`expected_by_opg()`](https://statmodels7.github.io/distributions7/reference/expected_by_opg.md)
for what `method = expect_opg()` computes.

## Examples

``` r
set.seed(1)
d <- gaussian1_distrib()
fit <- fit_distrib(d, distrib_rng(d, 400, list(mu = 1, sigma = 2)))

vcov(fit)
#>                  mu         sigma
#> mu     9.376678e-03 -4.213744e-19
#> sigma -4.213744e-19  4.688339e-03
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
#> [1] -4.213744e-19

# The two informations are different estimators of one quantity: close on
# a well-specified fit of this size, and not identical.
sqrt(diag(vcov(fit, information = observed_info())))
#>         mu      sigma 
#> 0.09683325 0.06847145 
sqrt(diag(vcov(fit, information = expected_info())))
#>         mu      sigma 
#> 0.09683325 0.06847145 
```
