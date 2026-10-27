# Information Matrices

The information of a likelihood is the matrix of second derivatives that
measures how sharply the log-likelihood \\\ell(\eta)\\ peaks at the
estimate \\\hat\eta\\ of the parameters on the link scale. Its inverse
is the estimated variance matrix of the estimates, so the choice of the
matrix is the choice of the standard errors. Two objects select it,
[`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md)
and
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md),
and both inherit from
[Information](https://statmodels7.github.io/distributions7/reference/Information.md).
A function that takes such an object, such as the
[`vcov()`](https://rdrr.io/r/stats/vcov.html) method of a fit (see
[`vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md)),
is told the matrix by the object and not by a character string. The
argument `information` of those functions is `NULL` by default, which
selects the matrix described under "The default" below.

## Value

No value. The topic describes the two objects, of classes
[ObservedInfo](https://statmodels7.github.io/distributions7/reference/Information.md)
and
[ExpectedInfo](https://statmodels7.github.io/distributions7/reference/Information.md),
that
[`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md)
and
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md)
return.

## Details

- [`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md):

  The negative Hessian of the log-likelihood at the data,
  \$\$-\frac{\partial^2
  \ell}{\partial\eta\\\partial\eta^\top}\bigg\|\_{\hat\eta}.\$\$ It is
  computed from the second derivatives that the family supplies in
  closed form, or from their numerical fallbacks when the family
  supplies only a density. It takes no argument. It is the curvature of
  the likelihood in the sample at hand, and every family has it. It
  describes the uncertainty only where the log-likelihood is twice
  differentiable in every parameter at the estimate.

- [`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md):

  The Fisher information, \$\$\mathbb{E}\left\[-\frac{\partial^2
  \ell}{\partial\eta\\\partial\eta^\top} \right\]\_{\hat\eta},\$\$ the
  average of the observed information over samples drawn from the model
  at the estimates. With `method = NULL`, the default, the closed form
  of the family is used where the family writes the expectation out (see
  [`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)),
  and
  [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
  otherwise. An object in `method` is used for every family, one of
  [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md),
  [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
  and
  [`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md);
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
  compares the three.

The two matrices agree asymptotically and differ in a finite sample. The
expected information is a function of the parameters alone and is
positive definite wherever the family is regular. The observed
information depends on the data, and it can be indefinite away from the
estimate or at an outlier. Standard errors from the observed information
describe the likelihood that was actually observed, which is the usual
recommendation (Efron and Hinkley, 1978). For a family that has no
closed-form expected information,
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
is a noisier estimate of the information than the observed matrix.

## The default

The argument `information` of
[`vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md),
[`confint.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md),
[`print.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/print.distrib_fit.md)
and
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
is `NULL`. It resolves to
[`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md)
when every parameter of the family has
[`params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.md)
of at least 2, that is, when the log-likelihood is twice differentiable
in every parameter, and to
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md)
with `method = NULL` otherwise, for the whole matrix. A family with a
kink in a parameter (the Laplace, the elastic net) has a singular or
undefined observed information in that parameter, and the closed-form
expected information gives finite standard errors. A fit states which
information it used whenever it is not the observed one.

## References

Efron, B. and Hinkley, D. V. (1978). Assessing the accuracy of the
maximum likelihood estimator: observed versus expected Fisher
information. *Biometrika*, 65, 457-483.

## See also

[`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md)
and
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md)
for the constructors;
[Information](https://statmodels7.github.io/distributions7/reference/Information.md)
for the classes;
[`vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md),
[`confint.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md)
and
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
for the functions that take the object;
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for the methods of the expected information.

## Examples

``` r
observed_info()
#> Observed information
expected_info()
#> Expected information (closed form where the family has one, otherwise the outer product of the scores)
expected_info(method = expect_integrate(identity = "score"))
#> Expected information (expectation by integration against the density (identity = score))

# The standard errors of one fit under the two matrices.
set.seed(1)
d <- gaussian1_distrib()
fit <- fit_distrib(d, distrib_rng(d, 300, list(mu = 1, sigma = 2)))
sqrt(diag(vcov(fit, information = observed_info())))
#>         mu      sigma 
#> 0.11109240 0.07855419 
sqrt(diag(vcov(fit, information = expected_info())))
#>         mu      sigma 
#> 0.11109240 0.07855419 
```
