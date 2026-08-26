# Expected Information of a Multivariate Distribution, by Sampling

The fallback for a multivariate family with no closed form: the
expectation of the observed Hessian, taken over a sample drawn from the
family. All four shipped families register their own closed forms, so
this runs only for a family written elsewhere.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object with no closed form of its own.

- y:

  An \\n \times p\\ numeric matrix of observations; only its row count
  is used, the expectation being over the law.

- theta:

  A named list of parameters, each component a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch.

- approx:

  One of `"bartlett"` (the default, equivalently `"opg"`) or `"mc"`.
  `"integrate"` is not available here and falls through to sampling.

- nsim:

  The Monte Carlo sample size. Defaults to `10000`. The error falls as
  its square root, so ten times the accuracy costs a hundred times the
  draws.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each vector constant.

## Why the one-dimensional routes do not survive

[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)'s
`"integrate"` splits a one-dimensional integral at quantiles, and there
are no quantiles here; `"bartlett"` would need the score's own higher
derivatives assembled over a \\p\\-dimensional partition sum. What is
left is sampling, and both admissible values of `approx` route to it.

## The result is a sample, not an integral

It carries Monte Carlo error of order `nsim^(-1/2)`, and two calls under
different seeds give two answers. Set a seed before calling if the
result must be reproducible, and expect a fit that inverts this matrix
to move with it. `"bartlett"` averages the outer product of the score
and `"mc"` averages the observed Hessian; the two agree in expectation
by the second Bartlett identity, and differ by sampling error at any
finite `nsim`.

## See also

[`distrib_expected_hessian.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MvGaussianDistrib.md)
and
[`distrib_expected_hessian.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MvStudentTDistrib.md)
for the closed forms,
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md),
which says which route a family takes, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
# Every shipped family overrides, so reach the fallback directly.
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)
base <- S7::method(distrib_expected_hessian, multivariate_distrib)

set.seed(1)
sampled <- base(d, matrix(0, 3, 2), distributions7:::align_theta(d, theta),
                approx = "mc", nsim = 4000)
exact <- distrib_expected_hessian(d, matrix(0, 3, 2), theta)
round(rbind(sampled = vapply(sampled, function(z) z[1], numeric(1)),
            closed = vapply(exact, function(z) z[1], numeric(1))), 3)
#>         mu1_mu1 mu2_mu2 sigma_log_L1_sigma_log_L1 sigma_log_L2_sigma_log_L2
#> sampled   -1.25      -1                    -2.415                    -2.008
#> closed    -1.25      -1                    -2.250                    -2.000
#>         sigma_L2.1_sigma_L2.1 mu1_mu2 mu1_sigma_log_L1 mu1_sigma_log_L2
#> sampled                -1.073     0.5           -0.009           -0.014
#> closed                 -1.000     0.5            0.000            0.000
#>         mu1_sigma_L2.1 mu2_sigma_log_L1 mu2_sigma_log_L2 mu2_sigma_L2.1
#> sampled          0.015            0.001            0.028         -0.001
#> closed           0.000            0.000            0.000          0.000
#>         sigma_log_L1_sigma_log_L2 sigma_log_L1_sigma_L2.1
#> sampled                    -0.003                   0.539
#> closed                      0.000                   0.500
#>         sigma_log_L2_sigma_L2.1
#> sampled                   0.006
#> closed                    0.000

# Two seeds give two answers, this being a sample.
set.seed(2)
again <- base(d, matrix(0, 3, 2), distributions7:::align_theta(d, theta),
              approx = "mc", nsim = 4000)
c(first = sampled$mu1_mu1[1], second = again$mu1_mu1[1])
#>  first second 
#>  -1.25  -1.25 
```
