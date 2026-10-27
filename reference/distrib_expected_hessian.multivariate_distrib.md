# Expected Information of a Multivariate Distribution, by Sampling

Computes the expected information of a multivariate family without a
closed form as an average over a sample drawn from the family. Every
shipped multivariate class registers its own closed form, so this method
runs only for a family written elsewhere.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object without a closed form of its own.

- y:

  An \\n \times p\\ numeric matrix of observations; only its row count
  is used, because the expectation is over the law.

- theta:

  A named list of parameters, in which each component is a single
  number.

- scale:

  One of `"parameter"` (the default) or `"link"`. The generic applies
  the link scale to the result of the method.

- method:

  `NULL` (the default), which selects
  [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md),
  or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).
  [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
  and
  [`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
  sample, and
  [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
  signals an error. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
in which each vector is constant.

## The routes

The `"integrate"` route of
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)
splits a one-dimensional integral at quantiles, and a multivariate
family has no quantiles, so `method = expect_integrate()` signals an
error. The other methods sample:
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
with `identity = "score"` average the negative outer product of the
score over simulated draws, and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
with `identity = "hessian"` averages the observed Hessian. The two
averages agree in expectation by the second Bartlett identity and differ
by sampling error at any finite `nsim`. The sample size is that of
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
and 10000 for
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md).
The routine
[`mv_expected_sampled()`](https://statmodels7.github.io/distributions7/reference/mv_expected_sampled.md)
does the sampling and also serves an explicit `method` for a
multivariate family that has a closed form.

## Monte Carlo error

The result carries Monte Carlo error of order `nsim^(-1/2)`, and two
calls under different seeds give different values. A seed set before the
call makes the result reproducible, and a fit that inverts this matrix
moves with it.

## See also

[`distrib_expected_hessian.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MvGaussianDistrib.md)
and
[`distrib_expected_hessian.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MvStudentTDistrib.md)
for two closed forms,
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md),
which reports whether a family has a closed form, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
# Every shipped family has its own method, so the fallback is called
# directly.
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)
base <- S7::method(distrib_expected_hessian, multivariate_distrib)

set.seed(1)
sampled <- base(d, matrix(0, 3, 2), distributions7:::align_theta(d, theta),
                method = expect_mc(nsim = 4000))
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

# Two seeds give two values, because the result is a sample.
set.seed(2)
again <- base(d, matrix(0, 3, 2), distributions7:::align_theta(d, theta),
              method = expect_mc(nsim = 4000))
c(first = sampled$mu1_mu1[1], second = again$mu1_mu1[1])
#>  first second 
#>  -1.25  -1.25 
```
