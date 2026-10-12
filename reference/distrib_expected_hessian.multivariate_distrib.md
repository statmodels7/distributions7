# Expected Information of a Multivariate Distribution, by Sampling

The fallback for a multivariate family without a closed form: an average
over a sample drawn from the family. Every shipped multivariate class
registers its own closed form, so this method runs only for a family
written elsewhere.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object without a closed form of its own.

- y:

  An \\n \times p\\ numeric matrix of observations; only its row count
  is used, the expectation being over the law.

- theta:

  A named list of parameters, each component a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`. The generic applies
  the link scale to the result of the method.

- approx:

  One of `"opg"` (the default), `"bartlett"` or `"mc"`; the first two
  are the same computation here. `"integrate"` signals an error.

- nsim:

  The Monte Carlo sample size. Defaults to `10000`. The error falls as
  its square root.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each vector constant.

## The routes

The `"integrate"` route of
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)
splits a one-dimensional integral at quantiles, and a multivariate
family has no quantiles, so `approx = "integrate"` signals an error. The
other values sample: `"opg"` and `"bartlett"` are the same computation
here, the average of the negative outer product of the score, and `"mc"`
averages the observed Hessian. The two averages agree in expectation by
the second Bartlett identity and differ by sampling error at any finite
`nsim`.

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

# Two seeds give two values, the result being a sample.
set.seed(2)
again <- base(d, matrix(0, 3, 2), distributions7:::align_theta(d, theta),
              approx = "mc", nsim = 4000)
c(first = sampled$mu1_mu1[1], second = again$mu1_mu1[1])
#>  first second 
#>  -1.25  -1.25 
```
