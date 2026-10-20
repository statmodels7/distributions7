# Kurtosis of a Multivariate Distribution

Signals an error, for the same reason as
[`skewness.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.multivariate_distrib.md):
a scalar kurtosis of a vector response is a choice among Mardia's,
Malkovich and Afifi's and the vector of the coordinatewise marginal
kurtoses, which are different numbers. A coordinatewise kurtosis is the
kurtosis of a marginal, where the marginal is a univariate object: a
Dirichlet's is a beta and has one, while an elliptical family's is a
one-dimensional object of its own class and signals the same error.

## Arguments

- x:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list of parameters. Not read: the error is signalled first.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

Does not return: it signals an error naming the family.

## See also

[`skewness.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.multivariate_distrib.md)
for the third moment,
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
for the coordinatewise route, and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
for the generic.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5, nu = 6)
try(kurtosis(d, theta))
#> Error : kurtosis() is not defined for 'multivariate student t [2d, sigma=log_cholesky]': a vector response has no single kurtosis: Mardia's, Malkovich and Afifi's and the vector of coordinatewise marginal kurtoses are different quantities. Use a univariate family instead.

# A Student t's marginal is a one-dimensional MvStudentT1Distrib, which
# signals the same error.
try(kurtosis(mv_marginal(d, theta, 1)$distrib,
             mv_marginal(d, theta, 1)$theta))
#> Error : kurtosis() is not defined for 'multivariate student t [1d, sigma=log_cholesky]': a vector response has no single kurtosis: Mardia's, Malkovich and Afifi's and the vector of coordinatewise marginal kurtoses are different quantities. Use a univariate family instead.

# A Dirichlet's marginal is a univariate beta and has a kurtosis.
b <- mv_marginal(dirichlet_distrib(3),
                 list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8), 1)
kurtosis(b$distrib, b$theta)
#> [1] -0.501495
```
