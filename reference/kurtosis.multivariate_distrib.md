# No Kurtosis Without Saying Which One

Signals an error, for the reason
[`skewness.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.multivariate_distrib.md)
gives: a scalar kurtosis for a vector response is a choice among
Mardia's, Malkovich-Afifi's and the vector of coordinatewise marginal
kurtoses, which are different numbers. A caller who wants a
coordinatewise one takes the marginal and asks it, where the marginal is
a genuinely univariate object: a Dirichlet's is a beta and answers, an
elliptical family's is a one-dimensional object of its own class and
inherits this refusal.

## Arguments

- x:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list of parameters. Not examined: the error is raised before
  it is read.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

Never returns: it always signals an error naming the family.

## See also

[`skewness.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.multivariate_distrib.md)
for the same refusal one moment down,
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
for the explicit route, and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
for the generic.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5, nu = 6)
try(kurtosis(d, theta))
#> Error : kurtosis() is not defined for 'multivariate student t [2d, sigma=log_cholesky]': a vector response has no single kurtosis: Mardia's, Malkovich-Afifi's and the vector of coordinatewise marginal kurtosises are different quantities. Ask a univariate family instead.

# A Student t's marginal is a one-dimensional MvStudentTDistrib, so it
# inherits the same refusal.
try(kurtosis(mv_marginal(d, theta, 1)$distrib,
             mv_marginal(d, theta, 1)$theta))
#> Error : kurtosis() is not defined for 'multivariate student t [1d, sigma=log_cholesky]': a vector response has no single kurtosis: Mardia's, Malkovich-Afifi's and the vector of coordinatewise marginal kurtosises are different quantities. Ask a univariate family instead.

# A Dirichlet's marginal is a genuine univariate beta and answers.
b <- mv_marginal(dirichlet_distrib(3),
                 list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8), 1)
kurtosis(b$distrib, b$theta)
#> [1] -0.501495
```
