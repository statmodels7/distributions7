# Skewness of a Multivariate Distribution

Signals an error. A scalar skewness of a vector response is a choice
among several quantities: Mardia's, Malkovich and Afifi's, and the
vector of the coordinatewise marginal skewnesses are different numbers,
and the method does not return one of them under the bare name.

## Arguments

- x:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list of parameters. It is not read, because the error is
  signalled first.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

Does not return: it signals an error naming the family.

## Details

[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
gives a coordinatewise skewness where the marginal is a univariate
object. A Dirichlet's marginal is a beta and has a skewness; an
elliptical family's marginal is a one-dimensional object of its own
class and signals the same error. Applying the generic to a marginal
specifies which quantity is computed.

## See also

[`kurtosis.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.multivariate_distrib.md)
for the fourth moment,
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
for the coordinatewise route, and
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)
try(skewness(d, theta))
#> Error : skewness() is not defined for 'multivariate gaussian [2d, sigma=log_cholesky]': a vector response has no single skewness: Mardia's, Malkovich and Afifi's and the vector of coordinatewise marginal skewnesses are different quantities. Use a univariate family instead.

# A Dirichlet's marginal is a univariate beta and has a skewness. An
# elliptical family's marginal is a one-dimensional object of its own class
# and signals the same error.
m <- mv_marginal(dirichlet_distrib(3),
                 list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8), 1)
skewness(m$distrib, m$theta)
#> [1] 0.1795466
try(skewness(mv_marginal(d, theta, 1)$distrib,
             mv_marginal(d, theta, 1)$theta))
#> Error : skewness() is not defined for 'multivariate gaussian [1d, sigma=log_cholesky]': a vector response has no single skewness: Mardia's, Malkovich and Afifi's and the vector of coordinatewise marginal skewnesses are different quantities. Use a univariate family instead.
```
