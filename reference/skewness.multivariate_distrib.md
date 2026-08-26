# No Skewness Without Saying Which One

Signals an error. A scalar skewness for a vector response is not one
quantity but a choice among several: Mardia's, Malkovich-Afifi's, and
the vector of coordinatewise marginal skewnesses are different numbers,
and they do not agree. Returning any one of them under the bare name
would make a choice the caller did not.

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

## Details

[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
is only a way round this where the marginal is a genuinely univariate
object. A Dirichlet's marginal is a beta and answers; an elliptical
family's marginal is a one-dimensional object of its own class, so it
inherits this refusal. Where it does answer, taking the marginal and
asking it is explicit about which quantity is meant, which is the point.

## See also

[`kurtosis.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.multivariate_distrib.md),
refused for the same reason,
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
for the explicit route, and
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
for the generic.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)
try(skewness(d, theta))
#> Error : skewness() is not defined for 'multivariate gaussian [2d, sigma=log_cholesky]': a vector response has no single skewness: Mardia's, Malkovich-Afifi's and the vector of coordinatewise marginal skewnesses are different quantities. Ask a univariate family instead.

# A Dirichlet's marginal is a genuine univariate beta and has one. An
# elliptical family's marginal is a one-dimensional object of its own
# class, so it inherits this refusal rather than answering zero.
m <- mv_marginal(dirichlet_distrib(3),
                 list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8), 1)
skewness(m$distrib, m$theta)
#> [1] 0.1795466
try(skewness(mv_marginal(d, theta, 1)$distrib,
             mv_marginal(d, theta, 1)$theta))
#> Error : skewness() is not defined for 'multivariate gaussian [1d, sigma=log_cholesky]': a vector response has no single skewness: Mardia's, Malkovich-Afifi's and the vector of coordinatewise marginal skewnesses are different quantities. Ask a univariate family instead.
```
