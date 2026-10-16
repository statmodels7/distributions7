# Distribution Function of a Multivariate Distribution

Signals an error. The distribution function of a multivariate law is
\\P(Y_1 \le q_1, \ldots, Y_p \le q_p)\\, an integral over an orthant. It
has no closed form in general, and the default of
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md),
which integrates the density over an interval of the line, has no
counterpart in several dimensions. A numerical orthant probability is a
separate computation with an accuracy question of its own, and the
package does not provide it.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- q:

  A numeric matrix of quantiles, one row per point. Not read: the error
  is signalled first.

- theta:

  A named list of parameters. It is not read.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

Does not return: it signals an error naming the family.

## Details

The method is registered on the class, so it applies to a family written
elsewhere as well as to the shipped ones. A family that can compute the
distribution function registers its own method.

## See also

[`distrib_quantile.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.multivariate_distrib.md)
for the quantile function,
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md),
which for the Dirichlet and the multinomial returns a univariate
distribution with a distribution function, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)

try(distrib_cdf(d, rbind(c(0, 0)), theta))
#> Error : distrib_cdf() is not defined for 'multivariate gaussian [2d, sigma=log_cholesky]': the distribution function is an integral over an orthant, with no closed form and no one-dimensional fallback.

# The marginal of a gaussian is a one-dimensional MvGaussian1Distrib, which
# signals the same error.
m <- mv_marginal(d, theta, 1)
class(m$distrib)[1]
#> [1] "distributions7::MvGaussian1Distrib"
try(distrib_cdf(m$distrib, 0, m$theta))
#> Error : distrib_cdf() is not defined for 'multivariate gaussian [1d, sigma=log_cholesky]': the distribution function is an integral over an orthant, with no closed form and no one-dimensional fallback.

# The marginal of a Dirichlet is a univariate beta, which has a distribution
# function.
b <- mv_marginal(dirichlet_distrib(3),
                 list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8), 1)
class(b$distrib)[1]
#> [1] "distributions7::Beta1Distrib"
distrib_cdf(b$distrib, 0.4, b$theta)
#> [1] 0.4565102
```
