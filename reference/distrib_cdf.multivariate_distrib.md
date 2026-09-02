# No Distribution Function in Several Dimensions

Signals an error. The distribution function of a multivariate law is
\\P(Y_1 \le q_1, \ldots, Y_p \le q_p)\\, an integral over an orthant.
There is no closed form for it in general, and the one-dimensional
fallback has no counterpart:
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)'s
default integrates the density over an interval of the line, and an
orthant is not one. A numerical orthant probability is a separate piece
of work with its own accuracy question, and this package does not
attempt it.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- q:

  A numeric matrix of quantiles, one row per point. Not examined: the
  error is raised before it is read.

- theta:

  A named list of parameters. Not examined.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

Never returns: it always signals an error naming the family.

## Details

The refusal is on the class, so it covers a family written elsewhere as
well as the four that ship. A family that CAN answer registers its own
method and the refusal never runs.

## See also

[`distrib_quantile.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.multivariate_distrib.md),
refused for a related reason,
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md),
which for some families returns a univariate distribution that does
answer, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)

try(distrib_cdf(d, rbind(c(0, 0)), theta))
#> Error : distrib_cdf() is not defined for 'multivariate gaussian [2d, sigma=log_cholesky]': the distribution function is an integral over an orthant, with no closed form and no one-dimensional fallback.

# A gaussian's marginal is a ONE-DIMENSIONAL MvGaussianDistrib, so it
# inherits the same refusal.
m <- mv_marginal(d, theta, 1)
class(m$distrib)[1]
#> [1] "distributions7::MvGaussian1Distrib"
try(distrib_cdf(m$distrib, 0, m$theta))
#> Error : distrib_cdf() is not defined for 'multivariate gaussian [1d, sigma=log_cholesky]': the distribution function is an integral over an orthant, with no closed form and no one-dimensional fallback.

# A Dirichlet's marginal is a genuine univariate beta, and answers.
b <- mv_marginal(dirichlet_distrib(3),
                 list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8), 1)
class(b$distrib)[1]
#> [1] "distributions7::Beta1Distrib"
distrib_cdf(b$distrib, 0.4, b$theta)
#> [1] 0.4565102
```
