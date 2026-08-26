# No Quantile Function in Several Dimensions

Signals an error. A quantile is defined by inverting a distribution
function along an ORDERING, and \\\mathbb{R}^p\\ has none for \\p \>
1\\. The obstruction is deeper than the one that stops
[`distrib_cdf.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.multivariate_distrib.md):
there a number exists and is hard to compute, here there is no number to
compute. Several quantities go by the name in the literature, and they
disagree.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- p:

  A numeric vector of probabilities. Not examined: the error is raised
  before it is read.

- theta:

  A named list of parameters. Not examined.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

Never returns: it always signals an error naming the family.

## See also

[`distrib_cdf.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.multivariate_distrib.md)
for the related refusal,
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md),
which for some families returns a univariate distribution that does
answer, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)

try(distrib_quantile(d, 0.5, theta))
#> Error : distrib_quantile() is not defined for 'multivariate gaussian [2d, sigma=log_cholesky]': a quantile inverts an ordering of the line, and there is none in several dimensions.

# A gaussian's marginal is a one-dimensional MvGaussianDistrib and refuses
# too; a Dirichlet's is a univariate beta and answers.
b <- mv_marginal(dirichlet_distrib(3),
                 list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8), 1)
distrib_quantile(b$distrib, c(0.025, 0.5, 0.975), b$theta)
#> [1] 0.1311312 0.4195581 0.7556417
```
