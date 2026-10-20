# Quantile Function of a Multivariate Distribution

Signals an error. A quantile inverts a distribution function along an
ordering of the line, and \\\mathbb{R}^p\\ has no ordering compatible
with its vector structure for \\p \> 1\\, so the multivariate quantile
is not one object: several quantities go by that name in the literature,
and they differ.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- p:

  A numeric vector of probabilities. It is not read, because the error
  is signalled first.

- theta:

  A named list of parameters. It is not read.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

Does not return: it signals an error naming the family.

## See also

[`distrib_cdf.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.multivariate_distrib.md)
for the distribution function,
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md),
which for the Dirichlet and the multinomial returns a univariate
distribution with a quantile function, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)

try(distrib_quantile(d, 0.5, theta))
#> Error : distrib_quantile() is not defined for 'multivariate gaussian [2d, sigma=log_cholesky]': a quantile inverts an ordering of the line, and there is none in several dimensions.

# The marginal of a gaussian signals the same error; the marginal of a
# Dirichlet is a univariate beta and has a quantile function.
b <- mv_marginal(dirichlet_distrib(3),
                 list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8), 1)
distrib_quantile(b$distrib, c(0.025, 0.5, 0.975), b$theta)
#> [1] 0.1311312 0.4195581 0.7556417
```
