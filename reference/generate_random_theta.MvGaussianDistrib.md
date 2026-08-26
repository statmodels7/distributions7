# Random Parameters for a Multivariate Gaussian

Draws a parameter vector for testing: each mean uniform on \\(-1, 1)\\
and each free value of the matrix parametrization uniform on \\(-0.4,
0.4)\\, which puts the matrix near the identity. The base class would
draw every free value from one wide range; on a log-Cholesky diagonal
that spans four orders of magnitude in the resulting variances, and a
covariance drawn that way is a starting point no fit recovers from. The
narrow band keeps
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
reproducible on this family.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of `distrib@n_params` single numbers, named and ordered as
`distrib@params`.

## See also

[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md),
which draws parameters this way, and
[`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md)
for the generic.

## Examples

``` r
d <- mvgaussian_distrib(2)

set.seed(11)
unlist(generate_random_theta(d))
#>          mu1          mu2 sigma_log_L1 sigma_log_L2   sigma_L2.1 
#> -0.445500412 -0.998963374  0.008486698 -0.388761673 -0.348248179 

# The matrix components stay inside (-0.4, 0.4) and the means inside
# (-1, 1), so the covariance drawn is never far from the identity.
set.seed(12)
many <- replicate(500, unlist(generate_random_theta(d)))
round(apply(many, 1, range), 2)
#>      mu1   mu2 sigma_log_L1 sigma_log_L2 sigma_L2.1
#> [1,]  -1 -1.00         -0.4         -0.4      -0.39
#> [2,]   1  0.98          0.4          0.4       0.40

# Every draw gives a positive definite covariance, the parametrization
# having no boundary to reach.
set.seed(13)
all(replicate(50, min(eigen(mv_sigma(d, generate_random_theta(d)),
                            only.values = TRUE)$values) > 0))
#> [1] TRUE
```
