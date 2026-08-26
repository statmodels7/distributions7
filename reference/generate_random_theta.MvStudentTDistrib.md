# Random Parameters for a Multivariate Student t

Draws a parameter vector for testing: each location uniform on \\(-1,
1)\\, each free value of the matrix parametrization uniform on \\(-0.4,
0.4)\\, and the degrees of freedom uniform on \\(3, 12)\\. The matrix
band puts the scale matrix near the identity, as it does for the
gaussian; the \\\nu\\ band is chosen so that the family is genuinely
heavy-tailed and every moment
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
compares still exists, the variance needing \\\nu \> 2\\ and the
kurtosis \\\nu \> 4\\.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of `distrib@n_params` single numbers, named and ordered as
`distrib@params`.

## See also

[`generate_random_theta.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.MvGaussianDistrib.md)
for the gaussian's bands,
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md),
which draws parameters this way, and
[`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md)
for the generic.

## Examples

``` r
d <- mvstudent_t_distrib(2)

set.seed(11)
unlist(generate_random_theta(d))
#>          mu1          mu2 sigma_log_L1 sigma_log_L2   sigma_L2.1           nu 
#> -0.445500412 -0.998963374  0.008486698 -0.388761673 -0.348248179 11.593643029 

# The three bands, over 400 draws.
set.seed(12)
round(apply(replicate(400, unlist(generate_random_theta(d))), 1, range), 2)
#>      mu1   mu2 sigma_log_L1 sigma_log_L2 sigma_L2.1    nu
#> [1,]  -1 -0.99         -0.4         -0.4       -0.4  3.03
#> [2,]   1  1.00          0.4          0.4        0.4 11.99

# Every draw keeps nu above 4, so the covariance and the kurtosis of every
# coordinate exist.
set.seed(13)
all(replicate(50, is.finite(variance(d, generate_random_theta(d))[1, 1])))
#> [1] TRUE
```
