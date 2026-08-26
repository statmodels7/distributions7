# Multinomial Covariance Matrix

Returns the covariance matrix of the counts, \$\$\operatorname{Cov}(Y_i,
Y_j) = n(\delta\_{ij}p_i - p_i p_j),\$\$ so the coordinate variances are
\\np_j(1-p_j)\\, the binomial ones, and every covariance is negative:
the counts share a fixed total, so one category can only grow at the
others' expense.

The matrix is **singular by construction**, the coordinates summing to
\\n\\, so the vector of ones is in its null space and the rank is
\\p-1\\. Its leading \\(p-1)\times(p-1)\\ block is the expected
information on the default chart; see
[`distrib_expected_hessian.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MultinomialDistrib.md).

## Arguments

- distrib:

  A `MultinomialDistrib` object, from
  [`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md).

- theta:

  A named list of the simplex's free values on the parameter scale, each
  of length 1.

## Value

A symmetric \\p \times p\\ numeric matrix of rank \\p-1\\.

## See also

[`mv_location.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_location.MultinomialDistrib.md)
for the mean,
[`variance.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.MultinomialDistrib.md),
which calls this,
[`mv_marginal.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.MultinomialDistrib.md)
for a coordinate's binomial law, and
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
for the generic.

## Examples

``` r
d <- multinomial_distrib(3, size = 5)
th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
S <- mv_sigma(d, th)
round(S, 6)
#>           [,1]      [,2]      [,3]
#> [1,]  1.222629 -0.550386 -0.672243
#> [2,] -0.550386  0.958122 -0.407736
#> [3,] -0.672243 -0.407736  1.079979

# Singular by construction: rank p - 1, with the ones vector in the null
# space.
c(rank = qr(S)$rank, dim = ncol(S))
#> rank  dim 
#>    2    3 
round(S %*% rep(1, 3), 12)
#>      [,1]
#> [1,]    0
#> [2,]    0
#> [3,]    0

# The diagonal is the binomial variance of each category.
pr <- mv_location(d, th) / 5
rbind(diagonal = diag(S), binomial = 5 * pr * (1 - pr))
#>              [,1]      [,2]     [,3]
#> diagonal 1.222629 0.9581222 1.079979
#> binomial 1.222629 0.9581222 1.079979
```
