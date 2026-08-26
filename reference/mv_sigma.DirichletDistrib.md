# Dirichlet Covariance Matrix

Returns the covariance matrix of the family, \$\$\operatorname{Cov}(Y_i,
Y_j) = \dfrac{\delta\_{ij}\mu_i - \mu_i\mu_j}{\phi + 1},\$\$ so the
coordinate variances are \\\mu_j(1-\mu_j)/(\phi+1)\\ and every
covariance is negative, a rise in one coordinate having to be paid for
by the others.

The matrix is **singular by construction**: the coordinates sum to one,
so the vector of ones is in its null space and the rank is \\p-1\\.
Anything that inverts a covariance must use the marginals or the free
vector instead. The concentration acts as a precision, every entry
falling as \\1/(\phi+1)\\.

## Arguments

- distrib:

  A `DirichletDistrib` object, from
  [`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md).

- theta:

  A named list of parameters on the parameter scale: the mean's free
  values followed by `phi`, each of length 1.

## Value

A symmetric \\p \times p\\ numeric matrix of rank \\p-1\\.

## See also

[`mv_location.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_location.DirichletDistrib.md)
for the mean,
[`variance.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.DirichletDistrib.md),
which calls this,
[`mv_marginal.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.DirichletDistrib.md)
for a coordinate's law, and
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
for the generic.

## Examples

``` r
d <- dirichlet_distrib(3)
th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
S <- mv_sigma(d, th)
round(S, 6)
#>           [,1]      [,2]      [,3]
#> [1,]  0.018810 -0.008467 -0.010342
#> [2,] -0.008467  0.014740 -0.006273
#> [3,] -0.010342 -0.006273  0.016615

# Singular by construction: the rank is p - 1 and the ones vector is in
# the null space.
c(rank = qr(S)$rank, dim = ncol(S))
#> rank  dim 
#>    2    3 
round(S %*% rep(1, 3), 12)
#>      [,1]
#> [1,]    0
#> [2,]    0
#> [3,]    0

# The concentration is a precision: every entry falls as 1 / (phi + 1).
vapply(c(2, 12, 100), function(p)
  mv_sigma(d, list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = p))[1, 1],
  numeric(1))
#> [1] 0.081508617 0.018809681 0.002421048
```
