# Dirichlet Mean Vector

Returns the point of the simplex the mean parameter carries at the free
values in `theta`. For this family that point **is** the mean,
\\\mathbb{E}\[Y_j\] = \mu_j\\, so
[`mean.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.DirichletDistrib.md)
delegates here. The value comes from
[`parameters7::param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.html)
on the object's own simplex chart, so it lies on the open simplex by
construction whatever the free values are.

## Arguments

- distrib:

  A `DirichletDistrib` object, from
  [`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md).

- theta:

  A named list of parameters on the parameter scale: the mean's free
  values followed by `phi`, each of length 1. The concentration is read
  but not used here.

## Value

A numeric vector of length \\p\\, strictly positive and summing to one.

## See also

[`mv_sigma.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.DirichletDistrib.md)
for the covariance,
[`mean.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.DirichletDistrib.md),
which calls this, and
[`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
for the generic.

## Examples

``` r
d <- dirichlet_distrib(3)
mu <- mv_location(d, list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12))
mu
#> [1] 0.4260125 0.2583897 0.3155978
sum(mu)
#> [1] 1

# Whatever the free values, the result is on the simplex.
m2 <- mv_location(d, list(mean_alr1 = -8, mean_alr2 = 5, phi = 12))
c(min = min(m2), sum = sum(m2))
#>          min          sum 
#> 2.245196e-06 1.000000e+00 
```
