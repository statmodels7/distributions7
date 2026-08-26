# Dirichlet Random Generation

Draws `n` independent Dirichlet vectors by the representation the family
is defined by: independent gamma variates with shapes \\\alpha_j =
\phi\mu_j\\ and a common rate, divided by their sum. The normalization
removes the rate, so
[`stats::rgamma()`](https://rdrr.io/r/stats/GammaDist.html)'s default of
1 is used, and every row of the result sums to one by construction. The
draws depend on `.Random.seed` in the usual way and consume \\p\\ of R's
streams per observation.

## Arguments

- distrib:

  A `DirichletDistrib` object, from
  [`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list of parameters on the parameter scale: the mean's free
  values followed by `phi`, each of length 1.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `n` rows and \\p\\ columns, each row strictly
positive and summing to one.

## See also

[`distrib_pdf.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.DirichletDistrib.md)
for the density,
[`mv_reference_draw.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.DirichletDistrib.md)
for the uniform proposal
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
integrates against,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- dirichlet_distrib(3)
th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)

# Every row is on the simplex.
set.seed(1)
Y <- distrib_rng(d, 4, th)
round(Y, 4)
#>        [,1]   [,2]   [,3]
#> [1,] 0.2777 0.5360 0.1863
#> [2,] 0.5536 0.2528 0.1936
#> [3,] 0.5921 0.1629 0.2450
#> [4,] 0.6241 0.1656 0.2103
rowSums(Y)
#> [1] 1 1 1 1

# The sample recovers the mean and the coordinate variances, which fall as
# 1 / (phi + 1).
set.seed(3)
Z <- distrib_rng(d, 3e5, th)
rbind(sample = c(mean(Z[, 1]), var(Z[, 1])),
      theoretical = c(mv_location(d, th)[1], mv_sigma(d, th)[1, 1]))
#>                  [,1]       [,2]
#> sample      0.4260478 0.01883675
#> theoretical 0.4260125 0.01880968
```
