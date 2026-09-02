# No Location Without a Family That Has One

Signals an error. Not every multivariate family has a location: a
Dirichlet is described by concentrations and a Wishart by a scale and a
count, and handing back the first \\p\\ parameters under the name of a
mean would be a wrong answer in the shape of a right one. A family that
HAS one registers its own method, and the two elliptical families
register
[`mv_leading_location()`](https://statmodels7.github.io/distributions7/reference/mv_leading_location.md).

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object with no method of its own.

- theta:

  A named list of parameters. Not examined: the error is raised before
  it is read.

## Value

Never returns: it always signals an error naming the family.

## See also

[`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
for the generic,
[`mv_leading_location()`](https://statmodels7.github.io/distributions7/reference/mv_leading_location.md)
for the implementation the elliptical families use, and
[`mv_location.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_location.MvGaussianDistrib.md)
for one that answers.

## Examples

``` r
# The Dirichlet registers no location.
try(mv_location(dirichlet_distrib(3),
                list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)))
#> [1] 0.4260125 0.2583897 0.3155978

# The gaussian does.
mv_location(mvgaussian1_distrib(2),
            list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0,
                 sigma_log_L2 = 0, sigma_L2.1 = 0))
#> v1 v2 
#>  1 -1 
```
