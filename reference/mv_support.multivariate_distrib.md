# Support of a Family Without an Enumerable One

Signals an error. A continuous family places no mass on a finite set of
points, and a discrete family whose support is infinite has no finite
matrix to return. A family with a finite support registers its own
method; among the shipped families the multinomial does.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object with no method of its own.

- theta:

  A named list of parameters. Not read: the error is signalled first.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

Does not return: it signals an error naming the family.

## See also

[`mv_support()`](https://statmodels7.github.io/distributions7/reference/mv_support.md)
for the generic and the multinomial's method,
[`has_mv_support()`](https://statmodels7.github.io/distributions7/reference/has_mv_support.md)
for the test that avoids this error, and
[`mv_reference_draw()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md),
the route that
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
takes for a continuous family.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0)
try(mv_support(d, theta))
#> Error : 'multivariate gaussian [2d, sigma=log_cholesky]' does not enumerate a support. A continuous family has no such set,
#>   and a discrete one whose support is infinite has no finite matrix to
#>   return; a family that does have one registers this generic.

# The test that check_distrib() applies first.
c(gaussian = distributions7:::has_mv_support(d),
  multinomial = distributions7:::has_mv_support(
    multinomial_distrib(3, size = 4)))
#>    gaussian multinomial 
#>       FALSE        TRUE 
```
