# No Enumerable Support

Signals an error. A continuous family places no mass on any finite set
of points, and a discrete family whose support is infinite has no finite
matrix to return, so there is nothing correct to hand back. A family
that DOES have a finite support registers its own method, and among the
four that ship only the multinomial does.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object with no method of its own.

- theta:

  A named list of parameters. Not examined: the error is raised before
  it is read.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

Never returns: it always signals an error naming the family.

## See also

[`mv_support()`](https://statmodels7.github.io/distributions7/reference/mv_support.md)
for the generic and the multinomial's answer,
[`has_mv_support()`](https://statmodels7.github.io/distributions7/reference/has_mv_support.md)
for the predicate that avoids this error, and
[`mv_reference_draw()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md),
the route a continuous family takes instead.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0)
try(mv_support(d, theta))
#> Error : 'multivariate gaussian [2d, sigma=log_cholesky]' does not enumerate a support. A continuous family has no such set,
#>   and a discrete one whose support is infinite has no finite matrix to
#>   return; a family that does have one registers this generic.

# The predicate a consumer asks first.
c(gaussian = distributions7:::has_mv_support(d),
  multinomial = distributions7:::has_mv_support(
    multinomial_distrib(3, size = 4)))
#>    gaussian multinomial 
#>       FALSE        TRUE 
```
