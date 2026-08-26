# Skewness of the Skew Normal in the Centered Parametrization

Returns \\\gamma_1\\, the third parameter, which in this parametrization
is the standardized third central moment by construction.
[`moment_const()`](https://statmodels7.github.io/distributions7/reference/moment_const.md)
recycles the result to the length the three parameters imply.

The value cannot leave \\(-0.9952717, 0.9952717)\\, the constructor
having bounded the parameter at the supremum the family reaches; see
[`sn_max_skew()`](https://statmodels7.github.io/distributions7/reference/sn_max_skew.md).

## Arguments

- x:

  A `SkewNormal2Distrib` object, from
  [`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md).

- theta:

  A named list with components `mu`, `sigma` and `gamma1`, in any order;
  it is aligned here.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of skewnesses, of the length the recycled parameters
imply.

## See also

[`sn_max_skew()`](https://statmodels7.github.io/distributions7/reference/sn_max_skew.md)
for the bound,
[`skewness.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.SkewNormal1Distrib.md)
for the same quantity computed from a shape, and
[`kurtosis.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.SkewNormal2Distrib.md),
which is not free.

## Examples

``` r
d <- skewnormal2_distrib()
skewness(d, list(mu = 3, sigma = 2, gamma1 = 0.6))
#> [1] 0.6

# It agrees with a quadrature of the standardized third central moment.
th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
c(parameter = skewness(d, th),
  quadrature = integrate(function(v) v^3 * distrib_pdf(d, v, th),
                         -Inf, Inf)$value)
#>  parameter quadrature 
#>        0.5        0.5 
```
