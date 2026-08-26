# Kurtosis of the Skew Normal in the Centered Parametrization

Returns the excess kurtosis, computed from
[`kurtosis.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.SkewNormal1Distrib.md)
at the implied direct parameters. It is the one moment this
parametrization does not name: fixing the first three uses up all three
parameters, so the fourth follows from them.

The consequence for use is that a skew normal cannot match an arbitrary
first four moments. At \\\gamma_1 = 0.5\\ the excess kurtosis is 0.347,
and it is determined by the skewness alone.

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

A numeric vector of excess kurtoses, of the length the recycled
parameters imply.

## See also

[`skewness.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.SkewNormal2Distrib.md),
which does name a parameter;
[`kurtosis.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.SkewNormal1Distrib.md)
for the closed form; and
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md),
whose degrees of freedom give the fourth moment a parameter of its own.

## Examples

``` r
d <- skewnormal2_distrib()

# The kurtosis follows from the skewness and from nothing else: changing
# the mean and the standard deviation leaves it alone.
c(a = kurtosis(d, list(mu = 0, sigma = 1, gamma1 = 0.5)),
  b = kurtosis(d, list(mu = 9, sigma = 4, gamma1 = 0.5)))
#>         a         b 
#> 0.3471199 0.3471199 

# It rises with the skewness, and is zero at symmetry.
vapply(c(0.001, 0.3, 0.6, 0.9),
       function(g) kurtosis(d, list(mu = 0, sigma = 1, gamma1 = g)), 0)
#> [1] 8.746873e-05 1.756633e-01 4.426439e-01 7.600512e-01
```
