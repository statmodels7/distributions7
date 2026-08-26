# Skewness of the Geometric Distribution

Closed form: \\\gamma_1 = (1+2\mu)/\sqrt{\mu(1+\mu)}\\. It is positive
at every mean and, unlike a Poisson's, it does not vanish as the counts
grow: it tends to 2 from above, the exponential's value, which is the
continuous limit of the family.

## Arguments

- x:

  A `GeometricDistrib`, from
  [`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).

- theta:

  A named list with one component, `mu` (positive), a numeric vector of
  length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, the length of `theta$mu`, above 2 throughout.

## See also

[`kurtosis.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.GeometricDistrib.md),
which tends to 6;
[`skewness.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.ExponentialDistrib.md),
the continuous limit;
[`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).

## Examples

``` r
d <- geometric_distrib()

# The published form, written out.
all.equal(skewness(d, list(mu = 3)), 7 / sqrt(12))
#> [1] TRUE

# It tends to the exponential's 2 as the mean grows.
skewness(d, list(mu = c(1, 10, 1000)))
#> [1] 2.121320 2.002271 2.000000
```
