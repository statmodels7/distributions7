# Excess Kurtosis of the Gaussian Distribution

Exactly zero at every parameter value. The raw fourth standardized
moment of a Gaussian is 3, and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
reports the excess over that number, so this family sits at the origin
of the scale by construction. Every other family's excess kurtosis is
read as a comparison with it.

## Arguments

- x:

  A `Gaussian1Distrib`, from
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

- theta:

  A named list with components `mu` and `sigma` (positive), each a
  numeric vector of length 1 or `n`. The values are not read, only their
  lengths.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of zeros, of length
`max(length(theta$mu), length(theta$sigma))`.

## See also

[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
for the convention;
[`kurtosis.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.LaplaceDistrib.md),
which is 3, and
[`kurtosis.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.LogisticDistrib.md),
which is 6/5;
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

## Examples

``` r
d <- gaussian1_distrib()

# Zero for every Gaussian, which is what the excess convention means.
kurtosis(d, list(mu = c(-2, 0, 5), sigma = c(1, 2, 3)))
#> [1] 0 0 0
```
