# Excess Kurtosis of a Distribution or Sample

Computes the excess kurtosis (fourth standardized moment minus 3). For
`distrib` objects it is evaluated numerically via
[`moment`](https://statmodels7.github.io/distributions7/reference/moment.md);
for numeric vectors the sample excess kurtosis (population denominator)
is returned.

## Usage

``` r
kurtosis(x, ...)
```

## Arguments

- x:

  An object inheriting from class `"distrib"`, or a numeric vector.

- ...:

  For `distrib` objects: `theta` and further arguments passed to
  [`moment`](https://statmodels7.github.io/distributions7/reference/moment.md).
  For numeric vectors: `na.rm`.

## Value

A numeric vector.

## Examples

``` r
kurtosis(gaussian_distrib(), list(mu = 0, sigma = 1))
#> [1] 2.876117e-10
kurtosis(gamma_distrib(), list(mu = 2, sigma2 = 1))
#> [1] 1.5
```
