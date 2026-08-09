# Standard Deviation of a Distribution or Sample

Computes the standard deviation as the square root of
[`variance`](https://statmodels7.github.io/distributions7/reference/variance.md).

## Usage

``` r
std_dev(x, ...)
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

## Details

\$\$\operatorname{sd}(Y) = \sqrt{\operatorname{Var}(Y)}.\$\$ For numeric
vectors the sample standard deviation
[`sd`](https://rdrr.io/r/stats/sd.html) is returned.

## See also

[`expectation`](https://statmodels7.github.io/distributions7/reference/expectation.md),
[`moment`](https://statmodels7.github.io/distributions7/reference/moment.md),
[`variance`](https://statmodels7.github.io/distributions7/reference/variance.md),
[`skewness`](https://statmodels7.github.io/distributions7/reference/skewness.md),
[`kurtosis`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)

## Examples

``` r
std_dev(gaussian1_distrib(), list(mu = 0, sigma = 2))
#> [1] 2
```
