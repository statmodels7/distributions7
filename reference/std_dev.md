# Standard Deviation of a Distribution or Sample

Computes the standard deviation as the square root of
[`variance`](https://statmodels7.github.io/distributions7/reference/variance.md).
For numeric vectors the sample standard deviation
[`sd`](https://rdrr.io/r/stats/sd.html) is returned.

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
