# How Many Observations a Response Holds

The number of observations in `y`: its length for a univariate
distribution, and the number of rows for a multivariate one.

## Usage

``` r
n_obs(distrib, y)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  The response.

## Value

A single integer.

## Details

Every place that used to write `length(y)` goes through this instead,
because for a matrix response `length(y)` counts entries rather than
observations, and the recycling checks built on it would ask for
parameters of length \\np\\.

## See also

[`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)

## Examples

``` r
n_obs(gaussian1_distrib(), c(1, 2, 3))
#> [1] 3
n_obs(mvgaussian_distrib(2), matrix(0, 5, 2))
#> [1] 5
```
