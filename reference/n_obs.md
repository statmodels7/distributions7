# How Many Observations a Response Holds

Returns the number of observations in `y`: its length for a univariate
distribution, and its number of rows for a multivariate one. For a
multivariate distribution a vector without a `dim` attribute counts as
one observation, or none when it is empty, whatever its length. The
functions of the package that count observations call this one, because
for a matrix response `length(y)` counts the \\np\\ entries, and a
recycling check built on it would require parameters of length \\np\\.

## Usage

``` r
n_obs(distrib, y)
```

## Arguments

- distrib:

  An object inheriting from `distrib`. Only its class is used, to decide
  which reading of `y` applies.

- y:

  The response: a numeric vector for a univariate distribution, or an
  \\n \times p\\ numeric matrix for a multivariate one.

## Value

A single integer.

## See also

[`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
and
[`as_mv_matrix()`](https://statmodels7.github.io/distributions7/reference/as_mv_matrix.md),
which puts a response in the shape that this function counts.

## Examples

``` r
# A univariate response is counted by length.
n_obs(gaussian1_distrib(), c(1.2, -0.4, 0.8))
#> [1] 3

# A multivariate one by rows, where length() gives the number of entries.
y <- matrix(0, 5, 2)
c(n_obs = n_obs(mvgaussian1_distrib(2), y), length = length(y))
#>  n_obs length 
#>      5     10 
```
