# Transpose and Simplify Parameter List Structure

Transposes a list structure (swapping "columns" and "rows") and
simplifies the inner elements into atomic vectors.

Turns a list of `k` equal-length vectors into a list of `n` vectors of
length `k`, one per observation, keeping the names.

## Usage

``` r
transpose_params(theta)
```

## Arguments

- theta:

  A list to be transposed.

## Value

A `list` where each element has been transposed and simplified to an
atomic vector.

## Examples

``` r
transpose_params(list(mu = c(0, 1), sigma = c(1, 2)))
#> [[1]]
#>    mu sigma 
#>     0     1 
#> 
#> [[2]]
#>    mu sigma 
#>     1     2 
#> 
```
