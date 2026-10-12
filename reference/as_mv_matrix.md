# Coerce a Multivariate Response to a Matrix

Puts `y` in the \\n \times p\\ form that every multivariate method
expects. A matrix of the right width is returned unchanged, and a
numeric vector of length \\p\\ is read as a single observation and
returned as a one-row matrix, which is the reading of a density at one
point. A vector of any other length is an error naming its length and
the dimension, and a response that is not numeric is an error as well.

## Usage

``` r
as_mv_matrix(distrib, y)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object, from which `n_dim` is read.

- y:

  A numeric matrix with `distrib@n_dim` columns, or a numeric vector of
  length `distrib@n_dim`. An empty vector gives a matrix with no rows. A
  vector of any other length, a matrix of the wrong width and a
  non-numeric `y` are errors.

## Value

A numeric matrix with `distrib@n_dim` columns.

## See also

[`n_obs()`](https://statmodels7.github.io/distributions7/reference/n_obs.md)
for the row count and
[`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
for the response convention.

## Examples

``` r
d <- mvgaussian1_distrib(2)

# A vector of length p is one observation.
dim(distributions7:::as_mv_matrix(d, c(1, -1)))
#> [1] 1 2

# A matrix of the right width is returned as it stands.
dim(distributions7:::as_mv_matrix(d, matrix(0, 4, 2)))
#> [1] 4 2

# Any other length is an error naming both numbers.
try(distributions7:::as_mv_matrix(d, c(1, 2, 3)))
#> Error : 'y' is a vector of length 3, which is read as one observation, but
#>   'multivariate gaussian [2d, sigma=log_cholesky]' has dimension 2. Supply a matrix with 2 columns for several
#>   observations.
```
