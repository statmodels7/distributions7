# How Many Observations a Response Holds

Returns the number of observations in `y`: its length for a univariate
distribution, and its number of ROWS for a multivariate one. Every place
in the package that would otherwise write `length(y)` goes through this,
because for a matrix response `length(y)` counts entries, and a
recycling check built on it would ask for parameters of length \\np\\.

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

## Details

A wrong length here is quiet. The multivariate gaussian's expected
information once built a zero vector with `length(y)` to stand for a
parameter that does not vary; that vector came out \\np\\ long, recycled
against the \\p\\-long components, and inflated every diagonal entry of
the information by a factor of \\p\\, so every standard error of a
multivariate fit was \\\sqrt{p}\\ too small. Nothing failed. Anywhere a
matrix response meets code written for a vector,
[`length()`](https://rdrr.io/r/base/length.html) is a defect and this is
the question to ask.

## See also

[`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
and
[`as_mv_matrix()`](https://statmodels7.github.io/distributions7/reference/as_mv_matrix.md),
which puts a response in the shape this counts.

## Examples

``` r
# A univariate response is counted by length.
n_obs(gaussian1_distrib(), c(1.2, -0.4, 0.8))
#> [1] 3

# A multivariate one by rows, where length() would give the entry count.
y <- matrix(0, 5, 2)
c(n_obs = n_obs(mvgaussian_distrib(2), y), length = length(y))
#>  n_obs length 
#>      5     10 
```
