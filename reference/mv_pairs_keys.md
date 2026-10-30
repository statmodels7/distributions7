# The Elements of a Panel Matrix and Their Colours

Returns the names of the elements that
[`mv_pairs_panels()`](https://statmodels7.github.io/distributions7/reference/mv_pairs_panels.md)
draws, as the names of a vector of colours: the density (its marginal
curves and its contours) and, when there are data, the kernel density
estimate of each marginal and the observations.

## Usage

``` r
mv_pairs_keys(data)
```

## Arguments

- data:

  The data matrix, or `NULL`.

## Value

A named character vector of colours, of length 1 without data and 3 with
data.
