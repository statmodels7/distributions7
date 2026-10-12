# Batched Series Summation with Rejection

Calls
[`numericals7::series_vec()`](https://statmodels7.github.io/numericals7/reference/series_vec.html)
with its warning muffled and promotes a row that did not converge to an
error that names the row.

## Usage

``` r
series_rows(term, from, n_rows)
```

## Arguments

- term:

  A function `term(k, i)` in the convention of `series_vec`.

- from:

  The first summation index.

- n_rows:

  The number of parameter rows.

## Value

A numeric vector of sums, one per row.
