# Batched Series Summation with Refusal

Calls
[`series_vec`](https://statmodels7.github.io/numericals7/reference/series_vec.html)
with its warning muffled and promotes a row that did not converge to an
error naming it: a series that does not converge is a failure of the
request, not a number.

## Usage

``` r
series_rows(term, from, n_rows)
```

## Arguments

- term:

  A function `term(k, i)` in `series_vec`'s contract.

- from:

  The first summation index.

- n_rows:

  The number of parameter rows.

## Value

A numeric vector of sums, one per row.
