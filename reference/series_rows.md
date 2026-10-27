# Batched Series Summation with Rejection

Calls
[`numericals7::series_vec()`](https://statmodels7.github.io/numericals7/reference/series_vec.html),
muffles its warning and promotes a row that did not converge to an error
that names the row. The tolerances and the term limit are those of the
[ExpectIntegrate](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)
object in force (see
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)).

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
