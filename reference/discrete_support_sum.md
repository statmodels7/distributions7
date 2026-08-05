# Summation over an Integer Support

Sums `term(k, i)` over the integers of `[from, to]` for each of `n_rows`
parameter rows at once. A finite support is summed directly in one
matrix evaluation; a support unbounded above goes through
[`series_vec`](https://statmodels7.github.io/numericals7/reference/series_vec.html);
one unbounded below is reflected; one unbounded on both sides is folded
around zero. A row whose series does not converge raises an error naming
it.

## Usage

``` r
discrete_support_sum(term, from, to, n_rows)
```

## Arguments

- term:

  A function `term(k, i)` of two equal-length integer vectors, returning
  the terms elementwise.

- from, to:

  The endpoints of the support, either possibly infinite.

- n_rows:

  The number of parameter rows.

## Value

A numeric vector of sums, one per row.
