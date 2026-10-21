# Batched Quadrature with Rejection

Calls
[`numericals7::quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.html),
muffles its warning and returns the per-row integrals, `NA` where the
accuracy was not reached. The calling function handles a failed row and
reports it.

## Usage

``` r
quad_rows(integrand, lower, upper)
```

## Arguments

- integrand:

  The integrand, in the matrix convention of `quad_vec`.

- lower, upper:

  Numeric vectors of panel endpoints.

## Value

A numeric vector of integrals, `NA` for failed rows.
