# Batched Quadrature with Rejection

Calls
[`quad_vec`](https://statmodels7.github.io/numericals7/reference/quad_vec.html)
with its warning muffled and returns the per-row integrals, `NA` where
the accuracy was not reached; the caller decides what a failed row means
and names it.

## Usage

``` r
quad_rows(integrand, lower, upper)
```

## Arguments

- integrand:

  The integrand, in `quad_vec`'s matrix contract.

- lower, upper:

  Numeric vectors of panel endpoints.

## Value

A numeric vector of integrals, `NA` for failed rows.
