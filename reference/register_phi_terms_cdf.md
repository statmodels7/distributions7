# Register the Four CDF Derivative Orders of a Normal-Tail Family

Turns a function returning the terms into the four methods.

## Usage

``` r
register_phi_terms_cdf(cls, term_fn)
```

## Arguments

- cls:

  The S7 class.

- term_fn:

  A function of `(distrib, q, theta)` returning the term list
  [`phi_terms_cdf_deriv_k`](https://statmodels7.github.io/distributions7/reference/phi_terms_cdf_deriv_k.md)
  consumes.

## Value

Invisibly `NULL`; called for the registration.
