# Register the Four CDF Derivative Orders of an Exponential Survival Family

Turns a function returning \\L\\ and its partial-derivative evaluator
into the four methods.

## Usage

``` r
register_surv_cdf(cls, pieces)
```

## Arguments

- cls:

  The S7 class.

- pieces:

  A function of `(distrib, q, theta)` returning a list with `Lval` and
  `Lderiv`.

## Value

Invisibly `NULL`; called for the registration.
