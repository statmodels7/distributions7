# Register the Third and Fourth Response Derivatives of a Family

Turns a function of `(distrib, y, theta, k)` into the two methods, so
that a family writes its rule once instead of twice.

## Usage

``` r
register_dy_k(cls, f)
```

## Arguments

- cls:

  The S7 class.

- f:

  The rule.

## Value

Invisibly `NULL`; called for the registration.
