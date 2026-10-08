# Does a Family Register a Method of Its Own for a Generic?

`TRUE` where the method of `generic` that dispatches on `x` is
registered on a class other than the four base classes, whose methods
are the numerical fallbacks.

## Usage

``` r
owns_method(x, generic)
```

## Arguments

- x:

  A distribution object.

- generic:

  An S7 generic dispatching on the distribution.

## Value

A single logical.

## See also

[`has_exact_deriv4()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv4.md),
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
