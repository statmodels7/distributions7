# The Order from Which the Fifth Derivative of the Information Is Differenced

Returns the highest order, 4, 3 or 2, at which the family registers its
own derivative of the expected information, and `NA` where it registers
none of them. For a wrapper (a family with a `parent_distrib` property)
the result is at most the one of the parent, since the wrapper's
derivatives of each order are built from the parent's of the same and
lower orders.

## Usage

``` r
d5expected_base(x)
```

## Arguments

- x:

  A distribution object.

## Value

An integer, 2 to 4, or `NA_integer_`.

## See also

[`numerical_d5expected()`](https://statmodels7.github.io/distributions7/reference/numerical_d5expected.md),
which differences that order.
