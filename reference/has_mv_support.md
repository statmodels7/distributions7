# Whether a Multivariate Family Enumerates Its Support

`TRUE` when the family registers
[`mv_support()`](https://statmodels7.github.io/distributions7/reference/mv_support.md),
which a discrete multivariate family does and a continuous one cannot.

## Usage

``` r
has_mv_support(x)
```

## Arguments

- x:

  An object inheriting from class
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md).

## Value

A single logical.

## Details

The test reads the registered method, because the multivariate branch
sits beside
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)
and
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md)
and not under either, so there is no class to test. The owning class of
a method is read through
[`is_class()`](https://statmodels7.github.io/distributions7/reference/is_class.md),
never with [`identical()`](https://rdrr.io/r/base/identical.html), which
is object identity and fails whenever the package's code is re-evaluated
in place of being loaded.

## See also

[`mv_support()`](https://statmodels7.github.io/distributions7/reference/mv_support.md),
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
