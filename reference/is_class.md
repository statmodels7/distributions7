# Whether an S7 Class Is the Given Base Class

Compares two S7 class objects, treating a class re-created from the same
definition as the same class.

## Usage

``` r
is_class(cls, base)
```

## Arguments

- cls:

  The S7 class recorded on a method.

- base:

  The base class to compare against.

## Value

A single logical.

## Details

Several functions here decide whether a method comes from the base class
or from a subclass that registered its own, by reading
`attr(m, "signature")[[1]]`, the documented S7 route. The comparison
that follows must not use
[`identical()`](https://rdrr.io/r/base/identical.html): on S7 class
objects that tests object identity, so it returns `FALSE` for a class
re-created from the same definition. A class is re-created whenever the
code of the package is re-evaluated instead of being loaded, as under
coverage instrumentation.

When [`identical()`](https://rdrr.io/r/base/identical.html) fails in
this way the comparison is silent, and a numerical fallback is taken for
an analytic method. The function therefore returns `TRUE` for identical
objects as a fast path, and otherwise compares the class name together
with the package name.

## See also

[`has_analytic_quantile()`](https://statmodels7.github.io/distributions7/reference/has_analytic_quantile.md),
[`has_exact_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/has_exact_cdf_deriv.md)
