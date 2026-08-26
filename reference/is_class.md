# Is an S7 Class the Given Base Class?

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

Several places here ask "did this method come from the base class, or
did the subclass register its own?", and answer it with the documented
S7 trick of reading `attr(m, "signature")[[1]]`. The comparison that
follows must not be
[`identical()`](https://rdrr.io/r/base/identical.html): on S7 class
objects that is object identity, so it returns `FALSE` for a class
re-created from the same definition. That happens whenever the package's
code is re-evaluated in place of being loaded, as it is under coverage
instrumentation.

The failure is silent and can be severe. In linkfunctions7 the same
mistake made a numerical fallback look like an analytic method and sent
the fourth derivative of the log link wrong by a factor of 900, under
coverage only, while every other check stayed green. Identity is kept as
a fast path and name-with-package as the answer.

## See also

[`has_analytic_quantile()`](https://statmodels7.github.io/distributions7/reference/has_analytic_quantile.md),
[`has_exact_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/has_exact_cdf_deriv.md)
