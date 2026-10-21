# Whether a Family Supplies Its Own Fourth Derivative

Reports whether the family supplies its own fourth-order derivatives,
rather than reaching one of the base classes' numerical fallbacks.

## Usage

``` r
has_exact_deriv4(x, ...)
```

## Arguments

- x:

  A distribution object.

- ...:

  Unused.

## Value

A single logical.

## Details

The numerical fifth order is one central difference of the fourth, so
its accuracy depends on the quantity that it differences. Where the
fourth order is the family's own, the result is an ordinary numerical
derivative of a quantity that the family computes directly. Where the
fourth order is the numerical fallback of the package, the fifth order
is a difference of a difference, and it is unreliable whatever the
tolerance.

A family in the second case still receives a value, since a distribution
needs only
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md),
and the value comes back without a warning. Its error can reach the
first significant digit, so it is not suitable for building on.
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
emits no order-5 row for such a family, because the check does not apply
to it.

The default method reads the class that owns the registered
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
method, the documented S7 route. A wrapper overrides it and returns the
answer of its parent instead: the fourth derivative of a wrapper is a
partition sum over the first four derivatives of the parent, so it is
the wrapper's own exactly when those of the parent are.

The predicate is one logical for the whole family. A family may own its
fourth-order method and still build part of it from single stencils on
analytic quantities; it returns `TRUE`, and the order-5 row of
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
is then emitted for it. That row appears only when `5` is among the
requested `orders`, which default to `1:4`.

## See also

[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md),
the consumer;
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md),
the predicate on which this is modeled.

## Examples

``` r
has_exact_deriv4(gaussian1_distrib())
#> [1] TRUE
has_exact_deriv4(truncated(gaussian1_distrib(), lower = 0))
#> [1] TRUE
```
