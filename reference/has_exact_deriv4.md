# Is a Family's Fourth Derivative Its Own

Whether the family supplies its own fourth-order derivatives, rather
than reaching one of the base classes' numerical fallbacks.

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

The fifth order is one central difference of the fourth, so what it is
worth depends entirely on what it differences. Where the fourth is the
family's own the result is an ordinary numerical derivative of a
quantity the family computes directly; where the fourth is the package's
own fallback the fifth is a difference of a difference, which this
package forbids everywhere else and which no tolerance can rescue.

A family reached that way still gets an answer – a distribution needs
only
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md),
and that promise is not withdrawn at the fifth order – but
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
emits no order-5 row for it, a check that does not apply being absent
rather than reported with a status every consumer would misread.

The default reads the class that owns the registered
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
method, the documented S7 route. A wrapper overrides it and asks its
parent instead: a wrapper's fourth derivative is a partition sum over
the parent's first four, so it is its own exactly when the parent's are.

A family may own its fourth-order method and still build part of it from
single stencils on analytic quantities, and the predicate is one logical
for the whole family.
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
is that case: fifteen of its thirty-five fourth-order components are
closed form and the twenty carrying `nu` are not, and it answers `TRUE`,
a model needing a fourth derivative of that family getting one. What
that costs is measured rather than asserted.
`check_distrib(skewt_distrib(), orders = 1:5)` emits an order-5 row
whose verdict depends on `nu`,
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
reading an absolute error and the stencil noise falling as `nu` grows:
swept over `nu` from 3 to 50 the row fails at `nu = 3` (3.7e-03 to
4.6e-03 against a tolerance of 1e-3), is marginal at 5, and passes from
8 upward with orders of margin. Read per component against its own scale
the fifth order of that family is untrustworthy at every `nu`, which is
why `orders` defaults to `1:4`.

## See also

[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md),
the consumer;
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md),
the predicate this is modelled on.

## Examples

``` r
has_exact_deriv4(gaussian1_distrib())
#> [1] TRUE
has_exact_deriv4(truncated(gaussian1_distrib(), lower = 0))
#> [1] TRUE
```
