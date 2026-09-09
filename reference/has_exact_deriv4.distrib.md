# Whether a Family's Fourth Derivative Is Analytic

The default: a
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
method owned by one of the base classes is the package's own numerical
fallback, and anything else is the family's own method.

## Arguments

- x:

  A distribution object.

- ...:

  Unused.

## Value

A single logical.

## Details

The comparison is by class name and package rather than with
[`identical()`](https://rdrr.io/r/base/identical.html), for the reason
recorded on
[`is_class()`](https://statmodels7.github.io/distributions7/reference/is_class.md):
[`identical()`](https://rdrr.io/r/base/identical.html) on an S7 class is
object identity and is false for a class re-created from the same
definition, as happens whenever a package's code is evaluated instead of
loaded.

A wrapper must not use this reading. Its fourth derivative is registered
on its own class, so the owner test says analytic, while the partition
sum it evaluates is over the **parent's** first four derivatives and is
analytic only when those are. Every wrapper therefore overrides the
generic and asks its parent, which is why the test is a generic at all.

## See also

[`has_exact_deriv4()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv4.md),
the generic;
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md),
which answers the same shape of question about the expected information.
