# Whether the Owner of the Method Settles the Question

The default: a method owned by one of the base classes is one of the
package's approximations, and anything else is taken to have written its
expected information out. It answers `FALSE` for a family that registers
nothing and `TRUE` for one that registers its own method.

The reading is right for 34 of the 40 univariate families and wrong for
two, which declare a method of their own for themselves. The
pseudo-Huber's method calls
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
and patches the two components that vanish by symmetry, and
`skewnormal2`'s chains onto `skewnormal1`, whose expected information is
the base class's quadrature. Both therefore **override this generic** in
place of relying on the owner test, which is why the test is a generic
at all.

## Arguments

- x:

  A distribution object.

- ...:

  Unused.

## Value

A single logical: `TRUE` where the family's expected information is a
closed-form expression, `FALSE` where it is a quadrature or a
simulation.

## See also

[`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md),
the predicate consumers call;
[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for what the approximations are;
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which rejects an `approx` where this is `TRUE`.
