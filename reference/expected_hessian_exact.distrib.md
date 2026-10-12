# Expected Information Read from the Class of Its Method

The default: a method owned by one of the base classes is one of the
package's approximations, and anything else is taken to have written its
expected information out. It answers `FALSE` for a family that registers
nothing and `TRUE` for one that registers its own method.

A family whose registered method chains onto a parent overrides this
generic and returns the answer of the parent, which is why the test is a
generic at all; the wrappers do so.

## Arguments

- x:

  A distribution object.

- ...:

  Unused.

## Value

A single logical: `TRUE` where the family computes its expected
information itself, in closed form or by a quadrature of its own,
`FALSE` where it reaches one of the approximations of the package.

## See also

[`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md),
the predicate consumers call;
[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for what the approximations are;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
which rejects an `approx` given to
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
where this is `TRUE`.
