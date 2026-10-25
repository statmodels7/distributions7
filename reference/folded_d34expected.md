# The Third and Fourth Derivatives of a Folded Family's Expected Information

Computes the derivatives of order `k` (3 or 4) of the expected Hessian
of a
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
family on the parameter scale, as integrals over \\y \> 0\\ against the
folded density, on the nodes of `fold_rule_cpp()` that
[`folded_expected()`](https://statmodels7.github.io/distributions7/reference/folded_expected.md)
uses.

## Usage

``` r
folded_d34expected(d, y, th, k, t = 1L)
```

## Arguments

- d:

  A folded family, from
  [`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md).

- y:

  The response, read for its length.

- th:

  The parameters.

- k:

  The order, `3L` or `4L`.

- t:

  Unused; it is the thread count of the registration.

## Value

A list of components keyed as
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
or
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md),
each of length `length(y)`.

## Details

With \\\ell\\ the folded log-density and \\E\_{ab} =
-E\[\ell_a\ell_b\]\\, the component for a multiset \\S\\ is minus the
integral of the polynomial of
[`d34e_plan()`](https://statmodels7.github.io/distributions7/reference/d34e_plan.md)
against the folded density. The polynomial reads the folded derivatives
of \\\ell\\ to order \\k + 1\\: those of orders one to four are the
family's own methods, and those of order five are formed by
`fold_deriv_k(5L)` from the parent's fifth derivatives, which the fourth
derivative of the information needs. The form \\E\[\ell\_{ab}\]\\ would
need the sixth derivative of the parent.

Where the parent's mass beyond the fold is below the square of the
machine epsilon, the components are the parent's own derivatives of the
expected information, as in
[`folded_expected()`](https://statmodels7.github.io/distributions7/reference/folded_expected.md).

A parent with a kink at its center (a parameter that is not smooth)
moves the point where the folded score jumps, and the integral then
gains boundary terms that involve the mixed derivatives of \\\ell\\ in
the response and the parameters beyond the first (the package has
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
only for one response derivative and one parameter). These parents, and
a parent with no scalar route, are served by
[`numerical_dexpected_higher()`](https://statmodels7.github.io/distributions7/reference/numerical_dexpected_higher.md),
one stencil on the second derivative.

## See also

[`folded_expected()`](https://statmodels7.github.io/distributions7/reference/folded_expected.md),
[`trunc_d34expected()`](https://statmodels7.github.io/distributions7/reference/trunc_d34expected.md)
