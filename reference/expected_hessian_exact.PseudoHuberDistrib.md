# The Pseudo-Huber Does Not Write Its Expected Information Out

Answers `FALSE`, declaring that this family's expected information is a
numerical approximation and not a formula, so that callers who branch on
the distinction branch correctly.

## Arguments

- x:

  A `PseudoHuberDistrib` object, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

`FALSE`, a logical of length 1.

## Details

The predicate's default reads the class a method is registered on, which
here would answer `TRUE` and be wrong.
[`distrib_expected_hessian.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.PseudoHuberDistrib.md)
is registered on this class, but what it does is call
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
and then replace the two components that vanish by symmetry: it
**improves** the approximation rather than replacing it.

The cost is the discriminator. Measured at 100 observations the method
takes about 11 seconds, where the families that do write their
information out answer in a median of 0.183 milliseconds. Two
consequences were live before the declaration:
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
rejected a legitimate `fisher_scoring(approx = )` here with a message
saying the family computes its expected information in closed form,
which is untrue; and its standard-error branch entered a multi-second
quadrature believing it a formula.

## See also

[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
for the generic and its default,
[`distrib_expected_hessian.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.PseudoHuberDistrib.md)
for the method this describes, and
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
for the consumer.
