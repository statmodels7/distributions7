# Default Expected Hessian for `distrib` Objects

This method is the fallback for a family that does not write its
expected information out. It rests on the second Bartlett identity,
\\\mathbb{E}\[\ell^{(ij)}\] = -\mathbb{E}\[\ell^{(i)}\ell^{(j)}\]\\,
which holds for a regular model and, unlike
\\\mathbb{E}\[\ell^{(ij)}\]\\ read directly, survives a log-likelihood
that is not differentiable in a parameter (the location of a Laplace,
where the observed Hessian is degenerate while the score variance is
still the information).

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  `"parameter"` or `"link"`, the scale on which the components are
  reported. The transformation is applied in the generic's body, so this
  method always returns the parameter scale.

- method:

  `NULL` (the default), which selects
  [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md),
  or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)
  that
  [`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
  applies. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of expected Hessian component vectors.

## Details

The argument `method` selects how the right-hand side is obtained, and
the methods differ in cost. With `method = NULL` the default is
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md),
which reads \\-\ell^{(i)}\ell^{(j)}\\ at each observation and takes no
expectation, so it costs one call to
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md);
`expect_integrate(identity = "score")` evaluates the expectation itself,
which is a sum over the support for a discrete family and a quadrature
for a continuous one, and is orders of magnitude dearer. See
[`expected_by_opg()`](https://statmodels7.github.io/distributions7/reference/expected_by_opg.md)
for the properties that the default loses and those that it retains.

The score is taken from
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
so it uses the analytical gradient where the family has one and finite
differences otherwise.

## See also

[`expected_by_opg()`](https://statmodels7.github.io/distributions7/reference/expected_by_opg.md)
and
[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md)
for the two readings of the identity, and
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
for the predicate that determines whether a family reaches this method
at all.
