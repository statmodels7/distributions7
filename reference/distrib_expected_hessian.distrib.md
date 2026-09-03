# Default Expected Hessian for `distrib` Objects

Fallback method, for a family that does not write its expected
information out. It rests on the second Bartlett identity,
\\\mathbb{E}\[\ell^{(ij)}\] = -\mathbb{E}\[\ell^{(i)}\ell^{(j)}\]\\,
which holds for a regular model and, unlike
\\\mathbb{E}\[\ell^{(ij)}\]\\ read directly, survives a log-likelihood
that is not differentiable in a parameter – the location of a Laplace,
where the observed Hessian is degenerate while the score variance is
still the information.

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  `"parameter"` or `"link"`, the scale the components are reported on.
  The transformation is applied in the generic's body, so this method
  always returns the parameter scale.

- approx:

  One of `"opg"`, `"bartlett"`, `"integrate"` or `"mc"`, the strategy
  [`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
  uses. Defaults to `"opg"`.

- nsim:

  Number of draws, read only by `approx = "mc"`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of expected Hessian component vectors.

## Details

`approx` says how the right-hand side is obtained, and the choice is a
choice of cost. The default `"opg"` reads \\-\ell^{(i)}\ell^{(j)}\\ at
each observation and takes no expectation, so it costs one call to
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md);
`"bartlett"` evaluates the expectation itself, which is a sum over the
support for a discrete family and a quadrature for a continuous one, and
is orders of magnitude dearer. See
[`expected_by_opg()`](https://statmodels7.github.io/distributions7/reference/expected_by_opg.md)
for what the default gives up and what it does not.

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
for the predicate that says whether a family reaches this method at all.
