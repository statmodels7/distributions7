# Expected Information of a Multivariate Distribution

Fallback for a multivariate family with no closed form: the expectation
is taken over draws from the distribution itself.

## Arguments

- distrib:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- y:

  An \\n \times p\\ matrix of observations; only its row count is used,
  the expectation being over the distribution rather than the data.

- theta:

  A named list of parameters.

- scale:

  Handled by the generic before dispatch.

- approx:

  One of `"bartlett"` (equivalently `"opg"`) or `"mc"`; `"integrate"` is
  rejected.

- nsim:

  Monte Carlo sample size.

- ...:

  Unused.

## Value

A named list keyed as
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md)`(distrib@params)`.

## Details

The one-dimensional routes do not survive the move to \\p\\ dimensions.
`"integrate"` builds its quadrature over an interval and is rejected
here; `"bartlett"` in the univariate package reaches
[`expectation`](https://statmodels7.github.io/distributions7/reference/expectation.md),
which is that same quadrature. What does generalize is sampling, so both
remaining routes draw from the family's own generator and differ in what
they average:

`"mc"` averages the observed Hessian, \\\mathbb{E}\[\ell^{(ij)}\]\\
directly. `"bartlett"` and `"opg"` average the outer product of the
score and negate it, which is the second Bartlett identity \\\mathcal{I}
= \mathbb{E}\[s s^\top\]\\; it needs no second derivative at all, and is
the only route that survives a family whose observed Hessian is
degenerate.

Both are Monte Carlo, so both carry an error of order
\\1/\sqrt{\texttt{nsim}}\\, and a fit that uses one is doing Fisher
scoring with a noisy information. That is a deliberate choice a caller
makes, which is why
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
rejects the argument for a family that has an exact expression.
