# Default Expected Hessian for `distrib` Objects

Fallback method: the expected Hessian is computed as the negative Fisher
information via the outer product of the score (gradient),
\\-\mathbb{E}\[\nabla\ell\\\nabla\ell^\top\]\\, obtained by numerical
[`expectation`](https://statmodels7.github.io/distributions7/reference/expectation.md).
This first-Bartlett form equals \\\mathbb{E}\[H\]\\ for regular
(twice-differentiable) models but, unlike \\\mathbb{E}\[H\]\\, remains
valid when the log-likelihood is non-differentiable in a parameter (e.g.
the location of a Laplace distribution), where the observed Hessian is
degenerate. The score is taken from
[`distrib_gradient`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
so it uses the analytical gradient when available and finite differences
otherwise.

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

## Value

A named list of expected Hessian component vectors.
