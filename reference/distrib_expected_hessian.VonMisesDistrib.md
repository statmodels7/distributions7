# von Mises Analytical Expected Hessian

\\\mathbb{E}\[\cos(Y-\mu)\] = A(\kappa)\\ and
\\\mathbb{E}\[\sin(Y-\mu)\] = 0\\ by symmetry, so
\$\$\mathbb{E}\[\ell^{(\mu\mu)}\] = -\kappa A(\kappa), \qquad
\mathbb{E}\[\ell^{(\mu\kappa)}\] = 0, \qquad
\mathbb{E}\[\ell^{(\kappa\kappa)}\] = -A'(\kappa)\$\$ The location and
the concentration are therefore orthogonal.

## Arguments

- distrib:

  A `VonMisesDistrib` object.

- y:

  A numeric vector of angles.

- theta:

  A list containing `mu` and `kappa`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is closed form.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list of expected second-derivative components.

## See also

[`vonmises_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises_distrib.md)
