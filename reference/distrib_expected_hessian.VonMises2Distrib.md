# von Mises Analytical Expected Hessian in the Resultant Length

Closed form. Since \\\mathbb{E}\[\cos(Y-\mu)\] = A\\ and
\\\mathbb{E}\[\sin(Y-\mu)\] = 0\\, the term in the second derivative of
the map drops and \$\$\mathbb{E}\[\ell^{(\mu\mu)}\] = -\kappa A, \qquad
\mathbb{E}\[\ell^{(\mu\rho)}\] = 0, \qquad
\mathbb{E}\[\ell^{(\rho\rho)}\] = -\dfrac{1}{A'(\kappa)}.\$\$ The last
is the inverse of the information in \\\kappa\\, which is what a
one-to-one reparametrization of a single parameter must give, and the
two parameters stay orthogonal.

## Arguments

- distrib:

  A `VonMises2Distrib` object.

- y:

  A numeric vector of angles.

- theta:

  A list with `mu` and `rho`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is closed form.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list of expected second derivatives.

## See also

[`vonmises2_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md)
