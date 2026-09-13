# Beta Derivatives of the Expected Information in Mean and Precision

The first and second derivatives of the expected information in \\(\mu,
\phi)\\, from a compiled kernel. With \\\alpha = \mu\phi\\, \\\beta =
(1-\mu)\phi\\ and polygamma functions evaluated there and at \\\phi\\,
\$\$\mathbb{E}\[\ell\_{\mu\mu}\] = -\phi^2(\psi'(\alpha) +
\psi'(\beta)),\quad \mathbb{E}\[\ell\_{\mu\phi}\] =
-\phi(\mu\psi'(\alpha) - (1-\mu)\psi'(\beta)),\$\$
\$\$\mathbb{E}\[\ell\_{\phi\phi}\] = \psi'(\phi) - \mu^2\psi'(\alpha) -
(1-\mu)^2\psi'(\beta),\$\$ and every component is the ordinary
derivative of one of the three through \\\partial\alpha/\partial\mu =
\phi\\, \\\partial\beta/\partial\mu = -\phi\\,
\\\partial\alpha/\partial\phi = \mu\\, \\\partial\beta/\partial\phi =
1-\mu\\. On the link scale the result is carried across by
[`dexpected_link()`](https://statmodels7.github.io/distributions7/reference/dexpected_link.md).

## Arguments

- distrib:

  A `Beta1Distrib` object.

- y:

  A numeric vector of observations, read for its length.

- theta:

  A named list with `mu` and `phi`.

- scale:

  `"parameter"` or `"link"`.

- approx, nsim:

  Unused.

- ...:

  Unused.

- threads:

  A single positive integer, the kernel's thread count.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).

## See also

[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md),
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
