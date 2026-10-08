# Derivatives of the Expected Information, von Mises

The first and second derivatives of the expected information in the
parameters, from the derivatives of the Bessel ratio \\A(\kappa) =
I_1(\kappa)/I_0(\kappa)\\ and of its inverse, computed by one compiled
kernel per order from numericals7's compiled Bessel ratio.

## Arguments

- distrib:

  A von Mises distribution object.

- y:

  A numeric vector of observations, read for its length.

- theta:

  A named list of parameters.

- scale:

  `"parameter"` or `"link"`.

- approx, nsim:

  Unused.

- ...:

  Unused.

- threads:

  A single positive integer, how many threads the kernel may use.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).

## Details

By the concentration, \\\mathbb{E}\[\ell\_{\mu\mu}\] = -\kappa
A(\kappa)\\ and \\\mathbb{E}\[\ell\_{\kappa\kappa}\] = -A'(\kappa)\\, so
the derivatives in \\\kappa\\ read \\A'\\, \\A''\\ and \\A'''\\. By the
mean resultant length \\\rho = A(\kappa)\\,
\\\mathbb{E}\[\ell\_{\mu\mu}\] = -\rho\\\kappa(\rho)\\ and
\\\mathbb{E}\[\ell\_{\rho\rho}\] = -\kappa'(\rho)\\, so they read the
inverse's derivatives \\\kappa'\\, \\\kappa''\\ and \\\kappa'''\\.
Nothing moves with the direction \\\mu\\. The first-order kernel reads
\\A\\ to \\A''\\ (or \\\kappa\\ to \\\kappa''\\), the second-order one
\\A'\\ to \\A'''\\ (or \\\kappa'\\ to \\\kappa'''\\).

## See also

[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md),
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
