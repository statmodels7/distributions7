# Gamma Derivatives of the Expected Information in Mean and Dispersion

The first and second derivatives of the expected information in \\(\mu,
\phi)\\, from a compiled kernel. With \\s = 1/\phi\\,
\\\mathbb{E}\[\ell\_{\mu\mu}\] = -s/\mu^2\\,
\\\mathbb{E}\[\ell\_{\mu\phi}\] = 0\\ and
\\\mathbb{E}\[\ell\_{\phi\phi}\] = q(s) = s^4(1/s - \psi'(s))\\; writing
\\f_2 = 1/s - \psi'(s)\\ and \\f_3\\, \\f_4\\ its derivatives,
\$\$\partial\_\phi \mathbb{E}\[\ell\_{\phi\phi}\] = -s^2 q'(s),\qquad
\partial\_{\phi\phi} \mathbb{E}\[\ell\_{\phi\phi}\] = s^4 q''(s) + 2 s^3
q'(s),\$\$ with \\q' = f_3 s^4 + 4 f_2 s^3\\ and \\q'' = f_4 s^4 + 8 f_3
s^3 + 12 f_2 s^2\\, each \\f_k\\ read in the form that forms no
cancellation at large \\s\\. The derivatives of \\-\mu^{-2}\phi^{-1}\\
are monomials. On the link scale the result is carried across by
[`dexpected_link()`](https://statmodels7.github.io/distributions7/reference/dexpected_link.md).

## Arguments

- distrib:

  A `Gamma1Distrib` object.

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
