# Negative Binomial Derivatives of the Expected Information, NB2

The first and second derivatives of the expected information in \\(\mu,
\theta)\\, from a compiled kernel. \\\mathbb{E}\[\ell\_{\mu\mu}\] =
-\theta/(\mu(\theta+\mu))\\ is differentiated in closed form and
\\\mathbb{E}\[\ell\_{\mu\theta}\] = 0\\.
\\\mathbb{E}\[\ell\_{\theta\theta}\] = S = \sum_k p_k U_k\\ is a sum
over the support whose mass depends on the parameters, so a derivative
moves the mass as well as the summand: \$\$\partial_x S = \sum_k p_k
(U\_{k,x} + U_k s_x(k)),\$\$ \$\$\partial\_{xy} S = \sum_k p_k
\big(U\_{k,xy} + U\_{k,x} s_y + U\_{k,y} s_x + U_k (s\_{xy} + s_x
s_y)\big),\$\$ with \\s\\ the score of the log-mass at \\k\\. Every
piece is an exact polynomial in \\k\\ or a sum
\\\sum\_{j\<k}(\theta+j)^{-r}\\ accumulated beside the mass, over the
recurrence and stopping rule the expected information itself uses. On
the link scale the result is carried across by
[`dexpected_link()`](https://statmodels7.github.io/distributions7/reference/dexpected_link.md).

## Arguments

- distrib:

  A `NegBin2Distrib` object.

- y:

  A numeric vector of observations, read for its length.

- theta:

  A named list with `mu` and `theta`.

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
