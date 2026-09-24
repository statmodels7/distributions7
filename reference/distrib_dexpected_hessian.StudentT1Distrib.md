# Derivatives of the Expected Information, Student t and Generalized Gamma

The first and second derivatives of the expected information in the
parameters, from compiled kernels.

## Arguments

- distrib:

  A distribution object of one of the classes above.

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

  A single positive integer, the kernel's thread count.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).

## Details

For the Student t, with \\w_1 = 1 + 1/\nu\\ and \\w_3 = 1 + 3/\nu\\,
\\\mathbb{E}\[\ell\_{\mu\mu}\] = -(w_1/w_3)/\sigma^2\\,
\\\mathbb{E}\[\ell\_{\sigma\sigma}\] = -(2/w_3)/\sigma^2\\,
\\\mathbb{E}\[\ell\_{\sigma\nu}\] = 2\nu^{-2}/(w_1 w_3 \sigma)\\ and
\$\$\mathbb{E}\[\ell\_{\nu\nu}\] = \tfrac14\\\psi'((\nu+1)/2) -
\psi'(\nu/2)\\ + \frac{\nu+5}{2\nu(\nu+1)(\nu+3)} = -\frac{7}{2\nu^4} +
\frac{13}{\nu^5} - \dots,\$\$ whose terms cancel from order
\\\nu^{-2}\\; its derivatives in \\\nu\\ are taken from the asymptotic
series above \\\nu = 30\\, the series' coefficients being exact integers
from the duplication identity \\\psi'((\nu+1)/2) - \psi'(\nu/2) =
4\psi'(\nu) - 2\psi'(\nu/2)\\.

For the generalized gamma by scale \\a\\ and shapes \\d, p\\, each
component is \\a^{\alpha} p^{\beta} F(k)\\ with \\k = d/p\\ and \\F\\ a
combination of \\\psi(k)\\, \\\psi(k+1)\\ and their derivatives,
differentiated through \\\partial k/\partial d = 1/p\\ and \\\partial
k/\partial p = -k/p\\.

## See also

[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md),
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
