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
\\\nu^{-2}\\; it and its derivatives in \\\nu\\ are taken from their
asymptotic series in \\1/\nu\\ above \\\nu = 20\\, derived from
Stirling's series, and the first and second derivatives come from one
kernel each.

For the generalized gamma by scale \\a\\ and shapes \\d, p\\, each
component is a power of \\1/a\\ times a closed form in \\(d, p)\\ and
the polygamma functions at \\k + 1\\, \\k = d/p\\, derived offline as
[`distrib_deriv3.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GenGamma1Distrib.md)
describes.

## See also

[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md),
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
