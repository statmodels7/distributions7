# Student's t Mixed Second-Response Derivatives

Closed form. With \\z = (y-\mu)/\sigma\\ and \\Q(z) = (\nu - z^2)/(\nu +
z^2)^2\\, the response Hessian is \\\ell^{(yy)} = -(\nu+1)Q/\sigma^2\\,
and with \\Q' = 2z(z^2 - 3\nu)/(\nu+z^2)^3\\ and \\\partial
Q/\partial\nu = (3z^2-\nu)/(\nu+z^2)^3\\, \$\$\partial\_\mu \ell^{(yy)}
= (\nu+1)Q'/\sigma^3, \quad \partial\_\sigma \ell^{(yy)} = (\nu+1)(2Q +
zQ')/\sigma^3, \quad \partial\_\nu \ell^{(yy)} = -(Q +
(\nu+1)\partial\_\nu Q)/\sigma^2.\$\$

## Arguments

- distrib:

  A `StudentT1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `sigma` and `nu`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.
