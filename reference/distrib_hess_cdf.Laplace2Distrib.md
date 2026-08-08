# Laplace Log-CDF Second Derivatives in Location and Rate

Closed form, exact away from \\q = \mu\\: with \\s =
\mathrm{sign}(q-\mu)\\ and \\a = \lvert q-\mu \rvert\\, \\\partial^{2}
F/\partial\mu^{2} = -\lambda s f\\, \\\partial^{2}
F/\partial\mu\\\partial\lambda = -(1/\lambda - a) f\\ and \\\partial^{2}
F/\partial\lambda^{2} = -s a^{2} f/\lambda\\. At the kink the second
derivative in \\\mu\\ does not exist, and the value returned is the
one-sided limit the sign convention picks out.

## Arguments

- distrib:

  A `Laplace2Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `lambda`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list keyed as
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## See also

[`laplace2_distrib`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
