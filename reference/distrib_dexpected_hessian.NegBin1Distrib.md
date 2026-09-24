# Derivatives of the Expected Information, Sums Over the Support

The first and second derivatives of the expected information in the
parameters for three families whose expected information is itself a sum
over the support.

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

  A single positive integer, passed to the kernels.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).

## Details

**NB1.** Every component is \\c_0 + c_1 G\\ with \\c_0, c_1\\ rational
in \\(\mu, \theta)\\ and \\G = \mathbb{E}\[\psi'(Y + r) - \psi'(r)\]\\,
\\r = \mu/\theta\\. The derivatives of \\G\\ are sums over the same mass
with the score of the mass beside the summand, the differences
\\\psi^{(n)}(r + y) - \psi^{(n)}(r)\\ taken by their exact recurrences,
from a compiled kernel. Toward the Poisson limit \\\theta \to 0\\ the
composition \\c_0 + c_1 G\\ cancels, which the expected information
itself already does (7.7e-07 relative at \\\theta = 0.005\\ against an
exact sum over the support, 1.6e-05 at \\10^{-3}\\), and each derivative
loses about one further factor of \\1/\theta\\: measured against the
exact sums, 1e-09 at \\\theta = 0.7\\, 1e-07 at 0.2, 1e-05 at 0.05.

**Beta-binomial.** The support \\\\0, \dots, n\\\\ is finite, so the
identities \$\$\partial_c \mathbb{E}\[\ell\_{ab}\] =
\mathbb{E}\[\ell\_{abc} + \ell\_{ab}\ell_c\]\$\$ and its second-order
counterpart, \$\$\partial\_{cd} \mathbb{E}\[\ell\_{ab}\] =
\mathbb{E}\[\ell\_{abcd} + \ell\_{abd}\ell_c + \ell\_{abc}\ell_d +
\ell\_{ab}\ell\_{cd} + \ell\_{ab}\ell_c\ell_d\],\$\$ are exact finite
sums against the mass. They are taken in the shapes by a compiled
kernel, where every derivative of the log-mass is a polygamma difference
at an integer shift and is therefore a sum of negative powers,
accumulated beside the mass with no polygamma called;
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
reads the same sums through its map to the shapes
([`betabinom1_dexpected()`](https://statmodels7.github.io/distributions7/reference/betabinom1_dexpected.md)).
[`support_dexpected()`](https://statmodels7.github.io/distributions7/reference/support_dexpected.md)
computes the same identities from the family's own observed derivatives,
sharing no arithmetic with the kernel, and the two agree to about
\\10^{-12}\\ relative; the kernel is 5 to 15 times faster at \\n =
4000\\.

## See also

[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md),
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
