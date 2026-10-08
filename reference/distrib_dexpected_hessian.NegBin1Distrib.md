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

**NB1.** The expected information is \\-A\\, \\W/\theta^2\\ and
\\-W/(\mu\theta)\\, with \\A\\ and \\W\\ the two sums over the support
of
[`distrib_expected_hessian.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.NegBin1Distrib.md).
Their derivatives are sums over the same mass with the score of the mass
beside the summand, \\\partial_c F = \mathbb{E}\[f_c + f s_c\]\\ and
\\\partial\_{ce} F = \mathbb{E}\[f\_{ce} + f_c s_e + f_e s_c + f(s_c
s_e + h\_{ce})\]\\, taken in \\(\mu, \theta)\\ by a compiled kernel,
with the scores \\s\\ and second derivatives \\h\\ of the log-mass
written so that no term of order \\\mu/\theta^2\\ is formed. Measured
against exact sums at 50 digits, the error relative to the largest
component of each order is about 2e-10 at \\\mu = 100\\, \\\theta =
0.05\\, and grows toward the Poisson limit, to 4e-08 at \\\theta =
0.005\\ and 5e-06 at \\\theta = 0.001\\ (\\\mu = 10\\), the fourth
derivative in \\\theta\\ being the worst component.

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
