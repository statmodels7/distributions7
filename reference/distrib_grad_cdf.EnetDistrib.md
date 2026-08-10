# Elastic-Net Log-CDF Derivatives

Closed form at every order. Each half of the distribution function is a
truncated Gaussian, so with \\z = q - \mu\\, \\s = \sqrt{c}\\ and \\x =
a/\sqrt{c}\\ it is \\e^{w}\Phi(X)\\ below the location and \\1 -
e^{w}\Phi(X)\\ above it, for \\X = \pm sz - x\\ and a weight \\w = -\log
M(x) + x^{2}/2 + \mathrm{const}\\ written through the Mills ratio the
family already carries. Both \\s\\ and \\x\\ are products of a function
of \\\lambda\\ and a function of \\\alpha\\, so their mixed partial
derivatives are products of one-variable ones.

The location is the non-regular direction, as in the Laplace the family
contains: the second derivative in \\\mu\\ carries a point mass at \\q =
\mu\\, and the formulas below hold on either side of it.

## Arguments

- distrib:

  An `EnetDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu`, `lambda` and `alpha`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

- ...:

  Unused.

## Value

A named list, one vector per component.

## See also

[`enet_distrib`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
