# Poisson-Inverse Gaussian Probability Mass Function

With \\c = 1 + 2\sigma\mu\\ and \\\alpha = \sqrt{c}/\sigma\\, \$\$P(Y =
y) = \sqrt{\dfrac{2\alpha}{\pi}}\\ \dfrac{\mu^y
e^{1/\sigma}}{(\alpha\sigma)^y\\ y!}\\ K\_{y-1/2}(\alpha),\$\$ evaluated
through the finite half-integer Bessel sum, in which the prefactors
cancel down to \\\ell(y) = y\log\mu - (y/2)\log c +
(1-\sqrt{c})/\sigma + \log S_y(\alpha) - \log y!\\. A non-integer or
negative \\y\\ has probability zero.

## Arguments

- distrib:

  A `Pig1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `sigma`.

- log:

  Logical; if `TRUE`, returns the log-probability.

## Value

A numeric vector of probability values.

## See also

[`pig1_distrib`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
