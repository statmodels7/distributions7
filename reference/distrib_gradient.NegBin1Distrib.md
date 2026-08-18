# NB1 Analytical Gradient

The chain rule through \\r = \mu/\theta\\. With \\P = \psi(y+r) -
\psi(r) - \log(1+\theta)\\, \$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{P}{\theta}\$\$ and the derivative in \\\theta\\ adds the terms in
which \\\theta\\ appears outside \\r\\.

## Arguments

- distrib:

  A `NegBin1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `theta`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

- threads:

  How many threads the kernel may use; below the measured internal
  threshold it stays sequential whatever the count says.

## Value

A named list with the `mu` and `theta` components.

## See also

[`negbin1_distrib`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
