# Beta-Binomial Analytical Gradient

The chain rule from the shapes, where every derivative is a difference
of digammas: \$\$\dfrac{\partial \ell}{\partial \alpha} =
\psi(y+\alpha) - \psi(\alpha) - \psi(n+S) + \psi(S), \quad S = \alpha +
\beta\$\$ and likewise in \\\beta\\, with \\\alpha = \mu/\sigma\\ and
\\\beta = (1-\mu)/\sigma\\.

## Arguments

- distrib:

  A `BetaBinom1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `sigma`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

- threads:

  How many threads the kernel may use; below the measured internal
  threshold it stays sequential whatever the count says.

## Value

A named list with the `mu` and `sigma` components.

## See also

[`betabinom1_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
