# Generalised Pareto Analytical Gradient

With \\z = y/\sigma\\, \\t = 1 + \xi z\\ and \\u = z/t\\,
\$\$\dfrac{\partial \ell}{\partial \sigma} = \dfrac{(\xi+1)u -
1}{\sigma}, \qquad \dfrac{\partial \ell}{\partial \xi} = \dfrac{\log
t}{\xi^2} - \left(1 + \dfrac{1}{\xi}\right)u\$\$ the second computed
through a series near \\\xi = 0\\, where its limit is \\z^2/2 - z\\.

## Arguments

- distrib:

  A `GPDDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `sigma` and `xi`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list with the `sigma` and `xi` components.

## See also

[`gpd_distrib`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
