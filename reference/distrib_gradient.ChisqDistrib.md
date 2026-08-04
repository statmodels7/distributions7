# Chi-Squared Analytical Gradient

\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{\log y - \log 2 -
\psi(\mu/2)}{2}\$\$ the only order that involves the response.

## Arguments

- distrib:

  A `ChisqDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list with the `mu` component.

## See also

[`chisq_distrib`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md)
