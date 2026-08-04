# Gumbel Analytical Observed Hessian

Closed-form second derivatives of the Gumbel log-density:
\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{w}{\sigma^2},
\qquad \dfrac{\partial^2 \ell}{\partial \mu \\ \partial \sigma} =
-\dfrac{1 - w + z w}{\sigma^2},\$\$ \$\$\dfrac{\partial^2 \ell}{\partial
\sigma^2} = \dfrac{1 - 2z + 2zw - z^2 w}{\sigma^2}.\$\$

## Arguments

- distrib:

  A `GumbelDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- scale:

  Either `"parameter"` or `"link"`.

- ...:

  Unused.

## Value

A named list of second derivatives.

## See also

[`gumbel_distrib`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
