# Chi-Squared Analytical Fourth-Order Derivative

\$\$\ell^{(\mu\mu\mu\mu)} = -\dfrac{\psi'''(\mu/2)}{16}\$\$ observed and
expected alike.

## Arguments

- distrib:

  A `ChisqDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

- expected:

  Logical; the two coincide, so it changes nothing.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is closed form.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list with the `mu_mu_mu_mu` component.

## See also

[`chisq_distrib`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md)
