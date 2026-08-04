# Chi-Squared Analytical Observed Hessian

\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} =
-\dfrac{\psi'(\mu/2)}{4}\$\$ which does not involve the response, so it
coincides with the expected information.

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

A named list with the `mu_mu` component.

## See also

[`chisq_distrib`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md)
