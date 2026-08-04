# Chi-Squared Analytical Expected Hessian

\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{\psi'(\mu/2)}{4}\$\$ identical to the observed Hessian on the
parameter scale, which the second derivative not depending on the
response makes exact rather than approximate. The two differ on the link
scale, where the chain rule adds a term proportional to the score; see
[`chisq_distrib`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

## Arguments

- distrib:

  A `ChisqDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is closed form.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list with the `mu_mu` component.

## See also

[`chisq_distrib`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md)
