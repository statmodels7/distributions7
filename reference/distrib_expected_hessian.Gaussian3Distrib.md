# Gaussian Analytical Expected Hessian in Mean and Precision

\$\$\mathbb{E}\[\ell^{(\mu\mu)}\] = -\tau, \qquad
\mathbb{E}\[\ell^{(\mu\tau)}\] = 0, \qquad
\mathbb{E}\[\ell^{(\tau\tau)}\] = -\dfrac{1}{2\tau^2}\$\$ Only the mixed
entry differs from the observed Hessian, which is what makes Fisher
scoring and Newton's method take the same step on the parameter scale
here.

## Arguments

- distrib:

  A `Gaussian3Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `tau`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is closed form.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list of expected second derivatives.

## See also

[`gaussian3_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md)
