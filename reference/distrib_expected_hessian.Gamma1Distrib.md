# Gamma Analytical Expected Hessian in Mean and Dispersion

Closed form: \$\$\mathbb{E}\[\ell^{(\mu\mu)}\] = -\dfrac{1}{\phi\mu^2},
\qquad \mathbb{E}\[\ell^{(\mu\phi)}\] = 0, \qquad
\mathbb{E}\[\ell^{(\phi\phi)}\] = s^4\left\\\dfrac{1}{s} -
\psi'(s)\right\\\$\$ The mean and the dispersion are orthogonal, which
is what makes this the natural parametrization for a generalized linear
model. The expectation uses \\\mathbb{E}\[\log(Y/\mu)\] = \psi(s) - \log
s\\.

## Arguments

- distrib:

  A `Gamma1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `phi`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is closed form.

- nsim:

  Ignored.

- ...:

  Unused.

- threads:

  How many threads the kernel may use; below the measured internal
  threshold it stays sequential whatever the count says.

## Value

A named list of expected second derivatives.

## See also

[`gamma1_distrib`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
