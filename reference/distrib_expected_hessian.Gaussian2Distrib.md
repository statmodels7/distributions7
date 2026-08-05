# Gaussian Analytical Expected Hessian in Mean and Variance

\$\$\mathbb{E}\[\ell^{(\mu\mu)}\] = -\dfrac{1}{v}, \qquad
\mathbb{E}\[\ell^{(\mu v)}\] = 0, \qquad \mathbb{E}\[\ell^{(vv)}\] =
-\dfrac{1}{2v^2}\$\$ The two parameters are orthogonal, as they are in
every parametrization of this family.

## Arguments

- distrib:

  A `Gaussian2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `sigma2`.

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

[`gaussian2_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
