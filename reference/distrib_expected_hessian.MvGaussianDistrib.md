# Multivariate Gaussian Expected Hessian

Closed form, and simpler than the observed one.
\$\$\mathbb{E}\[\ell^{(\mu_a \mu_b)}\] = -(\Sigma^{-1})\_{ab}, \qquad
\mathbb{E}\[\ell^{(\mu_a \eta_k)}\] = 0, \qquad
\mathbb{E}\[\ell^{(\eta_k \eta_l)}\] =
-\tfrac{1}{2}\mathrm{tr}(\Sigma^{-1} A_k \Sigma^{-1} A_l).\$\$

## Arguments

- distrib:

  A
  [`MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object.

- y:

  An \\n \times p\\ matrix of observations.

- theta:

  A named list of parameters.

- scale:

  Handled by the generic; the two scales coincide here.

- approx:

  Ignored: the expectation is exact.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list keyed as
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md)`(distrib@params)`.

## Details

The mixed block vanishes because \\\mathbb{E}\[w\] = 0\\, which is the
orthogonality of the mean and the covariance parameters that makes
Fisher scoring on this family so well behaved. The matrix block needs no
second derivative at all: the terms in \\A\_{kl}\\ cancel between the
log-determinant and the quadratic form.
