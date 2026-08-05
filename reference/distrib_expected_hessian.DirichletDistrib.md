# Dirichlet Analytical Expected Hessian

Closed form. \\\mathbb{E}\[\log y_j\] = \psi(\alpha_j) - \psi(\phi)\\,
so \\\mathbb{E}\[g_j\] = -\psi(\phi)\\ for every \\j\\; the columns of
\\A\\ and every second-derivative vector of the simplex sum to zero, so
each data-carrying term drops out and
\$\$\mathbb{E}\[\ell^{(\eta_k\eta_l)}\] = -\phi^2\sum_j t_j
A\_{jk}A\_{jl}, \qquad \mathbb{E}\[\ell^{(\eta_k\phi)}\] = -\phi\sum_j
t_j \mu_j A\_{jk}\$\$

## Arguments

- distrib:

  A `DirichletDistrib` object.

- y:

  A matrix with one row per observation.

- theta:

  A named list of parameters.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is closed form.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list of expected second-derivative components.

## See also

[`dirichlet_distrib`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
