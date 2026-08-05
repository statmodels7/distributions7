# Dirichlet Analytical Observed Hessian

The second derivatives of the same expressions, with \\t_j =
\psi'(\alpha_j)\\: \$\$\ell^{(\eta_k\eta_l)} = \phi\sum_j\left(-t_j\phi
A\_{jk}A\_{jl} + g_j B\_{j,kl}\right), \qquad \ell^{(\phi\phi)} =
\psi'(\phi) - \sum_j t_j \mu_j^2\$\$ the last free of the data, the
family being an exponential family in \\\log y\\.

## Arguments

- distrib:

  A `DirichletDistrib` object.

- y:

  A matrix with one row per observation.

- theta:

  A named list of parameters.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of second-derivative components.

## See also

[`dirichlet_distrib`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
