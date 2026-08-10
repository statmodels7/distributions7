# Dirichlet Analytical Third and Fourth Derivatives

Closed form. The log-density is \\\log\Gamma(\phi) - \sum_j
\log\Gamma(\alpha_j) + \sum_j(\alpha_j - 1)\log y_j\\ with \\\alpha =
\phi\mu(\eta)\\, so every term depends on one coordinate of the simplex
and the chain rule is one univariate partition sum per coordinate. The
concentration also enters directly, through \\\log\Gamma(\phi)\\, which
contributes only to the component all of whose indices name it.

## Arguments

- distrib:

  A `DirichletDistrib` object.

- y:

  A matrix with one row per observation.

- theta:

  A named list of parameters.

- expected:

  Whether to return expected derivatives.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  The approximation used when `expected` is `TRUE`.

- nsim:

  The Monte Carlo sample size.

- ...:

  Unused.

## Value

A named list of third-derivative components.

## See also

[`dirichlet_distrib`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
